import pysam
import numpy as np
import re, os, sys
import pandas as pd
from multiprocessing import Pool
import functools
from scipy import stats
import pysam
import tqdm
import logging
import warnings
from bam_filter.utils import is_debug, calc_chunksize, fast_flatten, initializer
import pyranges as pr
from collections import defaultdict

log = logging.getLogger("my_logger")

sys.setrecursionlimit(10 ** 6)

# Function to calculate evenness of coverage
def coverage_evenness(coverage):
    """
    Calculate the evenness of coverage
    """
    # get coverage evenness
    # covEvenness = (
    #     (breadth * 100) / (100.0 * expBreadth * expBreadth) if expBreadth > 0 else 0
    # )
    # C = mean(X)
    # D2 = X[X<=C]
    # N = len(X)
    # n = len(D2)
    # E = 1 - (n - sum(D2) / C) / N
    C = float(round(np.mean(coverage)))
    D2 = [x for x in coverage if x <= C]
    if len(D2) == 0:  # pragma: no cover
        covEvenness = 1.0
    else:
        covEvenness = 1.0 - (len(D2) - sum(D2) / C) / len(coverage)

    return covEvenness


# function to calculate GC content
def calc_gc_content(seq):
    """Calculate GC content of a sequence

    Args:
        seq (str): DNA sequence

    Returns:
        int: Number of GC content
    """
    gc = seq.count("G") + seq.count("C")
    return gc


def create_pyranges(reference, starts, ends, strands):
    """[summary]

    Args:
        reference ([type]): [description]
        starts ([type]): [description]
        ends ([type]): [description]
        strands ([type]): [description]
    """
    chromosomes = [reference] * len(starts)
    chromosomes = pd.Series(chromosomes).astype("category")
    starts = pd.Series(starts)
    ends = pd.Series(ends)
    strands = pd.Series(strands).astype("category")

    return pr.PyRanges(
        pd.DataFrame(
            {"Chromosome": chromosomes, "Start": starts, "End": ends, "Strand": strands}
        )
    )


def get_bam_stats(params, ref_lengths=None, scale=1e6):
    """
    Worker function per chromosome
    loop over a bam file and create tuple with lists containing metrics:
    two definitions of the edit distances to the reference genome scaled by aligned read length
    """
    bam, reference = params
    samfile = pysam.AlignmentFile(bam, "rb")
    edit_distances = []
    # edit_distances_md = []
    ani_nm = []
    # ani_md = []
    read_length = []
    read_aligned_length = []
    read_mapq = []
    read_aln_score = []
    read_names = []
    read_gc_content = []
    n_alns = 0
    if ref_lengths is None:
        reference_length = int(samfile.get_reference_length(reference))
        bam_reference_length = reference_length
    else:
        reference_length = int(ref_lengths.loc[reference, "length"])
        bam_reference_length = int(samfile.get_reference_length(reference))

    log.debug(f"Processing reference {reference}")
    log.debug(f"Reference length: {reference_length:,}")
    log.debug(f"BAM reference length: {bam_reference_length:,}")
    starts = []
    ends = []
    strands = []
    for aln in samfile.fetch(reference=reference, multiple_iterators=False):
        n_alns += 1
        if aln.has_tag("AS"):
            read_aln_score.append(aln.get_tag("AS"))
        else:
            read_aln_score.append(np.nan)

        if aln.has_tag("AS"):
            edit_distances.append(aln.get_tag("NM"))
            ani_nm.append((1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100)
        else:
            edit_distances.append(np.nan)
            ani_nm.append(np.nan)

        read_gc_content.append(calc_gc_content(aln.query_sequence))
        read_length.append(aln.infer_read_length())
        read_aligned_length.append(aln.query_alignment_length)
        read_mapq.append(aln.mapping_quality)
        read_names.append(aln.query_name)
        # check if strand is reverse
        if aln.is_reverse:
            strand = "-"
        else:
            strand = "+"
        starts.append(aln.reference_start)
        ends.append(aln.reference_end)
        strands.append(strand)

    # get bases covered by reads pileup
    cov_pos = [
        pileupcolumn.n
        for pileupcolumn in samfile.pileup(
            reference,
            start=None,
            stop=None,
            region=None,
            stepper="nofilter",
            min_mapping_quality=0,
        )
    ]
    samfile.close()
    # convert datafrane to pyranges
    ranges = create_pyranges(reference, starts, ends, strands)
    ranges = ranges.merge(strand=False).lengths().to_list()
    max_covered_bases = np.max(ranges)
    mean_covered_bases = np.mean(ranges)
    bases_covered = int(len(cov_pos))
    # get SD from covered bases
    cov_sd = np.std(cov_pos, ddof=1)
    # get average coverage
    mean_coverage = sum(cov_pos) / reference_length
    mean_coverage_covered = sum(cov_pos) / bases_covered

    breadth = bases_covered / reference_length
    exp_breadth = 1 - np.exp(-mean_coverage)
    breadth_exp_ratio = breadth / exp_breadth

    if breadth_exp_ratio > 1:
        breadth_exp_ratio = 1.0

    cov_evenness = coverage_evenness(cov_pos)
    gc_content = (np.sum(read_gc_content) / np.sum(read_length)) * 100
    c_v = cov_sd / mean_coverage
    read_mapq = [np.nan if x == 255 else x for x in read_mapq]

    tax_abund_aln = round((n_alns / reference_length) * scale)
    tax_abund_read = round((len(set(read_names)) / reference_length) * scale)

    log.debug(f"Number of reads: {len(set(read_names)):,}")
    log.debug(f"Number of alignments: {n_alns:,}")
    log.debug(f"Bases covered: {bases_covered:,}")
    log.debug(f"Mean coverage: {mean_coverage:.2f}")
    log.debug(f"Mean coverage covered: {mean_coverage_covered:.2f}")
    log.debug(f"Max covered bases: {max_covered_bases:,}")
    log.debug(f"Mean covered bases: {mean_covered_bases:.2f}")
    log.debug(f"SD: {cov_sd:.2f}")
    log.debug(f"Breadth: {breadth:.2f}")
    log.debug(f"Exp. breadth: {exp_breadth:.2f}")
    log.debug(f"Breadth/exp. ratio: {breadth_exp_ratio:.2f}")
    log.debug(f"Cov. evenness: {cov_evenness:.2f}")
    log.debug(f"C_v: {c_v:.2f}")
    log.debug(f"Mean mapq: {np.mean(read_mapq):.2f}")
    log.debug(f"GC content: {gc_content:.2f}")
    log.debug(f"Taxonomic abundance (alns): {tax_abund_aln:.2f}")
    log.debug(f"Taxonomic abundance (reads): {tax_abund_read:.2f}")
    data = BamAlignment(
        reference=reference,
        n_alns=n_alns,
        reference_length=reference_length,
        bam_reference_length=bam_reference_length,
        mean_coverage=mean_coverage,
        mean_coverage_covered=mean_coverage_covered,
        bases_covered=bases_covered,
        max_covered_bases=max_covered_bases,
        mean_covered_bases=mean_covered_bases,
        cov_evenness=cov_evenness,
        breadth=breadth,
        exp_breadth=exp_breadth,
        breadth_exp_ratio=breadth_exp_ratio,
        c_v=c_v,
        edit_distances=edit_distances,
        # edit_distances_md=edit_distances_md,
        ani_nm=ani_nm,
        # ani_md=ani_md,
        read_length=read_length,
        read_gc_content=read_gc_content,
        read_aligned_length=read_aligned_length,
        mapping_quality=read_mapq,
        read_names=set(read_names),
        read_aln_score=read_aln_score,
        tax_abund_aln=tax_abund_aln,
        tax_abund_read=tax_abund_read,
    )
    return data


class BamAlignment:
    """
    Class to store alignment information
    """

    def __init__(
        self,
        reference,
        n_alns,
        read_length,
        read_gc_content,
        read_aligned_length,
        mapping_quality,
        edit_distances,
        # edit_distances_md,
        ani_nm,
        # ani_md,
        bases_covered,
        max_covered_bases,
        mean_covered_bases,
        mean_coverage,
        mean_coverage_covered,
        reference_length,
        bam_reference_length,
        breadth,
        exp_breadth,
        breadth_exp_ratio,
        c_v,
        cov_evenness,
        read_names,
        read_aln_score,
        tax_abund_aln,
        tax_abund_read,
    ):
        self.reference = reference
        self.n_alns = n_alns
        self.read_length = read_length
        self.read_gc_content = read_gc_content
        self.read_aligned_length = read_aligned_length
        self.read_aln_score = read_aln_score
        self.mapping_quality = mapping_quality
        self.edit_distances = edit_distances
        # self.edit_distances_md = edit_distances_md
        self.ani_nm = ani_nm
        # self.ani_md = ani_md
        self.bases_covered = bases_covered
        self.max_covered_bases = max_covered_bases
        self.mean_covered_bases = mean_covered_bases
        self.mean_coverage = mean_coverage
        self.mean_coverage_covered = mean_coverage_covered
        self.reference_length = reference_length
        self.bam_reference_length = bam_reference_length
        self.breadth = breadth
        self.exp_breadth = exp_breadth
        self.breadth_exp_ratio = breadth_exp_ratio
        self.c_v = c_v
        self.cov_evenness = cov_evenness
        self.read_names = read_names
        self.tax_abund_aln = tax_abund_aln
        self.tax_abund_read = tax_abund_read
        # function to convert class to dict

    def as_dict(self):
        return {
            "reference": self.reference,
            "n_reads": self.read_names,
            "n_alns": self.n_alns,
            "read_length": self.read_length,
            "read_gc_content": self.read_gc_content,
            "read_aligned_length": self.read_aligned_length,
            "read_aln_score": self.read_aln_score,
            "mapping_quality": self.mapping_quality,
            "edit_distances": self.edit_distances,
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "ani_nm": self.ani_nm,
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "max_covered_bases": self.max_covered_bases,
            "mean_covered_bases": self.mean_covered_bases,
            "mean_coverage": self.mean_coverage,
            "mean_coverage_covered": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
            "tax_abund_read": self.tax_abund_read,
            "tax_abund_aln": self.tax_abund_aln,
        }

    def to_summary(self):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            read_length_mean = np.mean(self.read_length)
            read_length_median = np.median(self.read_length)
            read_length_std = np.std(self.read_length, ddof=1)
            read_length_max = np.max(self.read_length)
            read_length_min = np.min(self.read_length)
            read_length_mode = stats.mode(self.read_length, keepdims=True).mode[0]
            read_aligned_length = np.mean(self.read_aligned_length)
            read_aln_score = np.mean(self.read_aln_score)
            mapping_quality = np.mean(self.mapping_quality)
            edit_distances = np.mean(self.edit_distances)
            read_ani_mean = np.mean(self.ani_nm)
            read_ani_std = np.std(self.ani_nm, ddof=1)
            read_ani_median = np.median(self.ani_nm)
            gc_content = (np.sum(self.read_gc_content) / np.sum(self.read_length)) * 100
        return {
            "reference": self.reference,
            "n_reads": len(self.read_names),
            "n_alns": self.n_alns,
            "read_length_mean": read_length_mean,
            "read_length_std": read_length_std,
            "read_length_min": read_length_min,
            "read_length_max": read_length_max,
            "read_length_median": read_length_median,
            "read_length_mode": read_length_mode,
            "gc_content": gc_content,
            "read_aligned_length": read_aligned_length,
            "read_aln_score": read_aln_score,
            "mapping_quality": mapping_quality,
            "edit_distances": edit_distances,
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "read_ani_mean": read_ani_mean,
            "read_ani_std": read_ani_std,
            "read_ani_median": read_ani_median,
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "max_covered_bases": self.max_covered_bases,
            "mean_covered_bases": self.mean_covered_bases,
            "coverage_mean": self.mean_coverage,
            "coverage_covered_mean": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "bam_reference_length": self.bam_reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
            "tax_abund_read": self.tax_abund_read,
            "tax_abund_aln": self.tax_abund_aln,
        }

    def get_read_length_freqs(self):
        frags = {}
        lengths = pd.Series(self.read_length)
        lengths = lengths.value_counts().sort_index()
        freqs = list(lengths / sum(lengths))
        frags[self.reference] = {"length": list(lengths.index), "freq": freqs}
        return frags


# Inspired from https://gigabaseorgigabyte.wordpress.com/2017/04/14/getting-the-edit-distance-from-a-bam-alignment-a-journey/
def process_bam(bam, threads=1, reference_lengths=None, scale=1e6):
    """
    Processing function: calls pool of worker functions
    to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    logging.info(f"Loading BAM file")
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, "rb")

    references = samfile.references

    chr_lengths = []
    for chrom in samfile.references:
        chr_lengths.append(samfile.get_reference_length(chrom))
    max_chr_length = np.max(chr_lengths)
    pysam.set_verbosity(save)
    ref_lengths = None

    if reference_lengths is not None:
        ref_lengths = pd.read_csv(
            reference_lengths, sep="\t", index_col=0, names=["reference", "length"]
        )
        # check if the dataframe contains all the References in the BAM file
        if not set(references).issubset(set(ref_lengths.index)):
            logging.error(
                f"The BAM file contains references not found in the reference lengths file"
            )
            sys.exit(1)
        max_chr_length = np.max(ref_lengths["length"].tolist())

    if not samfile.has_index():
        logging.info(f"BAM index not found. Indexing...")
        if max_chr_length > 536870912:
            logging.info(f"A reference is longer than 2^29, indexing with csi")
            pysam.index(bam, "-c")
        else:
            pysam.index(bam)

        logging.info(f"Reloading BAM file")
        samfile = pysam.AlignmentFile(
            bam, "rb"
        )  # Need to reload the samfile after creating index

    total_refs = samfile.nreferences
    logging.info(f"Found {total_refs:,} reference sequences")
    # logging.info(f"Found {samfile.mapped:,} alignments")
    logging.info(f"Removing references without mappings...")
    # Remove references without mapped reads
    alns_in_ref = defaultdict(int)
    for aln in tqdm.tqdm(
        samfile.fetch(until_eof=True),
        total=samfile.mapped,
        leave=False,
        ncols=80,
        desc=f"Alignments processed",
    ):
        alns_in_ref[aln.reference_name] += 1

    references = [ref for ref, count in alns_in_ref.items() if count > 0]

    if len(references) == 0:
        logging.error("No reference sequences with alignments found in the BAM file")
        sys.exit(1)

    logging.info(f"Keeping {len(references):,} references")

    params = zip([bam] * len(references), references)
    try:
        logging.info(f"Getting stats for each reference...")

        if is_debug():
            data = list(
                map(
                    functools.partial(
                        get_bam_stats, ref_lengths=ref_lengths, scale=scale
                    ),
                    params,
                )
            )
        else:

            p = Pool(
                threads,
                initializer=initializer,
                initargs=([params, ref_lengths, scale],),
            )
            c_size = calc_chunksize(
                n_workers=threads, len_iterable=len(references), factor=4
            )
            data = list(
                tqdm.tqdm(
                    p.imap_unordered(
                        functools.partial(
                            get_bam_stats, ref_lengths=ref_lengths, scale=scale
                        ),
                        params,
                        chunksize=20,
                    ),
                    total=len(references),
                    leave=False,
                    ncols=80,
                    desc=f"References processed",
                )
            )

            p.close()
            p.join()

    except KeyboardInterrupt:
        logging.info(f"User canceled the operation. Terminating jobs.")
        p.terminate()
        p.join()
        sys.exit(0)
    return data


def filter_reference_BAM(
    bam,
    df,
    filter_conditions,
    threads,
    out_files,
    sort_memory,
    only_stats_filtered,
    sort_by_name=False,
):
    """Filter BAM based on certain conditions

    Args:
        bam (str): BAM file location
        data (pandas.DataFrame): Reference statistics
        filter_conditions (dict): A dictionary with the filter conditions to be used
        out_files (dict): Where to save the BAM files.
    """
    logging.info("Filtering stats...")
    logging.info(
        f"min_read_count >= {filter_conditions['min_read_count']} & min_read_length >= {filter_conditions['min_read_length']} & min_read_ani >= {filter_conditions['min_read_ani']} & min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} &  min_breadth >= {filter_conditions['min_breadth']} & min_coverage_evenness >= {filter_conditions['min_coverage_evenness']}"
    )
    df_filtered = df.loc[
        (df["n_reads"] >= filter_conditions["min_read_count"])
        & (df["read_length_mean"] >= filter_conditions["min_read_length"])
        & (df["read_ani_mean"] >= filter_conditions["min_read_ani"])
        & (df["breadth_exp_ratio"] >= filter_conditions["min_expected_breadth_ratio"])
        & (df["breadth"] >= filter_conditions["min_breadth"])
        & (df["cov_evenness"] >= filter_conditions["min_coverage_evenness"])
    ]
    if len(df_filtered.index) > 0:
        logging.info(f"Saving filtered stats...")
        df_filtered.to_csv(
            out_files["stats_filtered"], sep="\t", index=False, compression="gzip"
        )
        if only_stats_filtered:
            logging.info("Skipping saving filtered BAM file.")
        else:
            logging.info("Writing filtered BAM file... (be patient)")
            refs_dict = dict(
                zip(df_filtered["reference"], df_filtered["reference_length"])
            )
            (ref_names, ref_lengths) = zip(*refs_dict.items())

            out_bam_file = pysam.Samfile(
                out_files["bam_filtered_tmp"],
                "wb",
                referencenames=list(ref_names),
                referencelengths=list(ref_lengths),
                threads=threads,
            )
            header = pysam.AlignmentHeader.from_references(
                list(ref_names), list(ref_lengths)
            )
            references = df_filtered["reference"].values

            logging.info("Filtering BAM file...")
            samfile = pysam.AlignmentFile(bam, "rb")
            for reference in tqdm.tqdm(
                references,
                total=len(references),
                leave=False,
                ncols=80,
                desc="References processed",
            ):
                for aln in samfile.fetch(reference=reference, multiple_iterators=False):
                    out_bam_file.write(
                        pysam.AlignedSegment.fromstring(aln.to_string(), header=header)
                    )
            out_bam_file.close()
            if sort_by_name:
                logging.info("Sorting BAM file by read name...")
                pysam.sort(
                    "-n",
                    "-@",
                    str(threads),
                    "-m",
                    str(sort_memory),
                    "-o",
                    out_files["bam_filtered"],
                    out_files["bam_filtered_tmp"],
                )
            else:
                pysam.sort(
                    "-@",
                    str(threads),
                    "-m",
                    str(sort_memory),
                    "-o",
                    out_files["bam_filtered"],
                    out_files["bam_filtered_tmp"],
                )

                save = pysam.set_verbosity(0)
                samfile = pysam.AlignmentFile(out_files["bam_filtered"], "rb")
                chr_lengths = []
                for chrom in samfile.references:
                    chr_lengths.append(samfile.get_reference_length(chrom))
                max_chr_length = np.max(chr_lengths)
                pysam.set_verbosity(save)
                samfile.close()

                logging.info("BAM index not found. Indexing...")
                if max_chr_length > 536870912:
                    logging.info("A reference is longer than 2^29, indexing with csi")
                    pysam.index(out_files["bam_filtered"], "-c")
                else:
                    pysam.index(out_files["bam_filtered"])

            os.remove(out_files["bam_filtered_tmp"])
    else:
        logging.info("No references meet the filter conditions. Skipping...")


def get_alns(params):
    bam, reference = params
    samfile = pysam.AlignmentFile(bam, "rb")
    alns = []
    for aln in samfile.fetch(reference=reference, multiple_iterators=False):
        alns.append(aln.to_string())
    return alns
