import pysam
import numpy as np
import os, sys
import pandas as pd
from multiprocessing import Pool
import functools
from scipy import stats
import pysam
import tqdm
import logging
import warnings
from bam_filter.utils import is_debug, calc_chunksize, initializer
from bam_filter.entropy import entropy, norm_entropy, gini_coeff, norm_gini_coeff
from collections import Counter, defaultdict
import pyranges as pr

# import cProfile as profile
# import pstats

# import pstats

import matplotlib.pyplot as plt

log = logging.getLogger("my_logger")


sys.setrecursionlimit(10**6)


def get_tad(cov, trim_min=10, trim_max=90):
    """
    Get the TAD of a coverage
    """
    cov = cov[
        (cov >= np.percentile(cov, trim_min)) & (cov <= np.percentile(cov, trim_max))
    ]

    if np.sum(cov) == 0:
        return 0, 0
    else:
        return np.sum(cov) / len(cov), len(cov)


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
    C = float(np.rint(np.mean(coverage)))
    D2 = coverage[coverage <= C]
    # print(len(D2))
    # D2 = [x for x in coverage if x <= C]
    # print(len(D2))
    # exit()
    if len(D2) == 0:  # pragma: no cover
        covEvenness = 1.0
    else:
        if C > 0:
            covEvenness = 1.0 - (len(D2) - np.sum(D2) / C) / len(coverage)
        else:
            covEvenness = 0.0

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


def get_bam_stats(
    params,
    ref_lengths=None,
    min_read_ani=90.0,
    scale=1e6,
    trim_ends=0,
    trim_min=10,
    trim_max=90,
    plot=False,
    plots_dir="coverage-plots",
    read_length_freqs=False,
):
    """
    Worker function per chromosome
    loop over a bam file and create tuple with lists containing metrics:
    two definitions of the edit distances to the reference genome scaled by aligned read length
    """
    # prof = profile.Profile()
    # prof.enable()
    bam, references = params
    results = []
    samfile = pysam.AlignmentFile(bam, "rb")
    read_hits = defaultdict(int)
    for reference in references:
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
        cov_np = np.zeros(samfile.get_reference_length(reference), dtype=int)

        for aln in samfile.fetch(
            contig=reference, multiple_iterators=False, until_eof=False
        ):
            ani_read = (1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100
            if ani_read >= min_read_ani:
                n_alns += 1
                read_hits[aln.query_name] += 1
                if aln.has_tag("AS"):
                    read_aln_score.append(aln.get_tag("AS"))
                else:
                    read_aln_score.append(np.nan)

                if aln.has_tag("AS"):
                    edit_distances.append(aln.get_tag("NM"))
                    ani_nm.append(ani_read)
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

                cov_np[aln.reference_start : aln.reference_end] += 1

        # get bases covered by reads pileup

        # cov_pos_raw = [
        #     (pileupcolumn.pos, pileupcolumn.n)
        #     for pileupcolumn in samfile.pileup(
        #         reference,
        #         start=None,
        #         stop=None,
        #         region=None,
        #         stepper="nofilter",
        #         min_mapping_quality=0,
        #         max_depth=100000000,
        #     )
        # ]
        if trim_ends > len(cov_np) or trim_ends * 2 > len(cov_np):
            log.warning(
                f"Trimming ends ({trim_ends}) is larger than reference length ({len(cov_np)}). Disabling trimming."
            )
            trim_ends = 0

        if trim_ends > 0:
            cov_np = cov_np[trim_ends:-trim_ends]

        mean_coverage_trunc, mean_coverage_trunc_len = get_tad(
            cov_np,
            trim_min=10,
            trim_max=90,
        )

        cov_pos = cov_np[cov_np > 0]
        cov_positions = np.where(cov_np > 0)[0]
        # convert datafrane to pyranges
        ranges = create_pyranges(reference, starts, ends, strands)
        if ranges.df.shape[0] == 0:
            log.debug(f"No alignments found for {reference}")
            return None, None, read_hits
        ranges_raw = ranges.merge(strand=False)
        ranges = ranges_raw.lengths().to_list()

        max_covered_bases = np.max(ranges)
        mean_covered_bases = np.mean(ranges)
        bases_covered = int(len(cov_pos))
        # get SD from covered bases
        cov_sd = np.std(cov_pos, ddof=1)
        cov_var = np.var(cov_pos, ddof=1)
        # get average coverage
        mean_coverage = np.sum(cov_pos) / (reference_length - (2 * trim_ends))
        mean_coverage_covered = np.sum(cov_pos) / bases_covered

        breadth = bases_covered / (reference_length - (2 * trim_ends))
        exp_breadth = 1 - np.exp(-mean_coverage)
        breadth_exp_ratio = breadth / exp_breadth

        if breadth_exp_ratio > 1:
            breadth_exp_ratio = 1.0

        # fill vector with zeros to match reference length
        # cov_pos_zeroes = np.pad(
        #     cov_pos, (0, reference_length - len(cov_pos)), "constant"
        # )
        cov_evenness = coverage_evenness(cov_np)
        gc_content = (np.sum(read_gc_content) / np.sum(read_length)) * 100
        c_v = cov_sd / mean_coverage
        d_i = cov_var / mean_coverage

        read_mapq = [np.nan if x == 255 else x for x in read_mapq]

        tax_abund_aln = round((n_alns / reference_length) * scale)
        tax_abund_read = round((len(set(read_names)) / reference_length) * scale)

        # Using the trimmed mean to estimate number of reads that map to the reference
        # This is to avoid the issue of having a very high coverage region that
        # skews the mean coverage
        # C = LN / G
        # • C stands for coverage
        # • G is the haploid genome length
        # • L is the read length
        # • N is the number of reads
        #
        if mean_coverage_trunc_len > 0 and mean_coverage_trunc > 0:
            n_reads_tad = round(
                (reference_length * mean_coverage_trunc) / np.mean(read_length)
            )
            tax_abund_tad = round((n_reads_tad / mean_coverage_trunc_len) * scale)
        else:
            n_reads_tad = 0
            tax_abund_tad = 0
        # Analyse site distribution
        n_sites = len(cov_pos)
        genome_length = bam_reference_length
        # Site density (sites per thousand bp)
        site_density = 1000 * n_sites / genome_length

        counts, bins = np.histogram(
            cov_positions, bins="auto", range=(0, genome_length)
        )

        n_bins = len(bins)

        entr = entropy(counts)  # Positional entropy
        norm_entr = norm_entropy(counts)  # Normalized positional entropy
        gini = gini_coeff(counts)  # Gini coefficient
        norm_gini = norm_gini_coeff(counts)  # Normalized Gini coefficient

        log.debug(f"Number of reads: {len(set(read_names)):,}")
        log.debug(f"Number of alignments: {n_alns:,}")
        log.debug(f"Bases covered: {bases_covered:,}")
        log.debug(f"Mean coverage: {mean_coverage:.2f}")
        log.debug(f"Mean coverage (truncated): {mean_coverage_trunc:.2f}")
        log.debug(f"Reference length (truncated): {mean_coverage_trunc_len:.2f}")
        log.debug(f"Mean coverage covered: {mean_coverage_covered:.2f}")
        log.debug(f"Max covered bases: {max_covered_bases:,}")
        log.debug(f"Mean covered bases: {mean_covered_bases:.2f}")
        log.debug(f"SD: {cov_sd:.2f}")
        log.debug(f"Breadth: {breadth:.2f}")
        log.debug(f"Exp. breadth: {exp_breadth:.2f}")
        log.debug(f"Breadth/exp. ratio: {breadth_exp_ratio:.2f}")
        log.debug(f"Number of bins: {n_bins}")
        log.debug(f"Site density: {site_density:.2f}")
        log.debug(f"Entropy (H): {entr:.2f}")
        log.debug(f"Normalized entropy (H*): {norm_entr:.2f}")
        log.debug(f"Gini coefficient (G): {gini:.2f}")
        log.debug(f"Normalized Gini coefficient (G*): {norm_gini:.2f}")
        log.debug(f"Cov. evenness: {cov_evenness:.2f}")
        log.debug(f"C_v: {c_v:.2f}")
        log.debug(f"D_i: {d_i:.2f}")
        log.debug(f"Mean mapq: {np.mean(read_mapq):.2f}")
        log.debug(f"GC content: {gc_content:.2f}")
        log.debug(f"Taxonomic abundance (alns): {tax_abund_aln:.2f}")
        log.debug(f"Taxonomic abundance (reads): {tax_abund_read:.2f}")
        log.debug(f"Taxonomic abundance (TAD): {tax_abund_tad:.2f}")
        log.debug(f"Number of reads (TAD): {n_reads_tad:,}")
        if plot:
            fig, ax = plt.subplots(nrows=1, ncols=1)  # create figure & 1 axis
            # infer number of bins using Freedman-Diaconis rule
            positions_cov_zeros = pd.DataFrame(
                {"pos": range(1, bam_reference_length + 1)}
            )
            positions_cov = pd.DataFrame({"pos": cov_positions, "cov": cov_pos})
            positions_cov = positions_cov_zeros.merge(
                positions_cov, on="pos", how="left"
            )
            positions_cov["cov"] = positions_cov["cov"].fillna(0)
            positions_cov["cov"] = positions_cov["cov"].astype(int)
            positions_cov["cov_binary"] = positions_cov["cov"].apply(
                lambda x: 1 if x > 0 else 0
            )
            plt.plot(
                positions_cov["pos"],
                positions_cov["cov"],
                color="c",
                ms=0.5,
            )
            plt.suptitle(f"{reference}")
            plt.title(
                f"cov:{mean_coverage:.4f} b/e:{breadth_exp_ratio:.2f} cov_e:{cov_evenness:.2f} entropy:{norm_entr:.2f} gini:{norm_gini:.2f}"
            )
            fig.savefig(f"{plots_dir}/{reference}_coverage.png", dpi=300)
            plt.close(fig)

        data = BamAlignment(
            reference=reference,
            n_alns=n_alns,
            reference_length=reference_length,
            bam_reference_length=bam_reference_length,
            mean_coverage=mean_coverage,
            mean_coverage_trunc=mean_coverage_trunc,
            mean_coverage_trunc_len=mean_coverage_trunc_len,
            mean_coverage_covered=mean_coverage_covered,
            bases_covered=bases_covered,
            max_covered_bases=max_covered_bases,
            mean_covered_bases=mean_covered_bases,
            cov_evenness=cov_evenness,
            breadth=breadth,
            exp_breadth=exp_breadth,
            breadth_exp_ratio=breadth_exp_ratio,
            n_bins=n_bins,
            site_density=site_density,
            entropy=entr,
            norm_entropy=norm_entr,
            gini=gini,
            norm_gini=norm_gini,
            c_v=c_v,
            d_i=d_i,
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
            tax_abund_tad=tax_abund_tad,
            n_reads_tad=n_reads_tad,
        )
        results.append(data)
    samfile.close()
    # prof.disable()
    # # print profiling output
    # stats = pstats.Stats(prof).strip_dirs().sort_stats("tottime")
    # stats.print_stats(5)  # top 10 rows
    # exit()
    results = list(filter(None, results))
    data_df = pd.DataFrame([x.to_summary() for x in results])
    # read_hits = (
    #     pd.DataFrame.from_dict(read_hits, orient="index", columns=["count"])
    #     .rename_axis("read_id")
    #     .reset_index()
    # )
    if read_length_freqs:
        read_lens = [x.get_read_length_freqs() for x in results]
        return (data_df, read_lens, read_hits)
    else:
        return (
            data_df,
            None,
            read_hits,
        )


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
        mean_coverage_trunc,
        mean_coverage_trunc_len,
        mean_coverage_covered,
        reference_length,
        bam_reference_length,
        breadth,
        exp_breadth,
        breadth_exp_ratio,
        n_bins,
        site_density,
        entropy,
        norm_entropy,
        gini,
        norm_gini,
        c_v,
        d_i,
        cov_evenness,
        read_names,
        read_aln_score,
        tax_abund_aln,
        tax_abund_read,
        tax_abund_tad,
        n_reads_tad,
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
        self.mean_coverage_trunc = mean_coverage_trunc
        self.mean_coverage_trunc_len = mean_coverage_trunc_len
        self.mean_coverage_covered = mean_coverage_covered
        self.reference_length = reference_length
        self.bam_reference_length = bam_reference_length
        self.breadth = breadth
        self.exp_breadth = exp_breadth
        self.breadth_exp_ratio = breadth_exp_ratio
        self.n_bins = n_bins
        self.site_density = site_density
        self.entropy = entropy
        self.norm_entropy = norm_entropy
        self.gini = gini
        self.norm_gini = norm_gini
        self.c_v = c_v
        self.d_i = d_i
        self.cov_evenness = cov_evenness
        self.read_names = read_names
        self.tax_abund_aln = tax_abund_aln
        self.tax_abund_read = tax_abund_read
        self.tax_abund_tad = tax_abund_tad
        self.n_reads_tad = n_reads_tad
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
            "mean_coverage_trunc": self.mean_coverage_trunc,
            "mean_coverage_trunc_len": self.mean_coverage_trunc_len,
            "mean_coverage_covered": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "n_bins": self.n_bins,
            "site_density": self.site_density,
            "entropy": self.entropy,
            "norm_entropy": self.norm_entropy,
            "gini": self.gini,
            "norm_gini": self.norm_gini,
            "c_v": self.c_v,
            "d_i": self.d_i,
            "cov_evenness": self.cov_evenness,
            "tax_abund_read": self.tax_abund_read,
            "tax_abund_aln": self.tax_abund_aln,
            "tax_abund_tad": self.tax_abund_tad,
            "n_reads_tad": self.n_reads_tad,
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
            "coverage_mean_trunc": self.mean_coverage_trunc,
            "coverage_mean_trunc_len": self.mean_coverage_trunc_len,
            "coverage_covered_mean": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "bam_reference_length": self.bam_reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "n_bins": self.n_bins,
            "site_density": self.site_density,
            "entropy": self.entropy,
            "norm_entropy": self.norm_entropy,
            "gini": self.gini,
            "norm_gini": self.norm_gini,
            "c_v": self.c_v,
            "d_i": self.d_i,
            "cov_evenness": self.cov_evenness,
            "tax_abund_read": self.tax_abund_read,
            "tax_abund_aln": self.tax_abund_aln,
            "tax_abund_tad": self.tax_abund_tad,
            "n_reads_tad": self.n_reads_tad,
        }

    def get_read_length_freqs(self):
        frags = {}
        lengths = pd.Series(self.read_length)
        lengths = lengths.value_counts().sort_index()
        freqs = list(lengths / np.sum(lengths))
        frags[self.reference] = {"length": list(lengths.index), "freq": freqs}
        return frags


# Inspired from https://gigabaseorgigabyte.wordpress.com/2017/04/14/getting-the-edit-distance-from-a-bam-alignment-a-journey/
def process_bam(
    bam,
    threads=1,
    reference_lengths=None,
    min_read_ani=90.0,
    min_read_count=10,
    trim_ends=0,
    trim_min=10,
    trim_max=90,
    scale=1e6,
    sort_memory="1G",
    plot=False,
    plots_dir="coverage-plots",
    chunksize=None,
    read_length_freqs=False,
):
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

    # Check if BAM files is not sorted by coordinates, sort it by coordinates
    if not samfile.header["HD"]["SO"] == "coordinate":
        log.info("BAM file is not sorted by coordinates, sorting it...")
        sorted_bam = bam.replace(".bam", ".bf-sorted.bam")
        pysam.sort("-@", str(threads), "-m", str(sort_memory), "-o", sorted_bam, bam)
        bam = sorted_bam
        samfile = pysam.AlignmentFile(bam, "rb")

    if not samfile.has_index():
        logging.info("BAM index not found. Indexing...")
        if max_chr_length > 536870912:
            logging.info("A reference is longer than 2^29, indexing with csi")
            pysam.index("-c", "-@", str(threads), bam)
        else:
            pysam.index(
                "-@",
                str(threads),
                bam,
            )

        logging.info("Reloading BAM file")
        samfile = pysam.AlignmentFile(
            bam, "rb"
        )  # Need to reload the samfile after creating index

    total_refs = samfile.nreferences
    logging.info(f"Found {total_refs:,} reference sequences")
    # logging.info(f"Found {samfile.mapped:,} alignments")
    logging.info(
        f"Removing references without mappings or less than {min_read_count} reads..."
    )
    # Remove references without mapped reads
    # alns_in_ref = defaultdict(int)
    # for aln in tqdm.tqdm(
    #     samfile.fetch(until_eof=True),
    #     total=samfile.mapped,
    #     leave=False,
    #     ncols=80,
    #     desc=f"Alignments processed",
    # ):
    #     alns_in_ref[aln.reference_name] += 1

    if plot:
        # Check if image folder exists
        if not os.path.exists(plots_dir):
            os.makedirs(plots_dir)

    references = [
        chrom.contig
        for chrom in samfile.get_index_statistics()
        if chrom.mapped >= min_read_count
    ]

    if len(references) == 0:
        logging.error("No reference sequences with alignments found in the BAM file")
        sys.exit(1)

    logging.info(f"Keeping {len(references):,} references")

    if (chunksize is not None) and ((len(references) // chunksize) > threads):
        c_size = chunksize
    else:
        c_size = calc_chunksize(
            n_workers=threads, len_iterable=len(references), factor=4
        )
    ref_chunks = [references[i : i + c_size] for i in range(0, len(references), c_size)]
    params = zip([bam] * len(ref_chunks), ref_chunks)
    try:
        logging.info(
            f"Processing {len(ref_chunks):,} chunks of {c_size:,} references each"
        )
        if is_debug():
            data = list(
                map(
                    functools.partial(
                        get_bam_stats,
                        ref_lengths=ref_lengths,
                        min_read_ani=min_read_ani,
                        trim_ends=0,
                        trim_min=trim_min,
                        trim_max=trim_max,
                        scale=scale,
                        plot=plot,
                        plots_dir=plots_dir,
                        read_length_freqs=read_length_freqs,
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

            data = list(
                tqdm.tqdm(
                    p.imap_unordered(
                        functools.partial(
                            get_bam_stats,
                            ref_lengths=ref_lengths,
                            min_read_ani=min_read_ani,
                            trim_ends=0,
                            trim_min=trim_min,
                            trim_max=trim_max,
                            scale=scale,
                            plot=plot,
                            plots_dir=plots_dir,
                            read_length_freqs=read_length_freqs,
                        ),
                        params,
                        chunksize=1,
                    ),
                    total=len(ref_chunks),
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
    transform_cov_evenness=False,
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
    if "min_norm_entropy" in filter_conditions and "min_norm_gini" in filter_conditions:
        logging.info(
            f"min_read_count >= {filter_conditions['min_read_count']} & min_read_length >= {filter_conditions['min_read_length']} & min_avg_read_ani >= {filter_conditions['min_avg_read_ani']} & min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} & min_breadth >= {filter_conditions['min_breadth']} & min_coverage_evenness >= {filter_conditions['min_coverage_evenness']} & min_coverage_mean >= {filter_conditions['min_coverage_mean']} & min_norm_entropy >= {filter_conditions['min_norm_entropy']} & min_norm_gini <= {filter_conditions['min_norm_gini']}"
        )
        # We transform the coverage_evenenness to 1.0 where the coverage is smaller than 1
        if transform_cov_evenness is True:
            df["cov_evenness_tmp"] = df["cov_evenness"]
            df["cov_evenness_tmp"] = np.where(
                np.rint(df.coverage_mean) < 1.0, 1.0, df.cov_evenness_tmp
            )
        else:
            df["cov_evenness_tmp"] = df["cov_evenness"]
        df_filtered = df.loc[
            (df["n_reads"] >= filter_conditions["min_read_count"])
            & (df["read_length_mean"] >= filter_conditions["min_read_length"])
            & (df["read_ani_mean"] >= filter_conditions["min_avg_read_ani"])
            & (
                df["breadth_exp_ratio"]
                >= filter_conditions["min_expected_breadth_ratio"]
            )
            & (df["breadth"] >= filter_conditions["min_breadth"])
            & (df["cov_evenness_tmp"] >= filter_conditions["min_coverage_evenness"])
            & (df["coverage_mean"] >= filter_conditions["min_coverage_mean"])
            & (df["norm_entropy"] >= filter_conditions["min_norm_entropy"])
            & (df["norm_gini"] <= filter_conditions["min_norm_gini"])
        ]
    else:
        logging.info(
            f"min_read_count >= {filter_conditions['min_read_count']} & min_read_length >= {filter_conditions['min_read_length']} & min_avg_read_ani >= {filter_conditions['min_avg_read_ani']} & min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} & min_breadth >= {filter_conditions['min_breadth']} & min_coverage_evenness >= {filter_conditions['min_coverage_evenness']} & min_coverage_mean >= {filter_conditions['min_coverage_mean']}"
        )
        if transform_cov_evenness is True:
            df["cov_evenness_tmp"] = df["cov_evenness"]
            df["cov_evenness_tmp"] = np.where(
                np.rint(df.coverage_mean) < 1.0, 1.0, df.cov_evenness_tmp
            )
        else:
            df["cov_evenness_tmp"] = df["cov_evenness"]

        df_filtered = df.loc[
            (df["n_reads"] >= filter_conditions["min_read_count"])
            & (df["read_length_mean"] >= filter_conditions["min_read_length"])
            & (df["read_ani_mean"] >= filter_conditions["min_avg_read_ani"])
            & (
                df["breadth_exp_ratio"]
                >= filter_conditions["min_expected_breadth_ratio"]
            )
            & (df["breadth"] >= filter_conditions["min_breadth"])
            & (df["cov_evenness_tmp"] >= filter_conditions["min_coverage_evenness"])
            & (df["coverage_mean"] >= filter_conditions["min_coverage_mean"])
        ]

    del df_filtered["cov_evenness_tmp"]

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
                for aln in samfile.fetch(
                    reference=reference, multiple_iterators=False, until_eof=True
                ):
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
                    pysam.index(
                        "-c",
                        "-@",
                        str(threads),
                        out_files["bam_filtered"],
                    )
                else:
                    pysam.index(
                        "-@",
                        str(threads),
                        out_files["bam_filtered"],
                    )

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
