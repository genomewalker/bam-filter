import pysam
import numpy as np
import re, os, sys
import pandas as pd
from multiprocessing import Pool
from scipy import stats
import pysam
import tqdm
import logging
from bam_filter.utils import is_debug, calc_chunksize, fast_flatten

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


def get_bam_stats(params):
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
    n_alns = 0
    reference_length = int(samfile.get_reference_length(reference))

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

        read_length.append(aln.infer_read_length())
        read_aligned_length.append(aln.query_alignment_length)
        read_mapq.append(aln.mapping_quality)
        read_names.append(aln.query_name)
    if n_alns > 1:
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

        bases_covered = int(len(cov_pos))
        # get SD from covered bases
        cov_sd = np.std(cov_pos)
        # get average coverage
        mean_coverage = sum(cov_pos) / reference_length
        mean_coverage_covered = sum(cov_pos) / bases_covered

        breadth = bases_covered / reference_length
        exp_breadth = 1 - np.exp(-mean_coverage)
        breadth_exp_ratio = breadth / exp_breadth

        cov_evenness = coverage_evenness(cov_pos)

        c_v = cov_sd / mean_coverage
        read_mapq = [np.nan if x == 255 else x for x in read_mapq]
        data = BamAlignment(
            reference=reference,
            n_alns=n_alns,
            reference_length=reference_length,
            mean_coverage=mean_coverage,
            mean_coverage_covered=mean_coverage_covered,
            bases_covered=bases_covered,
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
            read_aligned_length=read_aligned_length,
            mapping_quality=read_mapq,
            read_names=read_names,
            read_aln_score=read_aln_score,
        )
    else:
        data = BamAlignment(
            reference=reference,
            n_alns=n_alns,
            reference_length=reference_length,
            mean_coverage=np.nan,
            mean_coverage_covered=np.nan,
            bases_covered=np.nan,
            cov_evenness=np.nan,
            breadth=np.nan,
            exp_breadth=np.nan,
            breadth_exp_ratio=np.nan,
            c_v=np.nan,
            edit_distances=np.nan,
            # edit_distances_md=np.nan,
            ani_nm=np.nan,
            # ani_md=np.nan,
            read_length=np.nan,
            read_aligned_length=np.nan,
            mapping_quality=np.nan,
            read_names=read_names,
            read_aln_score=read_aln_score,
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
        read_aligned_length,
        mapping_quality,
        edit_distances,
        # edit_distances_md,
        ani_nm,
        # ani_md,
        bases_covered,
        mean_coverage,
        mean_coverage_covered,
        reference_length,
        breadth,
        exp_breadth,
        breadth_exp_ratio,
        c_v,
        cov_evenness,
        read_names,
        read_aln_score,
    ):
        self.reference = reference
        self.n_alns = n_alns
        self.read_length = read_length
        self.read_aligned_length = read_aligned_length
        self.read_aln_score = read_aln_score
        self.mapping_quality = mapping_quality
        self.edit_distances = edit_distances
        # self.edit_distances_md = edit_distances_md
        self.ani_nm = ani_nm
        # self.ani_md = ani_md
        self.bases_covered = bases_covered
        self.mean_coverage = mean_coverage
        self.mean_coverage_covered = mean_coverage_covered
        self.reference_length = reference_length
        self.breadth = breadth
        self.exp_breadth = exp_breadth
        self.breadth_exp_ratio = breadth_exp_ratio
        self.c_v = c_v
        self.cov_evenness = cov_evenness
        self.read_names = read_names
        # function to convert class to dict

    def as_dict(self):
        return {
            "reference": self.reference,
            "n_reads": self.read_names,
            "n_alns": self.n_alns,
            "read_length": self.read_length,
            "read_aligned_length": self.read_aligned_length,
            "read_aln_score": self.read_aln_score,
            "mapping_quality": self.mapping_quality,
            "edit_distances": self.edit_distances,
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "ani_nm": self.ani_nm,
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "mean_coverage": self.mean_coverage,
            "mean_coverage_covered": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
        }

    def to_summary(self):

        return {
            "reference": self.reference,
            "n_reads": len(self.read_names),
            "n_alns": self.n_alns,
            "read_length_mean": np.mean(self.read_length),
            "read_length_median": np.median(self.read_length),
            "read_length_mode": stats.mode(self.read_length)[0][0],
            "read_aligned_length": np.mean(self.read_aligned_length),
            "read_aln_score": np.mean(self.read_aln_score),
            "mapping_quality": np.nanmean(self.mapping_quality),
            "edit_distances": np.mean(self.edit_distances),
            # "edit_distances_md": np.mean(self.edit_distances_md),
            "read_ani_mean": np.mean(self.ani_nm),
            # "ani_md": np.mean(self.ani_md),
            "bases_covered": self.bases_covered,
            "coverage_mean": self.mean_coverage,
            "coverage_covered_mean": self.mean_coverage_covered,
            "reference_length": self.reference_length,
            "breadth": self.breadth,
            "exp_breadth": self.exp_breadth,
            "breadth_exp_ratio": self.breadth_exp_ratio,
            "c_v": self.c_v,
            "cov_evenness": self.cov_evenness,
        }


def process_bam(bam, threads=1):
    """
    Processing function: calls pool of worker functions
    to extract from a bam file two definitions of the edit distances to the reference genome scaled by read length
    Returned in a pandas DataFrame
    """
    logging.info(f"Loading BAM file")
    save = pysam.set_verbosity(0)
    samfile = pysam.AlignmentFile(bam, "rb")
    pysam.set_verbosity(save)

    if not samfile.has_index():
        logging.info(f"BAM index not found. Indexing...")
        pysam.index(bam)

        logging.info(f"Reloading BAM file")
        samfile = pysam.AlignmentFile(
            bam, "rb"
        )  # Need to reload the samfile after creating index

    total_refs = samfile.nreferences
    logging.info(f"Found {total_refs:,} reference sequences")
    logging.info(f"Found {samfile.mapped:,} aligned reads")

    references = samfile.references[0:100]
    params = zip([bam] * len(references), references)
    try:
        logging.info(f"Getting stats for each reference...")

        if is_debug():
            data = list(map(get_bam_stats, params))
        else:

            p = Pool(threads)
            c_size = calc_chunksize(threads, len(references))
            data = list(
                tqdm.tqdm(
                    p.imap_unordered(get_bam_stats, params, chunksize=c_size),
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
        sys.exit()
    return data


def filter_reference_BAM(bam, df, filter_conditions, threads, out_files, sort_memory):
    """Filter BAM based on certain conditions

    Args:
        bam (str): BAM file location
        data (pandas.DataFrame): Reference statistics
        filter_conditions (dict): A dictionary with the filter conditions to be used
        out_files (dict): Where to save the BAM files.
    """
    logging.info(f"Filtering stats...")
    logging.info(
        f"min_read_count >= {filter_conditions['min_read_count']} & min_read_length >= {filter_conditions['min_read_length']} & min_read_ani >= {filter_conditions['min_read_ani']} & min_expected_breadth_ratio >= {filter_conditions['min_expected_breadth_ratio']} & min_coverage_evenness >= {filter_conditions['min_coverage_evenness']}"
    )

    df_filtered = df.loc[
        (df["n_reads"] >= filter_conditions["min_read_count"])
        & (df["read_length_mean"] >= filter_conditions["min_read_length"])
        & (df["read_ani_mean"] >= filter_conditions["min_read_ani"])
        & (df["breadth_exp_ratio"] >= filter_conditions["min_expected_breadth_ratio"])
        & (df["cov_evenness"] >= filter_conditions["min_coverage_evenness"])
    ]

    refs_dict = dict(zip(df_filtered["reference"], df_filtered["reference_length"]))
    (ref_names, ref_lengths) = zip(*refs_dict.items())

    out_bam_file = pysam.Samfile(
        out_files["bam_filtered_tmp"],
        "wb",
        referencenames=list(ref_names),
        referencelengths=list(ref_lengths),
        threads=threads,
    )
    header = pysam.AlignmentHeader.from_references(list(ref_names), list(ref_lengths))
    logging.info(f"Writing filtered BAM file... (be patient)")
    references = df_filtered["reference"].values
    params = zip([bam] * len(references), references)
    logging.info(f"Saving filtered stats...")
    df_filtered.to_csv(
        out_files["stats_filtered"], sep="\t", index=False, compression="gzip"
    )
    try:
        logging.info(f"Filtering BAM file...")

        if is_debug():
            alns = list(map(get_alns, params))
        else:

            p = Pool(threads)
            c_size = calc_chunksize(threads, len(references))
            alns = list(
                tqdm.tqdm(
                    p.imap_unordered(get_alns, params, chunksize=c_size),
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
        sys.exit()

    samfile = pysam.AlignmentFile(bam, "rb")
    for aln in fast_flatten(alns):
        out_bam_file.write(pysam.AlignedSegment.fromstring(aln, header=header))
    out_bam_file.close()
    pysam.sort(
        "-@",
        str(threads),
        "-m",
        str(sort_memory),
        "-o",
        out_files["bam_filtered"],
        out_files["bam_filtered_tmp"],
    )
    pysam.index(out_files["bam_filtered"])
    os.remove(out_files["bam_filtered_tmp"])


def get_alns(params):
    bam, reference = params
    samfile = pysam.AlignmentFile(bam, "rb")
    alns = []
    for aln in samfile.fetch(reference=reference, multiple_iterators=False):
        alns.append(aln.to_string())
    return alns
