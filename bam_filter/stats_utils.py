import numpy as np
import pandas as pd
import logging
import tqdm
import warnings
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import pysam
from bam_filter.entropy import entropy, norm_entropy, gini_coeff, norm_gini_coeff
from bam_filter.utils import is_debug

log = logging.getLogger("my_logger")


def get_tad(cov, trim_min=10, trim_max=90):
    """
    Get the TAD (Truncated Average Depth) of coverage more efficiently

    Args:
        cov: numpy array of coverage values
        trim_min: minimum percentile to include (default: 10)
        trim_max: maximum percentile to include (default: 90)

    Returns:
        tuple: (TAD value, number of positions used)
    """
    if len(cov) == 0:
        return (0, 0)

    # Calculate percentiles once
    min_val = np.percentile(cov, trim_min)
    max_val = np.percentile(cov, trim_max)

    # Use boolean indexing in one step
    mask = (cov >= min_val) & (cov <= max_val)
    filtered_cov = cov[mask]

    count = len(filtered_cov)
    return (np.sum(filtered_cov) / count, count) if count > 0 else (0, 0)


def coverage_evenness(coverage):
    """
    Calculate the evenness of coverage using a more efficient implementation.

    Formula: E = 1 - (|D2| - sum(D2)/C) / |D|
    where:
    - C is rounded mean coverage
    - D2 is subset of positions with coverage <= C
    - |D| is total number of positions

    Args:
        coverage: numpy array of coverage values

    Returns:
        float: Evenness score between 0 and 1
        - 1.0 indicates perfectly even coverage
        - 0.0 indicates completely uneven coverage or zero coverage
    """
    # Handle empty input
    if len(coverage) == 0:
        return 0.0

    # Calculate mean coverage once and round
    C = float(np.rint(np.mean(coverage)))

    # Early returns for edge cases
    if C == 0:
        return 0.0

    # Use boolean indexing for efficient filtering
    mask = coverage <= C
    D2 = coverage[mask]

    if len(D2) == 0:
        return 1.0

    # Calculate evenness score in one step
    return 1.0 - (len(D2) - np.sum(D2) / C) / len(coverage)


def calc_gc_content(seq):
    """Calculate GC content of a sequence"""
    return seq.count("G") + seq.count("C")


import numpy as np


def merge_intervals(starts, ends):
    """
    Merge overlapping intervals efficiently using numpy arrays.

    Args:
        starts: array-like of interval start positions
        ends: array-like of interval end positions

    Returns:
        tuple: (merged_starts, merged_ends) as lists

    Example:
        >>> merge_intervals([1,2,8], [5,6,10])
        ([1, 8], [6, 10])
    """
    if not starts:
        return [], []

    # Convert to numpy arrays if not already
    if not isinstance(starts, np.ndarray):
        starts = np.array(starts)
        ends = np.array(ends)

    # Sort by start positions
    idx = np.argsort(starts)
    starts = starts[idx]
    ends = ends[idx]

    # Find where new intervals start
    new_interval_mask = starts[1:] > ends[:-1]

    # Get indices where new intervals start
    split_points = np.where(new_interval_mask)[0] + 1

    if len(split_points) == 0:
        # All intervals overlap
        return [starts[0]], [ends.max()]

    # Split into groups of overlapping intervals
    start_groups = np.split(starts, split_points)
    end_groups = np.split(ends, split_points)

    # Take first start and max end from each group
    merged_starts = [group[0] for group in start_groups]
    merged_ends = [group.max() for group in end_groups]

    return merged_starts, merged_ends


class BamAlignment:
    """Class to store alignment information"""

    def __init__(
        self,
        reference,
        n_alns,
        read_length,
        read_gc_content,
        read_aligned_length,
        mapping_quality,
        edit_distances,
        ani_nm,
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
        self.ani_nm = ani_nm
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
            "read_ani_mean": read_ani_mean,
            "read_ani_std": read_ani_std,
            "read_ani_median": read_ani_median,
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
        """Calculate read length frequencies"""
        lengths = pd.Series(self.read_length)
        lengths = lengths.value_counts().sort_index()
        freqs = list(lengths / np.sum(lengths))
        return {self.reference: {"length": list(lengths.index), "freq": freqs}}


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
    threads=1,
):
    """Calculate comprehensive statistics from BAM file alignments"""
    bam, references = params
    results = []

    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
        bam_reference_lengths = {
            reference: np.int64(samfile.get_reference_length(reference))
            for reference in references
        }
        read_hits = defaultdict(int)

        for reference in references:
            edit_distances = []
            ani_nm = []
            read_length = []
            read_aligned_length = []
            read_mapq = []
            read_aln_score = []
            read_names = set()
            read_gc_content = []
            n_alns = 0

            if ref_lengths is None:
                reference_length = bam_reference_lengths[reference]
            else:
                reference_length = np.int64(ref_lengths.loc[reference, "length"])
            bam_reference_length = bam_reference_lengths[reference]

            log.debug(f"Processing reference {reference}")
            log.debug(f"Reference length: {reference_length:,}")
            log.debug(f"BAM reference length: {bam_reference_length:,}")

            starts = []
            ends = []
            strands = []
            cov_np = np.zeros(samfile.get_reference_length(reference), dtype=np.int32)

            for aln in samfile.fetch(reference, multiple_iterators=False):
                ani_read = (1 - ((aln.get_tag("NM") / aln.infer_query_length()))) * 100
                if ani_read >= min_read_ani:
                    n_alns += 1
                    read_hits[aln.query_name] += 1

                    # Collect alignment statistics
                    if aln.has_tag("AS"):
                        read_aln_score.append(aln.get_tag("AS"))
                        edit_distances.append(aln.get_tag("NM"))
                        ani_nm.append(ani_read)
                    else:
                        read_aln_score.append(np.nan)
                        edit_distances.append(np.nan)
                        ani_nm.append(np.nan)

                    read_gc_content.append(calc_gc_content(aln.query_sequence))
                    read_length.append(aln.infer_read_length())
                    read_aligned_length.append(aln.query_alignment_length)
                    read_mapq.append(aln.mapping_quality)
                    read_names.add(aln.query_name)

                    # Record alignment position info
                    strand = "-" if aln.is_reverse else "+"
                    starts.append(aln.reference_start)
                    ends.append(aln.reference_end)
                    strands.append(strand)

                    # Update coverage array
                    cov_np[aln.reference_start : aln.reference_end] += 1

            if n_alns > 0:
                if trim_ends > len(cov_np) or trim_ends * 2 > len(cov_np):
                    log.warning(
                        f"Trimming ends ({trim_ends}) is larger than reference length ({len(cov_np)}). Disabling trimming."
                    )
                    trim_ends = 0

                if trim_ends > 0:
                    cov_np = cov_np[trim_ends:-trim_ends]

                # Calculate coverage statistics
                cov_pos = cov_np[cov_np > 0]
                cov_positions = np.nonzero(cov_np)[0]

                # Use merge_intervals instead of PyRanges
                merged_starts, merged_ends = merge_intervals(starts, ends)
                ranges = np.array(
                    [end - start for start, end in zip(merged_starts, merged_ends)]
                )

                max_covered_bases = np.max(ranges)
                mean_covered_bases = np.mean(ranges)
                bases_covered = len(cov_pos)

                # Calculate coverage metrics
                mean_coverage_trunc, mean_coverage_trunc_len = get_tad(
                    cov_np,
                    trim_min=trim_min,
                    trim_max=trim_max,
                )

                # Coverage calculations
                mean_coverage = np.sum(cov_pos) / (reference_length - (2 * trim_ends))
                mean_coverage_covered = (
                    np.sum(cov_pos) / bases_covered if bases_covered > 0 else 0
                )

                breadth = bases_covered / (reference_length - (2 * trim_ends))
                exp_breadth = 1 - np.exp(-mean_coverage)
                breadth_exp_ratio = (
                    min(breadth / exp_breadth, 1.0) if exp_breadth > 0 else 0.0
                )

                # Calculate evenness and statistical metrics
                cov_evenness = coverage_evenness(cov_np)
                c_v = (
                    np.std(cov_pos, ddof=1) / mean_coverage if mean_coverage > 0 else 0
                )
                d_i = (
                    np.var(cov_pos, ddof=1) / mean_coverage if mean_coverage > 0 else 0
                )

                # Convert mapping quality 255 to NaN
                read_mapq = [np.nan if x == 255 else x for x in read_mapq]

                # Calculate abundance metrics
                tax_abund_aln = round((n_alns / reference_length) * scale)
                tax_abund_read = round((len(read_names) / reference_length) * scale)

                if mean_coverage_trunc_len > 0 and mean_coverage_trunc > 0:
                    n_reads_tad = round(
                        (reference_length * mean_coverage_trunc) / np.mean(read_length)
                    )
                    tax_abund_tad = round(
                        (n_reads_tad / mean_coverage_trunc_len) * scale
                    )
                else:
                    n_reads_tad = 0
                    tax_abund_tad = 0
                # Analyze site distribution
                counts, bins = np.histogram(
                    cov_positions, bins="auto", range=(0, bam_reference_length)
                )
                n_bins = len(bins)
                site_density = 1000 * len(cov_pos) / bam_reference_length

                # Calculate entropy metrics
                entr = entropy(counts)
                norm_entr = norm_entropy(counts)
                gini = gini_coeff(counts)
                norm_gini = norm_gini_coeff(counts)

                if plot:
                    _plot_coverage(
                        reference,
                        cov_positions,
                        cov_pos,
                        bam_reference_length,
                        mean_coverage,
                        breadth_exp_ratio,
                        cov_evenness,
                        norm_entr,
                        norm_gini,
                        plots_dir,
                    )

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
                    ani_nm=ani_nm,
                    read_length=read_length,
                    read_gc_content=read_gc_content,
                    read_aligned_length=read_aligned_length,
                    mapping_quality=read_mapq,
                    read_names=read_names,
                    read_aln_score=read_aln_score,
                    tax_abund_aln=tax_abund_aln,
                    tax_abund_read=tax_abund_read,
                    tax_abund_tad=tax_abund_tad,
                    n_reads_tad=n_reads_tad,
                )
                results.append(data)

    results = list(filter(None, results))
    data_df = pd.DataFrame([x.to_summary() for x in results])

    if read_length_freqs:
        read_lens = [x.get_read_length_freqs() for x in results]
        return (data_df, read_lens, read_hits)
    else:
        return (data_df, None, read_hits)


def _plot_coverage(
    reference,
    cov_positions,
    cov_pos,
    bam_reference_length,
    mean_coverage,
    breadth_exp_ratio,
    cov_evenness,
    norm_entr,
    norm_gini,
    plots_dir,
):
    """Create coverage plot for a reference sequence"""
    fig, ax = plt.subplots(nrows=1, ncols=1)
    positions_cov_zeros = pd.DataFrame({"pos": range(1, bam_reference_length + 1)})
    positions_cov = pd.DataFrame({"pos": cov_positions, "cov": cov_pos})
    positions_cov = positions_cov_zeros.merge(positions_cov, on="pos", how="left")
    positions_cov["cov"] = positions_cov["cov"].fillna(0)
    positions_cov["cov"] = positions_cov["cov"].astype(int)

    plt.plot(
        positions_cov["pos"],
        positions_cov["cov"],
        color="c",
        ms=0.5,
    )
    plt.suptitle(f"{reference}")
    plt.title(
        f"cov:{mean_coverage:.4f} b/e:{breadth_exp_ratio:.2f} "
        f"cov_e:{cov_evenness:.2f} entropy:{norm_entr:.2f} gini:{norm_gini:.2f}"
    )
    fig.savefig(f"{plots_dir}/{reference}_coverage.png", dpi=300)
    plt.close(fig)
