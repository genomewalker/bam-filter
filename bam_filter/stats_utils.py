import numpy as np
import pandas as pd
import logging
import warnings
from collections import defaultdict
from scipy import stats
import matplotlib.pyplot as plt
import pysam
from bam_filter.entropy import entropy, norm_entropy, gini_coeff, norm_gini_coeff
from numba import jit

log = logging.getLogger("my_logger")


def initializer_lengths(init_dict, init_dict2):
    global ref_lengths
    ref_lengths = init_dict
    global bam_reference_lengths
    bam_reference_lengths = init_dict2


# @jit(nopython=True, cache=True)
# def get_tad(cov, trim_min=10, trim_max=90):
#     """
#     Fast TAD calculation that maintains perfect accuracy.
#     Uses direct integer counting for exact results.
#     """
#     if len(cov) == 0:
#         return 0.0, 0

#     n = len(cov)
#     min_count = (n * trim_min) // 100
#     max_count = (n * trim_max) // 100

#     if min_count >= max_count:
#         return 0.0, 0

#     # For small arrays, just use sort directly
#     if n < 500:
#         sorted_cov = np.sort(cov)
#         filtered = sorted_cov[min_count:max_count]
#         return np.sum(filtered) / len(filtered), len(filtered)

#     # Find max and min values in one pass
#     max_val = cov[0]
#     for i in range(n):
#         if cov[i] > max_val:
#             max_val = cov[i]

#     # Use stable counting sort
#     counts = np.zeros(int(max_val) + 1, dtype=np.int32)
#     for i in range(n):
#         counts[int(cov[i])] += 1

#     # Calculate running sum and find trim points
#     curr_count = 0
#     start_val = -1
#     end_val = -1

#     # Find start value
#     for i in range(len(counts)):
#         curr_count += counts[i]
#         if curr_count > min_count:
#             start_val = i
#             break

#     # Reset and find end value
#     curr_count = 0
#     for i in range(len(counts)):
#         curr_count += counts[i]
#         if curr_count > max_count:
#             end_val = i
#             break

#     if start_val == -1 or end_val == -1:
#         return 0.0, 0

#     # Calculate exact sum and count
#     total = 0.0
#     count = 0
#     for i in range(start_val, end_val + 1):
#         total += i * counts[i]
#         count += counts[i]

#     if count == 0:
#         return 0.0, 0

#     return total / count, count


@jit(nopython=True, cache=True)
def get_tad(cov, trim_min=10, trim_max=90):
    """
    Ultrafast TAD calculation optimized for coverage data characteristics:
    - Non-negative integers
    - Often sparse (many zeros)
    - Values typically cluster in lower ranges
    """
    n = len(cov)
    if n == 0:
        return 0.0, 0

    # Early return for uniform coverage
    if cov[0] == cov[-1]:  # Check if first and last are same
        if cov[n // 2] == cov[0]:  # Check middle
            return float(cov[0]), n

    min_count = (n * trim_min) // 100
    max_count = (n * trim_max) // 100

    if min_count >= max_count:
        return 0.0, 0

    # For tiny arrays, direct sort
    if n < 32:  # Reduced threshold since we handle uniform case above
        sorted_cov = np.sort(cov)
        filtered = sorted_cov[min_count:max_count]
        return np.sum(filtered) / len(filtered), len(filtered)

    # Count zeros up front since they're common in coverage data
    zeros = 0
    for i in range(n):
        if cov[i] == 0:
            zeros += 1

    # Skip zeros array allocation if all values are above trimming point
    if zeros < min_count:
        counts = np.zeros(int(cov.max()) + 1, dtype=np.int32)
        start_idx = 1  # Skip zero bin
    else:
        counts = np.zeros(int(cov.max()) + 1, dtype=np.int32)
        counts[0] = zeros
        start_idx = 0

    # Count only non-zero values
    for i in range(n):
        val = cov[i]
        if val > 0:  # Skip zeros since we counted them
            counts[val] += 1

    # Direct trim point search optimized for coverage distribution
    start_val = -1
    end_val = -1
    curr_count = 0

    for i in range(start_idx, len(counts)):
        curr_count += counts[i]
        if start_val == -1 and curr_count > min_count:
            start_val = i
        if curr_count > max_count:
            end_val = i
            break

    if start_val == -1 or end_val == -1:
        return 0.0, 0

    # Fast sum calculation
    total = 0.0
    count = 0
    for i in range(start_val, end_val + 1):
        c = counts[i]
        total += i * c
        count += c

    if count == 0:
        return 0.0, 0

    return total / count, count


# @jit(nopython=True, cache=True)
# def coverage_evenness(coverage):
#     """
#     Optimized coverage evenness calculation that avoids nonzero()
#     """
#     if len(coverage) == 0:
#         return 0.0

#     # Calculate mean directly without np.mean()
#     total = 0.0
#     n = len(coverage)
#     for i in range(n):
#         total += coverage[i]
#     C = round(total / n)

#     if C == 0:
#         return 0.0

#     # Count values <= C and their sum directly
#     count_leq_C = 0
#     sum_leq_C = 0.0
#     for i in range(n):
#         val = coverage[i]
#         if val <= C:
#             count_leq_C += 1
#             sum_leq_C += val

#     return 1.0 - (count_leq_C - sum_leq_C / C) / n


@jit(nopython=True, cache=True, fastmath=True)
def get_coverage_stats(coverage):
    """
    Ultra-fast coverage stats using minimal operations
    Returns: (evenness, number_of_positions, positions)
    """
    n = len(coverage)
    positions = np.empty(n, dtype=np.int64)  # Preallocate max size
    pos_idx = 0
    total = 0.0

    # Single pass for both positions and total
    for i in range(n):
        val = coverage[i]
        if val > 0:
            positions[pos_idx] = i
            pos_idx += 1
            total += val

    # Early exit for empty
    if pos_idx == 0:
        return 0.0, 0, positions[:0]

    # Fast C calculation
    C = round(total / n)
    if C == 0:
        return 0.0, pos_idx, positions[:pos_idx]

    # Direct count and sum for values <= C
    count_leq_C = 0
    sum_leq_C = 0.0

    # Raw loop is faster than numpy operations for this case
    for i in range(n):
        val = coverage[i]
        if val <= C:
            count_leq_C += 1
            sum_leq_C += val

    # Final calculation
    evenness = 1.0 - (count_leq_C - sum_leq_C / C) / n
    return evenness, pos_idx, positions[:pos_idx]


def calc_gc_content(seq):
    """Calculate GC content of a sequence"""
    return seq.count("G") + seq.count("C")


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
        """Calculate read length frequencies using optimized methods"""
        unique_lengths, counts = np.unique(self.read_length, return_counts=True)
        freqs = counts / counts.sum()
        return {
            self.reference: {"length": unique_lengths.tolist(), "freq": freqs.tolist()}
        }


# import cProfile, pstats



def get_bam_stats(
    params,
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
    # profiler = cProfile.Profile()
    # profiler.enable()

    """Calculate comprehensive statistics from BAM file alignments"""
    bam, references = params
    results = []

    with pysam.AlignmentFile(bam, "rb", threads=threads) as samfile:
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
            min_read_ani /= 100

            reference_length = ref_lengths[reference]
            bam_reference_length = bam_reference_lengths[reference]

            log.debug(f"Processing reference {reference}")
            log.debug(f"Reference length: {reference_length:,}")
            log.debug(f"BAM reference length: {bam_reference_length:,}")

            starts = []
            ends = []
            strands = []
            cov_np = np.zeros(reference_length, dtype=np.int32)

            for aln in samfile.fetch(
                reference, multiple_iterators=False, until_eof=True
            ):
                ani_read = 1 - ((aln.get_tag("NM") / aln.infer_query_length()))
                if ani_read >= min_read_ani:
                    n_alns += 1
                    read_hits[aln.query_name] += 1

                    # Collect alignment statistics
                    if aln.has_tag("AS"):
                        read_aln_score.append(aln.get_tag("AS"))
                        edit_distances.append(aln.get_tag("NM"))
                        ani_nm.append(ani_read * 100)
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

                # Get all coverage statistics in one pass
                cov_evenness, bases_covered, cov_positions = get_coverage_stats(cov_np)

                # Create cov_pos array only once
                cov_pos = cov_np[cov_positions]

                # Use merge_intervals instead of PyRanges
                merged_starts, merged_ends = merge_intervals(starts, ends)
                ranges = np.array(
                    [end - start for start, end in zip(merged_starts, merged_ends)]
                )

                max_covered_bases = np.max(ranges)
                mean_covered_bases = np.mean(ranges)

                # Calculate coverage metrics
                mean_coverage_trunc, mean_coverage_trunc_len = get_tad(
                    cov_np,
                    trim_min=trim_min,
                    trim_max=trim_max,
                )

                # Coverage calculations using pre-computed values
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
                # cov_evenness = coverage_evenness(cov_np)
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

    # profiler.disable()
    # stats = pstats.Stats(profiler).sort_stats("totti")
    # stats.print_stats(25)

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
