# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
# -*- coding: utf-8 -*-

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free, realloc, calloc, qsort
from libc.math cimport log, log10, exp, sqrt, ceil

from bam_filter.stats cimport RefStats, RLEInterval, RLECoverage
from bam_filter.stats_helpers cimport compare_pairs, compare_int64

# Mathematical constants
cdef double LOG2E = 1.4426950408889634  # log2(e), for converting ln to log2

cdef extern from "time.h":
    cdef struct timespec:
        long tv_sec
        long tv_nsec
    int clock_gettime(int clk_id, timespec *tp) nogil
    int CLOCK_MONOTONIC


cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    double bf_monotonic_seconds() nogil


cdef inline double calculate_coverage_evenness_from_rle(
    int64_t ref_length,
    int64_t* starts, int64_t* ends, int32_t* depths, int64_t n_intervals
) noexcept nogil:
    """
    Calculate coverage evenness from RLE coverage intervals.
    """

    if ref_length <= 0:
        return 0.0

    # Step 1: Calculate total coverage efficiently from RLE
    cdef int64_t total_coverage = 0
    cdef int64_t total_covered_bases = 0
    cdef int64_t i

    cdef int64_t interval_len = 0
    for i in range(n_intervals):
        interval_len = ends[i] - starts[i]
        total_coverage += depths[i] * interval_len
        total_covered_bases += interval_len

    # Calculate mean and round (matches Python's np.rint(np.mean(coverage)))
    cdef double mean_coverage = <double>total_coverage / <double>ref_length
    # Clamp to int32_t range to prevent overflow
    cdef double rounded = mean_coverage + (0.5 if mean_coverage >= 0 else -0.5)
    cdef int32_t C
    if rounded > 2147483647.0:
        C = 2147483647  # INT32_MAX
    elif rounded < -2147483648.0:
        C = -2147483648  # INT32_MIN
    else:
        C = <int32_t>rounded

    # Step 2: Count positions and sum values <= C efficiently from RLE
    cdef int64_t n_leq_C = 0      # len(D2)
    cdef int64_t sum_leq_C = 0    # sum(D2)

    # Process covered intervals
    cdef double depth_val = 0.0
    for i in range(n_intervals):
        depth_val = <double>depths[i]
        if depth_val <= C:
            interval_len = ends[i] - starts[i]
            n_leq_C += interval_len
            sum_leq_C += depths[i] * interval_len

    # Add uncovered regions (depth = 0, always <= C if C >= 0)
    cdef int64_t uncovered_bases = ref_length - total_covered_bases
    if uncovered_bases > 0 and C >= 0.0:
        n_leq_C += uncovered_bases
        # sum_leq_C += 0 * uncovered_bases = 0 (no change)

    # Step 3: Apply Python's exact formula
    cdef double cov_evenness

    if n_leq_C == 0:  # if len(D2) == 0
        cov_evenness = 1.0
    else:
        if C > 0.0:
            # Python: 1.0 - (len(D2) - np.sum(D2) / C) / len(coverage)
            cov_evenness = 1.0 - (<double>n_leq_C - <double>sum_leq_C / C) / <double>ref_length
        else:
            cov_evenness = 0.0

    return cov_evenness


cdef inline double calculate_entropy(int32_t* counts, int64_t n_bins) nogil:
    """Calculate Shannon entropy."""
    cdef double entropy = 0.0
    cdef int64_t total = 0
    cdef int64_t i
    cdef double p

    for i in range(n_bins):
        total += counts[i]

    if total == 0:
        return 0.0

    for i in range(n_bins):
        if counts[i] > 0:
            p = counts[i] / <double>total
            entropy -= p * log(p)

    return entropy


cdef inline double calculate_spatial_entropy(int32_t* counts, int64_t n_bins) nogil:
    """Calculate spatial entropy: measures uniformity of covered position distribution.

    This measures how evenly distributed covered positions are across the reference,
    NOT the uniformity of coverage depths. Optimized for sparse histograms.

    Low spatial entropy = positions clustered (pileup)
    High spatial entropy = positions evenly spread
    """
    cdef double entropy = 0.0
    cdef int64_t total = 0
    cdef int64_t i
    cdef double p

    # Calculate total from non-zero bins only
    for i in range(n_bins):
        total += counts[i]  # This includes zeros but they don't contribute

    if total == 0:
        return 0.0

    # Only process non-zero bins (zeros contribute 0 * log(0) = 0)
    for i in range(n_bins):
        if counts[i] > 0:  # Skip zero bins - mathematically equivalent
            p = counts[i] / <double>total
            entropy -= p * log(p)

    return entropy


cdef inline double calculate_gini(int32_t* counts, int64_t n_bins) nogil:
    """Calculate Gini coefficient."""
    if n_bins <= 1:
        return 0.0

    cdef int64_t* sorted_counts = <int64_t*>malloc(n_bins * sizeof(int64_t))
    if sorted_counts == NULL:
        return 0.0

    cdef int64_t total = 0
    cdef int64_t i

    for i in range(n_bins):
        sorted_counts[i] = counts[i]
        total += counts[i]

    if total == 0:
        free(sorted_counts)
        return 0.0

    qsort(sorted_counts, n_bins, sizeof(int64_t), compare_pairs)

    cdef double gini_sum = 0.0
    for i in range(n_bins):
        gini_sum += (2 * (i + 1) - n_bins - 1) * sorted_counts[i]

    cdef double gini = gini_sum / (n_bins * total)

    free(sorted_counts)
    return gini


cdef inline double calculate_gini_smart(int32_t* counts, int64_t n_bins) nogil:
    """Mathematically identical Gini, optimized for sparse histograms."""
    if n_bins <= 1:
        return 0.0

    cdef int64_t total = 0
    cdef int64_t non_zero_count = 0
    cdef int64_t i, j, rank

    # Count non-zeros and total
    for i in range(n_bins):
        if counts[i] > 0:
            non_zero_count += 1
        total += counts[i]

    if total == 0:
        return 0.0

    cdef int64_t zero_count = n_bins - non_zero_count

    # Create minimal sorted array: [all zeros, then sorted non-zeros]
    cdef int64_t* non_zero_values = <int64_t*>malloc(non_zero_count * sizeof(int64_t))
    if non_zero_values == NULL:
        return 0.0

    # Extract non-zero values
    j = 0
    for i in range(n_bins):
        if counts[i] > 0:
            non_zero_values[j] = counts[i]
            j += 1

    # Sort only the non-zero values (much smaller array)
    qsort(non_zero_values, non_zero_count, sizeof(int64_t), compare_int64)

    # Calculate Gini: zeros come first in sorted order, then non-zeros
    cdef double gini_sum = 0.0

    # Contribution from zeros (rank 1 to zero_count)
    # Each zero contributes: (2*rank - n_bins - 1) * 0 = 0
    # So zeros contribute nothing to gini_sum

    # Contribution from non-zero values (rank zero_count+1 to n_bins)
    for i in range(non_zero_count):
        rank = zero_count + i + 1  # Position in full sorted array
        gini_sum += (2 * rank - n_bins - 1) * non_zero_values[i]

    cdef double gini = gini_sum / (n_bins * total)

    free(non_zero_values)
    return gini


cdef inline void get_tad_from_rle(
    int64_t ref_length,
    int64_t* starts, int64_t* ends, int32_t* depths, int64_t n_intervals,
    double* result_mean, int64_t* result_len,
    int trim_min=10, int trim_max=90
) noexcept nogil:
    """
    Fixed TAD calculation using histogram approach for outlier-resistant depth estimation.
    Trims coverage value distribution (not positional) to maintain stability.
    """
    result_mean[0] = 0.0
    result_len[0] = 0
    if ref_length <= 0 or n_intervals == 0:
        return

    # Calculate total covered positions and find max depth
    cdef int64_t total_covered_positions = 0
    cdef int32_t max_depth = 0
    cdef int64_t i
    for i in range(n_intervals):
        total_covered_positions += ends[i] - starts[i]
        if depths[i] > max_depth:
            max_depth = depths[i]

    if max_depth == 0 or total_covered_positions == 0:
        return

    # Create histogram: depth -> count of positions with that depth
    cdef int64_t* hist = <int64_t*>calloc(max_depth + 1, sizeof(int64_t))
    if hist == NULL:
        return

    # Fill histogram efficiently from RLE
    for i in range(n_intervals):
        hist[depths[i]] += ends[i] - starts[i]

    # CRITICAL FIX: Use total_covered_positions, not ref_length for percentiles
    cdef int64_t trim_start_pos = (total_covered_positions * trim_min) // 100
    cdef int64_t trim_end_pos = (total_covered_positions * trim_max) // 100

    # Ensure we have meaningful trim range
    if trim_end_pos <= trim_start_pos:
        trim_start_pos = 0
        trim_end_pos = total_covered_positions

    # Find depth thresholds using cumulative distribution
    cdef int64_t cumsum = 0
    cdef int32_t min_depth = 0, max_depth_thresh = max_depth
    cdef bint found_min = False

    # Find minimum depth threshold (skip bottom trim_min%)
    for i in range(max_depth + 1):
        cumsum += hist[i]
        if not found_min and cumsum > trim_start_pos:
            min_depth = <int32_t>i
            found_min = True
        if cumsum >= trim_end_pos:
            max_depth_thresh = <int32_t>i
            break

    # Calculate TAD from histogram within trimmed depth range
    cdef int64_t sum_coverage = 0, count_positions = 0
    for i in range(min_depth, max_depth_thresh + 1):
        if hist[i] > 0:
            sum_coverage += i * hist[i]  # depth * number_of_positions_with_that_depth
            count_positions += hist[i]

    if count_positions > 0:
        result_mean[0] = <double>sum_coverage / count_positions
        result_len[0] = count_positions

    free(hist)


cdef RLECoverage* initialize_rle_from_length(int64_t ref_length) noexcept nogil:
    """Initialize RLE coverage structure for a reference of given length."""
    cdef RLECoverage* rle = <RLECoverage*>malloc(sizeof(RLECoverage))
    if rle == NULL:
        return NULL

    rle.ref_length = ref_length
    rle.n_intervals = 0
    rle.capacity = 1000  # Start with reasonable capacity
    rle.intervals = <RLEInterval*>malloc(rle.capacity * sizeof(RLEInterval))

    if rle.intervals == NULL:
        free(rle)
        return NULL

    return rle


cdef void destroy_rle_coverage(RLECoverage* rle) noexcept nogil:
    """Free RLE coverage structure."""
    if rle != NULL:
        if rle.intervals != NULL:
            free(rle.intervals)
        free(rle)


cdef int add_coverage_interval(RLECoverage* rle, int64_t start, int64_t end, int32_t count) noexcept nogil:
    """Add a coverage interval to the RLE structure."""
    cdef int64_t new_capacity  # Declare at function start
    cdef RLEInterval* interval

    if rle == NULL or start >= end:
        return -1

    # Ensure we have capacity
    if rle.n_intervals >= rle.capacity:
        # Check for overflow before doubling (prevent wrapping)
        new_capacity = rle.capacity * 2
        if new_capacity < rle.capacity or new_capacity > 1000000000:
            # Capacity overflow or unreasonably large allocation
            return -1
        rle.capacity = new_capacity
        rle.intervals = <RLEInterval*>realloc(rle.intervals, rle.capacity * sizeof(RLEInterval))
        if rle.intervals == NULL:
            return -1

    # Add the interval (merge later if needed)
    interval = &rle.intervals[rle.n_intervals]
    interval.start = start
    interval.end = end
    interval.count = count
    rle.n_intervals += 1

    return 0


cdef inline double calculate_norm_entropy_inline(int32_t* counts, int64_t n_bins) nogil:
    """Calculate normalized entropy exactly matching Python implementation."""
    cdef int64_t total = 0
    cdef int64_t i
    cdef double entropy_val = 0.0
    cdef double p, max_entropy = 0.0, freq
    cdef int64_t quotient, remainder

    if n_bins <= 1:
        return 1.0

    for i in range(n_bins):
        total += counts[i]

    if total <= 1:
        return 1.0

    # Calculate actual entropy: H = -sum(p * log(p))
    for i in range(n_bins):
        if counts[i] > 0:
            p = counts[i] / <double>total
            entropy_val -= p * log(p)

    # Calculate maximum possible entropy with this distribution
    # Distribute counts as evenly as possible among bins
    quotient = total // n_bins
    remainder = total % n_bins

    # Bins with quotient items
    if quotient > 0:
        freq = quotient / <double>total
        max_entropy -= (n_bins - remainder) * freq * log(freq)

    # Bins with quotient + 1 items
    if remainder > 0 and (quotient + 1) > 0:
        freq = (quotient + 1) / <double>total
        max_entropy -= remainder * freq * log(freq)

    if max_entropy == 0.0:
        return 1.0

    return entropy_val / max_entropy


cdef inline double calculate_normalized_spatial_entropy(int32_t* counts, int64_t n_bins) nogil:
    """Calculate normalized spatial entropy: measures spatial distribution uniformity [0,1].

    Normalizes by the maximum achievable entropy given the constraints (total counts, n_bins),
    NOT the theoretical maximum (log(n_bins)). This accounts for discrete distribution constraints.

    Returns:
        0.0 = Maximally clustered (pileup)
        1.0 = Maximally uniform spatial distribution
    """
    if n_bins <= 1:
        return 1.0

    cdef int64_t total = 0
    cdef double entropy = 0.0
    cdef double p
    cdef int64_t i

    # Calculate total and entropy from non-zero bins only
    for i in range(n_bins):
        if counts[i] > 0:
            total += counts[i]

    if total <= 1:
        return 1.0

    for i in range(n_bins):
        if counts[i] > 0:
            p = counts[i] / <double>total
            entropy -= p * log(p)

    cdef double max_entropy = log(<double>n_bins)
    return entropy / max_entropy if max_entropy > 0.0 else 1.0


cdef inline double calculate_norm_gini_inline(int32_t* counts, int64_t n_bins) nogil:
    """Calculate normalized Gini coefficient exactly matching Python implementation."""
    cdef int64_t total = 0
    cdef int64_t i, j
    cdef double gini_sum = 0.0, gini_val, min_gini_sum = 0.0, min_gini, max_gini
    cdef int64_t quotient, remainder
    cdef int64_t* sorted_counts = NULL
    cdef int64_t* even_counts = NULL

    if n_bins <= 1:
        return 0.0

    for i in range(n_bins):
        total += counts[i]

    if total == 0:
        return 0.0

    # Calculate actual Gini coefficient using standard algorithm
    sorted_counts = <int64_t*>malloc(n_bins * sizeof(int64_t))
    if sorted_counts == NULL:
        return 0.0

    for i in range(n_bins):
        sorted_counts[i] = counts[i]

    # Sort counts for Gini calculation
    qsort(sorted_counts, n_bins, sizeof(int64_t), compare_int64)

    # Calculate Gini using standard formula
    for i in range(n_bins):
        gini_sum += (2 * (i + 1) - n_bins - 1) * sorted_counts[i]

    gini_val = gini_sum / (n_bins * total)

    # Calculate minimum possible Gini (most even distribution)
    quotient = total // n_bins
    remainder = total % n_bins

    # Create most even distribution
    even_counts = <int64_t*>malloc(n_bins * sizeof(int64_t))
    if even_counts == NULL:
        free(sorted_counts)
        return 0.0

    for i in range(n_bins - remainder):
        even_counts[i] = quotient
    for i in range(n_bins - remainder, n_bins):
        even_counts[i] = quotient + 1

    # Sort even distribution
    qsort(even_counts, n_bins, sizeof(int64_t), compare_int64)

    # Calculate min Gini
    for i in range(n_bins):
        min_gini_sum += (2 * (i + 1) - n_bins - 1) * even_counts[i]

    min_gini = min_gini_sum / (n_bins * total)

    # Calculate maximum possible Gini (most uneven: all in one bin)
    max_gini = <double>(n_bins - 1) / n_bins

    free(sorted_counts)
    free(even_counts)

    # Normalize: (actual - min) / (max - min)
    cdef double denominator = max_gini - min_gini
    if denominator == 0.0:
        # When max == min, all distributions have same Gini
        # Return 0.0 (perfectly even) if gini_val matches, else 0.5 (neutral)
        return 0.0 if abs(gini_val - min_gini) < 1e-10 else 0.5

    return (gini_val - min_gini) / denominator


cdef inline double calculate_norm_gini_smart(int32_t* counts, int64_t n_bins) nogil:
    """Mathematically identical normalized Gini, optimized for sparse data."""
    if n_bins <= 1:
        return 0.0

    # Use the smart Gini calculation
    cdef double actual_gini = calculate_gini_smart(counts, n_bins)

    # For normalization, we need min and max possible Gini
    cdef int64_t total = 0
    for i in range(n_bins):
        total += counts[i]

    if total == 0:
        return 0.0

    # Min Gini: perfectly even distribution
    # Max Gini: all mass in one bin = (n_bins - 1) / n_bins
    cdef double max_gini = <double>(n_bins - 1) / n_bins
    cdef double min_gini = 0.0  # Perfectly even case

    if max_gini <= min_gini:
        return 0.0

    return (actual_gini - min_gini) / (max_gini - min_gini)


# =============================================================================
# Weighted (mass-based) histogram functions for posterior-weighted coverage
# =============================================================================

cdef double calculate_weighted_spatial_entropy(double* mass, int64_t n_bins) noexcept nogil:
    """Calculate normalized spatial entropy from mass-weighted histogram.

    Args:
        mass: Histogram where mass[b] = sum of (segment_length * depth) for bin b
        n_bins: Number of histogram bins

    Returns:
        Normalized entropy in [0, 1] where 1 = perfectly uniform distribution
    """
    if n_bins <= 1:
        return 1.0

    cdef double total = 0.0
    cdef double entropy = 0.0
    cdef double p
    cdef int64_t i

    for i in range(n_bins):
        total += mass[i]

    if total <= 0.0:
        return 1.0

    for i in range(n_bins):
        if mass[i] > 0.0:
            p = mass[i] / total
            entropy -= p * log(p)

    cdef double max_entropy = log(<double>n_bins)
    return entropy / max_entropy if max_entropy > 0.0 else 1.0


cdef double calculate_weighted_gini(double* mass, int64_t n_bins) noexcept nogil:
    """Calculate Gini coefficient from mass-weighted histogram.

    Args:
        mass: Histogram where mass[b] = sum of coverage mass for bin b
        n_bins: Number of histogram bins

    Returns:
        Gini coefficient in [0, 1] where 0 = perfect equality
    """
    if n_bins <= 1:
        return 0.0

    cdef double total = 0.0
    cdef int64_t non_zero_count = 0
    cdef int64_t i, j, rank

    for i in range(n_bins):
        if mass[i] > 0.0:
            non_zero_count += 1
        total += mass[i]

    if total <= 0.0:
        return 0.0

    cdef int64_t zero_count = n_bins - non_zero_count

    cdef double* non_zero_values = <double*>malloc(non_zero_count * sizeof(double))
    if non_zero_values == NULL:
        return 0.0

    j = 0
    for i in range(n_bins):
        if mass[i] > 0.0:
            non_zero_values[j] = mass[i]
            j += 1

    qsort(non_zero_values, non_zero_count, sizeof(double), compare_double)

    cdef double gini_sum = 0.0
    for i in range(non_zero_count):
        rank = zero_count + i + 1
        gini_sum += (2 * rank - n_bins - 1) * non_zero_values[i]

    cdef double gini = gini_sum / (n_bins * total)

    free(non_zero_values)
    return gini


cdef int compare_double(const void* a, const void* b) noexcept nogil:
    """Comparison function for qsort on doubles."""
    cdef double va = (<double*>a)[0]
    cdef double vb = (<double*>b)[0]
    if va < vb:
        return -1
    elif va > vb:
        return 1
    return 0


cdef double calculate_norm_weighted_gini(double* mass, int64_t n_bins) noexcept nogil:
    """Calculate normalized Gini coefficient from mass-weighted histogram.

    Normalizes by max possible Gini = (n_bins - 1) / n_bins.

    Args:
        mass: Histogram where mass[b] = sum of coverage mass for bin b
        n_bins: Number of histogram bins

    Returns:
        Normalized Gini in [0, 1] where 0 = perfect equality, 1 = max inequality
    """
    if n_bins <= 1:
        return 0.0

    cdef double actual_gini = calculate_weighted_gini(mass, n_bins)
    cdef double max_gini = <double>(n_bins - 1) / n_bins

    if max_gini <= 0.0:
        return 0.0

    return actual_gini / max_gini


cdef int64_t estimate_histogram_bins_from_length(int64_t ref_length) noexcept nogil:
    """Estimate histogram bins from reference length using simple heuristic.

    Uses sqrt(ref_length / 100) with bounds [10, 1000].
    """
    cdef double bins_f = sqrt(<double>ref_length / 100.0)
    cdef int64_t n_bins = <int64_t>bins_f
    if n_bins < 10:
        n_bins = 10
    if n_bins > 1000:
        n_bins = 1000
    return n_bins


cdef inline int64_t estimate_histogram_bins(
    int64_t n_positions,
    int64_t* interval_starts,
    int64_t* interval_ends,
    int64_t n_intervals,
    double range_span
) noexcept nogil:
    """Estimate optimal number of histogram bins for spatial distribution analysis.

    Uses adaptive binning strategy:
    - Small datasets (n < 100): Sturges' rule (log2(n) + 1)
    - Large datasets (n >= 100): min(Freedman-Diaconis, Sturges)

    Freedman-Diaconis uses IQR for robustness to outliers.
    Sturges is simpler and works well for smaller datasets.

    Args:
        n_positions: Total number of covered positions
        interval_starts: RLE interval start positions
        interval_ends: RLE interval end positions
        n_intervals: Number of RLE intervals
        range_span: Total range to bin over (e.g., genome length)

    Returns:
        Optimal number of bins (>= 1)
    """
    if n_positions <= 1:
        return 1

    cdef double log_count = log(<double>n_positions)
    cdef int64_t q1_idx = n_positions // 4
    cdef int64_t q3_idx = (3 * n_positions) // 4
    cdef double pos_q1 = 0.0, pos_q3 = 0.0
    cdef double pos_first = 0.0, pos_last = 0.0
    cdef int64_t acc = 0, interval_len
    cdef int64_t k
    cdef double iqr, data_range, fd_bw, sturges_bw, bin_width
    cdef int64_t n_bins

    # Find Q1 position by scanning intervals
    for k in range(n_intervals):
        interval_len = interval_ends[k] - interval_starts[k]
        if acc + interval_len > q1_idx:
            pos_q1 = <double>(interval_starts[k] + (q1_idx - acc))
            break
        acc += interval_len

    # Find Q3 position
    acc = 0
    for k in range(n_intervals):
        interval_len = interval_ends[k] - interval_starts[k]
        if acc + interval_len > q3_idx:
            pos_q3 = <double>(interval_starts[k] + (q3_idx - acc))
            break
        acc += interval_len

    iqr = pos_q3 - pos_q1

    # Find first and last covered positions
    for k in range(n_intervals):
        interval_len = interval_ends[k] - interval_starts[k]
        if interval_len > 0:
            pos_first = <double>interval_starts[k]
            break

    for k in range(n_intervals):
        interval_len = interval_ends[k] - interval_starts[k]
        if interval_len > 0:
            pos_last = <double>(interval_ends[k] - 1)

    data_range = pos_last - pos_first

    # Adaptive binning strategy
    if n_positions < 100:
        # Small datasets: use Sturges' rule
        sturges_bw = data_range / (log_count * LOG2E + 1.0) if data_range > 0.0 else 0.0
        bin_width = sturges_bw if sturges_bw > 0.0 else 0.0
    else:
        # Large datasets: use Freedman-Diaconis, fallback to Sturges
        fd_bw = 2.0 * iqr * exp(-log_count / 3.0)
        if fd_bw > 0.0:
            sturges_bw = data_range / (log_count * LOG2E + 1.0) if data_range > 0.0 else 0.0
            bin_width = fd_bw if (sturges_bw <= 0.0 or fd_bw < sturges_bw) else sturges_bw
        else:
            sturges_bw = data_range / (log_count * LOG2E + 1.0) if data_range > 0.0 else 0.0
            bin_width = sturges_bw if sturges_bw > 0.0 else 0.0

    # Convert bin width to bin count
    if bin_width > 0.0:
        n_bins = <int64_t>ceil(range_span / bin_width)
    else:
        # Fallback: use Sturges directly on total range
        n_bins = <int64_t>(log_count * LOG2E + 1.0)

    # Ensure at least 1 bin
    return n_bins if n_bins >= 1 else 1


cdef void calculate_rle_coverage_stats(RLECoverage* rle, RefStats* stats, int trim_min, int trim_max) noexcept nogil:
    """Accurate RLE coverage stats with interval merging and depth calculation."""
    # Instrumentation: measure phases within coverage stats
    cdef double cov_start = bf_monotonic_seconds()
    cdef double cov_events_end, cov_merge_end, cov_tad_end
    # Declare all variables at the top
    cdef timespec ts_covstat_start, ts_covstat_end, ts_tad_start, ts_tad_end
    cdef double elapsed_covstat = 0.0, elapsed_tad = 0.0
    cdef int64_t n_events, i, j, rle_idx, curr_depth, total_coverage, total_bases_covered
    cdef int64_t max_covered_len, cov_len, current_interval_start, n_intervals, interval_len
    cdef int64_t uncovered_bases, n_covered_positions, genome_length, n_bins_calc
    cdef int64_t pos_count, pos_idx, interval_pos, q1_idx, q3_idx, bin_idx, total
    cdef int64_t trimmed_count, trim_idx
    cdef int64_t* starts = NULL
    cdef int64_t* ends = NULL
    cdef int32_t* depths = NULL
    cdef int32_t* hist_counts = NULL
    cdef int64_t* events = NULL
    cdef int64_t* covered_positions = NULL
    cdef int64_t* trimmed_positions = NULL
    cdef int64_t* sorted_counts = NULL
    cdef int64_t* even_counts = NULL
    cdef RLEInterval* iv
    cdef double mean_interval = 0.0, m2_interval = 0.0, delta = 0.0
    cdef double sum_interval_len_sq = 0.0  # For WCB calculation
    cdef double coverage_variance = 0.0, coverage_sd = 0.0, mean_cov = 0.0
    cdef double sum_squared_deviations = 0.0, depth_val = 0.0, deviation = 0.0, squared_dev = 0.0
    cdef double zero_deviation = 0.0, zero_squared_dev = 0.0
    cdef double bin_width = 0.0, first_edge, last_edge, data_range, iqr, fd_bw, sturges_bw, bin_size, bin_size_inv, log_count
    # Timing measurements for histogram pipeline
    cdef double iq_start, iq_end, bw_start, bw_end, fill_start, fill_end
    # Histogram summary variables
    cdef int32_t hist_max
    cdef int32_t hist_min
    cdef int64_t hist_nonzero
    cdef double hist_mean, hist_sd
    cdef int64_t hist_total
    cdef double hist_total_sq
    cdef double var_h
    cdef int32_t v
    cdef double entropy_val, p, max_entropy, freq, gini_sum, gini_val, min_gini_sum, min_gini, max_gini
    cdef double cube_root_factor, data_ptp, trimmed_ptp
    cdef double pos_q1, pos_q3, pos_first, pos_last
    cdef int64_t acc, k
    cdef int64_t iv_start, iv_end, bin_start_pos, bin_end_pos, add_count, full_bin_length, bin_boundary
    cdef double bin_boundary_d
    cdef int start_bin, end_bin, b
    cdef int64_t quotient, remainder

    if rle == NULL or stats == NULL or rle.n_intervals == 0:
        stats.bases_covered = 0
        stats.max_covered_bases = 0
        stats.mean_covered_bases = 0.0
        stats.mean_coverage = 0.0
        stats.breadth = 0.0
        stats.exp_breadth = 0.0
        stats.breadth_exp_ratio = 0.0
        stats.cov_evenness = 0.0
        stats.mean_coverage_trunc = 0.0
        stats.mean_coverage_trunc_len = 0
        stats.n_bins = 0
        stats.spatial_entropy = 0.0
        stats.norm_spatial_entropy = 0.0
        stats.gini = 0.0
        stats.norm_gini = 0.0
        # Ensure timing fields are zeroed when no intervals
        stats.cov_events_sec = 0.0
        stats.cov_merge_sec = 0.0
        stats.cov_tad_sec = 0.0
        stats.cov_total_sec = 0.0
        # Ensure histogram timing fields are zeroed when no intervals
    # Histogram per-phase timing removed; keep histogram summary fields only
        stats.mean_coverage_covered = 0.0
        stats.site_density = 0.0
        # New contamination detection metrics
        stats.n_intervals = 0
        stats.sum_interval_length_sq = 0.0
        stats.weighted_contiguity_breadth = 0.0
        stats.mega_genome_sparsity_index = 0.0
        stats.coverage_compressibility_ratio = 0.0
        stats.feature_space_clustering_score = 0.0
        # Authenticity metrics (p-value computed post-hoc across all refs)
        stats.authenticity_score = 0.0
        stats.authenticity_pvalue = 1.0
        return

    clock_gettime(CLOCK_MONOTONIC, &ts_covstat_start)

    # Gather all events (start: +count, end: -count)
    n_events = rle.n_intervals * 2
    events = <int64_t*>malloc(n_events * 2 * sizeof(int64_t))
    cov_events_end = bf_monotonic_seconds()
    if events == NULL:
        return

    j = 0
    for i in range(rle.n_intervals):
        iv = &rle.intervals[i]
        events[j*2] = iv.start
        events[j*2+1] = iv.count
        j += 1
        events[j*2] = iv.end
        events[j*2+1] = -iv.count
        j += 1

    # Sort events by position
    qsort(events, n_events, 2 * sizeof(int64_t), compare_pairs)

    # Initialize variables
    curr_depth = 0
    total_coverage = 0
    total_bases_covered = 0
    max_covered_len = 0
    current_interval_start = -1
    n_intervals = 0

    # Pre-allocate arrays for RLE TAD calculation and histogram
    starts = <int64_t*>malloc(n_events * sizeof(int64_t))
    ends = <int64_t*>malloc(n_events * sizeof(int64_t))
    depths = <int32_t*>malloc(n_events * sizeof(int32_t))
    if starts == NULL or ends == NULL or depths == NULL:
        free(events)
        if starts != NULL: free(starts)
        if ends != NULL: free(ends)
        if depths != NULL: free(depths)
        return
    rle_idx = 0

    # Process events to build merged intervals
    for i in range(n_events):
        curr_depth += events[i*2+1]
        if i+1 < n_events and events[(i+1)*2] != events[i*2]:
            cov_len = events[(i+1)*2] - events[i*2]
            if curr_depth > 0:
                if current_interval_start == -1:
                    current_interval_start = events[i*2]
                total_coverage += curr_depth * cov_len
                total_bases_covered += cov_len
                # Fill RLE arrays for TAD
                starts[rle_idx] = events[i*2]
                ends[rle_idx] = events[(i+1)*2]
                depths[rle_idx] = <int32_t>curr_depth
                rle_idx += 1;
            else:
                if current_interval_start != -1:
                    interval_len = events[i*2] - current_interval_start
                    n_intervals += 1
                    delta = (<double>interval_len) - mean_interval
                    mean_interval += delta / n_intervals
                    m2_interval += delta * ((<double>interval_len) - mean_interval)
                    sum_interval_len_sq += <double>interval_len * <double>interval_len  # WCB
                    if interval_len > max_covered_len:
                        max_covered_len = interval_len
                    current_interval_start = -1

    # Handle last interval if it ends at the last event
    if current_interval_start != -1:
        interval_len = events[(n_events-1)*2] - current_interval_start
        n_intervals += 1
        delta = (<double>interval_len) - mean_interval
        mean_interval += delta / n_intervals
        m2_interval += delta * ((<double>interval_len) - mean_interval)
        sum_interval_len_sq += <double>interval_len * <double>interval_len  # WCB
        if interval_len > max_covered_len:
            max_covered_len = interval_len

    # Set basic stats
    stats.bases_covered = total_bases_covered
    stats.total_coverage = total_coverage
    stats.max_covered_bases = max_covered_len
    stats.mean_covered_bases = mean_interval if n_intervals > 0 else 0.0
    stats.mean_coverage = <double>total_coverage / rle.ref_length
    stats.breadth = <double>total_bases_covered / rle.ref_length
    stats.exp_breadth = 1.0 - exp(-stats.mean_coverage) if stats.mean_coverage > 0 else 0.0
    stats.breadth_exp_ratio = (
        min(stats.breadth / stats.exp_breadth, 1.0) if stats.exp_breadth > 0 else 0.0
    )

    # Store interval stats and calculate WCB
    stats.n_intervals = n_intervals
    stats.sum_interval_length_sq = sum_interval_len_sq
    # WCB = sum(interval_len^2) / bases_covered^2 - Herfindahl-like concentration index
    # Range: 1/n_intervals (uniform) to 1.0 (single interval)
    # High WCB = few long intervals (good), Low WCB = many scattered fragments (bad)
    cdef double bases_cov_sq = <double>total_bases_covered * <double>total_bases_covered
    stats.weighted_contiguity_breadth = sum_interval_len_sq / bases_cov_sq if bases_cov_sq > 0.0 else 0.0

    # MGSI = Mega-Genome Sparsity Index
    # Compares observed breadth to expected breadth under Poisson model
    # Expected breadth: E[b] = 1 - exp(-n_reads * read_len / ref_length)
    # MGSI = log10(E[b] / b) - high values indicate suspiciously sparse coverage
    cdef double expected_breadth_mgsi = 0.0
    cdef double read_len_estimate = stats.read_length_mean if stats.read_length_mean > 0 else 50.0
    cdef double expected_cov = 0.0
    if rle.ref_length > 0 and stats.n_reads > 0:
        expected_cov = (<double>stats.n_reads * read_len_estimate) / <double>rle.ref_length
        expected_breadth_mgsi = 1.0 - exp(-expected_cov)
    if stats.breadth > 0.0 and expected_breadth_mgsi > 0.0:
        stats.mega_genome_sparsity_index = log10(expected_breadth_mgsi / stats.breadth)
    else:
        stats.mega_genome_sparsity_index = 0.0

    # CCR = Coverage Compressibility Ratio
    # Measures fragmentation: n_intervals / (bases_covered / read_len)
    # High CCR = many tiny scattered islands per "read's worth" of coverage = noise
    # Low CCR = contiguous coverage = real signal
    cdef double reads_worth = total_bases_covered / read_len_estimate if read_len_estimate > 0 else 1.0
    if reads_worth > 0.0:
        stats.coverage_compressibility_ratio = <double>n_intervals / reads_worth
    else:
        stats.coverage_compressibility_ratio = 0.0

    # FSCS = Feature-Space Clustering Score
    # Measures how tightly reads cluster in feature space (GC, complexity)
    # HIGH FSCS = tight clustering in feature space = reads hitting conserved niche = NOISE
    # LOW FSCS = reads sampling diverse genome regions = REAL signal
    # Formula: FSCS = 1 / (1 + combined_cv) where cv = coefficient of variation
    cdef double gc_cv = 0.0
    cdef double dust_cv = 0.0
    cdef double combined_cv = 0.0
    if stats.read_gc_content_mean > 0.0:
        gc_cv = stats.read_gc_content_std / stats.read_gc_content_mean
    if stats.dust_mean > 0.0:
        dust_cv = stats.dust_std / stats.dust_mean
    combined_cv = sqrt(gc_cv * gc_cv + dust_cv * dust_cv)
    stats.feature_space_clustering_score = 1.0 / (1.0 + combined_cv)

    # Calculate coverage standard deviation and variance for c_v and d_i (only covered positions)
    # Welford online variance (weighted by segment length) — avoids catastrophic cancellation.
    cdef int64_t n_cov = 0
    cdef double wf_mean_cov = 0.0
    cdef double wf_M2_cov = 0.0
    cdef double var_cov = 0.0
    cdef double sd_cov = 0.0
    cdef double wf_delta = 0.0
    cdef int64_t interval_len_cov = 0
    for i in range(rle_idx):
        depth_val = <double>depths[i]
        interval_len_cov = ends[i] - starts[i]
        if depth_val > 0.0:
            n_cov += interval_len_cov
            wf_delta = depth_val - wf_mean_cov
            wf_mean_cov += wf_delta * interval_len_cov / n_cov
            wf_M2_cov += wf_delta * (depth_val - wf_mean_cov) * interval_len_cov
    if n_cov > 1:
        var_cov = wf_M2_cov / (n_cov - 1)
        sd_cov = sqrt(var_cov)
    else:
        var_cov = 0.0
        sd_cov = 0.0
    # Calculate c_v and d_i using mean_coverage (per-base, including zeros)
    if stats.mean_coverage > 0.0:
        stats.c_v = sd_cov / stats.mean_coverage
        stats.d_i = var_cov / stats.mean_coverage
    else:
        stats.c_v = 0.0
        stats.d_i = 0.0

    # Calculate coverage evenness
    stats.cov_evenness = calculate_coverage_evenness_from_rle(
        rle.ref_length, starts, ends, depths, rle_idx
    )
    cov_merge_end = bf_monotonic_seconds()

    # Calculate mean coverage for covered positions
    if total_bases_covered > 0:
        stats.mean_coverage_covered = <double>total_coverage / total_bases_covered
    else:
        stats.mean_coverage_covered = 0.0

    # NumPy-compatible histogram calculation (interval-driven, cache-friendly)
    # Use merged intervals directly to compute percentile positions and build
    # histogram counts without materializing every covered position. This
    # reduces memory traffic and avoids random memory accesses that cause
    # cache misses for large genomes/coverage.
    n_covered_positions = 0
    for pos_idx in range(rle_idx):
        n_covered_positions += ends[pos_idx] - starts[pos_idx]

    genome_length = rle.ref_length
    n_bins_calc = 1

    if n_covered_positions > 0:
        # Range for histogram
        first_edge = 0.0
        last_edge = <double>genome_length

        # Trimmed count (positions already lie within [0, genome_length])
        trimmed_count = n_covered_positions

        # Estimate optimal number of bins using adaptive strategy
        n_bins_calc = estimate_histogram_bins(
            trimmed_count,
            starts,
            ends,
            rle_idx,
            last_edge - first_edge
        )

        # Build histogram with performance optimizations
        hist_counts = <int32_t*>calloc(n_bins_calc, sizeof(int32_t))
        if hist_counts != NULL:
            data_range = last_edge - first_edge
            # Guard both divisions: check n_bins_calc > 0 AND data_range > 0
            if n_bins_calc > 0 and data_range > 0.0:
                bin_size = data_range / n_bins_calc
                # Precompute inverse for multiplication instead of division
                bin_size_inv = n_bins_calc / data_range
            else:
                # Degenerate case: use safe defaults
                bin_size = max(data_range, 1.0)
                bin_size_inv = 0.0

            # Time the histogram filling
            fill_start = bf_monotonic_seconds()

            for i in range(rle_idx):
                iv_start = starts[i]
                iv_end = ends[i]
                if iv_end <= iv_start:
                    continue

                # Use multiplication instead of division for bin calculation
                start_bin = <int>(<double>iv_start * bin_size_inv)
                end_bin = <int>(<double>(iv_end - 1) * bin_size_inv)

                # Bounds checking
                if start_bin < 0:
                    start_bin = 0
                elif start_bin >= n_bins_calc:
                    start_bin = n_bins_calc - 1

                if end_bin < 0:
                    end_bin = 0
                elif end_bin >= n_bins_calc:
                    end_bin = n_bins_calc - 1

                if start_bin == end_bin:
                    # Most common case for sparse data - single bin
                    hist_counts[start_bin] += <int32_t>(iv_end - iv_start)
                else:
                    # Simplified multi-bin calculation
                    if bin_size_inv > 0.0:
                        # First partial bin: compute boundary as double then cast
                        bin_boundary_d = (<double>(start_bin + 1)) / bin_size_inv
                        bin_boundary = <int64_t>bin_boundary_d
                        add_count = bin_boundary - iv_start
                        if add_count > 0:
                            hist_counts[start_bin] += <int32_t>add_count

                        # Full middle bins - simplified calculation using while loop to avoid Python range
                        full_bin_length = <int64_t>(1.0 / bin_size_inv)
                        b = start_bin + 1
                        while b < end_bin:
                            hist_counts[b] += <int32_t>full_bin_length
                            b += 1

                        # Last partial bin
                        bin_boundary_d = (<double>end_bin) / bin_size_inv
                        bin_boundary = <int64_t>bin_boundary_d
                        add_count = iv_end - bin_boundary
                        if add_count > 0:
                            hist_counts[end_bin] += <int32_t>add_count
                    else:
                        # Fallback: if bin_size_inv is zero, attribute all to start_bin
                        hist_counts[start_bin] += <int32_t>(iv_end - iv_start)

            fill_end = bf_monotonic_seconds()

            # Compute histogram summary (min, max, mean, sd, nonzero)
            hist_max = -2147483648
            hist_min = 2147483647
            hist_nonzero = 0
            hist_mean = 0.0
            hist_sd = 0.0
            hist_total = 0
            hist_total_sq = 0.0
            for i in range(n_bins_calc):
                v = hist_counts[i]
                if v > 0:
                    hist_nonzero += 1
                if v > hist_max:
                    hist_max = v
                if v < hist_min:
                    hist_min = v
                hist_total += v
                hist_total_sq += (<double>v) * (<double>v)

            if hist_min == 2147483647:
                hist_min = 0

            if n_bins_calc > 0:
                hist_mean = <double>hist_total / <double>n_bins_calc
                if n_bins_calc > 1:
                    var_h = (hist_total_sq - (<double>hist_total * <double>hist_total) / <double>n_bins_calc) / (<double>(n_bins_calc - 1))
                    hist_sd = sqrt(var_h) if var_h > 0.0 else 0.0
                else:
                    hist_sd = 0.0

            # Calculate spatial entropy and Gini (no per-phase timing/logging)
            stats.spatial_entropy = calculate_spatial_entropy(hist_counts, n_bins_calc)
            stats.norm_spatial_entropy = calculate_normalized_spatial_entropy(hist_counts, n_bins_calc)
            stats.gini = calculate_gini_smart(hist_counts, n_bins_calc)
            stats.norm_gini = calculate_norm_gini_smart(hist_counts, n_bins_calc)
            # Cache histogram summary into stats (so caller can aggregate/print summaries later)
            stats.hist_min = hist_min
            stats.hist_max = hist_max
            stats.hist_nonzero = hist_nonzero
            stats.hist_mean = hist_mean
            stats.hist_sd = hist_sd

            free(hist_counts)
        else:
            # Memory allocation failed
            stats.spatial_entropy = 0.0
            stats.norm_spatial_entropy = 0.0
            stats.gini = 0.0
            stats.norm_gini = 0.0
    else:
        # No covered positions
        stats.spatial_entropy = 0.0
        stats.norm_spatial_entropy = 0.0
        stats.gini = 0.0
        stats.norm_gini = 0.0

    stats.n_bins = n_bins_calc

    # Calculate site_density (sites per 1000 bp)
    if genome_length > 0:
        stats.site_density = 1000.0 * <double>total_bases_covered / <double>genome_length
    else:
        stats.site_density = 0.0

    # Authenticity score = norm_spatial_entropy - norm_gini
    # Higher score = more even coverage = more likely authentic
    # P-value is computed post-hoc across all references (requires distribution)
    stats.authenticity_score = stats.norm_spatial_entropy - stats.norm_gini
    stats.authenticity_pvalue = 1.0  # Placeholder, computed after all refs processed

    # Calculate TAD
    clock_gettime(CLOCK_MONOTONIC, &ts_tad_start)
    get_tad_from_rle(rle.ref_length, starts, ends, depths, rle_idx,
                     &stats.mean_coverage_trunc, &stats.mean_coverage_trunc_len,
                     trim_min, trim_max)
    clock_gettime(CLOCK_MONOTONIC, &ts_tad_end)
    elapsed_tad = <double>(ts_tad_end.tv_sec - ts_tad_start.tv_sec) + <double>(ts_tad_end.tv_nsec - ts_tad_start.tv_nsec) / 1e9
    cov_tad_end = bf_monotonic_seconds()

    clock_gettime(CLOCK_MONOTONIC, &ts_covstat_end)
    elapsed_covstat = <double>(ts_covstat_end.tv_sec - ts_covstat_start.tv_sec) + <double>(ts_covstat_end.tv_nsec - ts_covstat_start.tv_nsec) / 1e9

    # Cleanup
    free(events)
    free(starts)
    free(ends)
    free(depths)
    # Verbose-only coverage stats phase durations (nogil-safe, level 2)
    bf_nogil_logf_verbose(2, NULL, "coverage: allocate+events took %.6f s\n", cov_events_end - cov_start)
    bf_nogil_logf_verbose(2, NULL, "coverage: merge+compute took %.6f s\n", cov_merge_end - cov_events_end)
    bf_nogil_logf_verbose(2, NULL, "coverage: TAD took %.6f s\n", cov_tad_end - cov_merge_end)
    # Store per-coverage-phase timings in stats for later aggregation/printing
    stats.cov_events_sec = cov_events_end - cov_start
    stats.cov_merge_sec = cov_merge_end - cov_events_end
    stats.cov_tad_sec = cov_tad_end - cov_merge_end
    stats.cov_total_sec = cov_tad_end - cov_start


cdef void calculate_abundance_metrics(
    RefStats* stats,
    int64_t n_alns,
    int64_t n_reads,
    double read_length_mean,
    int64_t scale
) noexcept nogil:
    """
    Calculate taxonomic abundance metrics from coverage and read statistics - FIXED VERSION.

    Parameters:
    -----------
    stats : RefStats*
        Statistics structure to update with abundance metrics
    n_alns : int64_t
        Number of alignments
    n_reads : int64_t
        Number of unique reads
    read_length_mean : double
        Mean read length
    scale : int64_t
        Scale factor (typically 1,000,000 for reads per million)
    """

    # Validate inputs
    if stats == NULL or stats.ref_length <= 0:
        stats.tax_abund_read = 0
        stats.tax_abund_aln = 0
        stats.tax_abund_tad = 0
        stats.n_reads_tad = 0
        return

    # Calculate tax_abund_aln: alignments per scale factor per reference length
    cdef double aln_rate = 0.0
    cdef double read_rate = 0.0
    if stats.ref_length > 0:
        aln_rate = <double>n_alns / <double>stats.ref_length
        stats.tax_abund_aln = <int64_t>(aln_rate * scale + (0.5 if aln_rate * scale >= 0 else -0.5))
    else:
        stats.tax_abund_aln = 0

    # Calculate tax_abund_read: unique reads per scale factor per reference length
    if stats.ref_length > 0:
        read_rate = <double>n_reads / <double>stats.ref_length
        stats.tax_abund_read = <int64_t>(read_rate * scale + (0.5 if read_rate * scale >= 0 else -0.5))
    else:
        stats.tax_abund_read = 0

    # FIXED TAD abundance calculation
    # Calculate n_reads_tad and tax_abund_tad using truncated coverage
    cdef double total_bases_tad = 0.0
    cdef double estimated_reads = 0.0
    cdef double tad_rate = 0.0

    if (stats.mean_coverage_trunc_len > 0 and stats.mean_coverage_trunc > 0 and
        read_length_mean > 0 and stats.ref_length > 0):

        # CONSERVATIVE APPROACH: Estimate reads based only on TAD regions
        # Total sequenced bases in regions that contributed to TAD calculation
        total_bases_tad = <double>stats.mean_coverage_trunc_len * stats.mean_coverage_trunc
        estimated_reads = total_bases_tad / read_length_mean
        stats.n_reads_tad = <int64_t>(estimated_reads + (0.5 if estimated_reads >= 0 else -0.5))

        # Calculate abundance per reference length (consistent with other abundance metrics)
        # This gives: "reads estimated from truncated depth per reference length"
        tad_rate = estimated_reads / <double>stats.ref_length
        stats.tax_abund_tad = <int64_t>(tad_rate * scale + (0.5 if tad_rate * scale >= 0 else -0.5))
    else:
        stats.n_reads_tad = 0
        stats.tax_abund_tad = 0
