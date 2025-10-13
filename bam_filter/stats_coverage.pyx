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
"""
Coverage-related helpers split out from stats.pyx to make compilation incremental.
"""

from libc.stdint cimport int32_t, int64_t
from libc.stdlib cimport malloc, free, calloc, realloc, qsort
from libc.math cimport sqrt, log2, log, exp, ceil

from bam_filter.stats_rle cimport rle_coverage_t, rle_interval_t
from bam_filter.stats_common cimport ref_stats_t

from bam_filter.stats_helpers cimport compare_pairs, compare_int64

# Minimal externs used by this module; declare locally to avoid depending on
# other .pxd files during incremental modular compilation.
cdef extern from "bam_filter/c_logging.h":
    double bf_monotonic_seconds() nogil
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil

cdef extern from "time.h":
    cdef struct timespec:
        long tv_sec
        long tv_nsec
    int clock_gettime(int clk_id, timespec *tp) nogil
    int CLOCK_MONOTONIC


cdef double compute_coverage_evenness_from_rle(
    int64_t ref_length,
    int32_t* starts, int32_t* ends, int32_t* depths, int64_t n_intervals
) noexcept nogil:
    if ref_length <= 0:
        return 0.0
    cdef int64_t total_coverage = 0
    cdef int64_t total_covered_bases = 0
    cdef int64_t i
    cdef int64_t interval_len = 0
    for i in range(n_intervals):
        interval_len = ends[i] - starts[i]
        total_coverage += depths[i] * interval_len
        total_covered_bases += interval_len
    cdef double mean_coverage = <double>total_coverage / <double>ref_length
    cdef int32_t C = <int32_t>(mean_coverage + (0.5 if mean_coverage >= 0 else -0.5))
    cdef int64_t n_leq_C = 0
    cdef int64_t sum_leq_C = 0
    cdef double depth_val = 0.0
    for i in range(n_intervals):
        depth_val = <double>depths[i]
        if depth_val <= C:
            interval_len = ends[i] - starts[i]
            n_leq_C += interval_len
            sum_leq_C += depths[i] * interval_len
    cdef int64_t uncovered_bases = ref_length - total_covered_bases
    if uncovered_bases > 0 and C >= 0.0:
        n_leq_C += uncovered_bases
    cdef double cov_evenness
    if n_leq_C == 0:
        cov_evenness = 1.0
    else:
        if C > 0.0:
            cov_evenness = 1.0 - (<double>n_leq_C - <double>sum_leq_C / C) / <double>ref_length
        else:
            cov_evenness = 0.0
    return cov_evenness


cdef inline void compute_tad_from_rle(
    int64_t ref_length,
    int32_t* starts, int32_t* ends, int32_t* depths, int64_t n_intervals,
    double* result_mean, int64_t* result_len,
    int trim_min, int trim_max
) noexcept nogil:
    result_mean[0] = 0.0
    result_len[0] = 0
    if ref_length <= 0 or n_intervals == 0:
        return
    cdef int64_t total_covered_positions = 0
    cdef int32_t max_depth = 0
    cdef int64_t i
    for i in range(n_intervals):
        total_covered_positions += ends[i] - starts[i]
        if depths[i] > max_depth:
            max_depth = depths[i]
    if max_depth == 0 or total_covered_positions == 0:
        return
    cdef int64_t* hist = <int64_t*>calloc(max_depth + 1, sizeof(int64_t))
    if hist == NULL:
        return
    for i in range(n_intervals):
        hist[depths[i]] += ends[i] - starts[i]
    cdef int64_t trim_start_pos = (total_covered_positions * trim_min) // 100
    cdef int64_t trim_end_pos = (total_covered_positions * trim_max) // 100
    if trim_end_pos <= trim_start_pos:
        trim_start_pos = 0
        trim_end_pos = total_covered_positions
    cdef int64_t cumsum = 0
    cdef int32_t min_depth = 0, max_depth_thresh = max_depth
    cdef bint found_min = False
    for i in range(max_depth + 1):
        cumsum += hist[i]
        if not found_min and cumsum > trim_start_pos:
            min_depth = <int32_t>i
            found_min = True
        if cumsum >= trim_end_pos:
            max_depth_thresh = <int32_t>i
            break
    cdef int64_t sum_coverage = 0, count_positions = 0
    for i in range(min_depth, max_depth_thresh + 1):
        if hist[i] > 0:
            sum_coverage += i * hist[i]
            count_positions += hist[i]
    if count_positions > 0:
        result_mean[0] = <double>sum_coverage / count_positions
        result_len[0] = count_positions
    free(hist)


cdef void compute_abundance_metrics(ref_stats_t* stats, int64_t n_alns, int64_t n_reads, double read_length_mean, int64_t scale) noexcept nogil:
    if stats == NULL or stats.ref_length <= 0:
        stats.tax_abund_read = 0
        stats.tax_abund_aln = 0
        stats.tax_abund_tad = 0
        stats.n_reads_tad = 0
        return
    cdef double aln_rate = 0.0
    cdef double read_rate = 0.0
    if stats.ref_length > 0:
        aln_rate = <double>n_alns / <double>stats.ref_length
        stats.tax_abund_aln = <int64_t>(aln_rate * scale + (0.5 if aln_rate * scale >= 0 else -0.5))
    else:
        stats.tax_abund_aln = 0
    if stats.ref_length > 0:
        read_rate = <double>n_reads / <double>stats.ref_length
        stats.tax_abund_read = <int64_t>(read_rate * scale + (0.5 if read_rate * scale >= 0 else -0.5))
    else:
        stats.tax_abund_read = 0
    cdef double total_bases_tad = 0.0
    cdef double estimated_reads = 0.0
    cdef double tad_rate = 0.0
    if (stats.mean_coverage_trunc_len > 0 and stats.mean_coverage_trunc > 0 and 
        read_length_mean > 0 and stats.ref_length > 0):
        total_bases_tad = <double>stats.mean_coverage_trunc_len * stats.mean_coverage_trunc
        estimated_reads = total_bases_tad / read_length_mean
        stats.n_reads_tad = <int64_t>(estimated_reads + (0.5 if estimated_reads >= 0 else -0.5))
        tad_rate = estimated_reads / <double>stats.ref_length
        stats.tax_abund_tad = <int64_t>(tad_rate * scale + (0.5 if tad_rate * scale >= 0 else -0.5))
    else:
        stats.n_reads_tad = 0
        stats.tax_abund_tad = 0


# Small helpers ported from stats.pyx: entropy and gini implementations used by coverage stats
cdef inline double compute_entropy(int32_t* counts, int64_t n_bins) except -1 nogil:
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


cdef inline double compute_gini(int32_t* counts, int64_t n_bins) except -1 nogil:
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


cdef void compute_rle_coverage_stats(rle_coverage_t* rle, ref_stats_t* stats, int trim_min, int trim_max) noexcept nogil:
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
    cdef int32_t* starts = NULL
    cdef int32_t* ends = NULL
    cdef int32_t* depths = NULL
    cdef int32_t* hist_counts = NULL
    cdef int64_t* events = NULL
    cdef int64_t* covered_positions = NULL
    cdef int64_t* trimmed_positions = NULL
    cdef int64_t* sorted_counts = NULL
    cdef int64_t* even_counts = NULL
    cdef rle_interval_t* iv
    cdef double mean_interval = 0.0, m2_interval = 0.0, delta = 0.0
    cdef double coverage_variance = 0.0, coverage_sd = 0.0, mean_cov = 0.0
    cdef double sum_squared_deviations = 0.0, depth_val = 0.0, deviation = 0.0, squared_dev = 0.0
    cdef double zero_deviation = 0.0, zero_squared_dev = 0.0
    cdef double bin_width = 0.0, first_edge, last_edge, data_range, iqr, fd_bw, sturges_bw, bin_size
    cdef double entropy_val, p, max_entropy, freq, gini_sum, gini_val, min_gini_sum, min_gini, max_gini
    cdef double cube_root_factor, data_ptp, trimmed_ptp
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
        stats.entropy = 0.0
        stats.norm_entropy = 0.0
        stats.gini = 0.0
        stats.norm_gini = 0.0
        stats.mean_coverage_covered = 0.0
        stats.site_density = 0.0
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
    starts = <int32_t*>malloc(n_events * sizeof(int32_t))
    ends = <int32_t*>malloc(n_events * sizeof(int32_t))
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
                starts[rle_idx] = <int32_t>events[i*2];
                ends[rle_idx] = <int32_t>events[(i+1)*2];
                depths[rle_idx] = <int32_t>curr_depth;
                rle_idx += 1;
            else:
                if current_interval_start != -1:
                    interval_len = events[i*2] - current_interval_start
                    n_intervals += 1
                    delta = (<double>interval_len) - mean_interval
                    mean_interval += delta / n_intervals
                    m2_interval += delta * ((<double>interval_len) - mean_interval)
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
    
    # Calculate coverage standard deviation and variance for c_v and d_i (only covered positions)
    cdef double sum_cov = 0.0
    cdef double sum_cov2 = 0.0
    cdef int64_t n_cov = 0
    cdef double mean_cov_cov = 0.0
    cdef double var_cov = 0.0
    cdef double sd_cov = 0.0
    for i in range(rle_idx):
        depth_val = <double>depths[i]
        interval_len = ends[i] - starts[i]
        if depth_val > 0.0:
            sum_cov += depth_val * interval_len
            sum_cov2 += depth_val * depth_val * interval_len
            n_cov += interval_len
    if n_cov > 1:
        mean_cov_cov = sum_cov / n_cov
        var_cov = (sum_cov2 - (sum_cov * sum_cov) / n_cov) / (n_cov - 1)
        sd_cov = sqrt(var_cov)
    else:
        mean_cov_cov = 0.0
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
    stats.cov_evenness = compute_coverage_evenness_from_rle(
        rle.ref_length, starts, ends, depths, rle_idx
    )
    cov_merge_end = bf_monotonic_seconds()
    
    # Calculate mean coverage for covered positions
    if total_bases_covered > 0:
        stats.mean_coverage_covered = <double>total_coverage / total_bases_covered
    else:
        stats.mean_coverage_covered = 0.0
    
    # NumPy-compatible histogram calculation
    n_covered_positions = 0
    #
    for pos_idx in range(rle_idx):
        #
        n_covered_positions += ends[pos_idx] - starts[pos_idx]
    genome_length = rle.ref_length
    n_bins_calc = 1
    
    if n_covered_positions > 0:
        # Collect positions that have coverage > 0 (matching Python's np.where(cov_np > 0)[0])
        covered_positions = <int64_t*>malloc(n_covered_positions * sizeof(int64_t))
        if covered_positions != NULL:
            pos_count = 0
            
            # Fill with positions that have coverage (same as Python cov_positions)
            for pos_idx in range(rle_idx):
                for interval_pos in range(starts[pos_idx], ends[pos_idx]):
                    if pos_count < n_covered_positions:
                        covered_positions[pos_count] = interval_pos
                        pos_count += 1

            # fprintf(stderr, "Covered positions: ");
            # for pos_idx in range(pos_count):
            #    fprintf(stderr, "%lld ", covered_positions[pos_idx])
            # fprintf(stderr, "\n");

            if pos_count > 0:
                # Apply NumPy's histogram_bin_edges with bins="auto" and range=(0, genome_length)
                # This exactly matches: np.histogram_bin_edges(cov_positions, bins="auto", range=(0, genome_length))
                
                first_edge = 0.0
                last_edge = <double>genome_length
                
                # NumPy's _get_bin_edges logic for bins="auto"
                if pos_count == 0:
                    n_bins_calc = 1
                else:
                    # Sort positions for percentile calculations
                    qsort(covered_positions, pos_count, sizeof(int64_t), compare_int64)
                    
                    # Apply range filtering (NumPy trims data to range before calculating bin width)
                    trimmed_count = 0
                    for pos_idx in range(pos_count):
                        if covered_positions[pos_idx] >= <int64_t>first_edge and covered_positions[pos_idx] <= <int64_t>last_edge:
                            trimmed_count += 1
                    
                    if trimmed_count <= 1:
                        n_bins_calc = 1
                    else:
                        # Create trimmed array for bin width calculation
                        trimmed_positions = <int64_t*>malloc(trimmed_count * sizeof(int64_t))
                        if trimmed_positions != NULL:
                            trim_idx = 0
                            for pos_idx in range(pos_count):
                                if covered_positions[pos_idx] >= <int64_t>first_edge and covered_positions[pos_idx] <= <int64_t>last_edge:
                                    trimmed_positions[trim_idx] = covered_positions[pos_idx]
                                    trim_idx += 1
                            
                            # Calculate bin width using NumPy's _hist_bin_auto logic
                            # FD rule: 2.0 * IQR * n^(-1/3)
                            q1_idx = trimmed_count // 4
                            q3_idx = (3 * trimmed_count) // 4
                            iqr = <double>(trimmed_positions[q3_idx] - trimmed_positions[q1_idx])
                            fd_bw = 2.0 * iqr * exp(-log(<double>trimmed_count) / 3.0)  # n^(-1/3)
                            
                            # Sturges rule: ptp(x) / (log2(n) + 1)  
                            # Note: NumPy uses ptp of TRIMMED data, not full range
                            trimmed_ptp = <double>(trimmed_positions[trimmed_count-1] - trimmed_positions[0])
                            sturges_bw = trimmed_ptp / (log2(<double>trimmed_count) + 1.0)
                            
                            # Auto rule: min(fd_bw, sturges_bw) if fd_bw > 0, else sturges_bw
                            if fd_bw > 0.0:
                                bin_width = fd_bw if fd_bw < sturges_bw else sturges_bw
                            else:
                                bin_width = sturges_bw
                            
                            # Convert bin width to number of bins using SPECIFIED range
                            data_ptp = last_edge - first_edge
                            if bin_width > 0.0:
                                n_bins_calc = <int64_t>ceil(data_ptp / bin_width)
                            else:
                                n_bins_calc = 1
                            
                            free(trimmed_positions)
                        else:
                            # Memory allocation failed - fallback
                            n_bins_calc = <int64_t>(log2(<double>pos_count) + 1.0)
                
                # Sanity check
                # n_bins_calc is no longer capped
                
                # Create histogram counts (matching np.histogram with range=(0, genome_length))
                hist_counts = <int32_t*>calloc(n_bins_calc, sizeof(int32_t))
                if hist_counts != NULL:
                    # Bin the positions using NumPy's logic
                    data_range = last_edge - first_edge  # Initialize data_range here
                    bin_size = data_range / n_bins_calc
                    
                    for pos_idx in range(pos_count):
                        # Use double precision for bin calculation, matching NumPy
                        pos_val = <double>covered_positions[pos_idx]
                        # Bin index: floor((pos_val - first_edge) / bin_size)
                        bin_idx = int((pos_val - first_edge) / bin_size)
                        # NumPy: rightmost edge is inclusive, all others exclusive
                        if pos_val == last_edge:
                            bin_idx = n_bins_calc - 1
                        # Clamp to valid range
                        if bin_idx < 0:
                            continue
                        if bin_idx >= n_bins_calc:
                            bin_idx = n_bins_calc - 1
                        hist_counts[bin_idx] += 1
                    
                    # Calculate entropy and Gini from histogram
                    stats.entropy = compute_entropy(hist_counts, n_bins_calc)
                    stats.norm_entropy = compute_normalized_entropy_inline(hist_counts, n_bins_calc)
                    stats.gini = compute_gini(hist_counts, n_bins_calc)
                    stats.norm_gini = compute_normalized_gini_inline(hist_counts, n_bins_calc)
                    
                    free(hist_counts)
                else:
                    # Memory allocation failed
                    stats.entropy = 0.0
                    stats.norm_entropy = 0.0
                    stats.gini = 0.0
                    stats.norm_gini = 0.0
            else:
                # No positions collected
                stats.entropy = 0.0
                stats.norm_entropy = 0.0
                stats.gini = 0.0
                stats.norm_gini = 0.0
            
            free(covered_positions)
        else:
            # Memory allocation failed - use fallback
            stats.entropy = 0.0
            stats.norm_entropy = 0.0
            stats.gini = 0.0
            stats.norm_gini = 0.0
    else:
        # No covered positions
        stats.entropy = 0.0
        stats.norm_entropy = 0.0
        stats.gini = 0.0
        stats.norm_gini = 0.0
    
    stats.n_bins = n_bins_calc
    
    # Calculate site_density (sites per 1000 bp)
    if genome_length > 0:
        stats.site_density = 1000.0 * <double>total_bases_covered / <double>genome_length
    else:
        stats.site_density = 0.0
    
    # Calculate TAD
    clock_gettime(CLOCK_MONOTONIC, &ts_tad_start)
    compute_tad_from_rle(rle.ref_length, starts, ends, depths, rle_idx, 
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

    # Inline entropy and Gini calculation functions
cdef inline double compute_normalized_entropy_inline(int32_t* counts, int64_t n_bins) noexcept nogil:
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


cdef inline double compute_normalized_gini_inline(int32_t* counts, int64_t n_bins) noexcept nogil:
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
    if max_gini - min_gini == 0.0:
        return 0.0
    
    return (gini_val - min_gini) / (max_gini - min_gini)
