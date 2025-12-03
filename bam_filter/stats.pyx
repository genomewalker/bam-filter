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

"""Reference statistics and filtering module.

Calculates comprehensive per-reference coverage statistics and applies
filtering criteria based on coverage evenness and quality metrics.
"""

#include <inttypes.h>

# Forward declare htslib types to ensure they're available for function pointers
# NOTE: htslib 1.3.x uses bam_hdr_t, newer versions (1.10+) use sam_hdr_t
# We use sam_hdr_t throughout the codebase, so add compatibility typedef and function wrappers
cdef extern from *:
    """
    #include "htslib/sam.h"

    /* Check htslib version and provide compatibility layer for old versions */
    #if !defined(HTS_VERSION) || HTS_VERSION < 101000
        /* Compatibility for htslib < 1.10 (uses bam_hdr_t instead of sam_hdr_t) */
        #ifndef sam_hdr_t
            typedef bam_hdr_t sam_hdr_t;
            #define sam_hdr_destroy bam_hdr_destroy
            #define sam_hdr_name2tid bam_name2id

            /* sam_hdr_tid2len doesn't exist, create inline wrapper */
            static inline uint32_t sam_hdr_tid2len(const sam_hdr_t *h, int tid) {
                return h->target_len[tid];
            }
        #endif
    #endif
    """
    ctypedef struct sam_hdr_t
    ctypedef struct samFile
    ctypedef struct hts_idx_t

    # Forward declare functions
    void sam_hdr_destroy(sam_hdr_t* h) nogil
    int64_t sam_hdr_tid2len(const sam_hdr_t* header, int tid) nogil
    int sam_hdr_name2tid(sam_hdr_t* header, const char* name) nogil

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint16_t, uint8_t, uint32_t, uint64_t, INT32_MAX
from libc.stdlib cimport malloc, free, realloc, calloc, qsort
from libc.string cimport memcpy, memset, strlen, strcmp, strcpy, strtok, strncmp
from libc.stdio cimport fprintf, FILE, fopen, fclose, stderr, snprintf
from libc.stdint cimport int64_t
from libc.math cimport sqrt, log2, log, exp, fabs, ceil, pow

from cython.parallel import prange
from cython cimport boundscheck, wraparound

from bam_filter.batch_utils cimport create_balanced_batches_greedy
from bam_filter.stats_io cimport write_output_files_complete

from bam_filter.reference_lengths cimport (
    TSVReferenceMap,
    load_tsv_reference_file,
    lookup_tsv_reference_length,
    free_tsv_reference_map,
    get_tsv_reference_count,
    print_tsv_reference_stats
)
from bam_filter.stats_bam_writer cimport (
    ReferenceFilter,
    create_reference_filter,
    destroy_reference_filter,
    write_filtered_bam_streaming,
)
from bam_filter.generic_filters cimport (
    GenericFilters,
    create_generic_filters,
    destroy_generic_filters,
    add_filter,
)
from bam_filter.stats_rle cimport (
    initialize_rle_from_length,
    destroy_rle_coverage,
    add_coverage_interval,
    calculate_rle_coverage_stats,
    calculate_abundance_metrics,
)
from bam_filter.stats_helpers cimport (
    count_gc_bases,
    count_reference_gc_bases,
    get_query_alignment_length,
    fnv1a_hash_read_id,
    compare_int32,
    compare_double,
    mode_from_sorted,
    calculate_dust_score,
)

from typing import Any
from bam_filter import logging as bf_logging

# Define timespec struct and clock_gettime manually for nogil timing
cdef extern from "time.h":
    cdef struct timespec:
        long tv_sec
        long tv_nsec
    int clock_gettime(int clk_id, timespec *tp) nogil
    int CLOCK_MONOTONIC

# nogil logging helpers
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf(const char* tag, const char* fmt, ...) nogil
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    void bf_nogil_log_fmt(const char* tag, const char* fmt, ...) nogil
    int bf_set_verbosity(int v) nogil
    int bf_get_verbosity() nogil
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    void bf_log_step_duration_notime(const char* tag, const char* step, double start, double end) nogil
    double bf_monotonic_seconds() nogil

cdef const char* STATS_TAG = b"STATS"


def _announce_stage(title: str, detail: str = "") -> None:
    """Emit a stage heading consistent with the processor pipeline output."""
    bf_logging.summary("")
    bf_logging.summary("┌─ %s", title)
    if detail:
        bf_logging.summary("│ %s", detail)
    bf_logging.summary("└─────────────────────────────────────────────────────────────")


def _stage_duration(stage: str, start_time: float) -> None:
    """Record a stage duration using the human-friendly summary channel."""
    duration = bf_monotonic_seconds() - start_time
    bf_logging.summary("stage=%s duration=%.2fs", stage, duration)


def _format_path(value: Any) -> str:
    """Render optional CLI paths for human-readable logging."""
    if value is None:
        return "none"
    if isinstance(value, (bytes, bytearray)):
        try:
            return value.decode("utf-8")
        except Exception:
            return "<binary-path>"
    return str(value)

cdef extern from "math.h":
    double log2(double x) nogil
    double ceil(double x) nogil

# External C library declarations
cdef extern from "seqid_khash.h":
    ctypedef long long khint64_t
    ctypedef int khint_t
    ctypedef khint_t khiter_t
    ctypedef struct kh_seqid_map_t:
        pass
    kh_seqid_map_t* kh_init_seqid_map() nogil
    void kh_destroy_seqid_map(kh_seqid_map_t*) nogil
    khint_t kh_put_seqid_map(kh_seqid_map_t*, khint64_t, int*) nogil
    khint_t kh_size(kh_seqid_map_t*) nogil

cdef extern from "taxonomy_khash.h":
    ctypedef struct kh_str_t:
        pass
    khint_t kh_get_str(const kh_str_t* h, const char* key) nogil
    khint_t kh_end(const kh_str_t* h) nogil

cdef extern from "lca_stats_khash.h":
    # Hash-based read hash → taxid hash (uint64_t → int32_t) - MEMORY EFFICIENT!
    ctypedef struct kh_read_hash_to_taxid_t:
        pass
    khint_t kh_get_read_hash_to_taxid(const kh_read_hash_to_taxid_t* h, uint64_t key) nogil
    khint_t kh_end_read_hash_to_taxid "kh_end" (const kh_read_hash_to_taxid_t* h) nogil

# Import htslib types and functions from centralized header
from bam_filter.processor_types cimport (
    bam1_t,
    bam1_core_t,
    samFile,
    sam_hdr_t,
    hts_idx_t,
    hts_itr_t,
    hts_open,
    hts_close,
    sam_hdr_read,
    sam_hdr_destroy,
    sam_hdr_tid2len,
    sam_hdr_name2tid,
    sam_hdr_init,
    sam_hdr_write,
    sam_hdr_add_line,
    sam_hdr_add_lines,
    sam_read1,
    sam_write1,
    bam_init1,
    bam_destroy1,
    sam_index_load,
    hts_idx_destroy,
    hts_idx_get_stat,
    sam_itr_queryi,
    sam_itr_next,
    hts_itr_destroy,
    bam_aux_get,
    bam_get_qname,
    bam_aux2i,
    bam_endpos,
    hts_set_threads,
)

cdef extern from "htslib/hts.h":
    int hts_set_opt(void *fp, int opt, ...) nogil
    cdef int HTS_OPT_CACHE_SIZE
    cdef int HTS_OPT_BLOCK_SIZE

# zlib streaming compression (optional gzip output)
cdef extern from "zlib.h":
    ctypedef void* gzFile
    gzFile gzopen(const char* path, const char* mode) nogil
    int gzclose(gzFile file) nogil
    int gzprintf(gzFile file, const char* format, ...) nogil
    int gzwrite(gzFile file, const void* buf, unsigned int len) nogil


# Fast endswith helper for .gz usable nogil
# Error codes
cdef enum BatchError:
    BATCH_OK = 0
    BATCH_FILE_OPEN_FAILED = 1
    BATCH_ITERATOR_FAILED = 2
    BATCH_ALLOC_FAILED = 3
    BATCH_BAM_INIT_FAILED = 4

# Data structures
cdef struct RefStats:
    # Basic counts
    int64_t n_reads
    int64_t n_alns
    
    # Read statistics
    double read_length_mean
    double read_length_std
    int read_length_median
    int min_read_length
    int max_read_length
    int read_length_mode
    
    # ANI and quality metrics
    double ani_mean
    double ani_std
    double ani_median
    double min_ani
    double max_ani
    # double gc_content  # Removed: only using gc_content_mean and gc_content_var
    double aligned_length_mean
    double aln_score_mean
    double aln_score_std
    double mapq_mean
    double mapq_std
    double edit_dist_mean
    double edit_dist_std
    double read_gc_content_mean         # Mean of per-read GC% (READ sequence)
    double read_gc_content_std          # Std of per-read GC% (READ sequence)
    double read_gc_content_total        # Overall GC content (total GC bases / total read length * 100) (READ)
    double ref_gc_content_mean          # Mean of per-read reference GC%
    double ref_gc_content_std           # Std of per-read reference GC%
    double ref_gc_content_total         # Overall reference GC content (total ref GC / total ref length * 100)
    
    # Coverage statistics
    int64_t bases_covered
    int64_t total_coverage
    double mean_coverage
    double mean_coverage_covered
    double breadth
    double exp_breadth
    double breadth_exp_ratio
    double cov_evenness
    double c_v
    double d_i
    
    # TAD (Truncated Average Depth) stats
    double mean_coverage_trunc
    int64_t mean_coverage_trunc_len
    int64_t n_reads_tad
    int64_t tax_abund_aln
    int64_t tax_abund_read
    int64_t tax_abund_tad
    
    # Coverage distribution stats
    double entropy
    double norm_entropy
    double gini
    double norm_gini
    int64_t n_bins
    double site_density
    # Histogram summary cached from calculate_rle_coverage_stats
    int32_t hist_min
    int32_t hist_max
    int64_t hist_nonzero
    double hist_mean
    double hist_sd

    # Per-reference processing time (seconds)
    double ref_seconds
    
    # Interval merging results
    int64_t max_covered_bases
    double mean_covered_bases
    
    # Reference lengths
    int64_t ref_length
    int64_t bam_ref_length

cdef struct BatchData:
    int64_t batch_id
    int64_t tid_start
    int64_t tid_end
    int error_code
    double batch_seconds
    int64_t n_refs_processed

cdef struct FilterConditions:
    int min_read_count
    
    # Filter values
    double min_avg_read_ani
    double min_expected_breadth_ratio
    double min_breadth
    double min_coverage_evenness
    double max_coeff_var
    double min_coverage_mean
    double min_norm_entropy
    double max_norm_gini
    
    # Filter enable flags - NEW
    bint enable_min_avg_read_ani
    bint enable_min_expected_breadth_ratio
    bint enable_min_breadth
    bint enable_min_coverage_evenness
    bint enable_max_coeff_var
    bint enable_min_coverage_mean
    bint enable_min_norm_entropy
    bint enable_max_norm_gini

cdef struct RefUniqueReads:
    kh_seqid_map_t* unique_reads_map

# RLE coverage data structures
cdef struct RLEInterval:
    int64_t start
    int64_t end
    int32_t count

cdef struct RLECoverage:
    RLEInterval* intervals
    int64_t n_intervals
    int64_t capacity
    int64_t ref_length

## Utility functions
cdef void sort_reference_pairs(int64_t* tids, int64_t* counts, int64_t n) noexcept nogil:
    """Quick sort for (tid, count) pairs by TID for memory locality.

    Parameters
    ----------
    tids : int64_t*
        Array of target IDs.
    counts : int64_t*
        Parallel array of counts associated with each tid.
    n : int64_t
        Number of elements in the arrays.
    """
    if n <= 1:
        return
    
    cdef int64_t pivot = tids[n // 2]
    cdef int64_t i = 0, j = n - 1
    cdef int64_t temp_tid, temp_count
    
    while i <= j:
        while tids[i] < pivot: i += 1
        while tids[j] > pivot: j -= 1
        
        if i <= j:
            # Swap both TID and count
            temp_tid = tids[i]
            temp_count = counts[i]
            tids[i] = tids[j]
            counts[i] = counts[j]
            tids[j] = temp_tid
            counts[j] = temp_count
            i += 1
            j -= 1
    
    if j > 0:
        sort_reference_pairs(tids, counts, j + 1)
    if i < n:
        sort_reference_pairs(&tids[i], &counts[i], n - i)



cdef void initialize_reference_stats(RefStats* stats, int64_t ref_length, int64_t bam_ref_length) noexcept nogil:
    """Initialize a RefStats structure with sensible defaults.

    Parameters
    ----------
    stats : RefStats*
        Pointer to an allocated RefStats structure (will be zeroed).
    ref_length : int64_t
        Reference length (from TSV if available, otherwise BAM-derived).
    bam_ref_length : int64_t
        Reference length read from the BAM header.

    Notes
    -----
    This function zeroes the structure and sets sentinel values for fields
    such as minimum read length so later accumulation is correct.
    """
    memset(stats, 0, sizeof(RefStats))
    stats.ref_length = ref_length
    stats.bam_ref_length = bam_ref_length
    stats.min_read_length = 2147483647
    stats.min_ani = 1.0
    stats.max_ani = 0.0

cdef BatchData* create_batch_data(int64_t batch_id, int64_t tid_start, int64_t tid_end) nogil:
    """Allocate and initialize a BatchData container used for processing groups of refs.

    Parameters
    ----------
    batch_id : int64_t
        Identifier for the created batch.
    tid_start : int64_t
        Inclusive start index of target IDs in this batch.
    tid_end : int64_t
        Exclusive end index of target IDs in this batch.

    Returns
    -------
    BatchData*
        Newly allocated BatchData or NULL on allocation failure.
    """
    cdef BatchData* batch = <BatchData*>malloc(sizeof(BatchData))
    if batch == NULL:
        return NULL

    batch.batch_id = batch_id
    batch.tid_start = tid_start
    batch.tid_end = tid_end
    batch.error_code = BATCH_OK
    return batch

cdef int calculate_reference_stats(
    samFile* htsfile,
    sam_hdr_t* header,
    hts_idx_t* idx,
    int64_t tid,
    int64_t num_alns,
    RefStats* stats,
    kh_seqid_map_t* unique_reads_map,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c,
    int64_t scale,
    int trim_ends,
    int trim_min,
    int trim_max,
    bint verbose,
    void* trusted_reads_hash_int
) nogil:
    """Compute detailed per-reference statistics from BAM alignments.

    This function iterates the alignments mapped to a single reference and
    computes coverage distributions (RLE), read-level summaries (lengths, GC,
    ANI), and per-reference aggregates used by filtering and reporting.

    Parameters
    ----------
    htsfile : samFile*
        Thread-local opened BAM file handle used for iteration.
    header : sam_hdr_t*
        BAM header pointer for resolving reference names/lengths.
    idx : hts_idx_t*
        Preloaded BAM index used for creating iterators.
    tid : int64_t
        Reference (target) ID to process.
    num_alns : int64_t
        Expected number of alignments (from index statistics).
    stats : RefStats*
        Output pointer to fill with computed statistics.
    unique_reads_map : kh_seqid_map_t*
        Hash set used to record unique read names for this reference.
    min_read_ani_c, min_read_length_c, max_read_length_c : numeric
        Filtering thresholds applied to reads before statistics aggregation.
    scale : int64_t
        Scaling factor for abundance calculations.
    trim_ends, trim_min, trim_max : int
        Parameters controlling TAD/coverage trimming.
    verbose : bint
        When True, emit additional progress information via nogil logging.
    trusted_reads_hash : kh_str_t*
        Optional hash set of read names to include; alignments whose read names
        are absent are skipped when provided.

    Returns
    -------
    int
        0 on success, negative error code on failure.
    """
    
    # Initialize iterator and BAM record
    cdef double ref_start = bf_monotonic_seconds()
    cdef hts_itr_t* iter = sam_itr_queryi(idx, tid, 0, 0x7fffffff)
    if iter == NULL:
        return -1
    
    cdef bam1_t* b = bam_init1()
    if b == NULL:
        hts_itr_destroy(iter)
        return -1
    
    # Initialize RLE coverage structure
    cdef RLECoverage* rle_coverage = initialize_rle_from_length(stats.ref_length)
    if rle_coverage == NULL:
        bam_destroy1(b)
        hts_itr_destroy(iter)
        return -1
    
    # Pre-allocate arrays with better initial estimates
    cdef int64_t capacity = max(num_alns, 1000)  # Better initial estimate
    cdef int32_t* read_lengths = <int32_t*>malloc(capacity * sizeof(int32_t))
    cdef double* ani_values = <double*>malloc(capacity * sizeof(double))
    cdef int32_t* new_read_lengths  # For safe realloc
    cdef double* new_ani_values     # For safe realloc

    if read_lengths == NULL or ani_values == NULL:
        if read_lengths != NULL: free(read_lengths)
        if ani_values != NULL: free(ani_values)
        destroy_rle_coverage(rle_coverage)
        bam_destroy1(b)
        hts_itr_destroy(iter)
        return -1
    
    # Counters and accumulators
    cdef int64_t n_alns = 0
    cdef int64_t total_read_length = 0
    cdef int64_t total_read_gc_bases = 0
    cdef int64_t total_ref_gc_bases = 0
    cdef int64_t total_ref_length = 0
    cdef double total_ani = 0.0
    
    # Min/max tracking
    cdef int32_t min_read_length = 2147483647  # INT32_MAX
    cdef int32_t max_read_length = 0
    
    # Welford's algorithm variables (combined for efficiency)
    cdef double qaln_mean = 0.0, qaln_M2 = 0.0
    cdef double as_mean = 0.0, as_M2 = 0.0
    cdef double nm_mean = 0.0, nm_M2 = 0.0
    cdef double mapq_mean = 0.0, mapq_M2 = 0.0
    cdef double ani_mean = 0.0, ani_M2 = 0.0
    cdef double dust_mean_acc = 0.0, dust_M2 = 0.0
    cdef double read_gc_mean = 0.0, read_gc_M2 = 0.0
    cdef double ref_gc_mean = 0.0, ref_gc_M2 = 0.0
    cdef int64_t qaln_count = 0, as_count = 0, nm_count = 0, mapq_count = 0, read_gc_count = 0, ref_gc_count = 0
    
    # Main processing loop
    cdef int ret = 0
    cdef int32_t read_length, read_gc_bases, ref_gc_bases, ref_length, nm_val, mapq_val
    cdef double ani, qaln_len, as_val, read_gc_percent, ref_gc_percent, dust_val
    cdef int64_t dust_count = 0
    cdef double delta  # Reused for all Welford calculations
    cdef int64_t start_pos, end_pos, i  # Add missing 'i' variable
    cdef char* qname
    cdef int64_t key
    cdef int ret_val
    cdef uint8_t* nm_tag
    cdef uint8_t* as_tag
    
    while True:
        ret = sam_itr_next(htsfile, iter, b)
        if ret < 0:
            break
            
        read_length = b.core.l_qseq
        if read_length < min_read_length_c or read_length > max_read_length_c:
            continue
        
        # Parse alignment tags
        # Parse NM tag once
        nm_tag = bam_aux_get(b, b"NM")
        nm_val = bam_aux2i(nm_tag) if nm_tag != NULL else -1
        
        # Calculate ANI immediately
        if nm_val >= 0 and read_length > 0:
            ani = (1.0 - (<double>nm_val / read_length)) * 100.0
        else:
            ani = 0.0
            
        if ani < min_read_ani_c:
            continue
        
        # Parse AS tag once
        as_val = 0.0
        as_tag = bam_aux_get(b, b"AS")
        if as_tag != NULL:
            if as_tag[0] == 102:  # 'f' - float
                memcpy(&as_val, <void*>(as_tag + 1), sizeof(float))
                as_count += 1
            else:  # integer
                as_val = <double>bam_aux2i(as_tag)
                as_count += 1
        
        # Get MAPQ once
        mapq_val = 255 if b.core.qual == 255 else b.core.qual
        
        # Analyze sequence data
        # Calculate READ GC content in one pass
        read_gc_bases = count_gc_bases(b)
        read_gc_percent = (<double>read_gc_bases / read_length) * 100.0

        # Calculate REFERENCE GC content
        ref_gc_bases = count_reference_gc_bases(b, &ref_length)
        if ref_length > 0:
            ref_gc_percent = (<double>ref_gc_bases / ref_length) * 100.0
        else:
            ref_gc_percent = read_gc_percent  # Fallback to read GC if no MD tag
            ref_length = read_length
            ref_gc_bases = read_gc_bases

        # Calculate DUST score (normalized 0-1)
        dust_val = calculate_dust_score(b)
        dust_count += 1
        delta = dust_val - dust_mean_acc
        dust_mean_acc += delta / dust_count
        dust_M2 += delta * (dust_val - dust_mean_acc)
        
        # Calculate query alignment length once
        qaln_len = get_query_alignment_length(b)
        
        # Get coverage positions once
        start_pos = <int64_t>b.core.pos
        end_pos = <int64_t>bam_endpos(b)
        
        # Get read name once for unique tracking
        qname = bam_get_qname(b)

        # Hash the read name for both unique tracking AND trusted filtering
        key = fnv1a_hash_read_id(qname)

        # Check if read is in trusted set (if filter is provided)
        # Memory-efficient hash-based lookup (70% less memory than string-based!)
        if trusted_reads_hash_int != NULL:
            if kh_get_read_hash_to_taxid(<kh_read_hash_to_taxid_t*>trusted_reads_hash_int, <uint64_t>key) == kh_end_read_hash_to_taxid(<kh_read_hash_to_taxid_t*>trusted_reads_hash_int):
                continue

        kh_put_seqid_map(unique_reads_map, key, &ret_val)
        
        # Update statistics
        # Resize arrays if needed (rarely)
        if n_alns >= capacity:
            capacity *= 2
            # Use temporary pointers to preserve originals on realloc failure
            new_read_lengths = <int32_t*>realloc(read_lengths, capacity * sizeof(int32_t))
            new_ani_values = <double*>realloc(ani_values, capacity * sizeof(double))
            if new_read_lengths == NULL or new_ani_values == NULL:
                # Cleanup original pointers on failure
                free(read_lengths)
                free(ani_values)
                destroy_rle_coverage(rle_coverage)
                bam_destroy1(b)
                hts_itr_destroy(iter)
                return -1
            read_lengths = new_read_lengths
            ani_values = new_ani_values
        
        # Store values for later median/mode calculation
        read_lengths[n_alns] = read_length
        ani_values[n_alns] = ani
        
        # Update coverage in single operation
        if end_pos > start_pos:
            add_coverage_interval(rle_coverage, start_pos, end_pos, 1)
        # print rle_coverage for debugging
        # fprintf(stderr, "Coverage updated: start=%lld, end=%lld, depth=1\n", start_pos, end_pos)
        # Update all running statistics in one pass
        total_read_length += read_length
        total_read_gc_bases += read_gc_bases
        total_ref_gc_bases += ref_gc_bases
        total_ref_length += ref_length
        total_ani += ani
        
        # Min/max tracking
        if read_length < min_read_length:
            min_read_length = read_length
        if read_length > max_read_length:
            max_read_length = read_length
        
        # Welford's algorithm updates (vectorized)
        n_alns += 1  # Increment once for all calculations
        
        # Query alignment length
        if qaln_len > 0:
            qaln_count += 1
            delta = qaln_len - qaln_mean
            qaln_mean += delta / qaln_count
            qaln_M2 += delta * (qaln_len - qaln_mean)
        
        # ANI (using pre-calculated value)
        delta = ani - ani_mean
        ani_mean += delta / n_alns
        ani_M2 += delta * (ani - ani_mean)
        
        # Alignment Score (using pre-parsed value)
        if as_count > 0:
            delta = as_val - as_mean
            as_mean += delta / as_count
            as_M2 += delta * (as_val - as_mean)
        
        # Edit distance (using pre-parsed value)
        if nm_val >= 0:
            nm_count += 1
            delta = nm_val - nm_mean
            nm_mean += delta / nm_count
            nm_M2 += delta * (nm_val - nm_mean)
        
        # MAPQ (using pre-parsed value)
        mapq_count += 1
        delta = mapq_val - mapq_mean
        mapq_mean += delta / mapq_count
        mapq_M2 += delta * (mapq_val - mapq_mean)

        # READ GC content (using pre-calculated value)
        read_gc_count += 1
        delta = read_gc_percent - read_gc_mean
        read_gc_mean += delta / read_gc_count
        read_gc_M2 += delta * (read_gc_percent - read_gc_mean)

        # REFERENCE GC content (using pre-calculated value)
        ref_gc_count += 1
        delta = ref_gc_percent - ref_gc_mean
        ref_gc_mean += delta / ref_gc_count
        ref_gc_M2 += delta * (ref_gc_percent - ref_gc_mean)
    
    # Post-processing and final calculations
    
    # Set basic counts
    stats.n_alns = n_alns
    stats.n_reads = kh_size(unique_reads_map)
    
    if n_alns == 0:
        # Handle empty case
        free(read_lengths)
        free(ani_values)
        destroy_rle_coverage(rle_coverage)
        bam_destroy1(b)
        hts_itr_destroy(iter)
        return 0
    
    # Calculate means and standard deviations
    stats.read_length_mean = <double>total_read_length / n_alns
    stats.min_read_length = min_read_length
    stats.max_read_length = max_read_length
    
    # Calculate read length std in single pass over stored values
    cdef double sum_sq_diff = 0.0
    cdef double diff
    for i in range(n_alns):
        diff = read_lengths[i] - stats.read_length_mean
        sum_sq_diff += diff * diff
    if n_alns > 1:
        stats.read_length_std = sqrt(sum_sq_diff / (n_alns - 1))
    else:
        stats.read_length_std = 0.0

    # Set all the Welford-calculated statistics
    stats.aligned_length_mean = qaln_mean
    stats.ani_mean = ani_mean
    if n_alns > 1:
        stats.ani_std = sqrt(ani_M2 / (n_alns - 1))
    else:
        stats.ani_std = 0.0
    stats.aln_score_mean = as_mean
    if as_count > 1:
        stats.aln_score_std = sqrt(as_M2 / (as_count - 1))
    else:
        stats.aln_score_std = 0.0
    stats.edit_dist_mean = nm_mean
    if nm_count > 1:
        stats.edit_dist_std = sqrt(nm_M2 / (nm_count - 1))
    else:
        stats.edit_dist_std = 0.0
    stats.mapq_mean = mapq_mean
    if mapq_count > 1:
        stats.mapq_std = sqrt(mapq_M2 / (mapq_count - 1))
    else:
        stats.mapq_std = 0.0
    stats.read_gc_content_mean = read_gc_mean
    if read_gc_count > 1:
        stats.read_gc_content_std = sqrt(read_gc_M2 / (read_gc_count - 1))
    else:
        stats.read_gc_content_std = 0.0
    if total_read_length > 0:
        stats.read_gc_content_total = (<double>total_read_gc_bases / <double>total_read_length) * 100.0
    else:
        stats.read_gc_content_total = 0.0
    stats.ref_gc_content_mean = ref_gc_mean
    if ref_gc_count > 1:
        stats.ref_gc_content_std = sqrt(ref_gc_M2 / (ref_gc_count - 1))
    else:
        stats.ref_gc_content_std = 0.0
    if total_ref_length > 0:
        stats.ref_gc_content_total = (<double>total_ref_gc_bases / <double>total_ref_length) * 100.0
    else:
        stats.ref_gc_content_total = 0.0
    if dust_count > 0:
        stats.dust_mean = dust_mean_acc
        if dust_count > 1:
            stats.dust_std = sqrt(dust_M2 / (dust_count - 1))
        else:
            stats.dust_std = 0.0
    else:
        stats.dust_mean = 0.0
        stats.dust_std = 0.0
    
    # Calculate median/mode (requires sorting - unavoidable)
    # Fix integer division for array indexing
    qsort(read_lengths, n_alns, sizeof(int32_t), compare_int32)
    cdef int64_t mid_idx = n_alns // 2
    if n_alns % 2 == 1:
        stats.read_length_median = read_lengths[mid_idx]
    else:
        # Prevent int32_t overflow by promoting to int64_t before addition
        stats.read_length_median = <int32_t>(((<int64_t>read_lengths[mid_idx - 1] + <int64_t>read_lengths[mid_idx]) // 2))
    stats.read_length_mode = mode_from_sorted(read_lengths, n_alns)
    
    qsort(ani_values, n_alns, sizeof(double), compare_double)
    if n_alns % 2 == 1:
        stats.ani_median = ani_values[mid_idx]
    else:
        stats.ani_median = (ani_values[mid_idx - 1] + ani_values[mid_idx]) / 2.0
    
    # Calculate coverage statistics
    # Coverage/TAD calculation (timed per-reference)
    cdef double cov_start_ref = bf_monotonic_seconds()
    calculate_rle_coverage_stats(rle_coverage, stats, trim_min, trim_max)
    cdef double cov_end_ref = bf_monotonic_seconds()
    bf_nogil_logf_verbose(2, STATS_TAG, "per-ref: coverage+TAD took %.6f s\n", cov_end_ref - cov_start_ref)
    
    # Calculate abundance metrics
    cdef double abund_start = bf_monotonic_seconds()
    calculate_abundance_metrics(
        stats,
        n_alns,
        kh_size(unique_reads_map),
        stats.read_length_mean,
        scale
    )
    cdef double abund_end = bf_monotonic_seconds()
    bf_nogil_logf_verbose(2, STATS_TAG, "per-ref: abundance took %.6f s\n", abund_end - abund_start)
    
    # Cleanup
    free(read_lengths)
    free(ani_values)
    destroy_rle_coverage(rle_coverage)
    bam_destroy1(b)
    hts_itr_destroy(iter)
    cdef double ref_end = bf_monotonic_seconds()
    bf_nogil_logf_verbose(2, STATS_TAG, "per-ref: total took %.6f s\n", ref_end - ref_start)
    # Store per-reference elapsed time for later reporting (caller may print summaries)
    stats.ref_seconds = ref_end - ref_start
    stats.abundance_seconds = abund_end - abund_start
    
    return 0


cdef int process_batches(
    samFile** thread_files,
    const char* bam_file_c,
    sam_hdr_t* header,
    hts_idx_t* idx,
    int64_t* tids_to_process,
    int64_t* tid_align_counts,
    int64_t* ref_lengths_array,
    int64_t* bam_ref_lengths_array,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int64_t n_batches,
    int c_num_threads,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c,
    int64_t scale,
    int trim_ends,
    int trim_min,
    int trim_max,
    bint verbose,
    bint show_progress,
    const char* output_c,
    const char* filtered_output_c,
    const char* filtered_bam_c,
    GenericFilters* gfilters
) nogil:
    
    # Declare all variables at the beginning
    cdef int n_refs = header.n_targets
    cdef RefStats* global_ref_stats = NULL
    cdef RefUniqueReads* ref_unique_reads = NULL
    cdef BatchData** batches = NULL
    cdef ReferenceFilter* ref_filter = NULL
    cdef int ret = 0
    cdef int i, j, batch_id, file_idx
    cdef double total_batch_seconds
    cdef int64_t total_refs_processed
    cdef int counted_batches
    cdef double avg_batch_sec
    cdef double avg_ref_sec
    cdef double avg_refs_per_batch
    cdef double batches_wall_start
    cdef double batches_wall_end
    cdef double batches_wall_seconds
    cdef double output_stage_start = 0.0
    cdef bint stage_output_started = False
    # Helpers for top-slowest selection and printing
    cdef int64_t* slow_pairs = NULL
    cdef int64_t n_candidates = 0
    cdef int64_t cand_idx = 0
    cdef int64_t top_k = 0
    cdef int64_t pair_idx
    cdef RefStats* rstat = NULL
    
    # Allocate global arrays
    global_ref_stats = <RefStats*>calloc(n_refs, sizeof(RefStats))
    ref_unique_reads = <RefUniqueReads*>calloc(n_refs, sizeof(RefUniqueReads))
    
    if global_ref_stats == NULL or ref_unique_reads == NULL:
        if global_ref_stats != NULL:
            free(global_ref_stats)
        if ref_unique_reads != NULL:
            free(ref_unique_reads)
        return -1
    
    # Initialize structures with TSV support
    for i in range(n_refs):
        initialize_reference_stats(&global_ref_stats[i], 
                      ref_lengths_array[i],
                      bam_ref_lengths_array[i])
        ref_unique_reads[i].unique_reads_map = kh_init_seqid_map()
        if ref_unique_reads[i].unique_reads_map == NULL:
            # Cleanup on failure
            for j in range(i):
                if ref_unique_reads[j].unique_reads_map != NULL:
                    kh_destroy_seqid_map(ref_unique_reads[j].unique_reads_map)
            free(ref_unique_reads)
            free(global_ref_stats)
            return -1
    
    # Create batch data structures
    batches = <BatchData**>malloc(n_batches * sizeof(BatchData*))
    if batches == NULL:
        for i in range(n_refs):
            if ref_unique_reads[i].unique_reads_map != NULL:
                kh_destroy_seqid_map(ref_unique_reads[i].unique_reads_map)
        free(ref_unique_reads)
        free(global_ref_stats)
        return -1
    for batch_id in range(n_batches):
        batches[batch_id] = create_batch_data(batch_id, batch_starts[batch_id], batch_ends[batch_id])
        if batches[batch_id] == NULL:
            for j in range(batch_id):
                if batches[j] != NULL:
                    free(batches[j])
            free(batches)
            for i in range(n_refs):
                if ref_unique_reads[i].unique_reads_map != NULL:
                    kh_destroy_seqid_map(ref_unique_reads[i].unique_reads_map)
            free(ref_unique_reads)
            free(global_ref_stats)
            return -1
    
    # Process batches with thread-specific file handles
    # Start wall-clock timer for batch processing (this measures real elapsed time
    # including parallel overlap). We will use this for the user-facing summary so
    # the reported time matches how long the processing actually took.
    batches_wall_start = bf_monotonic_seconds()

    if c_num_threads == 1:
        for batch_id in range(n_batches):
            ret = process_reference_batch(
                thread_files[0], header, idx, tids_to_process, tid_align_counts, batches[batch_id],
                global_ref_stats, ref_unique_reads,
                min_read_ani_c, min_read_length_c, max_read_length_c,
                scale, trim_ends, trim_min, trim_max,
                verbose,
                NULL  # No trusted reads filter for regular stats
            )
            if ret != 0:
                break
    else:
        # Parallelize batch processing using prange. Each iteration works on a disjoint
        # set of references (batches are non-overlapping), and each thread uses a
        # thread-local file handle indexed by threadid(), so this is safe nogil.
        from cython.parallel import prange, threadid

        for batch_id in prange(n_batches, schedule='guided', num_threads=c_num_threads, nogil=True):
            file_idx = threadid()
            # We intentionally ignore the return value here because process_reference_batch
            # will set batches[batch_id].error_code on failure; we will check them after
            # the parallel region (can't break from inside prange).
            process_reference_batch(
                thread_files[file_idx], header, idx, tids_to_process, tid_align_counts, batches[batch_id],
                global_ref_stats, ref_unique_reads,
                min_read_ani_c, min_read_length_c, max_read_length_c,
                scale, trim_ends, trim_min, trim_max,
                verbose,
                NULL  # No trusted reads filter for regular stats
            )

        # After the parallel region, inspect batch error codes and report first error if any.
        for batch_id in range(n_batches):
            if batches[batch_id].error_code != BATCH_OK:
                ret = -1
                break
        # Stop wall-clock timer here (after parallel region finished)
        batches_wall_end = bf_monotonic_seconds()
        batches_wall_seconds = batches_wall_end - batches_wall_start

        # Aggregate per-batch timings (sum of individual batch durations) for
        # diagnostics but do NOT use the sum as the primary user-facing time
        # because it double-counts work done in parallel.
        total_batch_seconds = 0.0
        total_refs_processed = 0
        counted_batches = 0
        for batch_id in range(n_batches):
            if batches[batch_id] != NULL and batches[batch_id].batch_seconds > 0.0:
                total_batch_seconds += batches[batch_id].batch_seconds
                total_refs_processed += batches[batch_id].n_refs_processed
                counted_batches += 1

        if counted_batches > 0:
            avg_batch_sec = total_batch_seconds / counted_batches
            # Use wall-clock time to compute avg seconds per reference (real elapsed)
            avg_ref_sec = batches_wall_seconds / total_refs_processed if total_refs_processed > 0 else 0.0
            avg_refs_per_batch = <double>total_refs_processed / counted_batches if counted_batches > 0 else 0.0
            # Print richer timing: wall-clock (primary), sum-of-batches (diagnostic), averages
            bf_nogil_logf_notime(STATS_TAG, "Batches: wall-clock %.3f s over %d batches (avg %.3f s/batch)\n",
                                  batches_wall_seconds, counted_batches, batches_wall_seconds / counted_batches)
            # Produce condensed summary: overall (already printed) plus top-N slowest references
            # NOTE: per-ref histogram timing/debugging has been removed; skip top-slowest-by-histogram reporting
    
    # Write TSV output files
    if ret == 0:
        if (output_c != NULL) or (filtered_output_c != NULL) or (filtered_bam_c != NULL):
            stage_output_started = True
            output_stage_start = bf_monotonic_seconds()
            with gil:
                stats_path = output_c if output_c != NULL else None
                filtered_stats_path = filtered_output_c if filtered_output_c != NULL else None
                filtered_bam_path = filtered_bam_c if filtered_bam_c != NULL else None
                _announce_stage("Output", "Writing statistics tables and optional filtered BAM")
                bf_logging.summary("Stats output: %s", _format_path(stats_path))
                bf_logging.summary("Filtered stats: %s", _format_path(filtered_stats_path))
                bf_logging.summary("Filtered BAM: %s", _format_path(filtered_bam_path))

        ret = write_output_files_complete(output_c, filtered_output_c, global_ref_stats,
                                          header, n_refs, NULL, gfilters)

        # Write filtered BAM if requested
        if ret == 0 and filtered_bam_c != NULL:
            bf_nogil_logf_notime(STATS_TAG, "Creating filtered BAM: %s\n", filtered_bam_c)

            # Create reference filter
            ref_filter = create_reference_filter(global_ref_stats, gfilters, n_refs)
            
            if ref_filter == NULL:
                bf_nogil_logf_notime(STATS_TAG, "Failed to create reference filter\n")
                ret = -1
            else:
                bf_nogil_logf_notime(STATS_TAG, "BAM filtering: %d/%d references pass criteria\n",
                              ref_filter.n_filtered_refs, n_refs)

                if ref_filter.n_filtered_refs > 0:
                    # Free large in-memory structures that are no longer needed
                    # before writing the filtered BAM to reduce peak memory usage.
                    # At this point the reference decision map (ref_filter) and
                    # the BAM header are sufficient for streaming the output.
                    if global_ref_stats != NULL:
                        free(global_ref_stats)
                        global_ref_stats = NULL

                    if ref_unique_reads != NULL:
                        # Individual maps should have been destroyed during processing,
                        # so it's safe to free the container array now.
                        free(ref_unique_reads)
                        ref_unique_reads = NULL

                    if batches != NULL:
                        for i in range(n_batches):
                            if batches[i] != NULL:
                                free(batches[i])
                        free(batches)
                        batches = NULL

                    if verbose:
                        bf_nogil_logf_notime(STATS_TAG, "[MEM] Freed global_ref_stats, ref_unique_reads and large temp arrays before BAM writing\n")

                    # Reuse an existing thread-local open samFile handle and the loaded index
                    # to avoid loading the BAM index twice and to reduce peak memory.
                    bf_nogil_logf_notime(STATS_TAG, "Writing filtered BAM to %s\n", filtered_bam_c)
                    ret = write_filtered_bam_streaming(
                        thread_files[0], idx, bam_file_c, filtered_bam_c, header, ref_filter, c_num_threads,
                        min_read_ani_c, min_read_length_c, max_read_length_c
                    )

                    if ret == 0:
                        bf_nogil_logf_notime(STATS_TAG, "Filtered BAM written successfully\n")
                    else:
                        bf_nogil_logf_notime(STATS_TAG, "Failed to write filtered BAM (error code %d)\n", ret)
                else:
                    bf_nogil_logf_notime(STATS_TAG, "No references pass filters; skipping filtered BAM output\n")
                    ret = 0  # Not an error condition
                
                destroy_reference_filter(ref_filter)

        if stage_output_started:
            with gil:
                _stage_duration("Output", output_stage_start)
    
    # Cleanup
    # Safe final cleanup: pointers may have been freed earlier to reduce peak memory.
    if batches != NULL:
        for batch_id in range(n_batches):
            if batches[batch_id] != NULL:
                free(batches[batch_id])
        free(batches)
        batches = NULL
    if ref_unique_reads != NULL:
        for i in range(n_refs):
            if ref_unique_reads[i].unique_reads_map != NULL:
                kh_destroy_seqid_map(ref_unique_reads[i].unique_reads_map)
        free(ref_unique_reads)
        ref_unique_reads = NULL
    if global_ref_stats != NULL:
        free(global_ref_stats)
        global_ref_stats = NULL
    
    return ret
cdef int process_reference_batch(
    samFile* htsfile,
    sam_hdr_t* header,
    hts_idx_t* idx,
    int64_t* tids_to_process,
    int64_t* tid_align_counts,
    BatchData* batch,
    RefStats* global_ref_stats,
    RefUniqueReads* ref_unique_reads,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c,
    int64_t scale,
    int trim_ends,
    int trim_min,
    int trim_max,
    bint verbose,
    void* trusted_reads_hash_int
) nogil:
    """Process a single batch of references."""
    # htsfile is now passed in, already opened for this thread
    
    cdef int64_t tid_idx, tid, num_alns
    cdef int ret = 0
    
    # Process each reference in the batch
    cdef double batch_start = bf_monotonic_seconds()
    cdef int64_t refs_processed = 0
    for tid_idx in range(batch.tid_start, batch.tid_end):
        tid = tids_to_process[tid_idx]
        num_alns = tid_align_counts[tid]
        ret = calculate_reference_stats(
            htsfile, header, idx, tid, num_alns,
            &global_ref_stats[tid],
            ref_unique_reads[tid].unique_reads_map,
            min_read_ani_c, min_read_length_c, max_read_length_c,
            scale, trim_ends, trim_min, trim_max,
            verbose,
            trusted_reads_hash_int
        )

        # Free unique_reads_map for this reference immediately after processing
        if ref_unique_reads[tid].unique_reads_map != NULL:
            kh_destroy_seqid_map(ref_unique_reads[tid].unique_reads_map)
            ref_unique_reads[tid].unique_reads_map = NULL

        if ret != 0:
            batch.error_code = BATCH_ITERATOR_FAILED
            break

        # Count references that were scanned (n_alns may be zero but we invoked processing)
        refs_processed += 1

    cdef double batch_end = bf_monotonic_seconds()
    batch.batch_seconds = batch_end - batch_start
    batch.n_refs_processed = refs_processed

    return ret


# Main Python interface function
# COMPLETE MODIFIED PYTHON INTERFACE FUNCTION - REPLACE YOUR EXISTING ONE

def compute_bam_stats(
    bam_file,
    batch_size_param=100,
    verbose=True,
    num_threads=1,
    show_progress=True,
    min_read_length=0,
    max_read_length=0x7fffffff,
    min_read_ani=0.0,
    min_read_count=1,
    output=None,
    filtered_output=None,
    filtered_bam=None,
    scale=1000000,
    trim_ends=0,
    trim_min=10,
    trim_max=90,
    reference_lengths_tsv=None,
    generic_filters=None,  # New: list of (column_index, min, max) tuples
    verbosity_level=None,
):
    """
    Compute comprehensive statistics for alignments in a BAM file.
    
    Parameters:
    -----------
    bam_file : str
        Path to the BAM file
    batch_size_param : int
        Batch size parameter for processing (default: 100)
    verbose : bool
        Enable verbose output (default: True)
    num_threads : int
        Number of threads to use (default: 1)
    show_progress : bool
        Show progress information (default: True)
    min_read_length : int
        Minimum read length filter (default: 0)
    max_read_length : int
        Maximum read length filter (default: 0x7fffffff)
    min_read_ani : float
        Minimum ANI filter as percentage (default: 90.0)
    min_read_count : int
        Minimum read count per reference (default: 1)
    output : str, optional
        Output file path for all statistics
    filtered_output : str, optional
        Output file path for filtered statistics
    filtered_bam : str, optional
        Output file path for filtered BAM containing alignments from passing references
    scale : int
        Scaling factor for abundance calculations (default: 1000000)
    trim_ends : int
        Number of bases to trim from reference ends (default: 0)
    trim_min : int
        Minimum percentile for TAD calculation (default: 10)
    trim_max : int
        Maximum percentile for TAD calculation (default: 90)
    reference_lengths_tsv : str, optional
        Path to TSV file with reference lengths (tab-separated: reference_name\tlength)
    filter_* : various
        Filter conditions for the filtered output file and filtered BAM
    verbosity_level : int, optional
        Explicit verbosity level propagated from the Python logging facade
    
    Returns:
    --------
    int : Return code (0 for success, negative for errors)
    """
    if verbosity_level is not None:
        bf_set_verbosity(int(verbosity_level))

    cdef double stage_input_start = 0.0
    cdef double stage_intake_start = 0.0
    cdef double stage_processing_start = 0.0
    cdef double stage_cleanup_start = 0.0

    bam_display = _format_path(bam_file)
    stats_display = _format_path(output)
    filtered_display = _format_path(filtered_output)
    filtered_bam_display = _format_path(filtered_bam)

    _announce_stage("Input validation", "Validating inputs, outputs, and filters")
    bf_logging.summary("Input BAM: %s", bam_display)
    bf_logging.summary(
        "Read filters: length %d-%d bp | ANI >= %.2f%%",
        min_read_length,
        max_read_length,
        min_read_ani,
    )
    bf_logging.summary(
        "Output targets: stats=%s | filtered=%s | filtered BAM=%s",
        stats_display,
        filtered_display,
        filtered_bam_display,
    )
    bf_logging.summary(
        "Threads: %d | Requested batches: %d | Scale=%d | Trim percentiles: %d-%d | Trim ends: %d",
        num_threads,
        batch_size_param,
        scale,
        trim_min,
        trim_max,
        trim_ends,
    )

    stage_input_start = bf_monotonic_seconds()

    # Input validation and encoding
    cdef bytes bam_file_bytes
    cdef bytes output_bytes = None
    cdef bytes filtered_output_bytes = None
    cdef bytes filtered_bam_bytes = None  # NEW
    cdef bytes tsv_file_bytes = None
    cdef bytes _ref_tsv_bytes = None
    cdef const char* _ref_tsv_c = NULL
    cdef bytes _tmp_tsv_bytes = None
    cdef bytes _tmp_filt_bam = None
    cdef const char* _tsv_log_c = NULL
    cdef const char* _filt_bam_c = NULL
    cdef bytes _out_bytes = None
    cdef const char* _out_c = NULL
    cdef bytes _fout_bytes = None
    cdef const char* _fout_c = NULL
    cdef bytes _fbam_bytes = None
    cdef const char* _fbam_c = NULL
    
    # TSV integration variables
    cdef TSVReferenceMap* tsv_map = NULL
    cdef const char* tsv_file_path_c = NULL
    
    # Timing: start
    cdef timespec ts_start, ts_end
    cdef timespec ts_global_start
    # Additional timespecs for TSV timing and other stages
    cdef timespec ts_tsv_start, ts_tsv_end
    cdef timespec ts_lookup_start, ts_lookup_end
    cdef timespec ts_scan_start, ts_scan_end
    cdef timespec ts_batch_start, ts_batch_end
    cdef timespec ts_proc_start, ts_proc_end
    cdef double elapsed_total = 0.0
    cdef double elapsed = 0.0
    cdef double tsv_load_sec = 0.0
    cdef double tsv_lookup_sec = 0.0
    cdef double scan_sec = 0.0
    cdef double batch_create_sec = 0.0
    cdef double proc_sec = 0.0
    cdef double elapsed_global = 0.0
    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    # Record a global start time that will not be overwritten by
    # subsequent sub-stage timing. We'll use this for the final
    # summary elapsed time (start -> end).
    clock_gettime(CLOCK_MONOTONIC, &ts_global_start)
    try:
        bam_file_bytes = bam_file.encode('utf-8') if isinstance(bam_file, str) else bam_file
    except (UnicodeError, AttributeError):
        bf_nogil_logf_notime(STATS_TAG, "Invalid BAM file path\n")
        return -1

    cdef const char* bam_file_c = bam_file_bytes

    # Load TSV reference file if provided - STRICT ERROR HANDLING
    if reference_lengths_tsv:
        import os
        if not os.path.exists(reference_lengths_tsv):
            # Prepare C string for logging
            _ref_tsv_bytes = reference_lengths_tsv.encode('utf-8') if isinstance(reference_lengths_tsv, str) else reference_lengths_tsv
            _ref_tsv_c = <const char*> _ref_tsv_bytes
            bf_nogil_logf_notime(STATS_TAG, "TSV file not found: %s\n", _ref_tsv_c)
            return -1
        
        # Keep a bytes object alive and extract C pointer for use in nogil calls
        tsv_file_bytes = reference_lengths_tsv.encode('utf-8') if isinstance(reference_lengths_tsv, str) else reference_lengths_tsv
        tsv_file_path_c = <const char*> tsv_file_bytes

        if verbose:
            _ref_tsv_c = <const char*> tsv_file_bytes
            bf_nogil_logf_verbose(1, STATS_TAG, "[TSV] Loading reference lengths from: %s\n", _ref_tsv_c)

        # Time the TSV load under nogil and use C-level fprintf for timing output
        clock_gettime(CLOCK_MONOTONIC, &ts_tsv_start)
        with nogil:
            tsv_map = load_tsv_reference_file(tsv_file_path_c)
        clock_gettime(CLOCK_MONOTONIC, &ts_tsv_end)
        tsv_load_sec = <double>(ts_tsv_end.tv_sec - ts_tsv_start.tv_sec) + <double>(ts_tsv_end.tv_nsec - ts_tsv_start.tv_nsec) / 1e9
        # Use nogil-safe logger for timing
        bf_nogil_logf_notime(STATS_TAG, "TSV load time: %.6f sec\n", tsv_load_sec)

        if not tsv_map:
            _ref_tsv_c = <const char*> tsv_file_bytes
            bf_nogil_logf_notime(STATS_TAG, "Failed to load TSV file: %s\n", _ref_tsv_c)
            return -1
        
        if verbose:
            bf_nogil_logf_verbose(1, STATS_TAG, "[TSV] Successfully loaded reference length overrides\n")
            # Call printing while holding the GIL to ensure safety
            print_tsv_reference_stats(tsv_map)

    clock_gettime(CLOCK_MONOTONIC, &ts_end)
    elapsed_total = <double>(ts_end.tv_sec - ts_start.tv_sec) + <double>(ts_end.tv_nsec - ts_start.tv_nsec) / 1e9
    if verbose:
        bf_nogil_logf_verbose(1, STATS_TAG, "Input validation + param setup: %.6f sec\n", elapsed)

    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    cdef const char* output_c = NULL
    cdef const char* filtered_output_c = NULL
    cdef const char* filtered_bam_c = NULL  # NEW
    
    if output is not None:
        try:
            output_bytes = output.encode('utf-8') if isinstance(output, str) else output
            output_c = output_bytes
        except (UnicodeError, AttributeError):
            bf_nogil_logf_notime(STATS_TAG, "Invalid output file path\n")
            if tsv_map:
                with nogil:
                    free_tsv_reference_map(tsv_map)
            return -1
    
    apply_filters = False
    if filtered_output is not None:
        try:
            filtered_output_bytes = filtered_output.encode('utf-8') if isinstance(filtered_output, str) else filtered_output
            filtered_output_c = filtered_output_bytes
            apply_filters = True  # Always apply filters if filtered_output is requested
        except (UnicodeError, AttributeError):
            bf_nogil_logf_notime(STATS_TAG, "Invalid filtered output file path\n")
            if tsv_map:
                with nogil:
                    free_tsv_reference_map(tsv_map)
            return -1

    # NEW: Handle filtered_bam parameter
    if filtered_bam is not None:
        try:
            filtered_bam_bytes = filtered_bam.encode('utf-8') if isinstance(filtered_bam, str) else filtered_bam
            filtered_bam_c = filtered_bam_bytes
            apply_filters = True  # Filtered BAM requires filtering to be enabled
            
            # Validate BAM extension
            if not (filtered_bam.endswith('.bam') or filtered_bam.endswith('.bam.gz')):
                bf_nogil_logf_verbose(1, STATS_TAG, "Filtered BAM output should use a .bam or .bam.gz extension\n")
                
        except (UnicodeError, AttributeError):
            bf_nogil_logf_notime(STATS_TAG, "Invalid filtered BAM file path\n")
            if tsv_map:
                with nogil:
                    free_tsv_reference_map(tsv_map)
            return -1

    # Validate that at least one output is specified
    if output_c == NULL and filtered_output_c == NULL and filtered_bam_c == NULL:
        bf_nogil_logf_notime(STATS_TAG, "At least one output target (output, filtered_output, filtered_bam) must be specified\n")
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -1

    # Validate that filters are provided if filtered BAM is requested
    if filtered_bam_c != NULL and (generic_filters is None or len(generic_filters) == 0):
        bf_nogil_logf_notime(STATS_TAG, "filtered_bam output requires at least one filter to be specified\n")
        bf_nogil_logf_notime(STATS_TAG, "  Use --filter 'column:min:max' to specify filters\n")
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -1

    clock_gettime(CLOCK_MONOTONIC, &ts_end)
    elapsed = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9
    if verbose:
        bf_nogil_logf_notime(STATS_TAG, "Output path validation: %.6f sec\n", elapsed)

    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    
    # Parameter validation and conversion
    cdef int c_num_threads = max(1, int(num_threads))
    if verbose:
        c_num_threads = 1  # Force single-threaded mode for timing
    cdef int c_min_read_length = max(0, int(min_read_length))
    cdef int c_max_read_length = min(0x7fffffff, max(c_min_read_length, int(max_read_length)))
    cdef double c_min_read_ani = max(0.0, min(100.0, float(min_read_ani)))
    cdef int64_t c_scale = max(1, int(scale))
    cdef int c_trim_ends = max(0, int(trim_ends))
    cdef int c_trim_min = max(0, min(100, int(trim_min)))
    cdef int c_trim_max = max(c_trim_min, min(100, int(trim_max)))
    cdef int c_ref_min_read_count = max(1, int(min_read_count))  # Minimum read count for including a reference

    # Set up generic filters
    cdef GenericFilters* gfilters = NULL
    cdef int n_filters = 0

    if generic_filters is not None and len(generic_filters) > 0:
        n_filters = len(generic_filters)
        gfilters = create_generic_filters(n_filters)
        if gfilters == NULL:
            raise MemoryError("Failed to allocate generic filters")

        # Add each filter from the list
        for col_idx, min_val, max_val in generic_filters:
            if add_filter(gfilters, col_idx, min_val, max_val) != 0:
                destroy_generic_filters(gfilters)
                raise MemoryError("Failed to add filter")

    clock_gettime(CLOCK_MONOTONIC, &ts_end)
    elapsed = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9
    if verbose:
        bf_nogil_logf_notime(STATS_TAG, "Parameter validation: %.6f sec\n", elapsed)

    _stage_duration("Input validation", stage_input_start)

    _announce_stage("BAM processing", "Opening BAM handles and scanning references")
    stage_intake_start = bf_monotonic_seconds()

    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    
    if verbose:
        bf_nogil_logf_notime(STATS_TAG, "Processing parameters:\n")
        bf_nogil_logf_notime(STATS_TAG, "  Threads: %d\n", c_num_threads)
        bf_nogil_logf_notime(STATS_TAG, "  Read length range: [%d, %d]\n", c_min_read_length, c_max_read_length)
        bf_nogil_logf_notime(STATS_TAG, "  Minimum ANI: %.2f%%\n", c_min_read_ani)
        bf_nogil_logf_notime(STATS_TAG, "  Min read count (reference/stats): %d\n", c_ref_min_read_count)
        bf_nogil_logf_notime(STATS_TAG, "  Scale factor: %lld\n", c_scale)
        bf_nogil_logf_notime(STATS_TAG, "  Trim ends: %d\n", c_trim_ends)
        bf_nogil_logf_notime(STATS_TAG, "  TAD percentiles: [%d, %d]\n", c_trim_min, c_trim_max)
    
    if tsv_map:
        # Use previously prepared tsv_file_bytes if available
        if tsv_file_bytes is not None:
            _tsv_log_c = <const char*> tsv_file_bytes
        else:
            _tmp_tsv_bytes = reference_lengths_tsv.encode('utf-8') if isinstance(reference_lengths_tsv, str) else reference_lengths_tsv
            _tsv_log_c = <const char*> _tmp_tsv_bytes
        bf_nogil_logf_notime(STATS_TAG, "  TSV reference file: %s\n", _tsv_log_c)
    
    if filtered_bam_c != NULL:
        # filtered_bam may be Python str or bytes; log safely
        # Use filtered_bam_bytes prepared earlier if present
        if filtered_bam_bytes is not None:
            _filt_bam_c = <const char*> filtered_bam_bytes
        else:
            _tmp_filt_bam = filtered_bam.encode('utf-8') if isinstance(filtered_bam, str) else filtered_bam
            _filt_bam_c = <const char*> _tmp_filt_bam
        bf_nogil_logf_notime(STATS_TAG, "  Filtered BAM output: %s\n", _filt_bam_c)
    
    # Open and validate BAM file
    bf_nogil_logf_notime(STATS_TAG, "Opening BAM file\n")
    
    # Pre-open BAM file handles for each thread
    cdef double preopen_start = bf_monotonic_seconds()
    bf_nogil_logf_notime(STATS_TAG, "Pre-opening %d BAM file handles\n", c_num_threads)
    cdef samFile** thread_files = <samFile**>malloc(c_num_threads * sizeof(samFile*))
    if thread_files == NULL:
        bf_nogil_logf_notime(STATS_TAG, "Failed to allocate file handle array\n")
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -1
    cdef int thread_idx
    for thread_idx in range(c_num_threads):
        thread_files[thread_idx] = hts_open(bam_file_c, b"r")
        if thread_files[thread_idx] == NULL:
            bf_nogil_logf_notime(STATS_TAG, "Failed to open BAM file handle %d\n", thread_idx)
            for j in range(thread_idx):
                if thread_files[j] != NULL:
                    hts_close(thread_files[j])
            free(thread_files)
            if tsv_map:
                with nogil:
                    free_tsv_reference_map(tsv_map)
            return -1
    # Log duration for pre-opening handles (always visible)
    cdef double preopen_end = bf_monotonic_seconds()
    bf_log_step_duration_notime(STATS_TAG, "Pre-open BAM handles", preopen_start, preopen_end)

    # Use the first handle for header/index reading
    cdef samFile* htsfile = thread_files[0]
    if htsfile == NULL:
        bf_nogil_logf_notime(STATS_TAG, "No valid BAM file handle for header/index\n")
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -1
    
    cdef sam_hdr_t* header = sam_hdr_read(htsfile)

    clock_gettime(CLOCK_MONOTONIC, &ts_end)
    elapsed = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9
    # Always print this timing (use no-time logger to avoid duplicating timestamps)
    bf_nogil_logf_notime(STATS_TAG, "BAM open + header read: %.6f sec\n", elapsed)

    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    if header == NULL:
        bf_nogil_logf_notime(STATS_TAG, "Failed to read BAM header\n")
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -2
    
    cdef hts_idx_t* idx = sam_index_load(htsfile, bam_file_c)
    if idx == NULL:
        bf_nogil_logf_notime(STATS_TAG, "Failed to load BAM index\n")
        sam_hdr_destroy(header)
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -3

    clock_gettime(CLOCK_MONOTONIC, &ts_end)
    elapsed = (ts_end.tv_sec - ts_start.tv_sec) + (ts_end.tv_nsec - ts_start.tv_nsec) / 1e9
    # Always print this timing
    bf_nogil_logf_notime(STATS_TAG, "BAM index load: %.6f sec\n", elapsed)

    # Also print number of references immediately after index load so it's visible
    if header != NULL:
        bf_nogil_logf_notime(STATS_TAG, "BAM references: %d\n", header.n_targets)

    clock_gettime(CLOCK_MONOTONIC, &ts_start)
    
    cdef int64_t n_refs = header.n_targets
    
    bf_nogil_logf_verbose(1, STATS_TAG, "Found %lld references in BAM file\n", n_refs)
    
    if n_refs <= 0:
        bf_nogil_logf_notime(STATS_TAG, "No references found in BAM file\n")
        
        # Perform cleanup
        hts_idx_destroy(idx)
        sam_hdr_destroy(header)
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        
        return -4
    
    # Allocate and initialize arrays
    bf_nogil_logf_verbose(1, STATS_TAG, "Initializing reference arrays\n")
    
    cdef int64_t* tid_align_counts = <int64_t*>calloc(n_refs, sizeof(int64_t))
    cdef int64_t* tids_to_process = <int64_t*>calloc(n_refs, sizeof(int64_t))
    cdef int64_t* ref_lengths_array = <int64_t*>calloc(n_refs, sizeof(int64_t))
    cdef int64_t* bam_ref_lengths_array = <int64_t*>calloc(n_refs, sizeof(int64_t))
    
    if (tid_align_counts == NULL or tids_to_process == NULL or 
        ref_lengths_array == NULL or bam_ref_lengths_array == NULL):
        bf_nogil_logf_notime(STATS_TAG, "Failed to allocate reference arrays\n")
        # Cleanup
        if tid_align_counts != NULL: free(tid_align_counts)
        if tids_to_process != NULL: free(tids_to_process)
        if ref_lengths_array != NULL: free(ref_lengths_array)
        if bam_ref_lengths_array != NULL: free(bam_ref_lengths_array)
        hts_idx_destroy(idx)
        sam_hdr_destroy(header)
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -5
    
    # Populate reference length arrays with TSV support
    cdef int64_t tid_idx
    cdef const char* ref_name
    cdef int64_t tsv_length, bam_length
    cdef int64_t i_proc, proc_tid
    cdef int32_t tsv_used_count = 0, bam_used_count = 0
    cdef int32_t n_tsv_entries
    cdef int32_t i_entry
    cdef const char* entry_name
    cdef int64_t entry_length
    cdef int32_t entry_tid
    cdef bint has_tsv_map = (tsv_map != NULL)

    # First, fill BAM lengths for all references under nogil (fast C loop)
    with nogil:
        for tid_idx in range(n_refs):
            bam_length = sam_hdr_tid2len(header, tid_idx)
            bam_ref_lengths_array[tid_idx] = bam_length
            # default to BAM length for the TSV-backed array for now
            ref_lengths_array[tid_idx] = bam_length if bam_length > 0 else 1
    
    # Scan references for mapped reads
    bf_nogil_logf_verbose(1, STATS_TAG, "Scanning references for mapped reads\n")
    
    cdef int64_t n_tids_to_process = 0
    cdef uint64_t mapped, unmapped
    cdef int64_t total_mapped_reads = 0
    
    clock_gettime(CLOCK_MONOTONIC, &ts_scan_start)
    for tid_idx in range(n_refs):
        if hts_idx_get_stat(idx, tid_idx, &mapped, &unmapped) == 0:
            tid_align_counts[tid_idx] = mapped
            total_mapped_reads += mapped
            # Reference-level filtering
            if tid_align_counts[tid_idx] >= c_ref_min_read_count:
                tids_to_process[n_tids_to_process] = tid_idx
                n_tids_to_process += 1
    clock_gettime(CLOCK_MONOTONIC, &ts_scan_end)
    scan_sec = <double>(ts_scan_end.tv_sec - ts_scan_start.tv_sec) + <double>(ts_scan_end.tv_nsec - ts_scan_start.tv_nsec) / 1e9
    bf_nogil_logf_notime(STATS_TAG, "Reference scan (hts_idx_get_stat) time: %.6f sec\n", scan_sec)
    

    bf_nogil_logf_verbose(1, STATS_TAG, "Large-scale processing configuration:\n")
    bf_nogil_logf_verbose(1, STATS_TAG, "  References: %lld to process (%lld total)\n", n_tids_to_process, n_refs)
    bf_nogil_logf_verbose(1, STATS_TAG, "  Total mapped reads: %lld\n", total_mapped_reads)
    bf_nogil_logf_verbose(1, STATS_TAG, "  Processing threads: %d\n", c_num_threads)

    if n_tids_to_process == 0:
        bf_nogil_logf_notime(STATS_TAG, "No references meet the minimum read count threshold\n")
        free(bam_ref_lengths_array)
        free(ref_lengths_array)
        free(tids_to_process)
        free(tid_align_counts)
        hts_idx_destroy(idx)
        sam_hdr_destroy(header)
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return 0

    # Apply TSV overrides only to references that will be processed
    if has_tsv_map and n_tids_to_process > 0:
        clock_gettime(CLOCK_MONOTONIC, &ts_lookup_start)
        with nogil:
            n_tsv_entries = get_tsv_reference_count(tsv_map)
            for i_entry in range(n_tsv_entries):
                entry_name = tsv_map.entries[i_entry].reference_name
                entry_length = tsv_map.entries[i_entry].reference_length
                if entry_name == NULL:
                    continue
                entry_tid = sam_hdr_name2tid(header, entry_name)
                if entry_tid < 0 or entry_tid >= n_refs:
                    continue
                if tid_align_counts[entry_tid] >= c_ref_min_read_count:
                    if entry_length > 0:
                        ref_lengths_array[entry_tid] = entry_length
                        tsv_used_count += 1
        clock_gettime(CLOCK_MONOTONIC, &ts_lookup_end)
        tsv_lookup_sec = <double>(ts_lookup_end.tv_sec - ts_lookup_start.tv_sec) + <double>(ts_lookup_end.tv_nsec - ts_lookup_start.tv_nsec) / 1e9
        bf_nogil_logf_notime(STATS_TAG, "TSV lookup time for %lld processed refs (scanned %d entries): %.6f sec\n", n_tids_to_process, n_tsv_entries, tsv_lookup_sec)

        bam_used_count = n_tids_to_process - tsv_used_count
        bf_nogil_logf_notime(STATS_TAG, "[TSV] overrides used: %d; BAM defaults used: %d; Processed refs: %lld/%lld\n",
                tsv_used_count, bam_used_count, n_tids_to_process, n_refs)
    else:
        bam_used_count = n_tids_to_process if n_tids_to_process >= 0 else 0

    if verbose and has_tsv_map:
        bf_nogil_logf_notime(STATS_TAG, "[TSV INTEGRATION] Reference lengths:\n")
        bf_nogil_logf_notime(STATS_TAG, "  TSV overrides used: %d\n", tsv_used_count)
        bf_nogil_logf_notime(STATS_TAG, "  BAM defaults used: %d\n", bam_used_count)
        bf_nogil_logf_notime(STATS_TAG, "  Processed references: %lld / %lld\n", n_tids_to_process, n_refs)
    
    # Create batches with load balancing
    cdef int64_t n_batches = min(c_num_threads * 4, n_tids_to_process)
    n_batches = max(1, min(n_batches, n_tids_to_process))
    
    cdef int64_t max_possible_batches = 512

    bf_nogil_logf_notime(STATS_TAG, "Initial batch allocation:\n")
    bf_nogil_logf_notime(STATS_TAG, "  Requested batches: %lld\n", n_batches)
    bf_nogil_logf_notime(STATS_TAG, "  Max possible batches: %lld\n", max_possible_batches)
    
    cdef int64_t* batch_starts = <int64_t*>calloc(max_possible_batches, sizeof(int64_t))
    cdef int64_t* batch_ends = <int64_t*>calloc(max_possible_batches, sizeof(int64_t))
    
    if batch_starts == NULL or batch_ends == NULL:
        bf_nogil_logf_notime(STATS_TAG, "Failed to allocate batch arrays\n")
        if batch_starts != NULL: free(batch_starts)
        if batch_ends != NULL: free(batch_ends)
        free(bam_ref_lengths_array)
        free(ref_lengths_array)
        free(tids_to_process)
        free(tid_align_counts)
        hts_idx_destroy(idx)
        sam_hdr_destroy(header)
        for thread_idx in range(c_num_threads):
            if thread_files[thread_idx] != NULL:
                hts_close(thread_files[thread_idx])
        free(thread_files)
        if tsv_map:
            with nogil:
                free_tsv_reference_map(tsv_map)
        return -7
    
    # Create batches using the corrected function
    cdef int ret_batch = create_balanced_batches_greedy(
        tids_to_process, tid_align_counts, n_tids_to_process, 
        max_possible_batches, batch_starts, batch_ends, num_threads, verbose
    )

    # Count actual batches created
    cdef int64_t actual_batches = 0
    cdef int64_t i
    for i in range(max_possible_batches):
        if batch_starts[i] < batch_ends[i]:
            actual_batches = i + 1
        else:
            break
    
    n_batches = actual_batches
    

    bf_nogil_logf_notime(STATS_TAG, "Final batch configuration:\n")
    bf_nogil_logf_notime(STATS_TAG, "  Created %lld batches\n", n_batches)
    bf_nogil_logf_notime(STATS_TAG, "  Average: %lld refs per batch\n", n_tids_to_process // n_batches)
    
    # Process batches
    cdef int ret
    # Time batch creation reporting
    clock_gettime(CLOCK_MONOTONIC, &ts_batch_start)
    # (No-op: batches already created above) - measure small overhead
    clock_gettime(CLOCK_MONOTONIC, &ts_batch_end)
    batch_create_sec = <double>(ts_batch_end.tv_sec - ts_batch_start.tv_sec) + <double>(ts_batch_end.tv_nsec - ts_batch_start.tv_nsec) / 1e9
    if verbose:
        bf_nogil_logf_notime(STATS_TAG, "Batch creation overhead: %.6f sec\n", batch_create_sec)

    bf_logging.summary(
        "References eligible: %d / %d | Mapped reads: %d",
        int(n_tids_to_process),
        header.n_targets,
        int(total_mapped_reads),
    )
    bf_logging.summary(
        "Batch layout: %d batches | Avg refs/batch: %d | Threads: %d",
        int(n_batches),
        int(n_tids_to_process // n_batches) if n_batches > 0 else 0,
        c_num_threads,
    )
    _stage_duration("BAM processing", stage_intake_start)

    _announce_stage("Statistics calculation", "Computing per-reference coverage and abundance metrics")
    bf_logging.summary("Processing %d batches with %d threads", int(n_batches), c_num_threads)
    stage_processing_start = bf_monotonic_seconds()

    # Time the main processing step
    clock_gettime(CLOCK_MONOTONIC, &ts_proc_start)
    bf_nogil_logf_notime(STATS_TAG, "Starting batch processing\n")
    ret = process_batches(
        thread_files,
        bam_file_c, header, idx, tids_to_process, tid_align_counts,
        ref_lengths_array, bam_ref_lengths_array,
        batch_starts, batch_ends, n_batches, c_num_threads,
        c_min_read_ani, c_min_read_length, c_max_read_length,
        c_scale, c_trim_ends, c_trim_min, c_trim_max,
        verbose, show_progress,
        output_c, filtered_output_c, filtered_bam_c,  # Pass filtered_bam_c
        gfilters
    )
    clock_gettime(CLOCK_MONOTONIC, &ts_proc_end)
    proc_sec = <double>(ts_proc_end.tv_sec - ts_proc_start.tv_sec) + <double>(ts_proc_end.tv_nsec - ts_proc_start.tv_nsec) / 1e9

    bf_logging.summary("Batch processing wall time: %.2fs", proc_sec)
    _stage_duration("Statistics calculation", stage_processing_start)

    _announce_stage("Cleanup", "Releasing temporary buffers and closing file handles")
    stage_cleanup_start = bf_monotonic_seconds()

    # Cleanup
    bf_nogil_logf_notime(STATS_TAG, "Releasing temporary resources\n")

    # Cleanup TSV map if loaded
    if tsv_map:
        with nogil:
            free_tsv_reference_map(tsv_map)

    # Cleanup generic filters
    if gfilters != NULL:
        destroy_generic_filters(gfilters)

    free(batch_starts)
    free(batch_ends)
    free(bam_ref_lengths_array)
    free(ref_lengths_array)
    free(tids_to_process)
    free(tid_align_counts)
    hts_idx_destroy(idx)
    sam_hdr_destroy(header)
    bf_nogil_logf_notime(STATS_TAG, "Closing pre-opened BAM handles\n")
    cdef double close_start = bf_monotonic_seconds()
    for thread_idx in range(c_num_threads):
        if thread_files[thread_idx] != NULL:
            hts_close(thread_files[thread_idx])
    free(thread_files)
    cdef double close_end = bf_monotonic_seconds()
    bf_log_step_duration_notime(STATS_TAG, "Close pre-opened handles", close_start, close_end)

    _stage_duration("Cleanup", stage_cleanup_start)

    # Report results
    if ret == 0:
        # Compute overall elapsed time from the global start timestamp
        clock_gettime(CLOCK_MONOTONIC, &ts_end)
        elapsed_global = <double>(ts_end.tv_sec - ts_global_start.tv_sec) + <double>(ts_end.tv_nsec - ts_global_start.tv_nsec) / 1e9
        # Single concise completion line with overall timing
        bf_nogil_logf_notime(STATS_TAG, "Filtering finished in %.6f s\n", elapsed_global)
        bf_nogil_logf_notime(STATS_TAG, "Reference filtering complete\n")
        if verbose:
            if output_c != NULL:
                _out_bytes = output.encode('utf-8') if isinstance(output, str) else output
                _out_c = <const char*> _out_bytes
                bf_nogil_logf_notime(STATS_TAG, "Summary statistics written to %s\n", _out_c)
            if filtered_output_c != NULL:
                _fout_bytes = filtered_output.encode('utf-8') if isinstance(filtered_output, str) else filtered_output
                _fout_c = <const char*> _fout_bytes
                bf_nogil_logf_notime(STATS_TAG, "Filtered statistics written to %s\n", _fout_c)
            if filtered_bam_c != NULL:
                _fbam_bytes = filtered_bam.encode('utf-8') if isinstance(filtered_bam, str) else filtered_bam
                _fbam_c = <const char*> _fbam_bytes
                bf_nogil_logf_notime(STATS_TAG, "Filtered BAM saved to %s\n", _fbam_c)
    else:
        bf_nogil_logf_notime(STATS_TAG, "Processing failed with code %d\n", ret)

    return ret
