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

"""Reference connectivity graph analysis.

This module provides the high-performance Cython implementation for
analyzing read-reference connectivity patterns, computing per-reference
statistics, building co-mapping graphs, and exporting TSV summaries.

Key functionality
-----------------
- Build efficient reference -> read indices in parallel
- Compute neighbor/connectivity metrics (exact and approximate)
- Optionally construct an igraph/filtered weighted graph for clustering
- Export per-reference TSVs and summary statistics

Implementation notes
--------------------
- Most heavy functions are ``nogil`` and designed for parallel execution
    using ``prange``. Be careful when editing to preserve nogil/noexcept
    annotations and avoid Python-level operations inside nogil regions.
"""

from libc.stdlib cimport malloc, free, calloc, realloc, qsort
from libc.string cimport memset, memcpy, strlen
from libc.math cimport sqrt as libc_sqrt, log, fabs, exp, fmax, fmin, log2
from libc.stdio cimport FILE, fopen, fclose
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t
from cython.parallel cimport prange, threadid
from bam_filter.processor cimport min_int32, max_int32, min_int64, max_int64
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor_graph cimport ReadRefsIndex, MultPairExtended, RefReadPair, SparseConnectivity, GraphConfig
from bam_filter.processor_graph_tsv cimport write_graph_tsv_c

# Import graph operations from processor_graph_ops
from bam_filter.processor_graph_ops cimport (
    WeightedGraph,
    build_weighted_graph_from_alignments,
    prune_low_weight_edges,
    calculate_graph_statistics,
    destroy_weighted_graph
)

# Import igraph functions for graph building and statistics

from bam_filter.processor_graph_ops cimport build_igraph_from_read_index
from bam_filter.processor_graph_ops cimport build_igraph_direct_from_read_index
from bam_filter.processor_igraph cimport *
from bam_filter.processor_graph_ops cimport create_weighted_graph

cdef extern from "zlib.h":
    ctypedef void* gzFile
    gzFile gzopen(const char* path, const char* mode) nogil
    int gzclose(gzFile file) nogil
    int gzprintf(gzFile file, const char* format, ...) nogil
    int gzwrite(gzFile file, const void* buf, unsigned int len) nogil

cdef extern from "time.h":
    cdef struct timespec:
        long tv_sec
        long tv_nsec
    int clock_gettime(int clk_id, timespec *tp) nogil
    int CLOCK_MONOTONIC

cdef extern from "bam_filter/c_logging.h":
    double bf_monotonic_seconds() nogil
    void bf_nogil_log(const char* tag, const char* msg) nogil
    void bf_nogil_log_fmt(const char* tag, const char* fmt_msg, long v) nogil
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from *:
    uint32_t __sync_fetch_and_add(uint32_t *ptr, uint32_t value) nogil

cdef inline bint c_str_endswith_gz(const char* s) nogil:
    """Check whether a C-style string ends with the suffix ".gz".

    Parameters
    ----------
    s : const char*
        NUL-terminated C string (may be NULL).

    Returns
    -------
    bint
        True (1) if the string ends with ".gz"; False (0) otherwise.
    """
    if s == NULL:
        return 0
    cdef int L = <int>strlen(s)
    if L < 3:
        return 0

    if s[L-3] == 46 and s[L-2] == 103 and s[L-1] == 122:
        return 1
    return 0


cdef uint32_t UINT32_MAX = <uint32_t>0xFFFFFFFF


cdef inline igraph_real_t get_vector_element(igraph_vector_t *v, igraph_integer_t i) nogil:
    """Access element from an ``igraph_vector_t`` by raw pointer.

    This bypasses some deprecated accessor functions in the C API and
    reads the underlying double buffer directly. Use only in nogil
    contexts where the igraph vector is known to be initialized.
    """
    return (<igraph_real_t**>v)[0][i]


# Helper to access igraph_vector_int_t elements (integers)
cdef inline igraph_integer_t get_vector_int_element(igraph_vector_int_t *v, igraph_integer_t i) nogil:
    """Access integer element from an ``igraph_vector_int_t`` by raw pointer.

    Parameters
    ----------
    v : igraph_vector_int_t*
        Pointer to igraph integer vector.
    i : igraph_integer_t
        Index to read.

    Returns
    -------
    igraph_integer_t
        Value at index ``i``.
    """
    return (<igraph_integer_t**>v)[0][i]


cdef void count_alignments_chunk(uint64_t start_aln, uint64_t end_aln,
                                MemoryPool* pool, uint64_t* thread_alignments_per_ref) noexcept nogil:
    """Count alignments per-reference in a half-open alignment index range.

    This function increments the provided per-reference counter array for each
    alignment in the interval [start_aln, end_aln). It is designed to run in
    nogil parallel regions and performs only C-level memory operations.

    Parameters
    ----------
    start_aln : uint64_t
        Inclusive start index of the alignment range.
    end_aln : uint64_t
        Exclusive end index of the alignment range.
    pool : MemoryPool*
        Pointer to the MemoryPool containing alignment records.
    thread_alignments_per_ref : uint64_t*
        Preallocated per-reference counter array (length pool.reference_count).

    Notes
    -----
    The caller must ensure ``thread_alignments_per_ref`` points to a zeroed
    array with at least ``pool.reference_count`` entries.
    """
    cdef uint64_t aln_i
    cdef uint32_t aln_ref

    for aln_i in range(start_aln, end_aln):
        aln_ref = pool.alignments[aln_i].reference_index
        if aln_ref < pool.reference_count:
            thread_alignments_per_ref[aln_ref] += 1


cdef void count_reads_chunk(uint32_t start_read, uint32_t end_read, MemoryPool* pool,
                                  uint32_t* thread_total_reads, uint32_t* thread_multimap_reads) noexcept nogil:
    """For reads in a range, count per-reference read occurrences and multimapping reads.

    Iterates reads in ``[start_read, end_read)`` and, for each read, determines the
    set of unique references the read maps to. It increments per-reference totals
    and per-reference multimapping counters (reads mapping to >1 references).

    Parameters
    ----------
    start_read : uint32_t
        Inclusive start read index.
    end_read : uint32_t
        Exclusive end read index.
    pool : MemoryPool*
        MemoryPool with per-read alignment start/count arrays and alignments.
    thread_total_reads : uint32_t*
        Per-reference array to increment when a read maps to that reference.
    thread_multimap_reads : uint32_t*
        Per-reference array to increment when a read maps to >1 references.

    Notes
    -----
    This function allocates temporary buffers; callers should prefer a
    preallocated scratch when calling many times to reduce allocation cost.
    """
    cdef uint32_t read_idx_c, aln_i, aln_ref, ref_count_in_read
    cdef uint64_t start_pos, end_pos
    cdef uint32_t n_alignments
    cdef char* ref_seen = NULL
    cdef uint32_t* unique_refs = NULL
    cdef uint32_t unique_capacity = 256
    cdef uint32_t* new_unique_refs
    cdef uint32_t j

    ref_seen = <char*>calloc(pool.reference_count, sizeof(char))
    if not ref_seen:
        return

    unique_refs = <uint32_t*>malloc(unique_capacity * sizeof(uint32_t))
    if not unique_refs:
        free(ref_seen)
        return

    for read_idx_c in range(start_read, end_read):
        start_pos = pool.read_alignment_starts[read_idx_c]
        end_pos = start_pos + pool.read_alignment_counts[read_idx_c]
        n_alignments = <uint32_t>(end_pos - start_pos)

        if n_alignments == 0:
            continue

        ref_count_in_read = 0

        for aln_i in range(start_pos, end_pos):
            aln_ref = pool.alignments[aln_i].reference_index
            if aln_ref >= pool.reference_count:
                continue

            if not ref_seen[aln_ref]:
                ref_seen[aln_ref] = 1

                if ref_count_in_read >= unique_capacity:
                    unique_capacity = unique_capacity * 2
                    new_unique_refs = <uint32_t*>realloc(unique_refs, unique_capacity * sizeof(uint32_t))
                    if not new_unique_refs:
                        break
                    unique_refs = new_unique_refs

                unique_refs[ref_count_in_read] = aln_ref
                ref_count_in_read += 1

        for j in range(ref_count_in_read):
            aln_ref = unique_refs[j]
            thread_total_reads[aln_ref] += 1
            if ref_count_in_read > 1:
                thread_multimap_reads[aln_ref] += 1

        for j in range(ref_count_in_read):
            ref_seen[unique_refs[j]] = 0

    free(ref_seen)
    free(unique_refs)


cdef void reduce_thread_arrays(uint32_t** thread_arrays, uint32_t* final_array,
                              uint32_t array_size, int num_threads) noexcept nogil:
    """Sum per-thread arrays into a single final array.

    Parameters
    ----------
    thread_arrays : uint32_t**
        Array of pointers where each pointer references a per-thread array.
    final_array : uint32_t*
        Destination array where the element-wise sums will be accumulated.
    array_size : uint32_t
        Length of each per-thread array.
    num_threads : int
        Number of threads / pointers in ``thread_arrays``.

    Notes
    -----
    The caller must ensure ``final_array`` is zeroed or contains the desired
    initial accumulators before calling this function.
    """

    cdef uint32_t ref_idx
    cdef int thread_id

    for thread_id in range(num_threads):
        if thread_arrays[thread_id]:
            for ref_idx in range(array_size):
                final_array[ref_idx] += thread_arrays[thread_id][ref_idx]


cdef void insertion_sort_uint32(uint32_t* arr, uint32_t n) noexcept nogil:
    """Insertion sort for small arrays of 32-bit unsigned integers.

    This implementation is intentionally simple and optimized for very small
    arrays (used as a fast path before falling back to qsort for larger arrays).

    Parameters
    ----------
    arr : uint32_t*
        Pointer to the array to sort in-place.
    n : uint32_t
        Number of elements in the array.
    """
    cdef uint32_t i, j, key
    for i in range(1, n):
        key = arr[i]
        j = i
        while j > 0 and arr[j-1] > key:
            arr[j] = arr[j-1]
            j -= 1
        arr[j] = key


cdef void calculate_dataset_summary_stats(MemoryPool* pool, DatasetSummaryStats* stats) noexcept nogil:
    # Compute basic dataset summary statistics from the MemoryPool and populate `stats`.
    # Args:
    #     pool: MemoryPool containing alignment/read data
    #     stats: pointer to DatasetSummaryStats struct to fill
    memset(stats, 0, sizeof(DatasetSummaryStats))
    stats.pmd_enabled = pool.pmd_enabled_for_output
    stats.min_score = 1e30
    stats.max_score = -1e30

    cdef uint32_t read_idx, alignment_count
    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count > 0:
            stats.total_reads_processed += 1
            if alignment_count == 1:
                stats.unique_reads += 1
            else:
                stats.multimapping_reads += 1

    stats.total_alignments_filtered = pool.alignment_count
    stats.total_alignments_input = pool.original_alignment_count

    cdef int64_t i
    cdef double score, delta, M2 = 0.0
    cdef double mean = 0.0
    cdef uint64_t n = 0
    cdef double pmd_mean = 0.0, pmd_M2 = 0.0
    cdef uint64_t pmd_n = 0
    cdef float pmd_score
    cdef double std_error = 0.0

    for i in range(pool.alignment_count):
        score = pool.alignments[i].alignment_score

        n += 1
        delta = score - mean
        mean += delta / n
        M2 += delta * (score - mean)

        if score < stats.min_score:
            stats.min_score = score
        if score > stats.max_score:
            stats.max_score = score

        if stats.pmd_enabled:
            pmd_score = pool.alignments[i].pmd_score

            pmd_n += 1
            delta = pmd_score - pmd_mean
            pmd_mean += delta / pmd_n
            pmd_M2 += delta * (pmd_score - pmd_mean)

    stats.score_count = n
    stats.mean_score = mean
    if n > 1:
        stats.score_variance = M2 / (n - 1)

        std_error = libc_sqrt(stats.score_variance / n)
        stats.score_ci_lower_95 = mean - 1.96 * std_error
        stats.score_ci_upper_95 = mean + 1.96 * std_error

    if stats.pmd_enabled and pmd_n > 1:
        stats.mean_pmd_score = pmd_mean
        stats.pmd_variance = pmd_M2 / (pmd_n - 1)
        stats.pmd_count = pmd_n
    elif stats.pmd_enabled and pmd_n == 1:
        stats.mean_pmd_score = pmd_mean
        stats.pmd_variance = 0.0
        stats.pmd_count = pmd_n

    stats.total_references = pool.reference_count
    for i in range(pool.reference_count):
        if pool.reference_lengths and pool.reference_lengths[i] > 0:
            stats.references_with_alignments += 1


cdef void calculate_reference_stats(MemoryPool* pool, sam_hdr_t* bam_header,
                                    ReferenceStats* ref_stats) noexcept nogil:
    # Compute per-reference statistics such as read classification counts and score statistics.
    # Args:
    #     pool: MemoryPool pointer with alignment and read data
    #     bam_header: BAM header pointer (used only for context if needed)
    #     ref_stats: preallocated array of ReferenceStats (length reference_count)
    cdef uint32_t ref_idx, read_idx, alignment_idx, current_ref
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count

    cdef double score_val, pmd_val, delta, delta2
    cdef uint64_t n_score
    cdef bint pmd_enabled = pool.pmd_enabled_for_output

    cdef uint32_t* read_refs = NULL
    cdef uint32_t* ref_counts = NULL
    cdef uint32_t ref_capacity = 0
    cdef uint32_t unique_ref_count
    cdef uint32_t i
    cdef bint found
    cdef uint32_t alignments_to_this_ref

    for ref_idx in range(pool.reference_count):
        ref_stats[ref_idx].total_reads = 0
        ref_stats[ref_idx].unique_reads = 0
        ref_stats[ref_idx].repeat_reads = 0
        ref_stats[ref_idx].shared_reads = 0

        ref_stats[ref_idx].alignment_count = 0
        ref_stats[ref_idx].score_mean = 0.0
        ref_stats[ref_idx].score_variance = 0.0
        ref_stats[ref_idx].score_min = 1e30
        ref_stats[ref_idx].score_max = -1e30
        ref_stats[ref_idx].score_std = 0.0

        ref_stats[ref_idx].pmd_available = pmd_enabled
        ref_stats[ref_idx].pmd_mean = 0.0
        ref_stats[ref_idx].pmd_variance = 0.0
        ref_stats[ref_idx].pmd_min = 1e30
        ref_stats[ref_idx].pmd_max = -1e30
        ref_stats[ref_idx].pmd_std = 0.0
        ref_stats[ref_idx].pmd_nonzero_count = 0

    bf_nogil_logf_notime(NULL, b"reference_stats: computing mutually exclusive read categories")

    for alignment_idx in range(pool.alignment_count):
        current_ref = pool.alignments[alignment_idx].reference_index
        if current_ref >= pool.reference_count:
            continue

        ref_stats[current_ref].alignment_count += 1

        score_val = <double>pool.alignments[alignment_idx].alignment_score
        n_score = ref_stats[current_ref].alignment_count

        if score_val < ref_stats[current_ref].score_min:
            ref_stats[current_ref].score_min = score_val
        if score_val > ref_stats[current_ref].score_max:
            ref_stats[current_ref].score_max = score_val

        if n_score == 1:
            ref_stats[current_ref].score_mean = score_val
            ref_stats[current_ref].score_variance = 0.0
        else:
            delta = score_val - ref_stats[current_ref].score_mean
            ref_stats[current_ref].score_mean += delta / n_score
            delta2 = score_val - ref_stats[current_ref].score_mean
            ref_stats[current_ref].score_variance += delta * delta2

        if pmd_enabled:
            pmd_val = <double>pool.alignments[alignment_idx].pmd_score

            if pmd_val != 0.0:
                ref_stats[current_ref].pmd_nonzero_count += 1

            if pmd_val < ref_stats[current_ref].pmd_min:
                ref_stats[current_ref].pmd_min = pmd_val
            if pmd_val > ref_stats[current_ref].pmd_max:
                ref_stats[current_ref].pmd_max = pmd_val

            if n_score == 1:
                ref_stats[current_ref].pmd_mean = pmd_val
                ref_stats[current_ref].pmd_variance = 0.0
            else:
                delta = pmd_val - ref_stats[current_ref].pmd_mean
                ref_stats[current_ref].pmd_mean += delta / n_score
                delta2 = pmd_val - ref_stats[current_ref].pmd_mean
                ref_stats[current_ref].pmd_variance += delta * delta2

    for ref_idx in range(pool.reference_count):
        if ref_stats[ref_idx].alignment_count > 1:
            ref_stats[ref_idx].score_variance /= (ref_stats[ref_idx].alignment_count - 1)
            ref_stats[ref_idx].score_std = libc_sqrt(ref_stats[ref_idx].score_variance)

            if pmd_enabled:
                ref_stats[ref_idx].pmd_variance /= (ref_stats[ref_idx].alignment_count - 1)
                ref_stats[ref_idx].pmd_std = libc_sqrt(ref_stats[ref_idx].pmd_variance)
        else:
            ref_stats[ref_idx].score_std = 0.0
            if pmd_enabled:
                ref_stats[ref_idx].pmd_std = 0.0

        if ref_stats[ref_idx].alignment_count == 0:
            ref_stats[ref_idx].score_min = 0.0
            if pmd_enabled:
                ref_stats[ref_idx].pmd_min = 0.0

    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        unique_ref_count = 0

        if alignment_count > ref_capacity:

            if read_refs: free(read_refs)
            if ref_counts: free(ref_counts)

            ref_capacity = alignment_count
            read_refs = <uint32_t*>malloc(ref_capacity * sizeof(uint32_t))
            ref_counts = <uint32_t*>malloc(ref_capacity * sizeof(uint32_t))

            if not read_refs or not ref_counts:
                bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate arrays for %u references per read\n", ref_capacity)
                if read_refs: free(read_refs)
                if ref_counts: free(ref_counts)
                return

        for alignment_idx in range(start_pos, end_pos):
            current_ref = pool.alignments[alignment_idx].reference_index
            if current_ref >= pool.reference_count:
                continue

            found = False
            i = <uint32_t>0
            while i < unique_ref_count:
                if read_refs[i] == current_ref:
                    ref_counts[i] += 1
                    found = True
                    break
                i += 1

            if not found:
                read_refs[unique_ref_count] = current_ref
                ref_counts[unique_ref_count] = 1
                unique_ref_count += 1

        # Classify this read into one of three mutually exclusive categories per reference:
        # - unique_reads: Read maps to ONLY this ref with exactly 1 alignment
        # - repeat_reads: Read maps to ONLY this ref with multiple alignments
        # - shared_reads: Read maps to MULTIPLE refs (regardless of alignment count to each)
        # These categories partition reads: total_reads = unique_reads + repeat_reads + shared_reads
        i = <uint32_t>0
        while i < unique_ref_count:
            current_ref = read_refs[i]
            alignments_to_this_ref = ref_counts[i]

            if current_ref < pool.reference_count:

                ref_stats[current_ref].total_reads += 1

                if unique_ref_count == 1:
                    # Read maps to ONLY this reference
                    if alignments_to_this_ref == 1:
                        ref_stats[current_ref].unique_reads += 1
                    else:
                        ref_stats[current_ref].repeat_reads += 1
                else:
                    # Read maps to MULTIPLE references (shared)
                    # Note: Even if this read has multiple alignments to this ref,
                    # it's still classified as "shared" not "repeat" because the
                    # categorization is based on reference count, not alignment count
                    ref_stats[current_ref].shared_reads += 1
            i += 1

    for ref_idx in range(pool.reference_count):
        if ref_stats[ref_idx].total_reads > 0:
            ref_stats[ref_idx].unique_percentage = (100.0 * ref_stats[ref_idx].unique_reads) / ref_stats[ref_idx].total_reads
            ref_stats[ref_idx].repeat_percentage = (100.0 * ref_stats[ref_idx].repeat_reads) / ref_stats[ref_idx].total_reads
            ref_stats[ref_idx].shared_percentage = (100.0 * ref_stats[ref_idx].shared_reads) / ref_stats[ref_idx].total_reads
            
            # Validate that categories sum correctly (mutually exclusive and exhaustive)
            category_sum = ref_stats[ref_idx].unique_reads + ref_stats[ref_idx].repeat_reads + ref_stats[ref_idx].shared_reads
            if category_sum != ref_stats[ref_idx].total_reads:
                bf_nogil_logf_notime(
                    NULL,
                    "ERROR: Ref %u category sum mismatch: total=%u vs unique+repeat+shared=%u\n",
                    ref_idx,
                    ref_stats[ref_idx].total_reads,
                    category_sum,
                )

    bf_nogil_log(b"GRAPH", b"reference_stats: read category tallies verified")

cdef int _multpair_extended_cmp(const void* a, const void* b) noexcept nogil:
    # Comparison for MultPairExtended used for sorting by neighbor_count, then read_count
    cdef MultPairExtended* A = <MultPairExtended*>a
    cdef MultPairExtended* B = <MultPairExtended*>b

    # Primary sort: neighbor_count (descending)
    if A.neighbor_count != B.neighbor_count:
        return -1 if A.neighbor_count > B.neighbor_count else 1

    if A.neighbor_count != B.neighbor_count:
        return -1 if A.neighbor_count > B.neighbor_count else 1

    if A.read_count != B.read_count:
        return -1 if A.read_count > B.read_count else 1

    return 0

cdef double calculate_elapsed_seconds(timespec* start, timespec* end) noexcept nogil:
    # Compute elapsed wall-clock seconds between two timespec points.
    # Args:
    #     start, end: pointers to timespec structures
    # Returns:
    #     elapsed time in seconds as a double
    cdef long sec_diff = end.tv_sec - start.tv_sec
    cdef long nsec_diff = end.tv_nsec - start.tv_nsec

    if nsec_diff < 0:
        sec_diff -= 1
        nsec_diff += 1000000000

    return <double>sec_diff + <double>nsec_diff / 1000000000.0


cdef void calculate_unique_read_counts_selective(MemoryPool* pool, ReferencePattern* pattern_data,
                                                  uint32_t* total_reads, 
                                                  int32_t low_coverage_threshold) noexcept nogil:
    """
    Calculate unique_read_count ONLY for references that survived probability filtering.
    This is more efficient than calculating for all references.
    
    A read is considered "unique" to a reference if it maps ONLY to that reference
    (not shared with any other reference).
    
    Args:
        pool: MemoryPool with alignment data
        pattern_data: Array to populate unique_read_count field
        total_reads: Total read counts per reference (from graph analysis)
        low_coverage_threshold: Minimum reads required (min_read_count)
    """
    cdef uint32_t read_idx, ref_idx, current_ref
    cdef uint64_t start_pos, end_pos, alignment_idx
    cdef uint32_t alignment_count
    
    # Track which references each read maps to
    cdef uint32_t* read_refs = NULL
    cdef uint32_t ref_capacity = 0
    cdef uint32_t unique_ref_count
    cdef uint32_t i
    cdef bint found
    cdef uint32_t surviving_refs = 0
    
    # Count how many references we're processing
    for ref_idx in range(pool.reference_count):
        if total_reads[ref_idx] >= <uint32_t>low_coverage_threshold:
            surviving_refs += 1
    
    bf_nogil_logf_notime(NULL, "  Processing unique read counts for %u references (out of %u total)\n",
                        surviving_refs, pool.reference_count)
    
    # Initialize unique_read_count to 0 for all references
    for ref_idx in range(pool.reference_count):
        pattern_data[ref_idx].unique_read_count = 0
    
    # Iterate through each read
    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue
        
        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count
        
        # Allocate/reallocate buffer if needed
        if alignment_count > ref_capacity:
            if read_refs:
                free(read_refs)
            
            ref_capacity = alignment_count
            read_refs = <uint32_t*>malloc(ref_capacity * sizeof(uint32_t))
            if not read_refs:
                bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate read_refs buffer\n")
                return
        
        # Find unique references for this read
        unique_ref_count = 0
        for alignment_idx in range(start_pos, end_pos):
            current_ref = pool.alignments[alignment_idx].reference_index
            if current_ref >= pool.reference_count:
                continue
            
            # Check if we've seen this reference already for this read
            found = False
            for i in range(unique_ref_count):
                if read_refs[i] == current_ref:
                    found = True
                    break
            
            if not found:
                read_refs[unique_ref_count] = current_ref
                unique_ref_count += 1
        
        # If this read maps to only ONE reference, increment that reference's unique count
        # BUT ONLY if that reference is in our surviving set
        if unique_ref_count == 1:
            current_ref = read_refs[0]
            if current_ref < pool.reference_count:
                if total_reads[current_ref] >= <uint32_t>low_coverage_threshold:
                    pattern_data[current_ref].unique_read_count += 1
    
    if read_refs:
        free(read_refs)


cdef void print_pattern_summary(MemoryPool* memory_pool, sam_hdr_t* bam_header,
                                 ReferenceMapping* mapping, ReferencePattern* pattern_data,
                                 GraphConfig* gconfig, ReferenceStats* ref_stats) noexcept nogil:
    """Emit a human-readable summary of co-mapping pattern statistics.

    Parameters
    ----------
    memory_pool : MemoryPool*
        Memory pool containing dataset statistics.
    bam_header : sam_hdr_t*
        BAM header for reference name/length lookups (may be NULL).
    mapping : ReferenceMapping*
        Optional mapping from compacted retained reference IDs to original TIDs.
    pattern_data : ReferencePattern*
        Array of per-reference pattern metrics.
    gconfig : GraphConfig*
        Optional configuration controlling displayed output; may be NULL.
    ref_stats : ReferenceStats*
        Optional per-reference statistics used to enrich the summary.

    Notes
    -----
    This function is nogil-safe and writes to stdout / nogil logging helpers.
    It does not allocate large temporary buffers.
    """
    cdef DatasetSummaryStats stats
    cdef uint32_t ref_idx
    cdef uint32_t display_count
    cdef const char* ref_name
    cdef int64_t ref_len
    cdef uint32_t original_tid
    cdef int64_t i, j, valid_count = 0
    cdef double penalty_percentile, total_percentage
    cdef int32_t low_sharing = 0, symmetric = 0, asymmetric = 0, hub = 0, isolated = 0
    cdef MultPairExtended* sorted_refs
    cdef int min_read_threshold = 1 

    cdef double min_penalty = 1e30, max_penalty = -1e30 
    cdef double penalty_sum = 0.0, penalty_sum_sq = 0.0
    cdef uint32_t penalty_bins[20]
    cdef uint32_t bin_idx

    cdef double range_starts[5]
    cdef double range_ends[5] 
    cdef uint32_t range_counts[5]
    cdef uint32_t samples_per_range
    cdef uint32_t range_idx, sample_count, total_samples = 0
    cdef double range_min, range_max, penalty_range
    cdef uint32_t refs_in_range, sample_step, sample_idx

    cdef MultPairExtended* range_refs = NULL
    cdef uint32_t range_count, actual_samples

    cdef uint64_t global_max_co_mappings = 0
    cdef uint64_t total_connections = 0
    cdef uint32_t refs_with_high_connections = 0

    cdef uint32_t global_unique_reads = 0, global_repeat_reads = 0, global_shared_reads = 0
    cdef uint32_t read_idx, alignment_count, unique_ref_count
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_idx
    cdef uint32_t* read_refs = NULL
    cdef char* ref_seen = NULL
    cdef bint found

    cdef double base_multimap, multimap_component, connection_penalty
    cdef double self_multimap, neighbor_influence, neighbor_penalty
    cdef double penalty_val, penalty_mean, penalty_variance, penalty_std
    cdef uint32_t neighbor_count_local
    cdef double neighbor_multimap_local, neighbor_connections_local

    memset(penalty_bins, 0, 20 * sizeof(uint32_t))
    memset(range_counts, 0, 5 * sizeof(uint32_t))

    if gconfig:
        display_count = gconfig.max_display_count if gconfig.max_display_count > 0 else 100
        samples_per_range = gconfig.samples_per_stratum if gconfig.samples_per_stratum > 0 else (display_count / 5)
    else:
        display_count = 100
        samples_per_range = display_count / 5

    calculate_dataset_summary_stats(memory_pool, &stats)

    if not ref_stats:
        bf_nogil_logf_notime(NULL, "ERROR: ref_stats not provided to print_pattern_summary\n")
        return

    sorted_refs = <MultPairExtended*>malloc(memory_pool.reference_count * sizeof(MultPairExtended))
    if not sorted_refs:
        bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate sorted reference array\n")
        return

    bf_nogil_logf_notime(NULL, "=== GLOBAL READ CLASSIFICATION ===\n")

    for read_idx in range(memory_pool.unique_read_count):
        alignment_count = memory_pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = memory_pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        if alignment_count == 1:
            global_unique_reads += 1
        else:

            read_refs = <uint32_t*>malloc(alignment_count * sizeof(uint32_t))
            ref_seen = <char*>calloc(memory_pool.reference_count, sizeof(char))

            if read_refs and ref_seen:
                unique_ref_count = 0

                for alignment_idx in range(start_pos, end_pos):
                    ref_idx = memory_pool.alignments[alignment_idx].reference_index
                    if ref_idx < memory_pool.reference_count and not ref_seen[ref_idx]:
                        ref_seen[ref_idx] = 1
                        read_refs[unique_ref_count] = ref_idx
                        unique_ref_count += 1

                if unique_ref_count == 1:
                    global_repeat_reads += 1
                else:
                    global_shared_reads += 1

            if read_refs: free(read_refs)
            if ref_seen: free(ref_seen)
            read_refs = NULL
            ref_seen = NULL

    cdef uint32_t total_classified = global_unique_reads + global_repeat_reads + global_shared_reads

    bf_nogil_logf_notime(NULL, "Global Read Classification:\n")
    bf_nogil_logf_notime(NULL, "  Total reads: %u\n", total_classified)
    bf_nogil_logf_notime(NULL, "  Unique reads: %u (%.1f%%) - maps only to 1 reference, 1 alignment\n",
                   global_unique_reads, 100.0 * global_unique_reads / total_classified if total_classified > 0 else 0.0)
    bf_nogil_logf_notime(NULL, "  Repeat reads: %u (%.1f%%) - maps only to 1 reference, multiple alignments\n",
                   global_repeat_reads, 100.0 * global_repeat_reads / total_classified if total_classified > 0 else 0.0)
    bf_nogil_logf_notime(NULL, "  Shared reads: %u (%.1f%%) - maps to multiple references\n",
                   global_shared_reads, 100.0 * global_shared_reads / total_classified if total_classified > 0 else 0.0)

    if total_classified != stats.total_reads_processed:
        bf_nogil_logf_notime(NULL, "  WARNING: Classification total (%u) != processed reads (%lu)\n", 
                      total_classified, stats.total_reads_processed)

    bf_nogil_logf_notime(NULL, "\nCross-validation with existing stats:\n")
    bf_nogil_logf_notime(NULL, "  Existing unique reads: %lu\n", stats.unique_reads)
    bf_nogil_logf_notime(NULL, "  Existing multi-mapping reads: %lu\n", stats.multimapping_reads)
    bf_nogil_logf_notime(NULL, "  New repeat + shared reads: %u\n", global_repeat_reads + global_shared_reads)

    if stats.unique_reads != global_unique_reads:
        bf_nogil_logf_notime(
            NULL,
            "  NOTE: Unique read counts differ (existing=%lu, new=%u)\n",
            stats.unique_reads,
            global_unique_reads,
        )

    bf_nogil_logf_notime(NULL, "\n=== MULTI-MAPPING ANALYSIS RESULTS ===\n")
    bf_nogil_logf_notime(NULL, "Reads:\n")
    bf_nogil_logf_notime(NULL, "  Total reads processed: %lu\n", stats.total_reads_processed)
    bf_nogil_logf_notime(NULL, "  Unique mapping reads: %lu (%.1f%%)\n", 
             stats.unique_reads, 100.0 * stats.unique_reads / stats.total_reads_processed)
    bf_nogil_logf_notime(NULL, "  Multi-mapping reads: %lu (%.1f%%)\n", 
             stats.multimapping_reads, 100.0 * stats.multimapping_reads / stats.total_reads_processed)

    bf_nogil_logf_notime(NULL, "\nAlignments:\n")
    bf_nogil_logf_notime(NULL, "  Total input alignments: %lu\n", stats.total_alignments_input)
    bf_nogil_logf_notime(NULL, "  Filtered alignments: %lu\n", stats.total_alignments_filtered)
    bf_nogil_logf_notime(
        NULL,
        "  Filtering rate: %.1f%%\n",
        100.0 * (stats.total_alignments_input - stats.total_alignments_filtered) / stats.total_alignments_input,
    )

    bf_nogil_logf_notime(NULL, "\nReferences:\n")
    bf_nogil_logf_notime(NULL, "  Total references: %u\n", stats.total_references)
    bf_nogil_logf_notime(NULL, "  References with alignments: %u\n", stats.references_with_alignments)

    bf_nogil_logf_notime(NULL, "\nAlignment Scores:\n")
    bf_nogil_logf_notime(NULL, "  Count: %lu\n", stats.score_count)
    bf_nogil_logf_notime(NULL, "  Range: [%.3f, %.3f]\n", stats.min_score, stats.max_score)
    bf_nogil_logf_notime(
        NULL,
        "  Mean: %.3f (95%% CI: %.3f-%.3f)\n",
        stats.mean_score,
        stats.score_ci_lower_95,
        stats.score_ci_upper_95,
    )
    bf_nogil_logf_notime(NULL, "  Std deviation: %.3f\n", libc_sqrt(stats.score_variance))

    if stats.pmd_enabled and stats.pmd_count > 0:
        bf_nogil_logf_notime(NULL, "\nPMD Scores:\n")
        bf_nogil_logf_notime(NULL, "  Count: %lu\n", stats.pmd_count)
        bf_nogil_logf_notime(NULL, "  Mean: %.4f\n", stats.mean_pmd_score)
        bf_nogil_logf_notime(NULL, "  Std deviation: %.4f\n", libc_sqrt(stats.pmd_variance))

    for ref_idx in range(memory_pool.reference_count):
        # Include references that have any alignments/read counts even if
        # they have no co-mappings recorded. Previously we skipped when
        # reads_with_comappings == 0 which caused "No valid references
        # found" when TSV had per-reference stats but no co-mapping data.
        # Consider a reference valid if it has any total reads in ref_stats
        # or has recorded co-mappings.
        if pattern_data[ref_idx].graph.reads_with_comappings == 0 and ref_stats[ref_idx].total_reads == 0:
            continue

        if pattern_data[ref_idx].graph.max_comappings > global_max_co_mappings:
            global_max_co_mappings = pattern_data[ref_idx].graph.max_comappings

        total_connections += pattern_data[ref_idx].graph.max_comappings

        if pattern_data[ref_idx].graph.max_comappings > 10000:
            refs_with_high_connections += 1

        sorted_refs[valid_count].neighbor_count = pattern_data[ref_idx].graph.connection_count
        sorted_refs[valid_count].read_count = ref_stats[ref_idx].total_reads
        sorted_refs[valid_count].idx = ref_idx
        # Calculate multimap fraction on-the-fly from actual read stats
        if ref_stats[ref_idx].total_reads > 0:
            sorted_refs[valid_count].conc = <float>(ref_stats[ref_idx].repeat_reads + ref_stats[ref_idx].shared_reads) / <float>ref_stats[ref_idx].total_reads
        else:
            sorted_refs[valid_count].conc = 0.0
        sorted_refs[valid_count].net = pattern_data[ref_idx].graph.avg_comappings_per_read
        valid_count += 1

    if valid_count == 0:
        bf_nogil_logf_notime(NULL, "No valid references found\n")
        free(sorted_refs)
        return

    bf_nogil_logf_notime(NULL, "\n=== MULTI-MAPPING STATISTICS ===\n")
    bf_nogil_logf_notime(NULL, "Maximum connections observed: %lu (NO ARTIFICIAL LIMITS)\n", global_max_co_mappings)
    bf_nogil_logf_notime(NULL, "Average max connections per reference: %.1f\n", 
                  <double>total_connections / <double>valid_count)
    bf_nogil_logf_notime(NULL, "References with >10K connections: %u\n", refs_with_high_connections)

    # Sort references by neighbor count (descending) then read count (descending)
    if valid_count > 1:
        qsort(sorted_refs, valid_count, sizeof(MultPairExtended), _multpair_extended_cmp)

    bf_nogil_logf_notime(NULL, "\n=== TOP REFERENCES BY CONNECTIVITY ===\n")
    bf_nogil_logf_notime(NULL, "Showing top %d references by neighbor count\n", min(display_count, valid_count))
    bf_nogil_logf_notime(NULL, "Sorted by: neighbors (desc) -> reads (desc)\n\n")

    cdef uint32_t actual_display = min_int32(display_count, valid_count)
    for i in range(actual_display):
        ref_idx = sorted_refs[i].idx

        if mapping and ref_idx < mapping.n_retained_refs:
            original_tid = mapping.new_to_old_tid[ref_idx]
            ref_name = sam_hdr_tid2name(bam_header, original_tid)
            ref_len = sam_hdr_tid2len(bam_header, original_tid)
        else:
            original_tid = ref_idx
            ref_name = b"<unmapped>"
            ref_len = 0

        if not ref_name:
            ref_name = b"<unknown>"

        bf_nogil_logf_notime(
            NULL,
            "  #%d: %u->%u (%s) [%ld bp]\n",
            i + 1,
            ref_idx,
            original_tid,
            ref_name,
            ref_len,
        )

        bf_nogil_logf_notime(NULL, "    Reference Metrics:\n")

        bf_nogil_logf_notime(NULL, "    Read Classification:\n")
        bf_nogil_logf_notime(NULL, "      Total reads: %u\n", ref_stats[ref_idx].total_reads)
        bf_nogil_logf_notime(
            NULL,
            "      Unique reads: %u (%.1f%%) - maps only here, 1 alignment\n",
            ref_stats[ref_idx].unique_reads,
            ref_stats[ref_idx].unique_percentage,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Repeat reads: %u (%.1f%%) - maps only here, multiple alignments\n",
            ref_stats[ref_idx].repeat_reads,
            ref_stats[ref_idx].repeat_percentage,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Shared reads: %u (%.1f%%) - also maps to other references\n",
            ref_stats[ref_idx].shared_reads,
            ref_stats[ref_idx].shared_percentage,
        )

        bf_nogil_logf_notime(NULL, "    Alignment Score Statistics:\n")
        bf_nogil_logf_notime(NULL, "      Total alignments: %lu\n", ref_stats[ref_idx].alignment_count)
        if ref_stats[ref_idx].alignment_count > 0:
            bf_nogil_logf_notime(
                NULL,
                "      Score range: [%.3f, %.3f]\n",
                ref_stats[ref_idx].score_min,
                ref_stats[ref_idx].score_max,
            )
            bf_nogil_logf_notime(
                NULL,
                "      Score mean: %.3f +/- %.3f\n",
                ref_stats[ref_idx].score_mean,
                ref_stats[ref_idx].score_std,
            )

        if ref_stats[ref_idx].pmd_available and ref_stats[ref_idx].alignment_count > 0:
            bf_nogil_logf_notime(NULL, "    PMD Score Statistics:\n")
            bf_nogil_logf_notime(
                NULL,
                "      PMD range: [%.4f, %.4f]\n",
                ref_stats[ref_idx].pmd_min,
                ref_stats[ref_idx].pmd_max,
            )
            bf_nogil_logf_notime(
                NULL,
                "      PMD mean: %.4f +/- %.4f\n",
                ref_stats[ref_idx].pmd_mean,
                ref_stats[ref_idx].pmd_std,
            )
            bf_nogil_logf_notime(
                NULL,
                "      Non-zero PMD scores: %lu (%.1f%%)\n",
                ref_stats[ref_idx].pmd_nonzero_count,
                100.0 * ref_stats[ref_idx].pmd_nonzero_count / ref_stats[ref_idx].alignment_count,
            )

        bf_nogil_logf_notime(NULL, "    Network Complexity:\n")
        bf_nogil_logf_notime(
            NULL,
            "      Connected references: %u\n",
            pattern_data[ref_idx].graph.connection_count,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Avg co-mappings per read: %.1f\n",
            pattern_data[ref_idx].graph.avg_comappings_per_read,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Max co-mappings observed: %lu\n",
            pattern_data[ref_idx].graph.max_comappings,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Reads with co-mappings: %u\n",
            pattern_data[ref_idx].graph.reads_with_comappings,
        )

        bf_nogil_logf_notime(NULL, "    Neighbor Quality Analysis:\n")

        neighbor_count_local = pattern_data[ref_idx].graph.connection_count
        neighbor_multimap_local = ref_stats[ref_idx].shared_percentage / 100.0
        neighbor_connections_local = pattern_data[ref_idx].graph.avg_comappings_per_read

        bf_nogil_logf_notime(NULL, "      Connected neighbors: %u\n", neighbor_count_local)
        bf_nogil_logf_notime(
            NULL,
            "      Avg neighbor multimap rate: %.1f%% (vs %.1f%% self)\n",
            neighbor_multimap_local * 100.0,
            ref_stats[ref_idx].shared_percentage + ref_stats[ref_idx].repeat_percentage,
        )
        bf_nogil_logf_notime(
            NULL,
            "      Avg neighbor connections: %.1f (vs %u self)\n",
            neighbor_connections_local,
            pattern_data[ref_idx].graph.connection_count,
        )

        self_multimap = (ref_stats[ref_idx].shared_percentage + ref_stats[ref_idx].repeat_percentage) / 100.0
        neighbor_influence = neighbor_multimap_local - self_multimap
        bf_nogil_logf_notime(NULL, "\n")

    free(sorted_refs)


# MAIN FUNCTION (SEQUENTIAL IMPLEMENTATION)

cdef int count_unique_refs_thread_local(uint32_t read_idx, MemoryPool* pool, 
                                           uint32_t array_size, uint32_t* thread_counts,
                                           uint32_t* scratch, uint32_t scratch_capacity) noexcept nogil:
    # Count unique references for a single read and increment thread-local counters.
    # Args:
    #     read_idx: index of the read to process
    #     pool: MemoryPool pointer
    #     array_size: upper bound for valid reference indices
    #     thread_counts: per-reference counter array to increment
    #     scratch, scratch_capacity: optional scratch buffer and capacity
    # Returns:
    #     0 on success, -1 on allocation failure

    cdef uint32_t ref_count = pool.read_alignment_counts[read_idx]
    cdef uint64_t start_pos, end_pos
    cdef uint32_t* read_refs
    cdef uint32_t valid_count = 0
    cdef uint32_t alignment_idx, ref_idx
    cdef uint32_t prev_ref
    cdef uint32_t i
    cdef uint32_t slot, other_ref
    cdef bint used_malloc

    if ref_count < 1:
        return 0

    start_pos = pool.read_alignment_starts[read_idx]
    end_pos = start_pos + ref_count

    used_malloc = False
    if scratch and scratch_capacity >= ref_count:
        read_refs = scratch
    else:
        read_refs = <uint32_t*>malloc(ref_count * sizeof(uint32_t))
        if not read_refs:
            return -1
        used_malloc = True

    for alignment_idx in range(start_pos, end_pos):
        ref_idx = pool.alignments[alignment_idx].reference_index
        if ref_idx < array_size:
            read_refs[valid_count] = ref_idx
            valid_count += 1

    if valid_count == 0:
        if used_malloc:
            free(read_refs)
        return 0

    if valid_count <= 32:
        insertion_sort_uint32(read_refs, valid_count)
    else:
        qsort(read_refs, valid_count, sizeof(uint32_t), _uint32_compare)

    prev_ref = read_refs[0]
    thread_counts[prev_ref] += 1 

    for i in range(1, valid_count):
        if read_refs[i] != prev_ref:
            thread_counts[read_refs[i]] += 1 
            prev_ref = read_refs[i]

    if used_malloc:
        free(read_refs)
    return 0

cdef int fill_ref_to_reads_thread_local(uint32_t read_idx, MemoryPool* pool,
                                 uint32_t array_size, uint32_t** ref_to_reads, uint32_t* write_positions,
                                 uint32_t* scratch, uint32_t scratch_capacity) noexcept nogil:
    # For a single read, write the read index into per-reference lists (ref_to_reads)
    # using atomic slots from write_positions.
    # Args:
    #     read_idx: index of the read
    #     pool: MemoryPool pointer
    #     array_size: number of references
    #     ref_to_reads: array of pointers where each points to a block to store read indices
    #     write_positions: per-ref atomic counters used to claim write slots
    #     scratch, scratch_capacity: optional scratch buffer
    # Returns:
    #     0 on success, -1 on allocation failure

    cdef uint32_t ref_count = pool.read_alignment_counts[read_idx]
    cdef uint64_t start_pos, end_pos
    cdef uint32_t* read_refs
    cdef uint32_t valid_count = 0
    cdef uint32_t alignment_idx, ref_idx
    cdef uint32_t prev_ref
    cdef uint32_t i
    cdef bint used_malloc

    if ref_count < 1:
        return 0

    start_pos = pool.read_alignment_starts[read_idx]
    end_pos = start_pos + ref_count

    used_malloc = False
    if scratch and scratch_capacity >= ref_count:
        read_refs = scratch
    else:
        read_refs = <uint32_t*>malloc(ref_count * sizeof(uint32_t))
        if not read_refs:
            return -1
        used_malloc = True

    for alignment_idx in range(start_pos, end_pos):
        ref_idx = pool.alignments[alignment_idx].reference_index
        if ref_idx < array_size:
            read_refs[valid_count] = ref_idx
            valid_count += 1

    if valid_count > 0:
        if valid_count <= 32:
            insertion_sort_uint32(read_refs, valid_count)
        else:
            qsort(read_refs, valid_count, sizeof(uint32_t), _uint32_compare)

        prev_ref = read_refs[0]
        slot = __sync_fetch_and_add(&write_positions[prev_ref], 1)
        ref_to_reads[prev_ref][slot] = read_idx

        for i in range(1, valid_count):
            if read_refs[i] != prev_ref:
                other_ref = read_refs[i]
                slot = __sync_fetch_and_add(&write_positions[other_ref], 1)
                ref_to_reads[other_ref][slot] = read_idx
                prev_ref = read_refs[i]

    if used_malloc:
        free(read_refs)
    return 0

cdef ReadIndex* build_read_index_parallel(MemoryPool* pool, uint32_t array_size, int num_threads) noexcept nogil:
    # Build an index mapping references -> list of read indices that map to them using parallel work.
    # Args:
    #     pool: MemoryPool containing read->alignment mappings
    #     array_size: number of references to consider
    #     num_threads: number of parallel threads to use
    # Returns:
    #     pointer to newly allocated ReadIndex or NULL on error

    cdef ReadIndex* index
    cdef uint32_t** thread_counts
    cdef uint64_t total_mappings = 0
    cdef uint32_t i
    cdef int thread_id
    cdef uint32_t* read_buffer
    cdef uint32_t** ref_to_reads
    cdef uint32_t* write_positions
    cdef uint32_t buffer_pos = 0
    cdef uint32_t read_idx
    cdef uint32_t* final_counts
    cdef uint32_t** scratch_ptrs
    cdef uint32_t* seq_scratch

    index = <ReadIndex*>malloc(sizeof(ReadIndex))
    if not index:
        return NULL

    thread_counts = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    if not thread_counts:
        free(index)
        return NULL

    cdef int cleanup_id
    for thread_id in range(num_threads):
        thread_counts[thread_id] = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        if not thread_counts[thread_id]:

            for cleanup_id in range(thread_id):
                free(thread_counts[cleanup_id])
            free(thread_counts)
            free(index)
            return NULL

    cdef uint32_t max_ref_count = 0
    for i in range(pool.unique_read_count):
        if pool.read_alignment_counts[i] > max_ref_count:
            max_ref_count = pool.read_alignment_counts[i]

    scratch_ptrs = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    if not scratch_ptrs:
        for thread_id in range(num_threads):
            free(thread_counts[thread_id])
        free(thread_counts)
        free(index)
        return NULL

    for thread_id in range(num_threads):
        if max_ref_count > 0:
            scratch_ptrs[thread_id] = <uint32_t*>malloc(max_ref_count * sizeof(uint32_t))
            if not scratch_ptrs[thread_id]:

                for cleanup_id in range(thread_id):
                    if scratch_ptrs[cleanup_id]: free(scratch_ptrs[cleanup_id])
                free(scratch_ptrs)
                for thread_id in range(num_threads):
                    free(thread_counts[thread_id])
                free(thread_counts)
                free(index)
                return NULL
        else:
            scratch_ptrs[thread_id] = NULL

    for read_idx in prange(pool.unique_read_count, nogil=True, schedule='static', num_threads=num_threads):
        thread_id = threadid()
        if thread_id < num_threads:
            count_unique_refs_thread_local(read_idx, pool, array_size, thread_counts[thread_id], scratch_ptrs[thread_id], max_ref_count)

    final_counts = <uint32_t*>calloc(array_size, sizeof(uint32_t))
    if not final_counts:
        for thread_id in range(num_threads):
            free(thread_counts[thread_id])
        free(thread_counts)
        free(index)
        return NULL

    for thread_id in range(num_threads):
        for i in range(array_size):
            final_counts[i] += thread_counts[thread_id][i]
        free(thread_counts[thread_id])
    free(thread_counts)

    for i in range(array_size):
        total_mappings += final_counts[i]

    if total_mappings == 0:
        free(final_counts)
        free(index)
        return NULL

    read_buffer = <uint32_t*>malloc(total_mappings * sizeof(uint32_t))
    ref_to_reads = <uint32_t**>malloc(array_size * sizeof(uint32_t*))
    write_positions = <uint32_t*>calloc(array_size, sizeof(uint32_t))

    if not read_buffer or not ref_to_reads or not write_positions:
        if read_buffer: free(read_buffer)
        if ref_to_reads: free(ref_to_reads)
        if write_positions: free(write_positions)
        free(final_counts)
        free(index)
        return NULL

    for i in range(array_size):
        ref_to_reads[i] = read_buffer + buffer_pos
        buffer_pos += final_counts[i]

    if max_ref_count == 0:

        pass
    else:

        for read_idx in prange(pool.unique_read_count, nogil=True, schedule='dynamic', num_threads=num_threads):
            thread_id = threadid()
            if thread_id >= num_threads:
                thread_id = 0
            fill_ref_to_reads_thread_local(read_idx, pool, array_size, ref_to_reads, write_positions, scratch_ptrs[thread_id], max_ref_count)

    for thread_id in range(num_threads):
        if scratch_ptrs[thread_id]: free(scratch_ptrs[thread_id])
    free(scratch_ptrs)

    index.ref_to_reads = ref_to_reads
    index.ref_read_counts = final_counts
    index.read_buffer = read_buffer

    free(write_positions)
    return index


cdef void compute_reference_neighbor_metrics(uint32_t start_ref, uint32_t end_ref,
                                                    ReadIndex* index, ReadRefsIndex* read_refs_idx, 
                                                    uint32_t array_size, uint32_t* connection_counts, 
                                                    MemoryPool* pool, uint32_t* stamp,
                                                    double* neighbor_multimap_avg,
                                                    double* neighbor_connections_avg, 
                                                    uint32_t* neighbor_counts,
                                                    float* dataset_multimap_fractions,
                                                    uint32_t* dataset_connections,
                                                    uint32_t low_coverage_threshold) noexcept nogil:
    # For each reference in [start_ref, end_ref), compute neighbor/connectivity metrics:
    # - number of connected neighbor references (unique)
    # - average neighbor multimap fraction and connections
    # Args:
    #     index/read_refs_idx: optional precomputed indices to accelerate neighbor discovery
    #     connection_counts: output per-ref connected neighbor counts
    #     stamp: per-thread marker array used to deduplicate neighbor counting
    #     dataset_multimap_fractions/dataset_connections: per-ref dataset-level metrics used in averages
    #     low_coverage_threshold: skip refs with fewer reads than this

    cdef uint32_t target_ref, read_count, read_idx, i
    cdef uint32_t ref_count, other_ref
    cdef uint32_t connection_count
    cdef uint32_t gen
    cdef uint32_t j
    cdef uint32_t* refs_ptr
    cdef double neighbor_multimap_sum, neighbor_connections_sum
    cdef uint32_t neighbor_count
    cdef uint64_t start_pos 

    if not stamp:
        return

    for target_ref in range(start_ref, end_ref):
        if target_ref >= array_size:
            continue

        read_count = index.ref_read_counts[target_ref]
        if read_count < low_coverage_threshold:
            connection_counts[target_ref] = 0
            if neighbor_multimap_avg:
                neighbor_multimap_avg[target_ref] = 0.0
            if neighbor_connections_avg:
                neighbor_connections_avg[target_ref] = 0.0
            if neighbor_counts:
                neighbor_counts[target_ref] = 0
            continue

        connection_count = 0
        neighbor_multimap_sum = 0.0
        neighbor_connections_sum = 0.0
        gen = target_ref + 1 

        for i in range(read_count):
            read_idx = index.ref_to_reads[target_ref][i]
            if read_idx >= pool.unique_read_count:
                continue

            if read_refs_idx:
                ref_count = read_refs_idx.counts[read_idx]
                if ref_count < 1:
                    continue

                refs_ptr = read_refs_idx.read_ptrs[read_idx]
                for j in range(ref_count):
                    other_ref = refs_ptr[j]
                    if other_ref != target_ref and other_ref < array_size:

                        if stamp[other_ref] != gen:
                            stamp[other_ref] = gen
                            connection_count += 1

                            if neighbor_multimap_avg and dataset_multimap_fractions:
                                neighbor_multimap_sum += dataset_multimap_fractions[other_ref]
            else:

                ref_count = pool.read_alignment_counts[read_idx]
                if ref_count < 1:
                    continue

                start_pos = pool.read_alignment_starts[read_idx]
                for j in range(ref_count):
                    other_ref = pool.alignments[start_pos + j].reference_index
                    if other_ref != target_ref and other_ref < array_size:
                        if stamp[other_ref] != gen:
                            stamp[other_ref] = gen
                            connection_count += 1

                            if neighbor_multimap_avg and dataset_multimap_fractions:
                                neighbor_multimap_sum += dataset_multimap_fractions[other_ref]

        connection_counts[target_ref] = connection_count
        if neighbor_counts:
            neighbor_counts[target_ref] = connection_count

        if connection_count > 0:
            if neighbor_multimap_avg:
                neighbor_multimap_avg[target_ref] = neighbor_multimap_sum / connection_count
            if neighbor_connections_avg:
                neighbor_connections_avg[target_ref] = neighbor_connections_sum / connection_count
        else:
            if neighbor_multimap_avg:
                neighbor_multimap_avg[target_ref] = 0.0
            if neighbor_connections_avg:
                neighbor_connections_avg[target_ref] = 0.0

cdef int connection_counting_with_index(ReadIndex* read_index, ReadRefsIndex* read_refs_idx,
                                       uint64_t array_size, uint32_t* exact_connection_counts,
                                       int num_threads, MemoryPool* pool,
                                       double* neighbor_multimap_avg,
                                       double* neighbor_connections_avg, 
                                       uint32_t* neighbor_counts,
                                       float* dataset_multimap_fractions,
                                       uint32_t* dataset_connections,
                                       uint32_t low_coverage_threshold) nogil:
    # Parallel wrapper that splits the reference range and calls compute_reference_neighbor_metrics
    # Args:
    #     read_index/read_refs_idx: indexes built by build_read_index_parallel/build_read_refs_index
    #     array_size: number of references
    #     exact_connection_counts: output array to store per-ref neighbor counts
    #     num_threads: number of threads to use
    #     pool: MemoryPool pointer
    # Returns:
    #     0 on success, -1 on allocation failure

    cdef int thread_id, chunk_id
    cdef uint32_t chunk_size = (array_size + num_threads - 1) / num_threads
    cdef uint32_t start_ref, end_ref
    cdef uint32_t** thread_stamps = NULL

    thread_stamps = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    if not thread_stamps:
        return -1

    for thread_id in range(num_threads):
        thread_stamps[thread_id] = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        if not thread_stamps[thread_id]:

            for cleanup_id in range(thread_id):
                free(thread_stamps[cleanup_id])
            free(thread_stamps)
            return -1

    for chunk_id in prange(num_threads, nogil=True, num_threads=num_threads, schedule='static'):
        start_ref = chunk_id * chunk_size
        end_ref = min_int32(start_ref + chunk_size, array_size)

        compute_reference_neighbor_metrics(start_ref, end_ref, read_index, read_refs_idx,
                                                  array_size, exact_connection_counts, pool,
                                                  thread_stamps[chunk_id],
                                                  neighbor_multimap_avg, neighbor_connections_avg,
                                                  neighbor_counts, dataset_multimap_fractions,
                                                  dataset_connections, low_coverage_threshold)

    for thread_id in range(num_threads):
        free(thread_stamps[thread_id])
    free(thread_stamps)

    return 0


cdef void compute_neighbor_connections_avg_from_exact(uint32_t start_ref, uint32_t end_ref,
                                                     ReadIndex* index, ReadRefsIndex* read_refs_idx,
                                                     uint32_t array_size, double* neighbor_connections_avg,
                                                     uint32_t* exact_connection_counts, MemoryPool* pool,
                                                     uint32_t* stamp) noexcept nogil:
    # For each reference in [start_ref, end_ref), compute the average of exact_connection_counts
    # over its unique neighbors. This is run after `exact_connection_counts` has been fully populated.
    cdef uint32_t target_ref, read_count, read_idx, i
    cdef uint32_t ref_count, other_ref
    cdef uint32_t gen
    cdef uint32_t j
    cdef uint32_t neighbor_count
    cdef double neighbor_sum
    cdef uint32_t* refs_ptr

    if not stamp or not neighbor_connections_avg or not exact_connection_counts:
        return

    for target_ref in range(start_ref, end_ref):
        if target_ref >= array_size:
            continue

        read_count = index.ref_read_counts[target_ref]
        if read_count == 0:
            neighbor_connections_avg[target_ref] = 0.0
            continue

        neighbor_sum = 0.0
        neighbor_count = 0
        gen = target_ref + 1

        for i in range(read_count):
            read_idx = index.ref_to_reads[target_ref][i]
            if read_idx >= pool.unique_read_count:
                continue

            if read_refs_idx:
                ref_count = read_refs_idx.counts[read_idx]
                if ref_count < 1:
                    continue
                refs_ptr = read_refs_idx.read_ptrs[read_idx]
                for j in range(ref_count):
                    other_ref = refs_ptr[j]
                    if other_ref != target_ref and other_ref < array_size:
                        if stamp[other_ref] != gen:
                            stamp[other_ref] = gen
                            neighbor_sum += <double>exact_connection_counts[other_ref]
                            neighbor_count += 1
            else:
                ref_count = pool.read_alignment_counts[read_idx]
                if ref_count < 1:
                    continue
                for j in range(ref_count):
                    other_ref = pool.alignments[pool.read_alignment_starts[read_idx] + j].reference_index
                    if other_ref != target_ref and other_ref < array_size:
                        if stamp[other_ref] != gen:
                            stamp[other_ref] = gen
                            neighbor_sum += <double>exact_connection_counts[other_ref]
                            neighbor_count += 1

        if neighbor_count > 0:
            neighbor_connections_avg[target_ref] = neighbor_sum / <double>neighbor_count
        else:
            neighbor_connections_avg[target_ref] = 0.0


cdef int connection_neighbor_connections_pass(ReadIndex* read_index, ReadRefsIndex* read_refs_idx,
                                              uint64_t array_size, double* neighbor_connections_avg,
                                              int num_threads, uint32_t* exact_connection_counts,
                                              MemoryPool* pool) nogil:
    cdef uint32_t chunk_size = (array_size + num_threads - 1) / num_threads
    cdef uint32_t start_ref, end_ref
    cdef uint32_t** thread_stamps = NULL
    cdef int thread_id

    thread_stamps = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    if not thread_stamps:
        return -1

    for thread_id in range(num_threads):
        thread_stamps[thread_id] = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        if not thread_stamps[thread_id]:
            for cleanup_id in range(thread_id):
                free(thread_stamps[cleanup_id])
            free(thread_stamps)
            return -1

    for thread_id in prange(num_threads, nogil=True, num_threads=num_threads, schedule='static'):
        start_ref = thread_id * chunk_size
        end_ref = min_int32(start_ref + chunk_size, array_size)
        compute_neighbor_connections_avg_from_exact(start_ref, end_ref, read_index, read_refs_idx,
                                                    array_size, neighbor_connections_avg,
                                                    exact_connection_counts, pool,
                                                    thread_stamps[thread_id])

    for thread_id in range(num_threads):
        free(thread_stamps[thread_id])
    free(thread_stamps)

    return 0


cdef void destroy_read_index(ReadIndex* index) noexcept nogil:
    # Free resources allocated by build_read_index_parallel.
    if not index:
        return
    if index.read_buffer:
        free(index.read_buffer)
    if index.ref_to_reads:
        free(index.ref_to_reads)
    if index.ref_read_counts:
        free(index.ref_read_counts)
    free(index)


# Helper functions for per-read counting/filling (top-level, nogil)
cdef uint32_t count_unique_refs_for_read(MemoryPool* pool, uint32_t read_idx, uint32_t array_size) noexcept nogil:
    # Count the number of unique reference IDs a single read maps to (bounded by array_size).
    # Args:
    #     pool: MemoryPool pointer
    #     read_idx: index of the read
    #     array_size: upper bound of valid reference indices
    # Returns:
    #     deduplicated count of unique references for the read

    cdef uint32_t ref_count = pool.read_alignment_counts[read_idx]
    cdef uint64_t start_pos, end_pos
    cdef uint32_t* tmp = NULL
    cdef uint32_t k = 0
    cdef uint32_t i
    cdef uint32_t dedup_count

    if ref_count < 1:
        return ref_count

    start_pos = pool.read_alignment_starts[read_idx]
    end_pos = start_pos + ref_count

    tmp = <uint32_t*>malloc(ref_count * sizeof(uint32_t))
    if not tmp:
        return 0

    for i in range(start_pos, end_pos):
        if pool.alignments[i].reference_index < array_size:
            tmp[k] = pool.alignments[i].reference_index
            k += 1

    if k == 0:
        free(tmp)
        return 0

    qsort(tmp, k, sizeof(uint32_t), _uint32_compare)
    dedup_count = 1
    for i in range(1, k):
        if tmp[i] != tmp[dedup_count - 1]:
            tmp[dedup_count] = tmp[i]
            dedup_count += 1

    free(tmp)
    return dedup_count


cdef void fill_unique_refs_for_read(ReadRefsIndex* rri, MemoryPool* pool, uint32_t read_idx, uint32_t array_size) noexcept nogil:
    # Fill the ReadRefsIndex entry for a single read with the unique reference IDs it maps to.
    # Args:
    #     rri: ReadRefsIndex with preallocated buffers
    #     pool: MemoryPool pointer
    #     read_idx: index of the read to process
    #     array_size: upper bound of valid reference indices

    cdef uint32_t ref_count = pool.read_alignment_counts[read_idx]
    cdef uint64_t start_pos, end_pos
    cdef uint32_t* tmp = NULL
    cdef uint32_t k = 0
    cdef uint32_t i, dedup_count

    if rri.counts[read_idx] == 0:
        return

    if ref_count < 1:
        if ref_count == 1:
            start_pos = pool.read_alignment_starts[read_idx]
            i = pool.alignments[start_pos].reference_index
            if i < array_size:
                rri.read_ptrs[read_idx][0] = i
        return

    start_pos = pool.read_alignment_starts[read_idx]
    end_pos = start_pos + ref_count

    tmp = <uint32_t*>malloc(ref_count * sizeof(uint32_t))
    if not tmp:
        return

    for i in range(start_pos, end_pos):
        if pool.alignments[i].reference_index < array_size:
            tmp[k] = pool.alignments[i].reference_index
            k += 1

    if k == 0:
        free(tmp)
        return

    qsort(tmp, k, sizeof(uint32_t), _uint32_compare)

    dedup_count = 0
    for i in range(k):
        if dedup_count == 0 or tmp[i] != rri.read_ptrs[read_idx][dedup_count - 1]:
            rri.read_ptrs[read_idx][dedup_count] = tmp[i]
            dedup_count += 1

    if dedup_count < rri.counts[read_idx]:
        rri.counts[read_idx] = dedup_count

    free(tmp)

cdef ReadRefsIndex* build_read_refs_index(MemoryPool* pool, uint32_t array_size, int num_threads) nogil:
    # Build an index of unique references per read (deduplicated) in parallel.
    # Args:
    #     pool: MemoryPool pointer
    #     array_size: upper bound for valid reference indices
    #     num_threads: number of threads to use
    # Returns:
    #     pointer to allocated ReadRefsIndex or NULL on failure

    cdef uint32_t read_idx, alignment_idx, ref_idx
    cdef uint64_t start_pos, end_pos
    cdef uint32_t ref_count
    cdef uint32_t* unique_counts = NULL
    cdef uint32_t i
    cdef uint64_t total_unique = 0
    cdef ReadRefsIndex* rri = NULL
    cdef uint32_t* offsets = NULL
    cdef uint32_t* tmp = NULL
    cdef uint32_t k, dedup_count

    if not pool:
            return NULL

    rri = <ReadRefsIndex*>malloc(sizeof(ReadRefsIndex))
    if not rri:
        return NULL
    rri.buffer = NULL
    rri.read_ptrs = NULL
    rri.counts = NULL
    rri.read_count = pool.unique_read_count

    unique_counts = <uint32_t*>calloc(rri.read_count, sizeof(uint32_t))
    if not unique_counts:
        free(rri)
        return NULL

    cdef int nt = num_threads
    for read_idx in prange(rri.read_count, nogil=True, num_threads=nt):
        unique_counts[read_idx] = count_unique_refs_for_read(pool, read_idx, array_size)

    offsets = <uint32_t*>malloc(rri.read_count * sizeof(uint32_t))
    if not offsets:
        free(unique_counts)
        free(rri)
        return NULL

    total_unique = 0
    for i in range(rri.read_count):
        offsets[i] = <uint32_t>total_unique
        total_unique += unique_counts[i]

    if total_unique == 0:
        free(offsets)
        free(unique_counts)
        free(rri)
        return NULL

    rri.buffer = <uint32_t*>malloc(total_unique * sizeof(uint32_t))
    rri.read_ptrs = <uint32_t**>malloc(rri.read_count * sizeof(uint32_t*))
    rri.counts = <uint32_t*>malloc(rri.read_count * sizeof(uint32_t))
    if not rri.buffer or not rri.read_ptrs or not rri.counts:
        if rri.buffer: free(rri.buffer)
        if rri.read_ptrs: free(rri.read_ptrs)
        if rri.counts: free(rri.counts)
        free(offsets)
        free(unique_counts)
        free(rri)
        return NULL

    for i in range(rri.read_count):
        rri.counts[i] = unique_counts[i]
        rri.read_ptrs[i] = rri.buffer + offsets[i]

    for read_idx in prange(rri.read_count, nogil=True, num_threads=nt):
        fill_unique_refs_for_read(rri, pool, read_idx, array_size)

    free(offsets)
    free(unique_counts)

    return rri


cdef void destroy_read_refs_index(ReadRefsIndex* rri) noexcept nogil:
    # Free memory allocated for a ReadRefsIndex
    if not rri:
        return
    if rri.buffer: free(rri.buffer)
    if rri.read_ptrs: free(rri.read_ptrs)
    if rri.counts: free(rri.counts)
    free(rri)

cdef WeightedGraph* analyze_reference_graph(MemoryPool* pool, ReferencePattern* pattern_data,
                                           int32_t min_read_count, EMAlgorithmConfig* config,
                                           sam_hdr_t* bam_header, ReferenceMapping* mapping,
                                           bint verbose, bint build_igraph, const char* tsv_export_path,
                                           uint32_t graph_min_edge_weight, TaxonomyDB* taxonomy_db) noexcept nogil:
    """High-level reference graph analysis pipeline.

    Performs a multi-phase analysis that converts alignment data in ``pool`` into
    per-reference statistics, read/reference indices, and (optionally) a filtered
    weighted graph suitable for clustering and TSV export.

    This function is optimized for nogil parallel execution and carefully manages
    temporary memory to reduce peak usage on large datasets. It does not perform
    any I/O itself except when TSV export is requested via ``tsv_export_path``.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments and per-read indices.
    pattern_data : ReferencePattern*
        Preallocated array to be filled with per-reference graph metrics.
    min_read_count : int32_t
        Minimum reads required for a reference to be considered in the graph.
    config : EMAlgorithmConfig*
        Configuration controlling threading and algorithm options.
    bam_header : sam_hdr_t*
        Optional BAM header for reference name/length lookups (may be NULL).
    mapping : ReferenceMapping*
        Optional mapping from compacted retained reference IDs to original TIDs.
    verbose : bint
        If True, emit progress and diagnostic messages.
    build_igraph : bint
        If True, construct and return a filtered WeightedGraph for downstream use.
    tsv_export_path : const char*
        If non-NULL, path where per-reference TSV summaries will be written.
    graph_min_edge_weight : uint32_t
        Minimum shared-read count to keep an edge when building the filtered graph.

    Returns
    -------
    WeightedGraph*
        Pointer to a filtered WeightedGraph when ``build_igraph`` is True and
        construction succeeded. Returns NULL on failure or when no graph is built.
    """

    cdef int num_threads = config.thread_count if config.thread_count > 0 else 1
    cdef uint32_t ref_idx, chunk_id
    cdef int thread_id
    cdef ReferenceStats* ref_stats = NULL
    cdef uint32_t low_coverage_threshold = <uint32_t>min_read_count
    cdef uint32_t suspicious_count = 0
    cdef uint32_t max_ref_id = 0
    cdef int64_t i
    cdef uint32_t array_size

    cdef double* neighbor_multimap_avg = NULL
    cdef double* neighbor_connections_avg = NULL
    cdef uint32_t* neighbor_counts = NULL

    cdef float* dataset_multimap_fractions = NULL
    cdef double* dataset_avg_co_mappings = NULL
    cdef uint32_t* dataset_connections = NULL
    cdef uint32_t* dataset_read_counts = NULL

    cdef uint32_t* total_reads = NULL
    cdef uint32_t* multimap_reads = NULL
    cdef uint64_t* alignments_per_ref = NULL
    cdef uint32_t** thread_total_reads = NULL
    cdef uint32_t** thread_multimap_reads = NULL
    cdef uint64_t** thread_alignments_per_ref = NULL

    cdef double* co_mapping_sums = NULL
    cdef uint32_t* co_mapping_counts = NULL
    cdef uint64_t* max_co_mappings = NULL
    cdef double* co_mapping_squares = NULL
    cdef double** thread_co_mapping_sums = NULL
    cdef uint32_t** thread_co_mapping_counts = NULL
    cdef uint64_t** thread_max_co_mappings = NULL
    cdef double** thread_co_mapping_squares = NULL

    cdef uint64_t chunk_size_alignments, start_aln, end_aln
    cdef uint64_t chunk_size_reads
    cdef uint32_t start_read, end_read
    cdef uint32_t chunk
    cdef uint32_t tid
    cdef ReadIndex* read_index = NULL
    cdef ReadRefsIndex* read_refs_idx = NULL
    cdef uint32_t* exact_connection_counts = NULL

    cdef float multimap_fraction
    cdef GraphMetrics exact_metrics
    cdef double mean_co_mappings, variance, std_dev
    cdef uint32_t exact_connections
    cdef uint32_t patterns_computed = 0
    cdef double* co_mapping_averages = NULL
    cdef double dataset_median_connections

    cdef timespec ts_phase_start, ts_phase_end
    # Declarations required for optional igraph construction (must be at function scope)
    cdef ReferenceStats* ref_stats_for_graph = NULL
    cdef igraph_t ig_graph
    cdef igraph_vector_t ig_weights
    cdef WeightedGraph* filtered_graph = NULL
    cdef igraph_t* igraph_ptr = NULL
    cdef igraph_vector_t* weights_ptr = NULL
    # Connected nodes optimization variables
    cdef uint32_t* connected_node_ids = NULL
    cdef uint32_t n_connected = 0
    cdef uint32_t n_isolated = 0
    cdef igraph_vector_int_t degree_vec
    cdef igraph_integer_t node_degree
    cdef int degree_ret
    cdef int build_result
    # igraph component vectors (declare at function scope to satisfy Cython)
    cdef igraph_vector_int_t comp_membership
    cdef igraph_vector_int_t comp_sizes
    cdef igraph_integer_t n_components = 0
    # builder-specific component counts (returned by direct igraph builder)
    cdef igraph_integer_t builder_n_components = 0
    cdef igraph_integer_t builder_singletons = 0
    cdef int cret
    cdef long i_comp
    cdef long singletons
    # Track whether comp vectors were successfully initialized (can't use try/except nogil)
    cdef int comp_membership_inited = 0
    cdef int comp_sizes_inited = 0

    if pool.reference_count == 0:
        bf_nogil_logf_notime(NULL, "ERROR: No references present in memory pool (reference_count=0)\n")
        return NULL

    # By default use the MemoryPool reference_count. However, after
    # remapping (remap_alignment_reference_ids) the compacted reference
    # space lives in `mapping->n_retained_refs`. If a mapping was provided
    # and contains a smaller retained set, prefer that compact size so
    # subsequent indexes (ReadIndex / ReadRefsIndex) are allocated using
    # the correct array bounds and counts remain consistent with the
    # remapped alignment reference IDs.
    array_size = pool.reference_count
    if mapping is not NULL and mapping.n_retained_refs > 0:
        array_size = mapping.n_retained_refs

    cdef uint64_t invalid_refs = 0
    for i in range(pool.alignment_count):
        if pool.alignments[i].reference_index >= array_size:
            invalid_refs += 1

    if invalid_refs:
        bf_nogil_logf_notime(NULL, "ERROR: Found %lu alignments with invalid reference IDs (>= %u). Aborting.\n", invalid_refs, array_size)
        return NULL

    if verbose:
        bf_nogil_logf_notime(NULL, "NETWORK-AWARE ANALYSIS: Starting with %d threads (max_ref_id=%u)\n", 
                      num_threads, max_ref_id)
        # Only mention pruning / filtered-graph behavior when igraph/filtered graph will be
        # constructed (i.e., when clustering or igraph build is requested). For TSV-only
        # runs we compute statistics from the ReadIndex/ReadRefsIndex and do NOT prune the
        # underlying co-mapping graph used for the stats.
        if build_igraph:
            bf_nogil_logf_notime(NULL, "Graph min edge weight: %u (edges with < %u shared reads will be pruned)\n",
                          graph_min_edge_weight, graph_min_edge_weight)
            bf_nogil_logf_notime(NULL, "=== BUILDING FILTERED WEIGHTED GRAPH FIRST ===\n")
            bf_nogil_logf_notime(NULL, "All statistics will be calculated from this filtered graph\n")
        else:
            bf_nogil_logf_notime(NULL, "Graph min edge weight: %u (pruning will be skipped because igraph/building filtered graph is disabled)\n",
                          graph_min_edge_weight)
            bf_nogil_logf_notime(NULL, "=== BUILDING FILTERED WEIGHTED GRAPH SKIPPED (TSV/stats-only path) ===\n")
            bf_nogil_logf_notime(NULL, "All statistics will be calculated from the ReadIndex/ReadRefsIndex (no graph pruning)\n")

    total_reads = <uint32_t*>calloc(array_size, sizeof(uint32_t))
    multimap_reads = <uint32_t*>calloc(array_size, sizeof(uint32_t))
    alignments_per_ref = <uint64_t*>calloc(array_size, sizeof(uint64_t))
    neighbor_multimap_avg = <double*>calloc(array_size, sizeof(double))
    neighbor_connections_avg = <double*>calloc(array_size, sizeof(double))
    neighbor_counts = <uint32_t*>calloc(array_size, sizeof(uint32_t))

    thread_total_reads = <uint32_t**>calloc(num_threads, sizeof(uint32_t*))
    thread_multimap_reads = <uint32_t**>calloc(num_threads, sizeof(uint32_t*))
    thread_alignments_per_ref = <uint64_t**>calloc(num_threads, sizeof(uint64_t*))

    if (not total_reads or not multimap_reads or not alignments_per_ref or
        not neighbor_multimap_avg or not neighbor_connections_avg or not neighbor_counts or
        not thread_total_reads or not thread_multimap_reads or not thread_alignments_per_ref):
        bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate main arrays\n")

        if total_reads: free(total_reads)
        if multimap_reads: free(multimap_reads)
        if alignments_per_ref: free(alignments_per_ref)
        if neighbor_multimap_avg: free(neighbor_multimap_avg)
        if neighbor_connections_avg: free(neighbor_connections_avg)
        if neighbor_counts: free(neighbor_counts)
        if thread_total_reads: free(thread_total_reads)
        if thread_multimap_reads: free(thread_multimap_reads)
        if thread_alignments_per_ref: free(thread_alignments_per_ref)
        return NULL

    for thread_id in range(num_threads):
        thread_total_reads[thread_id] = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        thread_multimap_reads[thread_id] = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        thread_alignments_per_ref[thread_id] = <uint64_t*>calloc(array_size, sizeof(uint64_t))
        if (not thread_total_reads[thread_id] or not thread_multimap_reads[thread_id] or 
            not thread_alignments_per_ref[thread_id]):
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate thread arrays\n")

            for cleanup_id in range(num_threads):
                if thread_total_reads[cleanup_id]: free(thread_total_reads[cleanup_id])
                if thread_multimap_reads[cleanup_id]: free(thread_multimap_reads[cleanup_id])
                if thread_alignments_per_ref[cleanup_id]: free(thread_alignments_per_ref[cleanup_id])
            free(total_reads); free(multimap_reads); free(alignments_per_ref)
            free(neighbor_multimap_avg); free(neighbor_connections_avg); free(neighbor_counts)
            free(thread_total_reads); free(thread_multimap_reads); free(thread_alignments_per_ref)
            return NULL

    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 1: Parallel alignment counting...\n")
    clock_gettime(CLOCK_MONOTONIC, &ts_phase_start)
    chunk_size_alignments = (pool.alignment_count + num_threads - 1) / num_threads

    for chunk_id in prange(num_threads, nogil=True, num_threads=num_threads):
        start_aln = chunk_id * chunk_size_alignments
        end_aln = min_int64(start_aln + chunk_size_alignments, pool.alignment_count)
        count_alignments_chunk(start_aln, end_aln, pool, thread_alignments_per_ref[chunk_id])

    for thread_id in range(num_threads):
        for ref_idx in range(array_size):
            alignments_per_ref[ref_idx] += thread_alignments_per_ref[thread_id][ref_idx]

    clock_gettime(CLOCK_MONOTONIC, &ts_phase_end)
    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 1: Completed in %.3f seconds\n", calculate_elapsed_seconds(&ts_phase_start, &ts_phase_end))

    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 2: Parallel read processing...\n")
    clock_gettime(CLOCK_MONOTONIC, &ts_phase_start)
    chunk_size_reads = (pool.final_unique_reads + num_threads - 1) / num_threads

    for chunk_id in prange(num_threads, nogil=True, num_threads=num_threads):
        start_read = chunk_id * chunk_size_reads
        end_read = min_int32(start_read + chunk_size_reads, pool.final_unique_reads)
        count_reads_chunk(start_read, end_read, pool, thread_total_reads[chunk_id], 
                        thread_multimap_reads[chunk_id])

    reduce_thread_arrays(thread_total_reads, total_reads, array_size, num_threads)
    reduce_thread_arrays(thread_multimap_reads, multimap_reads, array_size, num_threads)

    clock_gettime(CLOCK_MONOTONIC, &ts_phase_end)
    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 2: Completed in %.3f seconds\n", calculate_elapsed_seconds(&ts_phase_start, &ts_phase_end))

    # PHASE 4 statistics accumulation was previously done via read-centric co-mapping
    # accumulators. Those temporary arrays used a lot of memory. We now compute all
    # neighbor/co-mapping statistics from the filtered igraph (see below). Skip the
    # pre-graph co-mapping accumulation to reduce peak memory.
    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 4: Skipping pre-graph co-mapping accumulation (computed from igraph instead)\n")

    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 5: Building read index and computing neighbor quality...\n")
    clock_gettime(CLOCK_MONOTONIC, &ts_phase_start)

    read_index = build_read_index_parallel(pool, array_size, num_threads)
    # NOTE: we avoid building ReadRefsIndex here to reduce peak memory.
    # read_refs_idx = build_read_refs_index(pool, array_size, num_threads)
    if not read_index:
        bf_nogil_logf_notime(NULL, "ERROR: Failed to build read index\n")
        return NULL

    # Build ReadRefsIndex (deduplicated per-read ref lists) for accurate co-mapping stats
    # This is memory-efficient compared to the older per-thread temporary accumulators
    read_refs_idx = build_read_refs_index(pool, array_size, num_threads)
    if not read_refs_idx:
        bf_nogil_logf_notime(NULL, "WARNING: Failed to build ReadRefsIndex; skipping exact co-mapping accumulation\n")
    else:
        # Allocate co-mapping accumulators
        co_mapping_sums = <double*>calloc(array_size, sizeof(double))
        co_mapping_counts = <uint32_t*>calloc(array_size, sizeof(uint32_t))
        max_co_mappings = <uint64_t*>calloc(array_size, sizeof(uint64_t))
        co_mapping_squares = <double*>calloc(array_size, sizeof(double))

        if (not co_mapping_sums or not co_mapping_counts or not max_co_mappings or not co_mapping_squares):
            bf_nogil_logf_notime(NULL, "WARNING: Failed to allocate co-mapping accumulators; skipping exact accumulation\n")
            if co_mapping_sums: free(co_mapping_sums); co_mapping_sums = NULL
            if co_mapping_counts: free(co_mapping_counts); co_mapping_counts = NULL
            if max_co_mappings: free(max_co_mappings); max_co_mappings = NULL
            if co_mapping_squares: free(co_mapping_squares); co_mapping_squares = NULL
        else:
            # Parallel accumulation from ReadRefsIndex
            chunk = (read_refs_idx.read_count + num_threads - 1) // num_threads
            for tid in prange(num_threads, nogil=True, num_threads=num_threads):
                start_read = tid * chunk
                end_read = min_int32(start_read + chunk, read_refs_idx.read_count)
                accumulate_co_mappings_from_readrefs(start_read, end_read, pool, read_refs_idx,
                                                    co_mapping_sums, co_mapping_counts,
                                                    max_co_mappings, co_mapping_squares, array_size)

            # Compute co_mapping_averages (average per read across all reads mapped to ref)
            co_mapping_averages = <double*>malloc(array_size * sizeof(double))
            if co_mapping_averages:
                for ref_idx in range(array_size):
                    if total_reads and total_reads[ref_idx] > 0:
                        co_mapping_averages[ref_idx] = co_mapping_sums[ref_idx] / <double>total_reads[ref_idx]
                    else:
                        co_mapping_averages[ref_idx] = 0.0

    exact_connection_counts = <uint32_t*>calloc(array_size, sizeof(uint32_t))
    if not exact_connection_counts:
        bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate connection counts array\n")
        destroy_read_index(read_index)
        return NULL

    dataset_multimap_fractions = <float*>malloc(array_size * sizeof(float))
    dataset_avg_co_mappings = <double*>malloc(array_size * sizeof(double))
    dataset_connections = <uint32_t*>malloc(array_size * sizeof(uint32_t))
    dataset_read_counts = <uint32_t*>malloc(array_size * sizeof(uint32_t))
    if (not dataset_multimap_fractions or not dataset_avg_co_mappings or 
        not dataset_connections or not dataset_read_counts):
        bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate dataset arrays\n")
        return NULL

    # Populate dataset-level arrays from read counts we already have. Average
    # co-mappings per ref will be computed later from the filtered igraph if
    # required; avoid allocating large co-mapping accumulators here.
    for i in range(array_size):
        if total_reads[i] > 0:
            dataset_multimap_fractions[i] = <float>multimap_reads[i] / <float>total_reads[i]
            dataset_avg_co_mappings[i] = 0.0
            dataset_connections[i] = 0 
            dataset_read_counts[i] = total_reads[i]
        else:
            dataset_multimap_fractions[i] = 0.0
            dataset_avg_co_mappings[i] = 0.0
            dataset_connections[i] = 0
            dataset_read_counts[i] = 0

    # Compute exact connection counts (neighbor counts) from the ReadIndex / ReadRefsIndex
    # This populates `exact_connection_counts` and neighbor arrays so TSV and later
    # logic has accurate per-ref neighbor metrics even when igraph is not built.
    if read_index:
        if connection_counting_with_index(read_index, read_refs_idx, array_size,
                                          exact_connection_counts, num_threads, pool,
                                          neighbor_multimap_avg, neighbor_connections_avg,
                                          neighbor_counts, dataset_multimap_fractions,
                                          dataset_connections, low_coverage_threshold) != 0:
            bf_nogil_logf_notime(NULL, "WARNING: connection_counting_with_index failed; neighbor metrics may be incomplete\n")
        else:
            # Compute neighbor connections averages in a dedicated second pass
            if neighbor_connections_avg:
                if connection_neighbor_connections_pass(read_index, read_refs_idx, array_size,
                                                        neighbor_connections_avg, num_threads,
                                                        exact_connection_counts, pool) != 0:
                    bf_nogil_logf_notime(NULL, "WARNING: connection_neighbor_connections_pass failed; neighbor connection averages incomplete\n")

    # Free large temporary buffers that are no longer needed before building igraph
    # This reduces peak memory pressure for large datasets.
    for thread_id in range(num_threads):
        if thread_total_reads and thread_total_reads[thread_id]:
            free(thread_total_reads[thread_id])
        if thread_multimap_reads and thread_multimap_reads[thread_id]:
            free(thread_multimap_reads[thread_id])
        if thread_alignments_per_ref and thread_alignments_per_ref[thread_id]:
            free(thread_alignments_per_ref[thread_id])
    if thread_total_reads:
        free(thread_total_reads)
        thread_total_reads = NULL
    if thread_multimap_reads:
        free(thread_multimap_reads)
        thread_multimap_reads = NULL
    if thread_alignments_per_ref:
        free(thread_alignments_per_ref)
        thread_alignments_per_ref = NULL

    # Note: co_mapping_* arrays are intentionally NOT freed here because they
    # may have been populated from ReadRefsIndex and are needed later by the TSV writer.
    # Freeing of these arrays happens at function cleanup near the end.

    if thread_co_mapping_sums:
        for thread_id in range(num_threads):
            if thread_co_mapping_sums[thread_id]: free(thread_co_mapping_sums[thread_id])
        free(thread_co_mapping_sums)
        thread_co_mapping_sums = NULL
    if thread_co_mapping_counts:
        for thread_id in range(num_threads):
            if thread_co_mapping_counts[thread_id]: free(thread_co_mapping_counts[thread_id])
        free(thread_co_mapping_counts)
        thread_co_mapping_counts = NULL
    if thread_max_co_mappings:
        for thread_id in range(num_threads):
            if thread_max_co_mappings[thread_id]: free(thread_max_co_mappings[thread_id])
        free(thread_max_co_mappings)
        thread_max_co_mappings = NULL
    if thread_co_mapping_squares:
        for thread_id in range(num_threads):
            if thread_co_mapping_squares[thread_id]: free(thread_co_mapping_squares[thread_id])
        free(thread_co_mapping_squares)
        thread_co_mapping_squares = NULL

    # If igraph build is not requested we can free read_refs_idx after producing co-mapping arrays
    if not build_igraph:
        if read_refs_idx:
            destroy_read_refs_index(read_refs_idx)
            read_refs_idx = NULL

    # Defer building the filtered/trimmed igraph until after we compute and
    # print the reference statistics summary. For now populate dataset_connections
    # from the exact_connection_counts (ReadIndex-derived) so summary uses the
    # unpruned graph metrics.
    for ref_idx in range(array_size):
        dataset_connections[ref_idx] = exact_connection_counts[ref_idx]

    clock_gettime(CLOCK_MONOTONIC, &ts_phase_end)
    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 5: Completed in %.3f seconds\n", calculate_elapsed_seconds(&ts_phase_start, &ts_phase_end))

    for ref_idx in range(pool.reference_count):
        # Initialize unique read count
        pattern_data[ref_idx].unique_read_count = 0

        # Initialize graph metrics
        memset(&pattern_data[ref_idx].graph, 0, sizeof(GraphMetrics))

        # Initialize topology
        pattern_data[ref_idx].component_id = ref_idx
        pattern_data[ref_idx].component_size = 1

        # Initialize Community fields with defaults (will be overwritten if clustering is run)
        pattern_data[ref_idx].community_id = 0
        pattern_data[ref_idx].community_cc = 0.0
        pattern_data[ref_idx].community_individual_cc = 0.0
        pattern_data[ref_idx].community_cc_threshold = 0.0
        pattern_data[ref_idx].community_keep_flag = 1  # default: keep all references

    dataset_median_connections = compute_dataset_median_connections(
        exact_connection_counts, total_reads, array_size, low_coverage_threshold)

    if verbose:
        bf_nogil_logf_notime(NULL, "PHASE 6: Computing network-aware patterns with neighbor quality...\n")
        bf_nogil_logf_notime(NULL, "Dataset median connections: %.3f\n", dataset_median_connections)

    for ref_idx in range(array_size):
        if ref_idx >= pool.reference_count:
            continue

        if total_reads[ref_idx] < low_coverage_threshold:
            continue

        # Get exact connection count computed earlier (from ReadIndex or igraph)
        exact_connections = exact_connection_counts[ref_idx]

        # Preserve any graph metrics populated by igraph; but ensure that the
        # basic connection_count and co-mapping metrics are reflected in the
        # pattern_data structure so on-screen summaries and downstream code
        # see consistent values regardless of the igraph path.
        pattern_data[ref_idx].graph.connection_count = exact_connections

        # Populate average co-mappings per read from precomputed array when present
        if co_mapping_averages:
            pattern_data[ref_idx].graph.avg_comappings_per_read = co_mapping_averages[ref_idx]

        # Populate max co-mappings observed when available
        if max_co_mappings:
            pattern_data[ref_idx].graph.max_comappings = max_co_mappings[ref_idx]

        # Number of reads that contributed to co-mapping stats
        if co_mapping_counts:
            pattern_data[ref_idx].graph.reads_with_comappings = co_mapping_counts[ref_idx]

        # unique_read_count will be populated from ref_stats after it's calculated
        patterns_computed += 1
    
    # Calculate full ref_stats only if needed for verbose/TSV (includes score stats, PMD, etc.)
    ref_stats = NULL
    if verbose or tsv_export_path:
        ref_stats = <ReferenceStats*>malloc(pool.reference_count * sizeof(ReferenceStats))
        if not ref_stats:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate reference statistics array\n")
        else:
            # Calculate all reference statistics including read categories and score stats
            calculate_reference_stats(pool, bam_header, ref_stats)
            
            # Populate pattern_data[].unique_read_count from ref_stats for filtering logic
            # This ensures consistency between TSV output and filtering decisions
            for ref_idx in range(pool.reference_count):
                pattern_data[ref_idx].unique_read_count = ref_stats[ref_idx].unique_reads

    if verbose:
        bf_nogil_logf_notime(NULL, "NETWORK-AWARE ANALYSIS: Completed %u patterns with neighbor quality metrics\n", patterns_computed)
        print_pattern_summary(pool, bam_header, mapping, pattern_data, NULL, ref_stats)

    # After printing the summary and computing stats from ReadIndex-derived
    # arrays, optionally build the trimmed/filtered igraph for clustering
    # and cache it for later Community operations. Building the igraph after the
    # summary ensures the reported metrics are not affected by pruning.
    if build_igraph:
        if verbose:
            bf_nogil_logf_notime(NULL, "Building igraph from ReadIndex (min_edge_weight=%u)...\n", graph_min_edge_weight)

        ref_stats_for_graph = <ReferenceStats*>calloc(array_size, sizeof(ReferenceStats))
        if not ref_stats_for_graph:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate ref_stats for graph\n")
            return NULL

        for ref_idx in range(array_size):
            ref_stats_for_graph[ref_idx].total_reads = total_reads[ref_idx]

        # Use direct builder to minimize peak memory (avoid WeightedGraph allocation)
        builder_n_components = 0
        builder_singletons = 0
        build_result = build_igraph_direct_from_read_index(
            &ig_graph, &ig_weights,
            pool, read_index, ref_stats_for_graph, array_size,
            min_read_count, graph_min_edge_weight, num_threads, verbose,
            &builder_n_components, &builder_singletons)

        if build_result != 0:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to build igraph\n")
            free(ref_stats_for_graph)
            return NULL

        if verbose:
            bf_nogil_logf_notime(NULL, "igraph built successfully\n")

        # Store igraph and its weight vector on the heap so it can be reused by Community
        filtered_graph = <WeightedGraph*>malloc(sizeof(WeightedGraph))
        if not filtered_graph:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate WeightedGraph wrapper\n")
            igraph_destroy(&ig_graph)
            igraph_vector_destroy(&ig_weights)
            free(ref_stats_for_graph)
            return NULL

        filtered_graph.nodes = NULL
        filtered_graph.num_nodes = array_size
        filtered_graph.total_weight = 0
        filtered_graph.num_edges = 0

        igraph_ptr = <igraph_t*>malloc(sizeof(igraph_t))
        if not igraph_ptr:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate igraph pointer\n")
            free(filtered_graph)
            igraph_destroy(&ig_graph)
            igraph_vector_destroy(&ig_weights)
            free(ref_stats_for_graph)
            return NULL

        igraph_ptr[0] = ig_graph
        filtered_graph.igraph_handle = <void*>igraph_ptr

        weights_ptr = <igraph_vector_t*>malloc(sizeof(igraph_vector_t))
        if not weights_ptr:
            bf_nogil_logf_notime(NULL, "ERROR: Failed to allocate weights pointer\n")
            free(igraph_ptr)
            free(filtered_graph)
            igraph_destroy(&ig_graph)
            igraph_vector_destroy(&ig_weights)
            free(ref_stats_for_graph)
            return NULL

        weights_ptr[0] = ig_weights
        filtered_graph.weights_handle = <void*>weights_ptr

        if verbose:
            bf_nogil_logf_notime(NULL, "Cached igraph and edge weights for Community clustering\n")

        # Extract connected nodes (degree > 0) for Community optimization
        # This avoids processing 88k+ isolated nodes that become singleton communities
        # Variables declared at function scope (lines 1934-1939)

        # Get all node degrees at once
        degree_ret = igraph_vector_int_init(&degree_vec, array_size)
        if degree_ret == IGRAPH_SUCCESS:
            degree_ret = igraph_degree(&ig_graph, &degree_vec, igraph_vss_all(),
                                      IGRAPH_ALL, False)  # loops=False (don't count self-loops)
            if degree_ret == IGRAPH_SUCCESS:
                # Count connected vs isolated nodes
                for ref_idx in range(array_size):
                    node_degree = get_vector_int_element(&degree_vec, ref_idx)
                    if node_degree > 0:
                        n_connected += 1
                    else:
                        n_isolated += 1

                if verbose:
                    bf_nogil_logf_notime(NULL, "Graph connectivity: %u connected nodes, %u isolated nodes\n",
                                        n_connected, n_isolated)

                # Allocate and fill connected node list
                if n_connected > 0:
                    connected_node_ids = <uint32_t*>malloc(n_connected * sizeof(uint32_t))
                    if not connected_node_ids:
                        bf_nogil_logf_notime(NULL, "WARNING: Failed to allocate connected nodes list; will process all nodes\n")
                        n_connected = 0
                    else:
                        n_connected = 0  # Reset for filling
                        for ref_idx in range(array_size):
                            node_degree = get_vector_int_element(&degree_vec, ref_idx)
                            if node_degree > 0:
                                connected_node_ids[n_connected] = ref_idx
                                n_connected += 1

            igraph_vector_int_destroy(&degree_vec)

        filtered_graph.connected_node_ids = connected_node_ids
        filtered_graph.n_connected_nodes = n_connected

        free(ref_stats_for_graph)

        # Connected components were computed inside the builder; report them
        if verbose:
            bf_nogil_logf_notime(NULL, "Connected components: %ld (singletons: %ld)\n", <long>builder_n_components, <long>builder_singletons)

    # Only write TSV here if clustering is disabled. If clustering is enabled,
    # the TSV will be written after Community clustering completes (from apply_cluster_aware_filtering)
    # so that Community results are included in the output.
    if tsv_export_path and not build_igraph:
        if ref_stats:
            # Write TSV with graph analysis results
            if write_graph_tsv(pool, bam_header, mapping,
                                      pattern_data, ref_stats,
                                      total_reads, multimap_reads, alignments_per_ref,
                                      exact_connection_counts,
                                      co_mapping_averages, max_co_mappings, co_mapping_counts,
                                      neighbor_multimap_avg, neighbor_connections_avg,
                                      neighbor_counts, array_size, dataset_median_connections,
                                      min_read_count, build_igraph, 0, tsv_export_path,
                                      taxonomy_db) != 0:
                bf_nogil_logf_notime(NULL, "ERROR: Failed to write graph analysis TSV\n")
        else:
            bf_nogil_logf_notime(NULL, "ERROR: Reference stats not calculated for TSV export\n")
    
    if ref_stats:
        free(ref_stats)

    destroy_read_index(read_index)
    if read_refs_idx:
        destroy_read_refs_index(read_refs_idx)

    # If build_igraph is true, store TSV data in filtered_graph for later use
    # after Community clustering completes. Otherwise free immediately.
    if build_igraph and filtered_graph:
        filtered_graph.tsv_total_reads = total_reads
        filtered_graph.tsv_multimap_reads = multimap_reads
        filtered_graph.tsv_alignments_per_ref = alignments_per_ref
        filtered_graph.tsv_exact_connection_counts = exact_connection_counts
        filtered_graph.tsv_co_mapping_averages = co_mapping_averages
        filtered_graph.tsv_max_co_mappings = max_co_mappings
        filtered_graph.tsv_co_mapping_counts = co_mapping_counts
        filtered_graph.tsv_neighbor_multimap_avg = neighbor_multimap_avg
        filtered_graph.tsv_neighbor_connections_avg = neighbor_connections_avg
        filtered_graph.tsv_neighbor_counts = neighbor_counts
        filtered_graph.tsv_dataset_median_connections = dataset_median_connections
        filtered_graph.tsv_array_size = array_size
        filtered_graph.tsv_min_read_count = min_read_count
        # Intermediate arrays are still freed below
    else:
        # Free TSV arrays immediately if not building igraph
        free(total_reads)
        free(multimap_reads)
        free(alignments_per_ref)
        free(exact_connection_counts)
        free(co_mapping_averages)
        free(max_co_mappings)
        free(co_mapping_counts)
        free(neighbor_multimap_avg)
        free(neighbor_connections_avg)
        free(neighbor_counts)

    # Free intermediate arrays that are not needed for TSV
    free(co_mapping_sums)
    free(co_mapping_squares)
    free(dataset_multimap_fractions)
    free(dataset_avg_co_mappings)
    free(dataset_connections)
    free(dataset_read_counts)

    if thread_co_mapping_sums:
        for thread_id in range(num_threads):
            if thread_co_mapping_sums[thread_id]: free(thread_co_mapping_sums[thread_id])
        free(thread_co_mapping_sums)
    if thread_co_mapping_counts:
        for thread_id in range(num_threads):
            if thread_co_mapping_counts[thread_id]: free(thread_co_mapping_counts[thread_id])
        free(thread_co_mapping_counts)
    if thread_max_co_mappings:
        for thread_id in range(num_threads):
            if thread_max_co_mappings[thread_id]: free(thread_max_co_mappings[thread_id])
        free(thread_max_co_mappings)
    if thread_co_mapping_squares:
        for thread_id in range(num_threads):
            if thread_co_mapping_squares[thread_id]: free(thread_co_mapping_squares[thread_id])
        free(thread_co_mapping_squares)
    if thread_total_reads:
        for thread_id in range(num_threads):
            if thread_total_reads[thread_id]: free(thread_total_reads[thread_id])
        free(thread_total_reads)
    if thread_multimap_reads:
        for thread_id in range(num_threads):
            if thread_multimap_reads[thread_id]: free(thread_multimap_reads[thread_id])
        free(thread_multimap_reads)
    if thread_alignments_per_ref:
        for thread_id in range(num_threads):
            if thread_alignments_per_ref[thread_id]: free(thread_alignments_per_ref[thread_id])
        free(thread_alignments_per_ref)

    # Return the filtered graph so it can be reused by Community clustering
    # Caller is responsible for calling destroy_weighted_graph() when done
    return filtered_graph


cdef void accumulate_co_mappings(uint32_t start_read, uint32_t end_read,
                                 MemoryPool* pool,
                                 uint32_t* read_refs_buffer,
                                 uint32_t* connection_counts,
                                 double* co_mapping_sums,
                                 uint32_t* co_mapping_counts,
                                 uint64_t* max_co_mappings, 
                                 double* co_mapping_squares,
                                 uint32_t array_size) noexcept nogil:
    # For a range of reads, accumulate co-mapping statistics per reference.
    # Args:
    #     start_read, end_read: read index range
    #     pool: MemoryPool pointer
    #     read_refs_buffer: optional buffer to reuse for per-read refs (not required)
    #     connection_counts: optional per-ref connection counter
    #     co_mapping_sums: per-ref sum of co-mappings observed across reads
    #     co_mapping_counts: per-ref number of reads contributing to the sum
    #     max_co_mappings: per-ref maximum co-mappings observed
    #     co_mapping_squares: per-ref sum-of-squares for variance calculation
    #     array_size: number of references

    cdef uint32_t read_idx, alignment_idx, ref_idx, ref_count
    cdef uint64_t start_pos, end_pos
    cdef uint32_t i, collected_count
    cdef uint64_t co_mappings 
    cdef double co_mappings_f
    cdef uint32_t* temp_refs = NULL
    cdef uint32_t max_possible_refs, alloc_size
    cdef char* ref_seen = NULL
    cdef uint32_t* new_refs

    for read_idx in range(start_read, end_read):
        if read_idx >= pool.unique_read_count:
            continue

        ref_count = pool.read_alignment_counts[read_idx]
        if ref_count < 1:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + ref_count

        if end_pos > <uint64_t>pool.alignment_count:
            end_pos = <uint64_t>pool.alignment_count

        max_possible_refs = <uint32_t>(end_pos - start_pos)
        if max_possible_refs == 0:
            continue

        alloc_size = max_int32(1024, max_possible_refs)
        temp_refs = <uint32_t*>malloc(alloc_size * sizeof(uint32_t))
        if not temp_refs:
            continue

        ref_seen = <char*>calloc(array_size, sizeof(char))
        if not ref_seen:
            free(temp_refs)
            continue

        collected_count = 0
        for alignment_idx in range(start_pos, end_pos):
            ref_idx = pool.alignments[alignment_idx].reference_index
            if (ref_idx < array_size and not ref_seen[ref_idx]):

                if collected_count >= alloc_size:
                    alloc_size = alloc_size * 2 
                    new_refs = <uint32_t*>realloc(temp_refs, alloc_size * sizeof(uint32_t))
                    if not new_refs:

                        break
                    temp_refs = new_refs

                ref_seen[ref_idx] = 1
                temp_refs[collected_count] = ref_idx
                collected_count += 1

        free(ref_seen)

        if collected_count <= 1:
            free(temp_refs)
            continue

        if collected_count > 1:
            qsort(temp_refs, collected_count, sizeof(uint32_t), _uint32_compare)

        co_mappings = <uint64_t>(collected_count - 1)
        co_mappings_f = <double>co_mappings

        for i in range(collected_count):
            ref_idx = temp_refs[i]
            if ref_idx < array_size:
                co_mapping_sums[ref_idx] += co_mappings_f
                co_mapping_counts[ref_idx] += 1
                co_mapping_squares[ref_idx] += co_mappings_f * co_mappings_f

                if co_mappings > max_co_mappings[ref_idx]:
                    max_co_mappings[ref_idx] = co_mappings

        free(temp_refs)


cdef void accumulate_co_mappings_from_readrefs(uint32_t start_read, uint32_t end_read,
                                              MemoryPool* pool,
                                              ReadRefsIndex* rri,
                                              double* co_mapping_sums,
                                              uint32_t* co_mapping_counts,
                                              uint64_t* max_co_mappings,
                                              double* co_mapping_squares,
                                              uint32_t array_size) noexcept nogil:
    # Accumulate co-mapping statistics using a ReadRefsIndex (deduplicated per-read ref lists).
    # This is lighter-weight than scanning raw alignments and avoids per-read dedup logic here.
    cdef uint32_t read_idx, i, ref_idx
    cdef uint32_t ref_count
    cdef double co_mappings_f
    cdef uint64_t co_mappings
    cdef uint32_t* refs_ptr

    if not rri:
        return

    for read_idx in range(start_read, end_read):
        if read_idx >= rri.read_count:
            continue

        ref_count = rri.counts[read_idx]
        if ref_count <= 1:
            continue

        refs_ptr = rri.read_ptrs[read_idx]
        co_mappings = <uint64_t>(ref_count - 1)
        co_mappings_f = <double>co_mappings

        for i in range(ref_count):
            ref_idx = refs_ptr[i]
            if ref_idx < array_size:
                co_mapping_sums[ref_idx] += co_mappings_f
                co_mapping_counts[ref_idx] += 1
                co_mapping_squares[ref_idx] += co_mappings_f * co_mappings_f
                if co_mappings > max_co_mappings[ref_idx]:
                    max_co_mappings[ref_idx] = co_mappings

# -------------------------------------------------------------------------------
# DATASET HELPERS
# Functions for computing dataset-level statistics
# -------------------------------------------------------------------------------
cdef double compute_dataset_median_connections(uint32_t* exact_connection_counts,
                                              uint32_t* ref_read_counts,
                                              uint32_t array_size,
                                              uint32_t min_read_threshold) noexcept nogil:
    # Compute the median number of connections across references that meet the read threshold.
    # Args:
    #     exact_connection_counts: per-ref neighbor counts
    #     ref_read_counts: per-ref read counts
    #     array_size: number of references
    #     min_read_threshold: minimum reads for a reference to be considered
    # Returns:
    #     median connections as a double (0.0 if no valid refs)

    cdef uint32_t valid_count = 0
    cdef uint32_t ref_idx
    for ref_idx in range(array_size):
        if ref_read_counts[ref_idx] >= min_read_threshold and exact_connection_counts[ref_idx] > 0:
            valid_count += 1

    if valid_count == 0:

        return 0.0

    cdef double* valid_connections = <double*>malloc(valid_count * sizeof(double))
    if not valid_connections:

        return 0.0

    cdef uint32_t i = 0
    for ref_idx in range(array_size):
        if ref_read_counts[ref_idx] >= min_read_threshold and exact_connection_counts[ref_idx] > 0:
            valid_connections[i] = <double>exact_connection_counts[ref_idx]
            i += 1

    qsort(valid_connections, valid_count, sizeof(double), _double_compare)

    cdef double median
    if valid_count % 2 == 0:

        median = (valid_connections[valid_count / 2 - 1] + valid_connections[valid_count / 2]) / 2.0
    else:

        median = valid_connections[valid_count / 2]

    free(valid_connections)
    return median

# -------------------------------------------------------------------------------
# SORTING HELPERS
# qsort comparison functions used across this module
# -------------------------------------------------------------------------------

cdef int _double_compare(const void* a, const void* b) noexcept nogil:
    # qsort comparison for doubles (ascending)
    cdef double va = (<double*>a)[0]
    cdef double vb = (<double*>b)[0]
    if va < vb:
        return -1
    elif va > vb:
        return 1
    else:
        return 0

# Add this comparison function if it doesn't exist
cdef int _uint32_compare(const void* a, const void* b) noexcept nogil:
    # qsort comparison for uint32_t values (ascending)

    cdef uint32_t va = (<uint32_t*>a)[0]
    cdef uint32_t vb = (<uint32_t*>b)[0]
    if va < vb:
        return -1
    elif va > vb:
        return 1
    else:
        return 0


cdef int write_graph_tsv(MemoryPool* pool, sam_hdr_t* bam_header,
                        ReferenceMapping* mapping, ReferencePattern* pattern_data,
                        ReferenceStats* ref_stats,
                        uint32_t* total_reads, uint32_t* multimap_reads, uint64_t* alignments_per_ref,
                        uint32_t* exact_connection_counts, double* co_mapping_averages,
                        uint64_t* max_co_mappings, uint32_t* co_mapping_counts,
                        double* neighbor_multimap_avg, double* neighbor_connections_avg,
                        uint32_t* neighbor_counts, uint32_t array_size,
                        double dataset_median_connections,
                        int32_t min_read_count,
                        bint include_clustering,
                        int outlier_method,
                        const char* tsv_path,
                        TaxonomyDB* taxonomy_db) noexcept nogil:
    # Forwarder kept for ABI compatibility: call the implementation in the
    # separate processor_graph_tsv extension which contains the full logic.
    return write_graph_tsv_c(pool, bam_header, mapping, pattern_data, ref_stats,
                             total_reads, multimap_reads, alignments_per_ref,
                             exact_connection_counts, co_mapping_averages,
                             max_co_mappings, co_mapping_counts,
                             neighbor_multimap_avg, neighbor_connections_avg,
                             neighbor_counts, array_size, dataset_median_connections,
                             min_read_count, include_clustering, outlier_method, tsv_path,
                             taxonomy_db)
