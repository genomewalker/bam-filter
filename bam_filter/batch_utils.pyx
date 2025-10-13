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

# batch_utils.pyx
# Shared batching and sorting utilities for BAM processing
from libc.stdint cimport int64_t
from libc.stdlib cimport malloc, free
from bam_filter import logging as bf_logging

LOG_TAG = "BATCH"


# Shared batching and sorting utilities for BAM processing
from libc.stdint cimport int64_t
from libc.stdlib cimport malloc, free
from bam_filter.processor_types cimport ProcessingError, CompactAlignment, ProcessingBatch
cdef int compare_pairs_desc(const void* a, const void* b) noexcept nogil:
    """Compare function for descending order sorting by read count."""
    cdef int64_t* pair_a = <int64_t*>a
    cdef int64_t* pair_b = <int64_t*>b
    if pair_a[0] > pair_b[0]:
        return -1
    elif pair_a[0] < pair_b[0]:
        return 1
    return 0

cdef void qsort_tid_pairs(int64_t* tids, int64_t* counts, int64_t n) noexcept nogil:
    """Quick sort for (tid, count) pairs by TID for memory locality"""
    if n <= 1:
        return
    cdef int64_t pivot = tids[n // 2]
    cdef int64_t i = 0, j = n - 1
    cdef int64_t temp_tid, temp_count
    while i <= j:
        while tids[i] < pivot: i += 1
        while tids[j] > pivot: j -= 1
        if i <= j:
            temp_tid = tids[i]
            temp_count = counts[i]
            tids[i] = tids[j]
            counts[i] = counts[j]
            tids[j] = temp_tid
            counts[j] = temp_count
            i += 1
            j -= 1
    if j > 0:
        qsort_tid_pairs(tids, counts, j + 1)
    if i < n:
        qsort_tid_pairs(&tids[i], &counts[i], n - i)

cdef int create_simple_reference_batches(
    int64_t* tids_to_process,
    int64_t* tid_align_counts,
    int64_t n_tids_to_process,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int target_refs_per_batch,
    int max_refs_per_batch,
    int64_t max_batches,
    bint verbose
) noexcept nogil:
    """
    Partition references into batches while preserving tid locality.
    """
    if n_tids_to_process == 0:
        return 0
    cdef int64_t* sorted_tids = <int64_t*>malloc(n_tids_to_process * sizeof(int64_t))
    cdef int64_t* sorted_counts = <int64_t*>malloc(n_tids_to_process * sizeof(int64_t))
    if sorted_tids == NULL or sorted_counts == NULL:
        if sorted_tids != NULL: free(sorted_tids)
        if sorted_counts != NULL: free(sorted_counts)
        return 0
    cdef int64_t i
    for i in range(n_tids_to_process):
        sorted_tids[i] = tids_to_process[i]
        sorted_counts[i] = tid_align_counts[tids_to_process[i]]
    
    # Sort TIDs so batches iterate in BAM order
    qsort_tid_pairs(sorted_tids, sorted_counts, n_tids_to_process)
    
    cdef int64_t total_alignments = 0
    for i in range(n_tids_to_process):
        total_alignments += sorted_counts[i]
    cdef int64_t estimated_batches = n_tids_to_process // target_refs_per_batch
    if estimated_batches > max_batches:
        estimated_batches = max_batches
        target_refs_per_batch = n_tids_to_process // estimated_batches
    if n_tids_to_process > 10000000:
        estimated_batches = min(max_batches, max(estimated_batches, n_tids_to_process // 5000))
        target_refs_per_batch = n_tids_to_process // estimated_batches
    cdef int64_t target_alignments_per_batch = total_alignments // estimated_batches if estimated_batches > 0 else total_alignments
    with gil:
        if verbose and bf_logging.should_log(bf_logging.LogLevel.TRACE):
            bf_logging.log(LOG_TAG, "create_simple_reference_batches: consecutive TID batching enabled", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  references={n_tids_to_process:,}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  alignments={total_alignments:,}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  tid_range={sorted_tids[0]}-{sorted_tids[n_tids_to_process-1]}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  estimated_batches={estimated_batches:,}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  target_refs_per_batch={target_refs_per_batch:,}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  target_alignments_per_batch={target_alignments_per_batch:,}", level=bf_logging.LogLevel.TRACE)

    for i in range(max_batches):
        batch_starts[i] = n_tids_to_process
        batch_ends[i] = n_tids_to_process
    cdef int64_t current_batch = 0
    cdef int64_t current_start = 0
    cdef int64_t current_batch_alignments = 0
    for i in range(n_tids_to_process):
        current_batch_alignments += sorted_counts[i]
        # If we've reached or exceeded the target alignments, close the batch
        if current_batch_alignments >= target_alignments_per_batch and current_batch < max_batches - 1:
            batch_starts[current_batch] = current_start
            batch_ends[current_batch] = i + 1
            current_batch += 1
            current_start = i + 1
            current_batch_alignments = 0
    # Final batch for any remaining TIDs
    if current_start < n_tids_to_process and current_batch < max_batches:
        batch_starts[current_batch] = current_start
        batch_ends[current_batch] = n_tids_to_process
        current_batch += 1
    cdef int64_t final_batches = current_batch
    
    # Copy sorted TIDs back to maintain contiguity
    for i in range(n_tids_to_process):
        tids_to_process[i] = sorted_tids[i]
        tid_align_counts[sorted_tids[i]] = sorted_counts[i]
    
    # Log batch details
    cdef int64_t batch_idx, batch_alns, batch_refs
    cdef int64_t total_check = 0, min_alns = 9223372036854775807, max_alns = 0
    cdef int64_t batch_min_tid, batch_max_tid
    cdef double parallelization_factor

    # Compute aggregates across all batches first (fix: avoid using only the displayed subset)
    for batch_idx in range(final_batches):
        batch_alns = 0
        for i in range(batch_starts[batch_idx], batch_ends[batch_idx]):
            batch_alns += sorted_counts[i]
        total_check += batch_alns
        if batch_alns < min_alns: min_alns = batch_alns
        if batch_alns > max_alns: max_alns = batch_alns

    avg_alns = <double>total_check / final_batches if final_batches > 0 else 0.0

    with gil:
        if bf_logging.should_log(bf_logging.LogLevel.INFO):
            bf_logging.log(
                LOG_TAG,
                f"Created {final_batches:,} batches (~{avg_alns:,.0f} alignments per batch) for {n_tids_to_process:,} references",
                level=bf_logging.LogLevel.INFO,
            )
        if verbose and bf_logging.should_log(bf_logging.LogLevel.TRACE):
            bf_logging.log(LOG_TAG, "create_simple_reference_batches: sample batch layout", level=bf_logging.LogLevel.TRACE)
            # Log a representative subset of batches
            for batch_idx in range(min(10, final_batches)):
                batch_alns = 0
                batch_refs = batch_ends[batch_idx] - batch_starts[batch_idx]
                for i in range(batch_starts[batch_idx], batch_ends[batch_idx]):
                    batch_alns += sorted_counts[i]
                batch_min_tid = sorted_tids[batch_starts[batch_idx]]
                batch_max_tid = sorted_tids[batch_ends[batch_idx]-1]
                bf_logging.log(LOG_TAG, f"  batch_{batch_idx}: refs={batch_refs:,}, alignments={batch_alns:,}, tids={batch_min_tid}-{batch_max_tid}", level=bf_logging.LogLevel.TRACE)
            if final_batches > 20:
                bf_logging.log(LOG_TAG, f"  showing first 10 of {final_batches:,} total batches", level=bf_logging.LogLevel.TRACE)
                for batch_idx in range(max(10, final_batches-5), final_batches):
                    batch_alns = 0
                    batch_refs = batch_ends[batch_idx] - batch_starts[batch_idx]
                    for i in range(batch_starts[batch_idx], batch_ends[batch_idx]):
                        batch_alns += sorted_counts[i]
                    batch_min_tid = sorted_tids[batch_starts[batch_idx]]
                    batch_max_tid = sorted_tids[batch_ends[batch_idx]-1]
                    bf_logging.log(LOG_TAG, f"  batch_{batch_idx}: refs={batch_refs:,}, alignments={batch_alns:,}, tids={batch_min_tid}-{batch_max_tid}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  min_alignments_per_batch={min_alns:,}", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, f"  max_alignments_per_batch={max_alns:,}", level=bf_logging.LogLevel.TRACE)
            parallelization_factor = <double>final_batches / 12.0
            bf_logging.log(LOG_TAG, f"  parallelization_factor_estimate={parallelization_factor:.2f}x threads", level=bf_logging.LogLevel.TRACE)
            bf_logging.log(LOG_TAG, "  memory_locality=consecutive_tid_access", level=bf_logging.LogLevel.TRACE)

        free(sorted_tids)
        free(sorted_counts)
    return final_batches

cdef int create_smart_batches_for_large_datasets(
    int64_t* tids_to_process,
    int64_t* tid_align_counts,
    int64_t n_tids_to_process,
    int64_t max_batches,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int num_threads,
    bint verbose
) noexcept nogil:
    """
    Determine target batch sizes for large reference sets.
    """
    cdef int64_t total_reads = 0
    cdef int64_t i, tid
    for i in range(n_tids_to_process):
        tid = tids_to_process[i]
        total_reads += tid_align_counts[tid]
    cdef double avg_reads_per_ref = <double>total_reads / n_tids_to_process if n_tids_to_process > 0 else 25.0
    cdef int target_refs_per_batch, max_refs_per_batch
    cdef int64_t target_batches
    cdef int64_t expected_alns_per_batch
    
    if n_tids_to_process > 20000000:
        target_batches = min(max_batches, max(10000, n_tids_to_process // 2000))
        target_refs_per_batch = n_tids_to_process // target_batches
        max_refs_per_batch = target_refs_per_batch * 2
        if verbose:
            with gil:
                bf_logging.log(LOG_TAG, "create_smart_batches: ultra-large dataset heuristics applied", level=bf_logging.LogLevel.DEBUG)
                bf_logging.log(LOG_TAG, f"  target_batches={target_batches:,}", level=bf_logging.LogLevel.DEBUG)
                bf_logging.log(LOG_TAG, f"  target_refs_per_batch={target_refs_per_batch:,}", level=bf_logging.LogLevel.DEBUG)
    elif n_tids_to_process > 10000000:
        target_batches = min(max_batches, max(5000, n_tids_to_process // 3000))
        target_refs_per_batch = n_tids_to_process // target_batches
        max_refs_per_batch = target_refs_per_batch * 2
    elif n_tids_to_process > 1000000:
        target_batches = min(max_batches, max(1000, n_tids_to_process // 5000))
        target_refs_per_batch = n_tids_to_process // target_batches
        max_refs_per_batch = target_refs_per_batch * 2
    else:
        target_batches = min(max_batches, max(num_threads * 8, n_tids_to_process // 10000))
        target_refs_per_batch = max(2000, n_tids_to_process // target_batches)
        max_refs_per_batch = target_refs_per_batch * 3
        
   #  target_refs_per_batch = max(100, min(10000, target_refs_per_batch))
    max_refs_per_batch = max(target_refs_per_batch, min(20000, max_refs_per_batch))
    target_batches = n_tids_to_process // target_refs_per_batch
    target_batches = min(max_batches, max(target_batches, num_threads * 8))

    expected_alns_per_batch = total_reads // target_batches
    with gil:
        if verbose:
            bf_logging.log(LOG_TAG, "create_smart_batches_for_large_datasets: planning summary", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  references={n_tids_to_process:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  alignments={total_reads:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  avg_alignments_per_ref={avg_reads_per_ref:.1f}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  threads={num_threads}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  target_batches={target_batches:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  target_refs_per_batch={target_refs_per_batch:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  max_refs_per_batch={max_refs_per_batch:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, f"  expected_alignments_per_batch~{expected_alns_per_batch:,}", level=bf_logging.LogLevel.DEBUG)
            bf_logging.log(LOG_TAG, "  memory_locality=consecutive TID access", level=bf_logging.LogLevel.DEBUG)
    cdef int actual_batches = create_simple_reference_batches(
        tids_to_process, tid_align_counts, n_tids_to_process,
        batch_starts, batch_ends,
        target_refs_per_batch, max_refs_per_batch, 
        max_batches, verbose
    )
    return actual_batches

cdef int create_balanced_batches_greedy(
    int64_t* reference_ids,
    int64_t* reference_alignment_counts,
    int64_t num_references,
    int64_t max_batches,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int num_threads,
    bint verbose
) noexcept nogil:
    """
    C-only: Create balanced batches using greedy algorithm for BAM references.
    Returns the number of batches created, with batch_starts and batch_ends filled in.
    """
    cdef int estimated_threads = max(1, min(32, <int>(max_batches / 4)))
    cdef int actual_batches = create_smart_batches_for_large_datasets(
        reference_ids, reference_alignment_counts, num_references,
        max_batches, batch_starts, batch_ends,
        num_threads, verbose
    )
    return actual_batches
