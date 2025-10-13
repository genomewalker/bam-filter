# cython: language_level=3
from libc.stdint cimport int64_t
# batch_utils.pxd
# Cython declarations for shared batching utilities

cdef int compare_pairs_desc(const void* a, const void* b) noexcept nogil
cdef void qsort_tid_pairs(int64_t* tids, int64_t* counts, int64_t n) noexcept nogil
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
) nogil
cdef int create_smart_batches_for_large_datasets(
    int64_t* tids_to_process,
    int64_t* tid_align_counts,
    int64_t n_tids_to_process,
    int64_t max_batches,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int num_threads,
    bint verbose
) nogil
cdef int create_balanced_batches_greedy(
    int64_t* reference_ids,
    int64_t* reference_alignment_counts,
    int64_t num_references,
    int64_t max_batches,
    int64_t* batch_starts,
    int64_t* batch_ends,
    int num_threads,
    bint verbose
) nogil
