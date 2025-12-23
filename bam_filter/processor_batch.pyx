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

"""Parallel batch processing for BAM alignments.

This module handles parallel reading and processing of BAM alignments using
balanced reference batches with thread-local storage for maximum throughput.
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t, uint8_t, uint16_t
from libc.stdlib cimport malloc, free, realloc
from libc.string cimport memcpy, memset
from libc.stddef cimport size_t
from libc.math cimport log, fmax

from bam_filter.processor_batch cimport *
from bam_filter.processor cimport (
    Alignment, MemoryPool, AlignmentScoringConfig,
    AlignmentCore, HierarchicalData, DamageCounts, BAMWriterAux,
)
from bam_filter.processor cimport min_int64, max_int64, INVALID_SEQUENTIAL_ID
from bam_filter.processor_sort cimport radix_sort_alignments_by_read_id, radix_sort_uint64, radix_sort_compact_by_position
from bam_filter.processor_hash cimport ThreadLocalHashMap, extract_read_hash_identifier
from bam_filter.processor_md_quality cimport (
    calculate_md_quality_score,
    calculate_md_quality_score_with_stats,
    alignment_passes_quality_filters,
    ANIStats,
)
from bam_filter.processor_types cimport ProcessingError, PROCESSING_SUCCESS, PROCESSING_ERROR_MEMORY_ALLOCATION

# Import PMD context, thread accumulator, and corrected ANI computation
from bam_filter.processor_pmd cimport (
    PMDGlobalContext,
    PMDStatsAccumulator,
    PMDCurve,
    ANISnapshot,
    compute_corrected_ani,
)
from cython.parallel cimport prange, threadid
from .processor_types cimport (
    BGZF,
    bam1_core_t,
    bam1_t,
    sam_hdr_t,
    hts_idx_t,
    hts_itr_t,
    htsFile,
    samFile,
    sam_read1,
    sam_write1,
    sam_hdr_read,
    sam_hdr_destroy,
    sam_hdr_str,
    sam_hdr_tid2name,
    sam_hdr_tid2len,
    sam_hdr_parse,
    sam_hdr_nref,
    bam_endpos,
    bam_init1,
    bam_destroy1,
    bam_dup1,
    bam_get_qname,
    bam_aux_get,
    bam_aux2i,
    bam_aux_del,
    bam_aux_append,
    sam_index_load,
    hts_idx_destroy,
    hts_idx_get_stat,
    sam_itr_queryi,
    sam_itr_next,
    hts_itr_destroy,
    sam_hdr_name2tid,
    hts_set_threads,
    bam_get_seq,
    bam_get_qual,
    bam_get_cigar,
    seq_nt16_str,
    bam_seqi_wrapper,
    hts_open,
    hts_close,
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil


cdef void* MAP_FAILED_PTR = <void*>-1

cdef ProcessingBatch* create_processing_batch(int64_t batch_id, int64_t ref_start,
                                               int64_t ref_end, int64_t expected_count) except NULL nogil:
    """Create and initialize a processing batch.

    Allocates a batch structure with intelligent initial capacity. Uses mmap
    for large batches to improve memory locality and reduce allocation overhead.

    Parameters
    ----------
    batch_id : int64_t
        Unique batch identifier
    ref_start : int64_t
        Starting reference index (inclusive)
    ref_end : int64_t
        Ending reference index (exclusive)
    expected_count : int64_t
        Expected number of alignments in batch

    Returns
    -------
    ProcessingBatch*
        Initialized batch structure, or NULL on allocation failure
    """
    cdef ProcessingBatch* batch = <ProcessingBatch*>malloc(sizeof(ProcessingBatch))
    if not batch:
        return NULL

    batch.batch_identifier = batch_id
    batch.reference_start_index = ref_start
    batch.reference_end_index = ref_end
    batch.expected_alignment_count = expected_count
    batch.actual_alignment_count = 0
    batch.error_status = PROCESSING_SUCCESS
    batch.processed_by_thread_id = -1
    batch.uses_mmap = False
    batch.mmap_size = 0

    cdef int64_t smart_initial = max_int64(5000, expected_count // 200)
    batch.batch_capacity = min_int64(50000, smart_initial)

    cdef size_t batch_memory_size = batch.batch_capacity * sizeof(BatchAlignment)
    cdef size_t mmap_threshold = 32 * 1024 * 1024

    if batch_memory_size > mmap_threshold:
        batch.batch_alignments = <BatchAlignment*>mmap(NULL, batch_memory_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
        if batch.batch_alignments == MAP_FAILED_PTR:
            batch.batch_alignments = <BatchAlignment*>malloc(batch_memory_size)
            if not batch.batch_alignments:
                free(batch)
                return NULL
            batch.uses_mmap = False
            batch.mmap_size = 0
        else:
            batch.uses_mmap = True
            batch.mmap_size = batch_memory_size
            madvise(batch.batch_alignments, batch_memory_size, MADV_SEQUENTIAL)
    else:
        batch.batch_alignments = <BatchAlignment*>malloc(batch_memory_size)
        if not batch.batch_alignments:
            free(batch)
            return NULL
        batch.uses_mmap = False
        batch.mmap_size = 0

    return batch


cdef int grow_batch_capacity(ProcessingBatch* batch) except -1 nogil:
    """Expand batch capacity by 25%.

    Reallocates alignment array when batch fills up. Handles both regular
    malloc and mmap allocations appropriately.

    Parameters
    ----------
    batch : ProcessingBatch*
        Batch to expand

    Returns
    -------
    int
        0 on success, -1 on allocation failure
    """
    cdef int64_t new_capacity = batch.batch_capacity + (batch.batch_capacity >> 2)
    cdef size_t old_size = batch.batch_capacity * sizeof(BatchAlignment)
    cdef size_t new_size = new_capacity * sizeof(BatchAlignment)
    cdef BatchAlignment* new_alignments = NULL
    cdef size_t mmap_threshold = 64 * 1024 * 1024

    if batch.uses_mmap:
        new_alignments = <BatchAlignment*>mmap(NULL, new_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
        if new_alignments == MAP_FAILED_PTR:
            bf_nogil_logf_notime(
                b"BATCH",
                "Smart growth %ld: mmap growth failed, trying malloc",
                batch.batch_identifier,
            )
            new_alignments = <BatchAlignment*>malloc(new_size)
            if not new_alignments:
                return -1
            memcpy(new_alignments, batch.batch_alignments, old_size)
            munmap(batch.batch_alignments, batch.mmap_size)
            batch.batch_alignments = new_alignments
            batch.uses_mmap = False
            batch.mmap_size = 0
        else:
            memcpy(new_alignments, batch.batch_alignments, old_size)
            munmap(batch.batch_alignments, batch.mmap_size)
            batch.batch_alignments = new_alignments
            batch.mmap_size = new_size
    else:
        if new_size <= mmap_threshold:
            new_alignments = <BatchAlignment*>realloc(batch.batch_alignments, new_size)
            if not new_alignments:
                return -1
            batch.batch_alignments = new_alignments
        else:
            new_alignments = <BatchAlignment*>mmap(NULL, new_size, PROT_READ | PROT_WRITE, MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
            if new_alignments == MAP_FAILED_PTR:
                bf_nogil_logf_notime(
                    b"BATCH",
                    "Smart growth %ld: mmap switch failed, using realloc",
                    batch.batch_identifier,
                )
                new_alignments = <BatchAlignment*>realloc(batch.batch_alignments, new_size)
                if not new_alignments:
                    return -1
                batch.batch_alignments = new_alignments
            else:
                memcpy(new_alignments, batch.batch_alignments, old_size)
                free(batch.batch_alignments)
                batch.batch_alignments = new_alignments
                batch.uses_mmap = True
                batch.mmap_size = new_size

    batch.batch_capacity = new_capacity
    return 0


cdef void destroy_processing_batch(ProcessingBatch* batch) noexcept nogil:
    """Free all memory associated with a processing batch.

    Parameters
    ----------
    batch : ProcessingBatch*
        Batch to destroy (safe to pass NULL)
    """
    if not batch:
        return
    if batch.batch_alignments:
        if batch.uses_mmap:
            munmap(batch.batch_alignments, batch.mmap_size)
        else:
            free(batch.batch_alignments)
        batch.batch_alignments = NULL
    free(batch)


cdef int process_batch_alignments(samFile* bam_file, sam_hdr_t* header,
                                           hts_idx_t* index, int64_t* reference_ids,
                                           ProcessingBatch* batch,
                                           AlignmentScoringConfig* scoring_config,
                                           ThreadLocalHashMap* thread_map,
                                           int32_t thread_id) except -1 nogil:
    """Process alignments for a batch of references.

    Reads BAM file for specified reference range, filters by quality, assigns
    read IDs using thread-local hash map, calculates alignment scores and
    optional PMD scores.

    Parameters
    ----------
    bam_file : samFile*
        Open BAM file handle
    header : sam_hdr_t*
        BAM header
    index : hts_idx_t*
        BAM index for random access
    reference_ids : int64_t*
        Array mapping batch indices to BAM reference IDs
    batch : ProcessingBatch*
        Batch structure to populate
    scoring_config : AlignmentScoringConfig*
        Scoring and filtering configuration
    thread_map : ThreadLocalHashMap*
        Thread-local hash map for read ID assignment
    thread_id : int32_t
        Thread identifier

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef hts_itr_t* iterator = NULL
    cdef bam1_t* bam_record = NULL
    cdef int64_t reference_id
    cdef int32_t result_code
    cdef int64_t alignment_position
    cdef double raw_score
    cdef uint64_t read_hash
    cdef uint32_t local_sequential_id
    cdef khint_t k
    cdef int ret
    cdef BatchAlignment temp_alignment
    cdef int64_t ref_idx
    cdef float temp_pmd_score = 0.0
    cdef float* pmd_ptr = NULL

    if scoring_config.calculate_pmd:
        pmd_ptr = &temp_pmd_score

    batch.processed_by_thread_id = thread_id

    bam_record = bam_init1()
    if not bam_record:
        batch.error_status = PROCESSING_ERROR_MEMORY_ALLOCATION
        return -1

    for ref_idx in range(batch.reference_start_index, batch.reference_end_index):
        reference_id = reference_ids[ref_idx]
        iterator = sam_itr_queryi(index, reference_id, 0, 0x7fffffff)
        if not iterator:
            continue

        alignment_position = 0
        while True:
            result_code = sam_itr_next(bam_file, iterator, bam_record)
            if result_code < 0:
                break

            if not alignment_passes_quality_filters(bam_record, scoring_config):
                alignment_position += 1
                continue

            read_hash = <uint64_t>extract_read_hash_identifier(bam_record)
            k = kh_get_seqid_map(thread_map.hash_to_id_map, read_hash)
            if k == kh_end_seqid_map(thread_map.hash_to_id_map):
                local_sequential_id = thread_map.next_local_id
                thread_map.next_local_id += 1
                k = kh_put_seqid_map(thread_map.hash_to_id_map, read_hash, &ret)
                if ret != -1:
                    kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0] = local_sequential_id
            else:
                local_sequential_id = kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0]

            temp_alignment.read_index = local_sequential_id
            temp_alignment.reference_index = <uint32_t>reference_id
            temp_alignment.alignment_position = <uint32_t>alignment_position
            raw_score = calculate_md_quality_score(bam_record, header, scoring_config, pmd_ptr)
            temp_alignment.alignment_score = raw_score
            temp_alignment.pmd_score = temp_pmd_score if scoring_config.calculate_pmd else 0.0

            if batch.actual_alignment_count >= batch.batch_capacity:
                if grow_batch_capacity(batch) != 0:
                    batch.error_status = PROCESSING_ERROR_MEMORY_ALLOCATION
                    hts_itr_destroy(iterator)
                    bam_destroy1(bam_record)
                    return -1

            batch.batch_alignments[batch.actual_alignment_count] = temp_alignment
            batch.actual_alignment_count += 1
            alignment_position += 1

        hts_itr_destroy(iterator)
        iterator = NULL

    bam_destroy1(bam_record)
    return 0


cdef int process_batch_alignments_with_pmd(samFile* bam_file, sam_hdr_t* header,
                                            hts_idx_t* index, int64_t* reference_ids,
                                            ProcessingBatch* batch,
                                            AlignmentScoringConfig* scoring_config,
                                            ThreadLocalHashMap* thread_map,
                                            int32_t thread_id,
                                            PMDGlobalContext* pmd_context) except -1 nogil:
    """Process alignments with PMD statistics collection.

    Extended version of process_batch_alignments that also collects PMD statistics
    into thread-local accumulators and stores ANI snapshot values for each alignment.

    Parameters
    ----------
    bam_file : samFile*
        Open BAM file handle
    header : sam_hdr_t*
        BAM header
    index : hts_idx_t*
        BAM index for random access
    reference_ids : int64_t*
        Array mapping batch indices to BAM reference IDs
    batch : ProcessingBatch*
        Batch structure to populate
    scoring_config : AlignmentScoringConfig*
        Scoring and filtering configuration
    thread_map : ThreadLocalHashMap*
        Thread-local hash map for read ID assignment
    thread_id : int32_t
        Thread identifier
    pmd_context : PMDGlobalContext*
        Global PMD context with thread-local accumulators (NULL to skip PMD stats)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef hts_itr_t* iterator = NULL
    cdef bam1_t* bam_record = NULL
    cdef int64_t reference_id
    cdef int32_t result_code
    cdef int64_t alignment_position
    cdef double raw_score
    cdef uint64_t read_hash
    cdef uint32_t local_sequential_id
    cdef khint_t k
    cdef int ret
    cdef BatchAlignment temp_alignment
    cdef int64_t ref_idx
    cdef float temp_pmd_score = 0.0
    cdef float* pmd_ptr = NULL
    cdef ANIStats ani_stats
    cdef ANIStats* ani_ptr = NULL
    cdef PMDStatsAccumulator* pmd_acc = NULL

    # Set up PMD output pointer if enabled
    if scoring_config.calculate_pmd:
        pmd_ptr = &temp_pmd_score

    # Get thread-local PMD accumulator if context provided and stats collection enabled
    if pmd_context != NULL and pmd_context.collect_stats:
        if thread_id >= 0 and thread_id < pmd_context.num_threads:
            pmd_acc = &pmd_context.thread_contexts[thread_id].stats
        ani_ptr = &ani_stats

    batch.processed_by_thread_id = thread_id

    bam_record = bam_init1()
    if not bam_record:
        batch.error_status = PROCESSING_ERROR_MEMORY_ALLOCATION
        return -1

    for ref_idx in range(batch.reference_start_index, batch.reference_end_index):
        reference_id = reference_ids[ref_idx]
        iterator = sam_itr_queryi(index, reference_id, 0, 0x7fffffff)
        if not iterator:
            continue

        alignment_position = 0
        while True:
            result_code = sam_itr_next(bam_file, iterator, bam_record)
            if result_code < 0:
                break

            if not alignment_passes_quality_filters(bam_record, scoring_config):
                alignment_position += 1
                continue

            read_hash = <uint64_t>extract_read_hash_identifier(bam_record)
            k = kh_get_seqid_map(thread_map.hash_to_id_map, read_hash)
            if k == kh_end_seqid_map(thread_map.hash_to_id_map):
                local_sequential_id = thread_map.next_local_id
                thread_map.next_local_id += 1
                k = kh_put_seqid_map(thread_map.hash_to_id_map, read_hash, &ret)
                if ret != -1:
                    kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0] = local_sequential_id
            else:
                local_sequential_id = kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0]

            temp_alignment.read_index = local_sequential_id
            temp_alignment.reference_index = <uint32_t>reference_id
            temp_alignment.alignment_position = <uint32_t>alignment_position

            # Use extended function that collects PMD stats and ANI values
            if pmd_acc != NULL or ani_ptr != NULL:
                # Reset ANI stats for this alignment
                if ani_ptr != NULL:
                    ani_stats.aligned_length = 0
                    ani_stats.match_count = 0
                    ani_stats.ct_5p_count = 0
                    ani_stats.ga_3p_count = 0
                    ani_stats.other_mm_count = 0
                    ani_stats.c_at_5p_count = 0
                    ani_stats.g_at_3p_count = 0
                    ani_stats.damage_llr = 0.0

                raw_score = calculate_md_quality_score_with_stats(
                    bam_record, header, scoring_config, pmd_ptr, ani_ptr, pmd_acc
                )

                # Store ANI values in batch alignment
                if ani_ptr != NULL:
                    temp_alignment.aligned_length = ani_stats.aligned_length
                    temp_alignment.match_count = ani_stats.match_count
                    temp_alignment.ct_5p_count = ani_stats.ct_5p_count
                    temp_alignment.ga_3p_count = ani_stats.ga_3p_count
                    temp_alignment.c_at_5p_count = ani_stats.c_at_5p_count
                    temp_alignment.g_at_3p_count = ani_stats.g_at_3p_count
                    temp_alignment.damage_llr = ani_stats.damage_llr
                else:
                    temp_alignment.aligned_length = 0
                    temp_alignment.match_count = 0
                    temp_alignment.ct_5p_count = 0
                    temp_alignment.ga_3p_count = 0
                    temp_alignment.c_at_5p_count = 0
                    temp_alignment.g_at_3p_count = 0
                    temp_alignment.damage_llr = 0.0
            else:
                # Fall back to original function when PMD stats not needed
                raw_score = calculate_md_quality_score(bam_record, header, scoring_config, pmd_ptr)
                temp_alignment.aligned_length = 0
                temp_alignment.match_count = 0
                temp_alignment.ct_5p_count = 0
                temp_alignment.ga_3p_count = 0
                temp_alignment.c_at_5p_count = 0
                temp_alignment.g_at_3p_count = 0
                temp_alignment.damage_llr = 0.0

            temp_alignment.alignment_score = raw_score
            temp_alignment.pmd_score = temp_pmd_score if scoring_config.calculate_pmd else 0.0

            if batch.actual_alignment_count >= batch.batch_capacity:
                if grow_batch_capacity(batch) != 0:
                    batch.error_status = PROCESSING_ERROR_MEMORY_ALLOCATION
                    hts_itr_destroy(iterator)
                    bam_destroy1(bam_record)
                    return -1

            batch.batch_alignments[batch.actual_alignment_count] = temp_alignment
            batch.actual_alignment_count += 1
            alignment_position += 1

        hts_itr_destroy(iterator)
        iterator = NULL

    bam_destroy1(bam_record)
    return 0


cdef int assign_global_sequential_ids_fast(ProcessingBatch** batches, int64_t batch_count,
                                          ThreadLocalHashMap** thread_maps, int num_threads) except -1 nogil:
    """Assign global sequential read IDs across all thread-local hash maps.

    Merges thread-local read IDs into globally consistent IDs. Sorts all unique
    hash values to ensure deterministic ordering, then creates mapping tables
    for each thread and updates all batch alignments.

    Parameters
    ----------
    batches : ProcessingBatch**
        Array of processing batches
    batch_count : int64_t
        Number of batches
    thread_maps : ThreadLocalHashMap**
        Array of thread-local hash maps
    num_threads : int
        Number of threads

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef kh_seqid_map_t* global_map = kh_init_seqid_map()
    cdef khint_t k, global_k
    cdef uint64_t hash_value
    cdef int ret, thread_id
    cdef uint32_t local_id, global_id
    cdef int64_t hash_capacity = batch_count * 1000
    cdef uint64_t* all_hashes = <uint64_t*>malloc(hash_capacity * sizeof(uint64_t))
    if not all_hashes:
        kh_destroy_seqid_map(global_map)
        return -1

    cdef int64_t total_unique_hashes = 0
    for thread_id in range(num_threads):
        if not thread_maps[thread_id]:
            continue
        k = 0
        while k < kh_end_seqid_map(thread_maps[thread_id].hash_to_id_map):
            if kh_exist_seqid_map(thread_maps[thread_id].hash_to_id_map, k):
                hash_value = kh_key_seqid_map(thread_maps[thread_id].hash_to_id_map, k)
                global_k = kh_get_seqid_map(global_map, hash_value)
                if global_k == kh_end_seqid_map(global_map):
                    global_k = kh_put_seqid_map(global_map, hash_value, &ret)
                    if ret != -1:
                        if total_unique_hashes >= hash_capacity:
                            hash_capacity *= 2
                            __tmp_hashes = <uint64_t*>realloc(all_hashes, hash_capacity * sizeof(uint64_t))
                            if not __tmp_hashes:
                                kh_destroy_seqid_map(global_map)
                                free(all_hashes)
                                return -1
                            all_hashes = __tmp_hashes
                        all_hashes[total_unique_hashes] = hash_value
                        total_unique_hashes += 1
            k += 1

    if total_unique_hashes > 1:
        radix_sort_uint64(all_hashes, total_unique_hashes)

    kh_destroy_seqid_map(global_map)
    global_map = kh_init_seqid_map()
    for hash_idx in range(total_unique_hashes):
        hash_value = all_hashes[hash_idx]
        global_k = kh_put_seqid_map(global_map, hash_value, &ret)
        if ret != -1:
            kh_val_seqid_map_wrap(global_map, global_k)[0] = <uint32_t>hash_idx

    free(all_hashes)

    cdef uint32_t** thread_local_to_global = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    cdef uint32_t* thread_max_ids = <uint32_t*>malloc(num_threads * sizeof(uint32_t))
    if not thread_local_to_global or not thread_max_ids:
        kh_destroy_seqid_map(global_map)
        if thread_local_to_global: free(thread_local_to_global)
        if thread_max_ids: free(thread_max_ids)
        return -1

    for thread_id in range(num_threads):
        thread_max_ids[thread_id] = 0
        if thread_maps[thread_id]:
            k = 0
            while k < kh_end_seqid_map(thread_maps[thread_id].hash_to_id_map):
                if kh_exist_seqid_map(thread_maps[thread_id].hash_to_id_map, k):
                    local_id = kh_val_seqid_map_wrap(thread_maps[thread_id].hash_to_id_map, k)[0]
                    if local_id > thread_max_ids[thread_id]:
                        thread_max_ids[thread_id] = local_id
                k += 1

    for thread_id in range(num_threads):
        if thread_max_ids[thread_id] > 0:
            thread_local_to_global[thread_id] = <uint32_t*>malloc((thread_max_ids[thread_id] + 1) * sizeof(uint32_t))
            if not thread_local_to_global[thread_id]:
                for t in range(thread_id):
                    if thread_local_to_global[t]:
                        free(thread_local_to_global[t])
                free(thread_local_to_global)
                free(thread_max_ids)
                kh_destroy_seqid_map(global_map)
                return -1
            memset(thread_local_to_global[thread_id], INVALID_SEQUENTIAL_ID, (thread_max_ids[thread_id] + 1) * sizeof(uint32_t))
        else:
            thread_local_to_global[thread_id] = NULL

    for thread_id in range(num_threads):
        if not thread_maps[thread_id] or not thread_local_to_global[thread_id]:
            continue
        k = 0
        while k < kh_end_seqid_map(thread_maps[thread_id].hash_to_id_map):
            if kh_exist_seqid_map(thread_maps[thread_id].hash_to_id_map, k):
                hash_value = kh_key_seqid_map(thread_maps[thread_id].hash_to_id_map, k)
                local_id = kh_val_seqid_map_wrap(thread_maps[thread_id].hash_to_id_map, k)[0]
                global_k = kh_get_seqid_map(global_map, hash_value)
                if global_k != kh_end_seqid_map(global_map):
                    global_id = kh_val_seqid_map_wrap(global_map, global_k)[0]
                    if local_id <= thread_max_ids[thread_id]:
                        thread_local_to_global[thread_id][local_id] = global_id
            k += 1

    cdef int64_t batch_idx, alignment_idx
    cdef uint32_t thread_local_id
    cdef int64_t remapped_count = 0
    cdef int64_t error_count = 0
    cdef int32_t actual_thread_id

    for batch_idx in range(batch_count):
        if not batches[batch_idx]:
            continue
        actual_thread_id = batches[batch_idx].processed_by_thread_id
        if actual_thread_id < 0 or actual_thread_id >= num_threads or not thread_local_to_global[actual_thread_id]:
            error_count += 1
            continue
        for alignment_idx in range(batches[batch_idx].actual_alignment_count):
            thread_local_id = batches[batch_idx].batch_alignments[alignment_idx].read_index
            if thread_local_id <= thread_max_ids[actual_thread_id]:
                global_id = thread_local_to_global[actual_thread_id][thread_local_id]
                if global_id != INVALID_SEQUENTIAL_ID:
                    batches[batch_idx].batch_alignments[alignment_idx].read_index = global_id
                    remapped_count += 1
                else:
                    error_count += 1
            else:
                error_count += 1

    for thread_id in range(num_threads):
        if thread_local_to_global[thread_id]:
            free(thread_local_to_global[thread_id])
    free(thread_local_to_global)
    free(thread_max_ids)
    kh_destroy_seqid_map(global_map)

    return 0 if error_count == 0 else -1


cdef int64_t count_actual_alignments(ProcessingBatch** batches, int64_t batch_count) except -1 nogil:
    """Count total alignments across all batches.

    Parameters
    ----------
    batches : ProcessingBatch**
        Array of processing batches
    batch_count : int64_t
        Number of batches

    Returns
    -------
    int64_t
        Total alignment count
    """
    cdef int64_t total_count = 0
    cdef int64_t batch_idx
    for batch_idx in range(batch_count):
        if batches[batch_idx]:
            total_count += batches[batch_idx].actual_alignment_count
    bf_nogil_logf_notime(
        b"BATCH",
        "Processed %ld alignments across %ld parallel batches",
        total_count,
        batch_count,
    )
    return total_count


cdef int64_t count_unique_refs_from_batches(ProcessingBatch** batches, int64_t batch_count) except -1 nogil:
    """Count unique references with alignments across all batches.

    Parameters
    ----------
    batches : ProcessingBatch**
        Array of processing batches
    batch_count : int64_t
        Number of batches

    Returns
    -------
    int64_t
        Number of unique references
    """
    cdef kh_seqid_map_t* ref_map = kh_init_seqid_map()
    cdef int64_t batch_idx, alignment_idx
    cdef uint32_t ref_id
    cdef int ret
    cdef int64_t unique_refs

    if not ref_map:
        return -1

    for batch_idx in range(batch_count):
        if not batches[batch_idx]:
            continue
        for alignment_idx in range(batches[batch_idx].actual_alignment_count):
            ref_id = batches[batch_idx].batch_alignments[alignment_idx].reference_index
            kh_put_seqid_map(ref_map, ref_id, &ret)

    unique_refs = kh_size(ref_map)
    kh_destroy_seqid_map(ref_map)
    bf_nogil_logf_notime(
        b"BATCH",
        "Identified %ld unique references with alignments",
        unique_refs,
    )
    return unique_refs


cdef void parallel_streaming_stats_optimized(ProcessingBatch** batches, int64_t batch_count,
                                            double* out_min, double* out_max, double* out_mean,
                                            double* out_variance, int64_t* out_count, int num_threads) noexcept nogil:
    """Calculate alignment score statistics in parallel using Welford's algorithm.

    Computes per-thread statistics then combines using parallel variance formula
    for numerical stability.

    Parameters
    ----------
    batches : ProcessingBatch**
        Array of processing batches
    batch_count : int64_t
        Number of batches
    out_min : double*
        Output minimum score
    out_max : double*
        Output maximum score
    out_mean : double*
        Output mean score
    out_variance : double*
        Output variance
    out_count : int64_t*
        Output total alignment count
    num_threads : int
        Number of threads for parallel processing
    """
    cdef double* thread_mins = <double*>malloc(num_threads * sizeof(double))
    cdef double* thread_maxs = <double*>malloc(num_threads * sizeof(double))
    cdef double* thread_means = <double*>malloc(num_threads * sizeof(double))
    cdef double* thread_M2s = <double*>malloc(num_threads * sizeof(double))
    cdef int64_t* thread_counts = <int64_t*>malloc(num_threads * sizeof(int64_t))
    cdef int thread_id, batch_idx, alignment_idx
    cdef float score
    cdef double delta, delta2

    if not thread_mins or not thread_maxs or not thread_means or not thread_M2s or not thread_counts:
        if thread_mins: free(thread_mins)
        if thread_maxs: free(thread_maxs)
        if thread_means: free(thread_means)
        if thread_M2s: free(thread_M2s)
        if thread_counts: free(thread_counts)
        return

    for thread_id in range(num_threads):
        thread_mins[thread_id] = 1e30
        thread_maxs[thread_id] = -1e30
        thread_means[thread_id] = 0.0
        thread_M2s[thread_id] = 0.0
        thread_counts[thread_id] = 0

    for batch_idx in prange(batch_count, num_threads=num_threads, schedule='static', nogil=True):
        thread_id = threadid()
        if not batches[batch_idx]:
            continue

        for alignment_idx in range(batches[batch_idx].actual_alignment_count):
            score = batches[batch_idx].batch_alignments[alignment_idx].alignment_score
            if score < thread_mins[thread_id]:
                thread_mins[thread_id] = score
            if score > thread_maxs[thread_id]:
                thread_maxs[thread_id] = score
            thread_counts[thread_id] += 1
            delta = score - thread_means[thread_id]
            thread_means[thread_id] += delta / thread_counts[thread_id]
            delta2 = score - thread_means[thread_id]
            thread_M2s[thread_id] += delta * delta2

    out_min[0] = thread_mins[0]
    out_max[0] = thread_maxs[0]
    out_count[0] = 0

    for thread_id in range(num_threads):
        if thread_counts[thread_id] > 0:
            if thread_mins[thread_id] < out_min[0]:
                out_min[0] = thread_mins[thread_id]
            if thread_maxs[thread_id] > out_max[0]:
                out_max[0] = thread_maxs[thread_id]
            out_count[0] += thread_counts[thread_id]

    cdef double combined_mean = 0.0
    cdef double combined_M2 = 0.0
    cdef int64_t total_count = 0
    cdef double delta_means
    cdef int64_t new_total

    for thread_id in range(num_threads):
        if thread_counts[thread_id] > 0:
            if total_count == 0:
                combined_mean = thread_means[thread_id]
                combined_M2 = thread_M2s[thread_id]
                total_count = thread_counts[thread_id]
            else:
                delta_means = thread_means[thread_id] - combined_mean
                new_total = total_count + thread_counts[thread_id]
                combined_mean = (total_count * combined_mean + thread_counts[thread_id] * thread_means[thread_id]) / new_total
                combined_M2 = combined_M2 + thread_M2s[thread_id] + delta_means * delta_means * total_count * thread_counts[thread_id] / new_total
                total_count = new_total

    out_mean[0] = combined_mean
    if total_count > 1:
        out_variance[0] = combined_M2 / (total_count - 1)
    else:
        out_variance[0] = 0.0

    if total_count == 0:
        out_min[0] = 0.0
        out_max[0] = 0.0
        out_mean[0] = 0.0
        out_variance[0] = 0.0

    free(thread_mins)
    free(thread_maxs)
    free(thread_means)
    free(thread_M2s)
    free(thread_counts)


cdef int populate_memory_pool_direct(MemoryPool* pool,
                                     ProcessingBatch** batches,
                                     int64_t batch_count,
                                     sam_hdr_t* header,
                                     int num_threads) except -1 nogil:
    """Transfer alignments from batches to memory pool and sort by read ID.

    Copies batch alignments directly to pool, destroys batches to free memory,
    then performs in-place radix sort and rebuilds read indexing structures.

    Parameters
    ----------
    pool : MemoryPool*
        Target memory pool
    batches : ProcessingBatch**
        Array of processing batches (will be destroyed)
    batch_count : int64_t
        Number of batches
    header : sam_hdr_t*
        BAM header (unused, retained for API compatibility)
    num_threads : int
        Number of threads (unused, retained for API compatibility)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    bf_nogil_logf_notime(b"BATCH", "memory_pool: strategy=direct_copy_sort")
    cdef int64_t batch_idx, local_idx
    cdef int64_t total_count = 0
    cdef int64_t* counts = <int64_t*>malloc(<size_t>batch_count * sizeof(int64_t))
    cdef int64_t* offsets = NULL
    cdef int64_t i
    cdef int64_t write_idx
    cdef BatchAlignment* src = NULL
    cdef Alignment* dst = NULL
    if not counts:
        write_idx = 0
        for batch_idx in range(batch_count):
            if not batches[batch_idx]:
                continue
            for local_idx in range(batches[batch_idx].actual_alignment_count):
                if write_idx >= pool.alignment_capacity:
                    return -1
                src = &batches[batch_idx].batch_alignments[local_idx]
                dst = &pool.alignments[write_idx]
                dst.read_index = src.read_index
                dst.reference_index = src.reference_index
                dst.alignment_position = src.alignment_position
                dst.alignment_score = src.alignment_score
                dst.pmd_score = src.pmd_score
                # Copy ANI snapshot fields for damage-corrected ANI computation
                dst.aligned_length = src.aligned_length
                dst.match_count = src.match_count
                dst.ct_5p_count = src.ct_5p_count
                dst.ga_3p_count = src.ga_3p_count
                dst.c_at_5p_count = src.c_at_5p_count
                dst.g_at_3p_count = src.g_at_3p_count
                dst.damage_llr = src.damage_llr
                # Initialize PMD correction fields (computed after PMD fitting)
                dst.corrected_ani = 0.0
                dst.passes_ani_filter = 1  # Assume passes until corrected
                write_idx += 1
            destroy_processing_batch(batches[batch_idx])
            batches[batch_idx] = NULL
        pool.alignment_count = write_idx
        radix_sort_alignments_by_read_id(pool.alignments, pool.alignment_count)
    else:
        for batch_idx in range(batch_count):
            if not batches[batch_idx]:
                counts[batch_idx] = 0
            else:
                counts[batch_idx] = batches[batch_idx].actual_alignment_count
            total_count += counts[batch_idx]

        bf_nogil_logf_notime(b"BATCH", "memory_pool: total_count=%lld capacity=%lld",
                             <long long>total_count, <long long>pool.alignment_capacity)

        if total_count > pool.alignment_capacity:
            bf_nogil_logf_notime(b"ERROR", "Insufficient capacity: need %lld but have %lld",
                                <long long>total_count, <long long>pool.alignment_capacity)
            free(counts)
            return -1

        offsets = <int64_t*>malloc(<size_t>batch_count * sizeof(int64_t))
        if not offsets:
            free(counts)
            return -1

        offsets[0] = 0
        for i in range(1, batch_count):
            offsets[i] = offsets[i-1] + counts[i-1]

        for batch_idx in range(batch_count):
            if not batches[batch_idx] or counts[batch_idx] == 0:
                if batches[batch_idx]:
                    destroy_processing_batch(batches[batch_idx])
                    batches[batch_idx] = NULL
                continue
            src = batches[batch_idx].batch_alignments
            dst = &pool.alignments[offsets[batch_idx]]
            for local_idx in range(counts[batch_idx]):
                dst[local_idx].read_index = src[local_idx].read_index
                dst[local_idx].reference_index = src[local_idx].reference_index
                dst[local_idx].alignment_position = src[local_idx].alignment_position
                dst[local_idx].alignment_score = src[local_idx].alignment_score
                dst[local_idx].pmd_score = src[local_idx].pmd_score
                # Copy ANI snapshot fields for damage-corrected ANI computation
                dst[local_idx].aligned_length = src[local_idx].aligned_length
                dst[local_idx].match_count = src[local_idx].match_count
                dst[local_idx].ct_5p_count = src[local_idx].ct_5p_count
                dst[local_idx].ga_3p_count = src[local_idx].ga_3p_count
                dst[local_idx].c_at_5p_count = src[local_idx].c_at_5p_count
                dst[local_idx].g_at_3p_count = src[local_idx].g_at_3p_count
                dst[local_idx].damage_llr = src[local_idx].damage_llr
                # Initialize PMD correction fields (computed after PMD fitting)
                dst[local_idx].corrected_ani = 0.0
                dst[local_idx].passes_ani_filter = 1  # Assume passes until corrected
            destroy_processing_batch(batches[batch_idx])
            batches[batch_idx] = NULL

        pool.alignment_count = total_count
        free(counts)
        free(offsets)
    radix_sort_alignments_by_read_id(pool.alignments, pool.alignment_count)

    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    cdef uint32_t read_id
    for i in range(pool.alignment_count):
        read_id = pool.alignments[i].read_index
        if read_id < pool.unique_read_count:
            pool.read_alignment_counts[read_id] += 1
        else:
            return -1

    pool.read_alignment_starts[0] = 0
    for i in range(1, pool.unique_read_count):
        pool.read_alignment_starts[i] = (pool.read_alignment_starts[i-1] + pool.read_alignment_counts[i-1])

    pool.final_unique_reads = 0
    for i in range(pool.unique_read_count):
        if pool.read_alignment_counts[i] > 0:
            pool.final_unique_reads += 1

    return 0


cdef int64_t count_unique_reads_from_thread_maps(ThreadLocalHashMap** thread_maps, int num_threads) noexcept nogil:
    """Count unique reads from thread-local hash maps.

    Merges all thread-local hash maps to count globally unique read hash values.

    Parameters
    ----------
    thread_maps : ThreadLocalHashMap**
        Array of thread-local hash maps
    num_threads : int
        Number of threads

    Returns
    -------
    int64_t
        Number of unique reads, or -1 on error
    """
    cdef kh_seqid_map_t* global_unique_map = kh_init_seqid_map()
    cdef int64_t unique_count = 0
    cdef khint_t k
    cdef uint64_t hash_value
    cdef int ret, thread_id

    if not global_unique_map:
        return -1

    for thread_id in range(num_threads):
        if not thread_maps[thread_id]:
            continue

        for k in range(kh_end_seqid_map(thread_maps[thread_id].hash_to_id_map)):
            if kh_exist_seqid_map(thread_maps[thread_id].hash_to_id_map, k):
                hash_value = kh_key_seqid_map(thread_maps[thread_id].hash_to_id_map, k)
                kh_put_seqid_map(global_unique_map, hash_value, &ret)

    unique_count = kh_size(global_unique_map)
    kh_destroy_seqid_map(global_unique_map)

    bf_nogil_logf_notime(
        b"BATCH",
        "Identified %ld unique reads from alignments",
        unique_count,
    )
    return unique_count


cdef int64_t count_alignments_passing_ani_filter(ProcessingBatch** batches,
                                                  int64_t batch_count,
                                                  PMDCurve* curve,
                                                  float min_ani_threshold,
                                                  float epsilon) noexcept nogil:
    """Count alignments that pass the corrected ANI filter.

    Iterates through all batches and computes corrected ANI for each alignment
    using the fitted PMD curve. Returns count of alignments meeting threshold.

    Parameters
    ----------
    batches : ProcessingBatch**
        Array of processing batches with alignments
    batch_count : int64_t
        Number of batches
    curve : PMDCurve*
        Fitted PMD damage curve (or NULL for raw ANI)
    min_ani_threshold : float
        Minimum ANI percentage (0-100) to pass filter
    epsilon : float
        Baseline sequencing error rate

    Returns
    -------
    int64_t
        Number of alignments passing the filter
    """
    cdef int64_t total_passing = 0
    cdef int64_t batch_idx, aln_idx
    cdef ProcessingBatch* batch
    cdef BatchAlignment* src
    cdef ANISnapshot snapshot
    cdef float corrected_ani, raw_ani

    for batch_idx in range(batch_count):
        batch = batches[batch_idx]
        if not batch:
            continue

        for aln_idx in range(batch.actual_alignment_count):
            src = &batch.batch_alignments[aln_idx]

            if curve != NULL:
                # Compute corrected ANI using PMD curve
                snapshot.aligned_length = src.aligned_length
                snapshot.match_count = src.match_count
                snapshot.ct_5p_count = src.ct_5p_count
                snapshot.ga_3p_count = src.ga_3p_count
                snapshot.c_at_5p_count = src.c_at_5p_count
                snapshot.g_at_3p_count = src.g_at_3p_count
                snapshot.other_mm_count = 0
                snapshot.flags = 0

                corrected_ani = compute_corrected_ani(&snapshot, curve, epsilon)
            else:
                # No PMD curve - use raw ANI
                if src.aligned_length > 0:
                    raw_ani = (<float>src.match_count / <float>src.aligned_length) * 100.0
                else:
                    raw_ani = 0.0
                corrected_ani = raw_ani

            if corrected_ani >= min_ani_threshold:
                total_passing += 1

    return total_passing


cdef int populate_memory_pool_filtered(MemoryPool* pool,
                                        ProcessingBatch** batches,
                                        int64_t batch_count,
                                        sam_hdr_t* header,
                                        PMDCurve* curve,
                                        float min_ani_threshold,
                                        float epsilon,
                                        int num_threads) except -1 nogil:
    """Transfer only ANI-filtered alignments to memory pool.

    Computes corrected ANI for each alignment using the fitted PMD curve,
    and only transfers alignments that pass the threshold. This reduces
    memory pool size compared to transferring all then filtering.

    Parameters
    ----------
    pool : MemoryPool*
        Target memory pool (must be sized for filtered count)
    batches : ProcessingBatch**
        Array of processing batches (will be destroyed)
    batch_count : int64_t
        Number of batches
    header : sam_hdr_t*
        BAM header (unused, retained for API compatibility)
    curve : PMDCurve*
        Fitted PMD damage curve (or NULL for raw ANI)
    min_ani_threshold : float
        Minimum ANI percentage (0-100) to pass filter
    epsilon : float
        Baseline sequencing error rate
    num_threads : int
        Number of threads (unused, retained for API compatibility)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    bf_nogil_logf_notime(b"BATCH", "memory_pool: strategy=filtered_copy_sort (ANI>=%.1f%%)",
                         min_ani_threshold)

    cdef int64_t batch_idx, aln_idx
    cdef int64_t write_idx = 0
    cdef ProcessingBatch* batch
    cdef BatchAlignment* src
    cdef Alignment* dst
    cdef ANISnapshot snapshot
    cdef float corrected_ani, raw_ani
    cdef int64_t passed_count = 0
    cdef int64_t total_count = 0

    # Single pass: compute ANI and copy passing alignments
    for batch_idx in range(batch_count):
        batch = batches[batch_idx]
        if not batch:
            continue

        for aln_idx in range(batch.actual_alignment_count):
            total_count += 1
            src = &batch.batch_alignments[aln_idx]

            if curve != NULL:
                # Compute corrected ANI using PMD curve
                snapshot.aligned_length = src.aligned_length
                snapshot.match_count = src.match_count
                snapshot.ct_5p_count = src.ct_5p_count
                snapshot.ga_3p_count = src.ga_3p_count
                snapshot.c_at_5p_count = src.c_at_5p_count
                snapshot.g_at_3p_count = src.g_at_3p_count
                snapshot.other_mm_count = 0
                snapshot.flags = 0

                corrected_ani = compute_corrected_ani(&snapshot, curve, epsilon)
            else:
                # No PMD curve - use raw ANI
                if src.aligned_length > 0:
                    raw_ani = (<float>src.match_count / <float>src.aligned_length) * 100.0
                else:
                    raw_ani = 0.0
                corrected_ani = raw_ani

            # Skip alignments below threshold
            if corrected_ani < min_ani_threshold:
                continue

            # Check capacity
            if write_idx >= pool.alignment_capacity:
                bf_nogil_logf_notime(b"ERROR",
                    "Filtered pool overflow: write_idx=%lld capacity=%lld",
                    <long long>write_idx, <long long>pool.alignment_capacity)
                return -1

            # Copy alignment to pool
            dst = &pool.alignments[write_idx]
            dst.read_index = src.read_index
            dst.reference_index = src.reference_index
            dst.alignment_position = src.alignment_position
            dst.alignment_score = src.alignment_score
            dst.pmd_score = src.pmd_score
            dst.aligned_length = src.aligned_length
            dst.match_count = src.match_count
            dst.ct_5p_count = src.ct_5p_count
            dst.ga_3p_count = src.ga_3p_count
            dst.c_at_5p_count = src.c_at_5p_count
            dst.g_at_3p_count = src.g_at_3p_count
            dst.damage_llr = src.damage_llr
            # Store computed ANI values
            dst.corrected_ani = corrected_ani
            dst.passes_ani_filter = 1  # All transferred alignments pass

            write_idx += 1
            passed_count += 1

        # Destroy batch after processing
        destroy_processing_batch(batch)
        batches[batch_idx] = NULL

    pool.alignment_count = write_idx

    bf_nogil_logf_notime(b"BATCH",
        "memory_pool: filtered %lld/%lld alignments (%.1f%% passed)",
        <long long>passed_count, <long long>total_count,
        100.0 * <double>passed_count / <double>total_count if total_count > 0 else 0.0)

    # Sort by read ID for efficient indexing
    radix_sort_alignments_by_read_id(pool.alignments, pool.alignment_count)

    # Rebuild read indexing structures
    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    cdef uint32_t read_id
    cdef int64_t i
    for i in range(pool.alignment_count):
        read_id = pool.alignments[i].read_index
        if read_id < pool.unique_read_count:
            pool.read_alignment_counts[read_id] += 1
        else:
            bf_nogil_logf_notime(b"ERROR", "Invalid read_index %u >= unique_read_count %u",
                                 read_id, pool.unique_read_count)
            return -1

    pool.read_alignment_starts[0] = 0
    for i in range(1, pool.unique_read_count):
        pool.read_alignment_starts[i] = pool.read_alignment_starts[i-1] + pool.read_alignment_counts[i-1]

    # Count unique reads with alignments after filtering
    pool.final_unique_reads = 0
    for i in range(pool.unique_read_count):
        if pool.read_alignment_counts[i] > 0:
            pool.final_unique_reads += 1

    bf_nogil_logf_notime(b"BATCH", "memory_pool: %lld unique reads with filtered alignments",
                         <long long>pool.final_unique_reads)

    return 0


cdef int populate_memory_pool_filtered_split(MemoryPool* pool,
                                              ProcessingBatch** batches,
                                              int64_t batch_count,
                                              sam_hdr_t* header,
                                              PMDCurve* curve,
                                              float min_ani_threshold,
                                              float epsilon,
                                              int num_threads) except -1 nogil:
    """Transfer ANI-filtered alignments to memory pool using split array storage.

    Memory-optimized version that populates split arrays (AlignmentCore, read_indices,
    optional HierarchicalData, DamageCounts, BAMWriterAux) instead of the legacy
    monolithic Alignment array.

    Parameters
    ----------
    pool : MemoryPool*
        Target memory pool with split arrays allocated (use_split_arrays=True)
    batches : ProcessingBatch**
        Array of processing batches (will be destroyed)
    batch_count : int64_t
        Number of batches
    header : sam_hdr_t*
        BAM header (unused, retained for API compatibility)
    curve : PMDCurve*
        Fitted PMD damage curve (or NULL for raw ANI)
    min_ani_threshold : float
        Minimum ANI percentage (0-100) to pass filter
    epsilon : float
        Baseline sequencing error rate
    num_threads : int
        Number of threads (unused, retained for API compatibility)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    bf_nogil_logf_notime(b"BATCH", "memory_pool: strategy=split_filtered (ANI>=%.1f%%)",
                         min_ani_threshold)

    cdef int64_t batch_idx, aln_idx
    cdef int64_t write_idx = 0
    cdef ProcessingBatch* batch
    cdef BatchAlignment* src
    cdef AlignmentCore* core
    cdef HierarchicalData* hier
    cdef DamageCounts* dmg
    cdef BAMWriterAux* aux
    cdef ANISnapshot snapshot
    cdef float corrected_ani, raw_ani
    cdef int64_t passed_count = 0
    cdef int64_t total_count = 0

    # Precompute log factors for hierarchical EM log-likelihoods
    # Uses D_avg from PMD curve (average of positions 1-8) for sample-fitted model
    cdef double D_avg = 0.1
    cdef double log_p_damage_anc, log_p_survive_anc, log_1m_eps, log_eps_over_3
    cdef double observed_damage, survived, other_mm, actual_opp
    cdef int pmd_pos

    if curve != NULL:
        D_avg = 0.0
        for pmd_pos in range(8):
            D_avg += curve.D_5p_noncpg[pmd_pos] + curve.D_3p_noncpg[pmd_pos]
        D_avg /= 16.0

    log_p_damage_anc = log(fmax(D_avg + (1.0 - D_avg) * epsilon / 3.0, 1e-15))
    log_p_survive_anc = log(fmax((1.0 - D_avg) * (1.0 - epsilon / 3.0), 1e-15))
    log_1m_eps = log(fmax(1.0 - epsilon, 1e-15))
    log_eps_over_3 = log(fmax(epsilon / 3.0, 1e-15))

    if pool.hierarchical != NULL:
        bf_nogil_logf_notime(b"BATCH",
            "log_L precompute: D_avg=%.4f log_dmg=%.3f log_surv=%.3f log_1m_eps=%.3f",
            D_avg, log_p_damage_anc, log_p_survive_anc, log_1m_eps)

    # Phase 1: Filter and populate split arrays
    for batch_idx in range(batch_count):
        batch = batches[batch_idx]
        if not batch:
            continue

        for aln_idx in range(batch.actual_alignment_count):
            total_count += 1
            src = &batch.batch_alignments[aln_idx]

            if curve != NULL:
                snapshot.aligned_length = src.aligned_length
                snapshot.match_count = src.match_count
                snapshot.ct_5p_count = src.ct_5p_count
                snapshot.ga_3p_count = src.ga_3p_count
                snapshot.c_at_5p_count = src.c_at_5p_count
                snapshot.g_at_3p_count = src.g_at_3p_count
                snapshot.other_mm_count = 0
                snapshot.flags = 0
                corrected_ani = compute_corrected_ani(&snapshot, curve, epsilon)
            else:
                if src.aligned_length > 0:
                    raw_ani = (<float>src.match_count / <float>src.aligned_length) * 100.0
                else:
                    raw_ani = 0.0
                corrected_ani = raw_ani

            if corrected_ani < min_ani_threshold:
                continue

            if write_idx >= pool.alignment_capacity:
                bf_nogil_logf_notime(b"ERROR",
                    "Split pool overflow: write_idx=%lld capacity=%lld",
                    <long long>write_idx, <long long>pool.alignment_capacity)
                return -1

            # Populate AlignmentCore (always)
            core = &pool.alignment_cores[write_idx]
            core.reference_index = src.reference_index
            core.alignment_score = src.alignment_score
            core.alignment_position = src.alignment_position
            core.aligned_length = src.aligned_length
            core.match_count = src.match_count

            # Populate read_indices (always)
            pool.read_indices[write_idx] = src.read_index

            # Populate HierarchicalData (optional)
            if pool.hierarchical != NULL:
                hier = &pool.hierarchical[write_idx]
                hier.damage_llr = src.damage_llr

                # Compute log-likelihoods using sample-fitted D_avg (not hardcoded)
                # log_L_anc = obs_damage * log_p_damage + survived * log_p_survive
                #           + matches * log(1-eps) + other_mm * log(eps/3)
                observed_damage = <double>(src.ct_5p_count + src.ga_3p_count)
                actual_opp = <double>(src.c_at_5p_count + src.g_at_3p_count)
                survived = fmax(0.0, actual_opp - observed_damage)
                other_mm = <double>(src.aligned_length - src.match_count) - observed_damage
                if other_mm < 0.0:
                    other_mm = 0.0

                hier.log_L_anc = <float>(
                    observed_damage * log_p_damage_anc +
                    survived * log_p_survive_anc +
                    <double>src.match_count * log_1m_eps +
                    other_mm * log_eps_over_3
                )

                # log_L_mod = matches * log(1-eps) + all_mm * log(eps/3)
                hier.log_L_mod = <float>(
                    <double>src.match_count * log_1m_eps +
                    <double>(src.aligned_length - src.match_count) * log_eps_over_3
                )

            # Populate DamageCounts (optional)
            if pool.damage_counts != NULL:
                dmg = &pool.damage_counts[write_idx]
                dmg.ct_5p_count = src.ct_5p_count
                dmg.ga_3p_count = src.ga_3p_count
                dmg.c_at_5p_count = src.c_at_5p_count
                dmg.g_at_3p_count = src.g_at_3p_count

            # Populate BAMWriterAux (optional)
            if pool.bam_aux != NULL:
                aux = &pool.bam_aux[write_idx]
                aux.pmd_score = src.pmd_score
                aux.corrected_ani = corrected_ani

            write_idx += 1
            passed_count += 1

        destroy_processing_batch(batch)
        batches[batch_idx] = NULL

    pool.alignment_count = write_idx

    bf_nogil_logf_notime(b"BATCH",
        "split_pool: filtered %lld/%lld alignments (%.1f%% passed)",
        <long long>passed_count, <long long>total_count,
        100.0 * <double>passed_count / <double>total_count if total_count > 0 else 0.0)

    # Phase 2: Sort by read_index using packed sort keys
    cdef uint64_t* sort_keys = NULL
    cdef int64_t i, j
    cdef uint32_t read_id, orig_idx
    cdef AlignmentCore temp_core
    cdef HierarchicalData temp_hier
    cdef DamageCounts temp_dmg
    cdef BAMWriterAux temp_aux
    cdef uint32_t temp_read_idx
    cdef uint64_t current_key
    cdef int64_t cycle_start, current, next_pos

    if pool.alignment_count > 1:
        sort_keys = <uint64_t*>malloc(pool.alignment_count * sizeof(uint64_t))
        if not sort_keys:
            bf_nogil_logf_notime(b"ERROR", "Failed to allocate sort keys")
            return -1

        for i in range(pool.alignment_count):
            sort_keys[i] = (<uint64_t>pool.read_indices[i] << 32) | <uint32_t>i

        radix_sort_uint64(sort_keys, pool.alignment_count)

        # Apply permutation using cycle-following for in-place reordering

        # Mark processed entries by setting high bit of original index
        for cycle_start in range(pool.alignment_count):
            orig_idx = <uint32_t>(sort_keys[cycle_start] & <uint64_t>0xFFFFFFFF)
            if orig_idx == cycle_start or (sort_keys[cycle_start] & <uint64_t>0x80000000):
                continue

            # Save starting element
            temp_core = pool.alignment_cores[cycle_start]
            temp_read_idx = pool.read_indices[cycle_start]
            if pool.hierarchical != NULL:
                temp_hier = pool.hierarchical[cycle_start]
            if pool.damage_counts != NULL:
                temp_dmg = pool.damage_counts[cycle_start]
            if pool.bam_aux != NULL:
                temp_aux = pool.bam_aux[cycle_start]

            current = cycle_start
            while True:
                next_pos = <int64_t>(sort_keys[current] & <uint64_t>0x7FFFFFFF)
                sort_keys[current] = sort_keys[current] | <uint64_t>0x80000000  # Mark as processed

                if next_pos == cycle_start:
                    pool.alignment_cores[current] = temp_core
                    pool.read_indices[current] = temp_read_idx
                    if pool.hierarchical != NULL:
                        pool.hierarchical[current] = temp_hier
                    if pool.damage_counts != NULL:
                        pool.damage_counts[current] = temp_dmg
                    if pool.bam_aux != NULL:
                        pool.bam_aux[current] = temp_aux
                    break

                pool.alignment_cores[current] = pool.alignment_cores[next_pos]
                pool.read_indices[current] = pool.read_indices[next_pos]
                if pool.hierarchical != NULL:
                    pool.hierarchical[current] = pool.hierarchical[next_pos]
                if pool.damage_counts != NULL:
                    pool.damage_counts[current] = pool.damage_counts[next_pos]
                if pool.bam_aux != NULL:
                    pool.bam_aux[current] = pool.bam_aux[next_pos]

                current = next_pos

        free(sort_keys)

    # Phase 3: Rebuild read indexing structures
    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    for i in range(pool.alignment_count):
        read_id = pool.read_indices[i]
        if read_id < pool.unique_read_count:
            pool.read_alignment_counts[read_id] += 1
        else:
            bf_nogil_logf_notime(b"ERROR", "Invalid read_index %u >= unique_read_count %u",
                                 read_id, pool.unique_read_count)
            return -1

    pool.read_alignment_starts[0] = 0
    for i in range(1, pool.unique_read_count):
        pool.read_alignment_starts[i] = pool.read_alignment_starts[i-1] + pool.read_alignment_counts[i-1]

    pool.final_unique_reads = 0
    cdef int64_t total_counts = 0
    cdef uint32_t min_read_id = 0xFFFFFFFF
    cdef uint32_t max_read_id = 0
    for i in range(pool.unique_read_count):
        if pool.read_alignment_counts[i] > 0:
            pool.final_unique_reads += 1
        total_counts += pool.read_alignment_counts[i]

    cdef uint32_t min_ref_id = 0xFFFFFFFF
    cdef uint32_t max_ref_id = 0
    for i in range(pool.alignment_count):
        if pool.read_indices[i] < min_read_id:
            min_read_id = pool.read_indices[i]
        if pool.read_indices[i] > max_read_id:
            max_read_id = pool.read_indices[i]
        if pool.alignment_cores[i].reference_index < min_ref_id:
            min_ref_id = pool.alignment_cores[i].reference_index
        if pool.alignment_cores[i].reference_index > max_ref_id:
            max_ref_id = pool.alignment_cores[i].reference_index

    cdef uint64_t last_end = pool.read_alignment_starts[pool.unique_read_count - 1] + pool.read_alignment_counts[pool.unique_read_count - 1]

    bf_nogil_logf_notime(b"BATCH", "split_pool: %lld unique reads with filtered alignments",
                         <long long>pool.final_unique_reads)
    bf_nogil_logf_notime(b"BATCH",
        "split_pool VERIFY: sum_counts=%lld aln_count=%lld read_id_range=[%u,%u] ref_id_range=[%u,%u] unique_read_count=%u ref_count=%u last_end=%llu",
        <long long>total_counts, <long long>pool.alignment_count,
        min_read_id, max_read_id, min_ref_id, max_ref_id, pool.unique_read_count, pool.reference_count, <unsigned long long>last_end)

    return 0
