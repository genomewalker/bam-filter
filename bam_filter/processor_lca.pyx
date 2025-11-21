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
High-performance Lowest Common Ancestor (LCA) computation for BAM alignments.

This module replaces the legacy Python implementation (`bam_filter.lca`) with
a Cython pipeline that reuses the existing batch-processing infrastructure.
It reads alignments directly through htslib, aggregates read-to-reference
relationships, and summarises taxonomy assignments using the internal
`TaxonomyDB` representation.

Key improvements:
    * Streaming BAM ingestion using the same batching helpers as the main
      processor (no Python-level `pysam` iteration).
    * Memory-efficient accumulation of read/reference associations with
      contiguous buffers.
    * metaDMG-compatible strict path intersection so each trusted read
      contributes exactly once to its best-scoring reference.
"""

LOG_TAG = "LCA"
cdef bytes LOG_TAG_B = b"LCA"

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t, uint8_t
from libc.stdlib cimport malloc, calloc, free, qsort, realloc
from libc.string cimport memcpy, strcmp, strlen, memcmp, memset
from libc.stdio cimport FILE, fopen, fclose, fprintf
from libc.math cimport exp
from libc.time cimport time, time_t

from cython.parallel cimport prange, threadid
from bam_filter.processor_em cimport PREFETCH_WRITE

# Cache line size for padding to avoid false sharing
DEF CACHE_LINE_SIZE = 64
DEF CACHE_LINE_DOUBLES = 8  # 64 bytes / 8 bytes per double

# Log level constants for C logging
DEF LOG_LEVEL_INFO = 2

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    double bf_monotonic_seconds() nogil

import gzip
import os
import operator
from pathlib import Path

from cpython.unicode cimport PyUnicode_FromString
from libc.string cimport strdup

from bam_filter.taxonomy_db cimport (
    TaxonomyDatabase,
    AccessionMapping,
    TaxonomyDB,
    TaxNode,
    build_lineage_string_nogil,
    get_rank_id,
    get_rank_name,
    get_name_at_rank_nogil,
    AccessionMap,
    compute_lca_nogil as c_compute_lca,
    compute_lca_for_array_nogil,
)
from bam_filter import logging as bf_logging
from bam_filter import taxonomy_db as taxonomy_py

from bam_filter.processor cimport AlignmentScoringConfig, INVALID_SEQUENTIAL_ID
from bam_filter.processor_types cimport (
    samFile,
    sam_hdr_t,
    hts_idx_t,
    hts_itr_t,
    bam1_t,
    hts_open,
    hts_close,
    hts_set_threads,
    sam_index_load,
    hts_idx_destroy,
    hts_idx_get_stat,
    sam_hdr_read,
    sam_hdr_destroy,
    sam_hdr_tid2name,
    sam_hdr_name2tid,
    sam_hdr_tid2len,
    sam_hdr_nref,
    sam_itr_queryi,
    sam_itr_next,
    hts_itr_destroy,
    bam_init1,
    bam_destroy1,
    bam_aux_get,
    bam_get_qname,
)

from bam_filter.batch_utils cimport create_balanced_batches_greedy

from bam_filter.processor_hash cimport (
    ThreadLocalHashMap,
    create_thread_local_hash_map,
    destroy_thread_local_hash_map,
    kh_seqid_map_t,
    khint_t,
    kh_init_seqid_map,
    kh_destroy_seqid_map,
    kh_get_seqid_map,
    kh_put_seqid_map,
    kh_end_seqid_map,
    kh_exist_seqid_map,
    kh_key_seqid_map,
    kh_val_seqid_map_wrap,
    kh_size,
    extract_read_hash_identifier,
    kh_seqid_name_map_t,
    kh_get_seqid_name_map,
    kh_put_seqid_name_map,
    kh_end_seqid_name_map,
    kh_exist_seqid_name_map,
    kh_val_seqid_name_map_wrap,
)

from bam_filter.processor_batch cimport (
    count_unique_reads_from_thread_maps,
)

from bam_filter.processor_sort cimport radix_sort_uint64, extract_byte_32

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_log(const char* tag, const char* msg) nogil
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from "zlib.h":
    ctypedef void* gzFile
    gzFile gzopen(const char* path, const char* mode) nogil
    int gzclose(gzFile file) nogil
    int gzprintf(gzFile file, const char* format, ...) nogil
    int gzwrite(gzFile file, const void* buf, unsigned int len) nogil

# ------------------------------------------------------------------
# LCA-specific batch processing with slim data structure
# ------------------------------------------------------------------

# Provide inline wrapper for bam_aux2f (not declared in our shared pxd)
cdef extern from *:
    """
    #include "htslib/sam.h"
    static inline float bam_aux2f_wrapper(const uint8_t *s) {
        return bam_aux2f(s);
    }
    """
    float bam_aux2f_wrapper(const uint8_t* s) nogil

from bam_filter.processor_md_quality cimport (
    calculate_md_quality_score,
    alignment_passes_quality_filters,
)

# Define C static string for ZS tag to avoid Python bytes in nogil regions
cdef extern from *:
    """
    static const char ZS_TAG_C[3] = {'Z','S','\0'};
    """
    const char ZS_TAG_C[3]

# Lightweight batch container for LCA processing using slim alignment structure
cdef struct LCAAlignment:
    uint32_t read_index        # 4 bytes - global read ID
    uint32_t reference_index   # 4 bytes - reference ID
    float    alignment_score   # 4 bytes - log-likelihood score
    # Total: 12 bytes (vs BatchAlignment's 20 bytes)

# Per-read LCA result for output
cdef struct PerReadLCA:
    uint32_t read_id          # 4 bytes - global read ID
    int32_t lca_taxid         # 4 bytes - assigned taxid
    uint32_t num_alignments   # 4 bytes - alignment count for this read
    uint32_t norm_ref_index   # 4 bytes - chosen reference index for normalization
    uint8_t is_trusted        # 1 byte - 1 if rank >= threshold, 0 otherwise
    # Total: 17 bytes (padding to 20 bytes)

# Simple memory pool for LCA (mirrors processor MemoryPool pattern)
cdef struct LCAPool:
    LCAAlignment* alignments
    int64_t alignment_count
    int64_t alignment_capacity
    int64_t unique_read_count
    uint32_t* read_alignment_starts   # offset into alignments array for each read
    uint32_t* read_alignment_counts   # number of alignments per read

cdef struct LCABatch:
    int64_t batch_identifier
    int64_t reference_start_index
    int64_t reference_end_index
    int64_t expected_alignment_count
    int64_t actual_alignment_count
    LCAAlignment* alignments
    int64_t batch_capacity
    int32_t error_status
    int32_t processed_by_thread_id

cdef LCABatch* create_lca_batch(int64_t batch_id, int64_t ref_start,
                                 int64_t ref_end, int64_t expected_count) noexcept nogil:
    """Create an LCA batch with slim 12-byte alignments."""
    cdef LCABatch* batch = <LCABatch*>malloc(sizeof(LCABatch))
    if batch == NULL:
        return NULL
    
    batch.batch_identifier = batch_id
    batch.reference_start_index = ref_start
    batch.reference_end_index = ref_end
    batch.expected_alignment_count = expected_count
    batch.actual_alignment_count = 0
    batch.error_status = 0
    batch.processed_by_thread_id = -1
    
    # Allocate with 20% buffer for alignment count estimation errors
    cdef int64_t initial_capacity = <int64_t>(expected_count * 1.2) + 1024
    batch.batch_capacity = initial_capacity
    batch.alignments = <LCAAlignment*>malloc(initial_capacity * sizeof(LCAAlignment))
    
    if batch.alignments == NULL:
        free(batch)
        return NULL
    
    return batch

cdef void destroy_lca_batch(LCABatch* batch) noexcept nogil:
    """Free LCA batch memory."""
    if batch != NULL:
        if batch.alignments != NULL:
            free(batch.alignments)
        free(batch)

cdef LCAPool* create_lca_pool(int64_t total_alignments, int64_t unique_reads) noexcept nogil:
    """Create LCA memory pool."""
    cdef LCAPool* pool = <LCAPool*>malloc(sizeof(LCAPool))
    if pool == NULL:
        return NULL
    
    pool.alignment_capacity = total_alignments
    pool.alignment_count = 0
    pool.unique_read_count = unique_reads
    
    pool.alignments = <LCAAlignment*>malloc(<size_t>total_alignments * sizeof(LCAAlignment))
    pool.read_alignment_starts = <uint32_t*>calloc(<size_t>unique_reads, sizeof(uint32_t))
    pool.read_alignment_counts = <uint32_t*>calloc(<size_t>unique_reads, sizeof(uint32_t))
    
    if pool.alignments == NULL or pool.read_alignment_starts == NULL or pool.read_alignment_counts == NULL:
        if pool.alignments: free(pool.alignments)
        if pool.read_alignment_starts: free(pool.read_alignment_starts)
        if pool.read_alignment_counts: free(pool.read_alignment_counts)
        free(pool)
        return NULL
    
    return pool

cdef void destroy_lca_pool(LCAPool* pool) noexcept nogil:
    """Free LCA pool memory."""
    if pool != NULL:
        if pool.alignments: free(pool.alignments)
        if pool.read_alignment_starts: free(pool.read_alignment_starts)
        if pool.read_alignment_counts: free(pool.read_alignment_counts)
        free(pool)

cdef int grow_lca_batch_capacity(LCABatch* batch) noexcept nogil:
    """Double the capacity of an LCA batch."""
    cdef int64_t new_capacity = batch.batch_capacity * 2
    cdef LCAAlignment* new_alignments = <LCAAlignment*>malloc(new_capacity * sizeof(LCAAlignment))
    
    if new_alignments == NULL:
        return -1
    
    memcpy(new_alignments, batch.alignments, batch.actual_alignment_count * sizeof(LCAAlignment))
    free(batch.alignments)
    batch.alignments = new_alignments
    batch.batch_capacity = new_capacity
    
    return 0

cdef inline int32_t _taxid_to_idx(TaxonomyDB* db, int32_t taxid) nogil:
    if taxid < 0 or taxid > db.max_taxid:
        return -1
    return db.taxid_to_idx[taxid]

# Climb a taxonomy node up to the nearest ancestor at the requested rank.
# If the exact rank is not present in the lineage, return the original taxid.
cdef inline int32_t _climb_to_rank(TaxonomyDB* db, int32_t taxid, int32_t target_rank_id) nogil:
    cdef int32_t idx = _taxid_to_idx(db, taxid)
    cdef int32_t parent
    if idx < 0:
        return taxid
    while idx >= 0:
        if db.nodes[idx].rank_id == target_rank_id:
            return db.nodes[idx].taxid
        parent = db.nodes[idx].parent_taxid
        # Stop if at root or invalid
        if parent == db.nodes[idx].taxid or parent < 0 or parent > db.max_taxid:
            break
        idx = _taxid_to_idx(db, parent)
    return taxid

cdef int populate_lca_pool(LCAPool* pool, LCABatch** batches, int64_t batch_count) noexcept nogil:
    """Copy batches to pool, destroy batches, radix-sort by read_index, build read index."""
    cdef int64_t batch_idx, local_idx, i
    cdef int64_t write_idx = 0
    cdef LCABatch* batch
    cdef LCAAlignment* src
    cdef LCAAlignment* dst
    
    # Copy all batches to pool
    for batch_idx in range(batch_count):
        batch = batches[batch_idx]
        if batch == NULL or batch.actual_alignment_count == 0:
            continue
        if write_idx + batch.actual_alignment_count > pool.alignment_capacity:
            return -1
        src = batch.alignments
        dst = &pool.alignments[write_idx]
        memcpy(dst, src, <size_t>batch.actual_alignment_count * sizeof(LCAAlignment))
        write_idx += batch.actual_alignment_count
        destroy_lca_batch(batch)
        batches[batch_idx] = NULL
    
    pool.alignment_count = write_idx
    
    # Radix-sort by read_index
    radix_sort_lca_by_read(pool.alignments, pool.alignment_count)
    
    # Build read index (starts and counts)
    memset(pool.read_alignment_counts, 0, <size_t>pool.unique_read_count * sizeof(uint32_t))
    cdef uint32_t read_id
    for i in range(pool.alignment_count):
        read_id = pool.alignments[i].read_index
        if read_id >= <uint32_t>pool.unique_read_count:
            return -1
        pool.read_alignment_counts[read_id] += 1
    
    pool.read_alignment_starts[0] = 0
    for i in range(1, pool.unique_read_count):
        pool.read_alignment_starts[i] = pool.read_alignment_starts[i-1] + pool.read_alignment_counts[i-1]
    
    return 0

cdef int process_lca_batch_alignments(
        samFile* bam_file,
        sam_hdr_t* header,
        hts_idx_t* index,
        int64_t* reference_ids,
        LCABatch* batch,
        AlignmentScoringConfig* scoring_config,
        ThreadLocalHashMap* thread_map,
        int32_t thread_id) noexcept nogil:
    """
    Process alignments for LCA with ZS:f TAG reuse optimization.
    
    Checks for pre-computed ZS:f TAG (from filter/reassign step) and reuses it
    if available. Falls back to calculating score from MD tag if TAG is missing.
    
    This uses the slim 12-byte LCAAlignment structure (40% memory savings vs BatchAlignment).
    """
    cdef int64_t reference_id
    cdef hts_itr_t* iterator = NULL
    cdef bam1_t* bam_record = NULL
    cdef int return_code
    cdef uint64_t read_hash
    cdef khint_t k
    cdef int ret
    cdef uint32_t local_sequential_id
    cdef LCAAlignment temp_alignment
    cdef float raw_score
    cdef uint8_t* zs_aux
    cdef int64_t ref_idx
    cdef const char* qname_ptr
    cdef khint_t kn
    
    # Record which thread processed this batch
    batch.processed_by_thread_id = thread_id

    bam_record = bam_init1()
    if bam_record == NULL:
        return -1
    
    # Process each reference in this batch
    for ref_idx in range(batch.reference_start_index, batch.reference_end_index):
        reference_id = reference_ids[ref_idx]
        iterator = sam_itr_queryi(index, reference_id, 0, 0x7fffffff)
        if iterator == NULL:
            continue
        
        while True:
            return_code = sam_itr_next(bam_file, iterator, bam_record)
            if return_code < 0:
                break
            
            # Apply filters (identity, length, etc.)
            if not alignment_passes_quality_filters(bam_record, scoring_config):
                continue
            
            # Get or create read ID and cache read name
            read_hash = <uint64_t>extract_read_hash_identifier(bam_record)
            k = kh_get_seqid_map(thread_map.hash_to_id_map, read_hash)
            
            if k == kh_end_seqid_map(thread_map.hash_to_id_map):
                local_sequential_id = thread_map.next_local_id
                thread_map.next_local_id += 1
                k = kh_put_seqid_map(thread_map.hash_to_id_map, read_hash, &ret)
                if ret != -1:
                    kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0] = local_sequential_id
                # store name if not present
                qname_ptr = bam_get_qname(bam_record)
                if qname_ptr != NULL:
                    kn = kh_get_seqid_name_map(thread_map.hash_to_name_map, read_hash)
                    if kn == kh_end_seqid_name_map(thread_map.hash_to_name_map):
                        kn = kh_put_seqid_name_map(thread_map.hash_to_name_map, read_hash, &ret)
                        if ret != -1:
                            kh_val_seqid_name_map_wrap(thread_map.hash_to_name_map, kn)[0] = strdup(<const char*>qname_ptr)
            else:
                local_sequential_id = kh_val_seqid_map_wrap(thread_map.hash_to_id_map, k)[0]
            
            # Populate slim LCA alignment structure
            temp_alignment.read_index = local_sequential_id
            temp_alignment.reference_index = <uint32_t>reference_id
            
            # Try to reuse ZS:f TAG from filter/reassign step
            zs_aux = bam_aux_get(bam_record, ZS_TAG_C)
            if zs_aux != NULL and not scoring_config.calculate_pmd:
                # Reuse pre-computed score from TAG (fast path)
                raw_score = bam_aux2f_wrapper(zs_aux)
            else:
                # Calculate score from MD tag (fallback or when PMD needed)
                raw_score = calculate_md_quality_score(bam_record, header, scoring_config, <float*>NULL)
            
            temp_alignment.alignment_score = raw_score
            
            # Grow batch if needed
            if batch.actual_alignment_count >= batch.batch_capacity:
                if grow_lca_batch_capacity(batch) != 0:
                    batch.error_status = -1
                    hts_itr_destroy(iterator)
                    bam_destroy1(bam_record)
                    return -1
            
            batch.alignments[batch.actual_alignment_count] = temp_alignment
            batch.actual_alignment_count += 1
        
        hts_itr_destroy(iterator)
        iterator = NULL
    
    bam_destroy1(bam_record)
    return 0

cdef inline int64_t _binary_search_uint64(uint64_t* arr, int64_t n, uint64_t target) nogil:
    """Binary search for target in sorted array."""
    cdef int64_t left = 0
    cdef int64_t right = n - 1
    cdef int64_t mid
    while left <= right:
        mid = (left + right) // 2
        if arr[mid] == target:
            return mid
        elif arr[mid] < target:
            left = mid + 1
        else:
            right = mid - 1
    return -1

cdef int assign_lca_global_ids_fast(LCABatch** batches, int64_t batch_count,
                                    ThreadLocalHashMap** thread_maps, int num_threads) noexcept nogil:
    """Assign global sequential read IDs using the same fast strategy as processing batches.

    This mirrors assign_global_sequential_ids_fast but operates on LCABatch/LCAAlignment.
    """
    cdef kh_seqid_map_t* global_map = kh_init_seqid_map()
    cdef khint_t k, global_k
    cdef uint64_t hash_value
    cdef int ret, thread_id
    cdef uint32_t local_id, global_id

    if global_map == NULL:
        return -1

    # Collect all unique hashes into a dynamically growing array
    cdef int64_t hash_capacity = batch_count * 1000 if batch_count > 0 else 1024
    cdef uint64_t* all_hashes = <uint64_t*>malloc(hash_capacity * sizeof(uint64_t))
    cdef uint64_t* __tmp_hashes
    if all_hashes == NULL:
        kh_destroy_seqid_map(global_map)
        return -1
    cdef int64_t total_unique_hashes = 0

    for thread_id in range(num_threads):
        if thread_maps[thread_id] == NULL:
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
                            if __tmp_hashes == NULL:
                                kh_destroy_seqid_map(global_map)
                                free(all_hashes)
                                return -1
                            all_hashes = __tmp_hashes
                        all_hashes[total_unique_hashes] = hash_value
                        total_unique_hashes += 1
            k += 1

    if total_unique_hashes > 1:
        radix_sort_uint64(all_hashes, total_unique_hashes)

    # Rebuild global map: hash -> global_id (sorted order)
    kh_destroy_seqid_map(global_map)
    global_map = kh_init_seqid_map()
    if global_map == NULL:
        free(all_hashes)
        return -1
    for k in range(total_unique_hashes):
        hash_value = all_hashes[k]
        global_k = kh_put_seqid_map(global_map, hash_value, &ret)
        if ret != -1:
            kh_val_seqid_map_wrap(global_map, global_k)[0] = <uint32_t>k

    free(all_hashes)

    # Build per-thread local_id -> global_id tables
    cdef uint32_t** thread_local_to_global = <uint32_t**>malloc(num_threads * sizeof(uint32_t*))
    cdef uint32_t* thread_max_ids = <uint32_t*>malloc(num_threads * sizeof(uint32_t))
    cdef int t
    if thread_local_to_global == NULL or thread_max_ids == NULL:
        if thread_local_to_global != NULL: free(thread_local_to_global)
        if thread_max_ids != NULL: free(thread_max_ids)
        kh_destroy_seqid_map(global_map)
        return -1

    for t in range(num_threads):
        thread_max_ids[t] = 0
        if thread_maps[t] != NULL:
            k = 0
            while k < kh_end_seqid_map(thread_maps[t].hash_to_id_map):
                if kh_exist_seqid_map(thread_maps[t].hash_to_id_map, k):
                    local_id = kh_val_seqid_map_wrap(thread_maps[t].hash_to_id_map, k)[0]
                    if local_id > thread_max_ids[t]:
                        thread_max_ids[t] = local_id
                k += 1

    for t in range(num_threads):
        if thread_max_ids[t] > 0:
            thread_local_to_global[t] = <uint32_t*>malloc((thread_max_ids[t] + 1) * sizeof(uint32_t))
            if thread_local_to_global[t] == NULL:
                for thread_id in range(t):
                    if thread_local_to_global[thread_id] != NULL:
                        free(thread_local_to_global[thread_id])
                free(thread_local_to_global)
                free(thread_max_ids)
                kh_destroy_seqid_map(global_map)
                return -1
            memset(thread_local_to_global[t], INVALID_SEQUENTIAL_ID, (thread_max_ids[t] + 1) * sizeof(uint32_t))
        else:
            thread_local_to_global[t] = NULL

    for t in range(num_threads):
        if thread_maps[t] == NULL or thread_local_to_global[t] == NULL:
            continue
        k = 0
        while k < kh_end_seqid_map(thread_maps[t].hash_to_id_map):
            if kh_exist_seqid_map(thread_maps[t].hash_to_id_map, k):
                hash_value = kh_key_seqid_map(thread_maps[t].hash_to_id_map, k)
                local_id = kh_val_seqid_map_wrap(thread_maps[t].hash_to_id_map, k)[0]
                global_k = kh_get_seqid_map(global_map, hash_value)
                if global_k != kh_end_seqid_map(global_map):
                    global_id = kh_val_seqid_map_wrap(global_map, global_k)[0]
                    if local_id <= thread_max_ids[t]:
                        thread_local_to_global[t][local_id] = global_id
            k += 1

    # Remap alignments
    cdef int64_t batch_idx, aln_idx
    cdef int32_t actual_thread_id
    cdef LCABatch* batch
    cdef LCAAlignment* aln

    for batch_idx in range(batch_count):
        batch = batches[batch_idx]
        if batch == NULL:
            continue
        actual_thread_id = batch.processed_by_thread_id
        if actual_thread_id < 0 or actual_thread_id >= num_threads or thread_local_to_global[actual_thread_id] == NULL:
            continue
        for aln_idx in range(batch.actual_alignment_count):
            aln = &batch.alignments[aln_idx]
            local_id = aln.read_index
            if local_id <= thread_max_ids[actual_thread_id]:
                global_id = thread_local_to_global[actual_thread_id][local_id]
                if global_id != INVALID_SEQUENTIAL_ID:
                    aln.read_index = global_id

    for t in range(num_threads):
        if thread_local_to_global[t] != NULL:
            free(thread_local_to_global[t])
    free(thread_local_to_global)
    free(thread_max_ids)
    kh_destroy_seqid_map(global_map)

    return 0

cdef int64_t count_lca_alignments(LCABatch** batches, int64_t batch_count) noexcept nogil:
    """Count total alignments across all LCA batches."""
    cdef int64_t total = 0
    cdef int64_t batch_idx
    for batch_idx in range(batch_count):
        total += batches[batch_idx].actual_alignment_count
    return total

# ------------------------------------------------------------------
# Module-level helpers for sorting children by name (nogil)
# ------------------------------------------------------------------
cdef struct ChildNamePair:
    int32_t idx
    const char* name

cdef int _cmp_child_name_pair(const void* a, const void* b) noexcept nogil:
    cdef const ChildNamePair* pa = <const ChildNamePair*>a
    cdef const ChildNamePair* pb = <const ChildNamePair*>b
    return strcmp(pa.name, pb.name)

# Comparator to sort LCA alignments by read_index so that per-read groups are contiguous
cdef int _cmp_lca_alignment_by_read(const void* a, const void* b) noexcept nogil:
    cdef const LCAAlignment* pa = <const LCAAlignment*>a
    cdef const LCAAlignment* pb = <const LCAAlignment*>b
    if pa.read_index < pb.read_index:
        return -1
    elif pa.read_index > pb.read_index:
        return 1
    else:
        return 0

# --- American-flag radix sort for LCAAlignment by read_index (uint32) ---
cdef inline void _swap_lca(LCAAlignment* a, LCAAlignment* b) noexcept nogil:
    cdef LCAAlignment tmp = a[0]
    a[0] = b[0]
    b[0] = tmp

cdef inline void _insertion_sort_lca_by_read(LCAAlignment* a, int64_t left, int64_t right) noexcept nogil:
    cdef int64_t i, j
    cdef LCAAlignment key
    cdef uint32_t key_read
    for i in range(left + 1, right + 1):
        key = a[i]
        key_read = key.read_index
        j = i - 1
        while j >= left and a[j].read_index > key_read:
            a[j + 1] = a[j]
            j -= 1
        a[j + 1] = key

cdef void _radix_sort_lca_range(LCAAlignment* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil:
    cdef int64_t n = hi - lo
    if n <= 32:
        _insertion_sort_lca_by_read(a, lo, hi - 1)
        return
    if shift_bits < 0:
        return

    cdef int64_t count[256]
    cdef int64_t starts[256]
    cdef int64_t ends[256]
    cdef int64_t i, s, e
    cdef int bucket
    cdef uint32_t b

    for bucket in range(256):
        count[bucket] = 0

    for i in range(lo, hi):
        count[extract_byte_32(a[i].read_index, shift_bits)] += 1

    starts[0] = lo
    for bucket in range(1, 256):
        starts[bucket] = starts[bucket - 1] + count[bucket - 1]
    for bucket in range(256):
        ends[bucket] = starts[bucket] + count[bucket]

    bucket = 0
    while bucket < 256:
        s = starts[bucket]
        e = ends[bucket]
        while s < e:
            b = extract_byte_32(a[s].read_index, shift_bits)
            if b == bucket:
                s += 1
                starts[bucket] = s
            else:
                _swap_lca(&a[s], &a[starts[b]])
                starts[b] += 1
        bucket += 1

    if shift_bits > 0:
        for bucket in range(256):
            s = ends[bucket] - count[bucket]
            e = ends[bucket]
            if e - s > 1:
                _radix_sort_lca_range(a, s, e, shift_bits - 8)

cdef void radix_sort_lca_by_read(LCAAlignment* a, int64_t n) noexcept nogil:
    if n > 1:
        _radix_sort_lca_range(a, 0, n, 24)

cdef inline int _process_one_read_lca(
        LCAAlignment* read_aligns,
        uint32_t count,
        int32_t* reference_taxids,
        int32_t rank_id,
        TaxonomyDB* taxdb_c,
        double* w_tls,
        double* r_tls,
    int32_t* out_lca_taxid,
    uint32_t* out_best_ref_index) noexcept nogil:
    """Process LCA for one read's alignments. Returns 0 if accepted, -1 if discarded.
    
    Parameters
    ----------
    out_lca_taxid : int32_t*
        Output parameter: if accepted, set to the computed LCA taxid; if discarded, set to -1
    """
    # Stack buffers for typical case (<256 alignments per read)
    cdef int32_t tax_buffer[256]
    cdef float score_buffer[256]
    cdef uint32_t ref_index_buffer[256]
    cdef int32_t* tax_heap = NULL
    cdef float* score_heap = NULL
    cdef uint32_t* ref_heap = NULL

    cdef int dedup_count = 0
    cdef int best_idx = -1
    cdef double max_score = -1e100
    cdef int32_t current_taxid, computed_lca_local, lca_rank_id, min_depth_required
    cdef uint32_t ref_index, best_ref_index
    cdef float score_val
    cdef int k, j
    cdef int32_t tmp_idx, parent_taxid, tmp_idx2, parent_taxid2

    # Use heap if count > 256
    if count > 256:
        tax_heap = <int32_t*>malloc(count * sizeof(int32_t))
        score_heap = <float*>malloc(count * sizeof(float))
        ref_heap = <uint32_t*>malloc(count * sizeof(uint32_t))
        if tax_heap == NULL or score_heap == NULL or ref_heap == NULL:
            if tax_heap: free(tax_heap)
            if score_heap: free(score_heap)
            if ref_heap: free(ref_heap)
            return -1

    # Deduplicate by taxid, keep best score per taxid
    for k in range(count):
        ref_index = read_aligns[k].reference_index
        current_taxid = reference_taxids[ref_index]
        if current_taxid < 0:
            continue
        # Skip alignments whose taxon has "no rank" to mimic metaDMG-cpp behavior.
        if get_node_rank_id_nogil(taxdb_c, current_taxid) == 0:
            continue
        score_val = read_aligns[k].alignment_score

        # Check if taxid already seen
        for j in range(dedup_count):
            if count > 256:
                if tax_heap[j] == current_taxid:
                    if score_val > score_heap[j]:
                        score_heap[j] = score_val
                        ref_heap[j] = ref_index
                    break
            else:
                if tax_buffer[j] == current_taxid:
                    if score_val > score_buffer[j]:
                        score_buffer[j] = score_val
                        ref_index_buffer[j] = ref_index
                    break
        else:
            # New taxid
            if count > 256:
                tax_heap[dedup_count] = current_taxid
                score_heap[dedup_count] = score_val
                ref_heap[dedup_count] = ref_index
            else:
                tax_buffer[dedup_count] = current_taxid
                score_buffer[dedup_count] = score_val
                ref_index_buffer[dedup_count] = ref_index
            if score_val > max_score:
                max_score = score_val
                best_idx = dedup_count
            dedup_count += 1

    if dedup_count == 0:
        if count > 256:
            free(tax_heap); free(score_heap); free(ref_heap)
        out_lca_taxid[0] = -1
        return -1
    
    # LCA computation and validation
    cdef int32_t lca_idx = -1
    cdef int32_t lca_depth = -1

    if dedup_count > 1:
        if count > 256:
            computed_lca_local = _compute_lca_for_buffer(taxdb_c, tax_heap, dedup_count)
        else:
            computed_lca_local = _compute_lca_for_buffer(taxdb_c, tax_buffer, dedup_count)
        if computed_lca_local < 0:
            if count > 256:
                free(tax_heap); free(score_heap); free(ref_heap)
            return -1
        # Convert taxid to index
        lca_idx = _taxid_to_idx(taxdb_c, computed_lca_local)
        if lca_idx < 0:
            if count > 256:
                free(tax_heap); free(score_heap); free(ref_heap)
            return -1
        lca_rank_id = get_node_rank_id_nogil(taxdb_c, computed_lca_local)
        lca_depth = taxdb_c.nodes[lca_idx].depth
    else:
        # Single taxid - it's its own LCA
        if count > 256:
            lca_idx = _taxid_to_idx(taxdb_c, tax_heap[0])
        else:
            lca_idx = _taxid_to_idx(taxdb_c, tax_buffer[0])
        if lca_idx < 0:
            if count > 256:
                free(tax_heap); free(score_heap); free(ref_heap)
            out_lca_taxid[0] = -1
            return -1
        lca_rank_id = get_node_rank_id_nogil(taxdb_c, lca_idx)
        lca_depth = taxdb_c.nodes[lca_idx].depth
    
    # Reject if LCA computation failed
    if lca_rank_id < 0:
        if count > 256:
            free(tax_heap); free(score_heap); free(ref_heap)
        out_lca_taxid[0] = -1
        return -1

    # Count read once at best reference (strict mode behavior)
    if count > 256:
        best_ref_index = ref_heap[best_idx]
        free(tax_heap); free(score_heap); free(ref_heap)
    else:
        best_ref_index = ref_index_buffer[best_idx]

    out_lca_taxid[0] = taxdb_c.nodes[lca_idx].taxid
    w_tls[best_ref_index] += 1.0
    r_tls[best_ref_index] += 1.0
    out_best_ref_index[0] = best_ref_index

    return 0

cdef void process_lca_from_pool_parallel(
        LCAPool* pool,
        int32_t* reference_taxids,
        int32_t rank_id,
        TaxonomyDB* taxdb_c,
        double* ref_weight_tls,
        double* ref_read_tls,
        uint64_t* accepted_tls,
        uint64_t* discarded_tls,
        int64_t refs_per_thread_padded,
        int num_threads,
    PerReadLCA** per_read_results_tls) noexcept nogil:
    """Process per-read LCA from sorted pool in parallel with TLS accumulators.
    
    If per_read_results_tls is not NULL, collect per-read LCA assignments for output.
    """
    cdef int64_t read_id
    cdef uint32_t start_offset, count
    cdef LCAAlignment* read_aligns
    cdef int64_t tid
    cdef double* w_tls
    cdef double* r_tls
    cdef int result
    cdef int32_t lca_taxid
    cdef uint32_t best_ref_idx_out
    cdef int32_t lca_rank_id_actual

    for read_id in prange(pool.unique_read_count, num_threads=num_threads, schedule='static', nogil=True):
        count = pool.read_alignment_counts[read_id]
        if count == 0:
            continue
        
        start_offset = pool.read_alignment_starts[read_id]
        read_aligns = &pool.alignments[start_offset]
        tid = threadid()
        w_tls = &ref_weight_tls[tid * refs_per_thread_padded]
        r_tls = &ref_read_tls[tid * refs_per_thread_padded]
        result = _process_one_read_lca(read_aligns, count, reference_taxids, rank_id, taxdb_c, w_tls, r_tls, &lca_taxid, &best_ref_idx_out)
        
        if result == 0:
            accepted_tls[tid * CACHE_LINE_DOUBLES] += 1
            # If per-read output requested, collect result
            if per_read_results_tls != NULL and per_read_results_tls[tid] != NULL:
                per_read_results_tls[tid][read_id].read_id = <uint32_t>read_id
                # Report actual LCA taxid (metaDMG-compatible: no rank climbing)
                # This preserves full taxonomic resolution for combining with metaDMG damage estimates
                per_read_results_tls[tid][read_id].lca_taxid = lca_taxid
                per_read_results_tls[tid][read_id].num_alignments = count
                per_read_results_tls[tid][read_id].norm_ref_index = best_ref_idx_out
                # Determine if this read's LCA rank meets the threshold (trusted for stats)
                # IMPORTANT: Lower rank_id = more specific (subspecies=1, species=2, genus=6)
                # A read is trusted if its rank is at or BELOW (<=) the threshold rank
                lca_rank_id_actual = get_node_rank_id_nogil(taxdb_c, lca_taxid)
                if lca_rank_id_actual <= rank_id and lca_rank_id_actual > 0:
                    per_read_results_tls[tid][read_id].is_trusted = 1
                else:
                    per_read_results_tls[tid][read_id].is_trusted = 0
        else:
            discarded_tls[tid * CACHE_LINE_DOUBLES] += 1
            # Also record discarded reads in per-read output (if requested)
            if per_read_results_tls != NULL and per_read_results_tls[tid] != NULL:
                per_read_results_tls[tid][read_id].read_id = <uint32_t>read_id
                per_read_results_tls[tid][read_id].lca_taxid = -1  # Invalid/discarded
                per_read_results_tls[tid][read_id].num_alignments = count
                per_read_results_tls[tid][read_id].norm_ref_index = 0
                per_read_results_tls[tid][read_id].is_trusted = 0

cdef extern from "stdio.h" nogil:
    int snprintf(char* s, size_t n, const char* format, ...)
    size_t fwrite(const void* ptr, size_t size, size_t nmemb, FILE* stream)
    int fprintf(FILE* stream, const char* format, ...)
    int fflush(FILE* stream)

cdef extern from "zlib.h" nogil:
    int gzwrite(gzFile file, const void* buf, unsigned int len)

# Buffer for writing output - write in 1MB chunks
cdef int _write_tree(TaxonomyDB* db,
                     int32_t idx_node,
                     double* agg_counts,
                     double* w_sums,
                     int32_t* offsets,
                     int32_t* kids,
                     bint use_gz,
                     gzFile gzfp,
                     FILE* fp,
                     char* lineage_buf) nogil:
    cdef int local_wrote = 0
    cdef TaxNode node = db.nodes[idx_node]
    cdef const char* nm_ptr = get_name_at_rank_nogil(db, node.taxid, node.rank_id)
    cdef const char* rank_ptr
    cdef int line_len
    
    if nm_ptr == NULL:
        nm_ptr = b""
    if build_lineage_string_nogil(db, node.taxid, lineage_buf, 4096) < 0:
        lineage_buf[0] = 0
    if db.rank_names != NULL:
        rank_ptr = db.rank_names[node.rank_id]
    else:
        rank_ptr = b""
    
    # Write directly
    if use_gz:
        gzprintf(gzfp, "%d\t%s\t%s\t%.0f\t%.6f\t%s\n", 
                 node.taxid, nm_ptr, rank_ptr, agg_counts[idx_node], w_sums[idx_node], lineage_buf)
    else:
        fprintf(fp, "%d\t%s\t%s\t%.0f\t%.6f\t%s\n",
                node.taxid, nm_ptr, rank_ptr, agg_counts[idx_node], w_sums[idx_node], lineage_buf)
    
    local_wrote += 1
    cdef int32_t s = offsets[idx_node]
    cdef int32_t e = offsets[idx_node + 1]
    cdef int32_t k
    for k in range(e - s):
        local_wrote += _write_tree(db, kids[s + k], agg_counts, w_sums, offsets, kids, use_gz, gzfp, fp, lineage_buf)
    return local_wrote

from bam_filter.reference_lengths cimport (
    TSVReferenceMap,
    load_tsv_reference_file,
    get_tsv_reference_count,
    lookup_reference_length,
    free_tsv_reference_map,
    print_tsv_reference_stats,
)

# -----------------------------------------------------------------------------
# Local helper data structures
# -----------------------------------------------------------------------------

# Slim 12-byte structure for LCA-specific alignment data (vs 20-byte BatchAlignment)
# Saves 40% memory: only stores what LCA needs (read, reference, score)
    


cdef class _LCAStatistics:
    """
    Lightweight container exposed to Python for unit testing.
    Stores per-node read counts and weights along with metadata.
    """
    cdef public dict rows
    cdef public dict metadata

    def __init__(self, dict rows, dict metadata):
        self.rows = rows
        self.metadata = metadata


cdef dict _load_stats_file(str stats_path):
    cdef dict result = {}
    if not stats_path:
        return None

    if stats_path.endswith(".gz"):
        fh = gzip.open(stats_path, "rt")
    else:
        fh = open(stats_path, "r")

    with fh as handle:
        header = None
        for line in handle:
            line = line.strip()
            if not line:
                continue
            header = line.split("\t")
            break

        if header is None:
            return {}

        idx_ref = header.index("reference") if "reference" in header else -1
        if idx_ref < 0:
            raise ValueError("Stats file missing 'reference' column")

        idx_reads = header.index("n_reads") if "n_reads" in header else -1
        idx_reads_tad = header.index("n_reads_tad") if "n_reads_tad" in header else -1
        idx_cov = header.index("coverage_mean_trunc_len") if "coverage_mean_trunc_len" in header else -1

        for line in handle:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= idx_ref:
                continue
            ref_name = parts[idx_ref]
            n_reads = float(parts[idx_reads]) if idx_reads >= 0 and len(parts) > idx_reads and parts[idx_reads] else 0.0
            n_reads_tad = float(parts[idx_reads_tad]) if idx_reads_tad >= 0 and len(parts) > idx_reads_tad and parts[idx_reads_tad] else 0.0
            coverage = float(parts[idx_cov]) if idx_cov >= 0 and len(parts) > idx_cov and parts[idx_cov] else 0.0
            result[ref_name] = (n_reads, n_reads_tad, coverage)

    return result


cdef inline int32_t _compute_lca_for_buffer(TaxonomyDB* db, int32_t* taxids, int32_t n) nogil:
    """
    Compute LCA for a buffer of taxids using the correct ngsLCA.cpp algorithm.

    This now calls compute_lca_for_array_nogil which implements the proper
    multi-way LCA computation instead of the INCORRECT iterative pairwise approach.
    """
    return compute_lca_for_array_nogil(db, taxids, n)


cdef inline double _safe_exp(double value) nogil:
    if value > 60.0:
        return exp(60.0)
    elif value < -60.0:
        return exp(-60.0)
    return exp(value)


cdef inline int32_t get_node_rank_id_nogil(TaxonomyDB* db, int32_t taxid) nogil:
    cdef int32_t idx = _taxid_to_idx(db, taxid)
    if idx < 0:
        return -1
    return db.nodes[idx].rank_id


cdef inline void _propagate_up(TaxonomyDB* db, int32_t taxid, double read_inc, double weight,
                               double* read_counts, double* weight_sums) noexcept nogil:
    """
    Propagate both read counts and weights up the taxonomy tree.
    
    Each read is counted at every node in its lineage (species, genus, family, etc.)
    so that parent nodes reflect the total reads assigned to them and their descendants.
    """
    cdef int32_t idx = _taxid_to_idx(db, taxid)
    cdef int32_t parent
    
    # Propagate both reads and weights up the entire lineage
    while idx >= 0:
        read_counts[idx] += read_inc
        weight_sums[idx] += weight
        parent = db.nodes[idx].parent_taxid
        if parent == db.nodes[idx].taxid:
            break
        idx = _taxid_to_idx(db, parent)


cdef inline void _propagate_length(TaxonomyDB* db, int32_t taxid, double length,
                                   double* length_sums) noexcept nogil:
    cdef int32_t idx = _taxid_to_idx(db, taxid)
    cdef int32_t parent
    while idx >= 0:
        length_sums[idx] += length
        parent = db.nodes[idx].parent_taxid
        if parent == db.nodes[idx].taxid:
            break
        idx = _taxid_to_idx(db, parent)


cdef inline int32_t _get_rank_depth(TaxonomyDB* db, int32_t taxid) nogil:
    """
    Count the number of taxonomic ranks from root to the given taxid.
    This is used to determine if the LCA reaches the minimum required rank depth.
    """
    cdef int32_t idx = _taxid_to_idx(db, taxid)
    cdef int32_t depth = 0
    cdef int32_t parent
    
    while idx >= 0:
        depth += 1
        parent = db.nodes[idx].parent_taxid
        if parent == db.nodes[idx].taxid:
            break
        idx = _taxid_to_idx(db, parent)
    
    return depth


cdef inline double _compute_best_probability(float* scores, int n) nogil:
    cdef int i
    cdef double max_score = -1e100
    for i in range(n):
        if scores[i] > max_score:
            max_score = scores[i]
    cdef double sum_exp = 0.0
    for i in range(n):
        sum_exp += _safe_exp(scores[i] - max_score)
    if sum_exp <= 0.0:
        return 1.0
    return _safe_exp(max_score - max_score) / sum_exp  # equals 1/sum_exp


# ------------------------------------------------------------------
# Top-level helpers (avoid closures inside cpdef)
# ------------------------------------------------------------------

cdef object _node_name_for_sort_py(TaxonomyDB* taxdb_c, int idx_node):
    """Return a Python string for a node's display name used for sorting."""
    cdef const char* nm = get_name_at_rank_nogil(taxdb_c, taxdb_c.nodes[idx_node].taxid, taxdb_c.nodes[idx_node].rank_id)
    if nm != NULL:
        return PyUnicode_FromString(nm)
    return ""

cdef int _write_node_recursive_py(object fp,
                                  TaxonomyDB* taxdb_c,
                                  double* agg_read_counts,
                                  double* weight_sums,
                                  dict parent_to_children,
                                  int idx_node) except -1:
    """Write one node and its subtree; return number of rows written."""
    cdef TaxNode node
    cdef const char* node_name_ptr
    cdef char lineage_buf[4096]
    cdef str name_str
    cdef str rank_str
    cdef str lineage_str
    cdef list children_list
    cdef list sort_pairs
    cdef int wrote = 0

    node = taxdb_c.nodes[idx_node]
    node_name_ptr = get_name_at_rank_nogil(taxdb_c, node.taxid, node.rank_id)
    if node_name_ptr != NULL:
        name_str = PyUnicode_FromString(node_name_ptr)
    else:
        name_str = ""
    rank_str = get_rank_name(node.rank_id)
    if build_lineage_string_nogil(taxdb_c, node.taxid, lineage_buf, 4096) < 0:
        lineage_buf[0] = 0
    lineage_str = PyUnicode_FromString(lineage_buf)
    fp.write(f"{node.taxid}\t{name_str}\t{rank_str}\t{agg_read_counts[idx_node]:.0f}\t{weight_sums[idx_node]:.6f}\t{lineage_str}\n")
    wrote += 1

    if idx_node in parent_to_children:
        children_list = parent_to_children[idx_node]
        # Build (name, idx) pairs to avoid lambdas/closures in cpdef context
        sort_pairs = [(_node_name_for_sort_py(taxdb_c, ch), ch) for ch in children_list]
        sort_pairs.sort()
        for _nm, ch in sort_pairs:
            wrote += _write_node_recursive_py(fp, taxdb_c, agg_read_counts, weight_sums, parent_to_children, ch)
    return wrote


# ------------------------------------------------------------------
# Python helper functions for logging and UI
# ------------------------------------------------------------------

def _announce_stage(title: str, detail: str = "") -> None:
    """Emit a human-readable stage separator for pipeline progress."""
    bf_logging.summary("")
    bf_logging.summary("┌─ %s", title)
    if detail:
        bf_logging.summary("│ %s", detail)
    bf_logging.summary("└─────────────────────────────────────────────────────────────")


def _log_phase_duration(phase: str, double start_time) -> None:
    """Log the duration of a completed phase."""
    cdef double duration = bf_monotonic_seconds() - start_time
    bf_logging.summary("│ Phase completed in %.2f seconds", duration)
    bf_logging.summary("")


# ------------------------------------------------------------------
# Main LCA processing function
# ------------------------------------------------------------------

cpdef tuple run_lca(
    str bam_path,
    str output_path,
    str rank="genus",
    bint custom_acc=False,
    int threads=1,
    double min_read_ani=0.0,
    int min_read_length=30,
    int min_read_count=1,
    long scale=1000000,
    bint verbose=True,
    bint calculate_pmd=False,
    str stats_path=None,
    str reference_lengths_tsv=None,
    object taxonomy_db_dir=None,
    object taxonomy_accession_map=None,
    object nodes_path=None,
    object names_path=None,
    object acc2taxid_path=None,
    str per_read_path=None
):
    """
    Execute Cython LCA pipeline and write summary results to TSV/GZ.

    Parameters
    ----------
    calculate_pmd : bool, optional
        If True, calculate PMD scores (slower but available for downstream use).
        If False, skip PMD calculation and use ZS:f TAG if available (faster).
        Default: False
    per_read_path : str, optional
        Path to write per-read LCA assignments (TSV format). If the path ends with
        '.gz', output will be gzip-compressed. Columns: read_name, taxid, rank, n_aln, tax_path.
        Default: None (per-read output disabled)
        Note: This feature is currently under development and will be available in a future release.
    
    Other parameters mirror the old Python implementation but the heavy lifting
    happens in compiled code.
    """
    cdef bytes bam_bytes = bam_path.encode("utf-8")
    cdef const char* bam_file_path = bam_bytes
    cdef bytes output_bytes = output_path.encode("utf-8")
    cdef const char* output_c = output_bytes

    if threads < 1:
        threads = 1

    cdef double scale_factor = <double>scale
    if scale_factor <= 0.0:
        bf_logging.warn("Invalid scale factor %.2f; defaulting to 1e6", scale_factor)
        scale_factor = 1000000.0

    # Per-read output setup
    cdef bint write_per_read = False
    cdef bytes per_read_bytes
    cdef const char* per_read_c = NULL
    cdef int cleanup_idx
    cdef int64_t tid_idx, read_id, k
    cdef uint64_t read_hash
    # per-read name mapping helpers (declare upfront for Cython)
    cdef object read_names = None
    cdef kh_seqid_map_t* global_map_names = NULL
    cdef khint_t k_g = 0
    cdef khint_t k_t = 0
    cdef khint_t k_n = 0
    cdef uint64_t hv = 0
    cdef int ret_g = 0
    cdef int64_t capacity = 0
    cdef uint64_t* all_hashes = NULL
    cdef uint64_t* tmp_hashes = NULL
    cdef int64_t total_hashes = 0
    cdef char* name_ptr = NULL
    cdef uint32_t gidx = 0
    # Buffer for building lineage strings when writing per-read output
    cdef char lineage_buf_perread[4096]
    if per_read_path is not None and per_read_path != "":
        write_per_read = True
        per_read_bytes = per_read_path.encode("utf-8")
        per_read_c = per_read_bytes
        if verbose:
            bf_logging.log(LOG_TAG, "Per-read LCA output will be written to: %s", per_read_path)

    cdef TSVReferenceMap* tsv_map = NULL
    cdef bytes ref_len_bytes
    cdef const char* ref_len_c = NULL
    cdef double tsv_start, tsv_end, tax_start, tax_end, bam_idx_start, bam_idx_end

    # Build detail message based on what's being loaded
    cdef str init_detail = "Loading taxonomy databases"
    if reference_lengths_tsv is not None and reference_lengths_tsv != "":
        init_detail = "Loading reference lengths and taxonomy databases"

    cdef double phase0_start = bf_monotonic_seconds()
    _announce_stage("Initialization", init_detail)

    if reference_lengths_tsv is not None and reference_lengths_tsv != "":
        if not os.path.exists(reference_lengths_tsv):
            bf_logging.warn("Reference lengths TSV not found: %s", reference_lengths_tsv)
        else:
            if verbose:
                bf_logging.log(LOG_TAG, "Loading reference length overrides from %s", reference_lengths_tsv)
            tsv_start = bf_monotonic_seconds()
            ref_len_bytes = reference_lengths_tsv.encode("utf-8")
            ref_len_c = ref_len_bytes
            with nogil:
                tsv_map = load_tsv_reference_file(ref_len_c)
            tsv_end = bf_monotonic_seconds()
            if tsv_map == NULL:
                bf_logging.warn("Failed to load reference length TSV; falling back to BAM header lengths")
            else:
                with nogil:
                    print_tsv_reference_stats(tsv_map)
                if verbose:
                    bf_logging.log(LOG_TAG, "TSV loading time: %.2f seconds", tsv_end - tsv_start)

    cdef str taxonomy_db_dir_str = None
    cdef str nodes_path_str = None
    cdef str names_path_str = None
    cdef str accession_map_path_str = None
    cdef TaxonomyDatabase tax_db_obj

    if taxonomy_db_dir is not None and taxonomy_db_dir != "":
        taxonomy_db_dir_str = os.fspath(taxonomy_db_dir)
        if verbose:
            bf_logging.log(LOG_TAG, "Loading taxonomy database from %s", taxonomy_db_dir_str)
        tax_start = bf_monotonic_seconds()
        tax_db_obj = TaxonomyDatabase.from_parquet(taxonomy_db_dir_str)
        tax_end = bf_monotonic_seconds()
        if verbose:
            bf_logging.log(LOG_TAG, "Taxonomy loading time: %.2f seconds", tax_end - tax_start)
            bf_logging.log(LOG_TAG, "About to process BAM file at time %.2fs into Phase 0", tax_end - phase0_start)
    else:
        if nodes_path is None or names_path is None:
            raise ValueError(
                "Either taxonomy_db_dir or both nodes_path and names_path must be provided"
            )
        nodes_path_str = os.fspath(nodes_path)
        names_path_str = os.fspath(names_path)
        if verbose:
            bf_logging.log(
                LOG_TAG,
                "Loading taxonomy database from nodes file %s and names file %s",
                nodes_path_str,
                names_path_str,
            )
        tax_db_obj = taxonomy_py.load_taxonomy_from_ncbi(
            nodes_path_str, names_path_str, num_threads=threads
        )

    if tax_db_obj is None:
        raise RuntimeError("Failed to load taxonomy database")

    if taxonomy_db_dir_str is not None:
        if taxonomy_accession_map is not None and taxonomy_accession_map != "":
            accession_map_path_str = os.fspath(taxonomy_accession_map)
        elif acc2taxid_path is not None and acc2taxid_path != "":
            accession_map_path_str = os.fspath(acc2taxid_path)
        else:
            candidate_map = os.path.join(taxonomy_db_dir_str, "accession_map.parquet")
            if os.path.exists(candidate_map):
                accession_map_path_str = candidate_map
            else:
                raise ValueError(
                    "No accession_map.parquet found in the taxonomy database directory."
                )
    else:
        if acc2taxid_path is None or acc2taxid_path == "":
            raise ValueError("acc2taxid_path is required when taxonomy_db_dir is not provided")
        accession_map_path_str = os.fspath(acc2taxid_path)

    # ------------------------------------------------------------------
    # BAM ingestion - open BAM first to get reference list for filtering
    # ------------------------------------------------------------------
    cdef samFile* bam_handle = NULL
    cdef sam_hdr_t* bam_header = NULL
    cdef hts_idx_t* bam_index = NULL
    cdef int64_t total_references = 0
    cdef int64_t* reference_alignment_counts = NULL
    cdef int64_t* reference_ids = NULL
    cdef int64_t* reference_lengths = NULL
    cdef int64_t* batch_starts = NULL
    cdef int64_t* batch_ends = NULL
    cdef LCABatch** lca_batches = NULL
    cdef ThreadLocalHashMap** thread_maps = NULL
    cdef int64_t batch_count = 0
    cdef int64_t batch_idx
    cdef int32_t thread_idx
    cdef int64_t i
    cdef const char* ref_name
    cdef bint reference_lengths_from_tsv = False
    cdef int32_t nentries, tid32
    cdef int64_t tsv_length
    
    cdef AlignmentScoringConfig scoring_config
    scoring_config.minimum_read_identity = min_read_ani
    scoring_config.minimum_read_length = min_read_length
    scoring_config.maximum_read_length = 100000
    scoring_config.global_min_score = 1e30
    scoring_config.global_max_score = -1e30
    scoring_config.calculate_pmd = calculate_pmd
    scoring_config.is_single_stranded = False

    bam_handle = hts_open(bam_file_path, b"r")
    if not bam_handle:
        raise RuntimeError(f"Failed to open BAM file: {bam_path}")

    if threads > 1:
        hts_set_threads(bam_handle, threads)

    bam_header = sam_hdr_read(bam_handle)
    if not bam_header:
        hts_close(bam_handle)
        raise RuntimeError("Failed to read BAM header")

    total_references = sam_hdr_nref(bam_header)
    if total_references <= 0:
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise RuntimeError("BAM header reports zero references")

    reference_alignment_counts = <int64_t*>calloc(total_references, sizeof(int64_t))
    reference_ids = <int64_t*>malloc(total_references * sizeof(int64_t))
    reference_lengths = <int64_t*>malloc(total_references * sizeof(int64_t))
    if not reference_alignment_counts or not reference_ids or not reference_lengths:
        if reference_alignment_counts:
            free(reference_alignment_counts)
        if reference_ids:
            free(reference_ids)
        if reference_lengths:
            free(reference_lengths)
        if tsv_map != NULL:
            free_tsv_reference_map(tsv_map)
            tsv_map = NULL
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise MemoryError("Failed to allocate reference alignment arrays")

    with nogil:
        for i in range(total_references):
            reference_lengths[i] = sam_hdr_tid2len(bam_header, i)
            if reference_lengths[i] <= 0:
                reference_lengths[i] = 1000
    # NOTE: We'll query reference names from bam_header when writing per-read output
    # to avoid 166k+ Python list.append() calls in Phase 0

    if tsv_map != NULL:
        reference_lengths_from_tsv = True
        # Fast approach: iterate TSV entries and lookup in BAM header (not the other way around!)
        # This is 100x faster than looping through 166k BAM refs and doing TSV hash lookups
        with nogil:
            nentries = get_tsv_reference_count(tsv_map)
            for i in range(nentries):
                ref_name = tsv_map.entries[i].reference_name
                tsv_length = tsv_map.entries[i].reference_length
                tid32 = sam_hdr_name2tid(bam_header, ref_name)
                if 0 <= tid32 < total_references and tsv_length > 0:
                    reference_lengths[tid32] = tsv_length
        free_tsv_reference_map(tsv_map)
        tsv_map = NULL

    bam_idx_start = bf_monotonic_seconds()
    if verbose:
        bf_logging.log(LOG_TAG, "Loading BAM index...")
    bam_index = sam_index_load(bam_handle, bam_file_path)
    bam_idx_end = bf_monotonic_seconds()
    if verbose:
        bf_logging.log(LOG_TAG, "BAM index loaded in %.2f seconds", bam_idx_end - bam_idx_start)
    if not bam_index:
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise RuntimeError("Failed to load BAM index")

    _log_phase_duration("Phase 0", phase0_start)

    # ------------------------------------------------------------------
    # Phase 1: Accession Mapping
    # ------------------------------------------------------------------
    cdef double phase1_start = bf_monotonic_seconds()
    
    # Extract reference accessions for filtered accession map loading
    _announce_stage("Accession Mapping", "Loading taxonomy mappings for BAM references")
    
    if verbose:
        bf_logging.log(LOG_TAG, "Extracting reference accessions for filtered loading...")
    
    reference_accessions = []
    for i in range(total_references):
        ref_name = sam_hdr_tid2name(bam_header, i)
        if ref_name != NULL:
            ref_name_str = ref_name.decode('utf-8')
            if len(ref_name_str) > 0:
                # Extract first token (accession) from reference name
                tokens = ref_name_str.split()
                if len(tokens) > 0:
                    reference_accessions.append(tokens[0])
    
    if verbose:
        bf_logging.log(LOG_TAG, "Loading accession map for %d references...", len(reference_accessions))

    # Load accession map with filtering
    cdef AccessionMapping acc_map_obj = taxonomy_py.load_accession_map_from_file(
        accession_map_path_str,
        accession_filter=reference_accessions
    )

    if acc_map_obj is None:
        raise RuntimeError("Failed to load accession mapping")

    cdef TaxonomyDB* taxdb_c = tax_db_obj.db

    if taxdb_c == NULL or acc_map_obj.amap == NULL:
        raise RuntimeError("Failed to initialise taxonomy structures")

    cdef dict stats_map = None
    if stats_path is not None and stats_path != "":
        if verbose:
            bf_logging.log(LOG_TAG, "Loading LCA stats file: %s", stats_path)
        stats_map = _load_stats_file(stats_path)

    cdef uint64_t mapped = 0
    cdef uint64_t unmapped = 0
    for i in range(total_references):
        hts_idx_get_stat(bam_index, i, &mapped, &unmapped)
        reference_alignment_counts[i] = mapped
        reference_ids[i] = i

    batch_starts = <int64_t*>malloc(total_references * sizeof(int64_t))
    batch_ends = <int64_t*>malloc(total_references * sizeof(int64_t))
    if not batch_starts or not batch_ends:
        if batch_starts:
            free(batch_starts)
        if batch_ends:
            free(batch_ends)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        hts_idx_destroy(bam_index)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise MemoryError("Failed to allocate batch offsets")

    batch_count = create_balanced_batches_greedy(
        reference_ids,
        reference_alignment_counts,
        total_references,
        total_references,
        batch_starts,
        batch_ends,
        threads,
        verbose
    )

    if batch_count <= 0:
        free(batch_starts)
        free(batch_ends)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        hts_idx_destroy(bam_index)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise RuntimeError("Failed to create processing batches")

    lca_batches = <LCABatch**>malloc(batch_count * sizeof(LCABatch*))
    if not lca_batches:
        free(batch_starts)
        free(batch_ends)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        hts_idx_destroy(bam_index)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        raise MemoryError("Failed to allocate batch pointer array")

    cdef int64_t expected = 0
    cdef int64_t ref_idx
    for batch_idx in range(batch_count):
        expected = 0
        for ref_idx in range(batch_starts[batch_idx], batch_ends[batch_idx]):
            expected += reference_alignment_counts[reference_ids[ref_idx]]
        lca_batches[batch_idx] = create_lca_batch(
            batch_idx,
            batch_starts[batch_idx],
            batch_ends[batch_idx],
            expected
        )
        if not lca_batches[batch_idx]:
            free(batch_starts)
            free(batch_ends)
            free(reference_alignment_counts)
            free(reference_ids)
            free(reference_lengths)
            hts_idx_destroy(bam_index)
            sam_hdr_destroy(bam_header)
            hts_close(bam_handle)
            for i in range(batch_idx):
                destroy_lca_batch(lca_batches[i])
            free(lca_batches)
            raise MemoryError("Failed to create LCA batch")

    free(batch_starts)
    free(batch_ends)

    cdef samFile** thread_handles = <samFile**>malloc(threads * sizeof(samFile*))
    if not thread_handles:
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        hts_idx_destroy(bam_index)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        for i in range(batch_count):
            destroy_lca_batch(lca_batches[i])
        free(lca_batches)
        raise MemoryError("Failed to allocate thread BAM handles")

    for thread_idx in range(threads):
        thread_handles[thread_idx] = NULL

    for thread_idx in range(threads):
        thread_handles[thread_idx] = hts_open(bam_file_path, b"r")
        if not thread_handles[thread_idx]:
            for i in range(thread_idx):
                if thread_handles[i]:
                    hts_close(thread_handles[i])
            free(thread_handles)
            free(reference_alignment_counts)
            free(reference_ids)
            free(reference_lengths)
            hts_idx_destroy(bam_index)
            sam_hdr_destroy(bam_header)
            hts_close(bam_handle)
            for i in range(batch_count):
                destroy_lca_batch(lca_batches[i])
            free(lca_batches)
            raise RuntimeError(f"Failed to open thread BAM handle {thread_idx}")

        if threads > 1:
            hts_set_threads(thread_handles[thread_idx], 1)

    thread_maps = <ThreadLocalHashMap**>malloc(threads * sizeof(ThreadLocalHashMap*))
    if not thread_maps:
        for thread_idx in range(threads):
            if thread_handles[thread_idx]:
                hts_close(thread_handles[thread_idx])
        free(thread_handles)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        hts_idx_destroy(bam_index)
        sam_hdr_destroy(bam_header)
        hts_close(bam_handle)
        for i in range(batch_count):
            destroy_lca_batch(lca_batches[i])
        free(lca_batches)
        raise MemoryError("Failed to allocate thread maps")

    for thread_idx in range(threads):
        thread_maps[thread_idx] = create_thread_local_hash_map()
        if not thread_maps[thread_idx]:
            for i in range(thread_idx):
                destroy_thread_local_hash_map(thread_maps[i])
            for i in range(threads):
                if thread_handles[i]:
                    hts_close(thread_handles[i])
            free(thread_maps)
            free(thread_handles)
            free(reference_alignment_counts)
            free(reference_ids)
            free(reference_lengths)
            hts_idx_destroy(bam_index)
            sam_hdr_destroy(bam_header)
            hts_close(bam_handle)
            for i in range(batch_count):
                destroy_lca_batch(lca_batches[i])
            free(lca_batches)
            raise MemoryError("Failed to create thread hash map")

    _log_phase_duration("Phase 1", phase1_start)

    # ------------------------------------------------------------------
    # Phase 2: Batch Streaming
    # ------------------------------------------------------------------
    cdef double phase2_start = bf_monotonic_seconds()
    
    _announce_stage("Batch Streaming", "Processing BAM alignments in parallel batches with global ID assignment")

    # Process batches (serial or parallel depending on threads)
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Starting parallel batch streaming (%lld batches, %d threads)",
        batch_count,
        threads,
    )
    cdef time_t batch_start = time(NULL)
    with nogil:
        for batch_idx in prange(batch_count, num_threads=threads, schedule='static'):
            thread_idx = threadid()
            process_lca_batch_alignments(
                thread_handles[thread_idx],
                bam_header,
                bam_index,
                reference_ids,
                lca_batches[batch_idx],
                &scoring_config,
                thread_maps[thread_idx],
                thread_idx
            )
    cdef time_t batch_end = time(NULL)
    cdef double batch_elapsed = <double>(batch_end - batch_start)
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Batch streaming complete in %.2f seconds (%.1f batches/sec)",
        batch_elapsed,
        batch_count / (batch_elapsed if batch_elapsed > 0 else 1.0),
    )

    for thread_idx in range(threads):
        if thread_handles[thread_idx]:
            hts_close(thread_handles[thread_idx])
            thread_handles[thread_idx] = NULL
    free(thread_handles)

    # htslib handle used for header/index can be closed now
    hts_close(bam_handle)
    hts_idx_destroy(bam_index)

    # Assign global sequential IDs
    bf_nogil_logf_notime(LOG_TAG_B, "Assigning global sequential read IDs")
    cdef time_t id_assign_start = time(NULL)
    with nogil:
        assign_lca_global_ids_fast(lca_batches, batch_count, thread_maps, threads)
    cdef time_t id_assign_end = time(NULL)
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Global ID assignment complete in %.2f seconds",
        <double>(id_assign_end - id_assign_start),
    )

    # Ensure per-batch alignments are grouped by read_index for correct streaming
    bf_nogil_logf_notime(LOG_TAG_B, "Sorting alignments by read index within each batch")
    cdef time_t sort_start = time(NULL)
    with nogil:
        for batch_idx in prange(batch_count, num_threads=threads, schedule='static'):
            if lca_batches[batch_idx] != NULL and lca_batches[batch_idx].actual_alignment_count > 1:
                radix_sort_lca_by_read(
                    lca_batches[batch_idx].alignments,
                    lca_batches[batch_idx].actual_alignment_count
                )
    cdef time_t sort_end = time(NULL)
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Batch sorting complete in %.2f seconds",
        <double>(sort_end - sort_start),
    )

    cdef int64_t total_alignments = 0
    cdef int64_t unique_reads = 0
    bf_nogil_logf_notime(LOG_TAG_B, "Counting alignments and unique reads")
    cdef time_t count_start = time(NULL)
    with nogil:
        total_alignments = count_lca_alignments(lca_batches, batch_count)
        unique_reads = count_unique_reads_from_thread_maps(thread_maps, threads)
    cdef time_t count_end = time(NULL)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Found %lld total alignments for %lld unique reads in %.2f seconds",
        total_alignments,
        unique_reads,
        <double>(count_end - count_start),
    )
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Average alignments per read: %.1f",
        <double>total_alignments / <double>unique_reads if unique_reads > 0 else 0.0,
    )

    if total_alignments == 0 or unique_reads == 0:
        for thread_idx in range(threads):
            destroy_thread_local_hash_map(thread_maps[thread_idx])
        free(thread_maps)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        for batch_idx in range(batch_count):
            destroy_lca_batch(lca_batches[batch_idx])
        free(lca_batches)
        raise RuntimeError("No valid alignments found in BAM")

    _log_phase_duration("Phase 2", phase2_start)

    # Stream through batches directly without building giant read_entries array
    # This avoids the ~17GB allocation that was causing OOM on large datasets
    cdef int64_t n_refs = total_references
    cdef int32_t* reference_taxids = NULL
    cdef int32_t taxid
    cdef double* reference_weight_sum = NULL
    cdef double* reference_read_sum = NULL
    cdef double* reference_tad_counts = NULL
    # Temp vars for accession token extraction
    cdef const char* ref_c
    cdef int acc_len
    cdef char* acc_key
    cdef int allocation_error = 0
    
    # If we plan per-read output, build read_id -> read_name mapping BEFORE destroying thread maps
    if write_per_read:
        # Build merged unique hash list
        global_map_names = kh_init_seqid_map()
        capacity = unique_reads if unique_reads > 0 else 1024
        all_hashes = <uint64_t*>malloc(capacity * sizeof(uint64_t))
        total_hashes = 0
        if global_map_names == NULL or all_hashes == NULL:
            if global_map_names != NULL: kh_destroy_seqid_map(global_map_names)
            if all_hashes != NULL: free(all_hashes)
            # Fall back: disable per-read to avoid crash
            write_per_read = False
        else:
            for thread_idx in range(threads):
                if thread_maps[thread_idx] == NULL:
                    continue
                k_t = 0
                while k_t < kh_end_seqid_map(thread_maps[thread_idx].hash_to_id_map):
                    if kh_exist_seqid_map(thread_maps[thread_idx].hash_to_id_map, k_t):
                        hv = kh_key_seqid_map(thread_maps[thread_idx].hash_to_id_map, k_t)
                        k_g = kh_get_seqid_map(global_map_names, hv)
                        if k_g == kh_end_seqid_map(global_map_names):
                            k_g = kh_put_seqid_map(global_map_names, hv, &ret_g)
                            if ret_g != -1:
                                if total_hashes >= capacity:
                                    capacity *= 2
                                    tmp_hashes = <uint64_t*>realloc(all_hashes, capacity * sizeof(uint64_t))
                                    if tmp_hashes == NULL:
                                        # disable per-read safely
                                        kh_destroy_seqid_map(global_map_names)
                                        free(all_hashes)
                                        global_map_names = NULL
                                        all_hashes = NULL
                                        write_per_read = False
                                        break
                                    all_hashes = tmp_hashes
                                all_hashes[total_hashes] = hv
                                total_hashes += 1
                    k_t += 1
                if not write_per_read:
                    break
            if write_per_read:
                if total_hashes > 1:
                    radix_sort_uint64(all_hashes, total_hashes)
                # Rebuild global map: hash -> global id
                kh_destroy_seqid_map(global_map_names)
                global_map_names = kh_init_seqid_map()
                if global_map_names == NULL:
                    free(all_hashes)
                    all_hashes = NULL
                    write_per_read = False
                else:
                    for i in range(total_hashes):
                        hv = all_hashes[i]
                        k_g = kh_put_seqid_map(global_map_names, hv, &ret_g)
                        if ret_g != -1:
                            kh_val_seqid_map_wrap(global_map_names, k_g)[0] = <uint32_t>i
                    # Build Python list of names indexed by global id
                    read_names = [None] * total_hashes
                    for thread_idx in range(threads):
                        if thread_maps[thread_idx] == NULL:
                            continue
                        k_t = 0
                        while k_t < kh_end_seqid_name_map(thread_maps[thread_idx].hash_to_name_map):
                            if kh_exist_seqid_name_map(thread_maps[thread_idx].hash_to_name_map, k_t):
                                # We need the hash key corresponding to this name entry. Iterate id_map for keys
                                # Use id_map to get key list; name_map shares the same keys set
                                # Fallback: iterate id_map to retrieve names via name_map lookup
                                pass
                            k_t += 1
                    # Since khash doesn't expose iteration over keys for name map directly here,
                    # iterate over id_map and fetch names from name_map with the same key
                    for thread_idx in range(threads):
                        if thread_maps[thread_idx] == NULL:
                            continue
                        k_t = 0
                        while k_t < kh_end_seqid_map(thread_maps[thread_idx].hash_to_id_map):
                            if kh_exist_seqid_map(thread_maps[thread_idx].hash_to_id_map, k_t):
                                hv = kh_key_seqid_map(thread_maps[thread_idx].hash_to_id_map, k_t)
                                k_g = kh_get_seqid_map(global_map_names, hv)
                                if k_g != kh_end_seqid_map(global_map_names):
                                    gidx = kh_val_seqid_map_wrap(global_map_names, k_g)[0]
                                    k_n = kh_get_seqid_name_map(thread_maps[thread_idx].hash_to_name_map, hv)
                                    if k_n != kh_end_seqid_name_map(thread_maps[thread_idx].hash_to_name_map):
                                        name_ptr = kh_val_seqid_name_map_wrap(thread_maps[thread_idx].hash_to_name_map, k_n)[0]
                                        if name_ptr != NULL and read_names[gidx] is None:
                                            read_names[gidx] = PyUnicode_FromString(name_ptr)
                            k_t += 1
                    kh_destroy_seqid_map(global_map_names)
                    free(all_hashes)

    # Now release thread maps
    for thread_idx in range(threads):
        destroy_thread_local_hash_map(thread_maps[thread_idx])
    free(thread_maps)
    
    if verbose:
        bf_logging.log(LOG_TAG, "Streaming LCA computation directly from %d batches (avoiding large memory allocation)", batch_count)
    
    # Precompute reference -> taxid mapping
    reference_taxids = <int32_t*>malloc(n_refs * sizeof(int32_t))
    if reference_taxids == NULL:
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        for batch_idx in range(batch_count):
            destroy_lca_batch(lca_batches[batch_idx])
        free(lca_batches)
        raise MemoryError("Failed to allocate reference taxid array")
    
    for i in range(n_refs):
        ref_name = sam_hdr_tid2name(bam_header, i)
        if ref_name == NULL:
            reference_taxids[i] = -1
            continue
        # Normalize reference name to accession token (first whitespace-delimited field)
        # so lookups match keys in accession_map.parquet (which are pure accessions).
        ref_c = ref_name
        acc_len = 0
        while ref_c[acc_len] != 0 and ref_c[acc_len] != 32 and ref_c[acc_len] != 9:
            acc_len += 1
        acc_key = <char*>malloc(acc_len + 1)
        if acc_key != NULL:
            memcpy(acc_key, ref_c, acc_len)
            acc_key[acc_len] = 0
            taxid = acc_map_obj._get_taxid_nogil(acc_key)
            free(acc_key)
        else:
            # Fallback to using the full reference name if allocation fails
            taxid = acc_map_obj._get_taxid_nogil(ref_name)
        reference_taxids[i] = taxid
    
    # Allocate reference weight/read buffers
    reference_weight_sum = <double*>calloc(n_refs, sizeof(double))
    reference_read_sum = <double*>calloc(n_refs, sizeof(double))
    reference_tad_counts = <double*>calloc(n_refs, sizeof(double))
    
    if (reference_weight_sum == NULL) or (reference_read_sum == NULL) or (reference_tad_counts == NULL):
        if reference_weight_sum != NULL:
            free(reference_weight_sum)
        if reference_read_sum != NULL:
            free(reference_read_sum)
        if reference_tad_counts != NULL:
            free(reference_tad_counts)
        free(reference_taxids)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        for batch_idx in range(batch_count):
            destroy_lca_batch(lca_batches[batch_idx])
        free(lca_batches)
        raise MemoryError("Failed to allocate reference weight buffers")
    
    # Initialize TAD counts - optimized to avoid slow Python dict lookups per reference
    cdef object stats_tuple
    cdef double val
    cdef bytes ref_name_bytes
    cdef dict stats_map_bytes
    cdef time_t start_time, end_time
    cdef double elapsed
    
    start_time = time(NULL)
    bf_logging.info("[%s] Initializing TAD counts...", LOG_TAG)
    if verbose:
        bf_logging.log(LOG_TAG, "Initializing TAD counts from stats file...")
    
    if stats_map is not None:
        # Pre-build a simple dict keyed by bytes to avoid repeated string conversions
        stats_map_bytes = {}
        for key, value in stats_map.items():
            if isinstance(key, str):
                stats_map_bytes[key.encode('utf-8')] = value
            else:
                stats_map_bytes[key] = value
        
        # Now populate TAD counts with fast byte-key lookups
        for i in range(n_refs):
            reference_tad_counts[i] = 0.0
            ref_name = sam_hdr_tid2name(bam_header, i)
            if ref_name != NULL:
                ref_name_bytes = ref_name  # C string auto-converts to bytes
                stats_tuple = stats_map_bytes.get(ref_name_bytes, None)
                if stats_tuple is not None:
                    try:
                        val = float(stats_tuple[1]) if len(stats_tuple) > 1 and stats_tuple[1] is not None else 0.0
                        if val > 0.0:
                            reference_tad_counts[i] = val
                    except Exception:
                        pass
    else:
        # No stats file: just zero-initialize
        for i in range(n_refs):
            reference_tad_counts[i] = 0.0

    end_time = time(NULL)
    elapsed = <double>(end_time - start_time)
    bf_logging.info("[%s] TAD counts initialized (%.2f seconds)", LOG_TAG, elapsed)
    if verbose:
        bf_logging.log(LOG_TAG, "TAD initialization details: completed")

    # Prepare storage for taxonomy summary
    cdef int32_t rank_id = get_rank_id(rank)
    if rank_id < 0:
        free(reference_weight_sum)
        free(reference_read_sum)
        free(reference_tad_counts)
        free(reference_taxids)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        raise ValueError(f"Unknown taxonomy rank '{rank}'")

    cdef int32_t n_nodes = taxdb_c.n_nodes
    cdef double* read_counts = <double*>calloc(n_nodes, sizeof(double))
    cdef double* weight_sums = <double*>calloc(n_nodes, sizeof(double))
    cdef double* length_sums = <double*>calloc(n_nodes, sizeof(double))
    # Aggregated read counts propagated to ancestors for hierarchical reporting
    cdef double* agg_read_counts = <double*>calloc(n_nodes, sizeof(double))
    if not read_counts or not weight_sums or not length_sums or not agg_read_counts:
        if read_counts:
            free(read_counts)
        if weight_sums:
            free(weight_sums)
        if length_sums:
            free(length_sums)
        if agg_read_counts:
            free(agg_read_counts)
        free(reference_taxids)
        free(reference_alignment_counts)
        free(reference_ids)
        free(reference_lengths)
        sam_hdr_destroy(bam_header)
        raise MemoryError("Failed to allocate taxonomy summary buffers")

    # ------------------------------------------------------------------
    # Parallel batch streaming - assign global IDs and collect alignments
    # (Phase 2 already announced and timed above)
    # ------------------------------------------------------------------
    cdef double* ref_weight_tls = NULL
    cdef double* ref_read_tls = NULL
    cdef uint64_t* accepted_tls = NULL
    cdef uint64_t* discarded_tls = NULL
    
    # Pad thread-local arrays to cache-line boundaries to avoid false sharing
    # Each thread gets n_refs + padding doubles to ensure cache-line alignment
    cdef int64_t refs_per_thread_padded = ((n_refs + CACHE_LINE_DOUBLES - 1) // CACHE_LINE_DOUBLES) * CACHE_LINE_DOUBLES

    start_time = time(NULL)
    ref_weight_tls = <double*>calloc(threads * refs_per_thread_padded, sizeof(double))
    ref_read_tls = <double*>calloc(threads * refs_per_thread_padded, sizeof(double))
    accepted_tls = <uint64_t*>calloc(threads * CACHE_LINE_DOUBLES, sizeof(uint64_t))  # pad scalars too
    discarded_tls = <uint64_t*>calloc(threads * CACHE_LINE_DOUBLES, sizeof(uint64_t))
    if (ref_weight_tls == NULL) or (ref_read_tls == NULL) or (accepted_tls == NULL) or (discarded_tls == NULL):
        if ref_weight_tls != NULL: free(ref_weight_tls)
        if ref_read_tls != NULL: free(ref_read_tls)
        if accepted_tls != NULL: free(accepted_tls)
        if discarded_tls != NULL: free(discarded_tls)
        free(read_counts); free(weight_sums); free(length_sums); free(agg_read_counts)
        free(reference_taxids); free(reference_weight_sum); free(reference_read_sum); free(reference_tad_counts)
        free(reference_alignment_counts); free(reference_ids); free(reference_lengths)
        sam_hdr_destroy(bam_header)
        for batch_idx in range(batch_count):
            destroy_lca_batch(lca_batches[batch_idx])
        free(lca_batches)
        raise MemoryError("Failed to allocate thread-local reduction buffers")

    # ------------------------------------------------------------------
    # Phase 3: Memory Pool
    # ------------------------------------------------------------------
    cdef double phase3_start = bf_monotonic_seconds()
    
    # ------------------------------------------------------------------
    # Memory pool creation - allocate global sorted pool like processor.py
    # ------------------------------------------------------------------
    _announce_stage("Memory Pool", "Creating global sorted pool for efficient per-read LCA computation")
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Creating global LCA pool for %lld alignments, %lld unique reads",
        total_alignments,
        unique_reads,
    )
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Pool memory: ~%.2f GB for alignments + ~%.2f MB for index",
        (total_alignments * 12) / 1e9,
        (unique_reads * 8) / 1e6,
    )
    
    # Create LCA pool (single global sorted array like processor.py MemoryPool)
    cdef LCAPool* lca_pool = NULL
    cdef time_t pool_create_start = time(NULL)
    with nogil:
        lca_pool = create_lca_pool(total_alignments, unique_reads)
    cdef time_t pool_create_end = time(NULL)
    cdef double pool_create_elapsed = <double>(pool_create_end - pool_create_start)
    
    if lca_pool == NULL:
        free(read_counts); free(weight_sums); free(length_sums); free(agg_read_counts)
        free(reference_taxids); free(reference_weight_sum); free(reference_read_sum); free(reference_tad_counts)
        free(reference_alignment_counts); free(reference_ids); free(reference_lengths)
        sam_hdr_destroy(bam_header)
        for batch_idx in range(batch_count):
            destroy_lca_batch(lca_batches[batch_idx])
        free(lca_batches)
        free(ref_weight_tls); free(ref_read_tls); free(accepted_tls); free(discarded_tls)
        raise MemoryError("Failed to create LCA pool")
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Pool created in %.2f seconds",
        pool_create_elapsed,
    )
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Copying %lld batches to pool",
        batch_count,
    )
    
    # Populate pool: copy batches, destroy them, radix-sort, build read index
    cdef time_t populate_start = time(NULL)
    cdef int populate_result
    with nogil:
        populate_result = populate_lca_pool(lca_pool, lca_batches, batch_count)
    cdef time_t populate_end = time(NULL)
    cdef double populate_elapsed = <double>(populate_end - populate_start)
    
    if populate_result != 0:
        destroy_lca_pool(lca_pool)
        free(read_counts); free(weight_sums); free(length_sums); free(agg_read_counts)
        free(reference_taxids); free(reference_weight_sum); free(reference_read_sum); free(reference_tad_counts)
        free(reference_alignment_counts); free(reference_ids); free(reference_lengths)
        sam_hdr_destroy(bam_header)
        free(lca_batches)
        free(ref_weight_tls); free(ref_read_tls); free(accepted_tls); free(discarded_tls)
        raise RuntimeError(f"Failed to populate LCA pool - alignment count exceeded capacity (truncation detected)")
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Pool populated in %.2f seconds (copied, sorted, indexed)",
        populate_elapsed,
    )
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Pool contains %lld alignments for %lld unique reads",
        lca_pool.alignment_count,
        lca_pool.unique_read_count,
    )
    
    free(lca_batches)
    lca_batches = NULL
    
    _log_phase_duration("Phase 3", phase3_start)

    # ------------------------------------------------------------------
    # Phase 4: Parallel LCA
    # ------------------------------------------------------------------
    cdef double phase4_start = bf_monotonic_seconds()
    
    # Per-read output setup (if requested)
    cdef PerReadLCA** per_read_results_tls = NULL
    cdef int64_t max_tid_idx
    
    if write_per_read:
        # Allocate per-thread result buffers
        per_read_results_tls = <PerReadLCA**>calloc(threads, sizeof(PerReadLCA*))
        if per_read_results_tls == NULL:
            bf_logging.error("Failed to allocate per-read results buffers")
            write_per_read = False
        else:
            for thread_idx in range(threads):
                per_read_results_tls[thread_idx] = <PerReadLCA*>calloc(unique_reads, sizeof(PerReadLCA))
                if per_read_results_tls[thread_idx] == NULL:
                    bf_logging.error("Failed to allocate per-read buffer for thread %d", thread_idx)
                    # Clean up allocated buffers
                    for cleanup_idx in range(thread_idx):
                        free(per_read_results_tls[cleanup_idx])
                    free(per_read_results_tls)
                    per_read_results_tls = NULL
                    write_per_read = False
                    break
    
    # ------------------------------------------------------------------
    # Parallel LCA processing - per-read dedup, LCA, mode-specific weighting
    # ------------------------------------------------------------------
    _announce_stage("Parallel LCA", "Computing taxonomic assignments with strict path intersection across %d threads" % threads)
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Starting parallel per-read LCA processing (%lld reads, %d threads, strict mode)",
        unique_reads,
        threads,
    )
    
    # Process per-read LCA from sorted pool in parallel
    cdef time_t lca_start = time(NULL)
    with nogil:
        process_lca_from_pool_parallel(
            lca_pool,
            reference_taxids,
            rank_id,
            taxdb_c,
            ref_weight_tls,
            ref_read_tls,
            accepted_tls,
            discarded_tls,
            refs_per_thread_padded,
            threads,
            per_read_results_tls
        )
    cdef time_t lca_end = time(NULL)
    cdef double lca_elapsed = <double>(lca_end - lca_start)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Parallel LCA processing complete in %.2f seconds (%.1f reads/sec)",
        lca_elapsed,
        unique_reads / (lca_elapsed if lca_elapsed > 0 else 1.0),
    )
    
    # Free pool immediately after processing
    cdef time_t pool_destroy_start = time(NULL)
    with nogil:
        destroy_lca_pool(lca_pool)
    cdef time_t pool_destroy_end = time(NULL)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Pool destroyed in %.2f seconds",
        <double>(pool_destroy_end - pool_destroy_start),
    )
    
    # ------------------------------------------------------------------
    # Thread-local reduction - merge TLS arrays from all threads
    # ------------------------------------------------------------------
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Reducing thread-local results from %d threads",
        threads,
    )
    cdef time_t reduce_start = time(NULL)

    # Reduce thread-local reference aggregates
    cdef int64_t ri
    cdef int32_t ti
    cdef uint64_t discarded_reads = 0
    cdef uint64_t accepted_reads = 0
    for ri in range(n_refs):
        reference_weight_sum[ri] = 0.0
        reference_read_sum[ri] = 0.0
    for ti in range(threads):
        for ri in range(n_refs):
            reference_weight_sum[ri] += ref_weight_tls[ti * refs_per_thread_padded + ri]
            reference_read_sum[ri] += ref_read_tls[ti * refs_per_thread_padded + ri]
        accepted_reads += accepted_tls[ti * CACHE_LINE_DOUBLES]
        discarded_reads += discarded_tls[ti * CACHE_LINE_DOUBLES]

    free(ref_weight_tls)
    free(ref_read_tls)
    free(accepted_tls)
    free(discarded_tls)
    
    cdef time_t reduce_end = time(NULL)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Reduction complete in %.2f seconds",
        <double>(reduce_end - reduce_start),
    )
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Read acceptance summary: %lld accepted (%.1f%%), %lld discarded (%.1f%%)",
        accepted_reads,
        100.0 * accepted_reads / (accepted_reads + discarded_reads) if (accepted_reads + discarded_reads) > 0 else 0.0,
        discarded_reads,
        100.0 * discarded_reads / (accepted_reads + discarded_reads) if (accepted_reads + discarded_reads) > 0 else 0.0,
    )
    
    # ------------------------------------------------------------------
    # Per-read output writing (if requested)
    # ------------------------------------------------------------------
    cdef const char* norm_ref_cstr = NULL
    if write_per_read and per_read_results_tls != NULL:
        if verbose:
            bf_logging.log(LOG_TAG, "Writing per-read LCA assignments to %s", per_read_path)

        # Merge per-thread results into a single list (Python side for simplicity)
        # Output ALL reads, including those without valid LCA (lca_taxid <= 0)
        per_read_list = []
        for thread_idx in range(threads):
            if per_read_results_tls[thread_idx] != NULL:
                for read_id in range(unique_reads):
                    # Output ALL reads, not just those with valid LCA
                    if per_read_results_tls[thread_idx][read_id].read_id == read_id:
                        per_read_list.append((
                            per_read_results_tls[thread_idx][read_id].read_id,
                            per_read_results_tls[thread_idx][read_id].lca_taxid,
                            per_read_results_tls[thread_idx][read_id].num_alignments,
                            per_read_results_tls[thread_idx][read_id].norm_ref_index,
                            per_read_results_tls[thread_idx][read_id].is_trusted
                        ))
        
        # Write to file (with gz support)
        use_gz = per_read_path.endswith('.gz')
        if use_gz:
            import gzip
            with gzip.open(per_read_path, 'wt') as f:
                f.write("read_name\tlca_taxid\tlca_rank\tn_aln\ttax_path\tnorm_ref\tnorm_ref_len\ttrusted\n")
                for read_id, lca_taxid, num_alignments, norm_ref_idx, is_trusted in per_read_list:
                    read_name = read_names[read_id] if read_names is not None else str(read_id)
                    rank_name = get_rank_name(tax_db_obj._get_rank_id_nogil(lca_taxid)) if lca_taxid > 0 else "unassigned"
                    if lca_taxid > 0 and build_lineage_string_nogil(taxdb_c, lca_taxid, lineage_buf_perread, 4096) >= 0:
                        lineage = PyUnicode_FromString(lineage_buf_perread)
                    else:
                        lineage = "unassigned"
                    # Query BAM header directly to avoid 166k+ Python list building in Phase 0
                    norm_ref_cstr = sam_hdr_tid2name(bam_header, norm_ref_idx) if norm_ref_idx >= 0 else NULL
                    norm_ref = norm_ref_cstr.decode('utf-8') if norm_ref_cstr != NULL else str(norm_ref_idx)
                    norm_ref_len = reference_lengths[norm_ref_idx] if norm_ref_idx >= 0 and norm_ref_idx < total_references else 0
                    f.write(f"{read_name}\t{lca_taxid}\t{rank_name}\t{num_alignments}\t{lineage}\t{norm_ref}\t{norm_ref_len}\t{is_trusted}\n")
        else:
            with open(per_read_path, 'w') as f:
                f.write("read_name\tlca_taxid\tlca_rank\tn_aln\ttax_path\tnorm_ref\tnorm_ref_len\ttrusted\n")
                for read_id, lca_taxid, num_alignments, norm_ref_idx2, is_trusted in per_read_list:
                    read_name = read_names[read_id] if read_names is not None else str(read_id)
                    rank_name = get_rank_name(tax_db_obj._get_rank_id_nogil(lca_taxid)) if lca_taxid > 0 else "unassigned"
                    if lca_taxid > 0 and build_lineage_string_nogil(taxdb_c, lca_taxid, lineage_buf_perread, 4096) >= 0:
                        lineage = PyUnicode_FromString(lineage_buf_perread)
                    else:
                        lineage = "unassigned"
                    # Query BAM header directly to avoid 166k+ Python list building in Phase 0
                    norm_ref_cstr = sam_hdr_tid2name(bam_header, norm_ref_idx2) if norm_ref_idx2 >= 0 else NULL
                    norm_ref2 = norm_ref_cstr.decode('utf-8') if norm_ref_cstr != NULL else str(norm_ref_idx2)
                    norm_ref_len2 = reference_lengths[norm_ref_idx2] if norm_ref_idx2 >= 0 and norm_ref_idx2 < total_references else 0
                    f.write(f"{read_name}\t{lca_taxid}\t{rank_name}\t{num_alignments}\t{lineage}\t{norm_ref2}\t{norm_ref_len2}\t{is_trusted}\n")
        
        if verbose:
            bf_logging.log(LOG_TAG, "Wrote %d per-read LCA assignments", len(per_read_list))
        
        # Free per-read buffers
        for thread_idx in range(threads):
            if per_read_results_tls[thread_idx] != NULL:
                free(per_read_results_tls[thread_idx])
        free(per_read_results_tls)

    _log_phase_duration("Phase 4", phase4_start)

    # ------------------------------------------------------------------
    # Phase 5: Taxonomy Aggregation
    # ------------------------------------------------------------------
    cdef double phase5_start = bf_monotonic_seconds()
    
    # ------------------------------------------------------------------
    # Taxonomy aggregation - propagate TAD-normalized weights up tree
    # ------------------------------------------------------------------
    _announce_stage("Taxonomy Aggregation", "Propagating abundances through taxonomy hierarchy with TAD normalization")
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Aggregating taxonomy assignments across %d references",
        n_refs,
    )
    cdef time_t agg_start = time(NULL)
    if verbose:
        bf_logging.log(LOG_TAG, "Applying TAD normalization and taxonomy aggregation")

    # Apply TAD-based normalisation and aggregate into taxonomy hierarchy
    for ref_idx in range(n_refs):
        if reference_weight_sum[ref_idx] <= 0.0 and reference_read_sum[ref_idx] <= 0.0:
            continue
        # Use TAD counts for more robust abundance estimates (trims outliers)
        # Fall back to raw counts if TAD not available
        if reference_tad_counts[ref_idx] > 0.0:
            weight = reference_tad_counts[ref_idx]
        else:
            weight = reference_weight_sum[ref_idx]
        lca_taxid = reference_taxids[ref_idx]
        if lca_taxid < 0:
            continue
        _propagate_up(taxdb_c, lca_taxid, reference_read_sum[ref_idx], weight, read_counts, weight_sums)

    cdef time_t agg_end = time(NULL)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Taxonomy aggregation complete in %.2f seconds",
        <double>(agg_end - agg_start),
    )
    
    # ------------------------------------------------------------------
    # Length aggregation - propagate reference lengths up taxonomy tree
    # ------------------------------------------------------------------
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Aggregating reference lengths across %d references",
        n_refs,
    )
    cdef time_t length_start = time(NULL)
    if verbose:
        bf_logging.log(LOG_TAG, "Aggregating reference lengths")

    # Aggregate total reference lengths per taxonomy node
    cdef double length_val
    for ref_idx in range(n_refs):
        lca_taxid = reference_taxids[ref_idx]
        if lca_taxid < 0:
            continue
        length_val = <double>reference_lengths[ref_idx]
        if length_val <= 0.0:
            continue
        _propagate_length(taxdb_c, lca_taxid, length_val, length_sums)
    
    cdef time_t length_end = time(NULL)
    
    bf_nogil_logf_notime(
        LOG_TAG_B,
        "Length aggregation complete in %.2f seconds",
        <double>(length_end - length_start),
    )

    _log_phase_duration("Phase 5", phase5_start)

    end_time = time(NULL)
    elapsed = <double>(end_time - start_time)
    bf_logging.info("[%s] Preparing output file (%.2f seconds)...", LOG_TAG, elapsed)
    if verbose:
        bf_logging.log(LOG_TAG, "Copying read counts for output")

    # Read counts are already at LCA level only - just copy them for output
    # Do NOT aggregate up the tree (that would inflate counts massively)
    cdef int32_t idx_node
    cdef int32_t parent_idx
    for idx_node in range(n_nodes):
        agg_read_counts[idx_node] = read_counts[idx_node]

    end_time = time(NULL)
    elapsed = <double>(end_time - start_time)
    bf_logging.info("[%s] Writing results (%.2f seconds)...", LOG_TAG, elapsed)
    if verbose:
        bf_logging.log(LOG_TAG, "Cleaning up intermediate structures")

    # Cleanup intermediate BAM structures (batches already freed earlier to reduce peak memory)
    free(reference_alignment_counts)
    free(reference_ids)
    free(reference_lengths)
    free(reference_weight_sum)
    free(reference_read_sum)
    free(reference_tad_counts)

    sam_hdr_destroy(bam_header)

    if verbose:
        bf_logging.log(LOG_TAG, "Cleanup complete, building output adjacency lists...")

    # Log filtering statistics
    if verbose:
        bf_logging.log(
            LOG_TAG, 
            "LCA rank filtering: %d reads accepted, %d reads discarded (below rank '%s')",
            accepted_reads,
            discarded_reads,
            rank
        )

    # ------------------------------------------------------------------
    # Write output file (pure C I/O, nogil)
    # ------------------------------------------------------------------
    if verbose:
        bf_logging.log(LOG_TAG, "Writing LCA summary to %s", output_path)

    cdef double total_weight_root = 0.0
    if taxdb_c != NULL and taxdb_c.root_idx >= 0 and taxdb_c.root_idx < n_nodes:
        total_weight_root = weight_sums[taxdb_c.root_idx]

    # ------------------------------------------------------------------
    # Ultra-fast cache-friendly file writing: hierarchical tree traversal with lexicographic sorting
    # ------------------------------------------------------------------
    start_time = time(NULL)
    bf_logging.info("[%s] Writing LCA summary to %s", LOG_TAG, output_path)
    
    # Pre-compute all formatted lines AND build parent-child adjacency for hierarchical output
    cdef int32_t i32
    cdef const char* node_name_ptr
    cdef const char* rank_ptr  
    cdef char lineage_buf_write[4096]
    cdef dict node_lines = {}  # idx -> (name, line_str) for fast lookup
    cdef dict parent_children = {}  # parent_idx -> list of child indices
    cdef TaxNode write_node
    cdef list children_list
    cdef list sort_pairs
    
    if verbose:
        bf_logging.log(LOG_TAG, "Pre-computing output lines for %d active nodes...", n_nodes)
    
    # Build formatted strings and parent-child relationships
    for i32 in range(n_nodes):
        # Write nodes that have either read counts OR weight/abundance (includes ancestors)
        if agg_read_counts[i32] <= 0.0 and weight_sums[i32] <= 0.0:
            continue
        
        write_node = taxdb_c.nodes[i32]
        node_name_ptr = get_name_at_rank_nogil(taxdb_c, write_node.taxid, write_node.rank_id)
        if node_name_ptr != NULL:
            name_str = PyUnicode_FromString(node_name_ptr)
        else:
            name_str = ""
        
        if taxdb_c.rank_names != NULL:
            rank_ptr = taxdb_c.rank_names[write_node.rank_id]
            rank_str = PyUnicode_FromString(rank_ptr) if rank_ptr != NULL else ""
        else:
            rank_str = ""
        
        if build_lineage_string_nogil(taxdb_c, write_node.taxid, lineage_buf_write, 4096) < 0:
            lineage_buf_write[0] = 0
        lineage_str = PyUnicode_FromString(lineage_buf_write)
        
        # Format line and cache it
        line_str = f"{write_node.taxid}\t{name_str}\t{rank_str}\t{agg_read_counts[i32]:.0f}\t{weight_sums[i32]:.6f}\t{lineage_str}\n"
        node_lines[i32] = (name_str, line_str)
        
        # Build parent-child relationships
        if write_node.parent_taxid != write_node.taxid:
            parent_idx = _taxid_to_idx(taxdb_c, write_node.parent_taxid)
            if parent_idx >= 0:
                if parent_idx not in parent_children:
                    parent_children[parent_idx] = []
                parent_children[parent_idx].append(i32)
    
    if verbose:
        bf_logging.log(LOG_TAG, "Pre-computed %d lines, sorting children lexicographically...", len(node_lines))
    
    # Sort children lexicographically at each level using operator.itemgetter
    for parent_idx in parent_children:
        children_list = parent_children[parent_idx]
        # Build (name, idx) pairs for sorting, then extract indices
        sort_pairs = [(node_lines[idx][0], idx) for idx in children_list if idx in node_lines]
        sort_pairs.sort(key=operator.itemgetter(0))
        parent_children[parent_idx] = [idx for _, idx in sort_pairs]
    
    # Open file
    if output_path.endswith('.gz'):
        fp = gzip.open(output_path, 'wt', compresslevel=6)
    else:
        fp = open(output_path, 'w')
    
    # Write header
    fp.write("taxid\tname\trank\tn_reads\tabundance\ttax_path\n")
    
    # Iterative tree traversal using a stack (avoid nested function/recursion)
    cdef int wrote = 0
    cdef list stack = []
    cdef int current_idx
    cdef list child_indices
    
    # Start from root
    if taxdb_c.root_idx >= 0 and taxdb_c.root_idx in node_lines:
        stack.append(taxdb_c.root_idx)
    
    # Process stack (depth-first traversal with children in reverse order for correct output)
    while len(stack) > 0:
        current_idx = stack.pop()
        if current_idx not in node_lines:
            continue
        
        # Write this node
        _, line_str = node_lines[current_idx]
        fp.write(line_str)
        wrote += 1
        
        # Add children to stack in REVERSE order (stack pops from end, so reverse gives correct order)
        if current_idx in parent_children:
            child_indices = parent_children[current_idx]
            for i in range(len(child_indices) - 1, -1, -1):
                stack.append(child_indices[i])
    
    fp.close()
    
    end_time = time(NULL)
    elapsed = <double>(end_time - start_time)
    bf_logging.info("[%s] File writing complete - %d nodes written (%.2f seconds)", LOG_TAG, wrote, elapsed)
    if verbose:
        bf_logging.log(LOG_TAG, "Output file successfully created")

    # Build dictionary for tests / downstream consumption
    rows = {}
    cdef char lineage_buf_dict[4096]
    cdef double total_length
    cdef double coverage_per_bp
    cdef double scaled_weight_per_bp
    cdef double node_scaled_value
    cdef double total_reference_length = 0.0
    for i in range(n_nodes):
        # Include nodes with either read counts OR weight/abundance
        if agg_read_counts[i] <= 0.0 and weight_sums[i] <= 0.0:
            continue
        node = taxdb_c.nodes[i]
        node_name_ptr = get_name_at_rank_nogil(taxdb_c, node.taxid, node.rank_id)
        if node_name_ptr != NULL:
            name_str = PyUnicode_FromString(node_name_ptr)
        else:
            name_str = ""
        rank_str = get_rank_name(node.rank_id)
        if build_lineage_string_nogil(taxdb_c, node.taxid, lineage_buf_dict, 4096) < 0:
            lineage_buf_dict[0] = 0
        lineage_str = PyUnicode_FromString(lineage_buf_dict)
        total_length = length_sums[i]
        if total_length > 0.0:
            coverage_per_bp = weight_sums[i] / total_length
        else:
            coverage_per_bp = 0.0
        if total_weight_root > 0.0:
            node_scaled_value = (weight_sums[i] / total_weight_root) * scale_factor
        else:
            node_scaled_value = 0.0
        if total_length > 0.0:
            scaled_weight_per_bp = node_scaled_value / total_length
        else:
            scaled_weight_per_bp = 0.0
        total_reference_length += total_length
        rows[node.taxid] = {
            "name": name_str,
            "rank": rank_str,
            # Expose aggregated hierarchical read counts in rows for downstream
            "reads": agg_read_counts[i],
            "norm_weight": weight_sums[i],
            "scaled_weight": node_scaled_value,
            "path": lineage_str,
            "length_bp": total_length,
            "norm_weight_per_bp": coverage_per_bp,
            "scaled_weight_per_bp": scaled_weight_per_bp,
        }

    free(read_counts)
    free(weight_sums)
    free(length_sums)
    free(agg_read_counts)
    free(reference_taxids)

    metadata = {
        "unique_reads": unique_reads,
        "total_alignments": total_alignments,
        "reference_lengths_source": "tsv" if reference_lengths_from_tsv else "bam",
        "reference_lengths_path": reference_lengths_tsv if reference_lengths_from_tsv else None,
        "total_reference_length_bp": total_reference_length,
        "scale_factor": scale_factor,
        "total_weight_root": total_weight_root,
        "taxonomy_db_dir": taxonomy_db_dir_str,
        "accession_map_path": accession_map_path_str,
    }
    if tsv_map != NULL:
        free_tsv_reference_map(tsv_map)
    return _LCAStatistics(rows, metadata), tax_db_obj
