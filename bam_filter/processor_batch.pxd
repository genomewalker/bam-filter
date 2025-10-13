# cython: language_level=3
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from libc.stddef cimport size_t

# Import shared types and HTSlib typedefs (use relative cimports so sibling pxds are found)
from .processor cimport Alignment, MemoryPool, AlignmentScoringConfig, INVALID_SEQUENTIAL_ID
from .processor_sort cimport radix_sort_uint64, radix_sort_alignments_by_read_id, radix_sort_compact_by_position
from .processor cimport hts_idx_t, hts_itr_t, samFile, sam_hdr_t
from .processor_types cimport CompactAlignment
from .processor_hash cimport ThreadLocalHashMap
from .processor_hash cimport (
    kh_seqid_map_t, 
    khint_t, 
    khint64_t, 
    kh_get_seqid_map, 
    kh_seqid_map_t, 
    kh_init_seqid_map,
    kh_end_seqid_map,
    kh_put_seqid_map,
    kh_exist_seqid_map,
    kh_destroy_seqid_map,
    kh_size,
    kh_key_seqid_map,
    kh_val_seqid_map_wrap
    )
from cython.parallel cimport prange

cdef struct BatchAlignment:
    uint32_t read_index
    uint32_t reference_index
    uint32_t alignment_position
    float    alignment_score
    float    pmd_score

cdef struct ProcessingBatch:
    int64_t batch_identifier
    int64_t reference_start_index
    int64_t reference_end_index
    int64_t expected_alignment_count
    int64_t actual_alignment_count
    BatchAlignment* batch_alignments
    int64_t batch_capacity
    int32_t error_status
    int32_t processed_by_thread_id
    bint uses_mmap
    size_t mmap_size

# Function prototypes (nogil where applicable)
cdef ProcessingBatch* create_processing_batch(int64_t batch_id, int64_t ref_start,
                                                                int64_t ref_end, int64_t expected_count) except NULL nogil 

cdef int grow_batch_capacity(ProcessingBatch* batch) except -1 nogil 

cdef void destroy_processing_batch(ProcessingBatch* batch) noexcept nogil

cdef int process_batch_alignments(samFile* bam_file, sam_hdr_t* header,
                                           hts_idx_t* index, int64_t* reference_ids,
                                           ProcessingBatch* batch,
                                           AlignmentScoringConfig* scoring_config,
                                           ThreadLocalHashMap* thread_map,
                                           int32_t thread_id) except -1 nogil 

cdef int assign_global_sequential_ids_fast(ProcessingBatch** batches, int64_t batch_count,
                                          ThreadLocalHashMap** thread_maps, int num_threads) except -1 nogil 

cdef int64_t count_actual_alignments(ProcessingBatch** batches, int64_t batch_count) except -1 nogil 

cdef int64_t count_unique_refs_from_batches(ProcessingBatch** batches, int64_t batch_count) except -1 nogil 

cdef void parallel_streaming_stats_optimized(ProcessingBatch** batches, int64_t batch_count,
                                            double* out_min, double* out_max, double* out_mean,
                                            double* out_variance, int64_t* out_count, int num_threads) noexcept nogil

cdef int populate_memory_pool_direct(MemoryPool* pool,
                                     ProcessingBatch** batches,
                                     int64_t batch_count,
                                     sam_hdr_t* header,
                                     int num_threads) except -1 nogil 


# Provide mmap/madvise prototypes if other modules cimport this pxd
cdef extern from "sys/mman.h" nogil:
    void* mmap(void* addr, size_t length, int prot, int flags, int fd, long offset) nogil
    int munmap(void* addr, size_t length) nogil
    int madvise(void* addr, size_t length, int advice) nogil
    int PROT_READ
    int PROT_WRITE
    int MAP_PRIVATE
    int MAP_ANONYMOUS
    int MADV_SEQUENTIAL

cdef int64_t count_unique_reads_from_thread_maps(ThreadLocalHashMap** thread_maps, int num_threads) noexcept nogil