# cython: language_level=3
from libc.stdint cimport int32_t, int64_t, uint32_t

from bam_filter.processor cimport MemoryPool

cdef MemoryPool* create_memory_pool(int64_t alignment_count,
                                    uint32_t reference_count,
                                    uint32_t unique_read_count,
                                    int64_t* filtered_ref_lengths,
                                    bint enable_pmd,
                                    int32_t max_threads) except NULL nogil

cdef MemoryPool* create_memory_pool_split(int64_t alignment_count,
                                          uint32_t reference_count,
                                          uint32_t unique_read_count,
                                          int64_t* filtered_ref_lengths,
                                          bint enable_pmd,
                                          bint enable_hierarchical,
                                          bint enable_damage_counts,
                                          int32_t max_threads) except NULL nogil

cdef void destroy_memory_pool(MemoryPool* pool) noexcept nogil

cdef int shrink_memory_pool(MemoryPool* pool) except -1 nogil

cdef void cleanup_presorted_memory() noexcept nogil

cdef void cleanup_em_intermediate_memory(MemoryPool* pool) noexcept nogil
