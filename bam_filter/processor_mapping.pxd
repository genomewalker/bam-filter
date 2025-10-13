# cython: language_level=3
from libc.stdint cimport uint32_t, int64_t
from bam_filter.processor cimport sam_hdr_t, MemoryPool

cdef struct ReferenceMapping:
    uint32_t* old_to_new_tid
    uint32_t* new_to_old_tid
    uint32_t n_retained_refs
    uint32_t n_original_refs

cdef ReferenceMapping* create_reference_mapping(MemoryPool* pool, sam_hdr_t* original_header) noexcept nogil
cdef void destroy_reference_mapping(ReferenceMapping* mapping) noexcept nogil
cdef int update_reference_mapping_after_filtering(ReferenceMapping* mapping, MemoryPool* pool) noexcept nogil
cdef sam_hdr_t* create_filtered_header_efficient(sam_hdr_t* original_header, ReferenceMapping* mapping) noexcept nogil
cdef int remap_alignment_reference_ids(MemoryPool* pool, ReferenceMapping* mapping) noexcept nogil
