# cython: language_level=3

from libc.stdint cimport int8_t, int16_t, int32_t, int64_t, uint16_t, uint8_t, uint32_t, uint64_t

# Lookup table for GC content calculation (A=0, C=1, G=1, T=0)
# The actual storage (definition) lives in stats_helpers.pyx; cimporters
# see the declaration here.
cdef int[16] GC_LOOKUP


# Declare minimal bam1_t type via extern from htslib header so we don't
# require a separate htslib Cython pxd file. The full definition exists in
# other modules that cdef extern from "htslib/sam.h"; here we only need an
# opaque type for pointer signatures.


# Use centralized htslib bindings from processor_types to avoid duplication
from bam_filter.processor_types cimport (
    bam1_core_t,
    bam1_t,
    bam_aux_get,
    bam_aux2i,
    bam_get_qname,
    bam_endpos,
)


cdef int32_t get_query_alignment_length(bam1_t *src) noexcept nogil
cdef int count_gc_bases(bam1_t* b) noexcept nogil
cdef int count_reference_gc_bases(bam1_t* b, int32_t* ref_length_out) noexcept nogil
cdef float compute_ani(bam1_t* b) noexcept nogil
cdef int extract_aux_int(uint8_t* aux) noexcept nogil
cdef int64_t fnv1a_hash_read_id(char* qname) noexcept nogil
cdef double calculate_dust_score(bam1_t* b) noexcept nogil

cdef int compare_int32(const void* a, const void* b) noexcept nogil
cdef int compare_double(const void* a, const void* b) noexcept nogil
cdef int compare_pairs(const void* a, const void* b) noexcept nogil
cdef int compare_int64(const void* a, const void* b) noexcept nogil

# Mode returns an int; declare exception value and place 'nogil' at end per Cython
cdef int mode_from_sorted(int32_t* sorted_vals, int64_t n) except -1 nogil
