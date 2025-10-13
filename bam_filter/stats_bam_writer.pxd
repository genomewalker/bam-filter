# cython: language_level=3
from libc.stdint cimport int32_t

from bam_filter.stats cimport RefStats, FilterConditions
from bam_filter.processor cimport samFile, sam_hdr_t

# Use centralized htslib bindings
from bam_filter.processor_types cimport (
    hts_idx_t,
    hts_itr_t,
    bam1_t,
)

cdef struct ReferenceFilter:
    int32_t* tid_mapping
    int32_t* reverse_mapping
    int32_t n_filtered_refs
    int32_t n_total_refs

cdef bint passes_filters(RefStats* stats, FilterConditions* filters) noexcept nogil

cdef ReferenceFilter* create_reference_filter(RefStats* global_ref_stats,
                                              FilterConditions* filters,
                                              int n_refs) noexcept nogil

cdef void destroy_reference_filter(ReferenceFilter* ref_filter) noexcept nogil

cdef int write_filtered_bam_streaming(
    samFile* input_bam,
    hts_idx_t* existing_index,
    const char* input_bam_path,
    const char* output_bam_path,
    sam_hdr_t* original_header,
    ReferenceFilter* ref_filter,
    int num_threads,
    double min_read_ani_c,
    int min_read_length_c,
    int max_read_length_c
) except -1 nogil
