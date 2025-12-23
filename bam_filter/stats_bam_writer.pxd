# cython: language_level=3
from libc.stdint cimport int32_t

from bam_filter.stats cimport RefStats
from bam_filter.generic_filters cimport GenericFilters
from bam_filter.processor_pmd cimport PMDCurve

# Use centralized htslib bindings
from bam_filter.processor_types cimport (
    hts_idx_t,
    hts_itr_t,
    bam1_t,
    samFile,
    sam_hdr_t,
)

cdef struct ReferenceFilter:
    int32_t* tid_mapping
    int32_t* reverse_mapping
    int32_t n_filtered_refs
    int32_t n_total_refs

cdef ReferenceFilter* create_reference_filter(RefStats* global_ref_stats,
                                              GenericFilters* gfilters,
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
    int max_read_length_c,
    PMDCurve* pmd_curve,
    float pmd_epsilon
) except -1 nogil
