# cython: language_level=3
from libc.stdint cimport uint8_t, int32_t, int64_t, uint64_t, uint16_t, uint32_t

# Import AlignmentScoringConfig from the main processor pxd so callers can pass the struct
from bam_filter.processor cimport AlignmentScoringConfig

# Use centralized htslib bindings from processor_types to avoid duplication
from bam_filter.processor_types cimport (
    bam1_core_t,
    bam1_t,
    sam_hdr_t,
    bam_endpos,
    bam_aux_get,
    bam_aux2i,
    bam_get_qname,
)

# Public C-visible functions provided by the MD/PMD quality module
cdef void initialize_quality_lookup_tables() noexcept nogil

cdef double calculate_md_quality_score(bam1_t* alignment,
                                       sam_hdr_t* header,
                                       AlignmentScoringConfig* config,
                                       float* pmd_result) noexcept nogil

cdef double calculate_md_score_fast_path(char* md_tag, uint8_t* qual_data, int32_t read_length) noexcept nogil

cdef double calculate_md_score_with_pmd_single_pass(char* md_tag, uint8_t* read_bases,
                                                    uint8_t* ref_sequence, uint8_t* qual_data,
                                                    int32_t read_length, bint is_single_stranded,
                                                    float* pmd_result) noexcept nogil

# Alignment filter helpers (implemented in the corresponding .pyx).
# Declared here so other modules can cimport them from
# bam_filter.processor_md_quality without Cython failing to locate
# per-symbol .pxd files during compilation.
cdef bint alignment_passes_quality_filters(bam1_t* alignment, AlignmentScoringConfig* config) noexcept nogil
cdef bint alignment_passes_quality_filters_with_ani(bam1_t* alignment, AlignmentScoringConfig* config, double* ani_out, int32_t* nm_out) noexcept nogil
