# cython: language_level=3
from libc.stdint cimport int64_t, int32_t
from .processor cimport ProcessingStats, MemoryPool

# Initialize statistics structure
cdef void init_processing_stats(ProcessingStats* stats) noexcept nogil

# Update statistics at each stage
cdef void update_initial_stats(ProcessingStats* stats, int64_t alignments, 
                               int64_t reads, int64_t references) noexcept nogil

cdef void update_quality_filter_stats(ProcessingStats* stats, int64_t alignments,
                                      int64_t reads, int64_t references) noexcept nogil

cdef void update_em_stats(ProcessingStats* stats, int32_t iterations,
                         bint converged, double likelihood) noexcept nogil

cdef void update_probability_filter_stats(ProcessingStats* stats, int64_t alignments,
                                          int64_t reads, int64_t references) noexcept nogil

cdef void update_graph_stats(ProcessingStats* stats, int64_t references,
                             int64_t patterns) noexcept nogil

cdef void update_unified_filter_stats(ProcessingStats* stats,
                                      int64_t alignments, int64_t reads, int64_t references,
                                      int64_t coverage_only, int64_t info_only, int64_t both,
                                      int64_t align_cov, int64_t align_info, int64_t align_both) noexcept nogil

cdef void update_final_output_stats(ProcessingStats* stats, int64_t alignments,
                                    int64_t reads, int64_t references) noexcept nogil

# Calculate summary metrics
cdef void calculate_summary_metrics(ProcessingStats* stats) noexcept nogil

# Print formatted statistics report
cdef void print_processing_stats(ProcessingStats* stats) noexcept nogil
