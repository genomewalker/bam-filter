from libc.stdint cimport int32_t, int64_t

from bam_filter.stats_rle cimport rle_interval_t, rle_coverage_t
from bam_filter.stats_structs cimport ref_stats_t

# Expose small, independent helpers that operate on RLE arrays and
# do not depend on other modules. Implementations are in
# `bam_filter/stats_coverage.pyx` and this .pxd is the canonical place
# to declare their signatures so other modules can cimport them.
cdef double compute_coverage_evenness_from_rle(int64_t ref_length, int32_t* starts, int32_t* ends, int32_t* depths, int64_t n_intervals) noexcept nogil

cdef void compute_tad_from_rle(
    int64_t ref_length,
    int32_t* starts, int32_t* ends, int32_t* depths, int64_t n_intervals,
    double* result_mean, int64_t* result_len,
    int trim_min, int trim_max
) noexcept nogil

cdef double compute_normalized_entropy_inline(int32_t* counts, int64_t n_bins) noexcept nogil

cdef double compute_normalized_gini_inline(int32_t* counts, int64_t n_bins) noexcept nogil

# Low-level helpers for entropy and Gini (used by other modules)
cdef double compute_entropy(int32_t* counts, int64_t n_bins) except -1 nogil

cdef double compute_gini(int32_t* counts, int64_t n_bins) except -1 nogil

cdef void compute_abundance_metrics(ref_stats_t* stats, int64_t n_alns, int64_t n_reads, double read_length_mean, int64_t scale) noexcept nogil

cdef void compute_rle_coverage_stats(rle_coverage_t* rle, ref_stats_t* stats, int trim_min, int trim_max) noexcept nogil
