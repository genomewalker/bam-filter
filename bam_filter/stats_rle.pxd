# cython: language_level=3
from libc.stdint cimport int32_t, int64_t

from bam_filter.stats cimport RefStats, RLEInterval, RLECoverage

cdef RLECoverage* initialize_rle_from_length(int64_t ref_length) noexcept nogil

cdef void destroy_rle_coverage(RLECoverage* rle) noexcept nogil

cdef int add_coverage_interval(RLECoverage* rle, int64_t start, int64_t end, int32_t count) noexcept nogil

cdef void calculate_rle_coverage_stats(RLECoverage* rle, RefStats* stats, int trim_min, int trim_max) noexcept nogil

cdef void calculate_abundance_metrics(
	RefStats* stats,
	int64_t n_alns,
	int64_t n_reads,
	double read_length_mean,
	int64_t scale
) noexcept nogil
