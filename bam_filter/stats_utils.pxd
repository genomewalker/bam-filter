# Header file for statistical utility functions
# cython: language_level=3

from libc.stdint cimport int32_t, int64_t
from stats_structures cimport ReferenceStats

cdef double calculate_entropy(int32_t* counts, int64_t n_bins) nogil
cdef double calculate_normalized_entropy(int32_t* counts, int64_t n_bins) nogil
cdef double calculate_gini_coefficient(int32_t* values, int64_t n_values) nogil
cdef double calculate_normalized_gini(int32_t* values, int64_t n_values) nogil
cdef void calculate_covered_regions(int32_t* coverage, int64_t ref_length, ReferenceStats* stats) nogil
cdef double calculate_tad(int32_t* coverage, int64_t ref_length, 
                          int trim_min, int trim_max, int64_t* tad_length) nogil
cdef double calculate_coverage_evenness(int32_t* coverage, int64_t ref_length) nogil
