# cython: language_level=3

cdef void safe_normalize_weights(double* weights, int dimension) noexcept nogil
# cython: language_level=3
cdef double stable_log_sum_exp(double log_a, double log_b) noexcept nogil
