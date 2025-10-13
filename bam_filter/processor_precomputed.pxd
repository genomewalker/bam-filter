# cython: language_level=3
from libc.stdint cimport uint32_t
from bam_filter.processor_types cimport PrecomputedWeights

# Declare noexcept to avoid automatic exception checks when called without GIL
cdef PrecomputedWeights* create_precomputed_weights(uint32_t n_weights) noexcept nogil

cdef void update_precomputed_weights(PrecomputedWeights* pw, double* weights) noexcept nogil

cdef void free_precomputed_weights(PrecomputedWeights* pw) noexcept nogil
