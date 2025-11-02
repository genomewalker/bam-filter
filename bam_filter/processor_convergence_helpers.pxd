# cython: language_level=3
# Cython declarations for convergence helpers
cdef double _median5(double* a) except -1.0 nogil
cdef int _count_filled(int iteration, int H) except -1 nogil
cdef void _mad_statistics(double* buf, int H, int filled, double* med_out, double* sigma_out) noexcept nogil
