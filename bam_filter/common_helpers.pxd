# cython: language_level=3
"""
Common inline helpers shared across Cython modules.

Place simple, nogil-safe inline functions here so other modules can
`cimport` them from a single canonical location.
"""
from libc.stdint cimport int32_t, int64_t, uint32_t
from libc.stddef cimport size_t

cdef extern from "unistd.h":
    int getpagesize() nogil

cdef inline uint32_t pack_position_length(uint32_t pos, uint32_t length) nogil:
    return (pos & 0xFFFFFF) | ((length & 0xFF) << 24)

cdef inline uint32_t extract_position(uint32_t packed) nogil:
    return packed & 0xFFFFFF

cdef inline uint32_t extract_length(uint32_t packed) nogil:
    return (packed >> 24) & 0xFF

cdef inline int64_t min_int64(int64_t a, int64_t b) noexcept nogil:
    return a if a < b else b

cdef inline int64_t max_int64(int64_t a, int64_t b) noexcept nogil:
    return a if a > b else b

cdef inline int32_t min_int32(int32_t a, int32_t b) noexcept nogil:
    return a if a < b else b

cdef inline int32_t max_int32(int32_t a, int32_t b) noexcept nogil:
    return a if a > b else b

cdef inline double min_double(double a, double b) noexcept nogil:
    return a if a < b else b

cdef inline double max_double(double a, double b) noexcept nogil:
    return a if a > b else b

cdef inline size_t page_size() nogil:
    cdef int ps = getpagesize()
    return <size_t>(ps if ps > 0 else 4096)
