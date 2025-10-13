# cython: language_level=3

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t
from bam_filter.processor cimport Alignment
from bam_filter.processor_bam_writer cimport CompactAlignment

cdef void swap_compact(CompactAlignment* a, CompactAlignment* b) noexcept nogil
cdef uint32_t extract_byte_32(uint32_t key, int shift_bits) noexcept nogil
cdef uint32_t extract_byte_64(uint64_t key, int shift_bits) noexcept nogil
cdef void radix_sort_uint64_range(uint64_t* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil
cdef void radix_sort_uint64(uint64_t* a, int64_t n) noexcept nogil

cdef void radix_sort_compact_range(CompactAlignment* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil
cdef void radix_sort_compact_by_position(CompactAlignment* a, int64_t n) noexcept nogil

cdef void swap_alignments(Alignment* a, Alignment* b) noexcept nogil
cdef void insertion_sort_alignments_by_read_id(Alignment* alignments, int64_t left, int64_t right) noexcept nogil
cdef uint32_t extract_byte_32(uint32_t key, int shift_bits) noexcept nogil
cdef void radix_sort_alignments_range(Alignment* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil
cdef void radix_sort_alignments_by_read_id(Alignment* alignments, int64_t count) noexcept nogil

cdef void heapify_alignments_range(Alignment* alignments, int64_t offset, int64_t n, int64_t i) noexcept nogil
cdef void heapsort_alignments_range(Alignment* alignments, int64_t left, int64_t right) noexcept nogil
cdef int64_t partition_alignments_inplace(Alignment* alignments, int64_t left, int64_t right) noexcept nogil
cdef void introsort_alignments(Alignment* alignments, int64_t left, int64_t right, int depth_limit, int num_threads) noexcept nogil
