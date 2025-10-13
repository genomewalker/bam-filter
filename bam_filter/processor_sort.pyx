# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
# -*- coding: utf-8 -*-

"""High-performance sorting algorithms for alignment data.

Provides American Flag radix sort and introsort implementations optimized
for sorting alignments by read ID or position with minimal memory overhead.
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t
from bam_filter.processor cimport Alignment
from cython.parallel cimport prange
from bam_filter.processor_bam_writer cimport CompactAlignment

cdef inline void swap_compact(CompactAlignment* a, CompactAlignment* b) noexcept nogil:
    """Swap two CompactAlignment structures in-place.

    Parameters
    ----------
    a : CompactAlignment*
        First structure
    b : CompactAlignment*
        Second structure
    """
    cdef CompactAlignment tmp = a[0]
    a[0] = b[0]
    b[0] = tmp

cdef inline uint32_t extract_byte_32(uint32_t key, int shift_bits) noexcept nogil:
    """Extract byte from 32-bit key at specified bit shift.

    Used by radix sort to extract individual bytes for bucket sorting.

    Parameters
    ----------
    key : uint32_t
        Integer key value
    shift_bits : int
        Bit shift amount (0, 8, 16, or 24)

    Returns
    -------
    uint32_t
        Extracted byte value (0-255)
    """
    return <uint32_t>((key >> shift_bits) & 0xFF)

cdef inline uint32_t extract_byte_64(uint64_t key, int shift_bits) noexcept nogil:
    """Extract byte from 64-bit key at specified bit shift.

    Used by radix sort to extract individual bytes for bucket sorting.

    Parameters
    ----------
    key : uint64_t
        Integer key value
    shift_bits : int
        Bit shift amount (0, 8, 16, 24, 32, 40, 48, or 56)

    Returns
    -------
    uint32_t
        Extracted byte value (0-255)
    """
    return <uint32_t>((key >> shift_bits) & 0xFF)

cdef void radix_sort_uint64_range(uint64_t* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil:
    """American Flag radix sort for uint64 array range.

    Recursive in-place radix sort that processes one byte at a time from
    most significant to least significant. Switches to insertion sort for
    small ranges (≤32 elements) or when all bytes are processed.

    American Flag algorithm uses in-place permutation rather than allocating
    auxiliary arrays, making it cache-efficient and memory-friendly.

    Parameters
    ----------
    a : uint64_t*
        Array to sort (modified in-place)
    lo : int64_t
        Start index (inclusive)
    hi : int64_t
        End index (exclusive)
    shift_bits : int
        Current byte position to sort by (56, 48, 40, ..., 0)
    """
    cdef int64_t n = hi - lo
    cdef int64_t count[256]
    cdef int64_t starts[256]
    cdef int64_t ends[256]
    cdef int64_t i, s, e
    cdef uint32_t b
    cdef int bucket
    cdef int64_t j
    cdef uint64_t key
    cdef uint64_t tmp
    if n <= 32 or shift_bits < 0:
        for i in range(lo+1, hi):
            key = a[i]; j = i - 1
            while j >= lo and a[j] > key:
                a[j+1] = a[j]; j -= 1
            a[j+1] = key
        return
    for bucket in range(256): count[bucket] = 0
    for i in range(lo, hi): count[extract_byte_64(a[i], shift_bits)] += 1
    starts[0] = lo
    for bucket in range(1,256): starts[bucket] = starts[bucket-1] + count[bucket-1]
    for bucket in range(256): ends[bucket] = starts[bucket] + count[bucket]
    bucket = 0
    while bucket < 256:
        s = starts[bucket]; e = ends[bucket]
        while s < e:
            b = extract_byte_64(a[s], shift_bits)
            if b == bucket:
                s += 1; starts[bucket] = s
            else:
                tmp = a[s]
                a[s] = a[starts[b]]; a[starts[b]] = tmp
                starts[b] += 1
        bucket += 1
    if shift_bits > 0:
        for bucket in range(256):
            s = ends[bucket] - count[bucket]; e = ends[bucket]
            if e - s > 1: radix_sort_uint64_range(a, s, e, shift_bits - 8)

cdef void radix_sort_uint64(uint64_t* a, int64_t n) noexcept nogil:
    """Sort uint64 array using American Flag radix sort.

    Parameters
    ----------
    a : uint64_t*
        Array to sort
    n : int64_t
        Array length
    """
    if n > 1: radix_sort_uint64_range(a, 0, n, 56)

cdef void radix_sort_compact_range(CompactAlignment* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil:
    """American Flag radix sort for CompactAlignment array range by position.

    Sorts CompactAlignment structures by their original_position field using
    in-place radix sort. Processes one byte at a time, switching to insertion
    sort for small ranges.

    Parameters
    ----------
    a : CompactAlignment*
        Array to sort (modified in-place)
    lo : int64_t
        Start index (inclusive)
    hi : int64_t
        End index (exclusive)
    shift_bits : int
        Current byte position to sort by (24, 16, 8, or 0)
    """
    cdef int64_t n = hi - lo
    cdef int64_t count[256]
    cdef int64_t starts[256]
    cdef int64_t ends[256]
    cdef int64_t i, s, e
    cdef uint32_t b
    cdef int bucket
    cdef int64_t j
    cdef CompactAlignment tmp
    cdef uint32_t key
    if n <= 32 or shift_bits < 0:
        for i in range(lo+1, hi):
            tmp = a[i]; key = tmp.original_position; j = i - 1
            while j >= lo and a[j].original_position > key:
                a[j+1] = a[j]; j -= 1
            a[j+1] = tmp
        return
    for bucket in range(256): count[bucket] = 0
    for i in range(lo, hi): count[extract_byte_32(a[i].original_position, shift_bits)] += 1
    starts[0] = lo
    for bucket in range(1,256): starts[bucket] = starts[bucket-1] + count[bucket-1]
    for bucket in range(256): ends[bucket] = starts[bucket] + count[bucket]
    bucket = 0
    while bucket < 256:
        s = starts[bucket]; e = ends[bucket]
        while s < e:
            b = extract_byte_32(a[s].original_position, shift_bits)
            if b == bucket:
                s += 1; starts[bucket] = s
            else:
                swap_compact(&a[s], &a[starts[b]])
                starts[b] += 1
        bucket += 1
    if shift_bits > 0:
        for bucket in range(256):
            s = ends[bucket] - count[bucket]; e = ends[bucket]
            if e - s > 1: radix_sort_compact_range(a, s, e, shift_bits - 8)

cdef void radix_sort_compact_by_position(CompactAlignment* a, int64_t n) noexcept nogil:
    """Sort CompactAlignment array by original_position using radix sort.

    Parameters
    ----------
    a : CompactAlignment*
        Array to sort
    n : int64_t
        Array length
    """
    if n > 1: radix_sort_compact_range(a, 0, n, 24)

cdef inline void swap_alignments(Alignment* a, Alignment* b) noexcept nogil:
    """Swap two Alignment structures in-place.

    Parameters
    ----------
    a : Alignment*
        First structure
    b : Alignment*
        Second structure
    """
    cdef Alignment temp = a[0]
    a[0] = b[0]
    b[0] = temp

cdef void insertion_sort_alignments_by_read_id(Alignment* alignments, int64_t left, int64_t right) noexcept nogil:
    """Insertion sort for small alignment ranges by read_index.

    Efficient for small arrays (typically ≤32 elements) where the low overhead
    of insertion sort outperforms more complex algorithms. Used as a fallback
    by radix sort and introsort.

    Parameters
    ----------
    alignments : Alignment*
        Array to sort (modified in-place)
    left : int64_t
        Start index (inclusive)
    right : int64_t
        End index (inclusive, note the difference from other sort functions)
    """
    cdef int64_t i, j
    cdef Alignment key
    cdef uint32_t key_read_index
    for i in range(left+1, right+1):
        key = alignments[i]
        key_read_index = key.read_index
        j = i - 1
        while j >= left and alignments[j].read_index > key_read_index:
            alignments[j+1] = alignments[j]
            j -= 1
        alignments[j+1] = key

cdef void radix_sort_alignments_range(Alignment* a, int64_t lo, int64_t hi, int shift_bits) noexcept nogil:
    """American Flag radix sort for Alignment array range by read_index.

    Core recursive function for in-place radix sorting of alignments. Processes
    one byte of the read_index field at a time, using bucket permutation to
    achieve in-place sorting with O(n) time complexity per byte.

    The American Flag algorithm computes bucket boundaries, then permutes elements
    in-place by cycling through buckets rather than using auxiliary arrays.

    Parameters
    ----------
    a : Alignment*
        Array to sort (modified in-place)
    lo : int64_t
        Start index (inclusive)
    hi : int64_t
        End index (exclusive)
    shift_bits : int
        Current byte position to sort by (24, 16, 8, or 0)
    """
    cdef int64_t count[256]
    cdef int64_t starts[256]
    cdef int64_t ends[256]
    cdef int64_t i, n = hi - lo
    cdef int bucket
    cdef int64_t s, e
    cdef uint32_t b

    if n <= 32:
        insertion_sort_alignments_by_read_id(a, lo, hi - 1)
        return
    if shift_bits < 0:
        return

    for bucket in range(256):
        count[bucket] = 0

    for i in range(lo, hi):
        count[extract_byte_32(a[i].read_index, shift_bits)] += 1

    starts[0] = lo
    for bucket in range(1, 256):
        starts[bucket] = starts[bucket - 1] + count[bucket - 1]

    for bucket in range(256):
        ends[bucket] = starts[bucket] + count[bucket]

    bucket = 0
    while bucket < 256:
        s = starts[bucket]
        e = ends[bucket]
        while s < e:
            b = extract_byte_32(a[s].read_index, shift_bits)
            if b == bucket:
                s += 1
                starts[bucket] = s
            else:
                swap_alignments(&a[s], &a[starts[b]])
                starts[b] += 1
        bucket += 1

    if shift_bits > 0:
        for bucket in range(256):
            s = ends[bucket] - count[bucket]
            e = ends[bucket]
            if e - s > 1:
                radix_sort_alignments_range(a, s, e, shift_bits - 8)

cdef void radix_sort_alignments_by_read_id(Alignment* alignments, int64_t count) noexcept nogil:
    """Sort Alignment array by read_index using American Flag radix sort.

    Parameters
    ----------
    alignments : Alignment*
        Array to sort in-place
    count : int64_t
        Array length
    """
    if count <= 1:
        return
    radix_sort_alignments_range(alignments, 0, count, 24)


cdef void heapify_alignments_range(Alignment* alignments, int64_t offset, int64_t n, int64_t i) noexcept nogil:
    """Heapify subtree for heapsort (used by introsort fallback).

    Maintains the max-heap property for a binary heap represented as an array.
    Called recursively to restore heap structure after swaps.

    Parameters
    ----------
    alignments : Alignment*
        Array containing the heap
    offset : int64_t
        Starting position of the heap in the array
    n : int64_t
        Size of the heap
    i : int64_t
        Index of root of subtree to heapify (relative to offset)
    """
    cdef int64_t largest = i
    cdef int64_t left_child = 2 * i + 1
    cdef int64_t right_child = 2 * i + 2
    if left_child < n and alignments[offset + left_child].read_index > alignments[offset + largest].read_index:
        largest = left_child
    if right_child < n and alignments[offset + right_child].read_index > alignments[offset + largest].read_index:
        largest = right_child
    if largest != i:
        swap_alignments(&alignments[offset + i], &alignments[offset + largest])
        heapify_alignments_range(alignments, offset, n, largest)

cdef void heapsort_alignments_range(Alignment* alignments, int64_t left, int64_t right) noexcept nogil:
    """Heapsort for alignment range (introsort fallback for deep recursion).

    O(n log n) worst-case sorting algorithm used as a fallback when quicksort
    recursion depth exceeds safe limits. Guarantees completion without risk
    of stack overflow.

    Parameters
    ----------
    alignments : Alignment*
        Array to sort (modified in-place)
    left : int64_t
        Start index (inclusive)
    right : int64_t
        End index (inclusive)
    """
    cdef int64_t n = right - left + 1
    cdef int64_t i
    for i in range(n // 2 - 1, -1, -1):
        heapify_alignments_range(alignments, left, n, i)
    for i in range(n - 1, 0, -1):
        swap_alignments(&alignments[left], &alignments[left + i])
        heapify_alignments_range(alignments, left, i, 0)

cdef int64_t partition_alignments_inplace(Alignment* alignments, int64_t left, int64_t right) noexcept nogil:
    """Partition alignments for quicksort (median-of-three pivot selection).

    Reorders array such that elements ≤ pivot are on the left and elements
    > pivot are on the right. Uses median-of-three pivot selection to avoid
    worst-case performance on partially sorted data.

    Parameters
    ----------
    alignments : Alignment*
        Array to partition (modified in-place)
    left : int64_t
        Start index (inclusive)
    right : int64_t
        End index (inclusive)

    Returns
    -------
    int64_t
        Final position of pivot element
    """
    cdef int64_t mid = left + (right - left) // 2
    cdef uint32_t pivot
    cdef int64_t i, j
    if alignments[left].read_index > alignments[mid].read_index:
        swap_alignments(&alignments[left], &alignments[mid])
    if alignments[mid].read_index > alignments[right].read_index:
        swap_alignments(&alignments[mid], &alignments[right])
    if alignments[left].read_index > alignments[mid].read_index:
        swap_alignments(&alignments[left], &alignments[mid])
    pivot = alignments[mid].read_index
    swap_alignments(&alignments[mid], &alignments[right])
    i = left - 1
    for j in range(left, right):
        if alignments[j].read_index <= pivot:
            i += 1
            swap_alignments(&alignments[i], &alignments[j])
    swap_alignments(&alignments[i + 1], &alignments[right])
    return i + 1

cdef void introsort_alignments(Alignment* alignments, int64_t left, int64_t right, int depth_limit, int num_threads) noexcept nogil:
    """Introsort for alignments (quicksort with heapsort fallback).

    Hybrid sorting algorithm using quicksort with heapsort fallback for
    deep recursion. Switches to insertion sort for small ranges. Supports
    optional parallelization for very large partitions.

    Parameters
    ----------
    alignments : Alignment*
        Array to sort by read_index
    left : int64_t
        Left boundary (inclusive)
    right : int64_t
        Right boundary (inclusive)
    depth_limit : int
        Recursion depth limit before switching to heapsort
    num_threads : int
        Number of threads for parallel partitioning (large ranges only)
    """
    cdef int64_t size = right - left + 1
    cdef int64_t pivot
    cdef int64_t left_size, right_size
    cdef int tid
    if size < 32:
        insertion_sort_alignments_by_read_id(alignments, left, right)
        return
    if depth_limit == 0:
        heapsort_alignments_range(alignments, left, right)
        return
    pivot = partition_alignments_inplace(alignments, left, right)
    left_size = pivot - left
    right_size = right - pivot
    if size > 100000 and num_threads > 1:
        for tid in prange(2, num_threads=2, schedule='static', nogil=True):
            if tid == 0 and left_size > 0:
                introsort_alignments(alignments, left, pivot - 1, depth_limit - 1, num_threads)
            elif tid == 1 and right_size > 0:
                introsort_alignments(alignments, pivot + 1, right, depth_limit - 1, num_threads)
    else:
        if left_size > 0:
            introsort_alignments(alignments, left, pivot - 1, depth_limit - 1, num_threads)
        if right_size > 0:
            introsort_alignments(alignments, pivot + 1, right, depth_limit - 1, num_threads)
