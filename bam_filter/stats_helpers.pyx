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
"""
Small stats helpers: comparison functions and a mode-from-sorted helper.

These are intentionally minimal and nogil-friendly so other Cython modules
can cimport and reuse them. We will progressively move more helpers here.
"""

from libc.stdint cimport int32_t, int64_t, uint8_t, uint32_t
from libc.stdlib cimport malloc, free
from libc.string cimport strlen, memset
from libc.math cimport sqrt

# minimal externs required by helpers
cdef extern from "htslib/sam.h":
    ctypedef struct bam1_t:
        pass
    uint8_t* bam_aux_get(bam1_t* b, const char* tag) nogil
    int bam_aux2i(const uint8_t *s) nogil
    char* bam_get_qname(bam1_t* b) nogil
    int32_t bam_endpos(bam1_t* b) nogil
    uint8_t* bam_get_seq(bam1_t* b) nogil
    int bam_seqi(const uint8_t* s, int i) nogil

# Lookup table for GC content calculation (A=0, C=1, G=1, T=0)
# Define it here so the C symbol is emitted once; other modules cimport
# the array from stats_helpers.pxd.
cdef int[16] GC_LOOKUP = [0, 0, 1, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]


cdef int compare_int32(const void* a, const void* b) noexcept nogil:
    cdef int32_t val_a = (<int32_t*>a)[0]
    cdef int32_t val_b = (<int32_t*>b)[0]
    if val_a < val_b:
        return -1
    elif val_a > val_b:
        return 1
    return 0


cdef int compare_double(const void* a, const void* b) noexcept nogil:
    cdef double val_a = (<double*>a)[0]
    cdef double val_b = (<double*>b)[0]
    if val_a < val_b:
        return -1
    elif val_a > val_b:
        return 1
    return 0


cdef int compare_pairs(const void* a, const void* b) noexcept nogil:
    cdef int64_t* pair_a = <int64_t*>a
    cdef int64_t* pair_b = <int64_t*>b
    if pair_a[0] < pair_b[0]:
        return -1
    elif pair_a[0] > pair_b[0]:
        return 1
    return 0


cdef int compare_int64(const void* a, const void* b) noexcept nogil:
    cdef int64_t val_a = (<int64_t*>a)[0]
    cdef int64_t val_b = (<int64_t*>b)[0]
    if val_a < val_b:
        return -1
    elif val_a > val_b:
        return 1
    return 0


cdef inline int mode_from_sorted(int32_t* sorted_vals, int64_t n) except -1 nogil:
    """Calculate mode from pre-sorted array - more efficient."""
    if n == 0:
        return 0

    cdef int mode = sorted_vals[0]
    cdef int max_count = 1
    cdef int current_count = 1
    cdef int64_t i

    for i in range(1, n):
        if sorted_vals[i] == sorted_vals[i - 1]:
            current_count += 1
        else:
            if current_count > max_count:
                max_count = current_count
                mode = sorted_vals[i - 1]
            current_count = 1

    if current_count > max_count:
        mode = sorted_vals[n - 1]

    return mode


cdef int32_t get_query_alignment_length(bam1_t *src) noexcept nogil:
    if src == NULL or src.core.n_cigar == 0:
        return src.core.l_qseq if src != NULL else 0
    cdef uint32_t *cigar_p = <uint32_t *>(src.data + src.core.l_qname)
    cdef uint32_t aligned_length = 0
    cdef int k
    cdef uint32_t op, op_len
    cdef int n_cigar = <int>src.core.n_cigar
    for k in range(n_cigar):
        op = cigar_p[k] & 0xf
        op_len = cigar_p[k] >> 4
        if op == 0 or op == 1 or op == 7 or op == 8:
            aligned_length += op_len
    return aligned_length


cdef int count_gc_bases(bam1_t* b) noexcept nogil:
    cdef int gc = 0, i, base1, base2
    cdef int l_qseq = b.core.l_qseq
    cdef uint8_t* seq = b.data + b.core.l_qname + (b.core.n_cigar * 4)
    cdef int pairs = l_qseq >> 1
    cdef uint8_t byte_val
    cdef int chunk_pairs = pairs & ~3
    for i in range(0, chunk_pairs, 4):
        byte_val = seq[i]
        gc += GC_LOOKUP[(byte_val >> 4) & 0xF] + GC_LOOKUP[byte_val & 0xF]
        byte_val = seq[i + 1]
        gc += GC_LOOKUP[(byte_val >> 4) & 0xF] + GC_LOOKUP[byte_val & 0xF]
        byte_val = seq[i + 2]
        gc += GC_LOOKUP[(byte_val >> 4) & 0xF] + GC_LOOKUP[byte_val & 0xF]
        byte_val = seq[i + 3]
        gc += GC_LOOKUP[(byte_val >> 4) & 0xF] + GC_LOOKUP[byte_val & 0xF]
    for i in range(chunk_pairs, pairs):
        byte_val = seq[i]
        base1 = (byte_val >> 4) & 0xF
        base2 = byte_val & 0xF
        gc += GC_LOOKUP[base1] + GC_LOOKUP[base2]
    if l_qseq & 1:
        byte_val = seq[pairs]
        base1 = (byte_val >> 4) & 0xF
        gc += GC_LOOKUP[base1]
    return gc


cdef int count_reference_gc_bases(bam1_t* b, int32_t* ref_length_out) noexcept nogil:
    """
    Count G and C bases in the REFERENCE sequence using MD tag.
    Efficient single-pass implementation without strlen.

    For damaged reads (C→T): reference has C, read has T.
    This function counts the original C in reference.
    """
    cdef int ref_gc = 0
    cdef int32_t ref_length = 0
    cdef int read_pos = 0
    cdef int i

    if b == NULL:
        ref_length_out[0] = 0
        return 0

    cdef uint8_t* seq = bam_get_seq(b)
    cdef int read_len = b.core.l_qseq

    cdef uint8_t* md_tag = bam_aux_get(b, b"MD")
    if md_tag == NULL:
        ref_length_out[0] = 0
        return 0

    # Verify MD tag is string type ('Z')
    if md_tag[0] != 90:
        ref_length_out[0] = 0
        return 0

    cdef char* md = <char*>(md_tag + 1)
    cdef char c
    cdef int num = 0
    cdef int base_code

    while True:
        c = md[0]
        if c == 0:
            break
        md += 1

        if c >= 48 and c <= 57:  # '0'-'9': matches
            num = c - 48
            while md[0] >= 48 and md[0] <= 57:
                num = num * 10 + (md[0] - 48)
                md += 1
            for i in range(num):
                if read_pos < read_len:
                    base_code = bam_seqi(seq, read_pos)
                    if base_code == 2 or base_code == 4:  # C or G
                        ref_gc += 1
                    ref_length += 1
                    read_pos += 1
        elif c == 94:  # '^': deletion in read
            while md[0] >= 65 and md[0] <= 90:  # A-Z
                if md[0] == 67 or md[0] == 71:  # 'C' or 'G'
                    ref_gc += 1
                ref_length += 1
                md += 1
        elif c >= 65 and c <= 90:  # A-Z: mismatch
            if c == 67 or c == 71:  # 'C' or 'G'
                ref_gc += 1
            ref_length += 1
            read_pos += 1

    ref_length_out[0] = ref_length
    return ref_gc


cdef float compute_ani(bam1_t* b) noexcept nogil:
    cdef uint8_t* aux = bam_aux_get(b, b"NM")
    cdef int nm = bam_aux2i(aux) if aux != NULL else -1
    cdef int l_qseq = b.core.l_qseq
    if nm < 0 or l_qseq == 0:
        return 0.0
    return (1.0 - (<float>nm / l_qseq)) * 100.0


cdef int extract_aux_int(uint8_t* aux) noexcept nogil:
    if aux == NULL:
        return -1
    cdef char tag_type = aux[0]
    cdef uint8_t* data = aux + 1
    if tag_type == 105:
        return (<int32_t*>data)[0]
    elif tag_type == 99:
        return (<int8_t*>data)[0]
    elif tag_type == 115:
        return (<int16_t*>data)[0]
    elif tag_type == 67:
        return (<uint8_t*>data)[0]
    elif tag_type == 83:
        return (<uint16_t*>data)[0]
    elif tag_type == 73:
        return <int>(<uint32_t*>data)[0]
    return -1


cdef int64_t fnv1a_hash_read_id(char* qname) noexcept nogil:
    cdef uint64_t hash_val = 14695981039346656037UL
    cdef uint64_t fnv_prime = 1099511628211UL
    cdef unsigned char c
    cdef int i = 0
    while i < 64 and qname[i] != 0:
        c = <unsigned char>qname[i]
        hash_val ^= c
        hash_val *= fnv_prime
        i += 1
    return <int64_t>hash_val


cdef inline int _base_to_idx(int b) noexcept nogil:
    if b == 1:
        return 0  # A
    elif b == 2:
        return 1  # C
    elif b == 4:
        return 2  # G
    elif b == 8:
        return 3  # T
    else:
        return -1

cdef double calculate_dust_score(bam1_t* b) noexcept nogil:
    """
    Compute a normalized DUST score (0-1) for a read.
    0 = high complexity, 1 = extremely repetitive.
    """
    if b == NULL or b.core.l_qseq < 3:
        return 0.0

    cdef int triplets_capacity = 64  # 4^3 combinations
    cdef int counts[64]
    cdef int i
    for i in range(64):
        counts[i] = 0

    cdef uint8_t* seq = bam_get_seq(b)
    cdef int base0, base1, base2
    cdef int trip_key
    cdef int valid_triplets = 0
    cdef int read_len = b.core.l_qseq

    for i in range(read_len - 2):
        base0 = _base_to_idx(bam_seqi(seq, i))
        base1 = _base_to_idx(bam_seqi(seq, i + 1))
        base2 = _base_to_idx(bam_seqi(seq, i + 2))
        if base0 < 0 or base1 < 0 or base2 < 0:
            continue  # skip ambiguous bases
        trip_key = base0 * 16 + base1 * 4 + base2
        counts[trip_key] += 1
        valid_triplets += 1

    if valid_triplets <= 1:
        return 0.0

    cdef double score = 0.0
    cdef double c
    for i in range(64):
        if counts[i] > 1:
            c = <double>counts[i]
            score += (c * (c - 1.0)) / 2.0

    cdef double denom = (<double>valid_triplets * (valid_triplets - 1.0)) / 2.0
    if denom <= 0.0:
        return 0.0
    return score / denom
