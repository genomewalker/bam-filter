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

"""MD tag parsing and alignment quality scoring.

This module parses the BAM `MD` tag (the read-wise mismatch/deletion string)
to compute alignment log-likelihood scores from base quality values.
It also integrates Post-Mortem Damage (PMD) likelihood components used for
ancient DNA damage-aware scoring (C->T and G->A transitions).

Notes
-----
- The identifier ``md`` in this module always refers to the BAM ``MD`` tag.
- Most functions operate ``nogil`` for performance; any edits must preserve
    the nogil/noexcept annotations where present.
"""

from libc.math cimport exp, log, fmax, fmin, pow as libc_pow
from libc.stdlib cimport malloc, free
from libc.stdint cimport uint8_t, int32_t

from bam_filter.processor cimport AlignmentScoringConfig

# HTSLib declarations (centralized in processor_types.pxd)
from .processor_types cimport (
    bam_get_seq, bam_get_qual, bam_aux_get, bam_aux2i,
    seq_nt16_str, bam_seqi_wrapper
)

cdef inline int parse_int(char** ptr_ref) nogil:
    """Fast integer parsing from character pointer.

    Parameters
    ----------
    ptr_ref : char**
        Pointer to character pointer (updated in-place as digits are consumed)

    Returns
    -------
    int
        Parsed integer value
    """
    cdef char* ptr = ptr_ref[0]
    cdef int num = 0
    if ptr[0] < 48 or ptr[0] > 57:
        return 0
    num = ptr[0] - 48
    ptr += 1
    if ptr[0] >= 48 and ptr[0] <= 57:
        num = num * 10 + (ptr[0] - 48)
        ptr += 1
        if ptr[0] >= 48 and ptr[0] <= 57:
            num = num * 10 + (ptr[0] - 48)
            ptr += 1
            while ptr[0] >= 48 and ptr[0] <= 57:
                num = num * 10 + (ptr[0] - 48)
                ptr += 1
    ptr_ref[0] = ptr
    return num


cdef double[94] PRECOMPUTED_LOG_P_CORRECT
cdef double[94] PRECOMPUTED_LOG_P_ERROR
cdef bint LOOKUP_TABLES_INITIALIZED = False

cdef double LOG_GAP_OPEN_PROB = -9.210340372
cdef double LOG_GAP_EXT_PROB  = -2.302585093
cdef double LOG_DELETION_PROB = -6.907755279
cdef double PHRED_TO_LN       = 0.23025850929940458

cdef double[94] EPSILON
cdef double[94] LOG_EPSILON
cdef double[94] LOG1M_EPSILON
cdef int SENTINEL_Q = 20

cdef double P_CONST = 0.3
cdef double C_CONST = 0.01
cdef double DECAY   = 0.7

DEF DMAX = 150
cdef double[DMAX] DZ
cdef double[DMAX] LOG_DZ
cdef double[DMAX] LOG1M_DZ

cdef unsigned char IS_A[256]
cdef unsigned char IS_C[256]
cdef unsigned char IS_G[256]
cdef unsigned char IS_T[256]
cdef unsigned char IS_UPPER[256]

cdef inline double lse2(double a, double b) nogil:
    """Numerically stable log-sum-exp for two values.

    Parameters
    ----------
    a : double
        First value in log-space.
    b : double
        Second value in log-space.

    Returns
    -------
    double
        log(exp(a) + exp(b)) computed in a numerically stable way.
    """
    cdef double m = a if a > b else b
    return m + log(exp(a - m) + exp(b - m))

cdef inline double lse3(double a, double b, double c) nogil:
    """Numerically stable log-sum-exp for three values.

    Parameters
    ----------
    a, b, c : double
        Values in log-space.

    Returns
    -------
    double
        log(exp(a) + exp(b) + exp(c)) computed stably.
    """
    cdef double m = a
    if b > m: m = b
    if c > m: m = c
    return m + log(exp(a - m) + exp(b - m) + exp(c - m))

cdef inline double lse4(double a, double b, double c, double d) nogil:
    """Numerically stable log-sum-exp for four values.

    Parameters
    ----------
    a, b, c, d : double
        Values in log-space.

    Returns
    -------
    double
        log(exp(a) + exp(b) + exp(c) + exp(d)) computed stably.
    """
    cdef double m = a
    if b > m: m = b
    if c > m: m = c
    if d > m: m = d
    return m + log(exp(a - m) + exp(b - m) + exp(c - m) + exp(d - m))

cdef inline void init_base_tables_once() noexcept nogil:
    """Initialize base lookup tables (A, C, G, T, uppercase) once.

    The tables are simple 256-entry lookup arrays mapping ASCII codes to
    boolean flags (0/1) for fast base-category checks in inner loops.
    Calling this repeatedly is safe because the function checks a sentinel
    byte in ``IS_UPPER`` before initializing.
    """
    if IS_UPPER[65]:
        return
    cdef int i
    for i in range(256):
        IS_A[i] = 1 if i == 65 else 0
        IS_C[i] = 1 if i == 67 else 0
        IS_G[i] = 1 if i == 71 else 0
        IS_T[i] = 1 if i == 84 else 0
        IS_UPPER[i] = 1 if (i >= 65 and i <= 90) else 0

cdef inline void get_D_terms(int z, double* logDz, double* log1mDz) noexcept nogil:
    """Get PMD damage probability terms for position z.

    Parameters
    ----------
    z : int
        Position from read end
    logDz : double*
        Output: log(D(z))
    log1mDz : double*
        Output: log(1 - D(z))
    """
    cdef double Dz
    if z <= DMAX:
        logDz[0]   = LOG_DZ[z-1]
        log1mDz[0] = LOG1M_DZ[z-1]
    else:
        Dz         = libc_pow(DECAY, z - 1) * P_CONST + C_CONST
        logDz[0]   = log(Dz)
        log1mDz[0] = log(1.0 - Dz)

cdef inline int decode_read_bases(uint8_t* seq_data, int32_t n, uint8_t* out) noexcept nogil:
    """Decode BAM sequence data to ASCII bases.

    Parameters
    ----------
    seq_data : uint8_t*
        Packed BAM sequence data (HTSlib internal format).
    n : int32_t
        Number of bases to decode.
    out : uint8_t*
        Output buffer (must be at least ``n`` bytes). Filled with ASCII
        characters (A/C/G/T) suitable for array-based comparisons.

    Returns
    -------
    int
        Zero on success.
    """
    cdef int i
    for i in range(n):
        out[i] = seq_nt16_str[bam_seqi_wrapper(seq_data, i)]
    return 0

cdef inline bint are_lookup_tables_initialized() nogil:
    """Return True if the lookup tables for quality/PMD have been initialized.

    This check is cheap and safe to call from nogil context.
    """
    return LOOKUP_TABLES_INITIALIZED

cdef void initialize_quality_lookup_tables() noexcept nogil:
    """Initialize all quality score and PMD lookup tables.

    Sets up:
    - Quality-to-error-rate conversion tables (EPSILON arrays)
    - PMD damage probability tables (DZ arrays)
    - Base character lookup tables (IS_A, IS_C, IS_G, IS_T)

    Thread-safe: Uses global flag to ensure single initialization.
    """
    global LOOKUP_TABLES_INITIALIZED
    cdef int qual, i
    cdef double p_error, p_correct
    cdef double powv

    if LOOKUP_TABLES_INITIALIZED:
        return

    init_base_tables_once()

    for qual in range(94):
        EPSILON[qual]       = 0.33333333 * exp(-qual * PHRED_TO_LN)
        LOG_EPSILON[qual]   = log(EPSILON[qual])
        LOG1M_EPSILON[qual] = log(1.0 - EPSILON[qual])

        p_error   = exp(-qual * PHRED_TO_LN)
        p_correct = 1.0 - p_error
        p_correct = fmax(fmin(p_correct, 0.9999), 0.0001)
        p_error   = fmax(p_error, 1e-10)
        PRECOMPUTED_LOG_P_CORRECT[qual] = log(p_correct)
        PRECOMPUTED_LOG_P_ERROR[qual]   = log(p_error)

    powv = 1.0
    for i in range(DMAX):
        if i == 0:
            powv = 1.0
        else:
            powv *= DECAY
        DZ[i]       = P_CONST * powv + C_CONST
        LOG_DZ[i]   = log(DZ[i])
        LOG1M_DZ[i] = log(1.0 - DZ[i])

    LOOKUP_TABLES_INITIALIZED = True

cdef inline int get_optimal_buffer_strategy(int32_t read_length) nogil:
    """Determine optimal buffer allocation strategy based on read length.

    Parameters
    ----------
    read_length : int32_t
        Length of read in bases

    Returns
    -------
    int
        1 for stack allocation (short reads), 0 for heap (long reads)
    """
    return 1 if read_length <= 1000 else 0

cdef inline uint8_t* allocate_sequence_buffers(int32_t read_length, uint8_t** ref_buffer_ptr) nogil:
    """Allocate buffers for read and reference sequences.

    Parameters
    ----------
    read_length : int32_t
        Length of read in bases
    ref_buffer_ptr : uint8_t**
        Output pointer for reference sequence buffer

    Returns
    -------
    uint8_t*
        Pointer to read sequence buffer, or NULL on allocation failure
    """
    cdef uint8_t* base_buffer = <uint8_t*>malloc(read_length * 2 * sizeof(uint8_t))
    if not base_buffer:
        return NULL
    ref_buffer_ptr[0] = base_buffer + read_length
    return base_buffer

cdef double calculate_md_quality_score(bam1_t* alignment,
                                       sam_hdr_t* header,
                                       AlignmentScoringConfig* config,
                                       float* pmd_result) noexcept nogil:
    """Calculate alignment quality score from MD tag with optional PMD calculation.

    Main entry point for MD-based scoring. Routes to fast path (no PMD) or
    full path (with PMD calculation) based on configuration.

    Parameters
    ----------
    alignment : bam1_t*
        BAM alignment record
    header : sam_hdr_t*
        BAM header (unused, retained for API compatibility)
    config : AlignmentScoringConfig*
        Scoring configuration
    pmd_result : float*
        Output for PMD score (NULL if not needed)

    Returns
    -------
    double
        Log-likelihood score, or -1e20 on error
    """
    if not LOOKUP_TABLES_INITIALIZED:
        initialize_quality_lookup_tables()
    if not alignment:
        return -1e20

    cdef int32_t read_length = alignment.core.l_qseq
    if read_length <= 0 or read_length > 50000:
        return -1e20

    cdef uint8_t* qual_data = bam_get_qual(alignment)
    cdef uint8_t* seq_data = bam_get_seq(alignment)
    cdef uint8_t* md_aux = bam_aux_get(alignment, b"MD")
    if not qual_data or not seq_data or not md_aux:
        return -1e20

    cdef bint calculate_pmd = (pmd_result != NULL and config.calculate_pmd)
    cdef bint is_single_stranded = config.is_single_stranded if calculate_pmd else False

    if not calculate_pmd:
        return calculate_md_score_fast_path(<char*>(md_aux + 1), qual_data, read_length)

    cdef uint8_t* ref_sequence = NULL
    cdef uint8_t* read_bases = allocate_sequence_buffers(read_length, &ref_sequence)
    if not read_bases:
        return -1e20

    cdef double result
    try:
        decode_read_bases(seq_data, read_length, read_bases)
        result = calculate_md_score_with_pmd_single_pass(
            <char*>(md_aux + 1), read_bases, ref_sequence, qual_data,
            read_length, is_single_stranded, pmd_result
        )
    finally:
        free(read_bases)

    return result

cdef double calculate_md_score_fast_path(char* md_tag, uint8_t* qual_data, int32_t read_length) noexcept nogil:
    """Calculate MD-based alignment score without PMD calculation (fast path).

    Parses MD tag to count matches, mismatches, and deletions, computing
    log-likelihood based on quality scores. Optimized for speed when PMD
    calculation is not needed.

    Parameters
    ----------
    md_tag : char*
        MD tag string from BAM record
    qual_data : uint8_t*
        Quality score array
    read_length : int32_t
        Length of read in bases

    Returns
    -------
    double
        Total log-likelihood score
    """
    cdef char* ptr = md_tag
    cdef int32_t read_pos = 0
    cdef int32_t match_count = 0, mismatch_count = 0, deletion_count = 0
    cdef double total_log_likelihood = 0.0
    cdef int32_t num, end_pos, i
    cdef uint8_t qual_score
    cdef int qual_idx
    cdef uint8_t* qptr
    cdef uint8_t* qend
    cdef int cnt
    cdef double match_accumulator

    while ptr[0] != 0 and read_pos < read_length:
        if ptr[0] >= 48 and ptr[0] <= 57:
            num = parse_int(&ptr)
            match_count += num
            end_pos = read_pos + num
            if end_pos > read_length:
                end_pos = read_length

            qptr = qual_data + read_pos
            qend = qual_data + end_pos
            match_accumulator = 0.0
            while qptr < qend:
                qual_score = qptr[0]
                qual_idx = qual_score if qual_score <= 93 else (SENTINEL_Q if qual_score == 255 else 93)
                cnt = 1
                while qptr + cnt < qend and qptr[cnt] == qual_score:
                    cnt += 1
                match_accumulator += cnt * PRECOMPUTED_LOG_P_CORRECT[qual_idx]
                qptr += cnt
            total_log_likelihood += match_accumulator
            read_pos = end_pos

        elif ptr[0] == 94:
            ptr += 1
            while ptr[0] != 0 and IS_UPPER[ptr[0]]:
                deletion_count += 1
                ptr += 1

        elif IS_UPPER[ptr[0]]:
            if read_pos < read_length:
                mismatch_count += 1
                qual_score = qual_data[read_pos]
                qual_idx = qual_score if qual_score <= 93 else (SENTINEL_Q if qual_score == 255 else 93)
                total_log_likelihood += PRECOMPUTED_LOG_P_ERROR[qual_idx]
                read_pos += 1
            ptr += 1
        else:
            ptr += 1

    total_log_likelihood += deletion_count * LOG_DELETION_PROB
    return total_log_likelihood

cdef double calculate_md_score_with_pmd_single_pass(char* md_tag, uint8_t* read_bases,
                                                    uint8_t* ref_sequence, uint8_t* qual_data,
                                                    int32_t read_length, bint is_single_stranded,
                                                    float* pmd_result) noexcept nogil:
    """Calculate alignment score with integrated PMD calculation (single pass).

    Parses MD tag while simultaneously computing PMD likelihood for C→T and G→A
    transitions characteristic of ancient DNA damage. Computes both alignment
    quality score and PMD score in a single pass for efficiency.

    Parameters
    ----------
    md_tag : char*
        MD tag string from BAM record
    read_bases : uint8_t*
        Decoded read sequence (ASCII)
    ref_sequence : uint8_t*
        Buffer for reconstructed reference sequence
    qual_data : uint8_t*
        Quality score array
    read_length : int32_t
        Length of read in bases
    is_single_stranded : bint
        True for single-stranded library (both ends damaged)
    pmd_result : float*
        Output for PMD log-likelihood ratio

    Returns
    -------
    double
        Total alignment log-likelihood score
    """
    cdef char* ptr = md_tag
    cdef int32_t read_pos = 0, match_count = 0, mismatch_count = 0, deletion_count = 0
    cdef int32_t pmd_corrected_mismatches = 0
    cdef double total_log_likelihood = 0.0
    cdef double log_pmd_likelihood = 0.0
    cdef double log_null_likelihood = 0.0

    cdef double log_pi = -6.907755278982137
    cdef double log_1_minus_pi = -0.001000500333583532
    cdef double log_C = -4.605170185988091
    cdef double log_1_minus_C = -0.010050335853501442

    cdef double* LOG_EPS_cached = LOG_EPSILON
    cdef double* LOG1M_EPS_cached = LOG1M_EPSILON
    cdef unsigned char* IS_C_cached = IS_C
    cdef unsigned char* IS_G_cached = IS_G
    cdef unsigned char* IS_UPPER_cached = IS_UPPER

    cdef uint8_t* rptr = read_bases
    cdef uint8_t* qptr = qual_data
    cdef int32_t rlen = read_length

    cdef int32_t num, end_pos, k
    cdef uint8_t ref_base, read_base, qual_score
    cdef int qual_idx
    cdef double log_eps, log_1meps
    cdef int32_t z_from_5prime, z_from_3prime
    cdef double tmp_logDz, tmp_log1mDz, tmp_logDy, tmp_log1mDy
    cdef double log_pmd_comp, log_null_comp
    cdef double log_p_err_pmd, log_p_err_null, lbf
    cdef unsigned char ch
    cdef bint is_c_base, is_g_base, is_ct_mismatch, is_ga_mismatch

    cdef double shared_term1, shared_term2, shared_1meps_1, shared_1meps_2
    cdef double shared_g1, shared_g2, base_term1, base_term2
    cdef double null_base1, null_base2, pmd_base1, pmd_base2
    cdef double null_shared1, null_shared2, at_comp, error_penalty

    while ptr[0] != 0 and read_pos < rlen:
        ch = <unsigned char>ptr[0]
        if ch >= 48 and ch <= 57:
            num = parse_int(&ptr)
            match_count += num
            end_pos = read_pos + num
            if end_pos > rlen:
                end_pos = rlen

            for k in range(read_pos, end_pos):
                read_base = rptr[k]
                ref_sequence[k] = read_base
                ref_base = read_base
                qual_score = qptr[k]
                qual_idx = qual_score if qual_score <= 93 else (SENTINEL_Q if qual_score == 255 else 93)
                is_c_base = IS_C_cached[ref_base]
                is_g_base = IS_G_cached[ref_base]
                if is_c_base or is_g_base:
                    log_eps = LOG_EPS_cached[qual_idx]
                    log_1meps = LOG1M_EPS_cached[qual_idx]
                    if is_c_base:
                        z_from_5prime = k + 1
                        if is_single_stranded:
                            z_from_3prime = rlen - k
                            get_D_terms(z_from_5prime, &tmp_logDz, &tmp_log1mDz)
                            get_D_terms(z_from_3prime, &tmp_logDy, &tmp_log1mDy)
                            shared_term1 = log_1_minus_pi + log_1meps
                            shared_term2 = log_pi + log_1meps
                            log_pmd_comp = lse2((shared_term1 + tmp_log1mDz + tmp_log1mDy),
                                                (shared_term2 + tmp_log1mDz + tmp_log1mDy))
                            log_null_comp = lse2((shared_term1 + log_1_minus_C),
                                                 (shared_term2 + log_1_minus_C))
                            total_log_likelihood += lse2(tmp_log1mDz + log_null_comp, tmp_logDz + log_pmd_comp)
                        else:
                            get_D_terms(z_from_5prime, &tmp_logDz, &tmp_log1mDz)
                            shared_1meps_1 = log_1_minus_pi + log_1meps
                            shared_1meps_2 = log_pi + log_1meps
                            log_pmd_comp = lse2((shared_1meps_1 + tmp_log1mDz),
                                                (shared_1meps_2 + tmp_log1mDz))
                            log_null_comp = lse2((shared_1meps_1 + log_1_minus_C),
                                                 (shared_1meps_2 + log_1_minus_C))
                            total_log_likelihood += lse2(tmp_log1mDz + log_null_comp, tmp_logDz + log_pmd_comp)
                    else:
                        if not is_single_stranded:
                            z_from_3prime = rlen - k
                            get_D_terms(z_from_3prime, &tmp_logDz, &tmp_log1mDz)
                            shared_g1 = log_1_minus_pi + log_1meps
                            shared_g2 = log_pi + log_1meps
                            log_pmd_comp = lse2((shared_g1 + tmp_log1mDz),
                                                (shared_g2 + tmp_log1mDz))
                            log_null_comp = lse2((shared_g1 + log_1_minus_C),
                                                 (shared_g2 + log_1_minus_C))
                            total_log_likelihood += lse2(tmp_log1mDz + log_null_comp, tmp_logDz + log_pmd_comp)
                        else:
                            log_pmd_comp = LOG1M_EPS_cached[qual_idx]
                            log_null_comp = lse2((log_1_minus_pi + log_1meps + log_1_minus_C),
                                                 (log_pi + log_1meps + log_1_minus_C))
                            total_log_likelihood += log_pmd_comp
                    log_pmd_likelihood += log_pmd_comp
                    log_null_likelihood += log_null_comp
                else:
                    at_comp = LOG1M_EPS_cached[qual_idx]
                    log_pmd_likelihood += at_comp
                    log_null_likelihood += at_comp
                    total_log_likelihood += at_comp
            read_pos = end_pos
        elif ch == 94:
            ptr += 1
            while ptr[0] != 0 and IS_UPPER_cached[<unsigned char>ptr[0]]:
                deletion_count += 1
                ptr += 1
        elif IS_UPPER_cached[ch]:
            if read_pos < rlen:
                mismatch_count += 1
                ref_base = <uint8_t>ch
                read_base = rptr[read_pos]
                ref_sequence[read_pos] = ref_base
                qual_score = qptr[read_pos]
                qual_idx = qual_score if qual_score <= 93 else (SENTINEL_Q if qual_score == 255 else 93)
                is_ct_mismatch = IS_C_cached[ref_base] and (read_base == 84)
                is_ga_mismatch = IS_G_cached[ref_base] and (read_base == 65)
                if is_ct_mismatch or is_ga_mismatch:
                    log_eps = LOG_EPS_cached[qual_idx]
                    log_1meps = LOG1M_EPS_cached[qual_idx]
                    if is_single_stranded and is_ct_mismatch:
                        z_from_5prime = read_pos + 1
                        z_from_3prime = rlen - read_pos
                        get_D_terms(z_from_5prime, &tmp_logDz, &tmp_log1mDz)
                        get_D_terms(z_from_3prime, &tmp_logDy, &tmp_log1mDy)
                        base_term1 = log_1_minus_pi + tmp_log1mDz + tmp_log1mDy
                        base_term2 = log_pi + tmp_log1mDz + tmp_log1mDy
                        log_p_err_pmd = lse4((base_term1 + log_eps),
                                             (log_1_minus_pi + log_1meps + tmp_logDz + tmp_log1mDy),
                                             (log_1_minus_pi + log_1meps + tmp_log1mDz + tmp_logDy),
                                             (base_term2 + log_eps))
                        null_base1 = log_1_minus_pi + log_1_minus_C
                        null_base2 = log_pi + log_1_minus_C
                        log_p_err_null = lse4((null_base1 + log_eps + log_1_minus_C),
                                              (log_1_minus_pi + log_1meps + log_C + log_1_minus_C),
                                              (log_1_minus_pi + log_1meps + log_1_minus_C + log_C),
                                              (null_base2 + log_eps + log_1_minus_C))
                        lbf = log_p_err_pmd - log_p_err_null
                        if lbf > 0.0:
                            pmd_corrected_mismatches += 1
                        log_pmd_likelihood += log_p_err_pmd
                        log_null_likelihood += log_p_err_null
                        total_log_likelihood += lse2(tmp_log1mDz + log_p_err_null, tmp_logDz + log_p_err_pmd)
                    else:
                        if is_ct_mismatch:
                            z_from_5prime = read_pos + 1
                            get_D_terms(z_from_5prime, &tmp_logDz, &tmp_log1mDz)
                        else:
                            z_from_3prime = rlen - read_pos
                            get_D_terms(z_from_3prime, &tmp_logDz, &tmp_log1mDz)
                        pmd_base1 = log_1_minus_pi + tmp_log1mDz
                        pmd_base2 = log_pi + tmp_log1mDz
                        log_p_err_pmd = lse3((pmd_base1 + log_eps),
                                             (log_1_minus_pi + log_1meps + tmp_logDz),
                                             (pmd_base2 + log_eps))
                        null_shared1 = log_1_minus_pi + log_1_minus_C
                        null_shared2 = log_pi + log_1_minus_C
                        log_p_err_null = lse3((null_shared1 + log_eps),
                                              (log_1_minus_pi + log_1meps + log_C),
                                              (null_shared2 + log_eps))
                        lbf = log_p_err_pmd - log_p_err_null
                        if lbf > 0.0:
                            pmd_corrected_mismatches += 1
                        log_pmd_likelihood += log_p_err_pmd
                        log_null_likelihood += log_p_err_null
                        total_log_likelihood += lse2(tmp_log1mDz + log_p_err_null, tmp_logDz + log_p_err_pmd)
                else:
                    error_penalty = LOG_EPS_cached[qual_idx]
                    total_log_likelihood += error_penalty
                    log_pmd_likelihood += error_penalty
                    log_null_likelihood += error_penalty
                read_pos += 1
            ptr += 1
        else:
            ptr += 1

    total_log_likelihood += deletion_count * LOG_DELETION_PROB
    if pmd_result:
        pmd_result[0] = <float>(log_pmd_likelihood - log_null_likelihood)
    return total_log_likelihood


cdef bint alignment_passes_quality_filters(bam1_t* alignment, AlignmentScoringConfig* config) noexcept nogil:
    """Check if alignment passes quality filters (length and ANI).

    Parameters
    ----------
    alignment : bam1_t*
        BAM alignment record
    config : AlignmentScoringConfig*
        Configuration with filter thresholds (NULL uses defaults)

    Returns
    -------
    bint
        True if alignment passes all filters
    """
    cdef int l_qseq
    cdef uint8_t* aux
    cdef int32_t nm_val
    cdef double pct_id
    cdef int nm_val2
    cdef double pct_id2
    cdef int32_t min_read_length = 30
    cdef int32_t max_read_length = 10000
    cdef double min_read_identity = 90.0

    if not alignment:
        return False
    if config:
        min_read_length = config.minimum_read_length
        max_read_length = config.maximum_read_length
        min_read_identity = config.minimum_read_identity

    l_qseq = alignment.core.l_qseq
    if l_qseq < min_read_length or l_qseq > max_read_length:
        return False

    aux = bam_aux_get(alignment, b"NM")
    if aux != NULL:
        if aux[0] == ord('i'):
            nm_val = (<int32_t*>(aux + 1))[0]
            if nm_val >= 0 and l_qseq > 0:
                pct_id = (1.0 - (nm_val / <double>l_qseq)) * 100.0
                if pct_id < min_read_identity:
                    return False
        else:
            nm_val2 = bam_aux2i(aux)
            if nm_val2 >= 0 and l_qseq > 0:
                pct_id2 = (1.0 - (nm_val2 / <double>l_qseq)) * 100.0
                if pct_id2 < min_read_identity:
                    return False

    return True


cdef bint alignment_passes_quality_filters_with_ani(bam1_t* alignment, AlignmentScoringConfig* config, double* ani_out, int32_t* nm_out) noexcept nogil:
    """Check quality filters and return ANI and NM values.

    Parameters
    ----------
    alignment : bam1_t*
        BAM alignment record
    config : AlignmentScoringConfig*
        Configuration with filter thresholds
    ani_out : double*
        Output for percent identity (NULL if not needed)
    nm_out : int32_t*
        Output for NM tag value (NULL if not needed)

    Returns
    -------
    bint
        True if alignment passes all filters
    """
    cdef int32_t l_qseq = alignment.core.l_qseq
    cdef uint8_t* aux
    cdef int32_t nm_val = -1
    cdef double pct_id

    if config.minimum_read_length > 0 and l_qseq < config.minimum_read_length:
        return False
    if config.maximum_read_length > 0 and l_qseq > config.maximum_read_length:
        return False

    aux = bam_aux_get(alignment, b"NM")
    if aux != NULL:
        if aux[0] == ord('i'):
            nm_val = (<int32_t*>(aux + 1))[0]
        else:
            nm_val = bam_aux2i(aux)
    else:
        nm_val = -1

    if nm_val >= 0 and l_qseq > 0:
        pct_id = (1.0 - (nm_val / <double>l_qseq)) * 100.0
    else:
        pct_id = 0.0

    if ani_out != NULL:
        ani_out[0] = pct_id
    if nm_out != NULL:
        nm_out[0] = nm_val

    if config.minimum_read_identity > 0.0 and pct_id < config.minimum_read_identity:
        return False

    return True
