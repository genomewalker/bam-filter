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

"""EM algorithm for read reassignment with SQUAREM acceleration."""

from cython.parallel import prange, threadid

from libc.math cimport exp, fabs, fmax, fmin, log, log2, sqrt as libc_sqrt, pow as libc_pow, INFINITY
from libc.stdint cimport int32_t, int64_t, uint16_t, uint32_t, uint64_t, uint8_t, uintptr_t
from libc.stdlib cimport calloc, free, malloc
from libc.string cimport memcpy, memset
from libc.float cimport DBL_EPSILON

from bam_filter.processor cimport MemoryPool, Alignment
from bam_filter.processor_fast_math cimport stable_log_sum_exp, safe_normalize_weights
from cpython.pycapsule cimport PyCapsule_GetPointer

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil


# =============================================================================
# Constants
# =============================================================================

cdef double NEG_INF = -1e20
cdef double LOG_ZERO = -1e10


# =============================================================================
# Helper Functions
# =============================================================================

cdef inline double* get_reference_weights(MemoryPool* pool) noexcept nogil:
    """Get pointer to reference weights array in unified buffer."""
    return pool.unified_buffer + pool.reference_weights_offset


cdef extern from *:
    """
    #ifdef __GNUC__
    #define PREFETCH_R(addr) __builtin_prefetch((addr), 0, 3)
    #define PREFETCH_W(addr) __builtin_prefetch((addr), 1, 3)
    #else
    #define PREFETCH_R(addr) ((void)0)
    #define PREFETCH_W(addr) ((void)0)
    #endif
    """
    void PREFETCH_R(void* addr) nogil
    void PREFETCH_W(void* addr) nogil


cdef inline void PREFETCH_READ(void* ptr) noexcept nogil:
    PREFETCH_R(ptr)


cdef inline void PREFETCH_WRITE(void* ptr) noexcept nogil:
    PREFETCH_W(ptr)


# =============================================================================
# State Management
# =============================================================================

cdef EMState* create_em_state(uint32_t n_refs, uint32_t n_reads,
                               EMConfig* config) noexcept nogil:
    """Allocate and initialize EMState structure."""
    cdef EMState* state = <EMState*>malloc(sizeof(EMState))
    if state == NULL:
        return NULL

    state.n_refs = n_refs
    state.n_reads = n_reads

    # Allocate phi weights and counts
    state.phi_weights = <double*>calloc(n_refs, sizeof(double))
    state.phi_counts = <double*>calloc(n_refs, sizeof(double))

    if state.phi_weights == NULL or state.phi_counts == NULL:
        free_em_state(state)
        return NULL

    # Initialize phi to uniform
    cdef double init_weight = 1.0 / <double>n_refs
    cdef uint32_t j
    for j in range(n_refs):
        state.phi_weights[j] = init_weight

    # Hierarchical arrays (if enabled)
    state.hierarchical_enabled = config.hierarchical_enabled
    state.gamma_values = NULL
    state.S_anc = NULL
    state.S_mod = NULL

    if config.hierarchical_enabled:
        state.gamma_values = <double*>calloc(n_refs, sizeof(double))
        state.S_anc = <double*>calloc(n_refs, sizeof(double))
        state.S_mod = <double*>calloc(n_refs, sizeof(double))

        if state.gamma_values == NULL or state.S_anc == NULL or state.S_mod == NULL:
            free_em_state(state)
            return NULL

        # Initialize gamma to 0.5 (uninformative)
        for j in range(n_refs):
            state.gamma_values[j] = 0.5

    # Unknown component
    state.unknown_enabled = config.unknown_enabled
    state.phi_unknown = 0.05 if config.unknown_enabled else 0.0
    state.S_unknown = 0.0

    # PMD priors (allocated but not initialized here - caller must set)
    state.omega_ancient = NULL
    if config.hierarchical_enabled:
        state.omega_ancient = <double*>calloc(n_reads, sizeof(double))
        if state.omega_ancient == NULL:
            free_em_state(state)
            return NULL
        # Default: uninformative prior (0.5)
        for j in range(n_reads):
            state.omega_ancient[j] = 0.5

    # Copy config values
    state.dirichlet_prior = config.dirichlet_prior
    state.gamma_prior = config.gamma_prior
    state.power_rho = config.power_rho
    state.unknown_margin = config.unknown_margin
    state.em_beta = 1.0  # Default temperature

    return state


cdef void free_em_state(EMState* state) noexcept nogil:
    """Free EMState and all its arrays."""
    if state == NULL:
        return

    if state.phi_weights != NULL:
        free(state.phi_weights)
    if state.phi_counts != NULL:
        free(state.phi_counts)
    if state.gamma_values != NULL:
        free(state.gamma_values)
    if state.S_anc != NULL:
        free(state.S_anc)
    if state.S_mod != NULL:
        free(state.S_mod)
    if state.omega_ancient != NULL:
        free(state.omega_ancient)

    free(state)


cdef void copy_em_state(EMState* dest, EMState* src) noexcept nogil:
    """Deep copy EMState (assumes dest is already allocated with same dimensions)."""
    if dest == NULL or src == NULL:
        return

    memcpy(dest.phi_weights, src.phi_weights, src.n_refs * sizeof(double))
    memcpy(dest.phi_counts, src.phi_counts, src.n_refs * sizeof(double))

    if src.hierarchical_enabled and dest.gamma_values != NULL:
        memcpy(dest.gamma_values, src.gamma_values, src.n_refs * sizeof(double))
        memcpy(dest.S_anc, src.S_anc, src.n_refs * sizeof(double))
        memcpy(dest.S_mod, src.S_mod, src.n_refs * sizeof(double))

    dest.phi_unknown = src.phi_unknown
    dest.S_unknown = src.S_unknown

    # omega_ancient is FIXED, no need to copy


cdef void reset_accumulators(EMState* state) noexcept nogil:
    """Reset E-step accumulators to zero."""
    if state == NULL:
        return

    memset(state.phi_counts, 0, state.n_refs * sizeof(double))

    if state.hierarchical_enabled:
        memset(state.S_anc, 0, state.n_refs * sizeof(double))
        memset(state.S_mod, 0, state.n_refs * sizeof(double))

    state.S_unknown = 0.0


cdef void normalize_phi(EMState* state) noexcept nogil:
    """Normalize phi weights to sum to 1 (including unknown if enabled)."""
    if state == NULL:
        return

    cdef double total = 0.0
    cdef uint32_t j

    for j in range(state.n_refs):
        state.phi_weights[j] = fmax(state.phi_weights[j], 1e-15)
        total += state.phi_weights[j]

    if state.unknown_enabled:
        state.phi_unknown = fmax(state.phi_unknown, 1e-15)
        total += state.phi_unknown

    if total > 0:
        for j in range(state.n_refs):
            state.phi_weights[j] /= total
        if state.unknown_enabled:
            state.phi_unknown /= total


# =============================================================================
# Likelihood Computation Helpers
# =============================================================================

cdef inline double compute_log_L_ancient(Alignment* aln,
                                          float D_avg_5p, float D_avg_3p,
                                          float epsilon) noexcept nogil:
    """
    Log-likelihood under ancient DNA model.

    For ancient DNA: expect damage at terminal positions (C->T at 5', G->A at 3').
    """
    cdef double log_L = 0.0
    cdef uint16_t aligned = aln.aligned_length
    cdef uint16_t matches = aln.match_count
    cdef uint8_t ct_5p = aln.ct_5p_count
    cdef uint8_t ga_3p = aln.ga_3p_count
    cdef uint8_t c_at_5p = aln.c_at_5p_count
    cdef uint8_t g_at_3p = aln.g_at_3p_count
    cdef int other_mm
    cdef double D_avg, actual_opp, survived, observed_damage

    if aligned == 0:
        return LOG_ZERO

    other_mm = aligned - matches - ct_5p - ga_3p
    if other_mm < 0:
        other_mm = 0

    D_avg = (D_avg_5p + D_avg_3p) / 2.0
    actual_opp = <double>(c_at_5p + g_at_3p)
    observed_damage = <double>(ct_5p + ga_3p)
    survived = fmax(0.0, actual_opp - observed_damage)

    cdef double p_damage_ancient = D_avg + (1.0 - D_avg) * epsilon / 3.0
    cdef double p_survive_ancient = (1.0 - D_avg) * (1.0 - epsilon / 3.0)
    cdef double log_p_damage = log(fmax(p_damage_ancient, 1e-10))
    cdef double log_p_survive = log(fmax(p_survive_ancient, 1e-10))
    cdef double log_p_match = log(fmax(1.0 - epsilon, 1e-10))
    cdef double log_p_error = log(fmax(epsilon / 3.0, 1e-10))

    log_L += observed_damage * log_p_damage
    log_L += survived * log_p_survive
    log_L += matches * log_p_match
    log_L += other_mm * log_p_error

    return log_L


cdef inline double compute_log_L_modern(Alignment* aln, float epsilon) noexcept nogil:
    """
    Log-likelihood under modern DNA model.

    For modern DNA: C->T and G->A are just sequencing errors.
    """
    cdef double log_L = 0.0
    cdef uint16_t aligned = aln.aligned_length
    cdef uint16_t matches = aln.match_count
    cdef uint8_t ct_5p = aln.ct_5p_count
    cdef uint8_t ga_3p = aln.ga_3p_count
    cdef uint8_t c_at_5p = aln.c_at_5p_count
    cdef uint8_t g_at_3p = aln.g_at_3p_count
    cdef int other_mm
    cdef double actual_opp, survived, observed_damage

    if aligned == 0:
        return LOG_ZERO

    other_mm = aligned - matches - ct_5p - ga_3p
    if other_mm < 0:
        other_mm = 0

    actual_opp = <double>(c_at_5p + g_at_3p)
    observed_damage = <double>(ct_5p + ga_3p)
    survived = fmax(0.0, actual_opp - observed_damage)

    cdef double log_p_error = log(fmax(epsilon / 3.0, 1e-10))
    cdef double log_p_survive = log(fmax(1.0 - epsilon / 3.0, 1e-10))
    cdef double log_p_match = log(fmax(1.0 - epsilon, 1e-10))

    log_L += observed_damage * log_p_error
    log_L += survived * log_p_survive
    log_L += matches * log_p_match
    log_L += other_mm * log_p_error

    return log_L


# =============================================================================
# E-Step with optimizations
# =============================================================================

# Maximum alignments per read for scratch buffer (stack allocation)
DEF MAX_SCRATCH_SIZE = 64


cdef inline void process_read_single(
    EMState* state, MemoryPool* pool, double* log_phi,
    uint32_t read_idx, uint64_t start_pos,
    double* accum_phi, double* accum_S_anc, double* accum_S_mod,
    float D_5p, float D_3p, float epsilon,
    bint use_hierarchical
) noexcept nogil:
    """Fast path for reads with single alignment."""
    cdef Alignment* aln = &pool.alignments[start_pos]
    cdef uint32_t ref_idx = aln.reference_index
    cdef double gamma_j, omega_anc, log_omega_anc, log_omega_mod
    cdef double log_gamma, log_1m_gamma, log_L_anc, log_L_mod, log_L_mix, p_anc

    if ref_idx >= state.n_refs:
        return

    accum_phi[ref_idx] += 1.0

    if use_hierarchical:
        gamma_j = state.gamma_values[ref_idx]
        omega_anc = 0.5
        if state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(1.0 - omega_anc, 1e-10))
        log_gamma = log(fmax(gamma_j, 1e-10))
        log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))
        log_L_anc = compute_log_L_ancient(aln, D_5p, D_3p, epsilon)
        log_L_mod = compute_log_L_modern(aln, epsilon)
        log_L_mix = stable_log_sum_exp(
            log_gamma + log_omega_anc + log_L_anc,
            log_1m_gamma + log_omega_mod + log_L_mod
        )
        p_anc = exp(log_gamma + log_omega_anc + log_L_anc - log_L_mix)
        p_anc = fmax(fmin(p_anc, 0.999), 0.001)
        accum_S_anc[ref_idx] += p_anc
        accum_S_mod[ref_idx] += 1.0 - p_anc


cdef inline void process_read_multi(
    EMState* state, MemoryPool* pool, double* log_phi, double log_phi_u,
    uint32_t read_idx, uint64_t start_pos, uint64_t end_pos, uint32_t alignment_count,
    double* accum_phi, double* accum_S_anc, double* accum_S_mod, double* accum_S_unknown,
    float D_5p, float D_3p, float epsilon, double unknown_margin,
    bint use_hierarchical, bint use_unknown
) noexcept nogil:
    """Process read with multiple alignments using scratch buffer."""
    cdef uint32_t n_refs = state.n_refs
    cdef uint32_t ref_idx, i, alignment_idx
    cdef double log_normalizer, log_weighted, posterior, log_max, s_max, s_unknown
    cdef double alignment_score, omega_anc, omega_mod, log_omega_anc, log_omega_mod
    cdef double gamma_j, log_gamma, log_1m_gamma, log_L_anc, log_L_mod, log_L_mix
    cdef double log_anc_term, p_anc, unknown_posterior
    cdef Alignment* aln

    # Scratch buffers (stack for small, heap for large)
    cdef double scratch_log_w[MAX_SCRATCH_SIZE]
    cdef double scratch_p_anc[MAX_SCRATCH_SIZE]
    cdef uint32_t scratch_refs[MAX_SCRATCH_SIZE]
    cdef double* heap_log_w = NULL
    cdef double* heap_p_anc = NULL
    cdef uint32_t* heap_refs = NULL
    cdef double* log_w_buf
    cdef double* p_anc_buf
    cdef uint32_t* ref_buf

    # Select buffer
    if alignment_count <= MAX_SCRATCH_SIZE:
        log_w_buf = scratch_log_w
        p_anc_buf = scratch_p_anc
        ref_buf = scratch_refs
    else:
        heap_log_w = <double*>malloc(alignment_count * sizeof(double))
        heap_refs = <uint32_t*>malloc(alignment_count * sizeof(uint32_t))
        if use_hierarchical:
            heap_p_anc = <double*>malloc(alignment_count * sizeof(double))
        if heap_log_w == NULL or heap_refs == NULL:
            if heap_log_w != NULL: free(heap_log_w)
            if heap_refs != NULL: free(heap_refs)
            if heap_p_anc != NULL: free(heap_p_anc)
            return
        log_w_buf = heap_log_w
        p_anc_buf = heap_p_anc
        ref_buf = heap_refs

    # Get PMD prior
    omega_anc = 0.5
    omega_mod = 0.5
    if use_hierarchical and state.omega_ancient != NULL:
        omega_anc = state.omega_ancient[read_idx]
        omega_mod = 1.0 - omega_anc
    log_omega_anc = log(fmax(omega_anc, 1e-10))
    log_omega_mod = log(fmax(omega_mod, 1e-10))

    # First pass: compute log weights and find max
    log_max = NEG_INF
    s_max = NEG_INF
    i = 0

    for alignment_idx in range(start_pos, end_pos):
        aln = &pool.alignments[alignment_idx]
        ref_idx = aln.reference_index

        # Prefetch next alignment
        if alignment_idx + 1 < end_pos:
            PREFETCH_READ(&pool.alignments[alignment_idx + 1])

        if ref_idx >= n_refs:
            continue

        ref_buf[i] = ref_idx
        alignment_score = <double>aln.alignment_score

        if alignment_score > s_max:
            s_max = alignment_score

        if use_hierarchical:
            gamma_j = state.gamma_values[ref_idx]
            log_gamma = log(fmax(gamma_j, 1e-10))
            log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))
            log_L_anc = compute_log_L_ancient(aln, D_5p, D_3p, epsilon)
            log_L_mod = compute_log_L_modern(aln, epsilon)
            log_anc_term = log_gamma + log_omega_anc + log_L_anc
            log_L_mix = stable_log_sum_exp(
                log_anc_term,
                log_1m_gamma + log_omega_mod + log_L_mod
            )
            log_weighted = log_phi[ref_idx] + log_L_mix
            p_anc_buf[i] = exp(log_anc_term - log_L_mix)
            p_anc_buf[i] = fmax(fmin(p_anc_buf[i], 0.999), 0.001)
        else:
            log_weighted = log_phi[ref_idx] + alignment_score

        log_w_buf[i] = log_weighted
        if log_weighted > log_max:
            log_max = log_weighted
        i += 1

    if i == 0:
        if heap_log_w != NULL: free(heap_log_w)
        if heap_refs != NULL: free(heap_refs)
        if heap_p_anc != NULL: free(heap_p_anc)
        return

    # Compute normalizer with max-subtraction
    log_normalizer = 0.0
    for alignment_idx in range(i):
        log_normalizer += exp(log_w_buf[alignment_idx] - log_max)

    if use_unknown:
        s_unknown = s_max - unknown_margin
        log_weighted = log_phi_u + s_unknown
        log_normalizer += exp(log_weighted - log_max)

    log_normalizer = log_max + log(log_normalizer)

    if log_normalizer <= NEG_INF + 1e10:
        if heap_log_w != NULL: free(heap_log_w)
        if heap_refs != NULL: free(heap_refs)
        if heap_p_anc != NULL: free(heap_p_anc)
        return

    # Second pass: compute posteriors from cached values
    for alignment_idx in range(i):
        ref_idx = ref_buf[alignment_idx]
        posterior = exp(log_w_buf[alignment_idx] - log_normalizer)
        posterior = fmax(fmin(posterior, 0.999), 1e-12)

        PREFETCH_WRITE(&accum_phi[ref_idx])
        accum_phi[ref_idx] += posterior

        if use_hierarchical:
            p_anc = p_anc_buf[alignment_idx]
            accum_S_anc[ref_idx] += posterior * p_anc
            accum_S_mod[ref_idx] += posterior * (1.0 - p_anc)

    if use_unknown:
        unknown_posterior = exp(log_phi_u + s_unknown - log_normalizer)
        unknown_posterior = fmax(fmin(unknown_posterior, 0.999), 1e-12)
        accum_S_unknown[0] += unknown_posterior

    if heap_log_w != NULL: free(heap_log_w)
    if heap_refs != NULL: free(heap_refs)
    if heap_p_anc != NULL: free(heap_p_anc)


cdef void e_step(EMState* state, void* pool_ptr,
                 EMConfig* config) noexcept nogil:
    """E-step computing responsibilities with optimizations for speed."""
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t read_idx, ref_idx
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count, n_refs = state.n_refs
    cdef uint32_t n_reads = pool.unique_read_count

    cdef float D_5p = config.D_avg_5p
    cdef float D_3p = config.D_avg_3p
    cdef float epsilon = config.epsilon_error
    cdef double unknown_margin = config.unknown_margin
    cdef bint use_unknown = state.unknown_enabled
    cdef bint use_hierarchical = state.hierarchical_enabled

    cdef double* log_phi = <double*>malloc(n_refs * sizeof(double))
    cdef double log_phi_u

    # Thread-local accumulators
    cdef int n_threads = config.thread_count if config.thread_count > 0 else 1
    cdef int tid
    cdef double* thread_phi_counts = NULL
    cdef double* thread_S_anc = NULL
    cdef double* thread_S_mod = NULL
    cdef double* thread_S_unknown = NULL
    cdef size_t thread_buf_size
    cdef double* my_phi
    cdef double* my_S_anc
    cdef double* my_S_mod
    cdef double* my_S_unknown

    if log_phi == NULL:
        return

    # Reset accumulators
    reset_accumulators(state)

    # Precompute log(phi_j)
    for ref_idx in range(n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    log_phi_u = log(fmax(state.phi_unknown, 1e-15)) if use_unknown else NEG_INF

    # Allocate thread-local accumulators
    if n_threads > 1:
        thread_buf_size = <size_t>n_threads * n_refs
        thread_phi_counts = <double*>calloc(thread_buf_size, sizeof(double))
        if use_hierarchical:
            thread_S_anc = <double*>calloc(thread_buf_size, sizeof(double))
            thread_S_mod = <double*>calloc(thread_buf_size, sizeof(double))
        if use_unknown:
            thread_S_unknown = <double*>calloc(n_threads, sizeof(double))

        if thread_phi_counts == NULL:
            n_threads = 1

    # Process reads in parallel
    for read_idx in prange(n_reads, nogil=True, num_threads=n_threads,
                           schedule='dynamic', chunksize=256):
        tid = threadid() if n_threads > 1 else 0
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Select accumulator buffers
        if n_threads > 1:
            my_phi = &thread_phi_counts[tid * n_refs]
            my_S_anc = &thread_S_anc[tid * n_refs] if use_hierarchical else NULL
            my_S_mod = &thread_S_mod[tid * n_refs] if use_hierarchical else NULL
            my_S_unknown = &thread_S_unknown[tid] if use_unknown else NULL
        else:
            my_phi = state.phi_counts
            my_S_anc = state.S_anc
            my_S_mod = state.S_mod
            my_S_unknown = &state.S_unknown

        # Fast path for single alignment
        if alignment_count == 1 and not use_unknown:
            process_read_single(
                state, pool, log_phi, read_idx, start_pos,
                my_phi, my_S_anc, my_S_mod,
                D_5p, D_3p, epsilon, use_hierarchical
            )
        else:
            process_read_multi(
                state, pool, log_phi, log_phi_u,
                read_idx, start_pos, end_pos, alignment_count,
                my_phi, my_S_anc, my_S_mod, my_S_unknown,
                D_5p, D_3p, epsilon, unknown_margin,
                use_hierarchical, use_unknown
            )

    # Reduce thread-local accumulators
    if n_threads > 1 and thread_phi_counts != NULL:
        for tid in range(n_threads):
            for ref_idx in range(n_refs):
                state.phi_counts[ref_idx] += thread_phi_counts[tid * n_refs + ref_idx]
                if use_hierarchical and thread_S_anc != NULL:
                    state.S_anc[ref_idx] += thread_S_anc[tid * n_refs + ref_idx]
                    state.S_mod[ref_idx] += thread_S_mod[tid * n_refs + ref_idx]
            if use_unknown and thread_S_unknown != NULL:
                state.S_unknown += thread_S_unknown[tid]

        free(thread_phi_counts)
        if thread_S_anc != NULL: free(thread_S_anc)
        if thread_S_mod != NULL: free(thread_S_mod)
        if thread_S_unknown != NULL: free(thread_S_unknown)

    free(log_phi)


# =============================================================================
# Unified M-Step (Standard EM in phi-space, NO power transform)
# =============================================================================

cdef void m_step(EMState* state, EMConfig* config) noexcept nogil:
    """
    Unified M-step updating all parameters using expected counts.

    CRITICAL: NO power transform here. We optimize in phi-space.

    phi_j^new = (E[c_j] + alpha) / (sum E[c_j'] + J*alpha + unknown terms)
    gamma_j^new = (S_anc_j + alpha_gamma) / (S_anc_j + S_mod_j + 2*alpha_gamma)
    """
    cdef double total_count = 0.0
    cdef double alpha = config.dirichlet_prior
    cdef double alpha_gamma = config.gamma_prior
    cdef uint32_t j
    cdef double denom

    # Compute total for normalization
    for j in range(state.n_refs):
        total_count += state.phi_counts[j] + alpha

    if state.unknown_enabled:
        total_count += state.S_unknown + alpha

    if total_count <= 0:
        total_count = 1.0

    # Update phi weights (MAP with Dirichlet prior)
    for j in range(state.n_refs):
        state.phi_weights[j] = (state.phi_counts[j] + alpha) / total_count

    # Update unknown weight
    if state.unknown_enabled:
        state.phi_unknown = (state.S_unknown + alpha) / total_count

    # Update gamma values (hierarchical)
    # Clamp to [eps, 1-eps] for numerical stability when taking logs
    cdef double gamma_eps = 1e-6
    if state.hierarchical_enabled:
        for j in range(state.n_refs):
            denom = state.S_anc[j] + state.S_mod[j] + 2.0 * alpha_gamma
            if denom > 0:
                state.gamma_values[j] = (state.S_anc[j] + alpha_gamma) / denom
                # Clamp for numerical safety
                state.gamma_values[j] = fmax(gamma_eps, fmin(1.0 - gamma_eps, state.gamma_values[j]))
            else:
                state.gamma_values[j] = 0.5  # Uninformative default

    # Ensure proper normalization
    normalize_phi(state)


# =============================================================================
# Log-Likelihood Computation
# =============================================================================

cdef double compute_log_likelihood(EMState* state, void* pool_ptr,
                                   EMConfig* config) noexcept nogil:
    """
    Compute total log-likelihood of current parameter values.

    log P(X | phi, gamma) = sum_i log sum_j phi_j × f(x_i | j, gamma_j)
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t read_idx, ref_idx, alignment_idx
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count
    cdef double total_ll = 0.0
    cdef double log_normalizer, log_weighted
    cdef double log_phi_j, alignment_score
    cdef int32_t valid_reads = 0

    # Hierarchical variables
    cdef double log_L_anc, log_L_mod, log_L_mix
    cdef double gamma_j, omega_anc, omega_mod
    cdef double log_gamma, log_1m_gamma, log_omega_anc, log_omega_mod
    cdef float D_5p = config.D_avg_5p
    cdef float D_3p = config.D_avg_3p
    cdef float epsilon = config.epsilon_error
    cdef bint use_hierarchical = state.hierarchical_enabled

    # Unknown variables
    cdef double log_phi_u, s_max, s_unknown
    cdef double unknown_margin = config.unknown_margin
    cdef bint use_unknown = state.unknown_enabled

    cdef Alignment* aln
    cdef double* log_phi = <double*>malloc(state.n_refs * sizeof(double))

    if log_phi == NULL:
        return NEG_INF

    # Precompute log(phi_j)
    for ref_idx in range(state.n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    log_phi_u = log(fmax(state.phi_unknown, 1e-15)) if use_unknown else NEG_INF

    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Get PMD prior
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # Find max score for unknown
        s_max = NEG_INF
        if use_unknown:
            for alignment_idx in range(start_pos, end_pos):
                aln = &pool.alignments[alignment_idx]
                if aln.alignment_score > s_max:
                    s_max = aln.alignment_score
            s_unknown = s_max - unknown_margin

        # Compute log P(read | phi)
        log_normalizer = NEG_INF

        for alignment_idx in range(start_pos, end_pos):
            aln = &pool.alignments[alignment_idx]
            ref_idx = aln.reference_index
            if ref_idx >= state.n_refs:
                continue

            alignment_score = <double>aln.alignment_score

            if use_hierarchical:
                gamma_j = state.gamma_values[ref_idx]
                log_gamma = log(fmax(gamma_j, 1e-10))
                log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                log_L_anc = compute_log_L_ancient(aln, D_5p, D_3p, epsilon)
                log_L_mod = compute_log_L_modern(aln, epsilon)

                log_L_mix = stable_log_sum_exp(
                    log_gamma + log_omega_anc + log_L_anc,
                    log_1m_gamma + log_omega_mod + log_L_mod
                )

                log_weighted = log_phi[ref_idx] + log_L_mix
            else:
                log_weighted = log_phi[ref_idx] + alignment_score

            log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted)

        if use_unknown:
            log_weighted = log_phi_u + s_unknown
            log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted)

        if log_normalizer > NEG_INF + 1e10:
            total_ll += log_normalizer
            valid_reads += 1

    free(log_phi)

    return total_ll / valid_reads if valid_reads > 0 else NEG_INF


# =============================================================================
# SQUAREM Acceleration
# =============================================================================

cdef SQUAREMState* create_squarem_state(uint32_t dimension) noexcept nogil:
    """Allocate SQUAREM working arrays."""
    cdef SQUAREMState* sq = <SQUAREMState*>malloc(sizeof(SQUAREMState))
    if sq == NULL:
        return NULL

    sq.dimension = dimension
    sq.allocated = False

    sq.theta_0 = <double*>malloc(dimension * sizeof(double))
    sq.theta_1 = <double*>malloc(dimension * sizeof(double))
    sq.theta_2 = <double*>malloc(dimension * sizeof(double))
    sq.r_vector = <double*>malloc(dimension * sizeof(double))
    sq.v_vector = <double*>malloc(dimension * sizeof(double))
    sq.theta_extrapolated = <double*>malloc(dimension * sizeof(double))

    if (sq.theta_0 == NULL or sq.theta_1 == NULL or sq.theta_2 == NULL or
        sq.r_vector == NULL or sq.v_vector == NULL or sq.theta_extrapolated == NULL):
        free_squarem_state(sq)
        return NULL

    sq.allocated = True
    return sq


cdef void free_squarem_state(SQUAREMState* sq) noexcept nogil:
    """Free SQUAREM state."""
    if sq == NULL:
        return

    if sq.theta_0 != NULL: free(sq.theta_0)
    if sq.theta_1 != NULL: free(sq.theta_1)
    if sq.theta_2 != NULL: free(sq.theta_2)
    if sq.r_vector != NULL: free(sq.r_vector)
    if sq.v_vector != NULL: free(sq.v_vector)
    if sq.theta_extrapolated != NULL: free(sq.theta_extrapolated)

    free(sq)


cdef void copy_phi_to_array(EMState* state, double* arr) noexcept nogil:
    """Copy phi weights to array."""
    memcpy(arr, state.phi_weights, state.n_refs * sizeof(double))


cdef void copy_array_to_phi(double* arr, EMState* state) noexcept nogil:
    """Copy array to phi weights."""
    memcpy(state.phi_weights, arr, state.n_refs * sizeof(double))


cdef bint squarem_step(EMState* state, void* pool_ptr, EMConfig* config,
                       SQUAREMState* sq, double* ll_out) noexcept nogil:
    """
    SQUAREM acceleration with monotonicity safeguard.

    1. Compute theta_1 = M(theta_0), theta_2 = M(theta_1)
    2. Extrapolate: theta_sq = theta_0 - 2*alpha*r + alpha^2*v
    3. Project to feasible set (simplex)
    4. If LL(theta_sq) >= LL(theta_1): accept
       Else: fall back to theta_1

    Returns True if accelerated step accepted.
    """
    cdef uint32_t i, n = state.n_refs
    cdef double r_norm = 0.0, v_norm = 0.0
    cdef double r_temp, v_temp
    cdef double alpha
    cdef double ll_0, ll_1, ll_sq
    cdef double total

    # Save theta_0
    copy_phi_to_array(state, sq.theta_0)
    ll_0 = compute_log_likelihood(state, pool_ptr, config)

    # theta_1 = M(theta_0)
    e_step(state, pool_ptr, config)
    m_step(state, config)
    copy_phi_to_array(state, sq.theta_1)
    ll_1 = compute_log_likelihood(state, pool_ptr, config)

    # theta_2 = M(theta_1)
    e_step(state, pool_ptr, config)
    m_step(state, config)
    copy_phi_to_array(state, sq.theta_2)

    # Compute r = theta_1 - theta_0, v = (theta_2 - theta_1) - r
    for i in range(n):
        r_temp = sq.theta_1[i] - sq.theta_0[i]
        sq.r_vector[i] = r_temp
        r_norm += r_temp * r_temp

        v_temp = (sq.theta_2[i] - sq.theta_1[i]) - r_temp
        sq.v_vector[i] = v_temp
        v_norm += v_temp * v_temp

    r_norm = libc_sqrt(r_norm)
    v_norm = libc_sqrt(v_norm)

    # Compute step length alpha = -||r|| / ||v||
    if v_norm > 1e-15:
        alpha = -r_norm / v_norm
    else:
        alpha = -1.0

    # Clamp alpha
    if alpha > -0.01:
        alpha = -0.01
    elif alpha < -50.0:
        alpha = -50.0

    # Extrapolate: theta_sq = theta_0 - 2*alpha*r + alpha^2*v
    total = 0.0
    for i in range(n):
        sq.theta_extrapolated[i] = (sq.theta_0[i]
                                    - 2.0 * alpha * sq.r_vector[i]
                                    + alpha * alpha * sq.v_vector[i])
        sq.theta_extrapolated[i] = fmax(sq.theta_extrapolated[i], 1e-15)
        total += sq.theta_extrapolated[i]

    # Normalize to simplex
    if total > 0:
        for i in range(n):
            sq.theta_extrapolated[i] /= total

    # Evaluate extrapolated point
    copy_array_to_phi(sq.theta_extrapolated, state)
    normalize_phi(state)
    ll_sq = compute_log_likelihood(state, pool_ptr, config)

    # SAFEGUARD: monotonicity check
    if ll_sq >= ll_1 - 1e-10:
        ll_out[0] = ll_sq
        bf_nogil_logf_notime(b"EM_UNIFIED", "SQUAREM accepted: alpha=%.3f LL=%.6f->%.6f",
                            alpha, ll_0, ll_sq)
        return True
    else:
        # Fall back to standard EM step theta_1
        copy_array_to_phi(sq.theta_1, state)
        normalize_phi(state)
        ll_out[0] = ll_1
        bf_nogil_logf_notime(b"EM_UNIFIED", "SQUAREM rejected (LL decreased): falling back to EM step")
        return False


# =============================================================================
# Output Transform (POST-PROCESSING ONLY)
# =============================================================================

cdef void transform_to_output(EMState* state, double* output_pi,
                              double power_rho) noexcept nogil:
    """
    Transform phi-space weights to reported pi weights.

    pi_j = phi_j^rho / sum_j' phi_j'^rho

    This is ONLY for reporting, not used in EM iteration.

    Power transform interpretation (temperature analogy):
    - rho < 1: FLATTENS distribution (counteracts rich-get-richer bias)
      - e.g., rho=0.7 makes concentrated distributions more uniform
    - rho = 1: No change (standard EM output)
    - rho > 1: SHARPENS distribution (amplifies differences)

    Note: Using phi^rho (NOT phi^(1/rho)) for flattening with rho < 1.
    This is analogous to temperature T = 1/rho > 1 which smooths distributions.
    """
    cdef double total = 0.0
    cdef uint32_t j

    for j in range(state.n_refs):
        # phi^rho: rho < 1 flattens, rho > 1 sharpens
        output_pi[j] = libc_pow(state.phi_weights[j], power_rho)
        total += output_pi[j]

    if total > 0:
        for j in range(state.n_refs):
            output_pi[j] /= total


# =============================================================================
# Convergence Tracking
# =============================================================================

cdef ConvergenceState* create_convergence_state(int32_t history_length) noexcept nogil:
    """Create convergence tracking state."""
    cdef ConvergenceState* conv = <ConvergenceState*>malloc(sizeof(ConvergenceState))
    if conv == NULL:
        return NULL

    conv.history_length = history_length
    conv.ll_history = <double*>calloc(history_length, sizeof(double))
    conv.param_history = <double*>calloc(history_length, sizeof(double))

    if conv.ll_history == NULL or conv.param_history == NULL:
        free_convergence_state(conv)
        return NULL

    conv.current_index = 0
    conv.filled_count = 0
    conv.current_ll = NEG_INF
    conv.prev_ll = NEG_INF

    return conv


cdef void free_convergence_state(ConvergenceState* conv) noexcept nogil:
    """Free convergence state."""
    if conv == NULL:
        return

    if conv.ll_history != NULL:
        free(conv.ll_history)
    if conv.param_history != NULL:
        free(conv.param_history)

    free(conv)


cdef void update_convergence_history(ConvergenceState* conv,
                                     double ll_change, double param_change) noexcept nogil:
    """Update convergence history with new values."""
    if conv == NULL:
        return

    conv.ll_history[conv.current_index] = ll_change
    conv.param_history[conv.current_index] = param_change

    conv.current_index = (conv.current_index + 1) % conv.history_length
    if conv.filled_count < conv.history_length:
        conv.filled_count += 1


cdef double compute_mad(double* values, int32_t n) noexcept nogil:
    """Compute Median Absolute Deviation."""
    if n <= 0:
        return 0.0

    # For small n, just use mean absolute deviation
    cdef double mean = 0.0
    cdef double mad = 0.0
    cdef int32_t i

    for i in range(n):
        mean += values[i]
    mean /= n

    for i in range(n):
        mad += fabs(values[i] - mean)
    mad /= n

    return mad * 1.4826  # Scale factor for consistency with standard deviation


cdef bint check_convergence_mad(ConvergenceState* conv,
                                double base_tolerance) noexcept nogil:
    """
    Check convergence using MAD-based criterion.

    Converged when:
    - MAD(LL changes) < tolerance
    - AND MAD(param changes) < tolerance
    """
    if conv == NULL or conv.filled_count < 3:
        return False

    cdef double ll_mad = compute_mad(conv.ll_history, conv.filled_count)
    cdef double param_mad = compute_mad(conv.param_history, conv.filled_count)

    cdef double ll_tol = fmax(base_tolerance, ll_mad)
    cdef double param_tol = fmax(base_tolerance, param_mad)

    # Check if recent changes are within tolerance
    cdef double recent_ll_change = conv.ll_history[(conv.current_index - 1 + conv.history_length) % conv.history_length]
    cdef double recent_param_change = conv.param_history[(conv.current_index - 1 + conv.history_length) % conv.history_length]

    return (recent_ll_change <= ll_tol) and (recent_param_change <= param_tol)


cdef double compute_param_change_norm(EMState* current, EMState* prev) noexcept nogil:
    """Compute L2 norm of parameter change."""
    cdef double norm_sq = 0.0
    cdef double diff
    cdef uint32_t j

    for j in range(current.n_refs):
        diff = current.phi_weights[j] - prev.phi_weights[j]
        norm_sq += diff * diff

    return libc_sqrt(norm_sq)


# =============================================================================
# Main Entry Point
# =============================================================================

cdef int run_em(void* pool_ptr, EMConfig* config,
                double* output_weights) noexcept nogil:
    """
    Run EM algorithm.

    Parameters
    ----------
    pool_ptr : void*
        Pointer to MemoryPool
    config : EMConfig*
        Algorithm configuration
    output_weights : double*
        Output array for final weights [n_refs]

    Returns
    -------
    int
        Number of iterations run, or -1 on error
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t n_reads = pool.unique_read_count

    # Create state
    cdef EMState* state = create_em_state(n_refs, n_reads, config)
    if state == NULL:
        bf_nogil_logf_notime(b"EM_UNIFIED", "ERROR: Failed to allocate EMState")
        return -1

    cdef EMState* prev_state = create_em_state(n_refs, n_reads, config)
    if prev_state == NULL:
        free_em_state(state)
        bf_nogil_logf_notime(b"EM_UNIFIED", "ERROR: Failed to allocate prev_state")
        return -1

    # Create convergence tracker
    cdef ConvergenceState* conv = create_convergence_state(config.history_length)
    if conv == NULL:
        free_em_state(state)
        free_em_state(prev_state)
        return -1

    # Create SQUAREM state if enabled
    cdef SQUAREMState* sq = NULL
    if config.squarem_enabled:
        sq = create_squarem_state(n_refs)
        if sq == NULL:
            bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: SQUAREM allocation failed, using standard EM")

    bf_nogil_logf_notime(
        b"EM_UNIFIED",
        "Starting: refs=%u reads=%u max_iter=%d tol=%.2e hierarchical=%s unknown=%s squarem=%s rho=%.2f",
        n_refs, n_reads, config.max_iterations, config.convergence_tolerance,
        b"true" if config.hierarchical_enabled else b"false",
        b"true" if config.unknown_enabled else b"false",
        b"true" if (config.squarem_enabled and sq != NULL) else b"false",
        config.power_rho
    )

    cdef int iteration
    cdef double current_ll, prev_ll = NEG_INF
    cdef double ll_change, param_change
    cdef double ll_new
    cdef bint converged = False
    cdef bint use_squarem

    for iteration in range(config.max_iterations):
        # Save previous state
        copy_em_state(prev_state, state)
        prev_ll = conv.current_ll if conv.current_ll > NEG_INF + 1e10 else compute_log_likelihood(state, pool_ptr, config)

        # Decide whether to use SQUAREM
        use_squarem = (config.squarem_enabled and sq != NULL and
                       iteration >= config.squarem_start_iter)

        if use_squarem:
            squarem_step(state, pool_ptr, config, sq, &current_ll)
        else:
            # Standard E-step and M-step
            e_step(state, pool_ptr, config)
            m_step(state, config)
            current_ll = compute_log_likelihood(state, pool_ptr, config)

        # Compute changes
        ll_change = fabs(current_ll - prev_ll) / fmax(fabs(current_ll), 1.0)
        param_change = compute_param_change_norm(state, prev_state)

        # Update convergence history
        update_convergence_history(conv, ll_change, param_change)
        conv.current_ll = current_ll
        conv.prev_ll = prev_ll

        # Log progress
        if iteration % 10 == 0 or iteration < 5:
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "Iter %d: LL=%.6f dLL=%.2e ||dphi||=%.2e",
                iteration, current_ll, ll_change, param_change
            )

        # Check convergence
        if check_convergence_mad(conv, config.convergence_tolerance):
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "CONVERGED at iter %d: LL=%.6f dLL=%.2e ||dphi||=%.2e",
                iteration, current_ll, ll_change, param_change
            )
            converged = True
            break

    if not converged:
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "MAX_ITER reached (%d): LL=%.6f dLL=%.2e ||dphi||=%.2e",
            config.max_iterations, current_ll, ll_change, param_change
        )

    # Transform to output (apply power rho POST-PROCESSING)
    transform_to_output(state, output_weights, config.power_rho)

    bf_nogil_logf_notime(
        b"EM_UNIFIED",
        "Output transform: power_rho=%.2f (post-processing only)",
        config.power_rho
    )

    # Cleanup
    if sq != NULL:
        free_squarem_state(sq)
    free_convergence_state(conv)
    free_em_state(prev_state)
    free_em_state(state)

    return iteration + 1


# =============================================================================
# Python-accessible wrapper
# =============================================================================

def run_em_python(pool_capsule, config_dict):
    """
    Python wrapper for EM algorithm.

    Parameters
    ----------
    pool_capsule : PyCapsule
        Capsule containing MemoryPool pointer
    config_dict : dict
        Configuration dictionary with keys:
        - max_iterations: int
        - convergence_tolerance: float
        - dirichlet_prior: float
        - gamma_prior: float
        - power_rho: float
        - unknown_enabled: bool
        - unknown_margin: float
        - hierarchical_enabled: bool
        - D_avg_5p, D_avg_3p, epsilon_error: float
        - squarem_enabled: bool
        - squarem_start_iter: int
        - thread_count: int

    Returns
    -------
    tuple
        (iterations, output_weights) or (error_code, None)
    """
    import numpy as np

    cdef void* pool_ptr = PyCapsule_GetPointer(pool_capsule, <char*>NULL)
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr

    # Build config
    cdef EMConfig config
    config.max_iterations = config_dict.get('max_iterations', 50)
    config.convergence_tolerance = config_dict.get('convergence_tolerance', 1e-6)
    config.dirichlet_prior = config_dict.get('dirichlet_prior', 0.01)
    config.gamma_prior = config_dict.get('gamma_prior', 1.0)
    config.power_rho = config_dict.get('power_rho', 1.0)
    config.unknown_enabled = config_dict.get('unknown_enabled', False)
    config.unknown_margin = config_dict.get('unknown_margin', 2.0)
    config.hierarchical_enabled = config_dict.get('hierarchical_enabled', False)
    config.D_avg_5p = config_dict.get('D_avg_5p', 0.02)
    config.D_avg_3p = config_dict.get('D_avg_3p', 0.02)
    config.epsilon_error = config_dict.get('epsilon_error', 0.01)
    config.squarem_enabled = config_dict.get('squarem_enabled', True)
    config.squarem_start_iter = config_dict.get('squarem_start_iter', 3)
    config.enable_globalization = config_dict.get('enable_globalization', True)
    config.backtrack_factor = config_dict.get('backtrack_factor', 0.5)
    config.max_backtrack_steps = config_dict.get('max_backtrack_steps', 5)
    config.thread_count = config_dict.get('thread_count', 1)
    config.history_length = config_dict.get('history_length', 5)

    # Allocate output
    cdef uint32_t n_refs = pool.reference_count
    output_weights = np.zeros(n_refs, dtype=np.float64)
    cdef double[::1] output_view = output_weights

    # Run EM
    cdef int result
    with nogil:
        result = run_em(pool_ptr, &config, &output_view[0])

    if result < 0:
        return (result, None)

    return (result, output_weights)


def execute_em_py(
    uintptr_t pool_ptr,
    int max_iterations,
    double convergence_tolerance,
    double dirichlet_prior,
    double power_rho,
    bint unknown_enabled,
    double unknown_margin,
    bint hierarchical_enabled,
    double D_avg_5p,
    double D_avg_3p,
    double epsilon_error,
    bint squarem_enabled,
    int squarem_start_iter,
    bint enable_globalization,
    double backtrack_factor,
    int max_backtrack_steps,
    int thread_count,
):
    """
    Python wrapper for EM algorithm.

    Updates the memory pool's reference weights in place.

    Parameters
    ----------
    pool_ptr : uintptr_t
        Pointer to MemoryPool
    max_iterations : int
        Maximum EM iterations
    convergence_tolerance : double
        Convergence threshold
    dirichlet_prior : double
        Dirichlet prior alpha
    power_rho : double
        Power transform exponent (1.0 = no transform, <1 = flattening)
    unknown_enabled : bool
        Enable unknown component
    unknown_margin : double
        Margin below best score for unknown
    hierarchical_enabled : bool
        Enable hierarchical ancient/modern
    D_avg_5p, D_avg_3p : double
        Average PMD damage rates
    epsilon_error : double
        Sequencing error rate
    squarem_enabled : bool
        Enable SQUAREM acceleration
    squarem_start_iter : int
        Start SQUAREM after this iteration
    enable_globalization : bool
        Enable backtracking
    backtrack_factor : double
        Backtracking step reduction
    max_backtrack_steps : int
        Maximum backtracking attempts
    thread_count : int
        Number of threads

    Returns
    -------
    int
        0 on success, -1 on failure
    """
    import numpy as np

    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t j

    # Build config
    cdef EMConfig config
    config.max_iterations = max_iterations
    config.convergence_tolerance = convergence_tolerance
    config.dirichlet_prior = dirichlet_prior
    config.gamma_prior = 1.0  # Default
    config.power_rho = power_rho
    config.unknown_enabled = unknown_enabled
    config.unknown_margin = unknown_margin
    config.hierarchical_enabled = hierarchical_enabled
    config.D_avg_5p = D_avg_5p
    config.D_avg_3p = D_avg_3p
    config.epsilon_error = epsilon_error
    config.squarem_enabled = squarem_enabled
    config.squarem_start_iter = squarem_start_iter
    config.enable_globalization = enable_globalization
    config.backtrack_factor = backtrack_factor
    config.max_backtrack_steps = max_backtrack_steps
    config.thread_count = thread_count
    config.history_length = 5  # Default

    # Allocate output buffer
    output_weights = np.zeros(n_refs, dtype=np.float64)
    cdef double[::1] output_view = output_weights

    # Run EM
    cdef int iterations
    cdef void* pool_void = <void*>pool_ptr
    cdef double* ref_weights
    with nogil:
        iterations = run_em(pool_void, &config, &output_view[0])

    if iterations < 0:
        return -1

    # Copy output weights back to pool's reference weights
    # This matches what execute_em_algorithm does
    ref_weights = get_reference_weights(pool)
    for j in range(n_refs):
        ref_weights[j] = output_view[j]

    # Update pool iteration count and convergence flag
    pool.iteration_count = iterations
    pool.algorithm_converged = (iterations < max_iterations)

    return 0
