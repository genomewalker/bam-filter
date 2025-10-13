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

"""Expectation-Maximization algorithm for read reassignment.

Implements the EM algorithm with SQUAREM acceleration (Varadhan & Roland 2008)
for reassigning multi-mapping reads to references based on alignment quality.
"""

from cython.parallel import prange, threadid

from libc.math cimport exp, fabs, fmax, fmin, log, log2, sqrt as libc_sqrt, pow as libc_pow, INFINITY
from libc.stdint cimport int32_t, int64_t, uint16_t, uint32_t, uint64_t, uint8_t, uintptr_t
from libc.stdio cimport sprintf
from libc.stdlib cimport calloc, free, malloc, realloc
from libc.string cimport memcpy, memmove, memset, strlen
from libc.time cimport clock, clock_t, CLOCKS_PER_SEC
from libc.float cimport DBL_EPSILON

from bam_filter.processor cimport (
    MemoryPool, Alignment, AlignmentScoringConfig,
    min_int32, max_int32, min_int64, max_int64, min_double, max_double
)
from bam_filter.processor_types cimport EMAlgorithmConfig
from bam_filter.processor_precomputed cimport (
    PrecomputedWeights, create_precomputed_weights, 
    free_precomputed_weights, update_precomputed_weights
)
from bam_filter.processor_fast_math cimport stable_log_sum_exp, safe_normalize_weights
from bam_filter.processor_convergence_helpers cimport _median5, _count_filled, _robust_sigma

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil

cdef inline double* get_reference_weights(MemoryPool* pool) noexcept nogil:
    """Get pointer to reference weights array in unified buffer."""
    return pool.unified_buffer + pool.reference_weights_offset

cdef inline double* get_temp_buffer_A(MemoryPool* pool) noexcept nogil:
    """Get pointer to temporary buffer A in unified buffer."""
    return pool.unified_buffer + pool.temp_buffer_A_offset

cdef inline double* get_temp_buffer_B(MemoryPool* pool) noexcept nogil:
    """Get pointer to temporary buffer B in unified buffer."""
    return pool.unified_buffer + pool.temp_buffer_B_offset

cdef inline void PREFETCH_READ(void* ptr) noexcept nogil:
    """Prefetch hint for read access (no-op, reserved for future optimization)."""
    pass

cdef inline void PREFETCH_WRITE(void* ptr) noexcept nogil:
    """Prefetch hint for write access (no-op, reserved for future optimization)."""
    pass

cdef double _dominance_last_base_strength = -1.0
cdef double _dominance_last_final_strength = -1.0
cdef double _dominance_last_concentration = -1.0
cdef double _dominance_last_entropy_threshold = -1.0
cdef double _dominance_last_auto_entropy = -1.0
cdef double _dominance_last_auto_density = -1.0
cdef double _dominance_last_auto_ref_factor = -1.0
cdef double _dominance_last_auto_clamped = -1.0
cdef bint _dominance_logged_mode = False
cdef bint _dominance_logged_adaptive = False
cdef bint _dominance_logged_skip = False
cdef bint _dominance_last_mode_manual = False
cdef int _dominance_last_iteration = -1
cdef int _dominance_iteration_context = -1

cdef double _dominance_iter_min_base = INFINITY
cdef double _dominance_iter_max_base = -INFINITY
cdef double _dominance_iter_last_base = -1.0
cdef double _dominance_iter_last_concentration = -1.0
cdef double _dominance_iter_last_confidence = -1.0
cdef double _dominance_iter_min_strength = INFINITY
cdef double _dominance_iter_max_strength = -INFINITY
cdef double _dominance_iter_last_strength = -1.0
cdef int _dominance_iter_update_count = 0
cdef bint _dominance_iter_adaptive = False


cdef inline void _dominance_reset_iteration_metrics() noexcept nogil:
    global _dominance_iter_min_base, _dominance_iter_max_base, _dominance_iter_last_base
    global _dominance_iter_last_concentration, _dominance_iter_last_confidence
    global _dominance_iter_min_strength, _dominance_iter_max_strength, _dominance_iter_last_strength
    global _dominance_iter_update_count, _dominance_iter_adaptive

    _dominance_iter_min_base = INFINITY
    _dominance_iter_max_base = -INFINITY
    _dominance_iter_last_base = -1.0
    _dominance_iter_last_concentration = -1.0
    _dominance_iter_last_confidence = -1.0
    _dominance_iter_min_strength = INFINITY
    _dominance_iter_max_strength = -INFINITY
    _dominance_iter_last_strength = -1.0
    _dominance_iter_update_count = 0
    _dominance_iter_adaptive = False


cdef inline void _dominance_record_iteration_metrics(double base_strength,
                                                     double concentration,
                                                     double confidence_factor,
                                                     double final_strength,
                                                     bint adaptive) noexcept nogil:
    global _dominance_iter_min_base, _dominance_iter_max_base, _dominance_iter_last_base
    global _dominance_iter_last_concentration, _dominance_iter_last_confidence
    global _dominance_iter_min_strength, _dominance_iter_max_strength, _dominance_iter_last_strength
    global _dominance_iter_update_count, _dominance_iter_adaptive

    _dominance_iter_adaptive = adaptive
    _dominance_iter_last_base = base_strength
    if base_strength < _dominance_iter_min_base:
        _dominance_iter_min_base = base_strength
    if base_strength > _dominance_iter_max_base:
        _dominance_iter_max_base = base_strength

    _dominance_iter_last_concentration = concentration
    _dominance_iter_last_confidence = confidence_factor

    if final_strength < _dominance_iter_min_strength:
        _dominance_iter_min_strength = final_strength
    if final_strength > _dominance_iter_max_strength:
        _dominance_iter_max_strength = final_strength
    _dominance_iter_last_strength = final_strength
    _dominance_iter_update_count += 1


cdef inline void _dominance_flush_iteration_metrics(int iteration) noexcept nogil:
    global _dominance_iter_min_base, _dominance_iter_max_base, _dominance_iter_last_base
    global _dominance_iter_last_concentration, _dominance_iter_last_confidence
    global _dominance_iter_min_strength, _dominance_iter_max_strength, _dominance_iter_last_strength
    global _dominance_iter_update_count, _dominance_iter_adaptive

    if iteration < 0 or _dominance_iter_update_count == 0:
        return

    cdef double base_span = fabs(_dominance_iter_max_base - _dominance_iter_min_base)
    cdef double strength_span = fabs(_dominance_iter_max_strength - _dominance_iter_min_strength)

    if _dominance_iter_adaptive:
        if strength_span < 5e-3 and base_span < 5e-3:
            bf_nogil_logf_notime(
                b"EM",
                "dominance_regularization: iteration=%d base=%.3f concentration=%.3f final_strength=%.3f updates=%d",
                iteration,
                _dominance_iter_last_base,
                _dominance_iter_last_concentration,
                _dominance_iter_last_strength,
                _dominance_iter_update_count,
            )
        else:
            bf_nogil_logf_notime(
                b"EM",
                "dominance_regularization: iteration=%d base_range=[%.3f, %.3f] concentration=%.3f final_strength_range=[%.3f, %.3f] updates=%d",
                iteration,
                _dominance_iter_min_base if _dominance_iter_min_base != INFINITY else _dominance_iter_last_base,
                _dominance_iter_max_base if _dominance_iter_max_base != -INFINITY else _dominance_iter_last_base,
                _dominance_iter_last_concentration,
                _dominance_iter_min_strength,
                _dominance_iter_max_strength,
                _dominance_iter_update_count,
            )
    else:
        if strength_span < 5e-3:
            bf_nogil_logf_notime(
                b"EM",
                "dominance_regularization: iteration=%d adaptive=false base=%.3f final_strength=%.3f updates=%d",
                iteration,
                _dominance_iter_last_base,
                _dominance_iter_last_strength,
                _dominance_iter_update_count,
            )
        else:
            bf_nogil_logf_notime(
                b"EM",
                "dominance_regularization: iteration=%d adaptive=false base_range=[%.3f, %.3f] final_strength_range=[%.3f, %.3f] updates=%d",
                iteration,
                _dominance_iter_min_base if _dominance_iter_min_base != INFINITY else _dominance_iter_last_base,
                _dominance_iter_max_base if _dominance_iter_max_base != -INFINITY else _dominance_iter_last_base,
                _dominance_iter_min_strength,
                _dominance_iter_max_strength,
                _dominance_iter_update_count,
            )

    _dominance_reset_iteration_metrics()

cdef double calculate_dataset_entropy(MemoryPool* pool) except -1.0 nogil:
    """Calculate normalized Shannon entropy of reference weight distribution.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing reference weights

    Returns
    -------
    double
        Normalized entropy [0, 1] where 0 is completely concentrated, 1 is uniform
    """
    cdef double entropy = 0.0
    cdef double* weights = get_reference_weights(pool)
    cdef uint32_t i
    
    for i in range(pool.reference_count):
        if weights[i] > 1e-15:
            entropy -= weights[i] * log(weights[i])
    
    cdef double max_entropy = log(<double>pool.reference_count)
    return entropy / max_entropy if max_entropy > 0.0 else 0.0

cdef double auto_tune_dominance_strength(MemoryPool* pool) except -1.0 nogil:
    """Automatically determine dominance regularization strength from dataset properties.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing reference weights and alignments

    Returns
    -------
    double
        Recommended dominance strength clamped to [0.1, 2.5]
    """
    global _dominance_last_auto_entropy, _dominance_last_auto_density
    global _dominance_last_auto_ref_factor, _dominance_last_auto_clamped
    cdef double entropy = calculate_dataset_entropy(pool)
    cdef double concentration = 1.0 - entropy
    
    cdef double avg_alignments_per_read = <double>pool.alignment_count / <double>pool.final_unique_reads
    cdef double alignment_density = fmin(avg_alignments_per_read / 20.0, 1.0)
    
    cdef double ref_count_factor = fmax(0.1, fmin(1.0, 50000.0 / <double>pool.reference_count))
    
    cdef double auto_strength = 0.2 + 1.0 * concentration + 0.3 * alignment_density
    auto_strength *= ref_count_factor
    
    cdef double final_strength = fmax(0.1, fmin(auto_strength, 2.5))
    
    cdef double tolerance = 1e-9
    if (not _dominance_logged_mode or
        fabs(entropy - _dominance_last_auto_entropy) > tolerance or
        fabs(alignment_density - _dominance_last_auto_density) > tolerance or
        fabs(ref_count_factor - _dominance_last_auto_ref_factor) > tolerance or
        fabs(final_strength - _dominance_last_auto_clamped) > tolerance):
        bf_nogil_logf_notime(
            b"EM",
            "dominance_regularization: auto_summary entropy=%.3f concentration=%.3f align_density=%.3f ref_factor=%.3f base=%.3f clamped=%.3f",
            entropy,
            concentration,
            alignment_density,
            ref_count_factor,
            auto_strength,
            final_strength,
        )
        _dominance_last_auto_entropy = entropy
        _dominance_last_auto_density = alignment_density
        _dominance_last_auto_ref_factor = ref_count_factor
        _dominance_last_auto_clamped = final_strength
    
    return final_strength

cdef double learn_entropy_threshold_from_data(MemoryPool* pool, EMAlgorithmConfig* config) except -1.0 nogil:
    """Learn appropriate entropy threshold from dataset characteristics.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool with alignment data
    config : EMAlgorithmConfig*
        Configuration (returns manual threshold if set)

    Returns
    -------
    double
        Learned or configured entropy threshold
    """
    if config.entropy_confidence_threshold > 0.0:
        return config.entropy_confidence_threshold
    
    cdef double n_refs = <double>pool.reference_count
    cdef double alignment_density = <double>pool.alignment_count / <double>pool.final_unique_reads
    
    cdef double base_threshold = 0.1
    cdef double density_adjustment = fmin(0.15, alignment_density / 30.0) 
    bf_nogil_logf_notime(
        b"EM",
        "dominance_regularization: learned_entropy_threshold=%.3f",
        base_threshold + density_adjustment,
    )
    return base_threshold + density_adjustment

cdef double calculate_adaptive_dominance_strength(MemoryPool* pool, EMAlgorithmConfig* config) except -1.0 nogil:
    global _dominance_iteration_context, _dominance_last_iteration
    global _dominance_logged_mode, _dominance_logged_adaptive, _dominance_logged_skip
    global _dominance_last_concentration, _dominance_last_entropy_threshold
    global _dominance_last_base_strength, _dominance_last_mode_manual
    global _dominance_last_final_strength
    cdef double base_strength
    cdef double dataset_entropy = calculate_dataset_entropy(pool)
    cdef double max_entropy = log(<double>pool.reference_count)
    cdef double final_strength_simple
    cdef double concentration = 1.0 - (dataset_entropy / max_entropy)
    cdef bint mode_manual
    cdef double adaptive_strength
    cdef double confidence_factor
    cdef double final_strength
    cdef double tolerance = 1e-9

    if _dominance_iteration_context != _dominance_last_iteration:
        _dominance_flush_iteration_metrics(_dominance_last_iteration)
        _dominance_logged_skip = False
        _dominance_last_iteration = _dominance_iteration_context

    if concentration < config.entropy_confidence_threshold:
        if (not _dominance_logged_skip or
            fabs(concentration - _dominance_last_concentration) > tolerance or
            fabs(config.entropy_confidence_threshold - _dominance_last_entropy_threshold) > tolerance):
            bf_nogil_logf_notime(
                b"EM",
                "dominance_regularization: skipped concentration=%.3f threshold=%.3f",
                concentration,
                config.entropy_confidence_threshold,
            )
            _dominance_logged_skip = True
            _dominance_last_concentration = concentration
            _dominance_last_entropy_threshold = config.entropy_confidence_threshold
        return 0.0
    
    _dominance_logged_skip = False
    
    if config.dominance_strength <= 0.0:
        base_strength = auto_tune_dominance_strength(pool)
        mode_manual = False
    else:
        base_strength = config.dominance_strength
        mode_manual = True

    if (not _dominance_logged_mode or
        fabs(base_strength - _dominance_last_base_strength) > tolerance or
        mode_manual != _dominance_last_mode_manual):
        bf_nogil_logf_notime(
            b"EM",
            "dominance_regularization: mode=%s base_strength=%.3f",
            b"manual" if mode_manual else b"auto",
            base_strength,
        )
        _dominance_logged_mode = True
        _dominance_last_base_strength = base_strength
        _dominance_last_mode_manual = mode_manual
    
    if not config.use_adaptive_dominance:
        final_strength_simple = fmax(config.min_penalty_strength,
                                     fmin(base_strength, config.max_penalty_strength))
        _dominance_record_iteration_metrics(
            base_strength,
            concentration,
            concentration,
            final_strength_simple,
            False,
        )
        return final_strength_simple
    
    adaptive_strength = base_strength * (1.0 + config.entropy_scaling_factor * concentration)
    confidence_factor = concentration
    final_strength = adaptive_strength * confidence_factor
    final_strength = fmax(config.min_penalty_strength, fmin(final_strength, config.max_penalty_strength))
    
    _dominance_record_iteration_metrics(
        base_strength,
        concentration,
        confidence_factor,
        final_strength,
        True,
    )
    
    return final_strength

cdef inline double compute_dominance_penalty(double pi_j, double penalty_strength) except -1.0 nogil:
    return exp(-penalty_strength * pi_j)

cdef void execute_em_expectation_step_vectorized(MemoryPool* pool,
                                                int32_t thread_count,
                                                PrecomputedWeights* precomp,
                                                EMAlgorithmConfig* config=NULL,
                                                int32_t current_iteration=0) noexcept nogil:
    """Execute E-step: calculate posterior probabilities and accumulate expected counts.

    Computes posterior probabilities for each alignment and accumulates expected
    reference counts. Uses parallel processing with thread-local accumulators.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    thread_count : int32_t
        Number of threads for parallel processing
    precomp : PrecomputedWeights*
        Precomputed log weights for efficiency
    config : EMAlgorithmConfig*
        Optional configuration (unused, retained for API compatibility)
    current_iteration : int32_t
        Current iteration number (unused, retained for API compatibility)
    """
    cdef int32_t read_idx
    cdef uint32_t ref_idx
    cdef uint64_t start_pos, end_pos, alignment_idx
    cdef uint32_t alignment_count
    cdef float alignment_score
    cdef double log_likelihood
    cdef double log_weighted_score, log_normalizer, posterior
    cdef double NEG_INF = -INFINITY
    cdef double uniform_weight
    cdef int thread_id, t, r
    cdef uint32_t i
    cdef int32_t ti
    cdef int max_threads = thread_count
    cdef double** thread_local_weights = NULL
    cdef int* thread_used = NULL
    cdef int32_t read_idx_c
    cdef int32_t unique_read_count_i = <int32_t>pool.unique_read_count

    cdef double* reference_weights = get_reference_weights(pool)
    cdef double* new_weights = get_temp_buffer_A(pool)

    update_precomputed_weights(precomp, reference_weights)

    memset(new_weights, 0, pool.reference_count * sizeof(double))

    thread_local_weights = <double**>calloc(max_threads, sizeof(double*))
    thread_used = <int*>calloc(max_threads, sizeof(int))

    if not thread_local_weights or not thread_used:
        for read_idx in range(unique_read_count_i):
            alignment_count = pool.read_alignment_counts[read_idx]
            if alignment_count == 0:
                continue

            start_pos = pool.read_alignment_starts[read_idx]
            end_pos = start_pos + alignment_count

            if alignment_count == 1:
                ref_idx = pool.alignments[start_pos].reference_index
                if ref_idx < pool.reference_count:
                    new_weights[ref_idx] += 1.0
                continue

            log_normalizer = NEG_INF

            for alignment_idx in range(start_pos, end_pos):
                ref_idx = pool.alignments[alignment_idx].reference_index
                if ref_idx < pool.reference_count:
                    alignment_score = pool.alignments[alignment_idx].alignment_score
                    log_likelihood = <double>alignment_score
                    log_weighted_score = precomp.log_weights[ref_idx] + log_likelihood
                    log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted_score)

            if log_normalizer > NEG_INF + 1e10:
                for alignment_idx in range(start_pos, end_pos):
                    ref_idx = pool.alignments[alignment_idx].reference_index
                    if ref_idx < pool.reference_count:
                        alignment_score = pool.alignments[alignment_idx].alignment_score
                        log_likelihood = <double>alignment_score
                        log_weighted_score = precomp.log_weights[ref_idx] + log_likelihood
                        posterior = exp(log_weighted_score - log_normalizer)
                        posterior = fmax(fmin(posterior, 0.999), 1e-12)
                        new_weights[ref_idx] += posterior
            else:
                uniform_weight = 1.0 / alignment_count if alignment_count > 0 else 0.0
                for alignment_idx in range(start_pos, end_pos):
                    ref_idx = pool.alignments[alignment_idx].reference_index
                    if ref_idx < pool.reference_count:
                        new_weights[ref_idx] += uniform_weight
        return

    for thread_id in range(max_threads):
        thread_local_weights[thread_id] = <double*>calloc(pool.reference_count, sizeof(double))
        if not thread_local_weights[thread_id]:
            for t in range(thread_id):
                if thread_local_weights[t]:
                    free(thread_local_weights[t])
            free(thread_local_weights)
            free(thread_used)
            execute_em_expectation_step_vectorized(pool, 1, precomp, config, current_iteration)
            return

    cdef uint32_t ref_count = pool.reference_count
    cdef double* log_weights_ptr = precomp.log_weights  # Standard log weights

    for read_idx_c in prange(unique_read_count_i, nogil=True, schedule='static', num_threads=thread_count):
        if read_idx_c + 1 < unique_read_count_i:
            PREFETCH_READ(&pool.read_alignment_starts[read_idx_c + 1])
            PREFETCH_READ(&pool.read_alignment_counts[read_idx_c + 1])
        thread_id = threadid()
        if thread_id >= max_threads:
            continue

        thread_used[thread_id] = 1
        read_idx = read_idx_c

        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        if alignment_count == 1:
            ref_idx = pool.alignments[start_pos].reference_index
            if ref_idx < ref_count:
                thread_local_weights[thread_id][ref_idx] += 1.0
            continue

        log_normalizer = NEG_INF

        for alignment_idx in range(start_pos, end_pos):
            if alignment_idx + 1 < end_pos:
                PREFETCH_READ(&pool.alignments[alignment_idx + 1].alignment_score)
            ref_idx = pool.alignments[alignment_idx].reference_index
            if ref_idx < ref_count:
                alignment_score = pool.alignments[alignment_idx].alignment_score
                log_likelihood = <double>alignment_score
                log_weighted_score = log_weights_ptr[ref_idx] + log_likelihood
                log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted_score)

        if log_normalizer > NEG_INF + 1e10:
            for alignment_idx in range(start_pos, end_pos):
                if alignment_idx + 1 < end_pos:
                    PREFETCH_READ(&pool.alignments[alignment_idx + 1].alignment_score)
                ref_idx = pool.alignments[alignment_idx].reference_index
                if ref_idx < ref_count:
                    alignment_score = pool.alignments[alignment_idx].alignment_score
                    log_likelihood = <double>alignment_score
                    log_weighted_score = log_weights_ptr[ref_idx] + log_likelihood
                    posterior = exp(log_weighted_score - log_normalizer)
                    posterior = fmax(fmin(posterior, 0.999), 1e-12)
                    thread_local_weights[thread_id][ref_idx] += posterior
        else:
            uniform_weight = 1.0 / alignment_count if alignment_count > 0 else 0.0
            for alignment_idx in range(start_pos, end_pos):
                ref_idx = pool.alignments[alignment_idx].reference_index
                if ref_idx < ref_count:
                    thread_local_weights[thread_id][ref_idx] += uniform_weight

    cdef uint32_t** thread_touched_ref_idx = <uint32_t**>calloc(max_threads, sizeof(uint32_t*))
    cdef int* thread_touched_count = <int*>calloc(max_threads, sizeof(int))
    for thread_id in range(max_threads):
        thread_touched_ref_idx[thread_id] = <uint32_t*>calloc(ref_count, sizeof(uint32_t))
        thread_touched_count[thread_id] = 0

    for thread_id in range(max_threads):
        if thread_used[thread_id]:
            for ref_idx in range(ref_count):
                if thread_local_weights[thread_id][ref_idx] != 0.0:
                    thread_touched_ref_idx[thread_id][thread_touched_count[thread_id]] = ref_idx;
                    thread_touched_count[thread_id] += 1

    for thread_id in range(max_threads):
        if thread_used[thread_id]:
            for ti in range(thread_touched_count[thread_id]):
                ref_idx = thread_touched_ref_idx[thread_id][ti]
                new_weights[ref_idx] += thread_local_weights[thread_id][ref_idx]

    for thread_id in range(max_threads):
        if thread_local_weights[thread_id]:
            free(thread_local_weights[thread_id])
        if thread_touched_ref_idx[thread_id]:
            free(thread_touched_ref_idx[thread_id])
    free(thread_local_weights)
    free(thread_touched_ref_idx)
    free(thread_touched_count)
    free(thread_used)

cdef void execute_em_maximization_step_vectorized(MemoryPool* pool,
                                                 double alpha_prior,
                                                 PrecomputedWeights* precomp,
                                                 EMAlgorithmConfig* config=NULL,
                                                 int current_iteration=0,
                                                 double current_ll=-1e20,
                                                 double prev_ll=-1e20) noexcept nogil:
    """Execute M-step: update reference weights from expected counts.

    Normalizes accumulated counts with optional Dirichlet prior and dominance
    regularization. Detects pathological convergence and activates emergency
    regularization if enabled.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing weights
    alpha_prior : double
        Dirichlet prior strength for regularization
    precomp : PrecomputedWeights*
        Precomputed weights structure (marked dirty after update)
    config : EMAlgorithmConfig*
        Optional configuration for regularization
    current_iteration : int
        Current iteration number
    current_ll : double
        Current log-likelihood
    prev_ll : double
        Previous log-likelihood
    """
    global _dominance_iteration_context
    cdef double* reference_weights = get_reference_weights(pool)
    cdef double* accumulated_weights = get_temp_buffer_A(pool)
    cdef double total_sum = 0.0
    cdef double effective_prior = fmax(alpha_prior, 0.01)
    cdef uint32_t i
    cdef double inv_total
    cdef double adaptive_strength
    cdef double current_weight
    cdef double penalty
    
    if (config and config.enable_emergency_regularization and 
        not config.dominance_regularization_active and current_iteration > 0):
        
        if detect_pathological_convergence(pool, config, current_iteration, current_ll, prev_ll):
            bf_nogil_logf_notime(b"EM", "EMERGENCY ACTIVATION: Dominance regularization enabled at iteration %d", 
                   current_iteration)
            config.dominance_regularization_active = True
            config.dominance_strength = 2.0
            config.entropy_scaling_factor = 2.0
            config.min_penalty_strength = 0.5
    
    for i in range(pool.reference_count):
        reference_weights[i] = accumulated_weights[i] + effective_prior
        total_sum += reference_weights[i]
    
    if total_sum > 1e-15:
        inv_total = 1.0 / total_sum
        for i in range(pool.reference_count):
            reference_weights[i] *= inv_total
    else:
        inv_total = 1.0 / pool.reference_count
        for i in range(pool.reference_count):
            reference_weights[i] = inv_total
        precomp.weights_dirty = True
        return
    
    if config and (config.enable_dominance_regularization or config.dominance_regularization_active):
        _dominance_iteration_context = current_iteration
        adaptive_strength = calculate_adaptive_dominance_strength(pool, config)
        
        if adaptive_strength > 0.0:
            total_sum = 0.0
            for i in range(pool.reference_count):
                current_weight = reference_weights[i]
                penalty = exp(-adaptive_strength * current_weight)
                reference_weights[i] *= penalty
                total_sum += reference_weights[i]
            
            if total_sum > 1e-15:
                inv_total = 1.0 / total_sum
                for i in range(pool.reference_count):
                    reference_weights[i] *= inv_total
            
    for i in range(pool.reference_count):
        reference_weights[i] = fmax(fmin(reference_weights[i], 0.999), 1e-12)
    
    precomp.weights_dirty = True

cdef double compute_log_likelihood_vectorized(MemoryPool* pool,
                                            PrecomputedWeights* precomp) except -1.0 nogil:
    """Calculate total log-likelihood of current weight assignment.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    precomp : PrecomputedWeights*
        Precomputed log weights

    Returns
    -------
    double
        Mean log-likelihood per read
    """
    cdef uint32_t read_idx, ref_idx
    cdef uint64_t start_pos, end_pos, alignment_idx
    cdef uint32_t alignment_count
    cdef float alignment_score
    cdef double log_likelihood
    cdef double log_weighted_score, log_normalizer
    cdef double total_log_likelihood = 0.0
    cdef double NEG_INF = -1e20
    cdef int32_t valid_reads = 0

    cdef double* reference_weights = get_reference_weights(pool)
    cdef double read_likelihood = NEG_INF

    update_precomputed_weights(precomp, reference_weights)

    for read_idx in range(pool.unique_read_count):
        PREFETCH_READ(&pool.read_alignment_starts[read_idx])
        PREFETCH_READ(&pool.read_alignment_counts[read_idx])
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        if alignment_count == 1:
            alignment_idx = start_pos
            ref_idx = pool.alignments[alignment_idx].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[alignment_idx].alignment_score                
                log_likelihood = <double>alignment_score
                log_weighted_score = precomp.log_weights[ref_idx] + log_likelihood
                read_likelihood = log_weighted_score
            
        else:
            log_normalizer = NEG_INF
            for alignment_idx in range(start_pos, end_pos):
                PREFETCH_READ(&pool.alignments[alignment_idx].alignment_score)
                ref_idx = pool.alignments[alignment_idx].reference_index
                if ref_idx < pool.reference_count:
                    alignment_score = pool.alignments[alignment_idx].alignment_score
                    
                    log_likelihood = <double>alignment_score
                    log_weighted_score = precomp.log_weights[ref_idx] + log_likelihood
                    log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted_score)

            read_likelihood = log_normalizer

        if read_likelihood > NEG_INF + 1e10:
            total_log_likelihood += read_likelihood
            valid_reads += 1

    return total_log_likelihood / valid_reads if valid_reads > 0 else NEG_INF

cdef bint detect_pathological_convergence(MemoryPool* pool, 
                                         EMAlgorithmConfig* config,
                                         int current_iteration, 
                                         double current_ll, 
                                         double prev_ll) noexcept nogil:

    cdef double* weights = get_reference_weights(pool)
    cdef double entropy = calculate_dataset_entropy(pool)
    cdef double max_weight = 0.0
    cdef uint32_t dominant_refs = 0
    cdef uint32_t i
    cdef double dominant_threshold
    
    if entropy < config.emergency_entropy_threshold:
        bf_nogil_logf_notime(
            b"EM",
            "emergency_trigger: reason=low_entropy entropy=%.3f threshold=%.3f",
            entropy,
            config.emergency_entropy_threshold,
        )
        return True
    
    dominant_threshold = fmax(0.01, 1.0 / (10.0 * libc_sqrt(<double>pool.reference_count)))
    
    for i in range(pool.reference_count):
        if weights[i] > max_weight:
            max_weight = weights[i]
        if weights[i] > dominant_threshold:
            dominant_refs += 1
    
    if max_weight > config.emergency_max_weight_threshold:
        bf_nogil_logf_notime(
            b"EM",
            "emergency_trigger: reason=dominant_weight max_weight=%.3f threshold=%.3f",
            max_weight,
            config.emergency_max_weight_threshold,
        )
        return True
    
    if pool.reference_count > 1000:
        if dominant_refs < config.emergency_min_dominant_refs:
            bf_nogil_logf_notime(
                b"EM",
                "emergency_trigger: reason=insufficient_dominant_refs detected=%u total=%u min_expected=%u threshold=%.4f",
                dominant_refs,
                pool.reference_count,
                config.emergency_min_dominant_refs,
                dominant_threshold,
            )
            return True
    
    cdef int stabilization_period = 3
    if (current_iteration > stabilization_period and 
        current_ll < prev_ll - config.emergency_likelihood_drop):
        bf_nogil_logf_notime(
            b"EM",
            "emergency_trigger: reason=likelihood_drop delta=%.6f threshold=%.6f stabilization_period=%d",
            prev_ll - current_ll,
            config.emergency_likelihood_drop,
            stabilization_period,
        )
        return True
    
    return False

cdef int execute_em_algorithm(MemoryPool* pool, EMAlgorithmConfig* config) except -1 nogil:
    """Execute EM algorithm with SQUAREM acceleration for read reassignment.

    Implements the Expectation-Maximization algorithm for reassigning multi-mapping
    reads to references. Includes SQUAREM acceleration (Varadhan & Roland 2008) with
    three steplength schemes (S1, S2, S3) and globalization with backtracking.

    Algorithm:
    1. E-step: Calculate expected read assignments based on current weights
    2. M-step: Update reference weights based on expected assignments
    3. SQUAREM: Accelerate convergence using extrapolation when enabled
    4. Globalization: Backtrack if extrapolation decreases likelihood

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments and weight arrays
    config : EMAlgorithmConfig*
        Algorithm configuration (iterations, tolerance, SQUAREM settings)

    Returns
    -------
    int
        0 on success, negative error code on failure

    Notes
    -----
    Reference: Varadhan, R. and Roland, C. (2008). Simple and Globally Convergent
    Methods for Accelerating the Convergence of Any EM Algorithm. Scandinavian
    Journal of Statistics, 35: 335-353.
    """
    cdef int32_t iteration = 0
    cdef double current_log_likelihood = -1e20
    cdef double prev_log_likelihood = -1e20
    cdef bint use_squarem = config.use_squarem_acceleration
    cdef double tolerance = fmax(config.convergence_tolerance, 1e-8)
    cdef int32_t max_iterations = config.maximum_iterations
    cdef double* ll_history = NULL
    cdef double* param_history = NULL
    cdef double* prev_weights = NULL
    cdef double alpha = -999.0
    cdef double* squarem_memory = NULL
    cdef double* theta_0 = NULL
    cdef double* theta_1 = NULL
    cdef double* theta_2 = NULL
    cdef double* r_vector = NULL
    cdef double* v_vector = NULL
    cdef double* theta_extrapolated = NULL
    cdef double r_norm, v_norm, r_temp, v_temp
    cdef uint32_t i
    cdef uint32_t ref_cnt
    cdef double r_temp0, r_temp1, r_temp2, r_temp3, r_temp4, r_temp5, r_temp6, r_temp7
    cdef double v_temp0, v_temp1, v_temp2, v_temp3, v_temp4, v_temp5, v_temp6, v_temp7
    cdef size_t array_size, total_size
    cdef double test_ll = 0.0
    cdef int backtrack_steps = 0
    cdef double* reference_weights = get_reference_weights(pool)
    cdef PrecomputedWeights* precomp = NULL
    cdef double backtrack_factor = 0.5
    cdef int total_backtrack_steps = 0
    cdef double ll_before_backtrack, ll_after_backtrack
    global _dominance_last_base_strength, _dominance_last_final_strength
    global _dominance_last_concentration, _dominance_last_entropy_threshold
    global _dominance_last_auto_entropy, _dominance_last_auto_density
    global _dominance_last_auto_ref_factor, _dominance_last_auto_clamped
    global _dominance_logged_mode, _dominance_logged_adaptive, _dominance_logged_skip
    global _dominance_last_mode_manual, _dominance_last_iteration
    global _dominance_iteration_context

    _dominance_last_base_strength = -1.0
    _dominance_last_final_strength = -1.0
    _dominance_last_concentration = -1.0
    _dominance_last_entropy_threshold = -1.0
    _dominance_last_auto_entropy = -1.0
    _dominance_last_auto_density = -1.0
    _dominance_last_auto_ref_factor = -1.0
    _dominance_last_auto_clamped = -1.0
    _dominance_logged_mode = False
    _dominance_logged_adaptive = False
    _dominance_logged_skip = False
    _dominance_last_mode_manual = False
    _dominance_last_iteration = -1
    _dominance_iteration_context = -1
    _dominance_reset_iteration_metrics()

    if config.enable_emergency_regularization:
        if config.emergency_entropy_threshold == 0.0:
            config.emergency_entropy_threshold = 0.15
        if config.emergency_max_weight_threshold == 0.0:
            config.emergency_max_weight_threshold = 0.3
        if config.emergency_min_dominant_refs == 0:
            config.emergency_min_dominant_refs = max_int32(10, pool.reference_count / 100)
        if config.emergency_likelihood_drop == 0.0:
            config.emergency_likelihood_drop = 1e-3

        config.dominance_regularization_active = False

        bf_nogil_logf_notime(
            b"EM",
            "em_session: mode=emergency_regularization entropy_threshold=%.3f max_weight=%.3f min_dominant_refs=%u ll_drop=%.6f",
            config.emergency_entropy_threshold,
            config.emergency_max_weight_threshold,
            config.emergency_min_dominant_refs,
            config.emergency_likelihood_drop,
        )
    elif config.enable_dominance_regularization:
        bf_nogil_logf_notime(
            b"EM",
            "em_session: mode=dominance_regularization base_strength=%.3f adaptive=%s entropy_scaling=%.3f",
            config.dominance_strength,
            b"true" if config.use_adaptive_dominance else b"false",
            config.entropy_scaling_factor,
        )
    else:
        bf_nogil_logf_notime(b"EM", "em_session: mode=standard")

    ll_history = <double*>calloc(5, sizeof(double))
    param_history = <double*>calloc(5, sizeof(double))
    if not ll_history or not param_history:
        if ll_history: free(ll_history)
        if param_history: free(param_history)
        return -1

    prev_weights = <double*>malloc(pool.reference_count * sizeof(double))
    if not prev_weights:
        free(ll_history)
        free(param_history)
        return -1

    precomp = create_precomputed_weights(pool.reference_count)
    if not precomp:
        free(prev_weights)
        free(ll_history)
        free(param_history)
        return -1

    if use_squarem:
        array_size = pool.reference_count * sizeof(double)
        total_size = array_size * 6
        squarem_memory = <double*>malloc(total_size)
        if not squarem_memory:
            free_precomputed_weights(precomp)
            free(prev_weights)
            free(ll_history)
            free(param_history)
            return -1
        theta_0 = squarem_memory
        theta_1 = squarem_memory + pool.reference_count
        theta_2 = squarem_memory + 2 * pool.reference_count
        r_vector = squarem_memory + 3 * pool.reference_count
        v_vector = squarem_memory + 4 * pool.reference_count
        theta_extrapolated = squarem_memory + 5 * pool.reference_count

    initialize_em_weights(pool, config)
    diagnose_initialization_quality(pool)

    memcpy(prev_weights, reference_weights, pool.reference_count * sizeof(double))

    if config.enable_emergency_regularization:
        bf_nogil_logf_notime(b"EM", "emergency_regularization: initial_state=inactive")
    elif config.enable_dominance_regularization:
        bf_nogil_logf_notime(
            b"EM",
            "dominance_regularization: initial_state=active base_strength=%.3f adaptive=%s entropy_scaling=%.3f range=[%.3f, %.3f]",
            config.dominance_strength,
            b"true" if config.use_adaptive_dominance else b"false",
            config.entropy_scaling_factor,
            config.min_penalty_strength,
            config.max_penalty_strength,
        )
    else:
        bf_nogil_logf_notime(b"EM", "regularization: mode=standard")

    bf_nogil_logf_notime(
        b"EM",
        "squarem: enabled=%s start_iteration=%d globalization=%s",
        b"true" if use_squarem else b"false",
        config.squarem_start_iter,
        b"true" if config.enable_globalization else b"false",
    )

    for iteration in range(max_iterations):
        prev_log_likelihood = current_log_likelihood
        memcpy(prev_weights, reference_weights, pool.reference_count * sizeof(double))

        if use_squarem and iteration >= config.squarem_start_iter:
            memcpy(theta_0, reference_weights, pool.reference_count * sizeof(double))

            execute_em_expectation_step_vectorized(pool, config.thread_count, precomp, config, iteration)
            execute_em_maximization_step_vectorized(pool, config.regularization_weight, precomp,
                                                   config, iteration, current_log_likelihood, prev_log_likelihood)
            memcpy(theta_1, reference_weights, pool.reference_count * sizeof(double))

            execute_em_expectation_step_vectorized(pool, config.thread_count, precomp, config, iteration)
            execute_em_maximization_step_vectorized(pool, config.regularization_weight, precomp,
                                                   config, iteration, current_log_likelihood, prev_log_likelihood)
            memcpy(theta_2, reference_weights, pool.reference_count * sizeof(double))

            r_norm = 0.0
            v_norm = 0.0
            ref_cnt = pool.reference_count
            i = 0
            while i + 8 <= ref_cnt:
                r_temp0 = theta_1[i]     - theta_0[i]
                r_temp1 = theta_1[i+1]   - theta_0[i+1]
                r_temp2 = theta_1[i+2]   - theta_0[i+2]
                r_temp3 = theta_1[i+3]   - theta_0[i+3]
                r_temp4 = theta_1[i+4]   - theta_0[i+4]
                r_temp5 = theta_1[i+5]   - theta_0[i+5]
                r_temp6 = theta_1[i+6]   - theta_0[i+6]
                r_temp7 = theta_1[i+7]   - theta_0[i+7]
                r_vector[i]   = r_temp0
                r_vector[i+1] = r_temp1
                r_vector[i+2] = r_temp2
                r_vector[i+3] = r_temp3
                r_vector[i+4] = r_temp4
                r_vector[i+5] = r_temp5
                r_vector[i+6] = r_temp6
                r_vector[i+7] = r_temp7
                v_temp0 = (theta_2[i]   - theta_1[i])   - r_temp0
                v_temp1 = (theta_2[i+1] - theta_1[i+1]) - r_temp1
                v_temp2 = (theta_2[i+2] - theta_1[i+2]) - r_temp2
                v_temp3 = (theta_2[i+3] - theta_1[i+3]) - r_temp3
                v_temp4 = (theta_2[i+4] - theta_1[i+4]) - r_temp4
                v_temp5 = (theta_2[i+5] - theta_1[i+5]) - r_temp5
                v_temp6 = (theta_2[i+6] - theta_1[i+6]) - r_temp6
                v_temp7 = (theta_2[i+7] - theta_1[i+7]) - r_temp7
                v_vector[i]   = v_temp0
                v_vector[i+1] = v_temp1
                v_vector[i+2] = v_temp2
                v_vector[i+3] = v_temp3
                v_vector[i+4] = v_temp4
                v_vector[i+5] = v_temp5
                v_vector[i+6] = v_temp6
                v_vector[i+7] = v_temp7
                r_norm += r_temp0*r_temp0 + r_temp1*r_temp1 + r_temp2*r_temp2 + r_temp3*r_temp3 + r_temp4*r_temp4 + r_temp5*r_temp5 + r_temp6*r_temp6 + r_temp7*r_temp7
                v_norm += v_temp0*v_temp0 + v_temp1*v_temp1 + v_temp2*v_temp2 + v_temp3*v_temp3 + v_temp4*v_temp4 + v_temp5*v_temp5 + v_temp6*v_temp6 + v_temp7*v_temp7
                i += 8

            while i < ref_cnt:
                r_temp = theta_1[i] - theta_0[i]
                r_vector[i] = r_temp
                v_temp = (theta_2[i] - theta_1[i]) - r_temp
                v_vector[i] = v_temp
                r_norm += r_temp * r_temp
                v_norm += v_temp * v_temp
                i += 1

            r_norm = libc_sqrt(r_norm)
            v_norm = libc_sqrt(v_norm)

            if v_norm > 1e-15:
                alpha = -r_norm / v_norm
            else:
                alpha = -1.0

            if alpha > -0.01:
                alpha = -0.01
            elif alpha < -50.0:
                alpha = -50.0

            for i in range(pool.reference_count):
                theta_extrapolated[i] = (theta_0[i] - 2.0 * alpha * r_vector[i] +
                                       alpha * alpha * v_vector[i])
                theta_extrapolated[i] = fmax(fmin(theta_extrapolated[i], 0.999), 1e-15)

            safe_normalize_weights(theta_extrapolated, pool.reference_count)

            backtrack_steps = 0
            if config.enable_globalization:
                memcpy(reference_weights, theta_extrapolated, pool.reference_count * sizeof(double))
                ll_before_backtrack = prev_log_likelihood
                test_ll = compute_log_likelihood_vectorized(pool, precomp)
                
                if test_ll < prev_log_likelihood - 1e-6:
                    bf_nogil_logf_notime(b"EM", "  BACKTRACKING: LL %.6f -> %.6f (alpha=%.3f triggered backtrack)",
                           ll_before_backtrack, test_ll, alpha)
                    
                    backtrack_factor = 0.5
                    while backtrack_steps < config.max_backtrack_steps and test_ll < prev_log_likelihood - 1e-6:
                        for i in range(pool.reference_count):
                            reference_weights[i] = (theta_0[i] + backtrack_factor * 
                                                   (theta_extrapolated[i] - theta_0[i]))
                            reference_weights[i] = fmax(fmin(reference_weights[i], 0.999), 1e-15)
                        
                        safe_normalize_weights(reference_weights, pool.reference_count)
                        test_ll = compute_log_likelihood_vectorized(pool, precomp)
                        
                        backtrack_factor *= config.backtrack_factor
                        backtrack_steps += 1
                    
                    ll_after_backtrack = test_ll
                    total_backtrack_steps += backtrack_steps
                    
                    bf_nogil_logf_notime(
                        b"EM",
                        "squarem_backtrack: steps=%d ll_before=%.6f ll_after=%.6f factor=%.3f",
                        backtrack_steps,
                        ll_before_backtrack,
                        ll_after_backtrack,
                        backtrack_factor,
                    )
                    
                    if test_ll < prev_log_likelihood - 1e-6:
                        bf_nogil_logf_notime(b"EM", "squarem_backtrack: status=fallback_to_standard_step")
                        memcpy(reference_weights, theta_1, pool.reference_count * sizeof(double))
                        alpha = -999.0
            else:
                memcpy(reference_weights, theta_extrapolated, pool.reference_count * sizeof(double))

            execute_em_expectation_step_vectorized(pool, config.thread_count, precomp, config, iteration)
            execute_em_maximization_step_vectorized(pool, config.regularization_weight, precomp,
                                                   config, iteration, current_log_likelihood, prev_log_likelihood)

        else:
            alpha = -999.0
            execute_em_expectation_step_vectorized(pool, config.thread_count, precomp, config, iteration)
            execute_em_maximization_step_vectorized(pool, config.regularization_weight, precomp,
                                                   config, iteration, current_log_likelihood, prev_log_likelihood)

        safe_normalize_weights(reference_weights, pool.reference_count)
        current_log_likelihood = compute_log_likelihood_vectorized(pool, precomp)

        if check_trend_based_convergence(current_log_likelihood, prev_log_likelihood,
                                        pool, prev_weights, tolerance, iteration,
                                        config, precomp, ll_history, param_history, alpha):
            pool.algorithm_converged = True
            bf_nogil_logf_notime(
                b"EM",
                "em_convergence: status=achieved iteration=%d",
                iteration + 1,
            )
            break

    if iteration >= max_iterations - 1:
        bf_nogil_logf_notime(
            b"EM",
            "em_convergence: status=max_iterations iterations=%d",
            max_iterations,
        )
        pool.algorithm_converged = False

    _dominance_flush_iteration_metrics(_dominance_last_iteration)

    if total_backtrack_steps > 0:
        bf_nogil_logf_notime(
            b"EM",
            "squarem_backtrack_summary: total_steps=%d iterations=%d avg_per_iteration=%.1f",
            total_backtrack_steps,
            iteration + 1,
            <double>total_backtrack_steps / <double>(iteration + 1),
        )

    pool.iteration_count = iteration + 1
    pool.final_log_likelihood = current_log_likelihood

    if squarem_memory: free(squarem_memory)
    free_precomputed_weights(precomp)
    free(prev_weights)
    free(ll_history)
    free(param_history)

    bf_nogil_logf_notime(
        b"EM",
        "em_summary: iterations=%d log_likelihood=%.6f converged=%s",
        pool.iteration_count,
        pool.final_log_likelihood,
        b"true" if pool.algorithm_converged else b"false",
    )
    
    if config.enable_emergency_regularization:
        bf_nogil_logf_notime(
            b"EM",
            "emergency_regularization: activated=%s",
            b"true" if config.dominance_regularization_active else b"false",
        )
    elif config.enable_dominance_regularization:
        bf_nogil_logf_notime(b"EM", "dominance_regularization: active_throughout=true")

    return 0

cdef int apply_probability_filtering_optimal(MemoryPool* pool, EMAlgorithmConfig* config) except -1 nogil:
    cdef uint32_t read_idx, ref_idx, rid
    cdef int64_t start_pos, end_pos, ai
    cdef uint32_t alignment_count
    cdef float alignment_score, uniform_zp
    cdef double log_lik 
    cdef double log_weighted, log_norm, posterior
    cdef double NEG_INF = -1e20
    cdef int64_t alignments_removed
    cdef int64_t write_idx = 0
    cdef bint keep_alignment
    cdef double min_threshold = fmax(config.minimum_probability_threshold, 1e-8)
    cdef double fraction_threshold = config.probability_fraction_filter
    cdef double* reference_weights = get_reference_weights(pool)
    cdef PrecomputedWeights* precomp = NULL
    cdef float* read_max_probs
    cdef int32_t* survivors_per_read
    cdef int64_t survivors_total = 0
    cdef uint64_t start_pos_rebuild
    cdef double thr
    cdef double p

    precomp = create_precomputed_weights(pool.reference_count)
    if not precomp:
        return -1
    update_precomputed_weights(precomp, reference_weights)

    bf_nogil_logf_notime(
        b"EM",
        "probability_filter_integrated: start pmd_output=%s total_alignments=%lld",
        b"enabled" if pool.pmd_enabled_for_output else b"disabled",
        <long long>pool.alignment_count,
    )

    if pool.scratch_read_max_probs == NULL or pool.scratch_survivors_per_read == NULL or pool.scratch_unique_read_count < <int32_t>pool.unique_read_count:
        if pool.scratch_read_max_probs != NULL:
            free(pool.scratch_read_max_probs)
        if pool.scratch_survivors_per_read != NULL:
            free(pool.scratch_survivors_per_read)
        pool.scratch_read_max_probs = <float*>calloc(pool.unique_read_count, sizeof(float))
        pool.scratch_survivors_per_read = <int32_t*>calloc(pool.unique_read_count, sizeof(int32_t))
        pool.scratch_unique_read_count = pool.unique_read_count
    else:
        memset(pool.scratch_read_max_probs, 0, pool.unique_read_count * sizeof(float))
        memset(pool.scratch_survivors_per_read, 0, pool.unique_read_count * sizeof(int32_t))
    read_max_probs = pool.scratch_read_max_probs
    survivors_per_read = pool.scratch_survivors_per_read

    for read_idx in prange(pool.unique_read_count, nogil=True, schedule='static', num_threads=config.thread_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue
        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count
        if alignment_count == 1:
            read_max_probs[read_idx] = 1.0
            if 1.0 >= min_threshold and (fraction_threshold == 0.0 or 1.0 >= fraction_threshold * 1.0):
                survivors_per_read[read_idx] = 1
            continue

        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                log_norm = stable_log_sum_exp(log_norm, log_weighted)
        if log_norm == NEG_INF:
            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            read_max_probs[read_idx] = uniform_zp
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                survivors_per_read[read_idx] = alignment_count
            continue

        p = 0.0
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                posterior = exp(log_weighted - log_norm)
                if posterior > p:
                    p = posterior
        read_max_probs[read_idx] = <float>p
        thr = fraction_threshold * p
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_lik = <double>alignment_score
                log_weighted = precomp.log_weights[ref_idx] + log_lik
                posterior = exp(log_weighted - log_norm)
                if posterior >= min_threshold and (fraction_threshold == 0.0 or posterior >= thr):
                    survivors_per_read[read_idx] += 1


    survivors_total = 0
    for read_idx in range(pool.unique_read_count):
        survivors_total += survivors_per_read[read_idx]

    bf_nogil_logf_notime(
        b"EM",
        "probability_filter_integrated: survivors=%lld total=%lld",
        <long long>survivors_total,
        <long long>pool.alignment_count,
    )

    if survivors_total <= 0:
        bf_nogil_logf_notime(b"WARN", "probability_filter_integrated: survivors=0 (aborting)")
        free_precomputed_weights(precomp)
        return -1


    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
    pool.precomputed_zp_values = <float*>malloc(survivors_total * sizeof(float))
    if not pool.precomputed_zp_values:
        free_precomputed_weights(precomp)
        return -1


    write_idx = 0
    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        if alignment_count == 1:
            p = 1.0
            if p >= min_threshold and (fraction_threshold == 0.0 or p >= fraction_threshold * 1.0):
                if write_idx != start_pos:
                    pool.alignments[write_idx] = pool.alignments[start_pos]
                pool.precomputed_zp_values[write_idx] = 1.0
                write_idx += 1
            continue


        log_norm = NEG_INF
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                log_norm = stable_log_sum_exp(log_norm, log_weighted)

        if log_norm == NEG_INF:

            uniform_zp = 1.0 / alignment_count if alignment_count > 0 else 0.0
            if uniform_zp >= min_threshold and (fraction_threshold == 0.0 or uniform_zp >= fraction_threshold * uniform_zp):
                for ai in range(start_pos, end_pos):
                    if write_idx != ai:
                        pool.alignments[write_idx] = pool.alignments[ai]
                    pool.precomputed_zp_values[write_idx] = uniform_zp
                    write_idx += 1
            continue

        thr = (<double>fraction_threshold) * (<double>read_max_probs[read_idx])
        for ai in range(start_pos, end_pos):
            ref_idx = pool.alignments[ai].reference_index
            if ref_idx < pool.reference_count:
                alignment_score = pool.alignments[ai].alignment_score
                log_weighted = precomp.log_weights[ref_idx] + alignment_score
                posterior = exp(log_weighted - log_norm)
                keep_alignment = (posterior >= min_threshold) and (fraction_threshold == 0.0 or posterior >= thr)
                if keep_alignment:
                    if write_idx != ai:
                        pool.alignments[write_idx] = pool.alignments[ai]
                    # clamp for stability
                    if posterior < 1e-12:
                        posterior = 1e-12
                    elif posterior > 0.999:
                        posterior = 0.999
                    pool.precomputed_zp_values[write_idx] = <float>posterior
                    write_idx += 1

    alignments_removed = pool.alignment_count - write_idx
    pool.alignment_count = write_idx
    pool.zp_values_computed = True

    bf_nogil_logf_notime(
        b"EM",
        "probability_filter_integrated: compacted_alignments=%llu removed=%lld",
        <unsigned long long>write_idx,
        <long long>alignments_removed,
    )


    memset(pool.read_alignment_counts, 0, pool.unique_read_count * sizeof(uint32_t))
    for ai in range(pool.alignment_count):
        rid = pool.alignments[ai].read_index
        if rid < pool.unique_read_count:
            pool.read_alignment_counts[rid] += 1
        else:
            bf_nogil_logf_notime(b"EM", "probability_filter_integrated_error: invalid_read_id=%u alignment=%llu", rid, <unsigned long long>ai)

    start_pos_rebuild = 0
    for rid in range(pool.unique_read_count):
        pool.read_alignment_starts[rid] = start_pos_rebuild
        start_pos_rebuild += pool.read_alignment_counts[rid]

    pool.final_unique_reads = 0
    for rid in range(pool.unique_read_count):
        if pool.read_alignment_counts[rid] > 0:
            pool.final_unique_reads += 1

    free_precomputed_weights(precomp)

    bf_nogil_logf_notime(b"EM", "probability_filter_integrated: completed")
    return 0


cdef int check_trend_based_convergence(double current_ll, double prev_ll,
                                      MemoryPool* pool, double* prev_weights,
                                      double base_tolerance, int iteration,
                                      EMAlgorithmConfig* config,
                                      PrecomputedWeights* precomp,
                                      double* ll_history, double* param_history,
                                      double alpha) noexcept nogil:
 
    cdef bint is_squarem
    cdef double ll_change, ll_rel_change, residual_norm
    cdef double dimension_scale, scaled_residual
    cdef int H, idx, filled
    cdef bint ll_ok_strict, param_ok_strict
    cdef bint have_models
    cdef double ll_tail, r_inf
    cdef int i0, i1, j0, j1, j2
    cdef double d_k, d_km1, rho
    cdef double r_k, r_km1, r_km2, dr, dr_prev, denom
    cdef bint pred_ll_ok, pred_r_ok
    cdef double ll_med, ll_sigma, pr_med, pr_sigma
    cdef double ll_eps_floor, param_eps_floor
    cdef double ll_tol_eff, param_tol_eff

    is_squarem = (alpha > -900.0)

    ll_change = fabs(current_ll - prev_ll)
    ll_rel_change = ll_change / fmax(fabs(current_ll), 1.0)
    residual_norm = compute_residual_norm_efficient(pool, prev_weights, config, precomp)

    dimension_scale = libc_sqrt(<double>pool.reference_count)
    scaled_residual = residual_norm / dimension_scale

    H = 5
    idx = iteration % H
    ll_history[idx] = ll_rel_change
    param_history[idx] = residual_norm
    filled = _count_filled(iteration, H)

    _robust_sigma(ll_history, H, filled, &ll_med, &ll_sigma)
    _robust_sigma(param_history, H, filled, &pr_med, &pr_sigma)

    ll_eps_floor    = 10.0 * DBL_EPSILON * fmax(1.0, fabs(current_ll))
    param_eps_floor = 10.0 * DBL_EPSILON

    ll_tol_eff    = fmax(base_tolerance, fmax(3.0 * ll_sigma, ll_eps_floor))
    param_tol_eff = fmax(base_tolerance, fmax(3.0 * (pr_sigma / dimension_scale), param_eps_floor))

    ll_ok_strict    = (ll_rel_change <= base_tolerance)
    param_ok_strict = (residual_norm <= base_tolerance) or (scaled_residual <= base_tolerance)

    if ll_ok_strict and param_ok_strict:
        bf_nogil_logf_notime(b"EM", "Iter %d: CONVERGED [strict] DLL=%.2e<=%.1e, ||F||=%.2e, ||F||sqrt(p)=%.2e<=%.1e",
               iteration + 1, ll_rel_change, base_tolerance, residual_norm, scaled_residual, base_tolerance)
        return True

    have_models = (filled >= 3)
    ll_tail = 1e300
    r_inf   = 1e300

    if have_models:
        i0 = idx
        i1 = (idx - 1 + H) % H
        d_k   = ll_history[i0]
        d_km1 = ll_history[i1]
        if d_km1 > 0.0:
            rho = d_k / d_km1
            if rho >= 0.0 and rho < 1.0 and d_k < d_km1:
                ll_tail = d_k * rho / (1.0 - rho)

        j0 = idx
        j1 = (idx - 1 + H) % H
        j2 = (idx - 2 + H) % H
        r_k   = param_history[j0]
        r_km1 = param_history[j1]
        r_km2 = param_history[j2]
        dr      = r_k   - r_km1
        dr_prev = r_km1 - r_km2
        denom   = (dr - dr_prev)
        if (dr < 0.0) and (denom != 0.0) and (denom < 0.0):
            r_inf = r_k - (dr * dr) / denom
        else:
            r_inf = r_k

    pred_ll_ok = have_models and ((ll_tail <= ll_tol_eff) or (ll_rel_change <= ll_tol_eff))
    pred_r_ok = have_models and ((r_inf / dimension_scale) <= param_tol_eff)

    if pred_ll_ok and pred_r_ok:
        bf_nogil_logf_notime(b"EM", "Iter %d: CONVERGED [predictive+robust] tail_LL=%.2e<=%.1e, Aitken(||F||)->%.2e (scaled=%.2e)<=%.1e",
               iteration + 1, ll_tail, ll_tol_eff, r_inf, r_inf / dimension_scale, param_tol_eff)
        return True

    if is_squarem:
        bf_nogil_logf_notime(
            b"EM",
            "Iter %d: DLL=%.2e (base=%.1e, eff=%.1e) ||F||=%.2e, ||F||sqrt(p)=%.2e (base=%.1e, eff=%.1e) alpha=%.3f "
            "[strict ok: LL=%s Param=%s | pred ok: LL=%s Param=%s | models=%s]\n",
            iteration + 1,
            ll_rel_change,
            base_tolerance,
            ll_tol_eff,
            residual_norm,
            scaled_residual,
            base_tolerance,
            param_tol_eff,
            alpha,
            b"Y" if ll_ok_strict else b"N",
            b"Y" if param_ok_strict else b"N",
            b"Y" if pred_ll_ok else b"N",
            b"Y" if pred_r_ok else b"N",
            b"Y" if have_models else b"N",
        )
    else:
        bf_nogil_logf_notime(
            b"EM",
            "Iter %d: DLL=%.2e (base=%.1e, eff=%.1e) ||F||=%.2e, ||F||sqrt(p)=%.2e (base=%.1e, eff=%.1e) "
            "[strict ok: LL=%s Param=%s | pred ok: LL=%s Param=%s | models=%s]\n",
            iteration + 1,
            ll_rel_change,
            base_tolerance,
            ll_tol_eff,
            residual_norm,
            scaled_residual,
            base_tolerance,
            param_tol_eff,
            b"Y" if ll_ok_strict else b"N",
            b"Y" if param_ok_strict else b"N",
            b"Y" if pred_ll_ok else b"N",
            b"Y" if pred_r_ok else b"N",
            b"Y" if have_models else b"N",
        )

    return False

cdef double compute_residual_norm_efficient(MemoryPool* pool, double* prev_weights,
                                           EMAlgorithmConfig* config,
                                           PrecomputedWeights* precomp) except -1.0 nogil:

    cdef double* current_weights = get_reference_weights(pool)
    cdef double* temp_weights = get_temp_buffer_B(pool)
    cdef double norm_squared = 0.0
    cdef uint32_t i
    cdef double diff

    memcpy(temp_weights, current_weights, pool.reference_count * sizeof(double))

    memcpy(current_weights, prev_weights, pool.reference_count * sizeof(double))

    cdef int32_t thread_cnt = config.thread_count
    execute_em_expectation_step_vectorized(pool, thread_cnt, precomp)
    execute_em_maximization_step_vectorized(pool, config.regularization_weight, precomp,
                                           config, 0, -1e20, -1e20)

    for i in range(pool.reference_count):
        diff = current_weights[i] - prev_weights[i]
        norm_squared += diff * diff

    memcpy(current_weights, temp_weights, pool.reference_count * sizeof(double))

    return libc_sqrt(norm_squared)


cdef void compute_initialization_stats(MemoryPool* pool, InitializationStats* stats) noexcept nogil:
    cdef int64_t i
    cdef float score
    cdef double sum_score = 0.0
    
    stats.min_score = 1e30
    stats.max_score = -1e30
    stats.total_alignments = pool.alignment_count
    
    for i in range(pool.alignment_count):
        score = pool.alignments[i].alignment_score
        if score < stats.min_score:
            stats.min_score = score
        if score > stats.max_score:
            stats.max_score = score
        sum_score += score
    
    stats.mean_score = sum_score / pool.alignment_count if pool.alignment_count > 0 else 0.0
    stats.score_range = stats.max_score - stats.min_score
    
    bf_nogil_logf_notime(
        b"EM",
        "initialization_stats: min_score=%.3f max_score=%.3f mean_score=%.3f range=%.3f",
        stats.min_score,
        stats.max_score,
        stats.mean_score,
        stats.score_range,
    )


cdef int64_t identify_unique_alignments(MemoryPool* pool, InitializationStats* stats,
                                       UniqueAlignment* unique_alignments, int64_t capacity) noexcept nogil:
    cdef int64_t unique_count = 0
    cdef UniqueAlignment* uniques = unique_alignments
    cdef uint32_t read_idx
    cdef uint64_t start_pos, end_pos, ai
    cdef uint32_t alignment_count
    cdef float best_score, second_best_score, score_gap, current_score
    cdef uint32_t best_ref_idx, second_best_ref_idx
    cdef double dynamic_threshold
    
    if not uniques:
        return -1
    
    dynamic_threshold = fmax(5.0, stats.score_range * 0.2)
    
    bf_nogil_logf_notime(
        b"EM",
        "unique_detection: gap_threshold=%.3f range_fraction=%.1f%%",
        dynamic_threshold,
        (dynamic_threshold / stats.score_range) * 100.0 if stats.score_range > 0 else 0.0,
    )
    
    for read_idx in range(pool.unique_read_count):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue
            
        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count
        
        if alignment_count == 1:
            ai = start_pos
            uniques[unique_count].read_index = read_idx
            uniques[unique_count].reference_index = pool.alignments[ai].reference_index
            uniques[unique_count].alignment_score = pool.alignments[ai].alignment_score
            uniques[unique_count].confidence_score = 100.0 
            unique_count += 1
            
        else:
            best_score = -1e30
            second_best_score = -1e30
            best_ref_idx = 0
            second_best_ref_idx = 0
            
            for ai in range(start_pos, end_pos):
                current_score = pool.alignments[ai].alignment_score
                if (current_score > best_score) or (current_score == best_score and pool.alignments[ai].reference_index < best_ref_idx):
                    second_best_score = best_score
                    second_best_ref_idx = best_ref_idx
                    best_score = current_score
                    best_ref_idx = pool.alignments[ai].reference_index
                elif (current_score > second_best_score) or (current_score == second_best_score and pool.alignments[ai].reference_index < second_best_ref_idx):
                    second_best_score = current_score
                    second_best_ref_idx = pool.alignments[ai].reference_index
            
            score_gap = best_score - second_best_score
            if score_gap >= dynamic_threshold:
                uniques[unique_count].read_index = read_idx
                uniques[unique_count].reference_index = best_ref_idx
                uniques[unique_count].alignment_score = best_score
                uniques[unique_count].confidence_score = score_gap
                unique_count += 1
    
    stats.unique_alignments = unique_count
    stats.uniqueness_ratio = <double>unique_count / <double>pool.final_unique_reads
    
    bf_nogil_logf_notime(
        b"EM",
        "unique_detection: unique_alignments=%lld total_reads=%u coverage=%.1f%%",
        <long long>unique_count,
        pool.final_unique_reads,
        stats.uniqueness_ratio * 100.0,
    )
    
    return unique_count

cdef void compute_quality_weights(MemoryPool* pool, InitializationStats* stats,
                                UniqueAlignment* unique_alignments, int64_t unique_count,
                                double* quality_weighted_counts) noexcept nogil:
    cdef int64_t i
    cdef uint32_t ref_idx
    cdef double score_percentile, quality_weight
    cdef double score_normalized, score_weight, confidence_weight
    
    memset(quality_weighted_counts, 0, pool.reference_count * sizeof(double))

    if unique_count <= 0:
        bf_nogil_logf_notime(b"EM", "quality_weighting: applied=false unique_alignments=0")
        return

    for i in range(unique_count):
        ref_idx = unique_alignments[i].reference_index
        if ref_idx >= pool.reference_count:
            continue
        
        if stats.score_range > 0.0:
            score_normalized = (unique_alignments[i].alignment_score - stats.min_score) / stats.score_range
        else:
            score_normalized = 0.5

        score_weight = 0.1 + 0.9 * score_normalized
        confidence_weight = fmin(1.0, unique_alignments[i].confidence_score / 10.0)
        quality_weight = 0.7 * score_weight + 0.3 * confidence_weight
        quality_weighted_counts[ref_idx] += quality_weight

    bf_nogil_logf_notime(
        b"EM",
        "quality_weighting: applied=true unique_alignments=%lld",
        <long long>unique_count,
    )

cdef void initialize_em_weights(MemoryPool* pool, EMAlgorithmConfig* config) noexcept nogil:
    """Initialize reference weights using quality-weighted unique alignments.

    Uses simplified initialization approach:
    - Quality-weighted counting from unique alignments
    - Dirichlet prior for regularization
    - No length normalization (inappropriate for short ancient DNA)
    - No penalty pre-correction (moved to post-EM filtering)

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool to initialize
    config : EMAlgorithmConfig*
        Configuration with prior strength parameter
    """
    cdef InitializationStats stats
    cdef UniqueAlignment* unique_alignments = NULL
    cdef int64_t unique_count
    cdef double* reference_weights = get_reference_weights(pool)
    cdef double* quality_counts = get_temp_buffer_A(pool)
    
    cdef uint32_t ref_idx
    cdef int64_t i
    cdef double uniform_weight = 1.0 / pool.reference_count
    cdef double alpha_prior = config.init_prior_strength  # Dirichlet concentration parameter from config
    cdef double total_sum = 0.0

    bf_nogil_logf_notime(b"EM", "initialization: stage=start")

    compute_initialization_stats(pool, &stats)
    
    unique_alignments = <UniqueAlignment*>malloc(pool.final_unique_reads * sizeof(UniqueAlignment))
    if not unique_alignments:
        bf_nogil_logf_notime(b"EM", "initialization_fallback: reason=allocation_failed strategy=uniform")
        for ref_idx in range(pool.reference_count):
            reference_weights[ref_idx] = uniform_weight
        return
    
    unique_count = identify_unique_alignments(pool, &stats, unique_alignments, pool.final_unique_reads)
    if unique_count <= 0:
        bf_nogil_logf_notime(b"EM", "initialization_fallback: reason=no_unique_alignments strategy=uniform")
        for ref_idx in range(pool.reference_count):
            reference_weights[ref_idx] = uniform_weight
        free(unique_alignments)
        return
    
    # Compute quality-weighted counts (keeps existing logic)
    compute_quality_weights(pool, &stats, unique_alignments, unique_count, quality_counts)

    total_sum = 0.0
    for ref_idx in range(pool.reference_count):
        reference_weights[ref_idx] = quality_counts[ref_idx] + alpha_prior
        total_sum += reference_weights[ref_idx]
    
    if total_sum > 0.0:
        for ref_idx in range(pool.reference_count):
            reference_weights[ref_idx] /= total_sum
    else:
        for ref_idx in range(pool.reference_count):
            reference_weights[ref_idx] = uniform_weight

    free(unique_alignments)
    cdef double min_weight = 1e30
    cdef double max_weight = 0.0
    cdef uint32_t refs_with_evidence = 0
    
    for ref_idx in range(pool.reference_count):
        if reference_weights[ref_idx] < min_weight:
            min_weight = reference_weights[ref_idx]
        if reference_weights[ref_idx] > max_weight:
            max_weight = reference_weights[ref_idx]
        if quality_counts[ref_idx] > 0.0:
            refs_with_evidence += 1
    
    bf_nogil_logf_notime(
        b"EM",
        "initialization_summary: unique_alignments=%lld references_with_signal=%u/%u weight_min=%.6f weight_max=%.6f prior=%.2f",
        <long long>unique_count,
        refs_with_evidence,
        pool.reference_count,
        min_weight,
        max_weight,
        alpha_prior,
    )

cdef void diagnose_initialization_quality(MemoryPool* pool) noexcept nogil:
    """
    Print diagnostics about the initial EM weight vector.

    This helper computes basic summary statistics (min/max, entropy,
    dominant reference counts) for the reference weight distribution and
    prints human-readable diagnostics to aid debugging of EM initialization.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing the initialized reference weights.
    """
    cdef double* weights = get_reference_weights(pool)
    cdef double entropy = 0.0
    cdef double max_weight = 0.0
    cdef double min_weight = 1e30
    cdef uint32_t dominant_refs = 0
    cdef uint32_t ref_idx
    cdef double w
    cdef double max_entropy = log(<double>pool.reference_count)
    cdef double entropy_ratio
    
    for ref_idx in range(pool.reference_count):
        w = weights[ref_idx]
        if w > max_weight:
            max_weight = w
        if w < min_weight:
            min_weight = w
        if w > 2.0 / pool.reference_count:
            dominant_refs += 1
        if w > 1e-15:
            entropy -= w * log(w)
    
    entropy_ratio = entropy / max_entropy if max_entropy > 0.0 else 0.0
    
    bf_nogil_logf_notime(
        b"EM",
        "initialization_diagnostics: weight_min=%.6f weight_max=%.6f ratio=%.1fx entropy=%.3f/%.3f entropy_pct=%.1f dominant_refs=%u/%u (%.1f%%)",
        min_weight,
        max_weight,
        max_weight / fmax(min_weight, 1e-15),
        entropy,
        max_entropy,
        entropy_ratio * 100.0,
        dominant_refs,
        pool.reference_count,
        (dominant_refs * 100.0) / pool.reference_count,
    )
