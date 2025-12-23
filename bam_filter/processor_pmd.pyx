# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# -*- coding: utf-8 -*-
"""PMD (Post-Mortem Damage) learning and damage-corrected ANI computation.

Implements:
1. Thread-safe PMD statistics accumulation during batch processing
2. Damage curve D(z) estimation with Beta prior regularization
3. Hierarchical ancient/modern reference classification
4. On-the-fly damage-corrected ANI computation
"""

from libc.math cimport exp, log, log1p, fabs, fmax, fmin, sqrt, pow as cpow
from libc.stdlib cimport malloc, calloc, free, realloc
from libc.string cimport memset, memcpy
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t, uintptr_t
from cython.parallel cimport prange

from bam_filter.processor_pmd cimport *
from bam_filter.processor cimport MemoryPool, Alignment

# Logging
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef const char* LOG_TAG = "pmd"


# =============================================================================
# Initialization and Cleanup
# =============================================================================

cdef PMDGlobalContext* create_pmd_context(int32_t num_threads,
                                           bint is_single_stranded,
                                           bint enable_hierarchical) noexcept nogil:
    """Create and initialize global PMD context.

    Parameters
    ----------
    num_threads : int32_t
        Number of threads for parallel processing
    is_single_stranded : bint
        True for single-stranded library, False for double-stranded
    enable_hierarchical : bint
        True to enable hierarchical ancient/modern classification

    Returns
    -------
    PMDGlobalContext*
        Initialized context, or NULL on allocation failure
    """
    cdef PMDGlobalContext* ctx = <PMDGlobalContext*>calloc(1, sizeof(PMDGlobalContext))
    if not ctx:
        return NULL

    cdef int32_t max_threads = num_threads if num_threads <= PMD_MAX_THREADS else PMD_MAX_THREADS

    ctx.thread_contexts = <PMDThreadContext*>calloc(max_threads, sizeof(PMDThreadContext))
    if not ctx.thread_contexts:
        free(ctx)
        return NULL

    ctx.num_threads = num_threads
    ctx.max_threads = max_threads
    ctx.is_single_stranded = is_single_stranded
    ctx.collect_stats = True
    ctx.enable_hierarchical = enable_hierarchical

    # Initialize thread contexts
    cdef int32_t i
    for i in range(max_threads):
        ctx.thread_contexts[i].thread_id = i
        ctx.thread_contexts[i].initialized = True
        # Stats are zero-initialized by calloc

    # Initialize model with default parameters
    ctx.model.params.P_prior_mean = 0.3
    ctx.model.params.P_prior_sd = 0.15
    ctx.model.params.lambda_prior_mean = 0.35  # Corresponds to DECAY=0.7
    ctx.model.params.lambda_prior_sd = 0.2
    ctx.model.params.C_prior_mean = 0.01
    ctx.model.params.C_prior_sd = 0.005
    ctx.model.params.omega_alpha = 0.5
    ctx.model.params.omega_beta = 5.0
    ctx.model.params.epsilon = 0.01

    ctx.model.curve.epsilon_baseline = 0.01
    ctx.model.curve.is_single_stranded = is_single_stranded

    ctx.initialized = True
    ctx.owns_memory = True

    bf_nogil_logf_notime(LOG_TAG, "PMD context created: threads=%d ss=%s hierarchical=%s",
                         num_threads,
                         b"yes" if is_single_stranded else b"no",
                         b"yes" if enable_hierarchical else b"no")

    return ctx


cdef void destroy_pmd_context(PMDGlobalContext* ctx) noexcept nogil:
    """Free all memory associated with PMD context."""
    if not ctx:
        return

    if ctx.thread_contexts:
        free(ctx.thread_contexts)

    if ctx.model.hierarchical:
        if ctx.model.hierarchical.ref_classes:
            free(ctx.model.hierarchical.ref_classes)
        free(ctx.model.hierarchical)

    if ctx.owns_memory:
        free(ctx)


cdef void reset_pmd_stats(PMDGlobalContext* ctx) noexcept nogil:
    """Reset all PMD statistics (thread-local and global)."""
    if not ctx:
        return

    cdef int32_t i
    for i in range(ctx.max_threads):
        memset(&ctx.thread_contexts[i].stats, 0, sizeof(PMDStatsAccumulator))

    memset(&ctx.model.stats, 0, sizeof(PMDStatsGlobal))
    ctx.model.stats_collected = False
    ctx.model.curve_fitted = False
    ctx.model.finalized = False


# =============================================================================
# Thread-Local Statistics Collection
# =============================================================================

cdef inline void pmd_record_position(PMDStatsAccumulator* acc,
                                      int end_type,
                                      int context_type,
                                      int position,
                                      bint is_mismatch) noexcept nogil:
    """Record a single position for PMD statistics.

    Called during MD tag parsing for each C (at 5' end) or G (at 3' end) position.

    Parameters
    ----------
    acc : PMDStatsAccumulator*
        Thread-local accumulator
    end_type : int
        PMD_END_5P (0) for 5' C→T, PMD_END_3P (1) for 3' G→A
    context_type : int
        PMD_CTX_NONCPG (0) or PMD_CTX_CPG (1)
    position : int
        1-based position from the relevant end (1 = first base)
    is_mismatch : bint
        True if this is a damage-consistent mismatch (C→T or G→A)
    """
    if position < 1 or position > PMD_MAX_POSITION:
        return

    cdef int idx = position - 1  # Convert to 0-based index

    if end_type == PMD_END_5P:
        if context_type == PMD_CTX_CPG:
            acc.n_5p_cpg[idx] += 1
            if is_mismatch:
                acc.k_5p_cpg[idx] += 1
        else:
            acc.n_5p_noncpg[idx] += 1
            if is_mismatch:
                acc.k_5p_noncpg[idx] += 1
    else:  # PMD_END_3P
        if context_type == PMD_CTX_CPG:
            acc.n_3p_cpg[idx] += 1
            if is_mismatch:
                acc.k_3p_cpg[idx] += 1
        else:
            acc.n_3p_noncpg[idx] += 1
            if is_mismatch:
                acc.k_3p_noncpg[idx] += 1


# =============================================================================
# Statistics Merging
# =============================================================================

cdef void merge_pmd_stats(PMDGlobalContext* ctx) noexcept nogil:
    """Merge thread-local PMD statistics into global counts.

    Called after all batch processing is complete, before curve fitting.
    """
    if not ctx or ctx.model.stats_collected:
        return

    cdef PMDStatsGlobal* global_stats = &ctx.model.stats
    cdef PMDStatsAccumulator* thread_stats
    cdef int32_t t, z

    # Reset global stats
    memset(global_stats, 0, sizeof(PMDStatsGlobal))

    # Merge from all threads
    for t in range(ctx.max_threads):
        thread_stats = &ctx.thread_contexts[t].stats

        for z in range(PMD_MAX_POSITION):
            global_stats.n_5p_noncpg[z] += thread_stats.n_5p_noncpg[z]
            global_stats.k_5p_noncpg[z] += thread_stats.k_5p_noncpg[z]
            global_stats.n_5p_cpg[z] += thread_stats.n_5p_cpg[z]
            global_stats.k_5p_cpg[z] += thread_stats.k_5p_cpg[z]
            global_stats.n_3p_noncpg[z] += thread_stats.n_3p_noncpg[z]
            global_stats.k_3p_noncpg[z] += thread_stats.k_3p_noncpg[z]
            global_stats.n_3p_cpg[z] += thread_stats.n_3p_cpg[z]
            global_stats.k_3p_cpg[z] += thread_stats.k_3p_cpg[z]

        global_stats.total_alignments += thread_stats.total_alignments
        global_stats.total_bases += thread_stats.total_bases

    ctx.model.stats_collected = True

    bf_nogil_logf_notime(LOG_TAG, "PMD stats merged: alignments=%llu bases=%llu",
                         global_stats.total_alignments, global_stats.total_bases)

    # Log summary of damage rates at position 1
    cdef double rate_5p = 0.0, rate_3p = 0.0
    if global_stats.n_5p_noncpg[0] + global_stats.n_5p_cpg[0] > 0:
        rate_5p = <double>(global_stats.k_5p_noncpg[0] + global_stats.k_5p_cpg[0]) / \
                  <double>(global_stats.n_5p_noncpg[0] + global_stats.n_5p_cpg[0])
    if global_stats.n_3p_noncpg[0] + global_stats.n_3p_cpg[0] > 0:
        rate_3p = <double>(global_stats.k_3p_noncpg[0] + global_stats.k_3p_cpg[0]) / \
                  <double>(global_stats.n_3p_noncpg[0] + global_stats.n_3p_cpg[0])

    bf_nogil_logf_notime(LOG_TAG, "Observed damage at pos 1: 5' C->T=%.3f%% 3' G->A=%.3f%%",
                         rate_5p * 100.0, rate_3p * 100.0)


# =============================================================================
# Damage Curve Fitting
# =============================================================================

cdef inline double _fit_omega_objective(double omega,
                                         uint64_t* n_arr, uint64_t* k_arr,
                                         double* f_arr,  # D_0(z) shape values
                                         double epsilon,
                                         double alpha, double beta,
                                         int n_positions) noexcept nogil:
    """Compute negative log-posterior for omega given fixed shape.

    Objective: -log p(omega | data) = -log L(data | omega) - log p(omega)

    where L uses binomial likelihood and p(omega) is Beta(alpha, beta).
    """
    cdef double log_lik = 0.0
    cdef double q_z, log_q, log_1mq
    cdef int z
    cdef uint64_t n, k

    for z in range(n_positions):
        n = n_arr[z]
        k = k_arr[z]
        if n == 0:
            continue

        # q(z; omega) = epsilon + omega * f(z)
        q_z = epsilon + omega * f_arr[z]
        q_z = fmax(1e-10, fmin(q_z, 1.0 - 1e-10))

        log_q = log(q_z)
        log_1mq = log(1.0 - q_z)

        log_lik += k * log_q + (n - k) * log_1mq

    # Beta prior: (alpha-1)*log(omega) + (beta-1)*log(1-omega)
    cdef double log_prior = 0.0
    if omega > 1e-10:
        log_prior += (alpha - 1.0) * log(omega + 1e-10)
    if omega < 1.0 - 1e-10:
        log_prior += (beta - 1.0) * log(1.0 - omega + 1e-10)

    return -(log_lik + log_prior)


cdef double _fit_omega(uint64_t* n_5p, uint64_t* k_5p,
                        uint64_t* n_3p, uint64_t* k_3p,
                        double* f_arr,
                        double epsilon,
                        double alpha, double beta) noexcept nogil:
    """Fit omega using golden section search.

    Combines 5' and 3' data for joint estimation.
    """
    cdef double a = 0.0, b = 1.0
    cdef double gr = 0.6180339887  # Golden ratio
    cdef double c = b - gr * (b - a)
    cdef double d = a + gr * (b - a)
    cdef double fc, fd
    cdef double tol = 1e-6
    cdef int max_iter = 50, i

    # Combine n and k arrays for joint fitting
    cdef uint64_t n_combined[40]
    cdef uint64_t k_combined[40]
    cdef double f_combined[40]

    for i in range(20):
        n_combined[i] = n_5p[i]
        k_combined[i] = k_5p[i]
        f_combined[i] = f_arr[i]
        n_combined[20 + i] = n_3p[i]
        k_combined[20 + i] = k_3p[i]
        f_combined[20 + i] = f_arr[i]

    for i in range(max_iter):
        fc = _fit_omega_objective(c, n_combined, k_combined, f_combined,
                                  epsilon, alpha, beta, 40)
        fd = _fit_omega_objective(d, n_combined, k_combined, f_combined,
                                  epsilon, alpha, beta, 40)

        if fc < fd:
            b = d
            d = c
            c = b - gr * (b - a)
        else:
            a = c
            c = d
            d = a + gr * (b - a)

        if fabs(b - a) < tol:
            break

    return (a + b) / 2.0


cdef double _estimate_lambda_from_data(
    uint64_t* n_5p, uint64_t* k_5p,
    uint64_t* n_3p, uint64_t* k_3p,
    double prior_mean, double prior_sd
) noexcept nogil:
    """Estimate lambda decay parameter from observed damage rates.

    Uses linear regression on log-transformed damage rates for positions 1-10.
    Baseline is estimated from positions 15-19.

    Returns lambda with shrinkage toward prior if data is insufficient.
    """
    cdef double baseline = 0.0
    cdef double baseline_n = 0.0
    cdef int z
    cdef double rate

    # Estimate baseline from positions 15-19 (0-indexed: 14-18)
    for z in range(14, 19):
        if n_5p[z] > 0:
            baseline += <double>k_5p[z] / <double>n_5p[z]
            baseline_n += 1.0
        if n_3p[z] > 0:
            baseline += <double>k_3p[z] / <double>n_3p[z]
            baseline_n += 1.0

    if baseline_n > 0:
        baseline = baseline / baseline_n
    else:
        baseline = 0.01

    # Linear regression: log(D[z] - baseline) = log(P) - lambda * (z-1)
    # We fit: y = a + b*x where y = log(D - baseline), x = z-1, b = -lambda
    cdef double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0, sum_xx = 0.0
    cdef double n_points = 0.0
    cdef double x, y, D_obs

    # Use positions 1-10 for fitting (0-indexed: 0-9)
    for z in range(10):
        # Combine 5' and 3' data
        if n_5p[z] > 100:  # Require at least 100 observations
            D_obs = <double>k_5p[z] / <double>n_5p[z]
            if D_obs > baseline + 0.005:  # Must be above baseline
                x = <double>z  # z-1 where z is 1-based, so z (0-indexed) = z-1
                y = log(D_obs - baseline)
                sum_x += x
                sum_y += y
                sum_xy += x * y
                sum_xx += x * x
                n_points += 1.0

        if n_3p[z] > 100:
            D_obs = <double>k_3p[z] / <double>n_3p[z]
            if D_obs > baseline + 0.005:
                x = <double>z
                y = log(D_obs - baseline)
                sum_x += x
                sum_y += y
                sum_xy += x * y
                sum_xx += x * x
                n_points += 1.0

    # Need at least 3 points for meaningful fit
    if n_points < 3:
        bf_nogil_logf_notime(LOG_TAG, "Insufficient data for lambda fit, using prior (n=%.0f)", n_points)
        return prior_mean

    # Linear regression: b = (n*sum_xy - sum_x*sum_y) / (n*sum_xx - sum_x^2)
    cdef double denom = n_points * sum_xx - sum_x * sum_x
    if fabs(denom) < 1e-10:
        return prior_mean

    cdef double b = (n_points * sum_xy - sum_x * sum_y) / denom
    cdef double lambda_fit = -b  # lambda = -slope

    # Clamp to reasonable range [0.1, 1.0]
    if lambda_fit < 0.1:
        lambda_fit = 0.1
    elif lambda_fit > 1.0:
        lambda_fit = 1.0

    # Shrink toward prior based on data quantity
    # weight = n_points / (n_points + prior_strength)
    cdef double prior_strength = 10.0
    cdef double weight = n_points / (n_points + prior_strength)
    cdef double lambda_final = weight * lambda_fit + (1.0 - weight) * prior_mean

    bf_nogil_logf_notime(LOG_TAG, "Lambda fit: raw=%.4f prior=%.4f final=%.4f (n=%.0f, baseline=%.4f)",
                         lambda_fit, prior_mean, lambda_final, n_points, baseline)

    return lambda_final


cdef int fit_pmd_curve(PMDGlobalContext* ctx, PMDCurveParams* params) noexcept nogil:
    """Fit PMD damage curve from collected statistics.

    Uses regularized MLE with Beta prior on omega for shrinkage toward
    no-damage when evidence is weak.

    Parameters
    ----------
    ctx : PMDGlobalContext*
        Global context with merged statistics
    params : PMDCurveParams*
        Fitting parameters (priors, epsilon)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    if not ctx or not ctx.model.stats_collected:
        return -1

    cdef PMDCurveParams* p = params if params else &ctx.model.params
    cdef PMDStatsGlobal* stats = &ctx.model.stats
    cdef PMDCurve* curve = &ctx.model.curve

    # Fit lambda from observed damage rates (with shrinkage to prior)
    cdef double lam = _estimate_lambda_from_data(
        stats.n_5p_noncpg, stats.k_5p_noncpg,
        stats.n_3p_noncpg, stats.k_3p_noncpg,
        p.lambda_prior_mean, p.lambda_prior_sd
    )

    # Use prior means for P and C (could optimize these too in future)
    cdef double P = p.P_prior_mean
    cdef double C = p.C_prior_mean
    cdef double epsilon = p.epsilon

    # Precompute D_0(z) shape: f(z) = P * exp(-lambda*(z-1)) + C
    cdef double f_arr[20]
    cdef int z
    for z in range(20):
        f_arr[z] = P * exp(-lam * z) + C

    # Fit omega for each end/context combination
    cdef double omega_5p_noncpg = _fit_omega(
        stats.n_5p_noncpg, stats.k_5p_noncpg,
        stats.n_3p_noncpg, stats.k_3p_noncpg,  # Use 3' for joint fit
        f_arr, epsilon, p.omega_alpha, p.omega_beta)

    # For simplicity, use single omega (can extend to per-context omega)
    cdef double omega = omega_5p_noncpg

    # Store fitted omega
    curve.omega = <float>omega

    # Compute final D(z) curves
    for z in range(20):
        curve.D_5p_noncpg[z] = <float>(omega * f_arr[z])
        curve.D_5p_cpg[z] = <float>(omega * f_arr[z] * 1.2)  # CpG slightly higher
        curve.D_3p_noncpg[z] = <float>(omega * f_arr[z])
        curve.D_3p_cpg[z] = <float>(omega * f_arr[z] * 1.2)

    # Store fitted parameters
    curve.P_5p_noncpg = <float>P
    curve.P_5p_cpg = <float>(P * 1.2)
    curve.P_3p_noncpg = <float>P
    curve.P_3p_cpg = <float>(P * 1.2)
    curve.lambda_decay = <float>lam
    curve.C_background = <float>C
    curve.epsilon_baseline = <float>epsilon

    ctx.model.curve_fitted = True

    # Compute D_avg (average of positions 0-7 for both ends) for comparison with kaiku
    cdef double D_avg_check = 0.0
    for z in range(8):
        D_avg_check += curve.D_5p_noncpg[z] + curve.D_3p_noncpg[z]
    D_avg_check /= 16.0

    bf_nogil_logf_notime(LOG_TAG, "PMD curve fitted: omega=%.4f lambda=%.4f D_avg=%.4f D(1)=%.3f",
                         omega, lam, D_avg_check, curve.D_5p_noncpg[0])

    return 0


cdef void finalize_pmd_model(PMDGlobalContext* ctx) noexcept nogil:
    """Finalize PMD model, making it read-only for EM and BAM writing."""
    if not ctx:
        return

    if not ctx.model.stats_collected:
        merge_pmd_stats(ctx)

    if not ctx.model.curve_fitted:
        fit_pmd_curve(ctx, NULL)

    ctx.model.finalized = True
    bf_nogil_logf_notime(LOG_TAG, "PMD model finalized (read-only)")


# =============================================================================
# ANI Computation
# =============================================================================

cdef inline float compute_raw_ani(ANISnapshot* snapshot) noexcept nogil:
    """Compute raw ANI from snapshot.

    Returns
    -------
    float
        ANI as percentage (0-100)
    """
    if snapshot.aligned_length == 0:
        return 0.0
    return 100.0 * <float>snapshot.match_count / <float>snapshot.aligned_length


cdef inline float compute_damage_weight(float D_z, float epsilon) noexcept nogil:
    """Compute posterior probability that a mismatch is damage-induced.

    w = D(z) / (D(z) + epsilon)

    Parameters
    ----------
    D_z : float
        Damage probability at position z
    epsilon : float
        Baseline error rate

    Returns
    -------
    float
        Probability in [0, 1]
    """
    cdef float denom = D_z + epsilon
    if denom < 1e-10:
        return 0.0
    return D_z / denom


cdef float compute_corrected_ani(ANISnapshot* snapshot,
                                  PMDCurve* curve,
                                  float epsilon) noexcept nogil:
    """Compute damage-corrected ANI from snapshot.

    Adjusts ANI by treating damage-eligible mismatches as partial matches:
    corrected_matches = matches + sum(w_i) for damage-eligible mismatches

    Parameters
    ----------
    snapshot : ANISnapshot*
        Compact ANI statistics
    curve : PMDCurve*
        Fitted damage curve
    epsilon : float
        Baseline error rate

    Returns
    -------
    float
        Damage-corrected ANI as percentage (0-100)
    """
    if snapshot.aligned_length == 0:
        return 0.0

    cdef float matches = <float>snapshot.match_count
    cdef float damage_contribution = 0.0

    # For C→T mismatches in first 8bp from 5' end
    # We don't have per-position breakdown, so use average D(z) for positions 1-8
    cdef int ct_count = snapshot.ct_5p_count
    cdef int ga_count = snapshot.ga_3p_count

    # Average damage weight for positions 1-8
    cdef float avg_D_5p = 0.0
    cdef float avg_D_3p = 0.0
    cdef int z
    for z in range(8):
        avg_D_5p += curve.D_5p_noncpg[z]
        avg_D_3p += curve.D_3p_noncpg[z]
    avg_D_5p /= 8.0
    avg_D_3p /= 8.0

    cdef float w_5p = compute_damage_weight(avg_D_5p, epsilon)
    cdef float w_3p = compute_damage_weight(avg_D_3p, epsilon)

    damage_contribution = ct_count * w_5p + ga_count * w_3p

    cdef float corrected_matches = matches + damage_contribution
    return 100.0 * corrected_matches / <float>snapshot.aligned_length


cdef ANIResult compute_ani_pair(ANISnapshot* snapshot,
                                 PMDCurve* curve,
                                 float epsilon) noexcept nogil:
    """Compute both raw and corrected ANI.

    Parameters
    ----------
    snapshot : ANISnapshot*
        Compact ANI statistics
    curve : PMDCurve*
        Fitted damage curve (can be NULL for raw-only)
    epsilon : float
        Baseline error rate

    Returns
    -------
    ANIResult
        Structure with raw_ani, corrected_ani, and damage_contribution
    """
    cdef ANIResult result
    result.raw_ani = compute_raw_ani(snapshot)

    if curve and curve.omega > 0.001:
        result.corrected_ani = compute_corrected_ani(snapshot, curve, epsilon)
        result.damage_contribution = result.corrected_ani - result.raw_ani
    else:
        result.corrected_ani = result.raw_ani
        result.damage_contribution = 0.0

    return result


# =============================================================================
# Hierarchical EM: Ancient/Modern Classification
# =============================================================================

cdef int init_hierarchical_em(PMDGlobalContext* ctx,
                               uint32_t reference_count) noexcept nogil:
    """Initialize hierarchical EM state for ancient/modern classification.

    Parameters
    ----------
    ctx : PMDGlobalContext*
        Global PMD context
    reference_count : uint32_t
        Number of references

    Returns
    -------
    int
        0 on success, -1 on allocation failure
    """
    if not ctx:
        return -1

    if ctx.model.hierarchical:
        # Already initialized, reset
        if ctx.model.hierarchical.ref_classes:
            memset(ctx.model.hierarchical.ref_classes, 0,
                   reference_count * sizeof(PMDReferenceClass))
        return 0

    ctx.model.hierarchical = <PMDHierarchicalEM*>calloc(1, sizeof(PMDHierarchicalEM))
    if not ctx.model.hierarchical:
        return -1

    ctx.model.hierarchical.ref_classes = <PMDReferenceClass*>calloc(
        reference_count, sizeof(PMDReferenceClass))
    if not ctx.model.hierarchical.ref_classes:
        free(ctx.model.hierarchical)
        ctx.model.hierarchical = NULL
        return -1

    ctx.model.hierarchical.reference_count = reference_count
    ctx.model.hierarchical.rho = 0.5  # Initial ancient fraction
    ctx.model.hierarchical.rho_prior_alpha = 1.0
    ctx.model.hierarchical.rho_prior_beta = 1.0
    ctx.model.hierarchical.iteration = 0
    ctx.model.hierarchical.converged = False

    ctx.model.hierarchical_enabled = True

    bf_nogil_logf_notime(LOG_TAG, "Hierarchical EM initialized: %u references",
                         reference_count)

    return 0


cdef inline void accumulate_ref_likelihood(PMDHierarchicalEM* hem,
                                            uint32_t ref_idx,
                                            float log_L_ancient,
                                            float log_L_modern) noexcept nogil:
    """Accumulate log-likelihoods for a reference from one read.

    Called during EM E-step for each read-reference assignment.
    """
    if not hem or ref_idx >= hem.reference_count:
        return

    hem.ref_classes[ref_idx].log_L_ancient += log_L_ancient
    hem.ref_classes[ref_idx].log_L_modern += log_L_modern
    hem.ref_classes[ref_idx].read_count += 1


cdef void update_gamma_and_rho(PMDHierarchicalEM* hem) noexcept nogil:
    """Update gamma (per-reference ancient posterior) and rho (global ancient fraction).

    Called at the end of each EM iteration.

    gamma_k = sigma(log(rho) + log_L_ancient_k - log(1-rho) - log_L_modern_k)
    rho_new = (alpha - 1 + sum(gamma_k)) / (alpha + beta - 2 + K)
    """
    if not hem:
        return

    cdef uint32_t k
    cdef float log_rho = log(hem.rho + 1e-10)
    cdef float log_1mrho = log(1.0 - hem.rho + 1e-10)
    cdef float log_odds, gamma_sum = 0.0
    cdef PMDReferenceClass* rc

    for k in range(hem.reference_count):
        rc = &hem.ref_classes[k]
        if rc.read_count == 0:
            rc.gamma = hem.rho  # Use prior for refs with no reads
            continue

        # log odds = log(rho) + log_L_ancient - log(1-rho) - log_L_modern
        log_odds = log_rho + rc.log_L_ancient - log_1mrho - rc.log_L_modern

        # gamma = sigma(log_odds) = 1 / (1 + exp(-log_odds))
        if log_odds > 20.0:
            rc.gamma = 1.0
        elif log_odds < -20.0:
            rc.gamma = 0.0
        else:
            rc.gamma = 1.0 / (1.0 + exp(-log_odds))

        gamma_sum += rc.gamma

    # Update rho with Beta prior
    cdef float K = <float>hem.reference_count
    hem.rho = (hem.rho_prior_alpha - 1.0 + gamma_sum) / \
              (hem.rho_prior_alpha + hem.rho_prior_beta - 2.0 + K)
    hem.rho = fmax(0.01, fmin(hem.rho, 0.99))

    hem.iteration += 1


cdef inline float get_reference_gamma(PMDHierarchicalEM* hem, uint32_t ref_idx) noexcept nogil:
    """Get posterior ancient probability for a reference."""
    if not hem or ref_idx >= hem.reference_count:
        return 0.5
    return hem.ref_classes[ref_idx].gamma


# =============================================================================
# Per-Position Likelihood Computation
# =============================================================================

cdef inline float compute_log_likelihood_ancient(uint8_t ref_base, uint8_t read_base,
                                                  float epsilon, float D_z) noexcept nogil:
    """Compute log P(read_base | ref_base, ancient model with damage D_z).

    For C→T or G→A eligible positions:
      P(T|C) = D_z*(1-epsilon) + (1-D_z)*(epsilon/3)
      P(C|C) = (1-D_z)*(1-epsilon) + D_z*(epsilon/3)

    For other positions: standard error model.
    """
    cdef float p, log_p

    # Check for C→T damage (ref=C, read=T)
    if ref_base == 67 and read_base == 84:  # 'C' and 'T'
        p = D_z * (1.0 - epsilon) + (1.0 - D_z) * (epsilon / 3.0)
        return log(fmax(p, 1e-30))

    # Check for G→A damage (ref=G, read=A)
    if ref_base == 71 and read_base == 65:  # 'G' and 'A'
        p = D_z * (1.0 - epsilon) + (1.0 - D_z) * (epsilon / 3.0)
        return log(fmax(p, 1e-30))

    # Match at C with possible back-mutation
    if ref_base == 67 and read_base == 67:  # C match
        p = (1.0 - D_z) * (1.0 - epsilon) + D_z * (epsilon / 3.0)
        return log(fmax(p, 1e-30))

    # Match at G with possible back-mutation
    if ref_base == 71 and read_base == 71:  # G match
        p = (1.0 - D_z) * (1.0 - epsilon) + D_z * (epsilon / 3.0)
        return log(fmax(p, 1e-30))

    # Standard match
    if ref_base == read_base:
        return log(1.0 - epsilon)

    # Standard mismatch
    return log(epsilon / 3.0)


cdef inline float compute_log_likelihood_modern(uint8_t ref_base, uint8_t read_base,
                                                 float epsilon) noexcept nogil:
    """Compute log P(read_base | ref_base, modern model - no damage).

    Standard sequencing error model:
      P(match) = 1 - epsilon
      P(specific mismatch) = epsilon / 3
    """
    if ref_base == read_base:
        return log(1.0 - epsilon)
    return log(epsilon / 3.0)


# =============================================================================
# Utility Functions
# =============================================================================

cdef inline float get_damage_prob(PMDCurve* curve, int end_type, int context_type,
                                   int position) noexcept nogil:
    """Get damage probability D(z) for a specific position.

    Parameters
    ----------
    curve : PMDCurve*
        Fitted damage curve
    end_type : int
        PMD_END_5P or PMD_END_3P
    context_type : int
        PMD_CTX_NONCPG or PMD_CTX_CPG
    position : int
        1-based position from end (clamped to 1-20)

    Returns
    -------
    float
        Damage probability D(z)
    """
    if not curve:
        return 0.0

    cdef int idx = position - 1
    if idx < 0:
        idx = 0
    elif idx >= PMD_MAX_POSITION:
        idx = PMD_MAX_POSITION - 1

    if end_type == PMD_END_5P:
        if context_type == PMD_CTX_CPG:
            return curve.D_5p_cpg[idx]
        return curve.D_5p_noncpg[idx]
    else:
        if context_type == PMD_CTX_CPG:
            return curve.D_3p_cpg[idx]
        return curve.D_3p_noncpg[idx]


# =============================================================================
# Apply PMD Corrections to All Alignments
# =============================================================================

cdef int64_t apply_pmd_corrections_to_pool(MemoryPool* pool,
                                            PMDCurve* curve,
                                            float min_ani_threshold,
                                            float epsilon,
                                            int num_threads) noexcept nogil:
    """Apply damage-corrected ANI to all alignments in memory pool.

    Computes damage-corrected ANI for each alignment using the fitted PMD curve
    and marks alignments that pass the ANI threshold.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments with ANI snapshots
    curve : PMDCurve*
        Fitted PMD damage curve
    min_ani_threshold : float
        Minimum ANI threshold (e.g., 90.0 for 90%)
    epsilon : float
        Baseline sequencing error rate (e.g., 0.01)
    num_threads : int
        Number of threads for parallel processing

    Returns
    -------
    int64_t
        Number of alignments that pass the corrected ANI threshold
    """
    if not pool or not pool.alignments or pool.alignment_count == 0:
        return 0

    cdef int64_t i
    cdef int64_t passed_count = 0
    cdef int64_t alignment_count = pool.alignment_count
    cdef Alignment* aln
    cdef ANISnapshot snapshot
    cdef float raw_ani, corrected_ani
    cdef float D_5p, D_3p
    cdef float damage_correction
    cdef int ct_pos, ga_pos

    # If no curve available, use raw ANI
    cdef bint has_curve = (curve != NULL and curve.omega > 0.001)

    # Precompute log probability factors once for hierarchical EM (hoisted from inner loop)
    cdef double D_avg = 0.0
    cdef double log_p_damage_anc, log_p_survive_anc
    cdef double log_1m_eps, log_eps_over_3
    cdef double observed_damage, survived, other_mm
    cdef int actual_opp

    if has_curve:
        D_avg = (curve.D_5p_noncpg[0] + curve.D_3p_noncpg[0]) / 2.0
    else:
        D_avg = 0.1  # default fallback

    # Precomputed log factors (computed ONCE, used for all alignments)
    log_p_damage_anc = log(fmax(D_avg + (1.0 - D_avg) * epsilon / 3.0, 1e-15))
    log_p_survive_anc = log(fmax((1.0 - D_avg) * (1.0 - epsilon / 3.0), 1e-15))
    log_1m_eps = log(fmax(1.0 - epsilon, 1e-15))
    log_eps_over_3 = log(fmax(epsilon / 3.0, 1e-15))

    bf_nogil_logf_notime(LOG_TAG, "Applying PMD corrections to %ld alignments (threshold=%.1f%%)",
                         alignment_count, min_ani_threshold)
    bf_nogil_logf_notime(LOG_TAG, "Precomputed EM factors: D_avg=%.4f log_p_dmg=%.4f log_p_surv=%.4f",
                         D_avg, log_p_damage_anc, log_p_survive_anc)

    # Process alignments in parallel
    for i in prange(alignment_count, nogil=True, num_threads=num_threads, schedule='static'):
        aln = &pool.alignments[i]

        # Compute raw ANI
        if aln.aligned_length > 0:
            raw_ani = (<float>aln.match_count / <float>aln.aligned_length) * 100.0
        else:
            raw_ani = 0.0
            aln.corrected_ani = 0.0
            aln.passes_ani_filter = 0
            aln.log_L_anc = -1000.0  # very low likelihood
            aln.log_L_mod = -1000.0
            continue

        if has_curve:
            # Compute damage-corrected ANI
            # Build ANI snapshot from alignment fields
            snapshot.aligned_length = aln.aligned_length
            snapshot.match_count = aln.match_count
            snapshot.ct_5p_count = aln.ct_5p_count
            snapshot.ga_3p_count = aln.ga_3p_count
            snapshot.other_mm_count = 0
            snapshot.flags = 0
            snapshot.c_at_5p_count = aln.c_at_5p_count
            snapshot.g_at_3p_count = aln.g_at_3p_count

            # Compute corrected ANI using the damage curve
            corrected_ani = compute_corrected_ani(&snapshot, curve, epsilon)
        else:
            # No curve - use raw ANI
            corrected_ani = raw_ani

        # Store corrected ANI
        aln.corrected_ani = corrected_ani

        # Check if passes threshold
        if corrected_ani >= min_ani_threshold:
            aln.passes_ani_filter = 1
        else:
            aln.passes_ani_filter = 0

        # Precompute log-likelihoods for hierarchical EM E-step
        # Using precomputed log factors - no log() calls in this hot path
        observed_damage = <double>(aln.ct_5p_count + aln.ga_3p_count)
        actual_opp = aln.c_at_5p_count + aln.g_at_3p_count
        survived = fmax(0.0, <double>actual_opp - observed_damage)
        other_mm = <double>(aln.aligned_length - aln.match_count) - observed_damage
        if other_mm < 0.0:
            other_mm = 0.0

        # log_L_anc = obs_damage * log_p_damage + survived * log_p_survive + matches * log(1-eps) + other_mm * log(eps/3)
        aln.log_L_anc = <float>(
            observed_damage * log_p_damage_anc +
            survived * log_p_survive_anc +
            <double>aln.match_count * log_1m_eps +
            other_mm * log_eps_over_3
        )

        # log_L_mod = matches * log(1-eps) + all_mm * log(eps/3) (all mismatches are errors)
        aln.log_L_mod = <float>(
            <double>aln.match_count * log_1m_eps +
            <double>(aln.aligned_length - aln.match_count) * log_eps_over_3
        )

    # Count passed alignments (serial for accuracy)
    for i in range(alignment_count):
        if pool.alignments[i].passes_ani_filter == 1:
            passed_count += 1

    bf_nogil_logf_notime(LOG_TAG, "PMD correction complete: %ld/%ld passed (%.1f%%)",
                         passed_count, alignment_count,
                         100.0 * <float>passed_count / <float>alignment_count if alignment_count > 0 else 0.0)

    return passed_count


cdef void update_alignment_scores_with_pmd(MemoryPool* pool,
                                            PMDCurve* curve,
                                            float damage_score_bonus,
                                            int num_threads) noexcept nogil:
    """Update alignment scores to incorporate damage correction.

    For ancient DNA, alignments with damage patterns should get a score boost
    because the mismatches are due to damage, not divergence.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments
    curve : PMDCurve*
        Fitted PMD damage curve
    damage_score_bonus : float
        Score bonus per damage-corrected mismatch (e.g., 0.5)
    num_threads : int
        Number of threads for parallel processing
    """
    if not pool or not pool.alignments or pool.alignment_count == 0:
        return
    if not curve or curve.omega < 0.001:
        return

    cdef int64_t i
    cdef int64_t alignment_count = pool.alignment_count
    cdef Alignment* aln
    cdef float score_adjustment
    cdef float D_5p_avg, D_3p_avg

    # Use average damage at position 1-3 for score adjustment
    D_5p_avg = (curve.D_5p_noncpg[0] + curve.D_5p_noncpg[1] + curve.D_5p_noncpg[2]) / 3.0
    D_3p_avg = (curve.D_3p_noncpg[0] + curve.D_3p_noncpg[1] + curve.D_3p_noncpg[2]) / 3.0

    bf_nogil_logf_notime(LOG_TAG, "Updating alignment scores with PMD bonus (D_5p=%.3f, D_3p=%.3f)",
                         D_5p_avg, D_3p_avg)

    for i in prange(alignment_count, nogil=True, num_threads=num_threads, schedule='static'):
        aln = &pool.alignments[i]

        # Calculate score adjustment based on damage-like mismatches
        # Each C->T at 5' or G->A at 3' gets partial credit back
        score_adjustment = 0.0

        if aln.ct_5p_count > 0:
            # C->T mismatches at 5' end: weight by damage probability
            score_adjustment += <float>aln.ct_5p_count * D_5p_avg * damage_score_bonus

        if aln.ga_3p_count > 0:
            # G->A mismatches at 3' end: weight by damage probability
            score_adjustment += <float>aln.ga_3p_count * D_3p_avg * damage_score_bonus

        # Add bonus to alignment score (scores are typically negative, so add positive)
        aln.alignment_score += score_adjustment


# =============================================================================
# Hierarchical EM: Per-alignment ancient/modern likelihood computation
# =============================================================================

cdef inline double compute_alignment_log_L_ancient(Alignment* aln,
                                                     float D_avg_5p,
                                                     float D_avg_3p,
                                                     float epsilon) noexcept nogil:
    """Compute log-likelihood of alignment under ancient model (with damage).

    Uses aggregated counts from ANI snapshot:
    - ct_5p_count: C→T mismatches at 5' end (positions 1-8)
    - ga_3p_count: G→A mismatches at 3' end (positions 1-8)
    - match_count: exact matches
    - aligned_length: total aligned bases

    Under ancient model:
    - C→T at 5' and G→A at 3' are expected (probability D_avg)
    - Other positions follow error model (probability ε)

    Parameters
    ----------
    aln : Alignment*
        Alignment with ANI snapshot fields
    D_avg_5p : float
        Average damage probability for 5' positions 1-8
    D_avg_3p : float
        Average damage probability for 3' positions 1-8
    epsilon : float
        Sequencing error rate (~0.01)

    Returns
    -------
    double
        Log-likelihood under ancient model
    """
    cdef double log_L = 0.0
    cdef uint16_t aligned = aln.aligned_length
    cdef uint16_t matches = aln.match_count
    cdef uint8_t ct_5p = aln.ct_5p_count
    cdef uint8_t ga_3p = aln.ga_3p_count
    cdef int other_mm

    if aligned == 0:
        return -1e10  # Invalid alignment

    # Other mismatches (non-damage)
    other_mm = aligned - matches - ct_5p - ga_3p
    if other_mm < 0:
        other_mm = 0

    # Ancient model likelihood:
    # - C→T at 5' end: P(T|C) = D × (1-ε) + (1-D) × (ε/3) ≈ D for typical D >> ε
    # - G→A at 3' end: same as above
    # - Matches: P(match) = 1 - ε
    # - Other mismatches: P(mm) = ε/3

    # Simplified: log P(damage mm) ≈ log(D_avg + ε/3) for numerical stability
    cdef double log_p_damage_5p = log(fmax(D_avg_5p + epsilon / 3.0, 1e-10))
    cdef double log_p_damage_3p = log(fmax(D_avg_3p + epsilon / 3.0, 1e-10))
    cdef double log_p_match = log(fmax(1.0 - epsilon, 1e-10))
    cdef double log_p_error = log(fmax(epsilon / 3.0, 1e-10))

    log_L += ct_5p * log_p_damage_5p      # 5' C→T damage
    log_L += ga_3p * log_p_damage_3p      # 3' G→A damage
    log_L += matches * log_p_match        # Exact matches
    log_L += other_mm * log_p_error       # Other mismatches (sequencing error)

    return log_L


cdef inline double compute_alignment_log_L_modern(Alignment* aln,
                                                    float epsilon) noexcept nogil:
    """Compute log-likelihood of alignment under modern model (no damage).

    Under modern model:
    - ALL mismatches are sequencing errors (probability ε/3)
    - Matches have probability 1-ε

    Parameters
    ----------
    aln : Alignment*
        Alignment with ANI snapshot fields
    epsilon : float
        Sequencing error rate (~0.01)

    Returns
    -------
    double
        Log-likelihood under modern model
    """
    cdef double log_L = 0.0
    cdef uint16_t aligned = aln.aligned_length
    cdef uint16_t matches = aln.match_count
    cdef int total_mm

    if aligned == 0:
        return -1e10  # Invalid alignment

    # Total mismatches under modern model
    total_mm = aligned - matches
    if total_mm < 0:
        total_mm = 0

    cdef double log_p_match = log(fmax(1.0 - epsilon, 1e-10))
    cdef double log_p_error = log(fmax(epsilon / 3.0, 1e-10))

    log_L += matches * log_p_match       # Exact matches
    log_L += total_mm * log_p_error      # ALL mismatches are errors

    return log_L


cdef void init_hierarchical_em_pool(MemoryPool* pool, PMDCurve* curve,
                                     float epsilon) noexcept nogil:
    """Initialize hierarchical EM state in memory pool.

    Allocates arrays for γ, η, S_anc, S_mod and precomputes D_avg values.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool to initialize
    curve : PMDCurve*
        Fitted PMD curve (for D_avg computation)
    epsilon : float
        Sequencing error rate
    """
    cdef uint32_t i
    cdef float D_sum_5p = 0.0, D_sum_3p = 0.0

    # Allocate arrays if not already allocated
    if pool.gamma_values == NULL:
        pool.gamma_values = <double*>calloc(pool.reference_count, sizeof(double))
    if pool.eta_values == NULL:
        pool.eta_values = <double*>calloc(pool.reference_count, sizeof(double))
    if pool.S_anc_accum == NULL:
        pool.S_anc_accum = <double*>calloc(pool.reference_count, sizeof(double))
    if pool.S_mod_accum == NULL:
        pool.S_mod_accum = <double*>calloc(pool.reference_count, sizeof(double))

    # Initialize γ to 0.5 (uninformative prior) and η to 0 (logit(0.5))
    for i in range(pool.reference_count):
        pool.gamma_values[i] = 0.5
        pool.eta_values[i] = 0.0  # logit(0.5) = 0

    # Initialize ρ to 0.5 and ζ to 0
    pool.rho_ancient = 0.5
    pool.zeta_ancient = 0.0  # logit(0.5) = 0
    pool.rho_prior_alpha = 1.0  # Flat Beta prior
    pool.rho_prior_beta = 1.0

    # Precompute average D(z) for positions 1-8
    if curve != NULL:
        for i in range(8):
            D_sum_5p += curve.D_5p_noncpg[i]
            D_sum_3p += curve.D_3p_noncpg[i]
        pool.D_avg_5p = D_sum_5p / 8.0
        pool.D_avg_3p = D_sum_3p / 8.0
    else:
        pool.D_avg_5p = 0.1  # Default if no curve
        pool.D_avg_3p = 0.1

    pool.epsilon_error = epsilon
    pool.hierarchical_em_enabled = True

    bf_nogil_logf_notime(LOG_TAG,
        "Hierarchical EM initialized: refs=%u D_avg_5p=%.3f D_avg_3p=%.3f epsilon=%.4f",
        pool.reference_count, pool.D_avg_5p, pool.D_avg_3p, epsilon)


cdef void free_hierarchical_em_pool(MemoryPool* pool) noexcept nogil:
    """Free hierarchical EM arrays from memory pool."""
    if pool.gamma_values != NULL:
        free(pool.gamma_values)
        pool.gamma_values = NULL
    if pool.eta_values != NULL:
        free(pool.eta_values)
        pool.eta_values = NULL
    if pool.S_anc_accum != NULL:
        free(pool.S_anc_accum)
        pool.S_anc_accum = NULL
    if pool.S_mod_accum != NULL:
        free(pool.S_mod_accum)
        pool.S_mod_accum = NULL
    pool.hierarchical_em_enabled = False


# =============================================================================
# Python-accessible functions
# =============================================================================

def init_hierarchical_em_py(uintptr_t pool_ptr, uintptr_t curve_ptr, float epsilon):
    """Python wrapper for initializing hierarchical EM state."""
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef PMDCurve* curve = <PMDCurve*>curve_ptr if curve_ptr != 0 else NULL
    with nogil:
        init_hierarchical_em_pool(pool, curve, epsilon)


def apply_pmd_corrections_py(uintptr_t pool_ptr, uintptr_t curve_ptr,
                              float min_ani_threshold, float epsilon, int num_threads):
    """Python wrapper for applying PMD corrections to alignments.

    Returns the number of alignments that pass the corrected ANI threshold.
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef PMDCurve* curve = <PMDCurve*>curve_ptr
    return apply_pmd_corrections_to_pool(pool, curve, min_ani_threshold, epsilon, num_threads)


def apply_raw_ani_filter_py(uintptr_t pool_ptr, float min_ani_threshold, int num_threads):
    """Apply raw ANI filtering when PMD is disabled.

    Sets passes_ani_filter based on raw ANI (match_count/aligned_length).
    Returns the number of alignments that pass the threshold.
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef int64_t i
    cdef int64_t passed_count = 0
    cdef int64_t alignment_count = pool.alignment_count
    cdef Alignment* aln
    cdef float raw_ani

    bf_nogil_logf_notime(LOG_TAG, "Applying raw ANI filter to %lld alignments (threshold=%.1f%%)",
                         <long long>alignment_count, min_ani_threshold)

    for i in prange(alignment_count, nogil=True, num_threads=num_threads, schedule='static'):
        aln = &pool.alignments[i]

        if aln.aligned_length > 0:
            raw_ani = (<float>aln.match_count / <float>aln.aligned_length) * 100.0
        else:
            raw_ani = 0.0

        aln.corrected_ani = raw_ani  # No correction, use raw

        if raw_ani >= min_ani_threshold:
            aln.passes_ani_filter = 1
        else:
            aln.passes_ani_filter = 0

    # Count passed (separate loop to avoid race condition)
    for i in range(alignment_count):
        if pool.alignments[i].passes_ani_filter == 1:
            passed_count += 1

    bf_nogil_logf_notime(LOG_TAG, "Raw ANI filter complete: %lld/%lld passed (%.1f%%)",
                         <long long>passed_count, <long long>alignment_count,
                         100.0 * <double>passed_count / <double>alignment_count if alignment_count > 0 else 0.0)

    return passed_count


def create_pmd_context_py(int num_threads, bint is_single_stranded, bint enable_hierarchical):
    """Python wrapper for creating PMD context."""
    cdef PMDGlobalContext* ctx = create_pmd_context(num_threads, is_single_stranded,
                                                     enable_hierarchical)
    if not ctx:
        raise MemoryError("Failed to allocate PMD context")
    return <uintptr_t>ctx


def destroy_pmd_context_py(uintptr_t ctx_ptr):
    """Python wrapper for destroying PMD context."""
    destroy_pmd_context(<PMDGlobalContext*>ctx_ptr)


def finalize_pmd_model_py(uintptr_t ctx_ptr):
    """Python wrapper for finalizing PMD model."""
    cdef PMDGlobalContext* ctx = <PMDGlobalContext*>ctx_ptr
    finalize_pmd_model(ctx)
    return {
        'omega': ctx.model.curve.omega,
        'D_1': ctx.model.curve.D_5p_noncpg[0],
        'D_5': ctx.model.curve.D_5p_noncpg[4],
        'D_10': ctx.model.curve.D_5p_noncpg[9],
        'stats_collected': ctx.model.stats_collected,
        'curve_fitted': ctx.model.curve_fitted,
        'total_alignments': ctx.model.stats.total_alignments,
    }


def get_pmd_stats_py(uintptr_t ctx_ptr):
    """Python wrapper for getting PMD statistics."""
    cdef PMDGlobalContext* ctx = <PMDGlobalContext*>ctx_ptr
    if not ctx:
        return None

    cdef PMDStatsGlobal* stats = &ctx.model.stats

    return {
        'total_alignments': stats.total_alignments,
        'total_bases': stats.total_bases,
        'n_5p_noncpg': [stats.n_5p_noncpg[i] for i in range(20)],
        'k_5p_noncpg': [stats.k_5p_noncpg[i] for i in range(20)],
        'n_3p_noncpg': [stats.n_3p_noncpg[i] for i in range(20)],
        'k_3p_noncpg': [stats.k_3p_noncpg[i] for i in range(20)],
    }
