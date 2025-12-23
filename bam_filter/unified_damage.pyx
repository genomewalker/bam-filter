# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# -*- coding: utf-8 -*-
"""Unified ancient DNA damage model.

Single model for all damage-related tasks:
- ANI correction
- Ancient/modern classification
- Taxonomic profiling

Model:
    p_{r,e}(z) = b_{r,e} + δ_{r,e} · exp(-(z-1)/τ_s)

Where:
- τ_s: sample-level decay (fitted with prior)
- δ_{r,e}: reference-level amplitude with Gamma shrinkage
- b_{r,e}: baseline from interior positions

Inference: Poisson EM with closed-form updates
"""

from libc.math cimport exp, log, log1p, fabs, fmax, fmin, sqrt, isnan, isinf
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memset
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t, uintptr_t

from cython.parallel cimport prange

from libc.string cimport memcpy

from bam_filter.unified_damage cimport *
from bam_filter.processor_pmd cimport RefDamageStats, PMDCurve
from bam_filter.processor cimport MemoryPool
from bam_filter.squarem cimport (
    SquaremState, SquaremConfig, SquaremBoxCtx,
    squarem_bytes, squarem_config_default, squarem_state_init,
    squarem_state_free, squarem_step, squarem_project_box
)

# Logging
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef const char* LOG_TAG = b"UNIFIED_DMG"

# Numerical constants
cdef double EPS = 1e-15
cdef double LOG_EPS = -34.5  # log(1e-15)
cdef double MAX_EXP_ARG = 50.0


# =============================================================================
# Decay function
# =============================================================================

cdef inline double decay_func(int z, double tau) noexcept nogil:
    """Compute exp(-(z-1)/tau) with numerical safeguards.

    z is 1-based position from end.
    """
    cdef double t = -(<double>(z - 1)) / tau
    if t < -MAX_EXP_ARG:
        return 0.0
    return exp(t)


cdef double compute_damage_prob(double delta, double tau, int position) noexcept nogil:
    """Compute damage probability at position (1-based)."""
    return delta * decay_func(position, tau)


# =============================================================================
# Context management
# =============================================================================

cdef UnifiedDamageContext* create_damage_context(
    uint32_t n_refs,
    bint is_single_stranded,
    bint mask_cpg
) noexcept nogil:
    """Create and initialize unified damage context."""
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>calloc(1, sizeof(UnifiedDamageContext))
    if ctx == NULL:
        return NULL

    # Allocate per-reference arrays
    ctx.ref_counts = <RefDamageCounts*>calloc(n_refs, sizeof(RefDamageCounts))
    ctx.ref_params = <RefDamageParams*>calloc(n_refs, sizeof(RefDamageParams))

    if ctx.ref_counts == NULL or ctx.ref_params == NULL:
        if ctx.ref_counts != NULL:
            free(ctx.ref_counts)
        if ctx.ref_params != NULL:
            free(ctx.ref_params)
        free(ctx)
        return NULL

    ctx.n_refs = n_refs
    ctx.is_single_stranded = is_single_stranded
    ctx.mask_cpg = mask_cpg
    ctx.use_controls = True
    ctx.owns_memory = True

    # Default convergence settings
    ctx.max_iterations = 50
    ctx.tol_tau = 1e-4
    ctx.tol_delta = 1e-3

    # Initialize sample parameters with priors
    ctx.sample.tau = 5.0
    ctx.sample.tau_prior_mean = log(5.0)  # log-normal prior
    ctx.sample.tau_prior_sd = 0.5
    ctx.sample.tau_min = 1.0
    ctx.sample.tau_max = 30.0
    ctx.sample.alpha = 10.0  # Gamma shrinkage strength
    ctx.sample.mu_5p = 0.02  # initial amplitude guess
    ctx.sample.mu_3p = 0.02
    ctx.sample.iteration = 0
    ctx.sample.converged = False

    bf_nogil_logf_notime(LOG_TAG, "Created context: refs=%u ss=%d cpg_mask=%d",
                         n_refs, is_single_stranded, mask_cpg)

    return ctx


cdef void destroy_damage_context(UnifiedDamageContext* ctx) noexcept nogil:
    """Free all memory associated with context."""
    if ctx == NULL:
        return

    if ctx.owns_memory:
        if ctx.ref_counts != NULL:
            free(ctx.ref_counts)
        if ctx.ref_params != NULL:
            free(ctx.ref_params)

    free(ctx)


cdef void reset_damage_counts(UnifiedDamageContext* ctx) noexcept nogil:
    """Reset all per-reference counts to zero."""
    if ctx == NULL or ctx.ref_counts == NULL:
        return

    memset(ctx.ref_counts, 0, ctx.n_refs * sizeof(RefDamageCounts))


cdef void set_tau_from_pmd(
    UnifiedDamageContext* ctx,
    double tau_value,
    bint fix_tau
) noexcept nogil:
    """Set tau from PMD-derived estimate.

    Parameters
    ----------
    ctx : UnifiedDamageContext*
        Damage model context
    tau_value : double
        Tau value (already converted from PMD lambda: tau = 1/lambda)
    fix_tau : bint
        If True, set very tight prior to prevent re-estimation
    """
    if ctx == NULL:
        return

    cdef double tau = tau_value

    # Validate and clamp
    if tau < 1.0:
        tau = 5.0  # Fallback to prior if invalid

    tau = fmax(ctx.sample.tau_min, fmin(ctx.sample.tau_max, tau))

    ctx.sample.tau = tau
    ctx.sample.tau_prior_mean = log(tau)

    if fix_tau:
        ctx.sample.tau_prior_sd = 0.05  # Tight but allows small adjustment
    else:
        ctx.sample.tau_prior_sd = 0.3

    bf_nogil_logf_notime(LOG_TAG, "Set tau from PMD: tau=%.2f (fix=%d, prior_sd=%.2f)",
                         tau, fix_tau, ctx.sample.tau_prior_sd)


# =============================================================================
# Count accumulation (called during batch processing)
# =============================================================================

cdef void accumulate_damage_position(
    RefDamageCounts* counts,
    int end_type,
    int is_control,
    int position,
    double weight,
    bint is_mismatch
) noexcept nogil:
    """Accumulate a single position observation.

    Parameters
    ----------
    counts : RefDamageCounts*
        Per-reference count accumulator
    end_type : int
        0 = 5' end, 1 = 3' end
    is_control : int
        0 = damage channel (C→T at 5', G→A at 3')
        1 = control channel (G→A at 5', C→T at 3')
    position : int
        1-based position from end (1-20)
    weight : double
        φ weight from EM (1.0 if unweighted)
    is_mismatch : bint
        True if damage-type mismatch observed
    """
    if position < 1 or position > DAMAGE_MAX_POSITION:
        return

    cdef int idx = position - 1  # convert to 0-based

    if end_type == 0:  # 5' end
        if is_control == 0:  # C→T channel
            counts.n_5p[idx] += weight
            if is_mismatch:
                counts.k_5p[idx] += weight
        else:  # G→A control
            counts.n_5p_ctrl[idx] += weight
            if is_mismatch:
                counts.k_5p_ctrl[idx] += weight
    else:  # 3' end
        if is_control == 0:  # G→A channel
            counts.n_3p[idx] += weight
            if is_mismatch:
                counts.k_3p[idx] += weight
        else:  # C→T control
            counts.n_3p_ctrl[idx] += weight
            if is_mismatch:
                counts.k_3p_ctrl[idx] += weight


# =============================================================================
# Baseline estimation from interior positions
# =============================================================================

cdef void estimate_baselines(UnifiedDamageContext* ctx) noexcept nogil:
    """Estimate baseline rates from interior positions (z >= 10).

    Also computes sample-level global baselines as fallback.
    """
    cdef uint32_t r
    cdef int z
    cdef double k_sum, n_sum
    cdef double k_global_5p = 0.0, n_global_5p = 0.0
    cdef double k_global_3p = 0.0, n_global_3p = 0.0
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        # 5' baseline from interior
        k_sum = 0.0
        n_sum = 0.0
        for z in range(DAMAGE_INTERIOR_START - 1, DAMAGE_MAX_POSITION):
            k_sum += counts.k_5p[z]
            n_sum += counts.n_5p[z]

        if n_sum > 0:
            params.baseline_5p = fmax(0.0, fmin(1.0, k_sum / n_sum))
        else:
            params.baseline_5p = 0.01  # default

        k_global_5p += k_sum
        n_global_5p += n_sum

        # 3' baseline from interior
        k_sum = 0.0
        n_sum = 0.0
        for z in range(DAMAGE_INTERIOR_START - 1, DAMAGE_MAX_POSITION):
            k_sum += counts.k_3p[z]
            n_sum += counts.n_3p[z]

        if n_sum > 0:
            params.baseline_3p = fmax(0.0, fmin(1.0, k_sum / n_sum))
        else:
            params.baseline_3p = 0.01

        k_global_3p += k_sum
        n_global_3p += n_sum

    # Global baselines (fallback for refs with no interior data)
    if n_global_5p > 0:
        ctx.sample.baseline_5p_global = k_global_5p / n_global_5p
    else:
        ctx.sample.baseline_5p_global = 0.01

    if n_global_3p > 0:
        ctx.sample.baseline_3p_global = k_global_3p / n_global_3p
    else:
        ctx.sample.baseline_3p_global = 0.01

    # Fill in refs with no interior data
    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        n_sum = 0.0
        for z in range(DAMAGE_INTERIOR_START - 1, DAMAGE_MAX_POSITION):
            n_sum += counts.n_5p[z]
        if n_sum < 1.0:
            params.baseline_5p = ctx.sample.baseline_5p_global

        n_sum = 0.0
        for z in range(DAMAGE_INTERIOR_START - 1, DAMAGE_MAX_POSITION):
            n_sum += counts.n_3p[z]
        if n_sum < 1.0:
            params.baseline_3p = ctx.sample.baseline_3p_global

    bf_nogil_logf_notime(LOG_TAG, "Baselines estimated: global_5p=%.4f global_3p=%.4f",
                         ctx.sample.baseline_5p_global, ctx.sample.baseline_3p_global)


# =============================================================================
# Parameter initialization
# =============================================================================

cdef void initialize_parameters(UnifiedDamageContext* ctx) noexcept nogil:
    """Initialize δ_r from data before EM iterations."""
    cdef uint32_t r
    cdef int z
    cdef double excess, weighted_opp, g
    cdef double delta_init
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params
    cdef double tau = ctx.sample.tau

    cdef double delta_sum_5p = 0.0
    cdef double delta_sum_3p = 0.0
    cdef uint32_t n_usable = 0
    cdef double total_opp

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        # Initialize δ_5p: excess mismatches over baseline
        excess = 0.0
        weighted_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            excess += fmax(0.0, counts.k_5p[z] - counts.n_5p[z] * params.baseline_5p)
            weighted_opp += counts.n_5p[z] * g

        if weighted_opp > EPS:
            params.delta_5p = fmax(0.0, excess / weighted_opp)
        else:
            params.delta_5p = 0.0

        # Initialize δ_3p
        excess = 0.0
        weighted_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            excess += fmax(0.0, counts.k_3p[z] - counts.n_3p[z] * params.baseline_3p)
            weighted_opp += counts.n_3p[z] * g

        if weighted_opp > EPS:
            params.delta_3p = fmax(0.0, excess / weighted_opp)
        else:
            params.delta_3p = 0.0

        # Track for global mean
        total_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            total_opp += counts.n_5p[z] + counts.n_3p[z]

        if total_opp >= 50:  # minimum for usable ref
            delta_sum_5p += params.delta_5p
            delta_sum_3p += params.delta_3p
            n_usable += 1

    # Set global means (shrinkage targets)
    if n_usable > 0:
        ctx.sample.mu_5p = delta_sum_5p / n_usable
        ctx.sample.mu_3p = delta_sum_3p / n_usable
    else:
        ctx.sample.mu_5p = 0.02
        ctx.sample.mu_3p = 0.02

    # Ensure positive
    ctx.sample.mu_5p = fmax(0.001, ctx.sample.mu_5p)
    ctx.sample.mu_3p = fmax(0.001, ctx.sample.mu_3p)

    bf_nogil_logf_notime(LOG_TAG, "Initialized: usable_refs=%u mu_5p=%.4f mu_3p=%.4f",
                         n_usable, ctx.sample.mu_5p, ctx.sample.mu_3p)


# =============================================================================
# EM E-step: compute expected damage counts
# =============================================================================

cdef void em_e_step(UnifiedDamageContext* ctx) noexcept nogil:
    """E-step: compute expected damage counts ŷ[z] for each reference."""
    cdef uint32_t r
    cdef int z
    cdef double g, w, denom
    cdef double delta, baseline, k, n
    cdef double Y_5p, Y_3p, B_5p, B_3p
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params
    cdef double tau = ctx.sample.tau

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        Y_5p = 0.0
        B_5p = 0.0
        Y_3p = 0.0
        B_3p = 0.0

        # 5' end
        delta = params.delta_5p
        baseline = params.baseline_5p
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            k = counts.k_5p[z]
            n = counts.n_5p[z]

            if n < EPS or k < EPS:
                continue

            denom = baseline + delta * g + EPS
            w = (delta * g) / denom

            Y_5p += k * w
            B_5p += k * (1.0 - w)

        # 3' end
        delta = params.delta_3p
        baseline = params.baseline_3p
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            k = counts.k_3p[z]
            n = counts.n_3p[z]

            if n < EPS or k < EPS:
                continue

            denom = baseline + delta * g + EPS
            w = (delta * g) / denom

            Y_3p += k * w
            B_3p += k * (1.0 - w)

        params.Y_5p = Y_5p
        params.Y_3p = Y_3p
        params.B_5p = B_5p
        params.B_3p = B_3p


# =============================================================================
# EM M-step: update δ_r with Gamma shrinkage
# =============================================================================

cdef void em_m_step_delta(UnifiedDamageContext* ctx) noexcept nogil:
    """M-step: update δ_r using posterior mean with coverage-aware Gamma shrinkage.

    Coverage-aware shrinkage: α_r = α_base / sqrt(1 + N_r / N_scale)
    - Low-coverage refs (small N_r): α_r ≈ α_base → strong shrinkage to global mean
    - High-coverage refs (large N_r): α_r → 0 → data dominates

    This empirical-Bayes approach prevents over-shrinking well-sampled references
    while maintaining regularization for sparse data.
    """
    cdef uint32_t r
    cdef int z
    cdef double A_r, g
    cdef double alpha_base = ctx.sample.alpha
    cdef double alpha_r  # per-reference adaptive alpha
    cdef double mu_5p = ctx.sample.mu_5p
    cdef double mu_3p = ctx.sample.mu_3p
    cdef double tau = ctx.sample.tau
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params

    cdef double delta_sum_5p = 0.0
    cdef double delta_sum_3p = 0.0
    cdef double weight_sum_5p = 0.0
    cdef double weight_sum_3p = 0.0
    cdef uint32_t n_counted = 0
    cdef double total_opp
    cdef double N_r  # effective coverage (total_weight)

    # N_scale: coverage at which shrinkage is halved (~sqrt(2) factor)
    # 100 reads is a reasonable threshold for aDNA damage estimation
    cdef double N_scale = 100.0

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        # Compute coverage-aware alpha
        N_r = counts.total_weight
        alpha_r = alpha_base / sqrt(1.0 + N_r / N_scale)

        # Update δ_5p: posterior mean = (α_r + Y) / (α_r/μ + A)
        A_r = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            A_r += counts.n_5p[z] * g

        if A_r > EPS:
            params.delta_5p = (alpha_r + params.Y_5p) / (alpha_r / mu_5p + A_r)
        else:
            params.delta_5p = 0.0

        # Update δ_3p
        A_r = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)
            A_r += counts.n_3p[z] * g

        if A_r > EPS:
            params.delta_3p = (alpha_r + params.Y_3p) / (alpha_r / mu_3p + A_r)
        else:
            params.delta_3p = 0.0

        # Clamp to valid range
        params.delta_5p = fmax(0.0, fmin(1.0, params.delta_5p))
        params.delta_3p = fmax(0.0, fmin(1.0, params.delta_3p))

        # Track for global mean update (coverage-weighted)
        total_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            total_opp += counts.n_5p[z] + counts.n_3p[z]

        if total_opp >= 10:
            # Weight by effective sample size for robust mean estimation
            delta_sum_5p += params.delta_5p * N_r
            delta_sum_3p += params.delta_3p * N_r
            weight_sum_5p += N_r
            weight_sum_3p += N_r
            n_counted += 1

    # Update global means (coverage-weighted empirical Bayes)
    if weight_sum_5p > EPS:
        ctx.sample.mu_5p = fmax(0.001, delta_sum_5p / weight_sum_5p)
    if weight_sum_3p > EPS:
        ctx.sample.mu_3p = fmax(0.001, delta_sum_3p / weight_sum_3p)


# =============================================================================
# EM M-step: update τ (1D optimization)
# =============================================================================

cdef double tau_objective(UnifiedDamageContext* ctx, double tau) noexcept nogil:
    """Compute negative log-posterior for τ (to minimize).

    Uses per-position weighted-least-squares-like objective that
    measures how well the decay curve fits the observed damage pattern.
    """
    cdef uint32_t r
    cdef int z
    cdef double g_new, g_old, log_ratio
    cdef double obj = 0.0
    cdef double old_tau = ctx.sample.tau
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params
    cdef double obs_rate, expected_rate
    cdef double weight, residual

    cdef double ref_opp

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        # Skip low-coverage refs (check directly)
        ref_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            ref_opp += counts.n_5p[z] + counts.n_3p[z]
        if ref_opp < 50:
            continue

        # 5' contribution: fit decay curve to observed rates
        for z in range(DAMAGE_FIT_POSITIONS):
            if counts.n_5p[z] < 10:
                continue

            g_new = decay_func(z + 1, tau)

            # Observed excess rate over baseline
            obs_rate = counts.k_5p[z] / counts.n_5p[z] - params.baseline_5p

            # Expected rate under current delta
            expected_rate = params.delta_5p * g_new

            # Weighted squared residual
            weight = counts.n_5p[z]
            if obs_rate > 0:
                residual = log(fmax(obs_rate, EPS)) - log(fmax(expected_rate, EPS))
                obj -= 0.5 * weight * residual * residual

        # 3' contribution
        for z in range(DAMAGE_FIT_POSITIONS):
            if counts.n_3p[z] < 10:
                continue

            g_new = decay_func(z + 1, tau)
            obs_rate = counts.k_3p[z] / counts.n_3p[z] - params.baseline_3p
            expected_rate = params.delta_3p * g_new

            weight = counts.n_3p[z]
            if obs_rate > 0:
                residual = log(fmax(obs_rate, EPS)) - log(fmax(expected_rate, EPS))
                obj -= 0.5 * weight * residual * residual

    # Log-normal prior on τ: -0.5 * ((log τ - log 5) / 0.5)^2
    cdef double log_tau = log(tau)
    cdef double prior = -0.5 * ((log_tau - ctx.sample.tau_prior_mean) / ctx.sample.tau_prior_sd) ** 2

    return -(obj + prior)  # return negative (we minimize)


cdef void em_m_step_tau(UnifiedDamageContext* ctx) noexcept nogil:
    """M-step: update τ using golden section search."""
    cdef double a = ctx.sample.tau_min
    cdef double b = ctx.sample.tau_max
    cdef double gr = 0.6180339887  # golden ratio
    cdef double c = b - gr * (b - a)
    cdef double d = a + gr * (b - a)
    cdef double fc, fd
    cdef double tol = 0.01
    cdef int max_iter = 30
    cdef int i

    for i in range(max_iter):
        fc = tau_objective(ctx, c)
        fd = tau_objective(ctx, d)

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

    ctx.sample.tau = (a + b) / 2.0


# =============================================================================
# Single EM iteration
# =============================================================================

cdef int run_em_iteration(UnifiedDamageContext* ctx) noexcept nogil:
    """Run one EM iteration. Returns 1 if converged, 0 otherwise."""
    cdef double old_tau = ctx.sample.tau
    cdef double max_delta_change = 0.0
    cdef uint32_t r
    cdef double old_delta_5p, old_delta_3p, change

    # Store old deltas
    cdef double* old_deltas_5p = <double*>malloc(ctx.n_refs * sizeof(double))
    cdef double* old_deltas_3p = <double*>malloc(ctx.n_refs * sizeof(double))

    if old_deltas_5p == NULL or old_deltas_3p == NULL:
        if old_deltas_5p != NULL:
            free(old_deltas_5p)
        if old_deltas_3p != NULL:
            free(old_deltas_3p)
        return 0

    for r in range(ctx.n_refs):
        old_deltas_5p[r] = ctx.ref_params[r].delta_5p
        old_deltas_3p[r] = ctx.ref_params[r].delta_3p

    # E-step
    em_e_step(ctx)

    # M-step: update deltas
    em_m_step_delta(ctx)

    # M-step: update tau (skip if prior is very tight, indicating fix_tau)
    # A tight prior (sd < 0.1) means we want to keep tau at the PMD-derived value
    if ctx.sample.tau_prior_sd >= 0.1:
        em_m_step_tau(ctx)

    ctx.sample.iteration += 1

    # Check convergence
    cdef double tau_change = fabs(log(ctx.sample.tau + EPS) - log(old_tau + EPS))

    for r in range(ctx.n_refs):
        change = fabs(log(ctx.ref_params[r].delta_5p + EPS) - log(old_deltas_5p[r] + EPS))
        if change > max_delta_change:
            max_delta_change = change

        change = fabs(log(ctx.ref_params[r].delta_3p + EPS) - log(old_deltas_3p[r] + EPS))
        if change > max_delta_change:
            max_delta_change = change

    free(old_deltas_5p)
    free(old_deltas_3p)

    if tau_change < ctx.tol_tau and max_delta_change < ctx.tol_delta:
        ctx.sample.converged = True
        return 1

    return 0


# =============================================================================
# SQUAREM integration for damage model
# =============================================================================

cdef inline uint32_t damage_param_count(UnifiedDamageContext* ctx) noexcept nogil:
    """Return total parameter count: tau + mu_5p + mu_3p + 2*n_refs deltas."""
    return 3 + 2 * ctx.n_refs


cdef void damage_pack_params(
    UnifiedDamageContext* ctx,
    double* theta
) noexcept nogil:
    """Pack damage model parameters into flat array.

    Layout: [tau, mu_5p, mu_3p, delta_5p[0..n_refs-1], delta_3p[0..n_refs-1]]
    """
    cdef uint32_t r
    theta[0] = ctx.sample.tau
    theta[1] = ctx.sample.mu_5p
    theta[2] = ctx.sample.mu_3p

    for r in range(ctx.n_refs):
        theta[3 + r] = ctx.ref_params[r].delta_5p
        theta[3 + ctx.n_refs + r] = ctx.ref_params[r].delta_3p


cdef void damage_unpack_params(
    double* theta,
    UnifiedDamageContext* ctx
) noexcept nogil:
    """Unpack flat array into damage model parameters."""
    cdef uint32_t r
    ctx.sample.tau = theta[0]
    ctx.sample.mu_5p = theta[1]
    ctx.sample.mu_3p = theta[2]

    for r in range(ctx.n_refs):
        ctx.ref_params[r].delta_5p = theta[3 + r]
        ctx.ref_params[r].delta_3p = theta[3 + ctx.n_refs + r]


cdef void damage_em_update(
    double* theta_in,
    double* theta_out,
    uint32_t n,
    void* ctx_ptr
) noexcept nogil:
    """EM fixed-point map for SQUAREM: theta_out = F(theta_in)."""
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr

    # Unpack input parameters
    damage_unpack_params(theta_in, ctx)

    # Run E-step
    em_e_step(ctx)

    # Run M-step for deltas
    em_m_step_delta(ctx)

    # Run M-step for tau (if not fixed)
    if ctx.sample.tau_prior_sd >= 0.1:
        em_m_step_tau(ctx)

    ctx.sample.iteration += 1

    # Pack output parameters
    damage_pack_params(ctx, theta_out)


cdef double damage_log_posterior(
    double* theta,
    uint32_t n,
    void* ctx_ptr
) noexcept nogil:
    """Compute log-posterior for SQUAREM acceptance check.

    Includes:
    - Poisson log-likelihood for damage counts
    - Log-normal prior on tau
    - Gamma shrinkage prior on deltas
    """
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr
    cdef uint32_t r
    cdef int z
    cdef double ll = 0.0
    cdef double tau, mu_5p, mu_3p, delta_5p, delta_3p
    cdef double g, p, k, nn
    cdef double baseline
    cdef RefDamageCounts* counts
    cdef double alpha_base = ctx.sample.alpha
    cdef double N_scale = 100.0
    cdef double alpha_r, N_r

    # Unpack parameters
    tau = theta[0]
    mu_5p = theta[1]
    mu_3p = theta[2]

    # Validate tau bounds
    if tau < ctx.sample.tau_min or tau > ctx.sample.tau_max:
        return -1e30

    # Validate mu bounds
    if mu_5p <= 0 or mu_5p > 1 or mu_3p <= 0 or mu_3p > 1:
        return -1e30

    # Log-normal prior on tau
    cdef double log_tau = log(tau)
    ll += -0.5 * ((log_tau - ctx.sample.tau_prior_mean) / ctx.sample.tau_prior_sd) ** 2

    # Per-reference contributions
    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        delta_5p = theta[3 + r]
        delta_3p = theta[3 + ctx.n_refs + r]

        # Validate delta bounds
        if delta_5p < 0 or delta_5p > 1 or delta_3p < 0 or delta_3p > 1:
            return -1e30

        N_r = counts.total_weight
        if N_r < 1:
            continue

        # Gamma prior on deltas (coverage-aware)
        alpha_r = alpha_base / sqrt(1.0 + N_r / N_scale)
        if alpha_r > 0 and mu_5p > 0:
            ll += (alpha_r - 1) * log(delta_5p + EPS) - alpha_r * delta_5p / mu_5p
        if alpha_r > 0 and mu_3p > 0:
            ll += (alpha_r - 1) * log(delta_3p + EPS) - alpha_r * delta_3p / mu_3p

        # Poisson log-likelihood for 5' end
        baseline = ctx.ref_params[r].baseline_5p
        for z in range(DAMAGE_FIT_POSITIONS):
            k = counts.k_5p[z]
            nn = counts.n_5p[z]
            if nn < 1:
                continue
            g = decay_func(z + 1, tau)
            p = baseline + delta_5p * g
            p = fmax(p, EPS)
            p = fmin(p, 1.0 - EPS)
            # Binomial approximation: k * log(p) + (n-k) * log(1-p)
            ll += k * log(p) + (nn - k) * log(1.0 - p)

        # Poisson log-likelihood for 3' end
        baseline = ctx.ref_params[r].baseline_3p
        for z in range(DAMAGE_FIT_POSITIONS):
            k = counts.k_3p[z]
            nn = counts.n_3p[z]
            if nn < 1:
                continue
            g = decay_func(z + 1, tau)
            p = baseline + delta_3p * g
            p = fmax(p, EPS)
            p = fmin(p, 1.0 - EPS)
            ll += k * log(p) + (nn - k) * log(1.0 - p)

    return ll


cdef void damage_project_params(
    double* theta,
    uint32_t n,
    void* ctx_ptr
) noexcept nogil:
    """Project damage parameters onto valid range.

    Constraints:
    - tau in [tau_min, tau_max]
    - mu in [EPS, 1]
    - delta in [0, 1]
    """
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr
    cdef uint32_t r

    # Clamp tau
    theta[0] = fmax(ctx.sample.tau_min, fmin(ctx.sample.tau_max, theta[0]))

    # Clamp mu
    theta[1] = fmax(EPS, fmin(1.0, theta[1]))
    theta[2] = fmax(EPS, fmin(1.0, theta[2]))

    # Clamp deltas
    for r in range(ctx.n_refs):
        theta[3 + r] = fmax(0.0, fmin(1.0, theta[3 + r]))
        theta[3 + ctx.n_refs + r] = fmax(0.0, fmin(1.0, theta[3 + ctx.n_refs + r]))


# =============================================================================
# Full model fitting
# =============================================================================

cdef int fit_damage_model(UnifiedDamageContext* ctx) noexcept nogil:
    """Fit the complete damage model using EM with SQUAREM acceleration."""
    cdef int i
    cdef int converged = 0
    cdef bint squarem_used
    cdef double ll_current = 0.0
    cdef double ll_prev = -1e30
    cdef double ll_change
    cdef double old_tau, max_delta_change, change
    cdef uint32_t r

    # SQUAREM state
    cdef uint32_t n_params = damage_param_count(ctx)
    cdef SquaremState sq
    cdef SquaremConfig sq_cfg
    cdef double* theta = NULL
    cdef bint squarem_enabled = True

    # Step 1: Estimate baselines from interior
    estimate_baselines(ctx)

    # Step 2: Initialize parameters
    initialize_parameters(ctx)

    bf_nogil_logf_notime(LOG_TAG, "Starting EM: tau=%.2f mu_5p=%.4f mu_3p=%.4f",
                         ctx.sample.tau, ctx.sample.mu_5p, ctx.sample.mu_3p)

    # Step 3: Initialize SQUAREM
    squarem_config_default(&sq_cfg)
    sq_cfg.enable = squarem_enabled
    sq_cfg.steplength_scheme = 3  # S3 recommended
    sq_cfg.max_backtracks = 4
    sq_cfg.enable_globalization = True

    if not squarem_state_init(&sq, n_params, NULL, 0):
        bf_nogil_logf_notime(LOG_TAG, "SQUAREM init failed, falling back to plain EM")
        squarem_enabled = False
    else:
        # Allocate parameter vector
        theta = <double*>malloc(n_params * sizeof(double))
        if theta == NULL:
            squarem_state_free(&sq)
            squarem_enabled = False

    # Step 4: Run EM iterations with SQUAREM
    if squarem_enabled:
        # Pack initial parameters
        damage_pack_params(ctx, theta)

        for i in range(ctx.max_iterations):
            old_tau = ctx.sample.tau

            # Store old deltas for convergence check
            squarem_used = squarem_step(
                &sq, theta, n_params,
                damage_em_update,
                damage_log_posterior,
                damage_project_params,
                <void*>ctx,  # update_ctx
                <void*>ctx,  # obj_ctx
                <void*>ctx,  # project_ctx
                &sq_cfg,
                &ll_current
            )

            # Unpack to check convergence
            damage_unpack_params(theta, ctx)

            # Check convergence based on log-likelihood change
            ll_change = fabs(ll_current - ll_prev)
            if ll_change < 1e-6 * fabs(ll_current) and i > 2:
                converged = 1
                break

            ll_prev = ll_current

        # Cleanup
        free(theta)
        squarem_state_free(&sq)

    else:
        # Fallback to plain EM
        for i in range(ctx.max_iterations):
            converged = run_em_iteration(ctx)
            if converged:
                break

    bf_nogil_logf_notime(LOG_TAG, "EM finished: iter=%d converged=%d tau=%.2f",
                         ctx.sample.iteration, converged, ctx.sample.tau)

    # Step 5: Compute output quantities
    compute_outputs(ctx)

    return converged


# =============================================================================
# Output computation
# =============================================================================

cdef void compute_outputs(UnifiedDamageContext* ctx) noexcept nogil:
    """Compute all output quantities from fitted model.

    Computes a unified authentication score (p_ancient) that combines:
    1. Damage authenticity: Y / (Y + B + 1)
    2. Coverage confidence: 1 - exp(-N / N_confident)
    3. Asymmetry consistency: genuine damage shows CT5/GA3 > GA5/CT3
    4. ANI quality: corrected ANI should be reasonable

    This multi-factor approach is more robust than damage-only scoring.
    """
    cdef uint32_t r
    cdef int z
    cdef double Y_total, B_total
    cdef double auth, p_anc
    cdef double log_bf
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params
    cdef double tau = ctx.sample.tau
    cdef double logit_p
    cdef double damage_signal, control_signal
    cdef double g, rate_damage, rate_baseline
    cdef double total_opp

    # Unified score parameters
    cdef double coverage_confidence, asymmetry_bonus, ani_penalty
    cdef double N_confident = 50.0  # coverage for ~63% confidence
    cdef double N_r

    # Final E-step to get clean Y/B values
    em_e_step(ctx)

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        Y_total = params.Y_5p + params.Y_3p
        B_total = params.B_5p + params.B_3p

        # Authenticity score: Y / (Y + B + 1)
        params.authenticity = Y_total / (Y_total + B_total + 1.0)

        # Coverage confidence: increases with sample size
        N_r = counts.total_weight
        coverage_confidence = 1.0 - exp(-N_r / N_confident)

        # Compute asymmetry early for bonus calculation
        damage_signal = 0.0
        control_signal = 0.0
        for z in range(5):  # positions 1-5
            damage_signal += counts.k_5p[z] + counts.k_3p[z]
            control_signal += counts.k_5p_ctrl[z] + counts.k_3p_ctrl[z]
        params.asymmetry = damage_signal - control_signal

        # Asymmetry bonus: genuine damage should show CT5+GA3 > GA5+CT3
        # sigmoid((asymmetry / max(1, total_signal)) - 0) gives 0.5 baseline, up to ~1.5x bonus
        if damage_signal + control_signal > EPS:
            asymmetry_bonus = 0.5 + 0.5 / (1.0 + exp(-params.asymmetry / (damage_signal + control_signal + 1.0)))
        else:
            asymmetry_bonus = 0.5

        # P(ancient) using enhanced logistic transform with multi-factor weighting
        # Base: authenticity and evidence strength
        logit_p = 8.0 * (params.authenticity - 0.5) + log(Y_total + 1.0)

        # Apply coverage confidence: low coverage → shrink toward 0.5
        logit_p = logit_p * coverage_confidence

        # Apply asymmetry bonus: consistent damage pattern → boost score
        logit_p = logit_p + log(fmax(asymmetry_bonus, 0.1))

        if logit_p > 20.0:
            params.p_ancient = 1.0 - EPS
        elif logit_p < -20.0:
            params.p_ancient = EPS
        else:
            params.p_ancient = 1.0 / (1.0 + exp(-logit_p))

        # Log Bayes factor (approximate)
        # log_BF ≈ Y * log(delta * g / baseline) summed over positions
        log_bf = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            g = decay_func(z + 1, tau)

            # 5' contribution
            if counts.n_5p[z] > EPS and params.baseline_5p > EPS:
                rate_damage = params.baseline_5p + params.delta_5p * g
                rate_baseline = params.baseline_5p
                if rate_damage > EPS:
                    log_bf += counts.k_5p[z] * log(rate_damage / rate_baseline)

            # 3' contribution
            if counts.n_3p[z] > EPS and params.baseline_3p > EPS:
                rate_damage = params.baseline_3p + params.delta_3p * g
                rate_baseline = params.baseline_3p
                if rate_damage > EPS:
                    log_bf += counts.k_3p[z] * log(rate_damage / rate_baseline)

        params.log_bf = log_bf

        # Quality flags
        params.has_damage_evidence = Y_total >= 3.0

        total_opp = 0.0
        for z in range(DAMAGE_FIT_POSITIONS):
            total_opp += counts.n_5p[z] + counts.n_3p[z]
        params.is_low_coverage = total_opp < 10.0

        # ANI computation
        if counts.total_aligned > 0:
            params.ani_raw = <double>counts.total_matches / <double>counts.total_aligned
        else:
            params.ani_raw = 0.0

        # Damage-corrected ANI
        params.damage_correction = compute_ani_correction(counts, params, tau)
        params.ani_corrected = params.ani_raw + params.damage_correction

    # Compute posterior predictive diagnostics
    compute_posterior_predictive(ctx)


cdef double compute_ani_correction(
    RefDamageCounts* counts,
    RefDamageParams* params,
    double tau
) noexcept nogil:
    """Compute ANI correction from expected damage mismatches."""
    cdef int z
    cdef double Y_corr = 0.0
    cdef double g, w, denom

    # Sum expected damage mismatches that should be "corrected"
    # (treated as matches rather than mismatches)

    # 5' end
    for z in range(DAMAGE_FIT_POSITIONS):
        g = decay_func(z + 1, tau)
        denom = params.baseline_5p + params.delta_5p * g + EPS
        w = (params.delta_5p * g) / denom
        Y_corr += counts.k_5p[z] * w

    # 3' end
    for z in range(DAMAGE_FIT_POSITIONS):
        g = decay_func(z + 1, tau)
        denom = params.baseline_3p + params.delta_3p * g + EPS
        w = (params.delta_3p * g) / denom
        Y_corr += counts.k_3p[z] * w

    # Correction = Y_corr / total_aligned
    if counts.total_aligned > 0:
        return Y_corr / <double>counts.total_aligned
    return 0.0


# =============================================================================
# Posterior predictive diagnostics
# =============================================================================

cdef inline double chi2_survival(double x, int df) noexcept nogil:
    """Approximate chi-squared survival function P(X > x) for df degrees of freedom.

    Uses Wilson-Hilferty approximation for df >= 3:
        Z = ((x/df)^(1/3) - (1 - 2/(9*df))) / sqrt(2/(9*df))
    Then use standard normal survival approximation.
    """
    cdef double z, t, p
    cdef double cube_root

    if df <= 0 or x <= 0:
        return 1.0

    if df == 1:
        # Chi-squared with df=1 is sqrt(chi2) ~ |N(0,1)|
        # P(X > x) = 2 * P(Z > sqrt(x)) for Z ~ N(0,1)
        z = sqrt(x)
        t = 1.0 / (1.0 + 0.2316419 * z)
        p = 0.3989422804 * exp(-0.5 * z * z) * t * (
            0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
        )
        return 2.0 * p
    elif df == 2:
        # Chi-squared with df=2 is exponential
        return exp(-0.5 * x)

    # Wilson-Hilferty approximation for df >= 3
    cube_root = exp(log(x / <double>df) / 3.0) if x > 0 else 0.0
    z = (cube_root - (1.0 - 2.0 / (9.0 * <double>df))) / sqrt(2.0 / (9.0 * <double>df))

    # Standard normal survival function approximation
    if z > 6.0:
        return 0.0
    elif z < -6.0:
        return 1.0

    # Horner's method for normal CDF approximation
    t = 1.0 / (1.0 + 0.2316419 * fabs(z))
    p = 0.3989422804 * exp(-0.5 * z * z) * t * (
        0.319381530 + t * (-0.356563782 + t * (1.781477937 + t * (-1.821255978 + t * 1.330274429)))
    )

    if z > 0:
        return p
    else:
        return 1.0 - p


cdef void compute_posterior_predictive(UnifiedDamageContext* ctx) noexcept nogil:
    """Compute posterior predictive diagnostics for model validation.

    For each reference, compute:
    1. Chi-squared statistic: Σ (k - E[k])² / E[k]
    2. Posterior predictive p-value from chi-squared distribution
    3. Dispersion factor: observed variance / expected variance
    4. Overdispersion flag when model doesn't fit

    These diagnostics help identify references where:
    - The exponential decay model is misspecified
    - Damage patterns are inconsistent with ancient DNA
    - Contamination mixes ancient and modern DNA
    """
    cdef uint32_t r
    cdef int z
    cdef double E_k, k_obs
    cdef double chi2, residual, sum_residual_sq
    cdef double g, rate
    cdef double tau = ctx.sample.tau
    cdef RefDamageCounts* counts
    cdef RefDamageParams* params
    cdef int df_count
    cdef double n_opp

    for r in range(ctx.n_refs):
        counts = &ctx.ref_counts[r]
        params = &ctx.ref_params[r]

        chi2 = 0.0
        sum_residual_sq = 0.0
        df_count = 0

        # 5' end positions
        for z in range(DAMAGE_FIT_POSITIONS):
            n_opp = counts.n_5p[z]
            k_obs = counts.k_5p[z]

            if n_opp < 1.0:
                continue

            # Expected count: E[k] = n * (baseline + delta * g(z))
            g = decay_func(z + 1, tau)
            rate = params.baseline_5p + params.delta_5p * g
            E_k = n_opp * rate

            if E_k > EPS:
                residual = k_obs - E_k
                chi2 += (residual * residual) / E_k
                sum_residual_sq += residual * residual
                df_count += 1

        # 3' end positions
        for z in range(DAMAGE_FIT_POSITIONS):
            n_opp = counts.n_3p[z]
            k_obs = counts.k_3p[z]

            if n_opp < 1.0:
                continue

            # Expected count
            g = decay_func(z + 1, tau)
            rate = params.baseline_3p + params.delta_3p * g
            E_k = n_opp * rate

            if E_k > EPS:
                residual = k_obs - E_k
                chi2 += (residual * residual) / E_k
                sum_residual_sq += residual * residual
                df_count += 1

        # Store results
        params.chi_squared = chi2

        # Degrees of freedom: positions counted minus parameters estimated (delta_5p, delta_3p, baselines)
        # Conservative: df = positions - 4 (2 deltas + 2 baselines)
        params.df = df_count - 4 if df_count > 4 else 1

        # Posterior predictive p-value: P(chi2 > observed | model)
        if params.df > 0:
            params.pp_pvalue = chi2_survival(chi2, params.df)
        else:
            params.pp_pvalue = 1.0

        # Dispersion factor: chi2/df should be ~1 if model fits
        if params.df > 0:
            params.dispersion_factor = chi2 / <double>params.df
        else:
            params.dispersion_factor = 1.0

        # Flag overdispersion: p < 0.01 or dispersion > 3
        params.is_overdispersed = (params.pp_pvalue < 0.01) or (params.dispersion_factor > 3.0)


# =============================================================================
# Python-accessible interface
# =============================================================================

def create_context_py(uint32_t n_refs, bint is_single_stranded=False, bint mask_cpg=True):
    """Create unified damage context.

    Parameters
    ----------
    n_refs : int
        Number of references
    is_single_stranded : bool
        True for single-stranded library
    mask_cpg : bool
        Mask CpG sites from counts

    Returns
    -------
    int
        Pointer to context (as integer)
    """
    cdef UnifiedDamageContext* ctx = create_damage_context(n_refs, is_single_stranded, mask_cpg)
    if ctx == NULL:
        raise MemoryError("Failed to allocate UnifiedDamageContext")
    return <uintptr_t>ctx


def destroy_context_py(uintptr_t ctx_ptr):
    """Free context memory."""
    destroy_damage_context(<UnifiedDamageContext*>ctx_ptr)


def fit_model_py(uintptr_t ctx_ptr):
    """Fit the damage model using EM.

    Returns
    -------
    dict
        Fitting results including tau, convergence status
    """
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr
    cdef int converged

    with nogil:
        converged = fit_damage_model(ctx)

    return {
        'converged': bool(converged),
        'iterations': ctx.sample.iteration,
        'tau': ctx.sample.tau,
        'mu_5p': ctx.sample.mu_5p,
        'mu_3p': ctx.sample.mu_3p,
        'baseline_5p_global': ctx.sample.baseline_5p_global,
        'baseline_3p_global': ctx.sample.baseline_3p_global,
    }


def get_reference_results_py(uintptr_t ctx_ptr, uint32_t ref_idx):
    """Get fitted parameters and outputs for a reference.

    Returns
    -------
    dict
        Per-reference damage model results
    """
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr

    if ref_idx >= ctx.n_refs:
        raise IndexError(f"Reference index {ref_idx} out of range (n_refs={ctx.n_refs})")

    cdef RefDamageParams* params = &ctx.ref_params[ref_idx]

    return {
        'delta_5p': params.delta_5p,
        'delta_3p': params.delta_3p,
        'baseline_5p': params.baseline_5p,
        'baseline_3p': params.baseline_3p,
        'Y_5p': params.Y_5p,
        'Y_3p': params.Y_3p,
        'B_5p': params.B_5p,
        'B_3p': params.B_3p,
        'authenticity': params.authenticity,
        'p_ancient': params.p_ancient,
        'log_bf': params.log_bf,
        'asymmetry': params.asymmetry,
        'ani_raw': params.ani_raw,
        'ani_corrected': params.ani_corrected,
        'damage_correction': params.damage_correction,
        'has_damage_evidence': bool(params.has_damage_evidence),
        'is_low_coverage': bool(params.is_low_coverage),
        'is_overdispersed': bool(params.is_overdispersed),
        'chi_squared': params.chi_squared,
        'pp_pvalue': params.pp_pvalue,
        'dispersion_factor': params.dispersion_factor,
        'df': params.df,
    }


def get_all_results_py(uintptr_t ctx_ptr):
    """Get results for all references as numpy arrays.

    Returns
    -------
    dict
        Arrays of per-reference results including diagnostics
    """
    import numpy as np

    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr
    cdef uint32_t n = ctx.n_refs
    cdef uint32_t r

    # Allocate arrays
    delta_5p = np.zeros(n, dtype=np.float64)
    delta_3p = np.zeros(n, dtype=np.float64)
    authenticity = np.zeros(n, dtype=np.float64)
    p_ancient = np.zeros(n, dtype=np.float64)
    log_bf = np.zeros(n, dtype=np.float64)
    ani_raw = np.zeros(n, dtype=np.float64)
    ani_corrected = np.zeros(n, dtype=np.float64)
    chi_squared = np.zeros(n, dtype=np.float64)
    pp_pvalue = np.zeros(n, dtype=np.float64)
    dispersion_factor = np.zeros(n, dtype=np.float64)
    is_overdispersed = np.zeros(n, dtype=np.uint8)

    cdef double[:] delta_5p_view = delta_5p
    cdef double[:] delta_3p_view = delta_3p
    cdef double[:] auth_view = authenticity
    cdef double[:] p_anc_view = p_ancient
    cdef double[:] log_bf_view = log_bf
    cdef double[:] ani_raw_view = ani_raw
    cdef double[:] ani_corr_view = ani_corrected
    cdef double[:] chi2_view = chi_squared
    cdef double[:] pp_view = pp_pvalue
    cdef double[:] disp_view = dispersion_factor
    cdef uint8_t[:] overdisp_view = is_overdispersed

    for r in range(n):
        delta_5p_view[r] = ctx.ref_params[r].delta_5p
        delta_3p_view[r] = ctx.ref_params[r].delta_3p
        auth_view[r] = ctx.ref_params[r].authenticity
        p_anc_view[r] = ctx.ref_params[r].p_ancient
        log_bf_view[r] = ctx.ref_params[r].log_bf
        ani_raw_view[r] = ctx.ref_params[r].ani_raw
        ani_corr_view[r] = ctx.ref_params[r].ani_corrected
        chi2_view[r] = ctx.ref_params[r].chi_squared
        pp_view[r] = ctx.ref_params[r].pp_pvalue
        disp_view[r] = ctx.ref_params[r].dispersion_factor
        overdisp_view[r] = ctx.ref_params[r].is_overdispersed

    return {
        'tau': ctx.sample.tau,
        'delta_5p': delta_5p,
        'delta_3p': delta_3p,
        'authenticity': authenticity,
        'p_ancient': p_ancient,
        'log_bf': log_bf,
        'ani_raw': ani_raw,
        'ani_corrected': ani_corrected,
        'chi_squared': chi_squared,
        'pp_pvalue': pp_pvalue,
        'dispersion_factor': dispersion_factor,
        'is_overdispersed': is_overdispersed,
    }


def set_counts_from_arrays_py(
    uintptr_t ctx_ptr,
    double[:, ::1] k_5p,
    double[:, ::1] n_5p,
    double[:, ::1] k_3p,
    double[:, ::1] n_3p,
):
    """Set damage counts from numpy arrays.

    Parameters
    ----------
    ctx_ptr : int
        Context pointer
    k_5p, n_5p : ndarray, shape (n_refs, 20)
        5' mismatch and opportunity counts
    k_3p, n_3p : ndarray, shape (n_refs, 20)
        3' mismatch and opportunity counts
    """
    cdef UnifiedDamageContext* ctx = <UnifiedDamageContext*>ctx_ptr
    cdef uint32_t r
    cdef int z

    if k_5p.shape[0] != ctx.n_refs or k_5p.shape[1] != 20:
        raise ValueError(f"Expected shape ({ctx.n_refs}, 20), got {k_5p.shape}")

    for r in range(ctx.n_refs):
        for z in range(20):
            ctx.ref_counts[r].k_5p[z] = k_5p[r, z]
            ctx.ref_counts[r].n_5p[z] = n_5p[r, z]
            ctx.ref_counts[r].k_3p[z] = k_3p[r, z]
            ctx.ref_counts[r].n_3p[z] = n_3p[r, z]


def fit_unified_damage_model_py(
    double[:, ::1] k_5p,
    double[:, ::1] n_5p,
    double[:, ::1] k_3p,
    double[:, ::1] n_3p,
    bint is_single_stranded=False,
    bint mask_cpg=False,
):
    """High-level function to fit the unified damage model.

    This is the main entry point for integration with the processing pipeline.

    Parameters
    ----------
    k_5p, n_5p : ndarray, shape (n_refs, 20)
        5' mismatch and opportunity counts per reference
    k_3p, n_3p : ndarray, shape (n_refs, 20)
        3' mismatch and opportunity counts per reference
    is_single_stranded : bool
        True for single-stranded library
    mask_cpg : bool
        Mask CpG sites from analysis

    Returns
    -------
    dict
        Results including per-reference metrics (authenticity, p_ancient, etc.)
        and sample-level parameters (tau, mu)
    """
    import numpy as np

    cdef uint32_t n_refs = k_5p.shape[0]
    cdef uint32_t r
    cdef int z
    cdef int converged
    cdef UnifiedDamageContext* ctx

    # Create context
    ctx = create_damage_context(n_refs, is_single_stranded, mask_cpg)
    if ctx == NULL:
        raise MemoryError("Failed to allocate UnifiedDamageContext")

    try:
        # Populate counts
        for r in range(n_refs):
            for z in range(20):
                ctx.ref_counts[r].k_5p[z] = k_5p[r, z]
                ctx.ref_counts[r].n_5p[z] = n_5p[r, z]
                ctx.ref_counts[r].k_3p[z] = k_3p[r, z]
                ctx.ref_counts[r].n_3p[z] = n_3p[r, z]

        # Fit model
        with nogil:
            converged = fit_damage_model(ctx)

        # Extract results
        delta_5p = np.zeros(n_refs, dtype=np.float64)
        delta_3p = np.zeros(n_refs, dtype=np.float64)
        baseline_5p = np.zeros(n_refs, dtype=np.float64)
        baseline_3p = np.zeros(n_refs, dtype=np.float64)
        authenticity = np.zeros(n_refs, dtype=np.float64)
        p_ancient = np.zeros(n_refs, dtype=np.float64)
        log_bf = np.zeros(n_refs, dtype=np.float64)
        ani_corrected = np.zeros(n_refs, dtype=np.float64)
        damage_correction = np.zeros(n_refs, dtype=np.float64)

        for r in range(n_refs):
            delta_5p[r] = ctx.ref_params[r].delta_5p
            delta_3p[r] = ctx.ref_params[r].delta_3p
            baseline_5p[r] = ctx.ref_params[r].baseline_5p
            baseline_3p[r] = ctx.ref_params[r].baseline_3p
            authenticity[r] = ctx.ref_params[r].authenticity
            p_ancient[r] = ctx.ref_params[r].p_ancient
            log_bf[r] = ctx.ref_params[r].log_bf
            ani_corrected[r] = ctx.ref_params[r].ani_corrected
            damage_correction[r] = ctx.ref_params[r].damage_correction

        # Summary statistics
        n_ancient = int(np.sum(p_ancient > 0.5))
        n_confident_ancient = int(np.sum(p_ancient > 0.9))
        n_confident_modern = int(np.sum(p_ancient < 0.1))

        return {
            'converged': bool(converged),
            'iterations': ctx.sample.iteration,
            'tau': ctx.sample.tau,
            'mu_5p': ctx.sample.mu_5p,
            'mu_3p': ctx.sample.mu_3p,
            'baseline_5p_global': ctx.sample.baseline_5p_global,
            'baseline_3p_global': ctx.sample.baseline_3p_global,
            # Per-reference arrays
            'delta_5p': delta_5p,
            'delta_3p': delta_3p,
            'baseline_5p': baseline_5p,
            'baseline_3p': baseline_3p,
            'authenticity': authenticity,
            'p_ancient': p_ancient,
            'log_bf': log_bf,
            'ani_corrected': ani_corrected,
            'damage_correction': damage_correction,
            # Summary
            'n_refs': int(n_refs),
            'n_ancient': n_ancient,
            'n_confident_ancient': n_confident_ancient,
            'n_confident_modern': n_confident_modern,
        }

    finally:
        destroy_damage_context(ctx)


def fit_unified_from_damage_stats_py(
    uintptr_t pool_ptr,
    uintptr_t stats_ptr,
    bint is_single_stranded,
):
    """Fit unified damage model from RefDamageStats and store results in MemoryPool.

    This is the integration point with the existing pipeline.

    Parameters
    ----------
    pool_ptr : uintptr_t
        Pointer to MemoryPool
    stats_ptr : uintptr_t
        Pointer to RefDamageStats array
    is_single_stranded : bool
        True for single-stranded library

    Returns
    -------
    dict
        Summary statistics
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef RefDamageStats* stats = <RefDamageStats*>stats_ptr
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t r
    cdef int z
    cdef int converged
    cdef UnifiedDamageContext* ctx

    cdef double sum_p_ancient = 0.0
    cdef double max_p_ancient = 0.0
    cdef double min_p_ancient = 1.0
    cdef uint32_t n_ancient = 0
    cdef uint32_t n_fitted = 0
    cdef PMDCurve* pmd_curve = NULL
    cdef double initial_tau = 5.0
    cdef uint32_t n_overdispersed = 0
    cdef double sum_dispersion = 0.0
    cdef double sum_pp_pvalue = 0.0

    ctx = create_damage_context(n_refs, is_single_stranded, 0)
    if ctx == NULL:
        raise MemoryError("Failed to allocate UnifiedDamageContext")

    # Get tau from PMD curve if available
    if pool.pmd_curve_ptr != NULL:
        pmd_curve = <PMDCurve*>pool.pmd_curve_ptr
        if pmd_curve.lambda_decay > 0.01:
            initial_tau = 1.0 / pmd_curve.lambda_decay
            set_tau_from_pmd(ctx, initial_tau, True)  # fix_tau=True

    try:
        # Copy counts from RefDamageStats to UnifiedDamageContext
        for r in range(n_refs):
            for z in range(20):
                ctx.ref_counts[r].k_5p[z] = stats[r].k_5p[z]
                ctx.ref_counts[r].n_5p[z] = stats[r].n_5p[z]
                ctx.ref_counts[r].k_3p[z] = stats[r].k_3p[z]
                ctx.ref_counts[r].n_3p[z] = stats[r].n_3p[z]
            ctx.ref_counts[r].total_weight = stats[r].total_weight

        # Fit model
        with nogil:
            converged = fit_damage_model(ctx)

        # Copy results to MemoryPool arrays
        if pool.damage_amplitude != NULL:
            for r in range(n_refs):
                pool.damage_amplitude[r] = ctx.ref_params[r].delta_5p
        if pool.damage_baseline != NULL:
            for r in range(n_refs):
                pool.damage_baseline[r] = ctx.ref_params[r].baseline_5p
        if pool.damage_log_bf != NULL:
            for r in range(n_refs):
                pool.damage_log_bf[r] = ctx.ref_params[r].log_bf

        # Summary statistics
        for r in range(n_refs):
            if stats[r].total_weight > 0.1:
                n_fitted += 1
                sum_p_ancient += ctx.ref_params[r].p_ancient
                sum_dispersion += ctx.ref_params[r].dispersion_factor
                sum_pp_pvalue += ctx.ref_params[r].pp_pvalue
                if ctx.ref_params[r].p_ancient > max_p_ancient:
                    max_p_ancient = ctx.ref_params[r].p_ancient
                if ctx.ref_params[r].p_ancient < min_p_ancient:
                    min_p_ancient = ctx.ref_params[r].p_ancient
                if ctx.ref_params[r].p_ancient > 0.5:
                    n_ancient += 1
                if ctx.ref_params[r].is_overdispersed:
                    n_overdispersed += 1

        bf_nogil_logf_notime(
            LOG_TAG,
            b"Unified model: tau=%.2f, mu_5p=%.4f, %u/%u ancient (p>0.5)",
            ctx.sample.tau, ctx.sample.mu_5p, n_ancient, n_fitted
        )
        bf_nogil_logf_notime(
            LOG_TAG,
            b"Diagnostics: %u/%u overdispersed (disp>3 or p<0.01), mean_disp=%.2f",
            n_overdispersed, n_fitted, sum_dispersion / n_fitted if n_fitted > 0 else 1.0
        )

        return {
            'converged': bool(converged),
            'iterations': ctx.sample.iteration,
            'tau': ctx.sample.tau,
            'mu_5p': ctx.sample.mu_5p,
            'mu_3p': ctx.sample.mu_3p,
            'baseline_5p_global': ctx.sample.baseline_5p_global,
            'n_refs': int(n_refs),
            'n_fitted': int(n_fitted),
            'n_ancient': int(n_ancient),
            'n_overdispersed': int(n_overdispersed),
            'mean_p_ancient': sum_p_ancient / n_fitted if n_fitted > 0 else 0.5,
            'min_p_ancient': float(min_p_ancient) if n_fitted > 0 else 0.5,
            'max_p_ancient': float(max_p_ancient) if n_fitted > 0 else 0.5,
            'mean_dispersion': sum_dispersion / n_fitted if n_fitted > 0 else 1.0,
            'mean_pp_pvalue': sum_pp_pvalue / n_fitted if n_fitted > 0 else 1.0,
        }

    finally:
        destroy_damage_context(ctx)
