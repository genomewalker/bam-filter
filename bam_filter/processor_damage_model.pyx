# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# -*- coding: utf-8 -*-
"""Bayesian damage model for ancient/modern reference classification.

Implements a hierarchical Beta-binomial model that:
1. Separates damage from divergence using interior baseline estimation
2. Estimates per-reference damage amplitude with shrinkage to global
3. Computes Bayes factor P(data|ancient) / P(data|modern)
4. Produces proper posterior P(ancient | data)

More principled than metaDMG's Z-score approach.
"""

from libc.math cimport exp, log, log1p, fabs, fmax, fmin, sqrt, pow as cpow, lgamma
from libc.stdlib cimport malloc, calloc, free
from libc.string cimport memset, memcpy
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t, uintptr_t

from cython.parallel cimport prange, parallel

from bam_filter.processor cimport MemoryPool, Alignment, AlignmentCore, DamageCounts
from bam_filter.processor_pmd cimport (
    RefDamageStats, DamageModelHyperparams, PMDCurve, ANISnapshot
)

# Logging
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef const char* LOG_TAG = b"DAMAGE_MODEL"

# Constants
cdef int INTERIOR_START = 15  # Positions 16-20 are "interior" (0-indexed: 15-19)
cdef int NUM_POSITIONS = 20
cdef double MIN_PROB = 1e-10
cdef double MAX_PROB = 1.0 - 1e-10


# =============================================================================
# Beta-Binomial Likelihood Functions
# =============================================================================

cdef inline double log_beta(double a, double b) noexcept nogil:
    """Log of Beta function: log(B(a,b)) = lgamma(a) + lgamma(b) - lgamma(a+b)"""
    return lgamma(a) + lgamma(b) - lgamma(a + b)


cdef inline double log_beta_binomial_pmf(double k, double n, double alpha, double beta) noexcept nogil:
    """Log probability mass function of Beta-Binomial distribution.

    P(k | n, α, β) = C(n,k) × B(k+α, n-k+β) / B(α, β)

    Parameters
    ----------
    k : double
        Number of successes (can be fractional for weighted counts)
    n : double
        Number of trials
    alpha, beta : double
        Beta distribution parameters

    Returns
    -------
    double
        Log probability
    """
    if n < 0.01:  # No data
        return 0.0

    # Handle edge cases
    if k < 0:
        k = 0
    if k > n:
        k = n

    # log C(n,k) = lgamma(n+1) - lgamma(k+1) - lgamma(n-k+1)
    # For fractional k,n we use the continuous extension
    cdef double log_binom = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)

    # log B(k+α, n-k+β) - log B(α, β)
    cdef double log_beta_ratio = log_beta(k + alpha, n - k + beta) - log_beta(alpha, beta)

    return log_binom + log_beta_ratio


cdef inline double beta_binomial_log_lik(
    double* k_arr, double* n_arr, int n_positions,
    double* mean_arr, double concentration
) noexcept nogil:
    """Compute total log-likelihood under Beta-Binomial model.

    For each position z:
        k_z ~ BetaBinom(n_z, α_z, β_z)
        where α_z = mean_z × c, β_z = (1 - mean_z) × c

    Parameters
    ----------
    k_arr : double*
        Observed mismatches at each position
    n_arr : double*
        Opportunities at each position
    n_positions : int
        Number of positions to include
    mean_arr : double*
        Expected mismatch rate at each position
    concentration : double
        Concentration parameter (higher = less overdispersion)

    Returns
    -------
    double
        Total log-likelihood
    """
    cdef double log_lik = 0.0
    cdef double mean_z, alpha_z, beta_z
    cdef int z

    for z in range(n_positions):
        if n_arr[z] < 0.01:  # Skip positions with no data
            continue

        mean_z = fmax(MIN_PROB, fmin(MAX_PROB, mean_arr[z]))
        alpha_z = mean_z * concentration
        beta_z = (1.0 - mean_z) * concentration

        # Ensure valid parameters
        if alpha_z < 0.01:
            alpha_z = 0.01
        if beta_z < 0.01:
            beta_z = 0.01

        log_lik += log_beta_binomial_pmf(k_arr[z], n_arr[z], alpha_z, beta_z)

    return log_lik


# =============================================================================
# Baseline Estimation from Interior Positions
# =============================================================================

cdef void estimate_baseline(
    RefDamageStats* stats,
    double prior_alpha,
    double prior_beta
) noexcept nogil:
    """Estimate baseline mismatch rate from interior positions.

    Uses Beta-Binomial conjugate update:
        baseline ~ Beta(α₀ + k_int, β₀ + n_int - k_int)

    Parameters
    ----------
    stats : RefDamageStats*
        Per-reference statistics (will be updated with baseline)
    prior_alpha, prior_beta : double
        Prior Beta parameters from global interior rate
    """
    # Accumulate interior counts (positions 16-20, i.e., indices 15-19)
    cdef double n_int = 0.0
    cdef double k_int = 0.0
    cdef int z

    for z in range(INTERIOR_START, NUM_POSITIONS):
        n_int += stats.n_5p[z] + stats.n_3p[z]
        k_int += stats.k_5p[z] + stats.k_3p[z]

    stats.n_interior = n_int
    stats.k_interior = k_int

    # Posterior parameters
    stats.baseline_alpha = prior_alpha + k_int
    stats.baseline_beta = prior_beta + (n_int - k_int)

    # Posterior mean as point estimate
    cdef double denom = stats.baseline_alpha + stats.baseline_beta
    if denom > 0:
        stats.baseline = stats.baseline_alpha / denom
    else:
        stats.baseline = prior_alpha / (prior_alpha + prior_beta)  # Fall back to prior


# =============================================================================
# Damage Amplitude Estimation and Bayes Factor
# =============================================================================

cdef void compute_mean_array_modern(
    double* mean_arr,
    double baseline,
    int n_positions
) noexcept nogil:
    """Fill mean array for modern model: p_z = baseline (constant)."""
    cdef int z
    for z in range(n_positions):
        mean_arr[z] = baseline


cdef void compute_mean_array_ancient(
    double* mean_arr,
    double baseline,
    double amplitude,
    double* D_shape,
    int n_positions
) noexcept nogil:
    """Fill mean array for ancient model: p_z = baseline + A × D_shape(z)."""
    cdef int z
    cdef double p_z
    for z in range(n_positions):
        p_z = baseline + amplitude * D_shape[z]
        mean_arr[z] = fmax(MIN_PROB, fmin(MAX_PROB, p_z))


cdef double estimate_amplitude_mle(
    RefDamageStats* stats,
    double* D_shape,
    int n_end_positions
) noexcept nogil:
    """Estimate damage amplitude using method of moments.

    A_hat = (observed_rate_at_ends - baseline) / mean(D_shape_at_ends)

    This is a fast approximation; could use MLE with optimization for more accuracy.
    """
    cdef double total_k = 0.0
    cdef double total_n = 0.0
    cdef double total_D = 0.0
    cdef int z

    # Use first 10 positions (damage zone)
    for z in range(n_end_positions):
        total_k += stats.k_5p[z] + stats.k_3p[z]
        total_n += stats.n_5p[z] + stats.n_3p[z]
        total_D += D_shape[z]

    if total_n < 1.0 or total_D < 0.01:
        return 0.0

    cdef double observed_rate = total_k / total_n
    cdef double excess = observed_rate - stats.baseline
    cdef double mean_D = total_D / n_end_positions

    if excess <= 0 or mean_D <= 0:
        return 0.0

    return excess / mean_D


cdef void compute_bayes_factor(
    RefDamageStats* stats,
    DamageModelHyperparams* hyper,
    int use_both_ends
) noexcept nogil:
    """Compute Bayes factor for ancient vs modern model.

    BF = P(data | ancient, A_hat) / P(data | modern)

    Uses point estimate for amplitude with prior regularization.
    """
    cdef double mean_modern[20]
    cdef double mean_ancient[20]
    cdef double k_combined[20]
    cdef double n_combined[20]
    cdef int z
    cdef int n_positions = 15  # Use positions 1-15 (indices 0-14) for damage signal

    # Combine 5' and 3' counts
    for z in range(NUM_POSITIONS):
        if use_both_ends:
            k_combined[z] = stats.k_5p[z] + stats.k_3p[z]
            n_combined[z] = stats.n_5p[z] + stats.n_3p[z]
        else:
            # Use 5' only (for single-stranded or asymmetric damage)
            k_combined[z] = stats.k_5p[z]
            n_combined[z] = stats.n_5p[z]

    # Estimate amplitude with shrinkage
    cdef double A_mle = estimate_amplitude_mle(stats, hyper.D_shape, 10)

    # Apply shrinkage toward prior mean
    # A_shrunk = (A_mle × n_eff + μ_A × prior_strength) / (n_eff + prior_strength)
    cdef double prior_strength = 1.0 / (hyper.amplitude_sigma * hyper.amplitude_sigma)
    cdef double A_prior_mean = exp(hyper.amplitude_mu)
    cdef double n_eff = stats.total_weight

    if n_eff + prior_strength > 0:
        stats.amplitude = (A_mle * n_eff + A_prior_mean * prior_strength) / (n_eff + prior_strength)
    else:
        stats.amplitude = A_prior_mean

    # Clamp amplitude to valid range
    stats.amplitude = fmax(0.0, fmin(1.0 - stats.baseline, stats.amplitude))

    # Compute mean arrays
    compute_mean_array_modern(mean_modern, stats.baseline, n_positions)
    compute_mean_array_ancient(mean_ancient, stats.baseline, stats.amplitude, hyper.D_shape, n_positions)

    # Use average concentration across positions
    cdef double avg_concentration = 0.0
    for z in range(n_positions):
        avg_concentration += hyper.concentration[z]
    avg_concentration /= n_positions
    if avg_concentration < 10.0:
        avg_concentration = 100.0  # Default

    # Compute log-likelihoods
    stats.log_lik_modern = beta_binomial_log_lik(
        k_combined, n_combined, n_positions, mean_modern, avg_concentration
    )
    stats.log_lik_ancient = beta_binomial_log_lik(
        k_combined, n_combined, n_positions, mean_ancient, avg_concentration
    )

    # Log Bayes factor
    stats.log_bf = stats.log_lik_ancient - stats.log_lik_modern

    # Posterior P(ancient | data) with prior rho
    # P(ancient | data) = BF × rho / (BF × rho + (1 - rho))
    # In log space: logit(P) = log_bf + logit(rho)
    cdef double log_odds_prior = log(hyper.rho) - log(1.0 - hyper.rho)
    cdef double log_odds_posterior = stats.log_bf + log_odds_prior

    # Convert to probability with numerical stability
    if log_odds_posterior > 20:
        stats.p_ancient = 1.0 - 1e-9
    elif log_odds_posterior < -20:
        stats.p_ancient = 1e-9
    else:
        stats.p_ancient = 1.0 / (1.0 + exp(-log_odds_posterior))


# =============================================================================
# Main Entry Points
# =============================================================================

cdef RefDamageStats* allocate_ref_damage_stats(uint32_t n_refs) noexcept nogil:
    """Allocate and initialize array of per-reference damage stats."""
    cdef RefDamageStats* stats = <RefDamageStats*>calloc(n_refs, sizeof(RefDamageStats))
    return stats


cdef void free_ref_damage_stats(RefDamageStats* stats) noexcept nogil:
    """Free per-reference damage stats array."""
    if stats != NULL:
        free(stats)


cdef DamageModelHyperparams* create_hyperparams(
    PMDCurve* curve,
    double global_baseline,
    double baseline_strength,
    double amplitude_mu,
    double amplitude_sigma,
    double rho,
    double concentration
) noexcept nogil:
    """Create and initialize hyperparameters from fitted PMD curve."""
    cdef DamageModelHyperparams* hyper = <DamageModelHyperparams*>calloc(1, sizeof(DamageModelHyperparams))
    if hyper == NULL:
        return NULL

    # Baseline prior from global estimate
    # Set α₀, β₀ such that mean = global_baseline, strength = baseline_strength
    hyper.baseline_alpha0 = global_baseline * baseline_strength
    hyper.baseline_beta0 = (1.0 - global_baseline) * baseline_strength

    # Amplitude prior (LogNormal)
    hyper.amplitude_mu = amplitude_mu
    hyper.amplitude_sigma = amplitude_sigma

    # Ancient prior
    hyper.rho = rho

    # Initialize D_shape from curve (normalized so D_shape[0] = 1.0)
    cdef double D_max = curve.D_5p_noncpg[0] if curve != NULL else 0.3
    cdef int z

    if D_max < 0.01:
        D_max = 0.3  # Default

    for z in range(NUM_POSITIONS):
        if curve != NULL:
            # Average of 5' and 3' damage
            hyper.D_shape[z] = (curve.D_5p_noncpg[z] + curve.D_3p_noncpg[z]) / (2.0 * D_max)
        else:
            # Default exponential decay
            hyper.D_shape[z] = cpow(0.7, z)

        # Set concentration (can be position-dependent if needed)
        hyper.concentration[z] = concentration

    return hyper


cdef void accumulate_alignment_damage(
    RefDamageStats* ref_stats,
    AlignmentCore* core,
    DamageCounts* dmg,
    double phi_weight,
    bint is_single_stranded
) noexcept nogil:
    """Accumulate damage counts from a single alignment weighted by phi.

    NOTE: This is a simplified version. Full implementation needs to extract
    per-position C→T and G→A counts from the MD tag, which requires access
    to the original sequence data or pre-computed per-position stats.

    For now, we use the summary stats (ct_5p_count, ga_3p_count) distributed
    across the first 8 positions.
    """
    # These are approximations - ideally we'd have per-position counts
    cdef int ct_5p = dmg.ct_5p_count
    cdef int ga_3p = dmg.ga_3p_count
    cdef int c_at_5p = dmg.c_at_5p_count
    cdef int g_at_3p = dmg.g_at_3p_count

    # Distribute counts across first 8 positions (approximation)
    # In reality, we'd want actual per-position counts from MD parsing
    cdef double avg_ct = <double>ct_5p / 8.0
    cdef double avg_ga = <double>ga_3p / 8.0
    cdef double avg_c = <double>c_at_5p / 8.0
    cdef double avg_g = <double>g_at_3p / 8.0

    cdef int z
    for z in range(8):
        ref_stats.n_5p[z] += avg_c * phi_weight
        ref_stats.k_5p[z] += avg_ct * phi_weight

        if not is_single_stranded:
            ref_stats.n_3p[z] += avg_g * phi_weight
            ref_stats.k_3p[z] += avg_ga * phi_weight

    ref_stats.total_weight += phi_weight
    ref_stats.n_alignments += 1


cdef void fit_reference_damage_model(
    RefDamageStats* stats,
    DamageModelHyperparams* hyper,
    bint use_both_ends
) noexcept nogil:
    """Fit damage model for a single reference: estimate baseline, amplitude, BF."""
    # Step 1: Estimate baseline from interior positions
    estimate_baseline(stats, hyper.baseline_alpha0, hyper.baseline_beta0)

    # Step 2: Compute Bayes factor
    compute_bayes_factor(stats, hyper, use_both_ends)


# =============================================================================
# Python-Accessible Functions
# =============================================================================

def allocate_damage_stats_py(uint32_t n_refs):
    """Allocate per-reference damage stats array."""
    cdef RefDamageStats* stats = allocate_ref_damage_stats(n_refs)
    if stats == NULL:
        raise MemoryError("Failed to allocate RefDamageStats array")
    return <uintptr_t>stats


def free_damage_stats_py(uintptr_t stats_ptr):
    """Free per-reference damage stats array."""
    free_ref_damage_stats(<RefDamageStats*>stats_ptr)


def create_hyperparams_py(
    uintptr_t curve_ptr,
    double global_baseline,
    double baseline_strength,
    double amplitude_mu,
    double amplitude_sigma,
    double rho,
    double concentration
):
    """Create hyperparameters from fitted PMD curve."""
    cdef PMDCurve* curve = <PMDCurve*>curve_ptr if curve_ptr != 0 else NULL
    cdef DamageModelHyperparams* hyper = create_hyperparams(
        curve, global_baseline, baseline_strength,
        amplitude_mu, amplitude_sigma, rho, concentration
    )
    if hyper == NULL:
        raise MemoryError("Failed to allocate DamageModelHyperparams")
    return <uintptr_t>hyper


def free_hyperparams_py(uintptr_t hyper_ptr):
    """Free hyperparameters."""
    if hyper_ptr != 0:
        free(<DamageModelHyperparams*>hyper_ptr)


def fit_all_references_py(
    uintptr_t pool_ptr,
    uintptr_t stats_ptr,
    uintptr_t hyper_ptr,
    bint use_both_ends,
    int num_threads
):
    """Fit damage model for all references in parallel.

    Returns dict with summary statistics.
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef RefDamageStats* stats = <RefDamageStats*>stats_ptr
    cdef DamageModelHyperparams* hyper = <DamageModelHyperparams*>hyper_ptr

    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t j

    # Fit each reference in parallel
    with nogil:
        for j in prange(n_refs, num_threads=num_threads, schedule='dynamic'):
            if stats[j].total_weight > 0.1:  # Only fit refs with some data
                fit_reference_damage_model(&stats[j], hyper, use_both_ends)

    # Copy damage model results to pool arrays for output
    # NOTE: Do NOT overwrite pool.gamma_values here - the EM hierarchical model
    # provides superior gamma estimates with shrinkage across references.
    # The damage model's per-reference p_ancient lacks this shrinkage.
    if pool.damage_amplitude != NULL:
        for j in range(n_refs):
            pool.damage_amplitude[j] = stats[j].amplitude
    if pool.damage_baseline != NULL:
        for j in range(n_refs):
            pool.damage_baseline[j] = stats[j].baseline
    if pool.damage_log_bf != NULL:
        for j in range(n_refs):
            pool.damage_log_bf[j] = stats[j].log_bf

    # Compute summary statistics
    cdef double sum_p_ancient = 0.0
    cdef double max_p_ancient = 0.0
    cdef double min_p_ancient = 1.0
    cdef uint32_t n_ancient = 0  # p > 0.5
    cdef uint32_t n_fitted = 0

    for j in range(n_refs):
        if stats[j].total_weight > 0.1:
            n_fitted += 1
            sum_p_ancient += stats[j].p_ancient
            if stats[j].p_ancient > max_p_ancient:
                max_p_ancient = stats[j].p_ancient
            if stats[j].p_ancient < min_p_ancient:
                min_p_ancient = stats[j].p_ancient
            if stats[j].p_ancient > 0.5:
                n_ancient += 1

    bf_nogil_logf_notime(
        LOG_TAG,
        "Fitted %u references: %u ancient (p>0.5), p_ancient range=[%.4f, %.4f]",
        n_fitted, n_ancient, min_p_ancient, max_p_ancient
    )

    return {
        'n_refs': n_refs,
        'n_fitted': n_fitted,
        'n_ancient': n_ancient,
        'mean_p_ancient': sum_p_ancient / n_fitted if n_fitted > 0 else 0.5,
        'min_p_ancient': min_p_ancient,
        'max_p_ancient': max_p_ancient,
    }


def accumulate_from_em_py(
    uintptr_t pool_ptr,
    uintptr_t stats_ptr,
    bint is_single_stranded,
    int num_threads
):
    """Accumulate damage counts from all alignments weighted by EM phi values.

    This should be called after EM converges to gather per-reference damage stats.
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef RefDamageStats* stats = <RefDamageStats*>stats_ptr

    cdef uint64_t n_alns = pool.alignment_count
    cdef uint32_t n_refs = pool.reference_count
    cdef AlignmentCore* core
    cdef DamageCounts* dmg
    cdef uint64_t i
    cdef uint32_t ref_idx
    cdef double phi
    cdef float* zp_values = NULL

    # Prefer true EM posterior probabilities if available.
    # `apply_probability_filtering()` computes/stores per-alignment posteriors in
    # `pool.precomputed_zp_values`, aligned with split arrays after compaction.
    if pool.zp_values_computed and pool.precomputed_zp_values != NULL:
        zp_values = pool.precomputed_zp_values

    bf_nogil_logf_notime(LOG_TAG, "Accumulating damage stats from %llu alignments", n_alns)

    # Clear existing stats
    with nogil:
        memset(stats, 0, n_refs * sizeof(RefDamageStats))

    # Accumulate (single-threaded for now to avoid race conditions)
    with nogil:
        for i in range(n_alns):
            core = &pool.alignment_cores[i]
            dmg = &pool.damage_counts[i]
            ref_idx = core.reference_index

            if ref_idx >= n_refs:
                continue

            # Get φ (assignment probability) from stored EM posteriors (ZP).
            # Fallbacks:
            # - If ZP values are unavailable, treat each surviving alignment equally (1.0).
            if zp_values != NULL:
                phi = <double>zp_values[i]
                if phi < 0.0:
                    phi = 0.0
                elif phi > 1.0:
                    phi = 1.0
            else:
                phi = 1.0

            accumulate_alignment_damage(&stats[ref_idx], core, dmg, phi, is_single_stranded)

    # Log summary
    cdef uint32_t n_with_data = 0
    cdef double total_weight = 0.0
    for ref_idx in range(n_refs):
        if stats[ref_idx].total_weight > 0:
            n_with_data += 1
            total_weight += stats[ref_idx].total_weight

    bf_nogil_logf_notime(
        LOG_TAG,
        "Accumulated: %u refs with data, total_weight=%.1f",
        n_with_data, total_weight
    )
