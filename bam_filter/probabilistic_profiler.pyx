# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
"""
Probabilistic Taxonomic Profiler - Cython Implementation

Novel approach combining:
- DCMS: Damage-Calibrated Mixture-of-Sources Model (per-reference)
- TBN-RN: Taxonomic Bayesian Network with belief propagation
- TAD-PUT: Probabilistic abundance estimation

No LCA required - uses soft probabilistic assignments throughout.
"""

from libc.math cimport log, exp, lgamma, sqrt, fabs
from libc.string cimport memset
from libc.stdlib cimport malloc, free, calloc
from libc.stdint cimport int32_t, int64_t, uint32_t
from cython.parallel cimport prange, parallel
from cpython.pycapsule cimport PyCapsule_GetPointer

# Import graph structures for direct GMRF integration
from bam_filter.processor_graph_ops cimport WeightedGraph, GraphNode
from bam_filter.processor_igraph cimport (
    igraph_t, igraph_vector_t, igraph_vector_int_t,
    igraph_integer_t, igraph_real_t, igraph_error_t,
    igraph_vcount, igraph_ecount, igraph_neighbors,
    igraph_vector_int_init, igraph_vector_int_destroy, igraph_vector_int_size,
    igraph_vector_size,
    IGRAPH_ALL, IGRAPH_NO_LOOPS,
)

from bam_filter.probabilistic_profiler cimport (
    ProfilerHyperparams, DamageHyperparams, PresenceHyperparams,
    EpochHyperparams, NoisyORHyperparams,
    RefDamageData, RefCoverageData, RefPosteriors, RefProfileData,
    TaxonBeliefs, TaxonNode,
)

# Constants
cdef double LOG_EPSILON = -700.0  # Minimum log value to avoid underflow
cdef double EPSILON = 1e-300


# =============================================================================
# CSR (COMPRESSED SPARSE ROW) UTILITIES FOR BELIEF PROPAGATION
# =============================================================================

cdef inline int32_t binary_search_int32(int32_t* arr, int32_t n, int32_t val) noexcept nogil:
    """Binary search for val in sorted arr. Returns index or -1 if not found."""
    cdef int32_t lo = 0, hi = n - 1, mid
    while lo <= hi:
        mid = (lo + hi) >> 1
        if arr[mid] == val:
            return mid
        elif arr[mid] < val:
            lo = mid + 1
        else:
            hi = mid - 1
    return -1


cdef struct CSRGraph:
    int32_t* indices     # Flattened child indices
    int32_t* ptr         # Start pointer for each node (size n_nodes + 1)
    int32_t n_nodes
    int32_t n_edges


cdef struct CSRRefs:
    double* p_present    # Flattened p_present values
    int32_t* ptr         # Start pointer for each taxon (size n_taxa + 1)
    int32_t n_taxa
    int32_t n_refs


# =============================================================================
# HYPERPARAMETER INITIALIZATION
# =============================================================================

cdef void init_default_hyperparams(ProfilerHyperparams* params) noexcept nogil:
    """Initialize hyperparameters with default values and precompute derived values."""
    # Damage model
    params.damage.p_0 = 0.40        # Ancient end damage mean
    params.damage.p_bg = 0.005      # Background error floor
    params.damage.tau_A = 5.0       # Decay length (bases)
    params.damage.kappa_A = 50.0    # Ancient Beta concentration
    params.damage.p_err = 0.005     # Modern error rate
    params.damage.kappa_M = 200.0   # Modern Beta concentration
    params.damage.pi_A_ref = 0.10   # Prior P(ancient) per reference

    # Precompute Beta params for damage model (avoids lgamma calls in hot loops)
    precompute_damage_params(&params.damage)

    # Presence model: P(present) = sigmoid(scale * authenticity + intercept)
    # where authenticity = norm_entropy - norm_gini
    params.presence.intercept = 0.0  # Centered at authenticity = 0
    params.presence.scale = 5.0      # Steepness of sigmoid

    # Epoch model
    params.epoch.pi_0 = 0.5         # Root prior P(ancient)
    params.epoch.rho_AA = 0.90      # P(child ancient | parent ancient)
    params.epoch.rho_MA = 0.01      # P(child ancient | parent modern)

    # Noisy-OR
    params.noisy_or.lambda_leak = 0.01
    params.noisy_or.q_fail = 0.10
    params.noisy_or.q_fail_ref = 0.10


# =============================================================================
# MATH UTILITIES
# =============================================================================

cdef inline double log_add_exp(double a, double b) noexcept nogil:
    """Numerically stable log(exp(a) + exp(b))."""
    cdef double max_val, diff
    if a > b:
        max_val = a
        diff = b - a
    else:
        max_val = b
        diff = a - b

    if diff < LOG_EPSILON:
        return max_val
    return max_val + log(1.0 + exp(diff))


cdef inline double sigmoid(double x) noexcept nogil:
    """Logistic sigmoid function."""
    if x > 700.0:
        return 1.0
    elif x < -700.0:
        return 0.0
    return 1.0 / (1.0 + exp(-x))


cdef double log_beta_function(double a, double b) noexcept nogil:
    """Compute log(B(a, b)) = log(Gamma(a)) + log(Gamma(b)) - log(Gamma(a+b))."""
    return lgamma(a) + lgamma(b) - lgamma(a + b)


cdef double log_binomial_coeff(int64_t n, int64_t k) noexcept nogil:
    """Compute log(C(n, k)) = log(n!) - log(k!) - log((n-k)!)."""
    if k < 0 or k > n:
        return LOG_EPSILON
    if k == 0 or k == n:
        return 0.0
    return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)


cdef double log_beta_binomial(int64_t k, int64_t n, double a, double b) noexcept nogil:
    """
    Compute log P(k | n, a, b) for Beta-Binomial distribution.

    P(k | n, a, b) = C(n, k) * B(k+a, n-k+b) / B(a, b)
    """
    if n < 0 or k < 0 or k > n:
        return LOG_EPSILON
    if n == 0:
        return 0.0

    cdef double log_binom = log_binomial_coeff(n, k)
    cdef double log_beta_post = log_beta_function(k + a, n - k + b)
    cdef double log_beta_prior = log_beta_function(a, b)

    return log_binom + log_beta_post - log_beta_prior


# =============================================================================
# DAMAGE MODEL (DCMS)
# =============================================================================

cdef void precompute_damage_params(DamageHyperparams* params) noexcept nogil:
    """
    Precompute Beta parameters and log-Beta-function values for damage model.
    Call this once after setting hyperparameters to avoid repeated lgamma calls.
    """
    cdef int j
    cdef double m

    # Modern parameters (constant across all positions)
    params.a_M = params.kappa_M * params.p_err
    params.b_M = params.kappa_M * (1.0 - params.p_err)
    params.log_beta_prior_M = log_beta_function(params.a_M, params.b_M)

    # Ancient parameters per position (1-indexed, stored 0-indexed)
    for j in range(20):
        # m_A[j] = p_bg + (p_0 - p_bg) * exp(-(j)/tau_A) for j=0..19 (position j+1)
        m = params.p_bg + (params.p_0 - params.p_bg) * exp(-(<double>j) / params.tau_A)
        params.a_A[j] = params.kappa_A * m
        params.b_A[j] = params.kappa_A * (1.0 - m)
        params.log_beta_prior_A[j] = log_beta_function(params.a_A[j], params.b_A[j])


def estimate_sample_background(
    double[:, ::1] damage_5p_n,
    double[:, ::1] damage_5p_k,
    double[:, ::1] damage_3p_n,
    double[:, ::1] damage_3p_k,
    double min_bg=0.005,
    double max_bg=0.10,
):
    """
    Estimate sample-specific background error rate from positions 15-20.

    The decay of ancient damage is negligible at positions 15-20 (exp(-14/5) ≈ 0.06),
    so the observed damage rate there approximates the background sequencing error
    plus any modern contamination signal.

    Parameters
    ----------
    damage_5p_n, damage_5p_k : 2D array (n_refs, 20)
        5' damage counts (n=total, k=damage events)
    damage_3p_n, damage_3p_k : 2D array (n_refs, 20)
        3' damage counts
    min_bg : float
        Minimum background rate (default 0.5%, prevents too permissive model)
    max_bg : float
        Maximum background rate (default 10%, prevents too strict model)

    Returns
    -------
    dict with:
        'p_bg': estimated background rate
        'p_err': estimated modern error rate (same as p_bg)
        'n_obs': total observations used for estimation
    """
    cdef int32_t n_refs = damage_5p_n.shape[0]
    cdef int32_t i, j
    cdef double total_n = 0.0, total_k = 0.0
    cdef double n_val, k_val
    cdef double bg_rate

    # Sum counts from positions 15-20 (indices 14-19)
    for i in range(n_refs):
        for j in range(14, 20):
            n_val = damage_5p_n[i, j]
            k_val = damage_5p_k[i, j]
            if n_val > 0:
                total_n += n_val
                total_k += k_val

            n_val = damage_3p_n[i, j]
            k_val = damage_3p_k[i, j]
            if n_val > 0:
                total_n += n_val
                total_k += k_val

    if total_n < 100:
        # Not enough data for reliable estimation, use default
        return {'p_bg': 0.005, 'p_err': 0.005, 'n_obs': int(total_n)}

    bg_rate = total_k / total_n

    # Clamp to reasonable bounds
    if bg_rate < min_bg:
        bg_rate = min_bg
    elif bg_rate > max_bg:
        bg_rate = max_bg

    return {'p_bg': bg_rate, 'p_err': bg_rate, 'n_obs': int(total_n)}


cdef inline double log_beta_binomial_precomputed(
    int64_t k, int64_t n, double a, double b, double log_beta_prior
) noexcept nogil:
    """
    Fast Beta-Binomial log-likelihood using precomputed log B(a,b).
    Avoids 2 lgamma calls per invocation.
    """
    if n < 0 or k < 0 or k > n:
        return LOG_EPSILON
    if n == 0:
        return 0.0

    cdef double log_binom = lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1)
    cdef double log_beta_post = lgamma(k + a) + lgamma(n - k + b) - lgamma(n + a + b)

    return log_binom + log_beta_post - log_beta_prior


cdef double compute_damage_log_bf(
    RefDamageData* damage,
    DamageHyperparams* params
) noexcept nogil:
    """
    Compute log Bayes factor for ancient vs modern damage pattern.
    Uses precomputed Beta parameters for efficiency.

    log_BF = sum_j [log P(k_5p[j] | n_5p[j], ancient) - log P(k_5p[j] | n_5p[j], modern)]
           + sum_j [log P(k_3p[j] | n_3p[j], ancient) - log P(k_3p[j] | n_3p[j], modern)]
    """
    cdef double log_bf = 0.0
    cdef double log_lik_ancient, log_lik_modern
    cdef int j
    cdef int64_t k, n

    # 5' end (C->T)
    for j in range(20):
        n = <int64_t>damage.n_5p[j]
        if n > 0:
            k = <int64_t>damage.k_5p[j]
            log_lik_ancient = log_beta_binomial_precomputed(
                k, n, params.a_A[j], params.b_A[j], params.log_beta_prior_A[j]
            )
            log_lik_modern = log_beta_binomial_precomputed(
                k, n, params.a_M, params.b_M, params.log_beta_prior_M
            )
            log_bf += log_lik_ancient - log_lik_modern

    # 3' end (G->A) - same exponential decay from the other end
    for j in range(20):
        n = <int64_t>damage.n_3p[j]
        if n > 0:
            k = <int64_t>damage.k_3p[j]
            log_lik_ancient = log_beta_binomial_precomputed(
                k, n, params.a_A[j], params.b_A[j], params.log_beta_prior_A[j]
            )
            log_lik_modern = log_beta_binomial_precomputed(
                k, n, params.a_M, params.b_M, params.log_beta_prior_M
            )
            log_bf += log_lik_ancient - log_lik_modern

    return log_bf


cdef double compute_p_ancient(double log_bf, double pi_A_ref) noexcept nogil:
    """
    Compute P(ancient | reference) from log Bayes factor and prior.

    P(ancient | r) = sigmoid(log_BF + prior_log_odds)
    where prior_log_odds = log(pi_A_ref / (1 - pi_A_ref))
    """
    cdef double prior_log_odds
    if pi_A_ref <= 0.0:
        prior_log_odds = -700.0
    elif pi_A_ref >= 1.0:
        prior_log_odds = 700.0
    else:
        prior_log_odds = log(pi_A_ref / (1.0 - pi_A_ref))

    return sigmoid(log_bf + prior_log_odds)


# =============================================================================
# PRESENCE MODEL
# =============================================================================

cdef double compute_p_present(
    RefCoverageData* coverage,
    PresenceHyperparams* params
) noexcept nogil:
    """
    Compute P(present | reference) using authenticity score.

    authenticity = norm_entropy - norm_gini
    P(present|r) = sigmoid(scale * authenticity + intercept)

    High authenticity (uniform coverage) -> high P(present)
    Low authenticity (clustered coverage) -> low P(present)
    """
    cdef double eta = params.scale * coverage.authenticity + params.intercept
    return sigmoid(eta)


# =============================================================================
# PER-REFERENCE COMPUTATION
# =============================================================================

cdef void compute_ref_posteriors(
    RefProfileData* ref,
    ProfilerHyperparams* params
) noexcept nogil:
    """Compute all posteriors for a reference."""
    # Damage model -> P(ancient|r)
    ref.posteriors.log_bf = compute_damage_log_bf(&ref.damage, &params.damage)
    ref.posteriors.p_ancient = compute_p_ancient(ref.posteriors.log_bf, params.damage.pi_A_ref)

    # Coverage model -> P(present|r)
    ref.posteriors.p_present = compute_p_present(&ref.coverage, &params.presence)

    # TAD split
    ref.posteriors.tad_ancient = ref.coverage.tad * ref.posteriors.p_ancient * ref.posteriors.p_present
    ref.posteriors.tad_modern = ref.coverage.tad * (1.0 - ref.posteriors.p_ancient) * ref.posteriors.p_present


# =============================================================================
# PARALLEL BATCH PROCESSING
# =============================================================================

cdef void compute_batch_posteriors_parallel(
    RefProfileData* refs,
    int32_t n_refs,
    ProfilerHyperparams* params,
    int num_threads
) noexcept nogil:
    """Compute posteriors for all references in parallel using OpenMP."""
    cdef int32_t i

    for i in prange(n_refs, nogil=True, num_threads=num_threads, schedule='static'):
        compute_ref_posteriors(&refs[i], params)


# =============================================================================
# BELIEF PROPAGATION
# =============================================================================

cdef void init_taxon_beliefs(TaxonBeliefs* beliefs) noexcept nogil:
    """Initialize taxon beliefs to neutral state."""
    beliefs.psi_ancient = 1.0
    beliefs.psi_modern = 1.0
    beliefs.n_refs = 0

    beliefs.b_ancient = 1.0
    beliefs.b_modern = 1.0

    beliefs.msg_up_ancient = 1.0
    beliefs.msg_up_modern = 1.0

    beliefs.msg_down_ancient = 1.0
    beliefs.msg_down_modern = 1.0

    beliefs.p_ancient = 0.5
    beliefs.p_modern = 0.5
    beliefs.p_present = 0.0

    beliefs.tad_ancient = 0.0
    beliefs.tad_modern = 0.0
    beliefs.tad_total = 0.0


cdef void accumulate_ref_observation(
    TaxonBeliefs* beliefs,
    double p_ancient_ref
) noexcept nogil:
    """
    Accumulate observation factor from a reference.

    Psi_t(A) = product of P(ancient|r) for refs in taxon
    Psi_t(M) = product of (1 - P(ancient|r)) for refs in taxon
    """
    # Clamp to avoid zeros
    cdef double p_anc = p_ancient_ref
    if p_anc < EPSILON:
        p_anc = EPSILON
    elif p_anc > 1.0 - EPSILON:
        p_anc = 1.0 - EPSILON

    beliefs.psi_ancient *= p_anc
    beliefs.psi_modern *= (1.0 - p_anc)
    beliefs.n_refs += 1


cdef void compute_local_belief(
    TaxonBeliefs* beliefs,
    EpochHyperparams* params,
    bint is_root
) noexcept nogil:
    """
    Compute local belief combining prior, observations, and child messages.

    b_t(e) = phi_t(e) * Psi_t(e) * product of child upward messages
    """
    cdef double phi_ancient, phi_modern

    if is_root:
        phi_ancient = params.pi_0
        phi_modern = 1.0 - params.pi_0
    else:
        phi_ancient = 1.0
        phi_modern = 1.0

    beliefs.b_ancient = phi_ancient * beliefs.psi_ancient * beliefs.msg_up_ancient
    beliefs.b_modern = phi_modern * beliefs.psi_modern * beliefs.msg_up_modern

    # Normalize to prevent underflow
    cdef double total = beliefs.b_ancient + beliefs.b_modern
    if total > 0.0:
        beliefs.b_ancient /= total
        beliefs.b_modern /= total
    else:
        beliefs.b_ancient = 0.5
        beliefs.b_modern = 0.5


cdef void compute_upward_message(
    TaxonBeliefs* child,
    double* msg_ancient,
    double* msg_modern,
    EpochHyperparams* params
) noexcept nogil:
    """
    Compute upward message from child to parent.

    mu_{c->p}(A) = rho_AA * b_c(A) + (1 - rho_AA) * b_c(M)
    mu_{c->p}(M) = rho_MA * b_c(A) + (1 - rho_MA) * b_c(M)
    """
    msg_ancient[0] = params.rho_AA * child.b_ancient + (1.0 - params.rho_AA) * child.b_modern
    msg_modern[0] = params.rho_MA * child.b_ancient + (1.0 - params.rho_MA) * child.b_modern

    # Normalize
    cdef double total = msg_ancient[0] + msg_modern[0]
    if total > 0.0:
        msg_ancient[0] /= total
        msg_modern[0] /= total
    else:
        msg_ancient[0] = 0.5
        msg_modern[0] = 0.5


cdef void receive_upward_message(
    TaxonBeliefs* parent,
    double msg_ancient,
    double msg_modern
) noexcept nogil:
    """Multiply incoming upward message into parent's accumulated messages."""
    parent.msg_up_ancient *= msg_ancient
    parent.msg_up_modern *= msg_modern

    # Periodic renormalization to prevent underflow
    cdef double total = parent.msg_up_ancient + parent.msg_up_modern
    if total > 0.0 and total < 1e-100:
        parent.msg_up_ancient /= total
        parent.msg_up_modern /= total


cdef void compute_downward_message(
    TaxonBeliefs* parent,
    TaxonBeliefs* child,
    double* msg_ancient,
    double* msg_modern,
    EpochHyperparams* params,
    bint is_root
) noexcept nogil:
    """
    Compute downward message from parent to child.

    First compute parent's context belief excluding this child:
    b_p^{\\c}(e) = phi_p(e) * Psi_p(e) * msg_from_parent * product of other children

    Then:
    mu_{p->c}(A) = rho_AA * b_p^{\\c}(A) + rho_MA * b_p^{\\c}(M)
    mu_{p->c}(M) = (1 - rho_AA) * b_p^{\\c}(A) + (1 - rho_MA) * b_p^{\\c}(M)
    """
    cdef double phi_ancient, phi_modern
    cdef double b_exc_ancient, b_exc_modern
    cdef double child_msg_ancient, child_msg_modern

    if is_root:
        phi_ancient = params.pi_0
        phi_modern = 1.0 - params.pi_0
    else:
        phi_ancient = 1.0
        phi_modern = 1.0

    # Get the upward message that this child contributed
    compute_upward_message(child, &child_msg_ancient, &child_msg_modern, params)

    # Parent's context excluding this child
    # We need to divide out this child's contribution from msg_up
    # b_p^{\c} = phi * psi * msg_from_parent * (msg_up / child_msg)
    cdef double msg_up_exc_ancient = parent.msg_up_ancient
    cdef double msg_up_exc_modern = parent.msg_up_modern

    if child_msg_ancient > EPSILON:
        msg_up_exc_ancient /= child_msg_ancient
    if child_msg_modern > EPSILON:
        msg_up_exc_modern /= child_msg_modern

    b_exc_ancient = phi_ancient * parent.psi_ancient * parent.msg_down_ancient * msg_up_exc_ancient
    b_exc_modern = phi_modern * parent.psi_modern * parent.msg_down_modern * msg_up_exc_modern

    # Normalize
    cdef double total = b_exc_ancient + b_exc_modern
    if total > 0.0:
        b_exc_ancient /= total
        b_exc_modern /= total
    else:
        b_exc_ancient = 0.5
        b_exc_modern = 0.5

    # Compute downward message
    msg_ancient[0] = params.rho_AA * b_exc_ancient + params.rho_MA * b_exc_modern
    msg_modern[0] = (1.0 - params.rho_AA) * b_exc_ancient + (1.0 - params.rho_MA) * b_exc_modern

    # Normalize
    total = msg_ancient[0] + msg_modern[0]
    if total > 0.0:
        msg_ancient[0] /= total
        msg_modern[0] /= total
    else:
        msg_ancient[0] = 0.5
        msg_modern[0] = 0.5


cdef void receive_downward_message(
    TaxonBeliefs* child,
    double msg_ancient,
    double msg_modern
) noexcept nogil:
    """Store incoming downward message."""
    child.msg_down_ancient = msg_ancient
    child.msg_down_modern = msg_modern


cdef void compute_marginals(
    TaxonBeliefs* beliefs,
    EpochHyperparams* params,
    bint is_root
) noexcept nogil:
    """
    Compute final marginals for P(E_t = ancient) and P(E_t = modern).

    P(E_t = e) proportional to phi_t(e) * Psi_t(e) * msg_from_parent * product of child msgs
    """
    cdef double phi_ancient, phi_modern
    cdef double unnorm_ancient, unnorm_modern, total

    if is_root:
        phi_ancient = params.pi_0
        phi_modern = 1.0 - params.pi_0
    else:
        phi_ancient = 1.0
        phi_modern = 1.0

    unnorm_ancient = phi_ancient * beliefs.psi_ancient * beliefs.msg_down_ancient * beliefs.msg_up_ancient
    unnorm_modern = phi_modern * beliefs.psi_modern * beliefs.msg_down_modern * beliefs.msg_up_modern

    total = unnorm_ancient + unnorm_modern
    if total > 0.0:
        beliefs.p_ancient = unnorm_ancient / total
        beliefs.p_modern = unnorm_modern / total
    else:
        beliefs.p_ancient = 0.5
        beliefs.p_modern = 0.5


# =============================================================================
# NOISY-OR PRESENCE
# =============================================================================

cdef double compute_noisy_or_presence_from_refs(
    double* p_present_refs,
    int32_t n_refs,
    NoisyORHyperparams* params
) noexcept nogil:
    """
    Compute P(present) for a leaf taxon from its references using Noisy-OR.

    P(present = 0) = (1 - lambda_leak) * product_r [1 - P(present|r) * (1 - q_fail_ref)]
    P(present = 1) = 1 - P(present = 0)
    """
    cdef double p_absent = 1.0 - params.lambda_leak
    cdef int i

    for i in range(n_refs):
        p_absent *= (1.0 - p_present_refs[i] * (1.0 - params.q_fail_ref))

    return 1.0 - p_absent


cdef double compute_noisy_or_presence_from_children(
    double* p_present_children,
    int32_t n_children,
    NoisyORHyperparams* params
) noexcept nogil:
    """
    Compute P(present) for an internal taxon from its children using Noisy-OR.

    P(present = 0) = (1 - lambda_leak) * product_c [1 - P(present_c) * (1 - q_fail)]
    P(present = 1) = 1 - P(present = 0)
    """
    cdef double p_absent = 1.0 - params.lambda_leak
    cdef int i

    for i in range(n_children):
        p_absent *= (1.0 - p_present_children[i] * (1.0 - params.q_fail))

    return 1.0 - p_absent


# =============================================================================
# TAD AGGREGATION
# =============================================================================

cdef void aggregate_tad_from_refs(
    TaxonBeliefs* beliefs,
    RefProfileData* refs,
    int32_t n_refs
) noexcept nogil:
    """Aggregate TAD from references."""
    cdef int i
    beliefs.tad_ancient = 0.0
    beliefs.tad_modern = 0.0

    for i in range(n_refs):
        beliefs.tad_ancient += refs[i].posteriors.tad_ancient
        beliefs.tad_modern += refs[i].posteriors.tad_modern

    beliefs.tad_total = beliefs.tad_ancient + beliefs.tad_modern


cdef void aggregate_tad_from_children(
    TaxonBeliefs* parent,
    TaxonBeliefs* children,
    int32_t n_children
) noexcept nogil:
    """Add children's TAD to parent's TAD."""
    cdef int i

    for i in range(n_children):
        parent.tad_ancient += children[i].tad_ancient
        parent.tad_modern += children[i].tad_modern

    parent.tad_total = parent.tad_ancient + parent.tad_modern


# =============================================================================
# BATCH COMPUTATION WITH NUMPY INTERFACE
# =============================================================================

def compute_batch_posteriors(
    double[:, ::1] damage_5p_n,
    double[:, ::1] damage_5p_k,
    double[:, ::1] damage_3p_n,
    double[:, ::1] damage_3p_k,
    double[::1] breadth,
    double[::1] mean_depth,
    double[::1] wcb,
    double[::1] norm_entropy,
    double[::1] norm_gini,
    long[::1] n_reads,
    double[::1] tad,
    int num_threads=4,
    dict hyperparams=None,
    bint adaptive_background=True,
):
    """
    Compute posteriors for batch of references using OpenMP parallelization.

    Parameters
    ----------
    damage_5p_n, damage_5p_k : 2D array (n_refs, 20)
        5' damage counts
    damage_3p_n, damage_3p_k : 2D array (n_refs, 20)
        3' damage counts
    breadth, mean_depth, wcb : 1D array (n_refs,)
        Coverage statistics
    norm_entropy, norm_gini : 1D array (n_refs,)
        Normalized spatial entropy and Gini coefficient [0, 1]
    n_reads : 1D array (n_refs,)
        Read counts
    tad : 1D array (n_refs,)
        Truncated average depth
    num_threads : int
        Number of OpenMP threads
    hyperparams : dict, optional
        Custom hyperparameters
    adaptive_background : bool
        If True (default), estimate p_bg and p_err from sample positions 15-20.
        This provides a statistically sound baseline for the modern hypothesis.
        The decay length tau remains fixed at 5.0 (literature value).

    Returns
    -------
    dict with arrays: 'log_bf', 'p_ancient', 'p_present', 'tad_ancient', 'tad_modern'
    Also includes 'estimated_background' if adaptive_background=True.
    """
    import numpy as np

    cdef int32_t n_refs = damage_5p_n.shape[0]
    cdef int32_t i, j
    cdef ProfilerHyperparams params
    cdef RefProfileData* refs = NULL
    estimated_bg = None

    init_default_hyperparams(&params)

    # Sample-adaptive background estimation (Empirical Bayes)
    # Estimates p_bg and p_err from positions 15-20 where ancient decay is negligible
    if adaptive_background:
        estimated_bg = estimate_sample_background(
            damage_5p_n, damage_5p_k, damage_3p_n, damage_3p_k
        )
        params.damage.p_bg = estimated_bg['p_bg']
        params.damage.p_err = estimated_bg['p_err']
        # tau remains at literature value (5.0) - NOT fitted from sample
        precompute_damage_params(&params.damage)

    if hyperparams is not None:
        if 'damage' in hyperparams:
            d = hyperparams['damage']
            if 'p_0' in d: params.damage.p_0 = d['p_0']
            if 'p_bg' in d: params.damage.p_bg = d['p_bg']
            if 'tau_A' in d: params.damage.tau_A = d['tau_A']
            if 'kappa_A' in d: params.damage.kappa_A = d['kappa_A']
            if 'p_err' in d: params.damage.p_err = d['p_err']
            if 'kappa_M' in d: params.damage.kappa_M = d['kappa_M']
            if 'pi_A_ref' in d: params.damage.pi_A_ref = d['pi_A_ref']
            # Recompute precomputed values after overriding damage params
            precompute_damage_params(&params.damage)
        if 'presence' in hyperparams:
            p = hyperparams['presence']
            if 'intercept' in p: params.presence.intercept = p['intercept']
            if 'scale' in p: params.presence.scale = p['scale']

    refs = <RefProfileData*>calloc(n_refs, sizeof(RefProfileData))
    if refs == NULL:
        raise MemoryError("Failed to allocate reference data")

    try:
        for i in range(n_refs):
            for j in range(20):
                refs[i].damage.n_5p[j] = damage_5p_n[i, j]
                refs[i].damage.k_5p[j] = damage_5p_k[i, j]
                refs[i].damage.n_3p[j] = damage_3p_n[i, j]
                refs[i].damage.k_3p[j] = damage_3p_k[i, j]

            refs[i].coverage.breadth = breadth[i]
            refs[i].coverage.mean_depth = mean_depth[i]
            refs[i].coverage.wcb = wcb[i]
            refs[i].coverage.norm_entropy = norm_entropy[i]
            refs[i].coverage.norm_gini = norm_gini[i]
            refs[i].coverage.authenticity = norm_entropy[i] - norm_gini[i]
            refs[i].coverage.n_reads = n_reads[i]
            refs[i].coverage.tad = tad[i]

        with nogil:
            compute_batch_posteriors_parallel(refs, n_refs, &params, num_threads)

        log_bf_out = np.empty(n_refs, dtype=np.float64)
        p_ancient_out = np.empty(n_refs, dtype=np.float64)
        p_present_out = np.empty(n_refs, dtype=np.float64)
        tad_ancient_out = np.empty(n_refs, dtype=np.float64)
        tad_modern_out = np.empty(n_refs, dtype=np.float64)

        for i in range(n_refs):
            log_bf_out[i] = refs[i].posteriors.log_bf
            p_ancient_out[i] = refs[i].posteriors.p_ancient
            p_present_out[i] = refs[i].posteriors.p_present
            tad_ancient_out[i] = refs[i].posteriors.tad_ancient
            tad_modern_out[i] = refs[i].posteriors.tad_modern

        result = {
            'log_bf': log_bf_out,
            'p_ancient': p_ancient_out,
            'p_present': p_present_out,
            'tad_ancient': tad_ancient_out,
            'tad_modern': tad_modern_out,
        }
        if estimated_bg is not None:
            result['estimated_background'] = estimated_bg
        return result
    finally:
        if refs != NULL:
            free(refs)


# =============================================================================
# PYTHON INTERFACE
# =============================================================================

def create_default_hyperparams():
    """Create default hyperparameters as a Python dictionary."""
    cdef ProfilerHyperparams params
    init_default_hyperparams(&params)

    return {
        'damage': {
            'p_0': params.damage.p_0,
            'p_bg': params.damage.p_bg,
            'tau_A': params.damage.tau_A,
            'kappa_A': params.damage.kappa_A,
            'p_err': params.damage.p_err,
            'kappa_M': params.damage.kappa_M,
            'pi_A_ref': params.damage.pi_A_ref,
        },
        'presence': {
            'intercept': params.presence.intercept,
            'scale': params.presence.scale,
        },
        'epoch': {
            'pi_0': params.epoch.pi_0,
            'rho_AA': params.epoch.rho_AA,
            'rho_MA': params.epoch.rho_MA,
        },
        'noisy_or': {
            'lambda_leak': params.noisy_or.lambda_leak,
            'q_fail': params.noisy_or.q_fail,
            'q_fail_ref': params.noisy_or.q_fail_ref,
        },
    }


def compute_reference_posteriors(
    dict damage_counts,
    dict coverage_stats,
    dict hyperparams=None,
):
    """
    Compute posteriors for a single reference.

    Parameters
    ----------
    damage_counts : dict
        Keys: 'n_5p', 'k_5p', 'n_3p', 'k_3p' (each a list of 20 floats)
    coverage_stats : dict
        Keys: 'breadth', 'mean_depth', 'wcb', 'norm_entropy', 'norm_gini', 'n_reads', 'tad'
    hyperparams : dict, optional
        Custom hyperparameters (uses defaults if not provided)

    Returns
    -------
    dict with keys: 'log_bf', 'p_ancient', 'p_present', 'tad_ancient', 'tad_modern'
    """
    cdef ProfilerHyperparams params
    cdef RefProfileData ref
    cdef int i

    init_default_hyperparams(&params)

    # Override with custom hyperparams if provided
    if hyperparams is not None:
        if 'damage' in hyperparams:
            d = hyperparams['damage']
            if 'p_0' in d: params.damage.p_0 = d['p_0']
            if 'p_bg' in d: params.damage.p_bg = d['p_bg']
            if 'tau_A' in d: params.damage.tau_A = d['tau_A']
            if 'kappa_A' in d: params.damage.kappa_A = d['kappa_A']
            if 'p_err' in d: params.damage.p_err = d['p_err']
            if 'kappa_M' in d: params.damage.kappa_M = d['kappa_M']
            if 'pi_A_ref' in d: params.damage.pi_A_ref = d['pi_A_ref']
            precompute_damage_params(&params.damage)

        if 'presence' in hyperparams:
            p = hyperparams['presence']
            if 'intercept' in p: params.presence.intercept = p['intercept']
            if 'scale' in p: params.presence.scale = p['scale']

    # Fill damage data
    n_5p = damage_counts.get('n_5p', [0.0] * 20)
    k_5p = damage_counts.get('k_5p', [0.0] * 20)
    n_3p = damage_counts.get('n_3p', [0.0] * 20)
    k_3p = damage_counts.get('k_3p', [0.0] * 20)

    for i in range(20):
        ref.damage.n_5p[i] = n_5p[i] if i < len(n_5p) else 0.0
        ref.damage.k_5p[i] = k_5p[i] if i < len(k_5p) else 0.0
        ref.damage.n_3p[i] = n_3p[i] if i < len(n_3p) else 0.0
        ref.damage.k_3p[i] = k_3p[i] if i < len(k_3p) else 0.0

    # Fill coverage data
    ref.coverage.breadth = coverage_stats.get('breadth', 0.0)
    ref.coverage.mean_depth = coverage_stats.get('mean_depth', 0.0)
    ref.coverage.wcb = coverage_stats.get('wcb', 0.0)
    ref.coverage.norm_entropy = coverage_stats.get('norm_entropy', coverage_stats.get('norm_spatial_entropy', 0.0))
    ref.coverage.norm_gini = coverage_stats.get('norm_gini', 0.0)
    ref.coverage.authenticity = ref.coverage.norm_entropy - ref.coverage.norm_gini
    ref.coverage.n_reads = coverage_stats.get('n_reads', 0)
    ref.coverage.tad = coverage_stats.get('tad', 0.0)

    # Compute posteriors
    compute_ref_posteriors(&ref, &params)

    return {
        'log_bf': ref.posteriors.log_bf,
        'p_ancient': ref.posteriors.p_ancient,
        'p_present': ref.posteriors.p_present,
        'tad_ancient': ref.posteriors.tad_ancient,
        'tad_modern': ref.posteriors.tad_modern,
    }


# =============================================================================
# CSR-BASED BELIEF PROPAGATION (pure C arrays, fully nogil)
# =============================================================================

cdef void _bp_bottom_up_csr(
    int32_t n_taxa,
    int32_t root_idx,
    int32_t* postorder,
    int32_t* parent_idx_arr,
    double* log_psi_anc,
    double* log_psi_mod,
    double* log_msg_up_anc,
    double* log_msg_up_mod,
    double log_pi_0,
    double log_1_minus_pi_0,
    double rho_AA,
    double rho_MA,
) noexcept nogil:
    """Bottom-up pass: compute upward messages in log-space."""
    cdef int32_t i, idx, parent_idx
    cdef double log_phi_ancient, log_phi_modern
    cdef double log_b_ancient, log_b_modern, log_total
    cdef double b_ancient, b_modern
    cdef double msg_ancient, msg_modern, total

    for i in range(n_taxa):
        idx = postorder[i]
        if idx == root_idx:
            log_phi_ancient = log_pi_0
            log_phi_modern = log_1_minus_pi_0
        else:
            log_phi_ancient = 0.0
            log_phi_modern = 0.0

        log_b_ancient = log_phi_ancient + log_psi_anc[idx] + log_msg_up_anc[idx]
        log_b_modern = log_phi_modern + log_psi_mod[idx] + log_msg_up_mod[idx]

        log_total = log_add_exp(log_b_ancient, log_b_modern)
        log_b_ancient -= log_total
        log_b_modern -= log_total

        b_ancient = exp(log_b_ancient)
        b_modern = exp(log_b_modern)

        if idx != root_idx:
            parent_idx = parent_idx_arr[idx]
            if parent_idx >= 0:
                msg_ancient = rho_AA * b_ancient + (1.0 - rho_AA) * b_modern
                msg_modern = rho_MA * b_ancient + (1.0 - rho_MA) * b_modern

                total = msg_ancient + msg_modern
                if total > 0:
                    msg_ancient /= total
                    msg_modern /= total

                if msg_ancient > 0:
                    log_msg_up_anc[parent_idx] += log(msg_ancient)
                else:
                    log_msg_up_anc[parent_idx] += LOG_EPSILON
                if msg_modern > 0:
                    log_msg_up_mod[parent_idx] += log(msg_modern)
                else:
                    log_msg_up_mod[parent_idx] += LOG_EPSILON


cdef void _bp_compute_beliefs_csr(
    int32_t n_taxa,
    int32_t root_idx,
    int32_t* postorder,
    double* log_psi_anc,
    double* log_psi_mod,
    double* log_msg_up_anc,
    double* log_msg_up_mod,
    double* log_b_anc_view,
    double* log_b_mod_view,
    double log_pi_0,
    double log_1_minus_pi_0,
) noexcept nogil:
    """Compute normalized log-beliefs for each taxon."""
    cdef int32_t i, idx
    cdef double log_phi_ancient, log_phi_modern
    cdef double log_b_ancient, log_b_modern, log_total

    for i in range(n_taxa):
        idx = postorder[i]
        if idx == root_idx:
            log_phi_ancient = log_pi_0
            log_phi_modern = log_1_minus_pi_0
        else:
            log_phi_ancient = 0.0
            log_phi_modern = 0.0

        log_b_ancient = log_phi_ancient + log_psi_anc[idx] + log_msg_up_anc[idx]
        log_b_modern = log_phi_modern + log_psi_mod[idx] + log_msg_up_mod[idx]

        log_total = log_add_exp(log_b_ancient, log_b_modern)
        log_b_anc_view[idx] = log_b_ancient - log_total
        log_b_mod_view[idx] = log_b_modern - log_total


cdef void _bp_top_down_csr(
    int32_t n_taxa,
    int32_t root_idx,
    int32_t* preorder,
    int32_t* children_indices,
    int32_t* children_ptr,
    double* log_psi_anc,
    double* log_psi_mod,
    double* log_msg_up_anc,
    double* log_msg_up_mod,
    double* log_msg_down_anc,
    double* log_msg_down_mod,
    double* log_b_anc_view,
    double* log_b_mod_view,
    double log_pi_0,
    double log_1_minus_pi_0,
    double rho_AA,
    double rho_MA,
) noexcept nogil:
    """Top-down pass: compute downward messages using CSR children structure."""
    cdef int32_t i, idx, j, child_idx
    cdef double log_phi_ancient, log_phi_modern
    cdef double child_b_anc, child_b_mod
    cdef double child_msg_anc, child_msg_mod, total
    cdef double log_child_msg_anc, log_child_msg_mod
    cdef double log_msg_up_exc_anc, log_msg_up_exc_mod
    cdef double log_b_exc_anc, log_b_exc_mod, log_total
    cdef double b_exc_anc, b_exc_mod
    cdef double msg_ancient, msg_modern

    for i in range(n_taxa):
        idx = preorder[i]
        if idx == root_idx:
            log_phi_ancient = log_pi_0
            log_phi_modern = log_1_minus_pi_0
        else:
            log_phi_ancient = 0.0
            log_phi_modern = 0.0

        n_children = children_ptr[idx + 1] - children_ptr[idx]

        for j in range(children_ptr[idx], children_ptr[idx + 1]):
            child_idx = children_indices[j]

            child_b_anc = exp(log_b_anc_view[child_idx])
            child_b_mod = exp(log_b_mod_view[child_idx])

            child_msg_anc = rho_AA * child_b_anc + (1.0 - rho_AA) * child_b_mod
            child_msg_mod = rho_MA * child_b_anc + (1.0 - rho_MA) * child_b_mod

            total = child_msg_anc + child_msg_mod
            if total > 0:
                child_msg_anc /= total
                child_msg_mod /= total

            if child_msg_anc > 0:
                log_child_msg_anc = log(child_msg_anc)
            else:
                log_child_msg_anc = LOG_EPSILON
            if child_msg_mod > 0:
                log_child_msg_mod = log(child_msg_mod)
            else:
                log_child_msg_mod = LOG_EPSILON

            log_msg_up_exc_anc = log_msg_up_anc[idx] - log_child_msg_anc
            log_msg_up_exc_mod = log_msg_up_mod[idx] - log_child_msg_mod

            log_b_exc_anc = log_phi_ancient + log_psi_anc[idx] + log_msg_down_anc[idx] + log_msg_up_exc_anc
            log_b_exc_mod = log_phi_modern + log_psi_mod[idx] + log_msg_down_mod[idx] + log_msg_up_exc_mod

            log_total = log_add_exp(log_b_exc_anc, log_b_exc_mod)
            b_exc_anc = exp(log_b_exc_anc - log_total)
            b_exc_mod = exp(log_b_exc_mod - log_total)

            msg_ancient = rho_AA * b_exc_anc + rho_MA * b_exc_mod
            msg_modern = (1.0 - rho_AA) * b_exc_anc + (1.0 - rho_MA) * b_exc_mod

            total = msg_ancient + msg_modern
            if total > 0:
                msg_ancient /= total
                msg_modern /= total

            if msg_ancient > 0:
                log_msg_down_anc[child_idx] = log(msg_ancient)
            else:
                log_msg_down_anc[child_idx] = LOG_EPSILON
            if msg_modern > 0:
                log_msg_down_mod[child_idx] = log(msg_modern)
            else:
                log_msg_down_mod[child_idx] = LOG_EPSILON


cdef void _bp_compute_marginals_csr(
    int32_t n_taxa,
    int32_t root_idx,
    int32_t* postorder,
    int32_t* children_indices,
    int32_t* children_ptr,
    int32_t* taxids_arr,
    double* log_psi_anc,
    double* log_psi_mod,
    double* log_msg_up_anc,
    double* log_msg_up_mod,
    double* log_msg_down_anc,
    double* log_msg_down_mod,
    double* ref_p_present_data,
    int32_t* ref_ptr,
    int32_t* n_refs_tax,
    double* tad_anc,
    double* tad_mod,
    double* p_ancient_out,
    double* p_present_out,
    double* tad_ancient_out,
    double* tad_modern_out,
    double log_pi_0,
    double log_1_minus_pi_0,
    double lambda_leak,
    double q_fail,
    double q_fail_ref,
) noexcept nogil:
    """Compute final marginals, presence, and TAD aggregation."""
    cdef int32_t i, idx, j, child_idx
    cdef double log_phi_ancient, log_phi_modern
    cdef double log_unnorm_anc, log_unnorm_mod, log_total
    cdef double p_ancient_final
    cdef double p_absent, p_present_from_refs, p_present_from_children, p_present_val
    cdef double tad_ancient_val, tad_modern_val

    for i in range(n_taxa):
        idx = postorder[i]
        if idx == root_idx:
            log_phi_ancient = log_pi_0
            log_phi_modern = log_1_minus_pi_0
        else:
            log_phi_ancient = 0.0
            log_phi_modern = 0.0

        log_unnorm_anc = log_phi_ancient + log_psi_anc[idx] + log_msg_down_anc[idx] + log_msg_up_anc[idx]
        log_unnorm_mod = log_phi_modern + log_psi_mod[idx] + log_msg_down_mod[idx] + log_msg_up_mod[idx]

        log_total = log_add_exp(log_unnorm_anc, log_unnorm_mod)
        p_ancient_final = exp(log_unnorm_anc - log_total)

        # Noisy-OR from refs
        p_present_from_refs = 0.0
        if ref_ptr[idx + 1] > ref_ptr[idx]:
            p_absent = 1.0 - lambda_leak
            for j in range(ref_ptr[idx], ref_ptr[idx + 1]):
                p_absent *= (1.0 - ref_p_present_data[j] * (1.0 - q_fail_ref))
            p_present_from_refs = 1.0 - p_absent

        # Noisy-OR from children
        p_present_from_children = 0.0
        if children_ptr[idx + 1] > children_ptr[idx]:
            p_absent = 1.0 - lambda_leak
            for j in range(children_ptr[idx], children_ptr[idx + 1]):
                child_idx = children_indices[j]
                p_absent *= (1.0 - p_present_out[child_idx] * (1.0 - q_fail))
            p_present_from_children = 1.0 - p_absent

        if p_present_from_refs > p_present_from_children:
            p_present_val = p_present_from_refs
        else:
            p_present_val = p_present_from_children

        # TAD aggregation
        tad_ancient_val = tad_anc[idx]
        tad_modern_val = tad_mod[idx]

        for j in range(children_ptr[idx], children_ptr[idx + 1]):
            child_idx = children_indices[j]
            tad_ancient_val += tad_ancient_out[child_idx]
            tad_modern_val += tad_modern_out[child_idx]

        p_ancient_out[idx] = p_ancient_final
        p_present_out[idx] = p_present_val
        tad_ancient_out[idx] = tad_ancient_val
        tad_modern_out[idx] = tad_modern_val


# =============================================================================
# NATIVE BELIEF PROPAGATION (array-based, no Python dicts in hot paths)
# =============================================================================

def run_belief_propagation_native(
    int32_t[::1] taxids,
    int32_t[::1] parent_taxids,
    int32_t[::1] ref_taxids,
    double[::1] ref_p_ancient,
    double[::1] ref_p_present,
    double[::1] ref_tad_ancient,
    double[::1] ref_tad_modern,
    int32_t root_taxid,
    dict hyperparams=None,
    int num_threads=4,
):
    """
    Run belief propagation using native arrays (no Python dicts in hot paths).

    Parameters
    ----------
    taxids : 1D array
        All taxon IDs in the tree
    parent_taxids : 1D array
        Parent taxon ID for each taxon (same as taxid for root)
    ref_taxids : 1D array
        Taxon ID for each reference
    ref_p_ancient, ref_p_present, ref_tad_ancient, ref_tad_modern : 1D arrays
        Per-reference posteriors (from compute_batch_posteriors)
    root_taxid : int
        Root taxon ID
    hyperparams : dict, optional
        Custom hyperparameters
    num_threads : int
        Number of threads for parallel operations

    Returns
    -------
    dict: taxid -> {'p_ancient', 'p_present', 'tad_ancient', 'tad_modern', 'tad_total', 'n_refs'}
    """
    import numpy as np

    cdef ProfilerHyperparams params
    init_default_hyperparams(&params)

    if hyperparams is not None:
        if 'epoch' in hyperparams:
            e = hyperparams['epoch']
            if 'pi_0' in e: params.epoch.pi_0 = e['pi_0']
            if 'rho_AA' in e: params.epoch.rho_AA = e['rho_AA']
            if 'rho_MA' in e: params.epoch.rho_MA = e['rho_MA']
        if 'noisy_or' in hyperparams:
            n = hyperparams['noisy_or']
            if 'lambda_leak' in n: params.noisy_or.lambda_leak = n['lambda_leak']
            if 'q_fail' in n: params.noisy_or.q_fail = n['q_fail']
            if 'q_fail_ref' in n: params.noisy_or.q_fail_ref = n['q_fail_ref']

    cdef int32_t n_taxa = taxids.shape[0]
    cdef int32_t n_refs = ref_taxids.shape[0]
    cdef int32_t i, j, idx, parent_idx, child_idx
    cdef int32_t root_idx = -1

    taxid_to_idx = {}
    for i in range(n_taxa):
        taxid_to_idx[taxids[i]] = i
        if taxids[i] == root_taxid:
            root_idx = i

    if root_idx < 0:
        raise ValueError(f"Root taxid {root_taxid} not found in taxids array")

    # Build parent_idx array (direct index lookup instead of dict)
    parent_idx_arr = np.full(n_taxa, -1, dtype=np.int32)
    cdef int32_t[::1] parent_idx_view = parent_idx_arr

    # Count children per node for CSR allocation
    children_count = np.zeros(n_taxa, dtype=np.int32)
    cdef int32_t[::1] children_cnt = children_count

    for i in range(n_taxa):
        if taxids[i] != root_taxid:
            parent_t = parent_taxids[i]
            if parent_t in taxid_to_idx:
                parent_idx = taxid_to_idx[parent_t]
                parent_idx_view[i] = parent_idx
                children_cnt[parent_idx] += 1

    # Build CSR pointers for children
    children_ptr_arr = np.zeros(n_taxa + 1, dtype=np.int32)
    cdef int32_t[::1] children_ptr = children_ptr_arr
    for i in range(n_taxa):
        children_ptr[i + 1] = children_ptr[i] + children_cnt[i]

    cdef int32_t total_children = children_ptr[n_taxa]
    children_indices_arr = np.zeros(total_children, dtype=np.int32)
    cdef int32_t[::1] children_indices = children_indices_arr

    # Fill children indices (reset counts as insertion pointers)
    for i in range(n_taxa):
        children_cnt[i] = 0

    for i in range(n_taxa):
        if parent_idx_view[i] >= 0:
            parent_idx = parent_idx_view[i]
            children_indices[children_ptr[parent_idx] + children_cnt[parent_idx]] = i
            children_cnt[parent_idx] += 1

    # Count refs per taxon for CSR allocation
    refs_count = np.zeros(n_taxa, dtype=np.int32)
    cdef int32_t[::1] refs_cnt = refs_count

    for i in range(n_refs):
        ref_t = ref_taxids[i]
        if ref_t in taxid_to_idx:
            idx = taxid_to_idx[ref_t]
            refs_cnt[idx] += 1

    # Build CSR pointers for refs
    ref_ptr_arr = np.zeros(n_taxa + 1, dtype=np.int32)
    cdef int32_t[::1] ref_ptr = ref_ptr_arr
    for i in range(n_taxa):
        ref_ptr[i + 1] = ref_ptr[i] + refs_cnt[i]

    cdef int32_t total_refs_mapped = ref_ptr[n_taxa]
    ref_p_present_data_arr = np.zeros(total_refs_mapped, dtype=np.float64)
    cdef double[::1] ref_p_present_data = ref_p_present_data_arr

    # Use log-space accumulators to prevent underflow with many refs/deep trees
    log_psi_ancient = np.zeros(n_taxa, dtype=np.float64)
    log_psi_modern = np.zeros(n_taxa, dtype=np.float64)
    n_refs_per_taxon = np.zeros(n_taxa, dtype=np.int32)
    log_msg_up_ancient = np.zeros(n_taxa, dtype=np.float64)
    log_msg_up_modern = np.zeros(n_taxa, dtype=np.float64)
    log_msg_down_ancient = np.zeros(n_taxa, dtype=np.float64)
    log_msg_down_modern = np.zeros(n_taxa, dtype=np.float64)
    tad_ancient_arr = np.zeros(n_taxa, dtype=np.float64)
    tad_modern_arr = np.zeros(n_taxa, dtype=np.float64)

    cdef double[::1] log_psi_anc = log_psi_ancient
    cdef double[::1] log_psi_mod = log_psi_modern
    cdef int32_t[::1] n_refs_tax = n_refs_per_taxon
    cdef double[::1] log_msg_up_anc = log_msg_up_ancient
    cdef double[::1] log_msg_up_mod = log_msg_up_modern
    cdef double[::1] log_msg_down_anc = log_msg_down_ancient
    cdef double[::1] log_msg_down_mod = log_msg_down_modern
    cdef double[::1] tad_anc = tad_ancient_arr
    cdef double[::1] tad_mod = tad_modern_arr

    # Reset refs_cnt as insertion pointers
    for i in range(n_taxa):
        refs_cnt[i] = 0

    cdef double p_anc_val, p_pres_val
    for i in range(n_refs):
        ref_t = ref_taxids[i]
        if ref_t in taxid_to_idx:
            idx = taxid_to_idx[ref_t]

            p_anc_val = ref_p_ancient[i]
            if p_anc_val < 1e-10:
                p_anc_val = 1e-10
            elif p_anc_val > 1.0 - 1e-10:
                p_anc_val = 1.0 - 1e-10

            log_psi_anc[idx] += log(p_anc_val)
            log_psi_mod[idx] += log(1.0 - p_anc_val)
            n_refs_tax[idx] += 1

            ref_p_present_data[ref_ptr[idx] + refs_cnt[idx]] = ref_p_present[i]
            refs_cnt[idx] += 1

            tad_anc[idx] += ref_tad_ancient[i]
            tad_mod[idx] += ref_tad_modern[i]

    # Build traversal arrays using Python (setup phase), then convert to C arrays
    postorder_list = []
    visited = set()
    stack = [(root_idx, False)]
    while stack:
        node, done = stack.pop()
        if done:
            postorder_list.append(node)
        else:
            stack.append((node, True))
            for k in range(children_ptr[node], children_ptr[node + 1]):
                c = children_indices[k]
                if c not in visited:
                    stack.append((c, False))
            visited.add(node)

    preorder_list = []
    stack = [root_idx]
    while stack:
        node = stack.pop()
        preorder_list.append(node)
        for k in range(children_ptr[node + 1] - 1, children_ptr[node] - 1, -1):
            stack.append(children_indices[k])

    postorder_arr = np.array(postorder_list, dtype=np.int32)
    preorder_arr = np.array(preorder_list, dtype=np.int32)
    cdef int32_t[::1] postorder = postorder_arr
    cdef int32_t[::1] preorder = preorder_arr

    cdef double rho_AA = params.epoch.rho_AA
    cdef double rho_MA = params.epoch.rho_MA
    cdef double pi_0 = params.epoch.pi_0
    cdef double log_pi_0 = log(pi_0) if pi_0 > 0 else LOG_EPSILON
    cdef double log_1_minus_pi_0 = log(1.0 - pi_0) if pi_0 < 1.0 else LOG_EPSILON

    # Beliefs storage
    log_b_ancient_arr = np.empty(n_taxa, dtype=np.float64)
    log_b_modern_arr = np.empty(n_taxa, dtype=np.float64)
    cdef double[::1] log_b_anc_view = log_b_ancient_arr
    cdef double[::1] log_b_mod_view = log_b_modern_arr

    # Output arrays
    p_ancient_out_arr = np.zeros(n_taxa, dtype=np.float64)
    p_present_out_arr = np.zeros(n_taxa, dtype=np.float64)
    tad_ancient_out_arr = np.zeros(n_taxa, dtype=np.float64)
    tad_modern_out_arr = np.zeros(n_taxa, dtype=np.float64)
    cdef double[::1] p_ancient_out = p_ancient_out_arr
    cdef double[::1] p_present_out = p_present_out_arr
    cdef double[::1] tad_ancient_out = tad_ancient_out_arr
    cdef double[::1] tad_modern_out = tad_modern_out_arr

    # Run BP passes using CSR nogil functions
    _bp_bottom_up_csr(
        n_taxa, root_idx,
        &postorder[0], &parent_idx_view[0],
        &log_psi_anc[0], &log_psi_mod[0],
        &log_msg_up_anc[0], &log_msg_up_mod[0],
        log_pi_0, log_1_minus_pi_0, rho_AA, rho_MA
    )

    _bp_compute_beliefs_csr(
        n_taxa, root_idx, &postorder[0],
        &log_psi_anc[0], &log_psi_mod[0],
        &log_msg_up_anc[0], &log_msg_up_mod[0],
        &log_b_anc_view[0], &log_b_mod_view[0],
        log_pi_0, log_1_minus_pi_0
    )

    _bp_top_down_csr(
        n_taxa, root_idx,
        &preorder[0], &children_indices[0], &children_ptr[0],
        &log_psi_anc[0], &log_psi_mod[0],
        &log_msg_up_anc[0], &log_msg_up_mod[0],
        &log_msg_down_anc[0], &log_msg_down_mod[0],
        &log_b_anc_view[0], &log_b_mod_view[0],
        log_pi_0, log_1_minus_pi_0, rho_AA, rho_MA
    )

    cdef double lambda_leak = params.noisy_or.lambda_leak
    cdef double q_fail = params.noisy_or.q_fail
    cdef double q_fail_ref = params.noisy_or.q_fail_ref

    _bp_compute_marginals_csr(
        n_taxa, root_idx,
        &postorder[0], &children_indices[0], &children_ptr[0],
        &taxids[0],
        &log_psi_anc[0], &log_psi_mod[0],
        &log_msg_up_anc[0], &log_msg_up_mod[0],
        &log_msg_down_anc[0], &log_msg_down_mod[0],
        &ref_p_present_data[0], &ref_ptr[0],
        &n_refs_tax[0], &tad_anc[0], &tad_mod[0],
        &p_ancient_out[0], &p_present_out[0],
        &tad_ancient_out[0], &tad_modern_out[0],
        log_pi_0, log_1_minus_pi_0,
        lambda_leak, q_fail, q_fail_ref
    )

    # Build results dict (only at the end)
    results = {}
    for i in range(n_taxa):
        taxid = taxids[i]
        results[taxid] = {
            'p_ancient': p_ancient_out[i],
            'p_present': p_present_out[i],
            'tad_ancient': tad_ancient_out[i],
            'tad_modern': tad_modern_out[i],
            'tad_total': tad_ancient_out[i] + tad_modern_out[i],
            'n_refs': int(n_refs_tax[i]),
        }

    return results


def run_belief_propagation(
    dict ref_data_by_taxid,
    dict taxonomy_tree,
    int32_t root_taxid,
    dict hyperparams=None,
):
    """
    Run belief propagation on taxonomy tree to compute taxon-level posteriors.

    Parameters
    ----------
    ref_data_by_taxid : dict
        taxid -> list of ref dicts, each with 'p_ancient', 'p_present', 'tad_ancient', 'tad_modern'
    taxonomy_tree : dict
        taxid -> {'parent': parent_taxid or None, 'children': [child_taxids]}
    root_taxid : int
        Root taxon ID
    hyperparams : dict, optional
        Custom hyperparameters

    Returns
    -------
    dict: taxid -> {'p_ancient', 'p_present', 'tad_ancient', 'tad_modern'}
    """
    cdef ProfilerHyperparams params
    init_default_hyperparams(&params)

    if hyperparams is not None:
        if 'epoch' in hyperparams:
            e = hyperparams['epoch']
            if 'pi_0' in e: params.epoch.pi_0 = e['pi_0']
            if 'rho_AA' in e: params.epoch.rho_AA = e['rho_AA']
            if 'rho_MA' in e: params.epoch.rho_MA = e['rho_MA']

        if 'noisy_or' in hyperparams:
            n = hyperparams['noisy_or']
            if 'lambda_leak' in n: params.noisy_or.lambda_leak = n['lambda_leak']
            if 'q_fail' in n: params.noisy_or.q_fail = n['q_fail']
            if 'q_fail_ref' in n: params.noisy_or.q_fail_ref = n['q_fail_ref']

    # Initialize beliefs for all taxa
    cdef dict beliefs = {}
    for taxid in taxonomy_tree:
        beliefs[taxid] = {
            'psi_ancient': 1.0,
            'psi_modern': 1.0,
            'n_refs': 0,
            'msg_up_ancient': 1.0,
            'msg_up_modern': 1.0,
            'msg_down_ancient': 1.0,
            'msg_down_modern': 1.0,
            'tad_ancient': 0.0,
            'tad_modern': 0.0,
            'p_present_refs': [],
            # For TAD-weighted authenticity aggregation
            'auth_weighted_sum': 0.0,
            'tad_sum': 0.0,
            # For TAD-weighted read exclusivity aggregation
            'exclusivity_weighted_sum': 0.0,
        }

    # Accumulate reference observations
    for taxid, refs in ref_data_by_taxid.items():
        if taxid not in beliefs:
            continue

        for ref in refs:
            p_anc = ref.get('p_ancient', 0.5)
            p_pres = ref.get('p_present', 0.0)

            # Clamp
            if p_anc < 1e-10:
                p_anc = 1e-10
            elif p_anc > 1.0 - 1e-10:
                p_anc = 1.0 - 1e-10

            beliefs[taxid]['psi_ancient'] *= p_anc
            beliefs[taxid]['psi_modern'] *= (1.0 - p_anc)
            beliefs[taxid]['n_refs'] += 1
            beliefs[taxid]['p_present_refs'].append(p_pres)

            beliefs[taxid]['tad_ancient'] += ref.get('tad_ancient', 0.0)
            beliefs[taxid]['tad_modern'] += ref.get('tad_modern', 0.0)

            # Accumulate TAD-weighted authenticity and exclusivity
            tad_val = ref.get('tad', 0.0)
            auth_val = ref.get('authenticity_score', 0.0)
            excl_val = ref.get('read_exclusivity', 1.0)
            beliefs[taxid]['auth_weighted_sum'] += auth_val * tad_val
            beliefs[taxid]['exclusivity_weighted_sum'] += excl_val * tad_val
            beliefs[taxid]['tad_sum'] += tad_val

    # Build traversal orders
    def get_postorder(tree, root):
        """Post-order traversal (children before parents)."""
        order = []
        stack = [(root, False)]
        while stack:
            node, visited = stack.pop()
            if visited:
                order.append(node)
            else:
                stack.append((node, True))
                for child in tree.get(node, {}).get('children', []):
                    stack.append((child, False))
        return order

    def get_preorder(tree, root):
        """Pre-order traversal (parents before children)."""
        order = []
        stack = [root]
        while stack:
            node = stack.pop()
            order.append(node)
            for child in reversed(tree.get(node, {}).get('children', [])):
                stack.append(child)
        return order

    postorder = get_postorder(taxonomy_tree, root_taxid)
    preorder = get_preorder(taxonomy_tree, root_taxid)

    cdef double rho_AA = params.epoch.rho_AA
    cdef double rho_MA = params.epoch.rho_MA
    cdef double pi_0 = params.epoch.pi_0

    # Bottom-up pass
    for taxid in postorder:
        b = beliefs[taxid]
        children = taxonomy_tree[taxid].get('children', [])
        is_root = (taxid == root_taxid)

        # Local belief = prior * observation * child messages
        phi_ancient = pi_0 if is_root else 1.0
        phi_modern = (1.0 - pi_0) if is_root else 1.0

        b_ancient = phi_ancient * b['psi_ancient'] * b['msg_up_ancient']
        b_modern = phi_modern * b['psi_modern'] * b['msg_up_modern']

        total = b_ancient + b_modern
        if total > 0:
            b_ancient /= total
            b_modern /= total
        else:
            b_ancient = 0.5
            b_modern = 0.5

        b['b_ancient'] = b_ancient
        b['b_modern'] = b_modern

        # Send upward message to parent
        parent_taxid = taxonomy_tree[taxid].get('parent')
        if parent_taxid is not None and parent_taxid in beliefs:
            msg_ancient = rho_AA * b_ancient + (1.0 - rho_AA) * b_modern
            msg_modern = rho_MA * b_ancient + (1.0 - rho_MA) * b_modern

            total = msg_ancient + msg_modern
            if total > 0:
                msg_ancient /= total
                msg_modern /= total

            beliefs[parent_taxid]['msg_up_ancient'] *= msg_ancient
            beliefs[parent_taxid]['msg_up_modern'] *= msg_modern

    # Top-down pass
    for taxid in preorder:
        b = beliefs[taxid]
        children = taxonomy_tree[taxid].get('children', [])
        is_root = (taxid == root_taxid)

        for child_taxid in children:
            if child_taxid not in beliefs:
                continue

            child_b = beliefs[child_taxid]

            # Get child's upward message
            child_b_anc = child_b.get('b_ancient', 0.5)
            child_b_mod = child_b.get('b_modern', 0.5)
            child_msg_anc = rho_AA * child_b_anc + (1.0 - rho_AA) * child_b_mod
            child_msg_mod = rho_MA * child_b_anc + (1.0 - rho_MA) * child_b_mod

            total = child_msg_anc + child_msg_mod
            if total > 0:
                child_msg_anc /= total
                child_msg_mod /= total

            # Parent context excluding this child
            phi_ancient = pi_0 if is_root else 1.0
            phi_modern = (1.0 - pi_0) if is_root else 1.0

            msg_up_exc_anc = b['msg_up_ancient']
            msg_up_exc_mod = b['msg_up_modern']

            if child_msg_anc > 1e-300:
                msg_up_exc_anc /= child_msg_anc
            if child_msg_mod > 1e-300:
                msg_up_exc_mod /= child_msg_mod

            b_exc_anc = phi_ancient * b['psi_ancient'] * b['msg_down_ancient'] * msg_up_exc_anc
            b_exc_mod = phi_modern * b['psi_modern'] * b['msg_down_modern'] * msg_up_exc_mod

            total = b_exc_anc + b_exc_mod
            if total > 0:
                b_exc_anc /= total
                b_exc_mod /= total
            else:
                b_exc_anc = 0.5
                b_exc_mod = 0.5

            # Downward message
            msg_down_anc = rho_AA * b_exc_anc + rho_MA * b_exc_mod
            msg_down_mod = (1.0 - rho_AA) * b_exc_anc + (1.0 - rho_MA) * b_exc_mod

            total = msg_down_anc + msg_down_mod
            if total > 0:
                msg_down_anc /= total
                msg_down_mod /= total

            child_b['msg_down_ancient'] = msg_down_anc
            child_b['msg_down_modern'] = msg_down_mod

    # Compute marginals and aggregate TAD
    cdef double lambda_leak = params.noisy_or.lambda_leak
    cdef double q_fail = params.noisy_or.q_fail
    cdef double q_fail_ref = params.noisy_or.q_fail_ref

    results = {}

    for taxid in postorder:
        b = beliefs[taxid]
        children = taxonomy_tree[taxid].get('children', [])
        is_root = (taxid == root_taxid)

        # Final marginal
        phi_ancient = pi_0 if is_root else 1.0
        phi_modern = (1.0 - pi_0) if is_root else 1.0

        unnorm_anc = phi_ancient * b['psi_ancient'] * b['msg_down_ancient'] * b['msg_up_ancient']
        unnorm_mod = phi_modern * b['psi_modern'] * b['msg_down_modern'] * b['msg_up_modern']

        total = unnorm_anc + unnorm_mod
        if total > 0:
            p_ancient = unnorm_anc / total
        else:
            p_ancient = 0.5

        # Noisy-OR presence from references
        p_present_refs = b.get('p_present_refs', [])
        if p_present_refs:
            p_absent = 1.0 - lambda_leak
            for p_ref in p_present_refs:
                p_absent *= (1.0 - p_ref * (1.0 - q_fail_ref))
            p_present_from_refs = 1.0 - p_absent
        else:
            p_present_from_refs = 0.0

        # Noisy-OR presence from children
        child_p_presents = []
        for child_taxid in children:
            if child_taxid in results:
                child_p_presents.append(results[child_taxid]['p_present'])

        if child_p_presents:
            p_absent = 1.0 - lambda_leak
            for p_child in child_p_presents:
                p_absent *= (1.0 - p_child * (1.0 - q_fail))
            p_present_from_children = 1.0 - p_absent
        else:
            p_present_from_children = 0.0

        # Combine: taxon is present if refs OR children indicate presence
        p_present = max(p_present_from_refs, p_present_from_children)

        # Aggregate TAD, authenticity, and exclusivity from children
        tad_ancient = b['tad_ancient']
        tad_modern = b['tad_modern']
        auth_weighted_sum = b['auth_weighted_sum']
        excl_weighted_sum = b['exclusivity_weighted_sum']
        tad_sum = b['tad_sum']

        for child_taxid in children:
            if child_taxid in results:
                tad_ancient += results[child_taxid]['tad_ancient']
                tad_modern += results[child_taxid]['tad_modern']
                auth_weighted_sum += results[child_taxid].get('auth_weighted_sum', 0.0)
                excl_weighted_sum += results[child_taxid].get('excl_weighted_sum', 0.0)
                tad_sum += results[child_taxid].get('tad_sum', 0.0)

        # Compute TAD-weighted averages
        authenticity_score = auth_weighted_sum / tad_sum if tad_sum > 0 else 0.0
        read_exclusivity = excl_weighted_sum / tad_sum if tad_sum > 0 else 1.0

        results[taxid] = {
            'p_ancient': p_ancient,
            'p_present': p_present,
            'tad_ancient': tad_ancient,
            'tad_modern': tad_modern,
            'tad_total': tad_ancient + tad_modern,
            'n_refs': b['n_refs'],
            'authenticity_score': authenticity_score,
            'read_exclusivity': read_exclusivity,
            'auth_weighted_sum': auth_weighted_sum,
            'excl_weighted_sum': excl_weighted_sum,
            'tad_sum': tad_sum,
        }

    return results


# =============================================================================
# GMRF (GAUSSIAN MARKOV RANDOM FIELD) IMPLEMENTATION
# =============================================================================
#
# Phase 1 of the unified Bayesian profiler: Graph-GMRF smoothing for
# reference-level ancientness estimates.
#
# The key insight: references sharing reads should have correlated ancientness
# values. We use a Gaussian Markov Random Field with the read-sharing graph
# Laplacian as the precision matrix.
#
# Model:
#   P(η) ∝ exp(-τ/2 · ηᵀ L η)
#        = exp(-τ/2 · Σ_{r~s} w_rs (η_r - η_s)²)
#
# where:
#   η_r = logit ancientness for reference r
#   L = D - W (graph Laplacian)
#   D = diagonal matrix of weighted degrees
#   W = adjacency matrix with w_rs = normalized shared reads
#   τ = smoothing precision
#
# Edge normalization (from design doc):
#   w_rs = shared_reads_rs / sqrt(n_reads_r × n_reads_s)
#
# This ensures the GMRF is stable across references with vastly different
# coverage depths.
# =============================================================================

from bam_filter.probabilistic_profiler cimport (
    GMRFHyperparams, GMRFState,
)


cdef void init_gmrf_hyperparams(GMRFHyperparams* params) noexcept nogil:
    """Initialize GMRF hyperparameters with defaults."""
    params.tau = 1.0           # Moderate smoothing
    params.mu_eta = 0.0        # Neutral prior (50-50)
    params.epsilon = 0.01      # Small regularization for isolated nodes
    params.max_iter = 100      # Maximum iterations
    params.tol = 1e-6          # Convergence tolerance


cdef GMRFState* create_gmrf_state(int32_t n_refs) noexcept nogil:
    """Allocate and initialize GMRF state for n_refs references."""
    cdef GMRFState* state = <GMRFState*>malloc(sizeof(GMRFState))
    if state == NULL:
        return NULL

    state.n_refs = n_refs
    state.converged = False

    state.eta = <double*>calloc(n_refs, sizeof(double))
    state.eta_data = <double*>calloc(n_refs, sizeof(double))
    state.gamma = <double*>calloc(n_refs, sizeof(double))
    state.degree = <double*>calloc(n_refs, sizeof(double))

    if (state.eta == NULL or state.eta_data == NULL or
            state.gamma == NULL or state.degree == NULL):
        destroy_gmrf_state(state)
        return NULL

    # Initialize gamma to 0.5 (neutral)
    cdef int32_t i
    for i in range(n_refs):
        state.gamma[i] = 0.5

    return state


cdef void destroy_gmrf_state(GMRFState* state) noexcept nogil:
    """Free GMRF state memory."""
    if state == NULL:
        return

    if state.eta != NULL:
        free(state.eta)
    if state.eta_data != NULL:
        free(state.eta_data)
    if state.gamma != NULL:
        free(state.gamma)
    if state.degree != NULL:
        free(state.degree)

    free(state)


cdef void init_eta_from_damage(
    GMRFState* state,
    double* log_bf,
    double prior_log_odds,
    int32_t n_refs
) noexcept nogil:
    """
    Initialize ancientness logits from damage model log Bayes factors.

    η_r = log_BF_r + prior_log_odds

    This gives us the data-driven initialization before smoothing.
    """
    cdef int32_t i
    cdef double eta_val

    for i in range(n_refs):
        eta_val = log_bf[i] + prior_log_odds

        # Clamp to reasonable range to avoid numerical issues
        if eta_val > 20.0:
            eta_val = 20.0
        elif eta_val < -20.0:
            eta_val = -20.0

        state.eta_data[i] = eta_val
        state.eta[i] = eta_val


cdef void build_gmrf_degrees(
    GMRFState* state,
    uint32_t** neighbors,
    uint32_t** weights,
    uint32_t* degrees,
    int32_t n_refs
) noexcept nogil:
    """
    Build weighted degree vector from graph structure.

    Uses normalized edge weights: w_rs / sqrt(n_r × n_s) summed over neighbors.
    The degree[r] = Σ_s w_rs (normalized).
    """
    cdef int32_t i, j
    cdef uint32_t neighbor_idx
    cdef double weight_sum

    for i in range(n_refs):
        weight_sum = 0.0
        for j in range(<int32_t>degrees[i]):
            # Sum raw edge weights for now
            # Normalization happens in the smoothing step
            weight_sum += <double>weights[i][j]
        state.degree[i] = weight_sum


cdef int run_gmrf_smoothing(
    GMRFState* state,
    uint32_t** neighbors,
    uint32_t** weights,
    uint32_t* degrees,
    uint32_t* n_reads,
    GMRFHyperparams* params,
    int num_threads
) noexcept nogil:
    """
    Run GMRF coordinate descent smoothing.

    For each reference r, we optimize:
        η_r = argmin_η [ data_fit + τ · Σ_s w_rs (η_r - η_s)² ]

    The optimal update is:
        η_r = (η_data_r + τ · Σ_s w_rs_norm · η_s) / (1 + τ · d_r_norm)

    where:
        w_rs_norm = shared_reads_rs / sqrt(n_reads_r × n_reads_s)
        d_r_norm = Σ_s w_rs_norm

    Returns 0 on success, -1 on error.
    """
    cdef int32_t iteration, i, j
    cdef int32_t n_refs = state.n_refs
    cdef double max_change, change
    cdef double eta_old, eta_new
    cdef double neighbor_sum, degree_norm, w_norm
    cdef uint32_t neighbor_idx
    cdef double sqrt_nr, sqrt_ns
    cdef double tau = params.tau
    cdef double epsilon = params.epsilon

    state.converged = False

    for iteration in range(params.max_iter):
        max_change = 0.0

        # Coordinate descent: update each η_r in sequence
        for i in range(n_refs):
            eta_old = state.eta[i]

            # Skip isolated nodes (no neighbors) - use data-driven value
            if degrees[i] == 0:
                state.eta[i] = state.eta_data[i]
                continue

            # Compute normalized neighbor sum and degree
            neighbor_sum = 0.0
            degree_norm = 0.0
            sqrt_nr = sqrt(<double>n_reads[i]) if n_reads[i] > 0 else 1.0

            for j in range(<int32_t>degrees[i]):
                neighbor_idx = neighbors[i][j]

                # Normalized edge weight: w_rs / sqrt(n_r × n_s)
                sqrt_ns = sqrt(<double>n_reads[neighbor_idx]) if n_reads[neighbor_idx] > 0 else 1.0
                w_norm = <double>weights[i][j] / (sqrt_nr * sqrt_ns + epsilon)

                neighbor_sum += w_norm * state.eta[neighbor_idx]
                degree_norm += w_norm

            # Optimal update: weighted average of data and neighbors
            # η_r = (η_data_r + τ · neighbor_weighted_sum) / (1 + τ · degree_norm)
            eta_new = (state.eta_data[i] + tau * neighbor_sum) / (1.0 + tau * degree_norm + epsilon)

            # Clamp to reasonable range
            if eta_new > 20.0:
                eta_new = 20.0
            elif eta_new < -20.0:
                eta_new = -20.0

            state.eta[i] = eta_new

            change = fabs(eta_new - eta_old)
            if change > max_change:
                max_change = change

        # Check convergence
        if max_change < params.tol:
            state.converged = True
            break

    return 0


cdef void compute_gamma_from_smoothed_eta(GMRFState* state) noexcept nogil:
    """Convert smoothed eta (logits) to gamma (probabilities) via sigmoid."""
    cdef int32_t i

    for i in range(state.n_refs):
        state.gamma[i] = sigmoid(state.eta[i])


# =============================================================================
# PYTHON INTERFACE FOR GMRF
# =============================================================================

def create_gmrf_hyperparams(
    double tau=1.0,
    double mu_eta=0.0,
    double epsilon=0.01,
    int max_iter=100,
    double tol=1e-6,
):
    """Create GMRF hyperparameters dictionary."""
    return {
        'tau': tau,
        'mu_eta': mu_eta,
        'epsilon': epsilon,
        'max_iter': max_iter,
        'tol': tol,
    }


def apply_gmrf_smoothing(
    double[::1] log_bf,
    double[::1] p_ancient_raw,
    long[::1] n_reads,
    list adjacency_list,
    list weight_list,
    dict gmrf_params=None,
    int num_threads=4,
):
    """
    Apply GMRF smoothing to reference-level ancientness estimates.

    Parameters
    ----------
    log_bf : 1D array (n_refs,)
        Log Bayes factors from damage model
    p_ancient_raw : 1D array (n_refs,)
        Raw P(ancient) from damage model (for prior log-odds extraction)
    n_reads : 1D array (n_refs,)
        Read counts per reference (for edge normalization)
    adjacency_list : list of arrays
        For each reference, array of neighbor indices
    weight_list : list of arrays
        For each reference, array of edge weights (shared read counts)
    gmrf_params : dict, optional
        GMRF hyperparameters (uses defaults if None)
    num_threads : int
        Number of threads (for future parallelization)

    Returns
    -------
    dict with:
        'eta_smoothed': smoothed ancientness logits
        'gamma_smoothed': smoothed P(ancient) values
        'converged': whether optimization converged
        'eta_data': data-driven initialization
    """
    import numpy as np

    cdef int32_t n_refs = log_bf.shape[0]
    cdef int32_t i, j
    cdef GMRFHyperparams params
    cdef GMRFState* state = NULL

    # Initialize hyperparameters
    init_gmrf_hyperparams(&params)
    if gmrf_params is not None:
        if 'tau' in gmrf_params: params.tau = gmrf_params['tau']
        if 'mu_eta' in gmrf_params: params.mu_eta = gmrf_params['mu_eta']
        if 'epsilon' in gmrf_params: params.epsilon = gmrf_params['epsilon']
        if 'max_iter' in gmrf_params: params.max_iter = gmrf_params['max_iter']
        if 'tol' in gmrf_params: params.tol = gmrf_params['tol']

    # Allocate state
    state = create_gmrf_state(n_refs)
    if state == NULL:
        raise MemoryError("Failed to allocate GMRF state")

    # Build C arrays for adjacency
    cdef uint32_t** neighbors_c = <uint32_t**>malloc(n_refs * sizeof(uint32_t*))
    cdef uint32_t** weights_c = <uint32_t**>malloc(n_refs * sizeof(uint32_t*))
    cdef uint32_t* degrees_c = <uint32_t*>malloc(n_refs * sizeof(uint32_t))
    cdef uint32_t* n_reads_c = <uint32_t*>malloc(n_refs * sizeof(uint32_t))

    cdef double mean_p_anc = 0.0
    cdef double prior_log_odds = 0.0

    if neighbors_c == NULL or weights_c == NULL or degrees_c == NULL or n_reads_c == NULL:
        if neighbors_c != NULL: free(neighbors_c)
        if weights_c != NULL: free(weights_c)
        if degrees_c != NULL: free(degrees_c)
        if n_reads_c != NULL: free(n_reads_c)
        destroy_gmrf_state(state)
        raise MemoryError("Failed to allocate adjacency arrays")

    try:
        # Copy adjacency to C arrays
        for i in range(n_refs):
            n_reads_c[i] = <uint32_t>n_reads[i]
            if adjacency_list[i] is not None and len(adjacency_list[i]) > 0:
                adj_arr = np.asarray(adjacency_list[i], dtype=np.uint32)
                wt_arr = np.asarray(weight_list[i], dtype=np.uint32)
                degrees_c[i] = <uint32_t>len(adj_arr)
                neighbors_c[i] = <uint32_t*>malloc(degrees_c[i] * sizeof(uint32_t))
                weights_c[i] = <uint32_t*>malloc(degrees_c[i] * sizeof(uint32_t))
                for j in range(<int32_t>degrees_c[i]):
                    neighbors_c[i][j] = adj_arr[j]
                    weights_c[i][j] = wt_arr[j]
            else:
                degrees_c[i] = 0
                neighbors_c[i] = NULL
                weights_c[i] = NULL

        # Compute prior log-odds from mean raw P(ancient)
        mean_p_anc = 0.0
        for i in range(n_refs):
            mean_p_anc += p_ancient_raw[i]
        mean_p_anc /= n_refs if n_refs > 0 else 1.0

        if mean_p_anc <= 0.001:
            prior_log_odds = -7.0
        elif mean_p_anc >= 0.999:
            prior_log_odds = 7.0
        else:
            prior_log_odds = log(mean_p_anc / (1.0 - mean_p_anc))

        # Initialize eta from damage model
        init_eta_from_damage(state, &log_bf[0], prior_log_odds, n_refs)

        # Build degrees
        build_gmrf_degrees(state, neighbors_c, weights_c, degrees_c, n_refs)

        # Run GMRF smoothing
        run_gmrf_smoothing(state, neighbors_c, weights_c, degrees_c,
                          n_reads_c, &params, num_threads)

        # Convert to gamma
        compute_gamma_from_smoothed_eta(state)

        # Copy results to numpy
        eta_smoothed = np.empty(n_refs, dtype=np.float64)
        gamma_smoothed = np.empty(n_refs, dtype=np.float64)
        eta_data = np.empty(n_refs, dtype=np.float64)

        for i in range(n_refs):
            eta_smoothed[i] = state.eta[i]
            gamma_smoothed[i] = state.gamma[i]
            eta_data[i] = state.eta_data[i]

        return {
            'eta_smoothed': eta_smoothed,
            'gamma_smoothed': gamma_smoothed,
            'converged': bool(state.converged),
            'eta_data': eta_data,
        }

    finally:
        # Free adjacency arrays
        for i in range(n_refs):
            if neighbors_c[i] != NULL:
                free(neighbors_c[i])
            if weights_c[i] != NULL:
                free(weights_c[i])
        free(neighbors_c)
        free(weights_c)
        free(degrees_c)
        free(n_reads_c)
        destroy_gmrf_state(state)


# =============================================================================
# CYTHON-LEVEL GMRF INTEGRATION WITH WEIGHTEDGRAPH
# =============================================================================
#
# This function integrates GMRF smoothing directly with the WeightedGraph
# structure, avoiding any Python/NumPy overhead or file I/O.
# =============================================================================

cdef int apply_gmrf_smoothing_to_graph(
    WeightedGraph* graph,
    double* gamma_values,
    double* damage_log_bf,
    uint32_t* n_reads,
    uint32_t n_refs,
    double tau,
    int max_iter,
    double tol,
    int verbose,
) noexcept nogil:
    """
    Apply GMRF smoothing directly using WeightedGraph adjacency structure.

    This is the pure-Cython integration point for GMRF smoothing, designed to be
    called from processor.pyx after graph analysis, without any file I/O.

    Parameters
    ----------
    graph : WeightedGraph*
        Graph containing adjacency lists (nodes[i].neighbors, weights, degree)
    gamma_values : double*
        P(ancient) per reference - modified in place with smoothed values
    damage_log_bf : double*
        Log Bayes factors from damage model (input only)
    n_reads : uint32_t*
        Read counts per reference (for edge normalization)
    n_refs : uint32_t
        Number of references
    tau : double
        GMRF smoothing precision (higher = more smoothing)
    max_iter : int
        Maximum coordinate descent iterations
    tol : double
        Convergence tolerance
    verbose : int
        Verbosity level

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef GMRFHyperparams params
    cdef GMRFState* state = NULL
    cdef int32_t i, j
    cdef double mean_gamma, prior_log_odds
    cdef GraphNode* node

    if graph == NULL or gamma_values == NULL or n_refs == 0:
        return -1

    # Initialize hyperparameters
    init_gmrf_hyperparams(&params)
    params.tau = tau
    params.max_iter = max_iter
    params.tol = tol

    # Allocate state
    state = create_gmrf_state(<int32_t>n_refs)
    if state == NULL:
        return -1

    # Compute prior log-odds from mean gamma (input P(ancient))
    mean_gamma = 0.0
    for i in range(<int32_t>n_refs):
        mean_gamma += gamma_values[i]
    mean_gamma /= <double>n_refs if n_refs > 0 else 1.0

    if mean_gamma <= 0.001:
        prior_log_odds = -7.0
    elif mean_gamma >= 0.999:
        prior_log_odds = 7.0
    else:
        prior_log_odds = log(mean_gamma / (1.0 - mean_gamma))

    # Initialize eta from damage log Bayes factors
    if damage_log_bf != NULL:
        init_eta_from_damage(state, damage_log_bf, prior_log_odds, <int32_t>n_refs)
    else:
        # Fall back to initializing from gamma values
        for i in range(<int32_t>n_refs):
            if gamma_values[i] <= 0.001:
                state.eta_data[i] = -7.0
            elif gamma_values[i] >= 0.999:
                state.eta_data[i] = 7.0
            else:
                state.eta_data[i] = log(gamma_values[i] / (1.0 - gamma_values[i]))
            state.eta[i] = state.eta_data[i]

    # Run GMRF smoothing using WeightedGraph directly
    # The run_gmrf_smoothing function expects uint32_t** but we have GraphNode*
    # We need to adapt by creating temporary pointer arrays pointing to graph data
    cdef uint32_t** neighbors_ptr = <uint32_t**>malloc(n_refs * sizeof(uint32_t*))
    cdef uint32_t** weights_ptr = <uint32_t**>malloc(n_refs * sizeof(uint32_t*))
    cdef uint32_t* degrees_ptr = <uint32_t*>malloc(n_refs * sizeof(uint32_t))

    if neighbors_ptr == NULL or weights_ptr == NULL or degrees_ptr == NULL:
        if neighbors_ptr != NULL: free(neighbors_ptr)
        if weights_ptr != NULL: free(weights_ptr)
        if degrees_ptr != NULL: free(degrees_ptr)
        destroy_gmrf_state(state)
        return -1

    # Point to existing graph data (no copy needed)
    for i in range(<int32_t>n_refs):
        if graph.nodes != NULL and i < <int32_t>graph.num_nodes:
            node = &graph.nodes[i]
            neighbors_ptr[i] = node.neighbors
            weights_ptr[i] = node.weights
            degrees_ptr[i] = node.degree
        else:
            neighbors_ptr[i] = NULL
            weights_ptr[i] = NULL
            degrees_ptr[i] = 0

    # Build degrees for normalization
    build_gmrf_degrees(state, neighbors_ptr, weights_ptr, degrees_ptr, <int32_t>n_refs)

    # Run GMRF smoothing
    run_gmrf_smoothing(state, neighbors_ptr, weights_ptr, degrees_ptr,
                       n_reads, &params, 1)

    # Convert smoothed eta to gamma
    compute_gamma_from_smoothed_eta(state)

    # Copy smoothed values back to gamma_values
    for i in range(<int32_t>n_refs):
        gamma_values[i] = state.gamma[i]

    # Cleanup (only free pointer arrays, not the data they point to)
    free(neighbors_ptr)
    free(weights_ptr)
    free(degrees_ptr)
    destroy_gmrf_state(state)

    return 0


# =============================================================================
# DIRECT REFSTATS -> POSTERIORS COMPUTATION (CYTHON-ONLY, NO PYTHON DICTS)
# =============================================================================
#
# This section provides functions to compute posteriors directly from RefStats
# structures, bypassing Python dict conversion entirely.
# =============================================================================

from bam_filter.stats cimport RefStats


cdef void populate_profile_from_refstats(
    RefStats* stats,
    RefProfileData* profile
) noexcept nogil:
    """
    Copy fields from RefStats to RefProfileData for posterior computation.

    This function bridges the stats module and the probabilistic profiler,
    avoiding any Python dict conversion.
    """
    cdef int j

    # Copy damage counts (n_5p, k_5p, n_3p, k_3p are double[20] arrays)
    for j in range(20):
        profile.damage.n_5p[j] = stats.n_5p[j]
        profile.damage.k_5p[j] = stats.k_5p[j]
        profile.damage.n_3p[j] = stats.n_3p[j]
        profile.damage.k_3p[j] = stats.k_3p[j]

    # Copy coverage data
    profile.coverage.breadth = stats.breadth
    profile.coverage.mean_depth = stats.mean_coverage
    profile.coverage.wcb = stats.weighted_contiguity_breadth
    profile.coverage.norm_entropy = stats.norm_spatial_entropy
    profile.coverage.norm_gini = stats.norm_gini
    profile.coverage.authenticity = stats.norm_spatial_entropy - stats.norm_gini
    profile.coverage.n_reads = stats.n_reads
    profile.coverage.tad = <double>stats.tax_abund_tad


cdef void compute_posteriors_from_refstats_array(
    RefStats* stats_array,
    int32_t n_refs,
    double* p_ancient_out,
    double* p_present_out,
    double* log_bf_out,
    double* tad_ancient_out,
    double* tad_modern_out,
    ProfilerHyperparams* params,
    int num_threads
) noexcept nogil:
    """
    Compute posteriors for all references from RefStats array.

    This is the core Cython function that computes damage-based p_ancient
    and coverage-based p_present directly from RefStats, with no Python overhead.

    Parameters
    ----------
    stats_array : RefStats*
        Array of per-reference statistics (from compute_bam_stats)
    n_refs : int32_t
        Number of references
    p_ancient_out, p_present_out, log_bf_out : double*
        Output arrays for posteriors (must be pre-allocated)
    tad_ancient_out, tad_modern_out : double*
        Output arrays for TAD splits (must be pre-allocated)
    params : ProfilerHyperparams*
        Hyperparameters for damage and presence models
    num_threads : int
        Number of OpenMP threads
    """
    cdef int32_t i
    cdef RefProfileData profile

    for i in prange(n_refs, nogil=True, num_threads=num_threads, schedule='static'):
        # Populate RefProfileData from RefStats
        populate_profile_from_refstats(&stats_array[i], &profile)

        # Compute posteriors
        compute_ref_posteriors(&profile, params)

        # Store results
        p_ancient_out[i] = profile.posteriors.p_ancient
        p_present_out[i] = profile.posteriors.p_present
        log_bf_out[i] = profile.posteriors.log_bf
        tad_ancient_out[i] = profile.posteriors.tad_ancient
        tad_modern_out[i] = profile.posteriors.tad_modern


def compute_posteriors_from_stats_array(
    stats_capsule,
    int n_refs,
    int num_threads=4,
    dict hyperparams=None,
):
    """
    Compute posteriors from RefStats array passed via PyCapsule.

    This is the Python entry point for computing posteriors directly from
    Cython-level RefStats, avoiding dict conversion overhead.

    Parameters
    ----------
    stats_capsule : PyCapsule
        Capsule containing RefStats* pointer
    n_refs : int
        Number of references
    num_threads : int
        Number of threads for parallel computation
    hyperparams : dict, optional
        Override default hyperparameters

    Returns
    -------
    dict with numpy arrays:
        'p_ancient', 'p_present', 'log_bf', 'tad_ancient', 'tad_modern'
    """
    import numpy as np

    # Keep reference to capsule to avoid dangling pointer
    cdef object capsule_ref = stats_capsule
    cdef RefStats* stats_ptr = <RefStats*>PyCapsule_GetPointer(capsule_ref, "RefStats")
    if stats_ptr == NULL:
        raise ValueError("Invalid stats capsule")

    cdef ProfilerHyperparams params
    init_default_hyperparams(&params)

    if hyperparams is not None:
        if 'damage' in hyperparams:
            d = hyperparams['damage']
            if 'p_0' in d: params.damage.p_0 = d['p_0']
            if 'p_bg' in d: params.damage.p_bg = d['p_bg']
            if 'tau_A' in d: params.damage.tau_A = d['tau_A']
            if 'kappa_A' in d: params.damage.kappa_A = d['kappa_A']
            if 'p_err' in d: params.damage.p_err = d['p_err']
            if 'kappa_M' in d: params.damage.kappa_M = d['kappa_M']
            if 'pi_A_ref' in d: params.damage.pi_A_ref = d['pi_A_ref']
            precompute_damage_params(&params.damage)
        if 'presence' in hyperparams:
            p = hyperparams['presence']
            if 'intercept' in p: params.presence.intercept = p['intercept']
            if 'scale' in p: params.presence.scale = p['scale']

    # Allocate output arrays
    p_ancient_arr = np.empty(n_refs, dtype=np.float64)
    p_present_arr = np.empty(n_refs, dtype=np.float64)
    log_bf_arr = np.empty(n_refs, dtype=np.float64)
    tad_ancient_arr = np.empty(n_refs, dtype=np.float64)
    tad_modern_arr = np.empty(n_refs, dtype=np.float64)

    cdef double[::1] p_ancient_view = p_ancient_arr
    cdef double[::1] p_present_view = p_present_arr
    cdef double[::1] log_bf_view = log_bf_arr
    cdef double[::1] tad_ancient_view = tad_ancient_arr
    cdef double[::1] tad_modern_view = tad_modern_arr

    with nogil:
        compute_posteriors_from_refstats_array(
            stats_ptr, <int32_t>n_refs,
            &p_ancient_view[0], &p_present_view[0], &log_bf_view[0],
            &tad_ancient_view[0], &tad_modern_view[0],
            &params, num_threads
        )

    return {
        'p_ancient': p_ancient_arr,
        'p_present': p_present_arr,
        'log_bf': log_bf_arr,
        'tad_ancient': tad_ancient_arr,
        'tad_modern': tad_modern_arr,
    }


# =============================================================================
# FULL CYTHON PIPELINE: STATS -> POSTERIORS -> TAXONOMY -> BELIEF PROPAGATION
# =============================================================================
#
# This section provides a single entry point that runs the entire profiler
# pipeline in Cython with minimal Python overhead.
# =============================================================================

from bam_filter.taxonomy_db cimport (
    TaxonomyDB, TaxNode, AccessionMap, TaxonomyDatabase, AccessionMapping,
    lookup_taxid_duckdb, get_rank_name,
)


def run_full_profile_from_capsule(
    stats_capsule,
    list ref_names,
    TaxonomyDatabase taxonomy_db,
    AccessionMapping acc_map,
    int num_threads=4,
    dict hyperparams=None,
    bint verbose=True,
):
    """
    Run full profiler pipeline from RefStats capsule to final profile.

    This is the main Cython entry point that does everything without going
    back to Python for intermediate steps.

    Parameters
    ----------
    stats_capsule : PyCapsule
        Capsule containing RefStats* array from compute_bam_stats
    ref_names : list of str
        Reference names (parallel to stats array)
    taxonomy_db : TaxonomyDatabase
        Cython taxonomy database object
    acc_map : AccessionMapping
        Cython accession-to-taxid mapping object
    num_threads : int
        Number of threads for parallel computation
    hyperparams : dict, optional
        Hyperparameter overrides
    verbose : bool
        Print progress messages

    Returns
    -------
    dict : taxid -> profile data with p_ancient, p_present, tad_ancient, tad_modern, etc.
    """
    import numpy as np

    cdef int32_t n_refs = len(ref_names)
    if n_refs == 0:
        return {}

    # Keep reference to capsule to avoid dangling pointer
    cdef object capsule_ref = stats_capsule
    cdef RefStats* stats_ptr = <RefStats*>PyCapsule_GetPointer(capsule_ref, "RefStats")
    if stats_ptr == NULL:
        raise ValueError("Invalid stats capsule")

    # Initialize hyperparameters
    cdef ProfilerHyperparams params
    init_default_hyperparams(&params)

    if hyperparams is not None:
        if 'damage' in hyperparams:
            d = hyperparams['damage']
            if 'p_0' in d: params.damage.p_0 = d['p_0']
            if 'p_bg' in d: params.damage.p_bg = d['p_bg']
            if 'tau_A' in d: params.damage.tau_A = d['tau_A']
            if 'kappa_A' in d: params.damage.kappa_A = d['kappa_A']
            if 'p_err' in d: params.damage.p_err = d['p_err']
            if 'kappa_M' in d: params.damage.kappa_M = d['kappa_M']
            if 'pi_A_ref' in d: params.damage.pi_A_ref = d['pi_A_ref']
            precompute_damage_params(&params.damage)
        if 'presence' in hyperparams:
            p = hyperparams['presence']
            if 'intercept' in p: params.presence.intercept = p['intercept']
            if 'scale' in p: params.presence.scale = p['scale']
        if 'epoch' in hyperparams:
            e = hyperparams['epoch']
            if 'pi_0' in e: params.epoch.pi_0 = e['pi_0']
            if 'rho_AA' in e: params.epoch.rho_AA = e['rho_AA']
            if 'rho_MA' in e: params.epoch.rho_MA = e['rho_MA']
        if 'noisy_or' in hyperparams:
            n = hyperparams['noisy_or']
            if 'lambda_leak' in n: params.noisy_or.lambda_leak = n['lambda_leak']
            if 'q_fail' in n: params.noisy_or.q_fail = n['q_fail']
            if 'q_fail_ref' in n: params.noisy_or.q_fail_ref = n['q_fail_ref']

    # Allocate posterior arrays
    cdef double* p_ancient = <double*>malloc(n_refs * sizeof(double))
    cdef double* p_present = <double*>malloc(n_refs * sizeof(double))
    cdef double* log_bf = <double*>malloc(n_refs * sizeof(double))
    cdef double* tad_ancient = <double*>malloc(n_refs * sizeof(double))
    cdef double* tad_modern = <double*>malloc(n_refs * sizeof(double))
    cdef int32_t* ref_taxids = <int32_t*>malloc(n_refs * sizeof(int32_t))

    if (p_ancient == NULL or p_present == NULL or log_bf == NULL or
        tad_ancient == NULL or tad_modern == NULL or ref_taxids == NULL):
        free(p_ancient); free(p_present); free(log_bf)
        free(tad_ancient); free(tad_modern); free(ref_taxids)
        raise MemoryError("Failed to allocate posterior arrays")

    cdef int32_t i, idx
    cdef bytes ref_name_bytes
    cdef const char* ref_name_cstr
    cdef int32_t taxid
    cdef int32_t valid_count = 0
    cdef TaxNode* node
    cdef int32_t rank_id, name_len
    cdef char* name_ptr

    # Step 1: Compute posteriors for all references
    if verbose:
        print(f"  Computing posteriors for {n_refs} references...")

    with nogil:
        compute_posteriors_from_refstats_array(
            stats_ptr, n_refs,
            p_ancient, p_present, log_bf,
            tad_ancient, tad_modern,
            &params, num_threads
        )

    # Step 2: Map reference names to taxids using AccessionMapping
    if verbose:
        print(f"  Mapping references to taxonomy...")

    for i in range(n_refs):
        ref_name_bytes = ref_names[i].encode('utf-8')
        ref_name_cstr = ref_name_bytes
        taxid = acc_map._get_taxid_nogil(ref_name_cstr)
        ref_taxids[i] = taxid
        if taxid > 0:
            valid_count += 1

    if verbose:
        print(f"  Mapped {valid_count}/{n_refs} references to valid taxids")

    if valid_count == 0:
        free(p_ancient); free(p_present); free(log_bf)
        free(tad_ancient); free(tad_modern); free(ref_taxids)
        return {}

    # Step 3: Build taxonomy tree structure for belief propagation
    # Collect all unique taxids and their parents
    cdef set taxid_set = set()
    cdef int32_t parent_taxid
    cdef TaxonomyDB* db = taxonomy_db.db

    for i in range(n_refs):
        taxid = ref_taxids[i]
        if taxid > 0:
            # Walk up to root, collecting all taxids in lineage
            while taxid > 0 and taxid <= db.max_taxid:
                taxid_set.add(taxid)
                idx = db.taxid_to_idx[taxid]
                if idx < 0 or idx >= db.n_nodes:
                    break
                parent_taxid = db.nodes[idx].parent_taxid
                if parent_taxid == taxid:
                    break  # Root reached
                taxid = parent_taxid

    if len(taxid_set) == 0:
        free(p_ancient); free(p_present); free(log_bf)
        free(tad_ancient); free(tad_modern); free(ref_taxids)
        return {}

    # Convert to numpy arrays for run_belief_propagation_native
    taxid_list = sorted(taxid_set)
    cdef int32_t n_taxa = len(taxid_list)

    taxids_arr = np.array(taxid_list, dtype=np.int32)
    parent_taxids_arr = np.zeros(n_taxa, dtype=np.int32)
    cdef int32_t[::1] taxids_view = taxids_arr
    cdef int32_t[::1] parent_taxids_view = parent_taxids_arr

    # Build taxid -> index map
    taxid_to_idx_map = {t: i for i, t in enumerate(taxid_list)}

    # Find root (taxid where parent == self)
    cdef int32_t root_taxid = 1
    for i in range(n_taxa):
        taxid = taxids_view[i]
        if taxid > 0 and taxid <= db.max_taxid:
            idx = db.taxid_to_idx[taxid]
            if idx >= 0 and idx < db.n_nodes:
                parent_taxid = db.nodes[idx].parent_taxid
                parent_taxids_view[i] = parent_taxid
                if parent_taxid == taxid:
                    root_taxid = taxid

    # Build ref_taxids numpy array (only valid ones)
    ref_taxids_valid_list = []
    ref_p_ancient_list = []
    ref_p_present_list = []
    ref_tad_ancient_list = []
    ref_tad_modern_list = []
    ref_names_mapped = []
    ref_coverage_list = []
    ref_authenticity_list = []
    ref_exclusivity_list = []
    ref_reads_list = []

    for i in range(n_refs):
        taxid = ref_taxids[i]
        if taxid > 0 and taxid in taxid_to_idx_map:
            ref_taxids_valid_list.append(taxid)
            ref_p_ancient_list.append(p_ancient[i])
            ref_p_present_list.append(p_present[i])
            ref_tad_ancient_list.append(tad_ancient[i])
            ref_tad_modern_list.append(tad_modern[i])
            ref_names_mapped.append(ref_names[i])
            ref_coverage_list.append(stats_ptr[i].breadth)
            ref_authenticity_list.append(stats_ptr[i].authenticity_score)
            ref_exclusivity_list.append(1.0)  # Will be computed if graph data available
            ref_reads_list.append(stats_ptr[i].n_reads)

    ref_taxids_arr = np.array(ref_taxids_valid_list, dtype=np.int32)
    ref_p_ancient_arr = np.array(ref_p_ancient_list, dtype=np.float64)
    ref_p_present_arr = np.array(ref_p_present_list, dtype=np.float64)
    ref_tad_ancient_arr = np.array(ref_tad_ancient_list, dtype=np.float64)
    ref_tad_modern_arr = np.array(ref_tad_modern_list, dtype=np.float64)

    # Free C arrays (data now in numpy)
    free(p_ancient); free(p_present); free(log_bf)
    free(tad_ancient); free(tad_modern); free(ref_taxids)

    if verbose:
        print(f"  Running belief propagation on {n_taxa} taxa...")

    # Step 4: Run belief propagation
    profile = run_belief_propagation_native(
        taxids_arr,
        parent_taxids_arr,
        ref_taxids_arr,
        ref_p_ancient_arr,
        ref_p_present_arr,
        ref_tad_ancient_arr,
        ref_tad_modern_arr,
        root_taxid,
        hyperparams,
        num_threads,
    )

    # Step 5: Add taxonomy metadata to profile
    for taxid_key in profile:
        taxid = int(taxid_key)
        if taxid > 0 and taxid <= db.max_taxid:
            idx = db.taxid_to_idx[taxid]
            if idx >= 0 and idx < db.n_nodes:
                node = &db.nodes[idx]
                rank_id = node.rank_id
                name_ptr = db.names_buffer + node.name_offset
                name_len = node.name_length
                profile[taxid]['rank'] = get_rank_name(rank_id)
                profile[taxid]['name'] = name_ptr[:name_len].decode('utf-8', errors='replace')
                profile[taxid]['depth'] = node.depth

    if verbose:
        print(f"  Profile complete: {len(profile)} taxa")

    return profile
