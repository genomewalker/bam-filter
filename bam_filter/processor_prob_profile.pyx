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

"""
Probabilistic Taxonomic Profiler - Bayesian model for ancient DNA.

Computes TWO distinct probability outputs per taxon:
1. P(present|taxon) - Detection: Is this taxon truly in the sample?
2. P(ancient|taxon) - Ancient: Is this taxon ancient DNA?

Pipeline: REASSIGN → FILTER → PROFILE
- Reads BAM from filter stage (with ZP tags containing EM posteriors φ)
- Computes ALL statistics de novo from alignments
- Uses φ-weighting for multi-mapped reads (not LCA hard assignment)

Three-pass algorithm:
1. Pass 1 (Bottom-up): Accumulate coverage, damage, and metrics from alignments
2. Pass 2 (Compute): Calculate detection metrics and fit damage model per taxid
3. Pass 3 (Top-down): Apply hierarchical shrinkage for P(ancient)
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t, uint8_t
from libc.stdlib cimport malloc, calloc, free, realloc
from libc.string cimport memset, strcmp, strcpy, strlen, memcpy
from libc.math cimport sqrt, log, exp, isnan, isinf, fabs, pow, log10

from bam_filter.stats cimport RefStats
from bam_filter.taxonomy_db cimport TaxonomyDB, TaxonomyDatabase, AccessionMap

from bam_filter import logging as bf_logging
from typing import TypedDict, Optional, Dict, List
import math

LOG_TAG = "PROB_PROFILE"

# =============================================================================
# Constants
# =============================================================================

DEF MAX_DAMAGE_POS = 20
DEF INTERIOR_START = 10  # Positions 10-20 used for baseline estimation


# =============================================================================
# Beta Distribution Functions
# =============================================================================

cdef inline double beta_mean(double alpha, double beta) noexcept nogil:
    """Compute mean of Beta(alpha, beta) distribution."""
    cdef double total = alpha + beta
    if total <= 0.0:
        return 0.5
    return alpha / total

cdef inline double beta_variance(double alpha, double beta) noexcept nogil:
    """Compute variance of Beta(alpha, beta) distribution."""
    cdef double total = alpha + beta
    if total <= 0.0:
        return 0.0
    return (alpha * beta) / (total * total * (total + 1.0))

cdef double _betacf(double a, double b, double x) noexcept nogil:
    """Continued fraction for incomplete beta function (Lentz's algorithm)."""
    cdef int MAXIT = 200
    cdef double EPS = 3.0e-12
    cdef double FPMIN = 1.0e-30

    cdef double qab = a + b
    cdef double qap = a + 1.0
    cdef double qam = a - 1.0
    cdef double c = 1.0
    cdef double d = 1.0 - qab * x / qap
    cdef double aa, del_val
    cdef int m, m2

    if fabs(d) < FPMIN:
        d = FPMIN
    d = 1.0 / d
    cdef double h = d

    for m in range(1, MAXIT + 1):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d = 1.0 + aa * d
        if fabs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if fabs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        h *= d * c

        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d = 1.0 + aa * d
        if fabs(d) < FPMIN:
            d = FPMIN
        c = 1.0 + aa / c
        if fabs(c) < FPMIN:
            c = FPMIN
        d = 1.0 / d
        del_val = d * c
        h *= del_val
        if fabs(del_val - 1.0) < EPS:
            break

    return h

cdef double _lngamma(double x) noexcept nogil:
    """Log-gamma function (Lanczos approximation)."""
    cdef double[6] cof
    cof[0] = 76.18009172947146
    cof[1] = -86.50532032941677
    cof[2] = 24.01409824083091
    cof[3] = -1.231739572450155
    cof[4] = 0.1208650973866179e-2
    cof[5] = -0.5395239384953e-5

    cdef double y = x
    cdef double tmp = x + 5.5
    tmp -= (x + 0.5) * log(tmp)
    cdef double ser = 1.000000000190015
    cdef int j

    for j in range(6):
        y += 1.0
        ser += cof[j] / y

    return -tmp + log(2.5066282746310005 * ser / x)

cdef double log_beta_function(double alpha, double beta) noexcept nogil:
    """Compute log(B(alpha, beta)) = log(Gamma(a)) + log(Gamma(b)) - log(Gamma(a+b))."""
    return _lngamma(alpha) + _lngamma(beta) - _lngamma(alpha + beta)

cdef double _betai(double a, double b, double x) noexcept nogil:
    """Incomplete beta function I_x(a, b)."""
    if x < 0.0 or x > 1.0:
        return 0.0
    if x == 0.0 or x == 1.0:
        return x

    cdef double bt = exp(_lngamma(a + b) - _lngamma(a) - _lngamma(b) +
                         a * log(x) + b * log(1.0 - x))

    if x < (a + 1.0) / (a + b + 2.0):
        return bt * _betacf(a, b, x) / a
    else:
        return 1.0 - bt * _betacf(b, a, 1.0 - x) / b

cdef double beta_quantile(double alpha, double beta, double p) noexcept nogil:
    """Compute quantile of Beta(alpha, beta) distribution using bisection."""
    if p <= 0.0:
        return 0.0
    if p >= 1.0:
        return 1.0
    if alpha <= 0.0 or beta <= 0.0:
        return 0.5

    cdef double lo = 0.0
    cdef double hi = 1.0
    cdef double mid, cdf
    cdef double tol = 1e-10
    cdef int max_iter = 100
    cdef int i

    for i in range(max_iter):
        mid = (lo + hi) / 2.0
        cdf = _betai(alpha, beta, mid)
        if fabs(cdf - p) < tol:
            return mid
        if cdf < p:
            lo = mid
        else:
            hi = mid

    return mid

cdef double beta_ci_lower_approx(double alpha, double beta) noexcept nogil:
    """Approximate 2.5% quantile using logit-normal approximation."""
    cdef double m = alpha / (alpha + beta)
    cdef double v = (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0))
    cdef double logit_m, sigma2, sigma, logit_low, p_low

    if m <= 0.0 or m >= 1.0 or v <= 0.0:
        return 0.0

    logit_m = log(m / (1.0 - m))
    sigma2 = v / (m * m * (1.0 - m) * (1.0 - m))
    if sigma2 <= 0.0:
        return m
    sigma = sqrt(sigma2)
    logit_low = logit_m - 1.96 * sigma
    p_low = 1.0 / (1.0 + exp(-logit_low))
    return p_low

cdef double beta_ci_upper_approx(double alpha, double beta) noexcept nogil:
    """Approximate 97.5% quantile using logit-normal approximation."""
    cdef double m = alpha / (alpha + beta)
    cdef double v = (alpha * beta) / ((alpha + beta) * (alpha + beta) * (alpha + beta + 1.0))
    cdef double logit_m, sigma2, sigma, logit_high, p_high

    if m <= 0.0 or m >= 1.0 or v <= 0.0:
        return 1.0

    logit_m = log(m / (1.0 - m))
    sigma2 = v / (m * m * (1.0 - m) * (1.0 - m))
    if sigma2 <= 0.0:
        return m
    sigma = sqrt(sigma2)
    logit_high = logit_m + 1.96 * sigma
    p_high = 1.0 / (1.0 + exp(-logit_high))
    return p_high


# =============================================================================
# Initialization Functions
# =============================================================================

cdef void init_prob_profile_params(ProbProfileParams* params) noexcept nogil:
    """Initialize hyperparameters with sensible defaults."""
    # P(ancient) Beta prior
    params.alpha0 = 1.0
    params.beta0 = 1.0
    params.lambda_shrink = 10.0

    # TAD normalization
    params.kappa_dirichlet = 2.0

    # Credible intervals
    params.ci_lower_quantile = 0.025
    params.ci_upper_quantile = 0.975

    # Detection model weights for P(present)
    params.w_wcb = 1.0        # High WCB = good (contiguous coverage)
    params.w_cpc = 1.0        # High CPC = good (complex regions covered)
    params.w_ori = -1.0       # High ORI = bad (excessive overlap)
    params.w_entropy = 1.0    # High entropy = good (even spread)
    params.w_gini = -1.0      # High gini = bad (uneven coverage)
    params.w_ani = 0.5        # High ANI = good (sequence identity)

    # Minimum thresholds
    params.min_reads = 1.0
    params.min_p_present = 0.1

cdef void init_taxid_damage_accum(TaxidDamageAccum* accum) noexcept nogil:
    """Initialize damage accumulator to zeros."""
    cdef int i
    for i in range(MAX_DAMAGE_POS):
        accum.n_5p[i] = 0.0
        accum.k_5p[i] = 0.0
        accum.n_3p[i] = 0.0
        accum.k_3p[i] = 0.0
    accum.n_interior = 0.0
    accum.k_interior = 0.0
    accum.total_weight = 0.0
    accum.n_alignments = 0

cdef void init_taxid_coverage_accum(TaxidCoverageAccum* accum) noexcept nogil:
    """Initialize coverage accumulator to zeros."""
    accum.total_ref_length = 0
    accum.n_refs = 0
    accum.bases_covered = 0.0
    accum.total_coverage = 0.0
    accum.total_coverage_sq = 0.0
    accum.tad_sum = 0.0
    accum.tad_sum_ancient = 0.0
    accum.tad_sum_modern = 0.0
    accum.sum_interval_length_sq = 0.0
    accum.n_intervals = 0
    accum.dust_sum = 0.0
    accum.dust_weight = 0.0
    accum.spatial_entropy_sum = 0.0
    accum.gini_sum = 0.0
    accum.spatial_weight = 0.0
    accum.ani_sum = 0.0
    accum.ani_weight = 0.0
    accum.mapq_sum = 0.0
    accum.total_aligned_bases = 0.0

cdef void init_taxid_profile_stats(TaxidProfileStats* stats, int32_t taxid) noexcept nogil:
    """Initialize a TaxidProfileStats structure."""
    stats.taxid = taxid
    stats.parent_taxid = -1
    stats.rank_id = 0
    stats.n_refs = 0

    # Coverage & abundance
    stats.total_ref_length = 0
    stats.n_reads = 0.0
    stats.n_alns = 0.0
    stats.tad_total = 0.0
    stats.tad_ancient = 0.0
    stats.tad_modern = 0.0
    stats.breadth = 0.0
    stats.mean_coverage = 0.0
    stats.mean_coverage_covered = 0.0

    # Detection metrics
    stats.wcb = 0.0
    stats.cpc = 0.0
    stats.ori = 0.0
    stats.norm_spatial_entropy = 0.0
    stats.norm_gini = 0.0
    stats.ani_mean = 0.0
    stats.ani_corrected = 0.0
    stats.mapq_mean = 0.0
    stats.dust_mean = 0.0

    # Damage
    init_taxid_damage_accum(&stats.damage)
    stats.damage_amplitude = 0.0
    stats.damage_baseline = 0.0
    stats.damage_decay = 0.0

    # Probability outputs
    stats.p_present = 0.0
    stats.p_present_ci_lower = 0.0
    stats.p_present_ci_upper = 1.0
    stats.p_ancient = 0.5
    stats.p_ancient_ci_lower = 0.0
    stats.p_ancient_ci_upper = 1.0
    stats.beta_alpha_post = 1.0
    stats.beta_beta_post = 1.0
    stats.damage_log_bf = 0.0
    stats.shrinkage_weight = 1.0
    stats.detection_score = 0.0


# =============================================================================
# Accumulation Functions (Pass 1: Bottom-up)
# =============================================================================

cdef void accumulate_alignment_to_taxid(
    TaxidProfileStats* stats,
    double phi,
    double gamma,
    int64_t ref_length,
    int64_t aln_start,
    int64_t aln_end,
    double ani,
    double mapq,
    double dust,
    const double* k_5p,
    const double* n_5p,
    const double* k_3p,
    const double* n_3p
) noexcept nogil:
    """Accumulate an alignment's contribution to a taxid (φ-weighted)."""
    cdef int i
    cdef int64_t aln_len = aln_end - aln_start

    # Basic counts
    stats.n_alns += phi
    stats.n_reads += phi

    # TAD channel split
    stats.tad_ancient += phi * gamma
    stats.tad_modern += phi * (1.0 - gamma)
    stats.tad_total += phi

    # Quality metrics (φ-weighted accumulation)
    if ani > 0:
        stats.ani_mean += phi * ani
    if mapq >= 0:
        stats.mapq_mean += phi * mapq
    if dust >= 0:
        stats.dust_mean += phi * dust

    # Damage accumulation (φ-weighted)
    if k_5p != NULL and n_5p != NULL:
        for i in range(MAX_DAMAGE_POS):
            stats.damage.k_5p[i] += phi * k_5p[i]
            stats.damage.n_5p[i] += phi * n_5p[i]

    if k_3p != NULL and n_3p != NULL:
        for i in range(MAX_DAMAGE_POS):
            stats.damage.k_3p[i] += phi * k_3p[i]
            stats.damage.n_3p[i] += phi * n_3p[i]

    # Interior baseline (positions 10-19)
    if n_5p != NULL and k_5p != NULL:
        for i in range(INTERIOR_START, MAX_DAMAGE_POS):
            stats.damage.n_interior += phi * n_5p[i]
            stats.damage.k_interior += phi * k_5p[i]

    stats.damage.total_weight += phi
    stats.damage.n_alignments += 1

cdef void accumulate_coverage_interval(
    TaxidProfileStats* stats,
    double phi,
    int64_t interval_start,
    int64_t interval_end,
    int32_t depth
) noexcept nogil:
    """Accumulate a coverage interval's contribution (φ-weighted)."""
    cdef int64_t interval_len = interval_end - interval_start
    cdef double weighted_len = phi * <double>interval_len
    cdef double weighted_depth = phi * <double>depth * <double>interval_len

    # For WCB: sum of squared interval lengths
    stats.wcb += phi * <double>(interval_len * interval_len)

    # Coverage totals
    stats.breadth += weighted_len
    stats.mean_coverage += weighted_depth


# =============================================================================
# Computation Functions (Pass 2: Compute metrics)
# =============================================================================

cdef void compute_detection_metrics(TaxidProfileStats* stats) noexcept nogil:
    """Compute detection metrics from accumulated values."""
    cdef double ref_len_sq
    cdef double n_weight = stats.n_alns

    if stats.total_ref_length <= 0 or n_weight <= 0:
        return

    # Normalize accumulated sums
    stats.ani_mean /= n_weight
    stats.mapq_mean /= n_weight
    stats.dust_mean /= n_weight

    # Breadth: fraction of genome covered
    stats.breadth /= <double>stats.total_ref_length

    # Mean coverage
    if stats.total_ref_length > 0:
        stats.mean_coverage /= <double>stats.total_ref_length

    # WCB: normalize by ref_length²
    ref_len_sq = <double>stats.total_ref_length * <double>stats.total_ref_length
    if ref_len_sq > 0:
        stats.wcb /= ref_len_sq

    # CPC: breadth × (1 - dust_mean)
    stats.cpc = stats.breadth * (1.0 - stats.dust_mean)

    # ORI: total_aligned_bases / bases_covered
    # (needs bases_covered which we should track separately)

cdef void fit_taxid_damage_model(TaxidProfileStats* stats) noexcept nogil:
    """Fit damage model to accumulated position counts."""
    cdef TaxidDamageAccum* d = &stats.damage
    cdef double baseline_rate, d1_rate, d1_excess
    cdef int i

    if d.total_weight <= 0 or d.n_interior <= 0:
        stats.damage_amplitude = 0.0
        stats.damage_baseline = 0.0
        stats.damage_log_bf = 0.0
        return

    # Estimate baseline from interior positions (where damage ≈ 0)
    baseline_rate = d.k_interior / d.n_interior if d.n_interior > 0 else 0.01
    stats.damage_baseline = baseline_rate

    # Estimate damage amplitude from position 1 (5' end)
    if d.n_5p[0] > 0:
        d1_rate = d.k_5p[0] / d.n_5p[0]
        d1_excess = d1_rate - baseline_rate
        stats.damage_amplitude = d1_excess if d1_excess > 0 else 0.0
    else:
        stats.damage_amplitude = 0.0

    # Simple log Bayes factor approximation
    # Compare binomial likelihood under ancient (with damage) vs modern (baseline only)
    if stats.damage_amplitude > 0.01 and d.n_5p[0] > 10:
        # Rough approximation: higher amplitude = more evidence for ancient
        stats.damage_log_bf = d.n_5p[0] * (
            d1_rate * log(d1_rate / (baseline_rate + 0.001) + 0.001) +
            (1 - d1_rate) * log((1 - d1_rate) / (1 - baseline_rate + 0.001) + 0.001)
        )
    else:
        stats.damage_log_bf = 0.0

cdef double compute_p_present(
    TaxidProfileStats* stats,
    ProbProfileParams* params
) noexcept nogil:
    """Compute P(present|taxon) using logistic combination of detection metrics."""
    cdef double logit_score = 0.0
    cdef double p_present

    # Combine metrics with learned weights
    # Each metric is transformed to contribute to log-odds
    logit_score += params.w_wcb * stats.wcb
    logit_score += params.w_cpc * stats.cpc
    logit_score += params.w_ori * stats.ori          # negative weight
    logit_score += params.w_entropy * stats.norm_spatial_entropy
    logit_score += params.w_gini * stats.norm_gini   # negative weight
    logit_score += params.w_ani * (stats.ani_mean / 100.0)  # Normalize ANI to [0,1]

    # Convert logit to probability
    p_present = 1.0 / (1.0 + exp(-logit_score))

    stats.p_present = p_present
    stats.detection_score = logit_score

    return p_present


# =============================================================================
# Hierarchical Shrinkage (Pass 3: Top-down)
# =============================================================================

cdef void apply_hierarchical_shrinkage(
    TaxidProfileStats* child,
    TaxidProfileStats* parent,
    ProbProfileParams* params
) noexcept nogil:
    """Apply hierarchical Beta shrinkage from parent to child."""
    cdef double mu_parent, alpha_prior, beta_prior
    cdef double c_t, d_t, n_t

    # Get parent's posterior mean for P(ancient)
    if parent != NULL:
        mu_parent = parent.p_ancient
    else:
        mu_parent = params.alpha0 / (params.alpha0 + params.beta0)

    # Child's ancient/modern counts
    c_t = child.tad_ancient
    d_t = child.tad_modern
    n_t = c_t + d_t

    # Compute prior from parent (hierarchical shrinkage)
    alpha_prior = params.alpha0 + params.lambda_shrink * mu_parent
    beta_prior = params.beta0 + params.lambda_shrink * (1.0 - mu_parent)

    # Posterior update
    child.beta_alpha_post = alpha_prior + c_t
    child.beta_beta_post = beta_prior + d_t

    # Shrinkage weight diagnostic
    if params.lambda_shrink + n_t > 0:
        child.shrinkage_weight = params.lambda_shrink / (params.lambda_shrink + n_t)
    else:
        child.shrinkage_weight = 1.0

cdef void compute_p_ancient_posterior(
    TaxidProfileStats* stats,
    ProbProfileParams* params
) noexcept nogil:
    cdef double alpha = stats.beta_alpha_post
    cdef double beta_param = stats.beta_beta_post
    cdef double total = alpha + beta_param
    cdef double m, v, logit_m, sigma2, sigma, logit_low, logit_high
    cdef double exp_low, exp_high, ci_low, ci_high

    if total > 0:
        stats.p_ancient = alpha / total
    else:
        stats.p_ancient = 0.5

    if total > 0:
        m = alpha / total
    else:
        m = 0.5

    if m <= 0.0 or m >= 1.0 or total <= 0:
        stats.p_ancient_ci_lower = 0.0
        stats.p_ancient_ci_upper = 1.0
        return

    v = (alpha * beta_param) / (total * total * (total + 1.0))
    if v <= 0:
        stats.p_ancient_ci_lower = m
        stats.p_ancient_ci_upper = m
        return

    logit_m = log(m / (1.0 - m))
    sigma2 = v / (m * m * (1.0 - m) * (1.0 - m))
    if sigma2 <= 0:
        stats.p_ancient_ci_lower = m
        stats.p_ancient_ci_upper = m
        return

    sigma = sqrt(sigma2)
    logit_low = logit_m - 1.96 * sigma
    logit_high = logit_m + 1.96 * sigma
    exp_low = exp(-logit_low)
    exp_high = exp(-logit_high)
    ci_low = 1.0 / (1.0 + exp_low)
    ci_high = 1.0 / (1.0 + exp_high)
    stats.p_ancient_ci_lower = ci_low
    stats.p_ancient_ci_upper = ci_high


# =============================================================================
# Python-level helpers
# =============================================================================

def _py_beta_ci_approx(alpha: float, beta: float) -> tuple:
    """Python-level CI approximation for Beta distribution.

    Returns (lower_2.5%, upper_97.5%) using logit-normal approximation.
    """
    total = alpha + beta
    if total <= 0:
        return (0.0, 1.0)

    m = alpha / total
    if m <= 0.0 or m >= 1.0:
        return (0.0, 1.0)

    v = (alpha * beta) / (total * total * (total + 1.0))
    if v <= 0:
        return (m, m)

    logit_m = math.log(m / (1.0 - m))
    sigma2 = v / (m * m * (1.0 - m) * (1.0 - m))
    if sigma2 <= 0:
        return (m, m)

    sigma = math.sqrt(sigma2)
    logit_low = logit_m - 1.96 * sigma
    logit_high = logit_m + 1.96 * sigma

    p_low = 1.0 / (1.0 + math.exp(-logit_low))
    p_high = 1.0 / (1.0 + math.exp(-logit_high))

    return (p_low, p_high)


# =============================================================================
# Python-level TypedDict for output
# =============================================================================

class TaxidProfileEntry(TypedDict):
    """Per-taxid probabilistic profile entry."""
    # Identity
    taxid: int
    parent_taxid: Optional[int]
    rank: str
    name: str

    # Abundance
    n_refs: int
    n_reads: float
    tad_total: float
    tad_ancient: float
    tad_modern: float

    # Detection metrics
    breadth: float
    wcb: float
    cpc: float
    ori: float
    norm_spatial_entropy: float
    norm_gini: float
    ani_mean: float
    dust_mean: float

    # Damage
    damage_amplitude: float
    damage_baseline: float
    damage_log_bf: float

    # Probability outputs
    p_present: float
    p_present_ci_lower: float
    p_present_ci_upper: float
    p_ancient: float
    p_ancient_ci_lower: float
    p_ancient_ci_upper: float

    # Diagnostics
    shrinkage_weight: float
    detection_score: float


def _create_profile_entry(
    taxid: int,
    parent_taxid: Optional[int],
    rank: str,
    name: str
) -> TaxidProfileEntry:
    """Create an initialized TaxidProfileEntry."""
    return TaxidProfileEntry(
        taxid=taxid,
        parent_taxid=parent_taxid,
        rank=rank,
        name=name,
        n_refs=0,
        n_reads=0.0,
        tad_total=0.0,
        tad_ancient=0.0,
        tad_modern=0.0,
        breadth=0.0,
        wcb=0.0,
        cpc=0.0,
        ori=0.0,
        norm_spatial_entropy=0.0,
        norm_gini=0.0,
        ani_mean=0.0,
        dust_mean=0.0,
        damage_amplitude=0.0,
        damage_baseline=0.0,
        damage_log_bf=0.0,
        p_present=0.0,
        p_present_ci_lower=0.0,
        p_present_ci_upper=1.0,
        p_ancient=0.5,
        p_ancient_ci_lower=0.0,
        p_ancient_ci_upper=1.0,
        shrinkage_weight=1.0,
        detection_score=0.0,
    )


def _ensure_profile_entry(
    aggregates: Dict[int, TaxidProfileEntry],
    taxid: int,
    taxonomy: TaxonomyDatabase
) -> TaxidProfileEntry:
    """Get or create a TaxidProfileEntry for a taxid."""
    entry = aggregates.get(taxid)
    if entry is None:
        parent_taxid = taxonomy.get_parent(taxid)
        rank = taxonomy.get_rank(taxid) or "no rank"
        name = taxonomy.get_name(taxid) or f"taxid_{taxid}"
        entry = _create_profile_entry(taxid, parent_taxid, rank, name)
        aggregates[taxid] = entry
    return entry


# =============================================================================
# Reference Statistics Accumulation (from pre-computed RefStats)
# =============================================================================

def _accumulate_ref_stats_to_taxid(
    entry: TaxidProfileEntry,
    ref_stats: dict,
    gamma: float,
    phi_sum: float
) -> None:
    """Accumulate pre-computed reference statistics into a taxid entry."""
    entry['n_refs'] += 1
    entry['n_reads'] += phi_sum

    # TAD channel split
    tad = ref_stats.get('tax_abund_tad', 0.0)
    entry['tad_total'] += tad
    entry['tad_ancient'] += tad * gamma
    entry['tad_modern'] += tad * (1.0 - gamma)

    # Detection metrics (weighted by φ)
    ref_len = ref_stats.get('ref_length', 1)
    weight = phi_sum

    # Accumulate for later normalization
    if '_accum' not in entry:
        entry['_accum'] = {
            'wcb_sum': 0.0,
            'breadth_sum': 0.0,
            'entropy_sum': 0.0,
            'gini_sum': 0.0,
            'dust_sum': 0.0,
            'ani_sum': 0.0,
            'ori_sum': 0.0,
            'weight_sum': 0.0,
            'ref_length_sum': 0,
        }

    accum = entry['_accum']

    # Weighted accumulation
    accum['wcb_sum'] += weight * ref_stats.get('weighted_contiguity_breadth', 0.0)
    accum['breadth_sum'] += weight * ref_stats.get('breadth', 0.0)
    accum['entropy_sum'] += weight * ref_stats.get('norm_spatial_entropy', 0.0)
    accum['gini_sum'] += weight * ref_stats.get('norm_gini', 0.0)
    accum['dust_sum'] += weight * ref_stats.get('dust_mean', 0.0)
    accum['ani_sum'] += weight * ref_stats.get('ani_mean', 0.0)
    accum['ori_sum'] += weight * ref_stats.get('overlap_redundancy_index', 0.0)
    accum['weight_sum'] += weight
    accum['ref_length_sum'] += ref_len

    # Damage metrics (from damage model output)
    if '_damage_accum' not in entry:
        entry['_damage_accum'] = {
            'amplitude_sum': 0.0,
            'baseline_sum': 0.0,
            'log_bf_sum': 0.0,
            'weight_sum': 0.0,
        }

    d_accum = entry['_damage_accum']
    d_accum['amplitude_sum'] += weight * ref_stats.get('damage_amplitude', 0.0)
    d_accum['baseline_sum'] += weight * ref_stats.get('damage_baseline', 0.0)
    d_accum['log_bf_sum'] += weight * ref_stats.get('damage_log_bf', 0.0)
    d_accum['weight_sum'] += weight


def _finalize_detection_metrics(entry: TaxidProfileEntry) -> None:
    """Finalize detection metrics by normalizing accumulated values."""
    accum = entry.get('_accum')
    if accum is None or accum['weight_sum'] <= 0:
        return

    w = accum['weight_sum']
    entry['wcb'] = accum['wcb_sum'] / w
    entry['breadth'] = accum['breadth_sum'] / w
    entry['norm_spatial_entropy'] = accum['entropy_sum'] / w
    entry['norm_gini'] = accum['gini_sum'] / w
    entry['dust_mean'] = accum['dust_sum'] / w
    entry['ani_mean'] = accum['ani_sum'] / w
    entry['ori'] = accum['ori_sum'] / w

    # CPC = breadth × (1 - dust)
    entry['cpc'] = entry['breadth'] * (1.0 - entry['dust_mean'])

    # Damage metrics
    d_accum = entry.get('_damage_accum')
    if d_accum and d_accum['weight_sum'] > 0:
        dw = d_accum['weight_sum']
        entry['damage_amplitude'] = d_accum['amplitude_sum'] / dw
        entry['damage_baseline'] = d_accum['baseline_sum'] / dw
        entry['damage_log_bf'] = d_accum['log_bf_sum'] / dw

    # Clean up temporary fields
    if '_accum' in entry:
        del entry['_accum']
    if '_damage_accum' in entry:
        del entry['_damage_accum']


def _compute_p_present_py(
    entry: TaxidProfileEntry,
    w_wcb: float = 1.0,
    w_cpc: float = 1.0,
    w_ori: float = -1.0,
    w_entropy: float = 1.0,
    w_gini: float = -1.0,
    w_ani: float = 0.5
) -> float:
    """Compute P(present) using logistic combination of detection metrics."""
    logit_score = 0.0

    logit_score += w_wcb * entry['wcb']
    logit_score += w_cpc * entry['cpc']
    logit_score += w_ori * entry['ori']
    logit_score += w_entropy * entry['norm_spatial_entropy']
    logit_score += w_gini * entry['norm_gini']
    logit_score += w_ani * (entry['ani_mean'] / 100.0)

    p_present = 1.0 / (1.0 + math.exp(-logit_score))

    entry['p_present'] = p_present
    entry['detection_score'] = logit_score

    return p_present


# =============================================================================
# Main Entry Points
# =============================================================================

def compute_probabilistic_profile(
    ref_stats_dict: Dict[str, dict],
    gamma_dict: Dict[str, float],
    taxonomy: TaxonomyDatabase,
    accession_to_taxid: Dict[str, int],
    alpha0: float = 1.0,
    beta0: float = 1.0,
    lambda_shrink: float = 10.0,
    kappa_dirichlet: float = 2.0,
    use_exact_quantiles: bool = False,
    verbose: bool = False
) -> Dict[int, TaxidProfileEntry]:
    """
    Compute probabilistic taxonomic profile.

    This version uses pre-computed reference statistics from the filter stage.
    For de novo computation from BAM, use compute_profile_from_bam().

    Args:
        ref_stats_dict: Per-reference statistics (from filter stage)
        gamma_dict: Per-reference P(ancient) from hierarchical EM
        taxonomy: TaxonomyDatabase instance
        accession_to_taxid: Mapping from reference accession to taxid
        alpha0: Global Beta prior alpha
        beta0: Global Beta prior beta
        lambda_shrink: Hierarchical shrinkage strength
        kappa_dirichlet: Dirichlet smoothing for TAD
        use_exact_quantiles: Use exact Beta quantiles
        verbose: Enable verbose logging

    Returns:
        Dictionary mapping taxid -> TaxidProfileEntry
    """
    if verbose:
        bf_logging.log(LOG_TAG, "Starting probabilistic profile computation")
        bf_logging.log(LOG_TAG, f"Hyperparameters: α0={alpha0}, β0={beta0}, λ={lambda_shrink}")

    aggregates: Dict[int, TaxidProfileEntry] = {}

    # Pass 1: Bottom-up accumulation
    for ref_name, stats in ref_stats_dict.items():
        taxid = accession_to_taxid.get(ref_name)
        if taxid is None or taxid < 0:
            continue

        gamma = gamma_dict.get(ref_name, 0.5)
        phi_sum = stats.get('em_posterior_sum', stats.get('total_reads', 0.0))

        if phi_sum <= 0:
            continue

        # Accumulate to leaf taxid
        entry = _ensure_profile_entry(aggregates, taxid, taxonomy)
        _accumulate_ref_stats_to_taxid(entry, stats, gamma, phi_sum)

        # Propagate to ancestors
        lineage = taxonomy.get_lineage(taxid)
        if lineage:
            for ancestor in lineage:
                if ancestor == taxid:
                    continue
                entry = _ensure_profile_entry(aggregates, ancestor, taxonomy)
                _accumulate_ref_stats_to_taxid(entry, stats, gamma, phi_sum)

    if not aggregates:
        if verbose:
            bf_logging.log(LOG_TAG, "No taxids accumulated")
        return {}

    if verbose:
        bf_logging.log(LOG_TAG, f"Pass 1: Accumulated stats for {len(aggregates)} taxids")

    # Pass 2: Finalize detection metrics and compute P(present)
    for entry in aggregates.values():
        _finalize_detection_metrics(entry)
        _compute_p_present_py(entry)

    if verbose:
        bf_logging.log(LOG_TAG, "Pass 2: Computed detection metrics and P(present)")

    # Pass 3: Hierarchical Beta shrinkage for P(ancient)
    ordered_taxids = sorted(aggregates.keys(),
                           key=lambda t: len(taxonomy.get_lineage(t) or []))

    for taxid in ordered_taxids:
        entry = aggregates[taxid]
        parent_id = entry['parent_taxid']

        # Get parent's P(ancient)
        if parent_id is None or parent_id not in aggregates:
            mu_parent = alpha0 / (alpha0 + beta0)
        else:
            mu_parent = aggregates[parent_id]['p_ancient']

        # Compute prior from parent
        alpha_prior = alpha0 + lambda_shrink * mu_parent
        beta_prior = beta0 + lambda_shrink * (1.0 - mu_parent)

        # Ancient/modern counts
        c_t = entry['tad_ancient']
        d_t = entry['tad_modern']
        n_t = c_t + d_t

        # Posterior
        alpha_post = alpha_prior + c_t
        beta_post = beta_prior + d_t
        total = alpha_post + beta_post

        if total > 0:
            entry['p_ancient'] = alpha_post / total
        else:
            entry['p_ancient'] = 0.5

        # Credible intervals
        if use_exact_quantiles:
            entry['p_ancient_ci_lower'] = beta_quantile(alpha_post, beta_post, 0.025)
            entry['p_ancient_ci_upper'] = beta_quantile(alpha_post, beta_post, 0.975)
        else:
            ci_low, ci_high = _py_beta_ci_approx(alpha_post, beta_post)
            entry['p_ancient_ci_lower'] = ci_low
            entry['p_ancient_ci_upper'] = ci_high

        # Shrinkage weight
        if lambda_shrink + n_t > 0:
            entry['shrinkage_weight'] = lambda_shrink / (lambda_shrink + n_t)

    if verbose:
        bf_logging.log(LOG_TAG, f"Pass 3: Applied hierarchical shrinkage for P(ancient)")
        bf_logging.log(LOG_TAG, f"Completed profile for {len(aggregates)} taxids")

    return aggregates


# =============================================================================
# Output Formatting
# =============================================================================

def compute_profile_from_taxon_stats(
    taxon_stats: Dict[int, dict],
    taxonomy: TaxonomyDatabase,
    alpha0: float = 1.0,
    beta0: float = 1.0,
    lambda_shrink: float = 10.0,
    kappa_dirichlet: float = 2.0,
    use_exact_quantiles: bool = False,
    verbose: bool = False
) -> Dict[int, TaxidProfileEntry]:
    """
    Compute probabilistic profile from pre-computed taxon-level statistics.

    This function takes the taxon stats dictionary (from LCA stats processing)
    and computes P(ancient) and P(present) for each taxon using Bayesian
    hierarchical shrinkage.

    Args:
        taxon_stats: Dict mapping taxid -> stats dict with keys:
            - total_reads, total_alns (or n_reads, n_alns)
            - breadth, coverage_mean, coverage_covered_mean
            - norm_spatial_entropy, norm_gini
            - weighted_contiguity_breadth, overlap_redundancy_index
            - dust_mean, read_ani_mean
            - authenticity_score, authenticity_pvalue
            - tax_abund_tad
        taxonomy: TaxonomyDatabase instance
        alpha0: Global Beta prior alpha
        beta0: Global Beta prior beta
        lambda_shrink: Hierarchical shrinkage strength
        kappa_dirichlet: Dirichlet smoothing for TAD
        use_exact_quantiles: Use exact Beta quantiles
        verbose: Enable verbose logging

    Returns:
        Dictionary mapping taxid -> TaxidProfileEntry
    """
    if verbose:
        bf_logging.log(LOG_TAG, "Computing probabilistic profile from %d taxids", len(taxon_stats))

    aggregates: Dict[int, TaxidProfileEntry] = {}

    # Pass 1: Create profile entries from taxon stats
    for taxid, stats in taxon_stats.items():
        n_reads = stats.get('total_reads', stats.get('n_reads', 0))
        if n_reads <= 0:
            continue

        parent_taxid = taxonomy.get_parent(taxid)
        rank = taxonomy.get_rank(taxid) or "no rank"
        name = taxonomy.get_name(taxid) or f"taxid_{taxid}"

        entry = _create_profile_entry(taxid, parent_taxid, rank, name)
        entry['n_refs'] = stats.get('n_refs', 1)
        entry['n_reads'] = n_reads

        # TAD abundance
        tad = stats.get('tax_abund_tad', 0.0)
        entry['tad_total'] = tad

        # Get ancient/modern split from damage model gamma_ancient
        # This is computed from C→T and G→A damage patterns
        gamma_est = stats.get('gamma_ancient', 0.5)

        # If gamma not available from damage model, fall back to authenticity_score
        if gamma_est == 0.5 and 'authenticity_score' in stats:
            auth_score = stats.get('authenticity_score', 0.0)
            if auth_score != 0:
                gamma_est = 1.0 / (1.0 + math.exp(-auth_score))

        entry['tad_ancient'] = tad * gamma_est
        entry['tad_modern'] = tad * (1.0 - gamma_est)

        # Detection metrics
        entry['breadth'] = stats.get('breadth', 0.0)
        entry['wcb'] = stats.get('weighted_contiguity_breadth', 0.0)
        entry['ori'] = stats.get('overlap_redundancy_index', 0.0)
        entry['norm_spatial_entropy'] = stats.get('norm_spatial_entropy', 0.0)
        entry['norm_gini'] = stats.get('norm_gini', 0.0)
        entry['ani_mean'] = stats.get('read_ani_mean', 0.0)
        entry['dust_mean'] = stats.get('dust_mean', 0.0)

        # CPC = breadth × (1 - dust)
        entry['cpc'] = entry['breadth'] * (1.0 - entry['dust_mean'])

        # Damage (if available from stats)
        entry['damage_amplitude'] = stats.get('damage_amplitude', 0.0)
        entry['damage_baseline'] = stats.get('damage_baseline', 0.0)
        entry['damage_log_bf'] = stats.get('damage_log_bf', 0.0)

        aggregates[taxid] = entry

    if not aggregates:
        if verbose:
            bf_logging.log(LOG_TAG, "No taxids to process")
        return {}

    if verbose:
        bf_logging.log(LOG_TAG, "Pass 1: Created %d profile entries", len(aggregates))

    # Pass 2: Compute P(present) for each taxid
    for entry in aggregates.values():
        _compute_p_present_py(entry)

    if verbose:
        bf_logging.log(LOG_TAG, "Pass 2: Computed P(present) for all taxids")

    # Pass 3: Hierarchical Beta shrinkage for P(ancient)
    ordered_taxids = sorted(aggregates.keys(),
                           key=lambda t: len(taxonomy.get_lineage(t) or []))

    for taxid in ordered_taxids:
        entry = aggregates[taxid]
        parent_id = entry['parent_taxid']

        # Get parent's P(ancient)
        if parent_id is None or parent_id not in aggregates:
            mu_parent = alpha0 / (alpha0 + beta0)
        else:
            mu_parent = aggregates[parent_id]['p_ancient']

        # Compute prior from parent
        alpha_prior = alpha0 + lambda_shrink * mu_parent
        beta_prior = beta0 + lambda_shrink * (1.0 - mu_parent)

        # Ancient/modern counts from TAD
        c_t = entry['tad_ancient']
        d_t = entry['tad_modern']
        n_t = c_t + d_t

        # Posterior
        alpha_post = alpha_prior + c_t
        beta_post = beta_prior + d_t
        total = alpha_post + beta_post

        if total > 0:
            entry['p_ancient'] = alpha_post / total
        else:
            entry['p_ancient'] = 0.5

        # Credible intervals
        if use_exact_quantiles:
            entry['p_ancient_ci_lower'] = beta_quantile(alpha_post, beta_post, 0.025)
            entry['p_ancient_ci_upper'] = beta_quantile(alpha_post, beta_post, 0.975)
        else:
            ci_low, ci_high = _py_beta_ci_approx(alpha_post, beta_post)
            entry['p_ancient_ci_lower'] = ci_low
            entry['p_ancient_ci_upper'] = ci_high

        # Shrinkage weight
        if lambda_shrink + n_t > 0:
            entry['shrinkage_weight'] = lambda_shrink / (lambda_shrink + n_t)

    if verbose:
        bf_logging.log(LOG_TAG, "Pass 3: Applied hierarchical shrinkage")
        bf_logging.log(LOG_TAG, "Profile complete for %d taxids", len(aggregates))

    return aggregates


def format_profile_tsv(
    aggregates: Dict[int, TaxidProfileEntry],
    taxonomy: TaxonomyDatabase,
    output_path: str,
    min_reads: float = 0.0,
    min_p_present: float = 0.0,
    verbose: bool = False
) -> int:
    """
    Write probabilistic profile to TSV file.

    Args:
        aggregates: TaxidProfileEntry dict
        taxonomy: TaxonomyDatabase for lineage strings
        output_path: Output TSV path
        min_reads: Minimum n_reads to include
        min_p_present: Minimum P(present) to include
        verbose: Enable verbose logging

    Returns:
        Number of taxids written
    """
    header = [
        "taxid", "parent_taxid", "rank", "name",
        "n_refs", "n_reads",
        "tad_total", "tad_ancient", "tad_modern",
        "breadth", "wcb", "cpc", "ori",
        "norm_spatial_entropy", "norm_gini",
        "ani_mean", "dust_mean",
        "damage_amplitude", "damage_baseline", "damage_log_bf",
        "p_present", "p_present_ci_lower", "p_present_ci_upper",
        "p_ancient", "p_ancient_ci_lower", "p_ancient_ci_upper",
        "shrinkage_weight", "detection_score"
    ]

    n_written = 0
    import gzip

    opener = gzip.open if output_path.endswith('.gz') else open

    with opener(output_path, 'wt') as f:
        f.write('\t'.join(header) + '\n')

        for taxid in sorted(aggregates.keys()):
            entry = aggregates[taxid]

            if entry['n_reads'] < min_reads:
                continue
            if entry['p_present'] < min_p_present:
                continue

            row = [
                str(entry['taxid']),
                str(entry['parent_taxid']) if entry['parent_taxid'] else '',
                entry['rank'],
                entry['name'],
                str(entry['n_refs']),
                f"{entry['n_reads']:.2f}",
                f"{entry['tad_total']:.6f}",
                f"{entry['tad_ancient']:.6f}",
                f"{entry['tad_modern']:.6f}",
                f"{entry['breadth']:.6f}",
                f"{entry['wcb']:.6f}",
                f"{entry['cpc']:.6f}",
                f"{entry['ori']:.4f}",
                f"{entry['norm_spatial_entropy']:.6f}",
                f"{entry['norm_gini']:.6f}",
                f"{entry['ani_mean']:.4f}",
                f"{entry['dust_mean']:.4f}",
                f"{entry['damage_amplitude']:.6f}",
                f"{entry['damage_baseline']:.6f}",
                f"{entry['damage_log_bf']:.4f}",
                f"{entry['p_present']:.6f}",
                f"{entry['p_present_ci_lower']:.6f}",
                f"{entry['p_present_ci_upper']:.6f}",
                f"{entry['p_ancient']:.6f}",
                f"{entry['p_ancient_ci_lower']:.6f}",
                f"{entry['p_ancient_ci_upper']:.6f}",
                f"{entry['shrinkage_weight']:.4f}",
                f"{entry['detection_score']:.4f}",
            ]
            f.write('\t'.join(row) + '\n')
            n_written += 1

    if verbose:
        bf_logging.log(LOG_TAG, f"Wrote {n_written} taxids to {output_path}")

    return n_written
