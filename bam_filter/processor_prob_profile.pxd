# cython: language_level=3
"""
C-level declarations for probabilistic taxonomic profiling.

This module implements a Bayesian model for ancient DNA profiling with TWO distinct
probability outputs per taxon:

1. P(present|taxon) - Detection probability: Is this taxon truly in the sample?
   - Computed from coverage quality metrics: WCB, CPC, ORI, spatial entropy, breadth
   - Distinguishes true positives from false positives/contaminants

2. P(ancient|taxon) - Ancient probability: Is this taxon ancient DNA?
   - Computed from damage model fitted at taxon level
   - Uses hierarchical Beta shrinkage for taxa with few reads

Pipeline: REASSIGN → FILTER → PROFILE
- Reads BAM from filter stage (with ZP tags containing EM posteriors φ)
- Computes ALL statistics de novo (like lca-stats)
- Uses φ-weighting for multi-mapped reads (not LCA hard assignment)

Three-pass algorithm:
1. Pass 1 (Bottom-up): Accumulate coverage, damage, and metrics from alignments
2. Pass 2 (Compute): Calculate detection metrics and fit damage model per taxid
3. Pass 3 (Top-down): Apply hierarchical shrinkage for P(ancient)
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t
from bam_filter.stats cimport RefStats, RLECoverage, RLEInterval
from bam_filter.taxonomy_db cimport TaxonomyDB, TaxonomyDatabase

# =============================================================================
# Constants
# =============================================================================

cdef enum ProfileConstants:
    PROFILE_MAX_DAMAGE_POS = 20    # Track damage up to 20bp from read ends
    PROFILE_N_ENTROPY_BINS = 100   # Bins for spatial entropy calculation


# =============================================================================
# Per-Taxid Damage Accumulator (φ-weighted)
# =============================================================================

# Per-taxid damage counts accumulated from alignments with φ-weighting.
# Sufficient statistics for fitting the damage model at taxon level.
# All counts are weighted by EM posterior φ = P(ref|read).
cdef struct TaxidDamageAccum:
    # 5' C→T damage (positions 1-20 from 5' end)
    double n_5p[20]               # φ-weighted count of C bases at each position
    double k_5p[20]               # φ-weighted count of C→T mismatches

    # 3' G→A damage (positions 1-20 from 3' end)
    double n_3p[20]               # φ-weighted count of G bases at each position
    double k_3p[20]               # φ-weighted count of G→A mismatches

    # Interior baseline (positions 10-20 where damage ≈ 0)
    double n_interior             # φ-weighted opportunities in interior
    double k_interior             # φ-weighted mismatches in interior

    # Totals
    double total_weight           # Sum of φ weights (effective read count)
    uint64_t n_alignments         # Raw alignment count


# =============================================================================
# Per-Taxid Coverage Accumulator (for RLE-based stats)
# =============================================================================

# Per-taxid coverage accumulation for computing detection metrics.
# We accumulate coverage intervals from all contributing references,
# weighted by φ, to compute WCB, CPC, ORI, spatial entropy, etc.
cdef struct TaxidCoverageAccum:
    # Reference geometry
    int64_t total_ref_length      # Sum of contributing reference lengths
    int64_t n_refs                # Number of contributing references

    # Coverage basics (φ-weighted)
    double bases_covered          # φ-weighted covered bases
    double total_coverage         # φ-weighted total depth
    double total_coverage_sq      # For variance calculation

    # For TAD (Truncated Average Depth)
    double tad_sum                # φ-weighted TAD contribution
    double tad_sum_ancient        # TAD × γ channel
    double tad_sum_modern         # TAD × (1-γ) channel

    # Interval statistics for WCB
    double sum_interval_length_sq # Sum of (interval_length)² for WCB
    int64_t n_intervals           # Number of coverage intervals

    # Complexity for CPC
    double dust_sum               # φ-weighted DUST complexity sum
    double dust_weight            # Weight for DUST averaging

    # Spatial distribution metrics
    double spatial_entropy_sum    # φ-weighted spatial entropy sum
    double gini_sum               # φ-weighted gini sum
    double spatial_weight         # Weight for spatial metrics averaging

    # Read quality metrics
    double ani_sum                # φ-weighted ANI sum
    double ani_weight             # Weight for ANI averaging
    double mapq_sum               # φ-weighted MAPQ sum

    # For ORI (Overlap Redundancy Index)
    double total_aligned_bases    # φ-weighted aligned bases


# =============================================================================
# Per-Taxid Profile Statistics (comprehensive)
# =============================================================================

# Comprehensive per-taxid statistics for probabilistic profiling.
# Combines:
# - Detection metrics for P(present): WCB, CPC, ORI, entropy, breadth
# - Damage accumulation for P(ancient): position counts
# - Abundance metrics: TAD split by ancient/modern channel
cdef struct TaxidProfileStats:
    # Identity
    int32_t taxid
    int32_t parent_taxid
    int32_t rank_id
    int32_t n_refs

    # =====================================================
    # COVERAGE & ABUNDANCE (for composition)
    # =====================================================
    int64_t total_ref_length      # Total reference length
    double n_reads                # φ-weighted read count (effective reads)
    double n_alns                 # φ-weighted alignment count

    # TAD abundance (split by channel)
    double tad_total              # Total TAD abundance
    double tad_ancient            # TAD × γ (ancient channel)
    double tad_modern             # TAD × (1-γ) (modern channel)

    # Coverage statistics
    double breadth                # Fraction of genome with coverage
    double mean_coverage          # Mean depth across genome
    double mean_coverage_covered  # Mean depth in covered regions

    # =====================================================
    # DETECTION METRICS (for P(present))
    # =====================================================
    # These distinguish true taxa from false positives/contaminants

    # WCB: Weighted Contiguity Breadth
    # High WCB = large contiguous blocks = likely true taxon
    # Low WCB = scattered fragments = likely false positive
    double wcb                    # sum(interval_len²) / ref_length²

    # CPC: Complexity-Penalized Coverage
    # Penalizes coverage in low-complexity (high DUST) regions
    # Contaminants cluster in repetitive regions
    double cpc                    # breadth × (1 - dust_mean)

    # ORI: Overlap Redundancy Index
    # High ORI = excessive overlap = suspicious (short repetitive fragments)
    double ori                    # total_aligned_bases / bases_covered

    # Spatial entropy: uniformity of coverage distribution
    # High entropy = even spread = true taxon
    # Low entropy = clustered = likely contaminant
    double norm_spatial_entropy   # Normalized to [0,1]

    # Gini coefficient: coverage inequality
    # High gini = uneven coverage = suspicious
    # Low gini = even coverage = likely true taxon
    double norm_gini              # Normalized gini [0,1]

    # Read quality
    double ani_mean               # Average nucleotide identity
    double ani_corrected          # Damage-corrected ANI
    double mapq_mean              # Mean mapping quality
    double dust_mean              # Mean DUST complexity (low = repetitive)

    # =====================================================
    # DAMAGE ACCUMULATION (for P(ancient))
    # =====================================================
    TaxidDamageAccum damage       # Position-specific damage counts

    # Fitted damage model (after Pass 2)
    double damage_amplitude       # D(1) - damage at position 1
    double damage_baseline        # Background mismatch rate
    double damage_decay           # Exponential decay rate

    # =====================================================
    # PROBABILITY OUTPUTS
    # =====================================================

    # P(present|taxon) - Detection probability
    double p_present              # Posterior P(taxon present in sample)
    double p_present_ci_lower     # 2.5% credible interval
    double p_present_ci_upper     # 97.5% credible interval

    # P(ancient|taxon) - Ancient DNA probability
    double p_ancient              # Posterior P(taxon is ancient)
    double p_ancient_ci_lower     # 2.5% credible interval
    double p_ancient_ci_upper     # 97.5% credible interval

    # Beta posterior parameters for P(ancient)
    double beta_alpha_post        # Posterior alpha
    double beta_beta_post         # Posterior beta

    # Log Bayes factor for damage (ancient vs modern)
    double damage_log_bf          # log[P(data|ancient) / P(data|modern)]

    # Diagnostics
    double shrinkage_weight       # Hierarchical shrinkage applied
    double detection_score        # Combined detection evidence


# =============================================================================
# Hyperparameters
# =============================================================================

# Hyperparameters for probabilistic profiling model.
cdef struct ProbProfileParams:
    # P(ancient) Beta prior
    double alpha0                 # Global Beta prior alpha (default: 1.0)
    double beta0                  # Global Beta prior beta (default: 1.0)
    double lambda_shrink          # Hierarchical shrinkage strength (default: 10.0)

    # TAD normalization
    double kappa_dirichlet        # Dirichlet smoothing for TAD (default: 2.0)

    # Credible interval quantiles
    double ci_lower_quantile      # CI lower bound (default: 0.025)
    double ci_upper_quantile      # CI upper bound (default: 0.975)

    # Detection model weights (for P(present) logistic combination)
    double w_wcb                  # Weight for WCB (default: 1.0, positive = good)
    double w_cpc                  # Weight for CPC (default: 1.0, positive = good)
    double w_ori                  # Weight for ORI (default: -1.0, negative = bad)
    double w_entropy              # Weight for spatial entropy (default: 1.0, positive = good)
    double w_gini                 # Weight for gini (default: -1.0, negative = bad)
    double w_ani                  # Weight for ANI (default: 0.5, positive = good)

    # Minimum thresholds
    double min_reads              # Minimum φ-weighted reads to report (default: 1.0)
    double min_p_present          # Minimum P(present) to report (default: 0.1)


# =============================================================================
# Beta Distribution Functions
# =============================================================================

cdef double beta_quantile(double alpha, double beta, double p) noexcept nogil
cdef double beta_mean(double alpha, double beta) noexcept nogil
cdef double beta_variance(double alpha, double beta) noexcept nogil
cdef double log_beta_function(double alpha, double beta) noexcept nogil

# CI approximations (needed for nogil calls)
cdef double beta_ci_lower_approx(double alpha, double beta) noexcept nogil
cdef double beta_ci_upper_approx(double alpha, double beta) noexcept nogil


# =============================================================================
# Initialization Functions
# =============================================================================

cdef void init_prob_profile_params(ProbProfileParams* params) noexcept nogil
cdef void init_taxid_profile_stats(TaxidProfileStats* stats, int32_t taxid) noexcept nogil
cdef void init_taxid_damage_accum(TaxidDamageAccum* accum) noexcept nogil
cdef void init_taxid_coverage_accum(TaxidCoverageAccum* accum) noexcept nogil


# =============================================================================
# Accumulation Functions (Pass 1: Bottom-up)
# =============================================================================

cdef void accumulate_alignment_to_taxid(
    TaxidProfileStats* stats,
    double phi,                   # EM posterior P(ref|read)
    double gamma,                 # P(ancient|ref)
    int64_t ref_length,
    int64_t aln_start,
    int64_t aln_end,
    double ani,
    double mapq,
    double dust,
    # Damage counts from this alignment
    const double* k_5p,           # C→T counts at positions 1-20
    const double* n_5p,           # C opportunities at positions 1-20
    const double* k_3p,           # G→A counts at positions 1-20
    const double* n_3p            # G opportunities at positions 1-20
) noexcept nogil


cdef void accumulate_coverage_interval(
    TaxidProfileStats* stats,
    double phi,                   # Weight
    int64_t interval_start,
    int64_t interval_end,
    int32_t depth
) noexcept nogil


# =============================================================================
# Computation Functions (Pass 2: Compute metrics)
# =============================================================================

cdef void compute_detection_metrics(TaxidProfileStats* stats) noexcept nogil
cdef void fit_taxid_damage_model(TaxidProfileStats* stats) noexcept nogil
cdef double compute_p_present(
    TaxidProfileStats* stats,
    ProbProfileParams* params
) noexcept nogil


# =============================================================================
# Hierarchical Shrinkage (Pass 3: Top-down)
# =============================================================================

cdef void apply_hierarchical_shrinkage(
    TaxidProfileStats* child,
    TaxidProfileStats* parent,
    ProbProfileParams* params
) noexcept nogil

cdef void compute_p_ancient_posterior(
    TaxidProfileStats* stats,
    ProbProfileParams* params
) noexcept nogil
