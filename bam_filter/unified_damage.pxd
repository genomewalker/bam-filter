# cython: language_level=3
# -*- coding: utf-8 -*-
"""Header file for unified ancient DNA damage model.

Single model for all damage-related tasks:
- ANI correction
- Ancient/modern classification
- Taxonomic profiling
"""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t

# Constants
cdef enum:
    DAMAGE_MAX_POSITION = 20
    DAMAGE_INTERIOR_START = 10
    DAMAGE_FIT_POSITIONS = 9


# Per-reference damage statistics (accumulated during batch processing)
cdef struct RefDamageCounts:
    # Mismatch counts k[z] for z=0..19 (0-indexed, position 1-20)
    double k_5p[20]      # C→T at 5' end
    double k_3p[20]      # G→A at 3' end
    double k_5p_ctrl[20] # G→A at 5' end (control - should be flat)
    double k_3p_ctrl[20] # C→T at 3' end (control - should be flat)

    # Opportunity counts n[z]
    double n_5p[20]      # C opportunities at 5' end
    double n_3p[20]      # G opportunities at 3' end
    double n_5p_ctrl[20] # G opportunities at 5' end
    double n_3p_ctrl[20] # C opportunities at 3' end

    # Totals for ANI
    uint64_t total_aligned
    uint64_t total_matches
    uint64_t total_mismatches

    # Coverage
    double total_weight  # φ-weighted read count
    uint32_t n_alignments


# Per-reference fitted parameters and outputs
cdef struct RefDamageParams:
    # Fitted amplitudes (separate for 5' and 3')
    double delta_5p
    double delta_3p

    # Baselines (from interior positions)
    double baseline_5p
    double baseline_3p

    # Expected damage counts (from E-step)
    double Y_5p  # expected damage mismatches at 5'
    double Y_3p  # expected damage mismatches at 3'
    double B_5p  # expected baseline mismatches at 5'
    double B_3p  # expected baseline mismatches at 3'

    # Output scores
    double authenticity     # Y / (Y + B + 1), bounded [0,1]
    double p_ancient        # posterior P(ancient)
    double log_bf           # log Bayes factor
    double asymmetry        # (CT5 + GA3) - (GA5 + CT3) signal

    # ANI
    double ani_raw
    double ani_corrected
    double damage_correction

    # Quality flags
    uint8_t has_damage_evidence  # Y > threshold
    uint8_t is_overdispersed     # diagnostic flag
    uint8_t is_low_coverage      # < min opportunities

    # Posterior predictive diagnostics
    double chi_squared           # goodness-of-fit statistic
    double pp_pvalue             # posterior predictive p-value
    double dispersion_factor     # observed/expected variance ratio
    int df                       # degrees of freedom


# Sample-level parameters
cdef struct SampleDamageParams:
    # Decay parameter
    double tau             # fitted decay length
    double tau_prior_mean  # prior mean (log scale)
    double tau_prior_sd    # prior sd (log scale)
    double tau_min         # lower bound
    double tau_max         # upper bound

    # Global amplitude (shrinkage target)
    double mu_5p           # mean amplitude 5' end
    double mu_3p           # mean amplitude 3' end
    double alpha           # Gamma shrinkage strength

    # Sample-level baselines (fallback)
    double baseline_5p_global
    double baseline_3p_global

    # EM state
    int iteration
    bint converged
    double log_likelihood


# Main context holding all state
cdef struct UnifiedDamageContext:
    # Sample parameters
    SampleDamageParams sample

    # Per-reference arrays (allocated separately)
    RefDamageCounts* ref_counts
    RefDamageParams* ref_params
    uint32_t n_refs

    # Configuration
    bint is_single_stranded
    bint mask_cpg           # exclude CpG from counts
    bint use_controls       # use GA5/CT3 as controls
    int max_iterations
    double tol_tau
    double tol_delta

    # Memory ownership
    bint owns_memory


# Core functions (nogil for performance)
cdef UnifiedDamageContext* create_damage_context(
    uint32_t n_refs,
    bint is_single_stranded,
    bint mask_cpg
) noexcept nogil

cdef void destroy_damage_context(UnifiedDamageContext* ctx) noexcept nogil

cdef void reset_damage_counts(UnifiedDamageContext* ctx) noexcept nogil

cdef void accumulate_damage_position(
    RefDamageCounts* counts,
    int end_type,      # 0=5', 1=3'
    int is_control,    # 0=damage channel, 1=control channel
    int position,      # 1-based position from end
    double weight,     # φ weight (1.0 if unweighted)
    bint is_mismatch   # True if mismatch observed
) noexcept nogil

cdef void estimate_baselines(UnifiedDamageContext* ctx) noexcept nogil

cdef void initialize_parameters(UnifiedDamageContext* ctx) noexcept nogil

cdef int run_em_iteration(UnifiedDamageContext* ctx) noexcept nogil

cdef int fit_damage_model(UnifiedDamageContext* ctx) noexcept nogil

cdef void compute_outputs(UnifiedDamageContext* ctx) noexcept nogil

cdef double compute_damage_prob(
    double delta,
    double tau,
    int position  # 1-based
) noexcept nogil

cdef double compute_ani_correction(
    RefDamageCounts* counts,
    RefDamageParams* params,
    double tau
) noexcept nogil

cdef void set_tau_from_pmd(
    UnifiedDamageContext* ctx,
    double tau_value,
    bint fix_tau
) noexcept nogil

cdef void compute_posterior_predictive(
    UnifiedDamageContext* ctx
) noexcept nogil
