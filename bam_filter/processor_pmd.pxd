# cython: language_level=3
# -*- coding: utf-8 -*-
# PMD (Post-Mortem Damage) learning and damage-corrected ANI computation.
#
# This module provides:
# 1. Thread-local PMD statistics accumulation during batch processing
# 2. Damage curve D(z) estimation with regularization
# 3. Hierarchical ancient/modern reference classification
# 4. Damage-corrected ANI computation for BAM output
#
# The PMD model uses an exponential decay: D(z) = omega * (P * exp(-lambda*(z-1)) + C)
# where omega is the library-specific damage presence (0-1).

from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t

# =============================================================================
# Constants (as enum for cimport compatibility)
# =============================================================================

cdef enum PMDConstants:
    PMD_MAX_POSITION = 20       # Track damage up to 20bp from read ends
    PMD_NUM_CONTEXTS = 2        # CpG vs non-CpG
    PMD_NUM_ENDS = 2            # 5' (C→T) vs 3' (G→A)
    PMD_MAX_THREADS = 64        # Maximum supported threads

cdef enum PMDEndType:
    PMD_END_5P = 0              # 5' end (C→T damage)
    PMD_END_3P = 1              # 3' end (G→A damage)

cdef enum PMDContextType:
    PMD_CTX_NONCPG = 0          # Non-CpG context
    PMD_CTX_CPG = 1             # CpG context


# =============================================================================
# PMD Statistics Accumulator (thread-local, collected during batch processing)
# Memory: 8 arrays × 20 positions × 8 bytes = 1.28 KB per thread
# =============================================================================

cdef struct PMDStatsAccumulator:
    # 5' end C→T damage counts (positions 1-20 from 5' end)
    uint64_t n_5p_noncpg[20]      # Eligible C sites, non-CpG context
    uint64_t k_5p_noncpg[20]      # C→T mismatches, non-CpG context
    uint64_t n_5p_cpg[20]         # Eligible C sites, CpG context
    uint64_t k_5p_cpg[20]         # C→T mismatches, CpG context

    # 3' end G→A damage counts (positions 1-20 from 3' end)
    uint64_t n_3p_noncpg[20]      # Eligible G sites, non-CpG context
    uint64_t k_3p_noncpg[20]      # G→A mismatches, non-CpG context
    uint64_t n_3p_cpg[20]         # Eligible G sites, CpG context
    uint64_t k_3p_cpg[20]         # G→A mismatches, CpG context

    # Metadata
    uint64_t total_alignments     # Number of alignments processed
    uint64_t total_bases          # Total bases examined
    int32_t thread_id             # Thread that owns this accumulator


# Global PMD statistics (merged from all threads)
cdef struct PMDStatsGlobal:
    uint64_t n_5p_noncpg[20]
    uint64_t k_5p_noncpg[20]
    uint64_t n_5p_cpg[20]
    uint64_t k_5p_cpg[20]
    uint64_t n_3p_noncpg[20]
    uint64_t k_3p_noncpg[20]
    uint64_t n_3p_cpg[20]
    uint64_t k_3p_cpg[20]
    uint64_t total_alignments
    uint64_t total_bases


# =============================================================================
# PMD Damage Curve Model
# Model: D(z) = omega * (P * exp(-lambda * (z-1)) + C)
# Separate curves for 5'/3' ends and CpG/non-CpG contexts.
# =============================================================================

cdef struct PMDCurve:
    # Damage presence (0 = no damage, 1 = full damage pattern)
    float omega

    # Damage probabilities at each position (precomputed from parameters)
    float D_5p_noncpg[20]         # D(z) for 5' C→T, non-CpG
    float D_5p_cpg[20]            # D(z) for 5' C→T, CpG
    float D_3p_noncpg[20]         # D(z) for 3' G→A, non-CpG
    float D_3p_cpg[20]            # D(z) for 3' G→A, CpG

    # Fitted parameters (for reference/debugging)
    float P_5p_noncpg             # Initial excess damage rate
    float P_5p_cpg
    float P_3p_noncpg
    float P_3p_cpg
    float lambda_decay            # Exponential decay rate (shared)
    float C_background            # Background damage rate (shared)

    # Library type
    bint is_single_stranded       # True for ss, False for ds

    # Baseline error rate (from quality scores)
    float epsilon_baseline        # ~0.01 for Q20


# Parameters for PMD curve fitting
cdef struct PMDCurveParams:
    # Prior means (from literature or previous runs)
    float P_prior_mean            # Prior for P (default: 0.3)
    float P_prior_sd              # Prior SD for P (default: 0.15)
    float lambda_prior_mean       # Prior for lambda (default: 0.35 = DECAY 0.7)
    float lambda_prior_sd         # Prior SD for lambda (default: 0.2)
    float C_prior_mean            # Prior for C (default: 0.01)
    float C_prior_sd              # Prior SD for C (default: 0.005)

    # Beta prior for omega (shrinks toward 0 for no-damage samples)
    float omega_alpha             # Beta alpha (default: 0.5)
    float omega_beta              # Beta beta (default: 5.0)

    # Baseline error rate
    float epsilon                 # Sequencing error rate (default: 0.01)


# =============================================================================
# Hierarchical EM: Ancient/Modern Reference Classification
# =============================================================================

# Per-reference ancient/modern classification state
cdef struct PMDReferenceClass:
    float log_L_ancient           # Accumulated log-likelihood under ancient model
    float log_L_modern            # Accumulated log-likelihood under modern model
    float gamma                   # Posterior P(ancient | reads) after EM
    uint32_t read_count           # Number of reads mapped to this reference


# State for hierarchical ancient/modern EM
cdef struct PMDHierarchicalEM:
    # Per-reference classification (array of size reference_count)
    PMDReferenceClass* ref_classes
    uint32_t reference_count

    # Global ancient fraction
    float rho                     # P(reference is ancient)
    float rho_prior_alpha         # Beta prior alpha (default: 1.0)
    float rho_prior_beta          # Beta prior beta (default: 1.0)

    # Convergence tracking
    int iteration
    float prev_log_likelihood
    bint converged


# =============================================================================
# ANI Snapshot (compact per-alignment storage)
# Memory: 8 bytes per alignment (vs re-parsing MD tag)
# =============================================================================

cdef struct ANISnapshot:
    uint16_t aligned_length       # Total aligned bases (max 65535)
    uint16_t match_count          # Exact matches
    uint8_t ct_5p_count           # C→T mismatches in first 8bp from 5' end
    uint8_t ga_3p_count           # G→A mismatches in last 8bp from 3' end
    uint8_t other_mm_count        # Other mismatches (capped at 255)
    uint8_t flags                 # Bit flags: 0x01=CpG context available
    # Damage opportunity counts (for hierarchical EM)
    uint8_t c_at_5p_count         # C bases in reference at first 8bp (5' damage zone)
    uint8_t g_at_3p_count         # G bases in reference at last 8bp (3' damage zone)


# Result of ANI computation
cdef struct ANIResult:
    float raw_ani                 # Raw ANI (matches / aligned_length * 100)
    float corrected_ani           # Damage-corrected ANI
    float damage_contribution     # Sum of w_i for damage-eligible mismatches


# =============================================================================
# Complete PMD Model (frozen after finalization)
# =============================================================================

cdef struct PMDModel:
    # Learned damage curve
    PMDCurve curve

    # Fitting parameters used
    PMDCurveParams params

    # Global statistics (merged from all threads)
    PMDStatsGlobal stats

    # Hierarchical EM state (optional, allocated if hierarchical EM enabled)
    PMDHierarchicalEM* hierarchical

    # Model state
    bint stats_collected          # True after batch processing complete
    bint curve_fitted             # True after D(z) estimation
    bint finalized                # True after model is frozen (read-only)
    bint hierarchical_enabled     # True if ancient/modern classification active


# =============================================================================
# Thread-Local PMD Context (for batch processing)
# =============================================================================

# Per-thread context for PMD statistics collection
cdef struct PMDThreadContext:
    PMDStatsAccumulator stats     # Thread-local accumulator
    int32_t thread_id
    bint initialized


# Global PMD context managing all thread contexts and model
cdef struct PMDGlobalContext:
    # Thread-local contexts
    PMDThreadContext* thread_contexts
    int32_t num_threads
    int32_t max_threads

    # Global model (finalized after batch processing)
    PMDModel model

    # Configuration
    bint is_single_stranded
    bint collect_stats            # Whether to collect PMD stats
    bint enable_hierarchical      # Whether to run hierarchical EM

    # Memory management
    bint initialized
    bint owns_memory


# =============================================================================
# Function declarations
# =============================================================================

# Initialization and cleanup
cdef PMDGlobalContext* create_pmd_context(int32_t num_threads,
                                           bint is_single_stranded,
                                           bint enable_hierarchical) noexcept nogil

cdef void destroy_pmd_context(PMDGlobalContext* ctx) noexcept nogil

cdef void reset_pmd_stats(PMDGlobalContext* ctx) noexcept nogil

# Thread-local stats collection (called during MD parsing)
cdef void pmd_record_position(PMDStatsAccumulator* acc,
                               int end_type,        # PMD_END_5P or PMD_END_3P
                               int context_type,    # PMD_CTX_NONCPG or PMD_CTX_CPG
                               int position,        # 1-based position from end
                               bint is_mismatch) noexcept nogil

# Stats merging and curve fitting
cdef void merge_pmd_stats(PMDGlobalContext* ctx) noexcept nogil

cdef int fit_pmd_curve(PMDGlobalContext* ctx, PMDCurveParams* params) noexcept nogil

cdef void finalize_pmd_model(PMDGlobalContext* ctx) noexcept nogil

# ANI computation
cdef float compute_raw_ani(ANISnapshot* snapshot) noexcept nogil

cdef float compute_corrected_ani(ANISnapshot* snapshot,
                                  PMDCurve* curve,
                                  float epsilon) noexcept nogil

cdef ANIResult compute_ani_pair(ANISnapshot* snapshot,
                                 PMDCurve* curve,
                                 float epsilon) noexcept nogil

# Hierarchical EM helpers
cdef int init_hierarchical_em(PMDGlobalContext* ctx,
                               uint32_t reference_count) noexcept nogil

cdef void accumulate_ref_likelihood(PMDHierarchicalEM* hem,
                                     uint32_t ref_idx,
                                     float log_L_ancient,
                                     float log_L_modern) noexcept nogil

cdef void update_gamma_and_rho(PMDHierarchicalEM* hem) noexcept nogil

cdef float get_reference_gamma(PMDHierarchicalEM* hem, uint32_t ref_idx) noexcept nogil

# Likelihood computation for ancient vs modern
cdef float compute_log_likelihood_ancient(uint8_t ref_base, uint8_t read_base,
                                           float epsilon, float D_z) noexcept nogil

cdef float compute_log_likelihood_modern(uint8_t ref_base, uint8_t read_base,
                                          float epsilon) noexcept nogil

# Utility functions
cdef float get_damage_prob(PMDCurve* curve, int end_type, int context_type,
                           int position) noexcept nogil

cdef float compute_damage_weight(float D_z, float epsilon) noexcept nogil

# =============================================================================
# Apply PMD Corrections to Alignments (requires MemoryPool from processor.pxd)
# Note: These functions are declared here but need processor.pxd imports in .pyx
# =============================================================================

# Forward declaration - actual import in .pyx file
# cdef int64_t apply_pmd_corrections_to_pool(MemoryPool* pool, PMDCurve* curve,
#                                             float min_ani_threshold, float epsilon,
#                                             int num_threads) noexcept nogil
#
# cdef void update_alignment_scores_with_pmd(MemoryPool* pool, PMDCurve* curve,
#                                             float damage_score_bonus,
#                                             int num_threads) noexcept nogil

# =============================================================================
# Hierarchical EM: Per-alignment ancient/modern likelihood computation
# Note: These require MemoryPool and Alignment from processor.pxd
# =============================================================================

# Forward declarations for hierarchical EM functions (implemented in processor_pmd.pyx)
# These are used by processor_em.pyx for joint ancient/modern classification

# cdef double compute_alignment_log_L_ancient(Alignment* aln, float D_avg_5p,
#                                              float D_avg_3p, float epsilon) noexcept nogil
#
# cdef double compute_alignment_log_L_modern(Alignment* aln, float epsilon) noexcept nogil
#
# cdef void init_hierarchical_em_pool(MemoryPool* pool, PMDCurve* curve,
#                                      float epsilon) noexcept nogil
#
# cdef void free_hierarchical_em_pool(MemoryPool* pool) noexcept nogil
