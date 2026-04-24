# cython: language_level=3
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

# Forward-declare HTSlib types here so pxd doesn't conflict with module-local externs.
# Individual .pyx files provide full definitions when they need to access fields.
from bam_filter.processor_types cimport (
    BGZF,
    htsFile,
    samFile,
    bam1_t,
    sam_hdr_t,
    sam_hdr_tid2name,
    sam_hdr_tid2len,
    hts_idx_t,
    hts_itr_t,
    sam_index_load,
    hts_idx_destroy,
    hts_idx_get_stat,
    sam_itr_queryi,
)


# ReferenceMapping is now defined in processor_mapping.pxd
from bam_filter.processor_mapping cimport ReferenceMapping



# Fast streaming BAM writer helpers (declared so other modules can use them)
cdef struct WriteBatch:
    bam1_t** records
    uint64_t* pool_indices
    uint32_t count
    uint32_t capacity
    uint32_t reference_id



# NOTE: write-batch helper functions were moved to
# `bam_filter.processor_bam_writer` to keep the core processor module small.
# Other modules should cimport the implementations from
# `bam_filter.processor_bam_writer.pxd` instead of relying on these prototypes.

# -----------------------------------------------------------------------------
# Export Alignment and MemoryPool structs so other Cython modules can access fields
# (graph_analysis requires direct access to MemoryPool fields).
# Keep this declaration in the .pxd so it becomes a complete type for cimports.
# -----------------------------------------------------------------------------
from libc.stddef cimport size_t
from libc.stdint cimport int64_t, uint16_t, uint8_t

# -----------------------------------------------------------------------------
# Memory-optimized alignment storage using split arrays
# Core struct (16 bytes) + optional arrays allocated only when needed
# -----------------------------------------------------------------------------

# Core alignment data - always allocated (16 bytes, cache-line friendly)
cdef struct AlignmentCore:
    uint32_t reference_index      # Reference sequence index
    float    alignment_score      # Log-likelihood alignment score
    uint32_t alignment_position   # Start position in reference
    uint16_t aligned_length       # Total aligned bases
    uint16_t match_count          # Exact matches (for AN tag)

# Hierarchical EM data - allocated only when hierarchical_em_enabled (12 bytes)
cdef struct HierarchicalData:
    float damage_llr              # log(L_ancient/L_modern) position-specific
    float log_L_anc               # log P(alignment | ancient DNA model)
    float log_L_mod               # log P(alignment | modern DNA model)

# Damage counts for gamma update - allocated only when needed (4 bytes packed)
cdef struct DamageCounts:
    uint8_t ct_5p_count           # C→T mismatches in first 8bp from 5' end
    uint8_t ga_3p_count           # G→A mismatches in last 8bp from 3' end
    uint8_t c_at_5p_count         # C bases in reference at first 8bp (damage zone)
    uint8_t g_at_3p_count         # G bases in reference at last 8bp (damage zone)

# BAM writer auxiliary data - allocated only for BAM output (8 bytes)
cdef struct BAMWriterAux:
    float    pmd_score            # PMD score for PM:f tag output
    float    corrected_ani        # Damage-corrected ANI percentage (0-100, for DA:f tag)

# Legacy Alignment struct for backwards compatibility during transition
# TODO: Remove after full migration to split arrays
cdef struct Alignment:
    uint32_t read_index
    uint32_t reference_index
    uint32_t alignment_position
    float    alignment_score
    float    pmd_score
    uint16_t aligned_length       # Total aligned bases
    uint16_t match_count          # Exact matches
    uint8_t  ct_5p_count          # C→T mismatches in first 8bp from 5' end
    uint8_t  ga_3p_count          # G→A mismatches in last 8bp from 3' end
    uint8_t  c_at_5p_count        # C bases in reference at first 8bp (5' damage zone)
    uint8_t  g_at_3p_count        # G bases in reference at last 8bp (3' damage zone)
    float    corrected_ani        # Damage-corrected ANI percentage (0-100)
    uint8_t  passes_ani_filter    # 1 if passes corrected ANI threshold, 0 otherwise
    float    damage_llr           # log(L_ancient/L_modern) using position-specific D(z)
    float    log_L_anc            # log P(alignment | ancient DNA model)
    float    log_L_mod            # log P(alignment | modern DNA model)


cdef struct MemoryPool:
    void* base
    size_t capacity
    size_t used

    # Reusable aligned SQUAREM scratch block
    double* squarem_block
    size_t squarem_block_capacity

    # Primary memory pool (single allocation for all data)
    void* memory_pool
    size_t pool_capacity
    size_t pool_utilized

    # Core alignment data (now includes PMD)
    Alignment* alignments
    int64_t alignment_count          # Current number of alignments stored
    int64_t alignment_capacity       # Allocated capacity for alignments array
    bint alignments_is_external
    int64_t original_alignment_count

    # Split array storage (memory-optimized path)
    AlignmentCore* alignment_cores   # Always allocated (16 bytes/alignment)
    uint32_t* read_indices           # Always allocated (4 bytes/alignment) - needed for sort/index
    HierarchicalData* hierarchical   # Optional: only when hierarchical_em_enabled
    DamageCounts* damage_counts      # Optional: only when damage gamma update needed
    BAMWriterAux* bam_aux            # Optional: only when writing BAM with AN/DA tags
    bint use_split_arrays            # True if using split arrays instead of Alignment*
    size_t read_indices_alloc_size   # mmap size for read_indices array (0 if malloc)
    size_t hierarchical_alloc_size   # mmap size for hierarchical array (0 if malloc)
    size_t damage_counts_alloc_size  # mmap size for damage_counts array (0 if malloc)
    size_t bam_aux_alloc_size        # mmap size for bam_aux array (0 if malloc)

    # Memory management flags
    bint hash_data_dumped
    size_t original_alignment_size

    # Read-based indexing
    uint64_t* read_alignment_starts
    uint32_t* read_alignment_counts
    uint32_t unique_read_count
    uint32_t final_unique_reads

    # Reference data
    int64_t* reference_lengths
    uint32_t reference_count

    # Single unified buffer instead of separate arrays
    double* unified_buffer
    size_t unified_buffer_size

    # Offsets into unified buffer
    size_t reference_weights_offset
    size_t temp_buffer_A_offset
    size_t temp_buffer_B_offset

    # Algorithm metadata
    int32_t iteration_count
    double final_log_likelihood
    bint algorithm_converged
    bint memory_owner
    size_t mmap_allocation_size

    # ZP values for fast BAM writing
    float* precomputed_zp_values
    bint zp_values_computed

    # PMD control flag and model for BAM writing
    bint pmd_enabled_for_output     # Whether to write PM/AN/DA tags to BAM
    void* pmd_curve_ptr             # Pointer to PMDCurve for damage-corrected ANI

    # Hierarchical EM: Ancient/Modern Reference Classification
    bint hierarchical_em_enabled    # Whether hierarchical EM is active
    double* gamma_values            # γ_k: P(ancient | ref k) per reference [0,1]
    double* damage_amplitude        # A_k: damage amplitude per reference (like metaDMG A_b)
    double* damage_baseline         # b_k: baseline divergence per reference
    double* damage_log_bf           # log Bayes factor: log[P(data|ancient)/P(data|modern)]
    double* eta_values              # η_k = logit(γ_k) for SQUAREM extrapolation
    double* S_anc_accum             # Accumulated weighted ancient posterior per ref
    double* S_mod_accum             # Accumulated weighted modern posterior per ref
    double rho_ancient              # ρ: Global fraction of ancient references [0,1]
    double zeta_ancient             # ζ = logit(ρ) for SQUAREM extrapolation
    double rho_prior_alpha          # Beta prior α for ρ (default 1.0)
    double rho_prior_beta           # Beta prior β for ρ (default 1.0)
    float D_avg_5p                   # Average D(z) for 5' positions 1-8 (precomputed)
    float D_avg_3p                   # Average D(z) for 3' positions 1-8 (precomputed)
    float epsilon_error              # Sequencing error rate (default 0.01)

    # Per-reference Bayesian damage model (computed after EM)
    void* ref_damage_stats          # RefDamageStats* array [reference_count]
    void* damage_hyperparams        # DamageModelHyperparams* (global hyperparameters)

    # Coverage-Weighted Reference Priors (CWRP)
    double* authenticity_scores     # Per-reference authenticity scores [0,1] (hard/unweighted - path A)
    bint cwrp_enabled               # Whether CWRP is active
    double cwrp_lambda              # CWRP weight parameter

    # Posterior-Weighted Coverage Authenticity (Path B - used by CWRP inside EM)
    double* authenticity_scores_post  # Posterior-weighted authenticity [0,1] for CWRP
    double* norm_entropy_post         # Posterior-weighted normalized spatial entropy [0,1]
    double* norm_gini_post            # Posterior-weighted normalized Gini [0,1]
    int32_t auth_update_interval_post # Update B every N iterations (default 3)
    double auth_scale_post            # Sigmoid scale for posterior authenticity (default 4.0)
    int32_t auth_lambda_ramp_iters    # Ramp lambda from 0 to target over N iterations (default 5)
    double sample_pi_override         # Sample-level P(ancient) gate; 0.0 = auto from PMD

    # Ancientness Field (Full Fix: Joint γ-authenticity model)
    double* eta_ancientness         # Latent ancientness η_j per reference
    double* anc_feat_entropy        # Normalized spatial entropy [0,1]
    double* anc_feat_gini           # Normalized Gini coefficient [0,1]
    double* anc_feat_damage_5p      # Damage score at 5' end [0,1]
    double* anc_feat_damage_3p      # Damage score at 3' end [0,1]
    double* anc_feat_short_frac     # Fraction of short fragments [0,1]
    double* anc_feat_mean_length    # Mean fragment length
    double* anc_feat_read_count     # Posterior-weighted read count per ref
    bint iterative_auth_enabled     # Update authenticity during EM iterations
    int32_t auth_update_interval    # Update every N iterations (default 5)
    double damage_weight            # Weight for damage in ancientness (default 1.0)
    int32_t low_cov_floor_reads     # Shrink authenticity below this threshold
    double low_cov_shrink_tau       # Strength of shrinkage to neutral

    # Pooled scratch arrays for filtering (pooled & reused)
    float* scratch_read_max_probs
    int32_t* scratch_survivors_per_read
    int32_t scratch_unique_read_count
    
    # Processing statistics tracker
    ProcessingStats* stats


# AlignmentScoringConfig used by bam_processor; expose here so other modules
# (like stats.pyx) can build and pass a config to the inline filter helper.
cdef struct AlignmentScoringConfig:
    double minimum_read_identity
    int32_t minimum_read_length
    int32_t maximum_read_length
    double global_min_score
    double global_max_score
    bint calculate_pmd
    bint is_single_stranded
    int32_t damage_window  # Window size for damage correction (default 8, max 15)


# Functions implemented in processor.pyx
# NOTE: get_reference_weights, get_temp_buffer_A, get_temp_buffer_B are now imported from processor_em

# NOTE: PREFETCH_READ and PREFETCH_WRITE are now imported from processor_em

# Processing statistics tracker - follows data through entire pipeline
cdef struct ProcessingStats:
    # Stage 1: Initial BAM reading
    int64_t initial_total_alignments      # Total alignments in input BAM
    int64_t initial_total_reads           # Total unique reads in input BAM
    int64_t initial_total_references      # Total references with alignments
    
    # Stage 2: Quality filtering (ANI, length)
    int64_t post_quality_alignments       # After ANI/length filters
    int64_t post_quality_reads            # Unique reads after quality filter
    int64_t post_quality_references       # References after quality filter
    int64_t filtered_quality_alignments   # Removed by quality filters
    
    # Stage 3: EM algorithm
    int32_t em_iterations                 # Number of EM iterations
    bint em_converged                     # Did EM converge?
    double em_final_likelihood            # Final log-likelihood
    
    # Stage 4: Probability filtering
    int64_t post_probability_alignments   # After removing low-probability alignments
    int64_t post_probability_reads        # Unique reads after probability filter
    int64_t post_probability_references   # References after probability filter
    int64_t filtered_probability_alignments  # Removed by probability filter
    
    # Stage 5: Graph analysis
    int64_t graph_analysis_references     # References analyzed in graph
    int64_t graph_patterns_computed       # Network patterns computed

    # Stage 6: Unified filtering (coverage + information)
    int64_t post_unified_alignments       # After unified filtering
    int64_t post_unified_reads            # Unique reads after unified filter
    int64_t post_unified_references       # References after unified filter

    # Unified filtering breakdown
    int64_t filtered_coverage_only        # Failed coverage only
    int64_t filtered_information_only     # Failed information only
    int64_t filtered_both_criteria        # Failed both criteria

    int64_t alignments_removed_coverage   # Alignments from coverage-failed refs
    int64_t alignments_removed_information  # Alignments from info-failed refs
    int64_t alignments_removed_both       # Alignments from both-failed refs

    # Taxonomy-informed filtering (if enabled)
    int32_t taxonomy_enabled              # 1 if taxonomy filtering configured, 0 otherwise
    int64_t taxonomy_strict_removed       # Removed by strict taxonomy filtering
    int64_t taxonomy_weighted_count       # References with weighted anomaly scores
    int64_t taxonomy_second_chance_restored  # Restored by second-chance validation

    # Edge removal + misannotation stats
    int64_t edge_removal_edges_found
    int64_t edge_removal_alignments_removed
    int64_t edge_removal_references_affected
    int64_t edge_removal_refs_lost_all_edges
    int64_t edge_removal_refs_lost_most_edges

    int64_t misannotation_confident
    int64_t misannotation_likely
    int64_t misannotation_warning
    int64_t misannotation_removed_total
    int64_t misannotation_review_total

    # Stage 7: Final output
    int64_t final_alignments_written      # Alignments written to output BAM
    int64_t final_reads_written           # Unique reads in output BAM
    int64_t final_references_written      # References in output BAM
    
    # Summary metrics
    double overall_alignment_retention    # Percentage of alignments kept
    double overall_read_retention         # Percentage of reads kept
    double overall_reference_retention    # Percentage of references kept


# Use centralized helper definitions (no runtime symbol export required)
from .common_helpers cimport (
    pack_position_length,
    extract_position,
    extract_length,
    min_int64,
    max_int64,
    min_int32,
    max_int32,
    min_double,
    max_double,
    page_size,
)

# Exported constants used by other modules
# NOTE: Sorting helpers were moved to `bam_filter.processor_sort` to avoid
# duplicating implementations in this module.  If you need to cimport the
# sorting functions (radix_sort_uint64, radix_sort_compact_by_position,
# radix_sort_alignments_by_read_id, etc.) do so from `bam_filter.processor_sort`.
cdef uint32_t INVALID_SEQUENTIAL_ID
