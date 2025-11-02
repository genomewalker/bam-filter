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

cdef struct Alignment:
    uint32_t read_index
    uint32_t reference_index
    uint32_t alignment_position
    float    alignment_score
    float    pmd_score


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

    # OPTIMIZED: Single unified buffer instead of separate arrays
    double* unified_buffer
    size_t unified_buffer_size

    # OPTIMIZED: Offsets into unified buffer
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

    # PMD control flag (for BAM writing)
    bint pmd_enabled_for_output     # Whether to write PM tags to BAM

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
