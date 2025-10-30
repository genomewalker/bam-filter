# cython: language_level=3

"""
Processor graph C-level interface declarations.

This header exposes C structs and function prototypes used by the
graph-analysis implementation in :mod:`processor_graph.pyx`.

Notes
-----
- The structs defined here are intentionally compact and laid out to
    be shared across multiple Cython modules. Keep changes backward
    compatible with other cimports unless you coordinate broader updates.
"""
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t

# Import WeightedGraph from processor_graph_ops
from bam_filter.processor_graph_ops cimport WeightedGraph

# ==============================================================================
# ==============================================================================
# The ReferencePattern and GraphMetrics structures have been streamlined.
#
# Structure organization:
# - GraphMetrics: connection_count, avg_comappings_per_read, max_comappings, etc.
# - ReferencePattern: graph (GraphMetrics), component_id, unique_read_count, leiden fields
# ==============================================================================


# Graph-based network metrics for each reference
cdef struct GraphMetrics:
    # Self network properties
    uint32_t connection_count          # number of connected neighbors
    double avg_comappings_per_read     # mean co-mappings per read with neighbors
    uint32_t max_comappings            # max co-mappings observed
    uint32_t reads_with_comappings     # reads that have neighbors

    # Neighbor aggregate properties
    double neighbor_avg_multimap       # mean multimap rate of neighbors
    double neighbor_avg_connections    # mean connection count of neighbors

# Per-reference network pattern (streamlined)
cdef struct ReferencePattern:
    # Network topology
    uint32_t component_id              # connected component ID
    uint32_t component_size            # size of connected component

    # Network metrics
    GraphMetrics graph                 # graph-based metrics
    uint32_t node_degree               # Node degree (number of edges) from igraph

    # Read classification
    uint32_t unique_read_count         # reads mapping only to this reference

    # Leiden clustering results
    uint32_t leiden_community_id       # Leiden community ID
    float leiden_community_cc          # Average clustering coefficient for this reference's community
    float leiden_individual_cc         # Individual clustering coefficient for this reference (Barrat's method)
    float leiden_cc_threshold          # Broken-stick threshold used for filtering this reference's community
    char leiden_keep_flag              # 1=keep, 0=remove from Leiden filtering
    float leiden_anomaly_score         # Anomaly score from multi-metric outlier detection (0=normal, 1=anomalous)
    float betweenness_centrality       # Betweenness centrality (bridge-ness metric from igraph)

    # Taxonomy information
    int32_t taxid                      # Taxonomy ID from accession mapping (-1 if not found)
    int32_t taxid_rank_id              # Rank ID of this taxid (for quick comparisons)
    int32_t taxid_depth                # Depth in taxonomy tree
    char taxonomy_flag                 # Flag indicating taxonomy-based anomaly: 0=normal, 1=potential_contamination, 2=cross_domain, 3=kingdom_mismatch


cdef struct DatasetSummaryStats:
    uint64_t total_reads_processed
    uint64_t unique_reads
    uint64_t multimapping_reads
    uint64_t total_alignments_input
    uint64_t total_alignments_filtered
    uint32_t total_references
    uint32_t references_with_alignments
    double min_score
    double max_score
    double mean_score
    double score_variance
    double score_ci_lower_95
    double score_ci_upper_95
    uint64_t score_count
    bint pmd_enabled
    double mean_pmd_score
    double pmd_variance
    uint64_t pmd_count

cdef struct StreamingStats:
    double running_mean
    double M2
    uint32_t observation_count
    uint64_t max_observed
    uint32_t nonzero_count

# Per-reference read statistics
cdef struct ReferenceStats:
    # Read category counts
    uint32_t total_reads
    uint32_t unique_reads
    uint32_t repeat_reads
    uint32_t shared_reads
    uint32_t connections              # number of connected references
    
    # Alignment scoring statistics
    uint64_t alignment_count
    double score_mean
    double score_std
    double score_min
    double score_max
    double score_variance             # Temporary: for compatibility
    
    # PMD scoring statistics (if available)
    bint pmd_available
    double pmd_mean
    double pmd_std
    double pmd_min
    double pmd_max
    uint64_t pmd_nonzero_count
    double pmd_variance               # Temporary: for compatibility
    
    # Calculated percentages (temporary: for compatibility)
    double unique_percentage
    double repeat_percentage
    double shared_percentage

cdef struct ReadIndex:
    uint32_t** ref_to_reads
    uint32_t* ref_read_counts
    uint32_t* ref_offsets       # offset into read_buffer for each reference (per-ref base index)
    uint32_t* read_buffer
    uint32_t* read_entry_starts   # parallel to read_buffer: alignment start for each entry
    uint32_t* read_entry_counts   # parallel to read_buffer: alignment count for each entry
    uint32_t* new_to_old         # mapping: new_index -> original reference id (for locality)
    uint32_t* old_to_new         # mapping: original reference id -> new_index


# Additional structs used internally by graph_analysis
cdef struct SparseConnectivity:
    uint32_t* connected_refs     # Dynamic array of connected reference IDs
    uint32_t capacity           # Allocated capacity
    uint32_t count              # Current count

cdef struct ReadRefsIndex:
    uint32_t* buffer            # single large buffer holding all per-read lists
    uint32_t** read_ptrs       # pointers into buffer for each read
    uint32_t* counts           # number of unique refs for each read
    uint32_t read_count        # number of reads

cdef struct MultPairExtended:
    uint32_t neighbor_count  # Connected references
    uint32_t read_count      # Total reads
    uint32_t idx            # Reference index
    float conc              # Sharing intensity
    double net              # Average co-mappings

cdef struct RefReadPair:
    uint32_t ref_idx
    uint32_t read_count

# Forward declarations for types defined elsewhere
from bam_filter.processor cimport MemoryPool, Alignment
from bam_filter.processor_mapping cimport ReferenceMapping
cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

# Use the shared EMAlgorithmConfig from bam_types to avoid duplicate
# incompatible declarations across modules. This ensures Cython emits
# a single compatible C struct for all cimports.
from bam_filter.processor_types cimport EMAlgorithmConfig

# Configuration struct for graph analysis parameters
cdef struct GraphConfig:
    uint32_t max_display_count
    uint32_t samples_per_stratum

    # Performance parameters
    uint32_t initial_array_capacity
    uint32_t sort_threshold
    uint32_t chunk_size_factor

    # Analysis thresholds
    uint64_t high_connection_warning
    uint64_t extreme_connection_warning

    # Penalty calculation weights
    double multimap_penalty_weight
    double connection_penalty_weight
    double neighbor_penalty_weight
    double comapping_penalty_weight

    # Statistical parameters
    uint32_t penalty_bins
    uint32_t penalty_ranges
    double confidence_interval_z

    # Memory parameters
    double array_growth_factor


# Function prototypes
# Mark these as noexcept nogil to match implementations and avoid exception checks
cdef void calculate_dataset_summary_stats(MemoryPool* pool, DatasetSummaryStats* stats) noexcept nogil
cdef void calculate_reference_stats(MemoryPool* pool, sam_hdr_t* bam_header,
                                   ReferenceStats* ref_stats) noexcept nogil
cdef void print_pattern_summary(MemoryPool* memory_pool, sam_hdr_t* bam_header,
                                ReferenceMapping* mapping, ReferencePattern* pattern_data,
                                GraphConfig* gconfig, ReferenceStats* ref_stats) noexcept nogil

cdef int count_unique_refs_thread_local(uint32_t read_idx, MemoryPool* pool,
                                           uint32_t array_size, uint32_t* thread_counts,
                                           uint32_t* scratch, uint32_t scratch_capacity) noexcept nogil

cdef int fill_ref_to_reads_thread_local(uint32_t read_idx, MemoryPool* pool,
                                 uint32_t array_size, uint32_t** ref_to_reads, uint32_t* write_positions,
                                 uint32_t* scratch, uint32_t scratch_capacity) nogil

cdef WeightedGraph* analyze_reference_graph(MemoryPool* pool, ReferencePattern* pattern_data,
                                           int32_t min_read_count, EMAlgorithmConfig* config,
                                           sam_hdr_t* bam_header, ReferenceMapping* mapping,
                                           bint verbose, bint build_igraph, const char* tsv_export_path,
                                           uint32_t graph_min_edge_weight) noexcept nogil

# Helper functions for cluster-aware filtering
cdef ReadIndex* build_read_index_parallel(MemoryPool* pool, uint32_t array_size, int num_threads) noexcept nogil
cdef void destroy_read_index(ReadIndex* index) noexcept nogil

# Information-theoretic filtering functions
cdef void accumulate_co_mappings(uint32_t start_read, uint32_t end_read,
                                 MemoryPool* pool, 
                                 uint32_t* read_refs_buffer, uint32_t* connection_counts,
                                 double* co_mapping_sums, uint32_t* co_mapping_counts,
                                 uint64_t* max_co_mappings, double* co_mapping_squares,
                                 uint32_t array_size) noexcept nogil

# TSV writing helper (exported so a small wrapper module can expose it)
cdef int write_graph_tsv(MemoryPool* pool, sam_hdr_t* bam_header,
                        ReferenceMapping* mapping, ReferencePattern* pattern_data,
                        ReferenceStats* ref_stats,
                        uint32_t* total_reads, uint32_t* multimap_reads, uint64_t* alignments_per_ref,
                        uint32_t* exact_connection_counts, double* co_mapping_averages,
                        uint64_t* max_co_mappings, uint32_t* co_mapping_counts,
                        double* neighbor_multimap_avg, double* neighbor_connections_avg,
                        uint32_t* neighbor_counts, uint32_t array_size,
                        double dataset_median_connections, int32_t min_read_count,
                        bint include_clustering, int outlier_method, const char* tsv_path) noexcept nogil


cdef int _uint32_compare(const void* a, const void* b) noexcept nogil
