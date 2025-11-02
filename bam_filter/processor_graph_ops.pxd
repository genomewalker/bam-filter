# cython: language_level=3

"""
Graph operations module - header file

This module provides graph data structures and operations for co-mapping analysis.
Separated from community detection (Leiden/union-find) to maintain clear separation of concerns.

Functions:
- create_weighted_graph: Create empty graph structure
- destroy_weighted_graph: Free graph memory
- add_edge: Add or increment edge weight
- build_weighted_graph_from_alignments: Build graph from alignment data
- prune_low_weight_edges: Remove edges below threshold
- calculate_graph_statistics: Calculate all graph metrics from filtered graph
"""

from libc.stdint cimport uint32_t, uint64_t, int32_t

from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReadIndex, ReferencePattern
from bam_filter.processor_igraph cimport (
    igraph_t,
    igraph_vector_t,
    igraph_vector_int_t,
    igraph_vs_t,
    igraph_integer_t,
)


# Graph data structures
cdef struct GraphNode:
    uint32_t* neighbors      # Sorted array of neighbor node IDs (after pruning)
    uint32_t* weights        # Corresponding edge weights (after pruning)
    uint32_t degree          # Number of edges (after pruning)
    uint32_t original_degree # Number of edges before pruning (for detecting pruned isolated nodes)
    uint32_t* original_neighbors  # Original neighbors before pruning (only saved if node becomes isolated)
    uint32_t capacity        # Allocated capacity
    uint64_t node_weight     # Sum of all edge weights for this node


cdef struct WeightedGraph:
    GraphNode* nodes         # Array of nodes
    uint32_t num_nodes       # Number of nodes
    uint64_t total_weight    # Sum of all edge weights (2× for undirected)
    uint64_t num_edges       # Number of edges (counted once per direction)
    void* igraph_handle      # igraph_t* pointer for reuse in community detection
    void* weights_handle     # igraph_vector_t* pointer for edge weights (cached with igraph)

    # Connected nodes optimization (for skipping isolated nodes in community detection)
    uint32_t* connected_node_ids  # Array of node IDs with degree > 0
    uint32_t n_connected_nodes    # Number of connected nodes

    # TSV data (cached for writing after community detection)
    uint32_t* tsv_total_reads
    uint32_t* tsv_multimap_reads
    uint64_t* tsv_alignments_per_ref
    uint32_t* tsv_exact_connection_counts
    double* tsv_co_mapping_averages
    uint64_t* tsv_max_co_mappings
    uint32_t* tsv_co_mapping_counts
    double* tsv_neighbor_multimap_avg
    double* tsv_neighbor_connections_avg
    uint32_t* tsv_neighbor_counts
    double tsv_dataset_median_connections
    uint32_t tsv_array_size
    int32_t tsv_min_read_count


# Graph creation and destruction
cdef WeightedGraph* create_weighted_graph(uint32_t num_nodes) except NULL nogil
cdef void destroy_weighted_graph(WeightedGraph* graph) noexcept nogil

# Graph operations
cdef int add_edge(WeightedGraph* graph, uint32_t from_node, uint32_t to_node,
                  uint32_t weight) noexcept nogil

# Graph building from alignment data
cdef WeightedGraph* build_weighted_graph_from_alignments(
    MemoryPool* pool,
    ReferenceStats* ref_stats,
    ReadIndex* read_index,
    uint32_t array_size,
    int32_t min_read_count
) except NULL nogil



# Build igraph from WeightedGraph (after pruning)
cdef int build_igraph_from_weighted_graph(
    igraph_t* ig_graph,
    igraph_vector_t* ig_weights,
    WeightedGraph* graph,
    bint verbose,
    igraph_integer_t* out_n_components,
    igraph_integer_t* out_singletons
) noexcept nogil

# Graph filtering
cdef void prune_low_weight_edges(WeightedGraph* graph, uint32_t min_edge_weight) noexcept nogil

# Graph statistics calculation
cdef void calculate_graph_statistics(WeightedGraph* graph,
                                     MemoryPool* pool,
                                     ReferencePattern* pattern_data,
                                     uint32_t* exact_connection_counts,
                                     double* neighbor_multimap_avg,
                                     double* neighbor_connections_avg,
                                     uint32_t* neighbor_counts,
                                     float* dataset_multimap_fractions,
                                     ReferenceStats* ref_stats,
                                     uint32_t num_refs) noexcept nogil

# Direct igraph builder (convenience wrapper that produces igraph_t and weights)
cdef int build_igraph_from_read_index(
    void* ig_graph,               # igraph_t* (opaque pointer)
    void* ig_weights,             # igraph_vector_t* (opaque pointer)
    MemoryPool* pool,
    ReadIndex* read_index,
    ReferenceStats* ref_stats,
    uint32_t num_refs,
    uint32_t min_read_count,
    uint32_t min_edge_weight,
    int num_threads,
    int verbose
) except -1 nogil

# Direct-to-igraph builder that avoids allocating full WeightedGraph
cdef int build_igraph_direct_from_read_index(
    void* ig_graph,               # igraph_t* (opaque pointer)
    void* ig_weights,             # igraph_vector_t* (opaque pointer)
    MemoryPool* pool,
    ReadIndex* read_index,
    ReferenceStats* ref_stats,
    uint32_t num_refs,
    uint32_t min_read_count,
    uint32_t min_edge_weight,
    int num_threads,
    int verbose,
    igraph_integer_t* out_n_components,
    igraph_integer_t* out_singletons
) except -1 nogil

# Broken-stick threshold picker (exported so other modules can call it nogil)
cdef uint32_t pick_min_edge_weight_broken_stick(MemoryPool* pool, ReadIndex* read_index, ReferenceStats* ref_stats, uint32_t num_refs, uint32_t min_read_count, double tol, int verbose, double tail_percentile, uint32_t min_tail_size) noexcept nogil

# Extract neighbor lists from igraph (for taxonomy analysis when nodes array is NULL)
cdef int extract_neighbors_from_igraph(WeightedGraph* graph, uint32_t num_refs,
                                        uint32_t*** out_neighbor_lists,
                                        uint32_t** out_neighbor_counts,
                                        int verbose) noexcept nogil
