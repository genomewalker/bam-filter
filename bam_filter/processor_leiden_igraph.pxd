# cython: language_level=3
"""
Header file for processor_leiden_igraph module.
Declares public functions that can be cimported by other Cython modules.
"""

from bam_filter.processor_graph_ops cimport WeightedGraph, ReferencePattern
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReadIndex
from bam_filter.processor_igraph cimport *
from libc.stdint cimport int32_t, uint32_t
from libc.stdlib cimport free

# --- Leiden public types & helpers ---
cdef struct LeidenResults:
    uint32_t num_nodes
    char* keep_flag              # 1 = keep, 0 = remove
    uint32_t* community_membership  # community ID for each reference
    uint32_t* component_membership  # connected component ID for each reference
    uint32_t* node_degree        # Node degree (number of edges) from igraph
    float* community_cc_values   # CC value for the community each reference belongs to
    float* individual_cc_values  # Individual CC value for each reference (Barrat's method)
    float* cc_threshold_values   # Broken-stick threshold used for each reference's community
    float* anomaly_scores        # Anomaly score for each reference (for multi-metric methods like Isolation Forest)
    float* betweenness_centrality  # Betweenness centrality for each reference (bridge-ness metric)

cdef struct CommunityStructure:
    # NOTE: CommunityStructure was previously declared here but is no longer
    # exported as part of the stable Cython API. Keep a minimal placeholder in
    # the .pxd for documentation purposes only. Consumers should use
    # `leiden_clustering` and the `LeidenResults` structure instead.
    uint32_t* node_to_community    # Map: node index -> community ID
    uint32_t* community_sizes      # Number of nodes in each community
    uint64_t* community_in_weights # Internal edge weights within communities
    uint64_t* community_tot_weights # Total weights (in + out) for communities
    uint32_t num_communities       # Current number of communities
    uint32_t max_communities       # Allocated capacity

# Main Leiden clustering function (exposed here so callers can cimport leiden_clustering)
cdef LeidenResults* leiden_clustering(
    WeightedGraph* graph,
    ReferenceStats* ref_stats,
    double resolution,
    int32_t max_iterations,
    bint verbose,
    int32_t thread_count,
    int outlier_method,  # 0 = MAD (default), 1 = IQR, 2 = IFOREST, 3 = LOF, 4 = ZSCORE
    int32_t iforest_n_trees,
    uint32_t iforest_subsample_size,
    double leiden_anomaly_threshold,
    uint32_t iforest_random_seed,
    uint32_t lof_k,
    double lof_contamination,
    double zscore_threshold,
    uint32_t* exact_connection_counts,
    double* co_mapping_averages,
    uint64_t* max_co_mappings,
    double* neighbor_multimap_avg,
    uint32_t array_size
) except NULL nogil
# Graph statistics using igraph (NEW - much faster than custom implementation)

