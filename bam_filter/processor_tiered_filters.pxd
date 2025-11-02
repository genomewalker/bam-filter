# cython: language_level=3
"""
Three-tier filtering approach:
- Tier 1: Structural role classification (PERIPHERAL, CORE, HUB, BRIDGE)
- Tier 2: Community coherence analysis
- Tier 3: Integrated decision matrix
"""

from libc.stdint cimport uint32_t, int32_t, uint8_t
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats
from bam_filter.taxonomy_db cimport TaxonomyDB


# Tier 1: Structural role enum
cpdef enum StructuralRole:
    ROLE_PERIPHERAL = 0  # Low degree, high CC - specific, likely clean
    ROLE_CORE = 1        # High degree, high CC - conserved, likely clean
    ROLE_HUB = 2         # High degree, low CC - promiscuous, suspicious
    ROLE_BRIDGE = 3      # High betweenness across communities - VERY suspicious

# Tier 3: Filtering decision enum
cpdef enum FilterDecision:
    DECISION_KEEP = 0      # Keep reference
    DECISION_REMOVE = 1    # Remove reference
    DECISION_REVIEW = 2    # Flag for manual review (not implemented yet)


# Tier 1: Classify structural role of a node
cdef StructuralRole classify_structural_role(
    float betweenness,
    float clustering_coefficient,
    uint32_t degree,
    uint32_t num_neighbor_communities,
    float betweenness_threshold,
    float cc_threshold,
    uint32_t hub_degree_threshold
) nogil


# Tier 2: Check taxonomy coherence of a community
cdef bint check_community_coherence(
    uint32_t* community_members,
    uint32_t n_members,
    ReferencePattern* pattern_data,
    TaxonomyDB* taxonomy_db,
    int32_t* lca_result_out
) nogil


# Tier 3: Make integrated filtering decision
cdef FilterDecision make_filtering_decision(
    StructuralRole role,
    uint8_t taxonomy_flag,
    bint community_coherent,
    bint strict_mode
) nogil


# Main entry point: Apply three-tier filtering
cdef int apply_tiered_filtering(
    ReferencePattern* pattern_data,
    ReferenceStats* ref_stats,
    uint32_t array_size,
    char* keep_flag,
    TaxonomyDB* taxonomy_db,
    float betweenness_threshold,
    float cc_threshold,
    uint32_t hub_degree_threshold,
    bint strict_mode,
    bint verbose,
    void* pool_handle,
    uint32_t** neighbor_lists,
    uint32_t* neighbor_counts,
    bint enable_edge_removal,
    bint flag_misannotations,
    char** out_alignment_keep_flags
) nogil
