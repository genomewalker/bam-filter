# cython: language_level=3
"""
Network QC metrics for evaluating EM quality based on graph structure.

Metrics:
- Dataset-level: N_eff, Gini, modularity, assortativity
- Per-reference: participation coefficient, neighbor entropy, z-score
"""

from libc.stdint cimport uint32_t, int32_t, uint64_t, int64_t, uint8_t
from libc.stddef cimport size_t


# ============================================================================
# Dataset-level QC metrics (computed once per EM run)
# ============================================================================

cdef struct NetworkQCDataset:
    # Diversity metrics
    double n_eff_weights           # Effective number of refs: 1/Σw_j²
    double n_eff_reads             # Effective number based on read counts
    double gini_weights            # Gini coefficient on EM weights
    double gini_reads              # Gini coefficient on read counts

    # Distribution metrics
    double entropy_weights         # Shannon entropy of weight distribution
    double max_weight              # Maximum single reference weight
    double top5_share              # % of total in top 5 references
    double top10_share             # % of total in top 10 references

    # Network structure metrics
    double modularity_taxonomy     # Modularity Q using taxonomy as partition
    double modularity_community    # Modularity Q using detected communities
    double assortativity_taxonomy  # Newman assortativity by taxonomy
    double global_clustering       # Global clustering coefficient

    # Mutual information metrics
    double mi_taxonomy_community   # MI between taxonomy and community partition
    double vi_taxonomy_community   # Variation of information

    # Read assignment quality
    double mean_posterior_entropy  # Mean per-read posterior entropy
    double high_entropy_fraction   # Fraction of reads with H > threshold
    double cross_community_mass    # Fraction of posterior mass crossing communities


# ============================================================================
# Per-reference QC metrics (computed for each reference)
# ============================================================================

cdef struct NetworkQCReference:
    # Neighbor taxonomic entropy: H_neighbor(j) = -Σ_a P_j(a) log P_j(a)
    # High entropy indicates taxonomically diverse neighbors
    # Note: HUB references are expected to have high entropy and are NOT flagged
    float neighbor_tax_entropy

    # Within-community degree z-score: z_j = (k_j^c - μ) / σ
    float within_community_zscore

    # Taxonomic ambiguity flag (replaces chimera detection)
    # Only set for non-HUB references with high neighbor entropy
    uint8_t tax_ambiguity_flag     # 0=clean, 1=biased, 2=mixed, 3=highly_mixed

    # Edge weight statistics
    float mean_edge_weight         # Mean edge weight to neighbors
    float max_edge_weight          # Max edge weight to any neighbor
    uint32_t cross_community_edges # Number of edges to other communities


# ============================================================================
# Configuration for QC computation
# ============================================================================

cdef struct NetworkQCConfig:
    # Entropy thresholds
    double high_entropy_threshold      # Threshold for "high entropy" reads (default: 2.0 bits)

    # Taxonomic ambiguity thresholds (replaces chimera detection)
    # These thresholds determine when neighbor entropy is considered "high"
    double entropy_biased_threshold    # Entropy above which ref is "biased" (default: 1.0 bits)
    double entropy_mixed_threshold     # Entropy above which ref is "mixed" (default: 2.0 bits)
    double entropy_highly_mixed_threshold  # Entropy above which ref is "highly_mixed" (default: 3.0 bits)

    # Taxonomy level for analysis
    int32_t taxonomy_rank_for_analysis # Rank ID for taxonomy grouping (default: genus=5)

    # Community detection
    bint use_taxonomy_as_community     # Use taxonomy instead of Leiden communities
    double community_resolution        # Resolution for Leiden (if used)


# ============================================================================
# Function declarations
# ============================================================================

# Dataset-level metrics
cdef NetworkQCDataset compute_dataset_qc_metrics(
    double* weights,
    uint32_t* read_counts,
    uint32_t n_refs,
    int32_t* community_ids,
    int32_t* taxonomy_ids,
    uint32_t* edge_list_src,
    uint32_t* edge_list_dst,
    double* edge_weights,
    uint64_t n_edges,
    NetworkQCConfig* config,
) noexcept nogil

# Per-reference metrics
cdef void compute_reference_qc_metrics(
    uint32_t ref_idx,
    NetworkQCReference* result,
    int32_t* community_ids,
    int32_t* taxonomy_ids,
    uint32_t* neighbors,
    double* neighbor_weights,
    uint32_t n_neighbors,
    double* community_strength_sums,  # Sum of edge weights per community
    double* community_strength_sq,    # Sum of squared edge weights per community
    uint32_t* community_counts,       # Number of refs per community
    uint32_t n_communities,
    int32_t ref_community,
    NetworkQCConfig* config,
) noexcept nogil

# Taxonomic ambiguity detection (replaces chimera detection)
# Note: structural_role array is used to skip HUB references (role=2)
cdef void compute_tax_ambiguity_flags(
    NetworkQCReference* ref_metrics,
    uint32_t n_refs,
    char* structural_roles,  # 0=PERIPHERAL, 1=CORE, 2=HUB, 3=BRIDGE
    NetworkQCConfig* config,
) noexcept nogil

# Utility functions
cdef double compute_gini(double* values, uint32_t n) noexcept nogil
cdef double compute_entropy(double* values, uint32_t n) noexcept nogil
cdef double compute_modularity(
    uint32_t* edge_src,
    uint32_t* edge_dst,
    double* edge_weights,
    uint64_t n_edges,
    int32_t* partition,
    uint32_t n_nodes,
) noexcept nogil
cdef double compute_assortativity_categorical(
    uint32_t* edge_src,
    uint32_t* edge_dst,
    double* edge_weights,
    uint64_t n_edges,
    int32_t* categories,
    uint32_t n_nodes,
) noexcept nogil
