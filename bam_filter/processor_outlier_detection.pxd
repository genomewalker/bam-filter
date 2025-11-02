# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

"""
Type definitions for multi-metric outlier detection methods.

This module provides various outlier detection algorithms for identifying
anomalous references in graph-based filtering, including:
- MAD (Median Absolute Deviation) - current default
- IQR (Interquartile Range) - current alternative
- Isolation Forest - tree-based multi-metric detection
- Local Outlier Factor (LOF) - density-based detection
"""

from libc.stdint cimport uint32_t, int32_t

# Outlier detection method enum
cdef enum OutlierMethod:
    OUTLIER_MAD = 0      # Median Absolute Deviation (MAD-based)
    OUTLIER_IQR = 1      # Interquartile Range (standard boxplot)
    OUTLIER_IFOREST = 2  # Isolation Forest (multi-metric)
    OUTLIER_LOF = 3      # Local Outlier Factor (density-based)
    OUTLIER_ZSCORE = 4   # Standard Z-score (parametric)

# Structure to hold multi-metric feature vector for a reference
cdef struct ReferenceFeatures:
    double clustering_coefficient     # Barrat weighted CC
    double node_degree               # Number of connections
    double neighbor_multimap_rate    # Average multimap of neighbors
    double max_comappings            # Maximum co-mapping intensity
    double avg_comappings_per_read   # Average co-mapping intensity
    double connected_neighbors       # Number of connected neighbors
    double multimap_percentage       # Reference's own multimap rate
    double betweenness_centrality    # Betweenness centrality (bridge-ness)
    uint32_t reference_idx           # Reference index (for tracking)

# Isolation Forest tree node structure
cdef struct IsolationTreeNode:
    int32_t left_child_idx    # Index of left child (-1 if leaf)
    int32_t right_child_idx   # Index of right child (-1 if leaf)
    int32_t split_feature     # Feature index to split on (-1 if leaf)
    double split_value        # Threshold value for split
    int32_t depth             # Depth in tree (for path length calculation)

# Isolation Forest model structure
cdef struct IsolationForestModel:
    IsolationTreeNode* trees          # Array of all tree nodes (flattened)
    int32_t* tree_root_indices        # Starting index for each tree
    int32_t n_trees                   # Number of trees in ensemble
    int32_t n_nodes_total             # Total number of nodes across all trees
    int32_t n_features                # Number of features used
    uint32_t subsample_size           # Subsample size used for training
    double contamination              # Expected contamination rate

# Local Outlier Factor (LOF) data structures
cdef struct LOFNeighbor:
    uint32_t neighbor_idx    # Index of neighbor
    double distance          # Distance to neighbor

cdef struct LOFScores:
    double* lrd_scores       # Local reachability density for each point
    double* lof_scores       # Local outlier factor for each point
    uint32_t n_points        # Number of points
    uint32_t k_neighbors     # Number of neighbors used

# Result structure for outlier detection
cdef struct OutlierDetectionResult:
    double* anomaly_scores   # Anomaly score for each reference
    bint* is_outlier         # Boolean flag: 1=outlier, 0=inlier
    double threshold         # Threshold used for classification
    uint32_t n_outliers      # Count of detected outliers
    uint32_t n_references    # Total number of references

# Public C-level API functions
cdef double calculate_mad_threshold_c(
    const double* values,
    uint32_t n_values,
    double z_threshold
) nogil

cdef double calculate_iqr_threshold_c(
    const double* values,
    uint32_t n_values,
    double iqr_multiplier
) nogil

cdef OutlierDetectionResult* detect_outliers_iforest_c(
    const ReferenceFeatures* features,
    uint32_t n_references,
    uint32_t n_trees,
    uint32_t subsample_size,
    double contamination,
    uint32_t random_seed
) nogil

cdef OutlierDetectionResult* detect_outliers_multivariate_c(
    const ReferenceFeatures* features,
    uint32_t n_references,
    OutlierMethod method,
    double contamination,
    uint32_t random_seed
) nogil

cdef void free_outlier_result(OutlierDetectionResult* result) nogil

# Helper functions for feature extraction
cdef void extract_features_from_graph_metrics(
    const double* cc_values,
    const double* degree_values,
    const double* neighbor_multimap,
    const double* max_comappings,
    const double* avg_comappings,
    uint32_t n_references,
    ReferenceFeatures* features_out
) nogil

# Normalization helpers (z-score within dataset)
cdef void normalize_features_zscore(
    ReferenceFeatures* features,
    uint32_t n_references,
    double* means_out,
    double* stds_out
) nogil
