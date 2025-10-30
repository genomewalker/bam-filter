# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# cython: nonecheck=False

"""
Multi-metric outlier detection for graph-based filtering.

This module implements various outlier detection algorithms optimized for
identifying anomalous references in BAM filtering based on graph topology metrics.

Algorithms implemented:
1. MAD (Median Absolute Deviation) - robust univariate outlier detection
2. IQR (Interquartile Range) - standard boxplot method
3. Isolation Forest - tree-based multi-metric anomaly detection
4. LOF (Local Outlier Factor) - density-based anomaly detection

The Isolation Forest implementation is lightweight and specifically tailored
for small to medium datasets (hundreds to thousands of references) typical
in metagenomic graph analysis.
"""

from libc.stdlib cimport malloc, free, rand, srand, qsort, RAND_MAX
from libc.math cimport fabs, log, sqrt, exp, isnan, isinf
from libc.string cimport memset, memcpy
from libc.stdint cimport uint32_t, int32_t
from libc.float cimport DBL_MAX
from cython.parallel cimport prange, parallel
cimport openmp

# Import type definitions
from bam_filter.processor_outlier_detection cimport (
    OutlierMethod,
    ReferenceFeatures,
    IsolationTreeNode,
    IsolationForestModel,
    LOFNeighbor,
    LOFScores,
    OutlierDetectionResult,
)

# Constants
DEF MAX_TREE_DEPTH = 20
DEF MIN_SAMPLES_SPLIT = 2
DEF EPSILON = 1e-10

# ============================================================================
# Utility Functions
# ============================================================================

cdef inline double median_nogil(double* arr, uint32_t n) nogil:
    """Calculate median (assumes array will be modified - sorts in place)."""
    if n == 0:
        return 0.0
    if n == 1:
        return arr[0]

    # Simple quickselect would be more efficient, but qsort is simpler for small n
    qsort(arr, n, sizeof(double), compare_doubles)

    if n % 2 == 0:
        return (arr[n // 2 - 1] + arr[n // 2]) / 2.0
    else:
        return arr[n // 2]

cdef int compare_doubles(const void* a, const void* b) noexcept nogil:
    """Comparison function for qsort."""
    cdef double da = (<double*>a)[0]
    cdef double db = (<double*>b)[0]
    if da < db:
        return -1
    elif da > db:
        return 1
    else:
        return 0

cdef inline double calculate_mean(const double* arr, uint32_t n) nogil:
    """Calculate arithmetic mean."""
    if n == 0:
        return 0.0
    cdef double sum_val = 0.0
    cdef uint32_t i
    for i in range(n):
        sum_val += arr[i]
    return sum_val / n

cdef inline double calculate_std(const double* arr, uint32_t n, double mean) nogil:
    """Calculate standard deviation given mean."""
    if n <= 1:
        return 0.0
    cdef double sum_sq = 0.0
    cdef double diff
    cdef uint32_t i
    for i in range(n):
        diff = arr[i] - mean
        sum_sq += diff * diff
    return sqrt(sum_sq / (n - 1))

# ============================================================================
# MAD (Median Absolute Deviation) - Current Default
# ============================================================================

cdef double calculate_mad_threshold_c(
    const double* values,
    uint32_t n_values,
    double z_threshold
) nogil:
    """
    Calculate outlier threshold using MAD (Median Absolute Deviation).

    This is the robust version of z-score outlier detection.
    Modified z-score = 0.6745 * (x - median) / MAD

    Returns:
        threshold: Values below this are considered outliers
    """
    if n_values < 3:
        return -DBL_MAX  # Keep everything if too few samples

    # Allocate working arrays
    cdef double* values_copy = <double*>malloc(n_values * sizeof(double))
    cdef double* deviations = <double*>malloc(n_values * sizeof(double))

    if values_copy == NULL or deviations == NULL:
        free(values_copy)
        free(deviations)
        return -DBL_MAX

    # Copy values (median calculation sorts in place)
    memcpy(values_copy, values, n_values * sizeof(double))

    # Calculate median
    cdef double median_val = median_nogil(values_copy, n_values)

    # Calculate absolute deviations from median
    cdef uint32_t i
    for i in range(n_values):
        deviations[i] = fabs(values[i] - median_val)

    # Calculate MAD (median of absolute deviations)
    cdef double mad = median_nogil(deviations, n_values)

    free(values_copy)
    free(deviations)

    # Handle zero MAD (all values identical)
    if mad < EPSILON:
        return -DBL_MAX  # Keep everything

    # Calculate threshold: median - (z_threshold * MAD / 0.6745)
    # 0.6745 is the constant to make MAD comparable to standard deviation
    cdef double threshold = median_val - (z_threshold * mad / 0.6745)

    return threshold

# ============================================================================
# IQR (Interquartile Range) - Standard Boxplot Method
# ============================================================================

cdef double calculate_iqr_threshold_c(
    const double* values,
    uint32_t n_values,
    double iqr_multiplier
) nogil:
    """
    Calculate outlier threshold using IQR (Interquartile Range).

    Standard boxplot rule: Q1 - iqr_multiplier * IQR

    Returns:
        threshold: Values below this are considered outliers
    """
    if n_values < 4:
        return -DBL_MAX  # Keep everything if too few samples

    # Allocate and copy (sorting required)
    cdef double* values_copy = <double*>malloc(n_values * sizeof(double))
    if values_copy == NULL:
        return -DBL_MAX

    memcpy(values_copy, values, n_values * sizeof(double))
    qsort(values_copy, n_values, sizeof(double), compare_doubles)

    # Calculate Q1 (25th percentile) and Q3 (75th percentile)
    cdef uint32_t q1_idx = n_values // 4
    cdef uint32_t q3_idx = 3 * n_values // 4
    cdef double q1 = values_copy[q1_idx]
    cdef double q3 = values_copy[q3_idx]

    free(values_copy)

    # Calculate IQR and threshold
    cdef double iqr = q3 - q1
    cdef double threshold = q1 - iqr_multiplier * iqr

    return threshold

# ============================================================================
# Isolation Forest - Tree-Based Multi-Metric Anomaly Detection
# ============================================================================

cdef inline double c_value_iforest(uint32_t n) nogil:
    """
    Expected path length for unsuccessful search in BST.
    This is used to normalize the anomaly score.
    """
    if n <= 1:
        return 0.0
    if n == 2:
        return 1.0

    # Harmonic number approximation: H(n-1)
    cdef double harmonic = log(n - 1.0) + 0.5772156649  # Euler-Mascheroni constant
    return 2.0 * harmonic - (2.0 * (n - 1.0) / n)

cdef void build_isolation_tree_recursive(
    const ReferenceFeatures* features,
    const uint32_t* sample_indices,
    uint32_t n_samples,
    uint32_t n_features,
    int32_t depth,
    int32_t max_depth,
    IsolationTreeNode* nodes,
    int32_t* node_count,
    uint32_t* rng_state
) nogil:
    """
    Recursively build an isolation tree.

    The tree isolates anomalies by random recursive partitioning.
    Anomalies are easier to isolate (shorter path length) than normal points.
    """
    cdef int32_t current_node_idx = node_count[0]
    node_count[0] += 1

    # Initialize node
    nodes[current_node_idx].depth = depth
    nodes[current_node_idx].left_child_idx = -1
    nodes[current_node_idx].right_child_idx = -1
    nodes[current_node_idx].split_feature = -1
    nodes[current_node_idx].split_value = 0.0

    # Stopping criteria: leaf node
    if depth >= max_depth or n_samples <= MIN_SAMPLES_SPLIT:
        return

    # Randomly select feature to split on
    cdef int32_t split_feature = xorshift32(rng_state) % n_features

    # Find min/max of selected feature in current sample
    cdef double min_val = DBL_MAX
    cdef double max_val = -DBL_MAX
    cdef uint32_t i, sample_idx
    cdef double feature_val

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        feature_val = get_feature_value(&features[sample_idx], split_feature)
        if feature_val < min_val:
            min_val = feature_val
        if feature_val > max_val:
            max_val = feature_val

    # If all values identical, make this a leaf
    if fabs(max_val - min_val) < EPSILON:
        return

    # Random split point between min and max
    cdef double split_value = min_val + (max_val - min_val) * (xorshift32(rng_state) / <double>0xFFFFFFFF)

    # Partition samples
    cdef uint32_t* left_indices = <uint32_t*>malloc(n_samples * sizeof(uint32_t))
    cdef uint32_t* right_indices = <uint32_t*>malloc(n_samples * sizeof(uint32_t))

    if left_indices == NULL or right_indices == NULL:
        free(left_indices)
        free(right_indices)
        return

    cdef uint32_t n_left = 0
    cdef uint32_t n_right = 0

    for i in range(n_samples):
        sample_idx = sample_indices[i]
        feature_val = get_feature_value(&features[sample_idx], split_feature)
        if feature_val < split_value:
            left_indices[n_left] = sample_idx
            n_left += 1
        else:
            right_indices[n_right] = sample_idx
            n_right += 1

    # If split didn't partition data, make this a leaf
    if n_left == 0 or n_right == 0:
        free(left_indices)
        free(right_indices)
        return

    # Store split information
    nodes[current_node_idx].split_feature = split_feature
    nodes[current_node_idx].split_value = split_value
    nodes[current_node_idx].left_child_idx = node_count[0]

    # Build left subtree
    build_isolation_tree_recursive(
        features, left_indices, n_left, n_features,
        depth + 1, max_depth, nodes, node_count, rng_state
    )

    nodes[current_node_idx].right_child_idx = node_count[0]

    # Build right subtree
    build_isolation_tree_recursive(
        features, right_indices, n_right, n_features,
        depth + 1, max_depth, nodes, node_count, rng_state
    )

    free(left_indices)
    free(right_indices)

cdef inline double get_feature_value(const ReferenceFeatures* feature, int32_t feature_idx) nogil:
    """Extract feature value by index."""
    if feature_idx == 0:
        return feature.clustering_coefficient
    elif feature_idx == 1:
        return feature.node_degree
    elif feature_idx == 2:
        return feature.neighbor_multimap_rate
    elif feature_idx == 3:
        return feature.max_comappings
    elif feature_idx == 4:
        return feature.avg_comappings_per_read
    elif feature_idx == 5:
        return feature.connected_neighbors
    elif feature_idx == 6:
        return feature.multimap_percentage
    elif feature_idx == 7:
        return feature.betweenness_centrality
    else:
        return 0.0

cdef inline uint32_t xorshift32(uint32_t* state) nogil:
    """Fast pseudorandom number generator."""
    cdef uint32_t x = state[0]
    x ^= x << 13
    x ^= x >> 17
    x ^= x << 5
    state[0] = x
    return x

cdef double compute_path_length(
    const ReferenceFeatures* feature,
    const IsolationTreeNode* nodes,
    int32_t root_idx
) nogil:
    """
    Compute path length from root to leaf for a given sample.
    This is the anomaly score for the sample in this tree.
    """
    cdef int32_t current_idx = root_idx
    cdef int32_t depth = 0
    cdef double feature_val

    while current_idx >= 0:
        if nodes[current_idx].split_feature < 0:
            # Leaf node
            return depth

        feature_val = get_feature_value(feature, nodes[current_idx].split_feature)

        if feature_val < nodes[current_idx].split_value:
            current_idx = nodes[current_idx].left_child_idx
        else:
            current_idx = nodes[current_idx].right_child_idx

        depth += 1

        # Safety: prevent infinite loops
        if depth > MAX_TREE_DEPTH * 2:
            break

    return depth

cdef OutlierDetectionResult* detect_outliers_iforest_c(
    const ReferenceFeatures* features,
    uint32_t n_references,
    uint32_t n_trees,
    uint32_t subsample_size,
    double contamination,
    uint32_t random_seed
) nogil:
    """
    Detect outliers using Isolation Forest.

    Parameters:
        features: Array of feature vectors (one per reference)
        n_references: Number of references
        n_trees: Number of trees in ensemble (default: 100)
        subsample_size: Subsample size for building each tree (default: 256 or n_references)
        contamination: Expected proportion of outliers (e.g., 0.1 for 10%)
        random_seed: Random seed for reproducibility

    Returns:
        OutlierDetectionResult with anomaly scores and outlier flags
    """
    if n_references < 3:
        return NULL

    # Allocate result
    cdef OutlierDetectionResult* result = <OutlierDetectionResult*>malloc(sizeof(OutlierDetectionResult))
    if result == NULL:
        return NULL

    result.anomaly_scores = <double*>malloc(n_references * sizeof(double))
    result.is_outlier = <bint*>malloc(n_references * sizeof(bint))

    if result.anomaly_scores == NULL or result.is_outlier == NULL:
        free_outlier_result(result)
        return NULL

    result.n_references = n_references
    memset(result.anomaly_scores, 0, n_references * sizeof(double))
    memset(result.is_outlier, 0, n_references * sizeof(bint))

    # Set subsample size (default: min(256, n_references))
    if subsample_size == 0 or subsample_size > n_references:
        subsample_size = n_references if n_references < 256 else 256

    # Calculate expected path length for normalization
    cdef double c_n = c_value_iforest(subsample_size)

    # Determine number of features
    cdef uint32_t n_features = 8  # All available features (including betweenness)

    # Build trees and compute average path length for each sample
    cdef uint32_t tree_idx, i, j
    cdef uint32_t rng_state = random_seed
    cdef int32_t max_depth = <int32_t>(log(subsample_size) / log(2.0) + 1)
    if max_depth > MAX_TREE_DEPTH:
        max_depth = MAX_TREE_DEPTH

    # Estimate maximum nodes per tree (conservative: 2^(depth+1) - 1)
    cdef int32_t max_nodes_per_tree = (1 << (max_depth + 1))

    # Allocate tree nodes (all trees in one array)
    cdef IsolationTreeNode* all_nodes = <IsolationTreeNode*>malloc(
        n_trees * max_nodes_per_tree * sizeof(IsolationTreeNode)
    )
    cdef int32_t* tree_roots = <int32_t*>malloc(n_trees * sizeof(int32_t))
    cdef uint32_t* sample_indices = <uint32_t*>malloc(subsample_size * sizeof(uint32_t))

    if all_nodes == NULL or tree_roots == NULL or sample_indices == NULL:
        free(all_nodes)
        free(tree_roots)
        free(sample_indices)
        free_outlier_result(result)
        return NULL

    # Build each tree
    cdef int32_t node_offset = 0
    cdef int32_t node_count[1]

    for tree_idx in range(n_trees):
        tree_roots[tree_idx] = node_offset
        node_count[0] = node_offset

        # Random subsample
        for i in range(subsample_size):
            sample_indices[i] = xorshift32(&rng_state) % n_references

        # Build tree
        build_isolation_tree_recursive(
            features, sample_indices, subsample_size, n_features,
            0, max_depth, all_nodes, node_count, &rng_state
        )

        node_offset = node_count[0]

    # Compute anomaly scores (average path length across all trees)
    cdef double path_length
    # Declarations for threshold selection / auto-thresholding
    cdef uint32_t threshold_idx, best_idx, fallback_idx, idx
    cdef double best_gap, gap
    cdef double* sc
    for i in range(n_references):
        path_length = <double>0.0
        for tree_idx in range(n_trees):
            path_length += compute_path_length(&features[i], all_nodes, tree_roots[tree_idx])

        # Average path length
        path_length /= n_trees

        # Normalized anomaly score: s = 2^(-E[h(x)] / c(n))
        # Where E[h(x)] is expected path length and c(n) is average path length for random data
        # Anomalies have shorter paths → higher score (closer to 1)
        # Normal points have longer paths → lower score (closer to 0)
        result.anomaly_scores[i] = exp(-path_length / c_n) if c_n > EPSILON else 0.5

    free(<void*>all_nodes)
    free(<void*>tree_roots)
    free(<void*>sample_indices)

    # Determine threshold based on contamination rate
    # Copy scores for sorting
    cdef double* scores_copy = <double*>malloc(<size_t>n_references * sizeof(double))
    if scores_copy == NULL:
        free_outlier_result(result)
        return <OutlierDetectionResult*>NULL

    memcpy(<void*>scores_copy, <const void*>result.anomaly_scores, <size_t>n_references * sizeof(double))
    qsort(<void*>scores_copy, <size_t>n_references, <size_t>sizeof(double), compare_doubles)

    # Threshold: percentile-based (default behavior)
    threshold_idx = <uint32_t>((1.0 - contamination) * n_references)
    if threshold_idx >= n_references:
        threshold_idx = n_references - 1

    result.threshold = (<double*>scores_copy)[threshold_idx]
    free(<void*>scores_copy)

    # Mark outliers (scores >= threshold)
    result.n_outliers = 0
    cdef double* pr = <double*>result.anomaly_scores
    cdef bint* po = <bint*>result.is_outlier
    for i in range(n_references):
        if pr[0] >= result.threshold:
            po[0] = True
            result.n_outliers += 1
        else:
            po[0] = False
        pr += 1
        po += 1

    return result

# ============================================================================
# Multi-variate Detection Wrapper
# ============================================================================

cdef OutlierDetectionResult* detect_outliers_multivariate_c(
    const ReferenceFeatures* features,
    uint32_t n_references,
    OutlierMethod method,
    double contamination,
    uint32_t random_seed
) nogil:
    """
    Unified interface for multi-metric outlier detection.

    Dispatches to appropriate algorithm based on method parameter.
    """
    if method == OUTLIER_IFOREST:
        return detect_outliers_iforest_c(
            features, n_references,
            100,  # n_trees
            0,    # subsample_size (auto)
            contamination,
            random_seed
        )
    elif method == OUTLIER_LOF:
        return detect_outliers_lof_c(features, n_references, 20, contamination)
    else:
        # Fallback: use MAD on clustering coefficient only
        return NULL

# ============================================================================
# LOF (Local Outlier Factor) - Stub for Future Implementation
# ============================================================================

cdef OutlierDetectionResult* detect_outliers_lof_c(
    const ReferenceFeatures* features,
    uint32_t n_references,
    uint32_t k_neighbors,
    double contamination
) nogil:
    """
    LOF implementation - placeholder for future development.

    LOF detects outliers based on local density deviation.
    Points in low-density regions (relative to neighbors) are outliers.
    """
    # TODO: Implement LOF algorithm
    return NULL

# ============================================================================
# Feature Extraction Helpers
# ============================================================================

cdef void extract_features_from_graph_metrics(
    const double* cc_values,
    const double* degree_values,
    const double* neighbor_multimap,
    const double* max_comappings,
    const double* avg_comappings,
    uint32_t n_references,
    ReferenceFeatures* features_out
) nogil:
    """
    Extract multi-metric feature vectors from graph analysis results.
    """
    cdef uint32_t i
    for i in range(n_references):
        features_out[i].reference_idx = i
        features_out[i].clustering_coefficient = cc_values[i] if cc_values != NULL else 0.0
        features_out[i].node_degree = degree_values[i] if degree_values != NULL else 0.0
        features_out[i].neighbor_multimap_rate = neighbor_multimap[i] if neighbor_multimap != NULL else 0.0
        features_out[i].max_comappings = max_comappings[i] if max_comappings != NULL else 0.0
        features_out[i].avg_comappings_per_read = avg_comappings[i] if avg_comappings != NULL else 0.0
        features_out[i].connected_neighbors = 0.0  # Populated elsewhere
        features_out[i].multimap_percentage = 0.0  # Populated elsewhere

cdef void normalize_features_zscore(
    ReferenceFeatures* features,
    uint32_t n_references,
    double* means_out,
    double* stds_out
) nogil:
    """
    Normalize features to z-scores (mean=0, std=1).

    This is important for Isolation Forest to avoid bias toward features
    with larger scales.
    """
    # Calculate means
    cdef uint32_t i
    cdef double mean_cc = 0.0, mean_degree = 0.0, mean_multimap = 0.0
    cdef double mean_max_co = 0.0, mean_avg_co = 0.0, mean_conn = 0.0, mean_self_mm = 0.0, mean_between = 0.0

    for i in range(n_references):
        mean_cc += features[i].clustering_coefficient
        mean_degree += features[i].node_degree
        mean_multimap += features[i].neighbor_multimap_rate
        mean_max_co += features[i].max_comappings
        mean_avg_co += features[i].avg_comappings_per_read
        mean_conn += features[i].connected_neighbors
        mean_self_mm += features[i].multimap_percentage
        mean_between += features[i].betweenness_centrality

    mean_cc /= n_references
    mean_degree /= n_references
    mean_multimap /= n_references
    mean_max_co /= n_references
    mean_avg_co /= n_references
    mean_conn /= n_references
    mean_self_mm /= n_references
    mean_between /= n_references

    # Calculate standard deviations
    cdef double std_cc = 0.0, std_degree = 0.0, std_multimap = 0.0
    cdef double std_max_co = 0.0, std_avg_co = 0.0, std_conn = 0.0, std_self_mm = 0.0, std_between = 0.0
    cdef double diff

    for i in range(n_references):
        diff = features[i].clustering_coefficient - mean_cc
        std_cc += diff * diff
        diff = features[i].node_degree - mean_degree
        std_degree += diff * diff
        diff = features[i].neighbor_multimap_rate - mean_multimap
        std_multimap += diff * diff
        diff = features[i].max_comappings - mean_max_co
        std_max_co += diff * diff
        diff = features[i].avg_comappings_per_read - mean_avg_co
        std_avg_co += diff * diff
        diff = features[i].connected_neighbors - mean_conn
        std_conn += diff * diff
        diff = features[i].multimap_percentage - mean_self_mm
        std_self_mm += diff * diff
        diff = features[i].betweenness_centrality - mean_between
        std_between += diff * diff

    std_cc = sqrt(std_cc / n_references)
    std_degree = sqrt(std_degree / n_references)
    std_multimap = sqrt(std_multimap / n_references)
    std_max_co = sqrt(std_max_co / n_references)
    std_avg_co = sqrt(std_avg_co / n_references)
    std_conn = sqrt(std_conn / n_references)
    std_self_mm = sqrt(std_self_mm / n_references)
    std_between = sqrt(std_between / n_references)

    # Normalize features
    for i in range(n_references):
        if std_cc > EPSILON:
            features[i].clustering_coefficient = (features[i].clustering_coefficient - mean_cc) / std_cc
        if std_degree > EPSILON:
            features[i].node_degree = (features[i].node_degree - mean_degree) / std_degree
        if std_multimap > EPSILON:
            features[i].neighbor_multimap_rate = (features[i].neighbor_multimap_rate - mean_multimap) / std_multimap
        if std_max_co > EPSILON:
            features[i].max_comappings = (features[i].max_comappings - mean_max_co) / std_max_co
        if std_avg_co > EPSILON:
            features[i].avg_comappings_per_read = (features[i].avg_comappings_per_read - mean_avg_co) / std_avg_co
        if std_conn > EPSILON:
            features[i].connected_neighbors = (features[i].connected_neighbors - mean_conn) / std_conn
        if std_self_mm > EPSILON:
            features[i].multimap_percentage = (features[i].multimap_percentage - mean_self_mm) / std_self_mm
        if std_between > EPSILON:
            features[i].betweenness_centrality = (features[i].betweenness_centrality - mean_between) / std_between

# ============================================================================
# Cleanup
# ============================================================================

cdef void free_outlier_result(OutlierDetectionResult* result) nogil:
    """Free memory allocated for outlier detection result."""
    if result != NULL:
        if result.anomaly_scores != NULL:
            free(result.anomaly_scores)
        if result.is_outlier != NULL:
            free(result.is_outlier)
        free(result)
