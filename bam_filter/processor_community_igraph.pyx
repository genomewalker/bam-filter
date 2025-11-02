# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
# -*- coding: utf-8 -*-

"""
Community detection (Leiden/union-find) on igraph graphs with clustering-coefficient based filtering.

WORKFLOW:
1. Identify connected components in the pre-built igraph representation.
2. Skip singletons and pairs (auto-keep, no filtering needed).
3. For each multi-node component, run the selected community detection algorithm (Leiden by default).
4. For each resulting community, calculate weighted clustering coefficients (Barrat).
5. Apply adaptive percentile thresholds to detect low-CC outliers.
6. KEEP references with CC >= threshold (cohesive/informative references).
7. REMOVE references with CC < threshold (hubs/stars, likely contamination).

FILTERING LOGIC:
- Uses adaptive percentile thresholds (5th-10th percentile) based on distribution shape.
- Default: 5th percentile (removes only bottom 5%, very conservative).
- High CC = reference's neighbors are connected to each other (clique/triangle pattern) -> KEEP.
- Low CC = reference connects unrelated references (hub/star pattern) -> REMOVE as contamination.
- Threshold represents the minimum "cohesiveness" required to be considered informative.
"""

from libc.stdlib cimport malloc, calloc, free, realloc, qsort
from libc.time cimport clock, clock_t, CLOCKS_PER_SEC
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t
from libc.string cimport memset
from cython.parallel cimport prange, threadid

# OpenMP for controlling thread count in igraph operations
cdef extern from "<omp.h>" nogil:
    void omp_set_num_threads(int num_threads)

from bam_filter.processor_graph cimport ReferenceStats
from bam_filter.processor_graph_ops cimport WeightedGraph, GraphNode
from bam_filter.processor_igraph cimport *

# Import new multi-metric outlier detection module
from bam_filter.processor_outlier_detection cimport (
    ReferenceFeatures,
    OutlierDetectionResult,
    OutlierMethod,
    OUTLIER_MAD,
    OUTLIER_IQR,
    OUTLIER_IFOREST,
    OUTLIER_LOF,
    OUTLIER_ZSCORE,
    detect_outliers_multivariate_c,
    detect_outliers_iforest_c,
    free_outlier_result,
    extract_features_from_graph_metrics,
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil
    int bf_should_log(int level) nogil
    enum:
        BF_LOG_LEVEL_TRACE


# ==============================================================================
# HELPER: Access igraph vector elements
# ==============================================================================

cdef inline igraph_integer_t get_vector_int_element(igraph_vector_int_t *v, igraph_integer_t i) nogil:
    """Access element from igraph_vector_int_t."""
    return (<igraph_integer_t**>v)[0][i]


cdef inline igraph_real_t get_vector_element(igraph_vector_t *v, igraph_integer_t i) nogil:
    """Access element from igraph_vector_t."""
    return (<igraph_real_t**>v)[0][i]


# ==============================================================================
# WEIGHTED CLUSTERING COEFFICIENTS (Barrat's Method)
# ==============================================================================


cdef int calculate_weighted_clustering_coefficients_for_vids(
    igraph_t* graph,
    igraph_vector_t* edge_weights,
    igraph_vector_int_t* vids_vector,
    float* clustering_coefficients,
    uint32_t n_members
) except -1 nogil:
    """
    Calculate Barrat weighted clustering coefficients for a subset of vertices
    specified by an igraph_vector_int_t (vertex ids). This uses igraph_vs_vector
    to create a vertex selector and calls igraph_transitivity_barrat on the full
    graph with that selector. The result is copied into `clustering_coefficients`
    (length n_members).
    """
    cdef int ret
    cdef igraph_vector_t result
    cdef igraph_vs_t vs
    cdef uint32_t i

    ret = igraph_vector_init(&result, <igraph_integer_t>n_members)
    if ret != 0:
        return -1

    ret = igraph_vs_vector(&vs, vids_vector)
    if ret != 0:
        igraph_vector_destroy(&result)
        return -1

    ret = igraph_transitivity_barrat(
        graph,
        &result,
        vs,
        edge_weights,
        IGRAPH_TRANSITIVITY_ZERO
    )

    igraph_vs_destroy(&vs)

    if ret != 0:
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "[ERROR] igraph_transitivity_barrat failed (ret=%d) for %u vertices\n",
            ret,
            n_members,
        )
        igraph_vector_destroy(&result)
        return -1

    for i in range(n_members):
        clustering_coefficients[i] = <float>get_vector_element(&result, <igraph_integer_t>i)

    igraph_vector_destroy(&result)
    return 0


# Helper: process a single community (runs nogil)
cdef int process_single_community_global(
    igraph_t* ig_graph,
    igraph_vector_t* ig_weights,
    uint32_t* member_list,
    uint32_t n_members,
    char* keep_flag,
    float* out_comm_cc,
    float* out_individual_ccs,
    float* out_threshold,
    uint32_t* out_node_degrees,
    uint32_t comm_id,
    bint verbose,
    int outlier_method,  # 0 = MAD (default), 1 = IQR (simplified approach)
    ReferenceStats* ref_stats,  # Per-reference statistics
    float* out_anomaly_scores,  # Output anomaly scores for each reference
    uint32_t* exact_connection_counts,  # Real graph metric: exact number of connected neighbors
    double* co_mapping_averages,         # Real graph metric: average co-mapping intensity per read
    uint64_t* max_co_mappings,           # Real graph metric: maximum co-mapping intensity
    double* neighbor_multimap_avg,       # Real graph metric: average multimap rate of neighbors
    float* betweenness_centrality        # Real graph metric: betweenness centrality from igraph
) nogil:
    """
    Compute Barrat clustering coefficients for a community (given as a member_list
    of global vertex ids), apply broken-stick threshold, and set keep_flag for
    kept vertices. Runs nogil and returns 0 on success.

    IMPORTANT: Extracts the induced subgraph for this community to calculate CC
    on the correct topology (only edges within the community).
    """
    cdef int ret
    cdef igraph_vector_int_t tmp_vec
    cdef uint32_t i
    cdef float* comm_clustering = NULL
    cdef float threshold
    cdef uint32_t kept_in_community
    cdef uint32_t best_i
    cdef float max_c
    cdef igraph_t subgraph
    cdef igraph_vector_t sub_weights

    if verbose and not bf_should_log(BF_LOG_LEVEL_TRACE):
        verbose = False

    if n_members == 0:
        return 0

    if n_members == 1:
        keep_flag[member_list[0]] = 1
        return 0

    if n_members == 2:
        keep_flag[member_list[0]] = 1
        keep_flag[member_list[1]] = 1
        return 0

    # Build vertex list for this community
    ret = igraph_vector_int_init(&tmp_vec, <igraph_integer_t>n_members)
    if ret != IGRAPH_SUCCESS:
        return -1

    for i in range(n_members):
        igraph_vector_int_set(&tmp_vec, <igraph_integer_t>i, <igraph_integer_t>member_list[i])

    # Extract induced subgraph for this community only
    # This ensures CC is calculated on the correct topology (intra-community edges only)
    ret = extract_subgraph_with_weights(ig_graph, ig_weights, &tmp_vec, &subgraph, &sub_weights)
    if ret != 0:
        if verbose:
            bf_nogil_logf_notime(
                b"COMMUNITY",
                "[WARNING] Failed to extract subgraph for community %u (size=%u), keeping all members\n",
                comm_id,
                n_members,
            )
        igraph_vector_int_destroy(&tmp_vec)
        # Auto-keep all members on error
        for i in range(n_members):
            keep_flag[member_list[i]] = 1
        return 0

    comm_clustering = <float*>malloc(n_members * sizeof(float))
    if not comm_clustering:
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        igraph_vector_int_destroy(&tmp_vec)
        return -1

    # Calculate CC on the subgraph (node IDs are now 0 to n_members-1)
    # Build a new vertex selector for the subgraph (all vertices)
    igraph_vector_int_destroy(&tmp_vec)
    ret = igraph_vector_int_init(&tmp_vec, <igraph_integer_t>n_members)
    if ret != IGRAPH_SUCCESS:
        free(comm_clustering)
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        return -1

    for i in range(n_members):
        igraph_vector_int_set(&tmp_vec, <igraph_integer_t>i, <igraph_integer_t>i)

    ret = calculate_weighted_clustering_coefficients_for_vids(&subgraph, &sub_weights, &tmp_vec, comm_clustering, n_members)

    if ret != 0:
        # CC calculation failed - cleanup and return
        igraph_vector_int_destroy(&tmp_vec)
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        if verbose:
            bf_nogil_logf_notime(
                b"COMMUNITY",
                "[WARNING] CC calculation failed for community %u (size=%u), keeping all members\n",
                comm_id,
                n_members,
            )
        for i in range(n_members):
            keep_flag[member_list[i]] = 1
        free(comm_clustering)
        return 0

    # Calculate node degrees in the subgraph (same topology as CC calculation)
    cdef igraph_vector_int_t subgraph_degrees
    cdef igraph_vs_t vs_sub_all
    if out_node_degrees is not NULL:
        ret = igraph_vector_int_init(&subgraph_degrees, 0)
        if ret == IGRAPH_SUCCESS:
            ret = igraph_vs_all(&vs_sub_all)
            if ret == IGRAPH_SUCCESS:
                ret = igraph_degree(&subgraph, &subgraph_degrees, vs_sub_all, IGRAPH_ALL, 0)
                igraph_vs_destroy(&vs_sub_all)
                if ret == IGRAPH_SUCCESS:
                    # Export subgraph degrees (map subgraph IDs back to full graph IDs)
                    for i in range(n_members):
                        out_node_degrees[member_list[i]] = <uint32_t>get_vector_int_element(&subgraph_degrees, <igraph_integer_t>i)
                igraph_vector_int_destroy(&subgraph_degrees)

    # Cleanup subgraph and vertex list
    igraph_vector_int_destroy(&tmp_vec)
    igraph_destroy(&subgraph)
    igraph_vector_destroy(&sub_weights)

    # Calculate statistical outlier threshold for CC values using simplified approach
    # Supports MAD (default) and IQR methods - both univariate on CC only
    # INTERPRETATION: CC >= threshold -> KEEP (cohesive), CC < threshold -> REMOVE (hub/contamination)
    # Note: Complex multi-metric methods removed - 3-tier enhanced filtering handles this better

    cdef uint32_t i_member

    if outlier_method == 1:
        # IQR method (univariate on CC only)
        threshold = calculate_statistical_outlier_threshold_iqr(comm_clustering, n_members, (verbose and bf_should_log(BF_LOG_LEVEL_TRACE)) and n_members > 50)
    else:
        # MAD method (default, outlier_method == 0, univariate on CC only)
        threshold = calculate_statistical_outlier_threshold_mad(comm_clustering, n_members, (verbose and bf_should_log(BF_LOG_LEVEL_TRACE)) and n_members > 50)

    # Store anomaly scores as CC values (for TSV export compatibility)
    if out_anomaly_scores != NULL:
        for i_member in range(n_members):
            ref_idx = member_list[i_member]
            # Use CC as anomaly score (lower CC = higher anomaly)
            out_anomaly_scores[ref_idx] = comm_clustering[i_member]

    # Audit logging: report per-community chosen threshold and community size.
    # printf is safe to call in nogil contexts.
    if verbose and bf_should_log(BF_LOG_LEVEL_TRACE):
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "[COMMUNITY] id=%u size=%u chosen_cc_threshold=%.6f\n",
            comm_id,
            n_members,
            threshold,
        )

    # Apply threshold: KEEP references with CC >= threshold (cohesive/informative)
    # REMOVE references with CC < threshold (hub/star pattern, likely contamination)
    kept_in_community = 0
    for i in range(n_members):
        if comm_clustering[i] >= threshold:
            keep_flag[member_list[i]] = 1
            kept_in_community += 1

    if kept_in_community == 0:
        # Safety: ensure at least one reference kept (highest CC = most cohesive)
        best_i = 0
        max_c = comm_clustering[0]
        for i in range(1, n_members):
            if comm_clustering[i] > max_c:
                max_c = comm_clustering[i]
                best_i = i
        keep_flag[member_list[best_i]] = 1
        kept_in_community = 1

    # Audit logging: report how many were kept/removed in this community
    if verbose and bf_should_log(BF_LOG_LEVEL_TRACE):
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "[COMMUNITY] id=%u size=%u threshold=%.6f kept=%u removed=%u\n",
            comm_id,
            n_members,
            threshold,
            kept_in_community,
            <uint32_t>(n_members - kept_in_community),
        )

    # compute and export average CC for this community
    cdef float avg_cc = 0.0
    for i in range(n_members):
        avg_cc += comm_clustering[i]
    avg_cc /= <float>n_members
    if out_comm_cc is not NULL:
        out_comm_cc[0] = avg_cc

    # Export individual CC values for each member (using full graph IDs)
    if out_individual_ccs is not NULL:
        for i in range(n_members):
            out_individual_ccs[member_list[i]] = comm_clustering[i]

    # Export threshold value for all members of this community
    if out_threshold is not NULL:
        out_threshold[0] = threshold

    free(comm_clustering)
    return 0


# ==============================================================================
# SUBGRAPH EXTRACTION
# ==============================================================================

cdef int extract_subgraph_with_weights(
    igraph_t* full_graph,
    igraph_vector_t* full_weights,
    igraph_vector_int_t* vertex_ids,
    igraph_t* subgraph_out,
    igraph_vector_t* subgraph_weights_out
) except -1 nogil:
    """
    Extract induced subgraph and corresponding edge weights.
    Returns 0 on success, -1 on failure.

    IMPORTANT: Uses igraph_induced_subgraph to create the subgraph, which
    preserves vertex order and creates a proper induced subgraph (all edges
    between the specified vertices). Then extracts edge weights by iterating
    over the subgraph edges and looking up their weights in the original graph.
    """
    cdef int ret
    cdef igraph_integer_t num_vertices = igraph_vector_int_size(vertex_ids)
    cdef igraph_integer_t i, subgraph_ecount
    cdef igraph_vs_t vs
    cdef igraph_integer_t from_id, to_id, eid
    cdef igraph_integer_t global_vid_from, global_vid_to

    if num_vertices == 0:
        return -1

    # Create vertex selector from vertex_ids
    ret = igraph_vs_vector(&vs, vertex_ids)
    if ret != IGRAPH_SUCCESS:
        return -1

    # Create induced subgraph (includes all edges between specified vertices)
    ret = igraph_induced_subgraph(full_graph, subgraph_out, vs, IGRAPH_SUBGRAPH_AUTO)
    igraph_vs_destroy(&vs)

    if ret != IGRAPH_SUCCESS:
        return -1

    # Get number of edges in subgraph
    subgraph_ecount = igraph_ecount(subgraph_out)

    # Initialize weights vector for subgraph
    ret = igraph_vector_init(subgraph_weights_out, subgraph_ecount)
    if ret != IGRAPH_SUCCESS:
        igraph_destroy(subgraph_out)
        return -1

    # For each edge in the subgraph, look up its weight in the original graph
    # The vertices in the subgraph are remapped (0..num_vertices-1), but we
    # can map them back to the original graph using vertex_ids
    for i in range(subgraph_ecount):
        # Get edge endpoints in subgraph (local vertex IDs: 0..num_vertices-1)
        ret = igraph_edge(subgraph_out, <igraph_integer_t>i, &from_id, &to_id)
        if ret != IGRAPH_SUCCESS:
            igraph_vector_destroy(subgraph_weights_out)
            igraph_destroy(subgraph_out)
            return -1

        # Map subgraph vertex IDs back to original graph vertex IDs
        global_vid_from = get_vector_int_element(vertex_ids, from_id)
        global_vid_to = get_vector_int_element(vertex_ids, to_id)

        # Get edge ID in original graph
        ret = igraph_get_eid(full_graph, &eid, global_vid_from, global_vid_to, IGRAPH_UNDIRECTED, False)
        if ret != IGRAPH_SUCCESS:
            # Edge doesn't exist in original graph (shouldn't happen for induced subgraph)
            igraph_vector_destroy(subgraph_weights_out)
            igraph_destroy(subgraph_out)
            return -1

        # Copy weight from original graph to subgraph weights vector
        (<igraph_real_t**>subgraph_weights_out)[0][i] = get_vector_element(full_weights, eid)

    return 0


# ==============================================================================
# BROKEN STICK MODEL
# ==============================================================================

cdef int _float_compare_ascending(const void* a, const void* b) noexcept nogil:
    """Compare function for qsort (ascending order)."""
    cdef float fa = (<float*>a)[0]
    cdef float fb = (<float*>b)[0]
    if fa < fb:
        return -1
    elif fa > fb:
        return 1
    else:
        return 0


cdef float calculate_statistical_outlier_threshold_mad(float* values, uint32_t n, bint verbose) nogil:
    """
    Calculate threshold for identifying low-CC hubs using Modified Z-Score (MAD method).

    APPROACH:
    Uses Median Absolute Deviation (MAD) for robust outlier detection.
    The modified z-score is: 0.6745 * (x - median) / MAD
    Values with |modified_z_score| > 3.5 are considered outliers.

    FILTERING INTERPRETATION:
    - CC >= threshold -> cohesive/well-connected -> informative reference -> KEEP
    - CC < threshold -> hub/star pattern (statistical outlier) -> promiscuous/contamination -> REMOVE

    High CC means the reference's neighbors are also connected to each other (clique/triangle).
    Low CC means the reference connects otherwise unrelated references (hub/star, likely contaminant).

    ADVANTAGES:
    - Only removes TRUE statistical outliers (may remove 0 references if data is clean!)
    - Robust to extreme outliers (unlike standard deviation)
    - Non-parametric (doesn't assume normal distribution)
    - Well-established in statistics literature

    Args:
        values: Array of clustering coefficients (range [0,1])
        n: Number of values
        verbose: Print diagnostic info

    Returns:
        Threshold value (REMOVE if CC < threshold, KEEP if CC >= threshold)
    """
    cdef uint32_t i
    cdef float* sorted_values = <float*>malloc(n * sizeof(float))
    cdef float* abs_deviations = NULL
    cdef float threshold
    cdef float median_cc, mad, modified_z_threshold
    cdef uint32_t median_idx
    cdef uint32_t num_outliers

    if not sorted_values:
        return 0.0  # No filtering if allocation fails

    # Copy and sort values ascending
    for i in range(n):
        sorted_values[i] = values[i]
    qsort(sorted_values, n, sizeof(float), _float_compare_ascending)

    # Calculate median
    median_idx = n / 2
    if n % 2 == 0 and n > 1:
        median_cc = (sorted_values[median_idx - 1] + sorted_values[median_idx]) / 2.0
    else:
        median_cc = sorted_values[median_idx]

    # Calculate absolute deviations from median
    abs_deviations = <float*>malloc(n * sizeof(float))
    if not abs_deviations:
        free(sorted_values)
        return 0.0

    for i in range(n):
        abs_deviations[i] = sorted_values[i] - median_cc
        if abs_deviations[i] < 0:
            abs_deviations[i] = -abs_deviations[i]

    # Sort absolute deviations to find MAD
    qsort(abs_deviations, n, sizeof(float), _float_compare_ascending)

    # MAD is the median of absolute deviations
    if n % 2 == 0 and n > 1:
        mad = (abs_deviations[median_idx - 1] + abs_deviations[median_idx]) / 2.0
    else:
        mad = abs_deviations[median_idx]

    free(abs_deviations)

    # Calculate key percentiles for diagnostics
    cdef float q25, q75, iqr
    cdef uint32_t q25_idx = n / 4
    cdef uint32_t q75_idx = 3 * n / 4
    q25 = sorted_values[q25_idx]
    q75 = sorted_values[q75_idx]
    iqr = q75 - q25

    if verbose:
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "    CC distribution: Q25=%.3f, Median=%.3f, Q75=%.3f, IQR=%.3f, MAD=%.3f\n",
            q25,
            median_cc,
            q75,
            iqr,
            mad,
        )
        bf_nogil_logf_notime(b"COMMUNITY", "    Sample size: %u\n", n)

    # Modified Z-score threshold: 3.5 is standard for outlier detection
    # We want to find LOW outliers (hubs), so we look at the lower tail
    # threshold = median - 3.5 * MAD / 0.6745
    if mad > 0.0:
        modified_z_threshold = 3.5
        threshold = median_cc - (modified_z_threshold * mad / 0.6745)

        # Ensure threshold is in valid range [0, 1]
        if threshold < 0.0:
            threshold = 0.0
        if threshold > median_cc:
            threshold = median_cc

        # Count how many values fall below threshold (outliers)
        num_outliers = 0
        for i in range(n):
            if sorted_values[i] < threshold:
                num_outliers += 1

        if verbose:
            bf_nogil_logf_notime(
                b"COMMUNITY",
                "    [MAD] Outlier threshold: %.6f (MAD-based, modified_z > 3.5)\n",
                threshold,
            )
            bf_nogil_logf_notime(
                b"COMMUNITY",
                "    [MAD] Will REMOVE %u outliers (%.1f%% of community)\n",
                num_outliers,
                100.0 * <float>num_outliers / <float>n,
            )
    else:
        # MAD is 0 (all values identical) - no outliers to remove
        threshold = median_cc
        if verbose:
            bf_nogil_logf_notime(
                b"COMMUNITY",
                "    [MAD] MAD=0 (all values identical), no outliers to remove\n",
            )

    free(sorted_values)
    return threshold


cdef float calculate_statistical_outlier_threshold_iqr(float* values, uint32_t n, bint verbose) nogil:
    """
    Calculate threshold for identifying low-CC hubs using IQR method.

    APPROACH:
    Uses Interquartile Range (IQR) for outlier detection.
    Lower bound = Q1 - 1.5 * IQR
    Upper bound = Q3 + 1.5 * IQR (not used here, we only care about low outliers)

    FILTERING INTERPRETATION:
    - CC >= threshold -> cohesive/well-connected -> informative reference -> KEEP
    - CC < threshold -> hub/star pattern (statistical outlier) -> promiscuous/contamination -> REMOVE

    ADVANTAGES:
    - Only removes TRUE statistical outliers (may remove 0 references if data is clean!)
    - Standard boxplot outlier method
    - Non-parametric (doesn't assume normal distribution)
    - Simple and interpretable

    Args:
        values: Array of clustering coefficients (range [0,1])
        n: Number of values
        verbose: Print diagnostic info

    Returns:
        Threshold value (REMOVE if CC < threshold, KEEP if CC >= threshold)
    """
    cdef uint32_t i
    cdef float* sorted_values = <float*>malloc(n * sizeof(float))
    cdef float threshold
    cdef float q25, q75, iqr, lower_bound
    cdef uint32_t q25_idx, q75_idx, median_idx
    cdef float median_cc
    cdef uint32_t num_outliers

    if not sorted_values:
        return 0.0  # No filtering if allocation fails

    # Copy and sort values ascending
    for i in range(n):
        sorted_values[i] = values[i]
    qsort(sorted_values, n, sizeof(float), _float_compare_ascending)

    # Calculate quartiles
    q25_idx = n / 4
    q75_idx = 3 * n / 4
    median_idx = n / 2

    q25 = sorted_values[q25_idx]
    q75 = sorted_values[q75_idx]
    median_cc = sorted_values[median_idx]
    iqr = q75 - q25

    # Lower outlier bound: Q1 - 1.5 * IQR
    lower_bound = q25 - 1.5 * iqr

    # Ensure threshold is in valid range [0, 1]
    threshold = lower_bound
    if threshold < 0.0:
        threshold = 0.0
    if threshold > 1.0:
        threshold = 1.0

    # Count how many values fall below threshold (outliers)
    num_outliers = 0
    for i in range(n):
        if sorted_values[i] < threshold:
            num_outliers += 1

    if verbose:
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "    CC distribution: Q25=%.3f, Median=%.3f, Q75=%.3f, IQR=%.3f\n",
            q25,
            median_cc,
            q75,
            iqr,
        )
        bf_nogil_logf_notime(b"COMMUNITY", "    Sample size: %u\n", n)
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "    [IQR] Outlier threshold: %.6f (Q1 - 1.5*IQR = %.3f - 1.5*%.3f)\n",
            threshold,
            q25,
            iqr,
        )
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "    [IQR] Will REMOVE %u outliers (%.1f%% of community)\n",
            num_outliers,
            100.0 * <float>num_outliers / <float>n,
        )

    free(sorted_values)
    return threshold


cdef float calculate_broken_stick_threshold(float* values, uint32_t n, bint verbose) nogil:
    """
    Statistical outlier detection using MAD (Median Absolute Deviation).

    Uses Modified Z-Score with MAD for robust outlier detection.
    Only removes TRUE statistical outliers (may remove 0 references if data is clean).

    This function maintains the old name for backward compatibility but uses
    the new MAD-based statistical approach instead of percentiles.
    """
    return calculate_statistical_outlier_threshold_mad(values, n, verbose)


# ==============================================================================
# COMPONENT INFO STRUCTURE
# ==============================================================================

cdef struct ComponentInfo:
    uint32_t component_id
    uint32_t* members
    uint32_t n_members
    uint32_t refs_removed
    double time_total
    double time_extract
    double time_barrat
    double time_detection
    double time_communities


# ==============================================================================
# COMMUNITY RESULTS STRUCTURE
# ==============================================================================

cdef struct CommunityResults:
    uint32_t num_nodes  # number of nodes processed
    char* keep_flag  # 1 = keep, 0 = remove
    uint32_t* community_membership  # community ID for each reference
    float* community_cc_values  # Average CC value for the community each reference belongs to
    float* individual_cc_values  # Individual CC value for each reference (Barrat's method)
    float* cc_threshold_values  # Broken-stick threshold used for each reference's community
    float* anomaly_scores  # Anomaly score for each reference (for multi-metric methods like Isolation Forest)
    uint32_t* node_degree  # Node degree (number of edges) - already in ReferencePattern, kept here for convenience
    uint32_t* num_neighbor_communities  # Number of distinct communities among neighbors (for Tier 1)


# ==============================================================================
# MAIN COMMUNITY DETECTION ENTRY POINT
# ==============================================================================

cdef CommunityResults* community_clustering(WeightedGraph* graph,
                             ReferenceStats* ref_stats,
                             double resolution,
                             int32_t max_iterations,
                             bint verbose,
                             int32_t thread_count,
                             int outlier_method,
                             uint32_t* exact_connection_counts,
                             double* co_mapping_averages,
                             uint64_t* max_co_mappings,
                             double* neighbor_multimap_avg,
                             uint32_t array_size) except NULL nogil:
    """
    Run Community clustering with clustering coefficient filtering.
    
    Returns CommunityResults struct containing:
    - keep_flag: 1 = keep, 0 = remove
    - community_membership: community ID for each reference
    - community_cc_values: CC value for the community each reference belongs to
    
    Workflow:
    1. Get cached igraph from WeightedGraph
    2. Find connected components
    3. Auto-keep all singletons (no processing needed)
    4. For each multi-node component, run Community + CC filtering (PARALLEL)
    5. Return CommunityResults struct
    """
    cdef igraph_t* ig_graph
    cdef igraph_vector_t* ig_weights
    cdef igraph_vector_int_t component_membership, component_sizes
    cdef igraph_integer_t num_components
    cdef char* keep_flag = NULL
    cdef ComponentInfo* components = NULL
    cdef uint32_t** component_members_lists = NULL
    cdef uint32_t* component_counts = NULL
    cdef uint32_t i, k, ref_idx, comm_id, ni
    cdef int ret
    cdef uint32_t idx
    cdef int effective_threads
    cdef uint32_t singletons_count = 0
    cdef uint32_t multi_node_count = 0
    cdef uint32_t* multi_node_indices = NULL
    cdef uint32_t multi_idx
    # Stats variables for multi-node component sizes
    cdef uint32_t min_size, max_size, sidx, cur_size
    cdef double sum_sizes
    # Timing aggregation temporaries (declare at function top for Cython)
    cdef double total_time, total_extract, total_barrat, total_detection, total_communities
    cdef double t_start, t_end
    cdef uint32_t ci
    cdef uint32_t key, key_size
    cdef int j
    # Additional C-level locals hoisted here to satisfy Cython's requirement
    # that all cdef declarations appear before executable statements.
    cdef igraph_vector_int_t global_membership
    cdef igraph_integer_t num_communities
    cdef igraph_vector_int_t community_sizes
    cdef igraph_vector_int_t community_membership_tmp
    cdef igraph_integer_t max_m
    cdef igraph_integer_t tmp_m
    cdef uint32_t num_comms_u
    cdef uint32_t comm_idx
    cdef uint32_t n_members_in_comm
    cdef uint32_t actual_num_comms
    cdef uint32_t* community_member_counts
    cdef uint32_t** community_member_lists
    cdef uint32_t member_pos
    cdef uint32_t large_comm_count
    cdef uint32_t singleton_comm_count
    cdef uint32_t pair_comm_count
    cdef uint32_t* active_communities
    cdef uint32_t ac_pos
    cdef uint32_t comms_to_process
    cdef igraph_vector_int_t vids_vec
    cdef igraph_vector_int_t tmp_vec
    cdef double t_community_total
    cdef igraph_integer_t v, m
    cdef igraph_integer_t gm

    if verbose and not bf_should_log(BF_LOG_LEVEL_TRACE):
        verbose = False
    cdef igraph_t subgraph
    cdef igraph_vector_t sub_weights
    cdef igraph_vector_int_t vids_to_detect
    cdef igraph_vector_int_t membership_sel
    cdef igraph_vs_t vs_sel
    # Variables for num_neighbor_communities calculation
    cdef igraph_vector_int_t neighbors_vec
    cdef igraph_integer_t neighbor_id
    cdef uint32_t neighbor_community
    cdef GraphNode* node  # For accessing cached neighbor lists
    # Auto-keep isolated nodes variables
    cdef uint8_t* is_connected = NULL
    cdef uint32_t auto_kept_count = 0
    cdef uint32_t pruned_isolated_count = 0
    cdef uint32_t node_id
    cdef uint32_t unique_communities[256]
    cdef uint32_t num_unique
    cdef uint32_t ui
    cdef bint found
    # Variables for betweenness normalization
    cdef float max_betweenness
    cdef float normalization_factor
    cdef igraph_integer_t n_selected
    cdef float* per_comm_avg
    cdef uint32_t n_members
    cdef float comm_cc
    cdef uint32_t total_assigned
    # Variables for hybrid bridge scoring
    cdef uint32_t inter_comm_edges
    cdef uint32_t my_comm, neighbor_comm
    cdef float comm_diversity, edge_score
    cdef uint32_t max_communities
    cdef float bridge_score
    # Variables for optimized betweenness calculation (bridge candidates)
    cdef uint32_t bridge_candidate_count
    cdef uint32_t max_bridge_nodes
    cdef igraph_vector_int_t bridge_candidates
    cdef igraph_vs_t vs_bridges
    cdef igraph_vector_t betweenness_vec_bridges
    cdef int effective_threads_btw
    cdef uint32_t chunk_size, tid, chunk_start, chunk_end, chunk_len
    cdef int thread_ret
    cdef igraph_vector_int_t chunk_candidates
    # Thread-local variables for parallel betweenness
    cdef uint32_t local_idx, local_i, local_ref_idx
    cdef igraph_vector_int_t local_candidates
    cdef igraph_vector_t local_betweenness
    cdef igraph_vs_t local_vs
    cdef int local_ret
    # Variables for algorithm selection (LPA vs Community)
    cdef bint use_lpa
    cdef const char* algorithm_name
    # Variables for igraph rebuilding
    cdef igraph_t* ig_graph_rebuilt
    cdef igraph_vector_t* ig_weights_rebuilt
    cdef igraph_vector_int_t edges_vec_rebuild
    cdef igraph_vector_t weights_vec_temp
    cdef uint32_t node_idx_rebuild, neighbor_idx_rebuild, neighbor_id_rebuild
    cdef uint32_t edge_count_rebuild

    # Note: We'll control parallelization manually for betweenness
    # Set OpenMP threads to 1 for igraph to avoid nested parallelism issues
    omp_set_num_threads(1)

    # ALWAYS log start (not conditional on verbose) for debugging
    bf_nogil_logf_notime(b"COMMUNITY", "community_clustering: ENTERED function nodes=%u\n", graph.num_nodes)

    # Start timing
    cdef clock_t t_phase_start = clock()
    cdef clock_t t_step_start, t_step_end
    cdef double elapsed_sec

    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "community_clustering: start resolution=%.3f max_iterations=%d threads=%d\n",
                           resolution, max_iterations, thread_count)

    # Use cached igraph from Phase 5
    # NOTE: graph.nodes has been freed, but the igraph structure is still valid
    bf_nogil_logf_notime(b"COMMUNITY", "Using cached igraph from Phase 5...\n")

    if not graph.igraph_handle:
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: No cached igraph found\n")
        return NULL

    if not graph.weights_handle:
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: No cached weights found\n")
        return NULL

    ig_graph = <igraph_t*>graph.igraph_handle
    ig_weights = <igraph_vector_t*>graph.weights_handle

    bf_nogil_logf_notime(b"COMMUNITY", "Cached igraph loaded successfully\n")

    bf_nogil_logf_notime(b"COMMUNITY", "Skipping graph stats (may hang), proceeding to component finding...\n")
    # SKIP: Graph stat calls hang - indicates potential corruption
    # if verbose:
    #     bf_nogil_logf_notime(b"COMMUNITY", "community_clustering: graph vertices=%ld edges=%ld\n",
    #                        <long>igraph_vcount(ig_graph), <long>igraph_ecount(ig_graph))

    # Find connected components
    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Initializing component_membership vector...\n")
    ret = igraph_vector_int_init(&component_membership, graph.num_nodes)
    if ret != IGRAPH_SUCCESS:
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to initialize component membership vector\n")
        return NULL
    
    ret = igraph_vector_int_init(&component_sizes, 0)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&component_membership)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to initialize component sizes vector\n")
        return NULL

    bf_nogil_logf_notime(b"COMMUNITY", "SKIPPING igraph_clusters (hangs on 101k nodes), running Community/LPA on full graph...\n")

    # WORKAROUND: igraph_clusters() hangs indefinitely on large graphs (101k nodes)
    # Skip component detection and treat entire graph as one component
    # Community/LPA will naturally separate disconnected parts into different communities
    num_components = 1

    # Set all nodes to component 0 (all in one big component)
    bf_nogil_logf_notime(b"COMMUNITY", "Setting all nodes to component 0...\n")
    for ref_idx in range(graph.num_nodes):
        igraph_vector_int_set(&component_membership, ref_idx, 0)

    bf_nogil_logf_notime(b"COMMUNITY", "Setting component size...\n")
    igraph_vector_int_push_back(&component_sizes, graph.num_nodes)
    ret = IGRAPH_SUCCESS

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Component setup complete, proceeding to Community/LPA... (%.2fs)\n", elapsed_sec)

    # OLD CODE (hangs):
    # ret = igraph_clusters(ig_graph, &component_membership, &component_sizes,
    #                      &num_components, IGRAPH_WEAK)

    bf_nogil_logf_notime(b"COMMUNITY", "Checking component detection result...\n")
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to find connected components\n")
        return NULL

    bf_nogil_logf_notime(b"COMMUNITY", "Component detection OK, components=%ld\n", <long>num_components)

    # Allocate CommunityResults struct
    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Allocating CommunityResults struct...\n")
    cdef CommunityResults* results = <CommunityResults*>malloc(sizeof(CommunityResults))
    if not results:
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate CommunityResults\n")
        return NULL
    
    results.num_nodes = graph.num_nodes
    
    # Allocate keep_flag
    bf_nogil_logf_notime(b"COMMUNITY", "Allocating keep_flag array (%u bytes)...\n", graph.num_nodes)
    results.keep_flag = <char*>calloc(graph.num_nodes, sizeof(char))
    if not results.keep_flag:
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate keep_flag\n")
        return NULL
    bf_nogil_logf_notime(b"COMMUNITY", "keep_flag allocated OK\n")

    # Allocate community membership array
    bf_nogil_logf_notime(b"COMMUNITY", "Allocating community_membership array (%u elements)...\n", graph.num_nodes)
    results.community_membership = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.community_membership:
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate community_membership\n")
        return NULL

    # Allocate component membership array
    results.component_membership = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.component_membership:
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate component_membership\n")
        return NULL

    # Allocate node degree array
    results.node_degree = <uint32_t*>malloc(graph.num_nodes * sizeof(uint32_t))
    if not results.node_degree:
        free(results.node_degree)
        free(results.component_membership)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate node_degree\n")
        return NULL

    # Allocate community CC values array
    results.community_cc_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.community_cc_values:
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate community_cc_values\n")
        return NULL

    # Allocate individual CC values array
    results.individual_cc_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.individual_cc_values:
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate individual_cc_values\n")
        return NULL

    # Allocate CC threshold values array
    results.cc_threshold_values = <float*>malloc(graph.num_nodes * sizeof(float))
    if not results.cc_threshold_values:
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate cc_threshold_values\n")
        return NULL

    # Allocate anomaly scores array (for multi-metric methods like Isolation Forest)
    results.anomaly_scores = <float*>calloc(graph.num_nodes, sizeof(float))
    if not results.anomaly_scores:
        free(results.cc_threshold_values)
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate anomaly_scores\n")
        return NULL

    # Allocate betweenness centrality array
    results.betweenness_centrality = <float*>calloc(graph.num_nodes, sizeof(float))
    if not results.betweenness_centrality:
        free(results.anomaly_scores)
        free(results.cc_threshold_values)
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate betweenness_centrality\n")
        return NULL

    # Allocate num_neighbor_communities array
    results.num_neighbor_communities = <uint32_t*>calloc(graph.num_nodes, sizeof(uint32_t))
    if not results.num_neighbor_communities:
        free(results.betweenness_centrality)
        free(results.anomaly_scores)
        free(results.cc_threshold_values)
        free(results.individual_cc_values)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to allocate num_neighbor_communities\n")
        return NULL

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Memory allocation complete (%.2fs)\n", elapsed_sec)

    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Setting component membership (all nodes in component 0)...\n")
    # Since we skipped igraph_clusters() and filled all nodes with component 0,
    # we can just set all values to 0 directly (much faster than get_vector_int_element)
    for ref_idx in range(graph.num_nodes):
        results.component_membership[ref_idx] = 0
    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Component membership set (%u nodes, %.2fs)\n", graph.num_nodes, elapsed_sec)

    # Initialize node degrees to 0 (will be filled during community processing)
    # Node degrees will be calculated from community subgraphs, not the full graph
    for ref_idx in range(graph.num_nodes):
        results.node_degree[ref_idx] = 0

    # DEFERRED: Betweenness calculation moved to AFTER Community clustering
    # We only need betweenness for bridge detection (nodes connecting multiple communities)
    # Calculating it on the full 100k+ node graph is O(n^3) and can take hours
    # Instead, we'll calculate it after Community and only for bridge candidates
    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "Deferring betweenness calculation until after Community clustering...\n")

    # Initialize all betweenness values to 0 for now
    for ref_idx in range(graph.num_nodes):
        results.betweenness_centrality[ref_idx] = 0.0

    # Initialize num_neighbor_communities to 0 (will be calculated AFTER community detection)
    for ref_idx in range(graph.num_nodes):
        results.num_neighbor_communities[ref_idx] = 0

    # Count members per component and identify singletons
    component_counts = <uint32_t*>calloc(num_components, sizeof(uint32_t))
    if not component_counts:
        free(results.node_degree)
        free(results.component_membership)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    t_step_start = clock()
    # Since all nodes are in component 0 (we skipped igraph_clusters),
    # component_counts[0] = graph.num_nodes, and there are no singletons
    component_counts[0] = graph.num_nodes
    singletons_count = 0
    multi_node_count = 1  # Only component 0, which has all nodes
    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Component counting: all %u nodes in single component 0 (%.2fs)\n", graph.num_nodes, elapsed_sec)

    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "  Singletons: %u (auto-keep, no processing)\n", singletons_count)
        bf_nogil_logf_notime(b"COMMUNITY", "  Multi-node components: %u (will process)\n", multi_node_count)

    # Auto-keep isolated nodes that were NEVER in the graph (< min_read_count)
    # Nodes that became isolated after edge pruning are likely contamination -> filter them
    # Variables declared at function scope (lines 842-844)
    if graph.connected_node_ids and graph.n_connected_nodes > 0:
        t_step_start = clock()
        # Build a fast lookup set for connected nodes
        is_connected = <uint8_t*>calloc(graph.num_nodes, sizeof(uint8_t))
        if is_connected:
            for i in range(graph.n_connected_nodes):
                node_id = graph.connected_node_ids[i]
                if node_id < graph.num_nodes:
                    is_connected[node_id] = 1

            # Handle isolated nodes:
            # 1. Nodes with < min_read_count: Auto-keep (never entered graph)
            # 2. Nodes isolated by edge pruning: Keep for now, will be evaluated by taxonomy later
            #    (important for aDNA: low coverage is expected, weak edges may be legitimate)
            auto_kept_count = 0
            pruned_isolated_count = 0
            for ref_idx in range(graph.num_nodes):
                if not is_connected[ref_idx]:
                    # Check if this node was filtered by min_read_count or by edge pruning
                    if ref_stats and ref_stats[ref_idx].total_reads < <uint32_t>graph.tsv_min_read_count:
                        # Never entered graph -> auto-keep
                        results.keep_flag[ref_idx] = 1
                        auto_kept_count += 1
                    else:
                        # Became isolated after edge pruning
                        # For aDNA, keep these nodes (low coverage expected) but mark for taxonomy evaluation
                        # The graph structure has original_degree saved to identify these later
                        results.keep_flag[ref_idx] = 1  # Conservative: keep by default
                        pruned_isolated_count += 1

            free(is_connected)
            t_step_end = clock()
            elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
            bf_nogil_logf_notime(b"COMMUNITY", "Auto-kept %u originally isolated nodes (< min_read_count), %u pruned isolated nodes (low coverage, kept for taxonomy eval, %.2fs)\n",
                                auto_kept_count, pruned_isolated_count, elapsed_sec)

    # If no multi-node components, we're done
    if multi_node_count == 0:
        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "No multi-node components to process!\n")
        free(component_counts)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return results
    
    # Build list of multi-node component indices
    multi_node_indices = <uint32_t*>malloc(multi_node_count * sizeof(uint32_t))
    if not multi_node_indices:
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    multi_idx = 0
    for k in range(num_components):
        if component_counts[k] > 1:
            multi_node_indices[multi_idx] = k
            multi_idx += 1

    # Compute and print basic size statistics for multi-node components
    # (min, max, mean) to help diagnostics when verbose is enabled.
    if verbose and multi_node_count > 0:
        # Use variables declared at function top (Cython requires cdef at top-level)
        sum_sizes = 0.0

        min_size = component_counts[multi_node_indices[0]]
        max_size = min_size

        for sidx in range(multi_node_count):
            cur_size = component_counts[multi_node_indices[sidx]]
            if cur_size < min_size:
                min_size = cur_size
            if cur_size > max_size:
                max_size = cur_size
            sum_sizes += <double>cur_size

        bf_nogil_logf_notime(
            b"COMMUNITY",
            "  Multi-node component sizes: min=%u, max=%u, mean=%.2f\n",
            min_size,
            max_size,
            sum_sizes / <double>multi_node_count,
        )
    
    # New approach: run Community once on the full graph, then process communities in parallel
    # Build membership vector for full graph using igraph_community_leiden

    # Allocate global membership vector (size = number of nodes) and set to -1
    # We'll run Community only on the non-singleton vertices (faster) and map
    # the resulting membership back into this full-size vector so downstream
    # code can remain unchanged.
    ret = igraph_vector_int_init(&global_membership, graph.num_nodes)
    if ret != IGRAPH_SUCCESS:
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    # initialize all entries to -1 (meaning: not assigned / singleton)
    for i in range(graph.num_nodes):
        igraph_vector_int_set(&global_membership, <igraph_integer_t>i, <igraph_integer_t>-1)

    # Build vertex list excluding singletons (component_counts[comm_id] > 1)

    ret = igraph_vector_int_init(&vids_to_detect, 0)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    t_step_start = clock()
    # OPTIMIZATION: Use pre-computed connected nodes list (only nodes with degree > 0)
    # This avoids processing ~88k isolated nodes that would become singleton communities
    if graph.connected_node_ids and graph.n_connected_nodes > 0:
        # Use connected nodes list from Phase 5
        bf_nogil_logf_notime(b"COMMUNITY", "Using pre-computed connected nodes list (%u nodes)\n", graph.n_connected_nodes)
        for i in range(graph.n_connected_nodes):
            node_id = graph.connected_node_ids[i]
            # Since all nodes are in component 0 and component_counts[0] = graph.num_nodes > 2,
            # all connected nodes pass the filter
            igraph_vector_int_push_back(&vids_to_detect, <igraph_integer_t>node_id)
    else:
        # Fallback: Select vertices for Community: exclude singletons and pairs (component size <= 2)
        bf_nogil_logf_notime(b"COMMUNITY", "No connected nodes list; selecting all nodes with component size > 2\n")
        for i in range(graph.num_nodes):
            comm_id = get_vector_int_element(&component_membership, i)
            if comm_id >= 0 and comm_id < num_components and component_counts[comm_id] > 2:
                igraph_vector_int_push_back(&vids_to_detect, <igraph_integer_t>i)

    n_selected = igraph_vector_int_size(&vids_to_detect)
    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC

    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "community_clustering: selected_vertices=%ld excluded_isolated=%u (%.2fs)\n",
                           <long>n_selected, graph.num_nodes - n_selected, elapsed_sec)

    if n_selected == 0:
        # Nothing to run Community on (all singletons); set num_communities=0
        num_communities = 0
        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "No non-singleton vertices - skipping Community\n")
        igraph_vector_int_destroy(&vids_to_detect)
    else:
        # Prepare membership vector for selected vertices
        ret = igraph_vector_int_init(&membership_sel, n_selected)
        if ret != IGRAPH_SUCCESS:
            igraph_vector_int_destroy(&vids_to_detect)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        t_step_start = clock()
        ret = extract_subgraph_with_weights(ig_graph, ig_weights, &vids_to_detect,
                                            &subgraph, &sub_weights)
        t_step_end = clock()
        elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
        bf_nogil_logf_notime(b"COMMUNITY", "Subgraph extraction complete (%.2fs)\n", elapsed_sec)

        if ret != 0:
            if verbose:
                bf_nogil_logf_notime(b"COMMUNITY", "ERROR: Failed to extract subgraph for Community subset\n")
            igraph_vector_int_destroy(&membership_sel)
            igraph_vector_int_destroy(&vids_to_detect)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        # Automatic algorithm selection based on graph size:
        # - Small graphs (<10k nodes): Use Community (better quality)
        # - Large graphs (≥10k nodes): Use Label Propagation Algorithm (LPA) (much faster, linear time)
        #
        # Both algorithms produce community assignments, which is all we need for:
        # - Tier 2 (community coherence checking)
        # - Tier 3 (filtering decision matrix)
        #
        # LPA is 10-100x faster on large graphs with minimal quality difference for filtering purposes.

        use_lpa = (graph.num_nodes >= 10000)

        if use_lpa:
            algorithm_name = b"LPA"
        else:
            algorithm_name = b"Community"

        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "community_detection: algorithm=%s nodes=%u\n",
                               algorithm_name, graph.num_nodes)

        # time community detection on the extracted subgraph
        t_start = clock()
        bf_nogil_logf_notime(b"COMMUNITY", "Starting community detection (%s)...\n", algorithm_name)

        if use_lpa:
            # Label Propagation Algorithm (LPA) - O(m) linear time
            # Fast and scalable for large graphs
            ret = igraph_community_label_propagation(&subgraph, &membership_sel,
                                                     IGRAPH_ALL,  # mode (edge direction)
                                                     &sub_weights,  # weights
                                                     NULL,  # initial membership (NULL = random)
                                                     NULL)  # fixed vertices (NULL = none)

            # LPA doesn't return num_communities, so we compute it from max membership
            if ret == IGRAPH_SUCCESS:
                num_communities = 0
                for i in range(n_selected):
                    tmp_m = get_vector_int_element(&membership_sel, i)
                    if tmp_m >= num_communities:
                        num_communities = tmp_m + 1
        else:
            # Community algorithm - O(m log n) with optimization
            # Better quality but slower
            ret = igraph_community_leiden(&subgraph, &sub_weights, NULL,
                                          <igraph_real_t>resolution, 0.01, False,
                                          <igraph_integer_t>max_iterations,
                                          &membership_sel, &num_communities, NULL)

        t_end = clock()
        t_community_total = (t_end - t_start) / <double>CLOCKS_PER_SEC
        bf_nogil_logf_notime(b"COMMUNITY", "Community detection complete: %s found %ld communities (%.2fs)\n",
                           algorithm_name, <long>num_communities, t_community_total)

        # cleanup subgraph resources
        t_step_start = clock()
        igraph_destroy(&subgraph)
        igraph_vector_destroy(&sub_weights)
        t_step_end = clock()
        elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
        bf_nogil_logf_notime(b"COMMUNITY", "Subgraph cleanup complete (%.2fs)\n", elapsed_sec)

        if ret != IGRAPH_SUCCESS:
            if verbose:
                bf_nogil_logf_notime(b"COMMUNITY", "ERROR: %s community detection failed on extracted subgraph\n",
                                   algorithm_name)
            igraph_vector_int_destroy(&membership_sel)
            igraph_vector_int_destroy(&vids_to_detect)
            igraph_vector_int_destroy(&global_membership)
            free(component_counts)
            free(results.community_cc_values)
            free(results.community_membership)
            free(results.keep_flag)
            free(results)
            igraph_vector_int_destroy(&component_membership)
            igraph_vector_int_destroy(&component_sizes)
            return NULL

        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "community_detection: algorithm=%s communities=%ld duration=%.3fs\n",
                               algorithm_name, <long>num_communities, t_community_total)

        # Map membership_sel (length n_selected) back to the full-size global_membership
        t_step_start = clock()
        for i in range(<igraph_integer_t>n_selected):
            v = get_vector_int_element(&vids_to_detect, i)
            m = get_vector_int_element(&membership_sel, i)
            igraph_vector_int_set(&global_membership, v, m)
        t_step_end = clock()
        elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
        bf_nogil_logf_notime(b"COMMUNITY", "Membership mapping complete (%ld vertices, %.2fs)\n", <long>n_selected, elapsed_sec)

        # Defensive check: recompute actual number of communities from membership_sel
        # Some igraph builds have returned unexpected num_communities values in the
        # wild (larger than the number of vertices selected). Recompute the
        # actual max membership id and adjust num_communities accordingly so
        # subsequent allocations use a safe, correct size.
        max_m = -1
        for i in range(<igraph_integer_t>n_selected):
            tmp_m = get_vector_int_element(&membership_sel, i)
            if tmp_m > max_m:
                max_m = tmp_m

        actual_num_comms = 0
        if max_m >= 0:
            actual_num_comms = <uint32_t>(max_m + 1)
        else:
            actual_num_comms = 0

        if <uint32_t>num_communities != actual_num_comms:
            if verbose:
                bf_nogil_logf_notime(
                    b"COMMUNITY",
                    "[LEIDEN] Adjusting num_communities: reported=%ld recomputed=%u\n",
                    <long>num_communities,
                    actual_num_comms,
                )
            # Use the recomputed value for allocations and downstream logic
            num_communities = <igraph_integer_t>actual_num_comms

        igraph_vector_int_destroy(&membership_sel)
        igraph_vector_int_destroy(&vids_to_detect)

    # Prepare lists for each community (only for communities with size > 0)
    # First count community sizes
    t_step_start = clock()
    community_member_counts = <uint32_t*>calloc(<size_t>num_communities, sizeof(uint32_t))
    if not community_member_counts:
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    bf_nogil_logf_notime(b"COMMUNITY", "Counting community members (%u nodes, %ld communities)...\n", graph.num_nodes, <long>num_communities)
    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&global_membership, ref_idx)
        if comm_id >= 0 and comm_id < num_communities:
            community_member_counts[comm_id] += 1
    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Community member counting complete (%.2fs)\n", elapsed_sec)

    # Diagnostic: report how many vertices were selected for Community
    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "[LEIDEN] Vertices selected for Community: %ld\n", <long>n_selected)

    # Allocate member lists
    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Allocating community member lists (%ld communities)...\n", <long>num_communities)
    community_member_lists = <uint32_t**>malloc(num_communities * sizeof(uint32_t*))
    if not community_member_lists:
        free(community_member_counts)
        igraph_vector_int_destroy(&global_membership)
        free(component_counts)
        free(results.community_cc_values)
        free(results.community_membership)
        free(results.keep_flag)
        free(results)
        igraph_vector_int_destroy(&component_membership)
        igraph_vector_int_destroy(&component_sizes)
        return NULL

    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] > 0:
            community_member_lists[comm_idx] = <uint32_t*>malloc(community_member_counts[comm_idx] * sizeof(uint32_t))
            if not community_member_lists[comm_idx]:
                # free previous
                for i in range(comm_idx):
                    if community_member_lists[i]: free(community_member_lists[i])
                free(community_member_lists)
                free(community_member_counts)
                igraph_vector_int_destroy(&global_membership)
                free(component_counts)
                free(results.community_cc_values)
                free(results.community_membership)
                free(results.keep_flag)
                free(results)
                igraph_vector_int_destroy(&component_membership)
                igraph_vector_int_destroy(&component_sizes)
                return NULL
        else:
            community_member_lists[comm_idx] = NULL
    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Community member list allocation complete (%.2fs)\n", elapsed_sec)

    # Fill member lists
    # reset counts to use as cursor
    for comm_idx in range(num_communities):
        community_member_counts[comm_idx] = 0

    for ref_idx in range(graph.num_nodes):
        comm_id = get_vector_int_element(&global_membership, ref_idx)
        if comm_id >= 0 and comm_id < num_communities:
            member_pos = community_member_counts[comm_id]
            community_member_lists[comm_id][member_pos] = ref_idx
            community_member_counts[comm_id] += 1
    
    # Initialize CC-related arrays to 0.0
    for ref_idx in range(graph.num_nodes):
        results.community_cc_values[ref_idx] = 0.0
        results.individual_cc_values[ref_idx] = 0.0
        results.cc_threshold_values[ref_idx] = 0.0
    # Build a compact list of communities that need full processing (size > 2)
    large_comm_count = 0
    singleton_comm_count = 0
    pair_comm_count = 0
    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] == 0:
            continue
        elif community_member_counts[comm_idx] == 1:
            singleton_comm_count += 1
        elif community_member_counts[comm_idx] == 2:
            pair_comm_count += 1
        else:
            large_comm_count += 1

    if verbose:
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "Community sizes: singletons=%u, pairs=%u, larger=%u (total communities=%ld)\n",
            singleton_comm_count,
            pair_comm_count,
            large_comm_count,
            <long>num_communities,
        )

    # Diagnostic: sum of all community_member_counts (should equal number of
    # nodes assigned to communities by Community; singletons/pairs may not be
    # included in the Community selection, so this is a sanity check).
    if verbose:
        total_assigned = 0
        for comm_idx in range(num_communities):
            total_assigned += community_member_counts[comm_idx]
        bf_nogil_logf_notime(b"COMMUNITY", "[LEIDEN] Total assigned to communities: %u\n", total_assigned)

    # Allocate per-community buffers (initialized to 0)
    cdef float* per_comm_thresholds = NULL
    if num_communities > 0:
        per_comm_avg = <float*>calloc(num_communities, sizeof(float))
        per_comm_thresholds = <float*>calloc(num_communities, sizeof(float))

    # Process ALL communities in parallel (process_single_community_global handles small ones efficiently)
    effective_threads = thread_count if thread_count > 0 else 1
    if effective_threads > <int>num_communities:
        effective_threads = <int>num_communities

    if verbose:
        bf_nogil_logf_notime(
            b"COMMUNITY",
            "Processing %u communities with %d threads...\n",
            <uint32_t>num_communities,
            effective_threads,
        )

    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Starting community processing (%ld communities, %d threads)...\n", <long>num_communities, effective_threads)

    # Process all communities (sizes 1, 2, and >2)
    # The process_single_community_global function handles all sizes efficiently
    for comm_idx in prange(num_communities, nogil=True, schedule='dynamic', num_threads=effective_threads):
        if community_member_counts[comm_idx] == 0:
            continue
        process_single_community_global(ig_graph, ig_weights,
                                        community_member_lists[comm_idx],
                                        community_member_counts[comm_idx],
                                        results.keep_flag,
                                        &per_comm_avg[comm_idx],
                                        results.individual_cc_values,
                                        &per_comm_thresholds[comm_idx],
                                        results.node_degree,
                                        comm_idx,
                                        verbose,
                                        outlier_method,
                                        ref_stats,
                                        results.anomaly_scores,
                                        exact_connection_counts,
                                        co_mapping_averages,
                                        max_co_mappings,
                                        neighbor_multimap_avg,
                                        results.betweenness_centrality)

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Community processing complete (%ld communities, %.2fs)\n", <long>num_communities, elapsed_sec)

    # Populate results.community_cc_values and cc_threshold_values from per_comm buffers
    t_step_start = clock()
    for comm_idx in range(num_communities):
        if community_member_counts[comm_idx] == 0:
            continue
        n_members = community_member_counts[comm_idx]
        if n_members <= 2:
            comm_cc = 0.0
        else:
            comm_cc = per_comm_avg[comm_idx]
        for i in range(n_members):
            ref_idx = community_member_lists[comm_idx][i]
            results.community_cc_values[ref_idx] = comm_cc
            # Store the threshold used for this community
            results.cc_threshold_values[ref_idx] = per_comm_thresholds[comm_idx]

    # Map global membership vector into results.community_membership for all nodes.
    # Nodes that were not part of the Community run (singletons) will have -1 in
    # global_membership; represent those as UINT32_MAX (0xFFFFFFFF) in the
    # public results so downstream TSV/BAM writers can encode 'no community'.
    for ref_idx in range(graph.num_nodes):
        gm = get_vector_int_element(&global_membership, <igraph_integer_t>ref_idx)
        if gm >= 0:
            results.community_membership[ref_idx] = <uint32_t>gm
        else:
            results.community_membership[ref_idx] = <uint32_t>0xFFFFFFFF

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Populating community membership (%u nodes, %.2fs)\n", graph.num_nodes, elapsed_sec)

    # NOW calculate num_neighbor_communities for each node (AFTER communities are assigned)
    # This is used for bridge detection (Tier 1 filtering)
    # OPTIMIZATION: Only calculate for connected nodes (isolated nodes have 0 neighbors)
    t_step_start = clock()
    cdef uint32_t nodes_to_check = graph.n_connected_nodes if graph.connected_node_ids else graph.num_nodes
    bf_nogil_logf_notime(b"COMMUNITY", "Calculating neighbor community counts for %u connected nodes...\n", nodes_to_check)

    # Use cached neighbor lists from graph structure (much faster!)
    if graph.connected_node_ids and graph.n_connected_nodes > 0 and graph.nodes:
        for i in range(graph.n_connected_nodes):
            ref_idx = graph.connected_node_ids[i]
            # Use pre-computed neighbors from WeightedGraph
            node = &graph.nodes[ref_idx]
            num_unique = 0

            # Iterate through cached neighbors
            for ni in range(node.degree):
                neighbor_id = node.neighbors[ni]
                if neighbor_id < graph.num_nodes:
                    neighbor_community = results.community_membership[neighbor_id]

                    # Check if this community is already in our unique list
                    found = False
                    for ui in range(num_unique):
                        if unique_communities[ui] == neighbor_community:
                            found = True
                            break

                    # If not found, add it
                    if not found and num_unique < 256:
                        unique_communities[num_unique] = neighbor_community
                        num_unique += 1

            results.num_neighbor_communities[ref_idx] = num_unique
        # Isolated nodes already have num_neighbor_communities = 0 (from calloc)
    else:
        # Fallback: process all nodes (old behavior)
        for ref_idx in range(graph.num_nodes):
            # Get neighbors of this node
            ret = igraph_vector_int_init(&neighbors_vec, 0)
            if ret == IGRAPH_SUCCESS:
                ret = igraph_neighbors(ig_graph, &neighbors_vec, ref_idx, IGRAPH_ALL)
                if ret == IGRAPH_SUCCESS:
                    # Count unique communities among neighbors
                    num_unique = 0
                    for ni in range(igraph_vector_int_size(&neighbors_vec)):
                        neighbor_id = get_vector_int_element(&neighbors_vec, ni)
                        if neighbor_id >= 0 and neighbor_id < <igraph_integer_t>graph.num_nodes:
                            neighbor_community = results.community_membership[neighbor_id]

                            # Check if this community is already in our unique list
                            found = False
                            for ui in range(num_unique):
                                if unique_communities[ui] == neighbor_community:
                                    found = True
                                    break

                            # If not found, add it
                            if not found and num_unique < 256:
                                unique_communities[num_unique] = neighbor_community
                                num_unique += 1

                    results.num_neighbor_communities[ref_idx] = num_unique
                else:
                    results.num_neighbor_communities[ref_idx] = 0

                igraph_vector_int_destroy(&neighbors_vec)
            else:
                results.num_neighbor_communities[ref_idx] = 0

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Neighbor community counts calculated (%u nodes, %.2fs)\n", graph.num_nodes, elapsed_sec)

    # Cleanup community lists
    for comm_idx in range(num_communities):
        if community_member_lists[comm_idx]:
            free(community_member_lists[comm_idx])
    free(community_member_lists)
    free(community_member_counts)
    if per_comm_avg:
        free(per_comm_avg)
    if per_comm_thresholds:
        free(per_comm_thresholds)
    free(component_counts)
    igraph_vector_int_destroy(&component_membership)
    igraph_vector_int_destroy(&component_sizes)
    
    if verbose:
        bf_nogil_logf_notime(b"COMMUNITY", "Community clustering with CC filtering completed\n")
        bf_nogil_logf_notime(b"COMMUNITY", "==========================================================\n\n")

    # NOW calculate betweenness, but ONLY for bridge candidates (nodes with neighbors in multiple communities)
    # This is vastly faster than global betweenness on 100k+ nodes
    t_step_start = clock()
    bf_nogil_logf_notime(b"COMMUNITY", "Identifying bridge candidates...\n")
    bridge_candidate_count = 0
    max_bridge_nodes = 10000  # Safety limit

    # Identify bridge candidates: nodes with num_neighbor_communities > 1
    # OPTIMIZATION: Only check connected nodes (isolated nodes have num_neighbor_communities=0)
    ret = igraph_vector_int_init(&bridge_candidates, 0)
    if ret == IGRAPH_SUCCESS:
        if graph.connected_node_ids and graph.n_connected_nodes > 0:
            # Use connected nodes list (much faster - skip isolated nodes)
            for i in range(graph.n_connected_nodes):
                ref_idx = graph.connected_node_ids[i]
                if ref_idx < graph.num_nodes and results.num_neighbor_communities[ref_idx] > 1:
                    igraph_vector_int_push_back(&bridge_candidates, <igraph_integer_t>ref_idx)
                    bridge_candidate_count += 1
        else:
            # Fallback: check all nodes
            for ref_idx in range(graph.num_nodes):
                if results.num_neighbor_communities[ref_idx] > 1:
                    igraph_vector_int_push_back(&bridge_candidates, <igraph_integer_t>ref_idx)
                    bridge_candidate_count += 1

        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: bridge_candidates=%u total_nodes=%u pct=%.1f%%\n",
                               bridge_candidate_count, graph.num_nodes,
                               (100.0 * bridge_candidate_count) / graph.num_nodes)

        # HYBRID APPROACH: Calculate betweenness only for bridge candidates
        # Step 1: Pre-filter using num_neighbor_communities >= 2 (already done above)
        # Step 2: Calculate betweenness only for these candidates (much faster!)
        #
        # Instead of calculating on full 101k graph (hours), we calculate on
        # a subgraph induced by bridge candidates and their immediate neighbors
        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "HYBRID: Using num_neighbor_communities pre-filter + targeted betweenness\n")

        if bridge_candidate_count > 0 and bridge_candidate_count < graph.num_nodes and bridge_candidate_count <= max_bridge_nodes:
            # HYBRID BETWEENNESS: Fast approximation for bridge candidates
            # Instead of exact betweenness on full graph (O(V*E)), use:
            # 1. Bridge score = num_neighbor_communities (already calculated)
            # 2. Normalize to [0, 1] range for compatibility
            # 3. Optional: Calculate exact betweenness if candidate count is very small

            if verbose:
                bf_nogil_logf_notime(b"COMMUNITY", "Computing bridge scores for %u candidates...\n", bridge_candidate_count)

            # Compute normalization factor based on max possible num_neighbor_communities
            max_communities = num_communities if num_communities > 0 else 1

            # For each bridge candidate, compute refined bridge score
            # Combines: num_neighbor_communities + inter-community edge count

            for i in range(bridge_candidate_count):
                ref_idx = <uint32_t>get_vector_int_element(&bridge_candidates, i)
                if ref_idx < graph.num_nodes:
                    # Get this node's community
                    my_comm = results.community_membership[ref_idx]

                    # Count inter-community edges (edges to different communities)
                    inter_comm_edges = 0
                    if graph.nodes and ref_idx < graph.num_nodes:
                        node = &graph.nodes[ref_idx]
                        for ni in range(node.degree):
                            neighbor_id = node.neighbors[ni]
                            if neighbor_id < graph.num_nodes:
                                neighbor_comm = results.community_membership[neighbor_id]
                                if neighbor_comm != my_comm:
                                    inter_comm_edges += 1

                    # Bridge score = combination of:
                    # - Community diversity (num_neighbor_communities)
                    # - Inter-community connectivity (inter_comm_edges)
                    # Normalize both to [0, 1] and combine
                    comm_diversity = <float>results.num_neighbor_communities[ref_idx] / <float>max_communities
                    edge_score = <float>inter_comm_edges / <float>(node.degree + 1)  # +1 to avoid div by 0

                    # Weighted combination: 60% community diversity, 40% edge connectivity
                    bridge_score = 0.6 * comm_diversity + 0.4 * edge_score
                    results.betweenness_centrality[ref_idx] = bridge_score

            if verbose:
                bf_nogil_logf_notime(b"COMMUNITY", "Bridge scores assigned (hybrid: community diversity + inter-comm edges)\n")

            # Optional: Calculate exact betweenness if very few candidates (< 100)
            if bridge_candidate_count < 100:
                if verbose:
                    bf_nogil_logf_notime(b"COMMUNITY", "Few candidates (%u) - calculating exact betweenness...\n", bridge_candidate_count)

                # Compute normalization factor for exact betweenness
                if graph.num_nodes > 2:
                    max_betweenness = <float>((graph.num_nodes - 1) * (graph.num_nodes - 2)) / 2.0
                    if max_betweenness > 0:
                        normalization_factor = 1.0 / max_betweenness
                    else:
                        normalization_factor = 1.0
                else:
                    normalization_factor = 1.0

                # Calculate exact betweenness for refinement
                if bridge_candidate_count >= 10:
                    if verbose:
                        bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: mode=batch_optimized candidates=%u\n",
                                           bridge_candidate_count)

                    # Single batched call for all bridge candidates
                    ret = igraph_vector_init(&betweenness_vec_bridges, 0)
                    if ret == IGRAPH_SUCCESS:
                        ret = igraph_vs_vector(&vs_bridges, &bridge_candidates)
                        if ret == IGRAPH_SUCCESS:
                            # Single call to igraph_betweenness for all candidates (most efficient)
                            ret = igraph_betweenness(ig_graph, &betweenness_vec_bridges, vs_bridges, 0, ig_weights)
                            igraph_vs_destroy(&vs_bridges)

                            if ret == IGRAPH_SUCCESS:
                                # Store results
                                for i in range(bridge_candidate_count):
                                    ref_idx = <uint32_t>get_vector_int_element(&bridge_candidates, i)
                                    if ref_idx < graph.num_nodes:
                                        results.betweenness_centrality[ref_idx] = <float>VECTOR(betweenness_vec_bridges)[i] * normalization_factor

                                if verbose:
                                    bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: batch_optimized completed\n")
                            else:
                                if verbose:
                                    bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: WARNING calculation_failed\n")

                            igraph_vector_destroy(&betweenness_vec_bridges)
                        else:
                            if verbose:
                                bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: WARNING vertex_selector_failed\n")
                    else:
                        if verbose:
                            bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: WARNING vector_init_failed\n")

                else:
                    # Sequential calculation for small bridge sets or single thread
                    if verbose:
                        bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: mode=sequential candidates=%u\n",
                                           bridge_candidate_count)

                    ret = igraph_vector_init(&betweenness_vec_bridges, 0)
                    if ret == IGRAPH_SUCCESS:
                        ret = igraph_vs_vector(&vs_bridges, &bridge_candidates)
                        if ret == IGRAPH_SUCCESS:
                            ret = igraph_betweenness(ig_graph, &betweenness_vec_bridges, vs_bridges, 0, ig_weights)
                            igraph_vs_destroy(&vs_bridges)

                            if ret == IGRAPH_SUCCESS:
                                for i in range(bridge_candidate_count):
                                    ref_idx = <uint32_t>get_vector_int_element(&bridge_candidates, i)
                                    if ref_idx < graph.num_nodes:
                                        results.betweenness_centrality[ref_idx] = <float>VECTOR(betweenness_vec_bridges)[i] * normalization_factor

                                if verbose:
                                    bf_nogil_logf_notime(b"COMMUNITY", "  OK Sequential betweenness calculation completed\n")
                            else:
                                if verbose:
                                    bf_nogil_logf_notime(b"COMMUNITY", "[WARNING] Betweenness calculation failed\n")

                            igraph_vector_destroy(&betweenness_vec_bridges)
                        else:
                            if verbose:
                                bf_nogil_logf_notime(b"COMMUNITY", "[WARNING] Failed to create vertex selector\n")
                    else:
                        if verbose:
                            bf_nogil_logf_notime(b"COMMUNITY", "[WARNING] Failed to initialize betweenness vector\n")
            else:
                if verbose:
                    bf_nogil_logf_notime(b"COMMUNITY", "  WARNING Warning: Too many bridge candidates (%u > %u limit)\n",
                                       bridge_candidate_count, max_bridge_nodes)
                    bf_nogil_logf_notime(b"COMMUNITY", "  -> Skipping betweenness calculation (would be too slow)\n")
                    bf_nogil_logf_notime(b"COMMUNITY", "  -> Bridge detection will use num_neighbor_communities only\n")
        else:
            if verbose:
                if bridge_candidate_count == 0:
                    bf_nogil_logf_notime(b"COMMUNITY", "  OK No bridge candidates found\n")
                    bf_nogil_logf_notime(b"COMMUNITY", "  -> All nodes in homogeneous communities (good!)\n")
                else:
                    bf_nogil_logf_notime(b"COMMUNITY", "  WARNING All nodes are bridge candidates\n")
                    bf_nogil_logf_notime(b"COMMUNITY", "  -> Skipping betweenness (pathological case)\n")

        igraph_vector_int_destroy(&bridge_candidates)
    else:
        if verbose:
            bf_nogil_logf_notime(b"COMMUNITY", "betweenness_centrality: ERROR failed_to_init_candidate_list\n")

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_step_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "Bridge detection and betweenness complete (%.2fs)\n", elapsed_sec)

    # NOTE: We don't destroy the cached igraph here - it belongs to WeightedGraph
    # and will be freed when destroy_weighted_graph() is called

    t_step_end = clock()
    elapsed_sec = (t_step_end - t_phase_start) / <double>CLOCKS_PER_SEC
    bf_nogil_logf_notime(b"COMMUNITY", "TOTAL community_clustering time: %.2fs\n", elapsed_sec)

    return results
