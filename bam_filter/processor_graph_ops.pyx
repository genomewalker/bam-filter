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
Graph operations module for co-mapping analysis.

This module provides graph data structures and operations separated from
the Leiden clustering algorithm to maintain clean separation of concerns.

RESPONSIBILITIES:
- Graph data structures (WeightedGraph, GraphNode)
- Graph creation and destruction
- Graph building from alignment data
- Edge weight pruning using tail-focused broken-stick model
- Graph statistics calculation

EDGE WEIGHT FILTERING:
- Uses tail-focused broken-stick model to detect threshold
- Focuses on upper tail of edge weight distribution (typically 75th-90th percentile)
- Detects where high-weight edges exceed harmonic expectations
- INTERPRETATION: Edge weight >= threshold → KEEP, < threshold → REMOVE
- Consistent with clustering coefficient filtering approach

PERFORMANCE OPTIMIZATIONS:
- Binary search for edge lookup: O(log d) instead of O(d)
- Sorted neighbor lists maintained during insertion
- Direct read iteration for O(R) graph building
- Memory-efficient incremental edge insertion
- Exact edge counting before allocation (no reallocation overhead)

COMPLEXITY:
- Graph building: O(R × k²) where R = unique reads, k = avg refs/read
- Edge pruning: O(E) where E = number of edges
- Statistics calculation: O(N × d_avg) where N = nodes, d_avg = average degree
"""

from libc.stdlib cimport malloc, calloc, free, realloc, qsort
from libc.string cimport memset
from libc.string cimport memcpy
from libc.stdint cimport uint32_t, uint64_t, int32_t, uint8_t

# Define UINT32_MAX constant
cdef uint32_t UINT32_MAX = 0xFFFFFFFF

from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReadIndex, ReferencePattern
from bam_filter.processor_graph_ops cimport GraphNode, WeightedGraph

# igraph C API (used for direct igraph creation)
from bam_filter.processor_igraph cimport *
from libc.math cimport ceil

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef inline igraph_integer_t get_vector_int_element_local(igraph_vector_int_t *v, igraph_integer_t i) nogil:
    """Fast local accessor for an igraph integer vector element.

    This function accesses the internal storage of an ``igraph_vector_int_t``
    directly and returns the element at index ``i``. It is intended for use
    in nogil code paths where calling the public igraph API would be too
    expensive. The caller must ensure the index is in bounds.

    Parameters
    ----------
    v : igraph_vector_int_t*
        Pointer to the igraph integer vector.
    i : igraph_integer_t
        Index to read.

    Returns
    -------
    igraph_integer_t
        The integer value at index ``i``.
    """
    return (<igraph_integer_t**>v)[0][i]


cdef void* safe_realloc(void* ptr, size_t old_bytes, size_t new_bytes) noexcept nogil:
    """Reallocate memory using malloc+memcpy+free as a safer alternative to realloc.

    This helper allocates a new buffer of ``new_bytes``, copies the lesser of
    ``old_bytes`` and ``new_bytes`` from the original pointer, frees the old
    buffer and returns the new pointer. It is used as a safer replacement for
    calling the platform ``realloc`` directly in tight nogil code paths.

    Parameters
    ----------
    ptr : void*
        Pointer to the previously-allocated buffer (may be ``NULL``).
    old_bytes : size_t
        Size in bytes of the original buffer.
    new_bytes : size_t
        Desired size in bytes for the new buffer.

    Returns
    -------
    void*
        Pointer to the newly-allocated buffer on success, or ``NULL`` if
        allocation fails. If ``ptr`` is ``NULL`` this behaves like ``malloc``.

    Notes
    -----
    - This function does not attempt to detect or protect against callers
      passing already-freed pointers. It exists primarily to avoid certain
      platform-specific realloc failure modes and to centralize allocation
      logging if needed. The function is ``nogil`` and uses only C library
      functions.
    """
    cdef void* newp
    if not ptr:
        newp = malloc(new_bytes)
        return newp

    newp = malloc(new_bytes)
    if not newp:
        return NULL
    # copy the lesser of old_bytes/new_bytes
    if old_bytes > new_bytes:
        memcpy(newp, ptr, new_bytes)
    else:
        memcpy(newp, ptr, old_bytes)
    free(ptr)
    return newp

# ==============================================================================
# WEIGHTED GRAPH MANAGEMENT
# ==============================================================================

cdef WeightedGraph* create_weighted_graph(uint32_t num_nodes) except NULL nogil:
    """Create an empty weighted graph with specified number of nodes.

    Parameters
    ----------
    num_nodes : uint32_t
        Number of vertices in the graph (typically the number of references).

    Returns
    -------
    WeightedGraph*
        Pointer to newly allocated WeightedGraph, or ``NULL`` on failure.
    """
    cdef WeightedGraph* graph = <WeightedGraph*>malloc(sizeof(WeightedGraph))
    cdef uint32_t i

    if not graph:
        return NULL

    graph.num_nodes = num_nodes
    graph.total_weight = 0
    graph.num_edges = 0

    graph.nodes = <GraphNode*>calloc(num_nodes, sizeof(GraphNode))
    if not graph.nodes:
        free(graph)
        return NULL

    # Initialize each node
    for i in range(num_nodes):
        graph.nodes[i].neighbors = NULL
        graph.nodes[i].weights = NULL
        graph.nodes[i].degree = 0
        graph.nodes[i].capacity = 0
        graph.nodes[i].node_weight = 0

    # Initialize TSV data fields to NULL
    graph.tsv_total_reads = NULL
    graph.tsv_multimap_reads = NULL
    graph.tsv_alignments_per_ref = NULL
    graph.tsv_exact_connection_counts = NULL
    graph.tsv_co_mapping_averages = NULL
    graph.tsv_max_co_mappings = NULL
    graph.tsv_co_mapping_counts = NULL
    graph.tsv_neighbor_multimap_avg = NULL
    graph.tsv_neighbor_connections_avg = NULL
    graph.tsv_neighbor_counts = NULL
    graph.tsv_dataset_median_connections = 0.0
    graph.tsv_array_size = 0
    graph.tsv_min_read_count = 0

    return graph


# ---------------------------------------------------------------------------
# Build igraph directly from ReadIndex (reuse WeightedGraph builder)
# ---------------------------------------------------------------------------
cdef int build_igraph_from_read_index(
    void* ig_graph_ptr,
    void* ig_weights_ptr,
    MemoryPool* pool,
    ReadIndex* read_index,
    ReferenceStats* ref_stats,
    uint32_t num_refs,
    uint32_t min_read_count,
    uint32_t min_edge_weight,
    int num_threads,
    int verbose
) except -1 nogil:
    """Build an igraph representation from a ReadIndex.

    This routine builds a full in-memory weighted graph from the provided
    ``ReadIndex`` (via the internal WeightedGraph builder), optionally prunes
    edges below ``min_edge_weight``, and converts the resulting structure into
    an ``igraph_t`` and an associated edge weight vector. The function
    initializes the opaque igraph pointers provided by the caller.

    Parameters
    ----------
    ig_graph_ptr : void*
        Pointer to an uninitialized ``igraph_t`` storage (caller-allocated).
    ig_weights_ptr : void*
        Pointer to an uninitialized ``igraph_vector_t`` storage for edge weights.
    pool : MemoryPool*
        Memory pool containing alignments and read indices.
    read_index : ReadIndex*
        Read-index mapping references -> read lists used to compute co-mappings.
    ref_stats : ReferenceStats*
        Optional per-reference stats used to filter low-coverage references.
    num_refs : uint32_t
        Number of references in the compacted namespace.
    min_read_count : uint32_t
        Minimum reads required for a reference to be considered.
    min_edge_weight : uint32_t
        Minimum shared-read weight to retain an edge (pruning threshold).
    num_threads : int
        Number of threads hint for internal builders (may be ignored).
    verbose : int
        Verbosity flag (non-zero prints progress messages).

    Returns
    -------
    int
        0 on success, -1 on failure. On success the caller-owned igraph
        containers pointed to by ``ig_graph_ptr`` and ``ig_weights_ptr`` are
        initialized and must be destroyed by the caller when no longer needed.

    Notes
    -----
    - The function allocates and frees an intermediate WeightedGraph; the
      created ``igraph_t`` is a separate object that the caller controls.
    - This function is ``nogil`` and must not perform Python API calls.
    """
    cdef WeightedGraph* graph = NULL
    cdef igraph_t* ig_graph = <igraph_t*>ig_graph_ptr
    cdef igraph_vector_t* ig_weights = <igraph_vector_t*>ig_weights_ptr
    cdef igraph_vector_int_t edges
    cdef int ret
    cdef uint64_t edge_i = 0
    cdef uint64_t out_idx = 0
    cdef uint32_t node_i, neigh

    # Build exact weighted graph from alignments (uses ReadIndex internally)
    graph = build_weighted_graph_from_alignments(pool, ref_stats, read_index, num_refs, <int32_t>min_read_count)
    if not graph:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: build_weighted_graph_from_alignments failed\n")
        return -1

    # Optionally prune low-weight edges
    if min_edge_weight > 1:
        prune_low_weight_edges(graph, min_edge_weight)

    # Number of undirected edges (graph.num_edges) should be exact
    if verbose:
        bf_nogil_logf_notime(b"IGRAPH OPS", "Converting WeightedGraph (%llu edges) to igraph...\n", <unsigned long long>graph.num_edges)

    # Initialize igraph vectors
    ret = igraph_vector_int_init(&edges, <uint64_t>graph.num_edges * 2)
    if ret != 0:
        destroy_weighted_graph(graph)
        return -1

    ret = igraph_vector_init(ig_weights, <uint64_t>graph.num_edges)
    if ret != 0:
        igraph_vector_int_destroy(&edges)
        destroy_weighted_graph(graph)
        return -1

    # Fill edges and weights, emitting each undirected edge once (neighbor > node)
    out_idx = 0
    for node_i in range(graph.num_nodes):
        for edge_i in range(graph.nodes[node_i].degree):
            neigh = graph.nodes[node_i].neighbors[edge_i]
            if neigh > node_i:
                igraph_vector_int_set(&edges, out_idx * 2, <igraph_integer_t>node_i)
                igraph_vector_int_set(&edges, out_idx * 2 + 1, <igraph_integer_t>neigh)
                igraph_vector_set(ig_weights, out_idx, <igraph_real_t>graph.nodes[node_i].weights[edge_i])
                out_idx += 1

    # Sanity: out_idx should equal graph.num_edges
    if out_idx != graph.num_edges:
        if verbose:
            bf_nogil_logf_notime(
                b"IGRAPH OPS",
                "WARNING: expected %lu edges but wrote %lu\n",
                <unsigned long>graph.num_edges,
                <unsigned long>out_idx,
            )

    # Create igraph
    ret = igraph_create(ig_graph, &edges, <igraph_integer_t>num_refs, 0)
    igraph_vector_int_destroy(&edges)

    if ret != 0:
        igraph_vector_destroy(ig_weights)
        destroy_weighted_graph(graph)
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: igraph_create failed\n")
        return -1

    if verbose:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "igraph created: %ld vertices, %ld edges\n",
            <long>igraph_vcount(ig_graph),
            <long>igraph_ecount(ig_graph),
        )

    # Cleanup weighted graph
    destroy_weighted_graph(graph)

    return 0


# ---------------------------------------------------------------------------
# Direct-to-igraph builder (avoid full WeightedGraph allocation)
# ---------------------------------------------------------------------------

# Add this counting function near the top of your file, after the SeenTracker functions:

cdef uint64_t count_edges_from_read_index(
    ReadIndex* read_index,
    ReferenceStats* ref_stats,
    uint32_t num_refs,
    uint32_t min_read_count,
    uint32_t min_edge_weight,
    MemoryPool* pool,
    int verbose
) nogil:
    """Count the exact number of undirected edges meeting a minimum weight.

    This function scans per-reference read lists and computes the number of
    unique undirected reference pairs whose shared read count meets or
    exceeds ``min_edge_weight``. It is intended to be used to pre-size
    igraph edge arrays so the direct builder can allocate exact capacity.

    Parameters
    ----------
    read_index : ReadIndex*
        ReadIndex mapping references -> read lists.
    ref_stats : ReferenceStats*
        Optional per-reference stats used to skip low-coverage references.
    num_refs : uint32_t
        Number of references in the compacted namespace.
    min_read_count : uint32_t
        Minimum reads required for a reference to participate in counting.
    min_edge_weight : uint32_t
        Minimum shared reads for an edge to be counted.
    pool : MemoryPool*
        Memory pool with alignments used to expand reads -> references per read.
    verbose : int
        Verbosity flag (non-zero prints progress messages).

    Returns
    -------
    uint64_t
        Exact number of undirected edges that meet the specified threshold.
        Returns 0 if allocation fails or if no edges meet the criteria.
    """
    cdef uint32_t ref_idx, other_ref_idx
    cdef uint64_t read_idx64
    cdef uint32_t read_i, read_count
    cdef uint64_t aln_start_i, aln_end_i, aln_idx64
    cdef uint64_t edge_count = 0
    cdef uint32_t* edge_weights = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    cdef uint32_t* neighbors = NULL
    cdef uint32_t neighbor_count = 0
    cdef uint32_t neighbor_capacity = 0
    cdef uint32_t kk
    cdef uint32_t* new_neighbors
    cdef uint64_t new_cap
    cdef uint64_t total_neighbors = 0  # DIAGNOSTIC
    cdef uint32_t max_neighbors = 0    # DIAGNOSTIC
    cdef uint8_t* read_seen_refs = NULL
    
    if not edge_weights:
        return 0
    
    if verbose:
        bf_nogil_logf_notime(b"IGRAPH OPS", "Counting edges from ReadIndex...\n")
    
    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue
        
        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0 or not read_index.ref_to_reads[ref_idx]:
            continue
        
        neighbor_count = 0

        # Accumulate weights AND track which neighbors exist
        # Use a temporary marker array to track which references this read has seen
        read_seen_refs = <uint8_t*>calloc(num_refs, sizeof(uint8_t))
        if not read_seen_refs:
            if neighbors:
                free(neighbors)
            free(edge_weights)
            return 0

        for read_i in range(read_count):
            read_idx64 = read_index.ref_to_reads[ref_idx][read_i]
            if read_idx64 >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx64]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx64]

            if aln_end_i > pool.alignment_count:
                continue

            # First pass: mark which references this read maps to (unique per read)
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                # Mark that this read has an alignment to other_ref_idx
                read_seen_refs[other_ref_idx] = 1

            # Second pass: increment edge weight once per unique reference for this read
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                # Only increment once per read (check and clear the marker)
                if read_seen_refs[other_ref_idx] == 1:
                    # First time seeing this neighbor reference across all reads?
                    if edge_weights[other_ref_idx] == 0:
                        if neighbor_count >= neighbor_capacity:
                            new_cap = neighbor_capacity * 2 if neighbor_capacity > 0 else 256
                            new_neighbors = <uint32_t*>safe_realloc(
                                neighbors,
                                neighbor_capacity * sizeof(uint32_t),
                                new_cap * sizeof(uint32_t)
                            )
                            if not new_neighbors:
                                if neighbors:
                                    free(neighbors)
                                free(edge_weights)
                                free(read_seen_refs)
                                return 0
                            neighbors = new_neighbors
                            neighbor_capacity = new_cap
                        neighbors[neighbor_count] = other_ref_idx
                        neighbor_count += 1

                    # Increment edge weight once per unique shared read
                    edge_weights[other_ref_idx] += 1
                    # Clear marker so we don't count this reference again for this read
                    read_seen_refs[other_ref_idx] = 0
        
        # DIAGNOSTIC: Track neighbor counts
        total_neighbors += neighbor_count
        if neighbor_count > max_neighbors:
            max_neighbors = neighbor_count

        # Count edges meeting threshold (ONLY from neighbor list)
        for kk in range(neighbor_count):
            other_ref_idx = neighbors[kk]
            if edge_weights[other_ref_idx] >= min_edge_weight and other_ref_idx > ref_idx:
                edge_count += 1

        # Reset weights AFTER counting (reset ALL neighbors)
        for kk in range(neighbor_count):
            edge_weights[neighbors[kk]] = 0

        # Free per-reference temporary array
        free(read_seen_refs)

    if neighbors:
        free(neighbors)
    free(edge_weights)
    
    if verbose:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "Exact edge count: %llu (%.1f MB for arrays)\n",
            edge_count,
            (edge_count * 12.0) / (1024 * 1024),
        )
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "Total neighbors seen: %llu, Max neighbors for single ref: %u\n",
            total_neighbors,
            max_neighbors,
        )
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "Avg neighbors per ref: %.1f\n",
            <double>total_neighbors / <double>num_refs,
        )
    
    return edge_count


# ---------------------------------------------------------------------------
# Broken-stick based min-edge-weight picker
# ---------------------------------------------------------------------------


cdef int _uint32_compare(const void* a, const void* b) noexcept nogil:
    """Compare two uint32_t values for qsort-style sorting.

    Parameters
    ----------
    a, b : const void*
        Pointers to 32-bit unsigned integer elements being compared.

    Returns
    -------
    int
        -1 if *a < *b, 1 if *a > *b, 0 if equal.
    """
    cdef uint32_t va = (<uint32_t*>a)[0]
    cdef uint32_t vb = (<uint32_t*>b)[0]
    if va < vb:
        return -1
    elif va > vb:
        return 1
    else:
        return 0


# Local float comparator for qsort (ascending) — defined here so this module
# can call qsort on float arrays without depending on other compilation units.
cdef int _float_compare_ascending(const void* a, const void* b) noexcept nogil:
    """Compare two floats for ascending qsort ordering.

    Parameters
    ----------
    a, b : const void*
        Pointers to float elements.

    Returns
    -------
    int
        -1 if *a < *b, 1 if *a > *b, 0 if equal.
    """
    cdef float fa = (<float*>a)[0]
    cdef float fb = (<float*>b)[0]
    if fa < fb:
        return -1
    elif fa > fb:
        return 1
    else:
        return 0


cdef uint32_t pick_min_edge_weight_broken_stick(
    MemoryPool* pool,
    ReadIndex* read_index,
    ReferenceStats* ref_stats,
    uint32_t num_refs,
    uint32_t min_read_count,
    double tol,
    int verbose,
    double tail_percentile,
    uint32_t min_tail_size
) noexcept nogil:
    """Suggest an integer min_edge_weight using a tail-focused broken-stick model.

    The routine analyzes the distribution of undirected edge weights (shared
    read counts) and returns a conservative integer threshold. It focuses on
    the upper tail of the distribution to avoid being dominated by the large
    number of low-weight noise edges.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignment data.
    read_index : ReadIndex*
        ReadIndex used to enumerate edges.
    ref_stats : ReferenceStats*
        Optional per-reference stats used to skip low-coverage refs.
    num_refs : uint32_t
        Number of references.
    min_read_count : uint32_t
        Minimum reads for a reference to be included.
    tol : double
        Tolerance factor used to detect deviation from harmonic expectation.
    verbose : int
        Verbosity flag.
    tail_percentile : double
        Fraction (0-1) indicating tail start (e.g., 0.75 for 75th percentile).
    min_tail_size : uint32_t
        Minimum number of tail samples required to run tail-focused detection.

    Returns
    -------
    uint32_t
        Suggested integer threshold; returns 0 when insufficient data is
        available (no filtering) or when allocations fail.
    """
    cdef uint64_t edge_count = 0
    cdef uint32_t* edge_weights = NULL
    cdef uint32_t* neighbors = NULL
    cdef uint32_t neighbor_count = 0
    cdef uint32_t neighbor_capacity = 0
    cdef uint64_t i, kk
    cdef uint64_t read_idx64
    cdef uint32_t read_i, read_count
    cdef uint64_t aln_start_i, aln_end_i, aln_idx64
    cdef uint32_t other_ref_idx
    cdef uint32_t ref_idx
    cdef uint32_t new_cap_i
    cdef uint32_t* new_neighbors
    cdef uint64_t written = 0
    cdef uint32_t* weights = NULL
    cdef float* sorted_vals = NULL
    cdef double sumv = 0.0
    cdef float* harmonic_suffix = NULL
    cdef uint32_t cutoff_rank = 0
    cdef uint32_t threshold_u = 0
    cdef float expected
    cdef uint64_t count_w1 = 0, count_w2 = 0, count_w3 = 0, count_w5 = 0, count_w10 = 0
    # Tail-focused analysis temporaries (declare at function scope)
    cdef uint32_t tail_start = 0
    cdef uint32_t tail_len = 0
    cdef double tail_sum = 0.0
    cdef uint32_t jj = 0
    cdef float* tail_harmonic = NULL
    cdef int tail_cutoff = -1
    cdef uint32_t threshold_tail = 0
    cdef uint32_t threshold_pct = 0
    cdef uint32_t chosen_tail_thr = 0
    cdef uint8_t* read_seen_refs = NULL

    # First, count edges to know exact allocation size
    edge_count = count_edges_from_read_index(read_index, ref_stats, num_refs, min_read_count, 1, pool, verbose)
    if edge_count == 0:
        return 0

    # Allocate arrays
    weights = <uint32_t*>malloc(edge_count * sizeof(uint32_t))
    if not weights:
        return 0

    # Temporary accumulator for neighbor weights
    edge_weights = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    if not edge_weights:
        free(weights)
        return 0

    # Iterate reads and collect undirected edge weights (only other_ref_idx > ref_idx)
    written = 0
    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue

        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0 or not read_index.ref_to_reads[ref_idx]:
            continue

        neighbor_count = 0
        neighbor_capacity = 0
        neighbors = NULL

        # Temporary marker array to track which references each read has seen
        read_seen_refs = <uint8_t*>calloc(num_refs, sizeof(uint8_t))
        if not read_seen_refs:
            free(edge_weights)
            free(weights)
            return 0

        # accumulate counts for neighbors of this ref (count each read once per reference)
        for read_i in range(read_count):
            read_idx64 = read_index.ref_to_reads[ref_idx][read_i]
            if read_idx64 >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx64]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx64]
            if aln_end_i > pool.alignment_count:
                continue

            # First pass: mark which references this read maps to
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1
                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue
                read_seen_refs[other_ref_idx] = 1

            # Second pass: increment edge weight once per unique reference for this read
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1
                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                if read_seen_refs[other_ref_idx] == 1:
                    if edge_weights[other_ref_idx] == 0:
                        # add to neighbor list
                        if neighbor_count >= neighbor_capacity:
                            new_cap_i = neighbor_capacity * 2 if neighbor_capacity > 0 else 256
                            new_neighbors = <uint32_t*>safe_realloc(neighbors, neighbor_capacity * sizeof(uint32_t), new_cap_i * sizeof(uint32_t))
                            if not new_neighbors:
                                if neighbors:
                                    free(neighbors)
                                free(edge_weights)
                                free(weights)
                                free(read_seen_refs)
                                return 0
                            neighbors = new_neighbors
                            neighbor_capacity = new_cap_i
                        neighbors[neighbor_count] = other_ref_idx
                        neighbor_count += 1

                    edge_weights[other_ref_idx] += 1
                    read_seen_refs[other_ref_idx] = 0  # Clear marker

        # For each neighbor > ref_idx add undirected edge weight
        for kk in range(neighbor_count):
            other_ref_idx = neighbors[kk]
            if other_ref_idx > ref_idx:
                weights[written] = edge_weights[other_ref_idx]
                written += 1

        # reset accumulator for seen neighbors
        for kk in range(neighbor_count):
            edge_weights[neighbors[kk]] = 0

        # Free per-reference temporary array
        free(read_seen_refs)

        if neighbors:
            free(neighbors)

    # written should equal edge_count (or less if min_edge_weight >1 in first pass)
    if written == 0:
        free(edge_weights)
        free(weights)
        return 0

    # convert to float and sort ascending
    sorted_vals = <float*>malloc(written * sizeof(float))
    if not sorted_vals:
        free(edge_weights)
        free(weights)
        return 0
    sumv = 0.0
    for i in range(written):
        sorted_vals[i] = <float>weights[i]
        sumv += sorted_vals[i]

    qsort(sorted_vals, written, sizeof(float), _float_compare_ascending)

    if verbose:
        bf_nogil_logf_notime(b"BROKEN-STICK", "Collected %llu edge weights for analysis\n", <unsigned long long>written)
        bf_nogil_logf_notime(b"BROKEN-STICK", "Edge weight statistics:\n")
        bf_nogil_logf_notime(b"BROKEN-STICK", "  Min: %.0f, Max: %.0f, Mean: %.2f\n", sorted_vals[0], sorted_vals[written-1], sumv / <double>written)
        bf_nogil_logf_notime(
            b"BROKEN-STICK",
            "  Median: %.0f, 75th percentile: %.0f, 90th percentile: %.0f\n",
            sorted_vals[written // 2],
            sorted_vals[<uint64_t>(written * 0.75)],
            sorted_vals[<uint64_t>(written * 0.9)],
        )
        bf_nogil_logf_notime(
            b"BROKEN-STICK",
            "  95th percentile: %.0f, 99th percentile: %.0f, 99.9th percentile: %.0f\n",
            sorted_vals[<uint64_t>(written * 0.95)],
            sorted_vals[<uint64_t>(written * 0.99)],
            sorted_vals[<uint64_t>(written * 0.999)],
        )
        bf_nogil_logf_notime(b"BROKEN-STICK", "  Total weight sum: %.0f\n", sumv)
        
        # Count edges at low weights for distribution insight
        for i in range(written):
            if sorted_vals[i] == 1:
                count_w1 += 1
            elif sorted_vals[i] == 2:
                count_w2 += 1
            elif sorted_vals[i] == 3:
                count_w3 += 1
            elif sorted_vals[i] <= 5:
                count_w5 += 1
            elif sorted_vals[i] <= 10:
                count_w10 += 1
        bf_nogil_logf_notime(b"BROKEN-STICK", "  Distribution (cumulative counts):\n")
        bf_nogil_logf_notime(
            b"BROKEN-STICK",
            "    Weight=1: %llu (%.1f%%), Weight<=2: %llu (%.1f%%), Weight<=3: %llu (%.1f%%)\n",
            count_w1,
            100.0 * count_w1 / written,
            count_w1 + count_w2,
            100.0 * (count_w1 + count_w2) / written,
            count_w1 + count_w2 + count_w3,
            100.0 * (count_w1 + count_w2 + count_w3) / written,
        )
        bf_nogil_logf_notime(
            b"BROKEN-STICK",
            "    Weight<=5: %llu (%.1f%%), Weight<=10: %llu (%.1f%%)\n",
            count_w5,
            100.0 * count_w5 / written,
            count_w10,
            100.0 * count_w10 / written,
        )

    # TAIL-FOCUSED BROKEN-STICK (consistent with CC filtering approach)
    # Run broken-stick detection only on the upper tail of the distribution
    # (above tail_percentile, typically 75th-90th percentile for edge weights).
    #
    # RATIONALE:
    # - Edge weight distributions are heavily skewed: massive number of low-weight edges (1-2 reads)
    # - These low-weight edges dominate global broken-stick calculations
    # - Tail-focused approach detects meaningful breaks in the high-weight signal region
    # - Ignores the noise-dominated lower distribution
    #
    # INTERPRETATION:
    # - Edges with weight >= threshold are true co-mapping signal → KEEP
    # - Edges with weight < threshold are spurious co-mapping → REMOVE
    if tail_percentile > 0.0 and tail_percentile < 1.0 and written >= min_tail_size:
        tail_start = <uint32_t>(written * tail_percentile)
        if tail_start >= written:
            tail_start = written - 1
        tail_len = <uint32_t>(written - tail_start)
        if tail_len >= min_tail_size:
            # compute tail sum
            tail_sum = 0.0
            for jj in range(tail_start, written):
                tail_sum += sorted_vals[jj]

            # allocate harmonic suffix for tail
            tail_harmonic = <float*>malloc(tail_len * sizeof(float))
            if tail_harmonic:
                for jj in range(tail_len - 1, -1, -1):
                    if jj == tail_len - 1:
                        tail_harmonic[jj] = 1.0 / (<float>(jj + 1))
                    else:
                        tail_harmonic[jj] = tail_harmonic[jj + 1] + 1.0 / (<float>(jj + 1))

                # Detect cutoff within tail: find where observed > expected * (1 + tolerance)
                # This identifies where edge weights exceed uniform harmonic expectations
                tail_cutoff = -1
                for jj in range(tail_len):
                    observed = sorted_vals[tail_start + jj]
                    expected_tail = tail_harmonic[jj] * (<float>(tail_sum) / <float>tail_len)
                    if observed > expected_tail * (1.0 + <float>tol):
                        tail_cutoff = jj
                        break

                threshold_tail = 0
                if tail_cutoff >= 0:
                    threshold_tail = <uint32_t>ceil(sorted_vals[tail_start + tail_cutoff])

                # Percentile candidate at tail start (conservative fallback)
                # If no break detected in tail, use the tail start percentile
                threshold_pct = <uint32_t>ceil(sorted_vals[tail_start])

                # Choose the MORE CONSERVATIVE (higher) threshold
                # Higher threshold = keep more edges (less aggressive filtering)
                chosen_tail_thr = threshold_pct
                if threshold_tail > chosen_tail_thr:
                    chosen_tail_thr = threshold_tail

                if verbose:
                    bf_nogil_logf_notime(
                        b"BROKEN-STICK",
                        "Tail-focused analysis: tail_start=%u tail_len=%u min_tail_size=%u\n",
                        tail_start,
                        tail_len,
                        min_tail_size,
                    )
                    if threshold_tail > 0:
                        bf_nogil_logf_notime(
                            b"BROKEN-STICK",
                            "Tail broken-stick cutoff observed: %u (tail_cutoff=%d)\n",
                            threshold_tail,
                            tail_cutoff,
                        )
                    else:
                        bf_nogil_logf_notime(
                            b"BROKEN-STICK",
                            "No tail broken-stick cutoff detected; using tail percentile candidate: %u\n",
                            threshold_pct,
                        )
                    bf_nogil_logf_notime(b"BROKEN-STICK", "Chosen tail-focused threshold: %u\n", chosen_tail_thr)

                free(tail_harmonic)

                # If tail produced a non-zero suggestion, use it immediately.
                if chosen_tail_thr > 0:
                    free(edge_weights)
                    free(weights)
                    free(sorted_vals)
                    return chosen_tail_thr

    # FALLBACK: Global broken-stick (if tail-focused didn't apply or succeed)
    # Only used when insufficient edges for tail analysis
    if written < 20:
        if verbose:
            bf_nogil_logf_notime(
                b"BROKEN-STICK",
                "Too few edges (%llu < 20), using 75th percentile fallback\n",
                <unsigned long long>written,
            )
        # Use 75th percentile as conservative fallback
        # This removes the bottom 75% of edges (keeps top 25% high-weight edges)
        cutoff_rank = <uint32_t>(written * 0.75)
        if cutoff_rank >= written:
            cutoff_rank = written - 1
        threshold_u = <uint32_t>ceil(sorted_vals[cutoff_rank])
        free(edge_weights)
        free(weights)
        free(sorted_vals)
        return threshold_u if threshold_u > 0 else 1

    # GLOBAL BROKEN-STICK (secondary fallback if tail-focused not used)
    # Apply classical broken-stick model to full distribution
    harmonic_suffix = <float*>malloc(written * sizeof(float))
    if not harmonic_suffix:
        free(edge_weights)
        free(weights)
        free(sorted_vals)
        return 0

    # Build harmonic suffix array for broken-stick calculation
    for i in range(written - 1, -1, -1):
        if i == written - 1:
            harmonic_suffix[i] = 1.0 / (<float>(i + 1))
        else:
            harmonic_suffix[i] = harmonic_suffix[i + 1] + 1.0 / (<float>(i + 1))

    # Detect where observed exceeds expected by tolerance
    cutoff_rank = written  # default: no cutoff found
    for i in range(written):
        expected = harmonic_suffix[i] * (sumv / <float>written)
        if sorted_vals[i] > expected * (1.0 + <float>tol):
            cutoff_rank = i
            break

    if cutoff_rank < written:
        # Broken-stick detected a cutoff
        threshold_u = <uint32_t>ceil(sorted_vals[cutoff_rank])
        if verbose:
            bf_nogil_logf_notime(
                b"BROKEN-STICK",
                "Global broken-stick detected cutoff at rank %u/%llu\n",
                cutoff_rank,
                <unsigned long long>written,
            )
            bf_nogil_logf_notime(
                b"BROKEN-STICK",
                "Observed weight: %.0f, Expected: %.2f (tol=%.2f)\n",
                sorted_vals[cutoff_rank],
                harmonic_suffix[cutoff_rank] * (sumv / <float>written),
                <float>tol,
            )
    else:
        # No break detected - use 75th percentile as fallback
        if verbose:
            bf_nogil_logf_notime(b"BROKEN-STICK", "No global cutoff found, using 75th percentile fallback\n")
        cutoff_rank = <uint32_t>(written * 0.75)
        if cutoff_rank >= written:
            cutoff_rank = written - 1
        threshold_u = <uint32_t>ceil(sorted_vals[cutoff_rank])

    # Ensure threshold is at least 1 (never return 0 which would disable filtering)
    if threshold_u == 0:
        threshold_u = 1

    if verbose:
        bf_nogil_logf_notime(
            b"BROKEN-STICK",
            "Final edge weight threshold: %u (edges < %u will be REMOVED)\n",
            threshold_u,
            threshold_u,
        )

    free(harmonic_suffix)
    free(edge_weights)
    free(weights)
    free(sorted_vals)

    return threshold_u





cdef int build_igraph_direct_from_read_index(
    void* ig_graph_ptr,
    void* ig_weights_ptr,
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
) except -1 nogil:
    """Build an igraph directly from ReadIndex using exact pre-allocation.

    This builder first counts the exact number of edges that will be emitted
    and then allocates exact-sized arrays for edge "from", "to" and weights,
    avoiding intermediate growth and large reallocations. It then fills the
    edge arrays and creates an ``igraph_t`` and associated weight vector.

    Parameters
    ----------
    ig_graph_ptr : void*
        Pointer to caller-provided storage for an ``igraph_t``.
    ig_weights_ptr : void*
        Pointer to caller-provided storage for an ``igraph_vector_t`` of weights.
    pool : MemoryPool*
        Memory pool with alignment/read data.
    read_index : ReadIndex*
        ReadIndex mapping references to read lists.
    ref_stats : ReferenceStats*
        Optional per-reference stats (used to skip low-coverage refs).
    num_refs : uint32_t
        Number of references.
    min_read_count : uint32_t
        Minimum reads per reference to be considered.
    min_edge_weight : uint32_t
        Minimum shared-read weight to keep an edge.
    num_threads : int
        Number of threads hint for counting/build steps.
    verbose : int
        Verbosity flag; non-zero prints progress messages.
    out_n_components : igraph_integer_t*
        Optional output pointer to receive number of connected components.
    out_singletons : igraph_integer_t*
        Optional output pointer to receive number of singleton components.

    Returns
    -------
    int
        0 on success, -1 on failure. On success, the igraph and weight vector
        are initialized in the provided caller-owned storage.
    """
    cdef igraph_t* ig_graph = <igraph_t*>ig_graph_ptr
    cdef igraph_vector_t* ig_weights = <igraph_vector_t*>ig_weights_ptr
    cdef uint32_t ref_idx, other_ref_idx
    cdef uint64_t read_idx64
    cdef uint32_t read_i, read_count
    cdef uint64_t aln_start_i, aln_end_i, aln_idx64
    cdef int ret
    cdef uint32_t tmp_i, kk
    cdef bint success = False
    
    # Edge storage - EXACT SIZE
    cdef uint32_t* edge_from = NULL
    cdef uint32_t* edge_to = NULL
    cdef uint32_t* edge_weight = NULL
    cdef uint64_t edge_capacity = 0
    cdef uint64_t edge_count = 0
    
    # Neighbor tracking
    cdef uint32_t* neighbors = NULL
    cdef uint32_t neighbor_count = 0
    cdef uint32_t neighbor_capacity = 50000
    cdef uint32_t* new_neighbors = NULL
    cdef uint64_t new_cap
    
    # Weight accumulation
    cdef uint32_t* ref_weights = NULL

    # Read marker for tracking unique references per read
    cdef uint8_t* read_seen_marker = NULL

    # Component analysis
    cdef igraph_vector_int_t comp_membership
    cdef igraph_vector_int_t comp_sizes
    cdef igraph_integer_t n_components = 0
    cdef long i_comp, singletons = 0
    cdef bint comp_initialized = False
    cdef bint comp_sizes_initialized = False

    # igraph building
    cdef igraph_vector_int_t empty_edges
    cdef igraph_vector_int_t edges_vec

    if verbose:
        bf_nogil_logf_notime(b"IGRAPH OPS", "Direct igraph builder with exact pre-allocation\n")

    # Validate inputs
    if not read_index or not read_index.ref_read_counts or not read_index.ref_to_reads:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Invalid read_index\n")
        return -1
    if not pool.read_alignment_starts or not pool.read_alignment_counts or not pool.alignments:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Invalid pool arrays\n")
        return -1

    # STEP 1: Count edges exactly
    edge_capacity = count_edges_from_read_index(
        read_index, ref_stats, num_refs, min_read_count, min_edge_weight, pool, verbose
    )
    
    if edge_capacity == 0:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "No edges found, creating empty graph\n")
        ret = igraph_vector_int_init(&empty_edges, 0)
        ret = igraph_vector_init(ig_weights, 0)
        ret = igraph_create(ig_graph, &empty_edges, <igraph_integer_t>num_refs, 0)
        igraph_vector_int_destroy(&empty_edges)
        if out_n_components:
            out_n_components[0] = 0
        if out_singletons:
            out_singletons[0] = 0
        return 0 if ret == 0 else -1

    # STEP 2: Allocate EXACT size (no growth needed)
    if verbose:
        bf_nogil_logf_notime(b"IGRAPH OPS", "Allocating exact capacity: %llu edges\n", edge_capacity)
    
    edge_from = <uint32_t*>malloc(edge_capacity * sizeof(uint32_t))
    edge_to = <uint32_t*>malloc(edge_capacity * sizeof(uint32_t))
    edge_weight = <uint32_t*>malloc(edge_capacity * sizeof(uint32_t))
    
    if not edge_from or not edge_to or not edge_weight:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Failed to allocate edge arrays\n")
        if edge_from:
            free(edge_from)
        if edge_to:
            free(edge_to)
        if edge_weight:
            free(edge_weight)
        return -1

    ref_weights = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    if not ref_weights:
        free(edge_from)
        free(edge_to)
        free(edge_weight)
        return -1

    neighbors = <uint32_t*>malloc(neighbor_capacity * sizeof(uint32_t))
    if not neighbors:
        free(edge_from)
        free(edge_to)
        free(edge_weight)
        free(ref_weights)
        return -1

    # STEP 3: Build edges (no reallocation)
    if verbose:
        bf_nogil_logf_notime(b"IGRAPH OPS", "Building edges...\n")
    
    edge_count = 0
    
    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue

        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0 or not read_index.ref_to_reads[ref_idx]:
            continue

        if verbose and (ref_idx & 0x3FFF) == 0:
            bf_nogil_logf_notime(
                b"IGRAPH OPS",
                "Processing ref %u/%u, edges: %llu/%llu (%.1f%%)\n",
                ref_idx,
                num_refs,
                edge_count,
                edge_capacity,
                100.0 * edge_count / edge_capacity,
            )

        neighbor_count = 0

        # Allocate temporary marker array for tracking unique references per read
        read_seen_marker = <uint8_t*>calloc(num_refs, sizeof(uint8_t))
        if not read_seen_marker:
            if verbose:
                bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Failed to allocate read marker array\n")
            free(neighbors)
            free(edge_from)
            free(edge_to)
            free(edge_weight)
            free(ref_weights)
            return -1

        # Accumulate weights (count each read once per unique reference)
        for read_i in range(read_count):
            read_idx64 = read_index.ref_to_reads[ref_idx][read_i]
            if read_idx64 >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx64]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx64]

            if aln_end_i > pool.alignment_count:
                continue

            # First pass: mark which references this read maps to
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                read_seen_marker[other_ref_idx] = 1

            # Second pass: increment edge weight once per unique reference for this read
            aln_idx64 = aln_start_i
            while aln_idx64 < aln_end_i:
                other_ref_idx = pool.alignments[aln_idx64].reference_index
                aln_idx64 += 1

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                # Only increment once per read
                if read_seen_marker[other_ref_idx] == 1:
                    # First time seeing this neighbor reference across all reads?
                    if ref_weights[other_ref_idx] == 0:
                        if neighbor_count >= neighbor_capacity:
                            # Grow neighbors (only list, not edge arrays)
                            new_cap = neighbor_capacity * 2
                            new_neighbors = <uint32_t*>safe_realloc(
                                neighbors,
                                neighbor_capacity * sizeof(uint32_t),
                                new_cap * sizeof(uint32_t)
                            )
                            if not new_neighbors:
                                if verbose:
                                    bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Failed to grow neighbors\n")
                                free(read_seen_marker)
                                free(neighbors)
                                free(edge_from)
                                free(edge_to)
                                free(edge_weight)
                                free(ref_weights)
                                return -1
                            neighbors = new_neighbors
                            neighbor_capacity = new_cap

                        neighbors[neighbor_count] = other_ref_idx
                        neighbor_count += 1

                    ref_weights[other_ref_idx] += 1
                    read_seen_marker[other_ref_idx] = 0

        # Emit edges (no realloc - we allocated exact size)
        for kk in range(neighbor_count):
            other_ref_idx = neighbors[kk]

            if ref_weights[other_ref_idx] >= min_edge_weight and other_ref_idx > ref_idx:
                if edge_count >= edge_capacity:
                    if verbose:
                        bf_nogil_logf_notime(
                            b"IGRAPH OPS",
                            "ERROR: Edge overflow at ref %u (count=%llu, cap=%llu)\n",
                            ref_idx,
                            edge_count,
                            edge_capacity,
                        )
                    break

                edge_from[edge_count] = ref_idx
                edge_to[edge_count] = other_ref_idx
                edge_weight[edge_count] = ref_weights[other_ref_idx]
                edge_count += 1

            ref_weights[other_ref_idx] = 0

        # Free marker array for this reference
        free(read_seen_marker)

    if verbose:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "Built %llu/%llu edges (%.1f%% utilization)\n",
            edge_count,
            edge_capacity,
            100.0 * edge_count / edge_capacity,
        )

    # STEP 4: Create igraph
    ret = igraph_vector_int_init(&edges_vec, edge_count * 2)
    if ret != 0:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Failed to init edges vector\n")
        free(neighbors)
        free(edge_from)
        free(edge_to)
        free(edge_weight)
        free(ref_weights)
        return -1
    
    ret = igraph_vector_init(ig_weights, edge_count)
    if ret != 0:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Failed to init weights vector\n")
        igraph_vector_int_destroy(&edges_vec)
        free(neighbors)
        free(edge_from)
        free(edge_to)
        free(edge_weight)
        free(ref_weights)
        return -1
    
    # Copy edges to igraph
    for tmp_i in range(edge_count):
        igraph_vector_int_set(&edges_vec, tmp_i * 2, <igraph_integer_t>edge_from[tmp_i])
        igraph_vector_int_set(&edges_vec, tmp_i * 2 + 1, <igraph_integer_t>edge_to[tmp_i])
        igraph_vector_set(ig_weights, tmp_i, <igraph_real_t>edge_weight[tmp_i])
    
    ret = igraph_create(ig_graph, &edges_vec, <igraph_integer_t>num_refs, 0)
    igraph_vector_int_destroy(&edges_vec)
    
    if ret != 0:
        if verbose:
            bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: igraph_create failed\n")
        igraph_vector_destroy(ig_weights)
        free(neighbors)
        free(edge_from)
        free(edge_to)
        free(edge_weight)
        free(ref_weights)
        return -1

    success = True
    
    if verbose:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "igraph created: %ld vertices, %ld edges\n",
            <long>igraph_vcount(ig_graph),
            <long>igraph_ecount(ig_graph),
        )

    # STEP 5: Compute components
    if out_n_components or out_singletons:
        ret = igraph_vector_int_init(&comp_membership, num_refs)
        if ret == 0:
            comp_initialized = True
            ret = igraph_vector_int_init(&comp_sizes, 0)
            if ret == 0:
                comp_sizes_initialized = True
                ret = igraph_clusters(ig_graph, &comp_membership, &comp_sizes, &n_components, IGRAPH_WEAK)
                if ret == 0:
                    for i_comp in range(igraph_vector_int_size(&comp_sizes)):
                        if get_vector_int_element_local(&comp_sizes, i_comp) == 1:
                            singletons += 1
                    if out_n_components:
                        out_n_components[0] = n_components
                    if out_singletons:
                        out_singletons[0] = <igraph_integer_t>singletons
        
        if comp_sizes_initialized:
            igraph_vector_int_destroy(&comp_sizes)
        if comp_initialized:
            igraph_vector_int_destroy(&comp_membership)

    # Cleanup
    free(neighbors)
    free(edge_from)
    free(edge_to)
    free(edge_weight)
    free(ref_weights)
    
    return 0 if success else -1

cdef void destroy_weighted_graph(WeightedGraph* graph) noexcept nogil:
    """Free all memory associated with a WeightedGraph.

    Parameters
    ----------
    graph : WeightedGraph*
        Pointer to the weighted graph to free (safe to pass NULL).

    Notes
    -----
    This function frees per-node neighbor/weight buffers as well as any TSV
    arrays attached to the graph. The caller must ensure the graph is no
    longer used after calling this function.
    """
    cdef uint32_t i

    if not graph:
        return

    if graph.nodes:
        for i in range(graph.num_nodes):
            if graph.nodes[i].neighbors:
                free(graph.nodes[i].neighbors)
            if graph.nodes[i].weights:
                free(graph.nodes[i].weights)
        free(graph.nodes)

    # Free TSV data if present
    if graph.tsv_total_reads:
        free(graph.tsv_total_reads)
    if graph.tsv_multimap_reads:
        free(graph.tsv_multimap_reads)
    if graph.tsv_alignments_per_ref:
        free(graph.tsv_alignments_per_ref)
    if graph.tsv_exact_connection_counts:
        free(graph.tsv_exact_connection_counts)
    if graph.tsv_co_mapping_averages:
        free(graph.tsv_co_mapping_averages)
    if graph.tsv_max_co_mappings:
        free(graph.tsv_max_co_mappings)
    if graph.tsv_co_mapping_counts:
        free(graph.tsv_co_mapping_counts)
    if graph.tsv_neighbor_multimap_avg:
        free(graph.tsv_neighbor_multimap_avg)
    if graph.tsv_neighbor_connections_avg:
        free(graph.tsv_neighbor_connections_avg)
    if graph.tsv_neighbor_counts:
        free(graph.tsv_neighbor_counts)

    free(graph)


cdef int add_edge(WeightedGraph* graph, uint32_t from_node, uint32_t to_node,
                  uint32_t weight) noexcept nogil:
    """Insert or increment an undirected edge in a WeightedGraph.

    The function maintains sorted neighbor lists for each node and uses a
    binary search to find existing edges. If the edge exists the stored
    weight is incremented; otherwise the edge is inserted at the correct
    position to keep neighbors sorted.

    Parameters
    ----------
    graph : WeightedGraph*
        Pointer to the graph to modify (must be non-NULL).
    from_node : uint32_t
        Source node index.
    to_node : uint32_t
        Destination node index.
    weight : uint32_t
        Weight to add for the edge.

    Returns
    -------
    int
        0 on success, -1 on allocation failure or invalid node indices.

    Notes
    -----
    - Self-loops are ignored (no change).
    - Both directions are updated; ``graph.num_edges`` increments by 1 for a
      newly inserted undirected edge and ``graph.total_weight`` increases by
      twice the provided weight to reflect the undirected representation.
    - This function is ``nogil`` and must not call Python APIs.
    """
    cdef GraphNode* node_from
    cdef GraphNode* node_to
    cdef uint32_t i, insert_pos
    cdef bint found
    cdef uint32_t new_capacity
    cdef uint32_t* new_neighbors
    cdef uint32_t* new_weights
    cdef int low, high, mid
    
    if from_node >= graph.num_nodes or to_node >= graph.num_nodes:
        return -1
    
    if from_node == to_node:
        return 0  # No self-loops
    
    # Add edge from -> to
    node_from = &graph.nodes[from_node]
    
    # Binary search for existing edge (neighbors are kept sorted)
    found = False
    insert_pos = node_from.degree
    low = 0
    high = <int>node_from.degree - 1
    
    while low <= high:
        mid = (low + high) // 2
        if node_from.neighbors[mid] == to_node:
            # Edge exists, increment weight
            node_from.weights[mid] += weight
            node_from.node_weight += weight
            found = True
            break
        elif node_from.neighbors[mid] < to_node:
            low = mid + 1
            insert_pos = low
        else:
            high = mid - 1
            insert_pos = low
    
    if not found:
        # Need to add new edge at insert_pos
        if node_from.degree >= node_from.capacity:
            # Expand capacity - use larger increments for efficiency
            new_capacity = node_from.capacity * 2 if node_from.capacity > 0 else 8
            
            new_neighbors = <uint32_t*>realloc(node_from.neighbors, 
                                              new_capacity * sizeof(uint32_t))
            new_weights = <uint32_t*>realloc(node_from.weights,
                                            new_capacity * sizeof(uint32_t))
            
            if not new_neighbors or not new_weights:
                return -1
            
            node_from.neighbors = new_neighbors
            node_from.weights = new_weights
            node_from.capacity = new_capacity
        
        # Shift elements to make room at insert_pos (maintain sorted order)
        for i in range(node_from.degree, insert_pos, -1):
            node_from.neighbors[i] = node_from.neighbors[i-1]
            node_from.weights[i] = node_from.weights[i-1]
        
        node_from.neighbors[insert_pos] = to_node
        node_from.weights[insert_pos] = weight
        node_from.degree += 1
        node_from.node_weight += weight
    
    # Add edge to -> from (undirected graph)
    node_to = &graph.nodes[to_node]
    
    # Binary search again for reverse edge
    found = False
    insert_pos = node_to.degree
    low = 0
    high = <int>node_to.degree - 1
    
    while low <= high:
        mid = (low + high) // 2
        if node_to.neighbors[mid] == from_node:
            node_to.weights[mid] += weight
            node_to.node_weight += weight
            found = True
            break
        elif node_to.neighbors[mid] < from_node:
            low = mid + 1
            insert_pos = low
        else:
            high = mid - 1
            insert_pos = low
    
    if not found:
        if node_to.degree >= node_to.capacity:
            new_capacity = node_to.capacity * 2 if node_to.capacity > 0 else 8
            
            new_neighbors = <uint32_t*>realloc(node_to.neighbors,
                                              new_capacity * sizeof(uint32_t))
            new_weights = <uint32_t*>realloc(node_to.weights,
                                            new_capacity * sizeof(uint32_t))
            
            if not new_neighbors or not new_weights:
                return -1
            
            node_to.neighbors = new_neighbors
            node_to.weights = new_weights
            node_to.capacity = new_capacity
        
        # Shift elements to maintain sorted order
        for i in range(node_to.degree, insert_pos, -1):
            node_to.neighbors[i] = node_to.neighbors[i-1]
            node_to.weights[i] = node_to.weights[i-1]
        
        node_to.neighbors[insert_pos] = from_node
        node_to.weights[insert_pos] = weight
        node_to.degree += 1
        node_to.node_weight += weight
    
    # Update total graph weight (edge counted twice for undirected)
    graph.total_weight += 2 * weight
    graph.num_edges += 1
    
    return 0


# ==============================================================================
# DIRECT GRAPH BUILDING FROM READINDEX (ULTRA-FAST, ZERO EXTRA MEMORY)
# ==============================================================================

# Helper structure for tracking seen references per read
# Used to avoid counting the same edge pair multiple times per read
cdef struct SeenTracker:
    char* seen_flags  # bitmap: 1 = already processed this ref for current read
    uint32_t capacity

cdef SeenTracker* create_seen_tracker(uint32_t capacity) except NULL nogil:
    """Create a SeenTracker used to mark references seen within a read.

    The SeenTracker provides a compact byte-array bitmap used to mark which
    reference indices have already been seen while scanning a single read.
    It is intended for reuse across many reads to avoid repeated allocations
    in tight, nogil code paths.

    Parameters
    ----------
    capacity : uint32_t
        Number of reference indices the tracker must support (typically
        equal to ``num_refs``).

    Returns
    -------
    SeenTracker*
        Pointer to an allocated SeenTracker, or ``NULL`` on allocation failure.
    """
    cdef SeenTracker* tracker = <SeenTracker*>malloc(sizeof(SeenTracker))
    if not tracker:
        return NULL
    
    tracker.capacity = capacity
    tracker.seen_flags = <char*>calloc(capacity, sizeof(char))
    if not tracker.seen_flags:
        free(tracker)
        return NULL
    
    return tracker

cdef void destroy_seen_tracker(SeenTracker* tracker) noexcept nogil:
    """Free resources held by a SeenTracker.

    Parameters
    ----------
    tracker : SeenTracker*
        Pointer to the tracker to free (safe to pass ``NULL``).

    Notes
    -----
    This frees the internal bitmap and the tracker struct itself.
    """
    if not tracker:
        return
    if tracker.seen_flags:
        free(tracker.seen_flags)
    free(tracker)

cdef inline void reset_seen_tracker(SeenTracker* tracker) noexcept nogil:
    """Reset all flags in the SeenTracker to zero.

    Parameters
    ----------
    tracker : SeenTracker*
        Tracker to reset. Must be non-NULL.
    """
    memset(tracker.seen_flags, 0, tracker.capacity * sizeof(char))


# ==============================================================================
# GRAPH BUILDING FROM ALIGNMENT DATA
# ==============================================================================

cdef WeightedGraph* build_weighted_graph_from_alignments(
    MemoryPool* pool,
    ReferenceStats* ref_stats,
    ReadIndex* read_index,
    uint32_t num_refs,
    int32_t min_read_count
) except NULL nogil:
    """Build a WeightedGraph directly from alignment data using ReadIndex.

    This routine constructs a full weighted, undirected graph of references
    where an edge weight is the number of unique reads shared between two
    references. The builder avoids large intermediate allocations by
    pre-counting unique neighbor counts and allocating exact per-node
    buffers, then populating neighbor lists in a second pass.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignment records and per-read indexing arrays.
    ref_stats : ReferenceStats*
        Optional per-reference statistics used to skip low-coverage
        references (may be ``NULL``).
    read_index : ReadIndex*
        ReadIndex mapping references to lists of read indices.
    num_refs : uint32_t
        Number of references.
    min_read_count : int32_t
        Minimum reads required for a reference to be considered.

    Returns
    -------
    WeightedGraph*
        Pointer to an allocated WeightedGraph on success, or ``NULL`` on
        allocation or validation failure. The returned graph must be freed
        with :func:`destroy_weighted_graph`.
    """
    cdef uint32_t ref_idx, other_ref_idx
    cdef uint32_t read_i, read_count
    cdef uint64_t read_idx
    cdef uint32_t* read_list
    cdef uint64_t aln_start_i, aln_end_i, aln_idx_i
    cdef WeightedGraph* graph = NULL
    cdef uint32_t* degree_counts = NULL
    cdef uint64_t reads_processed = 0
    cdef uint64_t total_weight = 0
    cdef uint64_t total_degree = 0
    bf_nogil_logf_notime(b"IGRAPH OPS", "Building weighted graph from ReadIndex (ZERO extra memory)...\n")
    bf_nogil_logf_notime(b"IGRAPH OPS", "  References: %u\n", num_refs)
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Unique reads: %u\n", pool.unique_read_count)

    # Phase 1: Pre-count ALL neighbor occurrences (not unique!)
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Phase 1: Pre-counting neighbor occurrences from ReadIndex...\n")
    degree_counts = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    if not degree_counts:
        return NULL
    
    # For each reference, count ALL neighbor occurrences across all reads
    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue
        
        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0:
            continue
        
        read_list = read_index.ref_to_reads[ref_idx]
        
        # For each read that maps to this reference
        for read_i in range(read_count):
            read_idx = read_list[read_i]
            if read_idx >= pool.unique_read_count:
                continue
            
            aln_start_i = pool.read_alignment_starts[read_idx]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx]
            
            # Count ALL other references this read maps to
            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index
                
                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue
                
                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue
                
                # Count every occurrence
                degree_counts[ref_idx] += 1
        
        reads_processed += read_count
    
    bf_nogil_logf_notime(b"IGRAPH OPS", "    Processed %llu read mappings\n", <unsigned long long>reads_processed)
    
    # Phase 2: Allocate exact graph capacity based on UNIQUE neighbors
    # We need to compress degree_counts to unique neighbors
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Phase 2: Calculating unique neighbors and allocating graph...\n")
    graph = create_weighted_graph(num_refs)
    if not graph:
        free(degree_counts)
        return NULL
    
    # We'll allocate during phase 3 as we discover unique neighbors
    # (degree_counts tells us max possible, but actual unique count is lower)
    
    # Phase 3: Fill graph with weighted edges
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Phase 3: Filling graph with weighted edges...\n")
    
    cdef uint32_t* edge_weights = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    if not edge_weights:
        destroy_weighted_graph(graph)
        free(degree_counts)
        return NULL
    
    cdef uint32_t* neighbor_list = NULL
    cdef uint32_t neighbor_len = 0
    cdef uint32_t neighbor_capacity = 0
    cdef uint32_t kk
    cdef uint32_t w
    cdef uint32_t pos, pos2
    # Growth-capacity temporaries declared at function scope (Cython forbids cdef in inner blocks)
    cdef uint32_t alloc_cap
    cdef uint32_t old_cap, new_cap
    cdef uint32_t old_cap2, new_cap2
    cdef uint32_t* write_pos = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    cdef uint32_t* unique_neighbor_counts = <uint32_t*>calloc(num_refs, sizeof(uint32_t))
    cdef uint8_t* read_seen_refs_global = NULL
    cdef uint8_t* read_seen_refs_local = NULL

    if not write_pos or not unique_neighbor_counts:
        free(edge_weights)
        edge_weights = NULL
        if write_pos:
            free(write_pos)
            write_pos = NULL
        if unique_neighbor_counts:
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
        destroy_weighted_graph(graph)
        free(degree_counts)
        degree_counts = NULL
        return NULL

    # First pass: count unique neighbors per reference (count each shared read once)
    read_seen_refs_global = <uint8_t*>calloc(num_refs, sizeof(uint8_t))
    if not read_seen_refs_global:
        free(edge_weights)
        edge_weights = NULL
        if write_pos:
            free(write_pos)
            write_pos = NULL
        if unique_neighbor_counts:
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
        destroy_weighted_graph(graph)
        free(degree_counts)
        degree_counts = NULL
        return NULL

    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue

        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0:
            continue

        read_list = read_index.ref_to_reads[ref_idx]

        for read_i in range(read_count):
            read_idx = read_list[read_i]
            if read_idx >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx]

            # Mark which references this read maps to
            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue

                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                read_seen_refs_global[other_ref_idx] = 1

            # Count and increment weights (once per read)
            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue

                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                if read_seen_refs_global[other_ref_idx] == 1:
                    if edge_weights[other_ref_idx] == 0:
                        unique_neighbor_counts[ref_idx] += 1
                    edge_weights[other_ref_idx] += 1
                    read_seen_refs_global[other_ref_idx] = 0

        # Clear edge_weights for this ref
        for read_i in range(read_count):
            read_idx = read_list[read_i]
            if read_idx >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx]

            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index
                if other_ref_idx < num_refs:
                    edge_weights[other_ref_idx] = 0

    free(read_seen_refs_global)
    
    # Allocate neighbor arrays. Ensure each node has a non-NULL buffer (capacity>=1)
    # so we never write to a NULL pointer when emitting the reverse edge.
    for ref_idx in range(num_refs):
        if unique_neighbor_counts[ref_idx] > 0:
            alloc_cap = unique_neighbor_counts[ref_idx]
        else:
            alloc_cap = 1
        graph.nodes[ref_idx].capacity = alloc_cap
        graph.nodes[ref_idx].neighbors = <uint32_t*>malloc(alloc_cap * sizeof(uint32_t))
        graph.nodes[ref_idx].weights = <uint32_t*>malloc(alloc_cap * sizeof(uint32_t))

        if not graph.nodes[ref_idx].neighbors or not graph.nodes[ref_idx].weights:
            free(edge_weights)
            edge_weights = NULL
            free(write_pos)
            write_pos = NULL
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
            destroy_weighted_graph(graph)
            free(degree_counts)
            degree_counts = NULL
            return NULL
    
    # Second pass: fill graph
    for ref_idx in range(num_refs):
        if ref_stats and ref_stats[ref_idx].total_reads < min_read_count:
            continue
        
        read_count = read_index.ref_read_counts[ref_idx]
        if read_count == 0:
            continue
        
        read_list = read_index.ref_to_reads[ref_idx]
        neighbor_len = 0
        
        if not neighbor_list:
            free(edge_weights)
            edge_weights = NULL
            free(write_pos)
            write_pos = NULL
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
            destroy_weighted_graph(graph)
            free(degree_counts)
            degree_counts = NULL
            return NULL

        # Temporary marker for this ref's reads
        read_seen_refs_local = <uint8_t*>calloc(num_refs, sizeof(uint8_t))
        if not read_seen_refs_local:
            if neighbor_list:
                free(neighbor_list)
            free(edge_weights)
            edge_weights = NULL
            free(write_pos)
            write_pos = NULL
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
            destroy_weighted_graph(graph)
            free(degree_counts)
            degree_counts = NULL
            return NULL

        # Accumulate weights (count each read once per reference)
        for read_i in range(read_count):
            read_idx = read_list[read_i]
            if read_idx >= pool.unique_read_count:
                continue

            aln_start_i = pool.read_alignment_starts[read_idx]
            aln_end_i = aln_start_i + pool.read_alignment_counts[read_idx]

            # Mark which references this read maps to
            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue

                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                read_seen_refs_local[other_ref_idx] = 1

            # Increment weights once per unique reference for this read
            for aln_idx_i in range(aln_start_i, aln_end_i):
                other_ref_idx = pool.alignments[aln_idx_i].reference_index

                if other_ref_idx >= num_refs or other_ref_idx == ref_idx:
                    continue

                if ref_stats and ref_stats[other_ref_idx].total_reads < min_read_count:
                    continue

                if read_seen_refs_local[other_ref_idx] == 1:
                    if edge_weights[other_ref_idx] == 0:
                        neighbor_list[neighbor_len] = other_ref_idx
                        neighbor_len += 1
                    edge_weights[other_ref_idx] += 1
                    read_seen_refs_local[other_ref_idx] = 0
        
        # Emit edges (with bounds checks and safe growth if necessary)
        for kk in range(neighbor_len):
            other_ref_idx = neighbor_list[kk]
            w = edge_weights[other_ref_idx]

            if w > 0 and other_ref_idx > ref_idx:
                # Ensure space for ref_idx
                pos = write_pos[ref_idx]
                if pos >= graph.nodes[ref_idx].capacity:
                    # grow capacity
                    old_cap = graph.nodes[ref_idx].capacity
                    new_cap = old_cap * 2 if old_cap > 0 else 2
                    graph.nodes[ref_idx].neighbors = <uint32_t*>safe_realloc(
                        graph.nodes[ref_idx].neighbors,
                        old_cap * sizeof(uint32_t),
                        new_cap * sizeof(uint32_t)
                    )
                    graph.nodes[ref_idx].weights = <uint32_t*>safe_realloc(
                        graph.nodes[ref_idx].weights,
                        old_cap * sizeof(uint32_t),
                        new_cap * sizeof(uint32_t)
                    )
                    if not graph.nodes[ref_idx].neighbors or not graph.nodes[ref_idx].weights:
                        if graph.nodes[ref_idx].neighbors:
                            free(graph.nodes[ref_idx].neighbors)
                            graph.nodes[ref_idx].neighbors = NULL
                        if graph.nodes[ref_idx].weights:
                            free(graph.nodes[ref_idx].weights)
                            graph.nodes[ref_idx].weights = NULL
                        free(neighbor_list)
                        neighbor_list = NULL
                        free(edge_weights)
                        edge_weights = NULL
                        free(write_pos)
                        write_pos = NULL
                        free(unique_neighbor_counts)
                        unique_neighbor_counts = NULL
                        destroy_weighted_graph(graph)
                        free(degree_counts)
                        degree_counts = NULL
                        return NULL
                    graph.nodes[ref_idx].capacity = new_cap

                graph.nodes[ref_idx].neighbors[pos] = other_ref_idx
                graph.nodes[ref_idx].weights[pos] = w
                graph.nodes[ref_idx].degree += 1
                graph.nodes[ref_idx].node_weight += w
                write_pos[ref_idx] = pos + 1

                # Ensure space for other_ref_idx (reverse entry)
                pos2 = write_pos[other_ref_idx]
                if pos2 >= graph.nodes[other_ref_idx].capacity:
                    old_cap2 = graph.nodes[other_ref_idx].capacity
                    new_cap2 = old_cap2 * 2 if old_cap2 > 0 else 2
                    graph.nodes[other_ref_idx].neighbors = <uint32_t*>safe_realloc(
                        graph.nodes[other_ref_idx].neighbors,
                        old_cap2 * sizeof(uint32_t),
                        new_cap2 * sizeof(uint32_t)
                    )
                    graph.nodes[other_ref_idx].weights = <uint32_t*>safe_realloc(
                        graph.nodes[other_ref_idx].weights,
                        old_cap2 * sizeof(uint32_t),
                        new_cap2 * sizeof(uint32_t)
                    )
                    if not graph.nodes[other_ref_idx].neighbors or not graph.nodes[other_ref_idx].weights:
                        if graph.nodes[other_ref_idx].neighbors:
                            free(graph.nodes[other_ref_idx].neighbors)
                            graph.nodes[other_ref_idx].neighbors = NULL
                        if graph.nodes[other_ref_idx].weights:
                            free(graph.nodes[other_ref_idx].weights)
                            graph.nodes[other_ref_idx].weights = NULL
                        free(neighbor_list)
                        neighbor_list = NULL
                        free(edge_weights)
                        edge_weights = NULL
                        free(write_pos)
                        write_pos = NULL
                        free(unique_neighbor_counts)
                        unique_neighbor_counts = NULL
                        destroy_weighted_graph(graph)
                        free(degree_counts)
                        degree_counts = NULL
                        return NULL
                    graph.nodes[other_ref_idx].capacity = new_cap2

                graph.nodes[other_ref_idx].neighbors[pos2] = ref_idx
                graph.nodes[other_ref_idx].weights[pos2] = w
                graph.nodes[other_ref_idx].degree += 1
                graph.nodes[other_ref_idx].node_weight += w
                write_pos[other_ref_idx] = pos2 + 1

                graph.num_edges += 1
                graph.total_weight += 2 * w

            # reset accumulator
            edge_weights[other_ref_idx] = 0

        # Free temporary marker array for this reference
        free(read_seen_refs_local)

        if neighbor_list:
            free(neighbor_list)
            neighbor_list = NULL
        
        if ((ref_idx + 1) % 1000 == 0):
            bf_nogil_logf_notime(
                b"IGRAPH OPS",
                "  Processed %u/%u refs, %llu edges added\n",
                ref_idx + 1,
                num_refs,
                <unsigned long long>graph.num_edges,
            )

    bf_nogil_logf_notime(b"IGRAPH OPS", "\n  Graph building complete: %llu edges\n", <unsigned long long>graph.num_edges)
    
    if edge_weights:
        free(edge_weights)
        edge_weights = NULL
    if write_pos:
        free(write_pos)
        write_pos = NULL
    if unique_neighbor_counts:
        free(unique_neighbor_counts)
        unique_neighbor_counts = NULL
    if degree_counts:
        free(degree_counts)
        degree_counts = NULL
    
    # Calculate statistics
    total_weight = 0
    total_degree = 0
    cdef uint32_t min_weight = UINT32_MAX
    cdef uint32_t max_weight = 0
    cdef uint64_t weight_sum = 0
    cdef uint64_t edge_count = 0
    cdef uint32_t edge_weight
    cdef uint32_t edge_i
    
    for ref_idx in range(num_refs):
        total_degree += graph.nodes[ref_idx].degree
        total_weight += graph.nodes[ref_idx].node_weight
        
        for edge_i in range(graph.nodes[ref_idx].degree):
            edge_weight = graph.nodes[ref_idx].weights[edge_i]
            if edge_weight < min_weight:
                min_weight = edge_weight
            if edge_weight > max_weight:
                max_weight = edge_weight
            weight_sum += edge_weight
            edge_count += 1
    
    graph.total_weight = total_weight / 2
    
    bf_nogil_logf_notime(b"IGRAPH OPS", "\n  === GRAPH STATISTICS (BEFORE PRUNING) ===\n")
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Total edges: %llu\n", <unsigned long long>graph.num_edges)
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Total edge weight: %llu\n", <unsigned long long>graph.total_weight)
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Average node degree: %.1f\n", <double>total_degree / <double>num_refs)

    if edge_count > 0:
        bf_nogil_logf_notime(b"IGRAPH OPS", "  Edge weight distribution:\n")
        bf_nogil_logf_notime(b"IGRAPH OPS", "    Min: %u\n", min_weight)
        bf_nogil_logf_notime(b"IGRAPH OPS", "    Max: %u\n", max_weight)
        bf_nogil_logf_notime(b"IGRAPH OPS", "    Mean: %.1f\n", <double>weight_sum / <double>edge_count)
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "    Total edges examined: %llu (x2 for undirected)\n",
            <unsigned long long>edge_count,
        )

    # Validation pass: ensure neighbor ids are within bounds and weights > 0
    cdef uint32_t validation_errors = 0
    cdef uint64_t degree_sum = 0
    for ref_idx in range(num_refs):
        degree_sum += graph.nodes[ref_idx].degree
        for edge_i in range(graph.nodes[ref_idx].degree):
            if graph.nodes[ref_idx].neighbors[edge_i] >= num_refs:
                validation_errors += 1
            if graph.nodes[ref_idx].weights[edge_i] == 0:
                validation_errors += 1

    if validation_errors > 0:
        bf_nogil_logf_notime(b"IGRAPH OPS", "ERROR: Graph validation failed with %u errors\n", validation_errors)
        destroy_weighted_graph(graph)
        if edge_weights:
            free(edge_weights)
            edge_weights = NULL
        if write_pos:
            free(write_pos)
            write_pos = NULL
        if unique_neighbor_counts:
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
        if degree_counts:
            free(degree_counts)
            degree_counts = NULL
        return NULL

    if degree_sum != 2 * graph.num_edges:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "ERROR: Degree sum mismatch: %llu != %llu (2*num_edges)\n",
            degree_sum,
            <unsigned long long>(2 * graph.num_edges),
        )
        destroy_weighted_graph(graph)
        if edge_weights:
            free(edge_weights)
            edge_weights = NULL
        if write_pos:
            free(write_pos)
            write_pos = NULL
        if unique_neighbor_counts:
            free(unique_neighbor_counts)
            unique_neighbor_counts = NULL
        if degree_counts:
            free(degree_counts)
            degree_counts = NULL
        return NULL

    return graph

# ==============================================================================
# GRAPH FILTERING
# ==============================================================================

cdef void prune_low_weight_edges(WeightedGraph* graph, uint32_t min_edge_weight) noexcept nogil:
    """Prune edges whose weight is below a minimum threshold.

    This function compacts each node's neighbor and weight arrays in-place,
    removing edges with weight < ``min_edge_weight``. It updates per-node
    degrees and node weights and recomputes ``graph.num_edges`` accordingly.

    Parameters
    ----------
    graph : WeightedGraph*
        Graph to prune. If ``graph`` is ``NULL`` the function returns
        immediately.
    min_edge_weight : uint32_t
        Exclusive lower bound for keeping an edge; edges with
        weight < ``min_edge_weight`` are removed. Values <= 1 are treated
        as a no-op (no pruning performed).

    Notes
    -----
    - The function modifies the graph in-place and is ``nogil``.
    - ``graph.num_edges`` is recomputed from the remaining edges and
      represents the number of undirected edges after pruning.
    """
    cdef uint32_t node_idx, i, new_degree
    cdef GraphNode* node
    cdef uint64_t edges_before = graph.num_edges
    cdef uint64_t edges_removed = 0
    
    if not graph or min_edge_weight <= 1:
        return
    
    bf_nogil_logf_notime(b"IGRAPH OPS", "  Pruning edges with weight < %u...\n", min_edge_weight)
    
    for node_idx in range(graph.num_nodes):
        node = &graph.nodes[node_idx]
        new_degree = 0
        
        # Compact the neighbor/weight arrays, keeping only edges >= min_weight
        for i in range(node.degree):
            if node.weights[i] >= min_edge_weight:
                if new_degree != i:
                    node.neighbors[new_degree] = node.neighbors[i]
                    node.weights[new_degree] = node.weights[i]
                new_degree += 1
            else:
                # Edge is being removed, update node weight
                node.node_weight -= node.weights[i]
        
        # Update degree
        edges_removed += (node.degree - new_degree)
        node.degree = new_degree

    graph.num_edges = edges_before - (edges_removed / 2)
    
    bf_nogil_logf_notime(
        b"IGRAPH OPS",
        "  Removed %lu edges, %lu edges remaining\n",
        edges_removed / 2,
        graph.num_edges,
    )


# ==============================================================================
# GRAPH STATISTICS CALCULATION
# ==============================================================================

cdef void calculate_graph_statistics(WeightedGraph* graph,
                                     MemoryPool* pool,
                                     ReferencePattern* pattern_data,
                                     uint32_t* exact_connection_counts,
                                     double* neighbor_multimap_avg,
                                     double* neighbor_connections_avg,
                                     uint32_t* neighbor_counts,
                                     float* dataset_multimap_fractions,
                                     ReferenceStats* ref_stats,
                                     uint32_t num_refs) noexcept nogil:
    """Compute per-reference and neighbor-level statistics from a graph.

    This routine computes a set of statistics derived from the provided
    weighted graph and populates the supplied output arrays. The function
    assumes the graph has already been filtered (for example via
    :func:`prune_low_weight_edges`) and computes statistics consistently from
    that filtered state.

    Parameters
    ----------
    graph : WeightedGraph*
        Filtered weighted graph to analyze.
    pool : MemoryPool*
        Memory pool containing alignment information used by some stats.
    pattern_data : ReferencePattern*
        Per-reference pattern accumulator used to compute co-mapping
        averages and maxima (writer-supplied).
    exact_connection_counts : uint32_t*
        Output array (length ``num_refs``) receiving exact connection counts.
    neighbor_multimap_avg : double*
        Output array receiving average neighbor multimap rates per reference.
    neighbor_connections_avg : double*
        Output array receiving average neighbor connection counts per ref.
    neighbor_counts : uint32_t*
        Output array receiving number of neighbors per reference.
    dataset_multimap_fractions : float*
        Optional per-dataset multimap fraction buffer used when computing
        neighbor-level multimap metrics (may be ``NULL``).
    ref_stats : ReferenceStats*
        Optional per-reference stats used to mask or scale statistics
        (may be ``NULL``).
    num_refs : uint32_t
        Number of references / length of output arrays.

    Notes
    -----
    - All outputs are written in-place to caller-provided buffers. The
      caller must allocate arrays with length at least ``num_refs``.
    - This function is ``nogil`` and must not call Python APIs.
    """
    cdef uint32_t ref_idx, i, neighbor_ref
    cdef uint32_t connection_count
    cdef uint64_t total_weight
    cdef double neighbor_multimap_sum, neighbor_connections_sum
    cdef uint32_t max_comapping
    
    bf_nogil_logf_notime(
        b"IGRAPH OPS",
        "Calculating graph statistics from filtered graph (%lu edges)...\n",
        graph.num_edges,
    )
    
    for ref_idx in range(num_refs):
        # Get connection count from graph degree
        connection_count = graph.nodes[ref_idx].degree
        exact_connection_counts[ref_idx] = connection_count
        pattern_data[ref_idx].graph.connection_count = connection_count
        
        if connection_count == 0:
            pattern_data[ref_idx].graph.avg_comappings_per_read = 0.0
            pattern_data[ref_idx].graph.max_comappings = 0
            pattern_data[ref_idx].graph.reads_with_comappings = 0
            neighbor_multimap_avg[ref_idx] = 0.0
            neighbor_connections_avg[ref_idx] = 0.0
            neighbor_counts[ref_idx] = 0
            continue
        
        # Calculate average co-mappings per read
        # This is the average degree (connections - 1) for reads mapping to this ref
        total_weight = 0
        max_comapping = 0
        for i in range(connection_count):
            # Weight = number of reads shared with this neighbor
            total_weight += graph.nodes[ref_idx].weights[i]
            if graph.nodes[ref_idx].weights[i] > max_comapping:
                max_comapping = graph.nodes[ref_idx].weights[i]
        
        # Average co-mappings = average number of OTHER refs each read maps to
        # This is approximately total_weight / num_reads_for_this_ref
        if ref_stats and ref_stats[ref_idx].total_reads > 0:
            pattern_data[ref_idx].graph.avg_comappings_per_read = \
                <double>total_weight / <double>ref_stats[ref_idx].total_reads
            pattern_data[ref_idx].graph.reads_with_comappings = ref_stats[ref_idx].total_reads
        else:
            pattern_data[ref_idx].graph.avg_comappings_per_read = <double>connection_count
            pattern_data[ref_idx].graph.reads_with_comappings = 0
        
        # max_comapping is the largest edge weight observed for this reference
        # (i.e. the maximum number of reads shared with a single neighbor).
        # Previously this was incorrectly set to `connection_count` which
        # inflated aggregate statistics. Use the computed `max_comapping`.
        pattern_data[ref_idx].graph.max_comappings = max_comapping
        
        # Calculate neighbor statistics
        neighbor_multimap_sum = 0.0
        neighbor_connections_sum = 0.0
        
        for i in range(connection_count):
            neighbor_ref = graph.nodes[ref_idx].neighbors[i]
            if neighbor_ref < num_refs:
                # Add neighbor's multimap fraction
                if dataset_multimap_fractions:
                    neighbor_multimap_sum += dataset_multimap_fractions[neighbor_ref]
                # Add neighbor's connection count
                neighbor_connections_sum += <double>graph.nodes[neighbor_ref].degree
        
        # Average neighbor stats
        if connection_count > 0:
            neighbor_multimap_avg[ref_idx] = neighbor_multimap_sum / <double>connection_count
            neighbor_connections_avg[ref_idx] = neighbor_connections_sum / <double>connection_count
        else:
            neighbor_multimap_avg[ref_idx] = 0.0
            neighbor_connections_avg[ref_idx] = 0.0
        
        neighbor_counts[ref_idx] = connection_count
    
    bf_nogil_logf_notime(b"IGRAPH OPS", "Graph statistics calculated for %u references\n", num_refs)


cdef int build_igraph_from_weighted_graph(
    igraph_t* ig_graph,
    igraph_vector_t* ig_weights,
    WeightedGraph* graph,
    bint verbose,
    igraph_integer_t* out_n_components,
    igraph_integer_t* out_singletons
) noexcept nogil:
    """Build an igraph representation from an existing WeightedGraph.

    This function converts the in-memory :class:`WeightedGraph` representation
    into an ``igraph_t`` and an associated weight vector. It emits each
    undirected edge once (i < neighbor) into the igraph edge list. Optionally
    computes the number of connected components and singleton components.

    Parameters
    ----------
    ig_graph : igraph_t*
        Caller-provided storage for the resulting igraph object (initialized
        by this function on success).
    ig_weights : igraph_vector_t*
        Caller-provided storage for the igraph edge weight vector.
    graph : WeightedGraph*
        Input weighted graph to convert (may be pruned already).
    verbose : bint
        When true, print progress and statistics to stdout.
    out_n_components : igraph_integer_t*
        Optional output pointer that receives the number of connected
        components (may be ``NULL``).
    out_singletons : igraph_integer_t*
        Optional output pointer that receives the number of singleton
        components (may be ``NULL``).

    Returns
    -------
    int
        0 on success, -1 on failure. On success the caller-owned igraph
        containers are initialized and must be destroyed by the caller.
    """
    cdef igraph_vector_int_t edges
    cdef uint32_t i, j, neighbor_idx
    cdef uint32_t neighbor_id, edge_weight
    cdef uint32_t edge_count = 0
    cdef int ret
    
    # Count edges
    for i in range(graph.num_nodes):
        edge_count += graph.nodes[i].degree
    edge_count = edge_count / 2  # Undirected
    
    if verbose:
        bf_nogil_logf_notime(
            b"IGRAPH OPS",
            "Building igraph from WeightedGraph: %u nodes, %u edges\n",
            graph.num_nodes,
            edge_count,
        )
    
    # Initialize edge list and weights
    ret = igraph_vector_int_init(&edges, edge_count * 2)
    if ret != IGRAPH_SUCCESS:
        return -1
    
    ret = igraph_vector_init(ig_weights, edge_count)
    if ret != IGRAPH_SUCCESS:
        igraph_vector_int_destroy(&edges)
        return -1
    
    # Populate edge list and weights
    cdef uint32_t edge_idx = 0
    for i in range(graph.num_nodes):
        for neighbor_idx in range(graph.nodes[i].degree):
            neighbor_id = graph.nodes[i].neighbors[neighbor_idx]
            edge_weight = graph.nodes[i].weights[neighbor_idx]
            
            # Add edge only once (i < neighbor_id)
            if i < neighbor_id:
                igraph_vector_int_set(&edges, edge_idx * 2, <igraph_integer_t>i)
                igraph_vector_int_set(&edges, edge_idx * 2 + 1, <igraph_integer_t>neighbor_id)
                igraph_vector_set(ig_weights, edge_idx, <igraph_real_t>edge_weight)
                edge_idx += 1
    
    # Create igraph
    ret = igraph_create(ig_graph, &edges, <igraph_integer_t>graph.num_nodes, IGRAPH_UNDIRECTED)
    igraph_vector_int_destroy(&edges)
    
    if ret != IGRAPH_SUCCESS:
        igraph_vector_destroy(ig_weights)
        return -1
    
    # Count connected components
    cdef igraph_vector_int_t membership
    cdef igraph_vector_int_t csize
    cdef igraph_integer_t no_components
    # Declare loop counters and accumulators at function scope to satisfy Cython
    cdef igraph_integer_t singletons
    cdef igraph_integer_t comp_idx
    
    # Compute components using igraph_components (replacement for deprecated igraph_clusters)
    ret = igraph_vector_int_init(&membership, graph.num_nodes)
    if ret == IGRAPH_SUCCESS:
        ret = igraph_vector_int_init(&csize, 0)
        if ret == IGRAPH_SUCCESS:
            # igraph_clusters (aka igraph_components in some igraph versions) returns
            # component sizes in csize and membership mapping. Use the older
            # `igraph_clusters` symbol which is available in the igraph C API
            # versions this code is built against.
            ret = igraph_clusters(ig_graph, &membership, &csize, &no_components, IGRAPH_WEAK)
            if ret == IGRAPH_SUCCESS:
                out_n_components[0] = no_components

                # Count singletons (components of size 1)
                singletons = 0
                for comp_idx in range(igraph_vector_int_size(&csize)):
                    if get_vector_int_element_local(&csize, comp_idx) == 1:
                        singletons += 1
                out_singletons[0] = singletons

            igraph_vector_int_destroy(&csize)
        igraph_vector_int_destroy(&membership)
    
    return 0
