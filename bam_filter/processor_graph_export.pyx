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
Graph export functionality for read-reference co-mapping networks.

Exports graphs to GraphML format with comprehensive node and edge attributes
for visualization and analysis in Cytoscape, igraph, Gephi, NetworkX, etc.
"""

from libc.stdio cimport FILE, fopen, fclose, fprintf, snprintf
from libc.string cimport strlen, strcpy
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReferencePattern, WeightedGraph
from bam_filter.processor_graph_ops cimport GraphNode
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor_igraph cimport (
    igraph_t, igraph_error_t, IGRAPH_SUCCESS,
    igraph_integer_t, igraph_vector_t, igraph_vector_int_t, igraph_real_t,
    igraph_vs_t, igraph_neimode_t, IGRAPH_ALL,
    igraph_ecount, igraph_edge
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t:
        int32_t n_targets
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

cdef extern from "igraph.h":
    ctypedef bint igraph_bool_t
    igraph_error_t igraph_write_graph_graphml(const igraph_t *graph, FILE *outstream,
                                              igraph_bool_t prefixattr) nogil

    # VECTOR macro for accessing vector elements
    igraph_real_t* VECTOR(igraph_vector_t v) nogil

    # Vector functions
    igraph_error_t igraph_vector_int_init(igraph_vector_int_t* v, igraph_integer_t size) nogil
    void igraph_vector_int_destroy(igraph_vector_int_t* v) nogil
    igraph_integer_t igraph_vector_int_size(const igraph_vector_int_t* v) nogil

    # Vertex selector functions
    igraph_error_t igraph_vs_all(igraph_vs_t* vs) nogil
    void igraph_vs_destroy(igraph_vs_t* vs) nogil

    # Graph query functions
    igraph_error_t igraph_degree(const igraph_t *graph, igraph_vector_int_t *res,
                                 igraph_vs_t vids, igraph_neimode_t mode,
                                 igraph_bool_t loops) nogil


cdef inline igraph_integer_t get_vector_int_element(igraph_vector_int_t *v, igraph_integer_t i) nogil:
    """Access element from igraph_vector_int_t."""
    return (<igraph_integer_t**>v)[0][i]


cdef uint32_t UINT32_MAX = 0xFFFFFFFF


cdef inline void xml_escape(const char* src, char* dst, uint32_t dst_size) noexcept nogil:
    """Escape XML special characters for safe GraphML output."""
    cdef uint32_t i = 0, j = 0
    cdef char c

    if not src or not dst or dst_size < 2:
        if dst and dst_size > 0:
            dst[0] = 0
        return

    while src[i] != 0 and j < dst_size - 6:  # Reserve space for longest escape sequence
        c = src[i]
        if c == 38:  # '&'
            dst[j] = 38; dst[j+1] = 97; dst[j+2] = 109; dst[j+3] = 112; dst[j+4] = 59  # &amp;
            j += 5
        elif c == 60:  # '<'
            dst[j] = 38; dst[j+1] = 108; dst[j+2] = 116; dst[j+3] = 59  # &lt;
            j += 4
        elif c == 62:  # '>'
            dst[j] = 38; dst[j+1] = 103; dst[j+2] = 116; dst[j+3] = 59  # &gt;
            j += 4
        elif c == 34:  # '"'
            dst[j] = 38; dst[j+1] = 113; dst[j+2] = 117; dst[j+3] = 111; dst[j+4] = 116; dst[j+5] = 59  # &quot;
            j += 6
        else:
            dst[j] = c
            j += 1
        i += 1

    dst[j] = 0


cdef int export_graph_graphml(
    WeightedGraph* graph,
    MemoryPool* pool,
    sam_hdr_t* bam_header,
    ReferenceMapping* mapping,
    ReferencePattern* pattern_data,
    ReferenceStats* ref_stats,
    const char* output_path,
    bint verbose
) noexcept nogil:
    """Export graph to GraphML format with comprehensive node and edge attributes.

    Uses igraph object for edge structure combined with custom Leiden clustering
    and reference statistics. Outputs XML format compatible with Cytoscape,
    igraph, Gephi, NetworkX, and yEd.

    Parameters
    ----------
    graph : WeightedGraph*
        Graph with igraph handle
    pool : MemoryPool*
        Memory pool for validation
    bam_header : sam_hdr_t*
        BAM header for reference names/lengths
    mapping : ReferenceMapping*
        Reference ID mapping
    pattern_data : ReferencePattern*
        Leiden clustering results
    ref_stats : ReferenceStats*
        Reference statistics
    output_path : const char*
        GraphML output file path
    verbose : bint
        Enable verbose output

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef FILE* f = NULL
    cdef igraph_t* ig_graph = NULL
    cdef uint32_t node_idx
    cdef const char* ref_name_c
    cdef int32_t tid32, original_tid
    cdef int64_t ref_len
    cdef char escaped_name[512]
    cdef uint32_t component_id, community_id
    cdef float community_cc, individual_cc, cc_threshold
    cdef char keep_flag
    cdef uint32_t total_reads, unique_reads, repeat_reads, shared_reads
    cdef uint64_t total_alignments
    cdef double multimap_pct, edge_weight
    cdef igraph_integer_t num_edges, edge_idx, from_node, to_node
    cdef igraph_vector_t* weights_vec = NULL

    if not graph or not output_path:
        return -1

    if verbose:
        bf_nogil_logf_notime(b"GRAPH", "\n=== EXPORTING GRAPH TO GRAPHML ===\n")
        bf_nogil_logf_notime(b"GRAPH", "Output: %s\n", output_path)
        bf_nogil_logf_notime(b"GRAPH", "Nodes: %u\n", graph.num_nodes)

    if not graph.igraph_handle:
        if verbose:
            bf_nogil_logf_notime(b"GRAPH", "WARNING: No igraph object available for export (igraph_handle is NULL)\n")
            bf_nogil_logf_notime(b"GRAPH", "         Graph export requires --clustering to be enabled\n")
        return -1

    ig_graph = <igraph_t*>graph.igraph_handle
    num_edges = igraph_ecount(ig_graph)

    if graph.weights_handle:
        weights_vec = <igraph_vector_t*>graph.weights_handle

    if verbose:
        bf_nogil_logf_notime(b"GRAPH", "Edges: %ld\n", <long>num_edges)

    f = fopen(output_path, b"w")
    if not f:
        bf_nogil_logf_notime(b"GRAPH", "ERROR: Failed to open %s for writing\n", output_path)
        return -1

    fprintf(f, b"<?xml version=\"1.0\" encoding=\"UTF-8\"?>\n")
    fprintf(f, b"<graphml xmlns=\"http://graphml.graphdrawing.org/xmlns\"\n")
    fprintf(f, b"         xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\"\n")
    fprintf(f, b"         xsi:schemaLocation=\"http://graphml.graphdrawing.org/xmlns\n")
    fprintf(f, b"         http://graphml.graphdrawing.org/xmlns/1.0/graphml.xsd\">\n\n")

    fprintf(f, b"  <!-- Node attributes -->\n")
    fprintf(f, b"  <key id=\"label\" for=\"node\" attr.name=\"label\" attr.type=\"string\"/>\n")
    fprintf(f, b"  <key id=\"length\" for=\"node\" attr.name=\"length\" attr.type=\"long\"/>\n")
    fprintf(f, b"  <key id=\"component_id\" for=\"node\" attr.name=\"component_id\" attr.type=\"string\"/>\n")
    fprintf(f, b"  <key id=\"community\" for=\"node\" attr.name=\"community\" attr.type=\"string\"/>\n")
    fprintf(f, b"  <key id=\"community_cc\" for=\"node\" attr.name=\"community_cc\" attr.type=\"double\"/>\n")
    fprintf(f, b"  <key id=\"individual_cc\" for=\"node\" attr.name=\"individual_cc\" attr.type=\"double\"/>\n")
    fprintf(f, b"  <key id=\"cc_threshold\" for=\"node\" attr.name=\"cc_threshold\" attr.type=\"double\"/>\n")
    fprintf(f, b"  <key id=\"keep_status\" for=\"node\" attr.name=\"keep_status\" attr.type=\"string\"/>\n")
    fprintf(f, b"  <key id=\"total_reads\" for=\"node\" attr.name=\"total_reads\" attr.type=\"int\"/>\n")
    fprintf(f, b"  <key id=\"unique_reads\" for=\"node\" attr.name=\"unique_reads\" attr.type=\"int\"/>\n")
    fprintf(f, b"  <key id=\"repeat_reads\" for=\"node\" attr.name=\"repeat_reads\" attr.type=\"int\"/>\n")
    fprintf(f, b"  <key id=\"shared_reads\" for=\"node\" attr.name=\"shared_reads\" attr.type=\"int\"/>\n")
    fprintf(f, b"  <key id=\"total_alignments\" for=\"node\" attr.name=\"total_alignments\" attr.type=\"long\"/>\n")
    fprintf(f, b"  <key id=\"multimap_pct\" for=\"node\" attr.name=\"multimap_pct\" attr.type=\"double\"/>\n\n")

    fprintf(f, b"  <!-- Edge attributes -->\n")
    fprintf(f, b"  <key id=\"weight\" for=\"edge\" attr.name=\"weight\" attr.type=\"double\"/>\n\n")

    fprintf(f, b"  <graph id=\"G\" edgedefault=\"undirected\">\n\n")

    fprintf(f, b"    <!-- Nodes (references) -->\n")

    cdef uint32_t exported_nodes = 0

    for node_idx in range(graph.num_nodes):
        if not pattern_data or node_idx >= pool.reference_count:
            continue

        node_degree = pattern_data[node_idx].node_degree

        if node_degree == 0:
            continue

        community_id = pattern_data[node_idx].leiden_community_id

        ref_name_c = b"unknown"
        ref_len = 0

        if mapping and mapping.new_to_old_tid and node_idx < mapping.n_retained_refs:
            tid32 = mapping.new_to_old_tid[node_idx]
            if tid32 >= 0 and bam_header:
                ref_name_c = sam_hdr_tid2name(bam_header, tid32)
                if ref_name_c:
                    ref_len = sam_hdr_tid2len(bam_header, tid32)
                else:
                    ref_name_c = b"unknown"

        if ref_name_c:
            xml_escape(ref_name_c, escaped_name, 512)
        else:
            escaped_name[0] = b'u'
            escaped_name[1] = b'n'
            escaped_name[2] = b'k'
            escaped_name[3] = b'n'
            escaped_name[4] = b'o'
            escaped_name[5] = b'w'
            escaped_name[6] = b'n'
            escaped_name[7] = 0

        component_id = pattern_data[node_idx].component_id
        community_id = pattern_data[node_idx].leiden_community_id
        community_cc = pattern_data[node_idx].leiden_community_cc
        individual_cc = pattern_data[node_idx].leiden_individual_cc
        cc_threshold = pattern_data[node_idx].leiden_cc_threshold
        keep_flag = pattern_data[node_idx].leiden_keep_flag

        total_reads = 0
        unique_reads = 0
        repeat_reads = 0
        shared_reads = 0
        total_alignments = 0
        multimap_pct = 0.0

        if ref_stats and node_idx < pool.reference_count:
            total_reads = ref_stats[node_idx].total_reads
            unique_reads = ref_stats[node_idx].unique_reads
            repeat_reads = ref_stats[node_idx].repeat_reads
            shared_reads = ref_stats[node_idx].shared_reads
            total_alignments = ref_stats[node_idx].alignment_count

            if total_reads > 0:
                multimap_pct = 100.0 * <double>(repeat_reads + shared_reads) / <double>total_reads

        fprintf(f, b"    <node id=\"n%u\">\n", node_idx)
        fprintf(f, b"      <data key=\"label\">%s</data>\n", escaped_name)
        fprintf(f, b"      <data key=\"length\">%ld</data>\n", ref_len)
        fprintf(f, b"      <data key=\"component_id\">comp_%u</data>\n", component_id)
        if community_id == UINT32_MAX:
            fprintf(f, b"      <data key=\"community\">singleton</data>\n")
        else:
            fprintf(f, b"      <data key=\"community\">comm_%u</data>\n", community_id)
        fprintf(f, b"      <data key=\"community_cc\">%.6f</data>\n", community_cc)
        fprintf(f, b"      <data key=\"individual_cc\">%.6f</data>\n", individual_cc)
        fprintf(f, b"      <data key=\"cc_threshold\">%.6f</data>\n", cc_threshold)
        if keep_flag:
            fprintf(f, b"      <data key=\"keep_status\">kept</data>\n")
        else:
            fprintf(f, b"      <data key=\"keep_status\">removed</data>\n")
        fprintf(f, b"      <data key=\"total_reads\">%u</data>\n", total_reads)
        fprintf(f, b"      <data key=\"unique_reads\">%u</data>\n", unique_reads)
        fprintf(f, b"      <data key=\"repeat_reads\">%u</data>\n", repeat_reads)
        fprintf(f, b"      <data key=\"shared_reads\">%u</data>\n", shared_reads)
        fprintf(f, b"      <data key=\"total_alignments\">%lu</data>\n", total_alignments)
        fprintf(f, b"      <data key=\"multimap_pct\">%.2f</data>\n", multimap_pct)
        fprintf(f, b"    </node>\n")

        exported_nodes += 1

    fprintf(f, b"\n    <!-- Edges (shared reads between references) -->\n")
    for edge_idx in range(num_edges):
        igraph_edge(ig_graph, edge_idx, &from_node, &to_node)

        edge_weight = 1.0
        if weights_vec:
            edge_weight = VECTOR(weights_vec[0])[edge_idx]

        fprintf(f, b"    <edge id=\"e%ld\" source=\"n%ld\" target=\"n%ld\">\n",
                <long>edge_idx, <long>from_node, <long>to_node)
        fprintf(f, b"      <data key=\"weight\">%.1f</data>\n", edge_weight)
        fprintf(f, b"    </edge>\n")

    fprintf(f, b"  </graph>\n")
    fprintf(f, b"</graphml>\n")

    fclose(f)

    if verbose:
        bf_nogil_logf_notime(b"GRAPH", "Graph exported successfully to %s\n", output_path)
        bf_nogil_logf_notime(
            b"GRAPH",
            "  Exported %u nodes (with degree > 0) and %ld edges\n",
            exported_nodes,
            <long>num_edges,
        )
        bf_nogil_logf_notime(
            b"GRAPH",
            "  Skipped %u isolated nodes (degree = 0)\n",
            graph.num_nodes - exported_nodes,
        )

    return 0
