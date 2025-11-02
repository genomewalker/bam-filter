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
from libc.string cimport strlen, strcpy, memcpy
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t
from libc.stddef cimport size_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReferencePattern, WeightedGraph
from bam_filter.processor_graph_ops cimport GraphNode
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.taxonomy_db cimport TaxonomyDB, TaxNode
from bam_filter.processor_igraph cimport (
    igraph_t, igraph_error_t, IGRAPH_SUCCESS,
    igraph_integer_t, igraph_vector_t, igraph_vector_int_t, igraph_real_t,
    igraph_vs_t, igraph_neimode_t, IGRAPH_ALL,
    igraph_ecount, igraph_edge
)

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

from bam_filter.processor cimport sam_hdr_t, sam_hdr_tid2name, sam_hdr_tid2len

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


# Lineage data structure for taxonomy information
cdef struct LineageData:
    const char* name_superkingdom
    const char* name_clade
    const char* name_kingdom
    const char* name_phylum
    const char* name_class
    const char* name_order
    const char* name_family
    const char* name_genus
    const char* name_species
    const char* name_subspecies
    char lineage_string[1024]


cdef inline int c_strcmp(const char* s1, const char* s2) nogil:
    """Fast C string comparison (returns 0 if equal)."""
    cdef int i = 0
    while s1[i] != 0 and s2[i] != 0:
        if s1[i] != s2[i]:
            return 1
        i += 1
    if s1[i] == s2[i]:
        return 0
    return 1


cdef void extract_lineage_data_graphml(TaxonomyDB* taxonomy_db, int32_t taxid, LineageData* out) nogil:
    """Extract lineage data for GraphML export (same logic as TSV version)."""
    # Initialize all to NULL/empty
    out.name_superkingdom = NULL
    out.name_clade = NULL
    out.name_kingdom = NULL
    out.name_phylum = NULL
    out.name_class = NULL
    out.name_order = NULL
    out.name_family = NULL
    out.name_genus = NULL
    out.name_species = NULL
    out.name_subspecies = NULL
    out.lineage_string[0] = 0

    if taxonomy_db == NULL or taxid <= 0 or taxid > taxonomy_db.max_taxid:
        return

    cdef int32_t current_idx = taxonomy_db.taxid_to_idx[taxid]
    if current_idx < 0:
        return

    cdef int32_t parent_taxid
    cdef int32_t rank_id
    cdef const char* node_name
    cdef const char* rank_name
    cdef TaxNode* node

    while current_idx >= 0:
        node = &taxonomy_db.nodes[current_idx]
        rank_id = node.rank_id
        node_name = taxonomy_db.names_buffer + node.name_offset

        if rank_id >= 0 and rank_id < taxonomy_db.n_ranks and taxonomy_db.rank_names != NULL:
            rank_name = taxonomy_db.rank_names[rank_id]
        else:
            rank_name = NULL

        if rank_name != NULL:
            if c_strcmp(rank_name, "superkingdom") == 0 or c_strcmp(rank_name, "domain") == 0:
                if out.name_superkingdom == NULL:
                    out.name_superkingdom = node_name
            elif c_strcmp(rank_name, "clade") == 0 or c_strcmp(rank_name, "lineage") == 0:
                if out.name_clade == NULL:
                    out.name_clade = node_name
            elif c_strcmp(rank_name, "kingdom") == 0:
                if out.name_kingdom == NULL:
                    out.name_kingdom = node_name
            elif c_strcmp(rank_name, "phylum") == 0:
                if out.name_phylum == NULL:
                    out.name_phylum = node_name
            elif c_strcmp(rank_name, "class") == 0:
                if out.name_class == NULL:
                    out.name_class = node_name
            elif c_strcmp(rank_name, "order") == 0:
                if out.name_order == NULL:
                    out.name_order = node_name
            elif c_strcmp(rank_name, "family") == 0:
                if out.name_family == NULL:
                    out.name_family = node_name
            elif c_strcmp(rank_name, "genus") == 0:
                if out.name_genus == NULL:
                    out.name_genus = node_name
            elif c_strcmp(rank_name, "species") == 0:
                if out.name_species == NULL:
                    out.name_species = node_name
            elif c_strcmp(rank_name, "subspecies") == 0 or c_strcmp(rank_name, "strain") == 0:
                if out.name_subspecies == NULL:
                    out.name_subspecies = node_name

        parent_taxid = node.parent_taxid
        if parent_taxid == node.taxid:
            break
        if parent_taxid < 0 or parent_taxid > taxonomy_db.max_taxid:
            break
        current_idx = taxonomy_db.taxid_to_idx[parent_taxid]

    # Build lineage string
    cdef char* ptr = out.lineage_string
    cdef size_t len

    if out.name_superkingdom != NULL:
        len = strlen(out.name_superkingdom)
        memcpy(ptr, out.name_superkingdom, len)
        ptr += len
    if out.name_kingdom != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_kingdom)
        memcpy(ptr, out.name_kingdom, len)
        ptr += len
    if out.name_phylum != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_phylum)
        memcpy(ptr, out.name_phylum, len)
        ptr += len
    if out.name_class != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_class)
        memcpy(ptr, out.name_class, len)
        ptr += len
    if out.name_order != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_order)
        memcpy(ptr, out.name_order, len)
        ptr += len
    if out.name_family != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_family)
        memcpy(ptr, out.name_family, len)
        ptr += len
    if out.name_genus != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_genus)
        memcpy(ptr, out.name_genus, len)
        ptr += len
    if out.name_species != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_species)
        memcpy(ptr, out.name_species, len)
        ptr += len
    if out.name_subspecies != NULL:
        if ptr != out.lineage_string:
            ptr[0] = 59
            ptr += 1
        len = strlen(out.name_subspecies)
        memcpy(ptr, out.name_subspecies, len)
        ptr += len
    ptr[0] = 0


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
    bint verbose,
    bint export_only_used_stats,
    TaxonomyDB* taxonomy_db
) noexcept nogil:
    """Export graph to GraphML format with comprehensive node and edge attributes.

    Uses igraph object for edge structure combined with custom Community clustering
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
        Community clustering results
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
    cdef float anomaly_score
    cdef float betweenness
    cdef uint32_t node_degree
    cdef uint32_t connected_neighbors
    cdef double avg_co_mappings

    # Tier result variables
    cdef char structural_role
    cdef uint32_t num_neighbor_communities
    cdef char community_coherent
    cdef char filter_decision
    cdef uint8_t taxonomy_outlier_score

    # Score variables (alignment and PMD)
    cdef double score_mean, score_std, score_min, score_max
    cdef uint32_t pmd_available
    cdef double pmd_mean, pmd_std, pmd_min, pmd_max
    cdef uint64_t pmd_nonzero_pct

    # String mappings for tier results
    cdef const char* role_names[4]
    role_names[0] = b"PERIPHERAL"
    role_names[1] = b"CORE"
    role_names[2] = b"HUB"
    role_names[3] = b"BRIDGE"

    cdef const char* decision_names[3]
    decision_names[0] = b"KEEP"
    decision_names[1] = b"REMOVE"
    decision_names[2] = b"REVIEW"

    cdef uint32_t max_co_map
    cdef double neighbor_mm_rate
    cdef uint32_t total_reads, unique_reads, repeat_reads, shared_reads
    cdef uint64_t total_alignments
    cdef double multimap_pct
    cdef double edge_weight
    cdef igraph_integer_t num_edges, edge_idx, from_node, to_node
    cdef igraph_vector_t* weights_vec = NULL
    # Taxonomy variables
    cdef int32_t taxid
    cdef LineageData lineage_data
    cdef char escaped_tax[512]

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
    fprintf(f, b"  <key id=\"individual_cc\" for=\"node\" attr.name=\"individual_cc\" attr.type=\"double\"/>\n\n")

    # Taxonomy attributes (conditional on taxonomy_db)
    if taxonomy_db != NULL:
        fprintf(f, b"  <key id=\"taxid\" for=\"node\" attr.name=\"taxid\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"taxid_rank_id\" for=\"node\" attr.name=\"taxid_rank_id\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"taxid_depth\" for=\"node\" attr.name=\"taxid_depth\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"taxonomy_flag\" for=\"node\" attr.name=\"taxonomy_flag\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_neighbors_total\" for=\"node\" attr.name=\"tax_neighbors_total\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_domain\" for=\"node\" attr.name=\"tax_mismatch_domain\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_kingdom\" for=\"node\" attr.name=\"tax_mismatch_kingdom\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_phylum\" for=\"node\" attr.name=\"tax_mismatch_phylum\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_class\" for=\"node\" attr.name=\"tax_mismatch_class\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_order\" for=\"node\" attr.name=\"tax_mismatch_order\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_mismatch_family\" for=\"node\" attr.name=\"tax_mismatch_family\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"tax_match_genus_below\" for=\"node\" attr.name=\"tax_match_genus_below\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"superkingdom\" for=\"node\" attr.name=\"superkingdom\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"clade\" for=\"node\" attr.name=\"clade\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"kingdom\" for=\"node\" attr.name=\"kingdom\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"phylum\" for=\"node\" attr.name=\"phylum\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"class\" for=\"node\" attr.name=\"class\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"order\" for=\"node\" attr.name=\"order\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"family\" for=\"node\" attr.name=\"family\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"genus\" for=\"node\" attr.name=\"genus\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"species\" for=\"node\" attr.name=\"species\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"subspecies\" for=\"node\" attr.name=\"subspecies\" attr.type=\"string\"/>\n\n")

    if not export_only_used_stats:
        fprintf(f, b"  <key id=\"betweenness_centrality\" for=\"node\" attr.name=\"betweenness_centrality\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"structural_role\" for=\"node\" attr.name=\"structural_role\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"num_neighbor_communities\" for=\"node\" attr.name=\"num_neighbor_communities\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"community_coherent\" for=\"node\" attr.name=\"community_coherent\" attr.type=\"boolean\"/>\n")
        fprintf(f, b"  <key id=\"filter_decision\" for=\"node\" attr.name=\"filter_decision\" attr.type=\"string\"/>\n")
        fprintf(f, b"  <key id=\"taxonomy_outlier_score\" for=\"node\" attr.name=\"taxonomy_outlier_score\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"node_degree\" for=\"node\" attr.name=\"node_degree\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"connected_neighbors\" for=\"node\" attr.name=\"connected_neighbors\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"avg_co_mappings_per_read\" for=\"node\" attr.name=\"avg_co_mappings_per_read\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"max_co_mappings_observed\" for=\"node\" attr.name=\"max_co_mappings_observed\" attr.type=\"long\"/>\n")
        fprintf(f, b"  <key id=\"neighbor_multimap_rate\" for=\"node\" attr.name=\"neighbor_multimap_rate\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"total_reads\" for=\"node\" attr.name=\"total_reads\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"unique_reads\" for=\"node\" attr.name=\"unique_reads\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"repeat_reads\" for=\"node\" attr.name=\"repeat_reads\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"shared_reads\" for=\"node\" attr.name=\"shared_reads\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"total_alignments\" for=\"node\" attr.name=\"total_alignments\" attr.type=\"long\"/>\n")
        fprintf(f, b"  <key id=\"multimap_pct\" for=\"node\" attr.name=\"multimap_pct\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"alignment_score_mean\" for=\"node\" attr.name=\"alignment_score_mean\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"alignment_score_std\" for=\"node\" attr.name=\"alignment_score_std\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"alignment_score_min\" for=\"node\" attr.name=\"alignment_score_min\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"alignment_score_max\" for=\"node\" attr.name=\"alignment_score_max\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"pmd_available\" for=\"node\" attr.name=\"pmd_available\" attr.type=\"int\"/>\n")
        fprintf(f, b"  <key id=\"pmd_score_mean\" for=\"node\" attr.name=\"pmd_score_mean\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"pmd_score_std\" for=\"node\" attr.name=\"pmd_score_std\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"pmd_score_min\" for=\"node\" attr.name=\"pmd_score_min\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"pmd_score_max\" for=\"node\" attr.name=\"pmd_score_max\" attr.type=\"double\"/>\n")
        fprintf(f, b"  <key id=\"pmd_nonzero_percentage\" for=\"node\" attr.name=\"pmd_nonzero_percentage\" attr.type=\"long\"/>\n\n")

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

        community_id = pattern_data[node_idx].community_id

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
        community_id = pattern_data[node_idx].community_id
        community_cc = pattern_data[node_idx].community_cc
        individual_cc = pattern_data[node_idx].community_individual_cc
        betweenness = pattern_data[node_idx].betweenness_centrality
        node_degree = pattern_data[node_idx].node_degree
        connected_neighbors = pattern_data[node_idx].graph.connection_count

        # Load tier results
        structural_role = pattern_data[node_idx].structural_role
        num_neighbor_communities = pattern_data[node_idx].num_neighbor_communities
        community_coherent = pattern_data[node_idx].community_coherent
        filter_decision = pattern_data[node_idx].filter_decision
        taxonomy_outlier_score = pattern_data[node_idx].taxonomy_outlier_score
        avg_co_mappings = pattern_data[node_idx].graph.avg_comappings_per_read
        max_co_map = pattern_data[node_idx].graph.max_comappings
        neighbor_mm_rate = pattern_data[node_idx].graph.neighbor_avg_multimap

        total_reads = 0
        unique_reads = 0
        repeat_reads = 0
        shared_reads = 0
        total_alignments = 0
        multimap_pct = 0.0
        score_mean = 0.0
        score_std = 0.0
        score_min = 0.0
        score_max = 0.0
        pmd_available = 0
        pmd_mean = 0.0
        pmd_std = 0.0
        pmd_min = 0.0
        pmd_max = 0.0
        pmd_nonzero_pct = 0

        if ref_stats != NULL and node_idx < pool.reference_count:
            total_reads = ref_stats[node_idx].total_reads
            unique_reads = ref_stats[node_idx].unique_reads
            repeat_reads = ref_stats[node_idx].repeat_reads
            shared_reads = ref_stats[node_idx].shared_reads
            total_alignments = ref_stats[node_idx].alignment_count
            score_mean = ref_stats[node_idx].score_mean
            score_std = ref_stats[node_idx].score_std
            score_min = ref_stats[node_idx].score_min
            score_max = ref_stats[node_idx].score_max
            pmd_available = 1 if ref_stats[node_idx].pmd_available else 0
            pmd_mean = ref_stats[node_idx].pmd_mean
            pmd_std = ref_stats[node_idx].pmd_std
            pmd_min = ref_stats[node_idx].pmd_min
            pmd_max = ref_stats[node_idx].pmd_max

            if total_reads > 0:
                multimap_pct = 100.0 * <double>(repeat_reads + shared_reads) / <double>total_reads
            if total_alignments > 0:
                pmd_nonzero_pct = <uint64_t>(100.0 * <double>ref_stats[node_idx].pmd_nonzero_count / <double>total_alignments)

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

        # Extract and write taxonomy information
        if taxonomy_db != NULL and pattern_data != NULL:
            taxid = pattern_data[node_idx].taxid

            if taxid > 0:
                extract_lineage_data_graphml(taxonomy_db, taxid, &lineage_data)

                fprintf(f, b"      <data key=\"taxid\">%d</data>\n", taxid)
                fprintf(f, b"      <data key=\"taxid_rank_id\">%d</data>\n", pattern_data[node_idx].taxid_rank_id)
                fprintf(f, b"      <data key=\"taxid_depth\">%d</data>\n", pattern_data[node_idx].taxid_depth)
                fprintf(f, b"      <data key=\"taxonomy_flag\">%d</data>\n", <int>pattern_data[node_idx].taxonomy_flag)
                fprintf(f, b"      <data key=\"tax_neighbors_total\">%u</data>\n", pattern_data[node_idx].tax_neighbors_total)
                fprintf(f, b"      <data key=\"tax_mismatch_domain\">%u</data>\n", pattern_data[node_idx].tax_mismatch_domain)
                fprintf(f, b"      <data key=\"tax_mismatch_kingdom\">%u</data>\n", pattern_data[node_idx].tax_mismatch_kingdom)
                fprintf(f, b"      <data key=\"tax_mismatch_phylum\">%u</data>\n", pattern_data[node_idx].tax_mismatch_phylum)
                fprintf(f, b"      <data key=\"tax_mismatch_class\">%u</data>\n", pattern_data[node_idx].tax_mismatch_class)
                fprintf(f, b"      <data key=\"tax_mismatch_order\">%u</data>\n", pattern_data[node_idx].tax_mismatch_order)
                fprintf(f, b"      <data key=\"tax_mismatch_family\">%u</data>\n", pattern_data[node_idx].tax_mismatch_family)
                fprintf(f, b"      <data key=\"tax_match_genus_below\">%u</data>\n", pattern_data[node_idx].tax_match_genus_below)

                if lineage_data.name_superkingdom != NULL:
                    xml_escape(lineage_data.name_superkingdom, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"superkingdom\">%s</data>\n", escaped_tax)
                if lineage_data.name_clade != NULL:
                    xml_escape(lineage_data.name_clade, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"clade\">%s</data>\n", escaped_tax)
                if lineage_data.name_kingdom != NULL:
                    xml_escape(lineage_data.name_kingdom, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"kingdom\">%s</data>\n", escaped_tax)
                if lineage_data.name_phylum != NULL:
                    xml_escape(lineage_data.name_phylum, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"phylum\">%s</data>\n", escaped_tax)
                if lineage_data.name_class != NULL:
                    xml_escape(lineage_data.name_class, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"class\">%s</data>\n", escaped_tax)
                if lineage_data.name_order != NULL:
                    xml_escape(lineage_data.name_order, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"order\">%s</data>\n", escaped_tax)
                if lineage_data.name_family != NULL:
                    xml_escape(lineage_data.name_family, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"family\">%s</data>\n", escaped_tax)
                if lineage_data.name_genus != NULL:
                    xml_escape(lineage_data.name_genus, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"genus\">%s</data>\n", escaped_tax)
                if lineage_data.name_species != NULL:
                    xml_escape(lineage_data.name_species, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"species\">%s</data>\n", escaped_tax)
                if lineage_data.name_subspecies != NULL:
                    xml_escape(lineage_data.name_subspecies, escaped_tax, 512)
                    fprintf(f, b"      <data key=\"subspecies\">%s</data>\n", escaped_tax)

        if not export_only_used_stats:
            fprintf(f, b"      <data key=\"betweenness_centrality\">%.6f</data>\n", betweenness)
            fprintf(f, b"      <data key=\"structural_role\">%s</data>\n", role_names[<int>structural_role])
            fprintf(f, b"      <data key=\"num_neighbor_communities\">%u</data>\n", num_neighbor_communities)
            fprintf(f, b"      <data key=\"community_coherent\">%s</data>\n", b"true" if community_coherent else b"false")
            fprintf(f, b"      <data key=\"filter_decision\">%s</data>\n", decision_names[<int>filter_decision])
            fprintf(f, b"      <data key=\"taxonomy_outlier_score\">%u</data>\n", taxonomy_outlier_score)
            fprintf(f, b"      <data key=\"node_degree\">%u</data>\n", node_degree)
            fprintf(f, b"      <data key=\"connected_neighbors\">%u</data>\n", connected_neighbors)
            fprintf(f, b"      <data key=\"avg_co_mappings_per_read\">%.6f</data>\n", avg_co_mappings)
            fprintf(f, b"      <data key=\"max_co_mappings_observed\">%lu</data>\n", max_co_map)
            fprintf(f, b"      <data key=\"neighbor_multimap_rate\">%.6f</data>\n", neighbor_mm_rate)
            fprintf(f, b"      <data key=\"total_reads\">%u</data>\n", total_reads)
            fprintf(f, b"      <data key=\"unique_reads\">%u</data>\n", unique_reads)
            fprintf(f, b"      <data key=\"repeat_reads\">%u</data>\n", repeat_reads)
            fprintf(f, b"      <data key=\"shared_reads\">%u</data>\n", shared_reads)
            fprintf(f, b"      <data key=\"total_alignments\">%lu</data>\n", total_alignments)
            fprintf(f, b"      <data key=\"multimap_pct\">%.2f</data>\n", multimap_pct)
            fprintf(f, b"      <data key=\"alignment_score_mean\">%.6f</data>\n", score_mean)
            fprintf(f, b"      <data key=\"alignment_score_std\">%.6f</data>\n", score_std)
            fprintf(f, b"      <data key=\"alignment_score_min\">%.6f</data>\n", score_min)
            fprintf(f, b"      <data key=\"alignment_score_max\">%.6f</data>\n", score_max)
            fprintf(f, b"      <data key=\"pmd_available\">%u</data>\n", pmd_available)
            fprintf(f, b"      <data key=\"pmd_score_mean\">%.6f</data>\n", pmd_mean)
            fprintf(f, b"      <data key=\"pmd_score_std\">%.6f</data>\n", pmd_std)
            fprintf(f, b"      <data key=\"pmd_score_min\">%.6f</data>\n", pmd_min)
            fprintf(f, b"      <data key=\"pmd_score_max\">%.6f</data>\n", pmd_max)
            fprintf(f, b"      <data key=\"pmd_nonzero_percentage\">%lu</data>\n", pmd_nonzero_pct)
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
