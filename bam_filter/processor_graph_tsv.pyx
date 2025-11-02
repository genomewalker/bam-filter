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

from cython.parallel cimport prange
# -*- coding: utf-8 -*-
"""
TSV writer extension for graph analysis output.

This module exposes a C-level function ``write_graph_tsv_c`` that other
Cython modules can cimport and call with low-level pointers. The function
writes a per-reference TSV (optionally gzipped) summarizing graph and
clustering metrics.
"""

from libc.stdlib cimport malloc, free, qsort
from libc.string cimport strlen, strcpy, memcpy
from libc.stdio cimport FILE, fopen, fclose, fprintf
from libc.stdint cimport uint32_t, uint64_t, int32_t, int64_t, uint8_t, uint16_t
from libc.stddef cimport size_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats, MultPairExtended
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.taxonomy_db cimport TaxonomyDB, TaxNode

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

cdef extern from "zlib.h":
    ctypedef void* gzFile
    gzFile gzopen(const char* path, const char* mode) nogil
    int gzclose(gzFile file) nogil
    int gzprintf(gzFile file, const char* format, ...) nogil

cdef extern from "stdio.h":
    int snprintf(char* s, size_t n, const char* format, ...) nogil


cdef inline bint c_str_endswith_gz(const char* s) nogil:
    """Check if C string ends with '.gz' extension."""
    if s == NULL:
        return 0
    cdef int L = <int>strlen(s)
    if L < 3:
        return 0
    if s[L-3] == 46 and s[L-2] == 103 and s[L-1] == 122:
        return 1
    return 0


cdef uint32_t UINT32_MAX = 0xFFFFFFFF

# Structure to hold precomputed lineage data for parallel processing
cdef struct LineageData:
    const char* name_superkingdom
    const char* name_clade           # Intermediate rank between superkingdom and kingdom (e.g., SAR)
    const char* name_kingdom
    const char* name_phylum
    const char* name_class
    const char* name_order
    const char* name_family
    const char* name_genus
    const char* name_species
    const char* name_subspecies

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


cdef void extract_lineage_data(TaxonomyDB* taxonomy_db, int32_t taxid, LineageData* out) nogil:
    """Extract lineage data for a single taxid in one tree traversal.

    Robust version that works with any taxonomy (NCBI, GTDB, custom):
    - Uses rank NAME comparison instead of hard-coded rank IDs
    - Handles variations in rank naming (superkingdom/domain, etc.)
    - Works even if rank IDs change between taxonomy versions
    """
    # Initialize all to NULL
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

    if taxonomy_db == NULL or taxid <= 0 or taxid > taxonomy_db.max_taxid:
        return

    # Single upward traversal collecting all rank names
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

        # Get rank name string for this node
        if rank_id >= 0 and rank_id < taxonomy_db.n_ranks and taxonomy_db.rank_names != NULL:
            rank_name = taxonomy_db.rank_names[rank_id]
        else:
            rank_name = NULL

        # Store name at appropriate rank using RANK NAME comparison
        # This is robust to taxonomy changes and custom taxonomies
        if rank_name != NULL:
            # Superkingdom / Domain (check multiple variations)
            if c_strcmp(rank_name, "superkingdom") == 0 or c_strcmp(rank_name, "domain") == 0:
                if out.name_superkingdom == NULL:  # Only set if not already set
                    out.name_superkingdom = node_name
            # Clade/Lineage (intermediate rank between superkingdom and kingdom, e.g., SAR or custom taxonomies)
            elif c_strcmp(rank_name, "clade") == 0 or c_strcmp(rank_name, "lineage") == 0:
                if out.name_clade == NULL:
                    out.name_clade = node_name
            # Kingdom
            elif c_strcmp(rank_name, "kingdom") == 0:
                if out.name_kingdom == NULL:
                    out.name_kingdom = node_name
            # Phylum
            elif c_strcmp(rank_name, "phylum") == 0:
                if out.name_phylum == NULL:
                    out.name_phylum = node_name
            # Class
            elif c_strcmp(rank_name, "class") == 0:
                if out.name_class == NULL:
                    out.name_class = node_name
            # Order
            elif c_strcmp(rank_name, "order") == 0:
                if out.name_order == NULL:
                    out.name_order = node_name
            # Family
            elif c_strcmp(rank_name, "family") == 0:
                if out.name_family == NULL:
                    out.name_family = node_name
            # Genus
            elif c_strcmp(rank_name, "genus") == 0:
                if out.name_genus == NULL:
                    out.name_genus = node_name
            # Species
            elif c_strcmp(rank_name, "species") == 0:
                if out.name_species == NULL:
                    out.name_species = node_name
            # Subspecies
            elif c_strcmp(rank_name, "subspecies") == 0 or c_strcmp(rank_name, "strain") == 0:
                if out.name_subspecies == NULL:
                    out.name_subspecies = node_name

        # Move to parent
        parent_taxid = node.parent_taxid

        # Check if we reached root (parent == self)
        if parent_taxid == node.taxid:
            break

        if parent_taxid < 0 or parent_taxid > taxonomy_db.max_taxid:
            break

        current_idx = taxonomy_db.taxid_to_idx[parent_taxid]

cdef int _multpair_extended_cmp_local(const void* a, const void* b) noexcept nogil:
    """Comparator for qsort: sort references by neighbor_count desc, then read_count desc."""
    cdef MultPairExtended* A = <MultPairExtended*>a
    cdef MultPairExtended* B = <MultPairExtended*>b

    if A.neighbor_count != B.neighbor_count:
        return -1 if A.neighbor_count > B.neighbor_count else 1

    if A.read_count != B.read_count:
        return -1 if A.read_count > B.read_count else 1

    return 0


cdef int write_graph_tsv_c(MemoryPool* pool, sam_hdr_t* bam_header,
                          ReferenceMapping* mapping, ReferencePattern* pattern_data,
                          ReferenceStats* ref_stats,
                          uint32_t* total_reads, uint32_t* multimap_reads, uint64_t* alignments_per_ref,
                          uint32_t* exact_connection_counts, double* co_mapping_averages,
                          uint64_t* max_co_mappings, uint32_t* co_mapping_counts,
                          double* neighbor_multimap_avg, double* neighbor_connections_avg,
                          uint32_t* neighbor_counts, uint32_t array_size,
                          double dataset_median_connections, int32_t min_read_count,
                          bint include_clustering, int outlier_method, const char* tsv_path,
                          TaxonomyDB* taxonomy_db) noexcept nogil:
    """
    Write a TSV (optionally gzipped) summarizing graph analysis results per
    reference.

    This C-level function implements the TSV export previously located in
    :file:`bam_filter/processor_graph.pyx::write_graph_tsv` and was moved here
    to reduce the size of the main graph module. It is a low-level, ``nogil``
    function and must not perform Python API calls.

    Parameters
    ----------
    pool : MemoryPool*
        Pointer to the shared memory pool describing retained references and
        counts. Used for sizing and reference_count lookup.
    bam_header : sam_hdr_t*
        HTSlib header pointer used to map tid -> reference name/length. May
        be NULL, in which case reference name will be reported as ``unknown``.
    mapping : ReferenceMapping*
        Optional mapping object that translates new tids to original tids.
    pattern_data : ReferencePattern*
        Optional per-reference pattern/graph statistics produced by the
        clustering/graph pipeline. Used to populate community detection fields.
    ref_stats : ReferenceStats*
        Optional array of per-reference numeric statistics (read counts,
        score means, PMD aggregates). May be NULL when statistics are
        unavailable.
    total_reads, multimap_reads : uint32_t*
        Arrays of per-reference read counts. If provided, used to compute
        multimap percentages and per-reference read totals.
    alignments_per_ref : uint64_t*
        Array of total alignment counts per reference (optional).
    exact_connection_counts : uint32_t*
        If provided, exact neighbor counts for each reference (connection
        degrees) used for sorting and node-degree columns.
    co_mapping_averages : double*
        Per-reference average co-mappings per read (optional).
    max_co_mappings : uint64_t*
        Maximum co-mapping observed per reference (optional).
    co_mapping_counts : uint32_t*
        Optional counts used for distribution reporting.
    neighbor_multimap_avg, neighbor_connections_avg : double*
        Optional per-reference neighbor-level aggregates.
    neighbor_counts : uint32_t*
        Per-reference neighbor counts array (optional).
    array_size : uint32_t
        Length of the arrays passed above (safety bound for indexing).
    dataset_median_connections : double
        Median connections value for the dataset (used for contextual output).
    min_read_count : int32_t
        Minimum reads threshold used when deciding which references to emit.
    include_clustering : bint
        If true, Community/community columns will be filled when pattern_data is
        available.
    tsv_path : const char*
        Filesystem path for the TSV output. If it ends with ``.gz`` the file
        will be written gzipped using zlib's ``gzopen``/``gzprintf`` APIs.

    Returns
    -------
    int
        0 on success, -1 on any error (file open failure or allocation
        failure). This function performs internal allocations with ``malloc``
        and frees them before returning.

    Notes
    -----
    - This function is ``nogil`` and ``noexcept``. Do not introduce Python
      calls or raise exceptions inside it.
    - The caller is responsible for ensuring the arrays passed in are at
      least ``array_size`` long and that ``pool`` describes the correct
      reference_count.
    """

    cdef FILE* f = NULL
    cdef gzFile gzfp = NULL
    cdef uint32_t ref_idx, i
    cdef uint32_t node_degree
    cdef const char* ref_name_c
    cdef int32_t tid32
    cdef int64_t ref_len
    cdef double multimap_pct
    cdef uint32_t treads
    cdef uint32_t repeat_reads_val
    cdef uint32_t shared_reads_val
    cdef uint32_t connected_neighbors
    cdef double avg_co_mappings
    cdef uint64_t max_co_map
    cdef double neighbor_mm_rate
    cdef uint32_t unique_reads
    cdef uint32_t repeat_reads
    cdef uint32_t shared_reads
    cdef uint64_t align_count
    cdef double score_mean
    cdef double score_std
    cdef double score_min
    cdef double score_max
    cdef uint32_t pmd_available
    cdef double pmd_mean
    cdef double pmd_std
    cdef double pmd_min
    cdef double pmd_max
    cdef uint64_t pmd_nonzero_pct
    cdef uint32_t component_id, community_id
    cdef float community_cc
    cdef float individual_cc
    cdef char keep_flag
    cdef float betweenness
    cdef int32_t taxid
    cdef int32_t taxid_rank_id
    cdef int32_t taxid_depth
    cdef char taxonomy_flag

    # Per-rank taxonomy mismatch counts
    cdef uint16_t tax_neighbors_total
    cdef uint16_t tax_mismatch_domain
    cdef uint16_t tax_mismatch_kingdom
    cdef uint16_t tax_mismatch_phylum
    cdef uint16_t tax_mismatch_class
    cdef uint16_t tax_mismatch_order
    cdef uint16_t tax_mismatch_family
    cdef uint16_t tax_match_genus_below

    # Tier result variables
    cdef char structural_role
    cdef uint32_t num_neighbor_communities
    cdef char community_coherent
    cdef char filter_decision
    cdef uint8_t taxonomy_outlier_score

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

    cdef uint32_t ref_idx_par  # For parallel lineage extraction
    cdef LineageData temp_lineage  # Temporary for on-demand lineage extraction
    cdef uint32_t cache_idx  # For lineage cache lookup
    cdef bint found_in_cache  # Cache hit flag
    # Lineage names for 10 standard ranks + clade
    cdef const char* name_superkingdom
    cdef const char* name_clade
    cdef const char* name_kingdom
    cdef const char* name_phylum
    cdef const char* name_class
    cdef const char* name_order
    cdef const char* name_family
    cdef const char* name_genus
    cdef const char* name_species
    cdef const char* name_subspecies
    cdef const char* node_name
    cdef MultPairExtended* sorted_refs = NULL
    cdef uint32_t valid_count = 0

    cdef uint32_t UINT32_MAX = <uint32_t>0xFFFFFFFF

    if not tsv_path:
        return -1

    if c_str_endswith_gz(tsv_path):
        gzfp = gzopen(tsv_path, "wb")
        if gzfp == NULL:
            return -1

        # Write header - conditional on taxonomy availability
        if taxonomy_db != NULL:
            gzprintf(gzfp, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tcommunity_id\tcommunity_cc\tcommunity_individual_cc\tbetweenness_centrality\tstructural_role\tnum_neighbor_communities\tcommunity_coherent\tfilter_decision\ttaxonomy_outlier_score\ttaxid\ttaxid_rank_id\ttaxid_depth\ttaxonomy_flag\ttax_neighbors_total\ttax_mismatch_domain\ttax_mismatch_kingdom\ttax_mismatch_phylum\ttax_mismatch_class\ttax_mismatch_order\ttax_mismatch_family\ttax_match_genus_below\tsuperkingdom\tclade\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\tmisannotation_flag\tmisannotation_confidence\tcross_domain_edges_before\tedges_after_removal\tcross_domain_fraction\n")
        else:
            gzprintf(gzfp, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tcommunity_id\tcommunity_cc\tcommunity_individual_cc\tbetweenness_centrality\tstructural_role\tnum_neighbor_communities\tcommunity_coherent\tfilter_decision\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\tmisannotation_flag\tmisannotation_confidence\tcross_domain_edges_before\tedges_after_removal\tcross_domain_fraction\n")
    else:
        f = fopen(tsv_path, b"w")
        if not f:
            return -1

        # Write header - conditional on taxonomy availability
        if taxonomy_db != NULL:
            fprintf(f, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tcommunity_id\tcommunity_cc\tcommunity_individual_cc\tbetweenness_centrality\tstructural_role\tnum_neighbor_communities\tcommunity_coherent\tfilter_decision\ttaxonomy_outlier_score\ttaxid\ttaxid_rank_id\ttaxid_depth\ttaxonomy_flag\ttax_neighbors_total\ttax_mismatch_domain\ttax_mismatch_kingdom\ttax_mismatch_phylum\ttax_mismatch_class\ttax_mismatch_order\ttax_mismatch_family\ttax_match_genus_below\tsuperkingdom\tclade\tkingdom\tphylum\tclass\torder\tfamily\tgenus\tspecies\tsubspecies\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\tmisannotation_flag\tmisannotation_confidence\tcross_domain_edges_before\tedges_after_removal\tcross_domain_fraction\n")
        else:
            fprintf(f, b"reference_name\treference_length_bp\tcomponent_id\tnode_degree\tcommunity_id\tcommunity_cc\tcommunity_individual_cc\tbetweenness_centrality\tstructural_role\tnum_neighbor_communities\tcommunity_coherent\tfilter_decision\ttotal_reads\tunique_reads\trepeat_reads\tshared_reads\ttotal_alignments\tmultimap_percentage\tconnected_neighbors\tavg_co_mappings_per_read\tmax_co_mappings_observed\tneighbor_multimap_rate\talignment_score_mean\talignment_score_std\talignment_score_min\talignment_score_max\tpmd_available\tpmd_score_mean\tpmd_score_std\tpmd_score_min\tpmd_score_max\tpmd_nonzero_percentage\tmisannotation_flag\tmisannotation_confidence\tcross_domain_edges_before\tedges_after_removal\tcross_domain_fraction\n")

    sorted_refs = <MultPairExtended*>malloc(pool.reference_count * sizeof(MultPairExtended))
    if not sorted_refs:
        if f: fclose(f)
        if gzfp != NULL: gzclose(gzfp)
        return -1

    # Cache for lineage data: map taxid -> LineageData
    # Use a simple linear search cache (small enough to be fast)
    cdef int32_t* cached_taxids = NULL
    cdef LineageData* cached_lineages = NULL
    cdef uint32_t cache_size = 0
    cdef uint32_t cache_capacity = 0

    if taxonomy_db != NULL:
        # Allocate initial cache capacity (will grow as needed)
        cache_capacity = 1024  # Start with 1024 unique taxids
        cached_taxids = <int32_t*>malloc(cache_capacity * sizeof(int32_t))
        cached_lineages = <LineageData*>malloc(cache_capacity * sizeof(LineageData))
        if not cached_taxids or not cached_lineages:
            if cached_taxids: free(cached_taxids)
            if cached_lineages: free(cached_lineages)
            cached_taxids = NULL
            cached_lineages = NULL

    for ref_idx in range(pool.reference_count):
        treads = 0
        if total_reads and ref_idx < pool.reference_count:
            treads = total_reads[ref_idx]
        if treads == 0:
            continue

        connected_neighbors = 0
        if ref_idx < array_size and exact_connection_counts:
            connected_neighbors = exact_connection_counts[ref_idx]

        sorted_refs[valid_count].neighbor_count = connected_neighbors
        sorted_refs[valid_count].read_count = treads
        sorted_refs[valid_count].idx = ref_idx

        if pattern_data and ref_idx < pool.reference_count:
            sorted_refs[valid_count].conc = <float>multimap_reads[ref_idx] / <float>total_reads[ref_idx] if total_reads[ref_idx] > 0 else 0.0
            sorted_refs[valid_count].net = pattern_data[ref_idx].graph.avg_comappings_per_read
        else:
            sorted_refs[valid_count].conc = 0.0
            sorted_refs[valid_count].net = 0.0

        valid_count += 1

    if valid_count > 1:
        qsort(sorted_refs, valid_count, sizeof(MultPairExtended), _multpair_extended_cmp_local)

    for i in range(valid_count):
        ref_idx = sorted_refs[i].idx

        ref_name_c = b"unknown"
        ref_len = 0

        if mapping and mapping.new_to_old_tid and ref_idx < mapping.n_retained_refs:
            tid32 = mapping.new_to_old_tid[ref_idx]
            if tid32 >= 0 and bam_header:
                ref_name_c = sam_hdr_tid2name(bam_header, tid32)
                ref_len = sam_hdr_tid2len(bam_header, tid32)

        if not ref_name_c:
            ref_name_c = b"unknown"

        treads = 0
        unique_reads = 0
        repeat_reads = 0
        shared_reads = 0
        if ref_stats and ref_idx < pool.reference_count:
            treads = ref_stats[ref_idx].total_reads
            unique_reads = ref_stats[ref_idx].unique_reads
            repeat_reads = ref_stats[ref_idx].repeat_reads
            shared_reads = ref_stats[ref_idx].shared_reads

        align_count = 0
        if ref_stats and ref_idx < pool.reference_count:
            align_count = ref_stats[ref_idx].alignment_count

        multimap_pct = 0.0
        if treads > 0:
            multimap_pct = 100.0 * <double>(repeat_reads + shared_reads) / <double>treads

        connected_neighbors = 0
        avg_co_mappings = 0.0
        max_co_map = 0
        if ref_idx < array_size:
            if exact_connection_counts:
                connected_neighbors = exact_connection_counts[ref_idx]
            if co_mapping_averages:
                avg_co_mappings = co_mapping_averages[ref_idx]
            elif pattern_data and ref_idx < pool.reference_count:
                avg_co_mappings = pattern_data[ref_idx].graph.avg_comappings_per_read

            if max_co_mappings:
                max_co_map = max_co_mappings[ref_idx]
            elif pattern_data and ref_idx < pool.reference_count:
                max_co_map = pattern_data[ref_idx].graph.max_comappings

        neighbor_mm_rate = 0.0
        if ref_idx < array_size and neighbor_multimap_avg:
            neighbor_mm_rate = neighbor_multimap_avg[ref_idx] * 100.0

        score_mean = 0.0
        score_std = 0.0
        score_min = 0.0
        score_max = 0.0
        if ref_stats and ref_idx < pool.reference_count:
            score_mean = ref_stats[ref_idx].score_mean
            score_std = ref_stats[ref_idx].score_std
            score_min = ref_stats[ref_idx].score_min
            score_max = ref_stats[ref_idx].score_max

        pmd_available = 0
        pmd_mean = 0.0
        pmd_std = 0.0
        pmd_min = 0.0
        pmd_max = 0.0
        pmd_nonzero_pct = 0
        if ref_stats and ref_idx < pool.reference_count:
            pmd_available = 1 if ref_stats[ref_idx].pmd_available else 0
            pmd_mean = ref_stats[ref_idx].pmd_mean
            pmd_std = ref_stats[ref_idx].pmd_std
            pmd_min = ref_stats[ref_idx].pmd_min
            pmd_max = ref_stats[ref_idx].pmd_max
            if align_count > 0:
                pmd_nonzero_pct = <uint64_t>(100.0 * <double>ref_stats[ref_idx].pmd_nonzero_count / <double>align_count)

        if avg_co_mappings != avg_co_mappings or avg_co_mappings > 1e9 or avg_co_mappings < -1e9:
            avg_co_mappings = 0.0

        component_id = <uint32_t>UINT32_MAX
        community_id = <uint32_t>UINT32_MAX
        community_cc = 0.0
        individual_cc = 0.0
        keep_flag = 0
        node_degree = 0
        betweenness = 0.0
        taxid = -1
        taxid_rank_id = -1
        taxid_depth = -1
        taxonomy_flag = 0

        # Tier result variables
        structural_role = 0
        num_neighbor_communities = 0
        community_coherent = 0
        filter_decision = 0
        taxonomy_outlier_score = 0

        # Misannotation detection variables
        misannotation_flag = 0
        misannotation_confidence = 0.0
        cross_domain_edges_before = 0
        edges_after_removal = 0
        cross_domain_fraction = 0.0

        # Initialize lineage names
        name_superkingdom = NULL
        name_clade = NULL
        name_kingdom = NULL
        name_phylum = NULL
        name_class = NULL
        name_order = NULL
        name_family = NULL
        name_genus = NULL
        name_species = NULL
        name_subspecies = NULL

        if pattern_data and ref_idx < pool.reference_count:
            component_id = pattern_data[ref_idx].component_id
            node_degree = pattern_data[ref_idx].node_degree
            community_id = pattern_data[ref_idx].community_id
            community_cc = pattern_data[ref_idx].community_cc
            individual_cc = pattern_data[ref_idx].community_individual_cc
            keep_flag = pattern_data[ref_idx].community_keep_flag
            betweenness = pattern_data[ref_idx].betweenness_centrality
            taxid = pattern_data[ref_idx].taxid
            taxid_rank_id = pattern_data[ref_idx].taxid_rank_id
            taxid_depth = pattern_data[ref_idx].taxid_depth
            taxonomy_flag = pattern_data[ref_idx].taxonomy_flag

            # Load per-rank taxonomy counts
            tax_neighbors_total = pattern_data[ref_idx].tax_neighbors_total
            tax_mismatch_domain = pattern_data[ref_idx].tax_mismatch_domain
            tax_mismatch_kingdom = pattern_data[ref_idx].tax_mismatch_kingdom
            tax_mismatch_phylum = pattern_data[ref_idx].tax_mismatch_phylum
            tax_mismatch_class = pattern_data[ref_idx].tax_mismatch_class
            tax_mismatch_order = pattern_data[ref_idx].tax_mismatch_order
            tax_mismatch_family = pattern_data[ref_idx].tax_mismatch_family
            tax_match_genus_below = pattern_data[ref_idx].tax_match_genus_below

            # Load tier results
            structural_role = pattern_data[ref_idx].structural_role
            num_neighbor_communities = pattern_data[ref_idx].num_neighbor_communities
            community_coherent = pattern_data[ref_idx].community_coherent
            filter_decision = pattern_data[ref_idx].filter_decision
            taxonomy_outlier_score = pattern_data[ref_idx].taxonomy_outlier_score

            # Load misannotation detection fields
            misannotation_flag = pattern_data[ref_idx].misannotation_flag
            misannotation_confidence = pattern_data[ref_idx].misannotation_confidence
            cross_domain_edges_before = pattern_data[ref_idx].cross_domain_edges_before
            edges_after_removal = pattern_data[ref_idx].edges_after_removal
            cross_domain_fraction = pattern_data[ref_idx].cross_domain_fraction

            # Extract lineage using cache (MUCH faster for repeated taxids)
            if taxonomy_db != NULL and taxid > 0 and cached_taxids != NULL:
                # Check cache first
                found_in_cache = 0
                for cache_idx in range(cache_size):
                    if cached_taxids[cache_idx] == taxid:
                        # Cache hit! Reuse the lineage data
                        name_superkingdom = cached_lineages[cache_idx].name_superkingdom
                        name_clade = cached_lineages[cache_idx].name_clade
                        name_kingdom = cached_lineages[cache_idx].name_kingdom
                        name_phylum = cached_lineages[cache_idx].name_phylum
                        name_class = cached_lineages[cache_idx].name_class
                        name_order = cached_lineages[cache_idx].name_order
                        name_family = cached_lineages[cache_idx].name_family
                        name_genus = cached_lineages[cache_idx].name_genus
                        name_species = cached_lineages[cache_idx].name_species
                        name_subspecies = cached_lineages[cache_idx].name_subspecies
                        found_in_cache = 1
                        break

                # Cache miss - extract and add to cache
                if not found_in_cache:
                    extract_lineage_data(taxonomy_db, taxid, &temp_lineage)
                    name_superkingdom = temp_lineage.name_superkingdom
                    name_clade = temp_lineage.name_clade
                    name_kingdom = temp_lineage.name_kingdom
                    name_phylum = temp_lineage.name_phylum
                    name_class = temp_lineage.name_class
                    name_order = temp_lineage.name_order
                    name_family = temp_lineage.name_family
                    name_genus = temp_lineage.name_genus
                    name_species = temp_lineage.name_species
                    name_subspecies = temp_lineage.name_subspecies

                    # Add to cache if space available
                    if cache_size < cache_capacity:
                        cached_taxids[cache_size] = taxid
                        cached_lineages[cache_size] = temp_lineage
                        cache_size += 1

        if gzfp != NULL:
            if taxonomy_db != NULL:
                # WITH TAXONOMY
                if community_id == UINT32_MAX:
                    if keep_flag:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            tax_neighbors_total, tax_mismatch_domain, tax_mismatch_kingdom, tax_mismatch_phylum, tax_mismatch_class, tax_mismatch_order, tax_mismatch_family, tax_match_genus_below,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            tax_neighbors_total, tax_mismatch_domain, tax_mismatch_kingdom, tax_mismatch_phylum, tax_mismatch_class, tax_mismatch_order, tax_mismatch_family, tax_match_genus_below,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                else:
                    # WITH TAXONOMY - community (non-singleton)
                    if keep_flag:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            tax_neighbors_total, tax_mismatch_domain, tax_mismatch_kingdom, tax_mismatch_phylum, tax_mismatch_class, tax_mismatch_order, tax_mismatch_family, tax_match_genus_below,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%u\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            tax_neighbors_total, tax_mismatch_domain, tax_mismatch_kingdom, tax_mismatch_phylum, tax_mismatch_class, tax_mismatch_order, tax_mismatch_family, tax_match_genus_below,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
            else:
                # WITHOUT TAXONOMY
                if community_id == UINT32_MAX:
                    if keep_flag:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                else:
                    # WITHOUT TAXONOMY - community (non-singleton)
                    if keep_flag:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        gzprintf(gzfp, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
        else:
            # NON-GZIPPED (fprintf)
            if taxonomy_db != NULL:
                # WITH TAXONOMY
                if community_id == UINT32_MAX:
                    if keep_flag:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                else:
                    # WITH TAXONOMY - community (non-singleton)
                    if keep_flag:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%d\t%d\t%d\t%d\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            taxonomy_outlier_score, taxid, taxid_rank_id, taxid_depth, <int>taxonomy_flag,
                            (name_superkingdom if name_superkingdom != NULL else <const char*>b""), (name_clade if name_clade != NULL else <const char*>b""), (name_kingdom if name_kingdom != NULL else <const char*>b""),
                            (name_phylum if name_phylum != NULL else <const char*>b""), (name_class if name_class != NULL else <const char*>b""),
                            (name_order if name_order != NULL else <const char*>b""), (name_family if name_family != NULL else <const char*>b""),
                            (name_genus if name_genus != NULL else <const char*>b""), (name_species if name_species != NULL else <const char*>b""),
                            (name_subspecies if name_subspecies != NULL else <const char*>b""),
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
            else:
                # WITHOUT TAXONOMY
                if community_id == UINT32_MAX:
                    if keep_flag:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tsingleton\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                else:
                    # WITHOUT TAXONOMY - community (non-singleton)
                    if keep_flag:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)
                    else:
                        fprintf(f, b"%s\t%ld\tcomp_%u\t%u\tcomm_%u\t%.6f\t%.6f\t%.6f\t%s\t%u\t%d\t%s\t%u\t%u\t%u\t%u\t%lu\t%.2f\t%u\t%.2f\t%lu\t%.2f\t%.6f\t%.6f\t%.6f\t%.6f\t%u\t%.6f\t%.6f\t%.6f\t%.6f\t%lu\t%u\t%.3f\t%u\t%u\t%.3f\n",
                            ref_name_c, ref_len, component_id, node_degree, community_id, community_cc, individual_cc, betweenness,
                            role_names[<int>structural_role], num_neighbor_communities, <int>community_coherent, decision_names[<int>filter_decision],
                            treads, unique_reads, repeat_reads, shared_reads, align_count, multimap_pct,
                            connected_neighbors, avg_co_mappings, max_co_map, neighbor_mm_rate,
                            score_mean, score_std, score_min, score_max, pmd_available,
                            pmd_mean, pmd_std, pmd_min, pmd_max, pmd_nonzero_pct,
                            misannotation_flag, misannotation_confidence,
                            cross_domain_edges_before, edges_after_removal, cross_domain_fraction)

    # Cleanup lineage cache
    if cached_taxids != NULL:
        free(cached_taxids)
    if cached_lineages != NULL:
        free(cached_lineages)

    free(sorted_refs)
    if gzfp != NULL:
        gzclose(gzfp)
    if f != NULL:
        fclose(f)
    return 0
