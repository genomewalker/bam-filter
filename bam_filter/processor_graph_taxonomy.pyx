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

"""
Taxonomy-aware graph analysis implementation.

This module enriches graph analysis with taxonomic information to identify:
- Contamination (references connecting across domains/kingdoms)
- Misassemblies (chimeric contigs connecting distantly related taxa)
- Multi-mapping artifacts (reads incorrectly mapped to distant relatives)
"""

from libc.stdlib cimport malloc, free, calloc
from libc.string cimport strcmp, strncpy, strlen, strchr, strstr
from libc.stdio cimport printf, fprintf, stderr, snprintf
from libc.stdint cimport int32_t, uint32_t, int64_t, uint64_t
from libc.math cimport sqrt as libc_sqrt, fabs

from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor_graph_taxonomy cimport (
    TaxonomyGraphConfig,
    extract_accession_from_refname,
    enrich_patterns_with_taxonomy,
    detect_taxonomy_anomalies,
    compute_lca_between_refs,
    get_taxonomy_flag_name
)
from bam_filter.processor_graph_ops cimport WeightedGraph, GraphNode
from bam_filter.taxonomy_db cimport (
    TaxonomyDB, AccessionMap,
    compute_lca_nogil, TaxonomyDatabase, AccessionMapping,
    lookup_taxid_duckdb, get_taxid_at_rank_nogil
)

# Rank IDs for domain/superkingdom
cdef int32_t RANK_SUPERKINGDOM = 24
cdef int32_t RANK_DOMAIN = 25

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

# khash types for accession lookup
cdef extern from "taxonomy_khash.h":
    ctypedef unsigned int khint_t
    ctypedef struct kh_str_t "khash_t(str)":
        khint_t n_buckets
        const char** keys
        int32_t* vals
    khint_t kh_get_str(kh_str_t* h, const char* key) nogil
    int kh_exist(kh_str_t* h, khint_t x) nogil
    khint_t kh_size(kh_str_t* h) nogil


# =============================================================================
# Accession extraction from reference names
# =============================================================================

cdef int32_t extract_accession_from_refname(const char* refname, char* accession_buf,
                                             int32_t buf_size) noexcept nogil:
    """
    Extract accession from reference name.

    Handles common formats:
    - GenBank: NC_000001.11
    - RefSeq: NZ_CP123456.1
    - EMBL: AL123456
    - Simple: just the accession
    - With descriptions: "NC_000001.11 Homo sapiens chromosome 1"
    - With pipes: "gi|123|ref|NC_000001.11|"

    Parameters
    ----------
    refname : const char*
        Reference name from BAM header
    accession_buf : char*
        Buffer to store extracted accession
    buf_size : int32_t
        Size of accession_buf

    Returns
    -------
    int32_t
        Length of extracted accession, or -1 if extraction failed
    """
    # Declare all variables at the top
    cdef int32_t len_refname
    cdef const char* pipe_pos
    cdef const char* last_pipe = NULL
    cdef const char* second_last_pipe = NULL
    cdef const char* current_pos
    cdef int32_t acc_len
    cdef const char* space_pos
    cdef const char* tab_pos
    cdef const char* end_pos

    if refname == NULL or accession_buf == NULL or buf_size < 2:
        return -1

    len_refname = <int32_t>strlen(refname)
    if len_refname == 0:
        return -1

    # Strategy 1: Look for pipe-delimited format (e.g., "gi|123|ref|NC_000001.11|")
    pipe_pos = strchr(refname, <int>b'|')
    current_pos = refname

    if pipe_pos != NULL:
        # Count pipes and find the accession field
        while current_pos != NULL and current_pos < refname + len_refname:
            pipe_pos = strchr(current_pos, <int>b'|')
            if pipe_pos == NULL:
                break
            second_last_pipe = last_pipe
            last_pipe = pipe_pos
            current_pos = pipe_pos + 1

        # Accession is typically between second-last and last pipe
        if second_last_pipe != NULL and last_pipe != NULL:
            acc_len = <int32_t>(last_pipe - second_last_pipe - 1)
            if acc_len > 0 and acc_len < buf_size:
                strncpy(accession_buf, second_last_pipe + 1, acc_len)
                accession_buf[acc_len] = 0
                return acc_len

    # Strategy 2: Extract first word (before space or tab)
    space_pos = strchr(refname, <int>b' ')
    tab_pos = strchr(refname, <int>b'\t')
    end_pos = NULL

    if space_pos != NULL and tab_pos != NULL:
        end_pos = space_pos if space_pos < tab_pos else tab_pos
    elif space_pos != NULL:
        end_pos = space_pos
    elif tab_pos != NULL:
        end_pos = tab_pos
    else:
        end_pos = refname + len_refname

    acc_len = <int32_t>(end_pos - refname)
    if acc_len > 0 and acc_len < buf_size:
        strncpy(accession_buf, refname, acc_len)
        accession_buf[acc_len] = 0
        return acc_len

    return -1


# =============================================================================
# Taxonomy enrichment
# =============================================================================

cdef int enrich_patterns_with_taxonomy(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    sam_hdr_t* bam_header,
    TaxonomyDB* taxdb,
    AccessionMap* accmap,
    ReferenceMapping* mapping,
    bint verbose
) noexcept nogil:
    """
    Enrich reference patterns with taxonomy information.

    For each reference:
    1. Extract accession from reference name
    2. Look up taxid from accession map
    3. Query taxonomy database for rank and depth
    4. Store in pattern_data

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns to enrich (modified in-place)
    n_refs : uint32_t
        Number of references
    bam_header : sam_hdr_t*
        BAM header for reference name lookup
    taxdb : TaxonomyDB*
        Taxonomy database
    accmap : AccessionMap*
        Accession to taxid mapping
    mapping : ReferenceMapping*
        Reference ID mapping (may be NULL)
    verbose : bint
        Enable verbose logging

    Returns
    -------
    int
        0 on success, -1 on error
    """
    if pattern_data == NULL or taxdb == NULL or accmap == NULL:
        return -1

    # Additional safety checks for taxonomy database structure
    if taxdb.taxid_to_idx == NULL or taxdb.nodes == NULL:
        if verbose:
            printf("ERROR: Taxonomy database has NULL internal structures\n")
        return -1

    cdef uint32_t ref_idx
    cdef int32_t tid32, orig_tid
    cdef const char* refname
    cdef char accession_buf[256]
    cdef int32_t acc_len
    cdef int32_t taxid
    cdef int32_t taxid_idx
    cdef int32_t rank_id
    cdef int32_t depth
    cdef uint32_t found_count = 0
    cdef uint32_t not_found_count = 0
    cdef kh_str_t* acc_hash = <kh_str_t*>accmap.acc_hash
    cdef khint_t k

    if verbose:
        printf("Enriching %u references with taxonomy information...\n", n_refs)
        printf("  Taxonomy DB: %d nodes, max_taxid=%d\n", taxdb.n_nodes, taxdb.max_taxid)

    for ref_idx in range(n_refs):
        # Initialize taxonomy fields to "not found"
        pattern_data[ref_idx].taxid = -1
        pattern_data[ref_idx].taxid_rank_id = -1
        pattern_data[ref_idx].taxid_depth = -1
        pattern_data[ref_idx].domain_taxid = -1
        pattern_data[ref_idx].taxonomy_flag = 0

        # Get reference name using the same approach as GraphML/TSV export
        refname = NULL
        if mapping and mapping.new_to_old_tid and ref_idx < mapping.n_retained_refs:
            tid32 = <int32_t>mapping.new_to_old_tid[ref_idx]
            if tid32 >= 0 and bam_header:
                refname = sam_hdr_tid2name(bam_header, tid32)
        elif bam_header:
            # Fallback if no mapping is provided
            refname = sam_hdr_tid2name(bam_header, <int32_t>ref_idx)

        if refname == NULL:
            not_found_count += 1
            continue

        # Extract accession
        acc_len = extract_accession_from_refname(refname, accession_buf, 256)
        if acc_len <= 0:
            not_found_count += 1
            continue

        # Look up taxid - check if using khash or DuckDB backend
        if acc_hash != NULL:
            # khash lookup
            k = kh_get_str(acc_hash, accession_buf)
            if not kh_exist(acc_hash, k):
                not_found_count += 1
                continue
            taxid = acc_hash.vals[k]
        elif accmap.duckdb_conn != NULL:
            # DuckDB on-demand lookup
            taxid = lookup_taxid_duckdb(accmap, accession_buf)
            if taxid < 0:
                not_found_count += 1
                continue
        else:
            # No valid backend
            not_found_count += 1
            continue
        if taxid < 0 or taxid > taxdb.max_taxid:
            not_found_count += 1
            continue

        taxid_idx = taxdb.taxid_to_idx[taxid]
        if taxid_idx < 0:
            not_found_count += 1
            continue

        # Store taxonomy information
        pattern_data[ref_idx].taxid = taxid
        pattern_data[ref_idx].taxid_rank_id = taxdb.nodes[taxid_idx].rank_id
        pattern_data[ref_idx].taxid_depth = taxdb.nodes[taxid_idx].depth
        # Get domain-level taxid for entropy computation (superkingdom rank = 24)
        pattern_data[ref_idx].domain_taxid = get_taxid_at_rank_nogil(taxdb, taxid, RANK_SUPERKINGDOM)
        found_count += 1

    if verbose:
        printf("  Taxonomy enrichment: %u found, %u not found (%.1f%% coverage)\n",
               found_count, not_found_count,
               100.0 * <float>found_count / <float>(found_count + not_found_count))

    return 0


# =============================================================================
# Taxonomy-based anomaly detection
# =============================================================================

cdef int detect_taxonomy_anomalies(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    MemoryPool* pool,
    TaxonomyDB* taxdb,
    TaxonomyGraphConfig* config,
    uint32_t** neighbor_lists,
    uint32_t* neighbor_counts,
    void* graph_handle,
    bint verbose
) noexcept nogil:
    """
    Detect taxonomy-based anomalies in graph connectivity.

    Flags references with suspicious connectivity patterns:
    - Cross-domain connections (e.g., Bacteria <-> Eukaryota)
    - Kingdom mismatches (e.g., Animalia <-> Fungi)
    - Genus-level mismatches with high connectivity

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array of reference patterns with taxonomy info
    n_refs : uint32_t
        Number of references
    pool : MemoryPool*
        Memory pool with alignment data
    taxdb : TaxonomyDB*
        Taxonomy database
    config : TaxonomyGraphConfig*
        Configuration for anomaly detection thresholds
    neighbor_lists : uint32_t**
        Per-reference neighbor lists (from graph analysis)
    neighbor_counts : uint32_t*
        Per-reference neighbor counts
    verbose : bint
        Enable verbose logging

    Returns
    -------
    int
        0 on success, -1 on error
    """
    if pattern_data == NULL or taxdb == NULL or config == NULL:
        return -1

    if not config.enabled:
        return 0

    cdef uint32_t ref_idx, neighbor_idx, neighbor_ref_idx
    cdef int32_t taxid1, taxid2, lca_taxid
    cdef int32_t rank_id1, rank_id2, lca_rank_id
    cdef int32_t depth1, depth2, lca_depth
    cdef uint32_t n_neighbors
    cdef uint32_t cross_domain_count = 0
    cdef uint32_t kingdom_mismatch_count = 0
    cdef uint32_t genus_mismatch_count = 0
    cdef uint32_t total_edges_analyzed = 0
    cdef float cross_domain_fraction
    cdef float kingdom_mismatch_fraction
    cdef float genus_mismatch_fraction
    cdef int32_t lca_idx
    cdef uint32_t flagged_count = 0
    cdef uint32_t flag_counts[4]

    # Rank IDs for comparison (from RANK_TO_ID in taxonomy_db.pyx)
    cdef int32_t RANK_DOMAIN = 25
    cdef int32_t RANK_SUPERKINGDOM = 24
    cdef int32_t RANK_KINGDOM = 23
    cdef int32_t RANK_PHYLUM = 20
    cdef int32_t RANK_CLASS = 17
    cdef int32_t RANK_ORDER = 13
    cdef int32_t RANK_FAMILY = 8
    cdef int32_t RANK_GENUS = 6

    if verbose:
        printf("Detecting taxonomy-based anomalies in graph connectivity...\n")
        printf("  n_refs=%u neighbor_lists=%p neighbor_counts=%p graph_handle=%p\n",
               n_refs, <void*>neighbor_lists, <void*>neighbor_counts, graph_handle)
        printf("  Configuration:\n")
        printf("    Min rank for comparison: %d\n", config.min_rank_id_for_comparison)
        printf("    Cross-domain threshold: %.2f%%\n", config.cross_domain_threshold * 100.0)
        printf("    Kingdom mismatch threshold: %.2f%%\n", config.kingdom_mismatch_threshold * 100.0)
        printf("    Genus mismatch threshold: %.2f%%\n", config.genus_mismatch_threshold * 100.0)

        # Sample a few neighbor counts to verify data
        if neighbor_counts != NULL and n_refs > 0:
            printf("  Sample neighbor_counts: [0]=%u [100]=%u [1000]=%u\n",
                   neighbor_counts[0],
                   neighbor_counts[100] if n_refs > 100 else 0,
                   neighbor_counts[1000] if n_refs > 1000 else 0)

    # Initialize per-rank counters for all references
    for ref_idx in range(n_refs):
        pattern_data[ref_idx].tax_neighbors_total = 0
        pattern_data[ref_idx].tax_mismatch_domain = 0
        pattern_data[ref_idx].tax_mismatch_kingdom = 0
        pattern_data[ref_idx].tax_mismatch_phylum = 0
        pattern_data[ref_idx].tax_mismatch_class = 0
        pattern_data[ref_idx].tax_mismatch_order = 0
        pattern_data[ref_idx].tax_mismatch_family = 0
        pattern_data[ref_idx].tax_match_genus_below = 0

    # Analyze each reference's neighbors
    for ref_idx in range(n_refs):
        taxid1 = pattern_data[ref_idx].taxid
        if taxid1 < 0:
            continue  # No taxonomy info for this reference

        rank_id1 = pattern_data[ref_idx].taxid_rank_id
        depth1 = pattern_data[ref_idx].taxid_depth

        if neighbor_lists == NULL or neighbor_counts == NULL:
            continue

        n_neighbors = neighbor_counts[ref_idx]
        if n_neighbors == 0:
            continue

        # Count taxonomic mismatches with neighbors (both old-style and per-rank)
        cross_domain_count = 0
        kingdom_mismatch_count = 0
        genus_mismatch_count = 0

        for neighbor_idx in range(n_neighbors):
            neighbor_ref_idx = neighbor_lists[ref_idx][neighbor_idx]
            if neighbor_ref_idx >= n_refs:
                continue

            taxid2 = pattern_data[neighbor_ref_idx].taxid
            if taxid2 < 0:
                continue  # No taxonomy info for neighbor

            rank_id2 = pattern_data[neighbor_ref_idx].taxid_rank_id
            depth2 = pattern_data[neighbor_ref_idx].taxid_depth

            # Compute LCA
            lca_taxid = compute_lca_nogil(taxdb, taxid1, taxid2)
            if lca_taxid < 0 or lca_taxid > taxdb.max_taxid:
                continue

            lca_idx = taxdb.taxid_to_idx[lca_taxid]
            if lca_idx < 0:
                continue

            lca_rank_id = taxdb.nodes[lca_idx].rank_id
            lca_depth = taxdb.nodes[lca_idx].depth

            total_edges_analyzed += 1
            pattern_data[ref_idx].tax_neighbors_total += 1

            # Classify by LCA rank (for per-rank analysis)
            # IMPORTANT: Rank IDs are HIGHER for LESS specific ranks (domain=25, genus=6, species=2)
            # So we use >= for high-level comparisons and <= for low-level comparisons
            if lca_rank_id >= RANK_SUPERKINGDOM:
                # Domain/Superkingdom level or above (most severe: LCA is at or above domain)
                pattern_data[ref_idx].tax_mismatch_domain += 1
                cross_domain_count += 1
            elif lca_rank_id == RANK_KINGDOM:
                # Kingdom level (LCA is exactly at kingdom)
                pattern_data[ref_idx].tax_mismatch_kingdom += 1
                kingdom_mismatch_count += 1
            elif lca_rank_id == RANK_PHYLUM:
                # Phylum level (LCA is exactly at phylum)
                pattern_data[ref_idx].tax_mismatch_phylum += 1
                genus_mismatch_count += 1
            elif lca_rank_id == RANK_CLASS:
                # Class level (LCA is exactly at class)
                pattern_data[ref_idx].tax_mismatch_class += 1
                genus_mismatch_count += 1
            elif lca_rank_id == RANK_ORDER:
                # Order level (LCA is exactly at order)
                pattern_data[ref_idx].tax_mismatch_order += 1
                genus_mismatch_count += 1
            elif lca_rank_id == RANK_FAMILY:
                # Family level (LCA is exactly at family)
                pattern_data[ref_idx].tax_mismatch_family += 1
                genus_mismatch_count += 1
            elif lca_rank_id <= RANK_GENUS:
                # Genus or below (good match: LCA is at genus, species, or subspecies)
                pattern_data[ref_idx].tax_match_genus_below += 1
            else:
                # Unclassified rank between family and genus
                genus_mismatch_count += 1

        # Calculate mismatch fractions
        if n_neighbors > 0:
            cross_domain_fraction = <float>cross_domain_count / <float>n_neighbors
            kingdom_mismatch_fraction = <float>kingdom_mismatch_count / <float>n_neighbors
            genus_mismatch_fraction = <float>genus_mismatch_count / <float>n_neighbors

            # Flag based on thresholds
            if cross_domain_fraction >= config.cross_domain_threshold:
                pattern_data[ref_idx].taxonomy_flag = 2  # cross_domain
            elif kingdom_mismatch_fraction >= config.kingdom_mismatch_threshold:
                pattern_data[ref_idx].taxonomy_flag = 3  # kingdom_mismatch
            elif genus_mismatch_fraction >= config.genus_mismatch_threshold:
                pattern_data[ref_idx].taxonomy_flag = 1  # potential_contamination

    # Check pruned isolated nodes using their original neighbors
    # These nodes have degree=0 now but had edges before pruning
    cdef WeightedGraph* graph = <WeightedGraph*>graph_handle
    cdef GraphNode* node
    cdef uint32_t orig_degree
    cdef uint32_t pruned_isolated_checked = 0

    if graph and graph.nodes:
        for ref_idx in range(n_refs):
            # Skip if already flagged or no taxonomy
            if pattern_data[ref_idx].taxonomy_flag > 0:
                continue

            taxid1 = pattern_data[ref_idx].taxid
            if taxid1 < 0:
                continue

            # Check if this is a pruned isolated node
            if ref_idx < graph.num_nodes:
                node = &graph.nodes[ref_idx]
                if node.degree == 0 and node.original_degree > 0 and node.original_neighbors:
                    # Pruned isolated node - check original neighbors
                    pruned_isolated_checked += 1
                    orig_degree = node.original_degree

                    # Count taxonomic mismatches with original neighbors
                    cross_domain_count = 0
                    kingdom_mismatch_count = 0
                    genus_mismatch_count = 0

                    rank_id1 = pattern_data[ref_idx].taxid_rank_id
                    depth1 = pattern_data[ref_idx].taxid_depth

                    for neighbor_idx in range(orig_degree):
                        neighbor_ref_idx = node.original_neighbors[neighbor_idx]
                        if neighbor_ref_idx >= n_refs:
                            continue

                        taxid2 = pattern_data[neighbor_ref_idx].taxid
                        if taxid2 < 0:
                            continue

                        rank_id2 = pattern_data[neighbor_ref_idx].taxid_rank_id
                        depth2 = pattern_data[neighbor_ref_idx].taxid_depth

                        # Compute LCA
                        lca_taxid = compute_lca_nogil(taxdb, taxid1, taxid2)
                        if lca_taxid < 0 or lca_taxid > taxdb.max_taxid:
                            continue

                        lca_idx = taxdb.taxid_to_idx[lca_taxid]
                        if lca_idx < 0:
                            continue

                        lca_rank_id = taxdb.nodes[lca_idx].rank_id

                        # Count per-rank for pruned isolated nodes
                        pattern_data[ref_idx].tax_neighbors_total += 1

                        # Classify by LCA rank
                        if lca_rank_id <= RANK_DOMAIN or lca_rank_id == RANK_SUPERKINGDOM:
                            pattern_data[ref_idx].tax_mismatch_domain += 1
                            cross_domain_count += 1
                        elif lca_rank_id == RANK_KINGDOM:
                            pattern_data[ref_idx].tax_mismatch_kingdom += 1
                            kingdom_mismatch_count += 1
                        elif lca_rank_id == RANK_PHYLUM:
                            pattern_data[ref_idx].tax_mismatch_phylum += 1
                            genus_mismatch_count += 1
                        elif lca_rank_id == RANK_CLASS:
                            pattern_data[ref_idx].tax_mismatch_class += 1
                            genus_mismatch_count += 1
                        elif lca_rank_id == RANK_ORDER:
                            pattern_data[ref_idx].tax_mismatch_order += 1
                            genus_mismatch_count += 1
                        elif lca_rank_id == RANK_FAMILY:
                            pattern_data[ref_idx].tax_mismatch_family += 1
                            genus_mismatch_count += 1
                        elif lca_rank_id >= RANK_GENUS:
                            pattern_data[ref_idx].tax_match_genus_below += 1
                        else:
                            genus_mismatch_count += 1

                    # Flag based on thresholds (using original neighbors)
                    if orig_degree > 0:
                        cross_domain_fraction = <float>cross_domain_count / <float>orig_degree
                        kingdom_mismatch_fraction = <float>kingdom_mismatch_count / <float>orig_degree
                        genus_mismatch_fraction = <float>genus_mismatch_count / <float>orig_degree

                        if cross_domain_fraction >= config.cross_domain_threshold:
                            pattern_data[ref_idx].taxonomy_flag = 2  # cross_domain
                        elif kingdom_mismatch_fraction >= config.kingdom_mismatch_threshold:
                            pattern_data[ref_idx].taxonomy_flag = 3  # kingdom_mismatch
                        elif genus_mismatch_fraction >= config.genus_mismatch_threshold:
                            pattern_data[ref_idx].taxonomy_flag = 1  # potential_contamination

        if verbose and pruned_isolated_checked > 0:
            printf("  Checked %u pruned isolated nodes using original neighbors\n", pruned_isolated_checked)

    if verbose:
        flagged_count = 0
        flag_counts[0] = 0
        flag_counts[1] = 0
        flag_counts[2] = 0
        flag_counts[3] = 0

        for ref_idx in range(n_refs):
            if pattern_data[ref_idx].taxonomy_flag > 0:
                flagged_count += 1
                flag_counts[<int>pattern_data[ref_idx].taxonomy_flag] += 1

        printf("  Taxonomy anomaly detection results:\n")
        printf("    Total edges analyzed: %u\n", total_edges_analyzed)
        printf("    Flagged references: %u / %u (%.1f%%)\n",
               flagged_count, n_refs, 100.0 * <float>flagged_count / <float>n_refs)
        printf("      Potential contamination: %u\n", flag_counts[1])
        printf("      Cross-domain: %u\n", flag_counts[2])
        printf("      Kingdom mismatch: %u\n", flag_counts[3])

    return 0


# =============================================================================
# Helper functions
# =============================================================================

cdef int32_t compute_lca_between_refs(
    int32_t taxid1,
    int32_t taxid2,
    TaxonomyDB* taxdb
) noexcept nogil:
    """Compute LCA between two references."""
    return compute_lca_nogil(taxdb, taxid1, taxid2)


cdef const char* get_taxonomy_flag_name(char flag) noexcept nogil:
    """Get human-readable name for taxonomy flag."""
    if flag == 0:
        return b"normal"
    elif flag == 1:
        return b"potential_contamination"
    elif flag == 2:
        return b"cross_domain"
    elif flag == 3:
        return b"kingdom_mismatch"
    else:
        return b"unknown"


# =============================================================================
# Python-level interface (for testing and standalone use)
# =============================================================================

def py_extract_accession(str refname):
    """
    Python wrapper for accession extraction (for testing).

    Parameters
    ----------
    refname : str
        Reference name

    Returns
    -------
    str or None
        Extracted accession, or None if extraction failed
    """
    cdef bytes refname_bytes = refname.encode('utf-8')
    cdef const char* refname_c = refname_bytes
    cdef char accession_buf[256]
    cdef int32_t acc_len

    with nogil:
        acc_len = extract_accession_from_refname(refname_c, accession_buf, 256)

    if acc_len > 0:
        return accession_buf[:acc_len].decode('utf-8')
    else:
        return None
