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
Three-stage filtering pipeline combining structural roles, taxonomy coherence,
and aggregated decision rules for graph-based contamination detection.
"""

from libc.stdio cimport printf
from libc.stdint cimport uint32_t, int32_t, uint8_t, uint16_t, int64_t, uint64_t
from libc.stdlib cimport malloc, free, calloc
from libc.string cimport memset
from libc.math cimport fabs

from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats
from bam_filter.taxonomy_db cimport TaxonomyDB, compute_lca_nogil
from bam_filter.processor cimport MemoryPool, Alignment

# Constant for invalid community ID
cdef uint32_t UINT32_MAX = 0xFFFFFFFF

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil


# ==============================================================================
# TIER 1: Structural Role Classification
# ==============================================================================

cdef StructuralRole classify_structural_role(
    float betweenness,
    float clustering_coefficient,
    uint32_t degree,
    uint32_t num_neighbor_communities,
    float betweenness_threshold,
    float cc_threshold,
    uint32_t hub_degree_threshold
) nogil:
    """
    Classify a node's structural role in the graph.

    Decision tree:
    1. High betweenness + crosses communities → BRIDGE (most suspicious)
    2. High degree + low CC → HUB (suspicious)
    3. High degree + high CC → CORE (conserved, likely clean)
    4. Low degree → PERIPHERAL (specific, likely clean)

    Parameters
    ----------
    betweenness : float
        Betweenness centrality (normalized)
    clustering_coefficient : float
        Local clustering coefficient
    degree : uint32_t
        Node degree (number of neighbors)
    num_neighbor_communities : uint32_t
        Number of distinct communities among neighbors
    betweenness_threshold : float
        Threshold for bridge detection (e.g., 0.01)
    cc_threshold : float
        Threshold for hub detection (e.g., 0.3)
    hub_degree_threshold : uint32_t
        Minimum degree for hub classification (e.g., 5)

    Returns
    -------
    StructuralRole
        Classified role
    """
    # Priority 1: look for bridges connecting multiple communities
    # Strategy:
    # 1. Require num_neighbor_communities >= 2
    # 2. Use betweenness when it is available to refine the decision

    if num_neighbor_communities >= 2:
        # Connects 2+ communities - candidate for BRIDGE
        if betweenness > betweenness_threshold:
            # Strong bridge signal (high betweenness confirms)
            return ROLE_BRIDGE
        elif betweenness > 0:
            # Moderate bridge signal (betweenness present but low)
            # This could be a weak bridge or boundary node
            if num_neighbor_communities >= 3:
                # Connects many communities even if betweenness low -> still BRIDGE
                return ROLE_BRIDGE
            # Connects only 2 communities with low betweenness -> check other criteria
        else:
            # No betweenness calculated (betweenness=0)
            # Fall back to num_neighbor_communities only
            return ROLE_BRIDGE

    # Priority 2: Check for HUB (promiscuous within community)
    if degree >= hub_degree_threshold and clustering_coefficient < cc_threshold:
        return ROLE_HUB

    # Priority 3: Check for CORE (highly connected but cohesive)
    if degree >= hub_degree_threshold and clustering_coefficient >= cc_threshold:
        return ROLE_CORE

    # Default: PERIPHERAL (low connectivity)
    return ROLE_PERIPHERAL


# ==============================================================================
# Taxonomy Outlier Score: Detect likely misannotations in databases
# ==============================================================================

cdef uint8_t compute_taxonomy_outlier_score(
    uint32_t node_degree,
    uint8_t taxonomy_flag,
    bint community_coherent,
    StructuralRole structural_role
) nogil:
    """
    Compute a score (0-100) indicating likelihood of taxonomic misannotation.

    High scores indicate a reference is likely misannotated in the database.

    Evidence for misannotation:
    - High connectivity (shares reads with many references)
    - Severe taxonomy mismatch with neighbors (taxonomy_flag = 2 or 3)
    - Community is incoherent (mixed taxonomy)
    - Structural role is BRIDGE or HUB (connects different groups)

    Score interpretation:
    - 0-20:   Low confidence misannotation (taxonomy seems correct)
    - 21-40:  Mild suspicion
    - 41-60:  Moderate confidence misannotation
    - 61-80:  High confidence misannotation
    - 81-100: Very high confidence (likely database error)

    Parameters
    ----------
    node_degree : uint32_t
        Number of connected neighbors
    taxonomy_flag : uint8_t
        0=normal, 1=genus mismatch, 2=cross-domain, 3=kingdom mismatch
    community_coherent : bint
        True if community has coherent taxonomy
    structural_role : StructuralRole
        PERIPHERAL, CORE, HUB, or BRIDGE

    Returns
    -------
    uint8_t
        Score from 0-100
    """
    cdef float score = 0.0

    # Component 1: Connectivity (max 30 points)
    # High connectivity + taxonomy problems = likely misannotation
    if node_degree >= 20:
        score += 30.0
    elif node_degree >= 10:
        score += 20.0
    elif node_degree >= 5:
        score += 10.0
    # Low degree (<5) gets 0 points - not enough evidence

    # Component 2: Taxonomy flag severity (max 40 points)
    # This is the strongest signal
    if taxonomy_flag == 3:  # Kingdom mismatch
        score += 40.0
    elif taxonomy_flag == 2:  # Cross-domain
        score += 40.0
    elif taxonomy_flag == 1:  # Genus mismatch
        score += 15.0
    # taxonomy_flag == 0 (normal) gets 0 points

    # Component 3: Community incoherence (max 20 points)
    if not community_coherent:
        score += 20.0

    # Component 4: Structural role (max 10 points)
    # BRIDGE/HUB with taxonomy problems = very suspicious
    if structural_role == ROLE_BRIDGE:
        score += 10.0
    elif structural_role == ROLE_HUB:
        score += 5.0

    # Cap at 100
    if score > 100.0:
        score = 100.0

    return <uint8_t>score


# ==============================================================================
# TIER 2: Community Coherence Analysis
# ==============================================================================

cdef bint check_community_coherence(
    uint32_t* community_members,
    uint32_t n_members,
    ReferencePattern* pattern_data,
    TaxonomyDB* taxonomy_db,
    int32_t* lca_result_out
) nogil:
    """
    Check if a community is taxonomically coherent.

    A community is coherent if all members share a reasonably close LCA
    (e.g., same phylum or closer).

    Example of COHERENT community:
        Member 1: Escherichia coli       (Bacteria;Proteobacteria;Gamma;Enterobacterales;...)
        Member 2: Salmonella enterica    (Bacteria;Proteobacteria;Gamma;Enterobacterales;...)
        Member 3: Klebsiella pneumoniae  (Bacteria;Proteobacteria;Gamma;Enterobacterales;...)
        LCA: Order level (Enterobacterales) → COHERENT ✓

    Example of INCOHERENT community:
        Member 1: Escherichia coli       (Bacteria;Proteobacteria;...)
        Member 2: Homo sapiens           (Eukaryota;Chordata;...)
        Member 3: Aspergillus fumigatus  (Eukaryota;Ascomycota;...)
        LCA: Root level → INCOHERENT ✗ (cross-domain contamination!)

    Parameters
    ----------
    community_members : uint32_t*
        Array of reference indices in this community
    n_members : uint32_t
        Number of members
    pattern_data : ReferencePattern*
        Reference patterns with taxonomy info
    taxonomy_db : TaxonomyDB*
        Taxonomy database for LCA computation
    lca_result_out : int32_t*
        Output: LCA taxid of the community (can be NULL)

    Returns
    -------
    bint
        True if coherent, False if incoherent
    """
    if taxonomy_db == NULL or n_members == 0:
        return True  # Can't check, assume coherent

    # Collect taxids for all members with valid taxonomy
    cdef int32_t taxids[256]  # Stack allocation, should be enough
    cdef uint32_t valid_taxid_count = 0
    cdef uint32_t i
    cdef int32_t taxid
    cdef uint32_t member_idx

    for i in range(n_members):
        if i >= 256:  # Safety check
            break
        member_idx = community_members[i]
        taxid = pattern_data[member_idx].taxid

        if taxid > 0:
            taxids[valid_taxid_count] = taxid
            valid_taxid_count += 1

    if valid_taxid_count == 0:
        return True  # No taxonomy info, assume coherent

    if valid_taxid_count == 1:
        if lca_result_out != NULL:
            lca_result_out[0] = taxids[0]
        return True  # Single taxid, trivially coherent

    # Compute LCA of all taxids
    cdef int32_t lca_taxid = taxids[0]
    for i in range(1, valid_taxid_count):
        lca_taxid = compute_lca_nogil(taxonomy_db, lca_taxid, taxids[i])
        if lca_taxid < 0:
            return False  # LCA computation failed, assume incoherent

    if lca_result_out != NULL:
        lca_result_out[0] = lca_taxid

    # Check LCA depth - if too shallow (near root), community is incoherent
    # Use ONLY depth-based checking to support custom taxonomies
    cdef int32_t lca_idx = taxonomy_db.taxid_to_idx[lca_taxid]
    if lca_idx < 0:
        return False

    cdef int32_t lca_depth = taxonomy_db.nodes[lca_idx].depth

    # Coherence criteria (depth-based only, works with any taxonomy):
    # - Depth ≥ 3: Below root/superkingdom level → COHERENT
    #   This allows families/orders/genera clustering together (expected biology)
    # - Depth < 3: Root or domain level → INCOHERENT (cross-domain contamination)
    #
    # Depth levels (typical):
    #   0 = root
    #   1 = superkingdom/domain (Bacteria, Eukaryota, Archaea)
    #   2 = kingdom/phylum (varies by taxonomy)
    #   3+ = more specific (class, order, family, genus, species)

    if lca_depth >= 3:
        return True  # Deep enough, coherent
    else:
        return False  # Too shallow, incoherent


# ==============================================================================
# TIER 3: Integrated Decision Matrix
# ==============================================================================

cdef FilterDecision make_filtering_decision(
    StructuralRole role,
    uint8_t taxonomy_flag,
    bint community_coherent,
    bint strict_mode
) nogil:
    """
    Make an integrated filtering decision combining all evidence.

    Decision Matrix:

    Role       | Taxonomy   | Community  | Strict | Decision
    -----------|------------|------------|--------|----------
    PERIPHERAL | Normal     | Coherent   | Any    | KEEP
    PERIPHERAL | Abnormal   | Coherent   | No     | KEEP (second chance)
    PERIPHERAL | Abnormal   | Coherent   | Yes    | REMOVE
    PERIPHERAL | Any        | Incoherent | Any    | REMOVE

    CORE       | Normal     | Coherent   | Any    | KEEP (conserved gene)
    CORE       | Abnormal   | Coherent   | No     | REVIEW (needs inspection)
    CORE       | Abnormal   | Coherent   | Yes    | REMOVE
    CORE       | Any        | Incoherent | Any    | REMOVE

    HUB        | Normal     | Coherent   | No     | REVIEW (may be conserved)
    HUB        | Normal     | Coherent   | Yes    | REMOVE
    HUB        | Abnormal   | Any        | Any    | REMOVE
    HUB        | Any        | Incoherent | Any    | REMOVE

    BRIDGE     | Normal     | Coherent   | No     | REVIEW (unusual topology)
    BRIDGE     | Normal     | Coherent   | Yes    | REMOVE
    BRIDGE     | Abnormal   | Any        | Any    | REMOVE (critical!)
    BRIDGE     | Any        | Incoherent | Any    | REMOVE (critical!)

    Taxonomy flags:
    - 0 = normal
    - 1 = genus mismatch
    - 2 = cross-domain
    - 3 = kingdom mismatch

    Parameters
    ----------
    role : StructuralRole
        Structural role from Tier 1
    taxonomy_flag : uint8_t
        Taxonomy anomaly flag
    community_coherent : bint
        Is community taxonomically coherent?
    strict_mode : bint
        Use strict filtering (less permissive)

    Returns
    -------
    FilterDecision
        KEEP, REMOVE, or REVIEW
    """
    cdef bint taxonomy_normal = (taxonomy_flag == 0)
    cdef bint taxonomy_severe = (taxonomy_flag == 2 or taxonomy_flag == 3)  # cross-domain or kingdom

    # BRIDGE: Most suspicious - nearly always remove
    if role == ROLE_BRIDGE:
        if not taxonomy_normal or not community_coherent:
            return DECISION_REMOVE  # Bridge with bad taxonomy/community = definite contamination
        if strict_mode:
            return DECISION_REMOVE  # In strict mode, don't trust bridges even with good taxonomy
        return DECISION_REVIEW  # Non-strict: flag for review

    # HUB: Promiscuous connector
    if role == ROLE_HUB:
        if not taxonomy_normal or not community_coherent:
            return DECISION_REMOVE  # Hub with bad taxonomy/community = contamination
        if strict_mode:
            return DECISION_REMOVE  # In strict mode, remove hubs
        return DECISION_REVIEW  # Non-strict: may be conserved gene, review

    # CORE: High connectivity but cohesive
    if role == ROLE_CORE:
        if not community_coherent:
            return DECISION_REMOVE  # Incoherent community
        if taxonomy_severe:
            return DECISION_REMOVE  # Severe taxonomy issues
        if not taxonomy_normal:
            # Minor taxonomy issues (genus mismatch)
            if strict_mode:
                return DECISION_REMOVE
            return DECISION_REVIEW
        return DECISION_KEEP  # Normal taxonomy, coherent = likely conserved gene

    # PERIPHERAL: Low connectivity, specific
    if role == ROLE_PERIPHERAL:
        if not community_coherent:
            return DECISION_REMOVE  # Even peripheral nodes in incoherent communities are suspicious
        if taxonomy_severe:
            return DECISION_REMOVE  # Severe taxonomy issues
        if not taxonomy_normal:
            # Minor taxonomy issues
            if strict_mode:
                return DECISION_REMOVE
            return DECISION_KEEP  # Second chance for peripheral with minor issues
        return DECISION_KEEP  # Normal case

    # Default: keep (shouldn't reach here)
    return DECISION_KEEP


# ==============================================================================
# Edge Removal for Cross-Domain Contamination
# ==============================================================================

cdef float compute_misannotation_score(
    uint32_t degree_before,
    uint32_t degree_after,
    uint32_t cross_domain_count,
    uint8_t taxonomy_flag,
    float clustering_coefficient
) noexcept nogil:
    """
    Compute confidence score for database misannotation.

    Signals prioritized by strength:
    1. Lost ALL edges (degree_after == 0): Very strong signal (0.6 points)
    2. Lost most edges (>75%): Strong signal (0.3 points)
    3. High cross-domain fraction (>75%): Strong contamination (0.3 points)
    4. Medium cross-domain fraction (>50%): Moderate contamination (0.2 points)
    5. Taxonomy flag (cross_domain=2): Detected in Phase 6 (0.1 points)
    6. Hub with low CC: Possible misassembly (0.1 points)

    Returns
    -------
    float
        Confidence score 0.0-1.0
    """
    cdef float score = 0.0
    cdef float cd_fraction

    if degree_before == 0:
        return 0.0  # No data

    # Signal 1: Lost ALL edges (strongest signal - indicates 100% cross-domain)
    if degree_after == 0 and degree_before > 0:
        score += 0.6

    # Signal 2: Lost most edges (>75%)
    elif <float>degree_after < <float>degree_before * 0.25:
        score += 0.3

    # Signal 3 & 4: Cross-domain fraction
    cd_fraction = <float>cross_domain_count / <float>degree_before
    if cd_fraction > 0.75:
        score += 0.3
    elif cd_fraction > 0.50:
        score += 0.2

    # Signal 5: Taxonomy flag from Phase 6 detection
    if taxonomy_flag == 2:  # cross_domain
        score += 0.1

    # Signal 6: Hub with low clustering coefficient (possible misassembly)
    if degree_before > 20 and clustering_coefficient < 0.3:
        score += 0.1

    # Cap at 1.0
    if score > 1.0:
        score = 1.0

    return score


cdef uint32_t count_remaining_edges_after_removal(
    uint32_t ref_idx,
    char* alignment_keep_flags,
    MemoryPool* pool,
    uint32_t** neighbor_lists,
    uint32_t* neighbor_counts
) noexcept nogil:
    """
    Count how many edges remain after edge removal for a reference.

    For each neighbor, check if at least one read still connects them
    after alignments have been marked for removal.

    Parameters
    ----------
    ref_idx : uint32_t
        Reference to check
    alignment_keep_flags : char*
        Flags indicating which alignments survive (1=keep, 0=remove)
    pool : MemoryPool*
        Memory pool with alignments
    neighbor_lists : uint32_t**
        Neighbor lists for each reference
    neighbor_counts : uint32_t*
        Neighbor counts for each reference

    Returns
    -------
    uint32_t
        Number of neighbors still connected after edge removal
    """
    if neighbor_counts == NULL or neighbor_lists == NULL or pool == NULL:
        return 0

    cdef uint32_t n_neighbors = neighbor_counts[ref_idx]
    if n_neighbors == 0:
        return 0

    cdef uint32_t remaining = 0
    cdef uint32_t neighbor_idx, neighbor_ref_idx
    cdef int64_t aln_idx
    cdef uint32_t aln_ref_idx, read_idx
    cdef uint64_t read_start, read_count, other_aln_idx
    cdef bint has_connection

    # For each neighbor, check if any reads still connect them
    for neighbor_idx in range(n_neighbors):
        neighbor_ref_idx = neighbor_lists[ref_idx][neighbor_idx]
        has_connection = False

        # Check all alignments to ref_idx
        for aln_idx in range(pool.alignment_count):
            aln_ref_idx = pool.alignments[aln_idx].reference_index

            if aln_ref_idx != ref_idx:
                continue

            # Check if this alignment survives edge removal
            if alignment_keep_flags != NULL and not alignment_keep_flags[aln_idx]:
                continue

            # Get the read and check its other alignments
            read_idx = pool.alignments[aln_idx].read_index
            read_start = pool.read_alignment_starts[read_idx]
            read_count = pool.read_alignment_counts[read_idx]

            other_aln_idx = read_start
            while other_aln_idx < read_start + read_count:
                if pool.alignments[other_aln_idx].reference_index == neighbor_ref_idx:
                    # Check if the neighbor alignment also survives
                    if alignment_keep_flags == NULL or alignment_keep_flags[other_aln_idx]:
                        has_connection = True
                        break
                other_aln_idx += 1

            if has_connection:
                break

        if has_connection:
            remaining += 1

    return remaining


cdef struct EdgeRemovalStats:
    uint32_t cross_domain_edges_found
    int64_t alignments_removed
    uint32_t references_affected
    uint32_t references_lost_all_edges
    uint32_t references_lost_most_edges    # Lost >75% of edges


cdef int remove_cross_domain_edges_for_reference(
    MemoryPool* pool,
    uint32_t ref_idx,
    uint32_t* neighbor_list,
    uint32_t n_neighbors,
    ReferencePattern* pattern_data,
    TaxonomyDB* taxonomy_db,
    char* alignment_keep_flags,
    EdgeRemovalStats* stats
) noexcept nogil:
    """
    Remove alignments between this reference and its cross-domain neighbors.

    This function is called for references in INCOHERENT communities only.
    It removes individual alignments (edges) instead of entire references.

    Algorithm:
    1. For each neighbor of this reference
    2. Compute LCA between ref_idx and neighbor
    3. If LCA is at domain level (lca_rank_id >= RANK_SUPERKINGDOM)
    4. Mark all alignments between these two references for removal

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing all alignments
    ref_idx : uint32_t
        Reference index to process
    neighbor_list : uint32_t*
        Array of neighbor reference indices
    n_neighbors : uint32_t
        Number of neighbors
    pattern_data : ReferencePattern*
        Reference patterns with taxonomy info
    taxonomy_db : TaxonomyDB*
        Taxonomy database for LCA computation
    alignment_keep_flags : char*
        Flags for each alignment (1=keep, 0=remove), modified in-place
    stats : EdgeRemovalStats*
        Statistics structure to update

    Returns
    -------
    int
        0 on success, -1 on error
    """
    if pool == NULL or neighbor_list == NULL or pattern_data == NULL or taxonomy_db == NULL:
        return -1

    if alignment_keep_flags == NULL or stats == NULL:
        return -1

    # Get taxonomy for this reference
    cdef int32_t taxid1 = pattern_data[ref_idx].taxid
    if taxid1 < 0:
        return 0  # No taxonomy, can't compute LCA

    # Rank IDs for comparison
    cdef int32_t RANK_SUPERKINGDOM = 24

    # Build set of cross-domain neighbors
    cdef uint32_t* cross_domain_neighbors = <uint32_t*>malloc(n_neighbors * sizeof(uint32_t))
    if cross_domain_neighbors == NULL:
        return -1

    cdef uint32_t n_cross_domain = 0
    cdef uint32_t neighbor_idx, neighbor_ref_idx
    cdef int32_t taxid2, lca_taxid, lca_idx, lca_rank_id

    # Identify cross-domain neighbors
    for neighbor_idx in range(n_neighbors):
        neighbor_ref_idx = neighbor_list[neighbor_idx]

        taxid2 = pattern_data[neighbor_ref_idx].taxid
        if taxid2 < 0:
            continue  # No taxonomy

        # Compute LCA
        lca_taxid = compute_lca_nogil(taxonomy_db, taxid1, taxid2)
        if lca_taxid < 0 or lca_taxid > taxonomy_db.max_taxid:
            continue

        lca_idx = taxonomy_db.taxid_to_idx[lca_taxid]
        if lca_idx < 0:
            continue

        lca_rank_id = taxonomy_db.nodes[lca_idx].rank_id

        # Cross-domain: LCA at or above superkingdom level
        if lca_rank_id >= RANK_SUPERKINGDOM:
            cross_domain_neighbors[n_cross_domain] = neighbor_ref_idx
            n_cross_domain += 1
            stats.cross_domain_edges_found += 1

    # If no cross-domain neighbors, nothing to remove
    if n_cross_domain == 0:
        free(cross_domain_neighbors)
        return 0

    # Now iterate through alignments and mark those to cross-domain neighbors
    cdef int64_t aln_idx
    cdef uint32_t aln_ref_idx
    cdef int64_t removed_count = 0
    cdef bint found
    cdef uint32_t read_idx
    cdef uint64_t read_start
    cdef uint32_t read_count
    cdef uint64_t other_aln_idx
    cdef uint32_t other_ref_idx
    cdef uint32_t cd_idx

    for aln_idx in range(pool.alignment_count):
        # Check if this alignment involves ref_idx
        aln_ref_idx = pool.alignments[aln_idx].reference_index

        if aln_ref_idx != ref_idx:
            continue  # Not this reference

        # Check if the READ mapped to this alignment also maps to a cross-domain neighbor
        # We need to check all alignments of this read to find cross-domain pairs
        read_idx = pool.alignments[aln_idx].read_index
        read_start = pool.read_alignment_starts[read_idx]
        read_count = pool.read_alignment_counts[read_idx]

        # Check this read's other alignments
        other_aln_idx = read_start
        while other_aln_idx < read_start + read_count:
            other_ref_idx = pool.alignments[other_aln_idx].reference_index

            # Check if other_ref_idx is in cross_domain_neighbors
            found = False
            for cd_idx in range(n_cross_domain):
                if cross_domain_neighbors[cd_idx] == other_ref_idx:
                    found = True
                    break

            if found:
                # This alignment and the cross-domain alignment form a contaminating edge
                # Mark BOTH for removal
                if alignment_keep_flags[aln_idx] == 1:
                    alignment_keep_flags[aln_idx] = 0
                    removed_count += 1
                if alignment_keep_flags[other_aln_idx] == 1:
                    alignment_keep_flags[other_aln_idx] = 0
                    removed_count += 1

            other_aln_idx += 1

    stats.alignments_removed += removed_count
    if removed_count > 0:
        stats.references_affected += 1

    free(cross_domain_neighbors)
    return 0


# ==============================================================================
# Main Entry Point
# ==============================================================================

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
) nogil:
    """
    Apply three-tier filtering to all references.

    This replaces the old CC-only filtering with a comprehensive approach:
    - Tier 1: Classify structural role (including bridge detection!)
    - Tier 2: Check community coherence
    - Tier 3: Make integrated decision

    Parameters
    ----------
    pattern_data : ReferencePattern*
        Array with Community clustering results and taxonomy
    ref_stats : ReferenceStats*
        Reference statistics (unused currently)
    array_size : uint32_t
        Number of references
    keep_flag : char*
        Keep flags (1=keep, 0=remove), modified in-place
    taxonomy_db : TaxonomyDB*
        Taxonomy database for coherence checking
    betweenness_threshold : float
        Threshold for bridge detection (default: 0.01)
    cc_threshold : float
        Threshold for hub detection (default: 0.3)
    hub_degree_threshold : uint32_t
        Minimum degree for hub classification (default: 5)
    strict_mode : bint
        Use strict filtering (default: False)
    verbose : bint
        Enable verbose logging

    Returns
    -------
    int
        0 on success, -1 on error
    """
    if pattern_data == NULL or keep_flag == NULL:
        return -1

    # Check taxonomy availability
    cdef bint taxonomy_available = (taxonomy_db != NULL)

    if verbose:
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "Applying tier-based filtering\\n"
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  Taxonomy: %s\\n",
            b"AVAILABLE" if taxonomy_available else b"NOT AVAILABLE (Tier 2 disabled)"
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  Betweenness threshold: %.4f\\n",
            betweenness_threshold
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  CC threshold: %.3f\\n",
            cc_threshold
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  Hub degree threshold: %u\\n",
            hub_degree_threshold
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  Strict mode: %s\\n\\n",
            b"YES" if strict_mode else b"NO"
        )

    # Statistics counters
    cdef uint32_t count_peripheral = 0
    cdef uint32_t count_core = 0
    cdef uint32_t count_hub = 0
    cdef uint32_t count_bridge = 0
    cdef uint32_t count_removed = 0
    cdef uint32_t count_kept = 0
    cdef uint32_t count_review = 0
    cdef uint32_t count_individual_override = 0  # Count of refs rescued by individual taxonomy override
    cdef uint32_t flagged_confident = 0
    cdef uint32_t flagged_likely = 0
    cdef uint32_t flagged_warning = 0
    cdef uint32_t degree_before = 0
    cdef uint32_t degree_after = 0
    cdef uint32_t cross_domain_count = 0
    cdef float confidence = 0.0
    cdef float cd_frac = 0.0
    cdef uint8_t tax_flag = 0
    cdef float individual_cc_local = 0.0
    cdef uint32_t* ref_alignment_counts = NULL
    cdef int64_t* ref_alignment_offsets = NULL
    cdef int64_t* ref_alignment_index_list = NULL
    cdef uint32_t* neighbor_marker = NULL
    cdef int64_t total_surviving_alignments = 0
    cdef bint compute_degree_after = 0
    cdef uint32_t n_neighbors = 0
    cdef uint32_t neighbor_idx = 0
    cdef uint32_t neighbor_ref_idx = 0
    cdef uint32_t other_ref_idx = 0
    cdef uint32_t marker_idx = 0
    cdef int64_t start = 0
    cdef int64_t end = 0
    cdef int64_t idx_pos = 0
    cdef bint original_removed = False

    # ========================================================================
    # TIER 2 PREPARATION: Build community membership and check coherence
    # ========================================================================

    # First pass: Count unique communities and find max community ID
    cdef uint32_t max_community_id = 0
    cdef uint32_t num_communities_found = 0
    cdef uint32_t ref_idx
    cdef uint32_t community_id
    cdef uint32_t comm_idx, member_count
    cdef uint32_t* members
    cdef uint32_t n_members
    cdef int32_t lca_result

    for ref_idx in range(array_size):
        if keep_flag[ref_idx] == 0:
            continue
        community_id = pattern_data[ref_idx].community_id
        if community_id != UINT32_MAX and community_id > max_community_id:
            max_community_id = community_id

    # Allocate arrays based on actual max community ID (+ 1 for 0-indexed)
    cdef uint32_t num_communities = max_community_id + 1
    cdef uint32_t max_members_per_community = 256  # Reasonable limit
    cdef uint32_t* community_member_lists = NULL
    cdef uint32_t* community_sizes = NULL
    cdef char* community_coherence_cache = NULL

    if num_communities > 0:
        community_member_lists = <uint32_t*>malloc(num_communities * max_members_per_community * sizeof(uint32_t))
        community_sizes = <uint32_t*>malloc(num_communities * sizeof(uint32_t))
        community_coherence_cache = <char*>malloc(num_communities * sizeof(char))

        if community_member_lists == NULL or community_sizes == NULL or community_coherence_cache == NULL:
            if community_member_lists: free(community_member_lists)
            if community_sizes: free(community_sizes)
            if community_coherence_cache: free(community_coherence_cache)
            return -1  # Memory allocation failed

        # Initialize
        for comm_idx in range(num_communities):
            community_sizes[comm_idx] = 0
            community_coherence_cache[comm_idx] = 1  # Default: coherent

    # Second pass: Build membership lists
    for ref_idx in range(array_size):
        if keep_flag[ref_idx] == 0:
            continue

        community_id = pattern_data[ref_idx].community_id
        if community_id == UINT32_MAX:
            continue  # Singleton, skip

        if community_id >= num_communities:
            continue  # Out of range, skip (shouldn't happen)

        # Add this ref to the community's member list
        member_count = community_sizes[community_id]
        if member_count < max_members_per_community:
            community_member_lists[community_id * max_members_per_community + member_count] = ref_idx
            community_sizes[community_id] = member_count + 1
            if member_count == 0:
                num_communities_found += 1

    # Third pass: Check coherence for each community (if taxonomy available)
    if taxonomy_available and num_communities > 0:
        for comm_idx in range(num_communities):
            if community_sizes[comm_idx] == 0:
                continue  # No members

            # Get pointer to member list for this community
            members = &community_member_lists[comm_idx * max_members_per_community]
            n_members = community_sizes[comm_idx]
            lca_result = -1

            # Check coherence
            community_coherence_cache[comm_idx] = <char>check_community_coherence(
                members, n_members, pattern_data, taxonomy_db, &lca_result
            )

    if verbose:
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  Found %u communities, checked coherence for all\\n",
            num_communities_found
        )

    # ========================================================================
    # EDGE REMOVAL: Remove cross-domain edges in incoherent communities
    # ========================================================================
    cdef MemoryPool* pool = NULL
    cdef char* alignment_keep_flags = NULL
    cdef EdgeRemovalStats edge_stats
    cdef int64_t aln_idx, new_count, old_idx
    edge_stats.cross_domain_edges_found = 0
    edge_stats.alignments_removed = 0
    edge_stats.references_affected = 0
    edge_stats.references_lost_all_edges = 0
    edge_stats.references_lost_most_edges = 0

    if enable_edge_removal and pool_handle != NULL and neighbor_lists != NULL and neighbor_counts != NULL:
        pool = <MemoryPool*>pool_handle

        if verbose:
            bf_nogil_logf_notime(
                b"TIERED_FILTER",
                "Removing cross-domain edges\\n"
            )

        # Allocate alignment keep flags (initialized to 1 = keep all)
        alignment_keep_flags = <char*>calloc(pool.alignment_count, sizeof(char))
        if alignment_keep_flags == NULL:
            if verbose:
                bf_nogil_logf_notime(b"TIERED_FILTER", "ERROR: Failed to allocate alignment flags\\n")
        else:
            # Initialize all to 1 (keep)
            for aln_idx in range(pool.alignment_count):
                alignment_keep_flags[aln_idx] = 1

            # Process each reference in an incoherent community
            for ref_idx in range(array_size):
                if keep_flag[ref_idx] == 0:
                    continue  # Already filtered

                community_id = pattern_data[ref_idx].community_id
                if community_id == UINT32_MAX:
                    continue  # Singleton

                # Check if community is incoherent
                if community_id < num_communities and community_coherence_cache[community_id] == 0:
                    # Incoherent community - apply edge removal
                    if neighbor_counts[ref_idx] > 0:
                        remove_cross_domain_edges_for_reference(
                            pool,
                            ref_idx,
                            neighbor_lists[ref_idx],
                            neighbor_counts[ref_idx],
                            pattern_data,
                            taxonomy_db,
                            alignment_keep_flags,
                            &edge_stats
                        )

            # Report edge removal statistics (compaction will happen later in processor_filters)
            if edge_stats.alignments_removed > 0:
                if verbose:
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  Cross-domain edges found: %u\\n",
                        edge_stats.cross_domain_edges_found
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  Alignments marked for removal: %ld\\n",
                        edge_stats.alignments_removed
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  References affected: %u\\n",
                        edge_stats.references_affected
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  (Compaction will happen after reference filtering)\\n\\n"
                    )

            # ========================================================================
            # MISANNOTATION DETECTION: Analyze refs that lost edges
            # ========================================================================
            if flag_misannotations:
                if verbose:
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "Detecting potential misannotations\\n"
                    )

                flagged_confident = 0
                flagged_likely = 0
                flagged_warning = 0
                compute_degree_after = (
                    alignment_keep_flags != NULL
                    and neighbor_lists != NULL
                    and neighbor_counts != NULL
                    and pool != NULL
                )

                if compute_degree_after:
                    ref_alignment_counts = <uint32_t*>calloc(array_size, sizeof(uint32_t))
                    if ref_alignment_counts == NULL:
                        compute_degree_after = 0

                if compute_degree_after:
                    for aln_idx in range(pool.alignment_count):
                        if alignment_keep_flags[aln_idx]:
                            ref_idx = pool.alignments[aln_idx].reference_index
                            if ref_idx < array_size:
                                ref_alignment_counts[ref_idx] += 1

                    ref_alignment_offsets = <int64_t*>malloc((array_size + 1) * sizeof(int64_t))
                    if ref_alignment_offsets == NULL:
                        compute_degree_after = 0

                if compute_degree_after:
                    total_surviving_alignments = 0
                    for ref_idx in range(array_size):
                        ref_alignment_offsets[ref_idx] = total_surviving_alignments
                        total_surviving_alignments += ref_alignment_counts[ref_idx]
                    ref_alignment_offsets[array_size] = total_surviving_alignments

                    if total_surviving_alignments > 0:
                        ref_alignment_index_list = <int64_t*>malloc(total_surviving_alignments * sizeof(int64_t))
                        if ref_alignment_index_list == NULL:
                            compute_degree_after = 0
                    else:
                        ref_alignment_index_list = NULL

                if compute_degree_after:
                    for ref_idx in range(array_size):
                        ref_alignment_counts[ref_idx] = 0

                    for aln_idx in range(pool.alignment_count):
                        if alignment_keep_flags[aln_idx]:
                            ref_idx = pool.alignments[aln_idx].reference_index
                            if ref_idx < array_size:
                                ref_alignment_index_list[
                                    ref_alignment_offsets[ref_idx] + ref_alignment_counts[ref_idx]
                                ] = aln_idx
                                ref_alignment_counts[ref_idx] += 1

                    neighbor_marker = <uint32_t*>calloc(pool.reference_count, sizeof(uint32_t))
                    if neighbor_marker == NULL:
                        compute_degree_after = 0
                if not compute_degree_after:
                    if ref_alignment_counts != NULL:
                        free(ref_alignment_counts)
                        ref_alignment_counts = NULL
                    if ref_alignment_offsets != NULL:
                        free(ref_alignment_offsets)
                        ref_alignment_offsets = NULL
                    if ref_alignment_index_list != NULL:
                        free(ref_alignment_index_list)
                        ref_alignment_index_list = NULL
                    if neighbor_marker != NULL:
                        free(neighbor_marker)
                        neighbor_marker = NULL

                for ref_idx in range(array_size):
                    if keep_flag[ref_idx] == 0:
                        continue  # Already filtered

                    # Get degree before edge removal
                    degree_before = pattern_data[ref_idx].node_degree
                    if degree_before == 0:
                        # Initialize to clean for refs with no edges
                        pattern_data[ref_idx].misannotation_flag = 0
                        pattern_data[ref_idx].misannotation_confidence = 0.0
                        pattern_data[ref_idx].cross_domain_edges_before = 0
                        pattern_data[ref_idx].edges_after_removal = 0
                        pattern_data[ref_idx].cross_domain_fraction = 0.0
                        continue

                    if compute_degree_after and ref_alignment_counts != NULL and ref_alignment_offsets != NULL and neighbor_marker != NULL:
                        if neighbor_counts[ref_idx] == 0 or neighbor_lists[ref_idx] == NULL:
                            degree_after = 0
                        elif ref_alignment_counts[ref_idx] == 0 or ref_alignment_index_list == NULL:
                            degree_after = 0
                        else:
                            n_neighbors = neighbor_counts[ref_idx]
                            neighbor_alive = <char*>malloc(n_neighbors * sizeof(char))
                            if neighbor_alive != NULL:
                                memset(neighbor_alive, 0, n_neighbors * sizeof(char))
                                # Build marker for fast lookup
                                for neighbor_idx in range(n_neighbors):
                                    neighbor_ref_idx = neighbor_lists[ref_idx][neighbor_idx]
                                    if neighbor_ref_idx < pool.reference_count:
                                        neighbor_marker[neighbor_ref_idx] = neighbor_idx + 1

                                start = ref_alignment_offsets[ref_idx]
                                end = start + ref_alignment_counts[ref_idx]
                                for idx_pos in range(start, end):
                                    aln_idx = ref_alignment_index_list[idx_pos]
                                    read_idx = pool.alignments[aln_idx].read_index
                                    read_start = pool.read_alignment_starts[read_idx]
                                    read_count = pool.read_alignment_counts[read_idx]
                                    other_aln_idx = read_start
                                    while other_aln_idx < read_start + read_count:
                                        if alignment_keep_flags[other_aln_idx]:
                                            other_ref_idx = pool.alignments[other_aln_idx].reference_index
                                            if other_ref_idx != ref_idx and other_ref_idx < pool.reference_count:
                                                marker_idx = neighbor_marker[other_ref_idx]
                                                if marker_idx != 0:
                                                    neighbor_alive[marker_idx - 1] = <char>1
                                                    # no break; continue scanning to mark other neighbors
                                        other_aln_idx += 1

                                degree_after = 0
                                for neighbor_idx in range(n_neighbors):
                                    if neighbor_alive[neighbor_idx]:
                                        degree_after += 1

                                # Reset markers for this reference
                                for neighbor_idx in range(n_neighbors):
                                    neighbor_ref_idx = neighbor_lists[ref_idx][neighbor_idx]
                                    if neighbor_ref_idx < pool.reference_count:
                                        neighbor_marker[neighbor_ref_idx] = 0
                                free(neighbor_alive)
                            else:
                                degree_after = pattern_data[ref_idx].node_degree
                    else:
                        degree_after = pattern_data[ref_idx].node_degree

                    # Get cross-domain count and other metrics
                    cross_domain_count = pattern_data[ref_idx].tax_mismatch_domain
                    tax_flag = <uint8_t>pattern_data[ref_idx].taxonomy_flag
                    individual_cc_local = pattern_data[ref_idx].community_individual_cc

                    # Store counts for TSV export
                    pattern_data[ref_idx].cross_domain_edges_before = cross_domain_count
                    pattern_data[ref_idx].edges_after_removal = degree_after
                    cd_frac = <float>cross_domain_count / <float>degree_before if degree_before > 0 else 0.0
                    pattern_data[ref_idx].cross_domain_fraction = cd_frac

                    # Compute misannotation confidence score
                    confidence = compute_misannotation_score(
                        degree_before,
                        degree_after,
                        cross_domain_count,
                        tax_flag,
                        individual_cc_local
                    )

                    pattern_data[ref_idx].misannotation_confidence = confidence

                    # Classify by confidence
                    if confidence >= 0.8:
                        pattern_data[ref_idx].misannotation_flag = 3  # CONFIDENT
                        flagged_confident += 1
                    elif confidence >= 0.6:
                        pattern_data[ref_idx].misannotation_flag = 2  # LIKELY
                        flagged_likely += 1
                    elif confidence >= 0.4:
                        pattern_data[ref_idx].misannotation_flag = 1  # WARNING
                        flagged_warning += 1
                    else:
                        pattern_data[ref_idx].misannotation_flag = 0  # CLEAN

                    # Track statistics
                    if degree_after == 0 and degree_before > 0:
                        edge_stats.references_lost_all_edges += 1
                    elif degree_after < degree_before * 0.25:
                        edge_stats.references_lost_most_edges += 1

                if verbose:
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  CONFIDENT misannotations (>0.8): %u\\n",
                        flagged_confident
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  LIKELY misannotations (0.6-0.8):  %u\\n",
                        flagged_likely
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  WARNING misannotations (0.4-0.6): %u\\n",
                        flagged_warning
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  Lost ALL edges: %u\\n",
                        edge_stats.references_lost_all_edges
                    )
                    bf_nogil_logf_notime(
                        b"TIERED_FILTER",
                        "  Lost >75%% edges: %u\\n\\n",
                        edge_stats.references_lost_most_edges
                    )
            else:
                # Initialize misannotation fields to clean when not flagging
                for ref_idx in range(array_size):
                    pattern_data[ref_idx].misannotation_flag = 0
                    pattern_data[ref_idx].misannotation_confidence = 0.0
                    pattern_data[ref_idx].cross_domain_edges_before = 0
                    pattern_data[ref_idx].edges_after_removal = 0
                    pattern_data[ref_idx].cross_domain_fraction = 0.0

            # DON'T free alignment_keep_flags - return it to caller for combined compaction
            # free(alignment_keep_flags)  // REMOVED

    if pool != NULL and pool.stats != NULL:
        pool.stats.edge_removal_edges_found = <int64_t>edge_stats.cross_domain_edges_found
        pool.stats.edge_removal_alignments_removed = edge_stats.alignments_removed
        pool.stats.edge_removal_references_affected = <int64_t>edge_stats.references_affected
        pool.stats.edge_removal_refs_lost_all_edges = <int64_t>edge_stats.references_lost_all_edges
        pool.stats.edge_removal_refs_lost_most_edges = <int64_t>edge_stats.references_lost_most_edges
        pool.stats.misannotation_confident = <int64_t>flagged_confident
        pool.stats.misannotation_likely = <int64_t>flagged_likely
        pool.stats.misannotation_warning = <int64_t>flagged_warning
        pool.stats.misannotation_removed_total = <int64_t>count_removed
        pool.stats.misannotation_review_total = <int64_t>count_review

    # Process each reference
    cdef StructuralRole role
    cdef FilterDecision decision
    cdef uint8_t taxonomy_flag
    cdef uint8_t effective_taxonomy_flag
    cdef uint8_t mis_flag
    cdef bint community_coherent
    cdef uint32_t num_neighbor_communities
    cdef bint individual_taxonomy_override
    cdef uint16_t tax_neighbors, genus_matches
    cdef float individual_cc_check
    cdef bint lost_all_edges

    # ref_idx and community_id already declared above for first pass
    for ref_idx in range(array_size):
        original_removed = False
        if keep_flag[ref_idx] == 0:
            if enable_edge_removal and pattern_data[ref_idx].community_keep_flag == 0:
                keep_flag[ref_idx] = 1
                original_removed = True
            else:
                continue

        # Get node properties from pattern_data
        betweenness = pattern_data[ref_idx].betweenness_centrality
        cc = pattern_data[ref_idx].community_individual_cc
        degree = pattern_data[ref_idx].node_degree
        taxonomy_flag = pattern_data[ref_idx].taxonomy_flag
        mis_flag = pattern_data[ref_idx].misannotation_flag
        community_id = pattern_data[ref_idx].community_id
        # Get actual num_neighbor_communities from graph analysis (computed in Community clustering)
        num_neighbor_communities = pattern_data[ref_idx].num_neighbor_communities
        lost_all_edges = (
            pattern_data[ref_idx].original_degree > 0
            and pattern_data[ref_idx].cross_domain_edges_before > 0
            and pattern_data[ref_idx].edges_after_removal == 0
        )

        # SPECIAL HANDLING: Pruned isolated nodes (degree=0 but had edges before pruning)
        # These are aDNA low-coverage refs that only had weak co-mapping edges
        # Be conservative: keep if no strong taxonomy signal against them
        if degree == 0 and pattern_data[ref_idx].original_degree > 0:
            # Pruned isolated node detected
            if taxonomy_available and taxonomy_flag >= 2:
                if enable_edge_removal:
                    # After edge removal, rely on misannotation flag rather than removing reference
                    keep_flag[ref_idx] = 1
                    count_review += 1
                    pattern_data[ref_idx].structural_role = <char>ROLE_PERIPHERAL
                    pattern_data[ref_idx].community_coherent = 0
                    pattern_data[ref_idx].filter_decision = <char>DECISION_REVIEW
                else:
                    # Has severe taxonomy problems (cross-domain or kingdom mismatch)
                    # Even with low coverage, this is suspicious -> REMOVE
                    keep_flag[ref_idx] = 0
                    count_removed += 1
                continue  # Skip normal processing
            else:
                # No severe taxonomy problems -> KEEP (conservative for aDNA)
                keep_flag[ref_idx] = 1
                count_kept += 1
                # Still classify for TSV export
                pattern_data[ref_idx].structural_role = <char>ROLE_PERIPHERAL
                pattern_data[ref_idx].community_coherent = 1
                pattern_data[ref_idx].filter_decision = <char>DECISION_KEEP
                pattern_data[ref_idx].taxonomy_outlier_score = 0
                continue  # Skip normal processing

        # TIER 1: Classify structural role
        role = classify_structural_role(
            betweenness, cc, degree, num_neighbor_communities,
            betweenness_threshold, cc_threshold, hub_degree_threshold
        )

        # Count roles
        if role == ROLE_PERIPHERAL:
            count_peripheral += 1
        elif role == ROLE_CORE:
            count_core += 1
        elif role == ROLE_HUB:
            count_hub += 1
        elif role == ROLE_BRIDGE:
            count_bridge += 1

        # TIER 2: Check community coherence (using cached results)
        if not taxonomy_available:
            community_coherent = True  # No taxonomy = assume coherent (conservative)
        elif community_id == UINT32_MAX:
            community_coherent = True  # Singleton = coherent by definition
        elif num_communities > 0 and community_id < num_communities:
            # Look up cached coherence result
            community_coherent = <bint>community_coherence_cache[community_id]
        else:
            community_coherent = True  # Out of range or no communities, assume coherent

        # TIER 2.5: Individual taxonomy override for incoherent communities
        # If community is globally incoherent (cross-domain mixing) BUT this specific reference
        # has 100% genus-level matches with its direct neighbors AND high CC, rescue it
        individual_taxonomy_override = False
        if not community_coherent and taxonomy_available:
            tax_neighbors = pattern_data[ref_idx].tax_neighbors_total
            genus_matches = pattern_data[ref_idx].tax_match_genus_below
            individual_cc_check = pattern_data[ref_idx].community_individual_cc

            # Criteria for override:
            # 1. Has neighbors with taxonomy info (tax_neighbors_total > 0)
            # 2. 100% of neighbors are genus-level matches (all neighbors taxonomically compatible)
            # 3. High clustering coefficient (well-integrated locally, CC >= 0.9)
            if tax_neighbors > 0 and genus_matches == tax_neighbors and individual_cc_check >= 0.9:
                individual_taxonomy_override = True
                community_coherent = True  # Override: treat as locally coherent
                count_individual_override += 1

        # TIER 3: Make decision
        # If no taxonomy, treat all taxonomy flags as normal
        effective_taxonomy_flag = taxonomy_flag if taxonomy_available else 0
        decision = make_filtering_decision(
            role, effective_taxonomy_flag, community_coherent, strict_mode
        )

        if enable_edge_removal and decision == DECISION_REMOVE:
            decision = DECISION_REVIEW
            if original_removed:
                # Ensure alignment removal statistics reflect restored reference
                pattern_data[ref_idx].community_keep_flag = 0

        if enable_edge_removal:
            if mis_flag >= 2 or lost_all_edges:
                decision = DECISION_REMOVE

        # Store tier results in pattern_data for TSV/GraphML export
        pattern_data[ref_idx].structural_role = <char>role
        pattern_data[ref_idx].community_coherent = <char>community_coherent
        pattern_data[ref_idx].filter_decision = <char>decision
        # Note: num_neighbor_communities already set from Community results, no need to store again

        # Compute taxonomy outlier score (0-100 indicating likelihood of database misannotation)
        pattern_data[ref_idx].taxonomy_outlier_score = compute_taxonomy_outlier_score(
            degree, effective_taxonomy_flag, community_coherent, role
        )

        # Apply decision
        if decision == DECISION_REMOVE:
            keep_flag[ref_idx] = 0
            count_removed += 1
        elif decision == DECISION_KEEP:
            keep_flag[ref_idx] = 1
            count_kept += 1
        elif decision == DECISION_REVIEW:
            # Treat REVIEW as KEEP (user can filter TSV by filter_decision column)
            keep_flag[ref_idx] = 1
            count_review += 1

    if verbose:
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "Structural role distribution:\\n"
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  PERIPHERAL: %u (specific, low connectivity)\\n",
            count_peripheral
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  CORE:       %u (conserved, cohesive)\\n",
            count_core
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  HUB:        %u (promiscuous, suspicious)\\n",
            count_hub
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  BRIDGE:     %u (inter-community connector)\\n",
            count_bridge
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "Filtering decisions:\\n"
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  KEPT:    %u\\n",
            count_kept
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  REMOVED: %u\\n",
            count_removed
        )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "  REVIEW:  %u (flagged for inspection)\\n",
            count_review
        )
        if count_individual_override > 0:
            bf_nogil_logf_notime(
                b"TIERED_FILTER",
                "Individual taxonomy rescue:\\n"
            )
            bf_nogil_logf_notime(
                b"TIERED_FILTER",
                "  Rescued %u refs from incoherent communities\\n",
                count_individual_override
            )
            bf_nogil_logf_notime(
                b"TIERED_FILTER",
                "  (100%% genus-level neighbor matches + CC >= 0.9)\\n"
            )
        bf_nogil_logf_notime(
            b"TIERED_FILTER",
            "\\n"
        )

    if ref_alignment_index_list != NULL:
        free(ref_alignment_index_list)
        ref_alignment_index_list = NULL
    if ref_alignment_offsets != NULL:
        free(ref_alignment_offsets)
        ref_alignment_offsets = NULL
    if ref_alignment_counts != NULL:
        free(ref_alignment_counts)
        ref_alignment_counts = NULL
    if neighbor_marker != NULL:
        free(neighbor_marker)
        neighbor_marker = NULL

    # Cleanup allocated memory
    if community_member_lists: free(community_member_lists)
    if community_sizes: free(community_sizes)
    if community_coherence_cache: free(community_coherence_cache)

    # Return alignment_keep_flags to caller for combined compaction
    if out_alignment_keep_flags != NULL:
        out_alignment_keep_flags[0] = alignment_keep_flags
    else:
        # Caller doesn't want flags, free them
        if alignment_keep_flags:
            free(alignment_keep_flags)

    return 0
