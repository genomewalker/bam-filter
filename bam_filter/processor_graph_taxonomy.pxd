# cython: language_level=3

"""
Taxonomy-aware graph analysis header.

This module provides functions for enriching graph analysis with taxonomic
information and detecting contamination/misassembly patterns based on
taxonomic incongruence.
"""

from libc.stdint cimport int32_t, uint32_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.taxonomy_db cimport TaxonomyDB, AccessionMap

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil


# Taxonomy-aware graph analysis configuration
cdef struct TaxonomyGraphConfig:
    bint enabled                        # Enable taxonomy-aware analysis
    int32_t min_rank_id_for_comparison  # Minimum taxonomic rank to use for comparison (e.g., genus=6)
    float cross_domain_threshold        # Threshold for flagging cross-domain connections
    float kingdom_mismatch_threshold    # Threshold for flagging kingdom mismatches
    float genus_mismatch_threshold      # Threshold for flagging genus mismatches


# Extract accession from reference name (handles various formats)
cdef int32_t extract_accession_from_refname(const char* refname, char* accession_buf,
                                             int32_t buf_size) noexcept nogil

# Enrich reference patterns with taxonomy information
cdef int enrich_patterns_with_taxonomy(
    ReferencePattern* pattern_data,
    uint32_t n_refs,
    sam_hdr_t* bam_header,
    TaxonomyDB* taxdb,
    AccessionMap* accmap,
    ReferenceMapping* mapping,
    bint verbose
) noexcept nogil

# Detect taxonomy-based anomalies in graph connectivity
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
) noexcept nogil

# Helper: Compute LCA between two references
cdef int32_t compute_lca_between_refs(
    int32_t taxid1,
    int32_t taxid2,
    TaxonomyDB* taxdb
) noexcept nogil

# Helper: Get rank name for logging
cdef const char* get_taxonomy_flag_name(char flag) noexcept nogil
