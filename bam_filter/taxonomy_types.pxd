# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False

"""
Type definitions for high-performance taxonomy database.
"""

from libc.stdint cimport int32_t, uint32_t, int64_t, uint64_t

# Taxonomy node structure (compact representation)
cdef struct TaxNode:
    int32_t taxid           # Taxonomy ID
    int32_t parent_taxid    # Parent taxonomy ID
    int32_t rank_id         # Rank ID (enumerated for fast comparison)
    int32_t name_offset     # Offset into names string array
    int32_t name_length     # Length of name string
    int32_t depth           # Depth in tree (root = 0)
    int32_t lca_table_idx   # Index into LCA precomputation table (-1 if not cached)

# Taxonomy database structure
cdef struct TaxonomyDB:
    TaxNode* nodes          # Array of taxonomy nodes (indexed by internal ID)
    int32_t* taxid_to_idx   # Hash table: taxid -> internal index
    char* names_buffer      # Concatenated name strings
    char** rank_names       # Array of rank name strings
    int32_t n_nodes         # Number of nodes
    int32_t n_ranks         # Number of distinct ranks
    int32_t max_taxid       # Maximum taxid value
    int32_t root_idx        # Index of root node

# Accession to taxid mapping structure
cdef struct AccessionMap:
    char** accessions       # Array of accession strings
    int32_t* taxids         # Corresponding taxid array
    int32_t n_entries       # Number of entries

# LCA precomputation cache structure
cdef struct LCACache:
    int32_t* lca_matrix     # Flattened matrix of precomputed LCAs
    int32_t* cached_taxids  # List of taxids in cache
    int32_t n_cached        # Number of cached taxids
    int32_t matrix_size     # Size of matrix (n_cached * n_cached)
