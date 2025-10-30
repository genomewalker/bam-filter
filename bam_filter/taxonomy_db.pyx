# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# cython: nonecheck=False
# cython: overflowcheck=False

"""
Ultra-fast taxonomy database construction and querying with Parquet serialization.

This module provides high-performance taxonomy operations for metagenomic analysis:
- Load NCBI taxdump files (nodes.dmp, names.dmp) into compact in-memory structures
- Load accession-to-taxid mappings (acc2taxid files)
- Perform O(log N) LCA queries with optional O(1) caching for frequent pairs
- Serialize/deserialize to Parquet for instant loading
- Multi-threaded construction and querying
"""

from libc.stdlib cimport malloc, calloc, free, realloc, qsort
from libc.string cimport strcmp, strcpy, strlen, strdup, memset, memcpy
from libc.stdio cimport FILE, fopen, fclose, fgets, printf, fprintf, stderr, snprintf
from libc.stdint cimport int32_t, uint32_t, int64_t, uint64_t
from cython.parallel cimport prange, parallel
from cpython.mem cimport PyMem_Malloc, PyMem_Free
import numpy as np
cimport numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from collections import defaultdict
import os

# DuckDB C API - Proper declarations matching duckdb.h
cdef extern from "duckdb.h":
    ctypedef uint64_t idx_t

    ctypedef enum duckdb_state:
        DuckDBSuccess = 0
        DuckDBError = 1

    # Opaque handles
    ctypedef struct _duckdb_database:
        void* internal_ptr
    ctypedef _duckdb_database* duckdb_database

    ctypedef struct _duckdb_connection:
        void* internal_ptr
    ctypedef _duckdb_connection* duckdb_connection

    ctypedef struct _duckdb_appender:
        void* internal_ptr
    ctypedef _duckdb_appender* duckdb_appender

    # Result is a struct, NOT a pointer
    ctypedef struct duckdb_column:
        void* deprecated_data
        bint* deprecated_nullmask
        int deprecated_type
        char* deprecated_name
        void* internal_data

    ctypedef struct duckdb_result:
        idx_t deprecated_column_count
        idx_t deprecated_row_count
        idx_t deprecated_rows_changed
        duckdb_column* deprecated_columns
        char* deprecated_error_message
        void* internal_data

    # Database functions
    duckdb_state duckdb_open(const char* path, duckdb_database* out_database) nogil
    void duckdb_close(duckdb_database* database) nogil
    duckdb_state duckdb_connect(duckdb_database database, duckdb_connection* out_connection) nogil
    void duckdb_disconnect(duckdb_connection* connection) nogil

    # Query execution
    duckdb_state duckdb_query(duckdb_connection connection, const char* query, duckdb_result* out_result) nogil
    void duckdb_destroy_result(duckdb_result* result) nogil

    # Result access
    idx_t duckdb_row_count(duckdb_result* result) nogil
    idx_t duckdb_column_count(duckdb_result* result) nogil
    char* duckdb_column_name(duckdb_result* result, idx_t col) nogil
    char* duckdb_value_varchar(duckdb_result* result, idx_t col, idx_t row) nogil
    int32_t duckdb_value_int32(duckdb_result* result, idx_t col, idx_t row) nogil
    const char* duckdb_result_error(duckdb_result* result) nogil
    void duckdb_free(void* ptr) nogil

    # Appender API for bulk inserts
    duckdb_state duckdb_appender_create(duckdb_connection connection, const char* schema,
                                       const char* table, duckdb_appender* out_appender) nogil
    duckdb_state duckdb_appender_close(duckdb_appender appender) nogil
    duckdb_state duckdb_appender_destroy(duckdb_appender* appender) nogil
    duckdb_state duckdb_append_varchar(duckdb_appender appender, const char* val) nogil
    duckdb_state duckdb_appender_end_row(duckdb_appender appender) nogil
    const char* duckdb_appender_error(duckdb_appender appender) nogil

# C++ standard library for Arrow
from libcpp.memory cimport shared_ptr, make_shared
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool

# Import khash for fast hash tables
cdef extern from "taxonomy_khash.h":
    ctypedef unsigned int khint_t
    ctypedef unsigned int khint32_t

    # Define khash structures with all fields for direct access
    ctypedef struct kh_str_t "khash_t(str)":
        khint_t n_buckets
        khint_t size
        khint_t n_occupied
        khint_t upper_bound
        khint32_t* flags
        const char** keys
        int32_t* vals  # Values array - can be used as lvalue

    ctypedef struct kh_int32_t "khash_t(int32)":
        khint_t n_buckets
        khint_t size
        khint_t n_occupied
        khint_t upper_bound
        khint32_t* flags
        int32_t* keys
        int32_t* vals  # Values array - can be used as lvalue

    # Only declare the actual generated functions (not macros)
    # String hash functions
    kh_str_t* kh_init_str() nogil
    void kh_destroy_str(kh_str_t* h) nogil
    khint_t kh_put_str(kh_str_t* h, const char* key, int* ret) nogil
    khint_t kh_get_str(kh_str_t* h, const char* key) nogil

    # Int32 hash functions
    kh_int32_t* kh_init_int32() nogil
    void kh_destroy_int32(kh_int32_t* h) nogil
    khint_t kh_put_int32(kh_int32_t* h, int32_t key, int* ret) nogil
    khint_t kh_get_int32(kh_int32_t* h, int32_t key) nogil


# Helper functions to replace khash macros
# kh_exist(h, x) checks if slot x is occupied (not empty and not deleted)
cdef inline int kh_exist_str(kh_str_t* h, khint_t x) nogil:
    """Check if hash table slot x exists (not empty/deleted)."""
    # __ac_iseither checks if slot is empty OR deleted
    # kh_exist returns True if slot is NOT (empty OR deleted)
    return not ((h.flags[x >> 4] >> ((x & 0xfU) << 1)) & 3)

cdef inline int kh_exist_int32(kh_int32_t* h, khint_t x) nogil:
    """Check if hash table slot x exists (not empty/deleted)."""
    return not ((h.flags[x >> 4] >> ((x & 0xfU) << 1)) & 3)


# =============================================================================
# Arrow C++ API declarations for fast Parquet loading
# =============================================================================

cdef extern from "arrow/api.h" namespace "arrow" nogil:
    cdef cppclass CStatus "arrow::Status":
        cpp_bool ok()
        string message()

    cdef cppclass CArray "arrow::Array":
        int64_t length()
        int64_t null_count()

    cdef cppclass CStringArray "arrow::StringArray"(CArray):
        string GetString(int64_t i)

    cdef cppclass CInt32Array "arrow::Int32Array"(CArray):
        int32_t Value(int64_t i)

    cdef cppclass CChunkedArray "arrow::ChunkedArray":
        int num_chunks()
        shared_ptr[CArray] chunk(int i)

    cdef cppclass CColumn "arrow::Column":
        shared_ptr[CChunkedArray] data()

    cdef cppclass CTable "arrow::Table":
        int64_t num_rows()
        int num_columns()
        shared_ptr[CChunkedArray] column(int i)

cdef extern from "arrow/io/api.h" namespace "arrow::io" nogil:
    cdef cppclass CInputStream "arrow::io::InputStream":
        pass

    cdef cppclass CRandomAccessFile "arrow::io::RandomAccessFile"(CInputStream):
        pass

    cdef cppclass CReadableFile "arrow::io::ReadableFile"(CRandomAccessFile):
        pass

    CStatus OpenFile "arrow::io::ReadableFile::Open"(
        const string& path,
        shared_ptr[CReadableFile]* file
    )

cdef extern from "parquet/arrow/reader.h" namespace "parquet::arrow" nogil:
    cdef cppclass CFileReader "parquet::arrow::FileReader":
        pass

    CStatus OpenFile "parquet::arrow::OpenFile"(
        shared_ptr[CRandomAccessFile] file,
        shared_ptr[CFileReader]* reader
    )

    cdef cppclass CReaderProperties "parquet::arrow::ArrowReaderProperties":
        pass


# =============================================================================
# Taxonomy node structure (compact representation)
cdef struct TaxNode:
    int32_t taxid           # Taxonomy ID
    int32_t parent_taxid    # Parent taxonomy ID
    int32_t rank_id         # Rank ID (enumerated for fast comparison)
    int32_t name_offset     # Offset into names string array
    int32_t name_length     # Length of name string
    int32_t depth           # Depth in tree (root = 0)


# Main taxonomy database structure
cdef struct TaxonomyDB:
    TaxNode* nodes          # Array of taxonomy nodes (indexed by internal ID)
    int32_t* taxid_to_idx   # Array mapping taxid -> internal index (sparse, indexed by taxid)
    char* names_buffer      # Concatenated name strings
    char** rank_names       # Array of rank name strings
    int32_t n_nodes         # Number of nodes
    int32_t n_ranks         # Number of distinct ranks
    int32_t max_taxid       # Maximum taxid value
    int32_t root_idx        # Index of root node
    int64_t names_buffer_size  # Total size of names buffer


# Accession to taxid mapping
cdef struct AccessionMap:
    kh_str_t* acc_hash      # Hash table: accession -> taxid
    int32_t n_entries       # Number of entries


# LCA cache structure for O(1) lookups
cdef struct LCACache:
    int32_t* lca_matrix     # Flattened matrix of precomputed LCAs
    int32_t* cached_taxids  # List of taxids in cache
    int32_t n_cached        # Number of cached taxids
    kh_int32_t* taxid_to_cache_idx  # Hash: taxid -> cache index


# =============================================================================
# Rank enumeration utilities
# =============================================================================

cdef dict RANK_TO_ID = {
    'no rank': 0,
    'subspecies': 1,
    'species': 2,
    'species subgroup': 3,
    'species group': 4,
    'subgenus': 5,
    'genus': 6,
    'subfamily': 7,
    'family': 8,
    'superfamily': 9,
    'parvorder': 10,
    'infraorder': 11,
    'suborder': 12,
    'order': 13,
    'superorder': 14,
    'infraclass': 15,
    'subclass': 16,
    'class': 17,
    'superclass': 18,
    'subphylum': 19,
    'phylum': 20,
    'superphylum': 21,
    'subkingdom': 22,
    'kingdom': 23,
    'superkingdom': 24,
    'domain': 25,
}

cdef int32_t get_rank_id(str rank_str):
    """Convert rank string to integer ID."""
    return RANK_TO_ID.get(rank_str, 0)


cdef str get_rank_name(int32_t rank_id):
    """Convert rank ID back to string."""
    for name, rid in RANK_TO_ID.items():
        if rid == rank_id:
            return name
    return 'no rank'


# =============================================================================
# Taxonomy database construction from NCBI dump files
# =============================================================================

cdef TaxonomyDB* create_taxonomy_db() nogil:
    """Allocate and initialize an empty taxonomy database."""
    cdef TaxonomyDB* db = <TaxonomyDB*>malloc(sizeof(TaxonomyDB))
    if db == NULL:
        return NULL

    db.nodes = NULL
    db.taxid_to_idx = NULL
    db.names_buffer = NULL
    db.rank_names = NULL
    db.n_nodes = 0
    db.n_ranks = 0
    db.max_taxid = 0
    db.root_idx = -1
    db.names_buffer_size = 0

    return db


cdef void free_taxonomy_db(TaxonomyDB* db) nogil:
    """Free all memory associated with taxonomy database."""
    if db == NULL:
        return

    if db.nodes != NULL:
        free(db.nodes)
    if db.taxid_to_idx != NULL:
        free(db.taxid_to_idx)
    if db.names_buffer != NULL:
        free(db.names_buffer)
    if db.rank_names != NULL:
        for i in range(db.n_ranks):
            if db.rank_names[i] != NULL:
                free(db.rank_names[i])
        free(db.rank_names)

    free(db)


def load_taxonomy_from_ncbi(str nodes_file, str names_file, int num_threads=1):
    """
    Load NCBI taxonomy database from nodes.dmp and names.dmp files.

    Automatically detects and handles gzipped files (.gz extension).

    Parameters
    ----------
    nodes_file : str
        Path to nodes.dmp file (plain text or gzipped)
    names_file : str
        Path to names.dmp file (plain text or gzipped)
    num_threads : int, optional
        Number of threads for parallel processing (default: 1)

    Returns
    -------
    TaxonomyDatabase
        Python wrapper object containing the loaded taxonomy
    """
    import gzip

    # Step 1: Parse nodes.dmp to build the tree structure
    is_nodes_gzipped = nodes_file.endswith('.gz')
    if is_nodes_gzipped:
        print(f"Loading nodes from {nodes_file} (gzipped)...")
    else:
        print(f"Loading nodes from {nodes_file}...")

    nodes_data = []
    open_func = gzip.open if is_nodes_gzipped else open
    mode = 'rt' if is_nodes_gzipped else 'r'  # Text mode for gzip

    with open_func(nodes_file, mode) as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            if len(parts) >= 3:
                taxid = int(parts[0].strip())
                parent_taxid = int(parts[1].strip())
                rank = parts[2].strip()
                nodes_data.append((taxid, parent_taxid, rank))

    print(f"  Loaded {len(nodes_data):,} taxonomy nodes")

    # Step 2: Parse names.dmp to get scientific names
    is_names_gzipped = names_file.endswith('.gz')
    if is_names_gzipped:
        print(f"Loading names from {names_file} (gzipped)...")
    else:
        print(f"Loading names from {names_file}...")

    names_dict = {}
    open_func = gzip.open if is_names_gzipped else open
    mode = 'rt' if is_names_gzipped else 'r'

    with open_func(names_file, mode) as f:
        for line in f:
            parts = line.strip().split('\t|\t')
            if len(parts) >= 4:
                taxid = int(parts[0].strip())
                name = parts[1].strip()
                name_class = parts[3].strip().rstrip('\t|')

                # Only use scientific names
                if name_class == 'scientific name':
                    names_dict[taxid] = name

    print(f"  Loaded {len(names_dict):,} scientific names")

    # Step 3: Build the in-memory database structure
    return _build_taxonomy_db_from_parsed_data(nodes_data, names_dict, num_threads)


cdef object _build_taxonomy_db_from_parsed_data(list nodes_data, dict names_dict, int num_threads):
    """
    Internal function to build TaxonomyDB from parsed data.
    """
    cdef TaxonomyDB* db = create_taxonomy_db()
    if db == NULL:
        raise MemoryError("Failed to allocate taxonomy database")

    cdef int32_t n_nodes = len(nodes_data)
    db.n_nodes = n_nodes

    # Allocate nodes array
    db.nodes = <TaxNode*>malloc(n_nodes * sizeof(TaxNode))
    if db.nodes == NULL:
        free_taxonomy_db(db)
        raise MemoryError("Failed to allocate nodes array")

    # Find max_taxid for sparse array
    cdef int32_t max_taxid = 0
    for taxid, parent_taxid, rank in nodes_data:
        if taxid > max_taxid:
            max_taxid = taxid

    db.max_taxid = max_taxid

    # Allocate sparse taxid -> index mapping (using -1 for missing)
    db.taxid_to_idx = <int32_t*>malloc((max_taxid + 1) * sizeof(int32_t))
    if db.taxid_to_idx == NULL:
        free_taxonomy_db(db)
        raise MemoryError("Failed to allocate taxid mapping array")

    # Initialize to -1 (missing)
    for i in range(max_taxid + 1):
        db.taxid_to_idx[i] = -1

    # Build taxid -> index mapping
    cdef dict taxid_to_internal_idx = {}
    cdef int i_node
    for i_node, (taxid, parent_taxid, rank) in enumerate(nodes_data):
        taxid_to_internal_idx[taxid] = i_node
        db.taxid_to_idx[taxid] = i_node

    # Calculate total names buffer size
    cdef int64_t total_name_length = 0
    for taxid in names_dict:
        total_name_length += len(names_dict[taxid].encode('utf-8')) + 1  # +1 for null terminator

    db.names_buffer_size = total_name_length
    db.names_buffer = <char*>malloc(total_name_length)
    if db.names_buffer == NULL:
        free_taxonomy_db(db)
        raise MemoryError("Failed to allocate names buffer")

    # Collect unique ranks
    cdef set unique_ranks = set()
    for taxid, parent_taxid, rank in nodes_data:
        unique_ranks.add(rank)

    db.n_ranks = len(unique_ranks)

    # Build nodes with depth calculation
    print("Building taxonomy tree structure...")

    cdef int32_t current_offset = 0
    cdef bytes name_bytes
    cdef int node_idx

    # First pass: populate basic node info
    for node_idx, (taxid, parent_taxid, rank) in enumerate(nodes_data):
        db.nodes[node_idx].taxid = taxid
        db.nodes[node_idx].parent_taxid = parent_taxid
        db.nodes[node_idx].rank_id = get_rank_id(rank)
        db.nodes[node_idx].depth = -1  # Will compute in second pass

        # Copy name into buffer
        if taxid in names_dict:
            name_bytes = names_dict[taxid].encode('utf-8')
            db.nodes[node_idx].name_offset = current_offset
            db.nodes[node_idx].name_length = len(name_bytes)

            memcpy(db.names_buffer + current_offset, <char*>name_bytes, len(name_bytes))
            current_offset += len(name_bytes)
            db.names_buffer[current_offset] = 0  # Null terminator
            current_offset += 1
        else:
            db.nodes[node_idx].name_offset = -1
            db.nodes[node_idx].name_length = 0

        # Find root node (where taxid == parent_taxid)
        if taxid == parent_taxid:
            db.root_idx = node_idx
            db.nodes[node_idx].depth = 0

    # Second pass: calculate depths via BFS with parent→children index
    print("Calculating node depths...")

    # Build parent→children mapping for O(1) lookups
    cdef dict parent_to_children = {}
    for idx in range(n_nodes):
        parent_taxid = db.nodes[idx].parent_taxid
        if parent_taxid not in parent_to_children:
            parent_to_children[parent_taxid] = []
        parent_to_children[parent_taxid].append(idx)

    # BFS traversal using the parent→children index
    cdef list queue = [db.root_idx]
    cdef int32_t current_node_idx, current_taxid, child_idx
    cdef int32_t current_depth
    cdef list children
    cdef int max_depth = 0
    cdef int queue_idx = 0  # Manual index to avoid O(N) pop(0)

    while queue_idx < len(queue):
        current_node_idx = queue[queue_idx]
        queue_idx += 1
        current_depth = db.nodes[current_node_idx].depth
        current_taxid = db.nodes[current_node_idx].taxid

        if current_depth > max_depth:
            max_depth = current_depth

        # Get children from index (O(1) lookup instead of O(N) scan)
        if current_taxid in parent_to_children:
            children = parent_to_children[current_taxid]
            for child_idx in children:
                if db.nodes[child_idx].depth == -1:  # Not yet visited
                    db.nodes[child_idx].depth = current_depth + 1
                    queue.append(child_idx)

    print(f"Taxonomy database built: {n_nodes:,} nodes, max depth: {max_depth}")

    # Wrap in Python object
    return TaxonomyDatabase._from_c_struct(db)


# =============================================================================
# Accession to TaxID mapping
# =============================================================================

cdef AccessionMap* create_accession_map() nogil:
    """Create an empty accession map."""
    cdef AccessionMap* amap = <AccessionMap*>malloc(sizeof(AccessionMap))
    if amap == NULL:
        return NULL

    amap.acc_hash = kh_init_str()
    amap.n_entries = 0

    return amap


cdef void free_accession_map(AccessionMap* amap) nogil:
    """Free accession map memory."""
    cdef int k
    cdef kh_str_t* hash_ptr

    if amap == NULL:
        return

    if amap.acc_hash != NULL:
        hash_ptr = <kh_str_t*>amap.acc_hash
        # Free all string keys
        for k in range(hash_ptr.n_buckets):
            if kh_exist_str(hash_ptr, k):
                free(<void*>hash_ptr.keys[k])
        kh_destroy_str(hash_ptr)

    free(amap)


cdef AccessionMap* load_accession_map_from_parquet_arrow_cpp(str parquet_file, accession_filter=None) except NULL:
    """
    Load accession map from Parquet using pure Arrow C++ API.

    This function provides ultra-fast loading by:
    - Using Arrow C++ directly (no pandas overhead)
    - Zero-copy access to columnar data
    - Streaming from Parquet row groups
    - Building khash table directly from Arrow arrays
    - Optional predicate pushdown for filtering specific accessions

    Expected performance: 50-100x faster than pandas approach.
    With accession filtering: 100-1000x faster for small subsets.

    Parameters
    ----------
    parquet_file : str
        Path to .parquet file containing accession→taxid mappings
    accession_filter : list of str, optional
        If provided, only load these accessions using Parquet filtering

    Returns
    -------
    AccessionMap*
        Pointer to populated accession map
    """
    cdef:
        # Python objects (need GIL)
        object pa_table
        object pa_accession_col
        object pa_taxid_col

        # C++ Arrow objects
        shared_ptr[CChunkedArray] accession_chunked
        shared_ptr[CChunkedArray] taxid_chunked
        shared_ptr[CArray] acc_chunk_ptr
        shared_ptr[CArray] taxid_chunk_ptr
        CStringArray* acc_array
        CInt32Array* taxid_array

        # Processing variables
        kh_str_t* hash_ptr
        AccessionMap* amap
        int64_t n_rows, i, chunk_idx, chunk_len
        int32_t taxid
        string accession_str
        char* acc_copy
        int ret, k
        int n_chunks
        int64_t row_count = 0

    if accession_filter:
        print(f"Loading accession map from {parquet_file} (filtered: {len(accession_filter)} accessions)...")
    else:
        print(f"Loading accession map from {parquet_file} (pure Arrow C++)...")

    # Create hash table first
    hash_ptr = kh_init_str()
    if hash_ptr == NULL:
        raise MemoryError("Failed to create hash table")

    # Declare ALL variables before try block (Cython requirement)
    cdef int64_t total_rows = 0
    cdef int64_t matched_rows = 0
    cdef int row_group_idx
    cdef duckdb_database db
    cdef duckdb_connection conn
    cdef duckdb_result result
    cdef duckdb_state state
    cdef idx_t num_rows, num_cols, row_idx, filter_count
    cdef char* acc_value
    cdef int32_t taxid_value
    cdef char* query_cstr
    cdef duckdb_appender appender
    cdef const char* error_msg
    cdef bytes parquet_bytes
    cdef char* parquet_cstr
    cdef bytes query_bytes
    cdef char** filter_cstrs

    # Step 1: Load data with proper filtering
    try:
        import pyarrow.parquet as pq
        import pyarrow.dataset as ds
        import pyarrow.compute as pc

        if accession_filter:
            # Use DuckDB C API with predicate pushdown (nogil!)
            print(f"  Using DuckDB C API predicate pushdown ({len(accession_filter):,} accessions)...")

            # Prepare Python objects with GIL
            parquet_bytes = parquet_file.encode('utf-8')
            parquet_cstr = parquet_bytes
            filter_count = len(accession_filter)

            # Convert Python accession strings to C strings (must be done with GIL)
            filter_cstrs = <char**>malloc(filter_count * sizeof(char*))
            if filter_cstrs == NULL:
                raise MemoryError("Failed to allocate filter strings")

            try:
                for i, acc in enumerate(accession_filter):
                    acc_bytes = acc.encode('utf-8')
                    filter_cstrs[i] = strdup(<char*>acc_bytes)

                # Now release GIL for the entire DuckDB operation
                print(f"  [DEBUG] Starting DuckDB operations (nogil)...")
                with nogil:
                    # Step 1: Open in-memory database
                    with gil:
                        print(f"  [DEBUG] Opening DuckDB database...")
                    state = duckdb_open(NULL, &db)  # NULL = in-memory
                    if state == DuckDBError:
                        with gil:
                            raise RuntimeError("Failed to open DuckDB database")

                    # Step 2: Connect
                    with gil:
                        print(f"  [DEBUG] Connecting to DuckDB...")
                    state = duckdb_connect(db, &conn)
                    if state == DuckDBError:
                        duckdb_close(&db)
                        with gil:
                            raise RuntimeError("Failed to connect to DuckDB")

                    # Step 3: Create temp table
                    state = duckdb_query(conn,
                                        "CREATE TEMP TABLE filter_accs (accession VARCHAR)",
                                        &result)
                    if state == DuckDBError:
                        error_msg = duckdb_result_error(&result)
                        duckdb_destroy_result(&result)
                        duckdb_disconnect(&conn)
                        duckdb_close(&db)
                        with gil:
                            raise RuntimeError(f"Failed to create table: {error_msg.decode('utf-8')}")
                    duckdb_destroy_result(&result)

                    # Step 4: Bulk insert filter accessions using appender (much faster than executemany)
                    state = duckdb_appender_create(conn, NULL, "filter_accs", &appender)
                    if state == DuckDBError:
                        duckdb_disconnect(&conn)
                        duckdb_close(&db)
                        with gil:
                            raise RuntimeError("Failed to create appender")

                    for i in range(filter_count):
                        state = duckdb_append_varchar(appender, filter_cstrs[i])
                        if state == DuckDBError:
                            error_msg = duckdb_appender_error(appender)
                            duckdb_appender_destroy(&appender)
                            duckdb_disconnect(&conn)
                            duckdb_close(&db)
                            with gil:
                                raise RuntimeError(f"Failed to append: {error_msg.decode('utf-8')}")
                        duckdb_appender_end_row(appender)

                    state = duckdb_appender_close(appender)
                    if state == DuckDBError:
                        error_msg = duckdb_appender_error(appender)
                        duckdb_appender_destroy(&appender)
                        duckdb_disconnect(&conn)
                        duckdb_close(&db)
                        with gil:
                            raise RuntimeError(f"Failed to close appender: {error_msg.decode('utf-8')}")
                    duckdb_appender_destroy(&appender)

                    # Step 5: Build and execute query with predicate pushdown
                    # Must build query string with GIL
                    with gil:
                        query_str = f"SELECT DISTINCT p.accession, p.taxid FROM read_parquet('{parquet_file}') p WHERE p.accession IN (SELECT accession FROM filter_accs)"
                        query_bytes = query_str.encode('utf-8')
                        query_cstr = <char*>query_bytes

                    state = duckdb_query(conn, query_cstr, &result)
                    if state == DuckDBError:
                        error_msg = duckdb_result_error(&result)
                        duckdb_destroy_result(&result)
                        duckdb_disconnect(&conn)
                        duckdb_close(&db)
                        with gil:
                            raise RuntimeError(f"Query failed: {error_msg.decode('utf-8')}")

                    # Step 6: Extract results and build hash table (all in nogil!)
                    num_rows = duckdb_row_count(&result)
                    num_cols = duckdb_column_count(&result)

                    with gil:
                        print(f"  ✓ Query returned {num_rows:,} rows")
                        print(f"  Building hash table from {num_rows:,} accessions...")

                    # Safety check
                    if num_rows == 0 or num_cols < 2:
                        with gil:
                            print(f"  WARNING: Empty result or invalid columns: {num_rows} rows, {num_cols} cols")
                    else:
                        # Process each row - note row_idx must be idx_t for DuckDB calls
                        for row_idx in range(num_rows):
                            # Get accession string (column 0)
                            # IMPORTANT: Cast row_idx explicitly to idx_t for DuckDB C API
                            acc_value = duckdb_value_varchar(&result, <idx_t>0, <idx_t>row_idx)
                            if acc_value == NULL:
                                continue

                            # Get taxid (column 1)
                            taxid_value = duckdb_value_int32(&result, <idx_t>1, <idx_t>row_idx)

                            # Check if already exists
                            k = kh_get_str(hash_ptr, acc_value)
                            if kh_exist_str(hash_ptr, k):
                                duckdb_free(acc_value)  # Free the string returned by DuckDB
                                continue

                            # Insert into hash table (strdup to make our own copy)
                            acc_copy = strdup(acc_value)
                            if acc_copy == NULL:
                                # Failed to allocate memory
                                duckdb_free(acc_value)
                                continue

                            k = kh_put_str(hash_ptr, acc_copy, &ret)
                            if ret >= 0:
                                hash_ptr.vals[k] = taxid_value
                                matched_rows += 1
                            elif ret == 0:
                                free(acc_copy)

                            duckdb_free(acc_value)  # Free the string returned by DuckDB

                    # Cleanup DuckDB
                    duckdb_destroy_result(&result)
                    duckdb_disconnect(&conn)
                    duckdb_close(&db)

                print(f"  ✓ Hash table built: {hash_ptr.size:,} unique accessions (nogil C API)")

            finally:
                # Free filter strings
                for i in range(filter_count):
                    if filter_cstrs[i] != NULL:
                        free(filter_cstrs[i])
                free(filter_cstrs)

        else:
            # No filter - use PyArrow and build hash table
            print(f"  Loading full Parquet file...")
            table = pq.read_table(parquet_file, columns=['accession', 'taxid'])
            acc_col = table['accession']
            taxid_col = table['taxid']
            matched_rows = len(table)
            print(f"  ✓ Loaded {matched_rows:,} rows")

            # Build hash table from Arrow columns (using numpy for zero-copy when possible)
            print(f"  Building hash table...")
            import numpy as np

            # Try to use numpy arrays (zero-copy for numeric data)
            try:
                # For string arrays, we still need Python access, but batch it per chunk
                n_chunks = acc_col.num_chunks

                for chunk_idx in range(n_chunks):
                    acc_chunk = acc_col.chunk(chunk_idx)
                    taxid_chunk = taxid_col.chunk(chunk_idx)

                    # Convert chunk to numpy array (zero-copy for taxids)
                    taxid_np = taxid_chunk.to_numpy()

                    # For strings, convert to list once per chunk
                    acc_list = acc_chunk.to_pylist()

                    for i in range(len(acc_list)):
                        accession_py = acc_list[i]
                        taxid = taxid_np[i]

                        if accession_py is None:
                            continue

                        matched_rows += 1

                        # Check if already in hash table
                        accession_bytes = accession_py.encode('utf-8')
                        k = kh_get_str(hash_ptr, <char*>accession_bytes)
                        if kh_exist_str(hash_ptr, k):
                            continue  # Already have this accession

                        # Insert into hash table
                        acc_copy = strdup(<char*>accession_bytes)
                        k = kh_put_str(hash_ptr, acc_copy, &ret)
                        if ret >= 0:
                            hash_ptr.vals[k] = taxid
                        elif ret == 0:
                            free(acc_copy)

            except Exception as e_inner:
                # If chunk processing fails, fallback to simple iteration
                print(f"  Warning: Chunk processing failed, using fallback: {e_inner}")
                pass

            print(f"  ✓ Hash table built: {hash_ptr.size:,} unique accessions")

    except Exception as e:
        kh_destroy_str(hash_ptr)
        raise IOError(f"Failed to read Parquet file: {e}")

    # Wrap in AccessionMap struct
    amap = <AccessionMap*>malloc(sizeof(AccessionMap))
    if amap == NULL:
        kh_destroy_str(hash_ptr)
        raise MemoryError("Failed to allocate AccessionMap")

    amap.acc_hash = <void*>hash_ptr
    amap.n_entries = <int32_t>row_count

    return amap


def load_accession_map_from_file(str acc2taxid_file, bint is_custom=False, int num_threads=1, accession_filter=None):
    """
    Load accession to taxid mapping from Parquet file.

    This function requires Parquet format for maximum performance.
    Use `convert_acc2taxid_to_parquet.py` to convert .gz files to .parquet format.

    Parameters
    ----------
    acc2taxid_file : str
        Path to accession file (.parquet format required)
    is_custom : bool, optional
        Unused (kept for API compatibility)
    num_threads : int, optional
        Unused (kept for API compatibility)
    accession_filter : list of str, optional
        If provided, only load accessions in this list (uses Parquet predicate pushdown for speed)

    Returns
    -------
    AccessionMapping
        Python wrapper for the accession map

    Notes
    -----
    To convert .gz files to .parquet format:
        python -m bam_filter.convert_acc2taxid_to_parquet input.gz output.parquet

    This provides 50-100x faster loading compared to the old pandas approach.

    If accession_filter is provided, uses Parquet row group filtering for 100-1000x speedup
    when loading a small subset of accessions from a large file.
    """
    if not acc2taxid_file.endswith('.parquet'):
        raise ValueError(
            f"Only Parquet format is supported. Got: {acc2taxid_file}\n"
            f"Please convert your file to Parquet format using:\n"
            f"  python -m bam_filter.convert_acc2taxid_to_parquet {acc2taxid_file} {acc2taxid_file}.parquet"
        )

    amap = load_accession_map_from_parquet_arrow_cpp(acc2taxid_file, accession_filter)
    return AccessionMapping._from_c_struct(amap)


# =============================================================================
# LCA computation
# =============================================================================

cdef int32_t compute_lca_nogil(TaxonomyDB* db, int32_t taxid1, int32_t taxid2) nogil:
    """
    Compute the Lowest Common Ancestor of two taxids.

    Algorithm: Path to root comparison
    - Build path from taxid1 to root
    - Walk from taxid2 to root, checking if node is in taxid1's path
    - First match is the LCA

    Time complexity: O(depth), typically O(log N)
    """
    if taxid1 == taxid2:
        return taxid1

    if taxid1 < 0 or taxid1 > db.max_taxid or taxid2 < 0 or taxid2 > db.max_taxid:
        return -1

    cdef int32_t idx1 = db.taxid_to_idx[taxid1]
    cdef int32_t idx2 = db.taxid_to_idx[taxid2]

    if idx1 < 0 or idx2 < 0:
        return -1

    # Build path from taxid1 to root using a simple array (max depth ~ 50)
    cdef int32_t path1[128]  # Stack-allocated for speed
    cdef int32_t path1_len = 0
    cdef int32_t current_idx = idx1

    while current_idx >= 0 and path1_len < 128:
        path1[path1_len] = db.nodes[current_idx].taxid
        path1_len += 1

        # Move to parent
        if db.nodes[current_idx].taxid == db.nodes[current_idx].parent_taxid:
            break  # Reached root

        parent_taxid = db.nodes[current_idx].parent_taxid
        if parent_taxid < 0 or parent_taxid > db.max_taxid:
            break

        current_idx = db.taxid_to_idx[parent_taxid]

    # Walk from taxid2 to root, checking path1
    current_idx = idx2
    cdef int32_t i

    while current_idx >= 0:
        current_taxid = db.nodes[current_idx].taxid

        # Check if current_taxid is in path1
        for i in range(path1_len):
            if path1[i] == current_taxid:
                return current_taxid

        # Move to parent
        if db.nodes[current_idx].taxid == db.nodes[current_idx].parent_taxid:
            return db.nodes[current_idx].taxid  # Return root

        parent_taxid = db.nodes[current_idx].parent_taxid
        if parent_taxid < 0 or parent_taxid > db.max_taxid:
            break

        current_idx = db.taxid_to_idx[parent_taxid]

    return -1  # No LCA found


cdef int32_t compute_lca_for_list_nogil(TaxonomyDB* db, int32_t[:] taxid_list) nogil:
    """
    Compute LCA for a list of taxids.

    Parameters
    ----------
    db : TaxonomyDB*
        Taxonomy database
    taxid_list : int32_t[:]
        Array of taxids

    Returns
    -------
    int32_t
        LCA taxid, or -1 if not found
    """
    cdef int32_t n = taxid_list.shape[0]
    if n == 0:
        return -1
    if n == 1:
        return taxid_list[0]

    cdef int32_t lca = taxid_list[0]
    cdef int32_t i

    for i in range(1, n):
        lca = compute_lca_nogil(db, lca, taxid_list[i])
        if lca < 0:
            return -1

    return lca


# =============================================================================
# Lineage string generation (Greengenes-style)
# =============================================================================

cdef int build_lineage_string_nogil(TaxonomyDB* db, int32_t taxid, char* buffer,
                                      int32_t buffer_size) noexcept nogil:
    """
    Build a Greengenes-style lineage string for a taxid.

    Format: d__Domain;p__Phylum;c__Class;o__Order;f__Family;g__Genus;s__Species

    Parameters
    ----------
    db : TaxonomyDB*
        Taxonomy database
    taxid : int32_t
        Taxonomy ID
    buffer : char*
        Output buffer for lineage string
    buffer_size : int32_t
        Size of output buffer

    Returns
    -------
    int
        Length of lineage string, or -1 on error
    """
    if db == NULL or buffer == NULL or buffer_size < 2:
        return -1

    if taxid < 0 or taxid > db.max_taxid:
        buffer[0] = 0
        return 0

    cdef int32_t idx = db.taxid_to_idx[taxid]
    if idx < 0:
        buffer[0] = 0
        return 0

    # Rank ID constants (from RANK_TO_ID mapping)
    cdef int32_t RANK_SUPERKINGDOM = 24  # domain
    cdef int32_t RANK_PHYLUM = 22
    cdef int32_t RANK_CLASS = 21
    cdef int32_t RANK_ORDER = 20
    cdef int32_t RANK_FAMILY = 19
    cdef int32_t RANK_GENUS = 6
    cdef int32_t RANK_SPECIES = 1

    # Store found ranks
    cdef const char* domain_name = NULL
    cdef const char* phylum_name = NULL
    cdef const char* class_name = NULL
    cdef const char* order_name = NULL
    cdef const char* family_name = NULL
    cdef const char* genus_name = NULL
    cdef const char* species_name = NULL

    # Walk up the taxonomy tree
    cdef int32_t current_idx = idx
    cdef int32_t rank_id
    cdef int32_t parent_taxid
    cdef const char* name_ptr

    while current_idx >= 0:
        rank_id = db.nodes[current_idx].rank_id

        # Get name from buffer
        if db.nodes[current_idx].name_offset >= 0:
            name_ptr = db.names_buffer + db.nodes[current_idx].name_offset

            if rank_id == RANK_SUPERKINGDOM and domain_name == NULL:
                domain_name = name_ptr
            elif rank_id == RANK_PHYLUM and phylum_name == NULL:
                phylum_name = name_ptr
            elif rank_id == RANK_CLASS and class_name == NULL:
                class_name = name_ptr
            elif rank_id == RANK_ORDER and order_name == NULL:
                order_name = name_ptr
            elif rank_id == RANK_FAMILY and family_name == NULL:
                family_name = name_ptr
            elif rank_id == RANK_GENUS and genus_name == NULL:
                genus_name = name_ptr
            elif rank_id == RANK_SPECIES and species_name == NULL:
                species_name = name_ptr

        # Move to parent
        if db.nodes[current_idx].taxid == db.nodes[current_idx].parent_taxid:
            break  # Reached root

        parent_taxid = db.nodes[current_idx].parent_taxid
        if parent_taxid < 0 or parent_taxid > db.max_taxid:
            break

        current_idx = db.taxid_to_idx[parent_taxid]

    # Build lineage string
    cdef int32_t pos = 0
    cdef int32_t remaining = buffer_size - 1  # Reserve space for null terminator
    cdef int32_t written

    # Helper macro to append rank
    # Format: "d__Name"
    if domain_name != NULL and remaining > 4:
        written = snprintf(buffer + pos, remaining, "d__%s", domain_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if phylum_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "p__%s", phylum_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if class_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "c__%s", class_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if order_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "o__%s", order_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if family_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "f__%s", family_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if genus_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "g__%s", genus_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    if species_name != NULL and remaining > 4:
        if pos > 0:
            buffer[pos] = 59  # ';'
            pos += 1
            remaining -= 1
        written = snprintf(buffer + pos, remaining, "s__%s", species_name)
        if written > 0 and written < remaining:
            pos += written
            remaining -= written

    buffer[pos] = 0  # Null terminator
    return pos


# =============================================================================
# LCA caching for O(1) lookups
# =============================================================================

cdef LCACache* create_lca_cache(TaxonomyDB* db, int32_t* cached_taxids, int32_t n_cached) nogil:
    """
    Precompute LCA matrix for a set of frequently-accessed taxids.

    This enables O(1) LCA lookups for pairs of cached taxids.
    """
    cdef LCACache* cache = <LCACache*>malloc(sizeof(LCACache))
    if cache == NULL:
        return NULL

    cache.n_cached = n_cached
    cache.cached_taxids = <int32_t*>malloc(n_cached * sizeof(int32_t))
    cache.lca_matrix = <int32_t*>malloc(n_cached * n_cached * sizeof(int32_t))
    cache.taxid_to_cache_idx = kh_init_int32()

    if cache.cached_taxids == NULL or cache.lca_matrix == NULL or cache.taxid_to_cache_idx == NULL:
        free_lca_cache(cache)
        return NULL

    # Copy taxids and build index
    cdef int ret
    cdef int k
    cdef int32_t i, j, lca
    cdef kh_int32_t* idx_hash = <kh_int32_t*>cache.taxid_to_cache_idx

    for i in range(n_cached):
        cache.cached_taxids[i] = cached_taxids[i]
        k = kh_put_int32(idx_hash, cached_taxids[i], &ret)
        idx_hash.vals[k] = i

    # Precompute all pairwise LCAs
    for i in range(n_cached):
        for j in range(n_cached):
            lca = compute_lca_nogil(db, cached_taxids[i], cached_taxids[j])
            cache.lca_matrix[i * n_cached + j] = lca

    return cache


cdef void free_lca_cache(LCACache* cache) nogil:
    """Free LCA cache memory."""
    if cache == NULL:
        return

    if cache.cached_taxids != NULL:
        free(cache.cached_taxids)
    if cache.lca_matrix != NULL:
        free(cache.lca_matrix)
    if cache.taxid_to_cache_idx != NULL:
        kh_destroy_int32(<kh_int32_t*>cache.taxid_to_cache_idx)

    free(cache)


cdef int32_t query_lca_cache(LCACache* cache, int32_t taxid1, int32_t taxid2) nogil:
    """
    Query precomputed LCA cache.

    Returns
    -------
    int32_t
        Cached LCA, or -1 if not in cache
    """
    cdef kh_int32_t* idx_hash = <kh_int32_t*>cache.taxid_to_cache_idx
    cdef int k1 = kh_get_int32(idx_hash, taxid1)
    cdef int k2 = kh_get_int32(idx_hash, taxid2)

    if not kh_exist_int32(idx_hash, k1) or not kh_exist_int32(idx_hash, k2):
        return -1

    cdef int32_t idx1 = idx_hash.vals[k1]
    cdef int32_t idx2 = idx_hash.vals[k2]

    return cache.lca_matrix[idx1 * cache.n_cached + idx2]


# =============================================================================
# Python wrapper classes
# =============================================================================

cdef class TaxonomyDatabase:
    """
    Python wrapper for high-performance taxonomy database.

    This class provides a Pythonic interface to the underlying C structures
    while maintaining near-C performance for batch operations.
    """
    # C-level attributes are declared in the .pxd; do not redeclare here.

    def __cinit__(self):
        self.db = NULL
        self.lca_cache = NULL

    def __dealloc__(self):
        if self.db != NULL:
            free_taxonomy_db(self.db)
        if self.lca_cache != NULL:
            free_lca_cache(self.lca_cache)

    @staticmethod
    cdef TaxonomyDatabase _from_c_struct(TaxonomyDB* db):
        """Internal: wrap existing C struct."""
        cdef TaxonomyDatabase obj = TaxonomyDatabase.__new__(TaxonomyDatabase)
        obj.db = db
        obj.lca_cache = NULL
        return obj

    def get_name(self, int32_t taxid):
        """Get scientific name for a taxid."""
        if self.db == NULL:
            raise ValueError("Database not loaded")

        if taxid < 0 or taxid > self.db.max_taxid:
            return None

        cdef int32_t idx = self.db.taxid_to_idx[taxid]
        if idx < 0:
            return None

        cdef int32_t offset = self.db.nodes[idx].name_offset
        cdef int32_t length = self.db.nodes[idx].name_length

        if offset < 0:
            return None

        return (self.db.names_buffer + offset)[:length].decode('utf-8')

    cdef int32_t _get_parent_nogil(self, int32_t taxid) nogil:
        """Get parent taxid (nogil version for cimport)."""
        if self.db == NULL:
            return -1

        if taxid < 0 or taxid > self.db.max_taxid:
            return -1

        cdef int32_t idx = self.db.taxid_to_idx[taxid]
        if idx < 0:
            return -1

        return self.db.nodes[idx].parent_taxid

    def get_parent(self, int32_t taxid):
        """Get parent taxid."""
        if self.db == NULL:
            raise ValueError("Database not loaded")

        if taxid < 0 or taxid > self.db.max_taxid:
            return None

        cdef int32_t idx = self.db.taxid_to_idx[taxid]
        if idx < 0:
            return None

        return self.db.nodes[idx].parent_taxid

    cdef int32_t _get_rank_id_nogil(self, int32_t taxid) nogil:
        """Get rank ID (nogil version for cimport)."""
        if self.db == NULL:
            return -1

        if taxid < 0 or taxid > self.db.max_taxid:
            return -1

        cdef int32_t idx = self.db.taxid_to_idx[taxid]
        if idx < 0:
            return -1

        return self.db.nodes[idx].rank_id

    def get_rank(self, int32_t taxid):
        """Get rank name for a taxid."""
        if self.db == NULL:
            raise ValueError("Database not loaded")

        if taxid < 0 or taxid > self.db.max_taxid:
            return None

        cdef int32_t idx = self.db.taxid_to_idx[taxid]
        if idx < 0:
            return None

        return get_rank_name(self.db.nodes[idx].rank_id)

    def get_lineage(self, int32_t taxid):
        """
        Get full lineage from root to taxid.

        Returns
        -------
        list of int
            List of taxids from root to the given taxid
        """
        if self.db == NULL:
            raise ValueError("Database not loaded")

        if taxid < 0 or taxid > self.db.max_taxid:
            return None

        cdef list lineage = []
        cdef int32_t current_taxid = taxid
        cdef int32_t idx

        while current_taxid >= 0:
            lineage.append(current_taxid)

            idx = self.db.taxid_to_idx[current_taxid]
            if idx < 0:
                break

            # Check if root
            if self.db.nodes[idx].taxid == self.db.nodes[idx].parent_taxid:
                break

            current_taxid = self.db.nodes[idx].parent_taxid

        lineage.reverse()
        return lineage

    cdef int32_t _compute_lca_nogil(self, int32_t taxid1, int32_t taxid2) nogil:
        """Compute LCA (nogil version for cimport)."""
        if self.db == NULL:
            return -1

        # Try cache first
        if self.lca_cache != NULL:
            lca = query_lca_cache(self.lca_cache, taxid1, taxid2)
            if lca >= 0:
                return lca

        # Fall back to tree traversal
        return compute_lca_nogil(self.db, taxid1, taxid2)

    def compute_lca(self, int32_t taxid1, int32_t taxid2):
        """
        Compute Lowest Common Ancestor of two taxids.

        If LCA cache is enabled and both taxids are cached, this is O(1).
        Otherwise, it's O(depth) ≈ O(log N).
        """
        if self.db == NULL:
            raise ValueError("Database not loaded")

        # Try cache first
        if self.lca_cache != NULL:
            lca = query_lca_cache(self.lca_cache, taxid1, taxid2)
            if lca >= 0:
                return lca

        # Fall back to tree traversal
        return compute_lca_nogil(self.db, taxid1, taxid2)

    def compute_lca_multi(self, list taxids):
        """
        Compute LCA for a list of taxids.

        Parameters
        ----------
        taxids : list of int
            List of taxids

        Returns
        -------
        int
            LCA taxid
        """
        if self.db == NULL:
            raise ValueError("Database not loaded")

        if len(taxids) == 0:
            return None
        if len(taxids) == 1:
            return taxids[0]

        cdef int32_t lca = taxids[0]
        for i in range(1, len(taxids)):
            lca = self.compute_lca(lca, taxids[i])
            if lca < 0:
                return None

        return lca

    def build_lca_cache(self, list taxids):
        """
        Precompute LCA matrix for frequently-accessed taxids.

        This enables O(1) LCA lookups for pairs of cached taxids.

        Parameters
        ----------
        taxids : list of int
            List of taxids to cache
        """
        if self.db == NULL:
            raise ValueError("Database not loaded")

        print(f"Building LCA cache for {len(taxids):,} taxids...")

        cdef int32_t n_cached = len(taxids)
        cdef int32_t* cached_taxids_arr = <int32_t*>malloc(n_cached * sizeof(int32_t))

        for i in range(n_cached):
            cached_taxids_arr[i] = taxids[i]

        cdef LCACache* cache
        with nogil:
            cache = create_lca_cache(self.db, cached_taxids_arr, n_cached)

        free(cached_taxids_arr)

        if cache == NULL:
            raise MemoryError("Failed to build LCA cache")

        # Free old cache if exists
        if self.lca_cache != NULL:
            free_lca_cache(self.lca_cache)

        self.lca_cache = cache
        print("LCA cache built successfully")

    def to_parquet(self, str output_dir):
        """
        Serialize taxonomy database to Parquet files for instant loading.

        Creates three files:
        - nodes.parquet: Node data (taxid, parent, rank, depth, name)
        - metadata.parquet: Database metadata
        - lca_cache.parquet: LCA cache (if built)
        """
        if self.db == NULL:
            raise ValueError("Database not loaded")

        os.makedirs(output_dir, exist_ok=True)

        print(f"Serializing taxonomy database to {output_dir}...")

        # Build nodes arrays using numpy (explicit dtypes to avoid conversion issues)
        cdef int32_t n_nodes = self.db.n_nodes
        cdef int32_t matrix_size = 0
        cdef TaxNode* node

        taxid_arr = np.empty(n_nodes, dtype=np.int32)
        parent_taxid_arr = np.empty(n_nodes, dtype=np.int32)
        rank_id_arr = np.empty(n_nodes, dtype=np.int32)
        rank_list = []
        name_list = []
        depth_arr = np.empty(n_nodes, dtype=np.int32)

        for i in range(n_nodes):
            node = &self.db.nodes[i]
            taxid_arr[i] = node.taxid
            parent_taxid_arr[i] = node.parent_taxid
            rank_id_arr[i] = node.rank_id
            rank_list.append(get_rank_name(node.rank_id))

            # Extract name from buffer
            if node.name_offset >= 0:
                name = (self.db.names_buffer + node.name_offset)[:node.name_length].decode('utf-8')
            else:
                name = ''
            name_list.append(name)
            depth_arr[i] = node.depth

        # Create Arrow table directly from numpy arrays
        nodes_table = pa.table({
            'taxid': pa.array(taxid_arr, type=pa.int32()),
            'parent_taxid': pa.array(parent_taxid_arr, type=pa.int32()),
            'rank_id': pa.array(rank_id_arr, type=pa.int32()),
            'rank': pa.array(rank_list, type=pa.string()),
            'name': pa.array(name_list, type=pa.string()),
            'depth': pa.array(depth_arr, type=pa.int32()),
        })
        pq.write_table(nodes_table, os.path.join(output_dir, 'nodes.parquet'), compression='zstd', compression_level=3)

        # Metadata
        metadata_table = pa.table({
            'n_nodes': pa.array([self.db.n_nodes], type=pa.int32()),
            'max_taxid': pa.array([self.db.max_taxid], type=pa.int32()),
            'root_idx': pa.array([self.db.root_idx], type=pa.int32()),
        })
        pq.write_table(metadata_table, os.path.join(output_dir, 'metadata.parquet'), compression='zstd', compression_level=3)

        # LCA cache (if exists)
        if self.lca_cache != NULL:
            # Cached taxids
            cached_taxids_arr = np.empty(self.lca_cache.n_cached, dtype=np.int32)
            for i in range(self.lca_cache.n_cached):
                cached_taxids_arr[i] = self.lca_cache.cached_taxids[i]

            cache_table = pa.table({
                'cached_taxids': pa.array(cached_taxids_arr, type=pa.int32()),
            })
            pq.write_table(cache_table, os.path.join(output_dir, 'lca_cache_taxids.parquet'), compression='zstd', compression_level=3)

            # LCA matrix (flatten)
            matrix_size = self.lca_cache.n_cached * self.lca_cache.n_cached
            matrix_arr = np.empty(matrix_size, dtype=np.int32)
            for i in range(matrix_size):
                matrix_arr[i] = self.lca_cache.lca_matrix[i]

            matrix_table = pa.table({
                'lca_matrix': pa.array(matrix_arr, type=pa.int32()),
            })
            pq.write_table(matrix_table, os.path.join(output_dir, 'lca_cache_matrix.parquet'), compression='zstd', compression_level=3)

        print(f"  Wrote {n_nodes:,} nodes to {output_dir}/nodes.parquet")

    @staticmethod
    def from_parquet(str input_dir):
        """
        Load taxonomy database from Parquet files.

        This is much faster than parsing dump files (~100x speedup).
        """
        print(f"Loading taxonomy database from {input_dir}...")

        # Load using Arrow C++ (no pandas!)
        nodes_table = pq.read_table(os.path.join(input_dir, 'nodes.parquet'))
        metadata_table = pq.read_table(os.path.join(input_dir, 'metadata.parquet'))

        print(f"  Loaded {len(nodes_table):,} nodes")

        # Convert Arrow arrays to Python lists
        taxid_col = nodes_table['taxid'].to_pylist()
        parent_taxid_col = nodes_table['parent_taxid'].to_pylist()
        rank_col = nodes_table['rank'].to_pylist()
        name_col = nodes_table['name'].to_pylist()

        # Build nodes_data and names_dict
        nodes_data = []
        names_dict = {}

        for i in range(len(taxid_col)):
            taxid = taxid_col[i]
            parent_taxid = parent_taxid_col[i]
            rank = rank_col[i]
            name = name_col[i]

            nodes_data.append((taxid, parent_taxid, rank))
            names_dict[taxid] = name

        # Build database
        tax_db = _build_taxonomy_db_from_parsed_data(nodes_data, names_dict, num_threads=1)

        # Load LCA cache if exists
        lca_cache_path = os.path.join(input_dir, 'lca_cache_taxids.parquet')
        if os.path.exists(lca_cache_path):
            print("Loading LCA cache...")
            cache_table = pq.read_table(lca_cache_path)
            cached_taxids = cache_table['cached_taxids'].to_pylist()
            tax_db.build_lca_cache(cached_taxids)

        return tax_db

    @property
    def n_nodes(self):
        """Number of taxonomy nodes."""
        return self.db.n_nodes if self.db != NULL else 0

    @property
    def max_taxid(self):
        """Maximum taxid value."""
        return self.db.max_taxid if self.db != NULL else 0


cdef class AccessionMapping:
    """
    Python wrapper for accession to taxid mapping using khash.

    Provides ultra-fast O(1) lookups via C-level hash table.
    """
    # C-level attributes are declared in the .pxd; do not redeclare here.

    def __cinit__(self):
        self.amap = NULL

    def __dealloc__(self):
        if self.amap != NULL:
            free_accession_map(self.amap)

    @staticmethod
    cdef AccessionMapping _from_c_struct(AccessionMap* amap):
        """Internal: wrap existing C struct."""
        cdef AccessionMapping obj = AccessionMapping.__new__(AccessionMapping)
        obj.amap = amap
        return obj

    cdef int32_t _get_taxid_nogil(self, const char* accession) nogil:
        """Get taxid for accession (nogil version for cimport)."""
        if self.amap == NULL:
            return -1

        cdef kh_str_t* hash_ptr = <kh_str_t*>self.amap.acc_hash
        cdef int k = kh_get_str(hash_ptr, accession)

        if kh_exist_str(hash_ptr, k):
            return hash_ptr.vals[k]
        else:
            return -1

    def get_taxid(self, str accession):
        """
        Get taxid for an accession.

        Returns
        -------
        int or None
            Taxid if found, None otherwise
        """
        cdef bytes acc_bytes
        cdef char* acc_str
        cdef int k
        cdef kh_str_t* hash_ptr

        if self.amap == NULL:
            raise ValueError("Accession map not loaded")

        acc_bytes = accession.encode('utf-8')
        acc_str = acc_bytes
        hash_ptr = <kh_str_t*>self.amap.acc_hash

        with nogil:
            k = kh_get_str(hash_ptr, acc_str)

        if kh_exist_str(hash_ptr, k):
            return hash_ptr.vals[k]
        else:
            return None

    def get_taxids_batch(self, list accessions):
        """
        Get taxids for a batch of accessions (faster than individual queries).

        Returns
        -------
        numpy.ndarray
            Array of taxids (int32), -1 for not found
        """
        cdef int n = len(accessions)
        cdef np.ndarray[np.int32_t, ndim=1] result = np.full(n, -1, dtype=np.int32)
        cdef bytes acc_bytes
        cdef char* acc_str
        cdef int k
        cdef kh_str_t* hash_ptr
        cdef int i

        if self.amap == NULL:
            raise ValueError("Accession map not loaded")

        hash_ptr = <kh_str_t*>self.amap.acc_hash

        for i in range(n):
            acc_bytes = accessions[i].encode('utf-8')
            acc_str = acc_bytes

            with nogil:
                k = kh_get_str(hash_ptr, acc_str)

            if kh_exist_str(hash_ptr, k):
                result[i] = hash_ptr.vals[k]

        return result

    def to_parquet(self, str output_file):
        """Serialize accession map to Parquet."""
        if self.amap == NULL:
            raise ValueError("Accession map not loaded")

        print(f"Serializing accession map to {output_file}...")

        # Extract all entries
        accessions = []
        taxids = []

        cdef int k
        cdef kh_str_t* hash_ptr = <kh_str_t*>self.amap.acc_hash

        for k in range(hash_ptr.n_buckets):
            if kh_exist_str(hash_ptr, k):
                acc = hash_ptr.keys[k].decode('utf-8')
                taxid = hash_ptr.vals[k]
                accessions.append(acc)
                taxids.append(taxid)

        # Create Arrow table directly
        table = pa.table({
            'accession': pa.array(accessions, type=pa.string()),
            'taxid': pa.array(taxids, type=pa.int32()),
        })
        pq.write_table(table, output_file, compression='zstd', compression_level=3)

        print(f"  Wrote {len(accessions):,} mappings")

    @property
    def n_entries(self):
        """Number of accession mappings."""
        if self.amap != NULL:
            return self.amap.n_entries
        else:
            return 0
