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

"""Read name hashing and thread-local ID management.

Provides FNV-1a hashing for read names and thread-local hash maps for
assigning sequential read IDs during parallel batch processing.
"""

from libc.stdint cimport int64_t, uint32_t
from libc.stdlib cimport malloc, free
from libc.string cimport strdup


cdef int64_t compute_read_name_hash(const char* read_name) noexcept nogil:
    """Compute FNV-1a hash for read name.

    Parameters
    ----------
    read_name : const char*
        Null-terminated read name string

    Returns
    -------
    int64_t
        Non-negative FNV-1a hash value
    """
    cdef int64_t hash_value = 14695981039346656037UL  # FNV offset basis
    cdef int i = 0
    while read_name[i] != 0:
        hash_value ^= <int64_t>read_name[i]
        hash_value *= 1099511628211UL  # FNV prime
        i += 1
    if hash_value < 0:
        hash_value = -hash_value
    return hash_value


cdef struct ThreadLocalHashMap:
    kh_seqid_map_t* hash_to_id_map
    kh_seqid_name_map_t* hash_to_name_map
    uint32_t next_local_id


cdef ThreadLocalHashMap* create_thread_local_hash_map() except NULL nogil:
    """Create thread-local hash map for read ID assignment.

    Returns
    -------
    ThreadLocalHashMap*
        Initialized hash map, or NULL on allocation failure
    """
    cdef ThreadLocalHashMap* map = <ThreadLocalHashMap*>malloc(sizeof(ThreadLocalHashMap))
    if not map:
        return NULL
    map.hash_to_id_map = kh_init_seqid_map()
    map.hash_to_name_map = kh_init_seqid_name_map()
    map.next_local_id = 0
    if not map.hash_to_id_map or not map.hash_to_name_map:
        free(map)
        return NULL
    return map


cdef void destroy_thread_local_hash_map(ThreadLocalHashMap* map) noexcept nogil:
    """Free thread-local hash map.

    Parameters
    ----------
    map : ThreadLocalHashMap*
        Hash map to destroy (safe to pass NULL)
    """
    cdef khint_t _k
    cdef char* _nm
    if not map:
        return
    if map.hash_to_id_map:
        kh_destroy_seqid_map(map.hash_to_id_map)
    # free stored names
    if map.hash_to_name_map:
        _k = 0
        while _k < kh_end_seqid_name_map(map.hash_to_name_map):
            if kh_exist_seqid_name_map(map.hash_to_name_map, _k):
                _nm = kh_val_seqid_name_map_wrap(map.hash_to_name_map, _k)[0]
                if _nm != NULL:
                    free(_nm)
            _k += 1
        kh_destroy_seqid_name_map(map.hash_to_name_map)
    free(map)


cdef int64_t extract_read_hash_identifier(bam1_t* alignment) except -1 nogil:
    """Extract read name hash from BAM alignment.

    Parameters
    ----------
    alignment : bam1_t*
        BAM alignment record

    Returns
    -------
    int64_t
        FNV-1a hash of read name, or -1 on error
    """
    cdef const char* query_name
    if not alignment:
        return -1
    query_name = bam_get_qname(alignment)
    if not query_name:
        return -1
    return compute_read_name_hash(<const char*>query_name)
