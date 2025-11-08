# cython: language_level=3
from libc.stdint cimport int64_t, uint32_t

# HTSlib types
from bam_filter.processor_types cimport (
    bam1_t,
    bam_get_qname,
)

cdef extern from *:
    """
    /* FNV-1a hash helper exported for use by other modules */
    """

cdef int64_t compute_read_name_hash(const char* read_name) noexcept nogil

# Thread-local hash map type (opaque mapping to seqid_khash)


cdef extern from "seqid_khash.h":
    ctypedef long long khint64_t
    ctypedef int khint_t
    ctypedef khint_t khiter_t

    ctypedef struct kh_seqid_map_t:
        pass
    ctypedef struct kh_seqid_name_map_t:
        pass

    kh_seqid_map_t* kh_init_seqid_map() nogil
    void kh_destroy_seqid_map(kh_seqid_map_t*) nogil
    khint_t kh_put_seqid_map(kh_seqid_map_t*, khint64_t, int*) nogil
    khint_t kh_get_seqid_map(kh_seqid_map_t*, khint64_t) nogil
    khint_t kh_end_seqid_map(kh_seqid_map_t*) nogil
    int kh_exist_seqid_map(kh_seqid_map_t*, khint_t) nogil
    int* kh_val_seqid_map_wrap(kh_seqid_map_t*, khint_t) nogil
    khint_t kh_size(kh_seqid_map_t*) nogil

    # name map (hash -> char*)
    kh_seqid_name_map_t* kh_init_seqid_name_map() nogil
    void kh_destroy_seqid_name_map(kh_seqid_name_map_t*) nogil
    khint_t kh_put_seqid_name_map(kh_seqid_name_map_t*, khint64_t, int*) nogil
    khint_t kh_get_seqid_name_map(kh_seqid_name_map_t*, khint64_t) nogil
    khint_t kh_end_seqid_name_map(kh_seqid_name_map_t*) nogil
    int kh_exist_seqid_name_map(kh_seqid_name_map_t*, khint_t) nogil
    char** kh_val_seqid_name_map_wrap(kh_seqid_name_map_t*, khint_t) nogil


cdef extern from *:
    """
    #define kh_key_seqid_map(h,k) kh_key(h,k)
    """
    khint64_t kh_key_seqid_map(kh_seqid_map_t*, khint_t) nogil


cdef struct ThreadLocalHashMap:
    kh_seqid_map_t* hash_to_id_map
    kh_seqid_name_map_t* hash_to_name_map
    uint32_t next_local_id

cdef ThreadLocalHashMap* create_thread_local_hash_map() except NULL nogil
cdef void destroy_thread_local_hash_map(ThreadLocalHashMap* map) noexcept nogil

cdef int64_t extract_read_hash_identifier(bam1_t* alignment) except -1 nogil
