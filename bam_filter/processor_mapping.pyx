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

"""Reference ID mapping and BAM header management.

Manages mapping between original and filtered reference IDs, and generates
BAM headers containing only references with retained alignments.
"""

from libc.stdint cimport uint32_t, int64_t, int32_t
from libc.stddef cimport size_t
from libc.stdlib cimport malloc, free, calloc
from libc.string cimport strcat, strncat, strlen, strstr, strchr, strncmp
from libc.stdio cimport sprintf

from bam_filter.processor cimport MemoryPool, AlignmentCore
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor cimport sam_hdr_t

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

# HTSlib declarations needed by mapping helpers
cdef extern from "htslib/sam.h":
    const char* sam_hdr_str(sam_hdr_t* header) nogil
    sam_hdr_t* sam_hdr_parse(size_t l_text, const char* text) nogil
    int32_t sam_hdr_nref(sam_hdr_t* header) nogil
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid) nogil
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid) nogil

cdef ReferenceMapping* create_reference_mapping(MemoryPool* pool,
                                                sam_hdr_t* original_header) noexcept nogil:
    """Create bidirectional reference ID mapping from pool alignments.

    Scans all alignments in pool to determine which original references are used,
    creates compact numbering, and builds bidirectional mapping arrays.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool with alignments (using original reference IDs)
    original_header : sam_hdr_t*
        Original BAM header with all references

    Returns
    -------
    ReferenceMapping*
        Mapping structure with old_to_new and new_to_old arrays, or NULL on error
    """
    cdef ReferenceMapping* mapping = <ReferenceMapping*>malloc(sizeof(ReferenceMapping))
    cdef char* ref_seen = NULL
    cdef uint32_t ref_id, new_tid = 0
    cdef int64_t alignment_idx

    if not mapping or not pool or not original_header:
        if mapping: free(mapping)
        return NULL

    mapping.n_original_refs = <uint32_t>sam_hdr_nref(original_header)
    mapping.n_retained_refs = 0

    mapping.old_to_new_tid = <uint32_t*>malloc(mapping.n_original_refs * sizeof(uint32_t))
    ref_seen = <char*>calloc(mapping.n_original_refs, sizeof(char))

    if not mapping.old_to_new_tid or not ref_seen:
        if mapping.old_to_new_tid: free(mapping.old_to_new_tid)
        if ref_seen: free(ref_seen)
        free(mapping)
        return NULL

    for ref_id in range(mapping.n_original_refs):
        mapping.old_to_new_tid[ref_id] = <uint32_t>(-1)

    for alignment_idx in range(pool.alignment_count):
        ref_id = pool.alignment_cores[alignment_idx].reference_index
        if ref_id < mapping.n_original_refs and not ref_seen[ref_id]:
            ref_seen[ref_id] = 1
            mapping.n_retained_refs += 1

    mapping.new_to_old_tid = <uint32_t*>malloc(mapping.n_retained_refs * sizeof(uint32_t))
    if not mapping.new_to_old_tid:
        free(mapping.old_to_new_tid)
        free(ref_seen)
        free(mapping)
        return NULL

    new_tid = 0
    for ref_id in range(mapping.n_original_refs):
        if ref_seen[ref_id]:
            mapping.old_to_new_tid[ref_id] = new_tid
            mapping.new_to_old_tid[new_tid] = ref_id
            new_tid += 1

    free(ref_seen)

    bf_nogil_logf_notime(
        b"MAPPING",
        "mapping_summary: original_refs=%u retained_refs=%u",
        mapping.n_original_refs,
        mapping.n_retained_refs,
    )

    return mapping


cdef void destroy_reference_mapping(ReferenceMapping* mapping) noexcept nogil:
    """Free reference mapping structure.

    Parameters
    ----------
    mapping : ReferenceMapping*
        Mapping to destroy (safe to pass NULL)
    """
    if not mapping:
        return
    if mapping.old_to_new_tid:
        free(mapping.old_to_new_tid)
    if mapping.new_to_old_tid:
        free(mapping.new_to_old_tid)
    free(mapping)


cdef int update_reference_mapping_after_filtering(ReferenceMapping* mapping,
                                                  MemoryPool* pool) noexcept nogil:
    """Update reference mapping after additional filtering removes references.

    Compacts reference numbering to eliminate gaps from filtered references,
    updates alignment reference IDs in pool, and rebuilds mapping arrays.

    Parameters
    ----------
    mapping : ReferenceMapping*
        Existing reference mapping
    pool : MemoryPool*
        Memory pool with alignments (reference IDs updated in-place)

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef char* ref_active = <char*>calloc(mapping.n_retained_refs, sizeof(char))
    cdef uint32_t* old_to_new_compact = <uint32_t*>malloc(mapping.n_retained_refs * sizeof(uint32_t))
    cdef uint32_t new_compact_count = 0
    cdef int64_t i

    if not ref_active or not old_to_new_compact:
        if ref_active: free(ref_active)
        if old_to_new_compact: free(old_to_new_compact)
        return -1

    for i in range(pool.alignment_count):
        compact_ref = pool.alignment_cores[i].reference_index
        if compact_ref < mapping.n_retained_refs:
            ref_active[compact_ref] = 1

    for i in range(mapping.n_retained_refs):
        if ref_active[i]:
            old_to_new_compact[i] = new_compact_count
            new_compact_count += 1
        else:
            old_to_new_compact[i] = <uint32_t>(-1)

    for i in range(pool.alignment_count):
        old_compact_ref = pool.alignment_cores[i].reference_index
        new_compact_ref = old_to_new_compact[old_compact_ref]
        pool.alignment_cores[i].reference_index = new_compact_ref

    cdef uint32_t* new_new_to_old = <uint32_t*>malloc(new_compact_count * sizeof(uint32_t))
    if not new_new_to_old:
        free(ref_active)
        free(old_to_new_compact)
        return -1

    cdef uint32_t new_idx = 0
    for old_compact in range(mapping.n_retained_refs):
        if ref_active[old_compact]:
            original_ref = mapping.new_to_old_tid[old_compact]
            new_new_to_old[new_idx] = original_ref
            new_idx += 1

    for original_ref in range(mapping.n_original_refs):
        old_compact = mapping.old_to_new_tid[original_ref]
        if old_compact != <uint32_t>(-1) and old_compact < mapping.n_retained_refs:
            new_compact = old_to_new_compact[old_compact]
            mapping.old_to_new_tid[original_ref] = new_compact

    free(mapping.new_to_old_tid)
    mapping.new_to_old_tid = new_new_to_old
    mapping.n_retained_refs = new_compact_count

    if pool is not NULL:
        pool.reference_count = mapping.n_retained_refs

    bf_nogil_logf_notime(
        b"MAPPING",
        "mapping_summary: compacted_refs_before=%u compacted_refs_after=%u",
        mapping.n_retained_refs,
        new_compact_count,
    )

    free(ref_active)
    free(old_to_new_compact)
    return 0


cdef sam_hdr_t* create_filtered_header_efficient(sam_hdr_t* original_header,
                                                 ReferenceMapping* mapping) noexcept nogil:
    """Create a new BAM header containing only retained references.

    Builds a new SAM header text by extracting @SQ lines for references
    retained in ``mapping`` and concatenating them with the non-@SQ
    portions of the original header. The resulting header is parsed
    with HTSlib and returned as a new ``sam_hdr_t*``. The caller is
    responsible for destroying the returned header with
    ``sam_hdr_destroy``.

    Parameters
    ----------
    original_header : sam_hdr_t*
        Original BAM header from which to extract metadata.
    mapping : ReferenceMapping*
        Mapping describing which original references are retained.

    Returns
    -------
    sam_hdr_t*
        Newly parsed header containing only retained references, or
        NULL on error.
    """
    cdef char* new_header_text = NULL
    cdef const char* original_text = NULL
    cdef size_t header_size = 0
    cdef char* sq_content = NULL
    cdef uint32_t new_tid, original_tid
    cdef const char* ref_name = NULL
    cdef int64_t ref_length
    cdef int bytes_written
    cdef size_t current_pos = 0
    cdef sam_hdr_t* filtered_header = NULL
    cdef const char* sq_start = NULL
    cdef const char* sq_end = NULL
    cdef const char* line_start = NULL
    cdef size_t pre_sq_len = 0
    cdef size_t post_sq_len = 0
    cdef size_t sq_content_len = 0

    if not original_header or not mapping:
        return NULL

    original_text = sam_hdr_str(original_header)
    if not original_text:
        return NULL

    sq_start = strstr(original_text, b"@SQ\t")
    sq_end = sq_start

    if sq_start:
        line_start = sq_start
        while line_start:
            line_start = strchr(line_start, 10)
            if line_start:
                line_start += 1
                if line_start[0] != <char>64 or strncmp(<char*>line_start, b"@SQ\t", 4) != 0:
                    sq_end = line_start
                    break
            else:
                sq_end = original_text + strlen(original_text)
                break
    else:
        sq_start = strstr(original_text, b"@HD\t")
        if sq_start:
            sq_start = strchr(sq_start, 10)
            if sq_start:
                sq_start += 1
            else:
                sq_start = original_text
        else:
            sq_start = original_text
        sq_end = sq_start

    pre_sq_len = sq_start - original_text
    post_sq_len = strlen(sq_end)

    sq_content_len = mapping.n_retained_refs * 100
    sq_content = <char*>malloc(sq_content_len)
    if not sq_content:
        return NULL

    sq_content[0] = 0
    current_pos = 0

    for new_tid in range(mapping.n_retained_refs):
        original_tid = mapping.new_to_old_tid[new_tid]
        ref_name = sam_hdr_tid2name(original_header, original_tid)
        ref_length = sam_hdr_tid2len(original_header, original_tid)

        if ref_name and current_pos + 100 <= sq_content_len:
            bytes_written = sprintf(sq_content + current_pos,
                                   "@SQ\tSN:%s\tLN:%ld\n",
                                   ref_name, ref_length)
            if bytes_written > 0:
                current_pos += bytes_written

    header_size = pre_sq_len + current_pos + post_sq_len + 1
    new_header_text = <char*>malloc(header_size)
    if not new_header_text:
        free(sq_content)
        return NULL

    new_header_text[0] = 0

    if pre_sq_len > 0:
        strncat(new_header_text, original_text, pre_sq_len)

    strcat(new_header_text, sq_content)

    if post_sq_len > 0:
        strcat(new_header_text, sq_end)

    filtered_header = sam_hdr_parse(strlen(new_header_text), new_header_text)

    free(sq_content)
    free(new_header_text)
    return filtered_header


cdef int remap_alignment_reference_ids(MemoryPool* pool,
                                      ReferenceMapping* mapping) noexcept nogil:
    """Remap alignment reference IDs in a MemoryPool using a ReferenceMapping.

    Scans all alignments stored in ``pool`` and replaces their original
    reference IDs with compact IDs defined by ``mapping.old_to_new_tid``.
    Returns 0 on success and -1 on error (for example when an alignment
    references an unknown original reference).

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments whose reference IDs will be updated.
    mapping : ReferenceMapping*
        Mapping from original reference IDs to compact (new) IDs.

    Returns
    -------
    int
        0 on success, -1 on error.
    """
    cdef int64_t alignment_idx
    cdef uint32_t old_tid, new_tid
    cdef int64_t remapped_count = 0

    if not pool or not mapping:
        return -1

    for alignment_idx in range(pool.alignment_count):
        old_tid = pool.alignment_cores[alignment_idx].reference_index

        if old_tid < mapping.n_original_refs:
            new_tid = mapping.old_to_new_tid[old_tid]
            if new_tid != <uint32_t>(-1):
                pool.alignment_cores[alignment_idx].reference_index = new_tid
                remapped_count += 1
            else:
                bf_nogil_logf_notime(
                    b"ERROR",
                    "Alignment references filtered reference %u",
                    old_tid,
                )
                return -1
        else:
            bf_nogil_logf_notime(
                b"ERROR",
                "Invalid reference ID %u (max %u)",
                old_tid,
                mapping.n_original_refs,
            )
            return -1

    bf_nogil_logf_notime(
        b"MAPPING",
        "Remapped %ld alignment reference IDs",
        remapped_count,
    )
    return 0
