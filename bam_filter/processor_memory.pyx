# cython: initializedcheck=False
# cython: embedsignature=False
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True

"""Memory pool management for alignment storage.

Provides centralized memory management for alignments with cache-efficient
layout and optional PMD (Post-Mortem Damage) score integration.
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t, uintptr_t
from libc.stddef cimport size_t
from libc.stdlib cimport free, realloc
from libc.string cimport memset

from bam_filter.processor cimport (
    MemoryPool, Alignment, INVALID_SEQUENTIAL_ID, ProcessingStats,
    AlignmentCore, HierarchicalData, DamageCounts, BAMWriterAux,
)
from bam_filter.common_helpers cimport page_size, max_int64
from bam_filter.processor_stats cimport init_processing_stats
from libc.stdlib cimport malloc, calloc, free, realloc

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from "stdlib.h":
    int posix_memalign(void **memptr, size_t alignment, size_t size) nogil

cdef extern from "unistd.h":
    int getpagesize() nogil

cdef extern from "sys/mman.h":
    void* mmap(void* addr, size_t length, int prot, int flags, int fd, long offset) nogil
    int munmap(void* addr, size_t length) nogil
    int madvise(void* addr, size_t length, int advice) nogil
    int PROT_READ
    int PROT_WRITE
    int MAP_PRIVATE
    int MAP_ANONYMOUS
    int MADV_SEQUENTIAL
    int MADV_DONTNEED

cdef void* MAP_FAILED_PTR = <void*>-1

cdef Alignment* global_temp_buffer = NULL
cdef int64_t global_temp_capacity = 0
cdef bint global_temp_in_use = False


cdef void cleanup_presorted_memory() noexcept nogil:
    """Clean up shared temporary buffer for presorted processing.

    Frees module-level temporary buffer and resets tracking variables.
    Safe to call from C-level cleanup paths without GIL.
    """
    global global_temp_buffer, global_temp_capacity, global_temp_in_use

    if global_temp_buffer:
        free(global_temp_buffer)
        global_temp_buffer = NULL
        global_temp_capacity = 0
        global_temp_in_use = False

    bf_nogil_logf_notime(b"CLEANUP", "Presorted cleanup completed (no sort buffers)")


cdef MemoryPool* create_memory_pool(int64_t alignment_count,
                                   uint32_t reference_count,
                                   uint32_t unique_read_count,
                                   int64_t* filtered_ref_lengths,
                                   bint enable_pmd,
                                   int32_t max_threads) except NULL nogil:
    """Create and initialize memory pool for alignment storage.

    Allocates unified memory pool with cache-efficient layout. Uses mmap for
    large allocations (>100MB) to improve memory locality. Validates and copies
    reference lengths from filtered array.

    Parameters
    ----------
    alignment_count : int64_t
        Total number of alignments to store
    reference_count : uint32_t
        Number of filtered references
    unique_read_count : uint32_t
        Number of unique reads
    filtered_ref_lengths : int64_t*
        Array of filtered reference lengths (must match reference_count)
    enable_pmd : bint
        Enable PMD score output
    max_threads : int32_t
        Maximum thread count (unused, retained for API compatibility)

    Returns
    -------
    MemoryPool*
        Initialized memory pool, or NULL on allocation failure

    Notes
    -----
    Uses progressive allocation for small datasets (<1M alignments) to avoid
    memory peaks. Large datasets allocate full capacity upfront.
    """
    cdef MemoryPool* pool
    cdef int64_t initial_alignment_count
    cdef size_t temp_size
    cdef size_t alignments_size
    cdef size_t reference_weights_size
    cdef size_t temp_buffers_size
    cdef size_t total_doubles
    cdef size_t read_starts_size
    cdef size_t read_counts_size
    cdef size_t ref_lengths_size
    cdef size_t unified_doubles_size
    cdef size_t unified_pool_size
    cdef char* memory_cursor
    cdef size_t mmap_threshold = 100 * 1024 * 1024
    cdef int64_t min_len, max_len, total_len
    cdef int invalid_count
    cdef uint32_t i
    cdef void* pool_ptr
    cdef void* align_ptr
    cdef void* unified_ptr

    if posix_memalign(&pool_ptr, 64, sizeof(MemoryPool)) != 0:
        return NULL
    pool = <MemoryPool*>pool_ptr
    if not pool:
        return NULL

    memset(pool, 0, sizeof(MemoryPool))

    # Always allocate full capacity to avoid buffer overflows
    # Progressive allocation was causing segfaults when actual count exceeded initial capacity
    initial_alignment_count = alignment_count
    bf_nogil_logf_notime(
        b"MEMORY",
        "alignment_allocation: mode=full alignments=%ld",
        alignment_count,
    )

    temp_size = max_int64(unique_read_count, reference_count)
    alignments_size = initial_alignment_count * sizeof(Alignment)

    if alignments_size > mmap_threshold:
        bf_nogil_logf_notime(
            b"MEMORY",
            "alignment_buffer: strategy=mmap requested_mb=%.1f",
            alignments_size / (1024.0 * 1024.0),
        )
        pool.alignments = <Alignment*>mmap(NULL, alignments_size, PROT_READ | PROT_WRITE,
                                           MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
        if pool.alignments == MAP_FAILED_PTR:
            bf_nogil_logf_notime(b"WARN", "memory: mmap allocation failed; falling back to malloc")
            if posix_memalign(&align_ptr, 64, alignments_size) != 0:
                free(pool)
                return NULL
            pool.alignments = <Alignment*>align_ptr
            pool.alignments_is_external = True
            pool.mmap_allocation_size = 0
        else:
            bf_nogil_logf_notime(
                b"MEMORY",
                "alignment_buffer: mmap_committed_mb=%.1f",
                alignments_size / (1024.0 * 1024.0),
            )
            pool.alignments_is_external = False
            pool.mmap_allocation_size = alignments_size
            madvise(pool.alignments, alignments_size, MADV_SEQUENTIAL)
    else:
        bf_nogil_logf_notime(
            b"MEMORY",
            "alignment_buffer: strategy=malloc size_mb=%.1f",
            alignments_size / (1024.0 * 1024.0),
        )
        if posix_memalign(&align_ptr, 64, alignments_size) != 0:
            free(pool)
            return NULL
        pool.alignments = <Alignment*>align_ptr
        pool.alignments_is_external = True
        pool.mmap_allocation_size = 0

    reference_weights_size = reference_count
    temp_buffers_size = 2 * temp_size
    total_doubles = reference_weights_size + temp_buffers_size
    read_starts_size = unique_read_count * sizeof(uint64_t)
    read_counts_size = unique_read_count * sizeof(uint32_t)
    ref_lengths_size = reference_count * sizeof(int64_t)
    unified_doubles_size = total_doubles * sizeof(double)

    unified_pool_size = (read_starts_size + read_counts_size + ref_lengths_size +
                         unified_doubles_size + 64)

    if posix_memalign(&unified_ptr, 64, unified_pool_size) != 0:
        if pool.mmap_allocation_size > 0:
            munmap(pool.alignments, pool.mmap_allocation_size)
        else:
            free(pool.alignments)
        free(pool)
        return NULL
    pool.memory_pool = unified_ptr
    pool.pool_capacity = unified_pool_size

    if pool.mmap_allocation_size == 0:
        pool.mmap_allocation_size = 0

    memory_cursor = <char*>pool.memory_pool

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.read_alignment_starts = <uint64_t*>memory_cursor
    memory_cursor += read_starts_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 3) & ~3)
    pool.read_alignment_counts = <uint32_t*>memory_cursor
    memory_cursor += read_counts_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.reference_lengths = <int64_t*>memory_cursor
    memory_cursor += ref_lengths_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.unified_buffer = <double*>memory_cursor
    pool.unified_buffer_size = total_doubles

    pool.reference_weights_offset = 0
    pool.temp_buffer_A_offset = reference_weights_size
    pool.temp_buffer_B_offset = reference_weights_size + temp_size

    pool.alignment_count = 0  # Initially no alignments stored
    pool.alignment_capacity = initial_alignment_count  # Allocated capacity
    pool.original_alignment_count = alignment_count
    pool.reference_count = reference_count
    pool.unique_read_count = unique_read_count
    pool.final_unique_reads = unique_read_count
    pool.memory_owner = True

    pool.pmd_enabled_for_output = enable_pmd

    pool.precomputed_zp_values = NULL
    pool.zp_values_computed = False

    if filtered_ref_lengths:
        for i in range(reference_count):
            pool.reference_lengths[i] = filtered_ref_lengths[i]

        bf_nogil_logf_notime(
            b"MEMORY",
            "reference_lengths: entries=%u source=filtered",
            reference_count,
        )

        if reference_count > 0:
            min_len = pool.reference_lengths[0]
            max_len = min_len
            total_len = 0
            invalid_count = 0

            for i in range(reference_count):
                if pool.reference_lengths[i] <= 0 or pool.reference_lengths[i] > 1000000000:
                    pool.reference_lengths[i] = 1000
                    invalid_count += 1

                if pool.reference_lengths[i] < min_len:
                    min_len = pool.reference_lengths[i]
                if pool.reference_lengths[i] > max_len:
                    max_len = pool.reference_lengths[i]
                total_len += pool.reference_lengths[i]

            bf_nogil_logf_notime(
                b"MEMORY",
                "reference_lengths_stats: min=%ld max=%ld total=%ld adjusted=%d",
                min_len,
                max_len,
                total_len,
                invalid_count,
            )
    else:
        for i in range(reference_count):
            pool.reference_lengths[i] = 1000
        bf_nogil_logf_notime(b"MEMORY", "reference_lengths: entries=%u source=default", reference_count)

    bf_nogil_logf_notime(
        b"MEMORY",
        "memory_pool: status=initialized alignments_mb=%.1f alignments=%ld pmd=%s references=%u storage=%s",
        alignments_size / (1024.0 * 1024.0),
        initial_alignment_count,
        b"enabled" if enable_pmd else b"disabled",
        reference_count,
        b"mmap" if pool.mmap_allocation_size > 0 else b"malloc",
    )

    pool.stats = <ProcessingStats*>calloc(1, sizeof(ProcessingStats))
    if pool.stats:
        init_processing_stats(pool.stats)
    else:
        bf_nogil_logf_notime(b"WARN", "memory: failed to allocate processing_stats")

    return pool


cdef MemoryPool* create_memory_pool_split(int64_t alignment_count,
                                          uint32_t reference_count,
                                          uint32_t unique_read_count,
                                          int64_t* filtered_ref_lengths,
                                          bint enable_pmd,
                                          bint enable_hierarchical,
                                          bint enable_damage_counts,
                                          int32_t max_threads) except NULL nogil:
    """Create memory pool with split array storage for reduced memory usage.

    Uses separate arrays for core alignment data and optional features:
    - AlignmentCore: 16 bytes/alignment (always allocated)
    - HierarchicalData: 12 bytes/alignment (only if enable_hierarchical)
    - DamageCounts: 4 bytes/alignment (only if enable_damage_counts)
    - BAMWriterAux: 8 bytes/alignment (only if enable_pmd for BAM writing)

    For 1 billion alignments:
    - Full Alignment struct: 46 GB
    - Split (no hierarchical): 16 GB (65% reduction)
    - Split (with hierarchical): 28 GB (39% reduction)
    - Split (all features): 40 GB (13% reduction)

    Parameters
    ----------
    alignment_count : int64_t
        Total number of alignments to store
    reference_count : uint32_t
        Number of filtered references
    unique_read_count : uint32_t
        Number of unique reads
    filtered_ref_lengths : int64_t*
        Array of filtered reference lengths
    enable_pmd : bint
        Enable PMD output and BAMWriterAux allocation
    enable_hierarchical : bint
        Enable hierarchical EM (allocates HierarchicalData)
    enable_damage_counts : bint
        Enable damage gamma update (allocates DamageCounts)
    max_threads : int32_t
        Maximum thread count

    Returns
    -------
    MemoryPool*
        Initialized memory pool using split arrays, or NULL on failure
    """
    cdef MemoryPool* pool
    cdef size_t core_size, hier_size, damage_size, aux_size
    cdef size_t temp_size
    cdef size_t reference_weights_size, temp_buffers_size, total_doubles
    cdef size_t read_starts_size, read_counts_size, ref_lengths_size
    cdef size_t unified_doubles_size, unified_pool_size
    cdef char* memory_cursor
    cdef size_t mmap_threshold = 100 * 1024 * 1024
    cdef int64_t min_len, max_len, total_len
    cdef int invalid_count
    cdef uint32_t i
    cdef void* pool_ptr
    cdef void* unified_ptr
    cdef void* ptr
    cdef double total_mb

    if posix_memalign(&pool_ptr, 64, sizeof(MemoryPool)) != 0:
        return NULL
    pool = <MemoryPool*>pool_ptr
    if not pool:
        return NULL

    memset(pool, 0, sizeof(MemoryPool))

    # Calculate sizes
    core_size = alignment_count * sizeof(AlignmentCore)
    cdef size_t indices_size = alignment_count * sizeof(uint32_t)  # read_indices always needed
    hier_size = alignment_count * sizeof(HierarchicalData) if enable_hierarchical else 0
    damage_size = alignment_count * sizeof(DamageCounts) if enable_damage_counts else 0
    aux_size = alignment_count * sizeof(BAMWriterAux) if enable_pmd else 0

    total_mb = (core_size + indices_size + hier_size + damage_size + aux_size) / (1024.0 * 1024.0)
    bf_nogil_logf_notime(
        b"MEMORY",
        "split_allocation: alignments=%ld core_mb=%.1f idx_mb=%.1f hier_mb=%.1f damage_mb=%.1f aux_mb=%.1f total_mb=%.1f",
        alignment_count,
        core_size / (1024.0 * 1024.0),
        indices_size / (1024.0 * 1024.0),
        hier_size / (1024.0 * 1024.0),
        damage_size / (1024.0 * 1024.0),
        aux_size / (1024.0 * 1024.0),
        total_mb,
    )

    # Allocate core array (always)
    if core_size > mmap_threshold:
        pool.alignment_cores = <AlignmentCore*>mmap(NULL, core_size, PROT_READ | PROT_WRITE,
                                                     MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
        if pool.alignment_cores == <AlignmentCore*>MAP_FAILED_PTR:
            bf_nogil_logf_notime(b"WARN", "split_memory: core mmap failed, trying malloc")
            if posix_memalign(&ptr, 64, core_size) != 0:
                free(pool)
                return NULL
            pool.alignment_cores = <AlignmentCore*>ptr
            pool.mmap_allocation_size = 0
        else:
            pool.mmap_allocation_size = core_size
            madvise(pool.alignment_cores, core_size, MADV_SEQUENTIAL)
    else:
        if posix_memalign(&ptr, 64, core_size) != 0:
            free(pool)
            return NULL
        pool.alignment_cores = <AlignmentCore*>ptr
        pool.mmap_allocation_size = 0

    # Allocate read_indices array (always needed for sorting/indexing)
    if indices_size > mmap_threshold:
        pool.read_indices = <uint32_t*>mmap(NULL, indices_size, PROT_READ | PROT_WRITE,
                                            MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
        if pool.read_indices == <uint32_t*>MAP_FAILED_PTR:
            bf_nogil_logf_notime(b"WARN", "split_memory: read_indices mmap failed, trying malloc")
            if posix_memalign(&ptr, 64, indices_size) != 0:
                _cleanup_split_partial(pool)
                return NULL
            pool.read_indices = <uint32_t*>ptr
            pool.read_indices_alloc_size = 0
        else:
            pool.read_indices_alloc_size = indices_size
            madvise(pool.read_indices, indices_size, MADV_SEQUENTIAL)
    else:
        if posix_memalign(&ptr, 64, indices_size) != 0:
            _cleanup_split_partial(pool)
            return NULL
        pool.read_indices = <uint32_t*>ptr
        pool.read_indices_alloc_size = 0

    # Allocate hierarchical array (optional)
    if enable_hierarchical and hier_size > 0:
        if hier_size > mmap_threshold:
            pool.hierarchical = <HierarchicalData*>mmap(NULL, hier_size, PROT_READ | PROT_WRITE,
                                                         MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
            if pool.hierarchical == <HierarchicalData*>MAP_FAILED_PTR:
                bf_nogil_logf_notime(b"WARN", "split_memory: hierarchical mmap failed, trying malloc")
                if posix_memalign(&ptr, 64, hier_size) != 0:
                    _cleanup_split_partial(pool)
                    return NULL
                pool.hierarchical = <HierarchicalData*>ptr
                pool.hierarchical_alloc_size = 0
            else:
                pool.hierarchical_alloc_size = hier_size
                madvise(pool.hierarchical, hier_size, MADV_SEQUENTIAL)
        else:
            if posix_memalign(&ptr, 64, hier_size) != 0:
                _cleanup_split_partial(pool)
                return NULL
            pool.hierarchical = <HierarchicalData*>ptr
            pool.hierarchical_alloc_size = 0

    # Allocate damage counts array (optional)
    if enable_damage_counts and damage_size > 0:
        if damage_size > mmap_threshold:
            pool.damage_counts = <DamageCounts*>mmap(NULL, damage_size, PROT_READ | PROT_WRITE,
                                                      MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
            if pool.damage_counts == <DamageCounts*>MAP_FAILED_PTR:
                bf_nogil_logf_notime(b"WARN", "split_memory: damage_counts mmap failed, trying malloc")
                if posix_memalign(&ptr, 64, damage_size) != 0:
                    _cleanup_split_partial(pool)
                    return NULL
                pool.damage_counts = <DamageCounts*>ptr
                pool.damage_counts_alloc_size = 0
            else:
                pool.damage_counts_alloc_size = damage_size
                madvise(pool.damage_counts, damage_size, MADV_SEQUENTIAL)
        else:
            if posix_memalign(&ptr, 64, damage_size) != 0:
                _cleanup_split_partial(pool)
                return NULL
            pool.damage_counts = <DamageCounts*>ptr
            pool.damage_counts_alloc_size = 0

    # Allocate BAM writer aux array (optional)
    if enable_pmd and aux_size > 0:
        if aux_size > mmap_threshold:
            pool.bam_aux = <BAMWriterAux*>mmap(NULL, aux_size, PROT_READ | PROT_WRITE,
                                                MAP_PRIVATE | MAP_ANONYMOUS, -1, 0)
            if pool.bam_aux == <BAMWriterAux*>MAP_FAILED_PTR:
                bf_nogil_logf_notime(b"WARN", "split_memory: bam_aux mmap failed, trying malloc")
                if posix_memalign(&ptr, 64, aux_size) != 0:
                    _cleanup_split_partial(pool)
                    return NULL
                pool.bam_aux = <BAMWriterAux*>ptr
                pool.bam_aux_alloc_size = 0
            else:
                pool.bam_aux_alloc_size = aux_size
                madvise(pool.bam_aux, aux_size, MADV_SEQUENTIAL)
        else:
            if posix_memalign(&ptr, 64, aux_size) != 0:
                _cleanup_split_partial(pool)
                return NULL
            pool.bam_aux = <BAMWriterAux*>ptr
            pool.bam_aux_alloc_size = 0

    pool.use_split_arrays = True
    pool.alignments = NULL  # Not using legacy array

    # Allocate unified buffer for read/reference metadata (same as original)
    temp_size = max_int64(unique_read_count, reference_count)
    reference_weights_size = reference_count
    temp_buffers_size = 2 * temp_size
    total_doubles = reference_weights_size + temp_buffers_size
    read_starts_size = unique_read_count * sizeof(uint64_t)
    read_counts_size = unique_read_count * sizeof(uint32_t)
    ref_lengths_size = reference_count * sizeof(int64_t)
    unified_doubles_size = total_doubles * sizeof(double)

    unified_pool_size = (read_starts_size + read_counts_size + ref_lengths_size +
                         unified_doubles_size + 64)

    if posix_memalign(&unified_ptr, 64, unified_pool_size) != 0:
        _cleanup_split_partial(pool)
        return NULL
    pool.memory_pool = unified_ptr
    pool.pool_capacity = unified_pool_size

    memory_cursor = <char*>pool.memory_pool

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.read_alignment_starts = <uint64_t*>memory_cursor
    memory_cursor += read_starts_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 3) & ~3)
    pool.read_alignment_counts = <uint32_t*>memory_cursor
    memory_cursor += read_counts_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.reference_lengths = <int64_t*>memory_cursor
    memory_cursor += ref_lengths_size

    memory_cursor = <char*>((<uintptr_t>memory_cursor + 7) & ~7)
    pool.unified_buffer = <double*>memory_cursor
    pool.unified_buffer_size = total_doubles

    pool.reference_weights_offset = 0
    pool.temp_buffer_A_offset = reference_weights_size
    pool.temp_buffer_B_offset = reference_weights_size + temp_size

    pool.alignment_count = 0
    pool.alignment_capacity = alignment_count
    pool.original_alignment_count = alignment_count
    pool.reference_count = reference_count
    pool.unique_read_count = unique_read_count
    pool.final_unique_reads = unique_read_count
    pool.memory_owner = True

    pool.pmd_enabled_for_output = enable_pmd
    pool.hierarchical_em_enabled = enable_hierarchical

    pool.precomputed_zp_values = NULL
    pool.zp_values_computed = False

    # Copy reference lengths
    if filtered_ref_lengths:
        for i in range(reference_count):
            pool.reference_lengths[i] = filtered_ref_lengths[i]

        if reference_count > 0:
            min_len = pool.reference_lengths[0]
            max_len = min_len
            total_len = 0
            invalid_count = 0

            for i in range(reference_count):
                if pool.reference_lengths[i] <= 0 or pool.reference_lengths[i] > 1000000000:
                    pool.reference_lengths[i] = 1000
                    invalid_count += 1

                if pool.reference_lengths[i] < min_len:
                    min_len = pool.reference_lengths[i]
                if pool.reference_lengths[i] > max_len:
                    max_len = pool.reference_lengths[i]
                total_len += pool.reference_lengths[i]
    else:
        for i in range(reference_count):
            pool.reference_lengths[i] = 1000

    bf_nogil_logf_notime(
        b"MEMORY",
        "split_pool: status=initialized total_mb=%.1f hier=%s damage=%s aux=%s refs=%u",
        total_mb,
        b"yes" if enable_hierarchical else b"no",
        b"yes" if enable_damage_counts else b"no",
        b"yes" if enable_pmd else b"no",
        reference_count,
    )

    pool.stats = <ProcessingStats*>calloc(1, sizeof(ProcessingStats))
    if pool.stats:
        init_processing_stats(pool.stats)

    return pool


cdef void _cleanup_split_partial(MemoryPool* pool) noexcept nogil:
    """Clean up partially allocated split arrays on failure."""
    if pool.alignment_cores:
        if pool.mmap_allocation_size > 0:
            munmap(pool.alignment_cores, pool.mmap_allocation_size)
        else:
            free(pool.alignment_cores)
    if pool.read_indices:
        if pool.read_indices_alloc_size > 0:
            munmap(pool.read_indices, pool.read_indices_alloc_size)
        else:
            free(pool.read_indices)
    if pool.hierarchical:
        if pool.hierarchical_alloc_size > 0:
            munmap(pool.hierarchical, pool.hierarchical_alloc_size)
        else:
            free(pool.hierarchical)
    if pool.damage_counts:
        if pool.damage_counts_alloc_size > 0:
            munmap(pool.damage_counts, pool.damage_counts_alloc_size)
        else:
            free(pool.damage_counts)
    if pool.bam_aux:
        if pool.bam_aux_alloc_size > 0:
            munmap(pool.bam_aux, pool.bam_aux_alloc_size)
        else:
            free(pool.bam_aux)
    free(pool)


cdef void destroy_memory_pool(MemoryPool* pool) noexcept nogil:
    """Free all memory associated with memory pool.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool to destroy (safe to pass NULL)
    """
    if not pool:
        return

    if pool.precomputed_zp_values:
        free(pool.precomputed_zp_values)
        pool.precomputed_zp_values = NULL

    if pool.stats:
        free(pool.stats)
        pool.stats = NULL

    if pool.scratch_read_max_probs != NULL:
        free(pool.scratch_read_max_probs)
        pool.scratch_read_max_probs = NULL
    if pool.scratch_survivors_per_read != NULL:
        free(pool.scratch_survivors_per_read)
        pool.scratch_survivors_per_read = NULL
    pool.scratch_unique_read_count = 0

    # Clean up split arrays if using split storage
    if pool.use_split_arrays:
        if pool.alignment_cores:
            if pool.mmap_allocation_size > 0:
                munmap(pool.alignment_cores, pool.mmap_allocation_size)
            else:
                free(pool.alignment_cores)
            pool.alignment_cores = NULL
        if pool.read_indices:
            if pool.read_indices_alloc_size > 0:
                munmap(pool.read_indices, pool.read_indices_alloc_size)
            else:
                free(pool.read_indices)
            pool.read_indices = NULL
        if pool.hierarchical:
            if pool.hierarchical_alloc_size > 0:
                munmap(pool.hierarchical, pool.hierarchical_alloc_size)
            else:
                free(pool.hierarchical)
            pool.hierarchical = NULL
        if pool.damage_counts:
            if pool.damage_counts_alloc_size > 0:
                munmap(pool.damage_counts, pool.damage_counts_alloc_size)
            else:
                free(pool.damage_counts)
            pool.damage_counts = NULL
        if pool.bam_aux:
            if pool.bam_aux_alloc_size > 0:
                munmap(pool.bam_aux, pool.bam_aux_alloc_size)
            else:
                free(pool.bam_aux)
            pool.bam_aux = NULL
    elif pool.alignments and not pool.alignments_is_external:
        # Clean up legacy alignment array
        if pool.mmap_allocation_size > 0:
            munmap(pool.alignments, pool.mmap_allocation_size)
        else:
            free(pool.alignments)
        pool.alignments = NULL

    if pool.memory_owner and pool.memory_pool:
        free(pool.memory_pool)
        pool.memory_pool = NULL

    free(pool)


cdef void cleanup_em_intermediate_memory(MemoryPool* pool) noexcept nogil:
    """Free EM algorithm intermediate memory after probability filtering.

    This function frees memory that is only needed during EM iterations,
    reducing peak memory usage before graph analysis stages. It preserves:
    - gamma_values: Needed for output (ancient/modern classification)
    - precomputed_zp_values: Needed for BAM writing
    - unified_buffer: Part of memory pool, contains reference lengths

    Memory freed:
    - squarem_block: SQUAREM acceleration scratch memory
    - eta_values: Hierarchical EM intermediate (logit of gamma)
    - S_anc_accum: Hierarchical EM ancient accumulator
    - S_mod_accum: Hierarchical EM modern accumulator
    - scratch_read_max_probs: Filtering scratch arrays
    - scratch_survivors_per_read: Filtering scratch arrays

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool to clean up (safe to pass NULL)
    """
    cdef size_t freed_bytes = 0

    if not pool:
        return

    # Free SQUAREM scratch memory
    if pool.squarem_block != NULL:
        freed_bytes += pool.squarem_block_capacity * sizeof(double)
        free(pool.squarem_block)
        pool.squarem_block = NULL
        pool.squarem_block_capacity = 0

    # Free hierarchical EM intermediate arrays (but NOT gamma_values - needed for output)
    if pool.eta_values != NULL:
        freed_bytes += pool.reference_count * sizeof(double)
        free(pool.eta_values)
        pool.eta_values = NULL

    if pool.S_anc_accum != NULL:
        freed_bytes += pool.reference_count * sizeof(double)
        free(pool.S_anc_accum)
        pool.S_anc_accum = NULL

    if pool.S_mod_accum != NULL:
        freed_bytes += pool.reference_count * sizeof(double)
        free(pool.S_mod_accum)
        pool.S_mod_accum = NULL

    # Free filtering scratch arrays
    if pool.scratch_read_max_probs != NULL:
        freed_bytes += pool.scratch_unique_read_count * sizeof(float)
        free(pool.scratch_read_max_probs)
        pool.scratch_read_max_probs = NULL

    if pool.scratch_survivors_per_read != NULL:
        freed_bytes += pool.scratch_unique_read_count * sizeof(int32_t)
        free(pool.scratch_survivors_per_read)
        pool.scratch_survivors_per_read = NULL

    pool.scratch_unique_read_count = 0

    if freed_bytes > 0:
        bf_nogil_logf_notime(b"MEMORY", "cleanup_em_intermediate: freed %.1f MB",
                            <double>freed_bytes / (1024.0 * 1024.0))


cdef int shrink_memory_pool(MemoryPool* pool) except -1 nogil:
    """Shrink memory pool to used capacity.

    Releases unused memory from alignment array. For mmap allocations, unmaps
    tail pages. For malloc allocations, reallocs to actual size.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool to shrink

    Returns
    -------
    int
        0 on success, -1 on error
    """
    cdef size_t used_bytes
    cdef size_t ps
    cdef size_t keep_bytes
    cdef size_t old_bytes
    cdef size_t tail_len
    cdef size_t target
    cdef char* base
    cdef void* tail
    cdef Alignment* newp

    if not pool or not pool.alignments:
        return 0

    used_bytes = <size_t>(pool.alignment_count) * sizeof(Alignment)

    if used_bytes == 0:
        if pool.mmap_allocation_size > 0:
            munmap(pool.alignments, pool.mmap_allocation_size)
            pool.mmap_allocation_size = 0
        else:
            free(pool.alignments)
        pool.alignments = NULL
        pool.original_alignment_count = 0
        return 0

    if pool.mmap_allocation_size > 0:
        ps = <size_t>(getpagesize() if getpagesize() > 0 else 4096)
        keep_bytes = (used_bytes + ps - 1) & ~(ps - 1)
        old_bytes = pool.mmap_allocation_size

        if keep_bytes < old_bytes:
            base = <char*>pool.alignments
            tail = <void*>(base + keep_bytes)
            tail_len = old_bytes - keep_bytes

            madvise(tail, tail_len, MADV_DONTNEED)
            if munmap(tail, tail_len) != 0:
                return -1

            pool.mmap_allocation_size = keep_bytes
            pool.original_alignment_count = <int64_t>(keep_bytes // sizeof(Alignment))
        return 0

    target = used_bytes if used_bytes > 0 else 1
    newp = <Alignment*>realloc(pool.alignments, target)
    if not newp and used_bytes > 0:
        return -1
    if newp:
        pool.alignments = newp
    pool.alignment_capacity = pool.alignment_count  # After shrink, capacity matches count
    pool.original_alignment_count = pool.alignment_count
    return 0
