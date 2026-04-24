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

"""EM algorithm for read reassignment with SQUAREM acceleration."""

from cython.parallel import prange, threadid

from libc.math cimport exp, fabs, fmax, fmin, log, log2, sqrt as libc_sqrt, pow as libc_pow, INFINITY
from libc.stdint cimport int32_t, int64_t, uint16_t, uint32_t, uint64_t, uint8_t, uintptr_t
from libc.stdlib cimport calloc, free, malloc
from libc.string cimport memcpy, memset
from libc.float cimport DBL_EPSILON

from bam_filter.processor cimport (
    MemoryPool, Alignment,
    AlignmentCore, HierarchicalData, DamageCounts, BAMWriterAux,
)
from bam_filter.processor_fast_math cimport stable_log_sum_exp, safe_normalize_weights
from bam_filter.processor_graph cimport (
    update_ancientness_field, ReferenceStats, init_ancientness_arrays,
    calculate_reference_coverage, calculate_reference_coverage_batched
)
from bam_filter.unified_damage cimport (
    UnifiedDamageContext, RefDamageCounts, RefDamageParams,
    create_damage_context, destroy_damage_context, reset_damage_counts,
    fit_damage_model, compute_outputs, DAMAGE_MAX_POSITION,
    estimate_baselines, initialize_parameters, set_tau_from_pmd
)
from bam_filter.processor_pmd cimport RefDamageStats, PMDCurve
from bam_filter.stats_rle cimport (
    calculate_weighted_spatial_entropy, calculate_norm_weighted_gini,
    estimate_histogram_bins_from_length
)
from cpython.pycapsule cimport PyCapsule_GetPointer

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef extern from "stdlib.h":
    void qsort(void* base, size_t nmemb, size_t size,
               int (*compar)(const void*, const void*)) nogil


# =============================================================================
# Constants
# =============================================================================

DEF MAX_SCRATCH_SIZE = 64
DEF POSTERIOR_BUF_INIT_CAP = 4096
DEF POSTERIOR_BUF_MAX_CAP = 10000000  # 10M events max per thread to prevent OOM
DEF STREAMING_HIST_BINS = 200  # Fixed bin count for streaming coverage histograms

cdef double NEG_INF = -1e20
cdef double LOG_ZERO = -1e10


# =============================================================================
# Split Array Accessor Functions
# =============================================================================
# Direct access to split arrays (AlignmentCore*, HierarchicalData*, etc.)

cdef inline uint32_t aln_ref_idx(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get reference index for alignment at idx."""
    return pool.alignment_cores[idx].reference_index


cdef inline float aln_score(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get alignment score for alignment at idx."""
    return pool.alignment_cores[idx].alignment_score


cdef inline uint32_t aln_position(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get alignment position for alignment at idx."""
    return pool.alignment_cores[idx].alignment_position


cdef inline uint16_t aln_length(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get aligned length for alignment at idx."""
    return pool.alignment_cores[idx].aligned_length


cdef inline float aln_damage_llr(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get damage log-likelihood ratio for alignment at idx."""
    return pool.hierarchical[idx].damage_llr


cdef inline float aln_log_L_anc(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get log P(alignment | ancient DNA model) for alignment at idx."""
    return pool.hierarchical[idx].log_L_anc


cdef inline float aln_log_L_mod(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get log P(alignment | modern DNA model) for alignment at idx."""
    return pool.hierarchical[idx].log_L_mod


cdef inline uint8_t aln_ct_5p(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get C->T count at 5' end for alignment at idx."""
    return pool.damage_counts[idx].ct_5p_count


cdef inline uint8_t aln_ga_3p(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get G->A count at 3' end for alignment at idx."""
    return pool.damage_counts[idx].ga_3p_count


cdef inline uint8_t aln_c_at_5p(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get C bases at 5' end for alignment at idx."""
    return pool.damage_counts[idx].c_at_5p_count


cdef inline uint8_t aln_g_at_3p(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get G bases at 3' end for alignment at idx."""
    return pool.damage_counts[idx].g_at_3p_count


cdef inline uint32_t aln_read_idx(MemoryPool* pool, int64_t idx) noexcept nogil:
    """Get read index for alignment at idx."""
    return pool.read_indices[idx]


# =============================================================================
# Posterior Event Buffer for Weighted Coverage Authenticity (Path B)
# =============================================================================

cdef struct PosteriorEvent:
    uint64_t aln_idx      # Index into pool.alignments
    float posterior       # Responsibility r_ij

cdef struct PosteriorEventBuffer:
    PosteriorEvent* events
    uint64_t len
    uint64_t cap


cdef inline int init_posterior_buffer(PosteriorEventBuffer* buf, uint64_t cap0) noexcept nogil:
    """Initialize a posterior event buffer with initial capacity."""
    buf.events = <PosteriorEvent*>malloc(cap0 * sizeof(PosteriorEvent))
    if buf.events == NULL:
        buf.len = 0
        buf.cap = 0
        return -1
    buf.len = 0
    buf.cap = cap0
    return 0


cdef inline int append_posterior(PosteriorEventBuffer* buf, uint64_t aln_idx, float posterior) noexcept nogil:
    """Append a posterior event to the buffer, growing if needed.

    Returns 0 on success, -1 on allocation failure, 1 if buffer is at max capacity (event skipped).
    """
    cdef uint64_t new_cap
    cdef PosteriorEvent* new_events

    if buf.len >= buf.cap:
        # Check max cap to prevent OOM on huge datasets
        if buf.cap >= POSTERIOR_BUF_MAX_CAP:
            return 1  # Buffer full, skip this event (not an error)
        new_cap = buf.cap * 2 if buf.cap > 0 else POSTERIOR_BUF_INIT_CAP
        if new_cap > POSTERIOR_BUF_MAX_CAP:
            new_cap = POSTERIOR_BUF_MAX_CAP
        new_events = <PosteriorEvent*>realloc(buf.events, new_cap * sizeof(PosteriorEvent))
        if new_events == NULL:
            return -1
        buf.events = new_events
        buf.cap = new_cap

    buf.events[buf.len].aln_idx = aln_idx
    buf.events[buf.len].posterior = posterior
    buf.len += 1
    return 0


cdef inline void free_posterior_buffer(PosteriorEventBuffer* buf) noexcept nogil:
    """Free the posterior event buffer."""
    if buf.events != NULL:
        free(buf.events)
    buf.events = NULL
    buf.len = 0
    buf.cap = 0


cdef extern from "stdlib.h":
    void* realloc(void* ptr, size_t size) nogil


# =============================================================================
# Streaming Coverage Histogram for Memory-Efficient Path B
# =============================================================================
#
# Instead of buffering O(n_alignments) events, we use O(n_refs × n_bins) memory.
# Each reference has a fixed-size histogram delta array. During read iteration,
# we accumulate +posterior at start_bin and -posterior at end_bin. After iteration,
# prefix-sum converts deltas to coverage histogram for entropy/Gini computation.

cdef struct StreamingCoverageHist:
    double* deltas     # [n_bins+1] delta encoding: +posterior at start, -posterior at end
    int64_t ref_length # Reference length for bin mapping
    double total_mass  # Sum of posteriors (for validation)


cdef inline int init_streaming_hist(StreamingCoverageHist* hist, int64_t ref_length) noexcept nogil:
    """Initialize streaming histogram for a reference."""
    hist.deltas = <double*>calloc(STREAMING_HIST_BINS + 1, sizeof(double))
    if hist.deltas == NULL:
        return -1
    hist.ref_length = ref_length
    hist.total_mass = 0.0
    return 0


cdef inline void free_streaming_hist(StreamingCoverageHist* hist) noexcept nogil:
    """Free streaming histogram."""
    if hist.deltas != NULL:
        free(hist.deltas)
        hist.deltas = NULL


cdef inline void reset_streaming_hist(StreamingCoverageHist* hist) noexcept nogil:
    """Reset histogram deltas to zero."""
    if hist.deltas != NULL:
        memset(hist.deltas, 0, (STREAMING_HIST_BINS + 1) * sizeof(double))
    hist.total_mass = 0.0


cdef inline void accumulate_alignment_to_hist(
    StreamingCoverageHist* hist,
    int64_t start_pos,
    int64_t end_pos,
    double posterior
) noexcept nogil:
    """Accumulate an alignment's posterior-weighted coverage to histogram.

    Uses delta encoding: +posterior at start_bin, -posterior at end_bin.
    After all alignments, prefix sum converts deltas to actual coverage.
    """
    cdef int64_t start_bin, end_bin
    cdef double bin_scale

    if hist.deltas == NULL or hist.ref_length <= 0:
        return

    # Clamp positions
    if start_pos < 0:
        start_pos = 0
    if end_pos > hist.ref_length:
        end_pos = hist.ref_length
    if start_pos >= end_pos:
        return

    # Map positions to bins
    bin_scale = <double>STREAMING_HIST_BINS / <double>hist.ref_length
    start_bin = <int64_t>(<double>start_pos * bin_scale)
    end_bin = <int64_t>(<double>end_pos * bin_scale)

    if start_bin >= STREAMING_HIST_BINS:
        start_bin = STREAMING_HIST_BINS - 1
    if end_bin > STREAMING_HIST_BINS:
        end_bin = STREAMING_HIST_BINS

    # Delta encoding
    hist.deltas[start_bin] += posterior
    hist.deltas[end_bin] -= posterior
    hist.total_mass += posterior


cdef inline void finalize_hist_and_compute_metrics(
    StreamingCoverageHist* hist,
    double* out_norm_entropy,
    double* out_norm_gini,
    double* out_auth_score,
    double scale
) noexcept nogil:
    """Convert delta encoding to coverage histogram and compute metrics.

    Performs prefix sum to get coverage, then computes entropy/Gini.
    """
    cdef double coverage[STREAMING_HIST_BINS]
    cdef double current_depth = 0.0
    cdef int i
    cdef double norm_entropy, norm_gini, raw_score

    if hist.deltas == NULL or hist.total_mass < 1e-10:
        out_norm_entropy[0] = 0.5
        out_norm_gini[0] = 0.5
        out_auth_score[0] = 0.5
        return

    # Prefix sum: deltas -> actual coverage per bin
    for i in range(STREAMING_HIST_BINS):
        current_depth += hist.deltas[i]
        coverage[i] = current_depth

    # Compute metrics using existing functions
    norm_entropy = calculate_weighted_spatial_entropy(coverage, STREAMING_HIST_BINS)
    norm_gini = calculate_norm_weighted_gini(coverage, STREAMING_HIST_BINS)

    # Authenticity score: sigmoid(scale * (entropy - gini))
    raw_score = norm_entropy - norm_gini
    out_norm_entropy[0] = norm_entropy
    out_norm_gini[0] = norm_gini
    out_auth_score[0] = 1.0 / (1.0 + exp(-scale * raw_score))


# Interval event for weighted coverage sweep (legacy event-based approach)
cdef struct IntervalEvent:
    int64_t pos
    double delta  # +weight at start, -weight at end


cdef int compare_interval_events(const void* a, const void* b) noexcept nogil:
    """Compare interval events by position for qsort."""
    cdef IntervalEvent* ea = <IntervalEvent*>a
    cdef IntervalEvent* eb = <IntervalEvent*>b
    if ea.pos < eb.pos:
        return -1
    elif ea.pos > eb.pos:
        return 1
    return 0


# =============================================================================
# Helper Functions
# =============================================================================

cdef inline double* get_reference_weights(MemoryPool* pool) noexcept nogil:
    """Get pointer to reference weights array in unified buffer."""
    return pool.unified_buffer + pool.reference_weights_offset


cdef extern from *:
    """
    #ifdef __GNUC__
    #define PREFETCH_R(addr) __builtin_prefetch((addr), 0, 3)
    #define PREFETCH_W(addr) __builtin_prefetch((addr), 1, 3)
    #else
    #define PREFETCH_R(addr) ((void)0)
    #define PREFETCH_W(addr) ((void)0)
    #endif
    """
    void PREFETCH_R(void* addr) nogil
    void PREFETCH_W(void* addr) nogil


cdef inline void PREFETCH_READ(void* ptr) noexcept nogil:
    PREFETCH_R(ptr)


cdef inline void PREFETCH_WRITE(void* ptr) noexcept nogil:
    PREFETCH_W(ptr)


cdef inline double stable_sigmoid(double x) noexcept nogil:
    """Numerically stable sigmoid function.

    For x >= 0: σ(x) = 1 / (1 + exp(-x))
    For x < 0:  σ(x) = exp(x) / (1 + exp(x))
    """
    cdef double ex
    if x >= 0:
        return 1.0 / (1.0 + exp(-x))
    else:
        ex = exp(x)
        return ex / (1.0 + ex)


cdef inline double stable_logit(double p) noexcept nogil:
    """Numerically stable logit function: log(p/(1-p)).

    Clamps p to (1e-10, 1-1e-10) to avoid log(0).
    """
    cdef double p_clamped = fmax(1e-10, fmin(1.0 - 1e-10, p))
    return log(p_clamped) - log(1.0 - p_clamped)


cdef inline double logodds_to_log_prob(double logodds) noexcept nogil:
    """Convert log-odds to log(probability) with numerical stability.

    Given logodds = log(p / (1-p)), compute log(p).

    For logodds >= 0: log(p) = -log(1 + exp(-logodds))
    For logodds < 0:  log(p) = logodds - log(1 + exp(logodds))
    """
    if logodds >= 0:
        return -log(1.0 + exp(-logodds))
    else:
        return logodds - log(1.0 + exp(logodds))


cdef inline double logodds_to_log_1mp(double logodds) noexcept nogil:
    """Convert log-odds to log(1-probability) with numerical stability.

    Given logodds = log(p / (1-p)), compute log(1-p).

    For logodds >= 0: log(1-p) = -logodds - log(1 + exp(-logodds))
    For logodds < 0:  log(1-p) = -log(1 + exp(logodds))
    """
    if logodds >= 0:
        return -logodds - log(1.0 + exp(-logodds))
    else:
        return -log(1.0 + exp(logodds))


# =============================================================================
# Sample-Level Ancientness Gate (π)
# =============================================================================
# π = P(sample is ancient) computed from global damage evidence.
# When π is low (modern sample), the hierarchical damage model is dampened.
# Formula: π = clamp(omega × D1 / (D1 + baseline), 0.01, 0.99)
# where D1 is terminal damage rate and baseline is sequencing error.

cdef double compute_sample_pi(PMDCurve* pmd_curve, double pi_override) noexcept nogil:
    """Compute sample-level P(ancient) from PMD curve parameters.

    This gates the hierarchical damage model based on sample-level evidence.
    For ancient samples (high omega, high D1): π ≈ 0.7-0.9
    For modern samples (low omega, D1 ≈ baseline): π ≈ 0.25-0.5

    Parameters
    ----------
    pmd_curve : PMDCurve*
        Fitted PMD curve with omega, D_5p_noncpg[0], and C_background
    pi_override : double
        If > 0, use this value instead of computing from PMD curve

    Returns
    -------
    double
        π in range [0.01, 0.99], or 0.5 if PMD curve unavailable
    """
    cdef double omega, D1, baseline, pi_raw, pi_clamped

    # Use override if specified
    if pi_override > 0.0:
        return fmax(0.01, fmin(0.99, pi_override))

    # Return neutral if no PMD curve
    if pmd_curve == NULL:
        return 0.5

    # Extract parameters from PMD curve
    omega = <double>pmd_curve.omega
    D1 = <double>pmd_curve.D_5p_noncpg[0]  # Terminal damage rate (position 1)
    baseline = <double>pmd_curve.C_background

    # Fallback to epsilon_baseline if C_background is zero
    if baseline < 1e-6:
        baseline = <double>pmd_curve.epsilon_baseline
    if baseline < 1e-6:
        baseline = 0.01  # Default sequencing error rate

    # Compute π = omega × D1 / (D1 + baseline)
    # This combines:
    # - omega: library-level damage presence (0.81 for MED, ~0.5 for uncertain)
    # - D1/(D1+baseline): signal-to-noise ratio at terminal position
    #
    # For MED sample: π = 0.81 × 0.252 / (0.252 + 0.0225) ≈ 0.74
    # For modern sample: π = 0.5 × 0.02 / (0.02 + 0.02) ≈ 0.25
    if D1 + baseline < 1e-6:
        pi_raw = omega * 0.5  # Neutral if no signal
    else:
        pi_raw = omega * (D1 / (D1 + baseline))

    # Clamp to valid range
    pi_clamped = fmax(0.01, fmin(0.99, pi_raw))

    return pi_clamped


cdef void apply_sample_pi_to_omega(EMState* state, double sample_pi,
                                    uint32_t n_reads) noexcept nogil:
    """Set all omega_ancient values to the sample-level π.

    This effectively replaces per-read omega with a sample-level gate,
    so the prior becomes: logit(p_prior) = logit(γ) + logit(π)
    """
    cdef uint32_t i

    if state == NULL or state.omega_ancient == NULL:
        return

    for i in range(n_reads):
        state.omega_ancient[i] = sample_pi


# =============================================================================
# Posterior-Weighted Coverage Authenticity (Path B)
# =============================================================================

cdef int init_posterior_auth_arrays(MemoryPool* pool) noexcept nogil:
    """Initialize posterior-weighted authenticity arrays in MemoryPool.

    Allocates:
    - authenticity_scores_post: sigmoid(scale * (entropy - gini))
    - norm_entropy_post: normalized spatial entropy from weighted coverage
    - norm_gini_post: normalized Gini from weighted coverage

    Returns 0 on success, -1 on allocation failure.
    """
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t j

    if pool == NULL or n_refs == 0:
        return -1

    pool.authenticity_scores_post = <double*>calloc(n_refs, sizeof(double))
    pool.norm_entropy_post = <double*>calloc(n_refs, sizeof(double))
    pool.norm_gini_post = <double*>calloc(n_refs, sizeof(double))

    if (pool.authenticity_scores_post == NULL or
        pool.norm_entropy_post == NULL or
        pool.norm_gini_post == NULL):
        if pool.authenticity_scores_post != NULL:
            free(pool.authenticity_scores_post)
        if pool.norm_entropy_post != NULL:
            free(pool.norm_entropy_post)
        if pool.norm_gini_post != NULL:
            free(pool.norm_gini_post)
        pool.authenticity_scores_post = NULL
        pool.norm_entropy_post = NULL
        pool.norm_gini_post = NULL
        return -1

    # Initialize to neutral (0.5) - will be updated after first E-step
    for j in range(n_refs):
        pool.authenticity_scores_post[j] = 0.5
        pool.norm_entropy_post[j] = 0.5
        pool.norm_gini_post[j] = 0.5

    return 0


cdef void compute_posterior_weighted_authenticity(
    MemoryPool* pool,
    PosteriorEventBuffer* thread_buffers,
    int n_threads,
    double scale
) noexcept nogil:
    """Compute posterior-weighted coverage authenticity for all references.

    This implements the event sweep algorithm from the plan:
    1. Flatten thread-local event buffers
    2. Bucket events by reference using counting sort
    3. For each reference, sweep to compute weighted coverage histogram
    4. Compute norm_entropy_post, norm_gini_post, authenticity_scores_post

    Args:
        pool: MemoryPool with alignments and output arrays
        thread_buffers: Array of per-thread PosteriorEventBuffer
        n_threads: Number of threads (length of thread_buffers array)
        scale: Sigmoid scale parameter (default 4.0)
    """
    cdef uint32_t n_refs = pool.reference_count
    cdef uint64_t total_events = 0
    cdef uint64_t* ref_counts = NULL
    cdef uint64_t* ref_offsets = NULL
    cdef uint64_t* event_indices = NULL
    cdef float* event_posteriors = NULL
    cdef int tid
    cdef uint64_t i, e, offset, ref_event_count
    cdef uint32_t ref_idx
    cdef uint64_t aln_idx
    cdef PosteriorEvent* evt

    # Count total events across all threads
    for tid in range(n_threads):
        total_events += thread_buffers[tid].len

    if total_events == 0:
        return

    # Allocate flat arrays for bucketing
    ref_counts = <uint64_t*>calloc(n_refs, sizeof(uint64_t))
    ref_offsets = <uint64_t*>calloc(n_refs + 1, sizeof(uint64_t))
    event_indices = <uint64_t*>malloc(total_events * sizeof(uint64_t))
    event_posteriors = <float*>malloc(total_events * sizeof(float))

    if (ref_counts == NULL or ref_offsets == NULL or
        event_indices == NULL or event_posteriors == NULL):
        if ref_counts != NULL: free(ref_counts)
        if ref_offsets != NULL: free(ref_offsets)
        if event_indices != NULL: free(event_indices)
        if event_posteriors != NULL: free(event_posteriors)
        return

    # First pass: count events per reference
    for tid in range(n_threads):
        for i in range(thread_buffers[tid].len):
            evt = &thread_buffers[tid].events[i]
            ref_idx = aln_ref_idx(pool, evt.aln_idx)
            if ref_idx < n_refs:
                ref_counts[ref_idx] += 1

    # Compute prefix sums for offsets
    ref_offsets[0] = 0
    for ref_idx in range(n_refs):
        ref_offsets[ref_idx + 1] = ref_offsets[ref_idx] + ref_counts[ref_idx]

    # Reset counts for use as insertion indices
    memset(ref_counts, 0, n_refs * sizeof(uint64_t))

    # Second pass: place events into buckets
    for tid in range(n_threads):
        for i in range(thread_buffers[tid].len):
            evt = &thread_buffers[tid].events[i]
            ref_idx = aln_ref_idx(pool, evt.aln_idx)
            if ref_idx < n_refs:
                offset = ref_offsets[ref_idx] + ref_counts[ref_idx]
                event_indices[offset] = evt.aln_idx
                event_posteriors[offset] = evt.posterior
                ref_counts[ref_idx] += 1

    # Process each reference
    cdef int64_t ref_length, n_bins
    cdef double* hist_mass = NULL
    cdef IntervalEvent* interval_events = NULL
    cdef int64_t n_interval_events
    cdef int64_t start, end, seg_start, seg_end, bin_idx
    cdef double weight, curr_depth, mass, bin_width, bin_width_inv
    cdef double norm_entropy, norm_gini, raw_score, auth_score

    for ref_idx in range(n_refs):
        ref_event_count = ref_counts[ref_idx]
        if ref_event_count == 0:
            # No events for this reference - leave at neutral
            continue

        ref_length = pool.reference_lengths[ref_idx]
        if ref_length <= 0:
            continue

        # Allocate interval events (2 per alignment event: start and end)
        n_interval_events = ref_event_count * 2
        interval_events = <IntervalEvent*>malloc(n_interval_events * sizeof(IntervalEvent))
        if interval_events == NULL:
            continue

        # Build interval events from alignment events
        offset = ref_offsets[ref_idx]
        for e in range(ref_event_count):
            aln_idx = event_indices[offset + e]
            weight = <double>event_posteriors[offset + e]
            start = <int64_t>aln_position(pool, aln_idx)
            end = start + <int64_t>aln_length(pool, aln_idx)

            # Clamp to reference bounds
            if start < 0:
                start = 0
            if end > ref_length:
                end = ref_length

            interval_events[e * 2].pos = start
            interval_events[e * 2].delta = weight
            interval_events[e * 2 + 1].pos = end
            interval_events[e * 2 + 1].delta = -weight

        # Sort interval events by position
        qsort(interval_events, n_interval_events, sizeof(IntervalEvent), compare_interval_events)

        # Determine number of histogram bins
        n_bins = estimate_histogram_bins_from_length(ref_length)

        # Allocate histogram
        hist_mass = <double*>calloc(n_bins, sizeof(double))
        if hist_mass == NULL:
            free(interval_events)
            continue

        bin_width = <double>ref_length / <double>n_bins
        bin_width_inv = <double>n_bins / <double>ref_length if ref_length > 0 else 0.0

        # Sweep to accumulate histogram mass
        curr_depth = 0.0
        seg_start = 0
        for e in range(n_interval_events):
            seg_end = interval_events[e].pos

            # Accumulate mass for segment [seg_start, seg_end) with curr_depth
            if curr_depth > 0.0 and seg_end > seg_start:
                mass = <double>(seg_end - seg_start) * curr_depth

                # Distribute mass to bins
                # Simple approach: attribute to bin containing segment midpoint
                # (More sophisticated: split across multiple bins)
                bin_idx = <int64_t>((<double>(seg_start + seg_end) / 2.0) * bin_width_inv)
                if bin_idx >= n_bins:
                    bin_idx = n_bins - 1
                if bin_idx < 0:
                    bin_idx = 0
                hist_mass[bin_idx] += mass

            # Update depth
            curr_depth += interval_events[e].delta
            seg_start = seg_end

        # Compute metrics from histogram
        norm_entropy = calculate_weighted_spatial_entropy(hist_mass, n_bins)
        norm_gini = calculate_norm_weighted_gini(hist_mass, n_bins)

        # Compute raw score and apply sigmoid
        raw_score = norm_entropy - norm_gini
        auth_score = 1.0 / (1.0 + exp(-scale * raw_score))

        # Store results
        pool.norm_entropy_post[ref_idx] = norm_entropy
        pool.norm_gini_post[ref_idx] = norm_gini
        pool.authenticity_scores_post[ref_idx] = auth_score

        free(hist_mass)
        free(interval_events)

    # Cleanup
    free(ref_counts)
    free(ref_offsets)
    free(event_indices)
    free(event_posteriors)


cdef void collect_posteriors_for_auth(
    EMState* state,
    MemoryPool* pool,
    PosteriorEventBuffer* thread_buffers,
    int n_threads,
    float posterior_threshold
) noexcept nogil:
    """Collect posterior responsibilities for Path B weighted coverage.

    Iterates over all reads and computes posteriors using current phi_weights.
    When hierarchical is enabled, uses log_L_mix with normalized priors.
    Only stores posteriors above threshold to reduce memory and computation.
    """
    cdef uint32_t n_refs = state.n_refs
    cdef uint32_t n_reads = pool.unique_read_count
    cdef uint32_t read_idx, alignment_count
    cdef uint64_t start_pos, end_pos, aln_idx
    cdef int tid
    cdef uint32_t ref_idx
    cdef double log_phi_j, log_sum, posterior, log_weighted
    cdef double log_max, log_normalizer
    cdef PosteriorEventBuffer* my_buf
    cdef bint use_hierarchical = state.hierarchical_enabled

    # Hierarchical variables
    cdef double gamma_j, log_gamma, log_1m_gamma
    cdef double omega_anc, omega_mod, log_omega_anc, log_omega_mod
    cdef double log_L_anc, log_L_mod, log_L_mix
    cdef double log_odds, log_p_prior_anc, log_p_prior_mod

    # Precompute log(phi) for each reference
    cdef double* log_phi = <double*>malloc(n_refs * sizeof(double))
    if log_phi == NULL:
        return

    for ref_idx in range(n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    # Process reads in parallel
    for read_idx in prange(n_reads, nogil=True, num_threads=n_threads,
                          schedule='dynamic', chunksize=256):
        tid = threadid() if n_threads > 1 else 0
        my_buf = &thread_buffers[tid]

        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Single alignment: posterior = 1.0
        if alignment_count == 1:
            if 1.0 >= posterior_threshold:
                append_posterior(my_buf, start_pos, 1.0)
            continue

        # Get per-read omega for hierarchical mode
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # Multi-alignment: compute posteriors
        # First pass: find max for numerical stability
        log_max = NEG_INF
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                if log_weighted > log_max:
                    log_max = log_weighted

        if log_max <= NEG_INF + 1e10:
            continue

        # Compute normalizer
        log_normalizer = 0.0
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                log_normalizer += exp(log_weighted - log_max)

        log_normalizer = log_max + log(log_normalizer)

        if log_normalizer <= NEG_INF + 1e10:
            continue

        # Second pass: compute posteriors and collect
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                posterior = exp(log_weighted - log_normalizer)
                posterior = fmax(fmin(posterior, 0.999), 1e-12)
                if posterior >= posterior_threshold:
                    append_posterior(my_buf, aln_idx, <float>posterior)

    free(log_phi)


# =============================================================================
# Streaming Posterior-Weighted Coverage Authenticity (Memory-Efficient Path B)
# =============================================================================

cdef void compute_streaming_posterior_auth(
    EMState* state,
    MemoryPool* pool,
    double scale,
    float posterior_threshold
) noexcept nogil:
    """Compute posterior-weighted coverage authenticity using streaming histograms.

    Memory-efficient alternative to event buffering. Uses O(n_refs × n_bins) memory
    instead of O(n_alignments). When hierarchical is enabled, uses log_L_mix with
    normalized priors for posterior computation.
    """
    cdef uint32_t n_refs = state.n_refs
    cdef uint32_t n_reads = pool.unique_read_count
    cdef uint32_t read_idx, ref_idx, alignment_count
    cdef uint64_t start_pos, end_pos, aln_idx
    cdef double posterior, log_phi_j, log_max, log_normalizer, log_weighted
    cdef int64_t aln_start, aln_end
    cdef StreamingCoverageHist* ref_hists = NULL
    cdef double* log_phi = NULL
    cdef uint32_t j
    cdef double norm_ent, norm_gini, auth_score
    cdef int64_t refs_with_data = 0
    cdef bint use_hierarchical = state.hierarchical_enabled

    # Hierarchical variables
    cdef double gamma_j, log_gamma, log_1m_gamma
    cdef double omega_anc, omega_mod, log_omega_anc, log_omega_mod
    cdef double log_L_anc, log_L_mod, log_L_mix
    cdef double log_odds, log_p_prior_anc, log_p_prior_mod

    # Allocate streaming histograms for all references
    ref_hists = <StreamingCoverageHist*>calloc(n_refs, sizeof(StreamingCoverageHist))
    if ref_hists == NULL:
        bf_nogil_logf_notime(b"EM_STREAMING", "ERROR: Failed to allocate streaming histograms")
        return

    # Initialize each histogram with reference length
    for ref_idx in range(n_refs):
        if init_streaming_hist(&ref_hists[ref_idx], pool.reference_lengths[ref_idx]) != 0:
            bf_nogil_logf_notime(b"EM_STREAMING", "ERROR: Failed to init histogram for ref %u", ref_idx)
            # Continue with partial allocation

    # Precompute log(phi) for each reference
    log_phi = <double*>malloc(n_refs * sizeof(double))
    if log_phi == NULL:
        for ref_idx in range(n_refs):
            free_streaming_hist(&ref_hists[ref_idx])
        free(ref_hists)
        return

    for ref_idx in range(n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    # Process all reads and accumulate to histograms
    for read_idx in range(n_reads):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Single alignment: posterior = 1.0
        if alignment_count == 1:
            ref_idx = aln_ref_idx(pool, start_pos)
            if ref_idx < n_refs and ref_hists[ref_idx].deltas != NULL:
                aln_start = <int64_t>aln_position(pool, start_pos)
                aln_end = aln_start + <int64_t>aln_length(pool, start_pos)
                accumulate_alignment_to_hist(&ref_hists[ref_idx], aln_start, aln_end, 1.0)
            continue

        # Get per-read omega for hierarchical mode
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # Multi-alignment: compute posteriors and accumulate
        # First pass: find max for numerical stability
        log_max = NEG_INF
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                if log_weighted > log_max:
                    log_max = log_weighted

        if log_max <= NEG_INF + 1e10:
            continue

        # Compute normalizer
        log_normalizer = 0.0
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                log_normalizer += exp(log_weighted - log_max)

        log_normalizer = log_max + log(log_normalizer)

        if log_normalizer <= NEG_INF + 1e10:
            continue

        # Second pass: compute posteriors and accumulate to histograms
        for aln_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, aln_idx)
            if ref_idx < n_refs:
                if use_hierarchical:
                    gamma_j = state.gamma_values[ref_idx]
                    log_gamma = log(fmax(gamma_j, 1e-10))
                    log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                    log_L_anc = <double>aln_log_L_anc(pool, aln_idx)
                    log_L_mod = <double>aln_log_L_mod(pool, aln_idx)

                    log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                    log_p_prior_anc = logodds_to_log_prob(log_odds)
                    log_p_prior_mod = logodds_to_log_1mp(log_odds)

                    log_L_mix = stable_log_sum_exp(
                        log_p_prior_anc + log_L_anc,
                        log_p_prior_mod + log_L_mod
                    )
                    log_weighted = log_phi[ref_idx] + log_L_mix
                else:
                    log_weighted = log_phi[ref_idx] + <double>aln_score(pool, aln_idx)

                posterior = exp(log_weighted - log_normalizer)
                posterior = fmax(fmin(posterior, 0.999), 1e-12)
                if posterior >= posterior_threshold and ref_hists[ref_idx].deltas != NULL:
                    aln_start = <int64_t>aln_position(pool, aln_idx)
                    aln_end = aln_start + <int64_t>aln_length(pool, aln_idx)
                    accumulate_alignment_to_hist(&ref_hists[ref_idx], aln_start, aln_end, posterior)

    # Finalize histograms and compute metrics for each reference
    for ref_idx in range(n_refs):
        if ref_hists[ref_idx].total_mass < 1e-10:
            # No data for this reference - leave at neutral
            continue

        refs_with_data += 1
        finalize_hist_and_compute_metrics(
            &ref_hists[ref_idx],
            &norm_ent, &norm_gini, &auth_score,
            scale
        )

        pool.norm_entropy_post[ref_idx] = norm_ent
        pool.norm_gini_post[ref_idx] = norm_gini
        pool.authenticity_scores_post[ref_idx] = auth_score

    bf_nogil_logf_notime(
        b"EM_STREAMING",
        "Streaming auth complete: %lld refs with data, memory=%lluKB",
        refs_with_data,
        <unsigned long long>(n_refs * (STREAMING_HIST_BINS + 1) * sizeof(double) / 1024)
    )

    # Cleanup
    free(log_phi)
    for ref_idx in range(n_refs):
        free_streaming_hist(&ref_hists[ref_idx])
    free(ref_hists)


# =============================================================================
# Unified Damage Model Integration
# =============================================================================

cdef void copy_damage_stats_to_context(
    UnifiedDamageContext* ctx,
    RefDamageStats* stats,
    uint32_t n_refs
) noexcept nogil:
    """Copy damage counts from RefDamageStats to UnifiedDamageContext.

    This transfers the pre-accumulated damage counts into the unified model
    for fitting during EM iterations.
    """
    cdef uint32_t r
    cdef int z
    cdef RefDamageCounts* counts
    cdef RefDamageStats* src

    if ctx == NULL or stats == NULL or ctx.ref_counts == NULL:
        return

    for r in range(n_refs):
        counts = &ctx.ref_counts[r]
        src = &stats[r]

        # Copy position-specific counts
        for z in range(DAMAGE_MAX_POSITION):
            counts.k_5p[z] = src.k_5p[z]
            counts.n_5p[z] = src.n_5p[z]
            counts.k_3p[z] = src.k_3p[z]
            counts.n_3p[z] = src.n_3p[z]

        # Copy totals
        counts.total_weight = src.total_weight
        counts.n_alignments = <uint32_t>src.total_weight


cdef void update_authenticity_from_damage(
    UnifiedDamageContext* ctx,
    double* authenticity_scores,
    uint32_t n_refs
) noexcept nogil:
    """Update authenticity scores from unified damage model results.

    After fitting the damage model, copy the authenticity scores to the
    array used by CWRP in the EM M-step.
    """
    cdef uint32_t r

    if ctx == NULL or ctx.ref_params == NULL or authenticity_scores == NULL:
        return

    for r in range(n_refs):
        authenticity_scores[r] = ctx.ref_params[r].authenticity


cdef void accumulate_damage_from_alignments(
    UnifiedDamageContext* ctx,
    MemoryPool* pool,
    bint is_single_stranded
) noexcept nogil:
    """Accumulate quality-weighted damage counts from all alignments.

    This scans all alignments in the pool and accumulates per-position damage
    counts into the context's ref_counts arrays. Called at EM initialization
    to bootstrap the damage model.

    Quality weighting: Each alignment's contribution is weighted by a sigmoid
    of its ZS score (alignment_score), down-weighting low-quality alignments.

    Note: Uses simplified distribution across first 8 positions since Alignment
    stores aggregate counts rather than per-position counts.
    """
    cdef uint64_t n_alns = pool.alignment_count
    cdef uint32_t n_refs = pool.reference_count
    cdef RefDamageCounts* counts
    cdef uint64_t i
    cdef uint32_t ref_idx
    cdef int z
    cdef double avg_ct, avg_ga, avg_c, avg_g
    cdef double alignment_score, quality_weight

    # Quality score threshold and scale for soft weighting
    cdef double score_threshold = -5.0
    cdef double score_scale = 2.0

    if ctx == NULL or ctx.ref_counts == NULL or pool.alignment_cores == NULL:
        return

    # Clear existing counts
    for ref_idx in range(n_refs):
        counts = &ctx.ref_counts[ref_idx]
        for z in range(DAMAGE_MAX_POSITION):
            counts.k_5p[z] = 0.0
            counts.n_5p[z] = 0.0
            counts.k_3p[z] = 0.0
            counts.n_3p[z] = 0.0
        counts.total_weight = 0.0
        counts.n_alignments = 0

    # Accumulate from all alignments with quality weighting
    for i in range(n_alns):
        ref_idx = aln_ref_idx(pool, i)

        if ref_idx >= n_refs:
            continue

        counts = &ctx.ref_counts[ref_idx]

        # Compute quality weight from ZS score
        alignment_score = <double>aln_score(pool, i)
        quality_weight = 1.0 / (1.0 + exp(-(alignment_score - score_threshold) / score_scale))

        # Extract damage counts from alignment
        avg_ct = <double>aln_ct_5p(pool, i) / 8.0 * quality_weight
        avg_c = <double>aln_c_at_5p(pool, i) / 8.0 * quality_weight
        avg_ga = <double>aln_ga_3p(pool, i) / 8.0 * quality_weight
        avg_g = <double>aln_g_at_3p(pool, i) / 8.0 * quality_weight

        # Distribute across first 8 positions (approximation)
        for z in range(8):
            counts.n_5p[z] += avg_c
            counts.k_5p[z] += avg_ct

            if not is_single_stranded:
                counts.n_3p[z] += avg_g
                counts.k_3p[z] += avg_ga

        counts.total_weight += quality_weight
        counts.n_alignments += 1


cdef void accumulate_damage_weighted(
    UnifiedDamageContext* ctx,
    EMState* state,
    MemoryPool* pool,
    EMConfig* config,
    bint is_single_stranded
) noexcept nogil:
    """Accumulate damage counts weighted by EM posteriors and alignment quality.

    This mirrors the E-step logic but accumulates damage counts instead of
    phi_counts. Each alignment's contribution is weighted by:
    1. Posterior probability of belonging to that reference (from EM)
    2. Alignment quality score (ZS score) - higher scores contribute more

    The quality weighting prevents low-quality alignments from corrupting
    damage estimates, which is especially important when mapQ is unreliable
    (e.g., bowtie2 -k 1000).
    """
    cdef uint32_t n_refs = state.n_refs
    cdef uint32_t n_reads = pool.unique_read_count
    cdef uint32_t read_idx, ref_idx, alignment_idx, i
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count
    cdef RefDamageCounts* counts
    cdef int z
    cdef uint64_t pool_idx

    cdef double* log_phi = <double*>malloc(n_refs * sizeof(double))
    if log_phi == NULL:
        return

    # Precompute log(phi_j)
    for ref_idx in range(n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    # Clear existing counts
    for ref_idx in range(n_refs):
        counts = &ctx.ref_counts[ref_idx]
        for z in range(DAMAGE_MAX_POSITION):
            counts.k_5p[z] = 0.0
            counts.n_5p[z] = 0.0
            counts.k_3p[z] = 0.0
            counts.n_3p[z] = 0.0
        counts.total_weight = 0.0
        counts.n_alignments = 0

    # Scratch buffers for posterior computation
    cdef double scratch_log_w[MAX_SCRATCH_SIZE]
    cdef uint32_t scratch_refs[MAX_SCRATCH_SIZE]
    cdef uint64_t scratch_aln_idx[MAX_SCRATCH_SIZE]
    cdef double* heap_log_w = NULL
    cdef uint32_t* heap_refs = NULL
    cdef uint64_t* heap_aln_idx = NULL
    cdef double* log_w_buf
    cdef uint32_t* ref_buf
    cdef uint64_t* aln_idx_buf

    cdef double log_max, log_normalizer, log_weighted, posterior
    cdef double alignment_score, gamma_j, log_gamma, log_1m_gamma
    cdef double omega_anc, omega_mod, log_omega_anc, log_omega_mod
    cdef double log_L_anc, log_L_mod, log_L_mix
    cdef double avg_ct, avg_ga, avg_c, avg_g, weighted_ct, weighted_ga, weighted_c, weighted_g
    cdef double quality_weight, combined_weight  # ZS score quality weighting

    cdef float D_5p = config.D_avg_5p
    cdef float D_3p = config.D_avg_3p
    cdef bint use_hierarchical = state.hierarchical_enabled

    # Quality score threshold and scale for soft weighting
    # Alignments below threshold contribute less to damage estimation
    # threshold = -5.0 corresponds to ~99.3% ANI with typical scoring
    cdef double score_threshold = -5.0
    cdef double score_scale = 2.0  # controls steepness of weighting curve

    # Process each read
    for read_idx in range(n_reads):
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Select buffer
        if alignment_count <= MAX_SCRATCH_SIZE:
            log_w_buf = scratch_log_w
            ref_buf = scratch_refs
            aln_idx_buf = scratch_aln_idx
        else:
            heap_log_w = <double*>malloc(alignment_count * sizeof(double))
            heap_refs = <uint32_t*>malloc(alignment_count * sizeof(uint32_t))
            heap_aln_idx = <uint64_t*>malloc(alignment_count * sizeof(uint64_t))
            if heap_log_w == NULL or heap_refs == NULL or heap_aln_idx == NULL:
                if heap_log_w != NULL: free(heap_log_w)
                if heap_refs != NULL: free(heap_refs)
                if heap_aln_idx != NULL: free(heap_aln_idx)
                continue
            log_w_buf = heap_log_w
            ref_buf = heap_refs
            aln_idx_buf = heap_aln_idx

        # Get per-read PMD prior
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # First pass: compute log weights
        log_max = NEG_INF
        i = 0
        for alignment_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, alignment_idx)
            if ref_idx >= n_refs:
                continue

            ref_buf[i] = ref_idx
            aln_idx_buf[i] = alignment_idx  # Store original pool index
            alignment_score = <double>aln_score(pool, alignment_idx)

            if use_hierarchical:
                gamma_j = state.gamma_values[ref_idx]
                log_gamma = log(fmax(gamma_j, 1e-10))
                log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))
                log_L_anc = <double>aln_log_L_anc(pool, alignment_idx)
                log_L_mod = <double>aln_log_L_mod(pool, alignment_idx)
                log_L_mix = stable_log_sum_exp(
                    log_gamma + log_omega_anc + log_L_anc,
                    log_1m_gamma + log_omega_mod + log_L_mod
                )
                log_weighted = log_phi[ref_idx] + log_L_mix
            else:
                log_weighted = log_phi[ref_idx] + alignment_score

            log_w_buf[i] = log_weighted
            if log_weighted > log_max:
                log_max = log_weighted
            i += 1

        if i == 0:
            if heap_log_w != NULL: free(heap_log_w)
            if heap_refs != NULL: free(heap_refs)
            if heap_aln_idx != NULL: free(heap_aln_idx)
            heap_log_w = NULL
            heap_refs = NULL
            heap_aln_idx = NULL
            continue

        # Compute normalizer
        log_normalizer = 0.0
        for alignment_idx in range(i):
            log_normalizer += exp(log_w_buf[alignment_idx] - log_max)
        log_normalizer = log_max + log(log_normalizer)

        if log_normalizer <= NEG_INF + 1e10:
            if heap_log_w != NULL: free(heap_log_w)
            if heap_refs != NULL: free(heap_refs)
            if heap_aln_idx != NULL: free(heap_aln_idx)
            heap_log_w = NULL
            heap_refs = NULL
            heap_aln_idx = NULL
            continue

        # Second pass: accumulate damage weighted by posterior and quality score
        alignment_idx = 0
        for alignment_idx in range(i):
            ref_idx = ref_buf[alignment_idx]
            posterior = exp(log_w_buf[alignment_idx] - log_normalizer)
            posterior = fmax(fmin(posterior, 0.999), 1e-12)

            # Get alignment damage counts using stored pool index
            pool_idx = aln_idx_buf[alignment_idx]
            avg_ct = <double>aln_ct_5p(pool, pool_idx) / 8.0
            avg_c = <double>aln_c_at_5p(pool, pool_idx) / 8.0
            avg_ga = <double>aln_ga_3p(pool, pool_idx) / 8.0
            avg_g = <double>aln_g_at_3p(pool, pool_idx) / 8.0

            # Compute quality weight from ZS score (alignment_score)
            alignment_score = <double>aln_score(pool, pool_idx)
            quality_weight = 1.0 / (1.0 + exp(-(alignment_score - score_threshold) / score_scale))
            combined_weight = posterior * quality_weight

            # Weight by combined posterior × quality
            weighted_ct = avg_ct * combined_weight
            weighted_c = avg_c * combined_weight
            weighted_ga = avg_ga * combined_weight
            weighted_g = avg_g * combined_weight

            counts = &ctx.ref_counts[ref_idx]
            for z in range(8):
                counts.n_5p[z] += weighted_c
                counts.k_5p[z] += weighted_ct
                if not is_single_stranded:
                    counts.n_3p[z] += weighted_g
                    counts.k_3p[z] += weighted_ga

            counts.total_weight += combined_weight
            counts.n_alignments += 1

        if heap_log_w != NULL: free(heap_log_w)
        if heap_refs != NULL: free(heap_refs)
        if heap_aln_idx != NULL: free(heap_aln_idx)
        heap_log_w = NULL
        heap_refs = NULL
        heap_aln_idx = NULL

    free(log_phi)


# =============================================================================
# State Management
# =============================================================================

cdef EMState* create_em_state(uint32_t n_refs, uint32_t n_reads,
                               EMConfig* config) noexcept nogil:
    """Allocate and initialize EMState structure."""
    cdef EMState* state = <EMState*>malloc(sizeof(EMState))
    if state == NULL:
        return NULL

    state.n_refs = n_refs
    state.n_reads = n_reads

    # Allocate phi weights and counts
    state.phi_weights = <double*>calloc(n_refs, sizeof(double))
    state.phi_counts = <double*>calloc(n_refs, sizeof(double))

    if state.phi_weights == NULL or state.phi_counts == NULL:
        free_em_state(state)
        return NULL

    # Initialize phi to uniform
    cdef double init_weight = 1.0 / <double>n_refs
    cdef uint32_t j
    for j in range(n_refs):
        state.phi_weights[j] = init_weight

    # Hierarchical arrays (if enabled)
    state.hierarchical_enabled = config.hierarchical_enabled
    state.gamma_values = NULL
    state.S_anc = NULL
    state.S_mod = NULL

    if config.hierarchical_enabled:
        state.gamma_values = <double*>calloc(n_refs, sizeof(double))
        state.S_anc = <double*>calloc(n_refs, sizeof(double))
        state.S_mod = <double*>calloc(n_refs, sizeof(double))

        if state.gamma_values == NULL or state.S_anc == NULL or state.S_mod == NULL:
            free_em_state(state)
            return NULL

        # Initialize gamma to 0.5 (uninformative)
        for j in range(n_refs):
            state.gamma_values[j] = 0.5

    # Unknown component
    state.unknown_enabled = config.unknown_enabled
    state.phi_unknown = 0.05 if config.unknown_enabled else 0.0
    state.S_unknown = 0.0

    # PMD priors (allocated but not initialized here - caller must set)
    state.omega_ancient = NULL
    if config.hierarchical_enabled:
        state.omega_ancient = <double*>calloc(n_reads, sizeof(double))
        if state.omega_ancient == NULL:
            free_em_state(state)
            return NULL
        # Default: uninformative prior (0.5)
        for j in range(n_reads):
            state.omega_ancient[j] = 0.5

    # Copy config values
    state.dirichlet_prior = config.dirichlet_prior
    state.gamma_prior = config.gamma_prior
    state.power_rho = config.power_rho
    state.unknown_margin = config.unknown_margin
    state.em_beta = 1.0  # Default temperature

    # Initialize unified damage model if enabled
    state.damage_ctx = NULL
    state.damage_model_enabled = config.unified_damage_enabled
    state.tau_update_interval = config.damage_update_interval if config.damage_update_interval > 0 else 5
    state.tau_current = config.initial_tau if config.initial_tau > 0 else 5.0

    if config.unified_damage_enabled:
        state.damage_ctx = create_damage_context(n_refs, config.is_single_stranded, 0)
        if state.damage_ctx == NULL:
            bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: Failed to create damage context")
            state.damage_model_enabled = False
        elif config.initial_tau > 0:
            set_tau_from_pmd(state.damage_ctx, config.initial_tau, config.fix_tau)

    return state


cdef void free_em_state(EMState* state) noexcept nogil:
    """Free EMState and all its arrays."""
    if state == NULL:
        return

    if state.phi_weights != NULL:
        free(state.phi_weights)
    if state.phi_counts != NULL:
        free(state.phi_counts)
    if state.gamma_values != NULL:
        free(state.gamma_values)
    if state.S_anc != NULL:
        free(state.S_anc)
    if state.S_mod != NULL:
        free(state.S_mod)
    if state.omega_ancient != NULL:
        free(state.omega_ancient)

    # Free unified damage context
    if state.damage_ctx != NULL:
        destroy_damage_context(state.damage_ctx)

    free(state)


cdef void copy_em_state(EMState* dest, EMState* src) noexcept nogil:
    """Deep copy EMState (assumes dest is already allocated with same dimensions)."""
    if dest == NULL or src == NULL:
        return

    memcpy(dest.phi_weights, src.phi_weights, src.n_refs * sizeof(double))
    memcpy(dest.phi_counts, src.phi_counts, src.n_refs * sizeof(double))

    if src.hierarchical_enabled and dest.gamma_values != NULL:
        memcpy(dest.gamma_values, src.gamma_values, src.n_refs * sizeof(double))
        memcpy(dest.S_anc, src.S_anc, src.n_refs * sizeof(double))
        memcpy(dest.S_mod, src.S_mod, src.n_refs * sizeof(double))

    dest.phi_unknown = src.phi_unknown
    dest.S_unknown = src.S_unknown

    # omega_ancient is FIXED, no need to copy


cdef void reset_accumulators(EMState* state) noexcept nogil:
    """Reset E-step accumulators to zero."""
    if state == NULL:
        return

    memset(state.phi_counts, 0, state.n_refs * sizeof(double))

    if state.hierarchical_enabled:
        memset(state.S_anc, 0, state.n_refs * sizeof(double))
        memset(state.S_mod, 0, state.n_refs * sizeof(double))

    state.S_unknown = 0.0


cdef void normalize_phi(EMState* state) noexcept nogil:
    """Normalize phi weights to sum to 1 (including unknown if enabled)."""
    if state == NULL:
        return

    cdef double total = 0.0
    cdef uint32_t j

    for j in range(state.n_refs):
        state.phi_weights[j] = fmax(state.phi_weights[j], 1e-15)
        total += state.phi_weights[j]

    if state.unknown_enabled:
        state.phi_unknown = fmax(state.phi_unknown, 1e-15)
        total += state.phi_unknown

    if total > 0:
        for j in range(state.n_refs):
            state.phi_weights[j] /= total
        if state.unknown_enabled:
            state.phi_unknown /= total


# =============================================================================
# E-Step with optimizations
# =============================================================================


cdef inline void process_read_single(
    EMState* state, MemoryPool* pool, double* log_phi,
    uint32_t read_idx, uint64_t start_pos,
    double* accum_phi, double* accum_S_anc, double* accum_S_mod,
    double log_omega_anc, double log_omega_mod,
    bint use_hierarchical
) noexcept nogil:
    """Fast path for reads with single alignment.

    For single-alignment reads, the posterior is 1.0 (no competition).
    We still compute p_anc using the unified likelihood model with normalized prior.
    """
    cdef uint32_t ref_idx = aln_ref_idx(pool, start_pos)
    cdef double gamma_j, log_gamma, log_1m_gamma
    cdef double log_L_anc, log_L_mod, log_L_mix, log_odds, p_anc
    cdef double log_p_prior_anc, log_p_prior_mod

    if ref_idx >= state.n_refs:
        return

    accum_phi[ref_idx] += 1.0

    if use_hierarchical:
        gamma_j = state.gamma_values[ref_idx]
        log_gamma = log(fmax(gamma_j, 1e-10))
        log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

        # Use sample-fitted log-likelihoods
        log_L_anc = <double>aln_log_L_anc(pool, start_pos)
        log_L_mod = <double>aln_log_L_mod(pool, start_pos)

        # Combined prior: logit(p_prior) = logit(gamma) + logit(omega)
        log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
        log_p_prior_anc = logodds_to_log_prob(log_odds)
        log_p_prior_mod = logodds_to_log_1mp(log_odds)

        # Compute log_L_mix for proper posterior
        log_L_mix = stable_log_sum_exp(
            log_p_prior_anc + log_L_anc,
            log_p_prior_mod + log_L_mod
        )

        # p_anc = prior_anc * L_anc / L_mix
        p_anc = exp(log_p_prior_anc + log_L_anc - log_L_mix)
        p_anc = fmax(fmin(p_anc, 0.999), 0.001)
        accum_S_anc[ref_idx] += p_anc
        accum_S_mod[ref_idx] += 1.0 - p_anc


cdef inline void process_read_multi(
    EMState* state, MemoryPool* pool, double* log_phi, double log_phi_u,
    uint32_t read_idx, uint64_t start_pos, uint64_t end_pos, uint32_t alignment_count,
    double* accum_phi, double* accum_S_anc, double* accum_S_mod, double* accum_S_unknown,
    double log_omega_anc, double log_omega_mod, double unknown_margin,
    bint use_hierarchical, bint use_unknown
) noexcept nogil:
    """Process read with multiple alignments using scratch buffer.

    When hierarchical is enabled, uses log_L_mix (unified model) for both
    reference assignment and p_anc computation. This matches compute_log_likelihood.

    The prior P(ancient) combines gamma (per-ref) and omega (per-read) via log-odds:
        logit(p_prior) = logit(gamma) + logit(omega)
    This is then normalized to proper probabilities for the mixture model.
    """
    cdef uint32_t n_refs = state.n_refs
    cdef uint32_t ref_idx, i, alignment_idx
    cdef double log_normalizer, log_weighted, posterior, log_max, s_max, s_unknown
    cdef double alignment_score, gamma_j, log_gamma, log_1m_gamma
    cdef double log_L_anc, log_L_mod, log_L_mix, log_odds, p_anc, unknown_posterior
    cdef double log_p_prior_anc, log_p_prior_mod

    # Scratch buffers (stack for small, heap for large)
    cdef double scratch_log_w[MAX_SCRATCH_SIZE]
    cdef double scratch_p_anc[MAX_SCRATCH_SIZE]
    cdef uint32_t scratch_refs[MAX_SCRATCH_SIZE]
    cdef double* heap_log_w = NULL
    cdef double* heap_p_anc = NULL
    cdef uint32_t* heap_refs = NULL
    cdef double* log_w_buf
    cdef double* p_anc_buf
    cdef uint32_t* ref_buf

    # Select buffer
    if alignment_count <= MAX_SCRATCH_SIZE:
        log_w_buf = scratch_log_w
        p_anc_buf = scratch_p_anc
        ref_buf = scratch_refs
    else:
        heap_log_w = <double*>malloc(alignment_count * sizeof(double))
        heap_refs = <uint32_t*>malloc(alignment_count * sizeof(uint32_t))
        if use_hierarchical:
            heap_p_anc = <double*>malloc(alignment_count * sizeof(double))
        if heap_log_w == NULL or heap_refs == NULL:
            if heap_log_w != NULL: free(heap_log_w)
            if heap_refs != NULL: free(heap_refs)
            if heap_p_anc != NULL: free(heap_p_anc)
            return
        log_w_buf = heap_log_w
        p_anc_buf = heap_p_anc
        ref_buf = heap_refs

    # First pass: compute log weights and find max
    log_max = NEG_INF
    s_max = NEG_INF
    i = 0

    for alignment_idx in range(start_pos, end_pos):
        ref_idx = aln_ref_idx(pool, alignment_idx)

        # Prefetch next alignment core
        if alignment_idx + 1 < end_pos:
            PREFETCH_READ(&pool.alignment_cores[alignment_idx + 1])

        if ref_idx >= n_refs:
            continue

        ref_buf[i] = ref_idx
        alignment_score = <double>aln_score(pool, alignment_idx)

        if use_hierarchical:
            gamma_j = state.gamma_values[ref_idx]
            log_gamma = log(fmax(gamma_j, 1e-10))
            log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

            # Get sample-fitted log-likelihoods
            log_L_anc = <double>aln_log_L_anc(pool, alignment_idx)
            log_L_mod = <double>aln_log_L_mod(pool, alignment_idx)

            # Combined prior: logit(p_prior) = logit(gamma) + logit(omega)
            # Convert to normalized log-probabilities
            log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
            log_p_prior_anc = logodds_to_log_prob(log_odds)
            log_p_prior_mod = logodds_to_log_1mp(log_odds)

            # Unified model with normalized prior
            log_L_mix = stable_log_sum_exp(
                log_p_prior_anc + log_L_anc,
                log_p_prior_mod + log_L_mod
            )
            log_weighted = log_phi[ref_idx] + log_L_mix

            # Track max log_L_mix for unknown component (same likelihood domain)
            if log_L_mix > s_max:
                s_max = log_L_mix

            # Compute p_anc using posterior = prior * likelihood / marginal
            # log(p_anc) = log_p_prior_anc + log_L_anc - log_L_mix
            p_anc_buf[i] = exp(log_p_prior_anc + log_L_anc - log_L_mix)
            p_anc_buf[i] = fmax(fmin(p_anc_buf[i], 0.999), 0.001)
        else:
            # Non-hierarchical: use raw alignment score
            log_weighted = log_phi[ref_idx] + alignment_score
            if alignment_score > s_max:
                s_max = alignment_score

        log_w_buf[i] = log_weighted
        if log_weighted > log_max:
            log_max = log_weighted
        i += 1

    if i == 0:
        if heap_log_w != NULL: free(heap_log_w)
        if heap_refs != NULL: free(heap_refs)
        if heap_p_anc != NULL: free(heap_p_anc)
        return

    # Compute normalizer with max-subtraction
    log_normalizer = 0.0
    for alignment_idx in range(i):
        log_normalizer += exp(log_w_buf[alignment_idx] - log_max)

    if use_unknown:
        s_unknown = s_max - unknown_margin
        log_weighted = log_phi_u + s_unknown
        log_normalizer += exp(log_weighted - log_max)

    log_normalizer = log_max + log(log_normalizer)

    if log_normalizer <= NEG_INF + 1e10:
        if heap_log_w != NULL: free(heap_log_w)
        if heap_refs != NULL: free(heap_refs)
        if heap_p_anc != NULL: free(heap_p_anc)
        return

    # Second pass: compute posteriors from cached values
    for alignment_idx in range(i):
        ref_idx = ref_buf[alignment_idx]
        posterior = exp(log_w_buf[alignment_idx] - log_normalizer)
        posterior = fmax(fmin(posterior, 0.999), 1e-12)

        PREFETCH_WRITE(&accum_phi[ref_idx])
        accum_phi[ref_idx] += posterior

        if use_hierarchical:
            p_anc = p_anc_buf[alignment_idx]
            accum_S_anc[ref_idx] += posterior * p_anc
            accum_S_mod[ref_idx] += posterior * (1.0 - p_anc)

    if use_unknown:
        unknown_posterior = exp(log_phi_u + s_unknown - log_normalizer)
        unknown_posterior = fmax(fmin(unknown_posterior, 0.999), 1e-12)
        accum_S_unknown[0] += unknown_posterior

    if heap_log_w != NULL: free(heap_log_w)
    if heap_refs != NULL: free(heap_refs)
    if heap_p_anc != NULL: free(heap_p_anc)


cdef void e_step(EMState* state, void* pool_ptr,
                 EMConfig* config) noexcept nogil:
    """E-step computing responsibilities with optimizations for speed.

    Uses unified likelihood model (log_L_mix) when hierarchical is enabled.
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t read_idx, ref_idx
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count, n_refs = state.n_refs
    cdef uint32_t n_reads = pool.unique_read_count

    cdef double unknown_margin = config.unknown_margin
    cdef bint use_unknown = state.unknown_enabled
    cdef bint use_hierarchical = state.hierarchical_enabled

    cdef double* log_phi = <double*>malloc(n_refs * sizeof(double))
    cdef double log_phi_u

    # Per-read omega values for unified model
    cdef double omega_anc, omega_mod, log_omega_anc, log_omega_mod

    # Thread-local accumulators
    cdef int n_threads = config.thread_count if config.thread_count > 0 else 1
    cdef int tid
    cdef double* thread_phi_counts = NULL
    cdef double* thread_S_anc = NULL
    cdef double* thread_S_mod = NULL
    cdef double* thread_S_unknown = NULL
    cdef size_t thread_buf_size
    cdef double* my_phi
    cdef double* my_S_anc
    cdef double* my_S_mod
    cdef double* my_S_unknown

    if log_phi == NULL:
        return

    # Reset accumulators
    reset_accumulators(state)

    # Precompute log(phi_j)
    for ref_idx in range(n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    log_phi_u = log(fmax(state.phi_unknown, 1e-15)) if use_unknown else NEG_INF

    # Allocate thread-local accumulators
    if n_threads > 1:
        thread_buf_size = <size_t>n_threads * n_refs
        thread_phi_counts = <double*>calloc(thread_buf_size, sizeof(double))
        if use_hierarchical:
            thread_S_anc = <double*>calloc(thread_buf_size, sizeof(double))
            thread_S_mod = <double*>calloc(thread_buf_size, sizeof(double))
        if use_unknown:
            thread_S_unknown = <double*>calloc(n_threads, sizeof(double))

        if thread_phi_counts == NULL:
            n_threads = 1

    # Process reads in parallel
    for read_idx in prange(n_reads, nogil=True, num_threads=n_threads,
                           schedule='dynamic', chunksize=256):
        tid = threadid() if n_threads > 1 else 0
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Compute per-read omega for unified model
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # Select accumulator buffers
        if n_threads > 1:
            my_phi = &thread_phi_counts[tid * n_refs]
            my_S_anc = &thread_S_anc[tid * n_refs] if use_hierarchical else NULL
            my_S_mod = &thread_S_mod[tid * n_refs] if use_hierarchical else NULL
            my_S_unknown = &thread_S_unknown[tid] if use_unknown else NULL
        else:
            my_phi = state.phi_counts
            my_S_anc = state.S_anc
            my_S_mod = state.S_mod
            my_S_unknown = &state.S_unknown

        # Fast path for single alignment
        if alignment_count == 1 and not use_unknown:
            process_read_single(
                state, pool, log_phi, read_idx, start_pos,
                my_phi, my_S_anc, my_S_mod,
                log_omega_anc, log_omega_mod, use_hierarchical
            )
        else:
            process_read_multi(
                state, pool, log_phi, log_phi_u,
                read_idx, start_pos, end_pos, alignment_count,
                my_phi, my_S_anc, my_S_mod, my_S_unknown,
                log_omega_anc, log_omega_mod, unknown_margin,
                use_hierarchical, use_unknown
            )

    # Reduce thread-local accumulators
    if n_threads > 1 and thread_phi_counts != NULL:
        for tid in range(n_threads):
            for ref_idx in range(n_refs):
                state.phi_counts[ref_idx] += thread_phi_counts[tid * n_refs + ref_idx]
                if use_hierarchical and thread_S_anc != NULL:
                    state.S_anc[ref_idx] += thread_S_anc[tid * n_refs + ref_idx]
                    state.S_mod[ref_idx] += thread_S_mod[tid * n_refs + ref_idx]
            if use_unknown and thread_S_unknown != NULL:
                state.S_unknown += thread_S_unknown[tid]

        free(thread_phi_counts)
        if thread_S_anc != NULL: free(thread_S_anc)
        if thread_S_mod != NULL: free(thread_S_mod)
        if thread_S_unknown != NULL: free(thread_S_unknown)

    free(log_phi)


# =============================================================================
# Unified M-Step (Standard EM in phi-space, NO power transform)
# =============================================================================

cdef void m_step(EMState* state, EMConfig* config) noexcept nogil:
    """
    Unified M-step updating all parameters using expected counts.

    CRITICAL: NO power transform here. We optimize in phi-space.

    phi_j^new = (E[c_j] + alpha_j) / (sum E[c_j'] + sum alpha_j' + unknown terms)
    gamma_j^new = (S_anc_j + alpha_gamma) / (S_anc_j + S_mod_j + 2*alpha_gamma)

    When CWRP enabled with softmax redistribution:
        alpha_j = alpha_0 * n_refs * softmax(lambda * authenticity_j)
    This keeps sum(alpha) constant while redistributing based on authenticity.
    Otherwise: alpha_j = alpha_0 (uniform Dirichlet)
    """
    cdef double total_count = 0.0
    cdef double alpha_0 = config.dirichlet_prior
    cdef double alpha_gamma = config.gamma_prior
    cdef double alpha_j
    cdef double cwrp_lambda = config.coverage_prior_lambda
    cdef bint use_cwrp = config.coverage_prior_enabled and config.authenticity_scores != NULL
    cdef uint32_t j
    cdef double denom

    # Softmax redistribution variables
    cdef double max_score = -INFINITY
    cdef double sum_exp = 0.0
    cdef double scaled_score
    cdef double total_prior_mass = alpha_0 * <double>state.n_refs

    # For softmax redistribution, first compute normalization constant
    if use_cwrp:
        # Find max for numerical stability
        for j in range(state.n_refs):
            scaled_score = cwrp_lambda * config.authenticity_scores[j]
            if scaled_score > max_score:
                max_score = scaled_score

        # Compute sum of exp(score - max) for softmax denominator
        for j in range(state.n_refs):
            scaled_score = cwrp_lambda * config.authenticity_scores[j]
            sum_exp += exp(scaled_score - max_score)

        # Avoid division by zero
        if sum_exp <= 0:
            sum_exp = 1.0

    # Compute total for normalization (with per-reference priors if CWRP enabled)
    for j in range(state.n_refs):
        if use_cwrp:
            # Softmax redistribution: alpha_j = total_prior_mass * softmax_j
            scaled_score = cwrp_lambda * config.authenticity_scores[j]
            alpha_j = total_prior_mass * exp(scaled_score - max_score) / sum_exp
        else:
            alpha_j = alpha_0
        total_count += state.phi_counts[j] + alpha_j

    if state.unknown_enabled:
        total_count += state.S_unknown + alpha_0

    if total_count <= 0:
        total_count = 1.0

    # Update phi weights (MAP with Dirichlet prior)
    for j in range(state.n_refs):
        if use_cwrp:
            # Softmax redistribution
            scaled_score = cwrp_lambda * config.authenticity_scores[j]
            alpha_j = total_prior_mass * exp(scaled_score - max_score) / sum_exp
        else:
            alpha_j = alpha_0
        state.phi_weights[j] = (state.phi_counts[j] + alpha_j) / total_count

    # Update unknown weight
    if state.unknown_enabled:
        state.phi_unknown = (state.S_unknown + alpha_0) / total_count

    # Update gamma values (hierarchical)
    # Clamp to [eps, 1-eps] for numerical stability when taking logs
    cdef double gamma_eps = 1e-6
    if state.hierarchical_enabled:
        for j in range(state.n_refs):
            denom = state.S_anc[j] + state.S_mod[j] + 2.0 * alpha_gamma
            if denom > 0:
                state.gamma_values[j] = (state.S_anc[j] + alpha_gamma) / denom
                # Clamp for numerical safety
                state.gamma_values[j] = fmax(gamma_eps, fmin(1.0 - gamma_eps, state.gamma_values[j]))
            else:
                state.gamma_values[j] = 0.5  # Uninformative default

    # Ensure proper normalization
    normalize_phi(state)


# =============================================================================
# Log-Likelihood Computation
# =============================================================================

cdef double compute_log_likelihood(EMState* state, void* pool_ptr,
                                   EMConfig* config) noexcept nogil:
    """
    Compute total log-likelihood of current parameter values (parallelized).

    log P(X | phi, gamma) = sum_i log sum_j phi_j × f(x_i | j, gamma_j)
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t read_idx, ref_idx, alignment_idx
    cdef uint64_t start_pos, end_pos
    cdef uint32_t alignment_count
    cdef double total_ll = 0.0
    cdef double log_normalizer, log_weighted
    cdef double alignment_score
    cdef int32_t valid_reads = 0

    # Hierarchical variables
    cdef double log_L_anc, log_L_mod, log_L_mix
    cdef double gamma_j, omega_anc, omega_mod
    cdef double log_gamma, log_1m_gamma, log_omega_anc, log_omega_mod
    cdef double log_odds, log_p_prior_anc, log_p_prior_mod
    cdef bint use_hierarchical = state.hierarchical_enabled

    # Unknown variables
    cdef double log_phi_u, s_max, s_unknown
    cdef double unknown_margin = config.unknown_margin
    cdef bint use_unknown = state.unknown_enabled

    cdef double* log_phi = <double*>malloc(state.n_refs * sizeof(double))

    if log_phi == NULL:
        return NEG_INF

    # Thread-local accumulators
    cdef int32_t num_threads = config.thread_count if config.thread_count > 0 else 1
    cdef double* thread_ll = <double*>calloc(num_threads, sizeof(double))
    cdef int32_t* thread_valid = <int32_t*>calloc(num_threads, sizeof(int32_t))
    cdef int tid

    if thread_ll == NULL or thread_valid == NULL:
        free(log_phi)
        if thread_ll != NULL:
            free(thread_ll)
        if thread_valid != NULL:
            free(thread_valid)
        return NEG_INF

    # Precompute log(phi_j)
    for ref_idx in range(state.n_refs):
        log_phi[ref_idx] = log(fmax(state.phi_weights[ref_idx], 1e-15))

    log_phi_u = log(fmax(state.phi_unknown, 1e-15)) if use_unknown else NEG_INF

    for read_idx in prange(pool.unique_read_count, nogil=True,
                           num_threads=num_threads, schedule='static'):
        tid = threadid()
        alignment_count = pool.read_alignment_counts[read_idx]
        if alignment_count == 0:
            continue

        start_pos = pool.read_alignment_starts[read_idx]
        end_pos = start_pos + alignment_count

        # Get PMD prior
        omega_anc = 0.5
        omega_mod = 0.5
        if use_hierarchical and state.omega_ancient != NULL:
            omega_anc = state.omega_ancient[read_idx]
            omega_mod = 1.0 - omega_anc
        log_omega_anc = log(fmax(omega_anc, 1e-10))
        log_omega_mod = log(fmax(omega_mod, 1e-10))

        # Compute log P(read | phi) and track max for unknown component
        log_normalizer = NEG_INF
        s_max = NEG_INF

        for alignment_idx in range(start_pos, end_pos):
            ref_idx = aln_ref_idx(pool, alignment_idx)
            if ref_idx >= state.n_refs:
                continue

            alignment_score = <double>aln_score(pool, alignment_idx)

            if use_hierarchical:
                gamma_j = state.gamma_values[ref_idx]
                log_gamma = log(fmax(gamma_j, 1e-10))
                log_1m_gamma = log(fmax(1.0 - gamma_j, 1e-10))

                # Use precomputed log-likelihoods
                log_L_anc = <double>aln_log_L_anc(pool, alignment_idx)
                log_L_mod = <double>aln_log_L_mod(pool, alignment_idx)

                # Combined prior: logit(p_prior) = logit(gamma) + logit(omega)
                log_odds = (log_gamma - log_1m_gamma) + (log_omega_anc - log_omega_mod)
                log_p_prior_anc = logodds_to_log_prob(log_odds)
                log_p_prior_mod = logodds_to_log_1mp(log_odds)

                # Unified model with normalized prior
                log_L_mix = stable_log_sum_exp(
                    log_p_prior_anc + log_L_anc,
                    log_p_prior_mod + log_L_mod
                )

                log_weighted = log_phi[ref_idx] + log_L_mix

                # Track max log_L_mix for unknown (same domain as E-step)
                if use_unknown and log_L_mix > s_max:
                    s_max = log_L_mix
            else:
                log_weighted = log_phi[ref_idx] + alignment_score
                # Non-hierarchical: use alignment_score for unknown
                if use_unknown and alignment_score > s_max:
                    s_max = alignment_score

            log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted)

        if use_unknown:
            s_unknown = s_max - unknown_margin
            log_weighted = log_phi_u + s_unknown
            log_normalizer = stable_log_sum_exp(log_normalizer, log_weighted)

        if log_normalizer > NEG_INF + 1e10:
            thread_ll[tid] += log_normalizer
            thread_valid[tid] += 1

    # Reduce thread-local results
    for tid in range(num_threads):
        total_ll += thread_ll[tid]
        valid_reads += thread_valid[tid]

    free(log_phi)
    free(thread_ll)
    free(thread_valid)

    return total_ll / valid_reads if valid_reads > 0 else NEG_INF


# =============================================================================
# SQUAREM Acceleration
# =============================================================================

cdef SQUAREMState* create_squarem_state(uint32_t dimension) noexcept nogil:
    """Allocate SQUAREM working arrays including gamma vectors for hierarchical mode."""
    cdef SQUAREMState* sq = <SQUAREMState*>malloc(sizeof(SQUAREMState))
    if sq == NULL:
        return NULL

    sq.dimension = dimension
    sq.allocated = False
    sq.gamma_allocated = False

    # Initialize all pointers to NULL for safe cleanup
    sq.theta_0 = NULL
    sq.theta_1 = NULL
    sq.theta_2 = NULL
    sq.r_vector = NULL
    sq.v_vector = NULL
    sq.theta_extrapolated = NULL
    sq.gamma_0 = NULL
    sq.gamma_1 = NULL
    sq.gamma_2 = NULL
    sq.r_gamma = NULL
    sq.v_gamma = NULL
    sq.gamma_extrapolated = NULL

    # Allocate phi vectors
    sq.theta_0 = <double*>malloc(dimension * sizeof(double))
    sq.theta_1 = <double*>malloc(dimension * sizeof(double))
    sq.theta_2 = <double*>malloc(dimension * sizeof(double))
    sq.r_vector = <double*>malloc(dimension * sizeof(double))
    sq.v_vector = <double*>malloc(dimension * sizeof(double))
    sq.theta_extrapolated = <double*>malloc(dimension * sizeof(double))

    if (sq.theta_0 == NULL or sq.theta_1 == NULL or sq.theta_2 == NULL or
        sq.r_vector == NULL or sq.v_vector == NULL or sq.theta_extrapolated == NULL):
        free_squarem_state(sq)
        return NULL

    sq.allocated = True

    # Allocate gamma vectors for hierarchical mode
    sq.gamma_0 = <double*>malloc(dimension * sizeof(double))
    sq.gamma_1 = <double*>malloc(dimension * sizeof(double))
    sq.gamma_2 = <double*>malloc(dimension * sizeof(double))
    sq.r_gamma = <double*>malloc(dimension * sizeof(double))
    sq.v_gamma = <double*>malloc(dimension * sizeof(double))
    sq.gamma_extrapolated = <double*>malloc(dimension * sizeof(double))

    if (sq.gamma_0 == NULL or sq.gamma_1 == NULL or sq.gamma_2 == NULL or
        sq.r_gamma == NULL or sq.v_gamma == NULL or sq.gamma_extrapolated == NULL):
        # Gamma allocation failed - continue without gamma extrapolation
        if sq.gamma_0 != NULL: free(sq.gamma_0)
        if sq.gamma_1 != NULL: free(sq.gamma_1)
        if sq.gamma_2 != NULL: free(sq.gamma_2)
        if sq.r_gamma != NULL: free(sq.r_gamma)
        if sq.v_gamma != NULL: free(sq.v_gamma)
        if sq.gamma_extrapolated != NULL: free(sq.gamma_extrapolated)
        sq.gamma_0 = NULL
        sq.gamma_1 = NULL
        sq.gamma_2 = NULL
        sq.r_gamma = NULL
        sq.v_gamma = NULL
        sq.gamma_extrapolated = NULL
        sq.gamma_allocated = False
    else:
        sq.gamma_allocated = True

    return sq


cdef void free_squarem_state(SQUAREMState* sq) noexcept nogil:
    """Free SQUAREM state."""
    if sq == NULL:
        return

    # Free phi vectors
    if sq.theta_0 != NULL: free(sq.theta_0)
    if sq.theta_1 != NULL: free(sq.theta_1)
    if sq.theta_2 != NULL: free(sq.theta_2)
    if sq.r_vector != NULL: free(sq.r_vector)
    if sq.v_vector != NULL: free(sq.v_vector)
    if sq.theta_extrapolated != NULL: free(sq.theta_extrapolated)

    # Free gamma vectors
    if sq.gamma_0 != NULL: free(sq.gamma_0)
    if sq.gamma_1 != NULL: free(sq.gamma_1)
    if sq.gamma_2 != NULL: free(sq.gamma_2)
    if sq.r_gamma != NULL: free(sq.r_gamma)
    if sq.v_gamma != NULL: free(sq.v_gamma)
    if sq.gamma_extrapolated != NULL: free(sq.gamma_extrapolated)

    free(sq)


cdef void copy_phi_to_array(EMState* state, double* arr) noexcept nogil:
    """Copy phi weights to array."""
    memcpy(arr, state.phi_weights, state.n_refs * sizeof(double))


cdef void copy_array_to_phi(double* arr, EMState* state) noexcept nogil:
    """Copy array to phi weights."""
    memcpy(state.phi_weights, arr, state.n_refs * sizeof(double))


cdef void copy_gamma_to_array(EMState* state, double* arr) noexcept nogil:
    """Copy gamma values to array."""
    if state.gamma_values != NULL and arr != NULL:
        memcpy(arr, state.gamma_values, state.n_refs * sizeof(double))


cdef void copy_array_to_gamma(double* arr, EMState* state) noexcept nogil:
    """Copy array to gamma values."""
    if arr != NULL and state.gamma_values != NULL:
        memcpy(state.gamma_values, arr, state.n_refs * sizeof(double))


cdef bint squarem_step(EMState* state, void* pool_ptr, EMConfig* config,
                       SQUAREMState* sq, double* ll_out) noexcept nogil:
    """
    SQUAREM acceleration with globalization (backtracking line search).

    Per Varadhan & Roland (2008):
    - S1: α = -||r|| / ||v||          (basic)
    - S2: α = -||r||² / (r·v)         (squared norm ratio)
    - S3: α = -(r·v) / ||v||²         (recommended, better stability)

    α interpretation:
    - α = -1: equivalent to theta_2 (two EM steps)
    - α < -1: extrapolate beyond theta_2 (accelerate)
    - -1 < α < 0: interpolate between theta_1 and theta_2 (conservative)

    Globalization backtracks toward α=-1 when extrapolation overshoots.

    Returns True if accelerated step accepted.
    """
    cdef uint32_t i, n = state.n_refs
    cdef double r_norm_sq = 0.0, v_norm_sq = 0.0, r_dot_v = 0.0
    cdef double r_temp, v_temp
    cdef double alpha, alpha_try
    cdef double ll_2, ll_sq
    cdef double total
    cdef int backtrack, max_backtracks
    cdef double rel_tol, bf

    # Save theta_2 auxiliary state for proper fallback restoration
    cdef double gamma_2_saved_phi_unknown = 0.0
    cdef double* gamma_2_saved = NULL
    cdef bint need_aux_restore = state.hierarchical_enabled or state.unknown_enabled

    if need_aux_restore and state.hierarchical_enabled and state.gamma_values != NULL:
        gamma_2_saved = <double*>malloc(n * sizeof(double))

    # Save theta_0 and gamma_0 (don't compute LL yet - expensive)
    copy_phi_to_array(state, sq.theta_0)
    if state.hierarchical_enabled and sq.gamma_allocated:
        copy_gamma_to_array(state, sq.gamma_0)

    # theta_1 = M(theta_0)
    e_step(state, pool_ptr, config)
    m_step(state, config)
    copy_phi_to_array(state, sq.theta_1)
    if state.hierarchical_enabled and sq.gamma_allocated:
        copy_gamma_to_array(state, sq.gamma_1)

    # theta_2 = M(theta_1)
    e_step(state, pool_ptr, config)
    m_step(state, config)
    copy_phi_to_array(state, sq.theta_2)
    if state.hierarchical_enabled and sq.gamma_allocated:
        copy_gamma_to_array(state, sq.gamma_2)

    # Save auxiliary state at theta_2 for proper fallback
    if state.unknown_enabled:
        gamma_2_saved_phi_unknown = state.phi_unknown
    if gamma_2_saved != NULL:
        memcpy(gamma_2_saved, state.gamma_values, n * sizeof(double))

    ll_2 = compute_log_likelihood(state, pool_ptr, config)

    # Compute r = theta_1 - theta_0, v = (theta_2 - theta_1) - r
    # Also compute squared norms and dot product for step schemes
    for i in range(n):
        r_temp = sq.theta_1[i] - sq.theta_0[i]
        sq.r_vector[i] = r_temp
        r_norm_sq += r_temp * r_temp

        v_temp = (sq.theta_2[i] - sq.theta_1[i]) - r_temp
        sq.v_vector[i] = v_temp
        v_norm_sq += v_temp * v_temp
        r_dot_v += r_temp * v_temp

    # Compute r_gamma = gamma_1 - gamma_0, v_gamma = (gamma_2 - gamma_1) - r_gamma
    # Include gamma in norm calculations for joint extrapolation
    if state.hierarchical_enabled and sq.gamma_allocated:
        for i in range(n):
            r_temp = sq.gamma_1[i] - sq.gamma_0[i]
            sq.r_gamma[i] = r_temp
            r_norm_sq += r_temp * r_temp

            v_temp = (sq.gamma_2[i] - sq.gamma_1[i]) - r_temp
            sq.v_gamma[i] = v_temp
            v_norm_sq += v_temp * v_temp
            r_dot_v += r_temp * v_temp

    # Compute step length alpha based on steplength_scheme
    # S1: α = -||r|| / ||v||
    # S2: α = -||r||² / (r·v)
    # S3: α = -(r·v) / ||v||²  (recommended)
    if config.steplength_scheme == 1:
        # S1: basic ratio of norms
        if v_norm_sq > 1e-30:
            alpha = -libc_sqrt(r_norm_sq / v_norm_sq)
        else:
            if gamma_2_saved != NULL:
                free(gamma_2_saved)
            ll_out[0] = ll_2
            return False
    elif config.steplength_scheme == 2:
        # S2: squared norm over dot product
        if fabs(r_dot_v) > 1e-30:
            alpha = -r_norm_sq / r_dot_v
        else:
            if gamma_2_saved != NULL:
                free(gamma_2_saved)
            ll_out[0] = ll_2
            return False
    else:
        # S3 (default, recommended): dot product over squared norm
        if v_norm_sq > 1e-30:
            alpha = -r_dot_v / v_norm_sq
        else:
            if gamma_2_saved != NULL:
                free(gamma_2_saved)
            ll_out[0] = ll_2
            return False

    # Clamp alpha to reasonable range
    # When α > -1: S3 doesn't suggest extrapolation, fall back to theta_2
    # When α < -300: prevent extreme extrapolation
    if alpha > -1.0:
        # S3 doesn't suggest going beyond theta_2, just accept theta_2
        if gamma_2_saved != NULL:
            free(gamma_2_saved)
        ll_out[0] = ll_2
        return False  # Accept theta_2 (already in state after the 2 EM steps)
    elif alpha < -300.0:
        alpha = -300.0  # Prevent extreme extrapolation

    # At this point α is in [-300, -1), meaning we're extrapolating beyond theta_2

    # Relative tolerance for acceptance
    rel_tol = 1e-6 * fabs(ll_2) if ll_2 != 0.0 else 1e-8

    # Globalization: backtracking line search toward alpha=-1
    if config.enable_globalization:
        max_backtracks = config.max_backtrack_steps if config.max_backtrack_steps > 0 else 4
        bf = config.backtrack_factor if config.backtrack_factor > 0 else 0.5
    else:
        max_backtracks = 1  # Just try once
        bf = 0.5

    for backtrack in range(max_backtracks):
        if backtrack == 0:
            alpha_try = alpha
        else:
            # Backtrack: move alpha toward -1 (theta_2 equivalent)
            # New alpha = -1 + bf^backtrack * (alpha - (-1))
            alpha_try = -1.0 + (bf ** backtrack) * (alpha + 1.0)

            # If we've backtracked to essentially theta_2, stop
            if alpha_try > -1.01:
                break

        # Reset phi_unknown to theta_2 value before each trial
        if state.unknown_enabled:
            state.phi_unknown = gamma_2_saved_phi_unknown

        # Extrapolate phi: theta_sq = theta_0 - 2*alpha*r + alpha^2*v
        total = 0.0
        for i in range(n):
            sq.theta_extrapolated[i] = (sq.theta_0[i]
                                        - 2.0 * alpha_try * sq.r_vector[i]
                                        + alpha_try * alpha_try * sq.v_vector[i])
            sq.theta_extrapolated[i] = fmax(sq.theta_extrapolated[i], 1e-15)
            total += sq.theta_extrapolated[i]

        # Normalize phi to simplex
        if total > 0:
            for i in range(n):
                sq.theta_extrapolated[i] /= total

        # Extrapolate gamma: gamma_sq = gamma_0 - 2*alpha*r_gamma + alpha^2*v_gamma
        # Clamp to [0.001, 0.999] to maintain valid probabilities
        if state.hierarchical_enabled and sq.gamma_allocated:
            for i in range(n):
                sq.gamma_extrapolated[i] = (sq.gamma_0[i]
                                            - 2.0 * alpha_try * sq.r_gamma[i]
                                            + alpha_try * alpha_try * sq.v_gamma[i])
                sq.gamma_extrapolated[i] = fmax(0.001, fmin(0.999, sq.gamma_extrapolated[i]))

        # Apply extrapolated phi AND gamma, then stabilize with one EM step
        # Per Varadhan & Roland (2008): stabilization improves robustness
        copy_array_to_phi(sq.theta_extrapolated, state)
        normalize_phi(state)
        if state.hierarchical_enabled and sq.gamma_allocated:
            copy_array_to_gamma(sq.gamma_extrapolated, state)
        e_step(state, pool_ptr, config)
        m_step(state, config)

        # Evaluate stabilized point
        ll_sq = compute_log_likelihood(state, pool_ptr, config)

        # Check if extrapolation improved over theta_2
        if ll_sq >= ll_2 - rel_tol:
            ll_out[0] = ll_sq
            if gamma_2_saved != NULL:
                free(gamma_2_saved)
            if backtrack == 0:
                bf_nogil_logf_notime(b"EM_UNIFIED",
                    "SQUAREM S%d accepted: alpha=%.2f LL=%.6f (gain=%.2e)",
                    config.steplength_scheme, alpha_try, ll_sq, ll_sq - ll_2)
            else:
                bf_nogil_logf_notime(b"EM_UNIFIED",
                    "SQUAREM S%d accepted after %d backtracks: alpha=%.2f LL=%.6f",
                    config.steplength_scheme, backtrack, alpha_try, ll_sq)
            return True

    # All backtracking attempts failed - fall back to theta_2
    # CRITICAL: Restore full state (phi, gamma, phi_unknown) to theta_2 values
    copy_array_to_phi(sq.theta_2, state)
    normalize_phi(state)

    # Restore auxiliary state that was saved at theta_2
    if state.unknown_enabled:
        state.phi_unknown = gamma_2_saved_phi_unknown
    if gamma_2_saved != NULL:
        memcpy(state.gamma_values, gamma_2_saved, n * sizeof(double))
        free(gamma_2_saved)

    bf_nogil_logf_notime(b"EM_UNIFIED",
        "SQUAREM S%d rejected: alpha=%.2f ll_sq=%.6f < ll_2=%.6f (diff=%.2e)",
        config.steplength_scheme, alpha, ll_sq, ll_2, ll_sq - ll_2)
    ll_out[0] = ll_2
    return False


# =============================================================================
# Output Transform (POST-PROCESSING ONLY)
# =============================================================================

cdef void transform_to_output(EMState* state, double* output_pi,
                              double power_rho) noexcept nogil:
    """
    Transform phi-space weights to reported pi weights.

    pi_j = phi_j^rho / sum_j' phi_j'^rho

    This is ONLY for reporting, not used in EM iteration.

    Power transform interpretation (temperature analogy):
    - rho < 1: FLATTENS distribution (counteracts rich-get-richer bias)
      - e.g., rho=0.7 makes concentrated distributions more uniform
    - rho = 1: No change (standard EM output)
    - rho > 1: SHARPENS distribution (amplifies differences)

    Note: Using phi^rho (NOT phi^(1/rho)) for flattening with rho < 1.
    This is analogous to temperature T = 1/rho > 1 which smooths distributions.
    """
    cdef double total = 0.0
    cdef uint32_t j

    for j in range(state.n_refs):
        # phi^rho: rho < 1 flattens, rho > 1 sharpens
        output_pi[j] = libc_pow(state.phi_weights[j], power_rho)
        total += output_pi[j]

    if total > 0:
        for j in range(state.n_refs):
            output_pi[j] /= total


# =============================================================================
# Convergence Tracking
# =============================================================================

cdef ConvergenceState* create_convergence_state(int32_t history_length) noexcept nogil:
    """Create convergence tracking state."""
    cdef ConvergenceState* conv = <ConvergenceState*>malloc(sizeof(ConvergenceState))
    if conv == NULL:
        return NULL

    conv.history_length = history_length
    conv.ll_history = <double*>calloc(history_length, sizeof(double))
    conv.param_history = <double*>calloc(history_length, sizeof(double))

    if conv.ll_history == NULL or conv.param_history == NULL:
        free_convergence_state(conv)
        return NULL

    conv.current_index = 0
    conv.filled_count = 0
    conv.current_ll = NEG_INF
    conv.prev_ll = NEG_INF

    return conv


cdef void free_convergence_state(ConvergenceState* conv) noexcept nogil:
    """Free convergence state."""
    if conv == NULL:
        return

    if conv.ll_history != NULL:
        free(conv.ll_history)
    if conv.param_history != NULL:
        free(conv.param_history)

    free(conv)


cdef void update_convergence_history(ConvergenceState* conv,
                                     double ll_change, double param_change) noexcept nogil:
    """Update convergence history with new values."""
    if conv == NULL:
        return

    conv.ll_history[conv.current_index] = ll_change
    conv.param_history[conv.current_index] = param_change

    conv.current_index = (conv.current_index + 1) % conv.history_length
    if conv.filled_count < conv.history_length:
        conv.filled_count += 1


cdef double compute_mad(double* values, int32_t n) noexcept nogil:
    """Compute Median Absolute Deviation."""
    if n <= 0:
        return 0.0

    # For small n, just use mean absolute deviation
    cdef double mean = 0.0
    cdef double mad = 0.0
    cdef int32_t i

    for i in range(n):
        mean += values[i]
    mean /= n

    for i in range(n):
        mad += fabs(values[i] - mean)
    mad /= n

    return mad * 1.4826  # Scale factor for consistency with standard deviation


cdef bint check_convergence_mad(ConvergenceState* conv,
                                double base_tolerance) noexcept nogil:
    """
    Check convergence using MAD-based criterion.

    Converged when:
    - At least 5 iterations completed (to allow hierarchical gamma to stabilize)
    - Recent LL change < base_tolerance (strict criterion for actual convergence)
    - AND Recent param change < sqrt(base_tolerance) (looser for params)
    """
    # Require minimum 5 iterations to allow gamma values to stabilize
    if conv == NULL or conv.filled_count < 5:
        return False

    # Check if recent changes are within tolerance (use actual tolerance, not MAD-scaled)
    cdef double recent_ll_change = conv.ll_history[(conv.current_index - 1 + conv.history_length) % conv.history_length]
    cdef double recent_param_change = conv.param_history[(conv.current_index - 1 + conv.history_length) % conv.history_length]

    # Strict convergence: actual changes must be small
    cdef double param_tol = libc_sqrt(base_tolerance)  # More lenient for params

    return (recent_ll_change <= base_tolerance) and (recent_param_change <= param_tol)


cdef double compute_param_change_norm(EMState* current, EMState* prev) noexcept nogil:
    """Compute L2 norm of parameter change."""
    cdef double norm_sq = 0.0
    cdef double diff
    cdef uint32_t j

    for j in range(current.n_refs):
        diff = current.phi_weights[j] - prev.phi_weights[j]
        norm_sq += diff * diff

    return libc_sqrt(norm_sq)


# =============================================================================
# Main Entry Point
# =============================================================================

cdef int run_em(void* pool_ptr, EMConfig* config,
                double* output_weights) noexcept nogil:
    """
    Run EM algorithm.

    Parameters
    ----------
    pool_ptr : void*
        Pointer to MemoryPool
    config : EMConfig*
        Algorithm configuration
    output_weights : double*
        Output array for final weights [n_refs]

    Returns
    -------
    int
        Number of iterations run, or -1 on error
    """
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t n_reads = pool.unique_read_count

    # Create state
    cdef EMState* state = create_em_state(n_refs, n_reads, config)
    if state == NULL:
        bf_nogil_logf_notime(b"EM_UNIFIED", "ERROR: Failed to allocate EMState")
        return -1

    cdef EMState* prev_state = create_em_state(n_refs, n_reads, config)
    if prev_state == NULL:
        free_em_state(state)
        bf_nogil_logf_notime(b"EM_UNIFIED", "ERROR: Failed to allocate prev_state")
        return -1

    # Create convergence tracker
    cdef ConvergenceState* conv = create_convergence_state(config.history_length)
    if conv == NULL:
        free_em_state(state)
        free_em_state(prev_state)
        return -1

    # Create SQUAREM state if enabled
    cdef SQUAREMState* sq = NULL

    # Debug variables for gamma analysis
    cdef double min_gamma_debug = 1.0
    cdef uint32_t min_idx_debug = 0
    cdef double s_tot_debug
    cdef double total_S_anc = 0.0
    cdef double total_S_mod = 0.0
    cdef double max_gamma = 0.0
    cdef double min_gamma = 1.0
    cdef uint32_t j_debug

    if config.squarem_enabled:
        sq = create_squarem_state(n_refs)
        if sq == NULL:
            bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: SQUAREM allocation failed, using standard EM")

    # Compute sample-level P(ancient) from PMD curve
    cdef PMDCurve* pmd_curve = NULL
    cdef double sample_pi = 0.5  # Default to neutral
    cdef double pmd_omega = 0.5
    cdef double pmd_D1 = 0.0
    cdef double pmd_baseline = 0.01

    if pool.pmd_curve_ptr != NULL:
        pmd_curve = <PMDCurve*>pool.pmd_curve_ptr
        sample_pi = compute_sample_pi(pmd_curve, config.sample_pi_override)
        pmd_omega = <double>pmd_curve.omega
        pmd_D1 = <double>pmd_curve.D_5p_noncpg[0]
        pmd_baseline = <double>pmd_curve.C_background
        if pmd_baseline < 1e-6:
            pmd_baseline = <double>pmd_curve.epsilon_baseline
    elif config.sample_pi_override > 0.0:
        sample_pi = fmax(0.01, fmin(0.99, config.sample_pi_override))

    # Store computed sample_pi in config for SQUAREM gating
    config.sample_pi_computed = sample_pi

    # Apply sample π to omega_ancient array (gates the hierarchical damage model)
    if config.hierarchical_enabled:
        apply_sample_pi_to_omega(state, sample_pi, n_reads)
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "Sample-level gate: pi=%.4f (omega=%.3f D1=%.4f baseline=%.4f)",
            sample_pi, pmd_omega, pmd_D1, pmd_baseline
        )

    # Warn if SQUAREM will be disabled due to low damage
    if config.squarem_enabled and sample_pi < config.squarem_pi_threshold:
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "SQUAREM disabled: pi=%.4f < threshold=%.2f (low damage causes phi-gamma desync)",
            sample_pi, config.squarem_pi_threshold
        )

    # Initialize unified damage model if enabled
    cdef bint damage_model_enabled = state.damage_model_enabled
    cdef RefDamageStats* ref_damage_stats = NULL
    cdef int32_t damage_update_interval = state.tau_update_interval
    cdef bint is_single_stranded = config.is_single_stranded

    if damage_model_enabled and state.damage_ctx != NULL:
        # Try to get pre-computed damage stats from pool first
        if pool.ref_damage_stats != NULL:
            ref_damage_stats = <RefDamageStats*>pool.ref_damage_stats
            copy_damage_stats_to_context(state.damage_ctx, ref_damage_stats, n_refs)
            bf_nogil_logf_notime(b"EM_UNIFIED", "Using pre-computed damage stats from pool")
        else:
            # Accumulate damage counts directly from alignments
            bf_nogil_logf_notime(b"EM_UNIFIED", "Accumulating damage stats from %llu alignments...",
                                 pool.alignment_count)
            accumulate_damage_from_alignments(state.damage_ctx, pool, is_single_stranded)

        # Initialize the model parameters
        estimate_baselines(state.damage_ctx)
        initialize_parameters(state.damage_ctx)

        # Initial fit
        fit_damage_model(state.damage_ctx)
        compute_outputs(state.damage_ctx)

        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "Unified damage model initialized: tau=%.2f mu_5p=%.4f mu_3p=%.4f update_interval=%d",
            state.damage_ctx.sample.tau,
            state.damage_ctx.sample.mu_5p,
            state.damage_ctx.sample.mu_3p,
            damage_update_interval
        )

    # Iterative ancientness field setup
    cdef bint iterative_auth_enabled = config.iterative_auth
    cdef int32_t auth_update_interval = config.auth_update_interval
    cdef ReferenceStats* ref_stats_for_auth = NULL

    if iterative_auth_enabled and config.coverage_prior_enabled:
        # Initialize ancientness arrays in pool
        if init_ancientness_arrays(pool) != 0:
            bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: Failed to init ancientness arrays, disabling iterative auth")
            iterative_auth_enabled = False
        else:
            # Allocate ReferenceStats for coverage computation
            ref_stats_for_auth = <ReferenceStats*>calloc(n_refs, sizeof(ReferenceStats))
            if ref_stats_for_auth == NULL:
                bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: Failed to alloc ref_stats, disabling iterative auth")
                iterative_auth_enabled = False
            else:
                bf_nogil_logf_notime(
                    b"EM_UNIFIED",
                    "Iterative ancientness enabled: update_interval=%d damage_weight=%.2f low_cov_floor=%d",
                    auth_update_interval, config.damage_weight, config.low_cov_floor_reads
                )

    # Posterior-Weighted Coverage Authenticity (Path B) setup - using streaming histograms
    cdef bint auth_post_enabled = config.auth_post_enabled and config.coverage_prior_enabled
    cdef int32_t auth_update_interval_post = config.auth_update_interval_post
    cdef double auth_scale_post = config.auth_scale_post
    cdef int32_t auth_lambda_ramp_iters = config.auth_lambda_ramp_iters
    cdef double effective_lambda

    if auth_post_enabled:
        # Initialize posterior auth arrays in pool (streaming approach - no event buffers needed)
        if init_posterior_auth_arrays(pool) != 0:
            bf_nogil_logf_notime(b"EM_UNIFIED", "WARNING: Failed to init posterior auth arrays, disabling Path B")
            auth_post_enabled = False
        else:
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "Posterior auth (Path B) enabled (streaming): update_interval=%d scale=%.2f lambda_ramp=%d bins=%d",
                auth_update_interval_post, auth_scale_post, auth_lambda_ramp_iters, STREAMING_HIST_BINS
            )

    bf_nogil_logf_notime(
        b"EM_UNIFIED",
        "Starting: refs=%u reads=%u max_iter=%d tol=%.2e hierarchical=%s unknown=%s squarem=%s rho=%.2f",
        n_refs, n_reads, config.max_iterations, config.convergence_tolerance,
        b"true" if config.hierarchical_enabled else b"false",
        b"true" if config.unknown_enabled else b"false",
        b"true" if (config.squarem_enabled and sq != NULL) else b"false",
        config.power_rho
    )

    # Debug: Check damage_llr distribution
    cdef double llr_sum = 0.0
    cdef double llr_min = 1e20
    cdef double llr_max = -1e20
    cdef int64_t llr_nonzero = 0
    cdef int64_t aln_idx
    cdef float cur_llr

    if config.hierarchical_enabled:
        for aln_idx in range(pool.alignment_count):
            cur_llr = aln_damage_llr(pool, aln_idx)
            llr_sum += cur_llr
            if cur_llr != 0.0:
                llr_nonzero += 1
            if cur_llr < llr_min:
                llr_min = cur_llr
            if cur_llr > llr_max:
                llr_max = cur_llr

        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "damage_llr stats: count=%lld nonzero=%lld (%.1f%%) min=%.4f max=%.4f mean=%.4f",
            <long long>pool.alignment_count, <long long>llr_nonzero,
            100.0 * <double>llr_nonzero / <double>pool.alignment_count if pool.alignment_count > 0 else 0.0,
            llr_min, llr_max, llr_sum / <double>pool.alignment_count if pool.alignment_count > 0 else 0.0
        )

    cdef int iteration
    cdef double current_ll, prev_ll = NEG_INF
    cdef double ll_change, param_change
    cdef double ll_new
    cdef bint converged = False
    cdef bint use_squarem

    for iteration in range(config.max_iterations):
        # Save previous state
        copy_em_state(prev_state, state)
        prev_ll = conv.current_ll if conv.current_ll > NEG_INF + 1e10 else compute_log_likelihood(state, pool_ptr, config)

        # Decide whether to use SQUAREM
        # Disable SQUAREM for low-damage samples (pi < threshold) to avoid phi-gamma
        # desynchronization that causes repeated step rejections
        use_squarem = (config.squarem_enabled and sq != NULL and
                       iteration >= config.squarem_start_iter and
                       config.sample_pi_computed >= config.squarem_pi_threshold)

        if use_squarem:
            squarem_step(state, pool_ptr, config, sq, &current_ll)
        else:
            # Standard E-step and M-step
            e_step(state, pool_ptr, config)

            # Path B: Posterior-weighted coverage authenticity update (streaming approach)
            if auth_post_enabled and (iteration + 1) % auth_update_interval_post == 0:
                # Compute weighted coverage authenticity using streaming histograms
                # Memory: O(n_refs × n_bins) instead of O(n_alignments)
                compute_streaming_posterior_auth(
                    state, pool, auth_scale_post, 0.05
                )

                # Lambda ramping: scale CWRP weight from 0 to target over ramp iterations
                if auth_lambda_ramp_iters > 0 and iteration < auth_lambda_ramp_iters:
                    effective_lambda = config.coverage_prior_lambda * (<double>(iteration + 1) / <double>auth_lambda_ramp_iters)
                else:
                    effective_lambda = config.coverage_prior_lambda

                # Use posterior-weighted authenticity for CWRP
                config.authenticity_scores = pool.authenticity_scores_post
                config.coverage_prior_lambda = effective_lambda

                if iteration % 10 == 0 or iteration < 5:
                    bf_nogil_logf_notime(
                        b"EM_UNIFIED",
                        "Path B auth update iter %d: lambda=%.3f (ramp %d/%d)",
                        iteration, effective_lambda, iteration + 1, auth_lambda_ramp_iters
                    )

            m_step(state, config)
            current_ll = compute_log_likelihood(state, pool_ptr, config)

        # Compute changes
        ll_change = fabs(current_ll - prev_ll) / fmax(fabs(current_ll), 1.0)
        param_change = compute_param_change_norm(state, prev_state)

        # Update convergence history
        update_convergence_history(conv, ll_change, param_change)
        conv.current_ll = current_ll
        conv.prev_ll = prev_ll

        # Log progress
        if iteration % 10 == 0 or iteration < 5:
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "Iter %d: LL=%.6f dLL=%.2e ||dphi||=%.2e",
                iteration, current_ll, ll_change, param_change
            )

        # Iterative ancientness field update
        if iterative_auth_enabled and (iteration + 1) % auth_update_interval == 0:
            # Recompute coverage metrics using current posterior weights
            # Use batched version for memory efficiency with large datasets
            calculate_reference_coverage_batched(pool, ref_stats_for_auth)

            # Update ancientness field from posterior-weighted features
            update_ancientness_field(
                pool, state.phi_weights, ref_stats_for_auth,
                config.damage_weight, config.low_cov_floor_reads, config.low_cov_shrink_tau
            )

            # Sync updated authenticity scores to config for next M-step
            config.authenticity_scores = pool.authenticity_scores

            if iteration % 10 == 0:
                bf_nogil_logf_notime(
                    b"EM_UNIFIED",
                    "Updated ancientness field at iter %d",
                    iteration
                )

        # Unified damage model update (periodic)
        if damage_model_enabled and state.damage_ctx != NULL and (iteration + 1) % damage_update_interval == 0:
            # Re-accumulate damage counts weighted by current posteriors
            # This makes the damage model responsive to changing reference assignments
            accumulate_damage_weighted(state.damage_ctx, state, pool, config, is_single_stranded)

            # Re-estimate baselines and parameters with new weighted counts
            estimate_baselines(state.damage_ctx)
            initialize_parameters(state.damage_ctx)

            # Re-fit the damage model
            fit_damage_model(state.damage_ctx)
            compute_outputs(state.damage_ctx)

            # Update tau in state for tracking
            state.tau_current = state.damage_ctx.sample.tau

            # Update authenticity_scores for CWRP if enabled
            if config.coverage_prior_enabled and pool.authenticity_scores != NULL:
                update_authenticity_from_damage(state.damage_ctx, pool.authenticity_scores, n_refs)
                config.authenticity_scores = pool.authenticity_scores

            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "Damage update iter %d: tau=%.2f mu_5p=%.4f (weighted re-accumulation)",
                iteration, state.damage_ctx.sample.tau, state.damage_ctx.sample.mu_5p
            )

        # Debug: After first iteration, print S_anc/S_mod and gamma stats
        if iteration == 0 and config.hierarchical_enabled:
            total_S_anc = 0.0
            total_S_mod = 0.0
            max_gamma = 0.0
            min_gamma = 1.0
            for j_debug in range(n_refs):
                total_S_anc += state.S_anc[j_debug]
                total_S_mod += state.S_mod[j_debug]
                if state.gamma_values[j_debug] > max_gamma:
                    max_gamma = state.gamma_values[j_debug]
                if state.gamma_values[j_debug] < min_gamma:
                    min_gamma = state.gamma_values[j_debug]
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "After iter 0: total_S_anc=%.2f total_S_mod=%.2f frac_anc=%.4f gamma_range=[%.4f, %.4f]",
                total_S_anc, total_S_mod, total_S_anc / (total_S_anc + total_S_mod) if (total_S_anc + total_S_mod) > 0 else 0,
                min_gamma, max_gamma
            )

        # Check convergence
        if check_convergence_mad(conv, config.convergence_tolerance):
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "CONVERGED at iter %d: LL=%.6f dLL=%.2e ||dphi||=%.2e",
                iteration, current_ll, ll_change, param_change
            )
            converged = True
            break

    if not converged:
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "MAX_ITER reached (%d): LL=%.6f dLL=%.2e ||dphi||=%.2e",
            config.max_iterations, current_ll, ll_change, param_change
        )

    # Transform to output (apply power rho POST-PROCESSING)
    transform_to_output(state, output_weights, config.power_rho)

    bf_nogil_logf_notime(
        b"EM_UNIFIED",
        "Output transform: power_rho=%.2f (post-processing only)",
        config.power_rho
    )

    # Copy hierarchical gamma values to pool for output
    cdef int gamma_copy_count = 0
    cdef double min_copy_gamma = 1.0
    cdef double max_copy_gamma = 0.0
    if config.hierarchical_enabled and state.gamma_values != NULL and pool.gamma_values != NULL:
        for j in range(n_refs):
            pool.gamma_values[j] = state.gamma_values[j]
            gamma_copy_count += 1
            if state.gamma_values[j] < min_copy_gamma:
                min_copy_gamma = state.gamma_values[j]
            if state.gamma_values[j] > max_copy_gamma:
                max_copy_gamma = state.gamma_values[j]
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "Copied %d gamma values to pool: range=[%.4f, %.4f]",
            gamma_copy_count, min_copy_gamma, max_copy_gamma
        )
    else:
        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "GAMMA COPY SKIPPED: hierarchical=%d state.gamma=%p pool.gamma=%p",
            <int>config.hierarchical_enabled,
            <void*>state.gamma_values,
            <void*>pool.gamma_values
        )

    # Copy S_anc and S_mod to pool for debugging/analysis
    if config.hierarchical_enabled and state.S_anc != NULL and state.S_mod != NULL:
        if pool.S_anc_accum != NULL and pool.S_mod_accum != NULL:
            for j in range(n_refs):
                pool.S_anc_accum[j] = state.S_anc[j]
                pool.S_mod_accum[j] = state.S_mod[j]

        # Debug: log ref with lowest gamma to check S_tot
        min_gamma_debug = 1.0
        min_idx_debug = 0
        for j in range(n_refs):
            if state.gamma_values[j] < min_gamma_debug and (state.S_anc[j] + state.S_mod[j]) > 0.01:
                min_gamma_debug = state.gamma_values[j]
                min_idx_debug = j

        if min_gamma_debug < 0.5:
            s_tot_debug = state.S_anc[min_idx_debug] + state.S_mod[min_idx_debug]
            bf_nogil_logf_notime(
                b"EM_UNIFIED",
                "DEBUG lowest gamma: ref=%u gamma=%.4f S_anc=%.2f S_mod=%.2f S_tot=%.2f frac_anc=%.3f",
                min_idx_debug, min_gamma_debug, state.S_anc[min_idx_debug], state.S_mod[min_idx_debug], s_tot_debug,
                state.S_anc[min_idx_debug] / s_tot_debug if s_tot_debug > 0 else 0.0
            )

    # Copy final unified damage model results to pool
    if damage_model_enabled and state.damage_ctx != NULL:
        # Final fit and output computation
        compute_outputs(state.damage_ctx)

        # Copy authenticity scores to pool
        if pool.authenticity_scores != NULL:
            update_authenticity_from_damage(state.damage_ctx, pool.authenticity_scores, n_refs)

        # Copy per-reference damage amplitudes to pool arrays
        if pool.damage_amplitude != NULL:
            for j in range(n_refs):
                # Use average of 5' and 3' amplitudes
                pool.damage_amplitude[j] = (state.damage_ctx.ref_params[j].delta_5p +
                                            state.damage_ctx.ref_params[j].delta_3p) * 0.5

        bf_nogil_logf_notime(
            b"EM_UNIFIED",
            "Damage model final: tau=%.2f mu_5p=%.4f mu_3p=%.4f converged=%d",
            state.damage_ctx.sample.tau,
            state.damage_ctx.sample.mu_5p,
            state.damage_ctx.sample.mu_3p,
            <int>state.damage_ctx.sample.converged
        )

    # Cleanup
    if ref_stats_for_auth != NULL:
        free(ref_stats_for_auth)
    if sq != NULL:
        free_squarem_state(sq)
    free_convergence_state(conv)
    free_em_state(prev_state)
    free_em_state(state)

    return iteration + 1


# =============================================================================
# Python-accessible wrapper
# =============================================================================

def run_em_python(pool_capsule, config_dict):
    """
    Python wrapper for EM algorithm.

    Parameters
    ----------
    pool_capsule : PyCapsule
        Capsule containing MemoryPool pointer
    config_dict : dict
        Configuration dictionary with keys:
        - max_iterations: int
        - convergence_tolerance: float
        - dirichlet_prior: float
        - gamma_prior: float
        - power_rho: float
        - unknown_enabled: bool
        - unknown_margin: float
        - hierarchical_enabled: bool
        - D_avg_5p, D_avg_3p, epsilon_error: float
        - squarem_enabled: bool
        - squarem_start_iter: int
        - thread_count: int

    Returns
    -------
    tuple
        (iterations, output_weights) or (error_code, None)
    """
    import numpy as np

    cdef void* pool_ptr = PyCapsule_GetPointer(pool_capsule, <char*>NULL)
    cdef MemoryPool* pool = <MemoryPool*>pool_ptr

    # Build config
    cdef EMConfig config
    config.max_iterations = config_dict.get('max_iterations', 50)
    config.convergence_tolerance = config_dict.get('convergence_tolerance', 1e-6)
    config.dirichlet_prior = config_dict.get('dirichlet_prior', 0.01)
    config.gamma_prior = config_dict.get('gamma_prior', 0.01)
    config.power_rho = config_dict.get('power_rho', 1.0)
    config.unknown_enabled = config_dict.get('unknown_enabled', False)
    config.unknown_margin = config_dict.get('unknown_margin', 2.0)
    config.hierarchical_enabled = config_dict.get('hierarchical_enabled', False)
    config.D_avg_5p = config_dict.get('D_avg_5p', 0.02)
    config.D_avg_3p = config_dict.get('D_avg_3p', 0.02)
    config.epsilon_error = config_dict.get('epsilon_error', 0.01)
    config.squarem_enabled = config_dict.get('squarem_enabled', True)
    config.squarem_start_iter = config_dict.get('squarem_start_iter', 3)
    config.enable_globalization = config_dict.get('enable_globalization', True)
    config.backtrack_factor = config_dict.get('backtrack_factor', 0.5)
    config.max_backtrack_steps = config_dict.get('max_backtrack_steps', 5)
    config.steplength_scheme = config_dict.get('steplength_scheme', 3)  # S3 recommended
    config.squarem_pi_threshold = config_dict.get('squarem_pi_threshold', 0.0)  # Gamma extrapolation fix enabled
    config.thread_count = config_dict.get('thread_count', 1)
    config.history_length = config_dict.get('history_length', 5)

    # CWRP defaults (disabled unless explicitly enabled)
    config.coverage_prior_enabled = False
    config.coverage_prior_lambda = 0.0
    config.authenticity_scores = NULL

    # Unified damage model defaults
    config.unified_damage_enabled = config_dict.get('unified_damage_enabled', False)
    config.damage_update_interval = config_dict.get('damage_update_interval', 5)
    config.is_single_stranded = config_dict.get('is_single_stranded', False)
    config.initial_tau = config_dict.get('initial_tau', 0.0)
    config.fix_tau = config_dict.get('fix_tau', False)

    # Iterative auth defaults
    config.iterative_auth = config_dict.get('iterative_auth', False)
    config.auth_update_interval = config_dict.get('auth_update_interval', 5)
    config.damage_weight = config_dict.get('damage_weight', 1.0)
    config.low_cov_floor_reads = config_dict.get('low_cov_floor_reads', 10)
    config.low_cov_shrink_tau = config_dict.get('low_cov_shrink_tau', 50.0)

    # Posterior-weighted auth (Path B) defaults
    config.auth_post_enabled = config_dict.get('auth_post_enabled', False)
    config.auth_update_interval_post = config_dict.get('auth_update_interval_post', 3)
    config.auth_scale_post = config_dict.get('auth_scale_post', 4.0)
    config.auth_lambda_ramp_iters = config_dict.get('auth_lambda_ramp_iters', 5)

    # Allocate output
    cdef uint32_t n_refs = pool.reference_count
    output_weights = np.zeros(n_refs, dtype=np.float64)
    cdef double[::1] output_view = output_weights

    # Run EM
    cdef int result
    with nogil:
        result = run_em(pool_ptr, &config, &output_view[0])

    if result < 0:
        return (result, None)

    return (result, output_weights)


def execute_em_py(
    uintptr_t pool_ptr,
    int max_iterations,
    double convergence_tolerance,
    double dirichlet_prior,
    double power_rho,
    bint unknown_enabled,
    double unknown_margin,
    bint hierarchical_enabled,
    double D_avg_5p,
    double D_avg_3p,
    double epsilon_error,
    bint squarem_enabled,
    int squarem_start_iter,
    bint enable_globalization,
    double backtrack_factor,
    int max_backtrack_steps,
    int steplength_scheme,
    int thread_count,
):
    """
    Python wrapper for EM algorithm.

    Updates the memory pool's reference weights in place.

    Parameters
    ----------
    pool_ptr : uintptr_t
        Pointer to MemoryPool
    max_iterations : int
        Maximum EM iterations
    convergence_tolerance : double
        Convergence threshold
    dirichlet_prior : double
        Dirichlet prior alpha
    power_rho : double
        Power transform exponent (1.0 = no transform, <1 = flattening)
    unknown_enabled : bool
        Enable unknown component
    unknown_margin : double
        Margin below best score for unknown
    hierarchical_enabled : bool
        Enable hierarchical ancient/modern
    D_avg_5p, D_avg_3p : double
        Average PMD damage rates
    epsilon_error : double
        Sequencing error rate
    squarem_enabled : bool
        Enable SQUAREM acceleration
    squarem_start_iter : int
        Start SQUAREM after this iteration
    enable_globalization : bool
        Enable backtracking
    backtrack_factor : double
        Backtracking step reduction
    max_backtrack_steps : int
        Maximum backtracking attempts
    thread_count : int
        Number of threads

    Returns
    -------
    int
        0 on success, -1 on failure
    """
    import numpy as np

    cdef MemoryPool* pool = <MemoryPool*>pool_ptr
    cdef uint32_t n_refs = pool.reference_count
    cdef uint32_t j

    # Build config
    cdef EMConfig config
    config.max_iterations = max_iterations
    config.convergence_tolerance = convergence_tolerance
    config.dirichlet_prior = dirichlet_prior
    # Gamma prior should be small relative to per-reference evidence (~0.01-0.1)
    # With 100k+ references, total evidence ~1e4, so per-ref ~0.1
    # Prior of 0.01 allows data to dominate while still regularizing
    config.gamma_prior = 0.01
    config.power_rho = power_rho
    config.unknown_enabled = unknown_enabled
    config.unknown_margin = unknown_margin
    config.hierarchical_enabled = hierarchical_enabled

    # Use pool's D_avg values when hierarchical enabled (set by init_hierarchical_em_py)
    if hierarchical_enabled and pool.hierarchical_em_enabled:
        config.D_avg_5p = pool.D_avg_5p
        config.D_avg_3p = pool.D_avg_3p
        config.epsilon_error = pool.epsilon_error
    else:
        config.D_avg_5p = D_avg_5p
        config.D_avg_3p = D_avg_3p
        config.epsilon_error = epsilon_error
    config.squarem_enabled = squarem_enabled
    config.squarem_start_iter = squarem_start_iter
    config.enable_globalization = enable_globalization
    config.backtrack_factor = backtrack_factor
    config.max_backtrack_steps = max_backtrack_steps
    config.steplength_scheme = steplength_scheme
    config.squarem_pi_threshold = 0.0  # Gamma extrapolation fix enabled
    config.thread_count = thread_count
    config.history_length = 5  # Default

    # CWRP: Use pool settings if available
    if pool.cwrp_enabled and pool.authenticity_scores != NULL:
        config.coverage_prior_enabled = True
        config.coverage_prior_lambda = pool.cwrp_lambda
        config.authenticity_scores = pool.authenticity_scores
    else:
        config.coverage_prior_enabled = False
        config.coverage_prior_lambda = 0.0
        config.authenticity_scores = NULL

    # Iterative ancientness field: Use pool settings
    config.iterative_auth = pool.iterative_auth_enabled
    config.auth_update_interval = pool.auth_update_interval
    config.damage_weight = pool.damage_weight
    config.low_cov_floor_reads = pool.low_cov_floor_reads
    config.low_cov_shrink_tau = pool.low_cov_shrink_tau

    # Posterior-weighted auth (Path B): Use pool settings
    config.auth_post_enabled = pool.auth_update_interval_post > 0
    config.auth_update_interval_post = pool.auth_update_interval_post if pool.auth_update_interval_post > 0 else 3
    config.auth_scale_post = pool.auth_scale_post if pool.auth_scale_post > 0 else 4.0
    config.auth_lambda_ramp_iters = pool.auth_lambda_ramp_iters if pool.auth_lambda_ramp_iters >= 0 else 5

    # Unified damage model: enabled when hierarchical mode is active (PMD processing)
    # The model will accumulate counts from alignments if pool.ref_damage_stats is NULL
    config.unified_damage_enabled = config.hierarchical_enabled and (pool.alignment_count > 0)
    config.damage_update_interval = 5  # Update every 5 iterations
    config.is_single_stranded = False  # Default to double-stranded (most common for ancient DNA)

    # Extract tau from PMD curve if available (uses per-position library-wide counts)
    cdef PMDCurve* pmd_curve = NULL
    config.initial_tau = 0.0
    config.fix_tau = False
    if pool.pmd_curve_ptr != NULL:
        pmd_curve = <PMDCurve*>pool.pmd_curve_ptr
        if pmd_curve.lambda_decay > 0.01:
            config.initial_tau = 1.0 / pmd_curve.lambda_decay
            config.fix_tau = True

    # Sample-level P(ancient) gate: 0.0 means compute from PMD curve automatically
    config.sample_pi_override = pool.sample_pi_override

    # Allocate output buffer
    output_weights = np.zeros(n_refs, dtype=np.float64)
    cdef double[::1] output_view = output_weights

    # Run EM
    cdef int iterations
    cdef void* pool_void = <void*>pool_ptr
    cdef double* ref_weights
    with nogil:
        iterations = run_em(pool_void, &config, &output_view[0])

    if iterations < 0:
        return -1

    # Copy output weights back to pool's reference weights
    # This matches what execute_em_algorithm does
    ref_weights = get_reference_weights(pool)
    for j in range(n_refs):
        ref_weights[j] = output_view[j]

    # Update pool iteration count and convergence flag
    pool.iteration_count = iterations
    pool.algorithm_converged = (iterations < max_iterations)

    return 0
