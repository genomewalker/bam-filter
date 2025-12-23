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

"""High-performance BAM file processing engine.

This module provides the core C/Cython implementation for BAM file processing,
including alignment scoring, EM algorithm execution, graph analysis, and filtered
BAM output generation.
"""

# Python imports
from cython.parallel import prange, threadid
import os

# C standard library imports
from libc.stdint cimport int32_t, int64_t, uint16_t, uint32_t, uint64_t, uint8_t, uintptr_t
from libc.math cimport exp, fabs, fmax, fmin, log, log2, sqrt as libc_sqrt, pow as libc_pow, INFINITY
from libc.stdlib cimport calloc, free, malloc, realloc
from libc.stdio cimport printf, sprintf
from libc.string cimport memcpy, memmove, memset, strlen
from libc.time cimport clock, clock_t, CLOCKS_PER_SEC

# Project module imports
from .batch_utils cimport create_balanced_batches_greedy
from .processor_memory cimport create_memory_pool, create_memory_pool_split, destroy_memory_pool, shrink_memory_pool, cleanup_presorted_memory, cleanup_em_intermediate_memory
from .processor_fast_math cimport stable_log_sum_exp, safe_normalize_weights
from .processor_types cimport PrecomputedWeights, EMAlgorithmConfig
from .processor_graph_ops cimport WeightedGraph, destroy_weighted_graph, pick_min_edge_weight_elbow
from .processor_igraph cimport igraph_t, igraph_vcount, igraph_ecount

from .processor_batch cimport (
    BatchAlignment,
    ProcessingBatch,
    create_processing_batch,
    destroy_processing_batch,
    process_batch_alignments,
    process_batch_alignments_with_pmd,
    assign_global_sequential_ids_fast,
    count_actual_alignments,
    count_unique_refs_from_batches,
    parallel_streaming_stats_optimized,
    populate_memory_pool_direct,
    populate_memory_pool_filtered,
    populate_memory_pool_filtered_split,
    count_alignments_passing_ani_filter,
    count_unique_reads_from_thread_maps,
)

from .processor_pmd cimport (
    PMDGlobalContext,
    PMDCurveParams,
    PMDCurve,
    create_pmd_context,
    destroy_pmd_context,
    merge_pmd_stats,
    fit_pmd_curve,
    finalize_pmd_model,
)
# Python wrapper for applying PMD corrections to alignments
from .processor_pmd import apply_pmd_corrections_py, apply_raw_ani_filter_py, init_hierarchical_em_py
# NOTE: processor_em.execute_em_py is imported lazily in run_reassignment_pipeline
# to avoid circular import (processor_em cimports from processor)
# NOTE: processor_damage_model is imported lazily in run_reassignment_pipeline
# to avoid circular import (it cimports from processor)

from .processor_filters cimport (
    apply_probability_filtering,
    apply_cluster_aware_filtering
)

from .processor_taxonomy_filters cimport (
    TaxonomyFilterStats
)

from .processor_graph_ops cimport (
    extract_neighbors_from_igraph
)

from .processor_stats cimport (
    init_processing_stats,
    update_initial_stats,
    update_quality_filter_stats,
    update_em_stats,
    update_probability_filter_stats,
    update_graph_stats,
    update_unified_filter_stats,
    update_final_output_stats,
    calculate_summary_metrics,
    print_processing_stats
)

from .processor_precomputed cimport (
    create_precomputed_weights,
    update_precomputed_weights,
    free_precomputed_weights,
)

from .processor_bam_writer cimport (
    CompactAlignment, LookupTable,
    create_lookup_table, destroy_lookup_table,
    WriteBatch, create_write_batch, destroy_write_batch,
    copy_bam_record, write_batch_to_bam,
    write_filtered_bam, write_bam_with_filtered_header,
)

from .processor_hash cimport (
    compute_read_name_hash,
    ThreadLocalHashMap,
    create_thread_local_hash_map,
    destroy_thread_local_hash_map,
    extract_read_hash_identifier,
)

from .processor_md_quality cimport (
    initialize_quality_lookup_tables,
    calculate_md_quality_score,
    calculate_md_score_fast_path,
    calculate_md_score_with_pmd_single_pass,
)

from .processor_mapping cimport (
    ReferenceMapping,
    create_reference_mapping,
    destroy_reference_mapping,
    remap_alignment_reference_ids,
    update_reference_mapping_after_filtering,
)

from .processor_graph cimport (
    ReferencePattern, ReferenceStats, ReadIndex,
    analyze_reference_graph,
    calculate_reference_stats,
    calculate_reference_coverage,
    compute_authenticity_scores,
    init_cwrp,
    build_read_index_parallel,
    destroy_read_index,
    write_graph_tsv
)

from .processor_graph_taxonomy cimport (
    TaxonomyGraphConfig,
    enrich_patterns_with_taxonomy,
    detect_taxonomy_anomalies
)

from .processor_taxonomy_filters cimport (
    TaxonomyFilterConfig
)

from .processor_network_qc cimport (
    NetworkQCConfig
)

from .probabilistic_profiler cimport (
    apply_gmrf_smoothing_to_graph,
)

from .taxonomy_db cimport (
    TaxonomyDB, AccessionMap,
    TaxonomyDatabase, AccessionMapping
)
from .taxonomy_db import load_accession_map_from_file

from .reference_lengths cimport (
    TSVReferenceMap,
    load_tsv_reference_file,
    get_tsv_reference_count,
    lookup_reference_length,
    free_tsv_reference_map,
    print_tsv_reference_stats,
)
cdef extern from "bam_filter/c_logging.h":
    double bf_monotonic_seconds() nogil

# Common helpers
from bam_filter import logging as bf_logging

LOG_TAG = "PROCESS"


def _info(message: str, *args) -> None:
    bf_logging.log(LOG_TAG, message, *args)


def _warn(message: str, *args) -> None:
    bf_logging.warn(message, *args)


def _error(message: str, *args) -> None:
    bf_logging.error(message, *args)


def _debug(level: int, message: str) -> None:
    bf_logging.verbose(level, LOG_TAG, "%s", message)


def _log_stage(stage: str, start_time: float, level: int = 0) -> None:
    duration = bf_monotonic_seconds() - start_time
    bf_logging.verbose(level, LOG_TAG, "stage=%s duration=%.2fs", stage, duration)


def _announce_stage(title: str, detail: str = "") -> None:
    """Emit a human-readable stage separator for pipeline progress."""
    bf_logging.summary("")
    bf_logging.summary("┌─ %s", title)
    if detail:
        bf_logging.summary("│ %s", detail)
    bf_logging.summary("└─────────────────────────────────────────────────────────────")


def _announce_stage_skip(title: str, reason: str = "") -> None:
    """Emit a stage separator marking a skipped pipeline phase."""
    bf_logging.summary("")
    bf_logging.summary("┌─ %s (skipped)", title)
    if reason:
        bf_logging.summary("│ Skipped: %s", reason)
    bf_logging.summary("└─────────────────────────────────────────────────────────────")

from .common_helpers cimport (
    pack_position_length,
    extract_position,
    extract_length,
    min_int64,
    max_int64,
    min_int32,
    max_int32,
    min_double,
    max_double,
    page_size,
)

# External C bindings (centralized in processor_types)
from .processor_types cimport (
    BGZF, bam1_core_t, bam1_t, sam_hdr_t, hts_idx_t, hts_itr_t,
    htsFile, samFile, sam_hdr_read, sam_hdr_destroy, sam_hdr_str,
    sam_hdr_tid2name, sam_hdr_tid2len, sam_hdr_parse, sam_hdr_nref,
    bam_endpos, sam_read1, sam_write1, sam_hdr_write, bam_init1,
    bam_destroy1, bam_dup1, bam_get_qname, bam_aux_get, bam_aux2i,
    bam_aux_del, bam_aux_append, sam_index_load, hts_idx_destroy,
    hts_idx_get_stat, sam_itr_queryi, sam_itr_next, hts_open, hts_close,
    sam_hdr_name2tid, hts_set_threads, bam_get_seq, bam_get_qual,
    bam_get_cigar, seq_nt16_str, getpagesize, bam_seqi_wrapper
)


# Module-level constants
cdef void* MAP_FAILED_PTR = <void*>-1



cdef enum ProcessingError:
    PROCESSING_SUCCESS = 0
    PROCESSING_ERROR_FILE_ACCESS = 1
    PROCESSING_ERROR_MEMORY_ALLOCATION = 2
    PROCESSING_ERROR_INVALID_DATA = 3
    PROCESSING_ERROR_ALGORITHM_FAILURE = 4




cdef void validate_pmd_scores(MemoryPool* pool, int64_t max_check=100) noexcept nogil:
    """Validate PMD score data in memory pool for debugging.

    Parameters
    ----------
    pool : MemoryPool*
        Memory pool containing alignments with PMD scores
    max_check : int64_t
        Maximum number of alignments to check
    """
    if not pool.pmd_enabled_for_output:
        printf(b"PMD validation: PMD output disabled\n")
        return

    cdef int64_t non_zero_pmd = 0
    cdef int64_t check_count = min_int64(pool.alignment_count, max_check)
    cdef float min_pmd = 1e30, max_pmd = -1e30
    cdef int i
    cdef float pmd_val

    for i in range(check_count):
        pmd_val = pool.alignments[i].pmd_score
        if pmd_val != 0.0:
            non_zero_pmd += 1
            if pmd_val < min_pmd:
                min_pmd = pmd_val
            if pmd_val > max_pmd:
                max_pmd = pmd_val

    printf(b"PMD validation: Checked %ld alignments\n", check_count)
    printf(b"  Non-zero PMD scores: %ld (%.1f%%)\n", non_zero_pmd,
           100.0 * non_zero_pmd / check_count)

    if non_zero_pmd > 0:
        printf(b"  PMD score range: %.4f to %.4f\n", min_pmd, max_pmd)

    if check_count > 0:
        printf(b"  Sample PMD values (first 3):\n")
        for i in range(min_int64(3, check_count)):
            printf(b"    [%d] PMD: %.4f\n", i, pool.alignments[i].pmd_score)

cdef uint32_t INVALID_SEQUENTIAL_ID = 0xFFFFFFFF


cdef score_alignments(
    bam_file,
    output_bam=None,
    num_threads=1,
    verbose=False,
    calculate_pmd=True,
    library_type="ds",
    hierarchical_pmd=False,
    reference_lengths_tsv=None,
    reference_stats_tsv=None,
    min_read_count=1,
    min_read_ani=0.0,
    min_read_length=30,
    max_read_length=10000,
    use_em=True,
    max_em_iterations=25,
    em_tolerance=1e-5,
    min_probability=1e-6,
    prob_fraction=0.3,
    prior_weight=0.01,
    lambda_scale=0.3,
    use_squarem_acceleration=True,
    enable_globalization=True,
    squarem_start_iter=2,
    backtrack_factor=0.5,
    max_backtrack_steps=5,
    steplength_scheme=3,
    enable_dominance_regularization=True,
    auto_tune_penalties=True,
    dominance_strength=None,
    use_adaptive_dominance=True,
    entropy_scaling_factor=1.0,
    min_penalty_strength=0.1,
    max_penalty_strength=5.0,
    entropy_confidence_threshold=0.0,
    enable_emergency_regularization=True,
    emergency_entropy_threshold=1.5,
    emergency_max_weight_threshold=0.9,
    emergency_likelihood_drop=1e-4,
    emergency_min_dominant_refs=0,
    init_prior_strength=0.1,
    information_threshold=-999.0,
    em_beta=1.0,
    em_length_correction=False,
    em_unknown_component=False,
    em_unknown_prior=0.05,
    em_unknown_score=-50.0,
    em_length_init=False,
    # === UNIFIED φ-SPACE EM PARAMETERS ===
    em_power_rho=1.0,              # ρ: 1.0 = standard, <1 reduces dominance
    em_unknown_adaptive=False,     # Enable adaptive unknown score (vs fixed s_unknown)
    em_unknown_margin=2.0,         # Δ: margin below best score for unknown
    em_length_output_exp=1.0,      # γ_len: 0=none, 1=per-base abundance
    em_length_prior_exp=1.0,       # η_prior: 0=flat, 1=length-proportional
    graph_min_edge_weight=2,
    clustering=False,
    community_resolution=1.0,
    community_max_iterations=10,
    outlier_method="otsu",
    graph_export=None,
    # Taxonomy parameters
    taxonomy_db=None,
    taxonomy_accession_map=None,
    taxonomy_min_rank=6,
    taxonomy_cross_domain_threshold=0.10,
    taxonomy_kingdom_threshold=0.25,
    taxonomy_genus_threshold=0.50,
    # Taxonomy-informed filtering parameters
    taxonomy_filter_enabled=False,
    taxonomy_strict_filter=True,
    taxonomy_strict_min_connections=5,
    taxonomy_weighted_outlier=True,
    taxonomy_anomaly_weight=2.0,
    taxonomy_second_chance=True,
    taxonomy_second_chance_cc=0.3,
    # Cross-domain removal options
    remove_cross_domain_alignments=False,
    remove_cross_domain_references=False,
    detect_misannotations=False,
    # Network QC & Taxonomic ambiguity detection parameters
    network_qc_filter=False,
    entropy_biased_threshold=1.0,
    entropy_mixed_threshold=2.0,
    entropy_highly_mixed_threshold=3.0,
    tax_ambiguity_removal_level=2,
    # Coverage-Weighted Reference Priors (CWRP)
    cwrp_lambda=0.0,
    iterative_auth=False,
    auth_update_interval=5,
    damage_weight=1.0,
    low_cov_floor=10,
    low_cov_shrink_tau=50.0,
    # Posterior-Weighted Coverage Authenticity (Path B)
    auth_post_enabled=False,
    auth_update_interval_post=3,
    auth_scale_post=4.0,
    auth_lambda_ramp_iters=5,
    # Sample-level P(ancient) gate
    sample_pi_override=0.0,
    # GMRF smoothing (profile-only)
    enable_gmrf=False,
    gmrf_tau=1.0,
):
    """Core BAM alignment scoring and filtering engine.

    Processes BAM file through complete pipeline: alignment scoring, EM algorithm,
    probability filtering, optional graph analysis and clustering, and filtered
    BAM output generation.

    Parameters
    ----------
    bam_file : str
        Input BAM file path
    output_bam : str, optional
        Output BAM file path
    reference_lengths_tsv : str, optional
        TSV file with reference length overrides (ref_name<tab>length)
    reference_stats_tsv : str, optional
        Output path for reference statistics TSV
    calculate_pmd : bool
        Calculate Post-Mortem Damage scores
    use_squarem_acceleration : bool
        Enable SQUAREM acceleration for EM (Varadhan & Roland 2008)
    clustering : bool
        Enable Community clustering for cluster-aware filtering
    graph_export : str, optional
        Export graph to GraphML format

    Returns
    -------
    dict
        Processing results with alignment counts, EM statistics, and metadata
    """

    # Initialize result variable
    result = {
        'summary': {
            'total_alignments': 0,
            'filtered_alignments': 0,
            'unique_reads': 0,
            'unique_references': 0,
            'em_iterations': 0,
            'em_converged': False,
            'final_likelihood': 0.0
        },
        'memory_pool_available': False,
        'memory_optimized': True,
        'pmd_enabled': calculate_pmd,
        'integrated_pmd': True,
    }

    if verbose:
        _info(f"Processing BAM file with reference length overrides")
        if reference_lengths_tsv:
            _info(f"  TSV file: {reference_lengths_tsv}")
        _info(f"  PMD calculation enabled" if calculate_pmd else "  PMD calculation disabled")

    # Configure parameters
    cdef AlignmentScoringConfig scoring_config
    scoring_config.minimum_read_identity = min_read_ani
    scoring_config.minimum_read_length = min_read_length
    scoring_config.maximum_read_length = max_read_length
    scoring_config.global_min_score = 1e30
    scoring_config.global_max_score = -1e30
    scoring_config.calculate_pmd = calculate_pmd
    scoring_config.is_single_stranded = (library_type == "ss")

    cdef EMAlgorithmConfig em_config
    em_config.maximum_iterations = max_em_iterations
    em_config.convergence_tolerance = em_tolerance
    em_config.minimum_probability_threshold = min_probability
    em_config.probability_fraction_filter = prob_fraction
    em_config.regularization_weight = prior_weight
    em_config.use_squarem_acceleration = use_squarem_acceleration
    em_config.enable_acceleration = use_squarem_acceleration
    em_config.thread_count = num_threads
    em_config.minimum_read_coverage = min_read_count
    em_config.squarem_start_iter = squarem_start_iter
    em_config.enable_globalization = enable_globalization
    em_config.backtrack_factor = backtrack_factor
    em_config.max_backtrack_steps = max_backtrack_steps
    em_config.steplength_scheme = steplength_scheme
    em_config.enable_dominance_regularization = enable_dominance_regularization
    em_config.use_adaptive_dominance = use_adaptive_dominance
    em_config.entropy_scaling_factor = entropy_scaling_factor
    em_config.min_penalty_strength = min_penalty_strength
    em_config.max_penalty_strength = max_penalty_strength
    em_config.entropy_confidence_threshold=entropy_confidence_threshold

    em_config.enable_emergency_regularization = enable_emergency_regularization
    em_config.dominance_regularization_active = False  # Start inactive
    em_config.emergency_entropy_threshold = emergency_entropy_threshold
    em_config.emergency_max_weight_threshold = emergency_max_weight_threshold
    em_config.emergency_likelihood_drop = emergency_likelihood_drop
    em_config.emergency_min_dominant_refs = 0  # Will be auto-calculated
    em_config.init_prior_strength = init_prior_strength
    em_config.information_threshold = information_threshold

    # Tempered EM: beta < 1 reduces rich-gets-richer bias
    em_config.em_beta = em_beta

    # Effective length correction: reduces length bias in EM
    em_config.em_length_correction = em_length_correction

    # Unknown component: absorbs reads not belonging to any reference
    em_config.em_unknown_component = em_unknown_component
    em_config.em_unknown_prior = em_unknown_prior
    em_config.em_unknown_score = em_unknown_score

    # Length-aware initialization
    em_config.em_length_init = em_length_init

    # === UNIFIED φ-SPACE EM PARAMETERS ===
    em_config.em_power_rho = em_power_rho
    em_config.em_unknown_adaptive = em_unknown_adaptive
    em_config.em_unknown_margin = em_unknown_margin
    em_config.em_length_output_exp = em_length_output_exp
    em_config.em_length_prior_exp = em_length_prior_exp

    # Dominance regularization is disabled - use em_power_rho for anti-dominance behavior
    em_config.enable_dominance_regularization = False
    em_config.dominance_strength = 0.0

    # Initialize C variables
    cdef samFile* bam_handle = NULL
    cdef sam_hdr_t* bam_header = NULL
    cdef hts_idx_t* bam_index = NULL
    cdef int64_t* reference_ids = NULL
    cdef int64_t* reference_alignment_counts = NULL
    cdef int64_t* reference_lengths = NULL  # All reference lengths from BAM header
    cdef int64_t* filtered_reference_lengths = NULL  # Filtered reference lengths
    cdef ProcessingBatch** processing_batches = NULL
    cdef MemoryPool* memory_pool = NULL
    cdef samFile** thread_handles = NULL
    cdef ThreadLocalHashMap** thread_maps = NULL

    # PMD context for stats collection
    cdef PMDGlobalContext* pmd_context = NULL
    cdef PMDCurveParams pmd_params
    cdef bint collect_pmd_stats = False
    cdef bint hierarchical_pmd_c = hierarchical_pmd

    # TSV integration variables
    cdef TSVReferenceMap* tsv_map = NULL
    cdef bytes tsv_file_bytes
    cdef const char* tsv_file_path_c = NULL

    cdef int num_threads_c
    cdef const char* bam_file_path
    cdef const char* output_bam_path
    cdef const char* graph_export_path_c = NULL
    cdef bytes bam_file_bytes
    cdef bytes output_bam_bytes
    cdef int64_t total_references, references_to_process = 0, total_alignments = 0
    cdef int64_t batch_count, reference_idx, batch_idx
    cdef uint64_t mapped_alignments, unmapped_alignments
    cdef int32_t thread_idx
    cdef int64_t total_valid_alignments = 0
    cdef int64_t unique_read_count = 0
    cdef int64_t unique_reference_count = 0
    cdef int ret
    cdef int64_t actual_total_alignments = 0
    cdef int64_t actual_unique_read_count = 0
    cdef int32_t bam_write_threads

    # Early reference length variables
    cdef int32_t tid32, dense, i
    cdef const char* ref_name = NULL
    cdef int64_t tsv_length, bam_length
    cdef int32_t used_tsv = 0, used_bam = 0
    cdef int32_t nentries = 0
    cdef char* ref_selected = NULL

    # Gamma remapping variables (for preserving EM gamma values after reference compaction)
    cdef double* old_gamma = NULL
    cdef double* new_gamma = NULL
    cdef uint32_t old_idx, new_idx

    # Filtered length validation variables
    cdef int64_t min_filt_len, max_filt_len, total_filt_len
    cdef int invalid_original_lengths = 0
    cdef int invalid_filtered_lengths = 0

    # Batch array declarations
    cdef int64_t* batch_starts = NULL
    cdef int64_t* batch_ends = NULL
    cdef int64_t max_batches
    cdef int64_t batch_start, batch_end
    cdef int64_t expected_alignments
    cdef int64_t ref_idx
    cdef ReferencePattern* pattern_data = NULL
    cdef ReferenceMapping* mapping = NULL

    # Streaming stats variables
    cdef double min_score = 0.0, max_score = 0.0, mean_score = 0.0, variance_score = 0.0
    cdef int64_t n_score = 0
    cdef double pipeline_start = 0.0
    cdef double stage_timer = 0.0

    # Taxonomy filtering statistics
    cdef TaxonomyFilterStats taxonomy_stats

    # Coverage-Weighted Reference Priors
    cdef double c_cwrp_lambda = cwrp_lambda
    cdef bint c_iterative_auth = iterative_auth
    cdef int32_t c_auth_update_interval = auth_update_interval
    cdef double c_damage_weight = damage_weight
    cdef int32_t c_low_cov_floor = low_cov_floor
    cdef double c_low_cov_shrink_tau = low_cov_shrink_tau

    # Posterior-Weighted Coverage Authenticity (Path B)
    cdef bint c_auth_post_enabled = auth_post_enabled
    cdef int32_t c_auth_update_interval_post = auth_update_interval_post
    cdef double c_auth_scale_post = auth_scale_post

    # Sample-level P(ancient) gate
    cdef double c_sample_pi_override = sample_pi_override

    # ANI filtering variables (for filtered memory pool creation)
    cdef float min_ani_threshold = 90.0
    cdef float c_epsilon = 0.01  # Sequencing error rate
    cdef PMDCurve* curve_ptr = NULL
    cdef int64_t filtered_alignment_count = 0
    cdef int32_t c_auth_lambda_ramp_iters = auth_lambda_ramp_iters

    # Taxonomy integration variables
    cdef TaxonomyDB* taxdb_c = NULL
    cdef AccessionMap* accmap_c = NULL
    cdef TaxonomyDatabase taxdb_obj
    cdef AccessionMapping accmap_obj
    cdef TaxonomyGraphConfig tax_config
    cdef uint32_t** neighbor_lists_tax = NULL
    cdef uint32_t* neighbor_counts_tax = NULL
    cdef bint neighbor_lists_allocated_from_igraph = False  # Track if we need deep free
    cdef double taxonomy_stage_timer = 0.0

    # Taxonomy-informed filtering variables
    cdef TaxonomyFilterConfig taxonomy_filter_config

    # Network QC filtering variables
    cdef NetworkQCConfig network_qc_config
    cdef uint8_t tax_ambiguity_removal_level_c

    # Cluster filtering variables
    cdef ReferenceStats* ref_stats = NULL
    cdef ReadIndex* read_index = NULL
    cdef WeightedGraph* filtered_graph = NULL
    cdef int cluster_result
    cdef int verbose_c  # C int version of verbose for nogil contexts
    cdef int use_community_c  # C int version for clustering algorithm choice (0=union-find, 1=Leiden)
    cdef double community_res_c  # C double for community_resolution
    cdef int community_parallel_c  # C int for community_parallel
    cdef int community_max_iter_c  # C int for community_max_iterations
    cdef int outlier_method_c  # C int for outlier_method (0=MAD, 1=IQR, 2=IFOREST, 3=LOF, 4=ZSCORE)
    cdef uint32_t graph_min_edge_weight_c  # C uint32_t for graph_min_edge_weight
    cdef uint32_t thr_c  # temporary holder for edge weight threshold
    cdef uint32_t min_read_count_c  # C copy of min_read_count for nogil calls
    cdef int64_t aligns_filtered_information = 0
    cdef int64_t alignments_before_filtering = 0  # Alignments before filtering

    # GMRF smoothing variables
    cdef double gmrf_stage_timer
    cdef int gmrf_result
    cdef int32_t gmrf_i
    cdef uint32_t* gmrf_n_reads = NULL
    cdef double gmrf_tau_c = gmrf_tau
    cdef bint enable_gmrf_c = enable_gmrf

    try:
        pipeline_start = bf_monotonic_seconds()
        num_threads_c = num_threads
        bam_file_bytes = bam_file.encode('utf-8')
        bam_file_path = bam_file_bytes

        # Copy python parameters into C types for nogil calls
        min_read_count_c = <uint32_t>min_read_count
        # C int version of verbose for use inside nogil regions
        verbose_c = 1 if verbose else 0
        # Convert outlier_method string to C int: 0=MAD (default), 1=IQR, 2=IFOREST, 3=LOF, 4=ZSCORE
        outlier_method_str = str(outlier_method) if outlier_method else "mad"
        method_lower = outlier_method_str.lower()
        if method_lower == "iqr":
            outlier_method_c = 1
        elif method_lower == "iforest":
            outlier_method_c = 2
        elif method_lower == "lof":
            outlier_method_c = 3
        elif method_lower == "zscore":
            outlier_method_c = 4
        else:  # "mad" or any other value defaults to MAD
            outlier_method_c = 0

        # Assign graph_min_edge_weight_c with special handling for auto/no filtering
        # Semantics:
        #   -1 => no filtering (we represent as 0 in C threshold to indicate "keep all")
        #    0 => auto: don't decide a threshold yet, compute later when a read_index is available
        #   >0 => explicit value provided by user
        if graph_min_edge_weight == -1:
            graph_min_edge_weight_c = 0  # No filtering
            if verbose:
                bf_logging.log("GRAPH", f"Mode chosen: none (no filtering). CLI value={graph_min_edge_weight}. threshold={graph_min_edge_weight_c}")
        elif graph_min_edge_weight == 0:
            # Auto: defer selection until we have a read_index; use 0 as a sentinel
            graph_min_edge_weight_c = 0
            if verbose:
                bf_logging.log("GRAPH", f"Mode chosen: auto (deferred until read index available). CLI value={graph_min_edge_weight}")
        else:
            graph_min_edge_weight_c = <uint32_t>graph_min_edge_weight  # Use specified value
            if verbose:
                bf_logging.log("GRAPH", f"Mode chosen: explicit value. CLI value={graph_min_edge_weight}. threshold={graph_min_edge_weight_c}")

        # Configure taxonomy-informed filtering
        taxonomy_filter_config.enabled = taxonomy_filter_enabled
        # Enable strict filtering if --remove-cross-domain-references is set (or if legacy strict_filter is True and no alignment removal)
        taxonomy_filter_config.enable_strict_filtering = remove_cross_domain_references or (taxonomy_strict_filter and not remove_cross_domain_alignments)
        taxonomy_filter_config.strict_min_connections = <uint32_t>taxonomy_strict_min_connections
        taxonomy_filter_config.enable_weighted_outlier_detection = taxonomy_weighted_outlier
        taxonomy_filter_config.taxonomy_anomaly_weight = <float>taxonomy_anomaly_weight
        taxonomy_filter_config.enable_second_chance = taxonomy_second_chance and taxonomy_filter_config.enable_strict_filtering
        taxonomy_filter_config.second_chance_cc_threshold = <float>taxonomy_second_chance_cc

        # Configure network QC filtering (taxonomic ambiguity detection)
        network_qc_config.high_entropy_threshold = 2.0  # Default
        network_qc_config.entropy_biased_threshold = entropy_biased_threshold
        network_qc_config.entropy_mixed_threshold = entropy_mixed_threshold
        network_qc_config.entropy_highly_mixed_threshold = entropy_highly_mixed_threshold
        network_qc_config.taxonomy_rank_for_analysis = 6  # genus level
        network_qc_config.use_taxonomy_as_community = False
        network_qc_config.community_resolution = community_resolution
        tax_ambiguity_removal_level_c = <uint8_t>tax_ambiguity_removal_level

        # We'll set taxonomy_enabled on stats once the memory pool exists

        if verbose:
            _info("Taxonomy filtering config: enabled=%s, strict=%s, weighted=%s, second_chance=%s",
                  taxonomy_filter_enabled, taxonomy_filter_config.enable_strict_filtering, taxonomy_weighted_outlier, taxonomy_second_chance)
            if remove_cross_domain_alignments:
                _info("Cross-domain alignment removal: ENABLED")
            if remove_cross_domain_references:
                _info("Cross-domain reference removal: ENABLED")
            if network_qc_filter:
                _info("Network QC filtering: ENABLED (entropy_thresholds=%.1f/%.1f/%.1f, removal_level=%d)",
                      entropy_biased_threshold, entropy_mixed_threshold, entropy_highly_mixed_threshold, tax_ambiguity_removal_level)

        if output_bam:
            output_bam_bytes = output_bam.encode('utf-8')
            output_bam_path = output_bam_bytes
        else:
            output_bam_path = NULL

        # Load TSV reference file if provided
        if reference_lengths_tsv:
            if not os.path.exists(reference_lengths_tsv):
                _warn(f"TSV file not found: {reference_lengths_tsv}, using BAM header only")
                reference_lengths_tsv = None
            else:
                tsv_file_bytes = reference_lengths_tsv.encode('utf-8')
                tsv_file_path_c = tsv_file_bytes

                _info(f"Loading reference lengths from: {reference_lengths_tsv}")

                with nogil:
                    tsv_map = load_tsv_reference_file(tsv_file_path_c)

                if not tsv_map:
                    _warn("Failed to load TSV file, using BAM header only")
                else:
                    _info(f"Successfully loaded reference length overrides")
                    with nogil:
                        print_tsv_reference_stats(tsv_map)

        if reference_stats_tsv is not None:
            tsv_file_bytes = reference_stats_tsv.encode('utf-8') if reference_stats_tsv != '' else b''
            tsv_file_path_c = tsv_file_bytes

        if graph_export is not None:
            graph_export_bytes = graph_export.encode('utf-8') if graph_export != '' else b''
            graph_export_path_c = graph_export_bytes

        # Open BAM file and initialize
        bam_handle = hts_open(bam_file_path, b"r")
        if not bam_handle:
            raise RuntimeError(f"Failed to open BAM file: {bam_file}")

        bam_header = sam_hdr_read(bam_handle)
        if not bam_header:
            raise RuntimeError("Failed to read BAM header")

        bam_index = sam_index_load(bam_handle, bam_file_path)
        if not bam_index:
            raise RuntimeError(f"Failed to load BAM index for: {bam_file}")

        total_references = bam_header.n_targets
        if total_references <= 0:
            raise RuntimeError("No references found in BAM header")

        # Resolve reference lengths for all references
        if verbose:
            _info(f"Resolving reference lengths for {total_references} references")

        # Allocate reference lengths array for ALL original references
        reference_lengths = <int64_t*>malloc(total_references * sizeof(int64_t))
        if not reference_lengths:
            raise MemoryError("Failed to allocate reference lengths array")

        # Initialize ALL reference lengths with BAM defaults first
        with nogil:
            for i in range(total_references):
                bam_length = sam_hdr_tid2len(bam_header, i)
                if bam_length > 0:
                    reference_lengths[i] = bam_length
                else:
                    reference_lengths[i] = 1000  # Consistent default
                    invalid_original_lengths += 1

        if tsv_map:
            if verbose:
                _info(f"Applying TSV overrides to original reference array...")

            nentries = get_tsv_reference_count(tsv_map)
            used_tsv = 0
            used_bam = 0

            for i in range(nentries):
                ref_name = tsv_map.entries[i].reference_name
                tsv_length = tsv_map.entries[i].reference_length

                tid32 = sam_hdr_name2tid(bam_header, ref_name)
                if 0 <= tid32 < total_references:
                    if tsv_length > 0:
                        reference_lengths[tid32] = tsv_length
                        used_tsv += 1
            if verbose:
                _info(f"TSV overrides applied to original array:")
                _info(f"  TSV overrides applied: {used_tsv}")
                _info(f"  BAM defaults used: {total_references - used_tsv}")

        reference_alignment_counts = <int64_t*>malloc(total_references * sizeof(int64_t))
        if not reference_alignment_counts:
            raise MemoryError("Failed to allocate reference alignment counts")

        for reference_idx in range(total_references):
            result_code = hts_idx_get_stat(bam_index, reference_idx, &mapped_alignments, &unmapped_alignments)
            reference_alignment_counts[reference_idx] = <int64_t>mapped_alignments if result_code == 0 else 0

        references_to_process = 0
        for reference_idx in range(total_references):
            if reference_alignment_counts[reference_idx] >= min_read_count:
                references_to_process += 1

        if references_to_process <= 0:
            raise RuntimeError(f"No references with >= {min_read_count} alignments found")

        reference_ids = <int64_t*>malloc(references_to_process * sizeof(int64_t))
        if not reference_ids:
            raise MemoryError("Failed to allocate reference IDs")

        reference_idx = 0
        for ref_loop in range(total_references):
            if reference_alignment_counts[ref_loop] >= min_read_count:
                reference_ids[reference_idx] = ref_loop
                reference_idx += 1

        # Create filtered reference lengths array
        filtered_reference_lengths = <int64_t*>malloc(references_to_process * sizeof(int64_t))
        if not filtered_reference_lengths:
            raise MemoryError("Failed to allocate filtered reference lengths")

        # Copy lengths for kept references only (maintains TSV overrides)
        reference_idx = 0
        for ref_loop in range(total_references):
            if reference_alignment_counts[ref_loop] >= min_read_count:
                # Copy from original array (which has TSV overrides applied)
                filtered_reference_lengths[reference_idx] = reference_lengths[ref_loop]
                reference_idx += 1

        # Validate filtered reference lengths
        with nogil:
            min_filt_len = filtered_reference_lengths[0] if references_to_process > 0 else 0
            max_filt_len = min_filt_len
            total_filt_len = 0
            
            for i in range(references_to_process):
                if filtered_reference_lengths[i] <= 0 or filtered_reference_lengths[i] > 1000000000:
                    filtered_reference_lengths[i] = 1000  # Fix garbage values
                    invalid_filtered_lengths += 1
                
                if filtered_reference_lengths[i] < min_filt_len:
                    min_filt_len = filtered_reference_lengths[i]
                if filtered_reference_lengths[i] > max_filt_len:
                    max_filt_len = filtered_reference_lengths[i]
                total_filt_len += filtered_reference_lengths[i]

        if verbose:
            _info(f"Created filtered reference lengths:")
            _info(f"  Original references: {total_references}")
            _info(f"  Filtered references: {references_to_process}")
            _info(f"  Length range: [{min_filt_len}, {max_filt_len}]")
            _info(f"  Total length: {total_filt_len}")
            _info(f"  Invalid original lengths fixed: {invalid_original_lengths}")
            _info(f"  Invalid filtered lengths fixed: {invalid_filtered_lengths}")

        # Calculate total alignment count
        total_alignments = 0
        for reference_idx in range(references_to_process):
            total_alignments += reference_alignment_counts[reference_ids[reference_idx]]

        if verbose:
            _info(f"Analysis complete:")
            _info(f"  Total references: {total_references}")
            _info(f"  References to process: {references_to_process}")
            _info(f"  Total alignments: {total_alignments}")
            _info(f"  Reference lengths: FILTERED array created with valid data")

        # Create balanced batches
        max_batches_temp = 1024
        if references_to_process < max_batches_temp:
            max_batches_temp = references_to_process
        max_batches = max_batches_temp
        batch_starts = <int64_t*>malloc(max_batches * sizeof(int64_t))
        batch_ends = <int64_t*>malloc(max_batches * sizeof(int64_t))
        if not batch_starts or not batch_ends:
            raise MemoryError("Failed to allocate batch boundary arrays")

        batch_count = create_balanced_batches_greedy(reference_ids, reference_alignment_counts, 
                                                   references_to_process, max_batches, 
                                                   batch_starts, batch_ends, num_threads, verbose)

        if batch_count <= 0:
            free(batch_starts)
            free(batch_ends)
            raise RuntimeError("No batches created by create_balanced_batches_greedy")

        # Create processing batches
        processing_batches = <ProcessingBatch**>malloc(batch_count * sizeof(ProcessingBatch*))
        if not processing_batches:
            free(batch_starts)
            free(batch_ends)
            raise MemoryError("Failed to allocate batch array")

        if verbose:
            _info(f"Creating {batch_count} batches for parallel processing")

        for batch_idx in range(batch_count):
            batch_start = batch_starts[batch_idx]
            batch_end = batch_ends[batch_idx]
            expected_alignments = 0
            for ref_idx in range(batch_start, batch_end):
                expected_alignments += reference_alignment_counts[reference_ids[ref_idx]]
            processing_batches[batch_idx] = create_processing_batch(
                batch_idx, batch_start, batch_end, expected_alignments)
            if not processing_batches[batch_idx]:
                free(batch_starts)
                free(batch_ends)
                raise MemoryError(f"Failed to create processing batch {batch_idx}")

        free(batch_starts)
        free(batch_ends)
        
        # Pre-open thread file handles
        thread_handles = <samFile**>malloc(num_threads_c * sizeof(samFile*))
        if not thread_handles:
            raise MemoryError("Failed to allocate thread file handles")

        for thread_idx in range(num_threads_c):
            thread_handles[thread_idx] = NULL

        for thread_idx in range(num_threads_c):
            thread_handles[thread_idx] = hts_open(bam_file_path, b"r")
            if not thread_handles[thread_idx]:
                raise RuntimeError(f"Failed to open BAM file handle for thread {thread_idx}")

        # Create thread-local hash maps
        thread_maps = <ThreadLocalHashMap**>malloc(num_threads_c * sizeof(ThreadLocalHashMap*))
        if not thread_maps:
            raise MemoryError("Failed to allocate thread hash map array")

        for thread_idx in range(num_threads_c):
            thread_maps[thread_idx] = create_thread_local_hash_map()
            if not thread_maps[thread_idx]:
                raise MemoryError(f"Failed to create thread hash map {thread_idx}")

        # Create PMD context for stats collection (if PMD scoring is enabled)
        collect_pmd_stats = scoring_config.calculate_pmd

        if collect_pmd_stats:
            with nogil:
                pmd_context = create_pmd_context(num_threads_c, scoring_config.is_single_stranded, hierarchical_pmd_c)
            if not pmd_context:
                raise MemoryError("Failed to create PMD context")
            if verbose:
                _info("PMD statistics collection enabled")

        # Process batches in parallel
        _announce_stage("Alignment Processing", "Loading alignments, applying quality filters, and computing alignment scores")
        stage_timer = bf_monotonic_seconds()
        if verbose:
            _info("Processing batches")

        with nogil:
            for batch_idx in prange(batch_count, num_threads=num_threads_c, schedule='static'):
                thread_idx = threadid()
                process_batch_alignments_with_pmd(
                    thread_handles[thread_idx], bam_header, bam_index,
                    reference_ids, processing_batches[batch_idx], &scoring_config,
                    thread_maps[thread_idx], thread_idx, pmd_context
                )

        _log_stage("Batch processing", stage_timer)

        # Close all file handles immediately
        for thread_idx in range(num_threads_c):
            if thread_handles[thread_idx]:
                hts_close(thread_handles[thread_idx])
                thread_handles[thread_idx] = NULL
        free(thread_handles)
        thread_handles = NULL

        # Close BAM index immediately
        if bam_index:
            hts_idx_destroy(bam_index)
            bam_index = NULL

        # Assign global sequential IDs
        if verbose:
            _info(f"Assigning global sequential IDs")

        with nogil:
            if assign_global_sequential_ids_fast(processing_batches, batch_count, thread_maps, num_threads_c) != 0:
                raise RuntimeError("Failed to assign global sequential IDs")

        # Count results BEFORE streaming starts
        if verbose:
            _info(f"Counting results before streaming cleanup")

        with nogil:
            actual_total_alignments = count_actual_alignments(processing_batches, batch_count)
            unique_read_count = count_unique_reads_from_thread_maps(thread_maps, num_threads_c)
            unique_reference_count = count_unique_refs_from_batches(processing_batches, batch_count)

            # Clean up thread maps immediately after counting
            for thread_idx in range(num_threads_c):
                if thread_maps[thread_idx]:
                    destroy_thread_local_hash_map(thread_maps[thread_idx])
                    thread_maps[thread_idx] = NULL
            free(thread_maps)
            thread_maps = NULL

        if actual_total_alignments == 0:
            raise RuntimeError("No valid alignments found after processing")

        # Quick streaming statistics calculation
        if verbose:
            _info(f"Calculating statistics before cleanup")

        with nogil:
            parallel_streaming_stats_optimized(processing_batches, batch_count, &min_score, &max_score,
                                              &mean_score, &variance_score, &n_score, num_threads_c)
        if verbose:
            _info(f"Statistics complete:")
            _info(f"  Valid alignments: {actual_total_alignments:,}")
            _info(f"  Unique reads: {unique_read_count:,}")
            _info(f"  Unique references: {unique_reference_count:,}")
            std_score = libc_sqrt(variance_score) if variance_score > 0.0 else 0.0
            _info(f"  Score stats: min={min_score:.4f}, max={max_score:.4f}, mean={mean_score:.4f}, std={std_score:.4f}")

        # Finalize PMD model if stats were collected
        if pmd_context != NULL:
            if verbose:
                _info("Finalizing PMD damage model from collected statistics")
            with nogil:
                # Merge thread-local stats into global
                merge_pmd_stats(pmd_context)
                # Fit D(z) curve with regularization (use default params)
                pmd_params.P_prior_mean = 0.3
                pmd_params.P_prior_sd = 0.15
                pmd_params.lambda_prior_mean = 0.35
                pmd_params.lambda_prior_sd = 0.2
                pmd_params.C_prior_mean = 0.01
                pmd_params.C_prior_sd = 0.005
                pmd_params.omega_alpha = 0.5
                pmd_params.omega_beta = 5.0
                pmd_params.epsilon = 0.01
                fit_pmd_curve(pmd_context, &pmd_params)
                # Freeze model for BAM writing
                finalize_pmd_model(pmd_context)
            if verbose:
                _info(f"  PMD omega (damage presence): {pmd_context.model.curve.omega:.4f}")
                _info(f"  PMD decay parameter: {pmd_context.model.curve.lambda_decay:.4f}")
                _info(f"  Alignments analyzed: {pmd_context.model.stats.total_alignments:,}")

        # Apply EM algorithm
        if use_em:
            # Determine ANI threshold from scoring config
            min_ani_threshold = scoring_config.minimum_read_identity if scoring_config.minimum_read_identity > 0 else 90.0

            # Get PMD curve pointer if available
            if pmd_context != NULL and pmd_context.model.finalized:
                curve_ptr = &pmd_context.model.curve

            # Count alignments passing ANI filter (uses PMD curve if available)
            if verbose:
                _info("Counting alignments passing ANI filter")
            _announce_stage("ANI Filter", "Computing corrected ANI and counting passing alignments")
            stage_timer = bf_monotonic_seconds()

            with nogil:
                filtered_alignment_count = count_alignments_passing_ani_filter(
                    processing_batches, batch_count,
                    curve_ptr, min_ani_threshold, c_epsilon
                )
            _log_stage("ANI filter count", stage_timer)

            if verbose:
                _info(f"  Alignments passing ANI >= {min_ani_threshold:.1f}%: {filtered_alignment_count:,} / {actual_total_alignments:,}")
                _info(f"  Memory savings: {100.0 * (1.0 - <double>filtered_alignment_count / <double>actual_total_alignments):.1f}%")

            if filtered_alignment_count == 0:
                raise RuntimeError("No alignments pass the ANI filter - check threshold setting")

            # Create memory pool sized for filtered alignments only
            # Use split arrays for reduced memory usage
            if verbose:
                _info(f"Creating memory pool for {filtered_alignment_count:,} filtered alignments (split arrays)")

            memory_pool = create_memory_pool_split(
                filtered_alignment_count, unique_reference_count, unique_read_count,
                filtered_reference_lengths, calculate_pmd,
                hierarchical_pmd,  # enable hierarchical arrays
                hierarchical_pmd,  # enable damage counts (needed for gamma update)
                num_threads_c)
            if not memory_pool:
                raise MemoryError("Failed to create memory pool")

            # Store PMD curve pointer for BAM writing
            if curve_ptr != NULL:
                memory_pool.pmd_curve_ptr = <void*>curve_ptr
            else:
                memory_pool.pmd_curve_ptr = NULL

            if memory_pool.stats != NULL:
                memory_pool.stats.taxonomy_enabled = 1 if taxonomy_filter_config.enabled else 0

            # Transfer only filtered alignments to pool
            if verbose:
                _info("Transferring filtered alignments to memory pool")

            _announce_stage("Memory Optimization", "Transferring ANI-filtered alignments to optimized memory structures")
            stage_timer = bf_monotonic_seconds()

            with nogil:
                if memory_pool.use_split_arrays:
                    if populate_memory_pool_filtered_split(memory_pool, processing_batches, batch_count,
                                                           bam_header, curve_ptr, min_ani_threshold,
                                                           c_epsilon, num_threads_c) != 0:
                        with gil:
                            raise RuntimeError("Failed to stream filtered batches to memory pool (split)")
                else:
                    if populate_memory_pool_filtered(memory_pool, processing_batches, batch_count,
                                                      bam_header, curve_ptr, min_ani_threshold,
                                                      c_epsilon, num_threads_c) != 0:
                        with gil:
                            raise RuntimeError("Failed to stream filtered batches to memory pool")
            _log_stage("Filtered pool streaming", stage_timer)

            if verbose:
                _info(f"  Final pool: {memory_pool.alignment_count:,} alignments, {memory_pool.final_unique_reads:,} unique reads")

            # CRITICAL: Remap reference indices to compact form BEFORE EM runs
            # Alignments from batches use original BAM header indices [0, n_original_refs-1]
            # But pool.reference_count is set to unique_reference_count (refs with alignments)
            # EM uses pool.reference_count as n_refs, so indices must be in [0, reference_count-1]
            #
            # We create the mapping once here and KEEP IT for later use (after probability filtering).
            # This ensures new_to_old_tid always maps to original BAM header indices.
            # After probability filtering, we use update_reference_mapping_after_filtering()
            # instead of creating a new mapping.
            mapping = create_reference_mapping(memory_pool, bam_header)
            if not mapping:
                raise RuntimeError("Failed to create initial reference mapping")

            if remap_alignment_reference_ids(memory_pool, mapping) != 0:
                destroy_reference_mapping(mapping)
                mapping = NULL
                raise RuntimeError("Failed to remap alignment reference IDs for EM")

            # Update reference count to match the compact mapping
            memory_pool.reference_count = mapping.n_retained_refs
            if verbose:
                _info(f"  Remapped {mapping.n_retained_refs} reference indices for EM")

            # Update initial stats (use total_references from BAM header, not just refs with alignments)
            with nogil:
                update_initial_stats(memory_pool.stats, memory_pool.alignment_count,
                                   unique_read_count, total_references)

            # Initialize hierarchical EM for ancient/modern classification if enabled and PMD available
            if pmd_context != NULL and pmd_context.model.finalized and hierarchical_pmd:
                if verbose:
                    _info("Initializing hierarchical EM for ancient/modern classification")
                init_hierarchical_em_py(
                    <uintptr_t>memory_pool,
                    <uintptr_t>&pmd_context.model.curve,
                    c_epsilon
                )

            # Update quality filter stats (now using ANI-filtered count)
            with nogil:
                update_quality_filter_stats(memory_pool.stats, memory_pool.alignment_count,
                                          unique_read_count, unique_reference_count)

            # Batch array cleanup (contents already freed by filtered streaming)
            if processing_batches:
                free(processing_batches)
                processing_batches = NULL

            # Initialize Coverage-Weighted Reference Priors if enabled
            if c_cwrp_lambda > 0.0:
                if verbose:
                    _info(f"Initializing CWRP with lambda={c_cwrp_lambda:.3f} iterative={c_iterative_auth}")
                with nogil:
                    if init_cwrp(memory_pool, c_cwrp_lambda, c_iterative_auth,
                                 c_auth_update_interval, c_damage_weight,
                                 c_low_cov_floor, c_low_cov_shrink_tau) != 0:
                        with gil:
                            raise RuntimeError("CWRP initialization failed")

            # Set Posterior-Weighted Coverage Authenticity (Path B) parameters on pool
            # These are read by execute_em_py to configure the EM loop
            with nogil:
                memory_pool.auth_update_interval_post = c_auth_update_interval_post if c_auth_post_enabled else 0
                memory_pool.auth_scale_post = c_auth_scale_post
                memory_pool.auth_lambda_ramp_iters = c_auth_lambda_ramp_iters
                memory_pool.sample_pi_override = c_sample_pi_override

            # Execute EM algorithm before graph analysis
            if verbose:
                _info("Running EM algorithm")

            # EM algorithm: phi-space optimization with SQUAREM acceleration
            # Lazy import to avoid circular dependency (processor_em cimports from processor)
            from .processor_em import execute_em_py
            _announce_stage("EM Optimization", "Running Expectation-Maximization algorithm with SQUAREM acceleration")
            stage_timer = bf_monotonic_seconds()
            if execute_em_py(
                <uintptr_t>memory_pool,
                max_em_iterations,
                em_tolerance,
                prior_weight,
                em_power_rho,
                em_unknown_component,
                em_unknown_margin,
                hierarchical_pmd,
                0.02,  # D_avg_5p (default)
                0.02,  # D_avg_3p (default)
                0.01,  # epsilon_error
                use_squarem_acceleration,
                squarem_start_iter,
                enable_globalization,
                backtrack_factor,
                max_backtrack_steps,
                steplength_scheme,
                num_threads,
            ) != 0:
                raise RuntimeError("EM algorithm execution failed")
            _log_stage("EM algorithm", stage_timer)

            # Update EM stats
            with nogil:
                update_em_stats(memory_pool.stats, memory_pool.iteration_count,
                              memory_pool.algorithm_converged, memory_pool.final_log_likelihood)

            if verbose:
                _info(f"EM converged after {memory_pool.iteration_count} iterations")
                _info("Applying probability filtering")

            _announce_stage("Confidence Filtering", "Removing low-confidence alignments based on posterior probabilities")
            stage_timer = bf_monotonic_seconds()
            if apply_probability_filtering(memory_pool, &em_config) != 0:
                raise RuntimeError("Probability filtering failed")
            _log_stage("Probability filtering", stage_timer)

            # Update probability filter stats
            with nogil:
                update_probability_filter_stats(memory_pool.stats, memory_pool.alignment_count,
                                              memory_pool.final_unique_reads, unique_reference_count)

            # Free EM intermediate memory before graph analysis to reduce peak memory
            with nogil:
                cleanup_em_intermediate_memory(memory_pool)

            if verbose:
                _info(f"Analyzing reference graph...")

            # Update the existing mapping to account for refs that lost all alignments during filtering.
            # This preserves the original BAM header index mapping while compacting to eliminate gaps.
            # The mapping was created before EM runs and contains original->compact1 indices.
            # update_reference_mapping_after_filtering updates it to original->compact2 indices.
            if mapping != NULL:
                if update_reference_mapping_after_filtering(mapping, memory_pool) != 0:
                    destroy_reference_mapping(mapping)
                    mapping = NULL
                    raise RuntimeError("Failed to update reference mapping after filtering")

                memory_pool.reference_count = mapping.n_retained_refs
                _debug(1, f"APPLY MAPPING: Updated memory_pool.reference_count -> {mapping.n_retained_refs}")
                if memory_pool.stats != NULL:
                    memory_pool.stats.post_probability_references = mapping.n_retained_refs
            else:
                raise RuntimeError("Reference mapping is NULL after probability filtering")

            # Run Bayesian damage model AFTER remapping to compute damage metrics per reference
            # Uses EM phi values as weights for damage count accumulation
            #
            # NOTE: The EM hierarchical gamma values are PRESERVED and remapped to new indices.
            # These provide superior estimates via information pooling across references.
            # The damage model computes amplitude, baseline, and log_bf for each reference,
            # but does NOT overwrite gamma_values.
            if calculate_pmd:
                _announce_stage("Damage Model", "Computing per-reference P(ancient) using unified damage model")
                stage_timer = bf_monotonic_seconds()
                from .processor_damage_model import (
                    allocate_damage_stats_py, free_damage_stats_py,
                    accumulate_from_em_py
                )
                from .unified_damage import fit_unified_from_damage_stats_py
                try:
                    # Remap EM gamma values to new reference indices (preserves hierarchical shrinkage)
                    old_gamma = memory_pool.gamma_values

                    if old_gamma != NULL and mapping != NULL:
                        new_gamma = <double*>calloc(memory_pool.reference_count, sizeof(double))
                        if new_gamma != NULL:
                            for new_idx in range(mapping.n_retained_refs):
                                old_idx = mapping.new_to_old_tid[new_idx]
                                if old_idx < mapping.n_original_refs:
                                    new_gamma[new_idx] = old_gamma[old_idx]
                            free(old_gamma)
                            memory_pool.gamma_values = new_gamma
                            _debug(1, f"Remapped {mapping.n_retained_refs} gamma values from EM")
                        else:
                            _warn("Failed to allocate remapped gamma array, keeping old")
                    elif old_gamma == NULL:
                        memory_pool.gamma_values = <double*>calloc(memory_pool.reference_count, sizeof(double))

                    # Allocate arrays for damage model results
                    if memory_pool.damage_amplitude != NULL:
                        free(memory_pool.damage_amplitude)
                    if memory_pool.damage_baseline != NULL:
                        free(memory_pool.damage_baseline)
                    if memory_pool.damage_log_bf != NULL:
                        free(memory_pool.damage_log_bf)
                    memory_pool.damage_amplitude = <double*>calloc(memory_pool.reference_count, sizeof(double))
                    memory_pool.damage_baseline = <double*>calloc(memory_pool.reference_count, sizeof(double))
                    memory_pool.damage_log_bf = <double*>calloc(memory_pool.reference_count, sizeof(double))
                    if memory_pool.gamma_values == NULL or memory_pool.damage_amplitude == NULL or memory_pool.damage_baseline == NULL or memory_pool.damage_log_bf == NULL:
                        raise MemoryError("Failed to allocate damage model arrays")

                    # Accumulate damage counts from alignments weighted by EM posteriors
                    damage_stats_ptr = allocate_damage_stats_py(memory_pool.reference_count)
                    is_ss = library_type == "ss"
                    accumulate_from_em_py(<uintptr_t>memory_pool, damage_stats_ptr, is_ss, num_threads)

                    # Fit unified damage model (data-driven tau, Gamma shrinkage for amplitudes)
                    damage_result = fit_unified_from_damage_stats_py(
                        <uintptr_t>memory_pool, damage_stats_ptr, is_ss
                    )
                    if verbose:
                        _info(f"Unified damage model: tau={damage_result['tau']:.2f}, {damage_result['n_ancient']}/{damage_result['n_fitted']} ancient (p>0.5)")
                    free_damage_stats_py(damage_stats_ptr)
                except Exception as e:
                    _warn(f"Damage model failed: {e}")
                _log_stage("Unified damage model", stage_timer)

            pattern_data = <ReferencePattern*>calloc(memory_pool.reference_count, sizeof(ReferencePattern))
            if not pattern_data:
                destroy_reference_mapping(mapping)
                raise MemoryError("Failed to allocate pattern_data array")

            if tsv_file_path_c:
                ref_stats = <ReferenceStats*>malloc(memory_pool.reference_count * sizeof(ReferenceStats))
                if not ref_stats:
                    free(pattern_data)
                    destroy_reference_mapping(mapping)
                    raise MemoryError("Failed to allocate ref_stats for graph analysis")

                with nogil:
                    calculate_reference_stats(memory_pool, bam_header, ref_stats)

                with nogil:
                    read_index = build_read_index_parallel(memory_pool, memory_pool.reference_count, num_threads_c)

                if not read_index:
                    free(ref_stats)
                    free(pattern_data)
                    destroy_reference_mapping(mapping)
                    raise RuntimeError("Failed to build read index for graph analysis")

                if graph_min_edge_weight == 0:
                    if verbose:
                        bf_logging.log("GRAPH", f"Auto edge threshold: running elbow detection (min_read_count={min_read_count})")
                    try:
                        with nogil:
                            thr_c = pick_min_edge_weight_elbow(<MemoryPool*>memory_pool, <ReadIndex*>read_index, <ReferenceStats*>ref_stats, <uint32_t>memory_pool.reference_count, min_read_count_c, 0.0, verbose_c, 0.0, 0)
                        graph_min_edge_weight_c = thr_c
                        if verbose:
                            bf_logging.log("GRAPH", f"Elbow threshold = {graph_min_edge_weight_c}")
                    except Exception as e:
                        if verbose:
                            bf_logging.warn(f"GRAPH: Edge threshold selection failed; proceeding with threshold={graph_min_edge_weight_c}. Error: {e}")

                _announce_stage("Connectivity Analysis", "Constructing reference connectivity graph and analyzing read categories")
                stage_timer = bf_monotonic_seconds()
                # Build graph if clustering is enabled OR GMRF smoothing is requested
                filtered_graph = analyze_reference_graph(memory_pool, pattern_data, em_config.minimum_read_coverage,
                                          &em_config, bam_header, mapping, verbose, clustering or enable_gmrf_c, tsv_file_path_c,
                                          graph_min_edge_weight_c, NULL, read_index)
                _log_stage("Reference graph analysis", stage_timer)

                # GMRF smoothing for reference-level ancientness estimates (profile-only)
                # Requires nodes array to be populated (done in analyze_reference_graph)
                if enable_gmrf_c and filtered_graph != NULL and memory_pool.gamma_values != NULL:
                    if filtered_graph.nodes != NULL and ref_stats != NULL:
                        _announce_stage("GMRF Smoothing", "Smoothing reference ancientness using graph structure")
                        gmrf_stage_timer = bf_monotonic_seconds()

                        gmrf_n_reads = <uint32_t*>malloc(
                            memory_pool.reference_count * sizeof(uint32_t))
                        if gmrf_n_reads != NULL:
                            for gmrf_i in range(<int32_t>memory_pool.reference_count):
                                gmrf_n_reads[gmrf_i] = ref_stats[gmrf_i].total_reads

                            with nogil:
                                gmrf_result = apply_gmrf_smoothing_to_graph(
                                    filtered_graph,
                                    memory_pool.gamma_values,
                                    memory_pool.damage_log_bf,
                                    gmrf_n_reads,
                                    <uint32_t>memory_pool.reference_count,
                                    gmrf_tau_c,
                                    100,   # max_iter
                                    1e-6,  # tol
                                    verbose_c,
                                )

                            free(gmrf_n_reads)
                            gmrf_n_reads = NULL

                            if gmrf_result == 0:
                                _log_stage("GMRF smoothing", gmrf_stage_timer)
                            else:
                                if verbose:
                                    bf_logging.warn("GMRF smoothing failed, using unsmoothed values")

                # Taxonomy-aware graph analysis (if databases provided)
                if taxonomy_db is not None and taxonomy_accession_map is not None:
                    _announce_stage("Taxonomy Enrichment", "Enriching graph with taxonomic information")
                    taxonomy_stage_timer = bf_monotonic_seconds()

                    from bam_filter.taxonomy_db import TaxonomyDatabase, load_accession_map_from_file

                    _info("Loading taxonomy database from %s", taxonomy_db)
                    try:
                        taxdb_obj = TaxonomyDatabase.from_parquet(taxonomy_db)
                        _info("  Loaded taxonomy: %d nodes", taxdb_obj.n_nodes)
                    except Exception as e:
                        _error("Failed to load taxonomy database: %s", str(e))
                        raise

                    taxdb_c = taxdb_obj.db

                    _info("Extracting reference accessions for filtered loading...")
                    reference_accessions = []
                    for ref_idx in range(memory_pool.reference_count):
                        if mapping != NULL and mapping.new_to_old_tid != NULL and ref_idx < mapping.n_retained_refs:
                            orig_tid = mapping.new_to_old_tid[ref_idx]
                        else:
                            orig_tid = ref_idx
                        ref_name = sam_hdr_tid2name(bam_header, orig_tid)
                        if ref_name != NULL:
                            ref_name_str = ref_name.decode('utf-8')
                            if len(ref_name_str) > 0:
                                tokens = ref_name_str.split()
                                if len(tokens) > 0:
                                    reference_accessions.append(tokens[0])

                    _info("Loading accession map for %d references...", len(reference_accessions))

                    try:
                        accmap_obj = load_accession_map_from_file(taxonomy_accession_map, accession_filter=reference_accessions)
                        _info("Accession map loaded successfully")
                    except Exception as e:
                        _error("Failed to load accession map: %s", str(e))
                        raise

                    accmap_c = accmap_obj.amap

                    tax_config.enabled = True
                    tax_config.min_rank_id_for_comparison = taxonomy_min_rank
                    tax_config.cross_domain_threshold = taxonomy_cross_domain_threshold
                    tax_config.kingdom_mismatch_threshold = taxonomy_kingdom_threshold
                    tax_config.genus_mismatch_threshold = taxonomy_genus_threshold

                    if verbose:
                        _info("Enriching references with taxonomy IDs...")

                    with nogil:
                        enrich_patterns_with_taxonomy(
                            pattern_data,
                            memory_pool.reference_count,
                            bam_header,
                            taxdb_c,
                            accmap_c,
                            mapping,
                            verbose_c  # Use C int version
                        )

                    _log_stage("Taxonomy ID enrichment", taxonomy_stage_timer)

                if clustering and not filtered_graph:
                    free(pattern_data)
                    destroy_reference_mapping(mapping)
                    raise RuntimeError("Reference pattern detection failed")

                if memory_pool.stats != NULL:
                    with nogil:
                        update_graph_stats(memory_pool.stats, memory_pool.reference_count, memory_pool.reference_count)

                verbose_c = 1 if verbose else 0
                alignments_before_filtering = memory_pool.alignment_count
            else:
                _announce_stage_skip("Connectivity Analysis", "Graph analysis disabled")
                filtered_graph = NULL
                verbose_c = 1 if verbose else 0
                alignments_before_filtering = memory_pool.alignment_count
                free(pattern_data)
                pattern_data = NULL
                if memory_pool.stats != NULL:
                    memory_pool.stats.graph_analysis_references = -1
                    memory_pool.stats.graph_patterns_computed = 0

            aligns_filtered_information = 0
            if clustering:
                if filtered_graph and filtered_graph.num_nodes >= 10000:
                    _announce_stage("Graph Filtering", "Applying cluster-aware refinement using Label Propagation Algorithm (fast, for large graphs)")
                else:
                    _announce_stage("Graph Filtering", "Applying cluster-aware refinement using Community algorithm")

                if not filtered_graph:
                    raise RuntimeError("filtered_graph is NULL at Phase 6 start")
                if not filtered_graph.igraph_handle:
                    raise RuntimeError("filtered_graph.igraph_handle is NULL at Phase 6 start")
                if not filtered_graph.weights_handle:
                    raise RuntimeError("filtered_graph.weights_handle is NULL at Phase 6 start")
                if not filtered_graph.tsv_exact_connection_counts:
                    raise RuntimeError("filtered_graph.tsv_exact_connection_counts is NULL")
                if not filtered_graph.tsv_co_mapping_averages:
                    raise RuntimeError("filtered_graph.tsv_co_mapping_averages is NULL")
                if not filtered_graph.tsv_max_co_mappings:
                    raise RuntimeError("filtered_graph.tsv_max_co_mappings is NULL")
                if not filtered_graph.tsv_neighbor_multimap_avg:
                    raise RuntimeError("filtered_graph.tsv_neighbor_multimap_avg is NULL")

                if not ref_stats:
                    ref_stats = <ReferenceStats*>malloc(memory_pool.reference_count * sizeof(ReferenceStats))
                    if not ref_stats:
                        free(pattern_data)
                        destroy_reference_mapping(mapping)
                        raise MemoryError("Failed to allocate ref_stats for cluster filtering")

                    with nogil:
                        calculate_reference_stats(memory_pool, bam_header, ref_stats)

                if not read_index:
                    with nogil:
                        read_index = build_read_index_parallel(memory_pool, memory_pool.reference_count, num_threads_c)

                    if not read_index:
                        if ref_stats:
                            free(ref_stats)
                            ref_stats = NULL
                        free(pattern_data)
                        destroy_reference_mapping(mapping)
                        raise RuntimeError("Failed to build read index for cluster filtering")

                # If user selected auto (0), ensure we have a concrete threshold.
                # Prefer the earlier selection (if we ran Otsu prior to igraph build).
                # Only run Otsu here if the threshold is still the 0 sentinel.
                if graph_min_edge_weight == 0:
                    if graph_min_edge_weight_c == 0:
                        if verbose:
                            bf_logging.log("GRAPH", f"Auto edge threshold: running elbow detection (min_read_count={min_read_count})")
                        try:
                            with nogil:
                                thr_c = pick_min_edge_weight_elbow(<MemoryPool*>memory_pool, <ReadIndex*>read_index, <ReferenceStats*>ref_stats, <uint32_t>memory_pool.reference_count, min_read_count_c, 0.0, verbose_c, 0.0, 0)
                            graph_min_edge_weight_c = thr_c
                            if verbose:
                                bf_logging.log("GRAPH", f"Elbow threshold = {graph_min_edge_weight_c}")
                        except Exception as e:
                            if verbose:
                                bf_logging.warn(f"GRAPH: Edge threshold selection failed; keeping threshold={graph_min_edge_weight_c}. Error: {e}")
                    else:
                        if verbose:
                            bf_logging.log("GRAPH", f"Using previously selected auto threshold = {graph_min_edge_weight_c}")

                use_community_c = 1 if clustering else 0
                community_res_c = community_resolution if community_resolution else 1.0
                community_parallel_c = 0
                community_max_iter_c = community_max_iterations if community_max_iterations else 10

                stage_timer = bf_monotonic_seconds()

                if not memory_pool:
                    raise RuntimeError("memory_pool is NULL before clustering")
                if not pattern_data:
                    raise RuntimeError("pattern_data is NULL before clustering")
                if not ref_stats:
                    raise RuntimeError("ref_stats is NULL before clustering")
                if not read_index:
                    raise RuntimeError("read_index is NULL before clustering")
                taxonomy_stats.strict_removed = 0
                taxonomy_stats.weighted_count = 0
                taxonomy_stats.second_chance_restored = 0

                remove_cross_domain_alignments_c = 1 if remove_cross_domain_alignments else 0
                detect_misannotations_c = 1 if detect_misannotations else 0

                cluster_result = apply_cluster_aware_filtering(
                    memory_pool, pattern_data, ref_stats, read_index,
                    memory_pool.reference_count, em_config.minimum_read_coverage,
                    verbose_c,
                    community_res_c, community_max_iter_c,
                    graph_min_edge_weight_c,
                    num_threads_c,
                    outlier_method_c,  # Outlier detection method (0=MAD, 1=IQR)
                    filtered_graph,  # Pass the pre-built filtered graph
                    bam_header,       # For TSV writing
                    mapping,          # For TSV writing
                    tsv_file_path_c,  # TSV will be written after Community completes
                    graph_export_path_c,  # GraphML export path
                    &taxonomy_filter_config,  # Taxonomy filter config
                    &taxonomy_stats,  # Output: taxonomy filtering statistics
                    taxdb_c,  # Taxonomy database for lineage information
                    0.01,  # betweenness_threshold: bridge detection
                    0.3,   # cc_threshold: hub detection
                    5,     # hub_degree_threshold: minimum degree for hub/core
                    False,  # strict_mode: non-strict by default
                    remove_cross_domain_alignments_c,  # Alignment removal toggle
                    detect_misannotations_c,  # Misannotation detection toggle
                    # Network QC filtering parameters
                    &network_qc_config if network_qc_filter else NULL,  # Network QC config (NULL = disabled)
                    tax_ambiguity_removal_level_c  # Tax ambiguity flag level for removal
                )

                with nogil:
                    destroy_read_index(read_index)
                    destroy_weighted_graph(filtered_graph)
                free(ref_stats)
                
                if cluster_result != 0:
                    free(pattern_data)
                    destroy_reference_mapping(mapping)
                    raise RuntimeError("Cluster-aware filtering failed")

                if memory_pool.stats != NULL:
                    memory_pool.stats.taxonomy_strict_removed = taxonomy_stats.strict_removed
                    memory_pool.stats.taxonomy_weighted_count = taxonomy_stats.weighted_count
                    memory_pool.stats.taxonomy_second_chance_restored = taxonomy_stats.second_chance_restored

                aligns_filtered_information = alignments_before_filtering - memory_pool.alignment_count

                free(pattern_data)
                pattern_data = NULL

                if update_reference_mapping_after_filtering(mapping, memory_pool) != 0:
                    destroy_reference_mapping(mapping)
                    raise RuntimeError("Failed to update reference mapping after filtering")

                if shrink_memory_pool(memory_pool) != 0:
                    raise RuntimeError("Failed to optimize memory pool after filtering")
                _log_stage("Cluster-aware filtering", stage_timer)
                
                if memory_pool.stats != NULL:
                    with nogil:
                        update_unified_filter_stats(memory_pool.stats, memory_pool.alignment_count,
                                                    memory_pool.final_unique_reads, mapping.n_retained_refs,
                                                    0, 0, 0,
                                                    0, aligns_filtered_information, 0)
            else:
                _announce_stage_skip("Advanced Filtering", "Cluster-aware filtering disabled")
                if filtered_graph:
                    with nogil:
                        destroy_weighted_graph(filtered_graph)
                    filtered_graph = NULL
                if read_index:
                    with nogil:
                        destroy_read_index(read_index)
                    read_index = NULL
                if ref_stats:
                    free(ref_stats)
                    ref_stats = NULL
                free(pattern_data)
                pattern_data = NULL
                if memory_pool.stats != NULL:
                    memory_pool.stats.post_unified_alignments = memory_pool.alignment_count
                    memory_pool.stats.post_unified_reads = memory_pool.final_unique_reads
                    memory_pool.stats.post_unified_references = mapping.n_retained_refs
                    memory_pool.stats.filtered_coverage_only = -1
                    memory_pool.stats.filtered_information_only = 0
                    memory_pool.stats.filtered_both_criteria = 0
                    memory_pool.stats.alignments_removed_coverage = 0
                    memory_pool.stats.alignments_removed_information = 0
                    memory_pool.stats.alignments_removed_both = 0

            if verbose:
                _info(f"EM algorithm completed:")
                _info(f"  Iterations: {memory_pool.iteration_count}")
                _info(f"  Converged: {memory_pool.algorithm_converged}")
                _info(f"  Final log-likelihood: {memory_pool.final_log_likelihood:.6f}")

            # Update final output stats with the compacted reference count
            with nogil:
                update_final_output_stats(memory_pool.stats, memory_pool.alignment_count,
                                        memory_pool.final_unique_reads, mapping.n_retained_refs)
            
            # Write output BAM if requested
            if output_bam:
                if verbose:
                    _info("Writing output BAM")
                    
                bam_write_threads = min_int32(num_threads_c, 4)
                _announce_stage("Output Generation", "Writing filtered alignments to optimized BAM file")
                stage_timer = bf_monotonic_seconds()
                if write_bam_with_filtered_header(memory_pool, bam_file_path, output_bam_path,
                                     bam_header, mapping, bam_write_threads) != 0:
                    raise RuntimeError("Failed to write filtered BAM file")
                _log_stage("BAM writing", stage_timer)

        # Calculate summary metrics and print comprehensive report
        if memory_pool and memory_pool.stats:
            with nogil:
                calculate_summary_metrics(memory_pool.stats)
                print_processing_stats(memory_pool.stats)
        
        # Update result
        if memory_pool:
            result = {
                'summary': {
                    'total_alignments': int(total_alignments),
                    'filtered_alignments': int(memory_pool.alignment_count),
                    'unique_reads': int(unique_read_count),
                    'final_unique_reads': int(memory_pool.final_unique_reads),
                    'unique_references': int(memory_pool.reference_count),  # Filtered count
                    'em_iterations': int(memory_pool.iteration_count),
                    'em_converged': bool(memory_pool.algorithm_converged),
                    'final_likelihood': float(memory_pool.final_log_likelihood)
                },
                'memory_optimized': True,
                'pmd_enabled': bool(memory_pool.pmd_enabled_for_output),
                'integrated_pmd': True,
                'reference_lengths_fixed': True,
                'tsv_overrides_applied': used_tsv if tsv_map else 0,
                'filtered_references': int(references_to_process),
                'library_type': 'ss' if scoring_config.is_single_stranded else 'ds',
                'zp_values_precomputed': bool(memory_pool.zp_values_computed)
            }
            
            # Add comprehensive processing stats if available
            if memory_pool.stats:
                result['processing_stats'] = {
                    'initial': {
                        'alignments': int(memory_pool.stats.initial_total_alignments),
                        'reads': int(memory_pool.stats.initial_total_reads),
                        'references': int(memory_pool.stats.initial_total_references)
                    },
                    'post_em': {
                        'iterations': int(memory_pool.stats.em_iterations),
                        'converged': bool(memory_pool.stats.em_converged),
                        'likelihood': float(memory_pool.stats.em_final_likelihood)
                    },
                    'post_probability_filter': {
                        'alignments': int(memory_pool.stats.post_probability_alignments),
                        'reads': int(memory_pool.stats.post_probability_reads),
                        'references': int(memory_pool.stats.post_probability_references),
                        'filtered': int(memory_pool.stats.filtered_probability_alignments)
                    },
                    'post_unified_filter': {
                        'alignments': int(memory_pool.stats.post_unified_alignments),
                        'reads': int(memory_pool.stats.post_unified_reads),
                        'references': int(memory_pool.stats.post_unified_references)
                    },
                    'final_output': {
                        'alignments': int(memory_pool.stats.final_alignments_written),
                        'reads': int(memory_pool.stats.final_reads_written),
                        'references': int(memory_pool.stats.final_references_written)
                    },
                    'retention': {
                        'alignment_retention': float(memory_pool.stats.overall_alignment_retention),
                        'read_retention': float(memory_pool.stats.overall_read_retention),
                        'reference_retention': float(memory_pool.stats.overall_reference_retention)
                    }
                }
            
            if tsv_file_path_c:
                result['reference_stats_tsv'] = tsv_file_bytes.decode('utf-8')

        return result

    finally:
        if thread_maps:
            for thread_idx in range(num_threads_c):
                if thread_maps[thread_idx]:
                    destroy_thread_local_hash_map(thread_maps[thread_idx])
            free(thread_maps)
        if processing_batches:
            for batch_idx in range(batch_count):
                if processing_batches[batch_idx]:
                    destroy_processing_batch(processing_batches[batch_idx])
            free(processing_batches)
        if reference_ids:
            free(reference_ids)
        if reference_alignment_counts:
            free(reference_alignment_counts)
        if reference_lengths:  # Original reference lengths
            free(reference_lengths)
        if filtered_reference_lengths:  # Filtered reference lengths
            free(filtered_reference_lengths)
        if tsv_map:  # Cleanup TSV map
            with nogil:
                free_tsv_reference_map(tsv_map)
        if bam_index:
            hts_idx_destroy(bam_index)
        if bam_header:
            sam_hdr_destroy(bam_header)
        if bam_handle:
            hts_close(bam_handle)
        # Clean up PMD context (must be after BAM writing completes)
        if pmd_context:
            with nogil:
                destroy_pmd_context(pmd_context)
        # Memory cleanup handled by caller
# ===============================================================================
# UPDATED PUBLIC INTERFACE FUNCTIONS
# ===============================================================================

def process_bam_with_em(
    # Core parameters
    bam_file,
    output_bam=None,
    num_threads=1,
    verbose=True,

    # PMD parameters
    calculate_pmd=False,
    library_type="ds",
    hierarchical_pmd=False,

    # TSV
    reference_lengths_tsv=None,
    # Graph analysis TSV export path (optional)
    reference_stats_tsv=None,

    # Read filtering
    min_read_count=1,
    min_read_length=30,
    max_read_length=10000,
    min_read_ani=0.0,

    # EM algorithm parameters
    max_em_iterations=25,
    em_tolerance=1e-5,
    min_probability=1e-6,
    prob_fraction=0.1,
    prior_weight=0.01,
    lambda_scale=0.3,
    use_squarem_acceleration=True,

    # SQUAREM parameters
    enable_globalization=True,
    squarem_start_iter=2,
    backtrack_factor=0.5,
    max_backtrack_steps=5,
    steplength_scheme=3,
    
    # Dominance regularization parameters
    enable_dominance_regularization=True,       # Master switch
    dominance_strength=None,           # Base penalty strength
    use_adaptive_dominance=True,           # Adapt to dataset entropy
    entropy_scaling_factor=1.0,             # How much entropy affects penalty
    min_penalty_strength=0.1,               # Minimum penalty strength
    max_penalty_strength=10.0,               # Maximum penalty strength
    entropy_confidence_threshold=0.0,        # Minimum entropy to apply adaptive
    enable_emergency_regularization=True,
    emergency_entropy_threshold=0.15,
    emergency_max_weight_threshold=0.3,
    emergency_likelihood_drop=1e-3,
    emergency_min_dominant_refs=0,
    init_prior_strength=0.1,
    information_threshold=-999.0,  # Information-theoretic filtering (-999.0 = disabled)

    # Tempered EM parameter (bias reduction)
    em_beta=1.0,  # Temperature: 1.0 = standard EM, 0.3-0.7 recommended for bias reduction

    # Effective length correction (reduces length bias)
    em_length_correction=False,  # Normalize by reference length: π_j ∝ counts_j / length_j

    # Unknown component (absorbs reads not belonging to any reference)
    em_unknown_component=False,  # Enable unknown/background component
    em_unknown_prior=0.05,       # Prior probability for unknown (5%)
    em_unknown_score=-50.0,      # Fixed log-likelihood score for unknown

    # Length-aware initialization
    em_length_init=False,  # Use length-weighted prior: π_j^(0) ∝ length_j

    # === UNIFIED φ-SPACE EM PARAMETERS (clean formulation) ===
    em_power_rho=1.0,              # ρ: 1.0 = standard, <1 reduces dominance (replaces dominance penalty)
    em_unknown_adaptive=False,     # Enable adaptive unknown score (vs fixed s_unknown)
    em_unknown_margin=2.0,         # Δ: margin below best score for unknown
    em_length_output_exp=1.0,      # γ_len: 0=none, 1=per-base abundance (π_j ∝ φ_j / L_j^γ_len)
    em_length_prior_exp=1.0,       # η_prior: 0=flat, 1=length-proportional (α_j = α0 × (L_j / mean_L)^η_prior)

    # Graph construction parameters (cluster-aware filtering always enabled)
    graph_min_edge_weight=0,  # Minimum edge weight to keep in graph (0=auto via elbow, -1=no filtering)

    # Clustering parameters
    clustering=False,  # Enable clustering/community detection (requires reference_stats_tsv)
    community_resolution=1.0,
    community_max_iterations=10,
    outlier_method="otsu",  # CC threshold method (Otsu bimodal separation)

    # Graph export
    graph_export=None,  # Export graph to GraphML format (optional)

    # Taxonomy parameters
    taxonomy_db=None,
    taxonomy_accession_map=None,
    taxonomy_min_rank=6,
    taxonomy_cross_domain_threshold=0.10,
    taxonomy_kingdom_threshold=0.25,
    taxonomy_genus_threshold=0.50,
    # Taxonomy-informed filtering parameters
    taxonomy_filter_enabled=False,
    taxonomy_strict_filter=True,
    taxonomy_strict_min_connections=5,
    taxonomy_weighted_outlier=True,
    taxonomy_anomaly_weight=2.0,
    taxonomy_second_chance=True,
    taxonomy_second_chance_cc=0.3,
    # Cross-domain removal options
    remove_cross_domain_alignments=False,
    remove_cross_domain_references=False,
    detect_misannotations=False,
    # Network QC & Taxonomic ambiguity detection parameters
    network_qc_filter=False,
    entropy_biased_threshold=1.0,
    entropy_mixed_threshold=2.0,
    entropy_highly_mixed_threshold=3.0,
    tax_ambiguity_removal_level=2,
    # Coverage-Weighted Reference Priors (CWRP)
    cwrp_lambda=0.0,
    iterative_auth=False,
    auth_update_interval=5,
    damage_weight=1.0,
    low_cov_floor=10,
    low_cov_shrink_tau=50.0,
    # Posterior-Weighted Coverage Authenticity (Path B)
    auth_post_enabled=False,
    auth_update_interval_post=3,
    auth_scale_post=4.0,
    auth_lambda_ramp_iters=5,
    # Sample-level P(ancient) gate
    sample_pi_override=0.0,
    # GMRF smoothing (profile-only)
    enable_gmrf=False,
    gmrf_tau=1.0,
):
    """
    High-level entry point to process a BAM file with the EM-based pipeline.

    Runs the full pipeline: optional PMD scoring, EM read reassignment (with
    SQUAREM acceleration), probability-based filtering, optional graph-based
    clustering, and writes a filtered BAM and optional TSV reports.

    Parameters
    ----------
    bam_file : str
        Path to input BAM file.
    output_bam : str, optional
        Path to write filtered BAM output. If None, filtered BAM is not written.
    num_threads : int
        Number of threads to use for parallel stages.
    calculate_pmd : bool
        If True, compute Post-Mortem Damage (PMD) scores and include them
        in scoring and output tags.
    clustering : bool
        Enable Community clustering and cluster-aware filtering when reference
        statistics are available.
    reference_lengths_tsv : str, optional
        TSV mapping of reference name -> length used to override header lengths.
    reference_stats_tsv : str, optional
        Path to write per-reference statistics TSV.

    Returns
    -------
    dict
        Summary information including alignment counts, EM statistics and
        processing metadata.
    """
    # Anti-dominance behavior is controlled via em_power_rho (φ-space EM)

    clock_start = bf_logging.start_timer()

    try:
    # Call score_alignments with dominance regularization adjustments
        result = score_alignments(
            bam_file=bam_file,
            output_bam=output_bam,
            num_threads=num_threads,
            verbose=verbose,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
            hierarchical_pmd=hierarchical_pmd,
            reference_lengths_tsv=reference_lengths_tsv,
            reference_stats_tsv=reference_stats_tsv,
            min_read_count=min_read_count,
            min_read_length=min_read_length,
            max_read_length=max_read_length,
            min_read_ani=min_read_ani,
            use_em=True,
            max_em_iterations=max_em_iterations,
            em_tolerance=em_tolerance,
            min_probability=min_probability,
            prob_fraction=prob_fraction,
            prior_weight=prior_weight,
            lambda_scale=lambda_scale,
            use_squarem_acceleration=use_squarem_acceleration,
            enable_globalization=enable_globalization,
            squarem_start_iter=squarem_start_iter,
            backtrack_factor=backtrack_factor,
            max_backtrack_steps=max_backtrack_steps,
            steplength_scheme=steplength_scheme,
            
            enable_dominance_regularization=enable_dominance_regularization,
            auto_tune_penalties=True,
            dominance_strength=dominance_strength,
            use_adaptive_dominance=use_adaptive_dominance,
            entropy_scaling_factor=entropy_scaling_factor,
            min_penalty_strength=min_penalty_strength,
            max_penalty_strength=max_penalty_strength,
            entropy_confidence_threshold=entropy_confidence_threshold,
            enable_emergency_regularization=enable_emergency_regularization,
            emergency_entropy_threshold=emergency_entropy_threshold,
            emergency_max_weight_threshold=emergency_max_weight_threshold,
            emergency_likelihood_drop=emergency_likelihood_drop,
            emergency_min_dominant_refs=emergency_min_dominant_refs,
            init_prior_strength=init_prior_strength,
            information_threshold=information_threshold,
            em_beta=em_beta,
            em_length_correction=em_length_correction,
            em_unknown_component=em_unknown_component,
            em_unknown_prior=em_unknown_prior,
            em_unknown_score=em_unknown_score,
            em_length_init=em_length_init,
            # === UNIFIED φ-SPACE EM PARAMETERS ===
            em_power_rho=em_power_rho,
            em_unknown_adaptive=em_unknown_adaptive,
            em_unknown_margin=em_unknown_margin,
            em_length_output_exp=em_length_output_exp,
            em_length_prior_exp=em_length_prior_exp,
            graph_min_edge_weight=graph_min_edge_weight,
            clustering=clustering,
            community_resolution=community_resolution,
            community_max_iterations=community_max_iterations,
            outlier_method=outlier_method,
            graph_export=graph_export,
            taxonomy_db=taxonomy_db,
            taxonomy_accession_map=taxonomy_accession_map,
            taxonomy_min_rank=taxonomy_min_rank,
            taxonomy_cross_domain_threshold=taxonomy_cross_domain_threshold,
            taxonomy_kingdom_threshold=taxonomy_kingdom_threshold,
            taxonomy_genus_threshold=taxonomy_genus_threshold,
            taxonomy_filter_enabled=taxonomy_filter_enabled,
            taxonomy_strict_filter=taxonomy_strict_filter,
            taxonomy_strict_min_connections=taxonomy_strict_min_connections,
            taxonomy_weighted_outlier=taxonomy_weighted_outlier,
            taxonomy_anomaly_weight=taxonomy_anomaly_weight,
            taxonomy_second_chance=taxonomy_second_chance,
            taxonomy_second_chance_cc=taxonomy_second_chance_cc,
            remove_cross_domain_alignments=remove_cross_domain_alignments,
            remove_cross_domain_references=remove_cross_domain_references,
            detect_misannotations=detect_misannotations,
            # Network QC & Taxonomic ambiguity detection
            network_qc_filter=network_qc_filter,
            entropy_biased_threshold=entropy_biased_threshold,
            entropy_mixed_threshold=entropy_mixed_threshold,
            entropy_highly_mixed_threshold=entropy_highly_mixed_threshold,
            tax_ambiguity_removal_level=tax_ambiguity_removal_level,
            # Coverage-Weighted Reference Priors
            cwrp_lambda=cwrp_lambda,
            iterative_auth=iterative_auth,
            auth_update_interval=auth_update_interval,
            damage_weight=damage_weight,
            low_cov_floor=low_cov_floor,
            low_cov_shrink_tau=low_cov_shrink_tau,
            # Posterior-Weighted Coverage Authenticity (Path B)
            auth_post_enabled=auth_post_enabled,
            auth_update_interval_post=auth_update_interval_post,
            auth_scale_post=auth_scale_post,
            auth_lambda_ramp_iters=auth_lambda_ramp_iters,
            # Sample-level P(ancient) gate
            sample_pi_override=sample_pi_override,
            # GMRF smoothing
            enable_gmrf=enable_gmrf,
            gmrf_tau=gmrf_tau,
        )

        processing_time = bf_monotonic_seconds() - clock_start

        # Enhanced result formatting
        summary = result.get('summary', {})
        formatted_result = {
            'success': True,
            'n_input_alignments': int(summary.get('total_alignments') or 0),
            'n_filtered_alignments': int(summary.get('filtered_alignments') or 0),
            'n_unique_reads': int(summary.get('unique_reads') or 0),
            'n_final_unique_reads': int(summary.get('final_unique_reads') or summary.get('unique_reads') or 0),
            'n_unique_references': int(summary.get('unique_references') or 0),
            'em_iterations': int(summary.get('em_iterations') or 0),
            'em_converged': bool(summary.get('em_converged') or False),
            'final_likelihood': float(summary.get('final_likelihood') or 0.0),
            'total_time': float(processing_time),
            'memory_pool_available': bool(result.get('memory_pool_available', False)),

            # PMD results
            'pmd_enabled': bool(result.get('pmd_enabled', False)),
            'library_type': result.get('library_type', 'ds'),
            
            # Existing features
            'zp_values_precomputed': bool(result.get('zp_values_precomputed', False)),
            'bam_writing_speedup': "100-1000x faster" if result.get('zp_values_precomputed', False) else "standard",
            'reference_lengths_used': bool(result.get('reference_lengths_used', False)),
            'squarem_used': bool(use_squarem_acceleration),
            'globalization_enabled': bool(enable_globalization),
            'steplength_scheme': f"S{steplength_scheme}",
            'paper_aligned': True,
            'implementation': 'Varadhan & Roland (2008) SQUAREM + PMD + Ultra-fast ZP'
        }

        if verbose:
            _info(f"PMD enabled: {formatted_result['pmd_enabled']}")
            _info(f"Implementation: {formatted_result['implementation']}")

        return formatted_result

    except Exception as e:
        processing_time = bf_monotonic_seconds() - clock_start

        if verbose:
            _error(f"Enhanced processing failed: {e}")

        return {
            'success': False,
            'error_message': str(e),
            'n_input_alignments': 0,
            'n_filtered_alignments': 0,
            'n_unique_reads': 0,
            'n_final_unique_reads': 0,
            'n_unique_references': 0,
            'em_iterations': 0,
            'em_converged': False,
            'final_likelihood': 0.0,
            'total_time': processing_time,
            'memory_pool_available': False,
            'pmd_enabled': calculate_pmd,
            'library_type': library_type,
            'zp_values_precomputed': False,
            'bam_writing_speedup': "failed",
            'reference_lengths_used': False,
            'squarem_used': use_squarem_acceleration,
            'globalization_enabled': enable_globalization,
            'steplength_scheme': f"S{steplength_scheme}",
            'paper_aligned': True,
            'implementation': 'Varadhan & Roland (2008) SQUAREM + PMD - FAILED'
        }
    finally:
        # Cleanup any remaining global memory
        with nogil:
            cleanup_presorted_memory()
