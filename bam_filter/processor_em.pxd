# cython: language_level=3
# -*- coding: utf-8 -*-

from libc.stdint cimport uint32_t, int32_t, int64_t, uint8_t
from libc.stddef cimport size_t

# Import core types from processor.pxd
from bam_filter.processor cimport MemoryPool, AlignmentScoringConfig
from bam_filter.processor_types cimport EMAlgorithmConfig
from bam_filter.processor_precomputed cimport PrecomputedWeights

# Helper functions for memory access
cdef double* get_reference_weights(MemoryPool* pool) noexcept nogil
cdef double* get_temp_buffer_A(MemoryPool* pool) noexcept nogil
cdef double* get_temp_buffer_B(MemoryPool* pool) noexcept nogil

# Prefetch optimization helpers
cdef void PREFETCH_READ(void* ptr) noexcept nogil
cdef void PREFETCH_WRITE(void* ptr) noexcept nogil

# EM algorithm functions
cdef double calculate_dataset_entropy(MemoryPool* pool) except -1.0 nogil
cdef double auto_tune_dominance_strength(MemoryPool* pool) except -1.0 nogil
cdef double learn_entropy_threshold_from_data(MemoryPool* pool, EMAlgorithmConfig* config) except -1.0 nogil
cdef double calculate_adaptive_dominance_strength(MemoryPool* pool, EMAlgorithmConfig* config) except -1.0 nogil
cdef double compute_dominance_penalty(double pi_j, double penalty_strength) except -1.0 nogil

cdef void execute_em_expectation_step_vectorized(MemoryPool* pool,
                                                int32_t thread_count,
                                                PrecomputedWeights* precomp,
                                                EMAlgorithmConfig* config=*,
                                                int32_t current_iteration=*) noexcept nogil

cdef void execute_em_maximization_step_vectorized(MemoryPool* pool,
                                                 double alpha_prior,
                                                 PrecomputedWeights* precomp,
                                                 EMAlgorithmConfig* config=*,
                                                 int current_iteration=*,
                                                 double current_ll=*,
                                                 double prev_ll=*) noexcept nogil

cdef double compute_log_likelihood_vectorized(MemoryPool* pool,
                                            PrecomputedWeights* precomp) except -1.0 nogil

cdef bint detect_pathological_convergence(MemoryPool* pool, 
                                         EMAlgorithmConfig* config,
                                         int current_iteration, 
                                         double current_ll, 
                                         double prev_ll) noexcept nogil

cdef int execute_em_algorithm(MemoryPool* pool, EMAlgorithmConfig* config) except -1 nogil

# NOTE: apply_probability_filtering moved to processor_filters.pyx for better organization

cdef int check_trend_based_convergence(double current_ll, double prev_ll,
                                      MemoryPool* pool, double* prev_weights,
                                      double base_tolerance, int iteration,
                                      EMAlgorithmConfig* config,
                                      PrecomputedWeights* precomp,
                                      double* ll_history, double* param_history,
                                      double alpha) noexcept nogil

cdef double compute_residual_norm_efficient(MemoryPool* pool, double* prev_weights,
                                           EMAlgorithmConfig* config,
                                           PrecomputedWeights* precomp) except -1.0 nogil

# Initialization functions  
cdef struct UniqueAlignment:
    uint32_t read_index
    uint32_t reference_index
    float alignment_score
    float confidence_score

cdef struct InitializationStats:
    double min_score
    double max_score
    double mean_score
    double score_range
    int64_t total_alignments
    int64_t unique_alignments
    double uniqueness_ratio

cdef void compute_initialization_stats(MemoryPool* pool, InitializationStats* stats) noexcept nogil
cdef int64_t identify_unique_alignments(MemoryPool* pool, InitializationStats* stats,
                                       UniqueAlignment* unique_alignments, int64_t capacity) noexcept nogil
cdef void compute_quality_weights(MemoryPool* pool, InitializationStats* stats,
                                UniqueAlignment* unique_alignments, int64_t unique_count,
                                double* quality_weights) noexcept nogil
cdef void initialize_em_weights(MemoryPool* pool, EMAlgorithmConfig* config) noexcept nogil
cdef void diagnose_initialization_quality(MemoryPool* pool) noexcept nogil
