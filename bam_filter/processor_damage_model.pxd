# cython: language_level=3
# -*- coding: utf-8 -*-
"""Declarations for Bayesian damage model (ancient/modern classification)."""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t

from bam_filter.processor cimport MemoryPool, Alignment, AlignmentCore, DamageCounts
from bam_filter.processor_pmd cimport (
    RefDamageStats, DamageModelHyperparams, PMDCurve
)

# Allocation and cleanup
cdef RefDamageStats* allocate_ref_damage_stats(uint32_t n_refs) noexcept nogil
cdef void free_ref_damage_stats(RefDamageStats* stats) noexcept nogil

# Hyperparameter creation
cdef DamageModelHyperparams* create_hyperparams(
    PMDCurve* curve,
    double global_baseline,
    double baseline_strength,
    double amplitude_mu,
    double amplitude_sigma,
    double rho,
    double concentration
) noexcept nogil

# Core computations
cdef void estimate_baseline(
    RefDamageStats* stats,
    double prior_alpha,
    double prior_beta
) noexcept nogil

cdef void compute_bayes_factor(
    RefDamageStats* stats,
    DamageModelHyperparams* hyper,
    int use_both_ends
) noexcept nogil

cdef void fit_reference_damage_model(
    RefDamageStats* stats,
    DamageModelHyperparams* hyper,
    bint use_both_ends
) noexcept nogil

# Accumulation
cdef void accumulate_alignment_damage(
    RefDamageStats* ref_stats,
    AlignmentCore* core,
    DamageCounts* dmg,
    double phi_weight,
    bint is_single_stranded
) noexcept nogil
