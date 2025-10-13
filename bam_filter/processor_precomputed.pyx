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

"""Precomputed weight transformations for EM algorithm.

Caches frequently-used transformations (log, inverse, sqrt) of reference weights
to avoid redundant calculations during EM iterations.
"""

from libc.stdlib cimport malloc, free
from libc.math cimport log, sqrt as libc_sqrt, fmax
from libc.stdint cimport uint32_t

from bam_filter.processor_types cimport PrecomputedWeights


cdef PrecomputedWeights* create_precomputed_weights(uint32_t n_weights) noexcept nogil:
    """Create precomputed weights structure.

    Allocates arrays for storing log, inverse, and square root transformations
    of reference weights. Marked as dirty initially to trigger computation.

    Parameters
    ----------
    n_weights : uint32_t
        Number of reference weights

    Returns
    -------
    PrecomputedWeights*
        Precomputed weights structure, or NULL on allocation failure
    """
    cdef PrecomputedWeights* pw = <PrecomputedWeights*>malloc(sizeof(PrecomputedWeights))
    if not pw:
        return NULL

    pw.log_weights = <double*>malloc(n_weights * sizeof(double))
    pw.inv_weights = <double*>malloc(n_weights * sizeof(double))
    pw.sqrt_weights = <double*>malloc(n_weights * sizeof(double))

    if not pw.log_weights or not pw.inv_weights or not pw.sqrt_weights:
        if pw.log_weights: free(pw.log_weights)
        if pw.inv_weights: free(pw.inv_weights)
        if pw.sqrt_weights: free(pw.sqrt_weights)
        free(pw)
        return NULL

    pw.n_weights = n_weights
    pw.weights_dirty = True
    return pw


cdef void update_precomputed_weights(PrecomputedWeights* pw, double* weights) noexcept nogil:
    """Update precomputed transformations from reference weights.

    Computes log, inverse, and square root of weights with numerical safety
    clamping. Skips computation if weights haven't changed (dirty flag).

    Parameters
    ----------
    pw : PrecomputedWeights*
        Precomputed weights structure
    weights : double*
        Array of reference weights
    """
    cdef uint32_t i
    cdef double w, MIN_WEIGHT = 1e-12

    if not pw.weights_dirty:
        return

    for i in range(pw.n_weights):
        w = fmax(weights[i], MIN_WEIGHT)
        pw.log_weights[i] = log(w)
        pw.inv_weights[i] = 1.0 / w
        pw.sqrt_weights[i] = libc_sqrt(w)

    pw.weights_dirty = False


cdef void free_precomputed_weights(PrecomputedWeights* pw) noexcept nogil:
    """Free precomputed weights structure.

    Parameters
    ----------
    pw : PrecomputedWeights*
        Precomputed weights to free (safe to pass NULL)
    """
    if not pw:
        return
    if pw.log_weights: free(pw.log_weights)
    if pw.inv_weights: free(pw.inv_weights)
    if pw.sqrt_weights: free(pw.sqrt_weights)
    free(pw)
