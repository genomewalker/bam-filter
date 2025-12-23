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

"""Fast numerical operations for EM algorithm.

Provides numerically-stable mathematical operations including log-sum-exp
and safe weight normalization to prevent numerical instability.
"""

from libc.math cimport exp, log, fabs, fmax, log1p

cdef double LOG2 = 0.6931471805599453  # log(2.0)


cdef double stable_log_sum_exp(double log_a, double log_b) noexcept nogil:
    """Numerically stable branchless log-sum-exp for two values.

    Computes log(exp(log_a) + exp(log_b)) without overflow/underflow issues
    by factoring out the larger exponent before computing the sum.

    Uses branchless fmax/fabs and fast approximation for large differences.
    When |log_a - log_b| > 20, the smaller term contributes < 2e-9 and is skipped.

    Parameters
    ----------
    log_a : double
        First log-space value
    log_b : double
        Second log-space value

    Returns
    -------
    double
        log(exp(log_a) + exp(log_b))
    """
    cdef double m = fmax(log_a, log_b)
    cdef double diff = fabs(log_a - log_b)

    # Fast path: if diff > 20, exp(-diff) < 2e-9, smaller term is negligible
    if diff > 20.0:
        return m

    # Use log1p for better numerical accuracy with small values
    return m + log1p(exp(-diff))


cdef void safe_normalize_weights(double* weights, int dimension) noexcept nogil:
    """Safely normalize weights array to sum to 1.0.

    Enforces bounds [1e-15, 0.999] on individual weights and normalizes to sum 1.0.
    Falls back to uniform distribution if no valid weights exist.

    Parameters
    ----------
    weights : double*
        Weight array to normalize (modified in-place)
    dimension : int
        Number of weights
    """
    cdef int i
    cdef double total = 0.0
    cdef double uniform_weight
    cdef bint has_valid_weights = False
    cdef double min_weight = 1e-15
    cdef double max_weight = 0.999
    cdef double inv_total

    for i in range(dimension):
        if weights[i] < min_weight:
            weights[i] = min_weight
        elif weights[i] > max_weight:
            weights[i] = max_weight

        total += weights[i]
        if weights[i] > 1e-12:
            has_valid_weights = True

    if total > min_weight and has_valid_weights:
        inv_total = 1.0 / total
        for i in range(dimension):
            weights[i] *= inv_total
            if weights[i] < min_weight:
                weights[i] = min_weight
            elif weights[i] > max_weight:
                weights[i] = max_weight
    else:
        uniform_weight = 1.0 / dimension
        for i in range(dimension):
            weights[i] = uniform_weight

