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

"""Convergence detection helpers for EM algorithm.

Provides robust statistical functions (median, MAD) for detecting EM
convergence anomalies and stagnation patterns.
"""

from libc.math cimport fabs

cdef double _median5(double* a) except -1.0 nogil:
    """Compute median of 5-element array using selection network.

    Parameters
    ----------
    a : double*
        Array with exactly 5 elements

    Returns
    -------
    double
        Median value
    """
    cdef double x0=a[0], x1=a[1], x2=a[2], x3=a[3], x4=a[4]
    cdef double t
    if x1 < x0: t=x0; x0=x1; x1=t
    if x2 < x1: t=x1; x1=x2; x2=t
    if x1 < x0: t=x0; x0=x1; x1=t
    if x3 < x2: t=x2; x2=x3; x3=t
    if x2 < x1: t=x1; x1=x2; x2=t
    if x1 < x0: t=x0; x0=x1; x1=t
    if x4 < x3: t=x3; x3=x4; x4=t
    if x3 < x2: t=x2; x2=x3; x3=t
    return x2

cdef int _count_filled(int iteration, int H) except -1 nogil:
    """Count filled entries in history buffer.

    Parameters
    ----------
    iteration : int
        Current iteration number
    H : int
        History buffer capacity

    Returns
    -------
    int
        Number of filled entries (capped at H)
    """
    cdef int filled = iteration + 1
    if filled > H: filled = H
    return filled

cdef void _robust_sigma(double* buf, int H, int filled, double* med_out, double* sigma_out) noexcept nogil:
    """Compute robust median and sigma (MAD-based) from history buffer.

    Uses Median Absolute Deviation (MAD) scaled by 1.4826 for consistency
    with normal distribution standard deviation. Handles partially-filled
    buffers by replicating last value.

    Parameters
    ----------
    buf : double*
        History buffer (size H)
    H : int
        Buffer capacity (assumes 5)
    filled : int
        Number of filled entries
    med_out : double*
        Output for median value
    sigma_out : double*
        Output for robust sigma (1.4826 * MAD)
    """
    cdef double tmp[5]
    cdef double dev[5]
    cdef int i
    for i in range(H):
        tmp[i] = buf[i]
    if filled < H:
        for i in range(filled, H):
            tmp[i] = tmp[filled-1]
    cdef double med = _median5(tmp)
    for i in range(H):
        dev[i] = fabs(tmp[i] - med)
    cdef double mad = _median5(dev)
    med_out[0] = med
    sigma_out[0] = 1.4826 * mad
