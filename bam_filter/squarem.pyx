# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: initializedcheck=False
# -*- coding: utf-8 -*-
"""Generic SQUAREM acceleration for EM algorithms.

SQUAREM (Squared Iterative Methods) accelerates EM convergence by
extrapolating along the direction of iteration.

The algorithm:
1. Run two EM steps: theta0 -> theta1 -> theta2
2. Compute r = theta1 - theta0, v = theta2 - 2*theta1 + theta0
3. Compute step length alpha using scheme S1/S2/S3
4. Extrapolate: theta_sq = theta0 - 2*alpha*r + alpha^2*v
5. Stabilize with one EM step
6. Accept if objective improves, else backtrack or fallback

Reference: Varadhan & Roland (2008)
"""

from libc.math cimport sqrt as libc_sqrt, fabs, fmax, fmin
from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, memset
from libc.stdint cimport uint32_t, int32_t

# Logging
cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_notime(const char* tag, const char* fmt, ...) nogil

cdef const char* LOG_TAG = b"SQUAREM"


# =============================================================================
# Memory management
# =============================================================================

cdef size_t squarem_bytes(uint32_t n) noexcept nogil:
    """Return bytes needed for SQUAREM state with n parameters.

    We need 6 arrays: theta0, theta1, theta2, r, v, theta_extrap
    """
    return 6 * n * sizeof(double)


cdef void squarem_config_default(SquaremConfig* cfg) noexcept nogil:
    """Initialize config with default values."""
    cfg.enable = True
    cfg.enable_globalization = True
    cfg.max_backtracks = 4
    cfg.backtrack_factor = 0.5
    cfg.steplength_scheme = 3  # S3 recommended
    cfg.alpha_min = -300.0
    cfg.alpha_max = -1.0
    cfg.rel_tol = 1e-6
    cfg.project_after_em = True


cdef bint squarem_state_init(
    SquaremState* sq,
    uint32_t n,
    double* external_block,
    size_t external_len
) noexcept nogil:
    """Initialize SQUAREM state."""
    cdef size_t needed = squarem_bytes(n)

    sq.n = n
    sq.block = NULL
    sq.owns_block = False

    if external_block != NULL and external_len >= needed:
        sq.block = external_block
        sq.block_len = external_len
        sq.owns_block = False
    else:
        sq.block = <double*>malloc(needed)
        if sq.block == NULL:
            return False
        sq.block_len = needed
        sq.owns_block = True

    # Partition the block into 6 arrays
    sq.theta0 = sq.block
    sq.theta1 = sq.block + n
    sq.theta2 = sq.block + 2 * n
    sq.r = sq.block + 3 * n
    sq.v = sq.block + 4 * n
    sq.theta_extrap = sq.block + 5 * n

    return True


cdef void squarem_state_free(SquaremState* sq) noexcept nogil:
    """Free SQUAREM state resources."""
    if sq.owns_block and sq.block != NULL:
        free(sq.block)
    sq.block = NULL
    sq.theta0 = NULL
    sq.theta1 = NULL
    sq.theta2 = NULL
    sq.r = NULL
    sq.v = NULL
    sq.theta_extrap = NULL
    sq.owns_block = False


# =============================================================================
# Core SQUAREM step
# =============================================================================

cdef bint squarem_step(
    SquaremState* sq,
    double* theta_io,
    uint32_t n,
    squarem_update_fn update_fn,
    squarem_objective_fn obj_fn,
    squarem_project_fn project_fn,
    void* update_ctx,
    void* obj_ctx,
    void* project_ctx,
    SquaremConfig* cfg,
    double* ll_out
) noexcept nogil:
    """Perform one SQUAREM step."""

    cdef uint32_t i
    cdef double r_norm_sq, v_norm_sq, r_dot_v
    cdef double alpha, alpha_try
    cdef double ll_2, ll_sq, rel_tol
    cdef int backtrack
    cdef double bf
    cdef bint accepted = False

    if not cfg.enable:
        # Just run one EM step
        update_fn(theta_io, theta_io, n, update_ctx)
        if project_fn != NULL and cfg.project_after_em:
            project_fn(theta_io, n, project_ctx)
        if ll_out != NULL and obj_fn != NULL:
            ll_out[0] = obj_fn(theta_io, n, obj_ctx)
        return False

    # Step 1: Save theta0
    memcpy(sq.theta0, theta_io, n * sizeof(double))

    # Step 2: First EM step -> theta1
    update_fn(sq.theta0, sq.theta1, n, update_ctx)
    if project_fn != NULL and cfg.project_after_em:
        project_fn(sq.theta1, n, project_ctx)

    # Step 3: Second EM step -> theta2
    update_fn(sq.theta1, sq.theta2, n, update_ctx)
    if project_fn != NULL and cfg.project_after_em:
        project_fn(sq.theta2, n, project_ctx)

    # Compute objective at theta2 for comparison
    if obj_fn != NULL:
        ll_2 = obj_fn(sq.theta2, n, obj_ctx)
    else:
        ll_2 = 0.0

    # Step 4: Compute r = theta1 - theta0, v = theta2 - 2*theta1 + theta0
    r_norm_sq = 0.0
    v_norm_sq = 0.0
    r_dot_v = 0.0

    for i in range(n):
        sq.r[i] = sq.theta1[i] - sq.theta0[i]
        sq.v[i] = sq.theta2[i] - 2.0 * sq.theta1[i] + sq.theta0[i]
        r_norm_sq += sq.r[i] * sq.r[i]
        v_norm_sq += sq.v[i] * sq.v[i]
        r_dot_v += sq.r[i] * sq.v[i]

    # Check for degenerate cases
    if r_norm_sq < 1e-30 or v_norm_sq < 1e-30:
        # Already converged or degenerate, use theta2
        memcpy(theta_io, sq.theta2, n * sizeof(double))
        if ll_out != NULL:
            ll_out[0] = ll_2
        return False

    # Step 5: Compute alpha using selected scheme
    if cfg.steplength_scheme == 1:
        # S1: alpha = -sqrt(||r||^2 / ||v||^2)
        alpha = -libc_sqrt(r_norm_sq / v_norm_sq)
    elif cfg.steplength_scheme == 2:
        # S2: alpha = -||r||^2 / (r . v)
        if fabs(r_dot_v) < 1e-30:
            memcpy(theta_io, sq.theta2, n * sizeof(double))
            if ll_out != NULL:
                ll_out[0] = ll_2
            return False
        alpha = -r_norm_sq / r_dot_v
    else:
        # S3 (default): alpha = -(r . v) / ||v||^2
        alpha = -r_dot_v / v_norm_sq

    # Clamp alpha to valid range
    if alpha > cfg.alpha_max:
        # No extrapolation suggested, use theta2
        memcpy(theta_io, sq.theta2, n * sizeof(double))
        if ll_out != NULL:
            ll_out[0] = ll_2
        return False

    if alpha < cfg.alpha_min:
        alpha = cfg.alpha_min

    # Relative tolerance for acceptance
    rel_tol = cfg.rel_tol * fabs(ll_2) if ll_2 != 0.0 else 1e-8

    # Step 6: Globalization with backtracking
    bf = cfg.backtrack_factor

    for backtrack in range(cfg.max_backtracks if cfg.enable_globalization else 1):
        if backtrack == 0:
            alpha_try = alpha
        else:
            # Backtrack: move alpha toward -1
            alpha_try = -1.0 + (bf ** backtrack) * (alpha + 1.0)
            if alpha_try > -1.01:
                break

        # Extrapolate: theta_sq = theta0 - 2*alpha*r + alpha^2*v
        for i in range(n):
            sq.theta_extrap[i] = (sq.theta0[i]
                                  - 2.0 * alpha_try * sq.r[i]
                                  + alpha_try * alpha_try * sq.v[i])

        # Project onto constraints
        if project_fn != NULL:
            project_fn(sq.theta_extrap, n, project_ctx)

        # Stabilize with one EM step
        update_fn(sq.theta_extrap, sq.theta_extrap, n, update_ctx)
        if project_fn != NULL and cfg.project_after_em:
            project_fn(sq.theta_extrap, n, project_ctx)

        # Evaluate objective
        if obj_fn != NULL:
            ll_sq = obj_fn(sq.theta_extrap, n, obj_ctx)

            # Accept if improved
            if ll_sq >= ll_2 - rel_tol:
                memcpy(theta_io, sq.theta_extrap, n * sizeof(double))
                if ll_out != NULL:
                    ll_out[0] = ll_sq
                accepted = True
                break
        else:
            # No objective function, accept unconditionally
            memcpy(theta_io, sq.theta_extrap, n * sizeof(double))
            accepted = True
            break

    if not accepted:
        # Fallback to theta2
        memcpy(theta_io, sq.theta2, n * sizeof(double))
        if ll_out != NULL:
            ll_out[0] = ll_2

    return accepted


# =============================================================================
# Built-in projection functions
# =============================================================================

cdef void squarem_project_simplex(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil:
    """Project theta onto probability simplex.

    Ensures: theta[i] >= min_value, sum(theta) = target_sum
    """
    cdef SquaremSimplexCtx* sctx = <SquaremSimplexCtx*>ctx
    cdef double min_val = 1e-15
    cdef double target = 1.0
    cdef double total = 0.0
    cdef uint32_t i

    if sctx != NULL:
        min_val = sctx.min_value
        target = sctx.target_sum

    # Clamp to minimum and compute sum
    for i in range(n):
        if theta[i] < min_val:
            theta[i] = min_val
        total += theta[i]

    # Normalize to target sum
    if total > 0:
        for i in range(n):
            theta[i] = theta[i] * target / total


cdef void squarem_project_box(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil:
    """Project theta onto box constraints.

    Clamps each component to [lower, upper].
    """
    cdef SquaremBoxCtx* bctx = <SquaremBoxCtx*>ctx
    cdef double lo, hi
    cdef uint32_t i

    if bctx == NULL:
        return

    for i in range(n):
        if bctx.lower != NULL:
            lo = bctx.lower[i]
        else:
            lo = bctx.lower_scalar

        if bctx.upper != NULL:
            hi = bctx.upper[i]
        else:
            hi = bctx.upper_scalar

        theta[i] = fmax(lo, fmin(hi, theta[i]))
