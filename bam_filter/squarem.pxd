# cython: language_level=3
# -*- coding: utf-8 -*-
"""Generic SQUAREM acceleration for EM algorithms.

SQUAREM (Squared Iterative Methods) accelerates EM convergence by
extrapolating along the direction of iteration. This module provides
a reusable implementation that works with any parameter vector.

Reference: Varadhan & Roland (2008) "Simple and Globally Convergent
Methods for Accelerating the Convergence of Any EM Algorithm"
"""

from libc.stdint cimport uint32_t, int32_t


# =============================================================================
# Function pointer types for caller-provided operations
# =============================================================================

# EM fixed-point map: theta_out = F(theta_in)
ctypedef void (*squarem_update_fn)(
    double* theta_in,
    double* theta_out,
    uint32_t n,
    void* ctx
) noexcept nogil

# Objective function (log-likelihood or log-posterior) to maximize
ctypedef double (*squarem_objective_fn)(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil

# Project theta onto constraint set (in-place)
ctypedef void (*squarem_project_fn)(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil


# =============================================================================
# SQUAREM state and configuration
# =============================================================================

cdef struct SquaremState:
    uint32_t n              # parameter vector dimension
    double* theta0          # starting point
    double* theta1          # after 1 EM step
    double* theta2          # after 2 EM steps
    double* r               # r = theta1 - theta0
    double* v               # v = theta2 - 2*theta1 + theta0
    double* theta_extrap    # extrapolated point
    double* block           # contiguous backing storage
    size_t block_len        # size of block in bytes
    bint owns_block         # True if we allocated block


cdef struct SquaremConfig:
    bint enable                 # master switch
    bint enable_globalization   # backtracking line search
    int32_t max_backtracks      # max backtracking steps (default 4)
    double backtrack_factor     # step reduction factor (default 0.5)
    int32_t steplength_scheme   # 1=S1, 2=S2, 3=S3 (default 3)
    double alpha_min            # minimum alpha (default -300)
    double alpha_max            # maximum alpha (default -1)
    double rel_tol              # relative tolerance for acceptance
    bint project_after_em       # apply projection after each EM step


# =============================================================================
# Projection context structs
# =============================================================================

# Context for box constraint projection
cdef struct SquaremBoxCtx:
    double* lower           # per-parameter lower bounds (or NULL for scalar)
    double* upper           # per-parameter upper bounds (or NULL for scalar)
    double lower_scalar     # scalar lower bound (used if lower is NULL)
    double upper_scalar     # scalar upper bound (used if upper is NULL)

# Context for simplex constraint projection
cdef struct SquaremSimplexCtx:
    double min_value        # minimum allowed value (e.g., 1e-15)
    double target_sum       # target sum (e.g., 1.0)


# =============================================================================
# Core functions
# =============================================================================

# Return bytes needed for SQUAREM state with n parameters
cdef size_t squarem_bytes(uint32_t n) noexcept nogil

# Initialize config with default values
cdef void squarem_config_default(SquaremConfig* cfg) noexcept nogil

# Initialize SQUAREM state. Returns True on success.
# If external_block is provided and large enough, uses it.
# Otherwise allocates internally.
cdef bint squarem_state_init(
    SquaremState* sq,
    uint32_t n,
    double* external_block,
    size_t external_len
) noexcept nogil

# Free SQUAREM state resources
cdef void squarem_state_free(SquaremState* sq) noexcept nogil

# Perform one SQUAREM step (equivalent to 2+ EM steps).
# Returns True if SQUAREM extrapolation was used, False if fell back to plain EM.
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
) noexcept nogil


# =============================================================================
# Built-in projection functions
# =============================================================================

# Project theta onto probability simplex
# ctx should be SquaremSimplexCtx* or NULL for defaults
cdef void squarem_project_simplex(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil

# Project theta onto box constraints
# ctx should be SquaremBoxCtx*
cdef void squarem_project_box(
    double* theta,
    uint32_t n,
    void* ctx
) noexcept nogil
