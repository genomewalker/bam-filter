# cython: language_level=3
"""EM algorithm type declarations."""

from libc.stdint cimport uint32_t, int32_t, uint64_t, int64_t, uint8_t
from libc.stddef cimport size_t

from bam_filter.processor cimport MemoryPool


# =============================================================================
# EMState: All EM parameters in a single coherent structure
# =============================================================================

cdef struct EMState:
    # --- Reference weights in phi-space (what EM optimizes) ---
    double* phi_weights          # phi_j for each reference [n_refs]
    double* phi_counts           # E[c_j] accumulated counts [n_refs]

    # --- Hierarchical: per-reference ancient/modern classification ---
    bint hierarchical_enabled
    double* gamma_values         # gamma_j = P(ancient | ref j) [n_refs]
    double* S_anc                # E[ancient reads] per ref [n_refs]
    double* S_mod                # E[modern reads] per ref [n_refs]

    # --- Unknown component ---
    bint unknown_enabled
    double phi_unknown           # phi_u weight for unknown
    double S_unknown             # E[unknown reads] accumulated

    # --- PMD-based per-read priors (FIXED throughout EM) ---
    # omega_ig = P(G=g | PMD_i) computed once from PMD scores
    # Stored per-read: omega_ancient[i], omega_modern[i] = 1 - omega_ancient[i]
    double* omega_ancient        # P(ancient | PMD_i) per read [n_reads]

    # --- Dimensions ---
    uint32_t n_refs              # Number of references
    uint32_t n_reads             # Number of unique reads

    # --- Configuration (copied for nogil access) ---
    double dirichlet_prior       # alpha_0 for phi regularization
    double gamma_prior           # alpha_gamma for gamma regularization
    double power_rho             # rho for output transform (post-processing only)
    double unknown_margin        # Delta for adaptive unknown score
    double em_beta               # Temperature parameter


# =============================================================================
# EMConfig: Algorithm configuration
# =============================================================================

cdef struct EMConfig:
    # Iteration control
    int32_t max_iterations
    double convergence_tolerance

    # Priors
    double dirichlet_prior       # alpha_0 for phi (Dirichlet prior strength)
    double gamma_prior           # alpha_gamma for gamma (Beta prior strength)

    # Power transform (POST-PROCESSING ONLY)
    double power_rho             # 1.0 = standard, <1 flattens output

    # Unknown component
    bint unknown_enabled
    double unknown_margin        # Delta below max score for unknown likelihood

    # Hierarchical (ancient/modern)
    bint hierarchical_enabled
    float D_avg_5p               # Average 5' damage rate
    float D_avg_3p               # Average 3' damage rate
    float epsilon_error          # Sequencing error rate

    # SQUAREM acceleration
    bint squarem_enabled
    int32_t squarem_start_iter   # Start SQUAREM after this many iterations
    bint enable_globalization    # Enable backtracking
    double backtrack_factor      # Factor for backtracking (0.5)
    int32_t max_backtrack_steps  # Max backtrack attempts

    # Threading
    int32_t thread_count

    # Convergence (MAD-based)
    int32_t history_length       # Window for MAD statistics (default 5)


# =============================================================================
# SQUAREM state for acceleration
# =============================================================================

cdef struct SQUAREMState:
    double* theta_0              # Starting point
    double* theta_1              # After first EM step
    double* theta_2              # After second EM step
    double* r_vector             # r = theta_1 - theta_0
    double* v_vector             # v = (theta_2 - theta_1) - r
    double* theta_extrapolated   # Proposed accelerated point
    uint32_t dimension           # Size of parameter vector
    bint allocated


# =============================================================================
# Convergence tracking
# =============================================================================

cdef struct ConvergenceState:
    double* ll_history           # Log-likelihood changes [history_length]
    double* param_history        # Parameter changes [history_length]
    int32_t history_length
    int32_t current_index
    int32_t filled_count
    double current_ll
    double prev_ll


# =============================================================================
# Function declarations
# =============================================================================

# --- State management ---
cdef EMState* create_em_state(uint32_t n_refs, uint32_t n_reads,
                               EMConfig* config) noexcept nogil
cdef void free_em_state(EMState* state) noexcept nogil
cdef void copy_em_state(EMState* dest, EMState* src) noexcept nogil
cdef void reset_accumulators(EMState* state) noexcept nogil

# --- Core EM operations ---
cdef void e_step(EMState* state, void* pool, EMConfig* config) noexcept nogil
cdef void m_step(EMState* state, EMConfig* config) noexcept nogil
cdef double compute_log_likelihood(EMState* state, void* pool,
                                   EMConfig* config) noexcept nogil

# --- SQUAREM acceleration ---
cdef SQUAREMState* create_squarem_state(uint32_t dimension) noexcept nogil
cdef void free_squarem_state(SQUAREMState* sq) noexcept nogil
cdef bint squarem_step(EMState* state, void* pool, EMConfig* config,
                       SQUAREMState* sq, double* ll_out) noexcept nogil

# --- Output transform (POST-PROCESSING) ---
cdef void transform_to_output(EMState* state, double* output_pi,
                              double power_rho) noexcept nogil

# --- Convergence checking ---
cdef ConvergenceState* create_convergence_state(int32_t history_length) noexcept nogil
cdef void free_convergence_state(ConvergenceState* conv) noexcept nogil
cdef void update_convergence_history(ConvergenceState* conv,
                                     double ll_change, double param_change) noexcept nogil
cdef bint check_convergence_mad(ConvergenceState* conv,
                                double base_tolerance) noexcept nogil

# --- Main entry point ---
cdef int run_em(void* pool, EMConfig* config, double* output_weights) noexcept nogil

# --- Utility ---
cdef double compute_param_change_norm(EMState* current, EMState* prev) noexcept nogil
cdef void normalize_phi(EMState* state) noexcept nogil

# --- Helpers ---
cdef double* get_reference_weights(MemoryPool* pool) noexcept nogil
cdef void PREFETCH_READ(void* ptr) noexcept nogil
cdef void PREFETCH_WRITE(void* ptr) noexcept nogil
