# cython: language_level=3
"""
Probabilistic Taxonomic Profiler - Cython Header

Novel approach combining:
- DCMS: Damage-Calibrated Mixture-of-Sources Model (per-reference)
- TBN-RN: Taxonomic Bayesian Network with belief propagation
- TAD-PUT: Probabilistic abundance estimation
"""

from libc.stdint cimport int32_t, int64_t, uint32_t, uint64_t


# =============================================================================
# HYPERPARAMETERS
# =============================================================================

cdef struct DamageHyperparams:
    # Ancient damage exponential decay parameters
    double p_0          # End damage mean (default: 0.40)
    double p_bg         # Background error floor (default: 0.005)
    double tau_A        # Decay length in bases (default: 5.0)
    double kappa_A      # Ancient Beta concentration (default: 50.0)

    # Modern flat error parameters
    double p_err        # Modern error rate (default: 0.005)
    double kappa_M      # Modern Beta concentration (default: 200.0)

    # Prior
    double pi_A_ref     # Prior P(ancient) per reference (default: 0.10)

    # Precomputed values (call precompute_damage_params to initialize)
    double a_M                  # Modern Beta alpha
    double b_M                  # Modern Beta beta
    double log_beta_prior_M     # log B(a_M, b_M)
    double a_A[20]              # Ancient Beta alpha per position
    double b_A[20]              # Ancient Beta beta per position
    double log_beta_prior_A[20] # log B(a_A[j], b_A[j]) per position


cdef struct PresenceHyperparams:
    # Authenticity sigmoid parameters
    # P(present) = sigmoid(scale * authenticity_score + intercept)
    # where authenticity_score = norm_entropy - norm_gini
    double intercept    # Sigmoid intercept (default: 0.0)
    double scale        # Sigmoid scale (default: 5.0)


cdef struct EpochHyperparams:
    # Root prior
    double pi_0         # P(E_root = ancient) (default: 0.5)

    # Transition probabilities
    double rho_AA       # P(child ancient | parent ancient) (default: 0.90)
    double rho_MA       # P(child ancient | parent modern) (default: 0.01)


cdef struct NoisyORHyperparams:
    # Noisy-OR parameters
    double lambda_leak  # Leak presence probability (default: 0.01)
    double q_fail       # Child->parent failure probability (default: 0.10)
    double q_fail_ref   # Reference->leaf failure probability (default: 0.10)


cdef struct ProfilerHyperparams:
    DamageHyperparams damage
    PresenceHyperparams presence
    EpochHyperparams epoch
    NoisyORHyperparams noisy_or


# =============================================================================
# PER-REFERENCE DATA
# =============================================================================

cdef struct RefDamageData:
    # Per-position damage counts (20 positions from each end)
    double n_5p[20]     # C opportunities at 5' positions
    double k_5p[20]     # C->T mismatches at 5' positions
    double n_3p[20]     # G opportunities at 3' positions
    double k_3p[20]     # G->A mismatches at 3' positions


cdef struct RefCoverageData:
    double breadth          # Fraction of reference covered
    double mean_depth       # Mean coverage depth
    double wcb              # Weighted Contiguity Breadth (Herfindahl)
    double norm_entropy     # Normalized spatial entropy [0, 1]
    double norm_gini        # Normalized Gini coefficient [0, 1]
    double authenticity     # norm_entropy - norm_gini [-1, 1]
    int64_t n_reads         # Number of mapped reads
    double tad              # Truncated average depth


cdef struct RefPosteriors:
    double log_bf           # Log Bayes factor (ancient vs modern)
    double p_ancient        # P(ancient | reference)
    double p_present        # P(present | reference)
    double tad_ancient      # TAD * P(ancient) * P(present)
    double tad_modern       # TAD * (1-P(ancient)) * P(present)


cdef struct RefProfileData:
    int64_t ref_id          # Reference ID
    int32_t taxid           # Taxonomic ID for this reference
    RefDamageData damage
    RefCoverageData coverage
    RefPosteriors posteriors


# =============================================================================
# PER-TAXON DATA (for belief propagation)
# =============================================================================

cdef struct TaxonBeliefs:
    # Observation factors from references
    double psi_ancient      # Product of P(ancient|r) for refs in taxon
    double psi_modern       # Product of (1-P(ancient|r)) for refs in taxon
    int32_t n_refs          # Number of references in this taxon

    # Local belief (before messages)
    double b_ancient        # Local belief for ancient state
    double b_modern         # Local belief for modern state

    # Messages from children (product)
    double msg_up_ancient   # Product of upward messages for ancient
    double msg_up_modern    # Product of upward messages for modern

    # Message from parent
    double msg_down_ancient # Downward message for ancient
    double msg_down_modern  # Downward message for modern

    # Final marginals
    double p_ancient        # P(E_t = ancient)
    double p_modern         # P(E_t = modern)

    # Presence (Noisy-OR)
    double p_present        # P(present_t)

    # Abundance
    double tad_ancient      # Aggregated ancient TAD
    double tad_modern       # Aggregated modern TAD
    double tad_total        # Total TAD


# =============================================================================
# TREE STRUCTURE FOR NATIVE BELIEF PROPAGATION
# =============================================================================

cdef struct TaxonNode:
    int32_t taxid
    int32_t parent_idx     # Index into nodes array (-1 for root)
    int32_t first_child    # Index of first child (-1 if leaf)
    int32_t next_sibling   # Index of next sibling (-1 if last)
    int32_t first_ref      # Index of first ref in refs array (-1 if none)
    int32_t n_refs         # Number of refs for this taxon
    int32_t n_children     # Number of children


# =============================================================================
# FUNCTION DECLARATIONS
# =============================================================================

# Hyperparameter initialization
cdef void init_default_hyperparams(ProfilerHyperparams* params) noexcept nogil

# Parallel batch processing
cdef void compute_batch_posteriors_parallel(
    RefProfileData* refs,
    int32_t n_refs,
    ProfilerHyperparams* params,
    int num_threads
) noexcept nogil

# Beta-Binomial computations
cdef double log_beta_function(double a, double b) noexcept nogil
cdef double log_beta_binomial(int64_t k, int64_t n, double a, double b) noexcept nogil

# Per-reference computations
cdef double compute_damage_log_bf(
    RefDamageData* damage,
    DamageHyperparams* params
) noexcept nogil

cdef double compute_p_ancient(double log_bf, double pi_A_ref) noexcept nogil

cdef double compute_p_present(
    RefCoverageData* coverage,
    PresenceHyperparams* params
) noexcept nogil

cdef void compute_ref_posteriors(
    RefProfileData* ref,
    ProfilerHyperparams* params
) noexcept nogil

# Belief propagation
cdef void init_taxon_beliefs(TaxonBeliefs* beliefs) noexcept nogil

cdef void accumulate_ref_observation(
    TaxonBeliefs* beliefs,
    double p_ancient_ref
) noexcept nogil

cdef void compute_upward_message(
    TaxonBeliefs* child,
    double* msg_ancient,
    double* msg_modern,
    EpochHyperparams* params
) noexcept nogil

cdef void receive_upward_message(
    TaxonBeliefs* parent,
    double msg_ancient,
    double msg_modern
) noexcept nogil

cdef void compute_downward_message(
    TaxonBeliefs* parent,
    TaxonBeliefs* child,
    double* msg_ancient,
    double* msg_modern,
    EpochHyperparams* params,
    bint is_root
) noexcept nogil

cdef void receive_downward_message(
    TaxonBeliefs* child,
    double msg_ancient,
    double msg_modern
) noexcept nogil

cdef void compute_marginals(
    TaxonBeliefs* beliefs,
    EpochHyperparams* params,
    bint is_root
) noexcept nogil

# Noisy-OR presence
cdef double compute_noisy_or_presence_from_refs(
    double* p_present_refs,
    int32_t n_refs,
    NoisyORHyperparams* params
) noexcept nogil

cdef double compute_noisy_or_presence_from_children(
    double* p_present_children,
    int32_t n_children,
    NoisyORHyperparams* params
) noexcept nogil

# TAD aggregation
cdef void aggregate_tad_from_refs(
    TaxonBeliefs* beliefs,
    RefProfileData* refs,
    int32_t n_refs
) noexcept nogil

cdef void aggregate_tad_from_children(
    TaxonBeliefs* parent,
    TaxonBeliefs* children,
    int32_t n_children
) noexcept nogil


# =============================================================================
# GMRF (GAUSSIAN MARKOV RANDOM FIELD) FOR REFERENCE-LEVEL SMOOTHING
# =============================================================================

cdef struct GMRFHyperparams:
    double tau           # Smoothing precision (higher = more smoothing)
    double mu_eta        # Global prior mean for logit ancientness
    double epsilon       # Regularization for isolated nodes
    int32_t max_iter     # Maximum coordinate descent iterations
    double tol           # Convergence tolerance


cdef struct GMRFState:
    double* eta                   # Ancientness logits for each reference
    double* eta_data              # Data-driven initialization (from damage model)
    double* gamma                 # Smoothed P(ancient) = sigmoid(eta)
    double* degree                # Weighted degree for each reference (sum of edge weights)
    int32_t n_refs                # Number of references
    bint converged                # Whether GMRF optimization converged


# GMRF initialization and destruction
cdef GMRFState* create_gmrf_state(int32_t n_refs) noexcept nogil
cdef void destroy_gmrf_state(GMRFState* state) noexcept nogil

# Initialize eta from damage model log Bayes factors
cdef void init_eta_from_damage(
    GMRFState* state,
    double* log_bf,
    double prior_log_odds,
    int32_t n_refs
) noexcept nogil

# Build weighted degrees from graph structure
cdef void build_gmrf_degrees(
    GMRFState* state,
    uint32_t** neighbors,
    uint32_t** weights,
    uint32_t* degrees,
    int32_t n_refs
) noexcept nogil

# Run GMRF coordinate descent smoothing
cdef int run_gmrf_smoothing(
    GMRFState* state,
    uint32_t** neighbors,
    uint32_t** weights,
    uint32_t* degrees,
    uint32_t* n_reads,
    GMRFHyperparams* params,
    int num_threads
) noexcept nogil

# Convert smoothed eta to gamma (P(ancient))
cdef void compute_gamma_from_smoothed_eta(GMRFState* state) noexcept nogil


# =============================================================================
# CYTHON-LEVEL GMRF INTEGRATION WITH WEIGHTEDGRAPH
# =============================================================================

# Import WeightedGraph for the integration function
from bam_filter.processor_graph_ops cimport WeightedGraph

# Direct integration with WeightedGraph (no file I/O, no Python overhead)
cdef int apply_gmrf_smoothing_to_graph(
    WeightedGraph* graph,
    double* gamma_values,
    double* damage_log_bf,
    uint32_t* n_reads,
    uint32_t n_refs,
    double tau,
    int max_iter,
    double tol,
    int verbose,
) noexcept nogil
