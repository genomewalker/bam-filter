# cython: language_level=3

"""
Shared C-level type declarations used across Cython modules in the
``bam_filter`` package.

Place small, widely-used struct and enum declarations here to avoid
duplicate incompatible definitions across modules. Keep this file
stable to minimize ABI churn.
"""
from libc.stdint cimport uint32_t

# Minimal shared type declarations for processor helper modules

cdef struct PrecomputedWeights:
    double* log_weights
    double* inv_weights
    double* sqrt_weights
    bint weights_dirty
    uint32_t n_weights

# -----------------------------------------------------------------------------
# Additional shared types moved here so processor modules can cimport a single
# central header instead of depending on `bam_types.pxd`.
# -----------------------------------------------------------------------------
from libc.stdint cimport int32_t, int64_t, uint16_t, uint8_t, uint32_t, uint64_t
from libc.stddef cimport size_t

cdef enum ProcessingError:
    PROCESSING_SUCCESS = 0
    PROCESSING_ERROR_FILE_ACCESS = 1
    PROCESSING_ERROR_MEMORY_ALLOCATION = 2
    PROCESSING_ERROR_INVALID_DATA = 3
    PROCESSING_ERROR_ALGORITHM_FAILURE = 4

cdef struct CompactAlignment:
    uint32_t read_index
    uint32_t reference_index
    uint32_t position_info
    float    alignment_score

# Pull small helpers from central header to avoid duplication
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

cdef struct ProcessingBatch:
    int64_t batch_identifier
    int64_t reference_start_index
    int64_t reference_end_index
    int64_t expected_alignment_count
    int64_t actual_alignment_count
    CompactAlignment* batch_alignments
    int64_t batch_capacity
    int32_t error_status

cdef struct SQUAREMAccelerator:
    double* x_prev2
    double* x_prev1
    double* x_curr
    double* temp_diff1
    double* temp_diff2
    double* temp_extrap
    int dimension
    int iteration_count
    double alpha_min
    double alpha_max
    bint has_history
    bint owns_memory

cdef struct EMAlgorithmConfig:
    int32_t maximum_iterations
    double convergence_tolerance
    double minimum_probability_threshold
    double probability_fraction_filter
    double regularization_weight
    double score_scaling_factor
    bint enable_acceleration
    bint use_squarem_acceleration
    int32_t thread_count
    int32_t minimum_read_coverage

    # Squarem control
    int32_t squarem_start_iter

    # Globalization / backtracking (stabilized SQUAREM control)
    bint enable_globalization
    double backtrack_factor
    int32_t max_backtrack_steps

    # Steplength scheme selection
    int32_t steplength_scheme

    # Dominance regularization (used by improved M-step)
    bint enable_dominance_regularization
    double dominance_strength
    bint use_adaptive_dominance
    double entropy_scaling_factor
    double min_penalty_strength
    double max_penalty_strength
    double entropy_confidence_threshold

    bint enable_emergency_regularization    # Allow automatic activation
    bint dominance_regularization_active    # Runtime activation flag
    double emergency_entropy_threshold      # Entropy below which to activate
    double emergency_max_weight_threshold   # Single ref weight above which to activate
    uint32_t emergency_min_dominant_refs    # Minimum dominant refs for large datasets
    double emergency_likelihood_drop        # LL drop threshold for activation
    double init_prior_strength              # Dirichlet prior strength for initialization
    
    # Information-theoretic filtering
    float information_threshold             # Threshold for information content filtering (0.0 = disabled)


cdef extern from "unistd.h":
    int getpagesize() nogil

cdef extern from "htslib/sam.h" nogil:
    ctypedef struct BGZF

    ctypedef struct bam1_core_t:
        int64_t tid
        int64_t pos
        uint64_t bin
        uint8_t qual
        uint8_t l_qname
        uint16_t flag
        uint32_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int64_t mpos
        int32_t isize

    ctypedef struct bam1_t:
        bam1_core_t core
        uint64_t id
        uint8_t* data
        int l_data
        uint32_t m_data

    ctypedef struct sam_hdr_t:
        char* text
        int n_targets
        char** target_name

    ctypedef struct hts_idx_t
    ctypedef struct hts_itr_t

    ctypedef union htsFile_fp:
        BGZF* bgzf

    ctypedef struct htsFile:
        htsFile_fp fp
    
    # Expose samFile alias used across the codebase
    ctypedef htsFile samFile

    # File operations
    samFile* hts_open(const char* fn, const char* mode)
    int hts_close(samFile* fp)

    # Header operations
    sam_hdr_t* sam_hdr_read(samFile* fp)
    void sam_hdr_destroy(sam_hdr_t* h)
    const char* sam_hdr_str(sam_hdr_t* header)
    const char* sam_hdr_tid2name(sam_hdr_t* header, int tid)
    int64_t sam_hdr_tid2len(sam_hdr_t* header, int tid)
    sam_hdr_t* sam_hdr_parse(size_t l_text, const char* text)
    int32_t sam_hdr_nref(sam_hdr_t* header)
    int32_t bam_endpos(bam1_t* b) nogil
    sam_hdr_t* sam_hdr_init()
    int sam_hdr_add_line(sam_hdr_t* h, const char* tag, ...)
    int sam_hdr_add_lines(sam_hdr_t* h, const char* text, int keep)

    # Read/Write
    int sam_read1(samFile* fp, sam_hdr_t* h, bam1_t* b)
    int sam_write1(samFile* fp, const sam_hdr_t* h, const bam1_t* b)
    int sam_hdr_write(samFile* fp, const sam_hdr_t* h)

    # BAM object management
    bam1_t* bam_init1()
    void bam_destroy1(bam1_t* b)
    bam1_t* bam_dup1(bam1_t* src)

    # BAM data access
    char* bam_get_qname(bam1_t* b)
    uint8_t* bam_aux_get(bam1_t* b, const char* tag)
    int bam_aux2i(const uint8_t* s)
    int bam_aux_del(bam1_t* b, uint8_t* s)
    int bam_aux_append(bam1_t* b, const char* tag, char type, int len, const uint8_t* data)

    # Index operations
    hts_idx_t* sam_index_load(samFile* fp, const char* fn)
    void hts_idx_destroy(hts_idx_t* idx)
    int hts_idx_get_stat(hts_idx_t* idx, int tid, uint64_t* mapped, uint64_t* unmapped)

    # Iterators
    hts_itr_t* sam_itr_queryi(const hts_idx_t* idx, int tid, int beg, int end)
    int sam_itr_next(samFile* fp, hts_itr_t* iter, bam1_t* b)
    void hts_itr_destroy(hts_itr_t* iter)
    int32_t sam_hdr_name2tid(sam_hdr_t* header, const char* name)

    # Threads
    int hts_set_threads(samFile* fp, int n)
    uint8_t* bam_get_seq(bam1_t *b)
    uint8_t* bam_get_qual(bam1_t *b)
    uint32_t* bam_get_cigar(bam1_t *b)


cdef extern from "htslib/hts.h" nogil:
    char* seq_nt16_str
    int hts_set_opt(void* fp, int opt, ...) nogil
    cdef int HTS_OPT_CACHE_SIZE

# Add a small inline wrapper for the bam_seqi macro so Cython can call it
cdef extern from *:
    """
    #include "htslib/sam.h"
    static inline int bam_seqi_wrapper(const uint8_t *s, int i) {
        return bam_seqi(s, i);
    }
    """
    int bam_seqi_wrapper(const uint8_t* s, int i) nogil
