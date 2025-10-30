# cython: language_level=3
# Declarations exported by processor_filters.pyx

from libc.stdint cimport int32_t, int64_t, uint32_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats, GraphMetrics, ReadIndex
from bam_filter.processor_types cimport EMAlgorithmConfig
from bam_filter.processor_graph_ops cimport WeightedGraph
from bam_filter.processor_mapping cimport ReferenceMapping

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t

# Probability filtering (alignment-level)
cdef int apply_probability_filtering(MemoryPool* pool, EMAlgorithmConfig* config) except -1 nogil

cdef int apply_cluster_aware_filtering(MemoryPool* pool,
                                       ReferencePattern* pattern_data,
                                       ReferenceStats* ref_stats,
                                       ReadIndex* read_index,
                                       uint32_t array_size,
                                       int32_t min_read_count,
                                       float score_threshold,
                                       bint verbose,
                                       int use_leiden,
                                       double leiden_resolution,
                                       bint leiden_parallel,
                                       int leiden_max_iterations,
                                       uint32_t graph_min_edge_weight,
                                       int32_t thread_count,
                                       int32_t iforest_n_trees,
                                       uint32_t iforest_subsample_size,
                                       double iforest_contamination,
                                       uint32_t iforest_random_seed,
                                       uint32_t lof_k,
                                       double lof_contamination,
                                       double zscore_threshold,
                                       WeightedGraph* existing_graph,
                                       sam_hdr_t* bam_header,
                                       ReferenceMapping* mapping,
                                       const char* tsv_export_path,
                                       const char* graph_export_path,
                                       int outlier_method) except -1 nogil
