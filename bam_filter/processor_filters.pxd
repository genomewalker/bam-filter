# cython: language_level=3
# Declarations exported by processor_filters.pyx

from libc.stdint cimport int32_t, int64_t, uint32_t
from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats, GraphMetrics, ReadIndex
from bam_filter.processor_types cimport EMAlgorithmConfig
from bam_filter.processor_graph_ops cimport WeightedGraph
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor_taxonomy_filters cimport TaxonomyFilterConfig, TaxonomyFilterStats
from bam_filter.taxonomy_db cimport TaxonomyDB

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
                                       bint verbose,
                                       double community_resolution,
                                       int community_max_iterations,
                                       uint32_t graph_min_edge_weight,
                                       int32_t thread_count,
                                       int outlier_method,
                                       WeightedGraph* existing_graph,
                                       sam_hdr_t* bam_header,
                                       ReferenceMapping* mapping,
                                       const char* tsv_export_path,
                                       const char* graph_export_path,
                                       TaxonomyFilterConfig* taxonomy_filter_config,
                                       TaxonomyFilterStats* taxonomy_stats_out,
                                       TaxonomyDB* taxonomy_db,
                                       float betweenness_threshold,
                                       float cc_threshold,
                                       uint32_t hub_degree_threshold,
                                       bint strict_mode,
                                       bint remove_cross_domain_edges,
                                       bint flag_misannotations) except -1 nogil
