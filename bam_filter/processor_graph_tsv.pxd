# cython: language_level=3
from libc.stdint cimport uint32_t, uint64_t, int32_t

from bam_filter.processor cimport MemoryPool
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.processor_graph cimport ReferencePattern, ReferenceStats

cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t

# C-level TSV writer forwarder (implemented in processor_graph_tsv.pyx)
cdef int write_graph_tsv_c(MemoryPool* pool, sam_hdr_t* bam_header,
                          ReferenceMapping* mapping, ReferencePattern* pattern_data,
                          ReferenceStats* ref_stats,
                          uint32_t* total_reads, uint32_t* multimap_reads, uint64_t* alignments_per_ref,
                          uint32_t* exact_connection_counts, double* co_mapping_averages,
                          uint64_t* max_co_mappings, uint32_t* co_mapping_counts,
                          double* neighbor_multimap_avg, double* neighbor_connections_avg,
                          uint32_t* neighbor_counts, uint32_t array_size,
                          double dataset_median_connections, int32_t min_read_count,
                          bint include_clustering, int outlier_method, const char* tsv_path) noexcept nogil
