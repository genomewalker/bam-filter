# cython: language_level=3
"""
Header file for graph export functionality.

Exports read-reference co-mapping graphs to standard network formats
(GraphML) for visualization in Cytoscape, igraph, Gephi, etc.
"""

from bam_filter.processor cimport MemoryPool
from bam_filter.processor_graph cimport ReferenceStats, ReferencePattern, WeightedGraph
from bam_filter.processor_mapping cimport ReferenceMapping
from bam_filter.taxonomy_db cimport TaxonomyDB

from bam_filter.processor cimport sam_hdr_t

# Main export function
cdef int export_graph_graphml(
    WeightedGraph* graph,
    MemoryPool* pool,
    sam_hdr_t* bam_header,
    ReferenceMapping* mapping,
    ReferencePattern* pattern_data,
    ReferenceStats* ref_stats,
    const char* output_path,
    bint verbose,
    bint export_only_used_stats,
    TaxonomyDB* taxonomy_db,
) nogil
