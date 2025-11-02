from libc.stdint cimport int32_t, uint32_t, int64_t, uint64_t

cdef extern from "igraph.h":
    # Core types
    ctypedef struct igraph_t:
        pass
    
    ctypedef struct igraph_vector_t:
        pass
    
    ctypedef struct igraph_vector_int_t:
        pass
    
    ctypedef struct igraph_vector_long_t:
        pass
    
    ctypedef struct igraph_vector_bool_t:
        pass

    # Vertex selector struct (declare here so it is a known type for signatures)
    ctypedef struct igraph_vs_t:
        pass
    
    # Modern igraph uses bool instead of enum
    ctypedef bint igraph_bool_t
    
    ctypedef double igraph_real_t
    ctypedef long int igraph_integer_t
    
    # Error type for newer igraph
    ctypedef int igraph_error_t
    
    # Vector functions - newer igraph API
    igraph_error_t igraph_vector_init(igraph_vector_t *v, igraph_integer_t size) nogil
    void igraph_vector_destroy(igraph_vector_t *v) nogil
    igraph_integer_t igraph_vector_size(const igraph_vector_t *v) nogil
    igraph_real_t* VECTOR(const igraph_vector_t v) nogil  # Macro - returns pointer to storage
    void igraph_vector_set(igraph_vector_t *v, igraph_integer_t pos, igraph_real_t value) nogil
    igraph_real_t* VECTOR_PTR(const igraph_vector_t v, igraph_integer_t pos) nogil
    
    igraph_error_t igraph_vector_int_init(igraph_vector_int_t *v, igraph_integer_t size) nogil
    void igraph_vector_int_destroy(igraph_vector_int_t *v) nogil
    igraph_integer_t igraph_vector_int_size(const igraph_vector_int_t *v) nogil
    void igraph_vector_int_set(igraph_vector_int_t *v, igraph_integer_t pos, igraph_integer_t value) nogil
    
    int igraph_vector_long_init(igraph_vector_long_t *v, long int size) nogil
    int igraph_vector_long_destroy(igraph_vector_long_t *v) nogil
    long int igraph_vector_long_size(const igraph_vector_long_t *v) nogil
    long int igraph_vector_long_e(const igraph_vector_long_t *v, long int pos) nogil
    void igraph_vector_long_set(igraph_vector_long_t *v, long int pos, long int value) nogil
    
    # Graph creation and destruction - newer API uses igraph_vector_int_t for edges
    igraph_error_t igraph_create(igraph_t *graph, const igraph_vector_int_t *edges, 
                                igraph_integer_t n, igraph_bool_t directed) nogil
    void igraph_destroy(igraph_t *graph) nogil
    igraph_integer_t igraph_vcount(const igraph_t *graph) nogil
    igraph_integer_t igraph_ecount(const igraph_t *graph) nogil
    
    # Graph construction
    int igraph_add_vertices(igraph_t *graph, igraph_integer_t nv, void *attr) nogil
    int igraph_add_edges(igraph_t *graph, const igraph_vector_t *edges, void *attr) nogil
    
    # Edge weights
    int igraph_es_all(igraph_es_t *es, igraph_edgeseq_type_t type) nogil
    
    # Community detection - Leiden algorithm (actual API signature)
    igraph_error_t igraph_community_leiden(const igraph_t *graph,
                                          const igraph_vector_t *edge_weights,
                                          const igraph_vector_t *node_weights,
                                          igraph_real_t resolution_parameter,
                                          igraph_real_t beta,
                                          igraph_bool_t start,
                                          igraph_integer_t n_iterations,
                                          igraph_vector_int_t *membership,
                                          igraph_integer_t *nb_clusters,
                                          igraph_real_t *quality) nogil

    # Community detection - Label Propagation Algorithm (LPA)
    # Fast O(m) algorithm for large graphs
    igraph_error_t igraph_community_label_propagation(const igraph_t *graph,
                                                      igraph_vector_int_t *membership,
                                                      igraph_neimode_t mode,
                                                      const igraph_vector_t *weights,
                                                      const igraph_vector_int_t *initial,
                                                      const igraph_vector_bool_t *fixed) nogil

    # Connected components / clusters (returns number of components in cno)
    ctypedef enum igraph_connectedness_t:
        IGRAPH_WEAK = 0
        IGRAPH_STRONG = 1

    igraph_error_t igraph_clusters(const igraph_t *graph,
                                   igraph_vector_int_t *membership,
                                   igraph_vector_int_t *csize,
                                   igraph_integer_t *cno,
                                   igraph_connectedness_t mode) nogil
    
    # Transitivity / clustering (Barrat weighted clustering)
    ctypedef enum igraph_transitivity_mode_t:
        IGRAPH_TRANSITIVITY_ZERO = 0
        IGRAPH_TRANSITIVITY_NA = 1

    igraph_error_t igraph_transitivity_barrat(const igraph_t *graph,
                                              igraph_vector_t *res,
                                              igraph_vs_t vids,
                                              const igraph_vector_t *weights,
                                              igraph_transitivity_mode_t mode) nogil

    # Vertex selector functions - THESE WERE MISSING
    igraph_error_t igraph_vs_all(igraph_vs_t *vs) nogil
    void igraph_vs_destroy(igraph_vs_t *vs) nogil
    igraph_vs_t igraph_vss_all() nogil
    igraph_error_t igraph_vs_vector(igraph_vs_t *vs, const igraph_vector_int_t *v) nogil

    # Graph analysis functions - THESE WERE MISSING
    ctypedef enum igraph_neimode_t:
        IGRAPH_OUT = 1
        IGRAPH_IN = 2
        IGRAPH_ALL = 3
    
    igraph_error_t igraph_degree(const igraph_t *graph,
                                 igraph_vector_int_t *res,
                                 igraph_vs_t vids,
                                 igraph_neimode_t mode,
                                 igraph_bool_t loops) nogil
    
    igraph_error_t igraph_strength(const igraph_t *graph,
                                   igraph_vector_t *res,
                                   igraph_vs_t vids,
                                   igraph_neimode_t mode,
                                   igraph_bool_t loops,
                                   const igraph_vector_t *weights) nogil
    
    igraph_error_t igraph_neighbors(const igraph_t *graph,
                                    igraph_vector_int_t *neis,
                                    igraph_integer_t vid,
                                    igraph_neimode_t mode) nogil
    
    igraph_error_t igraph_get_eid(const igraph_t *graph,
                                  igraph_integer_t *eid,
                                  igraph_integer_t from_,
                                  igraph_integer_t to,
                                  igraph_bool_t directed,
                                  igraph_bool_t error) nogil

    # Directedness constants (used as booleans in some igraph calls)
    ctypedef enum igraph_directedness_t:
        IGRAPH_UNDIRECTED = 0
        IGRAPH_DIRECTED = 1
    
    # Alternative: Multilevel (Louvain) - faster but less accurate
    int igraph_community_multilevel(const igraph_t *graph,
                                   const igraph_vector_t *weights,
                                   igraph_real_t resolution,
                                   igraph_vector_int_t *membership,
                                   igraph_matrix_t *memberships,
                                   igraph_vector_t *modularity) nogil
    
    # Modularity calculation
    int igraph_modularity(const igraph_t *graph,
                         const igraph_vector_int_t *membership,
                         const igraph_vector_t *weights,
                         igraph_real_t resolution,
                         igraph_bool_t directed,
                         igraph_real_t *modularity) nogil
    
    # Error handling - newer igraph uses igraph_error_t
    cdef enum:
        IGRAPH_SUCCESS = 0
        IGRAPH_FAILURE = 1
        IGRAPH_ENOMEM = 2
    
    # Matrix type (needed for some functions)
    ctypedef struct igraph_matrix_t:
        pass
    
    int igraph_matrix_init(igraph_matrix_t *m, long int nrow, long int ncol) nogil
    void igraph_matrix_destroy(igraph_matrix_t *m) nogil
    
    # Edge selector type
    ctypedef struct igraph_es_t:
        pass
    
    ctypedef enum igraph_edgeseq_type_t:
        IGRAPH_ES_ALL = 0

    # Helper to get edge endpoints by edge id
    igraph_error_t igraph_edge(const igraph_t *graph, igraph_integer_t eid,
                              igraph_integer_t *from_, igraph_integer_t *to) nogil

    # Vector push back for int vectors
    igraph_error_t igraph_vector_int_push_back(igraph_vector_int_t *v, igraph_integer_t val) nogil
    igraph_error_t igraph_vector_push_back(igraph_vector_t *v, igraph_real_t val) nogil

    # Convert vector_int to edge sequence
    igraph_error_t igraph_es_vector(igraph_es_t *es, const igraph_vector_int_t *v) nogil

    # Create subgraph from edge sequence
    igraph_error_t igraph_subgraph_from_edges(const igraph_t *graph, igraph_t *res,
                                             igraph_es_t es, igraph_bool_t delete_vertices) nogil
    
    # Induced subgraph
    ctypedef enum igraph_subgraph_implementation_t:
        IGRAPH_SUBGRAPH_AUTO = 0
        IGRAPH_SUBGRAPH_COPY_AND_DELETE = 1
        IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH = 2
    
    igraph_error_t igraph_induced_subgraph(const igraph_t *graph, igraph_t *res,
                                          igraph_vs_t vids,
                                          igraph_subgraph_implementation_t impl) nogil
    
    # Betweenness centrality
    igraph_error_t igraph_betweenness(const igraph_t *graph, igraph_vector_t *res,
                                     const igraph_vs_t vids, igraph_bool_t directed,
                                     const igraph_vector_t *weights) nogil