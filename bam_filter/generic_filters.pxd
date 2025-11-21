# cython: language_level=3
# Generic column-based filtering declarations

from bam_filter.stats cimport RefStats

# Maximum number of filters
cdef int MAX_FILTERS

# Structure to hold a single column filter rule
cdef struct ColumnFilter:
    int column_index      # Which column to filter (0-44 for the 45 columns)
    double min_value      # Minimum value (use -inf for no lower bound)
    double max_value      # Maximum value (use +inf for no upper bound)
    bint is_active        # Whether this filter is active

# Structure to hold all column filters
cdef struct GenericFilters:
    ColumnFilter* filters  # Array of filter rules
    int n_filters         # Number of active filters
    int capacity          # Allocated capacity

# Function declarations
cdef GenericFilters* create_generic_filters(int initial_capacity) noexcept nogil
cdef void destroy_generic_filters(GenericFilters* gf) noexcept nogil
cdef int add_filter(GenericFilters* gf, int column_index, double min_val, double max_val) noexcept nogil
cdef bint passes_generic_filters(RefStats* stats, GenericFilters* gf) noexcept nogil
cdef double get_column_value(RefStats* stats, int column_index) noexcept nogil
