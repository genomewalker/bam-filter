# Generic Column-Based Filtering System

## Overview

This document describes the new generic column-based filtering system for filterBAM that replaces the old hardcoded filter arguments with a flexible, user-friendly syntax.

## Design Goals

1. **User-Friendly**: Single intuitive syntax instead of remembering many different flags
2. **Self-Documenting**: Column names are explicit and descriptive
3. **Extensible**: Easy to filter on any of the 45+ exported TSV columns
4. **Type-Safe**: Validates data types and value ranges
5. **Flexible**: Support for min-only, max-only, or range filters

## Syntax

### Format
```
--filter "column_name:min:max"
```

### Rules
- **Both bounds**: `column:10:100` → 10 ≤ value ≤ 100
- **Only min**: `column:10:` → value ≥ 10
- **Only max**: `column::100` → value ≤ 100
- **Multiple filters**: Comma-separated: `col1:10:,col2::100,col3:5:50`

### Examples
```bash
# Simple: read_ani_mean >= 90
--filter "read_ani_mean:90:"

# Maximum bound: breadth <= 0.95
--filter "breadth::0.95"

# Range: 30 <= read_length_mean <= 10000
--filter "read_length_mean:30:10000"

# Multiple filters
--filter "read_ani_mean:90:,coverage_mean:5:,breadth:0.4:0.95,coverage_evenness:0.1:"
```

## Column Names (Improved)

All column names have been improved for clarity and consistency:

### Old → New Mappings

#### Reference & Counts
- `reference` → `reference_name`
- `n_reads` → `read_count`
- `n_alns` → `alignment_count`
- `n_bins` → `bin_count`
- `n_reads_tad` → `read_count_tad`

#### Alignment Metrics
- `read_aligned_length` → `read_aligned_length_mean`
- `read_aln_score` → `read_alignment_score_mean`
- `mapping_quality` → `mapping_quality_mean`
- `edit_distances` → `edit_distance_mean`

#### Coverage Metrics
- `max_covered_bases` → `bases_covered_max`
- `mean_covered_bases` → `bases_covered_mean`
- `coverage_mean_trunc` → `coverage_mean_trimmed`
- `coverage_mean_trunc_len` → `coverage_mean_trimmed_length`
- `coverage_covered_mean` → `coverage_mean_covered_only`
- `bam_reference_length` → `reference_length_bam`

#### Breadth Metrics
- `exp_breadth` → `breadth_expected`
- `breadth_exp_ratio` → `breadth_expected_ratio`

#### Distribution Metrics
- `norm_entropy` → `entropy_normalized`
- `gini` → `gini_coefficient`
- `norm_gini` → `gini_coefficient_normalized`
- `c_v` → `coefficient_of_variation`
- `d_i` → `diversity_index`
- `cov_evenness` → `coverage_evenness`

#### Abundance Metrics
- `tax_abund_read` → `abundance_read_based`
- `tax_abund_aln` → `abundance_alignment_based`
- `tax_abund_tad` → `abundance_tad`

## Data Type Constraints

Each column has type constraints and valid ranges:

| Column | Type | Min | Max | Description |
|--------|------|-----|-----|-------------|
| `read_count` | int | 0 | ∞ | Number of reads |
| `read_ani_mean` | float | 0 | 100 | Mean ANI percentage |
| `coverage_mean` | float | 0 | ∞ | Mean coverage depth |
| `breadth` | float | 0 | 1 | Breadth of coverage (proportion) |
| `entropy_normalized` | float | 0 | 1 | Normalized entropy |
| `gini_coefficient_normalized` | float | 0 | 1 | Normalized Gini coefficient |
| `coefficient_of_variation` | float | 0 | ∞ | CV = SD/mean |
| `diversity_index` | float | 0 | 1 | Diversity index |
| `coverage_evenness` | float | 0 | 1 | Coverage evenness |

*(See `bam_filter/filter_parser.py` for complete specifications)*

## Implementation Architecture

### 1. Python Layer

**File**: `bam_filter/filter_parser.py`

Key components:
- `ColumnSpec`: Defines column properties (name, type, valid range)
- `FilterSpec`: Represents a single filter rule
- `parse_filter_string()`: Parses user input into structured filters
- `export_filters_to_cython_format()`: Converts to Cython-compatible format

**Validation**:
- Column name existence (with typo suggestions)
- Data type matching
- Value range checking
- Min < max validation

### 2. Cython Layer

**Files**:
- `bam_filter/generic_filters.pxd` (declarations)
- `bam_filter/generic_filters.pyx` (implementation)

Key components:
- `ColumnFilter`: C struct for filter rule (column_index, min, max)
- `GenericFilters`: Container for all filters
- `passes_generic_filters()`: Fast nogil filtering function
- `get_column_value()`: Extracts column value from RefStats

**Performance**:
- All filtering done in nogil context
- Direct struct access (no Python overhead)
- Column values accessed by index for speed

### 3. Integration

**Modified Files**:
- `bam_filter/utils.py`: CLI arguments (replaced 11 args with 1)
- `bam_filter/filter.py`: Parse filters and pass to Cython
- `bam_filter/stats_io.pyx`: Updated TSV headers
- `bam_filter/stats.pyx`: Accept generic filters parameter

## CLI Interface

### New Arguments

```bash
filterBAM filter --bam input.bam \
  --stats output.tsv \
  --filter "read_ani_mean:90:,coverage_mean:5:"
```

### List Available Columns

```bash
filterBAM filter --list-columns
```

Output:
```
=== Filterable Columns ===

Format: column_name [type] (min, max)

Read Statistics:
  read_count                          [int] (min=0)
  alignment_count                     [int] (min=0)
  read_length_mean                    [float] (min=0, max=100000)
  ...

Coverage Metrics:
  coverage_mean                       [float] (min=0)
  breadth                             [float] (min=0, max=1)
  coverage_evenness                   [float] (min=0, max=1)
  ...

Example usage:
  --filter 'read_ani_mean:90:'
  --filter 'breadth::0.95,coverage_mean:5:'
```

### Comparison: Old vs New

**Old (removed)**:
```bash
filterBAM filter --bam input.bam \
  -A 90 -l 30 -L 10000 -n 1 \
  -b 0 -e 0 -g 1.0 -B 0 \
  -a 90 -c 0 -V inf -C 0
```

**New**:
```bash
filterBAM filter --bam input.bam \
  --filter "read_ani_mean:90:,read_count:1:,breadth_expected_ratio:0:,\
entropy_normalized:0:,gini_coefficient_normalized::1.0,coverage_mean:0:,\
coverage_evenness:0:"
```

## Error Handling

### Example Errors

**Unknown column**:
```
ERROR: Unknown column: 'read_ani'
Did you mean: read_ani_mean, read_ani_std, read_ani_median?
```

**Invalid syntax**:
```
ERROR: Invalid filter syntax: 'coverage_mean:5'
Expected format: 'column:min:max' (use empty for no bound)
Examples:
  'read_ani_mean:90:'      (>= 90)
  'breadth::0.95'          (<= 0.95)
  'coverage_mean:5:1000'   (between 5 and 1000)
```

**Type mismatch**:
```
ERROR: Invalid minimum value 'abc' for column 'coverage_mean'
Expected float type
```

**Range violation**:
```
ERROR: Minimum value 150 for 'read_ani_mean' is above valid range
(must be <= 100)
```

**Invalid range**:
```
ERROR: Invalid range for 'coverage_mean': min (100) must be < max (10)
```

## Migration Notes

### Breaking Changes

1. **Removed CLI arguments** (no backward compatibility):
   - `-A, --min-read-ani` (for reference filtering - read-level still exists)
   - `-n, --min-read-count`
   - `-b, --min-expected-breadth-ratio`
   - `-e, --min-normalized-entropy`
   - `-g, --min-normalized-gini`
   - `-B, --min-breadth`
   - `-a, --min-avg-read-ani`
   - `-c, --min-coverage-evenness`
   - `-V, --min-coeff-var`
   - `-C, --min-coverage-mean`
   - `--include-low-detection`

2. **Changed TSV column headers**: All output TSV files now use new column names

3. **Kept arguments** (read-level filters):
   - `-A, --min-read-ani`: Minimum read ANI (filters individual reads, not references)
   - `-l, --min-read-length`: Minimum read length
   - `-L, --max-read-length`: Maximum read length

### Migration Examples

```bash
# Old
filterBAM filter --bam in.bam -a 90 -c 0.1

# New
filterBAM filter --bam in.bam --filter "read_ani_mean:90:,coverage_evenness:0.1:"
```

```bash
# Old
filterBAM filter --bam in.bam -B 0.5 -C 5

# New
filterBAM filter --bam in.bam --filter "breadth:0.5:,coverage_mean:5:"
```

## Testing

### Unit Tests

Test file: `tests/test_filter_parser.py`

```python
def test_parse_single_filter():
    filters = parse_filter_string("read_ani_mean:90:")
    assert "read_ani_mean" in filters
    assert filters["read_ani_mean"].min == 90
    assert filters["read_ani_mean"].max is None

def test_parse_range_filter():
    filters = parse_filter_string("coverage_mean:5:100")
    assert filters["coverage_mean"].min == 5
    assert filters["coverage_mean"].max == 100

def test_unknown_column():
    with pytest.raises(ValueError, match="Unknown column"):
        parse_filter_string("invalid_col:10:")

def test_type_validation():
    with pytest.raises(ValueError, match="Expected int type"):
        parse_filter_string("read_count:10.5:")
```

### Integration Tests

```bash
# Test basic filtering
filterBAM filter --bam test.bam --stats output.tsv \
  --filter "read_ani_mean:90:,coverage_mean:5:"

# Test list columns
filterBAM filter --list-columns

# Test error handling
filterBAM filter --bam test.bam --stats out.tsv \
  --filter "invalid_column:10:"  # Should error with suggestion
```

## Future Enhancements

1. **Logical operators**: Support AND/OR combinations
   ```bash
   --filter "(breadth:0.5: AND coverage_mean:5:) OR read_count:1000:"
   ```

2. **Computed columns**: Filter on derived metrics
   ```bash
   --filter "breadth/breadth_expected:0.8:"  # Breadth ratio
   ```

3. **Percentile filters**: Filter by quantiles
   ```bash
   --filter "coverage_mean:p25:p75"  # Middle 50%
   ```

4. **Filter presets**: Named filter combinations
   ```bash
   --filter-preset high_quality  # Predefined quality filters
   ```

## Summary

The new generic filtering system provides:

✅ **Improved UX**: Descriptive column names, clear syntax
✅ **Flexibility**: Filter on any column with min/max/range
✅ **Type Safety**: Validates types and ranges
✅ **Discoverability**: `--list-columns` shows all options
✅ **Performance**: Fast Cython implementation with nogil
✅ **Maintainability**: Easy to add new filterable columns

**No backward compatibility**: Clean break for better long-term design.
