# BAM to Parquet Conversion

## Overview

The `to-parquet` subcommand converts coordinate-sorted BAM files to partitioned Parquet format optimized for fast querying with DuckDB. This is especially useful for metagenomic datasets with millions of references and billions of alignments.

## Features

- **Metagenomic-scale**: Handles billions of alignments with constant memory usage
- **Hash-partitioned**: 256 partitions for even distribution and parallel processing
- **DuckDB-optimized**: Sorted data, column statistics, and compression for fast queries
- **Lossless**: Preserves all alignment information with calculated metrics (ANI, scores)
- **Space-efficient**: 70-85% compression vs original BAM

## Basic Usage

```bash
# Convert BAM to Parquet (default settings)
filterBAM to-parquet \
    --bam input.bam \
    --output output.parquet

# With custom settings
filterBAM to-parquet \
    --bam input.bam \
    --output output.parquet \
    --num-partitions 256 \
    --compression zstd \
    --compression-level 5 \
    --threads 4 \
    --min-mapq 20
```

## Output Structure

```
output.parquet/
├── _metadata.json                 # Conversion metadata and statistics
├── references/
│   └── references.parquet         # Reference ID → name mapping
└── alignments/
    ├── ref_partition=0000/
    │   └── data.parquet          # Alignments for refs 0, 256, 512, ...
    ├── ref_partition=0001/
    │   └── data.parquet          # Alignments for refs 1, 257, 513, ...
    └── ...
```

## Schema

### Alignments Table
- `read_id` (uint32): Sequential read identifier
- `ref_id` (uint32): Reference ID (join with references table)
- `position` (int32): 0-based leftmost position
- `end_position` (int32): 0-based rightmost position
- `mapq` (uint8): Mapping quality (0-255)
- `flag` (uint16): SAM flags
- `ani` (float32): Alignment identity percentage
- `alignment_score` (float32): Log-likelihood score
- `pmd_score` (float32): PMD score (nullable)
- `num_mismatches` (uint16): Edit distance (NM tag)
- `alignment_length` (uint16): Alignment length
- `template_length` (int32): Template length (TLEN)
- `mate_ref_id` (int32): Mate reference ID
- `mate_position` (int32): Mate position
- `read_name` (string): Read name (optional)
- `cigar` (string): CIGAR string
- `sequence` (string): Read sequence (optional)
- `quality` (binary): Quality scores

### References Table
- `ref_id` (uint32): Reference identifier
- `ref_name` (string): Reference name
- `ref_length` (uint32): Reference length
- `ref_partition` (uint16): Partition number

## Querying with DuckDB

### Python

```python
import duckdb

# Connect to DuckDB
conn = duckdb.connect()

# Query 1: Coverage stats per reference
result = conn.execute("""
    SELECT 
        r.ref_name,
        COUNT(*) as num_alignments,
        AVG(a.ani) as mean_ani,
        AVG(a.mapq) as mean_mapq
    FROM 'output.parquet/alignments/**/*.parquet' a
    JOIN 'output.parquet/references/references.parquet' r 
      ON a.ref_id = r.ref_id
    WHERE a.mapq >= 20
    GROUP BY r.ref_name
    ORDER BY num_alignments DESC
    LIMIT 100
""").df()

# Query 2: Alignments in specific region
result = conn.execute("""
    SELECT a.*, r.ref_name
    FROM 'output.parquet/alignments/ref_partition=0042/*.parquet' a
    JOIN 'output.parquet/references/references.parquet' r USING (ref_id)
    WHERE a.ref_id = 10666
      AND a.position BETWEEN 1000000 AND 2000000
      AND a.ani >= 95
    ORDER BY a.position
""").df()

# Query 3: Multi-mapping reads
result = conn.execute("""
    SELECT 
        read_id,
        COUNT(DISTINCT ref_id) as num_refs,
        MAX(alignment_score) - MIN(alignment_score) as score_range
    FROM 'output.parquet/alignments/**/*.parquet'
    GROUP BY read_id
    HAVING COUNT(*) > 1
    ORDER BY num_refs DESC
    LIMIT 1000
""").df()
```

### DuckDB CLI

```sql
-- Load DuckDB CLI
duckdb

-- Query directly
SELECT 
    r.ref_name,
    COUNT(*) as coverage
FROM 'output.parquet/alignments/**/*.parquet' a
JOIN 'output.parquet/references/references.parquet' r USING (ref_id)
WHERE a.mapq >= 30 AND a.ani >= 95
GROUP BY r.ref_name
ORDER BY coverage DESC;

-- Export subset to CSV
COPY (
    SELECT * FROM 'output.parquet/alignments/**/*.parquet'
    WHERE ref_id IN (SELECT ref_id FROM 'output.parquet/references/references.parquet' 
                     WHERE ref_name LIKE 'Escherichia_coli%')
) TO 'ecoli_alignments.csv' WITH (HEADER);
```

## Command-line Options

### Output Options
- `--output`, `-o`: Output directory (required)

### Format Options
- `--num-partitions`: Number of hash partitions (default: 256)
- `--batch-size`: Alignments per batch (default: 100,000)
- `--compression`: Codec - snappy, zstd, gzip, none (default: zstd)
- `--compression-level`: Compression level for zstd (default: 3)
- `--include-read-names`: Include read names (increases size ~30-50%)
- `--no-sequences`: Exclude sequences (reduces size ~20-30%)

### Filter Options
- `--min-mapq`: Minimum mapping quality (default: 0)
- `--min-read-length`: Minimum read length (default: 0)
- `--max-read-length`: Maximum read length (default: 100,000)
- `--min-read-ani`: Minimum alignment identity % (default: 0.0)

### Common Options
- `--threads`, `-t`: Number of threads for BAM reading (default: 1)
- `--verbose`, `-v`: Increase verbosity

## Performance

### Typical Metagenomic Dataset
```
Input:
- 10 billion alignments
- 10 million references
- 500 GB BAM file

Output:
- ~100 GB Parquet (80% reduction)
- ~600 MB references table
- Conversion time: 2-4 hours (streaming)
- Memory usage: ~200 MB (constant)

Query Performance:
- Full scan: ~30 seconds (vs 4 hours for BAM)
- Partition-filtered: ~1 second
- Position range: <100ms (sorted + statistics)
```

## Tips

1. **Choose appropriate partitions**: 256 works well for 1M-100M references
2. **Adjust compression**: Use zstd level 3-5 for balance, level 7-9 for max compression
3. **Omit read names**: Save 30-50% space if you don't need them
4. **Keep sequences**: Usually worth keeping for analysis
5. **Filter on write**: Use `--min-mapq` to exclude low-quality alignments early
6. **Use DuckDB**: Much faster than loading Parquet into pandas/polars for large datasets

## Examples

### Example 1: Quick conversion
```bash
filterBAM to-parquet --bam input.bam --output output.parquet
```

### Example 2: High-quality alignments only
```bash
filterBAM to-parquet \
    --bam input.bam \
    --output output.parquet \
    --min-mapq 30 \
    --min-read-ani 95.0 \
    --compression zstd \
    --compression-level 5
```

### Example 3: Space-optimized (no read names or sequences)
```bash
filterBAM to-parquet \
    --bam input.bam \
    --output output.parquet \
    --no-sequences \
    --compression zstd \
    --compression-level 9
```

### Example 4: Fast processing with filtering
```bash
filterBAM to-parquet \
    --bam input.bam \
    --output output.parquet \
    --threads 8 \
    --batch-size 200000 \
    --compression snappy \
    --min-read-length 30
```

## Troubleshooting

**Q: Conversion is slow**
A: Increase `--batch-size` to 200,000 or use snappy compression

**Q: Output is too large**
A: Use `--no-sequences`, increase `--compression-level`, or filter more aggressively

**Q: Out of memory**
A: Reduce `--batch-size` to 50,000 (should use <100MB regardless)

**Q: Can't query specific region**
A: Make sure you're querying the correct partition (use modulo to find it)

## See Also

- [PyArrow Parquet documentation](https://arrow.apache.org/docs/python/parquet.html)
- [DuckDB documentation](https://duckdb.org/docs/)
- [filterBAM filter command](README.md) - Filter references by coverage
