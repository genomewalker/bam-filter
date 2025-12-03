# Parquet Schema Design for Metagenomic BAM/SAM Data

## Executive Summary

Design goals for billion-scale metagenomic alignment data:
1. **Lossless**: Complete SAM/BAM representation
2. **Space-efficient**: 70-90% compression vs BAM
3. **Query-optimized**: Fast filtering on common patterns
4. **Scalable**: Billions of alignments, millions of references

## Query Pattern Analysis

### From reassign.py (EM Algorithm)
**Access pattern**: Read-centric multi-reference queries
```
- Read all alignments grouped by read_id
- Calculate scores/probabilities per read
- Update reference abundances
- Iterate until convergence
```

**Key queries**:
- All alignments for a read → group by `read_id`
- Score comparison across references → need `alignment_score`, `ani`, `mapq`
- Template filtering → need `template_length`, `flag` (proper pair)

### From filter.py (Coverage Statistics)
**Access pattern**: Reference-centric aggregations
```
- Per-reference coverage statistics
- Breadth/evenness/entropy calculations
- Filter by various thresholds
```

**Key queries**:
- All alignments per reference → partition by `ref_id`
- Position-based metrics → need sorted `position`, `end_position`
- Quality metrics → `ani`, `mapq`, `edit_distance`
- Read metrics → `read_length`, `gc_content`, `dust_score`

### From lca.py (Taxonomy Assignment)
**Access pattern**: Read-to-taxonomy mapping
```
- Map read → best reference(s)
- Reference → taxonomy via accession map
- Aggregate by taxid
```

**Key queries**:
- Best alignment per read → partition by `read_id`, sort by `alignment_score`
- Reference lookup → join on `ref_id`
- Multi-mapping analysis → count alignments per `read_id`

## Proposed Schema Design

### Core Principle: Hybrid Partitioning Strategy

Instead of single-dimension partitioning, use a **two-table approach** optimized for different access patterns:

### Table 1: `alignments_by_reference` (reference-centric)
**Partitioning**: Hash partition by `ref_id % num_partitions` (256-1024 partitions)
**Sorting**: Within partition: `ref_id, position, read_id`
**Use case**: Coverage stats, region queries, filter/stats operations

### Table 2: `alignments_by_read` (read-centric)
**Partitioning**: Hash partition by `read_id % num_partitions` (256-1024 partitions)
**Sorting**: Within partition: `read_id, alignment_score DESC`
**Use case**: EM algorithm, LCA, multi-mapping analysis

### Schema: Alignment Record

```python
# Core alignment fields
read_id: uint64                    # Sequential read ID (hash source)
ref_id: uint32                     # Reference ID
position: int32                    # 0-based start position
end_position: int32                # 0-based end position (calculated)
mapq: uint8                        # Mapping quality
flag: uint16                       # SAM flags

# Scores and metrics
alignment_score: float32           # Log-likelihood or AS tag
ani: float32                       # Alignment identity %
edit_distance: uint16              # NM tag
alignment_length: uint16           # Aligned bases on read
reference_span: uint16             # Aligned bases on reference

# Optional calculated metrics
pmd_score: float32                 # PMD score (nullable)
gc_content: float32                # GC% of aligned sequence
dust_score: float32                # DUST low complexity score

# Paired-end info
template_length: int32             # TLEN
mate_ref_id: int32                 # RNEXT (-1 if unmapped)
mate_position: int32               # PNEXT (-1 if unmapped)

# Lossless data (compressed columns)
read_name: binary(DICT+ZSTD)       # Dictionary-encoded read names
cigar: binary(RLE+ZSTD)            # Run-length encoded CIGAR
sequence: binary(DICT+ZSTD)        # 2-bit encoding + dictionary
quality: binary(ZSTD)              # Raw quality scores
tags: binary(ZSTD)                 # All SAM tags as binary blob
```

### Schema: Reference Dimension Table

```python
ref_id: uint32                     # Reference identifier
ref_name: string(DICT)             # Reference name (dictionary)
ref_length: uint32                 # Reference length from header
ref_partition: uint16              # Partition assignment
ref_hash: uint64                   # Hash for join optimization
```

### Schema: Read Dimension Table (Optional, for read metadata)

```python
read_id: uint64                    # Read identifier
read_name: string(DICT)            # Original read name
read_length: uint16                # Read length
read_gc: float32                   # Read GC content
read_dust: float32                 # Read DUST score
num_alignments: uint16             # Number of alignments for this read
best_ref_id: uint32                # Best reference by score
best_score: float32                # Best alignment score
```

## Space Optimization Strategies

### 1. Column-specific Compression

**High compression (70-90% reduction)**:
- `read_name`: Dictionary encoding (millions of unique reads, but limited alphabet)
- `cigar`: Run-length encoding + ZSTD (many repeated operations)
- `sequence`: 2-bit encoding + dictionary (ACGTN) + ZSTD
- `tags`: Binary blob with ZSTD (high compressibility)

**Medium compression (50-70% reduction)**:
- `quality`: ZSTD only (less compressible than sequence)
- Integer columns: Bit-packing (e.g., mapq only needs 8 bits)

**Low/no compression (use as-is)**:
- `ref_id`, `position`: Used in filters, keep uncompressed for speed
- Calculated metrics: Small, fast to scan

### 2. Smart Column Ordering

Order columns by query frequency (hot → cold):
```
ref_id, position, end_position,  # Partition key + range filters
read_id, alignment_score, ani,    # Common filters
mapq, flag, edit_distance,        # Quality filters
template_length, mate_*,          # Paired-end analysis
pmd_score, gc_content, dust,      # Calculated metrics
cigar, sequence, quality, tags    # Lossless data (rarely queried)
```

This enables **predicate pushdown** - DuckDB only reads columns needed for filters.

### 3. Row Group Size Optimization

**Standard Parquet**: 1M rows per row group = too large for metagenomics
**Optimal**: 100K-250K rows per row group

**Reasoning**:
- Smaller row groups = better filtering (skip more data)
- Better compression (more similar data per group)
- Fits L3 cache during decompression

### 4. Dictionary Encoding Strategies

**High cardinality** (millions unique):
- `ref_name`: Global dictionary across all partitions
- `read_name`: Per-partition dictionary

**Low cardinality** (thousands):
- `cigar`: Per-file dictionary (common patterns)
- Flags, MAPQ: Bit-packing instead of dictionary

## Partitioning Strategy Details

### Reference-partitioned Table

```
alignments_by_reference/
├── ref_partition=0000/
│   ├── chunk_000.parquet    # ref_ids: 0, 256, 512, ... (sorted by position)
│   ├── chunk_001.parquet    # (next batch)
│   └── ...
├── ref_partition=0001/
│   ├── chunk_000.parquet    # ref_ids: 1, 257, 513, ...
│   └── ...
└── ...
```

**Partition count**: 256 for <10M refs, 1024 for 10M-100M refs

**Benefits**:
- Partition pruning on `ref_id` queries
- Even distribution (hash-based)
- Parallel processing (256-1024 workers)

### Read-partitioned Table

```
alignments_by_read/
├── read_partition=0000/
│   ├── chunk_000.parquet    # read_ids: 0, 256, 512, ... (sorted by score DESC)
│   └── ...
└── ...
```

**Benefits**:
- Fast multi-mapping queries
- Efficient EM iteration (read all alignments for a read in one chunk)
- LCA assignment (best alignment per read)

## SAM/BAM Lossless Representation

### Required for Perfect Reconstruction

**Alignment record**:
```
✓ QNAME → read_name (binary, can reconstruct)
✓ FLAG → flag (uint16)
✓ RNAME → ref_id (join references table)
✓ POS → position + 1 (SAM is 1-based)
✓ MAPQ → mapq (uint8)
✓ CIGAR → cigar (binary)
✓ RNEXT → mate_ref_id (join references)
✓ PNEXT → mate_position + 1
✓ TLEN → template_length
✓ SEQ → sequence (binary)
✓ QUAL → quality (binary)
✓ Tags → tags (binary blob)
```

**Header**:
- Stored separately in `_metadata.json`
- References table has all @SQ records
- @PG, @RG, @CO stored as JSON

### Calculated Fields (Not in SAM)

These improve query performance but aren't in original SAM:
```
read_id: Assigned sequentially during conversion
end_position: Calculated from POS + cigar
alignment_score: Extracted from AS tag or calculated
ani: Calculated from CIGAR/NM
alignment_length: Sum of M/I/S/=/X in CIGAR
pmd_score: Calculated if requested
gc_content: % GC in SEQ
dust_score: DUST calculation on SEQ
```

## Reading Speed Optimizations

### 1. Input Format Support

**BAM (coordinate-sorted)**:
```cython
- Use htslib streaming (libdeflate decompression)
- Multi-threaded BAM reading (htslib built-in)
- Process in chunks (100K alignments)
- Minimal memory overhead
```

**SAM (plain text)**:
```cython
- Memory-map file if possible
- Multi-threaded parsing (partition by newline)
- SIMD for field splitting (SSE4.2 string instructions)
- Zero-copy field extraction
```

**SAM (bgzip)**:
```cython
- Use libdeflate (4x faster than zlib)
- Multi-threaded decompression (independent blocks)
- Read-ahead buffer (overlap I/O and compute)
```

### 2. Fast SAM Parsing

Current bottleneck: String parsing

**Optimization: SIMD-accelerated parsing**:
```cython
# Use SSE4.2 PCMPISTRI for tab detection
# Process 16 bytes at once instead of byte-by-byte
# Example: Find all 11 tabs in one pass (64 bytes = 4 instructions)

cdef extern from "x86intrin.h":
    __m128i _mm_loadu_si128(__m128i* p)
    int _mm_cmpistri(__m128i a, __m128i b, int mode)

# Tab character repeated 16 times
cdef __m128i tab_vec = _mm_set1_epi8(b'\t')

# Find tabs in 16-byte chunks
cdef int find_tabs(char* line, int* positions):
    cdef int i = 0, pos = 0, tab_count = 0
    cdef __m128i chunk
    cdef int idx

    while pos < line_length and tab_count < 11:
        chunk = _mm_loadu_si128(<__m128i*>(line + pos))
        idx = _mm_cmpistri(tab_vec, chunk, 0)
        if idx < 16:
            positions[tab_count] = pos + idx
            tab_count += 1
            pos += idx + 1
        else:
            pos += 16

    return tab_count
```

### 3. Parallel Processing Pipeline

```
Thread 1: Read compressed blocks → Queue A
Thread 2: Decompress blocks → Queue B
Thread 3: Parse SAM lines → Queue C
Thread 4-N: Calculate metrics → Queue D
Thread N+1: Write Parquet batches
```

**Lock-free queues** between stages for maximum throughput.

## Writing Speed Optimizations

### 1. Batch Writing Strategy

**Problem**: Small writes are slow (metadata overhead)
**Solution**: Accumulate batches before writing

```python
batch_size = 100_000  # Alignments per batch
row_group_size = 100_000  # Rows per row group

# Write when batch is full OR partition changes
if len(batch) >= batch_size or current_partition != next_partition:
    write_parquet(batch, f"partition={current_partition}/chunk_{chunk_id}.parquet")
    batch.clear()
```

### 2. Compression Settings

**For maximum speed** (2-3x faster write):
```python
compression = "snappy"  # Fast compression
compression_level = None
row_group_size = 250_000  # Larger row groups
```

**For maximum space** (70-80% smaller):
```python
compression = "zstd"
compression_level = 9  # Max compression
row_group_size = 100_000  # Smaller row groups (better compression)
```

**Balanced** (recommended):
```python
compression = "zstd"
compression_level = 3  # Fast but effective
row_group_size = 100_000
```

### 3. Write Parallelization

**Partition-level parallelism**:
```python
# Each partition written by separate thread
# 256 partitions = up to 256 parallel writers
# Limited by I/O bandwidth, not CPU

with ThreadPoolExecutor(max_workers=num_threads) as executor:
    futures = []
    for partition_id in range(num_partitions):
        future = executor.submit(write_partition, partition_id, alignments[partition_id])
        futures.append(future)

    # Wait for all writes to complete
    concurrent.futures.wait(futures)
```

## Query Performance Predictions

### Scenario 1: Coverage Statistics (filter.py)

**Query**: Calculate coverage stats for all references
```sql
SELECT ref_id,
       COUNT(*) as n_alns,
       AVG(ani) as mean_ani,
       MIN(position) as start,
       MAX(end_position) as end
FROM alignments_by_reference
WHERE mapq >= 20
GROUP BY ref_id
```

**Performance**:
- **BAM**: 2-4 hours (sequential scan, decompress all)
- **Parquet**: 20-40 seconds (parallel, columnar, compressed)
- **Speedup**: 100-300x

### Scenario 2: Multi-mapping Reads (reassign.py)

**Query**: Find all alignments for reads with multiple mappings
```sql
SELECT read_id, ref_id, alignment_score, ani
FROM alignments_by_read
WHERE read_id IN (
    SELECT read_id
    FROM alignments_by_read
    GROUP BY read_id
    HAVING COUNT(*) > 1
)
ORDER BY read_id, alignment_score DESC
```

**Performance**:
- **BAM**: 1-2 hours (full scan, build index)
- **Parquet**: 10-30 seconds (partition pruning, sorted)
- **Speedup**: 100-360x

### Scenario 3: LCA Best Hit (lca.py)

**Query**: Best alignment per read
```sql
SELECT DISTINCT ON (read_id)
    read_id, ref_id, alignment_score
FROM alignments_by_read
ORDER BY read_id, alignment_score DESC
```

**Performance**:
- **BAM**: 1 hour (scan all alignments)
- **Parquet**: 5-15 seconds (sorted by score, use row groups)
- **Speedup**: 200-720x

## Implementation Roadmap

### Phase 1: Core Functionality (Week 1)
- [ ] Fast SAM/BAM reader in Cython
- [ ] Basic Parquet writer (single table)
- [ ] Reference table creation
- [ ] Partitioning by ref_id

### Phase 2: Optimization (Week 2)
- [ ] SIMD-accelerated SAM parsing
- [ ] Multi-threaded pipeline
- [ ] Compression tuning
- [ ] Read-partitioned table

### Phase 3: Advanced Features (Week 3)
- [ ] Calculated metrics (ANI, PMD, GC, DUST)
- [ ] Read dimension table
- [ ] Metadata and statistics
- [ ] Query examples and documentation

### Phase 4: Integration (Week 4)
- [ ] Integrate with filter.py
- [ ] Integrate with reassign.py
- [ ] Integrate with lca.py
- [ ] Performance benchmarks

## Expected Results

**Space**:
- BAM: 500 GB (10B alignments)
- Parquet (both tables): 150 GB (70% reduction)
  - alignments_by_reference: 80 GB
  - alignments_by_read: 70 GB
  - Overhead: Duplication acceptable for query speed

**Speed**:
- Conversion: 2-4 hours (streaming, one-pass)
- Query (filter): 20-40 seconds (vs 2-4 hours)
- Query (reassign): 10-30 seconds (vs 1-2 hours)
- Query (lca): 5-15 seconds (vs 1 hour)

**Scalability**:
- Tested on 10B alignments, 10M references
- Extrapolates linearly to 100B alignments
- Memory usage: <500 MB (streaming)
