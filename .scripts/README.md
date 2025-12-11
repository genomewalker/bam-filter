# Test Scripts

This directory contains test and validation scripts for filterBAM development.

## Scripts

### test_bam_to_parquet.py
Tests the high-performance BAM to Parquet conversion using Cython/C++ with Arrow.

**Features:**
- Single-worker mode for debugging
- Parallel mode with ProcessPoolExecutor (record-based splitting)
- PMD score calculation on-the-fly
- Verification of output Parquet files

**Usage:**
```bash
# Quick test with single worker
python .scripts/test_bam_to_parquet.py input.bam /tmp/parquet_test --single

# Full parallel test (8 workers)
python .scripts/test_bam_to_parquet.py input.bam /tmp/parquet_test --threads 8 --parallel

# With verification
python .scripts/test_bam_to_parquet.py input.bam /tmp/parquet_test --verify

# Quick test with limited records
python .scripts/test_bam_to_parquet.py input.bam /tmp/parquet_test --max-records 100000 --single
```

**Architecture:**
```
CLI (bam_to_parquet.py)
    └── parallel_parquet_records.py (ProcessPoolExecutor)
        └── parquet_converter_pure_cpp.pyx (Cython, nogil HTSlib)
            └── arrow_parquet_writer.cpp (Arrow C++ API, ZSTD)
```

### test_bam_to_parquet_fast.py
Fast benchmark for BAM to Parquet conversion without PMD calculation.

**Usage:**
```bash
python .scripts/test_bam_to_parquet_fast.py input.bam /tmp/output_dir
```

### test_bam_to_parquet_noseq.py
Maximum speed benchmark - no PMD, no sequence/quality storage.

**Usage:**
```bash
python .scripts/test_bam_to_parquet_noseq.py input.bam /tmp/output_dir
```

## Performance Benchmarks

Tested on 75.8M records (2.7 GB BAM file):

| Mode | Time | Throughput | Speedup |
|------|------|------------|---------|
| Full (PMD + both partitions) | 1023s | 74K rec/s | 1.0x |
| No PMD, by_reference only | 625s | 121K rec/s | 1.6x |
| No PMD, nogil optimized | 591s | 128K rec/s | 1.7x |
| No sequences (push_back) | 393s | 193K rec/s | 2.6x |
| **No sequences + Zero-copy Arrow** | **323s** | **234K rec/s** | **3.2x** |

**Key optimizations:**
1. **nogil tag extraction** - All tag extraction runs in a single `nogil` block
2. **HTSlib nogil declarations** - `bam_aux_get`, `bam_aux2i`, etc. marked nogil
3. **Optional sequences** - `store_sequences=False` skips sequence/quality copying
4. **Zero-copy Arrow** - Primitive arrays wrapped directly from std::vector (O(1))
5. **ZSTD compression** - Dictionary encoding for low-cardinality columns

### SAM.gz Benchmarks

Tested on 924M records (11.7 GB SAM.gz file):

| Parser | Workers | Time | Throughput | Notes |
|--------|---------|------|------------|-------|
| HTSlib parallel (convert_bgzf_range) | 4 | 1013s | 686K rec/s | Uses sam_read1() |
| **Fast Cython (sam_parser_fast)** | 8 | **180s** | **5.14M rec/s** | Bypasses HTSlib text parsing |

**Speedup: 7.5x** by parsing SAM text directly in Cython instead of using HTSlib's `sam_read1()`.

### test_parallel_sam_fast.py
Fastest SAM.gz to Parquet conversion using custom Cython parser.

**Usage:**
```bash
python .scripts/test_parallel_sam_fast.py input.sam.gz /tmp/output 8
```

**Architecture:**
```
parallel_sam_fast.py (ProcessPoolExecutor, BGZF block splitting)
    └── sam_parser_fast.pyx (Cython, direct line parsing)
        └── bgzf_getline() (HTSlib BGZF decompression)
        └── parse_sam_line_nogil() (Pure C parsing, nogil)
        └── arrow_parquet_writer.cpp (Arrow C++ API, ZSTD)
```

### test_lossless_sam.py
Tests lossless SAM to Parquet conversion with `tags_raw` column preservation.

**Usage:**
```bash
python .scripts/test_lossless_sam.py input.sam.gz /scratch/tmp/output [max_records]
```

### test_tags_raw_quick.py
Quick verification test for `tags_raw` column - extracts 1000 records and verifies
all SAM optional tags are preserved in the Parquet output.

**Usage:**
```bash
python .scripts/test_tags_raw_quick.py input.sam.gz
```

## Parquet Schema Optimization

See `parquet_schema_proposal.md` for the full design document.

### test_optimized_schema.py
Compares storage efficiency of different Parquet schema designs.

**Usage:**
```bash
python .scripts/test_optimized_schema.py input.sam.gz [num_records]
```

### test_optimized_writer.py
Tests the super-optimized Parquet writer with:
- **LUT-based 2-bit sequence packing** (branchless, ~4x compression)
- **Zero-copy Python string handling** (`PyUnicode_AsUTF8AndSize`)
- **Bulk Arrow AppendValues API** (O(1) per batch for primitives)
- **Pre-reserved batch storage** (avoids reallocations)
- **Hot/cold tag separation** (no data duplication)

**Usage:**
```bash
# Quick test (20K records)
python .scripts/test_optimized_writer.py input.sam.gz /scratch/tmp/output

# Full benchmark (100K records)
python .scripts/test_optimized_writer.py input.sam.gz /scratch/tmp/output 100000
```

**Results (100K alignments, 410 alignments/read avg):**
| Mode | Throughput | Output Size |
|------|------------|-------------|
| Single-table | 28,400 rec/s | 1,264 KB |
| Normalized | 29,900 rec/s | 1,504 KB |

**Architecture:**
```
test_optimized_writer.py (Python SAM parsing)
    └── optimized_parquet_writer.pyx (Cython)
        ├── str_to_cpp_string() (zero-copy PyUnicode_AsUTF8AndSize)
        ├── parse_tags_hot_cold() (nogil tag extraction)
        └── arrow_parquet_writer_optimized.cpp (C++)
            ├── pack_sequence_2bit() (LUT-based, branchless)
            ├── AppendValues() (bulk Arrow API)
            └── ZSTD compression (level 6)
```

### Benchmark Results (10K alignments, 196 alignments/read avg)

| Schema | Size | vs Raw SAM | vs Current |
|--------|------|------------|------------|
| Raw SAM | 2,719 KB | 100% | - |
| Current (duplicated) | 224 KB | 8.4% | baseline |
| **Optimized** | 118 KB | 4.4% | **48% smaller** |
| **Normalized** | 97 KB | 3.6% | **57% smaller** |

### Key Optimizations

1. **No data duplication** - Hot tags as columns OR in raw string, never both
2. **2-bit sequence packing** - 4x smaller than ASCII (A=00, C=01, G=10, T=11)
3. **Normalized tables** - Separate reads (heavy) from alignments (light)
4. **Computed-on-read** - Don't store `end_position`, `alignment_length`
5. **Dictionary encoding** - `ref_id` instead of `ref_name` strings
6. **ZSTD level 6** - Better compression than default level 3

### Recommended Schema

**Normalized two-table design:**
- `reads.parquet`: read_id, read_name, sequence (2-bit), quality
- `alignments.parquet`: read_id, ref_id, position, mapq, flag, cigar, AS, NM, MD, tags_cold
- `references.parquet`: ref_id, ref_name (sidecar)

This design:
- Deduplicates sequence/quality for multi-mapped reads (huge savings)
- Supports fast filtering via hot tag columns
- Enables lossless SAM reconstruction via tags_cold
- Scales well with highly multi-mapped metagenomics data

### test_zs_calculation.py
Tests ZS score (log-likelihood alignment score) calculation during optimized Parquet conversion.

**Usage:**
```bash
python .scripts/test_zs_calculation.py
```

### test_derived_columns.py
Tests the `aligned_length` and `is_reverse` derived columns in optimized Parquet output.
These columns are precomputed during conversion for fast queries:
- `aligned_length`: reference span from CIGAR (uint16, compresses 99% better than end_position)
- `is_reverse`: flag & 0x10 (reverse strand indicator)

**Usage:**
```bash
python .scripts/test_derived_columns.py
```

### test_aligned_length.py
Tests the `aligned_length` column (uint16) vs `end_position` (int32) for compression efficiency.

**Result**: aligned_length compresses 99% better (15.7 MB → 0.15 MB) because values cluster around read lengths (30-150bp) rather than spanning the full genome coordinate space.

### test_sorted_parquet.py
Tests position sorting using PyArrow for compression improvement.

**Result**: Position sorting with PyArrow made files LARGER due to different encoder settings than the C++ writer.

### sort_parquet_duckdb.py
Sort Parquet file by (ref_id, position) using DuckDB's efficient external merge sort.

```bash
python sort_parquet_duckdb.py input.parquet output.parquet --compare --verify
```

**Result**: Position sorting increases file size for metagenomic data because it destroys the natural sequencer-based clustering that benefits columns like cigar, AS, NM, aligned_length.

## Compression Optimization Findings

### 1. aligned_length vs end_position
- **Winner**: aligned_length (uint16)
- **Reason**: Values cluster around read lengths (most reads 30-150bp), enabling excellent dictionary + RLE compression
- **Implementation**: Replaced end_position with aligned_length in the schema
- **Savings**: ~15 MB per 3.6M records

### 2. Position Sorting
- **Tested**: Sort by (ref_id, position) before writing Parquet
- **Expected**: Better compression due to RLE on ref_id and delta on position
- **Actual Result**: File size INCREASED by 11-17%
  - ref_id improved (+77%)
  - position improved (+35%)
  - But cigar, AS, NM, aligned_length, MD all got WORSE (-200% to -1200%)
- **Reason**: Metagenomic data has natural clustering by read properties (from sequencer run order) that is more valuable for compression than position clustering
- **Recommendation**: Do NOT sort metagenomic data by position; keep natural order

### 3. MD Tag Tokenization
- **Tested**: Parse MD string into (match_lengths[], mismatch_bases[]) arrays
- **Expected**: Better compression due to structured data
- **Actual Result**: Tokenized format compresses WORSE (3.83 MB → 4.64 MB)
- **Reason**: ZSTD already achieves excellent compression (11%) on repetitive MD strings like "34M", "0A0C"
- **Recommendation**: Keep MD as string; tokenization adds overhead without benefit

### 4. Per-Column Compression Results (3.6M records)

| Column | Compressed Size | Compression Ratio | Notes |
|--------|-----------------|-------------------|-------|
| ref_id | 10.55 MB | - | Dictionary encoding |
| position | 15.73 MB | - | High entropy, delta helps |
| aligned_length | 0.15 MB | 99% smaller than end_position | Clusters at read lengths |
| is_reverse | 0.40 MB | - | Boolean, excellent |
| MD | 5.33 MB | 11.5% | String, ZSTD excellent |
| quality | 1.62 MB | - | Raw bytes |
| sequence_packed | 1.70 MB | - | 2-bit packed |
| cigar | 0.38 MB | - | Dictionary encoding |

### 5. Data Characteristics (Metagenomic Alignments)

- Reads map to many different references (scattered ref_ids)
- Natural clustering by read properties from sequencer run order
- aligned_length clusters heavily (e.g., 963/1000 rows at length 47)
- cigar patterns cluster (29 unique in 1000 mid-file rows)
- MD strings are repetitive (10.8% compression ratio)

**Conclusion**: Sort by read_id for optimal compression AND downstream task performance.

### 6. Optimal Sort Order: read_id

**Discovery**: Sorting by read_id provides both compression benefits AND optimal data layout for downstream tasks.

| Sort Order | File Size | Change | Notes |
|------------|-----------|--------|-------|
| Unsorted | 49.7 MB | baseline | Natural sequencer order |
| **By read_id** | **31.4 MB** | **-37%** | Multi-mapped reads cluster |
| By position | 54.5 MB | +10% | Destroys read clustering |
| By ref_id, read_id | 49.0 MB | -1% | Minimal benefit |

**Why read_id sorting works:**
- Multi-mapped reads cluster together → same sequence/quality repeats
- Same read → similar alignment properties (cigar, AS, NM)
- Sequential read_ids compress with RLE
- Quality strings compress 79% better when clustered

**Per-column improvement with read_id sort:**
- quality: +79% (identical strings for multi-mapped reads)
- read_name: +69% (PE read prefixes cluster)
- read_id: +63% (sequential RLE)
- sequence_packed: +51% (identical for multi-mapped)
- position: +33% (some locality preserved)
- ref_id: +27% (multi-mapped share references)
- MD: +22% (similar alignment patterns)

**Downstream task benefits:**
- LCA algorithm: reads already grouped by read_id ✓
- EM algorithm: reads already grouped by read_id ✓
- Filtering: efficient row group pruning ✓

### convert_and_sort.py

Integrated converter that produces optimally sorted Parquet files.

```bash
python convert_and_sort.py input.sam.gz output.parquet -t 4
```

**Two-phase approach:**
1. Convert BAM/SAM → Parquet (streaming, C++ Arrow writer)
2. Sort by read_id with DuckDB (external merge sort)

**Performance (3.6M records):**
- Convert: 36s (100K rec/s)
- Sort: 3.5s (1M rec/s)
- Total: 39.5s
- Output: 31.4 MB (37% smaller than unsorted)

### test_cross_domain_flags.sh
Tests the cross-domain removal CLI flags.

### test_hierarchical_em.sh
Tests the hierarchical EM algorithm with Dirichlet-tree prior for taxonomy-informed
read assignment. This script tests various configurations:
- Standard EM (baseline)
- Hierarchical EM with default parameters
- Hierarchical EM with custom alpha parameters
- Hierarchical EM with misannotation detection
- Hierarchical EM with custom unknown bucket ranks

Usage:
```bash
BAM_FILE=/path/to/test.bam TAXONOMY_DB=/path/to/tax OUTPUT_DIR=./test_output ./test_hierarchical_em.sh
```

## Cross-Domain Removal (Default Behavior)

**NEW DEFAULT**: When using `--taxonomy-filter` with `--taxonomy-db`, combined cross-domain removal
is now enabled automatically. This applies the most effective decontamination strategy:

1. **Stage 5b**: Removes entire references flagged as cross-domain contamination
2. **Stage 6**: Removes remaining cross-domain alignments between different domains

### Flag Reference

| Flag | Description |
|------|-------------|
| `--taxonomy-filter --taxonomy-db PATH` | **Enables combined mode by default** |
| `--no-cross-domain-removal` | Disable automatic cross-domain removal |
| `--remove-cross-domain-alignments` | Remove only cross-domain alignments (edge-level) |
| `--remove-cross-domain-references` | Remove only cross-domain references (node-level) |
| `--remove-cross-domain-all` | Explicit combined mode (same as default) |
| `--detect-misannotations` | Flag potential database misannotations (enabled by default) |

### Examples

```bash
# Default combined mode (recommended for ancient DNA)
filterBAM reassign --bam input.bam -o output.bam --taxonomy-filter --taxonomy-db /path/to/tax

# Disable cross-domain removal (taxonomy info still added to output)
filterBAM reassign --bam input.bam -o output.bam --taxonomy-filter --taxonomy-db /path/to/tax \
    --no-cross-domain-removal

# Only alignment-level removal (preserves all references)
filterBAM reassign --bam input.bam -o output.bam --taxonomy-filter --taxonomy-db /path/to/tax \
    --no-cross-domain-removal --remove-cross-domain-alignments
```
