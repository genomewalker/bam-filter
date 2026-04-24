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

## Tempered EM Testing

### test_tempered_em.sh
Tests the tempered EM parameter (`--em-beta`) for bias reduction in read assignment.

**Background**:
Standard EM algorithm suffers from "rich-gets-richer" bias where abundant references
accumulate even more reads unfairly. Tempered EM scales the log-likelihood by β < 1,
flattening posteriors and reducing this bias.

```
p(j|r) ∝ π_j × exp(β × score)
```

- β = 1.0: Standard EM behavior
- β = 0.3-0.7: Recommended for bias reduction
- Lower β = more uniform distribution, less winner-take-all

**Usage:**
```bash
# Run comparison test (β=1.0 vs β=0.5)
./test_tempered_em.sh /path/to/input.bam

# Default uses Viking test BAM
./test_tempered_em.sh
```

**Expected outcomes:**
- Standard EM (β=1.0): More concentrated assignments, potentially biased
- Tempered EM (β=0.5): More distributed assignments, reduced bias toward abundant refs

## Parallel SAM Conversion Analysis

### parallel_sam_analysis_report.md
Comprehensive architectural analysis of the current parallel SAM to Parquet conversion implementation:

- **Current Performance**: 686K records/sec with 4 BGZF workers
- **Architecture Overview**: BGZF block-based parallel processing with ProcessPoolExecutor
- **Bottleneck Analysis**: File I/O contention, HTSlib overhead, single-threaded workers
- **Scaling Recommendations**: Hybrid process/thread model, NUMA awareness, memory optimizations
- **Expected Improvements**: Up to 5.5M+ records/sec on 32 cores

Key findings:
1. Current architecture uses smart BGZF block boundaries for true parallelism
2. HTSlib parsing overhead can be reduced 40% with custom Cython parser
3. Storage I/O becomes bottleneck at 16+ workers
4. Memory hierarchy optimizations needed for 32+ cores

The analysis provides a roadmap for scaling from current 686K rec/s to 5.5M+ rec/s performance.


## EM Algorithm Quality Evaluation

### compare_em_quality.sh
Runs a side-by-side comparison of standard EM vs unified φ-space EM implementation.

**Features:**
- Runs both EM variants on the same BAM file
- Generates reference statistics for each
- Outputs comparison summary

**Usage:**
```bash
# Edit BAM and TAXDB paths in script first, then:
bash .scripts/compare_em_quality.sh
```

### analyze_em_quality.py
Analyzes and compares EM quality metrics from reference statistics files.

**Metrics calculated:**
- **Entropy**: Shannon entropy of weight distribution (higher = more diverse)
- **Gini coefficient**: Inequality measure (0=equal, 1=concentrated)
- **Max weight**: Largest single reference weight
- **Top-5 share**: Percentage of reads in top 5 references

**Usage:**
```bash
# After running compare_em_quality.sh:
python .scripts/analyze_em_quality.py
```

**Expected φ-space EM improvements:**
- Higher entropy (reduced rich-get-richer effect)
- Lower Gini coefficient (more equal distribution)
- Lower max weight (reduced single-reference dominance)
- Lower top-5 share (less concentration in few references)

### test_network_qc.py
Tests the network QC metrics module using the EM comparison data.

**Metrics computed:**
- **N_eff** (effective number): Diversity measure based on inverse Simpson index
- **Gini coefficient**: Inequality measure for read distribution
- **Entropy**: Shannon entropy of read distribution
- **Top-k shares**: Concentration in top references

**Usage:**
```bash
# First run the EM comparison
bash .scripts/compare_em_quality.sh

# Then analyze with network QC metrics
python .scripts/test_network_qc.py
```

**Expected output:**
The φ-space EM should show improvements across all metrics:
- Higher N_eff (more diverse)
- Lower Gini (more equal)
- Higher entropy (more diverse)
- Lower max reads (less dominance)
- Lower top-k shares (less concentration)

### test_network_qc_filtering.py
Tests the network QC-based chimera/contamination filtering functionality.

**Features:**
- Tests `--network-qc-filter` option for chimera detection
- Uses participation coefficient and neighbor taxonomy entropy
- Compares results with and without network QC filtering

**Usage:**
```bash
python .scripts/test_network_qc_filtering.py /path/to/input.bam /path/to/taxonomy_db
```

**CLI Options Added:**
- `--network-qc-filter`: Enable network QC-based chimera detection
- `--chimera-threshold FLOAT`: Combined chimera score threshold (default: 0.5)
- `--participation-threshold FLOAT`: Participation coefficient threshold (default: 0.5)
- `--neighbor-entropy-threshold FLOAT`: Neighbor taxonomy entropy threshold (default: 0.7)
- `--chimera-removal-level INT`: Removal level 0-3 (default: 2=likely)

**Chimera Detection Logic:**
References are flagged as chimeras based on:
1. **Participation coefficient (P)**: High P indicates connections spread across multiple communities
2. **Neighbor taxonomy entropy (H)**: High H indicates taxonomically diverse neighbors
3. **Combined score**: Weighted combination of P and H

**Chimera Flag Levels:**
- 0 = Clean (no chimera signal)
- 1 = Suspect (weak signal)
- 2 = Likely (moderate signal, default removal threshold)
- 3 = Confident (strong signal)

## Simplified CLI Testing

### test_simplified_cli.py
Tests the simplified CLI option mapping logic in reassign.py.

**Purpose:**
Validates that the new user-friendly CLI options correctly map to internal filterBAM parameters.

**Options Tested:**

| Option | Values | Description |
|--------|--------|-------------|
| `--filter-mode` | none, structural, taxonomy, strict | Filtering strategy |
| `--cross-domain-mode` | none, alignments, references, all | Cross-domain removal scope |
| `--sensitivity` | low, medium, high | Threshold aggressiveness |
| `--em-preset` | fast, balanced, thorough | EM algorithm configuration |
| `--em-unknown-component` | on/off | Unknown component for unmatched reads |
| `--hierarchical-pmd` | on/off | Hierarchical ancient/modern classification |

**Usage:**
```bash
python .scripts/test_simplified_cli.py
```

**Test Coverage:**
1. **Filter Mode Mapping**: Verifies each mode sets correct internal flags:
   - `none`: Pure EM, no graph filtering
   - `structural`: Graph topology only, no taxonomy
   - `taxonomy`: Graph + taxonomy with cross-domain removal
   - `strict`: Aggressive filtering with all features enabled

2. **Sensitivity Mapping**: Verifies threshold presets:
   - `low`: Conservative (0.15, 0.35, 0.60 thresholds)
   - `medium`: Balanced (0.10, 0.25, 0.50 thresholds)
   - `high`: Aggressive (0.05, 0.15, 0.40 thresholds)

3. **Cross-Domain Mode Mapping**: Verifies removal scope:
   - `none`: No cross-domain removal
   - `alignments`: Only remove cross-domain alignments
   - `references`: Only remove cross-domain references
   - `all`: Remove both alignments and references

4. **EM Preset Mapping**: Verifies EM algorithm presets:
   - `fast`: 25 iterations, 1e-5 tolerance (quick convergence)
   - `balanced`: 50 iterations, 1e-6 tolerance (default)
   - `thorough`: 100 iterations, 1e-8 tolerance (complex samples)

5. **EM Defaults**: Verifies ancient DNA optimized defaults:
   - `--hierarchical-pmd`: True (ancient/modern classification)
   - `--em-unknown-component`: True (absorbs unmapped reads)

6. **Hidden Options**: Verifies detailed EM options are hidden from basic help:
   - `--max-em-iterations`, `--em-tolerance`, `--em-beta`, etc.
   - All SQUAREM acceleration options

**Two-Tier Help System:**
- `--help`: Shows simplified options only (for most users)
- `--advanced-help`: Shows ALL options including hidden granular settings (for experts)

**Ancient DNA Optimized Defaults:**
The tool defaults are optimized for ancient DNA analysis:
- `--hierarchical-pmd=True`: Automatically classify references as ancient/modern based on PMD damage
- `--em-unknown-component=True`: Include unknown bucket to absorb reads that don't match well

Use `--no-hierarchical-pmd` and `--no-em-unknown-component` to disable these features.

**Example:**
```bash
# Basic help (simplified options)
filterBAM reassign --help

# Advanced help (all options including expert settings)
filterBAM reassign --advanced-help

# Ancient DNA processing (default settings)
filterBAM reassign --bam input.bam -o output.bam

# Faster processing with more iterations
filterBAM reassign --bam input.bam -o output.bam --em-preset fast

# Modern DNA (disable ancient-specific features)
filterBAM reassign --bam input.bam -o output.bam --no-hierarchical-pmd --no-em-unknown-component
```

## Unified EM Framework

### UNIFIED_EM_REFACTORING_PLAN.md

Comprehensive design document for the mathematically correct unified EM implementation.

**Key Design Principles (GPT-5.2-pro validated):**

1. **φ-space optimization**: EM optimizes in φ-space (true mixture weights)
2. **Power transform is POST-PROCESSING ONLY**: `π_j = φ_j^ρ / Σ φ_k^ρ` applied only at output
3. **γ_j in E-step**: Per-reference ancient fraction must influence E-step
4. **ω_ig is FIXED**: PMD-based per-read priors are fixed throughout EM
5. **Unknown in same normalization**: Unknown component competes with known refs
6. **SQUAREM + safeguard = valid GEM**: Monotonic likelihood guaranteed

**Corrected Formulas:**

```
E-step (hierarchical):
  r_{igj} ∝ φ_j × [γ_j × ω_{i,anc} × L_anc + (1-γ_j) × ω_{i,mod} × L_mod]

M-step (standard EM):
  φ_j = (Σ_i r_{ij} + α) / (Σ_k Σ_i r_{ik} + K×α)
  γ_j = (S_anc_j + α_γ) / (S_anc_j + S_mod_j + 2α_γ)

Output transform (POST-PROCESSING):
  π_j = φ_j^ρ / Σ_k φ_k^ρ  (ρ < 1 FLATTENS distribution)
```

**Power Transform Direction:**
- `ρ < 1` (e.g., 0.7): FLATTENS distribution (counteracts rich-get-richer)
- `ρ = 1`: No change (standard EM)
- `ρ > 1`: SHARPENS distribution (amplifies differences)

### processor_em_unified.pyx

New unified EM implementation in `bam_filter/processor_em_unified.pyx`.

**Key Features:**
- Single coherent EMState struct with all parameters
- Unified e_step_unified() computing all responsibilities in one pass
- Standard m_step_unified() with Dirichlet prior (no power transform)
- SQUAREM acceleration with monotonicity safeguard (valid GEM)
- Post-processing transform_to_output() for power rho

**Usage:**
```bash
# Enable unified EM with --em-unified flag
filterBAM reassign --bam input.bam -o output.bam --em-unified

# With power transform (ρ=0.7 flattens distribution)
filterBAM reassign --bam input.bam -o output.bam --em-unified --em-power-rho 0.7
```

### test_unified_em.py

Test script comparing legacy EM with unified EM implementation.

**Features:**
- Runs both EM variants on the same BAM file
- Extracts EM metrics (iterations, convergence, SQUAREM stats)
- Compares output alignment counts
- Mathematical correctness checks

**Usage:**
```bash
# Basic comparison test
python .scripts/test_unified_em.py /path/to/input.bam

# With taxonomy database
python .scripts/test_unified_em.py /path/to/input.bam --taxonomy-db /path/to/tax

# Keep output files for inspection
python .scripts/test_unified_em.py /path/to/input.bam --keep-output --output-dir ./em_comparison
```

**Expected Outcomes:**
- Unified EM should converge (CONVERGED message)
- SQUAREM acceleration should be active
- Output should contain alignments
- Power transform applied only at output (not during iterations)

### validate_unified_em_extended.py

Extended validation script implementing GPT-5.2-pro's three sanity checks.

**Checks Performed:**
1. **Log-likelihood monotonicity**: Verifies EM is doing proper coordinate ascent
2. **Lost read characteristics**: Analyzes reads that lost alignments
3. **Entropy comparison**: Compares posterior entropy between implementations

**Usage:**
```bash
# Quick mode (log-likelihood only)
python .scripts/validate_unified_em_extended.py /path/to/input.bam --skip-read-analysis

# Full mode with read analysis
python .scripts/validate_unified_em_extended.py /path/to/input.bam --keep-output
```

### analyze_em_with_tags.py

Efficient EM comparison using filterBAM-specific tags (ZP, ZS, PM, AN, DA).

**Tags Used:**
- **ZP**: Posterior probability (P(ref|read))
- **ZS**: Alignment score (log-likelihood)
- **PM**: PMD score
- **AN**: Raw ANI (%)
- **DA**: Damaged/corrected ANI (%)

**Features:**
- Samples alignments efficiently with `shuf`
- Computes per-read posterior entropy
- Analyzes characteristics of "lost" reads
- Compares entropy changes for common reads

**Usage:**
```bash
python .scripts/analyze_em_with_tags.py legacy.bam unified.bam --sample-size 50000
```

**Validation Results (MED-2022-10 dataset, 17M alignments):**

| Metric | Legacy EM | Unified EM |
|--------|-----------|------------|
| Iterations | 8 | 3 |
| Mean ZP | 0.093 | 0.098 |
| Mean entropy | 0.137 | 0.122 |

**Key Findings (GPT-5.2-pro validated):**

1. **Log-likelihood monotonicity**: PASS
   - Iter 0→3: -17.29 → -17.23 (strictly increasing)
   - Smooth decay of changes (+0.056, +0.006, +0.002)

2. **Entropy comparison**:
   - Unified EM has LOWER mean entropy (0.122 vs 0.137)
   - For common reads: 27.6% lower, 26.7% higher, 45.7% similar
   - **Interpretation**: Unified EM provides sharper assignments for genuinely informative reads

3. **Lost reads analysis**:
   - "Lost" reads (in legacy, not unified) have LOWER entropy (0.061 vs 0.137)
   - **Interpretation**: Legacy EM was OVER-CONFIDENT on these reads
   - The in-loop power transform artificially inflated confidence on marginal reads
   - Unified EM correctly refuses to assign these artificially confident reads

**GPT-5.2-pro Validation Quote:**
> "The pattern you see is textbook over-confidence: legacy pushes some moderately
> informative reads into artificially low-entropy, high-posterior categories that
> unified EM refuses to 'believe.'"

## EM Algorithm Technical Documentation

### EM_ACCELERATION_MECHANISMS.md

Comprehensive technical document describing the EM acceleration and convergence mechanisms:

**Contents:**
1. **SQUAREM Acceleration** (Varadhan & Roland 2008)
   - Quadratic extrapolation formula: θ* = θ₀ - 2α·r + α²·v
   - Three steplength schemes (S1, S2, S3)
   - Globalization with backtracking
   - Geometric interpretation as quasi-Newton step

2. **MAD-Based Convergence Detection**
   - Median Absolute Deviation (MAD) scaled by 1.4826
   - 5-element rolling history buffers
   - Adaptive tolerance: tol_eff = max(base_tol, 3×σ_MAD)
   - Limit-cycle detection for power-reweighted EM

3. **Hierarchical EM for Ancient/Modern Classification**
   - Per-reference γ_k = P(ancient | reference k)
   - PMD damage patterns (C→T at 5', G→A at 3')
   - Beta-Bernoulli conjugate update: γ_k = (S_anc + α)/(S_anc + S_mod + 2α)

4. **Full Algorithm Integration**
   - One-iteration walkthrough diagram
   - Control flow flowchart
   - Numerical example with actual calculations
   - Complete pseudocode (plain EM vs accelerated EM)

**Key Equations:**
```
SQUAREM extrapolation:
  r = θ₁ - θ₀ (first-order difference)
  v = θ₂ - θ₁ - r (second-order difference)
  α = -||r||/||v|| (steplength)
  θ* = θ₀ - 2α·r + α²·v

MAD-based tolerance:
  σ_MAD = 1.4826 × median(|x_i - median(x)|)
  tol_eff = max(base_tolerance, 3 × σ_MAD)

Hierarchical γ update:
  γ_k = (S_anc_k + α) / (S_anc_k + S_mod_k + 2α)
```

**Usage:**
```bash
# View the full document
cat .scripts/EM_ACCELERATION_MECHANISMS.md
```

### benchmark_em_optimization.py

Benchmarks the optimized E-step implementation against baseline timing.

**Optimizations Implemented:**
1. **Fast path for single-alignment reads**: Skip log/exp when alignment_count==1
2. **Per-thread accumulators with prange**: Parallel E-step with thread-local phi_counts
3. **Scratch buffer caching**: Store log weights in first pass, reuse in second
4. **Max-subtraction log-sum-exp**: Numerically stable normalization
5. **Software prefetching**: __builtin_prefetch for next alignment

**Usage:**
```bash
# Basic benchmark
python .scripts/benchmark_em_optimization.py /path/to/input.bam

# With reference BAM for numerical comparison
python .scripts/benchmark_em_optimization.py /path/to/input.bam --reference-bam /path/to/baseline.bam

# Custom baseline timing (from previous run)
python .scripts/benchmark_em_optimization.py /path/to/input.bam --baseline-em-time 334.28
```

**Metrics Reported:**
- Total execution time
- EM algorithm time
- EM iterations to convergence
- Alignments kept after probability filtering
- Speedup vs baseline (x factor)
- Numerical equivalence check (MD5 comparison)

## Hierarchical EM Gamma Analysis

### analyze_gamma_metadmg_correlation.py

Analyzes correlation between hierarchical EM gamma values and MetaDMG A_b damage estimates.

**Key Findings:**

1. **No correlation for high-read-count refs**: The gamma values do NOT correlate with MetaDMG A_b
   for references with many multi-mapping reads. This is because:
   - Multi-mapping reads have diffuse posteriors spread across many references
   - Each reference gets only a fraction of each read's contribution
   - The effective sample size (S_tot) per reference is small
   - The prior (alpha=1.0) dominates, pushing gamma → 0.5

2. **Good discrimination for unique-read refs**: Gamma values ARE discriminative for refs with unique reads:
   - High gamma (>0.7): ~1,100 refs with avg ~10 unique reads that show damage
   - Low gamma (<0.3): ~350 refs with avg ~12 unique reads that don't show damage

3. **MetaDMG refs have many multi-mapping reads**: The 118 MetaDMG refs all have enough reads for
   damage fitting, but these reads are often multi-mapped (avg 40-50% multimap rate), diluting
   the damage signal in the gamma estimate.

**Usage:**
```bash
python .scripts/analyze_gamma_metadmg_correlation.py \
    /path/to/metadmg_results.dfit.txt.gz \
    /path/to/reassign_stats.tsv.gz
```

**Output:**
- Pearson correlation between gamma and A_b
- Stats for high/low damage refs
- Stats for high/low gamma refs
- Multimap rate analysis

**Recommendations:**
1. For ancient/modern classification of multi-mapping references, use per-alignment damage metrics
   (ZP tags with PM/DA scores) rather than per-reference gamma
2. The gamma estimate is most reliable for references with many unique reads
3. Consider reducing the prior strength (alpha < 1.0) for samples with many multi-mapping reads

**Related Files:**
- `bam_filter/processor_em.pyx`: Hierarchical EM implementation
- `bam_filter/processor_pmd.pyx`: PMD damage calculation
- `EM_ACCELERATION_MECHANISMS.md`: Full EM algorithm documentation

## Bayesian Damage Model Validation

### validate_damage_model.sh

Validates the Bayesian damage model outputs against MetaDMG results.

**Damage Model Output Columns:**

| Column | Description |
|--------|-------------|
| `gamma_ancient` | P(ancient \| data) - Bayesian posterior probability [0-1] |
| `damage_amplitude` | A - damage rate at position 1 (like MetaDMG A) |
| `damage_baseline` | b - baseline divergence (constant component) |
| `damage_log_bf` | log[P(data\|ancient)/P(data\|modern)] - fit quality |

**Usage:**
```bash
bash .scripts/validate_damage_model.sh
```

**Key Findings:**

1. **Amplitude correlation with MetaDMG**: Pearson r=0.77 (strong correlation)
   - Our values are ~1.4x higher due to EM-weighted accumulation

2. **Log Bayes Factor distribution**:
   - 25% negative (modern hypothesis fits better)
   - 51% low (0-5, weak evidence for ancient)
   - 20% mid (5-20, moderate evidence)
   - 4% high (>=20, strong evidence for ancient)

3. **Detecting weird damage patterns**:
   - Low log_bf + high amplitude = poor fit to exponential decay model
   - These may indicate contamination, chimeric sequences, or unusual damage

**Interpretation:**
- **High gamma_ancient (>0.8) + high log_bf (>20)**: Confidently ancient
- **Low gamma_ancient (<0.4) + low log_bf (<5)**: Likely modern or inconclusive
- **High amplitude + low log_bf**: Weird damage pattern, investigate manually

### compare_damage_models.py

Validates filterBAM's hierarchical EM damage parameters against metaDMG's dfit output.

**Usage:**
```bash
python .scripts/compare_damage_models.py \
  /path/to/metadmg_results.dfit.txt.gz \
  /path/to/filterbam_stats.tsv
```

**Comparison Metrics:**
- Damage amplitude correlation (filterBAM vs metaDMG A parameter)
- Ancient classification concordance (gamma > 0.5 vs Zfit > 2)
- Confusion matrix and sensitivity/precision
- Disagreement analysis

**Validation Results (MED-2022-10 dataset):**

| Metric | Value |
|--------|-------|
| References matched | 118 |
| Amplitude correlation (Pearson r) | 0.73 |
| Classification concordance | 90.7% |
| Sensitivity | 100% |
| Precision | 90.7% |

**Key Findings:**
1. filterBAM's damage_amplitude correlates well with metaDMG's A parameter (r=0.73)
2. All 107 ancient references identified by metaDMG (Zfit > 2) are also classified as ancient by filterBAM (gamma > 0.5)
3. 11 borderline references (Zfit 1-2) are classified as ancient by filterBAM but modern by metaDMG - these are marginal cases
4. Full dataset classification: 42% ancient, 38% borderline, 19% modern

## Code Review Documentation

### REASSIGN_WORKFLOW_REVIEW.md

Comprehensive code review of the reassign pipeline including:

- **Bugs found**: Legacy graph builder code, gamma_values semantic overwrite
- **Performance analysis**: Bottlenecks identified and mitigations
- **Numerical stability**: Log-space arithmetic verified throughout
- **Concurrency**: Thread safety verified for parallel sections
- **Recommendations**: Immediate, short-term, and long-term fixes

**Key Findings:**
- Core EM/PMD/Bayesian pipeline is numerically sound
- 10x speedup achieved with precomputed log-likelihoods
- Legacy `build_weighted_graph_from_alignments` should be removed (dead code)
- Hash table sizes could be dynamic instead of fixed 1.5GB

**Usage:**
```bash
cat .scripts/REASSIGN_WORKFLOW_REVIEW.md
```

## Damage Correction Pipeline Testing

### test_damage_pipeline.sh

Full pipeline test for damage-corrected filtering: reassign → filter → lca → prob_profile.

**Usage:**
```bash
bash .scripts/test_damage_pipeline.sh
```

**Configuration (edit script):**
- `INPUT_BAM`: Input BAM file
- `TAXONOMY`: Taxonomy database directory
- `OUTDIR`: Output directory
- `PREFIX`: Output file prefix

**Pipeline Steps:**
1. **Reassign**: EM-based read reassignment with damage correction
2. **Index**: samtools index of reassigned BAM
3. **Filter**: Reference filtering with damage-corrected ANI (`--damage-correction`)
4. **Index**: samtools index of filtered BAM
5. **LCA**: Lowest common ancestor analysis
6. **Prob Profile**: Probabilistic profile generation

### verify_intervals.py

Manual verification of interval counting against BAM data.

**Note:** When comparing with `samtools depth`, remember that samtools depth **excludes secondary alignments by default** (flag 0x100). For reassigned BAMs where most alignments are secondary, use:
```bash
# Count all positions including secondary alignments
samtools view -F 0x4 input.bam "ref_name" | awk '...'
```

The stats calculation correctly includes ALL alignments (essential for reassigned BAMs).

## Probabilistic Profiler Optimization Testing

### test_profiler_optimizations.py

Tests the optimized Cython probabilistic profiler (`probabilistic_profiler.pyx`).

**Optimizations Tested:**
1. **OpenMP parallelization**: Batch posterior computation with prange
2. **Precomputed Beta params**: Avoids repeated lgamma calls in damage model
3. **Log-space belief propagation**: Prevents underflow with many refs/deep trees

**Usage:**
```bash
python .scripts/test_profiler_optimizations.py
```

**Test Coverage:**
1. **Single reference computation**: Validates damage model and presence model
2. **Batch posteriors**: Measures speedup across thread counts (1, 2, 4, 8)
3. **Belief propagation**: Tests log-space numerics on taxonomy tree
4. **Underflow test**: Verifies no NaN/Inf with 200 refs per taxon

**Results (10K refs, Dec 2025):**

| Threads | Time | Speedup |
|---------|------|---------|
| 1 | 115ms | 1.00x |
| 2 | 63ms | 1.83x |
| 4 | 61ms | 1.88x |
| 8 | 55ms | 2.09x |

**Log-space validation:**
- 50 taxa, 10K refs (200/taxon): No NaN, No Inf, No Zero
- Underflow prevention working correctly

### benchmark_csr_bp.py

Benchmarks the CSR (Compressed Sparse Row) optimized belief propagation implementation.

**Optimizations Implemented:**
1. **CSR children graph**: Replaces Python dict of lists with contiguous arrays
2. **CSR refs per taxon**: Replaces Python dict of lists with contiguous arrays
3. **C array traversal**: Replaces Python list for postorder/preorder with typed arrays
4. **Nogil hot loops**: All BP passes run in pure C without GIL

**Usage:**
```bash
python .scripts/benchmark_csr_bp.py
```

**Results (Dec 2025):**

| Taxa | Refs | Time (ms) | Refs/sec |
|------|------|-----------|----------|
| 100 | 500 | 0.21 | 2.3M |
| 1,000 | 5,000 | 2.3 | 2.1M |
| 10,000 | 50,000 | 27 | 1.9M |
| 50,000 | 100,000 | 104 | 957K |
| 100,000 | 200,000 | 267 | 748K |

**Key Performance Notes:**
- Scales linearly with number of taxa and references
- No numerical issues (NaN/Inf) even at large scales
- Hot loops run in nogil for maximum throughput
- CSR format provides cache-friendly memory access

## DCMS Probabilistic Profiler

### test_dcms_from_bam.py

Tests the end-to-end DCMS profiler computing ALL statistics from BAM.

**Purpose:**
Validates the complete DCMS workflow that computes everything from scratch:
1. Coverage statistics (breadth, entropy, gini, wcb) via `compute_bam_stats`
2. Graph analysis (unique_reads, shared_reads, connected_neighbors) via `process_bam_with_em`
3. Damage parameters (gamma_ancient, damage_log_bf) via hierarchical PMD

No pre-computed TSV files required - just the BAM and taxonomy database.

**Usage:**
```bash
python .scripts/test_dcms_from_bam.py
```

**Output:**
- Coverage statistics: /tmp/LV7008866922.filt_coverage_*.tsv.gz
- Graph statistics: /tmp/LV7008866922.filt_graph_*.tsv.gz
- DCMS profile: /tmp/dcms_from_bam.tsv

**Test Results (Dec 2025):**

| Taxon | p_present | Exclusivity | TAD |
|-------|-----------|-------------|-----|
| Homo sapiens | 1.000 | 71.3% | 103,945 |
| Staphylococcus aureus | 0.084 | 7.7% | 930 |
| Plasmodium vivax | 0.502 | 40.0% | 2,399 |
| Rattus norvegicus | 0.904 | 10.1% | 22 |

**Performance:**
- Total time: ~184s (34s coverage + 148s graph analysis + 21s profiling)
- Two-pass BAM processing for graph metrics
- Temporary files cleaned up automatically

### test_exclusivity_integration.py

Tests the read exclusivity penalty integration in the DCMS (Damage-Calibrated Mixture-of-Sources) profiler.

**Purpose:**
Validates that the Beta-Binomial exclusivity penalty correctly reduces p_present for taxa with
high cross-mapping (low read exclusivity).

**Key Features:**
1. **Beta-Binomial penalty model**: Bayesian approach handles uncertainty in low-count cases
2. **Neighbor penalty**: Taxa sharing reads with many neighbors are penalized more
3. **TAD-weighted aggregation**: Exclusivity propagated up taxonomy tree weighted by abundance

**Penalty Model Parameters:**
- `eta=0.05`: Minimum penalty floor (prevents complete zeroing)
- `alpha=0.7`: Power for exclusivity curve (lower = gentler penalty)
- `n0=100`: Neighbor count for half-penalty
- `prior_strength=1.0`: Beta-Binomial prior strength

**Usage:**
```bash
python .scripts/test_exclusivity_integration.py
```

**Expected Output:**
- Penalty function test cases showing behavior across exclusivity range
- Comparison of p_present with and without exclusivity penalty
- Top taxa with largest p_present reduction
- Taxa with low exclusivity but still high p_present (flagged for review)

**Test Results (Dec 2025):**

| Taxon | p_present (no) | p_present (excl) | Exclusivity | TAD |
|-------|----------------|------------------|-------------|-----|
| Rattus norvegicus | 1.000 | 0.387 | 0.8% | 5.7 |
| Staphylococcus aureus | 0.455 | 0.091 | 9.1% | 1,026 |
| Plasmodium vivax | 0.700 | 0.530 | 55.6% | 2,937 |
| Homo sapiens | 1.000 | 1.000 | 77.7% | 110,227 |

**Interpretation:**
- Taxa with very low exclusivity (Rattus: 0.8%) receive strong penalty
- Taxa with moderate exclusivity (Plasmodium: 56%) receive mild penalty
- Taxa with high exclusivity (Homo sapiens: 78%) are not penalized
- High-TAD taxa can still have high p_present despite moderate exclusivity if signal is strong

### analyze_graph_contaminants.py

Analyzes graph metrics for known contaminant taxa by joining accession map with graph analysis.

**Usage:**
```bash
python .scripts/analyze_graph_contaminants.py
```

**Output:**
- Per-reference metrics: total reads, unique reads, shared reads, % unique, neighbors, tax entropy
- Per-taxon summary: aggregated counts across all accessions

### test_graph_contaminants.py

Tests graph analysis output columns and searches for potential contaminants by name pattern.

**Usage:**
```bash
python .scripts/test_graph_contaminants.py
```

**Note:** Uses reference name pattern matching, which may not work for accession-based datasets.
For accession-based analysis, use `analyze_graph_contaminants.py` instead.

## Coverage-Weighted Reference Priors (CWRP)

### test_em_authenticity.py

Tests and evaluates the CWRP full fix implementation including:

1. **Authenticity calculation**: How eta (latent ancientness) is computed from features
2. **Softmax redistribution**: How Dirichlet priors are redistributed based on authenticity
3. **Low-coverage shrinkage**: How low-coverage references are protected from false signals

**Authenticity Formula:**

```
eta = cov_score + damage_weight * dmg_score + 0.5 * len_score + interactions

where:
  cov_score = entropy - gini           (uniform coverage = positive)
  dmg_score = 2 * (damage_5p + damage_3p) / 2  (ancient damage = positive)
  len_score = short_frac - mean_length / (mean_length + 50)  (short = positive)

  interactions:
    + 0.3 * cov_score * dmg_score    (coverage-damage synergy)
    + 0.2 * dmg_score * len_score    (damage-length synergy)

authenticity = sigmoid(eta) = 1 / (1 + exp(-eta))
```

**Softmax Redistribution (M-step):**

```
alpha_j = alpha_0 * n_refs * softmax(lambda * authenticity_j)

where softmax_j = exp(scaled_j - max_score) / sum_k(exp(scaled_k - max_score))

Key property: sum(alpha_j) = alpha_0 * n_refs (constant total prior mass)
```

**Low-Coverage Shrinkage:**

```
eta_shrunk = eta * n_reads / (n_reads + tau)

Below low_cov_floor: shrinks toward 0 (neutral)
Above low_cov_floor: unchanged (data-driven)
```

**Usage:**
```bash
python .scripts/test_em_authenticity.py
```

**CLI Options for CWRP:**

| Option | Default | Description |
|--------|---------|-------------|
| `--cwrp-lambda` | 0.0 | Weight for coverage-weighted priors (0=disabled) |
| `--iterative-auth` | False | Enable iterative ancientness updates during EM |
| `--auth-update-interval` | 5 | Iterations between ancientness updates |
| `--damage-weight` | 1.0 | Weight for damage features in ancientness |
| `--low-cov-floor` | 10 | Shrink ancientness toward neutral below this read count |
| `--low-cov-shrink-tau` | 50.0 | Shrinkage strength for low-coverage references |

**Example Usage:**
```bash
# Standard CWRP (static authenticity)
filterBAM reassign --bam input.bam -o output.bam --cwrp-lambda 2.0

# CWRP with iterative updates (recommended for complex samples)
filterBAM reassign --bam input.bam -o output.bam --cwrp-lambda 2.0 \
    --iterative-auth --auth-update-interval 3 --damage-weight 1.5
```

## EM Algorithm Parallelization

### benchmark_em_parallelization.sh

Benchmarks the EM algorithm parallelization optimizations.

**Optimizations Implemented:**

1. **Parallelized compute_log_likelihood**: The log-likelihood computation loop was sequential.
   Changed from `for read_idx in range(...)` to `prange(...)` with thread-local accumulators
   for `total_ll` and `valid_reads`, reducing thread contention.

2. **Branchless log-sum-exp**: Replaced branchy `if log_a > log_b` comparisons with
   `fmax(log_a, log_b)` and `fabs(log_a - log_b)`. Added fast path for `|diff| > 20`
   where the smaller term is negligible (contributes < 2e-9). Uses `log1p(exp(-diff))`
   for better numerical accuracy.

**Files Modified:**
- `bam_filter/processor_em.pyx`: `compute_log_likelihood()` parallelized with prange
- `bam_filter/processor_fast_math.pyx`: `stable_log_sum_exp()` optimized

**Usage:**
```bash
bash .scripts/benchmark_em_parallelization.sh
```

**Performance Results (MED-2022-10, 2.7GB, 75.8M alignments, 4 threads):**

| Stage | Before | After | Speedup |
|-------|--------|-------|---------|
| EM Algorithm | 218.8s | 61.86s | **3.54x** |
| Total Pipeline | 356s | 188s | **1.89x** |

**Key Metrics:**
- 1.65M unique reads
- 166,461 references
- 10 EM iterations to convergence
- 4.8GB peak memory

The EM algorithm was the primary bottleneck (61% of runtime before optimization).
Parallelizing `compute_log_likelihood` eliminated the sequential bottleneck that
was called twice per SQUAREM iteration.

## Authenticity P-Value Fix (Dec 2025)

### Bug Description

The `authenticity_pvalue` column had inverted semantics that made filtering impossible:

**Before (incorrect):**
- Low score (-0.93) → Low pvalue (0.000) - detected as "authentic"
- High score (0.93) → High pvalue (0.999) - detected as "contamination"

**After (correct):**
- Low score (-0.93) → High pvalue (0.999) - not significantly authentic
- High score (0.93) → Low pvalue (0.006) - significantly authentic

### Fix Applied

Changed pvalue calculation from `CDF(z)` to `1 - CDF(z)` in two files:

1. **bam_filter/stats.pyx** (line 1750):
```cython
# Before: ref_stats[i].authenticity_pvalue = _std_normal_cdf(z_score)
# After:
ref_stats[i].authenticity_pvalue = 1.0 - _std_normal_cdf(z_score)
```

2. **bam_filter/processor_lca_stats.pyx** (line 2578):
```python
# Before: pvalue = scipy_stats.norm.cdf(score, loc=mu, scale=std)
# After:
pvalue = 1.0 - scipy_stats.norm.cdf(score, loc=mu, scale=std)
```

### Interpretation

After the fix:
- **authenticity_pvalue <= 0.05**: Reference has significantly HIGH authenticity score (likely authentic)
- **authenticity_pvalue > 0.05**: Reference score not significantly different from random (may be contamination)

### Filter Example

```bash
# Filter for authentic references (high score, low pvalue, high ANI)
filterBAM filter --bam input.bam --bam-filtered output.bam \
    --stats stats.tsv --stats-filtered filtered.tsv \
    --filter 'authenticity_score:0.5:,authenticity_pvalue::0.05,read_ani_corrected_mean:95:'
```

### Test Results

Before fix: 0 references passed the filter combination
After fix: 2,100 references passed (2,048,236 alignments)

## Unified Bayesian Profiler Design (Dec 2025)

### UNIFIED_BAYESIAN_PROFILER_DESIGN.md

Comprehensive mathematical specification for a rigorous Bayesian taxonomic profiler.

**Problem Addressed:**
The current profiler multiplies reference-level P(ancient) values assuming independence:
```python
beliefs[taxid]['psi_ancient'] *= p_anc  # INCORRECT
```
When references share reads, this violates independence and produces overconfident posteriors.

**Correct Model Structure:**

1. **Read-level Factor Graph**:
   - Each read has latent origin Z_i ∈ {refs} (single true source)
   - Each read has ancient indicator A_i ∈ {0, 1}
   - Each read counted EXACTLY ONCE in likelihood

2. **Graph-GMRF Prior on Ancientness**:
   ```
   P(η) ∝ exp(-τ/2 · Σ_{r~s} W_rs (η_r - η_s)²)
   ```
   Where W_rs = shared read count between refs r and s
   References sharing many reads have correlated ancientness

3. **Taxonomy Tree Prior**:
   - Abundances propagate via Dirichlet-Multinomial
   - Ancientness propagates via Beta hierarchy
   - Parent-child smoothing with concentration κ

4. **Authenticity Integration**:
   - Coverage uniformity (entropy - gini) as presence signal
   - Separate from damage-based ancientness
   - P(authentic ancient) = P(ancient) × P(present)

**Key Mathematical Guarantees:**
- Each read contributes ONE likelihood term (no double-counting)
- Graph structure encoded in prior (not likelihood multiplication)
- Posterior calibration reflects true uncertainty
- Taxonomy hierarchy respected in aggregation

**Implementation Phases:**
1. Replace direct multiplication with GMRF-smoothed η
2. Restructure to read-level factor graph
3. Proper belief propagation on taxonomy
4. Unified authenticity model

See `.scripts/UNIFIED_BAYESIAN_PROFILER_DESIGN.md` for full mathematical specification.

### test_gmrf_smoothing.py

Tests the GMRF (Gaussian Markov Random Field) smoothing implementation for Phase 1 of the unified profiler.

**Purpose:**
Validates that GMRF smoothing correctly correlates reference-level ancientness estimates based on shared reads.

**Key Features:**
1. **Edge normalization**: `w_rs = shared_reads / sqrt(n_reads_r × n_reads_s)`
2. **Coordinate descent**: Iterative optimization with convergence checking
3. **Isolated node handling**: Nodes without neighbors retain data-driven values
4. **Tau parameter**: Controls smoothing strength (higher = more smoothing)

**Tests:**
1. **Linear chain**: Verifies smoothing on 0-1-2-3-4 connected graph
2. **Isolated nodes**: Verifies isolated nodes keep original values
3. **Tau effect**: Higher tau produces more variance reduction
4. **Edge normalization**: Different read counts handled correctly

**Usage:**
```bash
python .scripts/test_gmrf_smoothing.py
```

**Expected Output:**
- All tests should PASS
- Linear chain: ~62% variance reduction with tau=1.0
- Higher tau: progressively lower variance (0.1→10.0: 0.081→0.0006)

### PROFILER_GMRF_VALIDATION_REPORT.md

Comprehensive validation report for the GMRF smoothing integration with the probabilistic profiler.

**Summary:**
Documents the successful integration of GMRF smoothing into the pipeline, with validation on Mediterranean sapropel metagenome data.

**Key Findings:**
- GMRF smoothing runs in 0.01s on 2114 connected nodes
- Pipeline processes 2.3M alignments in ~30s total
- Validates expected taxonomic detections (Emiliania coccolithophores)

**Contents:**
- Pipeline stages and timing breakdown
- GMRF integration validation
- Alignment and output statistics
- Top taxa by presence and authenticity
- Haptophyte (coccolithophore) detection as biological validation

## Per-Reference Damage Analysis (Dec 2025)

### analyze_profile_by_rank.py

Analyzes DCMS profiler results by taxonomic rank to understand p_ancient and p_present distribution.

**Purpose:**
Validates that damage-based p_ancient values vary correctly across taxa and ranks.

**Usage:**
```bash
python .scripts/analyze_profile_by_rank.py
```

**Output:**
- Distribution by taxonomic rank (p_ancient mean±std, min-max)
- Detailed threshold analysis (>0.1, >0.5, >0.9, =1.0)
- High-confidence ancient taxa (p_ancient > 0.9 AND p_present > 0.9)
- Abundance distribution by rank (TAD sums)
- Sanity checks

### analyze_reference_level.py

Deep-dive analysis of reference-level statistics to understand authenticity score distribution.

**Purpose:**
Validates per-reference coverage metrics and damage distribution.

**Usage:**
```bash
python .scripts/analyze_reference_level.py
```

**Output:**
- Reference-level metric distributions (breadth, entropy, gini, authenticity, WCB, damage)
- Authenticity score histogram
- 5' and 3' damage distribution
- Extreme cases (lowest/highest authenticity, lowest damage)
- Correlation matrix between metrics

### verify_beta_binomial_shrinkage.py

Demonstrates how the Beta-Binomial model handles low-count references.

**Purpose:**
Shows that the model provides natural shrinkage toward prior for extreme observations.

**Usage:**
```bash
python .scripts/verify_beta_binomial_shrinkage.py
```

**Key Parameters:**
- `kappa_A = 50.0` (Ancient Beta concentration)
- `kappa_M = 200.0` (Modern Beta concentration)
- `p_0 = 0.40` (Expected ancient damage at position 1)
- `p_err = 0.005` (Expected modern error rate)

**Output:**
- Log Bayes factor and P(ancient) for test cases
- Effect of sample size on P(ancient)
- Tables for 40%, 100%, and 0% observed damage rates

### check_extreme_damage_refs.py

Checks extreme damage references and their taxon-level p_ancient values.

**Usage:**
```bash
python .scripts/check_extreme_damage_refs.py
```

**Output:**
- 100% damage references and their taxon p_ancient
- 0% damage references and their taxon p_ancient
- P(ancient) summary by damage category

**Note:**
The p_ancient values shown are taxon-level (aggregated via belief propagation).

### check_ref_damage_likelihood.py

Computes reference-level log Bayes factors and p_ancient BEFORE taxon aggregation.

**Purpose:**
Validates per-reference damage model behavior independently of taxonomy.

**Usage:**
```bash
python .scripts/check_ref_damage_likelihood.py
```

**Output:**
- Lowest p_ancient references (most modern-like)
- Highest p_ancient references (most ancient-like)
- Reference-level p_ancient distribution histogram
- Correlation between observed damage and p_ancient
- P(ancient) statistics by damage category

**Key Findings (Dec 2025):**
- 26 refs (1.2%) have p_ancient < 0.1 (modern-like)
- 2150 refs (98.6%) have p_ancient > 0.99 (ancient-like)
- Correlation between observed damage and p_ancient: 0.27

### identify_modern_contamination.py

Identifies potential modern contamination based on low p_ancient values.

**Usage:**
```bash
python .scripts/identify_modern_contamination.py
```

**Output:**
- Taxa with p_ancient < 0.9 AND p_present > 0.5
- Strongest modern signal (p_ancient < 0.5 AND p_present > 0.3)
- Cross-check with references having <10% 5' damage

### test_full_dcms_pipeline.py

Tests the full DCMS profiler pipeline with capsule mode.

**Usage:**
```bash
python .scripts/test_full_dcms_pipeline.py
```

**Steps:**
1. Compute BAM stats with `return_stats='capsule'`
2. Load taxonomy database
3. Load filtered accession map
4. Run full Cython profiler pipeline
5. Display results summary

## Per-Reference Damage Bug Fix (Dec 2025)

### Bug Description

All references were showing identical damage values (~34.5% ± 0.6%) despite having different
actual damage patterns. This caused p_ancient to be nearly 1.0 for all taxa.

### Root Cause

The global `pmd_acc` (PMDStatsAccumulator) was being copied to every RefStats instead of
per-reference damage. In `bam_filter/stats.pyx`:

```cython
# INCORRECT: Used global accumulator
if pmd_acc != NULL:
    for z in range(20):
        stats.n_5p[z] = <double>(pmd_acc.n_5p_noncpg[z] + pmd_acc.n_5p_cpg[z])
```

### Fix Applied

Added per-reference PMD accumulator in `bam_filter/stats.pyx`:

1. **Declare local accumulator** (around line 763):
```cython
cdef PMDStatsAccumulator ref_pmd_acc
memset(&ref_pmd_acc, 0, sizeof(PMDStatsAccumulator))
```

2. **Use in damage collection** (around line 821):
```cython
calculate_md_quality_score_with_stats(
    b, header, &scoring_config,
    NULL,       # pmd_result
    &ani_stats,
    &ref_pmd_acc if collect_damage_stats else NULL  # Per-reference accumulator
)
```

3. **Store in RefStats** (around lines 1099-1106):
```cython
for z in range(20):
    stats.n_5p[z] = <double>(ref_pmd_acc.n_5p_noncpg[z] + ref_pmd_acc.n_5p_cpg[z])
    stats.k_5p[z] = <double>(ref_pmd_acc.k_5p_noncpg[z] + ref_pmd_acc.k_5p_cpg[z])
    stats.n_3p[z] = <double>(ref_pmd_acc.n_3p_noncpg[z] + ref_pmd_acc.n_3p_cpg[z])
    stats.k_3p[z] = <double>(ref_pmd_acc.k_3p_noncpg[z] + ref_pmd_acc.k_3p_cpg[z])
```

### Verification

**Before fix:**
- 5' Damage std = 0.6% (all references same)
- p_ancient = 1.0 for all taxa

**After fix:**
- 5' Damage std = 13.6% (proper variation)
- p_ancient distribution: 1.2% modern-like, 98.6% ancient-like
- Correlation with observed damage: 0.27

## Sample-Adaptive Damage Model (Empirical Bayes)

### Background

The original damage model used hardcoded parameters for the "modern" (non-ancient) hypothesis:
- `p_err = 0.005` (0.5% background error rate)
- `p_bg = 0.005` (0.5% background damage)
- `tau = 5.0` (decay length in bases)

This made the model too permissive - ANY reference with damage above 0.5% showing decay pattern
would be classified as ancient, even when actual background error rates are 5x higher.

### Problem

Analysis showed:
- Sample background (positions 15-20): **2.5%** vs hardcoded **0.5%**
- Result: 98.6% of refs classified as ancient with p>0.99
- Many refs with low damage (~9%) incorrectly favored as "ancient"

### Solution: Sample-Adaptive Background

The `estimate_sample_background()` function estimates p_bg and p_err from positions 15-20
where ancient damage decay is negligible (exp(-14/5) ≈ 6%). This is the Empirical Bayes approach
used by metaDMG and other ancient DNA tools.

**Key insight**: The decay length `tau` should NOT be fitted from the sample because it's
determined by physics (overhang length, deamination chemistry). Only the background rate
should be sample-adaptive.

### Implementation

```python
from bam_filter.probabilistic_profiler import (
    estimate_sample_background,
    compute_batch_posteriors,
)

# Estimate background from sample
bg = estimate_sample_background(damage_5p_n, damage_5p_k, damage_3p_n, damage_3p_k)
print(f"Estimated p_err: {100*bg['p_err']:.2f}%")  # e.g., 2.53%

# Run profiler with adaptive background (now the default)
result = compute_batch_posteriors(
    damage_5p_n, damage_5p_k, damage_3p_n, damage_3p_k,
    breadth, mean_depth, wcb, norm_entropy, norm_gini, n_reads, tad,
    adaptive_background=True,  # Default
)
```

### Parameters

| Parameter | Default | Adaptive | Description |
|-----------|---------|----------|-------------|
| `p_bg` | 0.5% | Estimated | Background damage floor |
| `p_err` | 0.5% | Estimated | Modern error rate |
| `tau` | 5.0 | **Fixed** | Decay length (literature value) |
| `p_0` | 40% | Fixed | Peak damage at position 1 |
| `kappa_A` | 50 | Fixed | Ancient Beta concentration |
| `kappa_M` | 200 | Fixed | Modern Beta concentration |

### Results

**MED-2022-10 dataset (2177 refs):**

| Model | Mean p_anc | p<0.5 | p<0.1 | p>0.99 |
|-------|------------|-------|-------|--------|
| Hardcoded (p_err=0.5%) | 0.9879 | 26 | 25 | 2148 |
| Adaptive (p_err=2.5%) | 0.9364 | 138 | 125 | 2002 |

- 112 more refs correctly reclassified as likely modern
- High-damage refs (30-36%) still correctly classified as ancient

### Test Script

```bash
python .scripts/test_adaptive_background.py
```

### Analysis Scripts

- `design_adaptive_model.py`: Original design and comparison of approaches
- `investigate_model_problem.py`: Root cause analysis showing why full parameter fitting fails
- `check_damage_model_fit.py`: Compare model parameters with actual sample data

## Unified Damage Model EM Integration (Dec 2025)

### test_unified_damage_em.py

Tests the integration of the unified damage model INTO the EM algorithm iterations.

**Purpose:**
Validates that the damage model parameters (τ decay, δ amplitudes) are jointly optimized
with reference abundances during EM iterations, rather than just as post-processing.

**Key Integration Points:**

1. **EM Initialization**: Damage counts accumulated from all alignments before EM starts
2. **Periodic Updates**: Damage model re-fitted every 5 iterations (configurable)
3. **Authenticity Feedback**: Updated authenticity scores fed into CWRP priors
4. **Final Export**: Results copied to MemoryPool after EM converges

**Log Messages:**
```
EM_UNIFIED: Accumulating damage stats from N alignments...
EM_UNIFIED: Unified damage model initialized: tau=X.XX mu_5p=X.XXXX update_interval=5
EM_UNIFIED: Damage model update at iter N: tau=X.XX mu_5p=X.XXXX
EM_UNIFIED: Damage model final: tau=X.XX mu_5p=X.XXXX converged=1
```

**Usage:**
```bash
python .scripts/test_unified_damage_em.py
```

**Configuration:**
- Edit BAM path and taxonomy DB in script
- Default: hierarchical_pmd=True enables damage model integration

**Implementation Files:**
- `bam_filter/processor_em.pyx`: EM integration code
- `bam_filter/processor_em.pxd`: Type declarations
- `bam_filter/unified_damage.pyx`: Unified damage model

### test_pmd_tau_init.py

Tests the tau (decay) parameter initialization from PMD curve.

**Purpose:**
Validates that the decay parameter τ is properly initialized from the PMD curve's
`lambda_decay` parameter (fitted from per-position library-wide counts) instead of
trying to estimate it from aggregate per-reference counts.

**Key Points:**

1. **Problem**: Aggregate counts (sum over 8 positions) make τ unidentifiable
2. **Solution**: Use PMD curve's λ which is fitted from per-position data
3. **Conversion**: τ = 1/λ (PMD uses `exp(-λz)`, unified uses `exp(-z/τ)`)
4. **Fixed Tau**: When PMD is available, τ is fixed (prior_sd=0.01)

**Code Path:**
```
pool.pmd_curve_ptr->lambda_decay  (per-position fit)
    -> config.initial_tau = 1.0 / lambda_decay
    -> set_tau_from_pmd(damage_ctx, lambda_decay, fix_tau=True)
    -> Unified model uses fixed tau
```

**Usage:**
```bash
python .scripts/test_pmd_tau_init.py
```

### test_unified_damage_model.py

Tests the unified damage model independently (without EM integration).

**Features:**
1. **Basic Workflow**: Create context → set counts → fit model → get results
2. **Edge Cases**: Low coverage, no damage, mixed scenarios
3. **Real Data**: Tests with actual parquet stats if available

**Usage:**
```bash
python .scripts/test_unified_damage_model.py
```

**Key Model Parameters:**
- `tau`: Decay length in bases (typically 5-30)
- `mu_5p`, `mu_3p`: Global damage amplitudes
- `delta_{5p,3p}`: Per-reference amplitudes
- `authenticity`: Combined score for ancient classification

### test_pmd_tau_fix.py

Verifies the fix for tau jumping to upper bound (30) during damage model EM.

**Background:**
The unified damage model was allowing tau to jump from the PMD-derived value
(~2.76) to the upper bound (30) even with a "tight" prior (sd=0.05). This
occurred because:

1. The prior penalty was unscaled relative to millions of data points
2. The data term overwhelmed the prior in the tau M-step optimization

**Fix Applied:**
1. Skip `em_m_step_tau()` when `prior_sd < 0.1` (indicates fix_tau=True)
2. Pass PMD-derived tau to final damage model in `fit_unified_from_damage_stats_py`

**Expected Output:**
- Lambda fit shows estimated value around 0.35-0.40
- All "tau=" log lines should show PMD-derived value (~2.76)
- tau should NOT jump to 30.00

**Usage:**
```bash
python .scripts/test_pmd_tau_fix.py
```

## Damage Model Robustness Improvements (Dec 2025)

### test_damage_robustness.py

Tests the comprehensive robustness improvements to the unified damage model.

**Improvements Tested:**

1. **Coverage-aware Gamma shrinkage** (`unified_damage.pyx:em_m_step_delta`):
   - Adaptive α_r = α_base / √(1 + N_r / N_scale)
   - Low-coverage refs get strong shrinkage to global mean
   - High-coverage refs let data dominate
   - Coverage-weighted global mean updates

2. **ZS score (alignment_score) weighting** (`processor_em.pyx`):
   - Sigmoid-based quality weighting: `1 / (1 + exp(-(score - threshold) / scale))`
   - Applied in both initial and EM-weighted damage accumulation
   - Uses alignment_score (ZS tag) instead of unreliable mapQ
   - Threshold calibrated to ~99.3% ANI

3. **Unified authentication score** (`unified_damage.pyx:compute_outputs`):
   - Base authenticity: Y / (Y + B + 1)
   - Coverage confidence: 1 - exp(-N / N_confident)
   - Asymmetry bonus: CT5'+GA3' vs GA5'+CT3' ratio
   - Multi-factor combination more robust than damage-only

**Usage:**
```bash
python .scripts/test_damage_robustness.py
```

**Verification Checklist:**
- Tau should stay at PMD-derived value (~2.76)
- Coverage-weighted damage accumulation visible in logs
- Quality-weighted accumulation active
- Unified auth score incorporates coverage confidence

### test_posterior_predictive.py

Tests the posterior predictive diagnostics for damage model validation.

**Diagnostics Computed:**

1. **Chi-squared goodness-of-fit**: Σ (observed - expected)² / expected
2. **Posterior predictive p-value**: P(χ² > observed | model) using Wilson-Hilferty approximation
3. **Dispersion factor**: χ²/df (should be ~1 if model fits well)
4. **Overdispersion flag**: p < 0.01 or dispersion > 3

**Purpose:**
Identifies references where the exponential decay model is misspecified:
- Contamination mixing ancient + modern DNA
- Damage patterns inconsistent with exponential decay
- Reference-specific sequencing artifacts

**Usage:**
```bash
python .scripts/test_posterior_predictive.py
```

**Interpretation:**
- `mean_disp ≈ 1.0`: Model fits well
- `mean_disp < 1.0`: Underdispersion (common with regularization)
- `mean_disp > 3.0`: Significant model misspecification
- Low % overdispersed: Most references fit exponential decay model
- High % overdispersed: Investigate outlier references

**Output Fields:**
| Field | Description |
|-------|-------------|
| `chi_squared` | Goodness-of-fit statistic |
| `pp_pvalue` | Posterior predictive p-value |
| `dispersion_factor` | χ²/df ratio |
| `is_overdispersed` | Boolean flag (p<0.01 or disp>3) |
| `df` | Degrees of freedom |

**Log Output:**
```
Diagnostics: X/Y overdispersed (disp>3 or p<0.01), mean_disp=Z.ZZ
```

## Key Log Metrics Reference

This section defines metrics that appear in filterBAM log output.

### Damage Model Metrics

| Metric | Definition |
|--------|------------|
| `D_avg` | Average expected damage rate across positions 1-8, weighted by exp(-z/τ) decay. Computed as Σ D(z)/8 where D(z) = δ × exp(-z/τ). |
| `D(1)`, `D(5)`, `D(10)` | Expected damage probability at positions 1, 5, 10 from read ends. |
| `tau` (τ) | Decay length parameter in bases. Damage decays as exp(-z/τ). Lower τ = faster decay. |
| `omega` (ω) | Peak damage amplitude at position 1 (from PMD curve). |
| `mu_5p`, `mu_3p` | Global mean damage amplitudes for 5' C→T and 3' G→A across all usable references. |
| `usable_refs` | References with ≥50 total damage opportunities (n_5p + n_3p summed over positions 1-8). Used for reliable parameter estimation. |

### EM Algorithm Metrics

| Metric | Definition |
|--------|------------|
| `S_anc` | Accumulated ancient evidence: Σ (φ-weighted posterior) × ω_{i,anc} × P(damage|ancient) |
| `S_mod` | Accumulated modern evidence: Σ (φ-weighted posterior) × ω_{i,mod} × P(damage|modern) |
| `frac_anc` | Ancient fraction = S_anc / (S_anc + S_mod). Proportion of evidence supporting ancient origin. |
| `gamma` (γ) | Per-reference ancient probability = (S_anc + α) / (S_anc + S_mod + 2α). Shrunk toward 0.5 for low-evidence refs. |
| `LL` | Log-likelihood per alignment (average over reads), not total. |
| `dLL` | Change in LL between iterations. |
| `||dphi||` | L∞ norm of φ changes (max absolute change across all refs). |

### Processing Statistics

| Stage | "References" Column Meaning |
|-------|---------------------------|
| Stage 1 | References in BAM header (database total) |
| Stage 2 | Refs with ≥1 alignment passing quality filters |
| Stage 4 | Refs with ≥1 alignment passing probability filter |

### Diagnostic Metrics

| Metric | Definition | Good Value |
|--------|------------|------------|
| `dispersion` | χ²/df ratio from posterior predictive check | ~1.0 |
| `pp_pvalue` | P(χ² > observed \| model) via Wilson-Hilferty | >0.01 |
| `is_overdispersed` | Flag when p<0.01 OR dispersion>3 | False |
