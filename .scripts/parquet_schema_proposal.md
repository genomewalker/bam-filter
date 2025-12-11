# Optimized Parquet Schema for filterBAM

## Problem Statement

Current Parquet output is **larger than SAM.gz** due to:
1. Data duplication (parsed tags + tags_raw)
2. ASCII string storage for sequences/quality
3. Redundant computed columns (end_position, alignment_length)
4. String ref_names instead of dictionary IDs

## Proposed Solution: Normalized Two-Table Design

### Design Principles
1. **No data duplication** - Each piece of data exists in exactly one place
2. **Normalized tables** - Reads (heavy) separate from Alignments (light)
3. **Computed-on-read** - Derive cheap values instead of storing
4. **Dictionary encoding** - IDs instead of strings where possible
5. **2-bit sequence packing** - 4x compression for DNA sequences

---

## Schema Design

### Table 1: `reads.parquet` (One row per unique read)

| Column | Type | Description | Encoding |
|--------|------|-------------|----------|
| read_id | UINT64 | Primary key | RLE |
| read_name | STRING | Optional, for SAM reconstruction | DICT if kept |
| sequence | BINARY | 2-bit packed (A=00,C=01,G=10,T=11) | ZSTD |
| quality | BINARY | Raw Phred bytes (0-93) | ZSTD |
| read_group | STRING | Optional | DICT |

**Size estimate**: ~25 bytes/read (vs ~150 bytes ASCII)

### Table 2: `alignments.parquet` (One row per alignment)

| Column | Type | Description | Encoding |
|--------|------|-------------|----------|
| read_id | UINT64 | FK to reads table | RLE |
| ref_id | UINT32 | Reference index | DICT |
| position | INT32 | 0-based start | DELTA |
| mapq | UINT8 | Mapping quality | RLE |
| flag | UINT16 | SAM flags | BITPACK |
| cigar | STRING | CIGAR string | DICT |
| mate_ref_id | INT32 | Optional | DICT |
| mate_position | INT32 | Optional | DELTA |
| template_length | INT32 | Optional | DELTA |
| tags_cold | STRING | Rare tags only (not AS,NM,MD,XS) | ZSTD |

**Hot tags as columns** (stripped from tags_cold):
| Column | Type | Description |
|--------|------|-------------|
| AS | INT16 | Alignment score |
| NM | UINT16 | Edit distance |
| XS | INT16 | Secondary score |
| MD | STRING | Mismatch string |

**Pipeline results** (written by downstream stages):
| Column | Type | Stage |
|--------|------|-------|
| ani | FLOAT32 | filter |
| pmd_score | FLOAT32 | filter |
| filter_passed | BOOL | filter |
| lca_taxid | INT32 | lca |
| reassigned_ref_id | INT32 | reassign |
| zp_posterior | FLOAT32 | reassign |

**Size estimate**: ~40 bytes/alignment (vs ~200 bytes in SAM)

### Table 3: `references.parquet` (Sidecar, one row per reference)

| Column | Type | Description |
|--------|------|-------------|
| ref_id | UINT32 | Primary key |
| ref_name | STRING | Full reference name |
| ref_length | INT64 | Reference length |
| taxid | INT32 | Optional taxonomy ID |
| domain | STRING | Optional (Bacteria, Archaea, etc.) |

---

## Computed-on-Read (NOT stored)

These values are **derived** from stored columns during queries:

| Derived | Computed From |
|---------|---------------|
| end_position | position + cigar_to_ref_len(cigar) |
| alignment_length | cigar_to_aln_len(cigar) |
| ref_name | JOIN references ON ref_id |

---

## Space Savings Analysis

### Current Schema (31 columns, all data)
```
Per alignment: ~200 bytes compressed
- read_name: 30 bytes
- ref_name: 20 bytes
- sequence: 75 bytes (ASCII)
- quality: 75 bytes (ASCII)
- tags_raw: 50 bytes (duplicates parsed tags)
- parsed tags: 40 bytes (duplicates tags_raw)
- positions/flags: 20 bytes
```

### Proposed Schema (Normalized)
```
Per READ (stored once):
- read_id: 8 bytes
- sequence: 19 bytes (2-bit packed 75bp)
- quality: 40 bytes (ZSTD compressed)
= ~67 bytes per unique read

Per ALIGNMENT:
- read_id: 8 bytes (RLE compressed)
- ref_id: 4 bytes (DICT)
- position: 4 bytes (DELTA)
- mapq/flag: 3 bytes
- cigar: 10 bytes (DICT)
- hot tags: 12 bytes
- cold tags: 10 bytes
= ~51 bytes per alignment
```

### Multi-mapped Read Savings

| Scenario | Current | Proposed | Savings |
|----------|---------|----------|---------|
| 1 alignment/read | 200 bytes | 118 bytes | 41% |
| 5 alignments/read | 1000 bytes | 322 bytes | 68% |
| 20 alignments/read | 4000 bytes | 1087 bytes | 73% |

---

## Implementation Plan

### Phase 1: Quick Wins (No schema change)
1. Remove `tags_raw`, keep only parsed hot tags + `tags_cold`
2. Remove `end_position`, `alignment_length` (compute from CIGAR)
3. Remove `ref_name`, use `ref_id` + sidecar
4. Enable ZSTD level 6+ compression

### Phase 2: Sequence Optimization
1. Implement 2-bit sequence packing
2. Store quality as raw bytes (not ASCII)
3. Add option to drop sequence/quality entirely

### Phase 3: Normalization
1. Split into `reads.parquet` + `alignments.parquet`
2. Deduplicate sequence/quality for multi-mapped reads
3. Add join logic to downstream stages

---

## Query Patterns

### Filter stage (quality filtering)
```sql
SELECT * FROM alignments
WHERE mapq >= 30 AND AS >= -20 AND NM <= 5
```
Fast: Uses hot tag columns directly.

### LCA stage (taxonomy)
```sql
SELECT read_id, array_agg(ref_id) as refs
FROM alignments
GROUP BY read_id
```
Fast: Clustering by read_id.

### Reassign stage (EM)
```sql
SELECT a.*, r.sequence
FROM alignments a
JOIN reads r ON a.read_id = r.read_id
WHERE a.filter_passed = true
```
Needs join, but alignments table is much smaller.

### SAM Export (lossless)
```sql
SELECT r.read_name, r.sequence, r.quality,
       ref.ref_name, a.position, a.cigar,
       concat_tags(a.AS, a.NM, a.MD, a.tags_cold) as tags
FROM alignments a
JOIN reads r ON a.read_id = r.read_id
JOIN references ref ON a.ref_id = ref.ref_id
```
Reconstructs full SAM from normalized tables.

---

## Workflow Integration

### reassign command
- Reads: `alignments.parquet` (filter_passed rows)
- Writes: `reassigned_ref_id`, `zp_posterior` columns
- Needs: sequence from `reads.parquet` for PMD

### filter command
- Reads: `alignments.parquet`
- Writes: `ani`, `pmd_score`, `filter_passed` columns
- Needs: sequence/quality from `reads.parquet`

### lca command
- Reads: `alignments.parquet` (grouped by read_id)
- Writes: `lca_taxid` column
- No need for sequence

---

## Estimated Final Sizes

| Format | Size (1M reads, 5 alignments each) |
|--------|-----------------------------------|
| SAM.gz | ~500 MB |
| BAM | ~400 MB |
| Current Parquet | ~600 MB (worse!) |
| **Proposed Parquet** | **~200 MB** (60% smaller than BAM) |
