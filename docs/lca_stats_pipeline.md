# LCA Stats Architecture

This document explains how `filterBAM lca --stats` transforms trusted LCA
assignments and BAM alignments into per-taxid quality metrics, how the
taxonomy is traversed, and how each metric is aggregated. When `--stats`
is omitted the `--output` file only contains the minimal LCA summary
(taxid/name/rank/n_reads/abundance/tax_path). Supplying `--stats`
augments **the same output file** with the metrics described below—no
additional files are produced.

## 1. Inputs and Trusted Read Filter

```
┌────────────────────┐     ┌──────────────────────────┐
│ BAM (alignments)   │     │ LCA per-read TSV         │
│  • read name       │     │  • read name             │
│  • reference tid   │     │  • LCA taxid             │
│  • tags (NM, AS…)  │     │  • trusted flag          │
└─────────┬──────────┘     └────────────┬─────────────┘
          │                             │
          │                             ▼
          │                 ┌────────────────────────┐
          │                 │ load_lca_assignments() │
          │                 │  → khash(read→LCA)     │
          │                 │  (trusted rows only)   │
          │                 └─────────┬──────────────┘
          │                           │
          ▼                           │
┌───────────────────────────────────────────────────┐
│ compute_trusted_reference_stats()                 │
│  • iterate BAM per reference (htslib iterator)    │
│  • drop alignment if read name ∉ trusted hash     │
│  • reuse stats.pyx logic → RefStats per tid       │
│      (ANI, GC, MAPQ, coverage RLE, entropy, TAD,  │
│       tax_abundances, etc.)                       │
└───────────────────────────────────────────────────┘
```

Only trusted reads ever contribute to `RefStats`, so the downstream aggregation
is perfectly aligned with the LCA mass we propagate through the taxonomy.

## 2. Mapping References to Taxids

Each BAM reference name is resolved through the accession map produced by
`filterBAM build-taxonomy`. That gives us `ref_taxid`, the biological node that
owns the reference sequence.

```
tid ──► sam_hdr_tid2name ──► accession_map ──► ref_taxid
```

## 3. Aggregating Up the Taxonomy

For every reference we:

1. **Accumulate at the leaf** (ref_taxid).  
   * Alignment-weighted stats: read length, GC%, ANI, MAPQ, edit distance,
     aligned length, alignment score.  
   * Length-weighted stats: coverage mean, truncated coverage, breadth, entropy,
     gini, coverage evenness, TAD-based abundances, bases covered.  
   * Per-reference summaries: reference length mean/std/min/max, plus per-ref
     means/stds for coverage, breadth, truncated coverage, and coverage on
     covered bases.  
   * “Best reference” tracking (the one with the most trusted alignments) so we
     can later report medians/modes taken from a concrete contig.

2. **Walk up the lineage** (species → genus → … → root) returned by
   `TaxonomyDatabase.get_lineage` and add the same `RefStats` to every ancestor.

```
ref_stats ──► leaf taxid (species)
            └► ancestor (genus)
               └► ancestor (family)
                  ...
```

All accumulators live in `_accum` dicts so we can later convert them into means,
variances, and CVs.

### Metric Flow Cheat Sheets

**Read-weighted metrics (length, GC, ANI, MAPQ, edit distance, aligned length, AS)**

```
Trusted alignments
     │
     ▼
Per-reference stats (alignment-weighted)
     │
     ▼
Leaf taxid accumulator  ──► ancestor accumulators ──► … ──► root
```

*Each alignment contributes proportionally to the read-level means and pooled
variances. When we walk up the lineage we reuse the same weighted sums.*

**Coverage/TAD metrics (coverage mean, breadth, entropy/gini, cov_evenness, TAD)**

```
Trusted alignments ─► coverage RLE per reference
                        │
                        └─ multiply by reference length
Leaf taxid sums (length-weighted) ─► ancestors ─► root
```

*Coverage-based values are length-weighted so the genome-scale average matches
“bases covered / total length”.*

**Per-reference summaries (`*_per_ref`, reference-length stats)**

```
References under taxid
     │
     ├─ Σ value per reference
     └─ Σ value² per reference
finalize_taxid_entries:
     mean = Σ / n_refs
     std  = sqrt((Σx² - (Σx)²/n) / (n-1))
```

*These capture heterogeneity between descendant references, independent of their
lengths.*

## 4. Propagating LCA Read Counts

Unique read counts are derived directly from the LCA hash (each trusted read
increments the taxid found in the per-read TSV). These counts are then pushed up
the lineage so that:

* every ancestor’s `n_reads` reflects the sum of trusted reads observed anywhere
  below it, and
* the root row’s `n_reads` matches the total number of trusted reads exactly.

Alignments (`n_alns`) and abundance metrics (`tax_abund_*`) are **not** derived
from the LCA file—they come from the trusted-only `RefStats`.

## 5. Finalization

`finalize_taxid_entries()` turns the raw accumulators into the reported values:

* **Alignment-weighted means/std** for read-level metrics (length, GC%, ANI,
  MAPQ, edit distance, aligned length, alignment score). This answers “How do
  the trusted reads assigned to this node behave?”
* **Length-weighted means/std** for coverage metrics (coverage mean, truncated
  coverage, breadth, entropy/gini, coverage evenness, site density). This keeps
  genome-scale interpretations intact—large contigs count proportionally more
  than tiny ones.
* **Per-reference summaries** (`*_per_ref` columns) using simple averages/std
  across references, revealing heterogeneity when descendants behave differently.
* **Reference length stats**: total length (used for abundance) plus mean/std/
  min/max across references.
* **Best reference snapshots** (`read_length_median`, `read_length_mode`,
  `read_ani_median`) drawn from the reference with the most trusted alignments
  so there’s always a concrete exemplar.

## 6. Output

`write_taxid_stats` emits a TSV (or TSV.GZ) sorted lexicographically by
`tax_path`. Each row includes:

* taxonomy columns (`taxid`, `name`, `rank`, `tax_path`)
* trusted read/alignment counts (`n_reads`, `n_alns`)
* reference length totals + per-reference length stats
* alignment-weighted read metrics (length, GC, ANI, MAPQ, edit distance,
  aligned length, AS, DUST low-complexity score)
* length-weighted coverage/TAD metrics (coverage mean, truncated coverage,
  breadth, entropy/gini, coverage evenness, site density)
* per-reference coverage summaries (`coverage_mean_per_ref`, `breadth_per_ref`,
  etc.) to gauge heterogeneity
* abundance metrics (`tax_abund_read`, `tax_abund_aln`, `tax_abund_tad`,
  `n_reads_tad`)

Only taxids with at least one trusted read are written; internal nodes that
contribute coverage but have zero LCA support remain internal to the calculations
but are not shown in the final table.

## Metric Weighting Cheat Sheet

| Category                                | Weighting / Source                                      |
|-----------------------------------------|----------------------------------------------------------|
| `n_reads`                               | From LCA hash; propagated once per read up the lineage  |
| `n_alns`, `tax_abund_*`, `n_reads_tad`  | Sums of trusted-only `RefStats` values                  |
| Read-level means/std (length, GC, ANI, MAPQ, edit distance, AS, aligned length) | Alignment-weighted (each trusted alignment contributes proportionally) |
| Coverage/TAD means/breadth/entropy/gini/cov_evenness/site_density | Reference-length-weighted (or covered length / TAD span) |
| Per-reference coverage/breadth/TAD columns (`*_per_ref`) | Simple mean/std across references                      |
| Reference length totals                 | Sums across references                                  |
| Reference length summaries (`reference_length_mean/std/min/max`) | Simple per-reference stats                              |
| Medians/modes                           | Taken from the “best” reference (most trusted alignments) |

This combination preserves mass and genome-scale intuition at ancestors, while
exposing how uneven the underlying reference set can be. The trusted-read filter
ensures that every number aligns with the per-read LCA assignments.
