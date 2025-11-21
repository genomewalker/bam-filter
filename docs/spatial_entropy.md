# Spatial Entropy: Measuring Positional Distribution Uniformity

## Overview

**Spatial entropy** measures how evenly covered positions are distributed across a reference sequence. It is a key metric for detecting coverage pileups and assessing alignment quality in metagenomic analyses.

## What It Measures

Spatial entropy quantifies the **uniformity of the spatial distribution** of covered positions, **NOT** the uniformity of coverage depths.

### Example

Consider a reference of 10,000 bp:

**Scenario 1: High Spatial Entropy (Good)**
```
Reference: |====|====|====|====|====|
Coverage:  | ** | ** | ** | ** | ** |
```
- Covered positions evenly distributed across the reference
- `norm_spatial_entropy ≈ 0.9-1.0`
- Interpretation: Good, uniform coverage

**Scenario 2: Low Spatial Entropy (Pileup)**
```
Reference: |====|====|====|====|====|
Coverage:  |****|    |    |    |    |
```
- All covered positions clustered in one region
- `norm_spatial_entropy ≈ 0.2-0.3`
- Interpretation: Pileup detected, likely spurious hit

## Calculation Method

### 1. **Histogram Construction**

We create a histogram of **positions** (not depths) across the reference:

```python
# For each covered position (where depth > 0):
#   - Bin the position spatially across the reference
#   - Count how many positions fall in each bin
```

The number of bins is determined using NumPy's "auto" binning (Freedman-Diaconis or Sturges rule).

### 2. **Raw Spatial Entropy**

Calculate Shannon entropy from the position histogram:

```
H = -Σ p_i × log(p_i)
```

where `p_i` is the proportion of covered positions in bin `i`.

### 3. **Normalized Spatial Entropy**

We normalize by the **achievable maximum entropy** given constraints:

```
norm_spatial_entropy = H / H_max_achievable
```

where `H_max_achievable` accounts for:
- Total number of covered positions
- Number of histogram bins
- Discrete distribution constraints

This is more rigorous than normalizing by `log(n_bins)` (theoretical maximum), because with discrete counts, perfect uniformity may not be achievable.

## Comparison to Other Entropy Metrics

### **Unicorn's Depth Entropy** (Coverage Depth Distribution)

```c
// Unicorn measures: "How uniform are the DEPTH VALUES?"
// Example: depths [1x, 1x, 2x, 2x, 5x]
// Histogram: {1: 2 bases, 2: 2 bases, 5: 1 base}
```

**Use case:** Detecting PCR bias, uneven amplification

### **Our Spatial Entropy** (Position Distribution)

```python
# We measure: "How evenly are POSITIONS distributed across the reference?"
# Example: positions [10, 11, 50, 51, 100] on 1000bp reference
# Bins positions spatially across reference
```

**Use case:** Detecting pileups, partial alignments, spurious hits

## Use Cases for Pileup Detection

### Filtering Strategy

```python
# Good alignments should have:
if breadth > 0.01 and norm_spatial_entropy > 0.6 and n_reads > 5:
    # High confidence hit

# Flag as potential pileup:
if norm_spatial_entropy < 0.4 and n_reads < 10:
    # Likely spurious hit - reads clustered in one region
```

### Real-World Example

From actual data (`Lib_KapK12135_collapsed.stats.tsv.gz`):

| Reference | n_reads | bases_covered | norm_spatial_entropy | Interpretation |
|-----------|---------|---------------|----------------------|----------------|
| `CBDBUJ010002393.1` | 1 | 30 | **0.2485** | Single read, localized hit |
| `JAKDEW010000738.1` | 1 | 120 | **0.2838** | Still clustered |
| Multi-read alignment | 50 | 5000 | **0.85** | Good spatial distribution |

## Edge Cases

| Scenario | norm_spatial_entropy | Behavior |
|----------|---------------------|----------|
| **No coverage** | 0.0 | Set to 0 (no data) |
| **1 position covered** | 1.0 | Early exit (undefined) |
| **Single read** | Low (0.2-0.4) | Clustered by definition |
| **Uniform coverage** | High (0.8-1.0) | Ideal distribution |

## Implementation Details

- **File:** `bam_filter/stats_rle.pyx`
- **Functions:**
  - `calculate_spatial_entropy()` - Raw entropy calculation
  - `calculate_normalized_spatial_entropy()` - Normalized [0,1] metric
- **Struct field:** `RefStats.spatial_entropy`, `RefStats.norm_spatial_entropy`
- **Output columns:** `spatial_entropy`, `norm_spatial_entropy`

## Key Advantages

1. **Pileup Detection:** Directly identifies clustered coverage patterns
2. **Quality Assessment:** Distinguishes between spurious and genuine alignments
3. **Mathematically Rigorous:** Accounts for discrete distribution constraints
4. **Complementary:** Works alongside breadth and coverage depth metrics

## References

- Shannon, C.E. (1948). "A Mathematical Theory of Communication"
- NumPy histogram binning: Freedman-Diaconis and Sturges rules
- Unicorn entropy (depth-based): https://github.com/GeoGenetics/unicorn
