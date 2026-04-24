
# filterBAM


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/bam-filter?include_prereleases&label=version)](https://github.com/genomewalker/bam-filter/releases) [![bam-filter](https://github.com/genomewalker/bam-filter/workflows/filterBAM_ci/badge.svg)](https://github.com/genomewalker/bam-filter/actions) [![PyPI](https://img.shields.io/pypi/v/bam-filter)](https://pypi.org/project/bam-filter/) [![Conda](https://img.shields.io/conda/v/genomewalker/bam-filter)](https://anaconda.org/genomewalker/bam-filter)


A simple tool to process a BAM file, filter references with uneven coverages, and estimate taxonomic abundances. FilterBAM has three main goals:

1. **Reassign reads** to the reference they belong using an E-M algorithm that takes into account the alignment score. 

2. **Estimate several metrics** for each reference in the BAM file and filter those references that do not meet the defined criteria.

3. **Perform an LCA (Last Common Ancestor) analysis** using the reads that pass the filtering criteria and estimate genome normalized abundances at the rank specified.



# Installation

We recommend managing the native dependencies with [**micromamba/conda**](https://mamba.readthedocs.io/en/latest/installation.html); the repository ships an environment specification that provides the compilers, HTSlib, igraph, Arrow, and DuckDB libraries required to compile the Cython extensions.

### Latest release (PyPI)

```bash
micromamba create -n bam-filter -f https://raw.githubusercontent.com/genomewalker/bam-filter/master/environment.yml
micromamba activate bam-filter
pip install bam-filter
```

### Development version

Install directly from the repository:

```bash
micromamba create -n bam-filter-dev -f environment.yml
micromamba activate bam-filter-dev
pip install -e .[dev,test]
pytest --maxfail=1 --disable-warnings -vv
```

### Building distribution artifacts

After activating the environment:

```bash
python -m build
```

This produces wheels and source distributions in `dist/` using the PEP 517 backend defined in `pyproject.toml`.


# Usage

filterBAM only needs a BAM file. For a complete list of options:

```
$ filterBAM --help
usage: filterBAM [-h] [--version] [--debug] {reassign,filter,lca} ...

A simple tool to calculate metrics from a BAM file and filter with uneven coverage.

positional arguments:
  {reassign,filter,lca}
                        positional arguments
    reassign            Reassign reads to references using an EM algorithm
    filter              Filter references based on coverage and other metrics
    lca                 Calculate LCA for each read and estimate abundances at each rank

optional arguments:
  -h, --help            show this help message and exit
  --version             Print program version
  --debug               Print debug messages (default: False)
```


## Read reassignment

Full list of options:


```bash
$ filterBAM reassign --help
usage: filterBAM reassign [-h] --bam BAM [-p STR] [-r FILE] [-t INT] [-i INT] [-s FLOAT] [-A FLOAT] [-l INT] [-L INT] [-n INT] [--match-reward INT] [--mismatch-penalty INT] [--gap-open-penalty INT] [--gap-extension-penalty INT] [--squarem-min-improvement FLOAT] [--squarem-max-step-factor FLOAT] [--anderson-interval INT]
                         [-o [FILE]] [-m STR] [-M INT] [-N] [--tmp-dir DIR] [--disable-sort]

optional arguments:
  -h, --help            show this help message and exit
  -p STR, --prefix STR  Prefix used for the output files (default: None)
  -r FILE, --reference-lengths FILE
                        File with references lengths (default: None)
  -t INT, --threads INT
                        Number of threads to use (default: 1)

required arguments:
  --bam BAM             BAM file containing aligned reads (default: None)

Re-assign optional arguments:
  -i INT, --iters INT   Number of iterations for the EM algorithm (default: 25)
  -s FLOAT, --scale FLOAT
                        Scale to select the best weighting alignments (default: 0.9)
  -A FLOAT, --min-read-ani FLOAT
                        Minimum read ANI to keep a read (default: 90.0)
  -l INT, --min-read-length INT
                        Minimum read length (default: 30)
  -L INT, --max-read-length INT
                        Maximum read length (default: inf)
  -n INT, --min-read-count INT
                        Minimum read count (default: 3)
  --match-reward INT    Match reward for the alignment score (default: 1)
  --mismatch-penalty INT
                        Mismatch penalty for the alignment score (default: -2)
  --gap-open-penalty INT
                        Gap open penalty for alignment score computation (default: 5)
  --gap-extension-penalty INT
                        Gap extension penalty for the alignment score (default: 2)
  --squarem-min-improvement FLOAT
                        Minimum relative improvement for SQUAREM convergence (default: 0.0001)
  --squarem-max-step-factor FLOAT
                        Maximum step size multiplier for SQUAREM stability (default: 4.0)
  --anderson-interval INT
                        Apply Anderson acceleration every N EM iterations (default: 5; set to 0 to disable Anderson)
  -o [FILE], --out-bam [FILE]
                        Save a BAM file without multimapping reads (default: None)
  -m STR, --sort-memory STR
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default: 1G)
  -M INT, --max-memory INT
                        Maximum memory to use for the EM algorithm (default: None)
  -N, --sort-by-name    Sort by read names (default: False)
  --tmp-dir DIR         Temporary directory (default: None)

miscellaneous arguments:
  --disable-sort        Disable sorting of the filtered BAM file (default: False)
```


One would run filterBAM `reassign` as follows:

```bash
filterBAM reassign --bam c55d4e2df1.dedup.bam  --threads 10 --iters 0 --min-read-ani 92 --reference-lengths gtdb-r202.len.map --out-bam c55d4e2df1.reassigned.bam
```


**--bam**: BAM file to process
**--threads**: Number of threads to use
**--iters**: Number of iterations for the EM algorithm. If set to 0, the EM algorithm will run until there are no more reads to reassign
**--min-read-ani**: Minimum read ANI to keep a read
**--reference-lengths**: File with the lengths of the references in the BAM file. This is used when multiple contigs have been concatenated with Ns.
**--out-bam**: Save a BAM file without multimapping reads
**--anderson-interval**: Apply Anderson acceleration every N EM iterations (default: 5). Set to 0 to disable Anderson acceleration entirely. Increasing this value can reduce overhead for large datasets.


## BAM filtering

Full list of options:

```bash
$ filterBAM filter --help
usage: filterBAM filter [-h] --bam BAM [-p STR] [-r FILE] [-t INT] [--reference-trim-length INT] [--trim-min INT] [--trim-max INT] [-A FLOAT] [-l INT] [-L INT] [-n INT] [-b FLOAT] [-e FLOAT]
                        [-g FLOAT] [-B FLOAT] [-a FLOAT] [-c FLOAT] [-V FLOAT] [-C FLOAT] [--include-low-detection] [-m STR] [-N] [--disable-sort] [--scale STR] --stats FILE
                        [--stats-filtered FILE] [--bam-filtered FILE] [--read-length-freqs [FILE]] [--read-hits-count [FILE]] [--knee-plot [FILE]] [--coverage-plots [FILE]] [--tmp-dir DIR]
                        [--low-memory]

optional arguments:
  -h, --help            show this help message and exit
  -p STR, --prefix STR  Prefix used for the output files (default: None)
  -r FILE, --reference-lengths FILE
                        File with references lengths (default: None)
  -t INT, --threads INT
                        Number of threads to use (default: 1)

required arguments:
  --bam BAM             BAM file containing aligned reads (default: None)

Filter required arguments:
  --stats FILE          Save a TSV file with the statistics for each reference (default: None)

filtering arguments:
  -A FLOAT, --min-read-ani FLOAT
                        Minimum read ANI to keep a read (default: 90.0)
  -l INT, --min-read-length INT
                        Minimum read length (default: 30)
  -L INT, --max-read-length INT
                        Maximum read length (default: inf)
  -n INT, --min-read-count INT
                        Minimum read count (default: 3)
  -b FLOAT, --min-expected-breadth-ratio FLOAT
                        Minimum expected breadth ratio (default: 0)
  -e FLOAT, --min-normalized-entropy FLOAT
                        Minimum normalized entropy (default: 0)
  -g FLOAT, --min-normalized-gini FLOAT
                        Minimum normalized Gini coefficient (default: 1.0)
  -B FLOAT, --min-breadth FLOAT
                        Minimum breadth (default: 0)
  -a FLOAT, --min-avg-read-ani FLOAT
                        Minimum average read ANI (default: 90.0)
  -c FLOAT, --min-coverage-evenness FLOAT
                        Minimum coverage evenness (default: 0)
  -V FLOAT, --min-coeff-var FLOAT
                        Minimum coverage evenness calculated as SD/MEAN (default: inf)
  -C FLOAT, --min-coverage-mean FLOAT
                        Minimum coverage mean (default: 0)
  --include-low-detection
                        Include those references that fulfill all filtering criteria but the coverage evenness is 0 (default: False)

miscellaneous arguments:
  --reference-trim-length INT
                        Exclude n bases at the ends of the reference sequences (default: 0)
  --trim-min INT        Remove coverage that are below this percentile. Used for the Truncated Average Depth (TAD) calculation (default: 10)
  --trim-max INT        Remove coverage that are above this percentile. Used for the Truncated Average Depth (TAD) calculation (default: 90)
  -m STR, --sort-memory STR
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default: 1G)
  -N, --sort-by-name    Sort by read names (default: False)
  --disable-sort        Disable sorting of the filtered BAM file (default: False)
  --scale STR           Scale taxonomic abundance by this factor; suffix K/M recognized (default: 1000000.0)
  --tmp-dir DIR         Temporary directory (default: None)
  --low-memory          Activate the low memory mode (default: False)

output arguments:
  --stats-filtered FILE
                        Save a TSV file with the statistics for each reference after filtering (default: None)
  --bam-filtered FILE
                        Save a BAM file with the references that passed the filtering criteria (default: None)
  --read-length-freqs [FILE]
                        Save a JSON file with the read length frequencies mapped to each reference (default: None)
  --read-hits-count [FILE]
                        Save a TSV file with the read hits frequencies mapped to each reference (default: None)
  --knee-plot [FILE]    Plot knee plot (default: None)
  --coverage-plots [FILE]
                        Folder where to save genome coverage plots (default: None)
```

One would run filterBAM `filter` as follows:


```bash
filterBAM filter --bam c55d4e2df1.reassigned.bam --bam-filtered c55d4e2df1.dedup.filtered.bam --stats c55d4e2df1.dedup.stats.tsv.gz --stats-filtered c55d4e2df1.dedup.stats-filtered.tsv.gz --threads 10 --min-read-ani 92 --min-normalized-entropy 0.6
```

**--stats**: Save a TSV file with the statistics for each reference

**--min-read-count**: Minimum number of reads mapped to a reference in the BAM file

**--min-expected-breadth-ratio**: Minimum expected breadth ratio needed to keep a reference.

**--min-read-ani**: Minimum average read ANI that a reference has

**--sort-by-name**: Sort filtered BAM file by read name so it can be used in [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp)

**--sort-memory**: Memory used for each thread when sorting the filtered BAM file

**--reference-lengths**: File with the lengths of the references in the BAM file. This is used to calculate the coverage estimates of each reference when multiple contigs have been concatenad with Ns.

**--threads**: Number of threads


The program will produce two main outputs:
 - A BAM file where the references that are below the defined threshold have been filtered out
 - A TSV file with statistics for each reference, with the following columns:
    - **reference**: Reference name
    - **n_reads**: Number of reads mapped to the reference
    - **n_alns**: Number of alignments in the reference
    - **read_length_mean**: Mean read length mapped to the reference
    - **read_length_std**: Standard deviation of read lengths mapped to the reference
    - **read_length_min**: Minimum read length mapped to the reference
    - **read_length_max**: Maximum read length mapped to the reference
    - **read_length_median**: Medium read length mapped to the reference
    - **read_length_mode**: Modal read length mapped to the reference
    - **gc_content**: Average GC content of the reads mapped to the reference
    - **read_aligned_length**: Average aligned read length mapped to the reference
    - **read_aln_score**: Average alignment score of the reads mapped to the reference
    - **mapping_quality**: Average mapping quality of the reads mapped to the reference
    - **edit_distances**: Average edit distance of the reads mapped to the reference
    - **read_ani_mean**: Average ANI of the reads mapped to the reference
    - **read_ani_std**: Standard deviation of ANI of the reads mapped to the reference
    - **read_ani_median**: Median ANI of the reads mapped to the reference
    - **bases_covered**: Number of bases covered by the reference
    - **max_covered_bases**: Maximum number of bases covered in the reference
    - **mean_covered_bases**: Average number of bases covered in the reference
    - **coverage_mean**: Mean depth of the reference
    - **coverage_mean_trunc**: Mean depth of the reference after removing the 10% and 90% of the coverage values (default: TAD80 as calculated [here](https://sfamjournals.onlinelibrary.wiley.com/doi/10.1111/1462-2920.15112))
    - **coverage_mean_trunc_len**: Length of the reference after being truncated by the TAD(X) values
    - **coverage_covered_mean**: Mean depth of the reference only counting covered bases
    - **reference_length**: Real reference length
    - **bam_reference_length**: Length reported by the BAM file
    - **breadth**: Breadth of coverage 
    - **exp_breadth**: Expected breadth of coverage. Using the equation: _expected_breadth = 1 - e<sup>-coverage</sup>_
    - **breadth_exp_ratio**: Ration between the observed and the expected depth of coverage
    - **n_bins**: Number of bins used to calculate the read coverage distribution
    - **site_density**: Site density of the reference
    - **entropy**: Entropy of the read coverage distribution
    - **norm_entropy**: Normalized entropy of the read coverage distribution
    - **gini**: Gini coefficient of the read coverage distribution
    - **norm_gini**: Normalized Gini coefficient of the read coverage distribution
    - **c_v**: Coefficient of variation of the coverage
    - **d_i**: Dispersion index
    - **cov_evenness**: Eveness of coverage as calculated [here](https://www.nature.com/articles/jhg201621).
    - **tax_abund_read**: Counts estimated using the number of reads and normalized by the reference length.
    - **tax_abund_aln**: Counts estimated using the number of alignments and normalized by the reference length.
    - **tax_abund_tad**: Counts estimated using the estimated number of reads in the TAD region and normalized by the length of the TAD region
    - **n_reads_tad**: Number of reads estimated in the TAD region using the equation *C = LN / G*, where C stands for the TAD coverage, N for the length of the TAD region and L for the average read length mapped to the reference.
    - **n_intervals**: Number of distinct coverage intervals (contiguous regions with coverage > 0)
    - **weighted_contiguity_breadth**: Weighted Contiguity Breadth (WCB) metric for detecting scattered false positives
    - **complexity_penalized_coverage**: Complexity-Penalized Coverage (CPC) for filtering low-complexity false positives
    - **overlap_redundancy_index**: Overlap Redundancy Index (ORI) measuring read stacking density


## Statistics Equations and Descriptions

This section provides detailed mathematical definitions for all statistics calculated by filterBAM.

### Read-Level Statistics

| Statistic | Equation | Description |
|-----------|---------|-------------|
| **read_length_mean** | $\bar{L} = \frac{1}{n}\sum_{i=1}^{n} L_i$ | Mean length of reads mapped to the reference |
| **read_length_std** | $\sigma_L = \sqrt{\frac{1}{n-1}\sum_{i=1}^{n}(L_i - \bar{L})^2}$ | Standard deviation of read lengths |
| **read_gc_content_mean** | $\bar{GC} = \frac{1}{n}\sum_{i=1}^{n} \frac{G_i + C_i}{L_i}$ | Mean GC content (proportion of G+C bases) |
| **dust_mean** | $\bar{D} = \frac{1}{n}\sum_{i=1}^{n} D_i$ | Mean DUST score (low-complexity indicator; higher = more repetitive) |
| **read_ani_mean** | $\bar{ANI} = \frac{1}{n}\sum_{i=1}^{n} \left(1 - \frac{edit\_dist_i}{aligned\_len_i}\right)$ | Mean Average Nucleotide Identity of reads |
| **aligned_length_mean** | $\bar{A} = \frac{1}{n}\sum_{i=1}^{n} A_i$ | Mean aligned length (excluding soft-clipped bases) |

### Coverage Statistics

| Statistic | Equation | Description |
|-----------|---------|-------------|
| **coverage_mean** | $\bar{C} = \frac{\sum_{p=1}^{L} c_p}{L}$ | Mean depth across reference length $L$ |
| **coverage_mean_trunc** | TAD = mean of coverage values between $P_{10}$ and $P_{90}$ | Truncated Average Depth (robust to outliers) |
| **coverage_covered_mean** | $\bar{C}_{cov} = \frac{\sum_{p: c_p > 0} c_p}{\|{p: c_p > 0}\|}$ | Mean depth only over covered positions |
| **bases_covered** | $B = \|{p : c_p > 0}\|$ | Number of positions with coverage > 0 |
| **breadth** | $b = \frac{B}{L}$ | Fraction of reference covered |
| **exp_breadth** | $b_{exp} = 1 - e^{-\bar{C}}$ | Expected breadth given mean coverage (Lander-Waterman) |
| **breadth_exp_ratio** | $R = \min\left(\frac{b}{b_{exp}}, 1\right)$ | Ratio of observed to expected breadth |

### Coverage Distribution Metrics

| Statistic | Equation | Description |
|-----------|---------|-------------|
| **spatial_entropy** | $H = -\sum_{i=1}^{k} p_i \log_2(p_i)$ where $p_i = \frac{c_i}{\sum c}$ | Shannon entropy of binned coverage histogram |
| **norm_spatial_entropy** | $H_{norm} = \frac{H}{\log_2(k)}$ | Normalized entropy (0-1 scale); 1 = perfectly uniform |
| **gini** | $G = \frac{\sum_{i=1}^{k}\sum_{j=1}^{k}\|x_i - x_j\|}{2k^2\bar{x}}$ | Gini coefficient of coverage inequality |
| **norm_gini** | $G_{norm} = 1 - G$ | Normalized Gini (1 = perfectly uniform) |
| **c_v** | $CV = \frac{\sigma}{\mu}$ | Coefficient of variation of coverage |
| **d_i** | $D = \frac{\sigma^2}{\mu}$ | Dispersion index (variance-to-mean ratio) |
| **cov_evenness** | $E = 1 - G$ | Coverage evenness (Pielou's J adapted) |

### Contamination Detection Metrics

These metrics help identify false positive taxonomic assignments caused by low-complexity sequences mapping to scattered genomic positions.

| Statistic | Equation | Description |
|-----------|---------|-------------|
| **n_intervals** | $N_{int} = $ count of contiguous coverage regions | Number of distinct coverage intervals |
| **weighted_contiguity_breadth** | $WCB = \frac{\sum_{i=1}^{N_{int}} \ell_i^2}{L^2}$ | Rewards long contiguous intervals; penalizes scattered tiny hits |
| **complexity_penalized_coverage** | $CPC = b \times (1 - \bar{D})$ | Breadth penalized by sequence complexity |
| **overlap_redundancy_index** | $ORI = \frac{\sum_{i=1}^{n} A_i}{B}$ | Total aligned bases divided by bases covered |

#### Interpretation Guide for Contamination Detection

| Metric | False Positives (contamination) | True Positives (real signal) |
|--------|--------------------------------|------------------------------|
| **WCB** | ≈ 0 (scattered tiny hits) | > 0.001 (contiguous coverage) |
| **CPC** | < 0.01 (low complexity) | > 0.1 (high complexity reads) |
| **dust_mean** | > 0.05 (repetitive) | < 0.03 (complex sequences) |
| **breadth** | ≈ 0 | > 0.1 |

**Example filter to remove false positives:**
```bash
filterBAM filter --bam input.bam --stats stats.tsv.gz \
    --filter complexity_penalized_coverage:0.01:,weighted_contiguity_breadth:0.0001:
```

### Abundance Metrics

| Statistic | Equation | Description |
|-----------|---------|-------------|
| **tax_abund_read** | $A_r = \frac{n \times S}{L}$ | Read-count abundance (reads per million bp) |
| **tax_abund_aln** | $A_a = \frac{n_{aln} \times S}{L}$ | Alignment-count abundance |
| **tax_abund_tad** | $A_{tad} = \frac{n_{tad} \times S}{L_{tad}}$ | TAD-normalized abundance (robust to coverage bias) |

Where $S$ is the scale factor (default: 1,000,000 for RPM normalization).

## LCA

Full list of options:
  
```bash
$ filterBAM lca --help
usage: filterBAM lca [-h] --bam BAM --taxonomy-db DIR [-p STR] [-r FILE] [-t INT] [--lca-rank STR] [--lca-summary [FILE]] [--scale STR] [-m STR] [--custom] [--stats [FILE]]

optional arguments:
  -h, --help            show this help message and exit
  -p STR, --prefix STR  Prefix used for the output files (default: None)
  -r FILE, --reference-lengths FILE
                        File with references lengths (default: None)
  -t INT, --threads INT
                        Number of threads to use (default: 1)

required arguments:
  --bam BAM             BAM file containing aligned reads (default: None)
  --taxonomy-db DIR     Directory with Parquet taxonomy database (includes accession_map.parquet)

LCA optional arguments:
  --lca-rank STR        Rank to use for LCA calculation (default: species)
  --lca-summary [FILE]  Save a TSV file with the LCA summary (default: None)
  --scale STR           Scale taxonomic abundance by this factor; suffix K/M recognized (default: 1000000.0)
  -m STR, --sort-memory STR
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default: 1G)
  --custom              Use custom taxdump files (default: False)
  --stats [FILE]        A TSV file from the filter subcommand (default: None)
```


> If you use the `--stat` option the LCA will use, when the possible, the reads inferred after calculating the TAD abundances. This is useful when the reads are not uniformly distributed across the reference. The program will produce a TSV file with the LCA summary.

One would run filterBAM `lca` as follows:

```bash
filterBAM lca --bam c55d4e2df1.dedup.filtered.bam --taxonomy-db ./taxonomy_db --threads 10 --lca-rank genus
```

**--taxonomy-db**: Path to the Parquet taxonomy database produced by `filterBAM build-taxonomy` (must include `accession_map.parquet`)

**--rank-lca**: Rank to use for LCA calculation 

**--scale**: Scale taxonomic abundance by this factor; suffix K/M recognized 


# Read Reassignment Algorithm in FilterBAM


FilterBAM uses an Expectation-Maximization (EM) algorithm to probabilistically reassign multi-mapping reads based on alignment quality. Anderson acceleration can be enabled to speed up convergence, especially for large and complex datasets.

### EM Initialization, Scoring, and Algorithm

FilterBAM uses an Expectation-Maximization (EM) algorithm to reassign multi-mapping reads based on alignment quality. The process is as follows:

#### 1. Alignment Scoring
For each alignment, a raw score $S$ is calculated as:

$S = r_m \cdot M - p_m \cdot X - p_o \cdot G - p_e \cdot E$

where:
- $M$: number of matches
- $X$: number of mismatches
- $G$: number of gap opens
- $E$: number of gap extensions
- $r_m$, $p_m$, $p_o$, $p_e$: user-configurable parameters (see CLI options)

To compare alignments of different lengths, the score is normalized:

$S'_{ij} = \frac{S_{ij}}{L_{ij}}$

where $L_{ij}$ is the alignment length for read $i$ and reference $j$.

#### 2. Softmax Normalization (Initialization)
For each read, the set of relative scores $S'_{ij}$ (across all references $j$) is transformed into probabilities using a softmax normalization:

$P_{ij} = \frac{\exp\left(\frac{S'_{ij} - \max_k S'_{ik}}{\sigma_i}\right)}{\sum_{k} \exp\left(\frac{S'_{ik} - \max_k S'_{ik}}{\sigma_i}\right)}$

where $\sigma_i$ is the range of $S'_{ik}$ for read $i$ (or 1 if the range is 0), and $\max_k S'_{ik}$ is subtracted for numerical stability. This ensures the probabilities for all alignments of a read sum to 1. These softmax-normalized probabilities are used as the initial probabilities for the EM algorithm, reflecting the relative quality of each alignment for a given read.

#### 3. EM Algorithm Steps

- **Initialization:** Use the softmax-normalized probabilities as the starting point for the EM algorithm. The softmax is only used for initialization; subsequent EM steps use linear normalization.

- **E-step (Expectation):** For each alignment, update the probability (responsibility) that a read originated from a reference, using the current weights $w_j$ and the alignment scores $S_{ij}$:
  
  $P(r_i|g_j) = \frac{S_{ij} w_j}{\sum_k S_{ik} w_k}$
  
  This is a linear normalization (not a softmax) and is repeated for each EM iteration. For each read $i$, the probabilities are normalized so that $\sum_j P(r_i|g_j) = 1$.

- **M-step (Maximization):** Update the weights for each reference by summing the responsibilities over all reads:
  
  $w_j = \frac{\sum_i P(r_i|g_j)}{\sum_{i,j} P(r_i|g_j)}$

- **Anderson Acceleration (optional):** To speed up convergence, FilterBAM can use Anderson acceleration, which extrapolates a better solution for the weights using a history of previous EM steps. This is user-configurable and only accepted if it improves the likelihood and produces valid weights.

- **Convergence and Filtering:** The EM process repeats E and M steps (optionally with Anderson acceleration) until weights stabilize or a maximum number of iterations is reached. After convergence, alignments are filtered: only those with probability above a minimum threshold and a fraction of the maximum for that read are retained. Probabilities are written to the output BAM as a tag for downstream analysis.

#### Practical Notes

- Anderson acceleration can greatly reduce the number of EM iterations required for convergence, especially for large and complex datasets. However, it introduces some computational overhead, so the `--anderson-interval` option allows users to balance speed and resource usage.
- Score normalization ensures that all alignments are fairly compared, regardless of their length or raw score scale.
- The EM process logs the likelihood and convergence status at each iteration, helping users monitor progress and diagnose issues.

#### Anderson Acceleration (in practice)

Anderson acceleration works by combining several previous EM steps to extrapolate a better solution for the weights. In mathematical terms, suppose $w^{(t)}$ is the vector of weights at EM iteration $t$, and $F(w)$ is the EM update operator (i.e., $w^{(t+1)} = F(w^{(t)})$). Anderson acceleration computes the next iterate as a linear combination of the most recent $m$ EM iterates and their residuals:

$$
w^{(t+1)} = (1 - \beta) F(w^{(t)}) + \beta \sum_{j=0}^{m-1} \alpha_j F(w^{(t-j)})
$$

where the coefficients $\alpha_j$ and $\beta$ are chosen to minimize the norm of the residuals $F(w^{(t-j)}) - w^{(t-j)}$ over the last $m$ steps (subject to $\sum_j \alpha_j = 1$). In practice, this means Anderson acceleration extrapolates a new solution from a small window of previous EM steps, and only accepts it if it improves the likelihood and produces valid weights. The `--anderson-interval` option allows users to control how often it is applied. All steps are implemented efficiently for large-scale data, using Numba and memory-mapped arrays.

### Applications and recommendations

One of the main applications of **bam-filter** is to reliably identify which potential organisms are present in a metagenomic ancient sample, and get relatively accurate taxonomic abundances, even when they are present in very low abundances. The resulting BAM file then can be used as input for [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp). We rely on several measures to discriminate between noise and a potential signal, analyzing the mapping results at two different levels:

- Is the observed breadth aligned with the expected one?
- Are the reads spread evenly across the reference or they are clumped in a few regions?

To assess the first question we use the concepts defined [here](https://doi.org/10.1016/0888-7543(88)90007-9). We estimate the ratio between the observed and expected breadth as a function of the coverage. If we get a **breadth_exp_ratio** close to 1, it means that the coverage we observed is close to the one we expect based on the calculated coverage. While this measure is already a strong indicator of a potential signal, we complement it with the metrics that measure the **normalized positional entropy** and the **normalized distribution inequality** (Gini coefficient) of the positions in the coverage. For details on how are calculated check [here](https://www.frontiersin.org/articles/10.3389/fmicb.2022.918015/full). These two metrics will help to identify those cases where we get a high **breadth_exp_ratio** but the coverage is not evenly distributed across the reference but instead is clumped in a few regions. One thing to be aware of is that we need to bin the reference to calculate those metrics. In our case, we use the ability of [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html) to identify the [numbers of bins](https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges), either using the Sturges or the Freedman-Diaconis rule. Finally, we use the [knee point detection algorithm](https://github.com/arvkevi/kneed) to identify the optimal values where to filter the Gini coefficient as a function of the positional entropy.

## LCA genome normalized abundances

The taxonomic abundances for each rank are estimated by normalizing the number of reads by the length of the reference. The LCA approach ranks taxonomic paths based on likelihood, selecting the most probable reference. The program also uses the TAD (Truncated Average Depth) estimated reads for the LCA analysis, which helps mitigate the effect of uneven coverages across references.
