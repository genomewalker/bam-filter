
# filterBAM


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/bam-filter?include_prereleases&label=version)](https://github.com/genomewalker/bam-filter/releases) [![bam-filter](https://github.com/genomewalker/bam-filter/workflows/filterBAM_ci/badge.svg)](https://github.com/genomewalker/bam-filter/actions) [![PyPI](https://img.shields.io/pypi/v/bam-filter)](https://pypi.org/project/bam-filter/) [![Conda](https://img.shields.io/conda/v/genomewalker/bam-filter)](https://anaconda.org/genomewalker/bam-filter)


A simple tool to process a BAM file, filter references with uneven coverages, and estimate taxonomic abundances. FilterBAM has three main goals:

1. **Reassign reads** to the reference they belong using an E-M algorithm that takes into account the alignment score. 

2. **Estimate several metrics** for each reference in the BAM file and filter those references that do not meet the defined criteria.

3. **Perform an LCA (Last Common Ancestor) analysis** using the reads that pass the filtering criteria and estimate genome normalized abundances at the rank specified.



# Installation

We recommend having [**conda**](https://docs.conda.io/en/latest/) installed to manage the virtual environments

### Using pip

First, we create a conda virtual environment with:

```bash
wget https://raw.githubusercontent.com/genomewalker/bam-filter/master/environment.yml
conda env create -f environment.yml
```

Then we proceed to install using pip:

```bash
pip install bam-filter
```

### Install from source to use the development version

Using pip

```bash
pip install git+https://github.com/genomewalker/bam-filter.git
```

By cloning in a dedicated conda environment

```bash
git clone https://github.com/genomewalker/bam-filter.git
cd bam-filter
conda env create -f environment.yml
conda activate bam-filter
pip install -e .
```


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
usage: filterBAM reassign [-h] --bam BAM [-p STR] [-r FILE] [-t INT] [-i INT] [-s FLOAT] [-A FLOAT] [-l INT] [-L INT] [-n INT] [--match-reward INT] [--mismatch-penalty INT]
                          [--gap-open-penalty INT] [--gap-extension-penalty INT] [--lambda FLOAT] [-k FLOAT] [-o [FILE]] [-m STR] [-M INT] [-N] [--tmp-dir DIR] [--disable-sort]

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
                        Scale to select the best weithing alignments (default: 0.9)
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
  --lambda FLOAT        Lambda parameter for the alignment score (default: 1.33)
  -k FLOAT              K parameter for the alignment score algorithm (default: 0.621)
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
**--reference-lengths**: File with the lengths of the references in the BAM file. This is used when multiple contigs have been concatenad with Ns.
**--out-bam**: Save a BAM file without multimapping reads


## BAM filtering

Full list of options:

```bash
$ filterBAM filter --help
usage: filterBAM filter [-h] --bam BAM [-p STR] [-r FILE] [-t INT] [--reference-trim-length INT] [--trim-min INT] [--trim-max INT] [-A FLOAT] [-l INT] [-L INT] [-n INT] [-b FLOAT] [-e FLOAT]
                        [-g FLOAT] [-B FLOAT] [-a FLOAT] [-c FLOAT] [-V FLOAT] [-C FLOAT] [--include-low-detection] [-m STR] [-N] [--disable-sort] [--scale STR] --stats [FILE]
                        [--stats-filtered [FILE]] [--bam-filtered [FILE]] [--read-length-freqs [FILE]] [--read-hits-count [FILE]] [--knee-plot [FILE]] [--coverage-plots [FILE]] [--tmp-dir DIR]
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
  --stats [FILE]        Save a TSV file with the statistics for each reference (default: None)

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
  --stats-filtered [FILE]
                        Save a TSV file with the statistics for each reference after filtering (default: None)
  --bam-filtered [FILE]
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

## LCA

Full list of options:
  
```bash
$ filterBAM lca --help
usage: filterBAM lca [-h] --bam BAM [-p STR] [-r FILE] [-t INT] [--names FILE] [--nodes FILE] [--acc2taxid FILE] [--lca-rank STR] [--lca-summary [FILE]] [--scale STR] [-m STR] [--custom]
                     [--stats [FILE]]

optional arguments:
  -h, --help            show this help message and exit
  -p STR, --prefix STR  Prefix used for the output files (default: None)
  -r FILE, --reference-lengths FILE
                        File with references lengths (default: None)
  -t INT, --threads INT
                        Number of threads to use (default: 1)

required arguments:
  --bam BAM             BAM file containing aligned reads (default: None)

LCA optional arguments:
  --names FILE          Names dmp file from taxonomy (default: None)
  --nodes FILE          Nodes dmp file from taxonomy (default: None)
  --acc2taxid FILE      acc2taxid file from taxonomy (default: None)
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
filterBAM lca --bam c55d4e2df1.dedup.filtered.bam --names ./taxonomy/names.dmp --nodes ./taxonomy/nodes.dmp --acc2taxid ./taxonomy/acc2taxid.map.gz --threads 10 --lca-rank genus
```

**--names**: Names dmp file from taxonomy 

**--nodes**: Nodes dmp file from taxonomy 

**--acc2taxid**: acc2taxid file from taxonomy 

**--rank-lca**: Rank to use for LCA calculation 

**--scale**: Scale taxonomic abundance by this factor; suffix K/M recognized 

## How the reassignment process works

The read reassignment algorithm aims to resolve multi-mapping reads by iteratively refining alignment probabilities between reads and reference sequences. It begins by calculating an initial score $S$ for each alignment:

$$S = r_m M - p_m X - p_o G - p_e E$$

where $M$, $X$, $G$, and $E$ represent the number of matches, mismatches, gap openings, and gap extensions respectively. The terms $r_m$, $p_m$, $p_o$, and $p_e$ are the corresponding rewards or penalties.

These raw scores are then normalized in two steps. First, they are shifted to ensure they are positive:

$$S' = S - \min(S) + 1$$

Then, they are divided by the alignment length $L$:

$$S'' = \frac{S'}{L}$$

Using these normalized scores, we initialize the probability $P(r_i|g_j)$ of read $r_i$ originating from reference sequence $g_j$:

$$P(r_i|g_j) = \frac{S''\_{ij}}{\sum_{k} S''\_{ik}}$$

The algorithm then enters an iterative refinement phase. In each iteration, it first calculates subject weights:

$$W_j = \sum_i P(r_i|g_j)$$

$$W'_j = \frac{W_j}{L_j}$$

where $L_j$ is the length of reference sequence $g_j$. These weights are used to update the probabilities:

$$P'(r_i|g_j) = P(r_i|g_j) \cdot W'_j$$

$$P''(r_i|g_j) = \frac{P'(r_i|g_j)}{\sum_k P'(r_i|g_k)}$$

Following this update, the algorithm filters alignments. For each read $r_i$, it calculates the maximum probability across all references:

$$P_{max}(r_i) = \max_j P''(r_i|g_j)$$

It then applies a scaling factor $\alpha$ (a value < 1) and retains only alignments satisfying:

$$P''(r_i|g_j) \geq \alpha \cdot P_{max}(r_i)$$

This iterative process continues until no more alignments are removed or a maximum number of iterations is reached. The final read-to-reference assignments are implicit in the filtering process, with only the highest probability alignments retained.

### Applications and recommendations

One of the main applications of **bam-filter** is to reliably identify which potential organisms are present in a metagenomic ancient sample, and get relatively accurate taxonomic abundances, even when they are present in very low abundances. The resulting BAM file then can be used as input for [metaDMG](https://github.com/metaDMG-dev/metaDMG-cpp). We rely on several measures to discriminate between noise and a potential signal, analyzing the mapping results at two different levels:

- Is the observed breadth aligned with the expected one?
- Are the reads spread evenly across the reference or they are clumped in a few regions?

To assess the first question we use the concepts defined [here](https://doi.org/10.1016/0888-7543(88)90007-9). We estimate the ratio between the observed and expected breadth as a function of the coverage. If we get a **breadth_exp_ratio** close to 1, it means that the coverage we observed is close to the one we expect based on the calculated coverage. While this measure is already a strong indicator of a potential signal, we complement it with the metrics that measure the **normalized positional entropy** and the **normalized distribution inequality** (Gini coefficient) of the positions in the coverage. For details on how are calculated check [here](https://www.frontiersin.org/articles/10.3389/fmicb.2022.918015/full). These two metrics will help to identify those cases where we get a high **breadth_exp_ratio** but the coverage is not evenly distributed across the reference but instead is clumped in a few regions. One thing to be aware of is that we need to bin the reference to calculate those metrics. In our case, we use the ability of [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html) to identify the [numbers of bins](https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges), either using the Sturges or the Freedman-Diaconis rule. Finally, we use the [knee point detection algorithm](https://github.com/arvkevi/kneed) to identify the optimal values where to filter the Gini coefficient as a function of the positional entropy.

## LCA genome normalized abundances

The taxonomic abundances for each rank are estimated by normalizing the number of reads by the length of the reference. The LCA approach ranks taxonomic paths based on likelihood, selecting the most probable reference. The program also uses the TAD (Truncated Average Depth) estimated reads for the LCA analysis, which helps mitigate the effect of uneven coverages across references.