
# filterBAM


[![GitHub release (latest by date including pre-releases)](https://img.shields.io/github/v/release/genomewalker/bam-filter?include_prereleases&label=version)](https://github.com/genomewalker/bam-filter/releases) [![bam-filter](https://github.com/genomewalker/bam-filter/workflows/filterBAM_ci/badge.svg)](https://github.com/genomewalker/bam-filter/actions) [![PyPI](https://img.shields.io/pypi/v/bam-filter)](https://pypi.org/project/bam-filter/) [![Conda](https://img.shields.io/conda/v/genomewalker/bam-filter)](https://anaconda.org/genomewalker/bam-filter)


A simple tool to calculate metrics from a BAM file and filter references with uneven coverages.

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

### Using conda

```bash
conda install -c conda-forge -c bioconda -c genomewalker bam-filter
```

### Install from source to use the development version

Using pip

```bash
pip install git+ssh://git@github.com/genomewalker/bam-filter.git
```

By cloning in a dedicated conda environment

```bash
git clone git@github.com:genomewalker/bam-filter.git
cd bam-filter
conda env create -f environment.yml
conda activate bam-filter
pip install -e .
```


# Usage

filterBAM only needs a BAM file. For a complete list of options:

```
$ filterBAM --help

usage: filterBAM [-h] [-t THREADS] [--reference-trim-length TRIM_ENDS] [--trim-min TRIM_MIN]
                 [--trim-max TRIM_MAX] [-p PREFIX] [-A MIN_READ_ANI] [-l MIN_READ_LENGTH]
                 [-n MIN_READ_COUNT] [-b MIN_EXPECTED_BREADTH_RATIO] [-e MIN_NORM_ENTROPY]
                 [-g MIN_NORM_GINI] [-B MIN_BREADTH] [-a MIN_AVG_READ_ANI] [-c MIN_COVERAGE_EVENNESS]
                 [-C MIN_COVERAGE_MEAN] [--include-low-detection] [-m SORT_MEMORY] [-N]
                 [--scale SCALE] [-r REFERENCE_LENGTHS] [--read-length-freqs] [--read-hits-count]
                 [--only-stats] [--plot] [--only-stats-filtered] [--chunk-size CHUNK_SIZE] [--debug]
                 [--version]
                 bam

A simple tool to calculate metrics from a BAM file and filter with uneven coverage.

positional arguments:
  bam                   BAM file containing aligned reads

optional arguments:
  -h, --help            show this help message and exit
  -t THREADS, --threads THREADS
                        Number of threads to use (default: 1)
  -p PREFIX, --prefix PREFIX
                        Prefix used for the output files (default: None)
  -m SORT_MEMORY, --sort-memory SORT_MEMORY
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default:
                        1G)
  -N, --sort-by-name    Sort by read names (default: False)
  --scale SCALE         Scale taxonomic abundance by this factor; suffix K/M recognized (default:
                        1000000.0)
  -r REFERENCE_LENGTHS, --reference-lengths REFERENCE_LENGTHS
                        File with references lengths (default: None)
  --read-length-freqs   Save a JSON file with the read length frequencies mapped to each reference
                        (default: False)
  --read-hits-count     Save a TSV file with the read hits frequencies mapped to each reference
                        (default: False)
  --only-stats          Only produce statistics and skip filtering (default: False)
  --plot                Plot genome coverage plots (default: False)
  --only-stats-filtered
                        Only filter statistics and skip BAM filtering (default: False)
  --chunk-size CHUNK_SIZE
                        Chunk size for parallel processing (default: None)
  --debug               Print debug messages (default: False)
  --version             Print program version

filtering arguments:
  -A MIN_READ_ANI, --min-read-ani MIN_READ_ANI
                        Minimum read ANI to keep a read (default: 90.0)
  -l MIN_READ_LENGTH, --min-read-length MIN_READ_LENGTH
                        Minimum read length (default: 30)
  -n MIN_READ_COUNT, --min-read-count MIN_READ_COUNT
                        Minimum read count (default: 10)
  -b MIN_EXPECTED_BREADTH_RATIO, --min-expected-breadth-ratio MIN_EXPECTED_BREADTH_RATIO
                        Minimum expected breadth ratio (default: 0.5)
  -e MIN_NORM_ENTROPY, --min-normalized-entropy MIN_NORM_ENTROPY
                        Minimum normalized entropy (default: auto)
  -g MIN_NORM_GINI, --min-normalized-gini MIN_NORM_GINI
                        Minimum normalized Gini coefficient (default: None)
  -B MIN_BREADTH, --min-breadth MIN_BREADTH
                        Minimum breadth (default: 0)
  -a MIN_AVG_READ_ANI, --min-avg-read-ani MIN_AVG_READ_ANI
                        Minimum average read ANI (default: 90.0)
  -c MIN_COVERAGE_EVENNESS, --min-coverage-evenness MIN_COVERAGE_EVENNESS
                        Minimum coverage evenness (default: 0)
  -C MIN_COVERAGE_MEAN, --min-coverage-mean MIN_COVERAGE_MEAN
                        Minimum coverage mean (default: 0)
  --include-low-detection
                        Include those references that fullfi all filtering criteria but the coverage
                        evenness is 0 (default: False)

miscellaneous arguments:
  --reference-trim-length TRIM_ENDS
                        Exclude n bases at the ends of the reference sequences (default: 0)
  --trim-min TRIM_MIN   Remove coverage that are below this percentile. Used for the Truncated Average
                        Depth (TAD) calculation (default: 10)
  --trim-max TRIM_MAX   Remove coverage that are above this percentile. Used for the Truncated Average
                        Depth (TAD) calculation (default: 90)
```

One would run filterBAM as:

```bash
filterBAM --min-read-count 100 --min-expected-breadth-ratio 0.75 --min-read-ani 98 --sort-by-name --sort-memory 1G --reference-lengths gtdb-r202.len.map --threads 16  c55d4e2df1.dedup.bam 
```

**--min-read-count**: Minimum number of reads mapped to a reference in the BAM file

**--min-expected-breadth-ratio**: Minimum expected breadth ratio needed to keep a reference.

**--min-read-ani**: Minimum average read ANI that a reference has

**--sort-by-name**: Sort filtered BAM file by read name so it can be used in [metaDMG](https://metadmg-dev.github.io/metaDMG-core/)

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
    - **exp_breadth**: Expected breadth of coverage. Using the formula: _expected_breadth = 1 - e<sup>-coverage</sup>_
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
    - **n_reads_tad**: Number of reads estimated in the TAD region using the formula *C = LN / G*, where C stands for the TAD coverage, N for the length of the TAD region and L for the average read length mapped to the reference.

> 

## Applications and recommendations

One of the main applications of **bam-filter** is to reliably identify which potential organisms are present in a metagenomic ancient sample, and get relatively accurate taxonomic abundances, even when they are present in very low abundances. The resulting BAM file then can be used as input for [metaDMG](https://metadmg-dev.github.io/metaDMG-core/). We rely on several measures to discriminate between noise and a potential signal, analyzing the mapping results at two different levels:

- Is the observed breadth aligned with the expected one?
- Are the reads spread evenly across the reference or they are clumped in a few regions?

To assess the first question we use the concepts defined [here](https://doi.org/10.1016/0888-7543(88)90007-9). We estimate the ratio between the observed and expected breadth as a function of the coverage. If we get a **breadth_exp_ratio** close to 1, it means that the coverage we observed is close to the one we expect based on the calculated coverage. While this measure is already a strong indicator of a potential signal, we complement it with the metrics that measure the **normalized positional entropy** and the **normalized distribution inequality** (Gini coefficient) of the positions in the coverage. For details on how are calculated check [here](https://www.frontiersin.org/articles/10.3389/fmicb.2022.918015/full). These two metrics will help to identify those cases where we get a high **breadth_exp_ratio** but the coverage is not evenly distributed across the reference but instead is clumped in a few regions. One thing to be aware of is that we need to bin the reference to calculate those metrics. In our case, we use the ability of [numpy.histogram](https://numpy.org/doc/stable/reference/generated/numpy.histogram.html) to identify the [numbers of bins](https://numpy.org/doc/stable/reference/generated/numpy.histogram_bin_edges.html#numpy.histogram_bin_edges), either using the Sturges or the Freedman-Diaconis rule. Finally, we use the [knee point detection algorithm](https://github.com/arvkevi/kneed) to identify the optimal values where to filter the Gini coefficient as a function of the positional entropy.