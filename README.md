
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

usage: filterBAM [-h] [-t THREADS] [-p PREFIX] [-l MIN_READ_LENGTH] [-n MIN_READ_COUNT] [-b MIN_EXPECTED_BREADTH_RATIO] [-B MIN_BREADTH] [-a MIN_READ_ANI] [-c MIN_COVERAGE_EVENNESS] [-m SORT_MEMORY] [-N] [--scale SCALE] [-r REFERENCE_LENGTHS]
                 [--read-length-freqs] [--only-stats] [--only-stats-filtered] [--debug] [--version]
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
  -l MIN_READ_LENGTH, --min-read-length MIN_READ_LENGTH
                        Minimum read length (default: 30)
  -n MIN_READ_COUNT, --min-read-count MIN_READ_COUNT
                        Minimum read count (default: 10)
  -b MIN_EXPECTED_BREADTH_RATIO, --min-expected-breadth-ratio MIN_EXPECTED_BREADTH_RATIO
                        Minimum expected breadth ratio (default: 0.5)
  -B MIN_BREADTH, --min-breadth MIN_BREADTH
                        Minimum breadth (default: 0)
  -a MIN_READ_ANI, --min-read-ani MIN_READ_ANI
                        Minimum average read ANI (default: 90.0)
  -c MIN_COVERAGE_EVENNESS, --min-coverage-evenness MIN_COVERAGE_EVENNESS
                        Minimum coverage evenness (default: 0)
  -m SORT_MEMORY, --sort-memory SORT_MEMORY
                        Set maximum memory per thread for sorting; suffix K/M/G recognized (default: 1G)
  -N, --sort-by-name          Sort by read names (default: False)
  --scale SCALE         Scale taxonomic abundance by this factor; suffix K/M recognized (default: 1000000.0)
  -r REFERENCE_LENGTHS, --reference-lengths REFERENCE_LENGTHS
                        File with references lengths (default: None)
  --read-length-freqs   Save a JSON file with the read length frequencies mapped to each reference (default: False)
  --only-stats          Only produce statistics and skip filtering (default: False)
  --only-stats-filtered
                        Only filter statistics and skip BAM filtering (default: False)
  --debug               Print debug messages (default: False)
  --version             Print program version
```

One would run filterBAM as:

```bash
filterBAM --min-read-count 100 --min-expected-breadth-ratio 0.75 --min-read-ani 98 --sort-by-name --sort-memory 1G --reference-lengths gtdb-r202.len.map --threads 16  c55d4e2df1.dedup.bam 
```

**--min-read-count**: Minimum number of reads mapped to a reference in the BAM file

**--min-expected-breadth-ratio**: Minimum expected breadth ratio needed to keep a reference. This is based on the concepts defined [here](https://instrain.readthedocs.io/en/latest/important_concepts.html#detecting-organisms-in-metagenomic-data). It basically estimates the ratio between the observed and expected breadth, the closest to 1 the more evenly distributed the mapped reads are and we can be more confident that the genome was detected.

**--min-read-ani**: Minimum average read ANI that a reference has

**--sort-by-name**: Sort filtered BAM file by read name so it can be used in metaDMG

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
    - **coverage_covered_mean**: Mean depth of the reference only counting covered bases
    - **reference_length**: Real reference length
    - **bam_reference_length**: Length reported by the BAM file
    - **breadth**: Breadth of coverage 
    - **exp_breadth**: Expected breadth of coverage. Using the formula: _expected_breadth = 1 - e<sup>-coverage</sup>_
    - **breadth_exp_ratio**: Ration between the obsrved and the expected depth of coverage
    - **c_v**: Coefficient of variation of the coverage
    - **cov_evenness**: Eveness of coverage as calculated [here](https://www.nature.com/articles/jhg201621).
    - **tax_abund_read**: Counts estimated like Woltka but using number of reads.
    - **tax_abund_aln**: Counts estimated like Woltka.
