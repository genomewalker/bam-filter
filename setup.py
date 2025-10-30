from setuptools import setup
import versioneer

from setuptools import Extension
from Cython.Build import cythonize
import numpy
import os

requirements = [
    "Cython>=0.29.24",
    "pandas>=1.3.3",
    "scipy>=1.9.0",
    "tqdm>=4.64.1",
    "pysam>=0.17.0",
    "numpy>1.24",
    "pyrle>=0.0.31",
    "sorted-nearest>=0.0.38",
    "pyranges>=0.0.118",
    "kneed>=0.8.1",
    "matplotlib>=3.6.0",
    "taxopy>=0.12.0",
    "python-datatable>=1.1.3",
    "networkx>=3.2.1",
    "psutil>=5.9.8",
    "numba>=0.60.0",
    "pyarrow>=10.0.0",  # For Parquet support in taxonomy module
]

conda_prefix = os.environ.get("CONDA_PREFIX", "")

# Common compilation settings for all extensions
common_compile_args = ["-fopenmp", "-O3", "-ffast-math", "-funroll-loops"]
common_link_args = ["-fopenmp"]
common_include_dirs = [
    numpy.get_include(),
    os.path.join(conda_prefix, "include"),
    os.path.join(os.path.dirname(__file__), "bam_filter", "include"),
    "bam_filter/include",  # Ensure relative path is always present
]
common_library_dirs = [os.path.join(conda_prefix, "lib")]
common_libraries = ["hts", "z", "igraph"]
common_include_dirs.append(os.path.join(conda_prefix, "include", "igraph"))

# Define all Cython extensions - both legacy and modular
ext_modules = [
    # Legacy monolithic module (for backward compatibility)
    # Extension(
    #     "bam_filter.bam_cython",
    #     ["bam_filter/bam_cython.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # BAM Writer module for efficient writing
    # Extension(
    #     "bam_filter.bam_writer",
    #     ["bam_filter/bam_writer.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # BAM Statistics module for comprehensive per-reference statistics (monolithic - backward compatibility)
    # Disabled temporarily due to compilation errors - use modular components instead
    Extension(
        "bam_filter.stats",
        ["bam_filter/stats.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.stats_helpers",
        ["bam_filter/stats_helpers.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # RLE coverage helpers (modularized from stats.pyx)
    Extension(
        "bam_filter.stats_rle",
        ["bam_filter/stats_rle.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # RLE coverage helpers (modularized from stats.pyx)
    Extension(
        "bam_filter.processor_graph",
        ["bam_filter/processor_graph.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Graph operations module (extracted from processor_leiden for better separation)
    Extension(
        "bam_filter.processor_graph_ops",
        ["bam_filter/processor_graph_ops.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # Taxonomy-aware graph analysis module
    Extension(
        "bam_filter.processor_graph_taxonomy",
        ["bam_filter/processor_graph_taxonomy.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # # Lightweight nogil logging/timing utilities (available to cimport)
    # Extension(
    #     "bam_filter.nogil_log",
    #     ["bam_filter/nogil_log.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # EM Core module
    # Extension(
    #     "bam_filter.em_core",
    #     ["bam_filter/em_core.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # # Fast Kernels module
    # Extension(
    #     "bam_filter.fast_kernels",
    #     ["bam_filter/fast_kernels.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # # Anderson Acceleration module
    # Extension(
    #     "bam_filter.anderson_accel",
    #     ["bam_filter/anderson_accel.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # Unified high-performance BAM processor (NEWEST - replaces reassign.py pipeline)
    Extension(
        "bam_filter.processor",
        ["bam_filter/processor.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),  # Enable GNU extensions
        ],
    ),
    # Dedicated batch-processing module (refactored out of processor.pyx)
    Extension(
        "bam_filter.processor_batch",
        ["bam_filter/processor_batch.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.processor_memory",
        ["bam_filter/processor_memory.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    Extension(
        "bam_filter.processor_fast_math",
        ["bam_filter/processor_fast_math.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.processor_precomputed",
        ["bam_filter/processor_precomputed.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.processor_filters",
        ["bam_filter/processor_filters.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # TSV writer for graph analysis (separated from processor_graph)
    Extension(
        "bam_filter.processor_graph_tsv",
        ["bam_filter/processor_graph_tsv.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # GraphML export for graph visualization (Cytoscape, igraph, etc.)
    Extension(
        "bam_filter.processor_graph_export",
        ["bam_filter/processor_graph_export.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # Fast Leiden clustering using igraph C library (10-100x speedup)
    Extension(
        "bam_filter.processor_leiden_igraph",
        ["bam_filter/processor_leiden_igraph.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # Mapping helpers extracted from processor.pyx (keeps processor.pyx smaller)
    Extension(
        "bam_filter.processor_mapping",
        ["bam_filter/processor_mapping.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.processor_md_quality",
        ["bam_filter/processor_md_quality.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_compile_dirs if False else common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Small Python-facing wrapper around the internal PMD scoring functions
    Extension(
        "bam_filter._pmd_wrapper",
        ["bam_filter/_pmd_wrapper.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Tiny helper: read-name hashing (extracted from processor.pyx)
    Extension(
        "bam_filter.processor_hash",
        ["bam_filter/processor_hash.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # BAM writer helpers (moved out of processor.pyx into a modular file)
    Extension(
        "bam_filter.processor_bam_writer",
        ["bam_filter/processor_bam_writer.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Sorting helpers module extracted from processor.pyx
    Extension(
        "bam_filter.processor_sort",
        ["bam_filter/processor_sort.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Convergence helper functions extracted from processor.pyx
    Extension(
        "bam_filter.processor_convergence_helpers",
        ["bam_filter/processor_convergence_helpers.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.processor_em",
        ["bam_filter/processor_em.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Processing statistics tracker
    Extension(
        "bam_filter.processor_stats",
        ["bam_filter/processor_stats.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.batch_utils",
        ["bam_filter/batch_utils.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.reference_lengths",
        ["bam_filter/reference_lengths.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.stats_io",
        ["bam_filter/stats_io.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    Extension(
        "bam_filter.stats_bam_writer",
        ["bam_filter/stats_bam_writer.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Outlier detection module for multi-metric graph filtering
    Extension(
        "bam_filter.processor_outlier_detection",
        ["bam_filter/processor_outlier_detection.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # Ultra-fast taxonomy database for LCA queries and reference taxonomic assignment
    # Now with Arrow C++ support for ultra-fast Parquet loading (50-100x speedup)
    Extension(
        "bam_filter.taxonomy_db",
        ["bam_filter/taxonomy_db.pyx"],
        extra_compile_args=common_compile_args + ["-std=c++17"],  # C++17 required for Arrow
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["arrow", "parquet", "duckdb"],  # Add Arrow, Parquet, and DuckDB
        language="c++",  # Use C++ compiler for Arrow API
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # # Modular architecture - Core batch processing framework
    # Extension(
    #     "bam_filter.batch_processor",
    #     ["bam_filter/batch_processor.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # # Modular architecture - Alignment score processor
    # # Extension(
    # #     "bam_filter.score_processor",
    # #     ["bam_filter/score_processor.pyx"],
    # #     extra_compile_args=common_compile_args,
    # #     extra_link_args=common_link_args,
    # #     include_dirs=common_include_dirs,
    # #     library_dirs=common_library_dirs,
    # #     libraries=common_libraries,
    # # ),
    # # Modular architecture - Statistics processor
    # Extension(
    #     "bam_filter.stats_processor",
    #     ["bam_filter/stats_processor.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
]

setup(
    setup_requires=[
        # Setuptools 18.0 properly handles Cython extensions.
        "setuptools>=18.0",
        "Cython>=0.29.24",
    ],
    name="bam-filter",
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description="A simple tool to filter references from a BAM file using different filter types",
    license="GNUv3",
    author="Antonio Fernandez-Guerra",
    author_email="antonio@metagenomics.eu",
    url="https://github.com/genomewalker/bam-filter",
    packages=["bam_filter"],
    entry_points={"console_scripts": ["filterBAM=bam_filter.__main__:main"]},
    install_requires=requirements,
    ext_modules=cythonize(ext_modules),
    keywords="bam-filter",
    classifiers=[
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
    ],
)
