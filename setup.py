from setuptools import setup
import versioneer

from setuptools import Extension
from Cython.Build import cythonize
import numpy
import os

conda_prefix = os.environ.get("CONDA_PREFIX", "")

# Common compilation settings for all extensions
common_compile_args = ["-fopenmp", "-O3", "-ffast-math", "-funroll-loops"]
common_link_args = ["-fopenmp", "-lmvec", "-lm"]
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
        "bam_filter.generic_filters",
        ["bam_filter/generic_filters.pyx"],
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
    # Taxonomy-informed filtering module (combines graph + taxonomy for better filtering)
    Extension(
        "bam_filter.processor_taxonomy_filters",
        ["bam_filter/processor_taxonomy_filters.pyx"],
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
    Extension(
        "bam_filter.processor_tiered_filters",
        ["bam_filter/processor_tiered_filters.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
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
        "bam_filter.processor_community_igraph",
        ["bam_filter/processor_community_igraph.pyx"],
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
    # PMD (Post-Mortem Damage) learning and damage-corrected ANI computation
    # Needs libm and libmvec for vectorized math functions (exp, log) with -ffast-math
    Extension(
        "bam_filter.processor_pmd",
        ["bam_filter/processor_pmd.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args + ["-lm", "-lmvec"],
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["m"],
    ),
    # Bayesian damage model for ancient/modern reference classification
    Extension(
        "bam_filter.processor_damage_model",
        ["bam_filter/processor_damage_model.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args + ["-lm"],
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["m"],
    ),
    Extension(
        "bam_filter.processor_lca",
        ["bam_filter/processor_lca.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # LCA Stats - per-taxid statistics aggregation
    Extension(
        "bam_filter.processor_lca_stats",
        ["bam_filter/processor_lca_stats.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Probabilistic taxonomic profiler - Bayesian hierarchical model for ancient DNA
    Extension(
        "bam_filter.processor_prob_profile",
        ["bam_filter/processor_prob_profile.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Novel probabilistic profiler - DCMS + belief propagation (no LCA)
    Extension(
        "bam_filter.probabilistic_profiler",
        ["bam_filter/probabilistic_profiler.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Generic SQUAREM acceleration for EM algorithms
    Extension(
        "bam_filter.squarem",
        ["bam_filter/squarem.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args + ["-lm"],
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["m"],
    ),
    # Unified ancient DNA damage model - single model for ANI correction, classification, profiling
    Extension(
        "bam_filter.unified_damage",
        ["bam_filter/unified_damage.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args + ["-lm"],
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["m"],
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
    # Unified TSV writer with optional compression (replaces fast_gzip_writer + fast_tsv_writer)
    Extension(
        "bam_filter.tsv_writer",
        ["bam_filter/tsv_writer.pyx"],
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
    # Network QC metrics for evaluating EM quality based on graph structure
    Extension(
        "bam_filter.processor_network_qc",
        ["bam_filter/processor_network_qc.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries,
    ),
    # Ultra-fast taxonomy database for LCA queries and reference taxonomic assignment
    # Now with Arrow C++ support for ultra-fast Parquet loading (50-100x speedup)
    Extension(
        "bam_filter.taxonomy_db",
        ["bam_filter/taxonomy_db.pyx"],
        extra_compile_args=common_compile_args
        + ["-std=c++17"],  # C++17 required for Arrow
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries
        + ["arrow", "parquet", "duckdb"],  # Add Arrow, Parquet, and DuckDB
        language="c++",  # Use C++ compiler for Arrow API
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # BAM to Parquet converter for metagenomic-scale data analysis with DuckDB
    Extension(
        "bam_filter.processor_parquet_writer",
        ["bam_filter/processor_parquet_writer.pyx"],
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
    # DuckDB-based Parquet converter (SAM/SAM.gz/BAM support, dual-table design)
    Extension(
        "bam_filter.parquet_converter_duckdb",
        ["bam_filter/parquet_converter_duckdb.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["duckdb"],
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # Streaming Parquet converter (parallel, batched writes, individual tag columns)
    # Temporarily disabled - has syntax errors
    # Extension(
    #     "bam_filter.parquet_converter_streaming",
    #     ["bam_filter/parquet_converter_streaming.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries + ["duckdb"],
    #     define_macros=[
    #         ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    #         ("_GNU_SOURCE", None),
    #     ],
    # ),
    # Batched Parquet converter (smart buffering, single-pass, columnar batches)
    Extension(
        "bam_filter.parquet_converter_batched",
        ["bam_filter/parquet_converter_batched.pyx"],
        extra_compile_args=common_compile_args,
        extra_link_args=common_link_args,
        include_dirs=common_include_dirs,
        library_dirs=common_library_dirs,
        libraries=common_libraries + ["duckdb"],
        define_macros=[
            ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
            ("_GNU_SOURCE", None),
        ],
    ),
    # # Arrow C++ Parquet converter (fastest: HTSlib + Arrow C++ API directly, no DuckDB)
    # Extension(
    #     "bam_filter.parquet_converter_arrow_cpp",
    #     sources=[
    #         "bam_filter/parquet_converter_arrow_cpp.pyx",
    #         "bam_filter/arrow_parquet_writer.cpp",
    #     ],
    #     extra_compile_args=common_compile_args + ["-std=c++17"],
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries + ["arrow", "parquet"],
    #     language="c++",
    #     define_macros=[
    #         ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    #         ("_GNU_SOURCE", None),
    #     ],
    # ),
    # # PURE C++ Parquet converter (ZERO Python overhead - maximum performance)
    # Extension(
    #     "bam_filter.parquet_converter_pure_cpp",
    #     sources=[
    #         "bam_filter/parquet_converter_pure_cpp.pyx",
    #         "bam_filter/arrow_parquet_writer.cpp",
    #     ],
    #     extra_compile_args=common_compile_args + ["-std=c++17", "-fopenmp"],
    #     extra_link_args=common_link_args + ["-fopenmp"],
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries + ["arrow", "parquet"],
    #     language="c++",
    #     define_macros=[
    #         ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    #         ("_GNU_SOURCE", None),
    #     ],
    # ),
    # # MINIMAL schema Parquet converter (only 16 fields - MAXIMUM SPEED)
    # Extension(
    #     "bam_filter.parquet_converter_minimal",
    #     sources=[
    #         "bam_filter/parquet_converter_minimal.pyx",
    #         "bam_filter/arrow_parquet_writer_minimal.cpp",
    #     ],
    #     extra_compile_args=common_compile_args + ["-std=c++17"],
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries + ["arrow", "parquet"],
    #     language="c++",
    #     define_macros=[
    #         ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    #         ("_GNU_SOURCE", None),
    #     ],
    # ),
    # # Pure C/OpenMP Parquet converter (TRUE parallel threading, shared memory)
    # Extension(
    #     "bam_filter.parquet_converter_openmp",
    #     sources=[
    #         "bam_filter/parquet_converter_openmp.pyx",
    #         "bam_filter/arrow_parquet_writer.cpp",
    #     ],
    #     extra_compile_args=common_compile_args + ["-std=c++17", "-fopenmp"],
    #     extra_link_args=common_link_args + ["-fopenmp"],
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries + ["arrow", "parquet"],
    #     language="c++",
    #     define_macros=[
    #         ("NPY_NO_DEPRECATED_API", "NPY_1_7_API_VERSION"),
    #         ("_GNU_SOURCE", None),
    #     ],
    # ),
    # # Fast unified alignment reader (SAM/SAM.gz/BAM) - replaced by DuckDB converter
    # Extension(
    #     "bam_filter.fast_alignment_reader",
    #     ["bam_filter/fast_alignment_reader.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # # Fast dual-table Parquet writer (by_reference + by_read) - replaced by DuckDB converter
    # Extension(
    #     "bam_filter.fast_parquet_writer",
    #     ["bam_filter/fast_parquet_writer.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
    # # Fast BAM/SAM to Parquet conversion pipeline - replaced by DuckDB converter
    # Extension(
    #     "bam_filter.bam_to_parquet_fast",
    #     ["bam_filter/bam_to_parquet_fast.pyx"],
    #     extra_compile_args=common_compile_args,
    #     extra_link_args=common_link_args,
    #     include_dirs=common_include_dirs,
    #     library_dirs=common_library_dirs,
    #     libraries=common_libraries,
    # ),
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
    cmdclass=versioneer.get_cmdclass(),
    ext_modules=cythonize(ext_modules, language_level="3"),
)
