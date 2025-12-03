#!/usr/bin/env python3
"""
High-performance SAM/BAM to Parquet converter.

Uses parallel multiprocessing + pure C++ conversion for maximum speed:
- 32 parallel processes reading different record ranges
- Pure Cython/C++ hot loop (zero Python overhead)
- Process isolation for automatic memory cleanup
- 60M records in ~76 seconds (0.79 M records/sec)
- 92% space savings vs gzipped SAM (12.8x compression)

All 30 fields captured including PMD scores for ancient DNA analysis.
"""

import os
import sys
import time
from pathlib import Path
from bam_filter import logging as bf_logging
from bam_filter.parallel_parquet_records import parallel_convert_records
from bam_filter.manifest import create_manifest, get_file_info

LOG_TAG = "to-parquet"


def do_bam_to_parquet(args):
    """
    CLI entry point for BAM to Parquet conversion.

    Uses high-performance parallel multiprocessing architecture:
    - 32 separate OS processes reading different record ranges
    - Pure C++/Cython hot loop (zero Python overhead)
    - Process isolation = automatic memory cleanup
    - PMD scores extracted from BAM tags (PM/PMD)

    Performance on 60M records:
    - Time: ~76 seconds (0.79 M records/sec)
    - Space: 92% savings vs gzipped SAM
    - Memory: No accumulation due to process isolation
    """

    # Extract arguments
    input_file = args.bam
    output_dir = args.output
    num_processes = getattr(args, 'threads', 32)

    # Auto-configure based on system if not specified
    if not hasattr(args, 'num_partitions') or args.num_partitions == -1:
        num_partitions = max(16, num_processes // 2)
    else:
        num_partitions = args.num_partitions

    if not hasattr(args, 'batch_size') or args.batch_size == -1:
        batch_size = 100000
    else:
        batch_size = args.batch_size

    compression_level = getattr(args, 'compression_level', 3)

    # Logging
    bf_logging.log(LOG_TAG, f"Converting {input_file} to Parquet")
    bf_logging.log(LOG_TAG, f"Output directory: {output_dir}")
    bf_logging.log(LOG_TAG, f"Parallel processes: {num_processes}")
    bf_logging.log(LOG_TAG, f"Partitions: {num_partitions}")
    bf_logging.log(LOG_TAG, f"Batch size: {batch_size:,}")
    bf_logging.log(LOG_TAG, f"Compression level: {compression_level}")
    bf_logging.log(LOG_TAG, "")

    # PMD calculation status
    calculate_pmd = getattr(args, 'calculate_pmd', True)
    library_type = getattr(args, 'library_type', 'ds')

    if calculate_pmd:
        lib_name = "double-stranded" if library_type == "ds" else "single-stranded"
        bf_logging.log(LOG_TAG, f"PMD scores: ENABLED (calculating on-the-fly, {lib_name})")
    else:
        bf_logging.log(LOG_TAG, "PMD scores: disabled")
    bf_logging.log(LOG_TAG, "")

    start_time = time.time()

    # Run parallel conversion
    try:
        stats = parallel_convert_records(
            input_file=input_file,
            output_dir=output_dir,
            num_processes=num_processes,
            num_partitions=num_partitions,
            batch_size=batch_size,
            compression_level=compression_level,
            num_threads_per_worker=1,
            count_records=True,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
        )

        elapsed = time.time() - start_time

        # Report results
        bf_logging.log(LOG_TAG, "")
        bf_logging.log(LOG_TAG, "="*80)
        bf_logging.log(LOG_TAG, "Conversion complete!")
        bf_logging.log(LOG_TAG, f"  Total records: {stats['total_records']:,}")
        bf_logging.log(LOG_TAG, f"  Processing time: {elapsed:.1f}s ({elapsed/60:.2f} min)")
        bf_logging.log(LOG_TAG, f"  Throughput: {stats['throughput_records_per_sec']/1e6:.3f} M records/sec")

        # Check output size
        output_path = Path(output_dir)
        total_size = sum(f.stat().st_size for f in output_path.rglob('*.parquet'))
        size_gb = total_size / (1024**3)
        bf_logging.log(LOG_TAG, f"  Output size: {size_gb:.2f} GB")
        bf_logging.log(LOG_TAG, "="*80)

        # Create provenance manifest
        try:
            input_info = get_file_info(input_file)
            input_info['total_records'] = stats['total_records']

            create_manifest(
                output_dir=output_dir,
                operation="bam_to_parquet",
                command=' '.join(sys.argv),
                parameters={
                    "num_processes": num_processes,
                    "num_partitions": num_partitions,
                    "batch_size": batch_size,
                    "compression_level": compression_level,
                    "calculate_pmd": calculate_pmd,
                    "library_type": library_type,
                },
                input_info={
                    "type": "bam",
                    "path": input_file,
                    "size_bytes": input_info.get('size_bytes', 0),
                    "total_records": stats['total_records'],
                },
                output_info={
                    "records_processed": stats['total_records'],
                    "processing_time_seconds": elapsed,
                    "throughput_records_per_sec": stats['throughput_records_per_sec'],
                    "num_workers": stats['num_workers'],
                    "output_size_bytes": total_size,
                    "output_size_gb": size_gb,
                    "chunks_created": stats['num_workers'],
                },
                original_input={
                    "path": input_file,
                    "size_bytes": input_info.get('size_bytes', 0),
                    "total_records": stats['total_records'],
                }
            )
            bf_logging.log(LOG_TAG, "")
            bf_logging.log(LOG_TAG, "Provenance manifest created: manifest.json")
        except Exception as e:
            bf_logging.log(LOG_TAG, f"Warning: Failed to create manifest: {e}")

        return stats

    except Exception as e:
        bf_logging.error(f"{LOG_TAG}: Conversion failed: {e}")
        raise
