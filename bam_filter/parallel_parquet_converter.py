#!/usr/bin/env python3
"""
Multi-process parallel SAM/BAM to Parquet converter.

Strategy: Split references among worker processes. Each worker:
1. Reads the full SAM/BAM file sequentially (bgzip-compressed)
2. Only processes records from assigned references
3. Writes to partition files

This achieves parallelization without needing random access or indices.
"""

import sys
import time
import multiprocessing as mp
from pathlib import Path
from typing import List, Dict, Set
import subprocess

# Import the pure C++ converter
from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp


def get_references_from_header(sam_file: str) -> List[str]:
    """Extract all reference names from SAM/BAM header using samtools."""
    result = subprocess.run(
        ['samtools', 'view', '-H', sam_file],
        capture_output=True,
        text=True,
        check=True
    )

    refs = []
    for line in result.stdout.splitlines():
        if line.startswith('@SQ'):
            # @SQ	SN:ref_name	LN:length
            for field in line.split('\t'):
                if field.startswith('SN:'):
                    refs.append(field[3:])
                    break

    return refs


def worker_process(
    worker_id: int,
    input_file: str,
    output_dir: str,
    assigned_refs: Set[str],
    num_partitions: int,
    batch_size: int,
    compression_level: int,
    num_bgzip_threads: int,
) -> Dict:
    """
    Worker process that reads the SAM/BAM file and processes only assigned references.

    Each worker:
    - Opens the SAM/BAM file independently
    - Reads all records sequentially
    - Only processes records from assigned_refs
    - Writes to worker-specific output directory
    """
    import pysam
    from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

    # Create worker-specific output directory
    worker_dir = Path(output_dir) / f"worker_{worker_id:02d}"
    worker_dir.mkdir(parents=True, exist_ok=True)

    print(f"[Worker {worker_id}] Processing {len(assigned_refs)} references: {sorted(list(assigned_refs))[:5]}...")

    start_time = time.time()

    # Open SAM/BAM file
    sam = pysam.AlignmentFile(input_file, "r", threads=num_bgzip_threads)

    # Create reference name to ID mapping
    ref_to_tid = {name: tid for tid, name in enumerate(sam.references)}
    assigned_tids = {ref_to_tid[ref] for ref in assigned_refs if ref in ref_to_tid}

    # We'll use a modified approach: write a filtered SAM to temp file, then convert
    # This is simpler than modifying the C++ converter to skip records

    # Actually, better approach: just call the C++ converter with a filter
    # But that requires modifying the Cython code...

    # Simplest approach: Create a temporary filtered SAM file for this worker
    temp_sam = worker_dir / "temp_filtered.sam"
    temp_sam_gz = worker_dir / "temp_filtered.sam.gz"

    # Write header + filtered records
    records_written = 0
    with pysam.AlignmentFile(str(temp_sam), "w", header=sam.header) as out:
        for read in sam:
            if read.reference_id in assigned_tids:
                out.write(read)
                records_written += 1

                if records_written % 100000 == 0:
                    print(f"[Worker {worker_id}] Filtered {records_written:,} records...")

    sam.close()

    print(f"[Worker {worker_id}] Filtered {records_written:,} records in {time.time()-start_time:.1f}s")

    # Compress the filtered SAM
    subprocess.run(['bgzip', '-@', str(num_bgzip_threads), str(temp_sam)], check=True)

    # Now convert to Parquet using the fast C++ converter
    print(f"[Worker {worker_id}] Converting to Parquet...")
    convert_start = time.time()

    stats = convert_sam_bam_to_parquet_pure_cpp(
        str(temp_sam_gz),
        str(worker_dir),
        num_partitions=num_partitions,
        batch_size=batch_size,
        write_by_reference=True,
        write_by_read=True,
        compression_level=compression_level,
        num_threads=num_bgzip_threads,
    )

    convert_time = time.time() - convert_start
    total_time = time.time() - start_time

    # Clean up temp file
    temp_sam_gz.unlink()

    print(f"[Worker {worker_id}] DONE: {records_written:,} records in {total_time:.1f}s total " +
          f"({convert_time:.1f}s convert, {records_written/convert_time:.0f} rec/s)")

    return {
        'worker_id': worker_id,
        'records_processed': records_written,
        'total_time': total_time,
        'convert_time': convert_time,
    }


def parallel_convert_sam_to_parquet(
    input_file: str,
    output_dir: str,
    num_workers: int = 32,
    num_partitions: int = 128,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_bgzip_threads_per_worker: int = 2,
):
    """
    Convert SAM/BAM to Parquet using multiple worker processes in parallel.

    Args:
        input_file: Path to SAM/BAM/SAM.gz file
        output_dir: Output directory for Parquet files
        num_workers: Number of parallel worker processes
        num_partitions: Number of partition files per worker
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_bgzip_threads_per_worker: BGzip threads per worker

    Returns:
        dict with statistics
    """
    overall_start = time.time()

    print(f"Parallel SAM/BAM to Parquet converter")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Workers: {num_workers}")
    print(f"Partitions per worker: {num_partitions}")
    print(f"BGzip threads per worker: {num_bgzip_threads_per_worker}")
    print()

    # Get all references
    print("Reading reference list from header...")
    refs = get_references_from_header(input_file)
    print(f"Found {len(refs)} references")

    if len(refs) == 0:
        raise ValueError("No references found in header!")

    # Split references among workers
    refs_per_worker = len(refs) // num_workers
    remainder = len(refs) % num_workers

    worker_assignments = []
    ref_idx = 0
    for worker_id in range(num_workers):
        # Give extra refs to first 'remainder' workers
        count = refs_per_worker + (1 if worker_id < remainder else 0)
        if count > 0:
            assigned = set(refs[ref_idx:ref_idx + count])
            worker_assignments.append((worker_id, assigned))
            ref_idx += count

    print(f"Split {len(refs)} references among {len(worker_assignments)} workers")
    for wid, assigned in worker_assignments[:5]:
        print(f"  Worker {wid}: {len(assigned)} refs")
    if len(worker_assignments) > 5:
        print(f"  ... and {len(worker_assignments) - 5} more workers")
    print()

    # Launch worker processes
    print(f"Launching {len(worker_assignments)} worker processes...")
    with mp.Pool(processes=len(worker_assignments)) as pool:
        results = pool.starmap(
            worker_process,
            [
                (
                    worker_id,
                    input_file,
                    output_dir,
                    assigned_refs,
                    num_partitions,
                    batch_size,
                    compression_level,
                    num_bgzip_threads_per_worker,
                )
                for worker_id, assigned_refs in worker_assignments
            ]
        )

    total_time = time.time() - overall_start
    total_records = sum(r['records_processed'] for r in results)

    print()
    print("="*80)
    print("PARALLEL CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.1f} minutes)")
    print(f"Throughput: {total_records/total_time:.0f} records/sec")
    print(f"Throughput: {total_records/total_time/1e6:.2f} M records/sec")
    print()
    print(f"Output in: {output_dir}/")
    print(f"  {len(worker_assignments)} worker directories with {num_partitions} partitions each")

    return {
        'total_records': total_records,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time,
        'num_workers': len(worker_assignments),
        'worker_results': results,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python parallel_parquet_converter.py <input.sam.gz> <output_dir> [num_workers]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_workers = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    stats = parallel_convert_sam_to_parquet(
        input_file,
        output_dir,
        num_workers=num_workers,
        num_partitions=16,  # Fewer partitions per worker since we have many workers
        batch_size=100000,
        compression_level=3,
        num_bgzip_threads_per_worker=2,  # 2 threads * 32 workers = 64 threads total
    )
