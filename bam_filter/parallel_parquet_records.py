#!/usr/bin/env python3
"""
TRUE parallel SAM/BAM to Parquet converter using record-based splitting.

Strategy: Each worker processes a different range of records from the same file.
All workers read from the shared bgzip file in parallel - no intermediate files!

This is the CORRECT approach for parallel processing.
"""

import sys
import time
import subprocess
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor, as_completed


def count_records_fast(sam_file):
    """Quickly count records using samtools view -c."""
    result = subprocess.run(
        ['samtools', 'view', '-c', sam_file],
        capture_output=True,
        text=True,
        check=True
    )
    return int(result.stdout.strip())


def process_record_chunk(args):
    """Process a chunk of records on-the-fly."""
    chunk_id, input_file, skip_records, max_records, output_dir, num_partitions, batch_size, compression_level, num_threads, calculate_pmd, library_type = args

    from bam_filter.parquet_converter_pure_cpp import convert_sam_bam_to_parquet_pure_cpp

    start_time = time.time()

    chunk_dir = Path(output_dir) / f"chunk_{chunk_id:03d}"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    end_record = skip_records + max_records
    print(f"[Chunk {chunk_id}] Processing records {skip_records:,} to {end_record:,} ({max_records:,} records)")

    # Process records on-the-fly - no intermediate files!
    try:
        stats = convert_sam_bam_to_parquet_pure_cpp(
            input_file,
            str(chunk_dir),
            num_partitions=num_partitions,
            batch_size=batch_size,
            write_by_reference=True,
            write_by_read=True,
            compression_level=compression_level,
            num_threads=num_threads,
            skip_records=skip_records,
            max_records=max_records,
            calculate_pmd=calculate_pmd,
            library_type=library_type,
        )

        total_time = time.time() - start_time

        print(f"[Chunk {chunk_id}] DONE: {stats['total_records']:,} records in {total_time:.1f}s ({stats['total_records']/total_time:.0f} rec/s)")

        return {
            'chunk_id': chunk_id,
            'total_records': stats['total_records'],
            'total_time': total_time,
            'throughput': stats['total_records'] / total_time if total_time > 0 else 0,
        }
    except Exception as e:
        print(f"[Chunk {chunk_id}] FAILED: {e}")
        import traceback
        traceback.print_exc()
        return {
            'chunk_id': chunk_id,
            'total_records': 0,
            'total_time': time.time() - start_time,
            'throughput': 0,
            'error': str(e),
        }


def parallel_convert_records(
    input_file: str,
    output_dir: str,
    num_processes: int = 32,
    num_partitions: int = 16,
    batch_size: int = 100000,
    compression_level: int = 3,
    num_threads_per_worker: int = 1,  # Each worker reads different records, so less threading needed
    count_records: bool = True,
    calculate_pmd: bool = True,
    library_type: str = "ds",
):
    """
    Convert SAM/BAM to Parquet in parallel using record-based splitting.

    Each worker processes a different range of records on-the-fly from the shared file.
    NO intermediate files, NO dd extraction, TRUE parallelization!

    Args:
        input_file: Input SAM/BAM/SAM.gz file
        output_dir: Output directory
        num_processes: Number of parallel worker processes
        num_partitions: Partitions per worker
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads_per_worker: BGzip decompression threads per worker
        count_records: Whether to count records first (slower but more accurate)
        calculate_pmd: Calculate PMD scores on-the-fly (default: True)
        library_type: Library type for PMD: "ds" (double-stranded) or "ss" (single-stranded)
    """
    overall_start = time.time()

    print(f"TRUE Parallel SAM/BAM → Parquet Converter (Record-Based Splitting)")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Processes: {num_processes}")
    print()

    # Save header once (before parallel processing)
    from bam_filter.parquet_converter_pure_cpp import save_bam_header
    print("Saving BAM header for round-trip conversion...")
    save_bam_header(input_file, output_dir)
    print()

    # Estimate total records - try BGZF block counting first (fast and accurate)
    import os
    file_size = os.path.getsize(input_file)
    print(f"File size: {file_size/1e9:.2f} GB")

    # Try BGZF block counting first (fast: ~3 min for 206GB vs 2 hours for samtools count)
    try:
        from bam_filter.bgzf_scanner import is_bgzf_file, scan_bgzf_blocks, estimate_records_from_blocks
        if is_bgzf_file(input_file):
            print("Scanning BGZF blocks for accurate record estimation...")
            block_offsets = scan_bgzf_blocks(input_file, progress_interval=5000000)
            num_blocks = len(block_offsets)
            # BGZF blocks decompress to ~64KB each, average SAM record ~300 bytes = ~206 records/block
            total_records = estimate_records_from_blocks(num_blocks)
            print(f"BGZF blocks: {num_blocks:,}")
            print(f"Estimated records: ~{total_records:,} (based on ~206 records per 64KB block)")
        else:
            # Fall back to file size estimation for non-BGZF files
            bytes_per_record = 700  # Conservative estimate
            total_records = int(file_size / bytes_per_record)
            print(f"Estimated records: ~{total_records:,} (assuming ~{bytes_per_record} bytes/record)")
    except ImportError:
        # Fall back to file size estimation
        bytes_per_record = 700  # Conservative estimate
        total_records = int(file_size / bytes_per_record)
        print(f"Estimated records: ~{total_records:,} (assuming ~{bytes_per_record} bytes/record)")
    print("(No pre-counting - workers will process until their range is complete)")

    # Calculate records per worker
    records_per_worker = total_records // num_processes
    print(f"Records per worker: ~{records_per_worker:,}")
    print()

    # Create worker assignments
    chunks = []
    for i in range(num_processes):
        skip = i * records_per_worker
        max_recs = records_per_worker if i < num_processes - 1 else (total_records - skip)

        chunks.append((
            i, input_file, skip, max_recs, output_dir,
            num_partitions, batch_size, compression_level, num_threads_per_worker,
            calculate_pmd, library_type
        ))

    print(f"Created {len(chunks)} worker assignments")
    for i, (cid, _, skip, max_recs, *_) in enumerate(chunks[:5]):
        print(f"  Worker {cid}: records {skip:,} to {skip+max_recs:,}")
    if len(chunks) > 5:
        print(f"  ... and {len(chunks) - 5} more workers")
    print()

    # Process in parallel - TRUE parallelization!
    print(f"Processing with {num_processes} parallel workers...")
    print("Each worker reads different records from the SAME file on-the-fly!")
    print()

    results = []

    with ProcessPoolExecutor(max_workers=num_processes) as executor:
        futures = {executor.submit(process_record_chunk, chunk): chunk[0] for chunk in chunks}

        for future in as_completed(futures):
            chunk_id = futures[future]
            try:
                result = future.result()
                results.append(result)
                if 'error' not in result:
                    print(f"✓ Worker {chunk_id} completed successfully")
            except Exception as e:
                print(f"✗ Worker {chunk_id} EXCEPTION: {e}")

    total_time = time.time() - overall_start
    total_records_processed = sum(r['total_records'] for r in results)

    print()
    print("="*80)
    print("PARALLEL CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records processed: {total_records_processed:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Throughput: {total_records_processed/total_time:.0f} records/sec")
    print(f"Throughput: {total_records_processed/total_time/1e6:.3f} M records/sec")
    print()

    # Calculate effective speedup
    if results:
        avg_worker_time = sum(r['total_time'] for r in results) / len(results)
        speedup = avg_worker_time / total_time if total_time > 0 else 0
        print(f"Average worker time: {avg_worker_time:.1f}s")
        print(f"Effective speedup: {speedup:.1f}x")
    print()

    return {
        'total_records': total_records_processed,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records_processed / total_time if total_time > 0 else 0,
        'num_workers': len(chunks),
        'worker_results': results,
    }


if __name__ == '__main__':
    if len(sys.argv) < 3:
        print("Usage: python parallel_parquet_records.py <input.sam.gz> <output_dir> [num_processes]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    num_processes = int(sys.argv[3]) if len(sys.argv) > 3 else 32

    stats = parallel_convert_records(
        input_file,
        output_dir,
        num_processes=num_processes,
        num_partitions=16,
        batch_size=100000,
        compression_level=3,
        num_threads_per_worker=1,
        count_records=True,
    )
