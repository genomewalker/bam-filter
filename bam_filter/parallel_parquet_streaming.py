#!/usr/bin/env python3
"""
Streaming Parallel BAM/SAM to Parquet Converter.

This approach uses a producer-consumer pattern:
1. Single reader process reads BAM/SAM records sequentially
2. Multiple worker processes convert record batches to Parquet in parallel

This guarantees all records are captured and works with any compression format.
"""

import os
import sys
import time
import multiprocessing as mp
from multiprocessing import Queue, Process
from pathlib import Path
from typing import Dict, Any, List, Optional
from queue import Empty


POISON_PILL = None


def parquet_writer_worker(
    worker_id: int,
    input_queue: Queue,
    output_dir: str,
    num_partitions: int,
    compression_level: int,
    calculate_pmd: bool,
    library_type: str,
    result_queue: Queue,
):
    """
    Worker process that receives record batches and writes them to Parquet.
    """
    from bam_filter.parquet_converter_pure_cpp import write_batch_to_parquet

    records_processed = 0
    batches_processed = 0
    start_time = time.time()

    chunk_dir = Path(output_dir) / f"chunk_{worker_id:03d}"
    chunk_dir.mkdir(parents=True, exist_ok=True)

    while True:
        try:
            item = input_queue.get(timeout=60)

            if item is POISON_PILL:
                # Put poison pill back for other workers
                input_queue.put(POISON_PILL)
                break

            batch_id, records = item
            batch_size = len(records)

            # Write batch to Parquet
            write_batch_to_parquet(
                records,
                str(chunk_dir),
                num_partitions=num_partitions,
                compression_level=compression_level,
                calculate_pmd=calculate_pmd,
                library_type=library_type,
            )

            records_processed += batch_size
            batches_processed += 1

            if batches_processed % 10 == 0:
                elapsed = time.time() - start_time
                rate = records_processed / elapsed if elapsed > 0 else 0
                print(f"[Worker {worker_id}] {records_processed:,} records, {rate:,.0f} rec/sec")

        except Empty:
            continue
        except Exception as e:
            print(f"[Worker {worker_id}] ERROR: {e}")
            import traceback
            traceback.print_exc()
            break

    elapsed = time.time() - start_time
    result_queue.put({
        'worker_id': worker_id,
        'records_processed': records_processed,
        'batches_processed': batches_processed,
        'elapsed_seconds': elapsed,
    })


def streaming_convert(
    input_file: str,
    output_dir: str,
    num_workers: int = 32,
    batch_size: int = 100000,
    num_partitions: int = 16,
    compression_level: int = 3,
    calculate_pmd: bool = True,
    library_type: str = "ds",
    queue_size: int = 64,
) -> Dict[str, Any]:
    """
    Stream BAM/SAM file through parallel Parquet writers.

    Args:
        input_file: Path to BAM/SAM file
        output_dir: Output directory for Parquet files
        num_workers: Number of parallel writer workers
        batch_size: Records per batch sent to workers
        num_partitions: Partitions per chunk for Parquet
        compression_level: Parquet compression level
        calculate_pmd: Whether to calculate PMD scores
        library_type: "ds" (double-stranded) or "ss" (single-stranded)
        queue_size: Size of inter-process queue

    Returns:
        Dictionary with conversion statistics
    """
    import pysam

    print("="*80)
    print("Streaming Parallel BAM/SAM → Parquet Converter")
    print("="*80)
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Workers: {num_workers}")
    print(f"Batch size: {batch_size:,}")
    print()

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Save header
    print("Saving BAM header...")
    from bam_filter.parquet_converter_pure_cpp import save_bam_header
    save_bam_header(input_file, output_dir)
    print()

    # Create queues
    work_queue = mp.Queue(maxsize=queue_size)
    result_queue = mp.Queue()

    # Start worker processes
    workers = []
    for i in range(num_workers):
        p = Process(
            target=parquet_writer_worker,
            args=(i, work_queue, output_dir, num_partitions, compression_level,
                  calculate_pmd, library_type, result_queue)
        )
        p.start()
        workers.append(p)

    print(f"Started {num_workers} worker processes")
    print()

    # Read records and distribute to workers
    start_time = time.time()
    total_records = 0
    batch_id = 0
    current_batch = []

    print("Reading and distributing records...")
    with pysam.AlignmentFile(input_file, "r", check_sq=False) as samfile:
        for read in samfile:
            # Collect record data
            record = {
                'qname': read.query_name,
                'flag': read.flag,
                'rname': read.reference_name or "*",
                'pos': read.reference_start,
                'mapq': read.mapping_quality,
                'cigar': read.cigarstring or "*",
                'rnext': read.next_reference_name or "*",
                'pnext': read.next_reference_start,
                'tlen': read.template_length,
                'seq': read.query_sequence or "*",
                'qual': read.qual or "*",
            }

            # Add optional tags
            for tag, value in read.tags:
                record[tag] = value

            current_batch.append(record)

            if len(current_batch) >= batch_size:
                work_queue.put((batch_id, current_batch))
                total_records += len(current_batch)
                batch_id += 1
                current_batch = []

                if batch_id % 100 == 0:
                    elapsed = time.time() - start_time
                    rate = total_records / elapsed if elapsed > 0 else 0
                    print(f"Read {total_records:,} records ({rate:,.0f} rec/sec)")

    # Send remaining batch
    if current_batch:
        work_queue.put((batch_id, current_batch))
        total_records += len(current_batch)

    read_time = time.time() - start_time
    print(f"\nFinished reading: {total_records:,} records in {read_time:.1f}s")

    # Signal workers to stop
    work_queue.put(POISON_PILL)

    # Wait for workers
    print("Waiting for workers to complete...")
    for p in workers:
        p.join()

    # Collect results
    worker_results = []
    while not result_queue.empty():
        worker_results.append(result_queue.get())

    total_time = time.time() - start_time

    # Print summary
    print()
    print("="*80)
    print("STREAMING CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records: {total_records:,}")
    print(f"Total time: {total_time:.1f}s ({total_time/60:.2f} minutes)")
    print(f"Read throughput: {total_records/read_time:,.0f} rec/sec")
    print(f"Overall throughput: {total_records/total_time:,.0f} rec/sec")
    print()

    worker_records = sum(r['records_processed'] for r in worker_results)
    print(f"Records written by workers: {worker_records:,}")

    return {
        'total_records': total_records,
        'read_time_seconds': read_time,
        'total_time_seconds': total_time,
        'throughput_records_per_sec': total_records / total_time if total_time > 0 else 0,
        'num_workers': num_workers,
        'worker_results': worker_results,
    }


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description="Streaming Parallel Converter")
    parser.add_argument('input_file', help="Input BAM/SAM file")
    parser.add_argument('output_dir', help="Output directory")
    parser.add_argument('-t', '--threads', type=int, default=32, help="Number of workers")
    parser.add_argument('--batch-size', type=int, default=100000, help="Batch size")

    args = parser.parse_args()

    stats = streaming_convert(
        args.input_file,
        args.output_dir,
        num_workers=args.threads,
        batch_size=args.batch_size,
    )
