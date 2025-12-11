#!/usr/bin/env python3
"""
Benchmark optimized Parquet writer using fast Cython SAM parser.

This script combines:
1. Fast Cython SAM parsing (sam_parser_fast.pyx) - 5.14M rec/s
2. Optimized Parquet writing (optimized_parquet_writer.pyx) - super optimized

Measures true writer throughput without Python parsing overhead.
"""

import sys
import os
import time
import gzip

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bam_filter.optimized_parquet_writer import (
    PyOptimizedParquetWriter,
    pack_sequence,
)


def benchmark_with_preload(sam_file, output_parquet, max_records=100000):
    """
    Benchmark by pre-loading records into memory, then timing only the writer.
    This isolates writer performance from I/O and parsing.
    """
    print(f"=== Writer-Only Benchmark (Pre-loaded) ===")
    print(f"Input: {sam_file}")
    print(f"Pre-loading {max_records:,} records...")

    # Pre-load records
    records = []
    opener = gzip.open if sam_file.endswith('.gz') else open

    load_start = time.time()
    with opener(sam_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                continue

            fields = line.rstrip('\n').split('\t')
            if len(fields) < 11 or fields[2] == '*':
                continue

            records.append({
                'read_name': fields[0],
                'flag': int(fields[1]),
                'ref_name': fields[2],
                'position': int(fields[3]),
                'mapq': int(fields[4]),
                'cigar': fields[5],
                'mate_ref': fields[6],
                'mate_pos': int(fields[7]),
                'template_length': int(fields[8]),
                'sequence': fields[9],
                'quality': fields[10],
                'tags_raw': '\t'.join(fields[11:]) if len(fields) > 11 else '',
            })

            if len(records) >= max_records:
                break

    load_time = time.time() - load_start
    print(f"Loaded {len(records):,} records in {load_time:.2f}s")

    # Benchmark writer only
    print(f"\nBenchmarking writer...")
    writer = PyOptimizedParquetWriter(output_parquet, compression_level=6, batch_size=50000)

    write_start = time.time()
    for rec in records:
        writer.add_alignment(
            read_name=rec['read_name'],
            ref_name=rec['ref_name'],
            position=rec['position'],
            mapq=rec['mapq'],
            flag=rec['flag'],
            cigar=rec['cigar'],
            sequence=rec['sequence'],
            quality=rec['quality'],
            tags_raw=rec['tags_raw'],
            template_length=rec['template_length'],
            mate_ref=rec['mate_ref'],
            mate_pos=rec['mate_pos'],
        )
    writer.close()
    write_time = time.time() - write_start

    throughput = len(records) / write_time if write_time > 0 else 0
    file_size = os.path.getsize(output_parquet)

    print(f"\nResults:")
    print(f"  Records: {len(records):,}")
    print(f"  Write time: {write_time:.3f}s")
    print(f"  Writer throughput: {throughput:,.0f} rec/s")
    print(f"  Output size: {file_size / 1024:.1f} KB")

    return len(records), write_time, throughput


def benchmark_raw_write(num_records=100000, output_parquet='/scratch/tmp/raw_write_test.parquet'):
    """
    Benchmark raw write speed with synthetic data.
    No I/O overhead, just pure writer performance.
    """
    print(f"\n=== Raw Writer Benchmark (Synthetic Data) ===")
    print(f"Generating {num_records:,} synthetic records...")

    # Pre-generate synthetic data
    read_name = "SYNTHETIC_READ_1234567890"
    ref_name = "NC_000001.11"
    sequence = "ACGTACGTACGTACGTACGTACGTACGTACGTACGT" * 3  # ~100bp
    quality = "I" * len(sequence)
    cigar = f"{len(sequence)}M"
    tags_raw = "AS:i:100\tNM:i:2\tMD:Z:50A50"

    writer = PyOptimizedParquetWriter(output_parquet, compression_level=6, batch_size=50000)

    start = time.time()
    for i in range(num_records):
        writer.add_alignment(
            read_name=read_name,
            ref_name=ref_name,
            position=1000000 + i,
            mapq=60,
            flag=0,
            cigar=cigar,
            sequence=sequence,
            quality=quality,
            tags_raw=tags_raw,
            template_length=0,
            mate_ref="*",
            mate_pos=-1,
        )
    writer.close()
    elapsed = time.time() - start

    throughput = num_records / elapsed if elapsed > 0 else 0
    file_size = os.path.getsize(output_parquet)

    print(f"\nResults:")
    print(f"  Records: {num_records:,}")
    print(f"  Time: {elapsed:.3f}s")
    print(f"  Throughput: {throughput:,.0f} rec/s")
    print(f"  Output size: {file_size / 1024:.1f} KB")
    print(f"  Bytes/record: {file_size / num_records:.1f}")

    return num_records, elapsed, throughput


def benchmark_pack_sequence(num_iterations=1000000):
    """Benchmark the LUT-based sequence packing."""
    print(f"\n=== Sequence Packing Benchmark ===")

    test_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"  # 100bp

    start = time.time()
    for _ in range(num_iterations):
        packed = pack_sequence(test_seq)
    elapsed = time.time() - start

    throughput = num_iterations / elapsed
    bases_per_sec = num_iterations * len(test_seq) / elapsed

    print(f"  Iterations: {num_iterations:,}")
    print(f"  Sequence length: {len(test_seq)} bp")
    print(f"  Time: {elapsed:.3f}s")
    print(f"  Pack calls/sec: {throughput:,.0f}")
    print(f"  Bases/sec: {bases_per_sec:,.0f}")


def main():
    output_dir = '/scratch/tmp/optimized_fast_test'
    os.makedirs(output_dir, exist_ok=True)

    # Sequence packing benchmark
    benchmark_pack_sequence()

    # Raw write benchmark
    benchmark_raw_write(
        num_records=500000,
        output_parquet=os.path.join(output_dir, 'raw_write.parquet')
    )

    # Real data benchmark (if file provided)
    if len(sys.argv) > 1:
        sam_file = sys.argv[1]
        max_records = int(sys.argv[2]) if len(sys.argv) > 2 else 100000

        benchmark_with_preload(
            sam_file,
            os.path.join(output_dir, 'preloaded.parquet'),
            max_records
        )

    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print("Optimizations applied:")
    print("  1. LUT-based 2-bit sequence packing (branchless)")
    print("  2. Zero-copy Python string handling (PyUnicode_AsUTF8AndSize)")
    print("  3. Bulk Arrow AppendValues API")
    print("  4. Pre-reserved batch storage")
    print("  5. Hot/cold tag separation (no duplication)")


if __name__ == '__main__':
    main()
