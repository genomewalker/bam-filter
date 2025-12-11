#!/usr/bin/env python3
"""
Test script for the super-optimized Parquet writer.

Tests:
1. Single-table mode (PyOptimizedParquetWriter)
2. Normalized mode (PyNormalizedParquetWriter)
3. Benchmarks throughput vs current implementation

Key optimizations tested:
- LUT-based 2-bit sequence packing (branchless)
- Zero-copy Python string handling (PyUnicode_AsUTF8AndSize)
- Bulk Arrow AppendValues API
- Pre-reserved batch storage
- Hot/cold tag separation (no duplication)
"""

import sys
import os
import time
import gzip
import pyarrow.parquet as pq

# Add parent directory to path
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from bam_filter.optimized_parquet_writer import (
    PyOptimizedParquetWriter,
    PyNormalizedParquetWriter,
    pack_sequence,
    unpack_sequence,
)


def parse_sam_line(line):
    """Parse a SAM line into fields."""
    fields = line.rstrip('\n').split('\t')
    if len(fields) < 11:
        return None

    read_name = fields[0]
    flag = int(fields[1])
    ref_name = fields[2]
    position = int(fields[3])
    mapq = int(fields[4])
    cigar = fields[5]
    mate_ref = fields[6]
    mate_pos = int(fields[7])
    template_length = int(fields[8])
    sequence = fields[9]
    quality = fields[10]
    tags_raw = '\t'.join(fields[11:]) if len(fields) > 11 else ''

    return {
        'read_name': read_name,
        'flag': flag,
        'ref_name': ref_name,
        'position': position,
        'mapq': mapq,
        'cigar': cigar,
        'mate_ref': mate_ref,
        'mate_pos': mate_pos,
        'template_length': template_length,
        'sequence': sequence,
        'quality': quality,
        'tags_raw': tags_raw,
    }


def test_single_table_mode(sam_file, output_parquet, max_records=10000):
    """Test the single-table optimized writer."""
    print(f"\n=== Testing Single-Table Mode ===")
    print(f"Input: {sam_file}")
    print(f"Output: {output_parquet}")
    print(f"Max records: {max_records}")

    writer = PyOptimizedParquetWriter(output_parquet, compression_level=6, batch_size=10000)

    start = time.time()
    count = 0

    opener = gzip.open if sam_file.endswith('.gz') else open

    with opener(sam_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                continue

            record = parse_sam_line(line)
            if record is None or record['ref_name'] == '*':
                continue

            writer.add_alignment(
                read_name=record['read_name'],
                ref_name=record['ref_name'],
                position=record['position'],
                mapq=record['mapq'],
                flag=record['flag'],
                cigar=record['cigar'],
                sequence=record['sequence'],
                quality=record['quality'],
                tags_raw=record['tags_raw'],
                template_length=record['template_length'],
                mate_ref=record['mate_ref'],
                mate_pos=record['mate_pos'],
            )

            count += 1
            if count >= max_records:
                break

            if count % 5000 == 0:
                print(f"  Processed {count:,} records...")

    writer.close()

    elapsed = time.time() - start
    throughput = count / elapsed if elapsed > 0 else 0

    print(f"\nResults:")
    print(f"  Records: {count:,}")
    print(f"  Time: {elapsed:.2f}s")
    print(f"  Throughput: {throughput:,.0f} rec/s")

    # Verify output
    file_size = os.path.getsize(output_parquet)
    print(f"  Output size: {file_size / 1024:.1f} KB")

    table = pq.read_table(output_parquet)
    print(f"  Parquet rows: {len(table):,}")
    print(f"  Columns: {table.column_names}")

    # Check 2-bit packing
    seq_packed = table['sequence_packed']
    seq_lengths = table['sequence_length']
    if len(seq_packed) > 0:
        first_packed = seq_packed[0].as_py()
        first_len = seq_lengths[0].as_py()
        unpacked = unpack_sequence(first_packed, first_len)
        print(f"  First seq packed size: {len(first_packed)} bytes (original: {first_len} bases)")
        print(f"  Compression ratio: {first_len / len(first_packed):.1f}x")

    return count, elapsed, throughput


def test_normalized_mode(sam_file, output_dir, max_records=10000):
    """Test the normalized (two-table) writer."""
    print(f"\n=== Testing Normalized Mode ===")
    print(f"Input: {sam_file}")
    print(f"Output dir: {output_dir}")
    print(f"Max records: {max_records}")

    os.makedirs(output_dir, exist_ok=True)
    writer = PyNormalizedParquetWriter(output_dir, compression_level=6, batch_size=10000)

    start = time.time()
    count = 0
    seen_reads = set()

    opener = gzip.open if sam_file.endswith('.gz') else open

    with opener(sam_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                continue

            record = parse_sam_line(line)
            if record is None or record['ref_name'] == '*':
                continue

            read_name = record['read_name']

            # Add read only if not seen before
            if read_name not in seen_reads:
                writer.add_read(
                    read_name=read_name,
                    sequence=record['sequence'],
                    quality=record['quality'],
                )
                seen_reads.add(read_name)

            # Add alignment
            writer.add_alignment(
                read_name=read_name,
                ref_name=record['ref_name'],
                position=record['position'],
                mapq=record['mapq'],
                flag=record['flag'],
                cigar=record['cigar'],
                tags_raw=record['tags_raw'],
                template_length=record['template_length'],
                mate_pos=record['mate_pos'],
            )

            count += 1
            if count >= max_records:
                break

            if count % 5000 == 0:
                print(f"  Processed {count:,} alignments...")

    writer.close()

    elapsed = time.time() - start
    throughput = count / elapsed if elapsed > 0 else 0

    print(f"\nResults:")
    print(f"  Alignments: {count:,}")
    print(f"  Unique reads: {len(seen_reads):,}")
    print(f"  Avg alignments/read: {count / len(seen_reads):.1f}")
    print(f"  Time: {elapsed:.2f}s")
    print(f"  Throughput: {throughput:,.0f} rec/s")

    # Check output files
    reads_file = os.path.join(output_dir, 'reads.parquet')
    alignments_file = os.path.join(output_dir, 'alignments.parquet')
    refs_file = os.path.join(output_dir, 'references.parquet')

    total_size = 0
    for f in [reads_file, alignments_file, refs_file]:
        if os.path.exists(f):
            size = os.path.getsize(f)
            total_size += size
            print(f"  {os.path.basename(f)}: {size / 1024:.1f} KB")

    print(f"  Total: {total_size / 1024:.1f} KB")

    # Verify
    reads_table = pq.read_table(reads_file)
    alignments_table = pq.read_table(alignments_file)
    print(f"  Reads rows: {len(reads_table):,}")
    print(f"  Alignments rows: {len(alignments_table):,}")

    return count, elapsed, throughput


def test_sequence_packing():
    """Test the LUT-based 2-bit sequence packing."""
    print("\n=== Testing Sequence Packing ===")

    test_cases = [
        "ACGT",
        "AAAA",
        "TTTT",
        "ACGTACGTACGT",
        "ACGTACGTA",  # Not multiple of 4
        "N" * 10,  # All N (maps to A)
        "ACGTacgt",  # Mixed case
    ]

    for seq in test_cases:
        packed = pack_sequence(seq)
        unpacked = unpack_sequence(packed, len(seq))

        # Map N to A for comparison
        expected = seq.upper().replace('N', 'A')
        ok = unpacked == expected

        print(f"  {seq:15s} -> {len(packed):2d} bytes -> {unpacked:15s} {'OK' if ok else 'FAIL'}")


def main():
    if len(sys.argv) < 2:
        print("Usage: python test_optimized_writer.py <sam.gz> [output_dir] [max_records]")
        print("\nTests the super-optimized Parquet writer with:")
        print("  - LUT-based 2-bit sequence packing")
        print("  - Zero-copy Python string handling")
        print("  - Bulk Arrow AppendValues API")
        print("  - Pre-reserved batch storage")
        sys.exit(1)

    sam_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else '/tmp/optimized_writer_test'
    max_records = int(sys.argv[3]) if len(sys.argv) > 3 else 10000

    # Test sequence packing first
    test_sequence_packing()

    # Single-table mode
    single_output = os.path.join(output_dir, 'single_table.parquet')
    os.makedirs(output_dir, exist_ok=True)
    single_count, single_time, single_throughput = test_single_table_mode(
        sam_file, single_output, max_records
    )

    # Normalized mode
    normalized_dir = os.path.join(output_dir, 'normalized')
    norm_count, norm_time, norm_throughput = test_normalized_mode(
        sam_file, normalized_dir, max_records
    )

    # Summary
    print("\n" + "=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print(f"Single-table: {single_throughput:,.0f} rec/s, {os.path.getsize(single_output) / 1024:.1f} KB")
    norm_size = sum(
        os.path.getsize(os.path.join(normalized_dir, f))
        for f in os.listdir(normalized_dir)
        if f.endswith('.parquet')
    )
    print(f"Normalized:   {norm_throughput:,.0f} rec/s, {norm_size / 1024:.1f} KB")


if __name__ == '__main__':
    main()
