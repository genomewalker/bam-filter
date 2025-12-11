#!/usr/bin/env python3
"""
Test end_position and is_reverse derived columns in optimized Parquet output.
"""

import sys
sys.path.insert(0, '/maps/projects/fernandezguerra/apps/repos/bam-filter')

from bam_filter.parquet_converter_optimized import convert_bam_to_optimized_parquet
import pyarrow.parquet as pq
import os


def test_derived_columns():
    """Test end_position and is_reverse columns."""
    input_file = "/scratch/tmp/zs_test/test_small.sam.gz"
    output_file = "/scratch/tmp/zs_test/test_derived.parquet"

    # Remove existing output
    if os.path.exists(output_file):
        os.remove(output_file)

    print("Testing derived columns (end_position, is_reverse)...")
    stats = convert_bam_to_optimized_parquet(
        input_file=input_file,
        output_file=output_file,
        batch_size=10000,
        compression_level=3,
        num_threads=1,
        calculate_pmd=False,
        library_type="ds",
        store_sequences=True
    )

    print(f"\nProcessed {stats['total_records']} records in {stats['processing_time_seconds']:.2f}s")

    # Check schema and data
    table = pq.read_table(output_file)
    print(f"\nSchema columns: {table.schema.names}")

    # Check if derived columns exist
    assert 'end_position' in table.schema.names, "Missing end_position column"
    assert 'is_reverse' in table.schema.names, "Missing is_reverse column"
    print("Both derived columns present in schema.")

    # Get sample data
    positions = table['position'].to_pylist()[:10]
    end_positions = table['end_position'].to_pylist()[:10]
    flags = table['flag'].to_pylist()[:10]
    is_reverse = table['is_reverse'].to_pylist()[:10]
    cigars = table['cigar'].to_pylist()[:10]

    print(f"\nSample data (first 10 records):")
    print(f"{'pos':>10} {'end_pos':>10} {'diff':>6} {'flag':>6} {'reverse':>8} {'cigar':>20}")
    print("-" * 70)
    for i in range(min(10, len(positions))):
        diff = end_positions[i] - positions[i] if end_positions[i] and positions[i] else 0
        expected_reverse = (flags[i] & 0x10) != 0
        reverse_match = "OK" if is_reverse[i] == expected_reverse else "MISMATCH"
        print(f"{positions[i]:>10} {end_positions[i]:>10} {diff:>6} {flags[i]:>6} {is_reverse[i]:>8} {cigars[i][:20]:>20}")

    # Validate is_reverse
    print("\nValidating is_reverse column...")
    all_flags = table['flag'].to_pylist()
    all_is_reverse = table['is_reverse'].to_pylist()
    mismatches = 0
    for i, (flag, rev) in enumerate(zip(all_flags, all_is_reverse)):
        expected = (flag & 0x10) != 0
        if rev != expected:
            mismatches += 1
            if mismatches <= 5:
                print(f"  Mismatch at {i}: flag={flag}, is_reverse={rev}, expected={expected}")

    if mismatches == 0:
        print(f"  All {len(all_flags)} records have correct is_reverse values.")
    else:
        print(f"  WARNING: {mismatches} mismatches found!")

    # Check end_position sanity
    print("\nValidating end_position column...")
    all_end_positions = table['end_position'].to_pylist()
    all_positions = table['position'].to_pylist()
    invalid_count = 0
    for i, (pos, end_pos) in enumerate(zip(all_positions, all_end_positions)):
        if end_pos < pos:
            invalid_count += 1
            if invalid_count <= 5:
                print(f"  Invalid at {i}: pos={pos}, end_pos={end_pos}")

    if invalid_count == 0:
        print(f"  All {len(all_positions)} records have valid end_position >= position.")
    else:
        print(f"  WARNING: {invalid_count} invalid end_positions found!")

    # Calculate alignment length distribution
    lengths = [end - pos for pos, end in zip(all_positions, all_end_positions) if end >= pos]
    if lengths:
        print(f"\nAlignment length stats:")
        print(f"  Min: {min(lengths)}")
        print(f"  Max: {max(lengths)}")
        print(f"  Mean: {sum(lengths)/len(lengths):.1f}")

    # Reverse strand stats
    reverse_count = sum(1 for r in all_is_reverse if r)
    forward_count = len(all_is_reverse) - reverse_count
    print(f"\nStrand distribution:")
    print(f"  Forward: {forward_count} ({100*forward_count/len(all_is_reverse):.1f}%)")
    print(f"  Reverse: {reverse_count} ({100*reverse_count/len(all_is_reverse):.1f}%)")

    print("\nTest completed successfully!")


if __name__ == "__main__":
    test_derived_columns()
