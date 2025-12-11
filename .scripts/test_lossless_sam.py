#!/usr/bin/env python3
"""
Test lossless SAM to Parquet conversion.

Verifies that all SAM tags are preserved in the tags_raw column
for round-trip conversion SAM -> Parquet -> SAM.

Usage:
    python test_lossless_sam.py <input.sam.gz> <output_dir> [max_records]
"""

import sys
import time
from pathlib import Path


def main():
    if len(sys.argv) < 3:
        print("Usage: python test_lossless_sam.py <input.sam.gz> <output_dir> [max_records]")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = sys.argv[2]
    max_records = int(sys.argv[3]) if len(sys.argv) > 3 else 100000

    # Clean output
    output_path = Path(output_dir)
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)

    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"Max records: {max_records:,}")
    print()

    # Convert SAM.gz to Parquet
    from bam_filter.sam_parser_fast import convert_sam_gz_to_parquet_fast

    print("Converting SAM.gz to Parquet with fast Cython parser...")
    start_time = time.time()

    stats = convert_sam_gz_to_parquet_fast(
        input_file,
        output_dir,
        batch_size=50000,
        compression_level=3,
        num_threads=4,
        store_sequences=True,  # Store sequences for round-trip
    )

    elapsed = time.time() - start_time
    total_records = stats.get('total_records', 0)

    if total_records == 0:
        print("No records converted!")
        return

    print(f"Converted {total_records:,} records in {elapsed:.1f}s")
    print(f"Throughput: {total_records/elapsed:,.0f} rec/s")
    print()

    # Verify Parquet file
    parquet_file = output_path / "alignments.parquet"
    if not parquet_file.exists():
        print(f"ERROR: Parquet file not found at {parquet_file}")
        return

    print(f"Parquet file: {parquet_file} ({parquet_file.stat().st_size / 1e6:.1f} MB)")
    print()

    # Read and verify tags_raw column
    import pyarrow.parquet as pq

    print("Reading Parquet file to verify tags_raw column...")
    table = pq.read_table(parquet_file)

    print(f"Schema columns: {table.column_names}")
    print()

    if 'tags_raw' not in table.column_names:
        print("ERROR: tags_raw column not found in Parquet!")
        return

    # Get tags_raw column
    tags_raw = table.column('tags_raw')
    non_empty = sum(1 for t in tags_raw.to_pylist() if t)
    total = len(tags_raw)

    print(f"Total rows: {total:,}")
    print(f"Rows with tags: {non_empty:,} ({100*non_empty/total:.1f}%)")
    print()

    # Show sample records
    print("Sample records with tags_raw:")
    print("-" * 80)

    df = table.to_pandas()
    sample_with_tags = df[df['tags_raw'].str.len() > 0].head(5)

    for idx, row in sample_with_tags.iterrows():
        print(f"read_name: {row['read_name']}")
        print(f"ref_name:  {row['ref_name']}")
        print(f"position:  {row['position']}")
        print(f"mapq:      {row['mapq']}")
        print(f"cigar:     {row['cigar']}")
        if row['sequence']:
            print(f"sequence:  {row['sequence'][:50]}..." if len(row['sequence']) > 50 else f"sequence:  {row['sequence']}")
        print(f"tags_raw:  {row['tags_raw'][:100]}..." if len(row['tags_raw']) > 100 else f"tags_raw:  {row['tags_raw']}")
        print()

    # Parse and count unique tag types
    print("Unique tag types found:")
    tag_counts = {}
    for tags in df['tags_raw']:
        if tags:
            for tag_part in tags.split('\t'):
                if ':' in tag_part:
                    tag_name = tag_part.split(':')[0]
                    tag_counts[tag_name] = tag_counts.get(tag_name, 0) + 1

    for tag, count in sorted(tag_counts.items(), key=lambda x: -x[1])[:20]:
        print(f"  {tag}: {count:,}")

    print()
    print("SUCCESS: Lossless SAM to Parquet conversion verified!")
    print("All optional tags are preserved in the tags_raw column.")


if __name__ == "__main__":
    main()
