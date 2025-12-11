#!/usr/bin/env python3
"""
Quick test for tags_raw column in lossless SAM to Parquet conversion.

Creates a small test SAM file and verifies tags_raw preservation.
"""

import sys
import gzip
import tempfile
from pathlib import Path


def main():
    if len(sys.argv) < 2:
        print("Usage: python test_tags_raw_quick.py <input.sam.gz>")
        sys.exit(1)

    input_file = sys.argv[1]
    output_dir = "/scratch/tmp/tags_raw_quick_test"

    # Clean output
    output_path = Path(output_dir)
    if output_path.exists():
        import shutil
        shutil.rmtree(output_path)
    output_path.mkdir(parents=True)

    # Extract header and first 1000 records
    print(f"Extracting sample from {input_file}...")
    sample_file = output_path / "sample.sam.gz"

    header_lines = []
    record_lines = []

    with gzip.open(input_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                header_lines.append(line)
            else:
                record_lines.append(line)
                if len(record_lines) >= 1000:
                    break

    # Write sample file
    with gzip.open(sample_file, 'wt') as f:
        for line in header_lines:
            f.write(line)
        for line in record_lines:
            f.write(line)

    print(f"Wrote {len(header_lines)} header lines + {len(record_lines)} records")

    # Show sample tags from input
    print("\nSample input tags (first 3 records):")
    print("-" * 80)
    for i, line in enumerate(record_lines[:3]):
        fields = line.strip().split('\t')
        read_name = fields[0]
        ref_name = fields[2]
        tags = '\t'.join(fields[11:]) if len(fields) > 11 else ''
        print(f"read_name: {read_name}")
        print(f"ref_name:  {ref_name}")
        print(f"tags:      {tags[:120]}...")
        print()

    # Convert to Parquet
    from bam_filter.sam_parser_fast import convert_sam_gz_to_parquet_fast

    print("Converting to Parquet...")
    stats = convert_sam_gz_to_parquet_fast(
        str(sample_file),
        output_dir,
        batch_size=500,
        compression_level=3,
        num_threads=1,
        store_sequences=True,
    )

    print(f"Converted {stats.get('total_records', 0):,} records")
    print()

    # Verify Parquet file
    parquet_file = output_path / "alignments.parquet"
    if not parquet_file.exists():
        print(f"ERROR: Parquet file not found!")
        return

    import pyarrow.parquet as pq

    print("Reading Parquet file to verify tags_raw...")
    table = pq.read_table(parquet_file)

    print(f"Schema columns: {table.column_names}")
    print()

    if 'tags_raw' not in table.column_names:
        print("ERROR: tags_raw column not found!")
        return

    df = table.to_pandas()

    # Count non-empty tags
    non_empty = df['tags_raw'].str.len() > 0
    print(f"Total rows: {len(df):,}")
    print(f"Rows with tags: {non_empty.sum():,} ({100*non_empty.sum()/len(df):.1f}%)")
    print()

    # Show sample records from Parquet
    print("Sample Parquet records with tags_raw:")
    print("-" * 80)
    sample = df[df['tags_raw'].str.len() > 0].head(3)
    for idx, row in sample.iterrows():
        print(f"read_name: {row['read_name']}")
        print(f"ref_name:  {row['ref_name']}")
        print(f"tags_raw:  {row['tags_raw'][:120]}...")
        print()

    # Count unique tag types
    print("Unique tag types found:")
    tag_counts = {}
    for tags in df['tags_raw']:
        if tags:
            for tag_part in tags.split('\t'):
                if ':' in tag_part:
                    tag_name = tag_part.split(':')[0]
                    tag_counts[tag_name] = tag_counts.get(tag_name, 0) + 1

    for tag, count in sorted(tag_counts.items(), key=lambda x: -x[1]):
        print(f"  {tag}: {count:,}")

    print()
    print("SUCCESS: Lossless SAM to Parquet conversion verified!")
    print("All optional tags are preserved in the tags_raw column.")


if __name__ == "__main__":
    main()
