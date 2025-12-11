#!/usr/bin/env python3
"""
Test optimized Parquet schema vs current schema.

Compares storage efficiency of:
1. Current schema (31 columns with duplication)
2. Optimized schema (no duplication, computed-on-read)
3. Normalized schema (separate reads/alignments tables)
"""

import sys
import gzip
import time
from pathlib import Path
import pyarrow as pa
import pyarrow.parquet as pq


def parse_sam_line(line):
    """Parse a SAM line into components."""
    fields = line.strip().split('\t')
    if len(fields) < 11:
        return None

    # Core fields
    read_name = fields[0]
    flag = int(fields[1])
    ref_name = fields[2]
    position = int(fields[3]) - 1  # Convert to 0-based
    mapq = int(fields[4])
    cigar = fields[5]
    mate_ref = fields[6]
    mate_pos = int(fields[7]) - 1 if fields[7] != '0' else -1
    template_len = int(fields[8])
    sequence = fields[9]
    quality = fields[10]

    # Parse optional tags
    tags = {}
    tags_cold = []
    hot_tags = {'AS', 'NM', 'XS', 'MD', 'XN', 'XM', 'XO', 'XG', 'YT'}

    for tag_str in fields[11:]:
        parts = tag_str.split(':')
        if len(parts) >= 3:
            tag_name = parts[0]
            tag_type = parts[1]
            tag_value = ':'.join(parts[2:])

            if tag_name in hot_tags:
                if tag_type == 'i':
                    tags[tag_name] = int(tag_value)
                elif tag_type == 'f':
                    tags[tag_name] = float(tag_value)
                else:
                    tags[tag_name] = tag_value
            else:
                tags_cold.append(tag_str)

    return {
        'read_name': read_name,
        'flag': flag,
        'ref_name': ref_name,
        'position': position,
        'mapq': mapq,
        'cigar': cigar,
        'mate_ref': mate_ref,
        'mate_pos': mate_pos,
        'template_len': template_len,
        'sequence': sequence,
        'quality': quality,
        'tags': tags,
        'tags_cold': '\t'.join(tags_cold) if tags_cold else '',
        'tags_raw': '\t'.join(fields[11:]) if len(fields) > 11 else '',
    }


def pack_sequence_2bit(seq):
    """Pack DNA sequence to 2-bit encoding (4x compression)."""
    encoding = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 0}
    packed = bytearray()
    for i in range(0, len(seq), 4):
        byte = 0
        for j in range(4):
            if i + j < len(seq):
                byte |= encoding.get(seq[i + j], 0) << (6 - 2 * j)
        packed.append(byte)
    return bytes(packed)


def write_current_schema(records, output_path, compression_level=3):
    """Write using current schema (with duplication)."""
    schema = pa.schema([
        ('read_id', pa.uint64()),
        ('read_name', pa.string()),
        ('ref_name', pa.string()),
        ('position', pa.int32()),
        ('end_position', pa.int32()),
        ('mapq', pa.uint8()),
        ('flag', pa.uint16()),
        ('alignment_length', pa.uint16()),
        ('cigar', pa.string()),
        ('sequence', pa.string()),
        ('quality', pa.string()),
        # Parsed tags (duplicated)
        ('alignment_score', pa.int32()),
        ('edit_distance', pa.int32()),
        ('xs_score', pa.int32()),
        ('md_string', pa.string()),
        # Raw tags (duplicated)
        ('tags_raw', pa.string()),
    ])

    data = {name: [] for name in schema.names}

    for i, rec in enumerate(records):
        data['read_id'].append(i)
        data['read_name'].append(rec['read_name'])
        data['ref_name'].append(rec['ref_name'])
        data['position'].append(rec['position'])
        data['end_position'].append(rec['position'] + len(rec['sequence']))
        data['mapq'].append(rec['mapq'])
        data['flag'].append(rec['flag'])
        data['alignment_length'].append(len(rec['sequence']))
        data['cigar'].append(rec['cigar'])
        data['sequence'].append(rec['sequence'])
        data['quality'].append(rec['quality'])
        data['alignment_score'].append(rec['tags'].get('AS', -1))
        data['edit_distance'].append(rec['tags'].get('NM', -1))
        data['xs_score'].append(rec['tags'].get('XS', -1))
        data['md_string'].append(rec['tags'].get('MD', ''))
        data['tags_raw'].append(rec['tags_raw'])

    table = pa.Table.from_pydict(data, schema=schema)
    pq.write_table(
        table, output_path,
        compression='zstd',
        compression_level=compression_level,
    )
    return output_path.stat().st_size


def write_optimized_schema(records, output_path, compression_level=6):
    """Write using optimized schema (no duplication, computed-on-read)."""
    schema = pa.schema([
        ('read_id', pa.uint64()),
        ('read_name', pa.string()),
        ('ref_id', pa.uint32()),  # Dictionary ID instead of string
        ('position', pa.int32()),
        # end_position: COMPUTED from cigar
        # alignment_length: COMPUTED from cigar
        ('mapq', pa.uint8()),
        ('flag', pa.uint16()),
        ('cigar', pa.string()),
        ('sequence', pa.binary()),  # 2-bit packed
        ('quality', pa.binary()),   # Raw bytes
        # Hot tags as columns (NOT in tags_cold)
        ('AS', pa.int16()),
        ('NM', pa.uint16()),
        ('XS', pa.int16()),
        ('MD', pa.string()),
        # Cold tags only (no duplication)
        ('tags_cold', pa.string()),
    ])

    # Build ref_name -> ref_id mapping
    ref_names = sorted(set(r['ref_name'] for r in records))
    ref_to_id = {name: i for i, name in enumerate(ref_names)}

    data = {name: [] for name in schema.names}

    for i, rec in enumerate(records):
        data['read_id'].append(i)
        data['read_name'].append(rec['read_name'])
        data['ref_id'].append(ref_to_id[rec['ref_name']])
        data['position'].append(rec['position'])
        data['mapq'].append(rec['mapq'])
        data['flag'].append(rec['flag'])
        data['cigar'].append(rec['cigar'])
        data['sequence'].append(pack_sequence_2bit(rec['sequence']))
        data['quality'].append(rec['quality'].encode('ascii'))
        data['AS'].append(rec['tags'].get('AS', -32768))
        data['NM'].append(rec['tags'].get('NM', 65535))
        data['XS'].append(rec['tags'].get('XS', -32768))
        data['MD'].append(rec['tags'].get('MD', ''))
        data['tags_cold'].append(rec['tags_cold'])

    table = pa.Table.from_pydict(data, schema=schema)
    pq.write_table(
        table, output_path,
        compression='zstd',
        compression_level=compression_level,
        use_dictionary=['ref_id', 'cigar'],
    )

    # Also write reference sidecar
    ref_schema = pa.schema([
        ('ref_id', pa.uint32()),
        ('ref_name', pa.string()),
    ])
    ref_data = {
        'ref_id': list(range(len(ref_names))),
        'ref_name': ref_names,
    }
    ref_table = pa.Table.from_pydict(ref_data, schema=ref_schema)
    ref_path = output_path.parent / 'references.parquet'
    pq.write_table(ref_table, ref_path, compression='zstd')

    return output_path.stat().st_size + ref_path.stat().st_size


def write_normalized_schema(records, output_dir, compression_level=6):
    """Write using normalized schema (separate reads/alignments)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(exist_ok=True)

    # Group alignments by read_name to find multi-mapped
    from collections import defaultdict
    reads_by_name = defaultdict(list)
    for rec in records:
        reads_by_name[rec['read_name']].append(rec)

    # Assign read_ids
    read_name_to_id = {name: i for i, name in enumerate(reads_by_name.keys())}

    # Build ref mapping
    ref_names = sorted(set(r['ref_name'] for r in records))
    ref_to_id = {name: i for i, name in enumerate(ref_names)}

    # Reads table (one row per unique read)
    reads_schema = pa.schema([
        ('read_id', pa.uint64()),
        ('read_name', pa.string()),
        ('sequence', pa.binary()),  # 2-bit packed
        ('quality', pa.binary()),
    ])

    reads_data = {name: [] for name in reads_schema.names}
    for read_name, alns in reads_by_name.items():
        reads_data['read_id'].append(read_name_to_id[read_name])
        reads_data['read_name'].append(read_name)
        # Use sequence/quality from first alignment
        reads_data['sequence'].append(pack_sequence_2bit(alns[0]['sequence']))
        reads_data['quality'].append(alns[0]['quality'].encode('ascii'))

    reads_table = pa.Table.from_pydict(reads_data, schema=reads_schema)
    reads_path = output_dir / 'reads.parquet'
    pq.write_table(reads_table, reads_path, compression='zstd', compression_level=compression_level)

    # Alignments table (one row per alignment)
    alns_schema = pa.schema([
        ('read_id', pa.uint64()),
        ('ref_id', pa.uint32()),
        ('position', pa.int32()),
        ('mapq', pa.uint8()),
        ('flag', pa.uint16()),
        ('cigar', pa.string()),
        ('AS', pa.int16()),
        ('NM', pa.uint16()),
        ('XS', pa.int16()),
        ('MD', pa.string()),
        ('tags_cold', pa.string()),
    ])

    alns_data = {name: [] for name in alns_schema.names}
    for rec in records:
        alns_data['read_id'].append(read_name_to_id[rec['read_name']])
        alns_data['ref_id'].append(ref_to_id[rec['ref_name']])
        alns_data['position'].append(rec['position'])
        alns_data['mapq'].append(rec['mapq'])
        alns_data['flag'].append(rec['flag'])
        alns_data['cigar'].append(rec['cigar'])
        alns_data['AS'].append(rec['tags'].get('AS', -32768))
        alns_data['NM'].append(rec['tags'].get('NM', 65535))
        alns_data['XS'].append(rec['tags'].get('XS', -32768))
        alns_data['MD'].append(rec['tags'].get('MD', ''))
        alns_data['tags_cold'].append(rec['tags_cold'])

    alns_table = pa.Table.from_pydict(alns_data, schema=alns_schema)
    alns_path = output_dir / 'alignments.parquet'
    pq.write_table(alns_table, alns_path, compression='zstd', compression_level=compression_level,
                   use_dictionary=['ref_id', 'cigar'])

    # References sidecar
    ref_schema = pa.schema([('ref_id', pa.uint32()), ('ref_name', pa.string())])
    ref_data = {'ref_id': list(range(len(ref_names))), 'ref_name': ref_names}
    ref_path = output_dir / 'references.parquet'
    pq.write_table(pa.Table.from_pydict(ref_data, schema=ref_schema), ref_path, compression='zstd')

    total_size = reads_path.stat().st_size + alns_path.stat().st_size + ref_path.stat().st_size

    return total_size, len(reads_by_name), len(records)


def main():
    if len(sys.argv) < 2:
        print("Usage: python test_optimized_schema.py <input.sam.gz> [num_records]")
        sys.exit(1)

    input_file = sys.argv[1]
    num_records = int(sys.argv[2]) if len(sys.argv) > 2 else 10000
    output_dir = Path("/scratch/tmp/schema_comparison")

    # Clean output
    import shutil
    if output_dir.exists():
        shutil.rmtree(output_dir)
    output_dir.mkdir(parents=True)

    print(f"Input: {input_file}")
    print(f"Records: {num_records:,}")
    print()

    # Parse records
    print("Parsing SAM records...")
    records = []
    with gzip.open(input_file, 'rt') as f:
        for line in f:
            if line.startswith('@'):
                continue
            rec = parse_sam_line(line)
            if rec:
                records.append(rec)
                if len(records) >= num_records:
                    break

    print(f"Parsed {len(records):,} records")

    # Count unique reads
    unique_reads = len(set(r['read_name'] for r in records))
    avg_alignments = len(records) / unique_reads
    print(f"Unique reads: {unique_reads:,}")
    print(f"Avg alignments/read: {avg_alignments:.1f}")
    print()

    # Calculate raw SAM size
    raw_sam_size = sum(
        len(r['read_name']) + len(r['ref_name']) + len(r['cigar']) +
        len(r['sequence']) + len(r['quality']) + len(r['tags_raw']) + 50
        for r in records
    )
    print(f"Raw SAM size: {raw_sam_size:,} bytes ({raw_sam_size/1024:.1f} KB)")
    print()

    # Test each schema
    print("=" * 60)
    print("SCHEMA COMPARISON")
    print("=" * 60)
    print()

    # Current schema
    current_path = output_dir / "current.parquet"
    current_size = write_current_schema(records, current_path)
    print(f"1. CURRENT SCHEMA (with duplication)")
    print(f"   Size: {current_size:,} bytes ({current_size/1024:.1f} KB)")
    print(f"   vs Raw SAM: {100*current_size/raw_sam_size:.1f}%")
    print()

    # Optimized schema
    optimized_path = output_dir / "optimized.parquet"
    optimized_size = write_optimized_schema(records, optimized_path)
    print(f"2. OPTIMIZED SCHEMA (no duplication, 2-bit seq)")
    print(f"   Size: {optimized_size:,} bytes ({optimized_size/1024:.1f} KB)")
    print(f"   vs Raw SAM: {100*optimized_size/raw_sam_size:.1f}%")
    print(f"   vs Current: {100*optimized_size/current_size:.1f}% ({100*(1-optimized_size/current_size):.0f}% smaller)")
    print()

    # Normalized schema
    normalized_dir = output_dir / "normalized"
    normalized_size, num_reads, num_alns = write_normalized_schema(records, normalized_dir)
    print(f"3. NORMALIZED SCHEMA (separate reads/alignments)")
    print(f"   Size: {normalized_size:,} bytes ({normalized_size/1024:.1f} KB)")
    print(f"   vs Raw SAM: {100*normalized_size/raw_sam_size:.1f}%")
    print(f"   vs Current: {100*normalized_size/current_size:.1f}% ({100*(1-normalized_size/current_size):.0f}% smaller)")
    print(f"   Reads table: {num_reads:,} rows")
    print(f"   Alignments table: {num_alns:,} rows")
    print()

    # Summary
    print("=" * 60)
    print("SUMMARY")
    print("=" * 60)
    print()
    print(f"{'Schema':<25} {'Size':>10} {'vs SAM':>10} {'vs Current':>12}")
    print("-" * 60)
    print(f"{'Raw SAM':<25} {raw_sam_size:>10,} {'100%':>10} {'-':>12}")
    print(f"{'Current (duplicated)':<25} {current_size:>10,} {100*current_size/raw_sam_size:>9.0f}% {'-':>12}")
    print(f"{'Optimized':<25} {optimized_size:>10,} {100*optimized_size/raw_sam_size:>9.0f}% {100*(1-optimized_size/current_size):>10.0f}%")
    print(f"{'Normalized':<25} {normalized_size:>10,} {100*normalized_size/raw_sam_size:>9.0f}% {100*(1-normalized_size/current_size):>10.0f}%")
    print()

    if avg_alignments > 1.5:
        print(f"NOTE: With {avg_alignments:.1f} alignments/read, normalized schema")
        print(f"      saves significant space by deduplicating sequences.")


if __name__ == "__main__":
    main()
