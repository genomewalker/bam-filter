# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

"""
Ultra-fast SAM parser in pure Cython - bypasses HTSlib text parsing.

SAM format: TAB-separated, 11 mandatory columns + optional tags
This parser reads BGZF-compressed SAM directly and parses in C.

Performance advantage over HTSlib sam_read1():
- HTSlib: decompress -> parse SAM text -> build bam1_t -> extract fields
- This: decompress -> parse directly to columnar arrays (skip bam1_t)

Supports:
- Single-threaded full file parsing
- Range-based parsing for parallel processing (BGZF block boundaries)
"""

from libc.stdlib cimport malloc, free, atoi, strtol, atof
from libc.string cimport strlen, strchr, memcpy, strncmp, strcmp
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool

import time
from pathlib import Path

# BGZF from HTSlib for parallel block decompression
cdef extern from "htslib/bgzf.h":
    ctypedef struct BGZF:
        pass

    BGZF* bgzf_open(const char* path, const char* mode) nogil
    int bgzf_close(BGZF* fp) nogil
    int64_t bgzf_read(BGZF* fp, void* data, size_t length) nogil
    int64_t bgzf_seek(BGZF* fp, int64_t pos, int whence) nogil
    int64_t bgzf_tell(BGZF* fp) nogil
    int bgzf_getline(BGZF* fp, int delim, void* str) nogil
    int bgzf_mt(BGZF* fp, int n_threads, int n_sub_blks) nogil

# kstring for line reading
cdef extern from "htslib/kstring.h":
    ctypedef struct kstring_t:
        size_t l
        size_t m
        char* s

    void ks_free(kstring_t* s) nogil

# Arrow C++
cdef extern from "arrow_parquet_writer.hpp" namespace "bam_filter":
    cdef cppclass AlignmentBatch:
        AlignmentBatch() except +
        size_t size()
        void clear()
        void reserve(size_t)

        vector[uint64_t] read_ids
        vector[string] read_names
        vector[uint32_t] ref_ids
        vector[string] ref_names
        vector[int32_t] positions
        vector[int32_t] end_positions
        vector[uint8_t] mapqs
        vector[uint16_t] flags
        vector[uint16_t] alignment_lengths
        vector[int32_t] template_lengths
        vector[int32_t] mate_ref_ids
        vector[int32_t] mate_positions
        vector[int32_t] alignment_scores
        vector[int32_t] xs_scores
        vector[int32_t] edit_distances
        vector[int32_t] num_mismatches
        vector[int32_t] num_gap_opens
        vector[int32_t] num_gap_extensions
        vector[string] md_strings
        vector[float] anis
        vector[float] pmd_scores
        vector[float] zs_scores
        vector[float] zp_posteriors
        vector[int32_t] lca_taxids
        vector[int32_t] reassigned_ref_ids
        vector[cbool] filter_passed
        vector[string] read_groups
        vector[string] cigars
        vector[string] sequences
        vector[string] qualities
        vector[string] tags_raw

    cdef cppclass ParquetWriter:
        ParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const AlignmentBatch& batch) except +
        void Close() except +


# Parsed SAM record structure
cdef struct ParsedSamRecord:
    # Mandatory fields (columns 1-11)
    const char* qname
    int qname_len
    uint16_t flag
    const char* rname
    int rname_len
    int32_t pos
    uint8_t mapq
    const char* cigar
    int cigar_len
    const char* rnext
    int rnext_len
    int32_t pnext
    int32_t tlen
    const char* seq
    int seq_len
    const char* qual
    int qual_len
    # ALL optional tags as raw string (for lossless round-trip)
    const char* tags_raw
    int tags_len
    # Common tags (extracted for fast querying)
    int32_t alignment_score  # AS:i:
    int32_t xs_score         # XS:i:
    int32_t edit_distance    # NM:i:
    int32_t num_mismatches   # XM:i:
    int32_t num_gap_opens    # XO:i:
    int32_t num_gap_ext      # XG:i:
    const char* md_str       # MD:Z:
    int md_len
    const char* rg_str       # RG:Z:
    int rg_len
    float ani
    # Validity
    cbool valid


cdef inline const char* find_tab(const char* s) noexcept nogil:
    """Find next tab character."""
    while s[0] != b'\t' and s[0] != b'\0' and s[0] != b'\n':
        s += 1
    return s


cdef inline int32_t parse_int_field(const char* start, const char* end) noexcept nogil:
    """Parse integer from field."""
    cdef char buf[32]
    cdef int length = <int>(end - start)
    if length >= 31:
        length = 30
    memcpy(buf, start, length)
    buf[length] = 0
    return <int32_t>atoi(buf)


cdef inline void parse_sam_line_nogil(
    const char* line,
    int line_len,
    ParsedSamRecord* rec
) noexcept nogil:
    """
    Parse a SAM line into record structure - pure C, nogil.

    SAM columns (TAB-separated):
    1. QNAME  2. FLAG  3. RNAME  4. POS  5. MAPQ  6. CIGAR
    7. RNEXT  8. PNEXT  9. TLEN  10. SEQ  11. QUAL  12+ TAGS
    """
    cdef const char* p = line
    cdef const char* field_start
    cdef const char* field_end
    cdef const char* tag_start
    cdef const char* line_end
    cdef int field_num = 0
    cdef char tag[3]
    cdef char tag_type

    # Initialize defaults
    rec.valid = False
    rec.tags_raw = NULL
    rec.tags_len = 0
    rec.alignment_score = -1
    rec.xs_score = -1
    rec.edit_distance = -1
    rec.num_mismatches = -1
    rec.num_gap_opens = -1
    rec.num_gap_ext = -1
    rec.md_str = NULL
    rec.md_len = 0
    rec.rg_str = NULL
    rec.rg_len = 0
    rec.ani = -1.0

    # Skip header lines
    if line[0] == b'@':
        return

    # Parse 11 mandatory fields
    while field_num < 11 and p[0] != b'\0' and p[0] != b'\n':
        field_start = p
        field_end = find_tab(p)

        if field_num == 0:  # QNAME
            rec.qname = field_start
            rec.qname_len = <int>(field_end - field_start)
        elif field_num == 1:  # FLAG
            rec.flag = <uint16_t>parse_int_field(field_start, field_end)
        elif field_num == 2:  # RNAME
            rec.rname = field_start
            rec.rname_len = <int>(field_end - field_start)
        elif field_num == 3:  # POS
            rec.pos = parse_int_field(field_start, field_end) - 1  # Convert to 0-based
        elif field_num == 4:  # MAPQ
            rec.mapq = <uint8_t>parse_int_field(field_start, field_end)
        elif field_num == 5:  # CIGAR
            rec.cigar = field_start
            rec.cigar_len = <int>(field_end - field_start)
        elif field_num == 6:  # RNEXT
            rec.rnext = field_start
            rec.rnext_len = <int>(field_end - field_start)
        elif field_num == 7:  # PNEXT
            rec.pnext = parse_int_field(field_start, field_end) - 1
        elif field_num == 8:  # TLEN
            rec.tlen = parse_int_field(field_start, field_end)
        elif field_num == 9:  # SEQ
            rec.seq = field_start
            rec.seq_len = <int>(field_end - field_start)
        elif field_num == 10:  # QUAL
            rec.qual = field_start
            rec.qual_len = <int>(field_end - field_start)

        field_num += 1
        p = field_end
        if p[0] == b'\t':
            p += 1

    if field_num < 11:
        return  # Invalid line

    rec.valid = True

    # Capture ALL optional tags as raw string for lossless round-trip
    # p now points to start of optional tags (after tab following QUAL)
    if p[0] != b'\0' and p[0] != b'\n':
        rec.tags_raw = p
        # Find end of line to get total tags length
        line_end = p
        while line_end[0] != b'\0' and line_end[0] != b'\n':
            line_end += 1
        rec.tags_len = <int>(line_end - p)

    # Also parse commonly-used tags for fast querying
    while p[0] != b'\0' and p[0] != b'\n':
        field_start = p
        field_end = find_tab(p)

        # Tags are TAG:TYPE:VALUE format
        if field_end - field_start >= 5 and field_start[2] == b':' and field_start[4] == b':':
            tag[0] = field_start[0]
            tag[1] = field_start[1]
            tag[2] = 0
            tag_type = field_start[3]
            tag_start = field_start + 5

            # AS:i: - Alignment score
            if tag[0] == b'A' and tag[1] == b'S' and tag_type == b'i':
                rec.alignment_score = parse_int_field(tag_start, field_end)
            # XS:i: - Suboptimal score
            elif tag[0] == b'X' and tag[1] == b'S' and tag_type == b'i':
                rec.xs_score = parse_int_field(tag_start, field_end)
            # NM:i: - Edit distance
            elif tag[0] == b'N' and tag[1] == b'M' and tag_type == b'i':
                rec.edit_distance = parse_int_field(tag_start, field_end)
                # Calculate ANI
                if rec.edit_distance >= 0 and rec.seq_len > 0:
                    rec.ani = (1.0 - <float>rec.edit_distance / <float>rec.seq_len) * 100.0
            # XM:i: - Mismatches
            elif tag[0] == b'X' and tag[1] == b'M' and tag_type == b'i':
                rec.num_mismatches = parse_int_field(tag_start, field_end)
            # XO:i: - Gap opens
            elif tag[0] == b'X' and tag[1] == b'O' and tag_type == b'i':
                rec.num_gap_opens = parse_int_field(tag_start, field_end)
            # XG:i: - Gap extensions
            elif tag[0] == b'X' and tag[1] == b'G' and tag_type == b'i':
                rec.num_gap_ext = parse_int_field(tag_start, field_end)
            # MD:Z: - MD string
            elif tag[0] == b'M' and tag[1] == b'D' and tag_type == b'Z':
                rec.md_str = tag_start
                rec.md_len = <int>(field_end - tag_start)
            # RG:Z: - Read group
            elif tag[0] == b'R' and tag[1] == b'G' and tag_type == b'Z':
                rec.rg_str = tag_start
                rec.rg_len = <int>(field_end - tag_start)

        p = field_end
        if p[0] == b'\t':
            p += 1


def convert_sam_gz_to_parquet_fast(
    str input_file,
    str output_dir,
    int batch_size = 100000,
    int compression_level = 3,
    int num_threads = 8,
    bint store_sequences = False,
):
    """
    Ultra-fast SAM.gz to Parquet converter using custom Cython parser.

    Bypasses HTSlib's sam_read1() text parsing overhead.
    Uses BGZF multi-threaded decompression + direct line parsing.

    Args:
        input_file: Path to SAM.gz file (BGZF compressed)
        output_dir: Output directory for Parquet
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads: Threads for BGZF decompression
        store_sequences: Store SEQ/QUAL strings

    Returns:
        Dictionary with statistics
    """
    cdef BGZF* fp = NULL
    cdef kstring_t line
    cdef int ret
    cdef uint64_t total_records = 0
    cdef uint64_t header_lines = 0
    cdef double start_time = time.time()
    cdef double last_print = start_time
    cdef ParsedSamRecord rec
    cdef cbool do_store_seq = store_sequences

    # C++ batch and writer
    cdef ParquetWriter* writer = NULL
    cdef AlignmentBatch* batch = new AlignmentBatch()

    # Initialize kstring
    line.l = 0
    line.m = 0
    line.s = NULL

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Open BGZF file
    input_bytes = input_file.encode('utf-8')
    fp = bgzf_open(input_bytes, b"r")
    if fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    # Enable multi-threaded decompression
    if num_threads > 1:
        bgzf_mt(fp, num_threads, 256)

    # Create writer
    output_file = str(output_path / "alignments.parquet")
    writer = new ParquetWriter(string(output_file.encode('utf-8')), compression_level)

    print(f"Processing {input_file} with fast Cython parser...")
    print(f"Threads: {num_threads}, Batch size: {batch_size:,}")

    try:
        while True:
            # Read line using BGZF
            with nogil:
                ret = bgzf_getline(fp, b'\n', &line)

            if ret < 0:
                break

            # Skip header lines
            if line.s[0] == b'@':
                header_lines += 1
                continue

            # Parse line in nogil block
            with nogil:
                parse_sam_line_nogil(line.s, <int>line.l, &rec)

            if not rec.valid:
                continue

            # Skip unmapped (RNAME = "*" or FLAG & 4)
            if rec.rname_len == 1 and rec.rname[0] == b'*':
                continue
            if rec.flag & 4:
                continue

            total_records += 1

            # Append to batch
            batch.read_ids.push_back(total_records)
            batch.read_names.push_back(string(rec.qname, rec.qname_len))
            batch.ref_ids.push_back(0)  # No ref_id in SAM text
            batch.ref_names.push_back(string(rec.rname, rec.rname_len))
            batch.positions.push_back(rec.pos)
            batch.end_positions.push_back(rec.pos + rec.seq_len)  # Approximate
            batch.mapqs.push_back(rec.mapq)
            batch.flags.push_back(rec.flag)
            batch.alignment_lengths.push_back(rec.seq_len)
            batch.template_lengths.push_back(rec.tlen)
            batch.mate_ref_ids.push_back(-1)
            batch.mate_positions.push_back(rec.pnext if rec.pnext >= 0 else -1)
            batch.alignment_scores.push_back(rec.alignment_score)
            batch.xs_scores.push_back(rec.xs_score)
            batch.edit_distances.push_back(rec.edit_distance)
            batch.num_mismatches.push_back(rec.num_mismatches)
            batch.num_gap_opens.push_back(rec.num_gap_opens)
            batch.num_gap_extensions.push_back(rec.num_gap_ext)

            if rec.md_str != NULL:
                batch.md_strings.push_back(string(rec.md_str, rec.md_len))
            else:
                batch.md_strings.push_back(string())

            batch.anis.push_back(rec.ani)
            batch.pmd_scores.push_back(-1.0)
            batch.zs_scores.push_back(-1.0)
            batch.zp_posteriors.push_back(-1.0)
            batch.lca_taxids.push_back(-1)
            batch.reassigned_ref_ids.push_back(-1)
            batch.filter_passed.push_back(True)

            if rec.rg_str != NULL:
                batch.read_groups.push_back(string(rec.rg_str, rec.rg_len))
            else:
                batch.read_groups.push_back(string())

            batch.cigars.push_back(string(rec.cigar, rec.cigar_len))

            if do_store_seq:
                batch.sequences.push_back(string(rec.seq, rec.seq_len))
                batch.qualities.push_back(string(rec.qual, rec.qual_len))
            else:
                batch.sequences.push_back(string())
                batch.qualities.push_back(string())

            # Store ALL tags for lossless round-trip
            if rec.tags_raw != NULL and rec.tags_len > 0:
                batch.tags_raw.push_back(string(rec.tags_raw, rec.tags_len))
            else:
                batch.tags_raw.push_back(string())

            # Flush batch
            if batch.size() >= batch_size:
                writer.WriteBatch(batch[0])
                batch.clear()

            # Progress
            if time.time() - last_print >= 10:
                elapsed = time.time() - start_time
                rate = total_records / elapsed if elapsed > 0 else 0
                print(f"Progress: {total_records:,} records ({rate/1e6:.2f} M/s)", flush=True)
                last_print = time.time()

        # Flush remaining
        if batch.size() > 0:
            writer.WriteBatch(batch[0])

        writer.Close()

        duration = time.time() - start_time

        print(f"\nDone! Header lines: {header_lines:,}")
        print(f"Records: {total_records:,}")
        print(f"Time: {duration:.1f}s")
        print(f"Throughput: {total_records/duration:,.0f} rec/s ({total_records/duration/1e6:.3f} M/s)")

        return {
            'total_records': total_records,
            'header_lines': header_lines,
            'processing_time_seconds': duration,
        }

    finally:
        ks_free(&line)
        del batch
        if writer != NULL:
            del writer
        if fp != NULL:
            bgzf_close(fp)


def convert_sam_gz_range_to_parquet_fast(
    str input_file,
    str output_dir,
    int64_t start_offset,
    int64_t end_offset,
    int batch_size = 100000,
    int compression_level = 3,
    int num_threads = 1,
    bint store_sequences = False,
):
    """
    Convert a specific BGZF byte range of SAM.gz to Parquet.

    For parallel processing - each worker processes a different byte range.
    Uses the fast Cython parser that bypasses HTSlib text parsing.

    Args:
        input_file: Path to SAM.gz file (BGZF compressed)
        output_dir: Output directory for Parquet
        start_offset: Starting BGZF block offset in bytes
        end_offset: Ending BGZF block offset in bytes
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads: Threads for BGZF decompression (per worker)
        store_sequences: Store SEQ/QUAL strings

    Returns:
        Dictionary with statistics
    """
    cdef BGZF* fp = NULL
    cdef kstring_t line
    cdef int ret
    cdef uint64_t total_records = 0
    cdef uint64_t skipped_lines = 0
    cdef double start_time = time.time()
    cdef ParsedSamRecord rec
    cdef cbool do_store_seq = store_sequences
    cdef int64_t current_block_offset
    cdef int64_t virtual_offset

    # C++ batch and writer
    cdef ParquetWriter* writer = NULL
    cdef AlignmentBatch* batch = new AlignmentBatch()

    # Initialize kstring
    line.l = 0
    line.m = 0
    line.s = NULL

    # Create output directory
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    # Open BGZF file
    input_bytes = input_file.encode('utf-8')
    fp = bgzf_open(input_bytes, b"r")
    if fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    # Enable multi-threaded decompression
    if num_threads > 1:
        bgzf_mt(fp, num_threads, 256)

    # Seek to start offset (virtual offset = block_offset << 16)
    if start_offset > 0:
        virtual_offset = start_offset << 16
        if bgzf_seek(fp, virtual_offset, 0) < 0:
            bgzf_close(fp)
            raise IOError(f"Failed to seek to offset {start_offset}")

        # After seeking to a block boundary, we might be in the middle of a line
        # Read and discard the partial line (if any)
        with nogil:
            ret = bgzf_getline(fp, b'\n', &line)
        if ret >= 0:
            skipped_lines += 1

    # Create writer
    output_file = str(output_path / "alignments.parquet")
    writer = new ParquetWriter(string(output_file.encode('utf-8')), compression_level)

    try:
        while True:
            # Check current position before reading
            current_block_offset = bgzf_tell(fp) >> 16
            if current_block_offset >= end_offset:
                break

            # Read line using BGZF
            with nogil:
                ret = bgzf_getline(fp, b'\n', &line)

            if ret < 0:
                break

            # Skip header lines (shouldn't be many after seeking)
            if line.s[0] == b'@':
                skipped_lines += 1
                continue

            # Parse line in nogil block
            with nogil:
                parse_sam_line_nogil(line.s, <int>line.l, &rec)

            if not rec.valid:
                skipped_lines += 1
                continue

            # Skip unmapped
            if rec.rname_len == 1 and rec.rname[0] == b'*':
                continue
            if rec.flag & 4:
                continue

            total_records += 1

            # Append to batch
            batch.read_ids.push_back(total_records)
            batch.read_names.push_back(string(rec.qname, rec.qname_len))
            batch.ref_ids.push_back(0)
            batch.ref_names.push_back(string(rec.rname, rec.rname_len))
            batch.positions.push_back(rec.pos)
            batch.end_positions.push_back(rec.pos + rec.seq_len)
            batch.mapqs.push_back(rec.mapq)
            batch.flags.push_back(rec.flag)
            batch.alignment_lengths.push_back(rec.seq_len)
            batch.template_lengths.push_back(rec.tlen)
            batch.mate_ref_ids.push_back(-1)
            batch.mate_positions.push_back(rec.pnext if rec.pnext >= 0 else -1)
            batch.alignment_scores.push_back(rec.alignment_score)
            batch.xs_scores.push_back(rec.xs_score)
            batch.edit_distances.push_back(rec.edit_distance)
            batch.num_mismatches.push_back(rec.num_mismatches)
            batch.num_gap_opens.push_back(rec.num_gap_opens)
            batch.num_gap_extensions.push_back(rec.num_gap_ext)

            if rec.md_str != NULL:
                batch.md_strings.push_back(string(rec.md_str, rec.md_len))
            else:
                batch.md_strings.push_back(string())

            batch.anis.push_back(rec.ani)
            batch.pmd_scores.push_back(-1.0)
            batch.zs_scores.push_back(-1.0)
            batch.zp_posteriors.push_back(-1.0)
            batch.lca_taxids.push_back(-1)
            batch.reassigned_ref_ids.push_back(-1)
            batch.filter_passed.push_back(True)

            if rec.rg_str != NULL:
                batch.read_groups.push_back(string(rec.rg_str, rec.rg_len))
            else:
                batch.read_groups.push_back(string())

            batch.cigars.push_back(string(rec.cigar, rec.cigar_len))

            if do_store_seq:
                batch.sequences.push_back(string(rec.seq, rec.seq_len))
                batch.qualities.push_back(string(rec.qual, rec.qual_len))
            else:
                batch.sequences.push_back(string())
                batch.qualities.push_back(string())

            # Store ALL optional tags for lossless round-trip
            if rec.tags_raw != NULL:
                batch.tags_raw.push_back(string(rec.tags_raw, rec.tags_len))
            else:
                batch.tags_raw.push_back(string())

            # Flush batch
            if batch.size() >= batch_size:
                writer.WriteBatch(batch[0])
                batch.clear()

        # Flush remaining
        if batch.size() > 0:
            writer.WriteBatch(batch[0])

        writer.Close()

        return {
            'total_records': total_records,
            'skipped_lines': skipped_lines,
            'processing_time_seconds': time.time() - start_time,
        }

    finally:
        ks_free(&line)
        del batch
        if writer != NULL:
            del writer
        if fp != NULL:
            bgzf_close(fp)
