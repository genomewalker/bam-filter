# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

"""
Optimized Parquet converter using HTSlib + optimized C++ Arrow writer.

Key optimizations:
1. LUT-based 2-bit sequence packing (4x compression)
2. Hot/cold tag separation (no duplication)
3. Bulk Arrow AppendValues API
4. Zero-copy string handling
5. Pre-reserved batch storage
"""

from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memcpy, strcmp, strchr
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int16_t, int32_t, int64_t
from libc.stdio cimport sprintf
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cpython cimport bool as py_bool

import time
from pathlib import Path
import os

# Import PMD calculation from processor_md_quality
from bam_filter.processor_md_quality cimport (
    calculate_md_quality_score,
    initialize_quality_lookup_tables
)
from bam_filter.processor cimport AlignmentScoringConfig

# HTSlib declarations
cdef extern from "htslib/sam.h":
    ctypedef struct sam_hdr_t:
        pass

    ctypedef struct bam_hdr_t:
        pass

    ctypedef struct bam1_core_t:
        int32_t tid
        int32_t pos
        uint16_t bin
        uint8_t qual
        uint8_t l_qname
        uint16_t flag
        uint32_t l_extranul
        uint16_t n_cigar
        int32_t l_qseq
        int32_t mtid
        int32_t mpos
        int32_t isize

    ctypedef struct bam1_t:
        bam1_core_t core
        uint64_t id
        uint8_t *data
        int l_data
        uint32_t m_data

    ctypedef struct samFile:
        pass

    samFile* sam_open(const char *fn, const char *mode)
    int sam_close(samFile *fp)
    bam_hdr_t* sam_hdr_read(samFile *fp)
    void bam_hdr_destroy(bam_hdr_t *h)
    bam1_t* bam_init1()
    void bam_destroy1(bam1_t *b)
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b)

    uint8_t* bam_get_seq(bam1_t *b) nogil
    uint8_t* bam_get_qual(bam1_t *b) nogil
    uint32_t* bam_get_cigar(bam1_t *b) nogil
    char* bam_get_qname(bam1_t *b) nogil

    uint8_t bam_seqi(uint8_t *s, int i) nogil

    uint8_t* bam_aux_get(bam1_t *b, const char tag[2]) nogil
    int32_t bam_aux2i(uint8_t *s) nogil
    float bam_aux2f(uint8_t *s) nogil
    char* bam_aux2Z(uint8_t *s) nogil
    uint8_t* bam_aux_first(bam1_t *b) nogil
    uint8_t* bam_aux_next(bam1_t *b, uint8_t *s) nogil

    char* sam_hdr_tid2name(bam_hdr_t *h, int tid) nogil
    int bam_endpos(bam1_t *b) nogil

    int sam_hdr_write(samFile *fp, const bam_hdr_t *h)
    char* sam_hdr_str(const bam_hdr_t *h)
    int sam_hdr_nref(const bam_hdr_t *h)
    int hts_set_threads(samFile *fp, int n)

# BGZF support
cdef extern from "htslib/bgzf.h":
    ctypedef struct BGZF:
        pass
    int64_t bgzf_seek(BGZF *fp, int64_t pos, int whence) nogil
    int64_t bgzf_tell(BGZF *fp) nogil

cdef extern from "htslib/tbx.h":
    BGZF* hts_get_bgzfp(htsFile *fp)

cdef extern from "htslib/hts.h":
    ctypedef struct htsFile:
        pass

# Optimized Arrow C++ writer
cdef extern from "arrow_parquet_writer_optimized.hpp" namespace "bam_filter":
    # LUT-based sequence packing
    string pack_sequence_2bit(const char* seq, size_t length)
    string unpack_sequence_2bit(const string& packed, size_t original_len)

    cdef cppclass OptimizedAlignmentBatch:
        OptimizedAlignmentBatch() except +
        size_t size()
        void clear()
        void reserve(size_t)

        # Core fields
        vector[uint64_t] read_ids
        vector[string] read_names
        vector[uint32_t] ref_ids
        vector[int32_t] positions
        vector[uint8_t] mapqs
        vector[uint16_t] flags
        vector[string] cigars
        vector[int32_t] template_lengths
        vector[int32_t] mate_ref_ids
        vector[int32_t] mate_positions

        # Sequences - 2-bit packed
        vector[string] sequences_packed
        vector[uint16_t] sequence_lengths
        vector[string] qualities

        # Hot tags as columns
        vector[int16_t] AS
        vector[uint16_t] NM
        vector[int16_t] XS
        vector[string] MD

        # Cold tags (hot tags stripped)
        vector[string] tags_cold

        # Derived fields (precomputed for fast queries)
        vector[uint16_t] aligned_lengths  # reference span from CIGAR (compresses well)
        vector[cbool] is_reverse          # flag & 0x10

        # Pipeline results
        vector[float] ani
        vector[float] zs_score           # Log-likelihood alignment score (ZS:f tag)
        vector[float] pmd_score          # Post-mortem damage score (PM:f tag)
        vector[cbool] filter_passed
        vector[int32_t] lca_taxid
        vector[int32_t] reassigned_ref_id
        vector[float] zp_posterior

        # Optional
        vector[string] read_groups

    cdef cppclass OptimizedParquetWriter:
        OptimizedParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const OptimizedAlignmentBatch& batch) except +
        void Close() except +

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(const char*)
        string(const char*, size_t)


# DNA encoding lookup
cdef char seq_nt16_str[16]
seq_nt16_str = [b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
                b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N']


cdef inline void extract_sequence(bam1_t *aln, char* result) nogil:
    """Extract DNA sequence - pure C."""
    cdef uint8_t *seq = bam_get_seq(aln)
    cdef int i
    cdef int l_seq = aln.core.l_qseq

    for i in range(l_seq):
        result[i] = seq_nt16_str[bam_seqi(seq, i)]
    result[l_seq] = 0


cdef inline void extract_cigar(bam1_t *aln, char* result) nogil:
    """Extract CIGAR string - pure C."""
    cdef uint32_t *cigar = bam_get_cigar(aln)
    cdef int n_cigar = aln.core.n_cigar
    cdef int i, pos = 0
    cdef uint32_t op
    cdef char op_char

    if n_cigar == 0:
        result[0] = b'*'
        result[1] = 0
        return

    for i in range(n_cigar):
        op = cigar[i]
        pos += sprintf(&result[pos], "%d", op >> 4)
        op_char = "MIDNSHP=X"[op & 0xf]
        result[pos] = op_char
        pos += 1

    result[pos] = 0


# Pre-extracted record for nogil processing
cdef struct ExtractedRecord:
    int32_t tid
    int32_t pos
    uint8_t mapq
    uint16_t flag
    int32_t l_qseq
    int32_t isize
    int32_t mtid
    int32_t mpos
    # Hot tags
    int16_t alignment_score  # AS
    uint16_t edit_distance   # NM
    int16_t xs_score         # XS
    # Derived fields
    uint16_t aligned_length  # reference span from CIGAR (compresses better than end_position)
    cbool is_reverse         # flag & 0x10
    # Computed
    float ani
    float zs_score           # Log-likelihood alignment score (ZS:f)
    float pmd_score          # Post-mortem damage score (PM:f)
    int32_t lca_taxid
    int32_t reassigned_ref_id
    float zp_posterior
    # Pointers
    const char* qname
    const char* ref_name
    const char* md_str
    const char* rg_str


cdef inline void extract_record_optimized(
    bam1_t* aln,
    bam_hdr_t* header,
    ExtractedRecord* rec,
    char* seq_buf,
    char* cigar_buf,
    char* tags_cold_buf,
    AlignmentScoringConfig* pmd_config,
    cbool calculate_pmd,
    cbool extract_sequence_flag
) noexcept nogil:
    """Extract record data with hot/cold tag separation."""
    cdef uint8_t* aux
    cdef uint8_t* aux_iter
    cdef float pmd_result = -1.0
    cdef char tag_name[3]
    cdef int tags_pos = 0
    cdef char type_char
    cdef int32_t int_val
    cdef float float_val
    cdef const char* str_val
    cdef cbool is_hot_tag

    # Core fields
    rec.tid = aln.core.tid
    rec.pos = aln.core.pos
    rec.mapq = aln.core.qual
    rec.flag = aln.core.flag
    rec.l_qseq = aln.core.l_qseq
    rec.isize = aln.core.isize
    rec.mtid = aln.core.mtid
    rec.mpos = aln.core.mpos

    # Derived fields (precomputed for fast queries)
    rec.aligned_length = <uint16_t>(bam_endpos(aln) - rec.pos)  # reference span from CIGAR
    rec.is_reverse = (rec.flag & 0x10) != 0                     # reverse strand flag

    # Names
    rec.qname = bam_get_qname(aln)
    rec.ref_name = sam_hdr_tid2name(header, aln.core.tid)
    if rec.ref_name == NULL:
        rec.ref_name = b"*"

    # Extract sequence and CIGAR
    if extract_sequence_flag:
        extract_sequence(aln, seq_buf)
    extract_cigar(aln, cigar_buf)

    # Initialize hot tags with defaults
    rec.alignment_score = -1
    rec.edit_distance = 65535  # Sentinel for missing
    rec.xs_score = -1
    rec.md_str = ""
    rec.rg_str = ""
    rec.lca_taxid = -1
    rec.reassigned_ref_id = -1
    rec.zp_posterior = -1.0

    # Extract hot tags
    aux = bam_aux_get(aln, b"AS")
    if aux != NULL:
        rec.alignment_score = <int16_t>bam_aux2i(aux)

    aux = bam_aux_get(aln, b"NM")
    if aux != NULL:
        rec.edit_distance = <uint16_t>bam_aux2i(aux)

    aux = bam_aux_get(aln, b"XS")
    if aux != NULL:
        rec.xs_score = <int16_t>bam_aux2i(aux)

    aux = bam_aux_get(aln, b"MD")
    if aux != NULL:
        rec.md_str = bam_aux2Z(aux)

    aux = bam_aux_get(aln, b"RG")
    if aux != NULL:
        rec.rg_str = bam_aux2Z(aux)

    # Pipeline tags
    aux = bam_aux_get(aln, b"ZT")
    if aux != NULL:
        rec.lca_taxid = bam_aux2i(aux)

    aux = bam_aux_get(aln, b"ZR")
    if aux != NULL:
        rec.reassigned_ref_id = bam_aux2i(aux)

    aux = bam_aux_get(aln, b"ZP")
    if aux != NULL:
        rec.zp_posterior = bam_aux2f(aux)

    # ANI calculation
    if rec.edit_distance < 65535 and rec.l_qseq > 0:
        rec.ani = (1.0 - <float>rec.edit_distance / <float>rec.l_qseq) * 100.0
    else:
        rec.ani = -1.0

    # ZS score calculation (log-likelihood alignment score)
    # Try to read from existing ZS tag first, otherwise calculate from MD
    aux = bam_aux_get(aln, b"ZS")
    if aux != NULL:
        rec.zs_score = bam_aux2f(aux)
    elif rec.md_str[0] != 0:
        # Calculate ZS from MD tag using fast path (no PMD output)
        rec.zs_score = <float>calculate_md_quality_score(aln, <sam_hdr_t*>header, pmd_config, <float*>NULL)
    else:
        rec.zs_score = -1e20

    # PMD calculation (post-mortem damage score)
    if calculate_pmd and rec.md_str[0] != 0:
        calculate_md_quality_score(aln, <sam_hdr_t*>header, pmd_config, &pmd_result)
        rec.pmd_score = pmd_result
    else:
        aux = bam_aux_get(aln, b"PM")
        if aux == NULL:
            aux = bam_aux_get(aln, b"PMD")
        rec.pmd_score = bam_aux2f(aux) if aux != NULL else -1.0

    # Build cold tags string (all tags except AS, NM, XS, MD, RG and pipeline tags)
    tags_cold_buf[0] = 0
    tags_pos = 0
    aux_iter = bam_aux_first(aln)

    while aux_iter != NULL:
        # Get tag name (2 chars before the type byte)
        tag_name[0] = <char>(aux_iter[-2])
        tag_name[1] = <char>(aux_iter[-1])
        tag_name[2] = 0

        # Skip hot tags and pipeline tags
        is_hot_tag = (
            (tag_name[0] == b'A' and tag_name[1] == b'S') or
            (tag_name[0] == b'N' and tag_name[1] == b'M') or
            (tag_name[0] == b'X' and tag_name[1] == b'S') or
            (tag_name[0] == b'M' and tag_name[1] == b'D') or
            (tag_name[0] == b'R' and tag_name[1] == b'G') or
            (tag_name[0] == b'Z' and tag_name[1] == b'T') or
            (tag_name[0] == b'Z' and tag_name[1] == b'R') or
            (tag_name[0] == b'Z' and tag_name[1] == b'P') or
            (tag_name[0] == b'Z' and tag_name[1] == b'S') or  # ZS score
            (tag_name[0] == b'P' and tag_name[1] == b'M')
        )

        if not is_hot_tag:
            type_char = <char>aux_iter[0]

            if tags_pos > 0:
                tags_cold_buf[tags_pos] = b'\t'
                tags_pos += 1

            # Add tag name
            tags_cold_buf[tags_pos] = tag_name[0]
            tags_cold_buf[tags_pos + 1] = tag_name[1]
            tags_cold_buf[tags_pos + 2] = b':'
            tags_cold_buf[tags_pos + 3] = type_char
            tags_cold_buf[tags_pos + 4] = b':'
            tags_pos += 5

            # Add value based on type
            if type_char == b'i' or type_char == b'I' or type_char == b's' or type_char == b'S' or type_char == b'c' or type_char == b'C':
                int_val = bam_aux2i(aux_iter)
                tags_pos += sprintf(&tags_cold_buf[tags_pos], "%d", int_val)
            elif type_char == b'f' or type_char == b'd':
                float_val = bam_aux2f(aux_iter)
                tags_pos += sprintf(&tags_cold_buf[tags_pos], "%.6g", float_val)
            elif type_char == b'Z' or type_char == b'H':
                str_val = bam_aux2Z(aux_iter)
                while str_val[0] != 0:
                    tags_cold_buf[tags_pos] = str_val[0]
                    tags_pos += 1
                    str_val += 1
            elif type_char == b'A':
                tags_cold_buf[tags_pos] = <char>aux_iter[1]
                tags_pos += 1

        aux_iter = bam_aux_next(aln, aux_iter)

    tags_cold_buf[tags_pos] = 0


def convert_bam_to_optimized_parquet(
    str input_file,
    str output_file,
    int batch_size = 100000,
    int compression_level = 6,
    int num_threads = 4,
    py_bool calculate_pmd = False,
    str library_type = "ds",
    py_bool store_sequences = True,
):
    """
    Convert BAM/SAM to optimized Parquet format.

    Uses the storage-efficient schema with:
    - 2-bit packed sequences (4x compression)
    - Hot/cold tag separation (no duplication)
    - Bulk Arrow AppendValues API

    Args:
        input_file: Input BAM/SAM file
        output_file: Output Parquet file
        batch_size: Records per batch (default: 100000)
        compression_level: ZSTD compression level (default: 6)
        num_threads: HTSlib decompression threads
        calculate_pmd: Calculate PMD scores on-the-fly
        library_type: "ds" (double-stranded) or "ss" (single-stranded)
        store_sequences: Store sequences (default: True)

    Returns:
        Dictionary with total_records and processing_time_seconds
    """
    cdef samFile *sam_fp
    cdef bam_hdr_t *header
    cdef bam1_t *aln
    cdef int ret
    cdef uint64_t total_records = 0
    cdef double start_time = time.time()
    cdef double last_print = start_time
    cdef char seq_buf[100000]
    cdef char cigar_buf[10000]
    cdef char tags_cold_buf[50000]
    cdef ExtractedRecord rec
    cdef cbool do_pmd = calculate_pmd
    cdef cbool do_store_seq = store_sequences
    cdef string packed_seq

    # PMD config
    cdef AlignmentScoringConfig pmd_config
    pmd_config.calculate_pmd = calculate_pmd
    pmd_config.is_single_stranded = (library_type == "ss")
    pmd_config.minimum_read_identity = 0.0
    pmd_config.minimum_read_length = 0
    pmd_config.maximum_read_length = 100000
    pmd_config.global_min_score = -1e20
    pmd_config.global_max_score = 1e20

    # Always initialize lookup tables for ZS score calculation
    initialize_quality_lookup_tables()

    # Create output directory
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # C++ objects
    cdef OptimizedParquetWriter* writer = NULL
    cdef OptimizedAlignmentBatch* batch = new OptimizedAlignmentBatch()
    batch.reserve(batch_size)

    # Open input file
    input_bytes = input_file.encode('utf-8')
    sam_fp = sam_open(input_bytes, b"r")
    if sam_fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    if num_threads > 1:
        hts_set_threads(sam_fp, num_threads - 1)

    header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        raise IOError(f"Failed to read header from {input_file}")

    aln = bam_init1()

    # Create writer
    output_bytes = output_file.encode('utf-8')
    writer = new OptimizedParquetWriter(string(output_bytes), compression_level)

    print(f"Converting {input_file} to optimized Parquet...")
    print(f"Output: {output_file}")
    print(f"Batch size: {batch_size:,}, Compression: ZSTD level {compression_level}")

    try:
        while True:
            ret = sam_read1(sam_fp, header, aln)
            if ret < 0:
                break

            # Skip unmapped
            if aln.core.tid < 0:
                continue

            # Extract record with hot/cold separation
            with nogil:
                extract_record_optimized(
                    aln, header, &rec, seq_buf, cigar_buf, tags_cold_buf,
                    &pmd_config, do_pmd, do_store_seq
                )

            # Append to batch
            batch.read_ids.push_back(total_records)
            batch.read_names.push_back(string(rec.qname))
            batch.ref_ids.push_back(<uint32_t>rec.tid)
            batch.positions.push_back(rec.pos)
            batch.mapqs.push_back(rec.mapq)
            batch.flags.push_back(rec.flag)
            batch.cigars.push_back(string(<const char*>cigar_buf))
            batch.template_lengths.push_back(rec.isize)
            batch.mate_ref_ids.push_back(rec.mtid)
            batch.mate_positions.push_back(rec.mpos)

            # Sequences - 2-bit packed
            if do_store_seq:
                packed_seq = pack_sequence_2bit(<const char*>seq_buf, <size_t>rec.l_qseq)
                batch.sequences_packed.push_back(packed_seq)
                batch.sequence_lengths.push_back(<uint16_t>rec.l_qseq)
                batch.qualities.push_back(string(<char*>bam_get_qual(aln), rec.l_qseq))
            else:
                batch.sequences_packed.push_back(string())
                batch.sequence_lengths.push_back(0)
                batch.qualities.push_back(string())

            # Hot tags
            batch.AS.push_back(rec.alignment_score)
            batch.NM.push_back(rec.edit_distance)
            batch.XS.push_back(rec.xs_score)
            batch.MD.push_back(string(rec.md_str))

            # Cold tags (hot tags stripped)
            batch.tags_cold.push_back(string(<const char*>tags_cold_buf))

            # Derived fields (precomputed for fast queries)
            batch.aligned_lengths.push_back(rec.aligned_length)
            batch.is_reverse.push_back(rec.is_reverse)

            # Pipeline results
            batch.ani.push_back(rec.ani)
            batch.zs_score.push_back(rec.zs_score)
            batch.pmd_score.push_back(rec.pmd_score)
            batch.filter_passed.push_back(True)
            batch.lca_taxid.push_back(rec.lca_taxid)
            batch.reassigned_ref_id.push_back(rec.reassigned_ref_id)
            batch.zp_posterior.push_back(rec.zp_posterior)

            # Read group
            batch.read_groups.push_back(string(rec.rg_str))

            total_records += 1

            # Flush batch
            if batch.size() >= batch_size:
                writer.WriteBatch(batch[0])
                batch.clear()
                batch.reserve(batch_size)

            # Progress
            if time.time() - last_print >= 10:
                elapsed = time.time() - start_time
                rate = total_records / elapsed / 1e6
                print(f"\rProgress: {total_records:,} records ({rate:.2f} M/sec)", end='', flush=True)
                last_print = time.time()

        # Flush remaining
        if batch.size() > 0:
            writer.WriteBatch(batch[0])

        writer.Close()

        duration = time.time() - start_time
        rate = total_records / duration if duration > 0 else 0

        print(f"\n\nCompleted: {total_records:,} records in {duration:.1f}s ({rate/1e6:.2f} M/sec)")

        # File size
        file_size = os.path.getsize(output_file)
        bytes_per_record = file_size / total_records if total_records > 0 else 0
        print(f"Output size: {file_size / 1024 / 1024:.1f} MB ({bytes_per_record:.1f} bytes/record)")

        return {
            'total_records': total_records,
            'processing_time_seconds': duration,
            'output_size_bytes': file_size,
        }

    finally:
        del batch
        if writer != NULL:
            del writer
        bam_destroy1(aln)
        bam_hdr_destroy(header)
        sam_close(sam_fp)


def convert_bgzf_range_to_optimized_parquet(
    str input_file,
    str output_file,
    int64_t start_offset,
    int64_t end_offset,
    int batch_size = 100000,
    int compression_level = 6,
    int num_threads = 1,
    py_bool calculate_pmd = False,
    str library_type = "ds",
    py_bool store_sequences = True,
):
    """
    Convert a BGZF byte range to optimized Parquet.

    For parallel processing - each worker handles a byte range.

    Args:
        input_file: Input BGZF file (BAM/SAM.gz)
        output_file: Output Parquet file
        start_offset: Starting BGZF block offset
        end_offset: Ending BGZF block offset
        batch_size: Records per batch
        compression_level: ZSTD compression level
        num_threads: HTSlib threads
        calculate_pmd: Calculate PMD scores
        library_type: "ds" or "ss"
        store_sequences: Store sequences

    Returns:
        Dictionary with statistics
    """
    cdef samFile *sam_fp
    cdef BGZF *bgzf_fp
    cdef bam_hdr_t *header
    cdef bam1_t *aln
    cdef int ret
    cdef uint64_t total_records = 0
    cdef double start_time = time.time()
    cdef char seq_buf[100000]
    cdef char cigar_buf[10000]
    cdef char tags_cold_buf[50000]
    cdef int64_t current_pos
    cdef ExtractedRecord rec
    cdef cbool do_pmd = calculate_pmd
    cdef cbool do_store_seq = store_sequences
    cdef string packed_seq

    # PMD config
    cdef AlignmentScoringConfig pmd_config
    pmd_config.calculate_pmd = calculate_pmd
    pmd_config.is_single_stranded = (library_type == "ss")
    pmd_config.minimum_read_identity = 0.0
    pmd_config.minimum_read_length = 0
    pmd_config.maximum_read_length = 100000
    pmd_config.global_min_score = -1e20
    pmd_config.global_max_score = 1e20

    # Always initialize lookup tables for ZS score calculation
    initialize_quality_lookup_tables()

    # Create output directory
    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    # C++ objects
    cdef OptimizedParquetWriter* writer = NULL
    cdef OptimizedAlignmentBatch* batch = new OptimizedAlignmentBatch()
    batch.reserve(batch_size)

    # Open file
    input_bytes = input_file.encode('utf-8')
    sam_fp = sam_open(input_bytes, b"r")
    if sam_fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    if num_threads > 1:
        hts_set_threads(sam_fp, num_threads - 1)

    header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        raise IOError(f"Failed to read header from {input_file}")

    # Get BGZF and seek
    bgzf_fp = hts_get_bgzfp(<htsFile*>sam_fp)
    if bgzf_fp == NULL:
        bam_hdr_destroy(header)
        sam_close(sam_fp)
        raise IOError(f"Failed to get BGZF pointer")

    cdef int64_t virtual_offset = start_offset << 16
    if bgzf_seek(bgzf_fp, virtual_offset, 0) < 0:
        bam_hdr_destroy(header)
        sam_close(sam_fp)
        raise IOError(f"Failed to seek to offset {start_offset}")

    aln = bam_init1()

    # Create writer
    output_bytes = output_file.encode('utf-8')
    writer = new OptimizedParquetWriter(string(output_bytes), compression_level)

    try:
        while True:
            ret = sam_read1(sam_fp, header, aln)
            if ret < 0:
                break

            # Check if past end
            current_pos = bgzf_tell(bgzf_fp) >> 16
            if current_pos >= end_offset:
                break

            if aln.core.tid < 0:
                continue

            total_records += 1

            with nogil:
                extract_record_optimized(
                    aln, header, &rec, seq_buf, cigar_buf, tags_cold_buf,
                    &pmd_config, do_pmd, do_store_seq
                )

            # Append to batch
            batch.read_ids.push_back(total_records)
            batch.read_names.push_back(string(rec.qname))
            batch.ref_ids.push_back(<uint32_t>rec.tid)
            batch.positions.push_back(rec.pos)
            batch.mapqs.push_back(rec.mapq)
            batch.flags.push_back(rec.flag)
            batch.cigars.push_back(string(<const char*>cigar_buf))
            batch.template_lengths.push_back(rec.isize)
            batch.mate_ref_ids.push_back(rec.mtid)
            batch.mate_positions.push_back(rec.mpos)

            if do_store_seq:
                packed_seq = pack_sequence_2bit(<const char*>seq_buf, <size_t>rec.l_qseq)
                batch.sequences_packed.push_back(packed_seq)
                batch.sequence_lengths.push_back(<uint16_t>rec.l_qseq)
                batch.qualities.push_back(string(<char*>bam_get_qual(aln), rec.l_qseq))
            else:
                batch.sequences_packed.push_back(string())
                batch.sequence_lengths.push_back(0)
                batch.qualities.push_back(string())

            batch.AS.push_back(rec.alignment_score)
            batch.NM.push_back(rec.edit_distance)
            batch.XS.push_back(rec.xs_score)
            batch.MD.push_back(string(rec.md_str))
            batch.tags_cold.push_back(string(<const char*>tags_cold_buf))

            # Derived fields (precomputed for fast queries)
            batch.aligned_lengths.push_back(rec.aligned_length)
            batch.is_reverse.push_back(rec.is_reverse)

            # Pipeline results
            batch.ani.push_back(rec.ani)
            batch.zs_score.push_back(rec.zs_score)
            batch.pmd_score.push_back(rec.pmd_score)
            batch.filter_passed.push_back(True)
            batch.lca_taxid.push_back(rec.lca_taxid)
            batch.reassigned_ref_id.push_back(rec.reassigned_ref_id)
            batch.zp_posterior.push_back(rec.zp_posterior)
            batch.read_groups.push_back(string(rec.rg_str))

            if batch.size() >= batch_size:
                writer.WriteBatch(batch[0])
                batch.clear()
                batch.reserve(batch_size)

        if batch.size() > 0:
            writer.WriteBatch(batch[0])

        writer.Close()

        return {
            'total_records': total_records,
            'processing_time_seconds': time.time() - start_time,
        }

    finally:
        del batch
        if writer != NULL:
            del writer
        bam_destroy1(aln)
        bam_hdr_destroy(header)
        sam_close(sam_fp)
