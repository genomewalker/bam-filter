# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

"""
Fast Parquet converter using Arrow C++ API directly.

Architecture:
- HTSlib for SAM/BAM reading (C library)
- Arrow C++ for Parquet writing (no DuckDB overhead)
- Columnar batch processing (100K records per batch)
- Hash partitioning for parallel query performance
"""

from libc.stdlib cimport malloc, free, realloc
from libc.string cimport strdup, strlen, memcpy
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t
from cpython.bytes cimport PyBytes_FromStringAndSize
from cpython cimport bool as py_bool

import os
from pathlib import Path
import time

# HTSlib declarations
cdef extern from "htslib/sam.h":
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

    ctypedef struct hts_itr_t:
        pass

    samFile* sam_open(const char *fn, const char *mode)
    int sam_close(samFile *fp)
    bam_hdr_t* sam_hdr_read(samFile *fp)
    void bam_hdr_destroy(bam_hdr_t *h)
    bam1_t* bam_init1()
    void bam_destroy1(bam1_t *b)
    int sam_read1(samFile *fp, bam_hdr_t *h, bam1_t *b)

    uint8_t* bam_get_seq(bam1_t *b)
    uint8_t* bam_get_qual(bam1_t *b)
    uint32_t* bam_get_cigar(bam1_t *b)
    char* bam_get_qname(bam1_t *b)
    uint8_t* bam_get_aux(bam1_t *b)

    uint8_t bam_seqi(uint8_t *s, int i)

    uint8_t* bam_aux_get(bam1_t *b, const char tag[2])
    int32_t bam_aux2i(uint8_t *s)
    float bam_aux2f(uint8_t *s)
    char* bam_aux2Z(uint8_t *s)

    char* sam_hdr_tid2name(bam_hdr_t *h, int tid)

    int bam_cigar2qlen(int n_cigar, uint32_t *cigar)
    int bam_endpos(bam1_t *b)

# Arrow C++ declarations
cdef extern from "arrow_parquet_writer.hpp" namespace "bam_filter":
    cdef cppclass AlignmentBatch:
        AlignmentBatch() except +
        size_t size()
        void clear()
        void reserve(size_t)

        # Vectors - direct access
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
        vector[bint] filter_passed
        vector[string] read_groups
        vector[string] cigars
        vector[string] sequences
        vector[string] qualities

    cdef cppclass ParquetWriter:
        ParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const AlignmentBatch& batch) except +
        void Close() except +

cdef extern from "<vector>" namespace "std":
    cdef cppclass vector[T]:
        void push_back(T&)
        void clear()
        void reserve(size_t)
        size_t size()
        T& operator[](size_t)

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(const char*)
        string(const char*, size_t)
        const char* c_str()

# C++ batch wrapper
cdef class CppAlignmentBatch:
    """Cython wrapper for C++ AlignmentBatch."""
    cdef AlignmentBatch* batch

    def __cinit__(self):
        self.batch = new AlignmentBatch()

    def __dealloc__(self):
        del self.batch

    cdef size_t size(self):
        return self.batch.size()

    cdef void clear(self):
        self.batch.clear()


cdef char seq_nt16_str[16]
seq_nt16_str = [b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
                b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N']


cdef inline void extract_sequence(bam1_t *aln, char* result):
    """Extract DNA sequence from BAM record."""
    cdef uint8_t *seq = bam_get_seq(aln)
    cdef int i
    cdef int l_seq = aln.core.l_qseq

    for i in range(l_seq):
        result[i] = seq_nt16_str[bam_seqi(seq, i)]
    result[l_seq] = 0


cdef inline void extract_cigar(bam1_t *aln, char* result):
    """Extract CIGAR string from BAM record."""
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
        # Format: <length><op>
        pos += sprintf(&result[pos], "%d", op >> 4)
        op_char = "MIDNSHP=X"[op & 0xf]
        result[pos] = op_char
        pos += 1

    result[pos] = 0


cdef extern from "stdio.h":
    int sprintf(char *str, const char *format, ...)


def convert_sam_bam_to_parquet_cpp(
    str input_file,
    str output_dir,
    int num_partitions = 128,
    int batch_size = 100000,
    py_bool write_by_reference = True,
    py_bool write_by_read = True,
    int compression_level = 3,
):
    """
    Convert SAM/BAM to Parquet using Arrow C++ directly.

    Parameters:
    -----------
    input_file : str
        Path to SAM/BAM/CRAM file
    output_dir : str
        Output directory for Parquet files
    num_partitions : int
        Number of hash partitions (default: 128)
    batch_size : int
        Records to accumulate before writing (default: 100K)
    write_by_reference : bool
        Create by_reference partitions
    write_by_read : bool
        Create by_read partitions
    compression_level : int
        ZSTD compression level (default: 3)

    Returns:
    --------
    dict with statistics
    """
    cdef samFile *sam_fp
    cdef bam_hdr_t *header
    cdef bam1_t *aln
    cdef int ret
    cdef uint64_t total_records = 0
    cdef double start_time = time.time()
    cdef double last_print = start_time
    cdef int ref_partition, read_partition
    cdef uint8_t *aux
    cdef const char *ref_name
    cdef char seq_buf[100000]
    cdef char cigar_buf[10000]
    cdef int32_t edit_dist
    cdef float ani_val

    # Create output directories
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if write_by_reference:
        (output_path / "by_reference").mkdir(exist_ok=True)
    if write_by_read:
        (output_path / "by_read").mkdir(exist_ok=True)

    # Open SAM/BAM file
    input_bytes = input_file.encode('utf-8')
    sam_fp = sam_open(input_bytes, b"r")
    if sam_fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        raise IOError(f"Failed to read header from {input_file}")

    aln = bam_init1()

    # Partition writers (lazy initialization)
    ref_writers = {}
    read_writers = {}

    # Partition buffers
    ref_buffers = {}
    read_buffers = {}

    print(f"Processing {input_file}...")
    print(f"Writing to {output_dir}")
    print(f"Partitions: {num_partitions}, Batch size: {batch_size:,}")

    try:
        while True:
            ret = sam_read1(sam_fp, header, aln)
            if ret < 0:
                break

            # Skip unmapped
            if aln.core.tid < 0:
                continue

            # Get reference name
            ref_name = sam_hdr_tid2name(header, aln.core.tid)
            if ref_name == NULL:
                ref_name = b"*"

            # Extract sequence and CIGAR
            extract_sequence(aln, seq_buf)
            extract_cigar(aln, cigar_buf)

            # Create record dict
            record = {
                'read_id': total_records,
                'read_name': bam_get_qname(aln).decode('utf-8'),
                'ref_id': aln.core.tid,
                'ref_name': ref_name.decode('utf-8'),
                'position': aln.core.pos,
                'end_position': bam_endpos(aln),
                'mapq': aln.core.qual,
                'flag': aln.core.flag,
                'alignment_length': aln.core.l_qseq,
                'template_length': aln.core.isize,
                'mate_ref_id': aln.core.mtid if aln.core.mtid >= 0 else -1,
                'mate_position': aln.core.mpos if aln.core.mpos >= 0 else -1,
            }

            # Extract optional tags
            aux = bam_aux_get(aln, b"AS")
            record['alignment_score'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"XS")
            record['xs_score'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"NM")
            edit_dist = bam_aux2i(aux) if aux != NULL else -1
            record['edit_distance'] = edit_dist

            aux = bam_aux_get(aln, b"XM")
            record['num_mismatches'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"XO")
            record['num_gap_opens'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"XG")
            record['num_gap_extensions'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"MD")
            record['md_string'] = bam_aux2Z(aux).decode('utf-8') if aux != NULL else ""

            # Calculate ANI
            if edit_dist >= 0 and aln.core.l_qseq > 0:
                ani_val = (1.0 - <float>edit_dist / <float>aln.core.l_qseq) * 100.0
                record['ani'] = ani_val
            else:
                record['ani'] = -1.0

            aux = bam_aux_get(aln, b"PM")
            if aux == NULL:
                aux = bam_aux_get(aln, b"PMD")
            record['pmd_score'] = bam_aux2f(aux) if aux != NULL else -1.0

            aux = bam_aux_get(aln, b"ZS")
            record['zs_score'] = bam_aux2f(aux) if aux != NULL else -1.0

            aux = bam_aux_get(aln, b"ZP")
            record['zp_posterior'] = bam_aux2f(aux) if aux != NULL else -1.0

            aux = bam_aux_get(aln, b"ZT")
            record['lca_taxid'] = bam_aux2i(aux) if aux != NULL else -1

            aux = bam_aux_get(aln, b"ZR")
            record['reassigned_ref_id'] = bam_aux2i(aux) if aux != NULL else -1

            record['filter_passed'] = True

            aux = bam_aux_get(aln, b"RG")
            record['read_group'] = bam_aux2Z(aux).decode('utf-8') if aux != NULL else ""

            record['cigar'] = cigar_buf.decode('utf-8')
            record['sequence'] = seq_buf[:aln.core.l_qseq].decode('utf-8')

            # Quality scores
            qual_bytes = bytes(bam_get_qual(aln)[:aln.core.l_qseq])
            record['quality'] = qual_bytes

            # Hash partition
            ref_partition = aln.core.tid % num_partitions
            read_partition = total_records % num_partitions

            # Add to buffers
            if write_by_reference:
                if ref_partition not in ref_buffers:
                    ref_buffers[ref_partition] = PyAlignmentBatch()

                _append_to_batch(ref_buffers[ref_partition], record)

                # Flush if full
                if len(ref_buffers[ref_partition]) >= batch_size:
                    _flush_partition_python(
                        ref_buffers[ref_partition],
                        ref_partition,
                        ref_writers,
                        str(output_path / "by_reference"),
                        compression_level,
                        True
                    )
                    ref_buffers[ref_partition].clear()

            if write_by_read:
                if read_partition not in read_buffers:
                    read_buffers[read_partition] = PyAlignmentBatch()

                _append_to_batch(read_buffers[read_partition], record)

                if len(read_buffers[read_partition]) >= batch_size:
                    _flush_partition_python(
                        read_buffers[read_partition],
                        read_partition,
                        read_writers,
                        str(output_path / "by_read"),
                        compression_level,
                        False
                    )
                    read_buffers[read_partition].clear()

            total_records += 1

            # Progress
            if time.time() - last_print >= 10:
                elapsed = time.time() - start_time
                print(f"\rProgress: {total_records:,} records ({total_records/elapsed/1e6:.2f} M/sec, {elapsed:.0f}s)",
                      end='', flush=True)
                last_print = time.time()

        # Flush remaining buffers
        print(f"\n\nFlushing remaining buffers...")
        for partition, buffer in ref_buffers.items():
            if len(buffer) > 0:
                _flush_partition_python(
                    buffer,
                    partition,
                    ref_writers,
                    str(output_path / "by_reference"),
                    compression_level,
                    True
                )

        for partition, buffer in read_buffers.items():
            if len(buffer) > 0:
                _flush_partition_python(
                    buffer,
                    partition,
                    read_writers,
                    str(output_path / "by_read"),
                    compression_level,
                    False
                )

        # Close all writers
        for writer in ref_writers.values():
            writer.close()
        for writer in read_writers.values():
            writer.close()

        duration = time.time() - start_time

        return {
            'total_records': total_records,
            'processing_time_seconds': duration,
        }

    finally:
        bam_destroy1(aln)
        bam_hdr_destroy(header)
        sam_close(sam_fp)


cdef inline void _append_to_batch_cpp(
    AlignmentBatch* batch,
    uint64_t read_id,
    const char* read_name,
    uint32_t ref_id,
    const char* ref_name,
    int32_t position,
    int32_t end_position,
    uint8_t mapq,
    uint16_t flag,
    uint16_t alignment_length,
    int32_t template_length,
    int32_t mate_ref_id,
    int32_t mate_position,
    int32_t alignment_score,
    int32_t xs_score,
    int32_t edit_distance,
    int32_t num_mismatches,
    int32_t num_gap_opens,
    int32_t num_gap_extensions,
    const char* md_string,
    float ani,
    float pmd_score,
    float zs_score,
    float zp_posterior,
    int32_t lca_taxid,
    int32_t reassigned_ref_id,
    bint filter_passed,
    const char* read_group,
    const char* cigar,
    const char* sequence,
    int seq_len,
    const char* quality,
    int qual_len
) nogil:
    """Append a record to the C++ batch - pure C, no Python."""
    batch.read_ids.push_back(read_id)
    batch.read_names.push_back(string(read_name))
    batch.ref_ids.push_back(ref_id)
    batch.ref_names.push_back(string(ref_name))
    batch.positions.push_back(position)
    batch.end_positions.push_back(end_position)
    batch.mapqs.push_back(mapq)
    batch.flags.push_back(flag)
    batch.alignment_lengths.push_back(alignment_length)
    batch.template_lengths.push_back(template_length)
    batch.mate_ref_ids.push_back(mate_ref_id)
    batch.mate_positions.push_back(mate_position)
    batch.alignment_scores.push_back(alignment_score)
    batch.xs_scores.push_back(xs_score)
    batch.edit_distances.push_back(edit_distance)
    batch.num_mismatches.push_back(num_mismatches)
    batch.num_gap_opens.push_back(num_gap_opens)
    batch.num_gap_extensions.push_back(num_gap_extensions)
    batch.md_strings.push_back(string(md_string))
    batch.anis.push_back(ani)
    batch.pmd_scores.push_back(pmd_score)
    batch.zs_scores.push_back(zs_score)
    batch.zp_posteriors.push_back(zp_posterior)
    batch.lca_taxids.push_back(lca_taxid)
    batch.reassigned_ref_ids.push_back(reassigned_ref_id)
    batch.filter_passed.push_back(filter_passed)
    batch.read_groups.push_back(string(read_group))
    batch.cigars.push_back(string(cigar))
    batch.sequences.push_back(string(sequence, seq_len))
    batch.qualities.push_back(string(quality, qual_len))


cdef void _flush_partition_python(
    PyAlignmentBatch batch,
    int partition_id,
    dict writers,
    str output_dir,
    int compression_level,
    py_bool is_reference
):
    """Flush batch using PyArrow (temporary - will switch to C++ later)."""
    import pyarrow as pa
    import pyarrow.parquet as pq

    if len(batch) == 0:
        return

    # Define schema
    schema = pa.schema([
        ('read_id', pa.uint64()),
        ('read_name', pa.string()),
        ('ref_id', pa.uint32()),
        ('ref_name', pa.string()),
        ('position', pa.int32()),
        ('end_position', pa.int32()),
        ('mapq', pa.uint8()),
        ('flag', pa.uint16()),
        ('alignment_length', pa.uint16()),
        ('template_length', pa.int32()),
        ('mate_ref_id', pa.int32()),
        ('mate_position', pa.int32()),
        ('alignment_score', pa.int32()),
        ('xs_score', pa.int32()),
        ('edit_distance', pa.int32()),
        ('num_mismatches', pa.int32()),
        ('num_gap_opens', pa.int32()),
        ('num_gap_extensions', pa.int32()),
        ('md_string', pa.string()),
        ('ani', pa.float32()),
        ('pmd_score', pa.float32()),
        ('zs_score', pa.float32()),
        ('zp_posterior', pa.float32()),
        ('lca_taxid', pa.int32()),
        ('reassigned_ref_id', pa.int32()),
        ('filter_passed', pa.bool_()),
        ('read_group', pa.string()),
        ('cigar', pa.string()),
        ('sequence', pa.binary()),
        ('quality', pa.binary()),
    ])

    # Create Arrow arrays
    arrays = [
        pa.array(batch.read_ids, type=pa.uint64()),
        pa.array(batch.read_names, type=pa.string()),
        pa.array(batch.ref_ids, type=pa.uint32()),
        pa.array(batch.ref_names, type=pa.string()),
        pa.array(batch.positions, type=pa.int32()),
        pa.array(batch.end_positions, type=pa.int32()),
        pa.array(batch.mapqs, type=pa.uint8()),
        pa.array(batch.flags, type=pa.uint16()),
        pa.array(batch.alignment_lengths, type=pa.uint16()),
        pa.array(batch.template_lengths, type=pa.int32()),
        pa.array(batch.mate_ref_ids, type=pa.int32()),
        pa.array(batch.mate_positions, type=pa.int32()),
        pa.array(batch.alignment_scores, type=pa.int32()),
        pa.array(batch.xs_scores, type=pa.int32()),
        pa.array(batch.edit_distances, type=pa.int32()),
        pa.array(batch.num_mismatches, type=pa.int32()),
        pa.array(batch.num_gap_opens, type=pa.int32()),
        pa.array(batch.num_gap_extensions, type=pa.int32()),
        pa.array(batch.md_strings, type=pa.string()),
        pa.array(batch.anis, type=pa.float32()),
        pa.array(batch.pmd_scores, type=pa.float32()),
        pa.array(batch.zs_scores, type=pa.float32()),
        pa.array(batch.zp_posteriors, type=pa.float32()),
        pa.array(batch.lca_taxids, type=pa.int32()),
        pa.array(batch.reassigned_ref_ids, type=pa.int32()),
        pa.array(batch.filter_passed, type=pa.bool_()),
        pa.array(batch.read_groups, type=pa.string()),
        pa.array(batch.cigars, type=pa.string()),
        pa.array(batch.sequences, type=pa.binary()),
        pa.array(batch.qualities, type=pa.binary()),
    ]

    # Create table
    table = pa.Table.from_arrays(arrays, schema=schema)

    # Get or create writer
    if partition_id not in writers:
        prefix = "by_reference" if is_reference else "by_read"
        filename = Path(output_dir) / f"{prefix}_p{partition_id:04d}.parquet"
        writers[partition_id] = pq.ParquetWriter(
            str(filename),
            schema,
            compression='zstd',
            compression_level=compression_level,
            use_dictionary=True,
        )

    # Write table
    writers[partition_id].write_table(table)
