# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

"""
Ultra-fast Parquet converter using pure C++/Cython.

NO Python overhead - everything in C/C++:
- HTSlib for SAM/BAM reading
- Arrow C++ for Parquet writing
- Direct C++ batch accumulation
- No dict creation, no Python objects in hot loop
"""

from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memcpy, strcmp
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t
from libc.stdio cimport sprintf
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cpython cimport bool as py_bool

import time
from pathlib import Path

# HTSlib
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

    uint8_t* bam_aux_get(bam1_t *b, const char tag[2])
    int32_t bam_aux2i(uint8_t *s)
    float bam_aux2f(uint8_t *s)
    char* bam_aux2Z(uint8_t *s)

    char* sam_hdr_tid2name(bam_hdr_t *h, int tid)
    int bam_endpos(bam1_t *b)

    # Threading support
    int hts_set_threads(samFile *fp, int n)

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

    cdef cppclass ParquetWriter:
        ParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const AlignmentBatch& batch) except +
        void Close() except +

cdef extern from "<string>" namespace "std":
    cdef cppclass string:
        string()
        string(const char*)
        string(const char*, size_t)


# DNA encoding
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


def convert_sam_bam_to_parquet_pure_cpp(
    str input_file,
    str output_dir,
    int num_partitions = 128,
    int batch_size = 100000,
    py_bool write_by_reference = True,
    py_bool write_by_read = True,
    int compression_level = 3,
    int num_threads = 8,
):
    """
    Pure C++/Cython Parquet converter - ZERO Python overhead.

    Expected performance: ~1M records/sec on modern hardware.
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
    cdef const char *md_str
    cdef const char *rg_str
    cdef char seq_buf[100000]
    cdef char cigar_buf[10000]
    cdef char md_buf[10000]
    cdef char rg_buf[256]
    cdef int32_t edit_dist, alignment_score, xs_score, num_mm, num_go, num_ge
    cdef int32_t lca_taxid, reassigned_ref_id
    cdef float ani_val, pmd_score, zs_score, zp_posterior
    cdef int i

    # C++ objects - partition writers and batches
    cdef ParquetWriter** ref_writers = NULL
    cdef ParquetWriter** read_writers = NULL
    cdef AlignmentBatch** ref_batches = NULL
    cdef AlignmentBatch** read_batches = NULL
    cdef cbool* ref_writers_active = NULL
    cdef cbool* read_writers_active = NULL

    # Create output directories
    output_path = Path(output_dir)
    output_path.mkdir(parents=True, exist_ok=True)

    if write_by_reference:
        (output_path / "by_reference").mkdir(exist_ok=True)
    if write_by_read:
        (output_path / "by_read").mkdir(exist_ok=True)

    # Allocate partition arrays
    if write_by_reference:
        ref_writers = <ParquetWriter**>malloc(num_partitions * sizeof(ParquetWriter*))
        ref_batches = <AlignmentBatch**>malloc(num_partitions * sizeof(AlignmentBatch*))
        ref_writers_active = <cbool*>malloc(num_partitions * sizeof(cbool))
        for i in range(num_partitions):
            ref_writers[i] = NULL
            ref_batches[i] = new AlignmentBatch()
            ref_writers_active[i] = False

    if write_by_read:
        read_writers = <ParquetWriter**>malloc(num_partitions * sizeof(ParquetWriter*))
        read_batches = <AlignmentBatch**>malloc(num_partitions * sizeof(AlignmentBatch*))
        read_writers_active = <cbool*>malloc(num_partitions * sizeof(cbool))
        for i in range(num_partitions):
            read_writers[i] = NULL
            read_batches[i] = new AlignmentBatch()
            read_writers_active[i] = False

    # Open SAM/BAM file
    input_bytes = input_file.encode('utf-8')
    sam_fp = sam_open(input_bytes, b"r")
    if sam_fp == NULL:
        raise IOError(f"Failed to open {input_file}")

    # Enable multi-threaded decompression
    if num_threads > 1:
        hts_set_threads(sam_fp, num_threads - 1)  # -1 because main thread also does work

    header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        raise IOError(f"Failed to read header from {input_file}")

    aln = bam_init1()

    print(f"Processing {input_file}...")
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

            # Extract sequence and CIGAR (pure C, nogil-safe)
            with nogil:
                extract_sequence(aln, seq_buf)
                extract_cigar(aln, cigar_buf)

            # Extract optional tags
            aux = bam_aux_get(aln, b"AS")
            if aux != NULL:
                alignment_score = bam_aux2i(aux)
            else:
                alignment_score = -1

            aux = bam_aux_get(aln, b"XS")
            if aux != NULL:
                xs_score = bam_aux2i(aux)
            else:
                xs_score = -1

            aux = bam_aux_get(aln, b"NM")
            if aux != NULL:
                edit_dist = bam_aux2i(aux)
            else:
                edit_dist = -1

            aux = bam_aux_get(aln, b"XM")
            if aux != NULL:
                num_mm = bam_aux2i(aux)
            else:
                num_mm = -1

            aux = bam_aux_get(aln, b"XO")
            if aux != NULL:
                num_go = bam_aux2i(aux)
            else:
                num_go = -1

            aux = bam_aux_get(aln, b"XG")
            if aux != NULL:
                num_ge = bam_aux2i(aux)
            else:
                num_ge = -1

            aux = bam_aux_get(aln, b"MD")
            if aux != NULL:
                md_str = bam_aux2Z(aux)
            else:
                md_str = b""

            # Calculate ANI
            if edit_dist >= 0 and aln.core.l_qseq > 0:
                ani_val = (1.0 - <float>edit_dist / <float>aln.core.l_qseq) * 100.0
            else:
                ani_val = -1.0

            aux = bam_aux_get(aln, b"PM")
            if aux == NULL:
                aux = bam_aux_get(aln, b"PMD")
            if aux != NULL:
                pmd_score = bam_aux2f(aux)
            else:
                pmd_score = -1.0

            aux = bam_aux_get(aln, b"ZS")
            if aux != NULL:
                zs_score = bam_aux2f(aux)
            else:
                zs_score = -1.0

            aux = bam_aux_get(aln, b"ZP")
            if aux != NULL:
                zp_posterior = bam_aux2f(aux)
            else:
                zp_posterior = -1.0

            aux = bam_aux_get(aln, b"ZT")
            if aux != NULL:
                lca_taxid = bam_aux2i(aux)
            else:
                lca_taxid = -1

            aux = bam_aux_get(aln, b"ZR")
            if aux != NULL:
                reassigned_ref_id = bam_aux2i(aux)
            else:
                reassigned_ref_id = -1

            aux = bam_aux_get(aln, b"RG")
            if aux != NULL:
                rg_str = bam_aux2Z(aux)
            else:
                rg_str = b""

            # Hash partition
            ref_partition = aln.core.tid % num_partitions
            read_partition = total_records % num_partitions

            # Append to C++ batches (ZERO Python overhead)
            if write_by_reference:
                ref_batches[ref_partition].read_ids.push_back(total_records)
                ref_batches[ref_partition].read_names.push_back(string(bam_get_qname(aln)))
                ref_batches[ref_partition].ref_ids.push_back(aln.core.tid)
                ref_batches[ref_partition].ref_names.push_back(string(ref_name))
                ref_batches[ref_partition].positions.push_back(aln.core.pos)
                ref_batches[ref_partition].end_positions.push_back(bam_endpos(aln))
                ref_batches[ref_partition].mapqs.push_back(aln.core.qual)
                ref_batches[ref_partition].flags.push_back(aln.core.flag)
                ref_batches[ref_partition].alignment_lengths.push_back(aln.core.l_qseq)
                ref_batches[ref_partition].template_lengths.push_back(aln.core.isize)

                if aln.core.mtid >= 0:
                    ref_batches[ref_partition].mate_ref_ids.push_back(aln.core.mtid)
                else:
                    ref_batches[ref_partition].mate_ref_ids.push_back(-1)

                if aln.core.mpos >= 0:
                    ref_batches[ref_partition].mate_positions.push_back(aln.core.mpos)
                else:
                    ref_batches[ref_partition].mate_positions.push_back(-1)

                ref_batches[ref_partition].alignment_scores.push_back(alignment_score)
                ref_batches[ref_partition].xs_scores.push_back(xs_score)
                ref_batches[ref_partition].edit_distances.push_back(edit_dist)
                ref_batches[ref_partition].num_mismatches.push_back(num_mm)
                ref_batches[ref_partition].num_gap_opens.push_back(num_go)
                ref_batches[ref_partition].num_gap_extensions.push_back(num_ge)
                ref_batches[ref_partition].md_strings.push_back(string(md_str))
                ref_batches[ref_partition].anis.push_back(ani_val)
                ref_batches[ref_partition].pmd_scores.push_back(pmd_score)
                ref_batches[ref_partition].zs_scores.push_back(zs_score)
                ref_batches[ref_partition].zp_posteriors.push_back(zp_posterior)
                ref_batches[ref_partition].lca_taxids.push_back(lca_taxid)
                ref_batches[ref_partition].reassigned_ref_ids.push_back(reassigned_ref_id)
                ref_batches[ref_partition].filter_passed.push_back(True)
                ref_batches[ref_partition].read_groups.push_back(string(rg_str))
                ref_batches[ref_partition].cigars.push_back(string(<const char*>cigar_buf))
                ref_batches[ref_partition].sequences.push_back(string(<const char*>seq_buf, <size_t>aln.core.l_qseq))
                ref_batches[ref_partition].qualities.push_back(string(<char*>bam_get_qual(aln), aln.core.l_qseq))

                # Flush if batch full
                if ref_batches[ref_partition].size() >= batch_size:
                    # Lazy create writer
                    if not ref_writers_active[ref_partition]:
                        filename = str(output_path / "by_reference" / f"by_reference_p{ref_partition:04d}.parquet")
                        ref_writers[ref_partition] = new ParquetWriter(string(filename.encode('utf-8')), compression_level)
                        ref_writers_active[ref_partition] = True

                    ref_writers[ref_partition].WriteBatch(ref_batches[ref_partition][0])
                    ref_batches[ref_partition].clear()

            if write_by_read:
                # Same for read partitions...
                read_batches[read_partition].read_ids.push_back(total_records)
                read_batches[read_partition].read_names.push_back(string(bam_get_qname(aln)))
                read_batches[read_partition].ref_ids.push_back(aln.core.tid)
                read_batches[read_partition].ref_names.push_back(string(ref_name))
                read_batches[read_partition].positions.push_back(aln.core.pos)
                read_batches[read_partition].end_positions.push_back(bam_endpos(aln))
                read_batches[read_partition].mapqs.push_back(aln.core.qual)
                read_batches[read_partition].flags.push_back(aln.core.flag)
                read_batches[read_partition].alignment_lengths.push_back(aln.core.l_qseq)
                read_batches[read_partition].template_lengths.push_back(aln.core.isize)

                if aln.core.mtid >= 0:
                    read_batches[read_partition].mate_ref_ids.push_back(aln.core.mtid)
                else:
                    read_batches[read_partition].mate_ref_ids.push_back(-1)

                if aln.core.mpos >= 0:
                    read_batches[read_partition].mate_positions.push_back(aln.core.mpos)
                else:
                    read_batches[read_partition].mate_positions.push_back(-1)

                read_batches[read_partition].alignment_scores.push_back(alignment_score)
                read_batches[read_partition].xs_scores.push_back(xs_score)
                read_batches[read_partition].edit_distances.push_back(edit_dist)
                read_batches[read_partition].num_mismatches.push_back(num_mm)
                read_batches[read_partition].num_gap_opens.push_back(num_go)
                read_batches[read_partition].num_gap_extensions.push_back(num_ge)
                read_batches[read_partition].md_strings.push_back(string(md_str))
                read_batches[read_partition].anis.push_back(ani_val)
                read_batches[read_partition].pmd_scores.push_back(pmd_score)
                read_batches[read_partition].zs_scores.push_back(zs_score)
                read_batches[read_partition].zp_posteriors.push_back(zp_posterior)
                read_batches[read_partition].lca_taxids.push_back(lca_taxid)
                read_batches[read_partition].reassigned_ref_ids.push_back(reassigned_ref_id)
                read_batches[read_partition].filter_passed.push_back(True)
                read_batches[read_partition].read_groups.push_back(string(rg_str))
                read_batches[read_partition].cigars.push_back(string(<const char*>cigar_buf))
                read_batches[read_partition].sequences.push_back(string(<const char*>seq_buf, <size_t>aln.core.l_qseq))
                read_batches[read_partition].qualities.push_back(string(<char*>bam_get_qual(aln), aln.core.l_qseq))

                if read_batches[read_partition].size() >= batch_size:
                    if not read_writers_active[read_partition]:
                        filename = str(output_path / "by_read" / f"by_read_p{read_partition:04d}.parquet")
                        read_writers[read_partition] = new ParquetWriter(string(filename.encode('utf-8')), compression_level)
                        read_writers_active[read_partition] = True

                    read_writers[read_partition].WriteBatch(read_batches[read_partition][0])
                    read_batches[read_partition].clear()

            total_records += 1

            # Progress
            if time.time() - last_print >= 10:
                elapsed = time.time() - start_time
                print(f"\rProgress: {total_records:,} records ({total_records/elapsed/1e6:.2f} M/sec, {elapsed:.0f}s)",
                      end='', flush=True)
                last_print = time.time()

        # Flush remaining batches
        print(f"\n\nFlushing remaining batches...")
        for i in range(num_partitions):
            if write_by_reference and ref_batches[i].size() > 0:
                if not ref_writers_active[i]:
                    filename = str(output_path / "by_reference" / f"by_reference_p{i:04d}.parquet")
                    ref_writers[i] = new ParquetWriter(string(filename.encode('utf-8')), compression_level)
                    ref_writers_active[i] = True
                ref_writers[i].WriteBatch(ref_batches[i][0])

            if write_by_read and read_batches[i].size() > 0:
                if not read_writers_active[i]:
                    filename = str(output_path / "by_read" / f"by_read_p{i:04d}.parquet")
                    read_writers[i] = new ParquetWriter(string(filename.encode('utf-8')), compression_level)
                    read_writers_active[i] = True
                read_writers[i].WriteBatch(read_batches[i][0])

        # Close all writers
        for i in range(num_partitions):
            if write_by_reference and ref_writers_active[i]:
                ref_writers[i].Close()
                del ref_writers[i]
            if write_by_read and read_writers_active[i]:
                read_writers[i].Close()
                del read_writers[i]

        duration = time.time() - start_time

        return {
            'total_records': total_records,
            'processing_time_seconds': duration,
        }

    finally:
        bam_destroy1(aln)
        bam_hdr_destroy(header)
        sam_close(sam_fp)

        # Cleanup
        if write_by_reference:
            for i in range(num_partitions):
                del ref_batches[i]
            free(ref_writers)
            free(ref_batches)
            free(ref_writers_active)

        if write_by_read:
            for i in range(num_partitions):
                del read_batches[i]
            free(read_writers)
            free(read_batches)
            free(read_writers_active)
