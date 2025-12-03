# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False

"""
Pure C/OpenMP parallel Parquet converter - TRUE threading with shared memory.

NO Python multiprocessing, NO memory duplication.
All threading in C with OpenMP.
"""

from libc.stdlib cimport malloc, free
from libc.string cimport strdup, strlen, memcpy, strcmp
from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int32_t, int64_t
from libc.stdio cimport sprintf, printf
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cbool
from cpython cimport bool as py_bool
from cython.parallel cimport prange
from libc.time cimport time, time_t

import os
from pathlib import Path

# HTSlib declarations
cdef extern from "htslib/sam.h" nogil:
    ctypedef struct bam_hdr_t:
        int32_t n_targets
        char **target_name

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
        int m_data

    ctypedef struct samFile:
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
    uint8_t bam_seqi(uint8_t *s, int i)

    uint8_t* bam_aux_get(bam1_t *b, const char tag[2])
    int32_t bam_aux2i(uint8_t *s)
    float bam_aux2f(uint8_t *s)
    char* bam_aux2Z(uint8_t *s)

    char* sam_hdr_tid2name(bam_hdr_t *h, int tid)
    int bam_endpos(bam1_t *b)
    int hts_set_threads(samFile *fp, int n)

cdef extern from "htslib/sam.h":
    char* BAM_CIGAR_STR

# Arrow C++ writer
cdef extern from "arrow_parquet_writer.hpp" namespace "bam_filter":
    cdef cppclass AlignmentBatch:
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

        size_t size() nogil
        void clear() nogil
        void shrink_to_fit() nogil
        void reserve(size_t capacity) nogil

    cdef cppclass ParquetWriter:
        ParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const AlignmentBatch& batch) except + nogil
        void Close() except + nogil


# Sequence lookup table
cdef char seq_nt16_str[16]
seq_nt16_str = [b'=', b'A', b'C', b'M', b'G', b'R', b'S', b'V',
                b'T', b'W', b'Y', b'H', b'K', b'D', b'B', b'N']


cdef inline void extract_sequence(bam1_t *aln, char* result) nogil:
    cdef uint8_t *seq = bam_get_seq(aln)
    cdef int l_qseq = aln.core.l_qseq
    cdef int i
    for i in range(l_qseq):
        result[i] = seq_nt16_str[bam_seqi(seq, i)]
    result[l_qseq] = 0


cdef inline void extract_cigar(bam1_t *aln, char* result) nogil:
    cdef uint32_t *cigar = bam_get_cigar(aln)
    cdef int n_cigar = aln.core.n_cigar
    cdef int pos = 0
    cdef int i, op, op_len

    for i in range(n_cigar):
        op = cigar[i] & 0xf
        op_len = cigar[i] >> 4
        pos += sprintf(&result[pos], "%d%c", op_len, BAM_CIGAR_STR[op])
    result[pos] = 0


def convert_sam_bam_to_parquet_openmp(
    str input_file,
    str output_dir,
    int num_threads = 32,
    int num_partitions = 16,
    int batch_size = 100000,
    py_bool write_by_reference = True,
    py_bool write_by_read = True,
    int compression_level = 3,
):
    """
    Pure C/OpenMP parallel Parquet converter.

    All 30 fields, true C threading with shared memory.
    NO Python multiprocessing, NO memory duplication.

    Args:
        input_file: SAM/BAM/SAM.gz file
        output_dir: Output directory
        num_threads: OpenMP threads (default: 32)
        num_partitions: Partitions per thread
        batch_size: Records per batch
        write_by_reference: By-reference partitions
        write_by_read: By-read partitions
        compression_level: ZSTD compression

    Returns:
        Statistics dict
    """
    # Estimate records from file size
    cdef uint64_t file_size = os.path.getsize(input_file)
    cdef uint64_t total_records_est = file_size // 700
    cdef uint64_t records_per_thread = total_records_est // num_threads

    print(f"Pure C/OpenMP Parallel Parquet Converter")
    print(f"Input: {input_file}")
    print(f"Output: {output_dir}")
    print(f"OpenMP threads: {num_threads} (shared memory)")
    print(f"File size: {file_size/1e9:.2f} GB")
    print(f"Estimated records: ~{total_records_est:,}")
    print(f"Records/thread: ~{records_per_thread:,}")
    print()

    Path(output_dir).mkdir(parents=True, exist_ok=True)

    # Pre-create all chunk directories (Python, with GIL)
    for thread_idx in range(num_threads):
        chunk_dir = Path(output_dir) / f"chunk_{thread_idx:03d}"
        if write_by_reference:
            (chunk_dir / "by_reference").mkdir(parents=True, exist_ok=True)
        if write_by_read:
            (chunk_dir / "by_read").mkdir(parents=True, exist_ok=True)

    cdef int i

    cdef bytes input_bytes = input_file.encode('utf-8')
    cdef bytes output_bytes = output_dir.encode('utf-8')
    cdef const char* input_cstr = <const char*>input_bytes
    cdef const char* output_cstr = <const char*>output_bytes
    cdef cbool c_write_by_ref = write_by_reference
    cdef cbool c_write_by_read = write_by_read

    # Thread results (shared memory)
    cdef uint64_t *thread_records = <uint64_t*>malloc(num_threads * sizeof(uint64_t))

    cdef time_t start_time = time(NULL)

    # OpenMP parallel region - each thread processes its record range
    printf("Starting %d OpenMP threads...\n", num_threads)

    for i in prange(num_threads, nogil=True, schedule='static', num_threads=num_threads):
        thread_records[i] = process_thread_chunk(
            input_cstr, output_cstr, i,
            i * records_per_thread,
            records_per_thread if i < num_threads - 1 else (total_records_est - i * records_per_thread),
            num_partitions, batch_size, compression_level,
            c_write_by_ref, c_write_by_read
        )

    # Collect results
    cdef uint64_t total_processed = 0
    for i in range(num_threads):
        total_processed += thread_records[i]

    cdef double elapsed = <double>(time(NULL) - start_time)

    free(thread_records)

    print()
    print("="*80)
    print("OPENMP CONVERSION COMPLETE")
    print("="*80)
    print(f"Total records: {total_processed:,}")
    print(f"Total time: {elapsed:.1f}s ({elapsed/60:.2f} min)")
    print(f"Throughput: {total_processed/elapsed:.0f} rec/sec")
    print(f"Throughput: {total_processed/elapsed/1e6:.3f} M rec/sec")
    print()

    return {
        'total_records': total_processed,
        'total_time_seconds': elapsed,
        'throughput_records_per_sec': total_processed / elapsed if elapsed > 0 else 0,
    }


cdef uint64_t process_thread_chunk(
    const char* input_file,
    const char* output_dir,
    int thread_id,
    uint64_t skip_records,
    uint64_t max_records,
    int num_partitions,
    int batch_size,
    int compression_level,
    cbool write_by_reference,
    cbool write_by_read,
) nogil:
    """Process a chunk of records in a thread (pure C, nogil)."""

    cdef samFile *sam_fp = sam_open(input_file, b"r")
    if sam_fp == NULL:
        return 0

    cdef bam_hdr_t *header = sam_hdr_read(sam_fp)
    if header == NULL:
        sam_close(sam_fp)
        return 0

    cdef bam1_t *aln = bam_init1()
    cdef uint64_t records_skipped = 0
    cdef uint64_t records_processed = 0
    cdef uint64_t total_records = skip_records
    cdef int ret
    cdef int i

    # Allocate batches (thread-local) - NO partitioning, each thread writes to single file
    cdef AlignmentBatch ref_batch
    cdef AlignmentBatch read_batch
    cdef ParquetWriter* ref_writer = NULL
    cdef ParquetWriter* read_writer = NULL
    cdef cbool ref_writer_active = False
    cdef cbool read_writer_active = False

    # Buffers for sequence and CIGAR
    cdef char seq_buf[100000]
    cdef char cigar_buf[10000]

    # Variables for field extraction
    cdef uint8_t *aux
    cdef const char *ref_name
    cdef const char *md_str
    cdef const char *rg_str
    cdef int32_t alignment_score, xs_score, edit_dist, num_mm, num_go, num_ge
    cdef int32_t lca_taxid, reassigned_ref_id
    cdef float ani_val, pmd_score, zs_score, zp_posterior
    cdef int ref_partition, read_partition
    cdef char filename_buf[1024]
    cdef char empty_str[1]
    cdef char default_ref[2]
    empty_str[0] = 0
    default_ref[0] = <char>ord('*')
    default_ref[1] = 0

    # Skip to assigned range
    while records_skipped < skip_records:
        ret = sam_read1(sam_fp, header, aln)
        if ret < 0:
            break
        records_skipped += 1

    # Process assigned records with full 30-field extraction
    while records_processed < max_records:
        ret = sam_read1(sam_fp, header, aln)
        if ret < 0:
            break

        # Skip unmapped
        if aln.core.tid < 0:
            continue

        records_processed += 1
        total_records += 1

        # Get reference name
        ref_name = sam_hdr_tid2name(header, aln.core.tid)
        if ref_name == NULL:
            ref_name = default_ref

        # Extract sequence and CIGAR
        extract_sequence(aln, seq_buf)
        extract_cigar(aln, cigar_buf)

        # Extract optional tags
        aux = bam_aux_get(aln, b"AS")
        alignment_score = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"XS")
        xs_score = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"NM")
        edit_dist = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"XM")
        num_mm = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"XO")
        num_go = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"XG")
        num_ge = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"MD")
        md_str = bam_aux2Z(aux) if aux != NULL else empty_str

        # Calculate ANI
        if edit_dist >= 0 and aln.core.l_qseq > 0:
            ani_val = (1.0 - <float>edit_dist / <float>aln.core.l_qseq) * 100.0
        else:
            ani_val = -1.0

        aux = bam_aux_get(aln, b"PM")
        if aux == NULL:
            aux = bam_aux_get(aln, b"PMD")
        pmd_score = bam_aux2f(aux) if aux != NULL else -1.0

        aux = bam_aux_get(aln, b"ZS")
        zs_score = bam_aux2f(aux) if aux != NULL else -1.0

        aux = bam_aux_get(aln, b"ZP")
        zp_posterior = bam_aux2f(aux) if aux != NULL else -1.0

        aux = bam_aux_get(aln, b"ZT")
        lca_taxid = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"ZR")
        reassigned_ref_id = bam_aux2i(aux) if aux != NULL else -1

        aux = bam_aux_get(aln, b"RG")
        rg_str = bam_aux2Z(aux) if aux != NULL else empty_str

        # Append to batches - NO partitioning, thread IS the partition
        if write_by_reference:
            ref_batch.read_ids.push_back(total_records)
            ref_batch.read_names.push_back(string(bam_get_qname(aln)))
            ref_batch.ref_ids.push_back(aln.core.tid)
            ref_batch.ref_names.push_back(string(ref_name))
            ref_batch.positions.push_back(aln.core.pos)
            ref_batch.end_positions.push_back(bam_endpos(aln))
            ref_batch.mapqs.push_back(aln.core.qual)
            ref_batch.flags.push_back(aln.core.flag)
            ref_batch.alignment_lengths.push_back(aln.core.l_qseq)
            ref_batch.template_lengths.push_back(aln.core.isize)
            ref_batch.mate_ref_ids.push_back(aln.core.mtid if aln.core.mtid >= 0 else -1)
            ref_batch.mate_positions.push_back(aln.core.mpos if aln.core.mpos >= 0 else -1)
            ref_batch.alignment_scores.push_back(alignment_score)
            ref_batch.xs_scores.push_back(xs_score)
            ref_batch.edit_distances.push_back(edit_dist)
            ref_batch.num_mismatches.push_back(num_mm)
            ref_batch.num_gap_opens.push_back(num_go)
            ref_batch.num_gap_extensions.push_back(num_ge)
            ref_batch.md_strings.push_back(string(md_str))
            ref_batch.anis.push_back(ani_val)
            ref_batch.pmd_scores.push_back(pmd_score)
            ref_batch.zs_scores.push_back(zs_score)
            ref_batch.zp_posteriors.push_back(zp_posterior)
            ref_batch.lca_taxids.push_back(lca_taxid)
            ref_batch.reassigned_ref_ids.push_back(reassigned_ref_id)
            ref_batch.filter_passed.push_back(True)
            ref_batch.read_groups.push_back(string(rg_str))
            ref_batch.cigars.push_back(string(<const char*>cigar_buf))
            ref_batch.sequences.push_back(string(<const char*>seq_buf, <size_t>aln.core.l_qseq))
            ref_batch.qualities.push_back(string(<char*>bam_get_qual(aln), aln.core.l_qseq))

            # Flush if batch full
            if ref_batch.size() >= <size_t>batch_size:
                # Lazy create writer
                if not ref_writer_active:
                    sprintf(filename_buf, "%s/chunk_%03d/by_reference.parquet",
                            output_dir, thread_id)
                    with gil:
                        ref_writer = new ParquetWriter(string(<const char*>filename_buf), compression_level)
                    ref_writer_active = True

                ref_writer.WriteBatch(ref_batch)
                ref_batch.clear()
                ref_batch.shrink_to_fit()  # Force memory release

        if write_by_read:
            read_batch.read_ids.push_back(total_records)
            read_batch.read_names.push_back(string(bam_get_qname(aln)))
            read_batch.ref_ids.push_back(aln.core.tid)
            read_batch.ref_names.push_back(string(ref_name))
            read_batch.positions.push_back(aln.core.pos)
            read_batch.end_positions.push_back(bam_endpos(aln))
            read_batch.mapqs.push_back(aln.core.qual)
            read_batch.flags.push_back(aln.core.flag)
            read_batch.alignment_lengths.push_back(aln.core.l_qseq)
            read_batch.template_lengths.push_back(aln.core.isize)
            read_batch.mate_ref_ids.push_back(aln.core.mtid if aln.core.mtid >= 0 else -1)
            read_batch.mate_positions.push_back(aln.core.mpos if aln.core.mpos >= 0 else -1)
            read_batch.alignment_scores.push_back(alignment_score)
            read_batch.xs_scores.push_back(xs_score)
            read_batch.edit_distances.push_back(edit_dist)
            read_batch.num_mismatches.push_back(num_mm)
            read_batch.num_gap_opens.push_back(num_go)
            read_batch.num_gap_extensions.push_back(num_ge)
            read_batch.md_strings.push_back(string(md_str))
            read_batch.anis.push_back(ani_val)
            read_batch.pmd_scores.push_back(pmd_score)
            read_batch.zs_scores.push_back(zs_score)
            read_batch.zp_posteriors.push_back(zp_posterior)
            read_batch.lca_taxids.push_back(lca_taxid)
            read_batch.reassigned_ref_ids.push_back(reassigned_ref_id)
            read_batch.filter_passed.push_back(True)
            read_batch.read_groups.push_back(string(rg_str))
            read_batch.cigars.push_back(string(<const char*>cigar_buf))
            read_batch.sequences.push_back(string(<const char*>seq_buf, <size_t>aln.core.l_qseq))
            read_batch.qualities.push_back(string(<char*>bam_get_qual(aln), aln.core.l_qseq))

            # Flush if batch full
            if read_batch.size() >= <size_t>batch_size:
                # Lazy create writer
                if not read_writer_active:
                    sprintf(filename_buf, "%s/chunk_%03d/by_read.parquet",
                            output_dir, thread_id)
                    with gil:
                        read_writer = new ParquetWriter(string(<const char*>filename_buf), compression_level)
                    read_writer_active = True

                read_writer.WriteBatch(read_batch)
                read_batch.clear()
                read_batch.shrink_to_fit()  # Force memory release

    # Flush remaining batches and cleanup
    if write_by_reference:
        if ref_batch.size() > 0:
            if not ref_writer_active:
                sprintf(filename_buf, "%s/chunk_%03d/by_reference.parquet",
                        output_dir, thread_id)
                with gil:
                    ref_writer = new ParquetWriter(string(<const char*>filename_buf), compression_level)
                ref_writer_active = True
            ref_writer.WriteBatch(ref_batch)
            ref_batch.clear()
            ref_batch.shrink_to_fit()  # Force memory release

        if ref_writer_active:
            ref_writer.Close()
            del ref_writer

    if write_by_read:
        if read_batch.size() > 0:
            if not read_writer_active:
                sprintf(filename_buf, "%s/chunk_%03d/by_read.parquet",
                        output_dir, thread_id)
                with gil:
                    read_writer = new ParquetWriter(string(<const char*>filename_buf), compression_level)
                read_writer_active = True
            read_writer.WriteBatch(read_batch)
            read_batch.clear()
            read_batch.shrink_to_fit()  # Force memory release

        if read_writer_active:
            read_writer.Close()
            del read_writer

    bam_destroy1(aln)
    bam_hdr_destroy(header)
    sam_close(sam_fp)

    return records_processed
