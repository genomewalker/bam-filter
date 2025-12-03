# cython: initializedcheck=False
# cython: embedsignature=False
# cython: binding=True
# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# cython: nonecheck=False
# cython: infer_types=True
# distutils: language = c++
"""
Fast BAM/SAM to Parquet conversion pipeline.

Single-pass streaming conversion with:
- Automatic format detection (SAM/SAM.gz/BAM)
- Dual-table writing (by_reference + by_read)
- Efficient memory usage (<500MB)
- Multi-threaded compression
"""

from libc.stdlib cimport malloc, free
from libc.string cimport strcpy, strlen
from libc.stdint cimport uint64_t, uint32_t, uint16_t, uint8_t, int32_t

from bam_filter.fast_alignment_reader cimport (
    AlignmentReader, alignment_reader_open, alignment_reader_read,
    alignment_reader_close, extract_read_name, extract_sequence,
    extract_quality, extract_cigar_string, get_tag_int, get_tag_float,
    bam1_t, bam_hdr_t
)
from bam_filter.fast_parquet_writer cimport (
    ParquetWriter, ParquetWriterConfig, AlignmentRecord,
    parquet_writer_create, parquet_writer_add_record,
    parquet_writer_flush, parquet_writer_close,
    calculate_ani, calculate_gc_content
)

import time
from pathlib import Path

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil
    double bf_monotonic_seconds() nogil


# ============================================================================
# Main Conversion Function
# ============================================================================

def convert_sam_bam_to_parquet(
    str input_file,
    str output_dir,
    int num_partitions=256,
    int batch_size=100000,
    str compression="zstd",
    int compression_level=3,
    bint write_by_reference=True,
    bint write_by_read=True,
    bint include_read_names=True,
    bint include_sequences=True,
    bint calculate_pmd=False,
    int num_threads=1,
):
    """
    Convert SAM/SAM.gz/BAM to dual-table Parquet format.

    Parameters
    ----------
    input_file : str
        Input alignment file (SAM, SAM.gz, or BAM)
    output_dir : str
        Output directory for Parquet files
    num_partitions : int
        Number of hash partitions (default: 256)
    batch_size : int
        Records per batch before writing (default: 100000)
    compression : str
        Compression codec: zstd, snappy, gzip, none (default: zstd)
    compression_level : int
        Compression level for zstd (default: 3)
    write_by_reference : bool
        Write alignments_by_reference table (default: True)
    write_by_read : bool
        Write alignments_by_read table (default: True)
    include_read_names : bool
        Include read names (adds ~30-50% size) (default: True)
    include_sequences : bool
        Include sequences (default: True)
    calculate_pmd : bool
        Calculate PMD scores (default: False)
    num_threads : int
        Number of threads for compression (default: 1)

    Returns
    -------
    dict
        Statistics: total_records, records_by_ref, records_by_read, duration_seconds
    """
    cdef bytes input_file_bytes = input_file.encode('utf-8')
    cdef bytes output_dir_bytes = output_dir.encode('utf-8')
    cdef bytes compression_bytes = compression.encode('utf-8')

    cdef double start_time = time.perf_counter()

    # Call Cython implementation
    cdef dict stats = _convert_impl(
        input_file_bytes,
        output_dir_bytes,
        num_partitions,
        batch_size,
        compression_bytes,
        compression_level,
        write_by_reference,
        write_by_read,
        include_read_names,
        include_sequences,
        calculate_pmd,
        num_threads
    )

    cdef double duration = time.perf_counter() - start_time
    stats['duration_seconds'] = duration

    return stats


cdef dict _convert_impl(
    bytes input_file,
    bytes output_dir,
    int num_partitions,
    int batch_size,
    bytes compression,
    int compression_level,
    bint write_by_reference,
    bint write_by_read,
    bint include_read_names,
    bint include_sequences,
    bint calculate_pmd,
    int num_threads
):
    """Cython implementation of conversion (can release GIL)."""

    # Prepare config
    cdef ParquetWriterConfig config
    config.output_path = <const char*>output_dir
    config.num_partitions = num_partitions
    config.batch_size = batch_size
    config.compression = <const char*>compression
    config.compression_level = compression_level
    config.write_by_reference = write_by_reference
    config.write_by_read = write_by_read
    config.include_read_names = include_read_names
    config.include_sequences = include_sequences
    config.num_threads = num_threads

    # Create output directories
    output_path = Path(output_dir.decode('utf-8'))
    output_path.mkdir(parents=True, exist_ok=True)

    if write_by_reference:
        (output_path / "alignments_by_reference").mkdir(exist_ok=True)

    if write_by_read:
        (output_path / "alignments_by_read").mkdir(exist_ok=True)

    # Process in nogil section for maximum performance
    cdef uint64_t total_records = 0
    cdef uint64_t records_by_ref = 0
    cdef uint64_t records_by_read = 0
    cdef double processing_time = 0.0
    cdef const char* input_file_cstr = <const char*>input_file

    with nogil:
        _process_alignments_nogil(
            input_file_cstr,
            &config,
            calculate_pmd,
            &total_records,
            &records_by_ref,
            &records_by_read,
            &processing_time
        )

    # Build stats dict with GIL
    return {
        'total_records': total_records,
        'records_by_reference': records_by_ref,
        'records_by_read': records_by_read,
        'processing_time_seconds': processing_time,
    }


cdef void _process_alignments_nogil(
    const char* input_file,
    ParquetWriterConfig* config,
    bint calculate_pmd,
    uint64_t* out_total_records,
    uint64_t* out_records_by_ref,
    uint64_t* out_records_by_read,
    double* out_processing_time
) nogil:
    """Main processing loop (nogil for performance)."""

    cdef double start_time = bf_monotonic_seconds()

    # Open reader
    cdef AlignmentReader* reader = alignment_reader_open(input_file)
    if reader == NULL:
        with gil:
            raise IOError(f"Failed to open input file: {input_file.decode('utf-8')}")

    # Create writer
    cdef ParquetWriter* writer = parquet_writer_create(config)
    if writer == NULL:
        alignment_reader_close(reader)
        with gil:
            raise IOError("Failed to create Parquet writer")

    # Allocate buffers for variable-length data
    cdef char* read_name_buf = <char*>malloc(512)
    cdef char* cigar_buf = <char*>malloc(2048)
    cdef char* sequence_buf = <char*>malloc(65536)
    cdef char* quality_buf = <char*>malloc(65536)

    if read_name_buf == NULL or cigar_buf == NULL or \
       sequence_buf == NULL or quality_buf == NULL:
        # Cleanup on allocation failure
        if read_name_buf != NULL:
            free(read_name_buf)
        if cigar_buf != NULL:
            free(cigar_buf)
        if sequence_buf != NULL:
            free(sequence_buf)
        if quality_buf != NULL:
            free(quality_buf)
        parquet_writer_close(writer)
        alignment_reader_close(reader)
        with gil:
            raise MemoryError("Failed to allocate buffers")

    # Processing loop
    cdef bam1_t* record
    cdef int read_result
    cdef AlignmentRecord aln_record
    cdef uint64_t records_processed = 0
    cdef uint64_t records_written = 0
    cdef double last_log_time = start_time
    cdef double current_time

    bf_nogil_logf_verbose(2, "CONVERT", "Starting conversion: %s\n", input_file)

    while True:
        read_result = alignment_reader_read(reader, &record)

        if read_result == 0:
            # EOF
            break
        elif read_result < 0:
            # Error
            bf_nogil_logf_verbose(1, "CONVERT", "Error reading record\n")
            break

        # Skip unmapped reads
        if record.core.tid < 0:
            continue

        # Extract fields from BAM record and populate AlignmentRecord
        _populate_alignment_record(
            &aln_record, record, records_processed,
            read_name_buf, cigar_buf, sequence_buf, quality_buf,
            config, calculate_pmd
        )

        # Add record to writer
        if parquet_writer_add_record(writer, &aln_record) < 0:
            bf_nogil_logf_verbose(1, "CONVERT", "Error writing record\n")
            break

        records_processed += 1
        records_written += 1

        # Log progress every 10 seconds
        current_time = bf_monotonic_seconds()
        if current_time - last_log_time >= 10.0:
            bf_nogil_logf_verbose(2, "CONVERT",
                                 "Processed %llu records (%.1f M/sec)\n",
                                 records_processed,
                                 records_processed / (current_time - start_time) / 1000000.0)
            last_log_time = current_time

    # Flush and close
    parquet_writer_flush(writer)

    cdef uint64_t final_by_ref = writer.records_by_reference
    cdef uint64_t final_by_read = writer.records_by_read

    parquet_writer_close(writer)
    alignment_reader_close(reader)

    # Free buffers
    free(read_name_buf)
    free(cigar_buf)
    free(sequence_buf)
    free(quality_buf)

    cdef double end_time = bf_monotonic_seconds()
    cdef double duration = end_time - start_time

    bf_nogil_logf_verbose(2, "CONVERT",
                         "Conversion complete: %llu records in %.2f seconds (%.1f M/sec)\n",
                         records_written, duration,
                         records_written / duration / 1000000.0)

    # Set output parameters
    out_total_records[0] = records_written
    out_records_by_ref[0] = final_by_ref
    out_records_by_read[0] = final_by_read
    out_processing_time[0] = duration


cdef void _populate_alignment_record(
    AlignmentRecord* aln_record,
    bam1_t* bam_record,
    uint64_t read_id,
    char* read_name_buf,
    char* cigar_buf,
    char* sequence_buf,
    char* quality_buf,
    ParquetWriterConfig* config,
    bint calculate_pmd
) nogil:
    """Populate AlignmentRecord from BAM record."""

    # IDs
    aln_record.read_id = read_id
    aln_record.ref_id = bam_record.core.tid

    # Position
    aln_record.position = bam_record.core.pos

    # Calculate end position from CIGAR
    cdef int ref_span = 0
    cdef uint32_t* cigar = bam_get_cigar(bam_record)
    cdef int n_cigar = bam_record.core.n_cigar
    cdef int i
    cdef int op

    for i in range(n_cigar):
        op = cigar[i] & 0xf
        # Count M, D, N, =, X operations (consume reference)
        if op == 0 or op == 2 or op == 3 or op == 7 or op == 8:
            ref_span += cigar[i] >> 4

    aln_record.end_position = bam_record.core.pos + ref_span
    aln_record.reference_span = ref_span

    # Quality
    aln_record.mapq = bam_record.core.qual
    aln_record.flag = bam_record.core.flag

    # Get alignment score (AS tag or calculate)
    cdef int32_t alignment_score_int = get_tag_int(bam_record, "AS", 0)
    aln_record.alignment_score = <float>alignment_score_int

    # Get edit distance (NM tag)
    aln_record.edit_distance = <uint16_t>get_tag_int(bam_record, "NM", 0)

    # Calculate alignment length and ANI
    aln_record.alignment_length = bam_record.core.l_qseq
    aln_record.ani = calculate_ani(NULL, aln_record.edit_distance, aln_record.alignment_length)

    # PMD score (optional)
    if calculate_pmd:
        aln_record.pmd_score = get_tag_float(bam_record, "PMD", -1.0)
    else:
        aln_record.pmd_score = -1.0

    # Paired-end info
    aln_record.template_length = bam_record.core.isize
    aln_record.mate_ref_id = bam_record.core.mtid
    aln_record.mate_position = bam_record.core.mpos

    # Extract variable-length data
    if config.include_read_names:
        aln_record.read_name_len = extract_read_name(bam_record, read_name_buf, 512)
        aln_record.read_name = read_name_buf
    else:
        aln_record.read_name = NULL
        aln_record.read_name_len = 0

    aln_record.cigar_len = extract_cigar_string(bam_record, cigar_buf, 2048)
    aln_record.cigar = cigar_buf

    if config.include_sequences:
        aln_record.sequence_len = extract_sequence(bam_record, sequence_buf, 65536)
        aln_record.sequence = sequence_buf

        aln_record.quality_len = extract_quality(bam_record, quality_buf, 65536)
        aln_record.quality = quality_buf

        # Calculate GC content
        aln_record.gc_content = calculate_gc_content(sequence_buf, aln_record.sequence_len)
    else:
        aln_record.sequence = NULL
        aln_record.sequence_len = 0
        aln_record.quality = NULL
        aln_record.quality_len = 0
        aln_record.gc_content = -1.0

    # DUST score (not calculated by default)
    aln_record.dust_score = -1.0

    # Tags (raw binary)
    cdef uint8_t* aux_start = bam_get_aux(bam_record)
    cdef int aux_len = bam_record.l_data - (aux_start - bam_record.data)

    if aux_len > 0:
        aln_record.tags = aux_start
        aln_record.tags_len = aux_len
    else:
        aln_record.tags = NULL
        aln_record.tags_len = 0


# Helper function
cdef extern from "htslib/sam.h":
    uint32_t* bam_get_cigar(bam1_t* b) nogil
    uint8_t* bam_get_aux(bam1_t* b) nogil
