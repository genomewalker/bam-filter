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
Fast dual-table Parquet writer for metagenomic alignments.

Architecture:
1. Alignments partitioned by ref_id → alignments_by_reference/
2. Alignments partitioned by read_id → alignments_by_read/
3. Batched writing with configurable compression
4. Efficient memory usage via streaming
"""

from libc.stdlib cimport malloc, calloc, free, realloc
from libc.string cimport memcpy, memset, strlen
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int64_t
from libc.math cimport log, floor

# Python imports for Parquet writing
import pyarrow as pa
import pyarrow.parquet as pq
from pathlib import Path

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil


# ============================================================================
# Schema Definitions
# ============================================================================

def get_alignment_schema(bint include_read_names=True, bint include_sequences=True):
    """
    Get PyArrow schema for alignment records.
    Optimized column order: hot fields first, cold fields last.
    """
    fields = [
        # Hot fields (frequently filtered)
        pa.field("read_id", pa.uint64(), nullable=False),
        pa.field("ref_id", pa.uint32(), nullable=False),
        pa.field("position", pa.int32(), nullable=False),
        pa.field("end_position", pa.int32(), nullable=False),
        pa.field("alignment_score", pa.float32(), nullable=False),
        pa.field("ani", pa.float32(), nullable=False),
        pa.field("mapq", pa.uint8(), nullable=False),
        pa.field("flag", pa.uint16(), nullable=False),

        # Medium-hot fields (quality metrics)
        pa.field("edit_distance", pa.uint16(), nullable=False),
        pa.field("alignment_length", pa.uint16(), nullable=False),
        pa.field("reference_span", pa.uint16(), nullable=False),

        # Calculated metrics (nullable)
        pa.field("pmd_score", pa.float32(), nullable=True),
        pa.field("gc_content", pa.float32(), nullable=True),
        pa.field("dust_score", pa.float32(), nullable=True),

        # Paired-end info
        pa.field("template_length", pa.int32(), nullable=False),
        pa.field("mate_ref_id", pa.int32(), nullable=False),
        pa.field("mate_position", pa.int32(), nullable=False),
    ]

    # Optional: read names (adds ~30-50% size)
    if include_read_names:
        fields.append(pa.field("read_name", pa.binary(), nullable=False))

    # Cold fields (lossless data - large but rarely queried)
    fields.extend([
        pa.field("cigar", pa.binary(), nullable=False),
    ])

    if include_sequences:
        fields.extend([
            pa.field("sequence", pa.binary(), nullable=False),
            pa.field("quality", pa.binary(), nullable=False),
        ])

    fields.append(pa.field("tags", pa.binary(), nullable=True))

    return pa.schema(fields)


def get_reference_schema():
    """Get PyArrow schema for reference dimension table."""
    return pa.schema([
        pa.field("ref_id", pa.uint32(), nullable=False),
        pa.field("ref_name", pa.string(), nullable=False),
        pa.field("ref_length", pa.uint32(), nullable=False),
        pa.field("ref_partition", pa.uint16(), nullable=False),
    ])


# ============================================================================
# Batch Management
# ============================================================================

cdef ParquetBatch* create_batch(int capacity, int partition_id) except NULL nogil:
    """Create a new batch for accumulating records."""
    cdef ParquetBatch* batch = <ParquetBatch*>malloc(sizeof(ParquetBatch))
    if batch == NULL:
        return NULL

    batch.capacity = capacity
    batch.count = 0
    batch.partition_id = partition_id

    batch.records = <AlignmentRecord*>calloc(capacity, sizeof(AlignmentRecord))
    if batch.records == NULL:
        free(batch)
        return NULL

    return batch


cdef int batch_add_record(ParquetBatch* batch, AlignmentRecord* record) except -1 nogil:
    """
    Add a record to the batch.
    Note: This copies the record structure, but NOT the variable-length data.
    Caller must ensure variable-length data remains valid until batch is written.
    """
    cdef int64_t new_capacity
    cdef AlignmentRecord* new_records

    if batch == NULL or record == NULL:
        return -1

    if batch.count >= batch.capacity:
        # Grow batch if needed
        new_capacity = batch.capacity * 2
        new_records = <AlignmentRecord*>realloc(
            batch.records, new_capacity * sizeof(AlignmentRecord))

        if new_records == NULL:
            return -1

        batch.records = new_records
        batch.capacity = new_capacity

    # Copy record structure
    memcpy(&batch.records[batch.count], record, sizeof(AlignmentRecord))
    batch.count += 1

    return 0


cdef void batch_clear(ParquetBatch* batch) nogil:
    """Clear batch for reuse (doesn't free memory)."""
    if batch != NULL:
        batch.count = 0


cdef void batch_destroy(ParquetBatch* batch) nogil:
    """Destroy batch and free all memory."""
    if batch != NULL:
        if batch.records != NULL:
            free(batch.records)
        free(batch)


# ============================================================================
# Helper Functions
# ============================================================================

cdef int calculate_partition_id(uint64_t id_value, int num_partitions) nogil:
    """Calculate partition ID using hash function."""
    # Simple modulo hash - distributes evenly
    return <int>(id_value % <uint64_t>num_partitions)


cdef float calculate_ani(const char* cigar, int edit_distance, int alignment_length) nogil:
    """
    Calculate Alignment Identity (ANI) percentage.
    ANI = (1 - edit_distance / alignment_length) * 100
    """
    if alignment_length <= 0:
        return 0.0

    cdef float identity = 1.0 - (<float>edit_distance / <float>alignment_length)
    return identity * 100.0


cdef float calculate_gc_content(const char* sequence, int seq_len) nogil:
    """Calculate GC content percentage."""
    if sequence == NULL or seq_len <= 0:
        return 0.0

    cdef int gc_count = 0
    cdef int i

    for i in range(seq_len):
        if sequence[i] == b'G' or sequence[i] == b'C' or \
           sequence[i] == b'g' or sequence[i] == b'c':
            gc_count += 1

    return (<float>gc_count / <float>seq_len) * 100.0


# ============================================================================
# Batch Writing to Parquet
# ============================================================================

cdef int batch_write_to_parquet(ParquetBatch* batch, const char* output_file,
                                 ParquetWriterConfig* config, bint by_read) except -1 nogil:
    """
    Write batch to Parquet file.
    This function releases the GIL to call Python code.
    """
    if batch == NULL or batch.count == 0:
        return 0

    # We need to release nogil and call Python
    # This is done in a with gil block
    with gil:
        return _write_batch_python(batch, output_file, config, by_read)


cdef int _write_batch_python(ParquetBatch* batch, const char* output_file,
                              ParquetWriterConfig* config, bint by_read) except -1 with gil:
    """
    Write batch to Parquet file using PyArrow (requires GIL).
    """
    try:
        # Build PyArrow arrays from batch
        arrays = _build_pyarrow_arrays(batch, config)

        # Create table
        schema = get_alignment_schema(
            include_read_names=config.include_read_names,
            include_sequences=config.include_sequences
        )
        table = pa.Table.from_arrays(arrays, schema=schema)

        # Sort table
        if by_read:
            # Sort by read_id, then by alignment_score DESC for LCA queries
            table = table.sort_by([
                ("read_id", "ascending"),
                ("alignment_score", "descending")
            ])
        else:
            # Sort by ref_id, position for coverage queries
            table = table.sort_by([
                ("ref_id", "ascending"),
                ("position", "ascending"),
                ("read_id", "ascending")
            ])

        # Write Parquet file
        output_path = Path(output_file.decode('utf-8'))
        output_path.parent.mkdir(parents=True, exist_ok=True)

        compression_str = config.compression.decode('utf-8') if config.compression else 'zstd'

        pq.write_table(
            table,
            output_path,
            compression=compression_str,
            compression_level=config.compression_level,
            row_group_size=min(100000, batch.count),  # 100K rows per row group
            use_dictionary=True,  # Enable dictionary encoding
            write_statistics=True,  # Write column statistics for filtering
        )

        return 0

    except Exception as e:
        import sys
        print(f"Error writing Parquet: {e}", file=sys.stderr)
        return -1


def _build_pyarrow_arrays(ParquetBatch* batch, ParquetWriterConfig* config):
    """Build PyArrow arrays from batch records (called with GIL)."""
    cdef int64_t i
    cdef AlignmentRecord* rec

    # Pre-allocate Python lists
    read_ids = []
    ref_ids = []
    positions = []
    end_positions = []
    alignment_scores = []
    anis = []
    mapqs = []
    flags = []
    edit_distances = []
    alignment_lengths = []
    reference_spans = []
    pmd_scores = []
    gc_contents = []
    dust_scores = []
    template_lengths = []
    mate_ref_ids = []
    mate_positions = []
    read_names = [] if config.include_read_names else None
    cigars = []
    sequences = [] if config.include_sequences else None
    qualities = [] if config.include_sequences else None
    tags = []

    # Extract data from C structs
    for i in range(batch.count):
        rec = &batch.records[i]

        read_ids.append(rec.read_id)
        ref_ids.append(rec.ref_id)
        positions.append(rec.position)
        end_positions.append(rec.end_position)
        alignment_scores.append(rec.alignment_score)
        anis.append(rec.ani)
        mapqs.append(rec.mapq)
        flags.append(rec.flag)
        edit_distances.append(rec.edit_distance)
        alignment_lengths.append(rec.alignment_length)
        reference_spans.append(rec.reference_span)

        # Nullable fields
        pmd_scores.append(rec.pmd_score if rec.pmd_score >= 0 else None)
        gc_contents.append(rec.gc_content if rec.gc_content >= 0 else None)
        dust_scores.append(rec.dust_score if rec.dust_score >= 0 else None)

        template_lengths.append(rec.template_length)
        mate_ref_ids.append(rec.mate_ref_id)
        mate_positions.append(rec.mate_position)

        # Variable-length fields
        if config.include_read_names and rec.read_name != NULL:
            read_names.append(bytes(rec.read_name[:rec.read_name_len]))

        if rec.cigar != NULL:
            cigars.append(bytes(rec.cigar[:rec.cigar_len]))

        if config.include_sequences:
            if rec.sequence != NULL:
                sequences.append(bytes(rec.sequence[:rec.sequence_len]))
            if rec.quality != NULL:
                qualities.append(bytes(rec.quality[:rec.quality_len]))

        if rec.tags != NULL and rec.tags_len > 0:
            tags.append(bytes(rec.tags[:rec.tags_len]))
        else:
            tags.append(None)

    # Build PyArrow arrays
    arrays = [
        pa.array(read_ids, type=pa.uint64()),
        pa.array(ref_ids, type=pa.uint32()),
        pa.array(positions, type=pa.int32()),
        pa.array(end_positions, type=pa.int32()),
        pa.array(alignment_scores, type=pa.float32()),
        pa.array(anis, type=pa.float32()),
        pa.array(mapqs, type=pa.uint8()),
        pa.array(flags, type=pa.uint16()),
        pa.array(edit_distances, type=pa.uint16()),
        pa.array(alignment_lengths, type=pa.uint16()),
        pa.array(reference_spans, type=pa.uint16()),
        pa.array(pmd_scores, type=pa.float32()),
        pa.array(gc_contents, type=pa.float32()),
        pa.array(dust_scores, type=pa.float32()),
        pa.array(template_lengths, type=pa.int32()),
        pa.array(mate_ref_ids, type=pa.int32()),
        pa.array(mate_positions, type=pa.int32()),
    ]

    if config.include_read_names:
        arrays.append(pa.array(read_names, type=pa.binary()))

    arrays.append(pa.array(cigars, type=pa.binary()))

    if config.include_sequences:
        arrays.append(pa.array(sequences, type=pa.binary()))
        arrays.append(pa.array(qualities, type=pa.binary()))

    arrays.append(pa.array(tags, type=pa.binary()))

    return arrays


# ============================================================================
# Writer API
# ============================================================================

cdef ParquetWriter* parquet_writer_create(ParquetWriterConfig* config) except NULL nogil:
    """Create a new Parquet writer with dual-table support."""
    cdef ParquetWriter* writer = <ParquetWriter*>malloc(sizeof(ParquetWriter))
    if writer == NULL:
        return NULL

    # Copy config
    memcpy(&writer.config, config, sizeof(ParquetWriterConfig))

    # Initialize statistics
    writer.total_records_written = 0
    writer.total_bytes_written = 0
    writer.records_by_reference = 0
    writer.records_by_read = 0

    writer.num_references = 0
    writer.reference_names = NULL
    writer.reference_lengths = NULL

    # Allocate partition batches
    cdef int num_partitions = config.num_partitions
    cdef int i

    if config.write_by_reference:
        writer.ref_batches = <ParquetBatch**>calloc(num_partitions, sizeof(ParquetBatch*))
        if writer.ref_batches == NULL:
            free(writer)
            return NULL

        for i in range(num_partitions):
            writer.ref_batches[i] = create_batch(config.batch_size, i)
            if writer.ref_batches[i] == NULL:
                # Cleanup on failure
                for j in range(i):
                    batch_destroy(writer.ref_batches[j])
                free(writer.ref_batches)
                free(writer)
                return NULL
    else:
        writer.ref_batches = NULL

    if config.write_by_read:
        writer.read_batches = <ParquetBatch**>calloc(num_partitions, sizeof(ParquetBatch*))
        if writer.read_batches == NULL:
            if writer.ref_batches != NULL:
                for i in range(num_partitions):
                    batch_destroy(writer.ref_batches[i])
                free(writer.ref_batches)
            free(writer)
            return NULL

        for i in range(num_partitions):
            writer.read_batches[i] = create_batch(config.batch_size, i)
            if writer.read_batches[i] == NULL:
                # Cleanup on failure
                for j in range(i):
                    batch_destroy(writer.read_batches[j])
                free(writer.read_batches)
                if writer.ref_batches != NULL:
                    for j in range(num_partitions):
                        batch_destroy(writer.ref_batches[j])
                    free(writer.ref_batches)
                free(writer)
                return NULL
    else:
        writer.read_batches = NULL

    bf_nogil_logf_verbose(2, "PARQUET", "Created writer: %d partitions, batch_size=%d\n",
                         num_partitions, config.batch_size)

    return writer


cdef int parquet_writer_add_record(ParquetWriter* writer, AlignmentRecord* record) except -1 nogil:
    """Add a record to the writer (will be added to appropriate partition batches)."""
    if writer == NULL or record == NULL:
        return -1

    cdef int ref_partition_id
    cdef int read_partition_id
    cdef int result

    # Add to reference-partitioned table
    if writer.config.write_by_reference and writer.ref_batches != NULL:
        ref_partition_id = calculate_partition_id(record.ref_id, writer.config.num_partitions)
        result = batch_add_record(writer.ref_batches[ref_partition_id], record)
        if result < 0:
            return -1

        # Flush if batch is full
        if writer.ref_batches[ref_partition_id].count >= writer.config.batch_size:
            result = _flush_ref_batch(writer, ref_partition_id)
            if result < 0:
                return -1

    # Add to read-partitioned table
    if writer.config.write_by_read and writer.read_batches != NULL:
        read_partition_id = calculate_partition_id(record.read_id, writer.config.num_partitions)
        result = batch_add_record(writer.read_batches[read_partition_id], record)
        if result < 0:
            return -1

        # Flush if batch is full
        if writer.read_batches[read_partition_id].count >= writer.config.batch_size:
            result = _flush_read_batch(writer, read_partition_id)
            if result < 0:
                return -1

    writer.total_records_written += 1
    return 0


cdef int _flush_ref_batch(ParquetWriter* writer, int partition_id) except -1 nogil:
    """Flush a reference-partitioned batch to disk."""
    if writer == NULL or writer.ref_batches == NULL:
        return -1

    cdef ParquetBatch* batch = writer.ref_batches[partition_id]
    if batch == NULL or batch.count == 0:
        return 0

    # Build output filename
    cdef char output_file[512]
    cdef int chunk_id = <int>(writer.records_by_reference / writer.config.batch_size)

    snprintf(output_file, sizeof(output_file),
             "%s/alignments_by_reference/ref_partition=%04d/chunk_%06d.parquet",
             writer.config.output_path, partition_id, chunk_id)

    # Write to Parquet
    cdef int result = batch_write_to_parquet(batch, output_file, &writer.config, False)
    if result < 0:
        return -1

    writer.records_by_reference += batch.count
    batch_clear(batch)

    return 0


cdef int _flush_read_batch(ParquetWriter* writer, int partition_id) except -1 nogil:
    """Flush a read-partitioned batch to disk."""
    if writer == NULL or writer.read_batches == NULL:
        return -1

    cdef ParquetBatch* batch = writer.read_batches[partition_id]
    if batch == NULL or batch.count == 0:
        return 0

    # Build output filename
    cdef char output_file[512]
    cdef int chunk_id = <int>(writer.records_by_read / writer.config.batch_size)

    snprintf(output_file, sizeof(output_file),
             "%s/alignments_by_read/read_partition=%04d/chunk_%06d.parquet",
             writer.config.output_path, partition_id, chunk_id)

    # Write to Parquet
    cdef int result = batch_write_to_parquet(batch, output_file, &writer.config, True)
    if result < 0:
        return -1

    writer.records_by_read += batch.count
    batch_clear(batch)

    return 0


cdef int parquet_writer_flush(ParquetWriter* writer) except -1 nogil:
    """Flush all pending batches to disk."""
    if writer == NULL:
        return -1

    cdef int i
    cdef int result

    # Flush all reference batches
    if writer.ref_batches != NULL:
        for i in range(writer.config.num_partitions):
            result = _flush_ref_batch(writer, i)
            if result < 0:
                return -1

    # Flush all read batches
    if writer.read_batches != NULL:
        for i in range(writer.config.num_partitions):
            result = _flush_read_batch(writer, i)
            if result < 0:
                return -1

    bf_nogil_logf_verbose(2, "PARQUET",
                         "Flushed all batches: %llu by_reference, %llu by_read\n",
                         writer.records_by_reference, writer.records_by_read)

    return 0


cdef int parquet_writer_close(ParquetWriter* writer) except -1 nogil:
    """Close writer and free all resources."""
    if writer == NULL:
        return 0

    # Flush any pending data
    parquet_writer_flush(writer)

    # Free batches
    cdef int i
    if writer.ref_batches != NULL:
        for i in range(writer.config.num_partitions):
            batch_destroy(writer.ref_batches[i])
        free(writer.ref_batches)

    if writer.read_batches != NULL:
        for i in range(writer.config.num_partitions):
            batch_destroy(writer.read_batches[i])
        free(writer.read_batches)

    # Free reference metadata
    if writer.reference_names != NULL:
        for i in range(writer.num_references):
            if writer.reference_names[i] != NULL:
                free(writer.reference_names[i])
        free(writer.reference_names)

    if writer.reference_lengths != NULL:
        free(writer.reference_lengths)

    bf_nogil_logf_verbose(2, "PARQUET", "Closed writer: %llu total records\n",
                         writer.total_records_written)

    free(writer)
    return 0


# C function for snprintf
cdef extern from "stdio.h":
    int snprintf(char* s, size_t n, const char* format, ...) nogil
