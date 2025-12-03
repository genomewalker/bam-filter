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
Fast unified alignment reader for SAM, SAM.gz, and BAM formats.

Features:
- Automatic format detection
- Unified API for all formats
- Zero-copy field extraction
- HTSlib-based for maximum compatibility
"""

from libc.stdlib cimport malloc, free
from libc.string cimport memcpy, strlen, strcmp, strstr
from libc.stdint cimport uint8_t, uint16_t, uint32_t, int32_t, uint64_t, int8_t, int16_t
from libc.stdio cimport FILE, fopen, fclose, fread

cdef extern from "bam_filter/c_logging.h":
    void bf_nogil_logf_verbose(int level, const char* tag, const char* fmt, ...) nogil


# ============================================================================
# Format Detection
# ============================================================================

cdef AlignmentFormat detect_format(const char* path) nogil:
    """
    Detect alignment file format from path and magic bytes.

    Detection order:
    1. Check file extension (.sam, .sam.gz, .bam)
    2. Check magic bytes (gzip: 0x1f8b, BAM: 0x1f8b0408 + "BAM")
    """
    cdef FILE* fp = NULL
    cdef unsigned char magic[4]
    cdef size_t bytes_read
    cdef AlignmentFormat format = FORMAT_UNKNOWN

    # Check extension first
    cdef const char* ext = strstr(path, ".bam")
    if ext != NULL and strcmp(ext, ".bam") == 0:
        return FORMAT_BAM

    ext = strstr(path, ".sam.gz")
    if ext != NULL and strcmp(ext, ".sam.gz") == 0:
        return FORMAT_SAM_GZ

    ext = strstr(path, ".sam")
    if ext != NULL and strcmp(ext, ".sam") == 0:
        return FORMAT_SAM

    # Check magic bytes for gzip/BAM
    fp = fopen(path, "rb")
    if fp == NULL:
        return FORMAT_UNKNOWN

    bytes_read = fread(magic, 1, 4, fp)
    fclose(fp)

    if bytes_read < 4:
        return FORMAT_UNKNOWN

    # Check for gzip magic (0x1f8b)
    if magic[0] == 0x1f and magic[1] == 0x8b:
        # Could be BAM or SAM.gz - assume BAM if no .sam.gz extension
        # HTSlib will handle it correctly either way
        return FORMAT_BAM

    # Default to SAM if no magic bytes match
    return FORMAT_SAM


# ============================================================================
# Reader API
# ============================================================================

cdef AlignmentReader* alignment_reader_open(const char* path) except NULL nogil:
    """
    Open an alignment file for reading.
    Automatically detects format and uses appropriate reader.

    Returns NULL on error.
    """
    cdef AlignmentReader* reader = <AlignmentReader*>malloc(sizeof(AlignmentReader))
    if reader == NULL:
        return NULL

    # Initialize
    reader.format = FORMAT_UNKNOWN
    reader.hts_file = NULL
    reader.header = NULL
    reader.current_record = NULL
    reader.records_read = 0
    reader.is_open = False

    # Detect format
    reader.format = detect_format(path)
    if reader.format == FORMAT_UNKNOWN:
        bf_nogil_logf_verbose(1, "READER", "Unknown format for: %s\n", path)
        free(reader)
        return NULL

    # Open with HTSlib
    # HTSlib automatically handles SAM/SAM.gz/BAM based on content
    reader.hts_file = sam_open(path, "r")
    if reader.hts_file == NULL:
        bf_nogil_logf_verbose(1, "READER", "Failed to open: %s\n", path)
        free(reader)
        return NULL

    # Read header
    reader.header = sam_hdr_read(reader.hts_file)
    if reader.header == NULL:
        bf_nogil_logf_verbose(1, "READER", "Failed to read header: %s\n", path)
        sam_close(reader.hts_file)
        free(reader)
        return NULL

    # Allocate record buffer
    reader.current_record = bam_init1()
    if reader.current_record == NULL:
        bf_nogil_logf_verbose(1, "READER", "Failed to allocate record buffer\n")
        bam_hdr_destroy(reader.header)
        sam_close(reader.hts_file)
        free(reader)
        return NULL

    reader.is_open = True

    cdef const char* format_name = "UNKNOWN"
    if reader.format == FORMAT_SAM:
        format_name = "SAM"
    elif reader.format == FORMAT_SAM_GZ:
        format_name = "SAM.gz"
    elif reader.format == FORMAT_BAM:
        format_name = "BAM"

    bf_nogil_logf_verbose(2, "READER", "Opened %s file: %s (%d references)\n",
                         format_name, path, reader.header.n_targets)

    return reader


cdef int alignment_reader_read(AlignmentReader* reader, bam1_t** record) except -1 nogil:
    """
    Read next alignment record.

    Returns:
        1 if record read successfully (*record points to valid data)
        0 if EOF
        -1 on error
    """
    if reader == NULL or not reader.is_open:
        return -1

    cdef int ret = sam_read1(reader.hts_file, reader.header, reader.current_record)

    if ret >= 0:
        # Success - record read
        record[0] = reader.current_record
        reader.records_read += 1
        return 1
    elif ret == -1:
        # EOF
        record[0] = NULL
        return 0
    else:
        # Error
        bf_nogil_logf_verbose(1, "READER", "Error reading record (code %d)\n", ret)
        record[0] = NULL
        return -1


cdef int alignment_reader_close(AlignmentReader* reader) except -1 nogil:
    """Close reader and free all resources."""
    if reader == NULL:
        return 0

    if reader.is_open:
        bf_nogil_logf_verbose(2, "READER", "Closing reader (read %llu records)\n",
                             reader.records_read)

        if reader.current_record != NULL:
            bam_destroy1(reader.current_record)

        if reader.header != NULL:
            bam_hdr_destroy(reader.header)

        if reader.hts_file != NULL:
            sam_close(reader.hts_file)

        reader.is_open = False

    free(reader)
    return 0


# ============================================================================
# Field Extraction Helpers
# ============================================================================

cdef int extract_read_name(bam1_t* record, char* buffer, int max_len) nogil:
    """
    Extract read name (QNAME) from BAM record.
    Returns length of name copied, or -1 on error.
    """
    if record == NULL or buffer == NULL or max_len <= 0:
        return -1

    cdef char* qname = bam_get_qname(record)
    if qname == NULL:
        return -1

    cdef int qname_len = record.core.l_qname - 1  # Exclude null terminator
    if qname_len >= max_len:
        qname_len = max_len - 1

    memcpy(buffer, qname, qname_len)
    buffer[qname_len] = 0

    return qname_len


cdef int extract_sequence(bam1_t* record, char* buffer, int max_len) nogil:
    """
    Extract DNA sequence from BAM record.
    Converts 4-bit encoding to ASCII (ACGTN).
    Returns sequence length, or -1 on error.
    """
    if record == NULL or buffer == NULL or max_len <= 0:
        return -1

    cdef int seq_len = record.core.l_qseq
    if seq_len >= max_len:
        seq_len = max_len - 1

    cdef uint8_t* seq = bam_get_seq(record)
    if seq == NULL:
        return -1

    # Decode 4-bit sequence to ASCII
    # BAM encoding: =ACMGRSVTWYHKDBN
    cdef const char* seq_nt16_str = "=ACMGRSVTWYHKDBN"
    cdef int i

    for i in range(seq_len):
        buffer[i] = seq_nt16_str[seq[i >> 1] >> ((~i & 1) << 2) & 0xf]

    buffer[seq_len] = 0
    return seq_len


cdef int extract_quality(bam1_t* record, char* buffer, int max_len) nogil:
    """
    Extract quality scores from BAM record.
    Converts to ASCII (Phred+33).
    Returns quality length, or -1 on error.
    """
    if record == NULL or buffer == NULL or max_len <= 0:
        return -1

    cdef int qual_len = record.core.l_qseq
    if qual_len >= max_len:
        qual_len = max_len - 1

    cdef uint8_t* qual = bam_get_qual(record)
    if qual == NULL or qual[0] == 0xff:
        # No quality scores
        buffer[0] = b'*'
        buffer[1] = 0
        return 1

    # Convert to Phred+33 ASCII
    cdef int i
    for i in range(qual_len):
        buffer[i] = qual[i] + 33

    buffer[qual_len] = 0
    return qual_len


cdef int extract_cigar_string(bam1_t* record, char* buffer, int max_len) nogil:
    """
    Extract CIGAR string from BAM record.
    Converts binary CIGAR to string representation.
    Returns CIGAR string length, or -1 on error.
    """
    if record == NULL or buffer == NULL or max_len <= 0:
        return -1

    cdef uint32_t* cigar = bam_get_cigar(record)
    cdef int n_cigar = record.core.n_cigar

    if n_cigar == 0:
        buffer[0] = b'*'
        buffer[1] = 0
        return 1

    cdef const char* cigar_ops = "MIDNSHP=X"
    cdef int pos = 0
    cdef int i
    cdef uint32_t op_len
    cdef int op
    cdef char num_buf[16]
    cdef int num_len

    for i in range(n_cigar):
        if pos >= max_len - 20:  # Safety margin
            break

        op_len = cigar[i] >> 4
        op = cigar[i] & 0xf

        # Convert number to string (simple integer to ASCII)
        num_len = snprintf(num_buf, sizeof(num_buf), "%u", op_len)

        if pos + num_len + 1 >= max_len:
            break

        memcpy(buffer + pos, num_buf, num_len)
        pos += num_len

        buffer[pos] = cigar_ops[op]
        pos += 1

    buffer[pos] = 0
    return pos


cdef uint8_t* find_tag(bam1_t* record, const char* tag) nogil:
    """
    Find auxiliary tag in BAM record.
    Returns pointer to tag data, or NULL if not found.
    """
    cdef uint8_t* aux
    cdef uint8_t* end
    cdef uint8_t tag_type
    cdef uint32_t array_len
    cdef uint8_t elem_type
    cdef int elem_size

    if record == NULL or tag == NULL:
        return NULL

    aux = bam_get_aux(record)
    end = record.data + record.l_data

    while aux < end:
        # Check tag name (2 bytes)
        if aux[0] == tag[0] and aux[1] == tag[1]:
            return aux + 2  # Return pointer to tag type

        # Skip to next tag
        aux += 2  # Tag name
        tag_type = aux[0]
        aux += 1  # Tag type

        # Skip tag data based on type
        if tag_type == b'A' or tag_type == b'c' or tag_type == b'C':
            aux += 1
        elif tag_type == b's' or tag_type == b'S':
            aux += 2
        elif tag_type == b'i' or tag_type == b'I' or tag_type == b'f':
            aux += 4
        elif tag_type == b'Z' or tag_type == b'H':
            # Null-terminated string
            while aux < end and aux[0] != 0:
                aux += 1
            aux += 1  # Skip null
        elif tag_type == b'B':
            # Array
            elem_type = aux[1]  # Element type (next byte)
            aux += 2  # Skip element type
            array_len = (<uint32_t*>aux)[0]
            aux += 4
            if elem_type == b'c' or elem_type == b'C':
                aux += array_len
            elif elem_type == b's' or elem_type == b'S':
                aux += array_len * 2
            elif elem_type == b'i' or elem_type == b'I' or elem_type == b'f':
                aux += array_len * 4

    return NULL


cdef int32_t get_tag_int(bam1_t* record, const char* tag, int32_t default_val) nogil:
    """
    Get integer value from auxiliary tag.
    Returns default_val if tag not found or wrong type.
    """
    cdef uint8_t* tag_data = find_tag(record, tag)
    if tag_data == NULL:
        return default_val

    cdef uint8_t tag_type = tag_data[0]
    tag_data += 1  # Skip type byte

    if tag_type == b'c':
        return <int32_t>(<int8_t*>tag_data)[0]
    elif tag_type == b'C':
        return <int32_t>(<uint8_t*>tag_data)[0]
    elif tag_type == b's':
        return <int32_t>(<int16_t*>tag_data)[0]
    elif tag_type == b'S':
        return <int32_t>(<uint16_t*>tag_data)[0]
    elif tag_type == b'i':
        return (<int32_t*>tag_data)[0]
    elif tag_type == b'I':
        return <int32_t>(<uint32_t*>tag_data)[0]

    return default_val


cdef float get_tag_float(bam1_t* record, const char* tag, float default_val) nogil:
    """
    Get float value from auxiliary tag.
    Returns default_val if tag not found or wrong type.
    """
    cdef uint8_t* tag_data = find_tag(record, tag)
    if tag_data == NULL:
        return default_val

    cdef uint8_t tag_type = tag_data[0]
    tag_data += 1  # Skip type byte

    if tag_type == b'f':
        return (<float*>tag_data)[0]
    elif tag_type == b'i':
        return <float>(<int32_t*>tag_data)[0]
    elif tag_type == b'I':
        return <float>(<uint32_t*>tag_data)[0]

    return default_val


# C function for snprintf
cdef extern from "stdio.h":
    int snprintf(char* s, size_t n, const char* format, ...) nogil
