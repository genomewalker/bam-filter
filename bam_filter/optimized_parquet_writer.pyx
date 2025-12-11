# cython: language_level=3
# cython: boundscheck=False
# cython: wraparound=False
# cython: cdivision=True
# distutils: language = c++
"""
Optimized Parquet Writer - Cython bindings for storage-efficient SAM/BAM to Parquet.

Key features:
- 2-bit sequence packing (4x compression)
- Hot/cold tag separation (no duplication)
- Normalized output (separate reads/alignments tables)
- Zero-copy Python string access (no intermediate bytes objects)
"""

from libc.stdint cimport uint8_t, uint16_t, uint32_t, uint64_t, int16_t, int32_t, int64_t
from libc.string cimport memcpy, strlen
from libcpp.string cimport string
from libcpp.vector cimport vector
from libcpp cimport bool as cpp_bool
from cpython.ref cimport PyObject

import os


# ============================================================================
# Zero-copy Python string access (avoids creating intermediate bytes objects)
# ============================================================================

cdef extern from "Python.h":
    const char* PyUnicode_AsUTF8AndSize(object o, Py_ssize_t *size) except NULL

cdef inline string str_to_cpp_string(str py_str) noexcept:
    """Zero-copy conversion from Python str to C++ string.
    Uses PyUnicode_AsUTF8AndSize to get pointer directly without creating bytes object.
    """
    cdef:
        const char* ptr
        Py_ssize_t length
    ptr = PyUnicode_AsUTF8AndSize(py_str, &length)
    return string(ptr, length)


# ============================================================================
# C++ declarations
# ============================================================================

cdef extern from "arrow_parquet_writer_optimized.hpp" namespace "bam_filter":
    # Sequence packing functions
    string pack_sequence_2bit(const char* seq, size_t length) nogil
    string unpack_sequence_2bit(const string& packed, size_t original_len) nogil

    # Optimized batch structure
    cdef cppclass OptimizedAlignmentBatch:
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
        vector[string] sequences_packed
        vector[uint16_t] sequence_lengths
        vector[string] qualities
        vector[int16_t] AS
        vector[uint16_t] NM
        vector[int16_t] XS
        vector[string] MD
        vector[string] tags_cold
        vector[float] ani
        vector[float] pmd_score
        vector[cpp_bool] filter_passed
        vector[int32_t] lca_taxid
        vector[int32_t] reassigned_ref_id
        vector[float] zp_posterior
        vector[string] read_groups

        size_t size()
        void clear()
        void reserve(size_t capacity)

    # Read batch for normalized schema
    cdef cppclass ReadBatch:
        vector[uint64_t] read_ids
        vector[string] read_names
        vector[string] sequences_packed
        vector[uint16_t] sequence_lengths
        vector[string] qualities
        vector[string] read_groups

        size_t size()
        void clear()
        void reserve(size_t capacity)

    # Alignment batch for normalized schema
    cdef cppclass AlignmentOnlyBatch:
        vector[uint64_t] read_ids
        vector[uint32_t] ref_ids
        vector[int32_t] positions
        vector[uint8_t] mapqs
        vector[uint16_t] flags
        vector[string] cigars
        vector[int32_t] template_lengths
        vector[int32_t] mate_ref_ids
        vector[int32_t] mate_positions
        vector[int16_t] AS
        vector[uint16_t] NM
        vector[int16_t] XS
        vector[string] MD
        vector[string] tags_cold
        vector[float] ani
        vector[float] pmd_score
        vector[cpp_bool] filter_passed
        vector[int32_t] lca_taxid
        vector[int32_t] reassigned_ref_id
        vector[float] zp_posterior

        size_t size()
        void clear()
        void reserve(size_t capacity)

    # Reference info
    cdef cppclass ReferenceInfo:
        uint32_t ref_id
        string ref_name
        int64_t ref_length
        int32_t taxid

    # Optimized writer (single table)
    cdef cppclass OptimizedParquetWriter:
        OptimizedParquetWriter(const string& filename, int compression_level) except +
        void WriteBatch(const OptimizedAlignmentBatch& batch) except +
        void Close() except +

    # Normalized writer (separate tables)
    cdef cppclass NormalizedParquetWriter:
        NormalizedParquetWriter(const string& output_dir, int compression_level) except +
        void WriteReadBatch(const ReadBatch& batch) except +
        void WriteAlignmentBatch(const AlignmentOnlyBatch& batch) except +
        void WriteReferences(const vector[ReferenceInfo]& refs) except +
        void Close() except +


# ============================================================================
# Hot tags to extract (stripped from tags_cold)
# ============================================================================

HOT_TAGS = {b'AS', b'NM', b'XS', b'MD', b'XN', b'XM', b'XO', b'XG', b'YT'}


cdef inline void parse_tags_hot_cold(
    const char* tags_raw, int tags_len,
    int16_t* AS_out, uint16_t* NM_out, int16_t* XS_out,
    string* MD_out, string* tags_cold_out
) noexcept nogil:
    """
    Parse tags string into hot tags (columns) and cold tags (remaining string).
    Hot tags are removed from the cold string to avoid duplication.
    """
    cdef:
        const char* p = tags_raw
        const char* end = tags_raw + tags_len
        const char* tag_start
        const char* tag_end
        char tag[3]
        char tag_type
        int int_val
        float float_val
        string cold_parts
        # ASCII constants (for nogil compatibility)
        char TAB = 9      # '\t'
        char NL = 10      # '\n'
        char COLON = 58   # ':'
        char MINUS = 45   # '-'
        char ZERO = 48    # '0'
        char NINE = 57    # '9'
        char A = 65
        char D = 68
        char G = 71
        char M = 77
        char N = 78
        char O = 79
        char S = 83
        char T = 84
        char X = 88
        char Y = 89
        char Z = 90
        char i_char = 105  # 'i'

    # Initialize outputs
    AS_out[0] = -32768  # Missing value sentinel
    NM_out[0] = 65535
    XS_out[0] = -32768
    MD_out[0] = string()
    tags_cold_out[0] = string()

    while p < end:
        # Skip leading tab
        if p[0] == TAB:
            p += 1
            continue

        tag_start = p

        # Read tag name (2 chars)
        if p + 2 >= end:
            break
        tag[0] = p[0]
        tag[1] = p[1]
        tag[2] = 0
        p += 2

        # Read colon and type
        if p >= end or p[0] != COLON:
            break
        p += 1
        if p >= end:
            break
        tag_type = p[0]
        p += 1
        if p >= end or p[0] != COLON:
            break
        p += 1

        # Find end of tag value
        tag_end = p
        while tag_end < end and tag_end[0] != TAB and tag_end[0] != NL:
            tag_end += 1

        # Check if this is a hot tag
        if tag[0] == A and tag[1] == S and tag_type == i_char:
            # Parse AS:i:value
            int_val = 0
            if p < tag_end:
                if p[0] == MINUS:
                    p += 1
                    while p < tag_end and p[0] >= ZERO and p[0] <= NINE:
                        int_val = int_val * 10 + (p[0] - ZERO)
                        p += 1
                    int_val = -int_val
                else:
                    while p < tag_end and p[0] >= ZERO and p[0] <= NINE:
                        int_val = int_val * 10 + (p[0] - ZERO)
                        p += 1
            AS_out[0] = <int16_t>int_val

        elif tag[0] == N and tag[1] == M and tag_type == i_char:
            # Parse NM:i:value
            int_val = 0
            while p < tag_end and p[0] >= ZERO and p[0] <= NINE:
                int_val = int_val * 10 + (p[0] - ZERO)
                p += 1
            NM_out[0] = <uint16_t>int_val

        elif tag[0] == X and tag[1] == S and tag_type == i_char:
            # Parse XS:i:value
            int_val = 0
            if p < tag_end:
                if p[0] == MINUS:
                    p += 1
                    while p < tag_end and p[0] >= ZERO and p[0] <= NINE:
                        int_val = int_val * 10 + (p[0] - ZERO)
                        p += 1
                    int_val = -int_val
                else:
                    while p < tag_end and p[0] >= ZERO and p[0] <= NINE:
                        int_val = int_val * 10 + (p[0] - ZERO)
                        p += 1
            XS_out[0] = <int16_t>int_val

        elif tag[0] == M and tag[1] == D and tag_type == Z:
            # Store MD string
            MD_out[0] = string(p, tag_end - p)

        elif (tag[0] == X and (tag[1] == N or tag[1] == M or tag[1] == O or tag[1] == G)) or \
             (tag[0] == Y and tag[1] == T):
            # Skip other hot tags (XN, XM, XO, XG, YT) - don't add to cold
            pass

        else:
            # Cold tag - add to output
            if cold_parts.size() > 0:
                cold_parts.push_back(TAB)
            cold_parts.append(tag_start, tag_end - tag_start)

        p = tag_end

    tags_cold_out[0] = cold_parts


# ============================================================================
# Python wrapper for optimized single-table writer
# ============================================================================

cdef class PyOptimizedParquetWriter:
    """
    Python wrapper for optimized Parquet writer (single table mode).

    Usage:
        writer = PyOptimizedParquetWriter("/path/to/output.parquet", compression_level=6)
        for batch in parse_sam(...):
            writer.write_batch(batch)
        writer.close()
    """
    cdef OptimizedParquetWriter* writer
    cdef OptimizedAlignmentBatch batch
    cdef dict ref_name_to_id
    cdef list ref_names
    cdef uint64_t next_read_id

    def __cinit__(self, str filename, int compression_level=6, int batch_size=50000):
        cdef string cpp_filename = str_to_cpp_string(filename)
        self.writer = new OptimizedParquetWriter(cpp_filename, compression_level)
        self.ref_name_to_id = {}
        self.ref_names = []
        self.next_read_id = 0
        self._batch_size = batch_size
        self._batch_reserved = False

    cdef int _batch_size
    cdef bint _batch_reserved

    def __dealloc__(self):
        if self.writer != NULL:
            del self.writer

    cdef inline void _ensure_batch_reserved(self) noexcept:
        """Pre-allocate batch memory to avoid reallocations during append."""
        if not self._batch_reserved:
            self.batch.reserve(self._batch_size)
            self._batch_reserved = True

    def get_ref_id(self, str ref_name):
        """Get or create ref_id for a reference name."""
        if ref_name not in self.ref_name_to_id:
            ref_id = len(self.ref_names)
            self.ref_name_to_id[ref_name] = ref_id
            self.ref_names.append(ref_name)
        return self.ref_name_to_id[ref_name]

    def add_alignment(
        self,
        str read_name,
        str ref_name,
        int position,
        int mapq,
        int flag,
        str cigar,
        str sequence,
        str quality,
        str tags_raw,
        int template_length=0,
        str mate_ref="*",
        int mate_pos=-1,
    ):
        """Add a single alignment to the batch (zero-copy string handling)."""
        cdef:
            # Zero-copy string conversion - no intermediate Python bytes objects
            string cpp_read_name = str_to_cpp_string(read_name)
            string cpp_cigar = str_to_cpp_string(cigar)
            string cpp_sequence = str_to_cpp_string(sequence)
            string cpp_quality = str_to_cpp_string(quality)
            string cpp_tags_raw = str_to_cpp_string(tags_raw)
            string packed_seq
            int16_t AS_val
            uint16_t NM_val
            int16_t XS_val
            string MD_val
            string tags_cold_val
            uint32_t ref_id
            Py_ssize_t seq_len = len(sequence)

        # Pre-reserve batch storage
        self._ensure_batch_reserved()

        # Get ref_id
        ref_id = self.get_ref_id(ref_name)

        # Pack sequence (uses LUT-based packing in C++)
        packed_seq = pack_sequence_2bit(cpp_sequence.c_str(), cpp_sequence.size())

        # Parse hot/cold tags (nogil)
        parse_tags_hot_cold(
            cpp_tags_raw.c_str(), cpp_tags_raw.size(),
            &AS_val, &NM_val, &XS_val, &MD_val, &tags_cold_val
        )

        # Add to batch - all push_back operations
        self.batch.read_ids.push_back(self.next_read_id)
        self.batch.read_names.push_back(cpp_read_name)
        self.batch.ref_ids.push_back(ref_id)
        self.batch.positions.push_back(position)
        self.batch.mapqs.push_back(<uint8_t>mapq)
        self.batch.flags.push_back(<uint16_t>flag)
        self.batch.cigars.push_back(cpp_cigar)
        self.batch.template_lengths.push_back(template_length)
        self.batch.mate_ref_ids.push_back(-1)  # TODO: parse mate ref
        self.batch.mate_positions.push_back(mate_pos)
        self.batch.sequences_packed.push_back(packed_seq)
        self.batch.sequence_lengths.push_back(<uint16_t>seq_len)
        self.batch.qualities.push_back(cpp_quality)
        self.batch.AS.push_back(AS_val)
        self.batch.NM.push_back(NM_val)
        self.batch.XS.push_back(XS_val)
        self.batch.MD.push_back(MD_val)
        self.batch.tags_cold.push_back(tags_cold_val)

        # Initialize result columns
        self.batch.ani.push_back(0.0)
        self.batch.pmd_score.push_back(0.0)
        self.batch.filter_passed.push_back(False)
        self.batch.lca_taxid.push_back(-1)
        self.batch.reassigned_ref_id.push_back(-1)
        self.batch.zp_posterior.push_back(0.0)
        self.batch.read_groups.push_back(string())

        self.next_read_id += 1

    def flush_batch(self):
        """Flush current batch to Parquet."""
        if self.batch.size() > 0:
            self.writer.WriteBatch(self.batch)
            self.batch.clear()
            self._batch_reserved = False  # Reset for next batch

    def close(self):
        """Close the writer."""
        self.flush_batch()
        self.writer.Close()

    @property
    def references(self):
        """Return list of reference names in order of ref_id."""
        return self.ref_names.copy()


# ============================================================================
# Python wrapper for normalized two-table writer
# ============================================================================

cdef class PyNormalizedParquetWriter:
    """
    Python wrapper for normalized Parquet writer (separate reads/alignments tables).

    This is the most storage-efficient format for multi-mapped reads.

    Usage:
        writer = PyNormalizedParquetWriter("/path/to/output_dir", compression_level=6)

        for read_name, alignments in group_by_read(parse_sam(...)):
            # Write read once
            writer.add_read(read_name, sequence, quality)

            # Write each alignment
            for aln in alignments:
                writer.add_alignment(read_name, ref_name, position, ...)

        writer.close()
    """
    cdef NormalizedParquetWriter* writer
    cdef ReadBatch read_batch
    cdef AlignmentOnlyBatch alignment_batch
    cdef dict read_name_to_id
    cdef dict ref_name_to_id
    cdef list ref_names
    cdef uint64_t next_read_id
    cdef int batch_size

    def __cinit__(self, str output_dir, int compression_level=6, int batch_size=50000):
        cdef string cpp_dir = str_to_cpp_string(output_dir)
        self.writer = new NormalizedParquetWriter(cpp_dir, compression_level)
        self.read_name_to_id = {}
        self.ref_name_to_id = {}
        self.ref_names = []
        self.next_read_id = 0
        self.batch_size = batch_size
        self._reads_reserved = False
        self._alignments_reserved = False

    cdef bint _reads_reserved
    cdef bint _alignments_reserved

    def __dealloc__(self):
        if self.writer != NULL:
            del self.writer

    cdef inline void _ensure_reads_reserved(self) noexcept:
        """Pre-allocate read batch memory."""
        if not self._reads_reserved:
            self.read_batch.reserve(self.batch_size)
            self._reads_reserved = True

    cdef inline void _ensure_alignments_reserved(self) noexcept:
        """Pre-allocate alignment batch memory."""
        if not self._alignments_reserved:
            self.alignment_batch.reserve(self.batch_size)
            self._alignments_reserved = True

    def get_ref_id(self, str ref_name):
        """Get or create ref_id for a reference name."""
        if ref_name not in self.ref_name_to_id:
            ref_id = len(self.ref_names)
            self.ref_name_to_id[ref_name] = ref_id
            self.ref_names.append(ref_name)
        return self.ref_name_to_id[ref_name]

    def get_read_id(self, str read_name):
        """Get or create read_id for a read name."""
        if read_name not in self.read_name_to_id:
            read_id = self.next_read_id
            self.read_name_to_id[read_name] = read_id
            self.next_read_id += 1
            return read_id, True  # New read
        return self.read_name_to_id[read_name], False  # Existing read

    def add_read(self, str read_name, str sequence, str quality, str read_group=""):
        """Add a unique read to the reads table (zero-copy)."""
        cdef:
            # Zero-copy string conversion
            string cpp_read_name = str_to_cpp_string(read_name)
            string cpp_sequence = str_to_cpp_string(sequence)
            string cpp_quality = str_to_cpp_string(quality)
            string cpp_read_group = str_to_cpp_string(read_group)
            string packed_seq
            uint64_t read_id
            Py_ssize_t seq_len = len(sequence)

        read_id, is_new = self.get_read_id(read_name)
        if not is_new:
            return  # Already added

        self._ensure_reads_reserved()
        packed_seq = pack_sequence_2bit(cpp_sequence.c_str(), cpp_sequence.size())

        self.read_batch.read_ids.push_back(read_id)
        self.read_batch.read_names.push_back(cpp_read_name)
        self.read_batch.sequences_packed.push_back(packed_seq)
        self.read_batch.sequence_lengths.push_back(<uint16_t>seq_len)
        self.read_batch.qualities.push_back(cpp_quality)
        self.read_batch.read_groups.push_back(cpp_read_group)

        if self.read_batch.size() >= self.batch_size:
            self.flush_reads()

    def add_alignment(
        self,
        str read_name,
        str ref_name,
        int position,
        int mapq,
        int flag,
        str cigar,
        str tags_raw,
        int template_length=0,
        int mate_ref_id=-1,
        int mate_pos=-1,
    ):
        """Add an alignment to the alignments table (zero-copy)."""
        cdef:
            # Zero-copy string conversion
            string cpp_cigar = str_to_cpp_string(cigar)
            string cpp_tags_raw = str_to_cpp_string(tags_raw)
            int16_t AS_val
            uint16_t NM_val
            int16_t XS_val
            string MD_val
            string tags_cold_val
            uint64_t read_id
            uint32_t ref_id

        self._ensure_alignments_reserved()

        # Get IDs
        read_id = self.read_name_to_id.get(read_name, self.next_read_id)
        ref_id = self.get_ref_id(ref_name)

        # Parse hot/cold tags (nogil)
        parse_tags_hot_cold(
            cpp_tags_raw.c_str(), cpp_tags_raw.size(),
            &AS_val, &NM_val, &XS_val, &MD_val, &tags_cold_val
        )

        # Add to batch
        self.alignment_batch.read_ids.push_back(read_id)
        self.alignment_batch.ref_ids.push_back(ref_id)
        self.alignment_batch.positions.push_back(position)
        self.alignment_batch.mapqs.push_back(<uint8_t>mapq)
        self.alignment_batch.flags.push_back(<uint16_t>flag)
        self.alignment_batch.cigars.push_back(cpp_cigar)
        self.alignment_batch.template_lengths.push_back(template_length)
        self.alignment_batch.mate_ref_ids.push_back(mate_ref_id)
        self.alignment_batch.mate_positions.push_back(mate_pos)
        self.alignment_batch.AS.push_back(AS_val)
        self.alignment_batch.NM.push_back(NM_val)
        self.alignment_batch.XS.push_back(XS_val)
        self.alignment_batch.MD.push_back(MD_val)
        self.alignment_batch.tags_cold.push_back(tags_cold_val)

        # Initialize result columns
        self.alignment_batch.ani.push_back(0.0)
        self.alignment_batch.pmd_score.push_back(0.0)
        self.alignment_batch.filter_passed.push_back(False)
        self.alignment_batch.lca_taxid.push_back(-1)
        self.alignment_batch.reassigned_ref_id.push_back(-1)
        self.alignment_batch.zp_posterior.push_back(0.0)

        if self.alignment_batch.size() >= self.batch_size:
            self.flush_alignments()

    def flush_reads(self):
        """Flush reads batch to Parquet."""
        if self.read_batch.size() > 0:
            self.writer.WriteReadBatch(self.read_batch)
            self.read_batch.clear()
            self._reads_reserved = False

    def flush_alignments(self):
        """Flush alignments batch to Parquet."""
        if self.alignment_batch.size() > 0:
            self.writer.WriteAlignmentBatch(self.alignment_batch)
            self.alignment_batch.clear()
            self._alignments_reserved = False

    def write_references(self):
        """Write reference sidecar file."""
        cdef vector[ReferenceInfo] refs
        cdef ReferenceInfo ref_info

        for i, ref_name in enumerate(self.ref_names):
            ref_info.ref_id = i
            ref_info.ref_name = str_to_cpp_string(ref_name)
            ref_info.ref_length = 0  # Unknown
            ref_info.taxid = -1  # Unknown
            refs.push_back(ref_info)

        self.writer.WriteReferences(refs)

    def close(self):
        """Close the writer and write sidecar files."""
        self.flush_reads()
        self.flush_alignments()
        self.write_references()
        self.writer.Close()

    @property
    def references(self):
        """Return list of reference names in order of ref_id."""
        return self.ref_names.copy()

    @property
    def num_reads(self):
        """Return number of unique reads."""
        return self.next_read_id

    @property
    def num_alignments(self):
        """Return number of alignments written."""
        return len(self.read_name_to_id)


# ============================================================================
# Convenience functions
# ============================================================================

def pack_sequence(str sequence):
    """Pack a DNA sequence to 2-bit encoding (zero-copy)."""
    cdef string cpp_seq = str_to_cpp_string(sequence)
    cdef string packed = pack_sequence_2bit(cpp_seq.c_str(), cpp_seq.size())
    return bytes(packed)


def unpack_sequence(bytes packed, int original_length):
    """Unpack a 2-bit encoded sequence back to ASCII."""
    cdef string cpp_packed = packed
    cdef string unpacked = unpack_sequence_2bit(cpp_packed, original_length)
    return unpacked.decode('utf-8')
