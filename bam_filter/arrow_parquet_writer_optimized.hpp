/*
 * Optimized Arrow Parquet Writer - Storage-efficient schema for SAM/BAM to Parquet
 *
 * Key optimizations:
 * 1. No data duplication (hot tags as columns, cold tags in string)
 * 2. 2-bit packed sequences (4x compression)
 * 3. Raw binary quality scores
 * 4. ref_id dictionary instead of ref_name strings
 * 5. Removed computed columns (end_position, alignment_length)
 * 6. Support for normalized output (separate reads/alignments tables)
 */

#ifndef ARROW_PARQUET_WRITER_OPTIMIZED_HPP
#define ARROW_PARQUET_WRITER_OPTIMIZED_HPP

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <memory>
#include <string>
#include <vector>
#include <unordered_map>

namespace bam_filter {

// LUT for base to 2-bit encoding (A=00, C=01, G=10, T=11, N=00)
// Pre-computed for all 256 byte values - branchless lookup
alignas(64) static constexpr uint8_t BASE_TO_2BIT[256] = {
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 0-15
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 16-31
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 32-47
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 48-63
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0, // 64-79:  A=65->0, C=67->1, G=71->2
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0, // 80-95:  T=84->3
    0,0,0,1,0,0,0,2,0,0,0,0,0,0,0,0, // 96-111: a=97->0, c=99->1, g=103->2
    0,0,0,0,3,0,0,0,0,0,0,0,0,0,0,0, // 112-127: t=116->3
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 128-143
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 144-159
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 160-175
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 176-191
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 192-207
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 208-223
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 224-239
    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, // 240-255
};

// Pack DNA sequence to 2-bit encoding using LUT (branchless, cache-friendly)
inline std::string pack_sequence_2bit(const char* seq, size_t len) {
    std::string packed;
    const size_t packed_len = (len + 3) / 4;
    packed.resize(packed_len);

    // Process 4 bases at a time (unrolled, branchless)
    size_t full_bytes = len / 4;
    const uint8_t* s = reinterpret_cast<const uint8_t*>(seq);

    for (size_t i = 0; i < full_bytes; ++i) {
        const size_t j = i * 4;
        packed[i] = static_cast<char>(
            (BASE_TO_2BIT[s[j]]     << 6) |
            (BASE_TO_2BIT[s[j + 1]] << 4) |
            (BASE_TO_2BIT[s[j + 2]] << 2) |
            (BASE_TO_2BIT[s[j + 3]])
        );
    }

    // Handle remaining 1-3 bases
    size_t remaining = len % 4;
    if (remaining > 0) {
        uint8_t byte = 0;
        size_t j = full_bytes * 4;
        for (size_t k = 0; k < remaining; ++k) {
            byte |= (BASE_TO_2BIT[s[j + k]] << (6 - 2 * k));
        }
        packed[full_bytes] = static_cast<char>(byte);
    }

    return packed;
}

// Unpack 2-bit sequence back to ASCII
inline std::string unpack_sequence_2bit(const std::string& packed, size_t original_len) {
    static const char bases[] = "ACGT";
    std::string seq;
    seq.reserve(original_len);

    for (size_t i = 0; i < packed.size() && seq.size() < original_len; ++i) {
        uint8_t byte = static_cast<uint8_t>(packed[i]);
        for (int j = 0; j < 4 && seq.size() < original_len; ++j) {
            uint8_t bits = (byte >> (6 - 2 * j)) & 0x03;
            seq.push_back(bases[bits]);
        }
    }
    return seq;
}

// Optimized alignment batch - no duplication, efficient storage
struct OptimizedAlignmentBatch {
    // Core fields
    std::vector<uint64_t> read_ids;
    std::vector<std::string> read_names;
    std::vector<uint32_t> ref_ids;           // Dictionary ID, not string
    std::vector<int32_t> positions;
    std::vector<uint8_t> mapqs;
    std::vector<uint16_t> flags;
    std::vector<std::string> cigars;
    std::vector<int32_t> template_lengths;
    std::vector<int32_t> mate_ref_ids;
    std::vector<int32_t> mate_positions;

    // Sequences - 2-bit packed binary
    std::vector<std::string> sequences_packed;  // 2-bit packed
    std::vector<uint16_t> sequence_lengths;     // Original length for unpacking
    std::vector<std::string> qualities;         // Raw bytes (not ASCII)

    // Hot tags as columns (NOT in tags_cold)
    std::vector<int16_t> AS;                    // Alignment score
    std::vector<uint16_t> NM;                   // Edit distance
    std::vector<int16_t> XS;                    // Secondary score
    std::vector<std::string> MD;                // Mismatch string

    // Cold tags (rare tags only, hot tags stripped)
    std::vector<std::string> tags_cold;

    // Derived fields (precomputed for fast queries)
    std::vector<uint16_t> aligned_lengths; // reference span from CIGAR (uint16 compresses better than end_position)
    std::vector<bool> is_reverse;          // flag & 0x10 (for strand filtering)

    // Pipeline results (written by downstream stages)
    std::vector<float> ani;
    std::vector<float> zs_score;          // Log-likelihood alignment score (ZS:f tag)
    std::vector<float> pmd_score;         // Post-mortem damage score (PM:f tag)
    std::vector<bool> filter_passed;
    std::vector<int32_t> lca_taxid;
    std::vector<int32_t> reassigned_ref_id;
    std::vector<float> zp_posterior;

    // Read group (optional)
    std::vector<std::string> read_groups;

    size_t size() const { return read_ids.size(); }

    void clear() {
        read_ids.clear();
        read_names.clear();
        ref_ids.clear();
        positions.clear();
        mapqs.clear();
        flags.clear();
        cigars.clear();
        template_lengths.clear();
        mate_ref_ids.clear();
        mate_positions.clear();
        sequences_packed.clear();
        sequence_lengths.clear();
        qualities.clear();
        AS.clear();
        NM.clear();
        XS.clear();
        MD.clear();
        tags_cold.clear();
        aligned_lengths.clear();
        is_reverse.clear();
        ani.clear();
        zs_score.clear();
        pmd_score.clear();
        filter_passed.clear();
        lca_taxid.clear();
        reassigned_ref_id.clear();
        zp_posterior.clear();
        read_groups.clear();
    }

    void reserve(size_t capacity) {
        read_ids.reserve(capacity);
        read_names.reserve(capacity);
        ref_ids.reserve(capacity);
        positions.reserve(capacity);
        mapqs.reserve(capacity);
        flags.reserve(capacity);
        cigars.reserve(capacity);
        template_lengths.reserve(capacity);
        mate_ref_ids.reserve(capacity);
        mate_positions.reserve(capacity);
        sequences_packed.reserve(capacity);
        sequence_lengths.reserve(capacity);
        qualities.reserve(capacity);
        AS.reserve(capacity);
        NM.reserve(capacity);
        XS.reserve(capacity);
        MD.reserve(capacity);
        tags_cold.reserve(capacity);
        aligned_lengths.reserve(capacity);
        is_reverse.reserve(capacity);
        ani.reserve(capacity);
        zs_score.reserve(capacity);
        pmd_score.reserve(capacity);
        filter_passed.reserve(capacity);
        lca_taxid.reserve(capacity);
        reassigned_ref_id.reserve(capacity);
        zp_posterior.reserve(capacity);
        read_groups.reserve(capacity);
    }
};

// Read batch for normalized schema (one row per unique read)
struct ReadBatch {
    std::vector<uint64_t> read_ids;
    std::vector<std::string> read_names;
    std::vector<std::string> sequences_packed;  // 2-bit packed
    std::vector<uint16_t> sequence_lengths;
    std::vector<std::string> qualities;         // Raw bytes
    std::vector<std::string> read_groups;

    size_t size() const { return read_ids.size(); }

    void clear() {
        read_ids.clear();
        read_names.clear();
        sequences_packed.clear();
        sequence_lengths.clear();
        qualities.clear();
        read_groups.clear();
    }

    void reserve(size_t capacity) {
        read_ids.reserve(capacity);
        read_names.reserve(capacity);
        sequences_packed.reserve(capacity);
        sequence_lengths.reserve(capacity);
        qualities.reserve(capacity);
        read_groups.reserve(capacity);
    }
};

// Alignment batch for normalized schema (no sequence/quality)
struct AlignmentOnlyBatch {
    std::vector<uint64_t> read_ids;             // FK to reads table
    std::vector<uint32_t> ref_ids;
    std::vector<int32_t> positions;
    std::vector<uint8_t> mapqs;
    std::vector<uint16_t> flags;
    std::vector<std::string> cigars;
    std::vector<int32_t> template_lengths;
    std::vector<int32_t> mate_ref_ids;
    std::vector<int32_t> mate_positions;

    // Hot tags
    std::vector<int16_t> AS;
    std::vector<uint16_t> NM;
    std::vector<int16_t> XS;
    std::vector<std::string> MD;
    std::vector<std::string> tags_cold;

    // Derived fields
    std::vector<uint16_t> aligned_lengths;
    std::vector<bool> is_reverse;

    // Pipeline results
    std::vector<float> ani;
    std::vector<float> zs_score;          // Log-likelihood alignment score (ZS:f tag)
    std::vector<float> pmd_score;         // Post-mortem damage score (PM:f tag)
    std::vector<bool> filter_passed;
    std::vector<int32_t> lca_taxid;
    std::vector<int32_t> reassigned_ref_id;
    std::vector<float> zp_posterior;

    size_t size() const { return read_ids.size(); }

    void clear() {
        read_ids.clear();
        ref_ids.clear();
        positions.clear();
        mapqs.clear();
        flags.clear();
        cigars.clear();
        template_lengths.clear();
        mate_ref_ids.clear();
        mate_positions.clear();
        AS.clear();
        NM.clear();
        XS.clear();
        MD.clear();
        tags_cold.clear();
        aligned_lengths.clear();
        is_reverse.clear();
        ani.clear();
        zs_score.clear();
        pmd_score.clear();
        filter_passed.clear();
        lca_taxid.clear();
        reassigned_ref_id.clear();
        zp_posterior.clear();
    }

    void reserve(size_t capacity) {
        read_ids.reserve(capacity);
        ref_ids.reserve(capacity);
        positions.reserve(capacity);
        mapqs.reserve(capacity);
        flags.reserve(capacity);
        cigars.reserve(capacity);
        template_lengths.reserve(capacity);
        mate_ref_ids.reserve(capacity);
        mate_positions.reserve(capacity);
        AS.reserve(capacity);
        NM.reserve(capacity);
        XS.reserve(capacity);
        MD.reserve(capacity);
        tags_cold.reserve(capacity);
        aligned_lengths.reserve(capacity);
        is_reverse.reserve(capacity);
        ani.reserve(capacity);
        zs_score.reserve(capacity);
        pmd_score.reserve(capacity);
        filter_passed.reserve(capacity);
        lca_taxid.reserve(capacity);
        reassigned_ref_id.reserve(capacity);
        zp_posterior.reserve(capacity);
    }
};

// Reference info for sidecar table
struct ReferenceInfo {
    uint32_t ref_id;
    std::string ref_name;
    int64_t ref_length;
    int32_t taxid;
};

// Optimized Parquet writer - single table mode
class OptimizedParquetWriter {
public:
    OptimizedParquetWriter(const std::string& filename, int compression_level = 6);
    ~OptimizedParquetWriter();

    arrow::Status WriteBatch(const OptimizedAlignmentBatch& batch);
    arrow::Status Close();

private:
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    bool closed_;

    static std::shared_ptr<arrow::Schema> CreateSchema();
    arrow::Result<std::shared_ptr<arrow::Table>> BatchToTable(const OptimizedAlignmentBatch& batch);
};

// Normalized Parquet writer - separate reads/alignments tables
class NormalizedParquetWriter {
public:
    NormalizedParquetWriter(const std::string& output_dir, int compression_level = 6);
    ~NormalizedParquetWriter();

    // Write reads batch
    arrow::Status WriteReadBatch(const ReadBatch& batch);

    // Write alignments batch
    arrow::Status WriteAlignmentBatch(const AlignmentOnlyBatch& batch);

    // Write reference sidecar
    arrow::Status WriteReferences(const std::vector<ReferenceInfo>& refs);

    arrow::Status Close();

private:
    std::string output_dir_;
    int compression_level_;

    std::shared_ptr<arrow::Schema> reads_schema_;
    std::shared_ptr<arrow::Schema> alignments_schema_;
    std::shared_ptr<arrow::Schema> references_schema_;

    std::shared_ptr<arrow::io::FileOutputStream> reads_outfile_;
    std::shared_ptr<arrow::io::FileOutputStream> alignments_outfile_;

    std::unique_ptr<parquet::arrow::FileWriter> reads_writer_;
    std::unique_ptr<parquet::arrow::FileWriter> alignments_writer_;

    bool closed_;

    static std::shared_ptr<arrow::Schema> CreateReadsSchema();
    static std::shared_ptr<arrow::Schema> CreateAlignmentsSchema();
    static std::shared_ptr<arrow::Schema> CreateReferencesSchema();

    arrow::Result<std::shared_ptr<arrow::Table>> ReadBatchToTable(const ReadBatch& batch);
    arrow::Result<std::shared_ptr<arrow::Table>> AlignmentBatchToTable(const AlignmentOnlyBatch& batch);
};

} // namespace bam_filter

#endif // ARROW_PARQUET_WRITER_OPTIMIZED_HPP
