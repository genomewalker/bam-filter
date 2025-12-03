/*
 * Arrow Parquet Writer - Minimal Schema Version
 *
 * Optimized for maximum speed with only essential SAM fields:
 * - 8 core SAM fields: qname, flag, rname, pos, mapq, cigar, seq, qual
 * - 8 aux tags: AS, XN, XM, XO, XG, NM, MD, YT
 *
 * Total: 16 fields (down from 30)
 */

#ifndef ARROW_PARQUET_WRITER_MINIMAL_HPP
#define ARROW_PARQUET_WRITER_MINIMAL_HPP

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <memory>
#include <string>
#include <vector>

namespace bam_filter {

// Minimal columnar batch - only essential fields
struct MinimalAlignmentBatch {
    // Core SAM fields
    std::vector<std::string> read_names;   // qname
    std::vector<uint16_t> flags;            // flag
    std::vector<std::string> ref_names;     // rname
    std::vector<int32_t> positions;         // pos
    std::vector<uint8_t> mapqs;             // mapq
    std::vector<std::string> cigars;        // cigar
    std::vector<std::string> sequences;     // seq
    std::vector<std::string> qualities;     // qual

    // Standard aux tags (integers)
    std::vector<int32_t> tag_AS;  // Alignment score
    std::vector<int32_t> tag_XN;  // Number of ambiguous bases
    std::vector<int32_t> tag_XM;  // Number of mismatches
    std::vector<int32_t> tag_XO;  // Number of gap opens
    std::vector<int32_t> tag_XG;  // Number of gap extensions
    std::vector<int32_t> tag_NM;  // Edit distance

    // String aux tags
    std::vector<std::string> tag_MD;  // Mismatch string
    std::vector<std::string> tag_YT;  // Alignment type (UU, CP, DP, UP)

    size_t size() const { return read_names.size(); }

    void clear() {
        read_names.clear();
        flags.clear();
        ref_names.clear();
        positions.clear();
        mapqs.clear();
        cigars.clear();
        sequences.clear();
        qualities.clear();
        tag_AS.clear();
        tag_XN.clear();
        tag_XM.clear();
        tag_XO.clear();
        tag_XG.clear();
        tag_NM.clear();
        tag_MD.clear();
        tag_YT.clear();
    }

    void reserve(size_t capacity) {
        read_names.reserve(capacity);
        flags.reserve(capacity);
        ref_names.reserve(capacity);
        positions.reserve(capacity);
        mapqs.reserve(capacity);
        cigars.reserve(capacity);
        sequences.reserve(capacity);
        qualities.reserve(capacity);
        tag_AS.reserve(capacity);
        tag_XN.reserve(capacity);
        tag_XM.reserve(capacity);
        tag_XO.reserve(capacity);
        tag_XG.reserve(capacity);
        tag_NM.reserve(capacity);
        tag_MD.reserve(capacity);
        tag_YT.reserve(capacity);
    }
};

// Minimal Parquet writer
class MinimalParquetWriter {
public:
    MinimalParquetWriter(const std::string& filename, int compression_level = 3);
    ~MinimalParquetWriter();

    // Write a batch of alignments
    arrow::Status WriteBatch(const MinimalAlignmentBatch& batch);

    // Close the file
    arrow::Status Close();

private:
    std::shared_ptr<arrow::Schema> schema_;
    std::shared_ptr<arrow::io::FileOutputStream> outfile_;
    std::unique_ptr<parquet::arrow::FileWriter> writer_;
    bool closed_;

    // Create Arrow schema
    static std::shared_ptr<arrow::Schema> CreateSchema();

    // Convert batch to Arrow table
    arrow::Result<std::shared_ptr<arrow::Table>> BatchToTable(const MinimalAlignmentBatch& batch);
};

} // namespace bam_filter

#endif // ARROW_PARQUET_WRITER_MINIMAL_HPP
