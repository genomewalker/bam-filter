/*
 * Arrow Parquet Writer - High-performance C++ wrapper for writing SAM/BAM to Parquet
 *
 * Uses Arrow C++ API directly for maximum performance.
 * No DuckDB, no Python overhead - pure C++ with Cython bindings.
 */

#ifndef ARROW_PARQUET_WRITER_HPP
#define ARROW_PARQUET_WRITER_HPP

#include <arrow/api.h>
#include <arrow/io/api.h>
#include <parquet/arrow/writer.h>
#include <memory>
#include <string>
#include <vector>

namespace bam_filter {

// Columnar batch buffer for accumulating records before writing
struct AlignmentBatch {
    std::vector<uint64_t> read_ids;
    std::vector<std::string> read_names;
    std::vector<uint32_t> ref_ids;
    std::vector<std::string> ref_names;
    std::vector<int32_t> positions;
    std::vector<int32_t> end_positions;
    std::vector<uint8_t> mapqs;
    std::vector<uint16_t> flags;
    std::vector<uint16_t> alignment_lengths;
    std::vector<int32_t> template_lengths;
    std::vector<int32_t> mate_ref_ids;
    std::vector<int32_t> mate_positions;

    // Tags (use -1 for NULL/missing values)
    std::vector<int32_t> alignment_scores;
    std::vector<int32_t> xs_scores;
    std::vector<int32_t> edit_distances;
    std::vector<int32_t> num_mismatches;
    std::vector<int32_t> num_gap_opens;
    std::vector<int32_t> num_gap_extensions;
    std::vector<std::string> md_strings;
    std::vector<float> anis;
    std::vector<float> pmd_scores;
    std::vector<float> zs_scores;
    std::vector<float> zp_posteriors;
    std::vector<int32_t> lca_taxids;
    std::vector<int32_t> reassigned_ref_ids;
    std::vector<bool> filter_passed;
    std::vector<std::string> read_groups;
    std::vector<std::string> cigars;
    std::vector<std::string> sequences;
    std::vector<std::string> qualities;

    size_t size() const { return read_ids.size(); }

    void clear() {
        read_ids.clear();
        read_names.clear();
        ref_ids.clear();
        ref_names.clear();
        positions.clear();
        end_positions.clear();
        mapqs.clear();
        flags.clear();
        alignment_lengths.clear();
        template_lengths.clear();
        mate_ref_ids.clear();
        mate_positions.clear();
        alignment_scores.clear();
        xs_scores.clear();
        edit_distances.clear();
        num_mismatches.clear();
        num_gap_opens.clear();
        num_gap_extensions.clear();
        md_strings.clear();
        anis.clear();
        pmd_scores.clear();
        zs_scores.clear();
        zp_posteriors.clear();
        lca_taxids.clear();
        reassigned_ref_ids.clear();
        filter_passed.clear();
        read_groups.clear();
        cigars.clear();
        sequences.clear();
        qualities.clear();
    }

    void shrink_to_fit() {
        read_ids.shrink_to_fit();
        read_names.shrink_to_fit();
        ref_ids.shrink_to_fit();
        ref_names.shrink_to_fit();
        positions.shrink_to_fit();
        end_positions.shrink_to_fit();
        mapqs.shrink_to_fit();
        flags.shrink_to_fit();
        alignment_lengths.shrink_to_fit();
        template_lengths.shrink_to_fit();
        mate_ref_ids.shrink_to_fit();
        mate_positions.shrink_to_fit();
        alignment_scores.shrink_to_fit();
        xs_scores.shrink_to_fit();
        edit_distances.shrink_to_fit();
        num_mismatches.shrink_to_fit();
        num_gap_opens.shrink_to_fit();
        num_gap_extensions.shrink_to_fit();
        md_strings.shrink_to_fit();
        anis.shrink_to_fit();
        pmd_scores.shrink_to_fit();
        zs_scores.shrink_to_fit();
        zp_posteriors.shrink_to_fit();
        lca_taxids.shrink_to_fit();
        reassigned_ref_ids.shrink_to_fit();
        filter_passed.shrink_to_fit();
        read_groups.shrink_to_fit();
        cigars.shrink_to_fit();
        sequences.shrink_to_fit();
        qualities.shrink_to_fit();
    }

    void reserve(size_t capacity) {
        read_ids.reserve(capacity);
        read_names.reserve(capacity);
        ref_ids.reserve(capacity);
        ref_names.reserve(capacity);
        positions.reserve(capacity);
        end_positions.reserve(capacity);
        mapqs.reserve(capacity);
        flags.reserve(capacity);
        alignment_lengths.reserve(capacity);
        template_lengths.reserve(capacity);
        mate_ref_ids.reserve(capacity);
        mate_positions.reserve(capacity);
        alignment_scores.reserve(capacity);
        xs_scores.reserve(capacity);
        edit_distances.reserve(capacity);
        num_mismatches.reserve(capacity);
        num_gap_opens.reserve(capacity);
        num_gap_extensions.reserve(capacity);
        md_strings.reserve(capacity);
        anis.reserve(capacity);
        pmd_scores.reserve(capacity);
        zs_scores.reserve(capacity);
        zp_posteriors.reserve(capacity);
        lca_taxids.reserve(capacity);
        reassigned_ref_ids.reserve(capacity);
        filter_passed.reserve(capacity);
        read_groups.reserve(capacity);
        cigars.reserve(capacity);
        sequences.reserve(capacity);
        qualities.reserve(capacity);
    }
};

// Parquet writer class - wraps Arrow's ParquetFileWriter
class ParquetWriter {
public:
    ParquetWriter(const std::string& filename, int compression_level = 3);
    ~ParquetWriter();

    // Write a batch of alignments
    arrow::Status WriteBatch(const AlignmentBatch& batch);

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
    arrow::Result<std::shared_ptr<arrow::Table>> BatchToTable(const AlignmentBatch& batch);
};

} // namespace bam_filter

#endif // ARROW_PARQUET_WRITER_HPP
