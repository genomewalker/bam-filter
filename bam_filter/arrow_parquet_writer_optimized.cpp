/*
 * Optimized Arrow Parquet Writer Implementation
 *
 * Storage-efficient schema for SAM/BAM to Parquet conversion.
 */

#include "arrow_parquet_writer_optimized.hpp"

#include <arrow/builder.h>
#include <arrow/table.h>
#include <arrow/type.h>
#include <parquet/arrow/writer.h>
#include <parquet/properties.h>

#include <stdexcept>
#include <iostream>
#include <filesystem>

namespace bam_filter {

#define ARROW_CHECK(expr)                                           \
    do {                                                            \
        arrow::Status _s = (expr);                                  \
        if (!_s.ok()) {                                             \
            throw std::runtime_error("Arrow error: " + _s.ToString()); \
        }                                                           \
    } while (0)

// ============================================================================
// OptimizedParquetWriter - Single table mode
// ============================================================================

std::shared_ptr<arrow::Schema> OptimizedParquetWriter::CreateSchema() {
    return arrow::schema({
        // Core fields
        arrow::field("read_id", arrow::uint64()),
        arrow::field("read_name", arrow::utf8()),
        arrow::field("ref_id", arrow::uint32()),        // Dictionary ID
        arrow::field("position", arrow::int32()),
        arrow::field("mapq", arrow::uint8()),
        arrow::field("flag", arrow::uint16()),
        arrow::field("cigar", arrow::utf8()),
        arrow::field("template_length", arrow::int32()),
        arrow::field("mate_ref_id", arrow::int32()),
        arrow::field("mate_position", arrow::int32()),

        // Sequences - 2-bit packed
        arrow::field("sequence_packed", arrow::binary()),
        arrow::field("sequence_length", arrow::uint16()),
        arrow::field("quality", arrow::binary()),

        // Hot tags as columns
        arrow::field("AS", arrow::int16()),             // Alignment score
        arrow::field("NM", arrow::uint16()),            // Edit distance
        arrow::field("XS", arrow::int16()),             // Secondary score
        arrow::field("MD", arrow::utf8()),              // Mismatch string

        // Cold tags (hot tags stripped)
        arrow::field("tags_cold", arrow::utf8()),

        // Derived fields (precomputed for fast queries)
        arrow::field("aligned_length", arrow::uint16()),  // reference span from CIGAR (compresses well)
        arrow::field("is_reverse", arrow::boolean()),     // flag & 0x10

        // Pipeline results
        arrow::field("ani", arrow::float32()),
        arrow::field("zs_score", arrow::float32()),     // Log-likelihood alignment score (ZS:f)
        arrow::field("pmd_score", arrow::float32()),    // Post-mortem damage score (PM:f)
        arrow::field("filter_passed", arrow::boolean()),
        arrow::field("lca_taxid", arrow::int32()),
        arrow::field("reassigned_ref_id", arrow::int32()),
        arrow::field("zp_posterior", arrow::float32()),

        // Optional
        arrow::field("read_group", arrow::utf8()),
    });
}

OptimizedParquetWriter::OptimizedParquetWriter(const std::string& filename, int compression_level)
    : closed_(false) {

    schema_ = CreateSchema();

    auto maybe_outfile = arrow::io::FileOutputStream::Open(filename);
    if (!maybe_outfile.ok()) {
        throw std::runtime_error("Failed to open output file: " + filename);
    }
    outfile_ = *maybe_outfile;

    auto builder = parquet::WriterProperties::Builder();
    builder.compression(parquet::Compression::ZSTD);
    builder.compression_level(compression_level);

    // Dictionary encoding for low-cardinality columns
    builder.enable_dictionary("ref_id");
    builder.enable_dictionary("cigar");
    builder.enable_dictionary("read_group");

    // Disable dictionary for high-cardinality columns
    builder.disable_dictionary("read_name");
    builder.disable_dictionary("sequence_packed");
    builder.disable_dictionary("quality");
    builder.disable_dictionary("MD");
    builder.disable_dictionary("tags_cold");

    auto writer_properties = builder.build();

    auto arrow_properties = parquet::ArrowWriterProperties::Builder()
        .store_schema()
        ->build();

    auto result = parquet::arrow::FileWriter::Open(
        *schema_,
        arrow::default_memory_pool(),
        outfile_,
        writer_properties,
        arrow_properties
    );

    if (!result.ok()) {
        throw std::runtime_error("Failed to create Parquet writer: " + result.status().ToString());
    }
    writer_ = std::move(*result);
}

OptimizedParquetWriter::~OptimizedParquetWriter() {
    if (!closed_) {
        try { Close(); } catch (...) {}
    }
}

// Helper to bulk-append strings to StringBuilder
static inline arrow::Status AppendStrings(arrow::StringBuilder& builder,
                                          const std::vector<std::string>& values) {
    for (const auto& v : values) {
        ARROW_RETURN_NOT_OK(builder.Append(v));
    }
    return arrow::Status::OK();
}

// Helper to bulk-append binary data to BinaryBuilder
static inline arrow::Status AppendBinary(arrow::BinaryBuilder& builder,
                                         const std::vector<std::string>& values) {
    for (const auto& v : values) {
        ARROW_RETURN_NOT_OK(builder.Append(v));
    }
    return arrow::Status::OK();
}

arrow::Result<std::shared_ptr<arrow::Table>> OptimizedParquetWriter::BatchToTable(
    const OptimizedAlignmentBatch& batch) {

    if (batch.size() == 0) {
        return arrow::Table::MakeEmpty(schema_);
    }

    const size_t n = batch.size();

    // Builders for all columns
    arrow::UInt64Builder read_id_builder;
    arrow::StringBuilder read_name_builder;
    arrow::UInt32Builder ref_id_builder;
    arrow::Int32Builder position_builder;
    arrow::UInt8Builder mapq_builder;
    arrow::UInt16Builder flag_builder;
    arrow::StringBuilder cigar_builder;
    arrow::Int32Builder template_length_builder;
    arrow::Int32Builder mate_ref_id_builder;
    arrow::Int32Builder mate_position_builder;
    arrow::BinaryBuilder sequence_packed_builder;
    arrow::UInt16Builder sequence_length_builder;
    arrow::BinaryBuilder quality_builder;
    arrow::Int16Builder AS_builder;
    arrow::UInt16Builder NM_builder;
    arrow::Int16Builder XS_builder;
    arrow::StringBuilder MD_builder;
    arrow::StringBuilder tags_cold_builder;
    arrow::UInt16Builder aligned_length_builder;
    arrow::BooleanBuilder is_reverse_builder;
    arrow::FloatBuilder ani_builder;
    arrow::FloatBuilder zs_score_builder;
    arrow::FloatBuilder pmd_score_builder;
    arrow::BooleanBuilder filter_passed_builder;
    arrow::Int32Builder lca_taxid_builder;
    arrow::Int32Builder reassigned_ref_id_builder;
    arrow::FloatBuilder zp_posterior_builder;
    arrow::StringBuilder read_group_builder;

    // BULK APPEND for primitive types - use AppendValues with raw pointers (O(1) copy)
    // This is 10-50x faster than per-row Append() calls

    // uint64 arrays - direct pointer bulk append
    ARROW_CHECK(read_id_builder.AppendValues(batch.read_ids.data(), n));

    // uint32 arrays
    ARROW_CHECK(ref_id_builder.AppendValues(batch.ref_ids.data(), n));

    // int32 arrays
    ARROW_CHECK(position_builder.AppendValues(batch.positions.data(), n));
    ARROW_CHECK(template_length_builder.AppendValues(batch.template_lengths.data(), n));
    ARROW_CHECK(mate_ref_id_builder.AppendValues(batch.mate_ref_ids.data(), n));
    ARROW_CHECK(mate_position_builder.AppendValues(batch.mate_positions.data(), n));
    ARROW_CHECK(lca_taxid_builder.AppendValues(batch.lca_taxid.data(), n));
    ARROW_CHECK(reassigned_ref_id_builder.AppendValues(batch.reassigned_ref_id.data(), n));

    // uint16 arrays
    ARROW_CHECK(flag_builder.AppendValues(batch.flags.data(), n));
    ARROW_CHECK(sequence_length_builder.AppendValues(batch.sequence_lengths.data(), n));
    ARROW_CHECK(NM_builder.AppendValues(batch.NM.data(), n));

    // int16 arrays
    ARROW_CHECK(AS_builder.AppendValues(batch.AS.data(), n));
    ARROW_CHECK(XS_builder.AppendValues(batch.XS.data(), n));

    // uint8 arrays
    ARROW_CHECK(mapq_builder.AppendValues(batch.mapqs.data(), n));

    // Derived fields - aligned_length
    ARROW_CHECK(aligned_length_builder.AppendValues(batch.aligned_lengths.data(), n));

    // float arrays
    ARROW_CHECK(ani_builder.AppendValues(batch.ani.data(), n));
    ARROW_CHECK(zs_score_builder.AppendValues(batch.zs_score.data(), n));
    ARROW_CHECK(pmd_score_builder.AppendValues(batch.pmd_score.data(), n));
    ARROW_CHECK(zp_posterior_builder.AppendValues(batch.zp_posterior.data(), n));

    // boolean arrays - need to convert vector<bool> to uint8_t array
    {
        std::vector<uint8_t> bool_vals(n);
        for (size_t i = 0; i < n; ++i) {
            bool_vals[i] = batch.is_reverse[i] ? 1 : 0;
        }
        ARROW_CHECK(is_reverse_builder.AppendValues(bool_vals.data(), n));
    }
    {
        std::vector<uint8_t> bool_vals(n);
        for (size_t i = 0; i < n; ++i) {
            bool_vals[i] = batch.filter_passed[i] ? 1 : 0;
        }
        ARROW_CHECK(filter_passed_builder.AppendValues(bool_vals.data(), n));
    }

    // String arrays - still need per-element append but with reserved capacity
    ARROW_CHECK(read_name_builder.Reserve(n));
    ARROW_CHECK(cigar_builder.Reserve(n));
    ARROW_CHECK(MD_builder.Reserve(n));
    ARROW_CHECK(tags_cold_builder.Reserve(n));
    ARROW_CHECK(read_group_builder.Reserve(n));

    ARROW_CHECK(AppendStrings(read_name_builder, batch.read_names));
    ARROW_CHECK(AppendStrings(cigar_builder, batch.cigars));
    ARROW_CHECK(AppendStrings(MD_builder, batch.MD));
    ARROW_CHECK(AppendStrings(tags_cold_builder, batch.tags_cold));
    ARROW_CHECK(AppendStrings(read_group_builder, batch.read_groups));

    // Binary arrays
    ARROW_CHECK(sequence_packed_builder.Reserve(n));
    ARROW_CHECK(quality_builder.Reserve(n));

    ARROW_CHECK(AppendBinary(sequence_packed_builder, batch.sequences_packed));
    ARROW_CHECK(AppendBinary(quality_builder, batch.qualities));

    // Finish arrays
    std::shared_ptr<arrow::Array> read_id_array, read_name_array, ref_id_array;
    std::shared_ptr<arrow::Array> position_array, mapq_array, flag_array, cigar_array;
    std::shared_ptr<arrow::Array> template_length_array, mate_ref_id_array, mate_position_array;
    std::shared_ptr<arrow::Array> sequence_packed_array, sequence_length_array, quality_array;
    std::shared_ptr<arrow::Array> AS_array, NM_array, XS_array, MD_array, tags_cold_array;
    std::shared_ptr<arrow::Array> aligned_length_array, is_reverse_array;
    std::shared_ptr<arrow::Array> ani_array, zs_score_array, pmd_score_array, filter_passed_array;
    std::shared_ptr<arrow::Array> lca_taxid_array, reassigned_ref_id_array, zp_posterior_array;
    std::shared_ptr<arrow::Array> read_group_array;

    ARROW_CHECK(read_id_builder.Finish(&read_id_array));
    ARROW_CHECK(read_name_builder.Finish(&read_name_array));
    ARROW_CHECK(ref_id_builder.Finish(&ref_id_array));
    ARROW_CHECK(position_builder.Finish(&position_array));
    ARROW_CHECK(mapq_builder.Finish(&mapq_array));
    ARROW_CHECK(flag_builder.Finish(&flag_array));
    ARROW_CHECK(cigar_builder.Finish(&cigar_array));
    ARROW_CHECK(template_length_builder.Finish(&template_length_array));
    ARROW_CHECK(mate_ref_id_builder.Finish(&mate_ref_id_array));
    ARROW_CHECK(mate_position_builder.Finish(&mate_position_array));
    ARROW_CHECK(sequence_packed_builder.Finish(&sequence_packed_array));
    ARROW_CHECK(sequence_length_builder.Finish(&sequence_length_array));
    ARROW_CHECK(quality_builder.Finish(&quality_array));
    ARROW_CHECK(AS_builder.Finish(&AS_array));
    ARROW_CHECK(NM_builder.Finish(&NM_array));
    ARROW_CHECK(XS_builder.Finish(&XS_array));
    ARROW_CHECK(MD_builder.Finish(&MD_array));
    ARROW_CHECK(tags_cold_builder.Finish(&tags_cold_array));
    ARROW_CHECK(aligned_length_builder.Finish(&aligned_length_array));
    ARROW_CHECK(is_reverse_builder.Finish(&is_reverse_array));
    ARROW_CHECK(ani_builder.Finish(&ani_array));
    ARROW_CHECK(zs_score_builder.Finish(&zs_score_array));
    ARROW_CHECK(pmd_score_builder.Finish(&pmd_score_array));
    ARROW_CHECK(filter_passed_builder.Finish(&filter_passed_array));
    ARROW_CHECK(lca_taxid_builder.Finish(&lca_taxid_array));
    ARROW_CHECK(reassigned_ref_id_builder.Finish(&reassigned_ref_id_array));
    ARROW_CHECK(zp_posterior_builder.Finish(&zp_posterior_array));
    ARROW_CHECK(read_group_builder.Finish(&read_group_array));

    return arrow::Table::Make(schema_, {
        read_id_array, read_name_array, ref_id_array, position_array,
        mapq_array, flag_array, cigar_array, template_length_array,
        mate_ref_id_array, mate_position_array,
        sequence_packed_array, sequence_length_array, quality_array,
        AS_array, NM_array, XS_array, MD_array, tags_cold_array,
        aligned_length_array, is_reverse_array,
        ani_array, zs_score_array, pmd_score_array, filter_passed_array,
        lca_taxid_array, reassigned_ref_id_array, zp_posterior_array,
        read_group_array
    });
}

arrow::Status OptimizedParquetWriter::WriteBatch(const OptimizedAlignmentBatch& batch) {
    if (closed_) return arrow::Status::Invalid("Writer is closed");
    if (batch.size() == 0) return arrow::Status::OK();

    auto table_result = BatchToTable(batch);
    if (!table_result.ok()) return table_result.status();

    return writer_->WriteTable(*(*table_result), batch.size());
}

arrow::Status OptimizedParquetWriter::Close() {
    if (closed_) return arrow::Status::OK();
    closed_ = true;
    auto status = writer_->Close();
    if (!status.ok()) return status;
    return outfile_->Close();
}

// ============================================================================
// NormalizedParquetWriter - Separate reads/alignments tables
// ============================================================================

std::shared_ptr<arrow::Schema> NormalizedParquetWriter::CreateReadsSchema() {
    return arrow::schema({
        arrow::field("read_id", arrow::uint64()),
        arrow::field("read_name", arrow::utf8()),
        arrow::field("sequence_packed", arrow::binary()),
        arrow::field("sequence_length", arrow::uint16()),
        arrow::field("quality", arrow::binary()),
        arrow::field("read_group", arrow::utf8()),
    });
}

std::shared_ptr<arrow::Schema> NormalizedParquetWriter::CreateAlignmentsSchema() {
    return arrow::schema({
        // Foreign key to reads
        arrow::field("read_id", arrow::uint64()),
        arrow::field("ref_id", arrow::uint32()),
        arrow::field("position", arrow::int32()),
        arrow::field("mapq", arrow::uint8()),
        arrow::field("flag", arrow::uint16()),
        arrow::field("cigar", arrow::utf8()),
        arrow::field("template_length", arrow::int32()),
        arrow::field("mate_ref_id", arrow::int32()),
        arrow::field("mate_position", arrow::int32()),

        // Hot tags
        arrow::field("AS", arrow::int16()),
        arrow::field("NM", arrow::uint16()),
        arrow::field("XS", arrow::int16()),
        arrow::field("MD", arrow::utf8()),
        arrow::field("tags_cold", arrow::utf8()),

        // Derived fields (precomputed for fast queries)
        arrow::field("aligned_length", arrow::uint16()),  // reference span from CIGAR (compresses well)
        arrow::field("is_reverse", arrow::boolean()),     // flag & 0x10

        // Pipeline results
        arrow::field("ani", arrow::float32()),
        arrow::field("zs_score", arrow::float32()),     // Log-likelihood alignment score (ZS:f)
        arrow::field("pmd_score", arrow::float32()),    // Post-mortem damage score (PM:f)
        arrow::field("filter_passed", arrow::boolean()),
        arrow::field("lca_taxid", arrow::int32()),
        arrow::field("reassigned_ref_id", arrow::int32()),
        arrow::field("zp_posterior", arrow::float32()),
    });
}

std::shared_ptr<arrow::Schema> NormalizedParquetWriter::CreateReferencesSchema() {
    return arrow::schema({
        arrow::field("ref_id", arrow::uint32()),
        arrow::field("ref_name", arrow::utf8()),
        arrow::field("ref_length", arrow::int64()),
        arrow::field("taxid", arrow::int32()),
    });
}

NormalizedParquetWriter::NormalizedParquetWriter(const std::string& output_dir, int compression_level)
    : output_dir_(output_dir), compression_level_(compression_level), closed_(false) {

    // Create output directory
    std::filesystem::create_directories(output_dir);

    reads_schema_ = CreateReadsSchema();
    alignments_schema_ = CreateAlignmentsSchema();
    references_schema_ = CreateReferencesSchema();

    // Writer properties
    auto builder = parquet::WriterProperties::Builder();
    builder.compression(parquet::Compression::ZSTD);
    builder.compression_level(compression_level);
    builder.enable_dictionary("ref_id");
    builder.enable_dictionary("cigar");
    builder.disable_dictionary("read_name");
    builder.disable_dictionary("sequence_packed");
    builder.disable_dictionary("quality");
    auto writer_properties = builder.build();

    auto arrow_properties = parquet::ArrowWriterProperties::Builder()
        .store_schema()
        ->build();

    // Open reads file
    std::string reads_path = output_dir + "/reads.parquet";
    auto reads_outfile = arrow::io::FileOutputStream::Open(reads_path);
    if (!reads_outfile.ok()) {
        throw std::runtime_error("Failed to open reads file: " + reads_path);
    }
    reads_outfile_ = *reads_outfile;

    auto reads_result = parquet::arrow::FileWriter::Open(
        *reads_schema_, arrow::default_memory_pool(), reads_outfile_,
        writer_properties, arrow_properties);
    if (!reads_result.ok()) {
        throw std::runtime_error("Failed to create reads writer");
    }
    reads_writer_ = std::move(*reads_result);

    // Open alignments file
    std::string alignments_path = output_dir + "/alignments.parquet";
    auto alignments_outfile = arrow::io::FileOutputStream::Open(alignments_path);
    if (!alignments_outfile.ok()) {
        throw std::runtime_error("Failed to open alignments file: " + alignments_path);
    }
    alignments_outfile_ = *alignments_outfile;

    auto alignments_result = parquet::arrow::FileWriter::Open(
        *alignments_schema_, arrow::default_memory_pool(), alignments_outfile_,
        writer_properties, arrow_properties);
    if (!alignments_result.ok()) {
        throw std::runtime_error("Failed to create alignments writer");
    }
    alignments_writer_ = std::move(*alignments_result);
}

NormalizedParquetWriter::~NormalizedParquetWriter() {
    if (!closed_) {
        try { Close(); } catch (...) {}
    }
}

arrow::Result<std::shared_ptr<arrow::Table>> NormalizedParquetWriter::ReadBatchToTable(
    const ReadBatch& batch) {

    if (batch.size() == 0) {
        return arrow::Table::MakeEmpty(reads_schema_);
    }

    size_t n = batch.size();

    arrow::UInt64Builder read_id_builder;
    arrow::StringBuilder read_name_builder;
    arrow::BinaryBuilder sequence_packed_builder;
    arrow::UInt16Builder sequence_length_builder;
    arrow::BinaryBuilder quality_builder;
    arrow::StringBuilder read_group_builder;

    ARROW_CHECK(read_id_builder.Reserve(n));
    ARROW_CHECK(read_name_builder.Reserve(n));
    ARROW_CHECK(sequence_packed_builder.Reserve(n));
    ARROW_CHECK(sequence_length_builder.Reserve(n));
    ARROW_CHECK(quality_builder.Reserve(n));
    ARROW_CHECK(read_group_builder.Reserve(n));

    for (size_t i = 0; i < n; ++i) {
        ARROW_CHECK(read_id_builder.Append(batch.read_ids[i]));
        ARROW_CHECK(read_name_builder.Append(batch.read_names[i]));
        ARROW_CHECK(sequence_packed_builder.Append(batch.sequences_packed[i]));
        ARROW_CHECK(sequence_length_builder.Append(batch.sequence_lengths[i]));
        ARROW_CHECK(quality_builder.Append(batch.qualities[i]));
        ARROW_CHECK(read_group_builder.Append(batch.read_groups[i]));
    }

    std::shared_ptr<arrow::Array> read_id_array, read_name_array;
    std::shared_ptr<arrow::Array> sequence_packed_array, sequence_length_array;
    std::shared_ptr<arrow::Array> quality_array, read_group_array;

    ARROW_CHECK(read_id_builder.Finish(&read_id_array));
    ARROW_CHECK(read_name_builder.Finish(&read_name_array));
    ARROW_CHECK(sequence_packed_builder.Finish(&sequence_packed_array));
    ARROW_CHECK(sequence_length_builder.Finish(&sequence_length_array));
    ARROW_CHECK(quality_builder.Finish(&quality_array));
    ARROW_CHECK(read_group_builder.Finish(&read_group_array));

    return arrow::Table::Make(reads_schema_, {
        read_id_array, read_name_array, sequence_packed_array,
        sequence_length_array, quality_array, read_group_array
    });
}

arrow::Result<std::shared_ptr<arrow::Table>> NormalizedParquetWriter::AlignmentBatchToTable(
    const AlignmentOnlyBatch& batch) {

    if (batch.size() == 0) {
        return arrow::Table::MakeEmpty(alignments_schema_);
    }

    const size_t n = batch.size();

    arrow::UInt64Builder read_id_builder;
    arrow::UInt32Builder ref_id_builder;
    arrow::Int32Builder position_builder;
    arrow::UInt8Builder mapq_builder;
    arrow::UInt16Builder flag_builder;
    arrow::StringBuilder cigar_builder;
    arrow::Int32Builder template_length_builder;
    arrow::Int32Builder mate_ref_id_builder;
    arrow::Int32Builder mate_position_builder;
    arrow::Int16Builder AS_builder;
    arrow::UInt16Builder NM_builder;
    arrow::Int16Builder XS_builder;
    arrow::StringBuilder MD_builder;
    arrow::StringBuilder tags_cold_builder;
    arrow::UInt16Builder aligned_length_builder;
    arrow::BooleanBuilder is_reverse_builder;
    arrow::FloatBuilder ani_builder;
    arrow::FloatBuilder zs_score_builder;
    arrow::FloatBuilder pmd_score_builder;
    arrow::BooleanBuilder filter_passed_builder;
    arrow::Int32Builder lca_taxid_builder;
    arrow::Int32Builder reassigned_ref_id_builder;
    arrow::FloatBuilder zp_posterior_builder;

    // BULK APPEND for primitive types - raw pointer O(1) copy
    ARROW_CHECK(read_id_builder.AppendValues(batch.read_ids.data(), n));
    ARROW_CHECK(ref_id_builder.AppendValues(batch.ref_ids.data(), n));
    ARROW_CHECK(position_builder.AppendValues(batch.positions.data(), n));
    ARROW_CHECK(mapq_builder.AppendValues(batch.mapqs.data(), n));
    ARROW_CHECK(flag_builder.AppendValues(batch.flags.data(), n));
    ARROW_CHECK(template_length_builder.AppendValues(batch.template_lengths.data(), n));
    ARROW_CHECK(mate_ref_id_builder.AppendValues(batch.mate_ref_ids.data(), n));
    ARROW_CHECK(mate_position_builder.AppendValues(batch.mate_positions.data(), n));
    ARROW_CHECK(AS_builder.AppendValues(batch.AS.data(), n));
    ARROW_CHECK(NM_builder.AppendValues(batch.NM.data(), n));
    ARROW_CHECK(XS_builder.AppendValues(batch.XS.data(), n));
    ARROW_CHECK(aligned_length_builder.AppendValues(batch.aligned_lengths.data(), n));
    ARROW_CHECK(ani_builder.AppendValues(batch.ani.data(), n));
    ARROW_CHECK(zs_score_builder.AppendValues(batch.zs_score.data(), n));
    ARROW_CHECK(pmd_score_builder.AppendValues(batch.pmd_score.data(), n));
    ARROW_CHECK(lca_taxid_builder.AppendValues(batch.lca_taxid.data(), n));
    ARROW_CHECK(reassigned_ref_id_builder.AppendValues(batch.reassigned_ref_id.data(), n));
    ARROW_CHECK(zp_posterior_builder.AppendValues(batch.zp_posterior.data(), n));

    // Boolean arrays conversion
    {
        std::vector<uint8_t> bool_vals(n);
        for (size_t i = 0; i < n; ++i) {
            bool_vals[i] = batch.is_reverse[i] ? 1 : 0;
        }
        ARROW_CHECK(is_reverse_builder.AppendValues(bool_vals.data(), n));
    }
    {
        std::vector<uint8_t> bool_vals(n);
        for (size_t i = 0; i < n; ++i) {
            bool_vals[i] = batch.filter_passed[i] ? 1 : 0;
        }
        ARROW_CHECK(filter_passed_builder.AppendValues(bool_vals.data(), n));
    }

    // String arrays with reserved capacity
    ARROW_CHECK(cigar_builder.Reserve(n));
    ARROW_CHECK(MD_builder.Reserve(n));
    ARROW_CHECK(tags_cold_builder.Reserve(n));

    ARROW_CHECK(AppendStrings(cigar_builder, batch.cigars));
    ARROW_CHECK(AppendStrings(MD_builder, batch.MD));
    ARROW_CHECK(AppendStrings(tags_cold_builder, batch.tags_cold));

    std::shared_ptr<arrow::Array> read_id_array, ref_id_array, position_array;
    std::shared_ptr<arrow::Array> mapq_array, flag_array, cigar_array;
    std::shared_ptr<arrow::Array> template_length_array, mate_ref_id_array, mate_position_array;
    std::shared_ptr<arrow::Array> AS_array, NM_array, XS_array, MD_array, tags_cold_array;
    std::shared_ptr<arrow::Array> aligned_length_array, is_reverse_array;
    std::shared_ptr<arrow::Array> ani_array, zs_score_array, pmd_score_array, filter_passed_array;
    std::shared_ptr<arrow::Array> lca_taxid_array, reassigned_ref_id_array, zp_posterior_array;

    ARROW_CHECK(read_id_builder.Finish(&read_id_array));
    ARROW_CHECK(ref_id_builder.Finish(&ref_id_array));
    ARROW_CHECK(position_builder.Finish(&position_array));
    ARROW_CHECK(mapq_builder.Finish(&mapq_array));
    ARROW_CHECK(flag_builder.Finish(&flag_array));
    ARROW_CHECK(cigar_builder.Finish(&cigar_array));
    ARROW_CHECK(template_length_builder.Finish(&template_length_array));
    ARROW_CHECK(mate_ref_id_builder.Finish(&mate_ref_id_array));
    ARROW_CHECK(mate_position_builder.Finish(&mate_position_array));
    ARROW_CHECK(AS_builder.Finish(&AS_array));
    ARROW_CHECK(NM_builder.Finish(&NM_array));
    ARROW_CHECK(XS_builder.Finish(&XS_array));
    ARROW_CHECK(MD_builder.Finish(&MD_array));
    ARROW_CHECK(tags_cold_builder.Finish(&tags_cold_array));
    ARROW_CHECK(aligned_length_builder.Finish(&aligned_length_array));
    ARROW_CHECK(is_reverse_builder.Finish(&is_reverse_array));
    ARROW_CHECK(ani_builder.Finish(&ani_array));
    ARROW_CHECK(zs_score_builder.Finish(&zs_score_array));
    ARROW_CHECK(pmd_score_builder.Finish(&pmd_score_array));
    ARROW_CHECK(filter_passed_builder.Finish(&filter_passed_array));
    ARROW_CHECK(lca_taxid_builder.Finish(&lca_taxid_array));
    ARROW_CHECK(reassigned_ref_id_builder.Finish(&reassigned_ref_id_array));
    ARROW_CHECK(zp_posterior_builder.Finish(&zp_posterior_array));

    return arrow::Table::Make(alignments_schema_, {
        read_id_array, ref_id_array, position_array,
        mapq_array, flag_array, cigar_array,
        template_length_array, mate_ref_id_array, mate_position_array,
        AS_array, NM_array, XS_array, MD_array, tags_cold_array,
        aligned_length_array, is_reverse_array,
        ani_array, zs_score_array, pmd_score_array, filter_passed_array,
        lca_taxid_array, reassigned_ref_id_array, zp_posterior_array
    });
}

arrow::Status NormalizedParquetWriter::WriteReadBatch(const ReadBatch& batch) {
    if (closed_) return arrow::Status::Invalid("Writer is closed");
    if (batch.size() == 0) return arrow::Status::OK();

    auto table_result = ReadBatchToTable(batch);
    if (!table_result.ok()) return table_result.status();

    return reads_writer_->WriteTable(*(*table_result), batch.size());
}

arrow::Status NormalizedParquetWriter::WriteAlignmentBatch(const AlignmentOnlyBatch& batch) {
    if (closed_) return arrow::Status::Invalid("Writer is closed");
    if (batch.size() == 0) return arrow::Status::OK();

    auto table_result = AlignmentBatchToTable(batch);
    if (!table_result.ok()) return table_result.status();

    return alignments_writer_->WriteTable(*(*table_result), batch.size());
}

arrow::Status NormalizedParquetWriter::WriteReferences(const std::vector<ReferenceInfo>& refs) {
    if (refs.empty()) return arrow::Status::OK();

    arrow::UInt32Builder ref_id_builder;
    arrow::StringBuilder ref_name_builder;
    arrow::Int64Builder ref_length_builder;
    arrow::Int32Builder taxid_builder;

    ARROW_CHECK(ref_id_builder.Reserve(refs.size()));
    ARROW_CHECK(ref_name_builder.Reserve(refs.size()));
    ARROW_CHECK(ref_length_builder.Reserve(refs.size()));
    ARROW_CHECK(taxid_builder.Reserve(refs.size()));

    for (const auto& ref : refs) {
        ARROW_CHECK(ref_id_builder.Append(ref.ref_id));
        ARROW_CHECK(ref_name_builder.Append(ref.ref_name));
        ARROW_CHECK(ref_length_builder.Append(ref.ref_length));
        ARROW_CHECK(taxid_builder.Append(ref.taxid));
    }

    std::shared_ptr<arrow::Array> ref_id_array, ref_name_array, ref_length_array, taxid_array;
    ARROW_CHECK(ref_id_builder.Finish(&ref_id_array));
    ARROW_CHECK(ref_name_builder.Finish(&ref_name_array));
    ARROW_CHECK(ref_length_builder.Finish(&ref_length_array));
    ARROW_CHECK(taxid_builder.Finish(&taxid_array));

    auto table = arrow::Table::Make(references_schema_, {
        ref_id_array, ref_name_array, ref_length_array, taxid_array
    });

    std::string ref_path = output_dir_ + "/references.parquet";
    auto outfile = arrow::io::FileOutputStream::Open(ref_path);
    if (!outfile.ok()) {
        return arrow::Status::IOError("Failed to open references file");
    }

    auto builder = parquet::WriterProperties::Builder();
    builder.compression(parquet::Compression::ZSTD);
    builder.compression_level(compression_level_);
    auto props = builder.build();

    return parquet::arrow::WriteTable(*table, arrow::default_memory_pool(), *outfile, 1000000, props);
}

arrow::Status NormalizedParquetWriter::Close() {
    if (closed_) return arrow::Status::OK();
    closed_ = true;

    auto status = reads_writer_->Close();
    if (!status.ok()) return status;
    status = reads_outfile_->Close();
    if (!status.ok()) return status;

    status = alignments_writer_->Close();
    if (!status.ok()) return status;
    return alignments_outfile_->Close();
}

} // namespace bam_filter
