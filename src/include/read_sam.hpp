#pragma once
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include <optional>
#include <unordered_map>

namespace duckdb {
class ReadSAMTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sam_paths;
		std::optional<std::unordered_map<std::string, uint64_t>> reference_lengths;
		bool include_filepath;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
		std::vector<miint::SAMRecordField> fields;

		explicit Data(const std::vector<std::string> &paths,
		              const std::optional<std::unordered_map<std::string, uint64_t>> &ref_lengths, bool include_fp)
		    : sam_paths(paths), reference_lengths(ref_lengths), include_filepath(include_fp),
		      names({"read_id",       "flags",           "reference", "position", "mapq",   "cigar",  "mate_reference",
		             "mate_position", "template_length", "tag_as",    "tag_xs",   "tag_ys", "tag_xn", "tag_xm",
		             "tag_xo",        "tag_xg",          "tag_nm",    "tag_yt",   "tag_md", "tag_sa"}),
		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::BIGINT,    // position
		             LogicalType::UTINYINT,  // mapq
		             LogicalType::VARCHAR,   // cigar
		             LogicalType::VARCHAR,   // mate_reference
		             LogicalType::BIGINT,    // mate_position
		             LogicalType::BIGINT,    // template_length
		             LogicalType::BIGINT,    // tag_as
		             LogicalType::BIGINT,    // tag_xs
		             LogicalType::BIGINT,    // tag_ys
		             LogicalType::BIGINT,    // tag_xn
		             LogicalType::BIGINT,    // tag_xm
		             LogicalType::BIGINT,    // tag_xo
		             LogicalType::BIGINT,    // tag_xg
		             LogicalType::BIGINT,    // tag_nm
		             LogicalType::VARCHAR,   // tag_yt
		             LogicalType::VARCHAR,   // tag_md
		             LogicalType::VARCHAR}), // tag_sa
		      fields({miint::SAMRecordField::READ_ID,
		              miint::SAMRecordField::FLAGS,
		              miint::SAMRecordField::REFERENCE,
		              miint::SAMRecordField::POSITION,
		              miint::SAMRecordField::MAPQ,
		              miint::SAMRecordField::CIGAR,
		              miint::SAMRecordField::MATE_REFERENCE,
		              miint::SAMRecordField::MATE_POSITION,
		              miint::SAMRecordField::TEMPLATE_LENGTH,
		              miint::SAMRecordField::TAG_AS,
		              miint::SAMRecordField::TAG_XS,
		              miint::SAMRecordField::TAG_YS,
		              miint::SAMRecordField::TAG_XN,
		              miint::SAMRecordField::TAG_XM,
		              miint::SAMRecordField::TAG_XO,
		              miint::SAMRecordField::TAG_XG,
		              miint::SAMRecordField::TAG_NM,
		              miint::SAMRecordField::TAG_YT,
		              miint::SAMRecordField::TAG_MD,
		              miint::SAMRecordField::TAG_SA}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::vector<std::unique_ptr<miint::SAMReader>> readers;
		std::vector<std::string> filepaths;
		size_t current_file_idx;
		bool finished;

		idx_t MaxThreads() const override {
			return 4;
		}

		GlobalState(const std::vector<std::string> &paths,
		            const std::optional<std::unordered_map<std::string, uint64_t>> &ref_lengths)
		    : current_file_idx(0), finished(false) {
			filepaths = paths;
			for (const auto &path : paths) {
				if (ref_lengths.has_value()) {
					readers.push_back(std::make_unique<miint::SAMReader>(path, ref_lengths.value()));
				} else {
					readers.push_back(std::make_unique<miint::SAMReader>(path));
				}
			}
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void SetResultVector(Vector &result_vector, const miint::SAMRecordField &field,
	                            const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorString(Vector &result_vector, const miint::SAMRecordField &field,
	                                  const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorStringNullable(Vector &result_vector, const miint::SAMRecordField &field,
	                                          const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorUInt8(Vector &result_vector, const miint::SAMRecordField &field,
	                                 const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorUInt16(Vector &result_vector, const miint::SAMRecordField &field,
	                                  const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorInt64(Vector &result_vector, const miint::SAMRecordField &field,
	                                 const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorInt64Nullable(Vector &result_vector, const miint::SAMRecordField &field,
	                                         const std::vector<miint::SAMRecord> &records);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath, size_t num_records);

	static TableFunction GetFunction();
};
}; // namespace duckdb
