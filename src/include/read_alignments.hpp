#pragma once
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <optional>
#include <thread>
#include <unordered_map>
#include <vector>

namespace duckdb {
class ReadAlignmentsTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> sam_paths;
		std::optional<std::string> reference_lengths_table;
		bool include_filepath;
		bool include_seq_qual;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
		std::vector<miint::SAMRecordField> fields;

		explicit Data(const std::vector<std::string> &paths, const std::optional<std::string> &ref_table,
		              bool include_fp, bool include_sq)
		    : sam_paths(paths), reference_lengths_table(ref_table), include_filepath(include_fp),
		      include_seq_qual(include_sq),
		      names({"read_id", "flags",          "reference",     "position",        "stop_position", "mapq",
		             "cigar",   "mate_reference", "mate_position", "template_length", "tag_as",        "tag_xs",
		             "tag_ys",  "tag_xn",         "tag_xm",        "tag_xo",          "tag_xg",        "tag_nm",
		             "tag_yt",  "tag_md",         "tag_sa"}),
		      types({LogicalType::VARCHAR,   // read_id
		             LogicalType::USMALLINT, // flags
		             LogicalType::VARCHAR,   // reference
		             LogicalType::BIGINT,    // position
		             LogicalType::BIGINT,    // stop_position
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
		      fields({miint::SAMRecordField::READ_ID,       miint::SAMRecordField::FLAGS,
		              miint::SAMRecordField::REFERENCE,     miint::SAMRecordField::POSITION,
		              miint::SAMRecordField::STOP_POSITION, miint::SAMRecordField::MAPQ,
		              miint::SAMRecordField::CIGAR,         miint::SAMRecordField::MATE_REFERENCE,
		              miint::SAMRecordField::MATE_POSITION, miint::SAMRecordField::TEMPLATE_LENGTH,
		              miint::SAMRecordField::TAG_AS,        miint::SAMRecordField::TAG_XS,
		              miint::SAMRecordField::TAG_YS,        miint::SAMRecordField::TAG_XN,
		              miint::SAMRecordField::TAG_XM,        miint::SAMRecordField::TAG_XO,
		              miint::SAMRecordField::TAG_XG,        miint::SAMRecordField::TAG_NM,
		              miint::SAMRecordField::TAG_YT,        miint::SAMRecordField::TAG_MD,
		              miint::SAMRecordField::TAG_SA}) {
			if (include_seq_qual) {
				names.emplace_back("sequence");
				types.emplace_back(LogicalType::VARCHAR);
				names.emplace_back("qual");
				types.emplace_back(LogicalType::LIST(LogicalType::UTINYINT));
			}
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
		size_t next_file_idx;

		idx_t MaxThreads() const override {
			return std::min<idx_t>(readers.size(), std::thread::hardware_concurrency());
		}

		GlobalState(const std::vector<std::string> &paths,
		            std::optional<std::unordered_map<std::string, uint64_t>> ref_lengths, bool include_seq_qual)
		    : next_file_idx(0) {
			filepaths = paths;
			for (const auto &path : paths) {
				if (ref_lengths.has_value()) {
					readers.push_back(std::make_unique<miint::SAMReader>(path, ref_lengths.value(), include_seq_qual));
				} else {
					readers.push_back(std::make_unique<miint::SAMReader>(path, include_seq_qual));
				}
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;

		LocalState() : current_file_idx(0), has_file(false) {
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                      GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values);
	static void SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values);
	static void SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values);
	static void SetResultVectorInt64Nullable(Vector &result_vector, const std::vector<int64_t> &values);
	static void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
}; // namespace duckdb
