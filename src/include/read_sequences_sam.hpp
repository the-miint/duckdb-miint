#pragma once
#include "SAMReader.hpp"
#include "QualScore.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <atomic>
#include <thread>
#include <vector>

namespace duckdb {
class ReadSequencesSamTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;
		bool uses_stdin;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp, bool stdin_used)
		    : file_paths(paths), include_filepath(include_fp), uses_stdin(stdin_used),
		      names({"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"}),
		      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
		             LogicalType::VARCHAR, LogicalType::LIST(LogicalType::UTINYINT),
		             LogicalType::LIST(LogicalType::UTINYINT)}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		};
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::vector<std::unique_ptr<miint::SAMReader>> readers;
		std::vector<std::string> filepaths;
		size_t next_file_idx;
		bool uses_stdin;
		std::vector<uint64_t> file_sequence_counters;

		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			auto hw_threads = std::thread::hardware_concurrency();
			if (hw_threads == 0) {
				hw_threads = 1;
			}
			return std::min<idx_t>(readers.size(), std::min<idx_t>(8, hw_threads));
		}

		GlobalState(const std::vector<std::string> &paths, bool stdin_used) : next_file_idx(0), uses_stdin(stdin_used) {
			filepaths = paths;
			for (size_t i = 0; i < paths.size(); i++) {
				readers.push_back(std::make_unique<miint::SAMReader>(paths[i], /*include_seq_qual=*/true,
				                                                     /*require_references=*/false));
				file_sequence_counters.emplace_back(1);
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

	static void SetResultVectorNull(Vector &result_vector);
	static void SetResultVectorString(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorStringNullable(Vector &result_vector, const std::vector<std::string> &values);
	static void SetResultVectorListUInt8(Vector &result_vector, const std::vector<miint::QualScore> &values);
	static void SetResultVectorFilepath(Vector &result_vector, const std::string &filepath);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
} // namespace duckdb
