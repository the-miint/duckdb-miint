#pragma once
#include "SFFReader.hpp"
#include "QualScore.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <thread>
#include <vector>

namespace duckdb {
class ReadSequencesSFFTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;
		bool trim;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp, bool do_trim)
		    : file_paths(paths), include_filepath(include_fp), trim(do_trim),
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
		std::vector<std::unique_ptr<miint::SFFReader>> readers; // Lazily opened (nullptr until claimed)
		std::vector<std::string> filepaths;
		bool trim;
		size_t next_file_idx;
		std::vector<uint64_t>
		    file_sequence_counters; // Per-file sequence counters (no atomic needed - file access is exclusive)

		idx_t MaxThreads() const override {
			auto hw_threads = std::thread::hardware_concurrency();
			if (hw_threads == 0) {
				hw_threads = 1;
			}
			return std::min<idx_t>(filepaths.size(), std::min<idx_t>(8, hw_threads));
		}

		GlobalState(const std::vector<std::string> &paths, bool do_trim)
		    : readers(paths.size()), filepaths(paths), trim(do_trim), next_file_idx(0),
		      file_sequence_counters(paths.size(), 1) {
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

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};
} // namespace duckdb
