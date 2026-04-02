#pragma once

#include "read_newick.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

class ReadJplaceNewickTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp);
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::vector<std::string> file_paths;
		size_t next_file_idx;
		FileSystem &fs;

		idx_t MaxThreads() const override {
			return std::min<idx_t>(file_paths.size(), std::thread::hardware_concurrency());
		}

		GlobalState(const std::vector<std::string> &paths, FileSystem &fs_ref)
		    : file_paths(paths), next_file_idx(0), fs(fs_ref) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		bool has_file = false;
		std::vector<ReadNewickTableFunction::NodeRow> current_rows;
		size_t current_row_idx = 0;
		std::string current_filepath;
	};

	// Extract the "tree" field value from jplace JSON content.
	static std::string ExtractTreeFromJplaceContent(const std::string &content, const std::string &path);

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
