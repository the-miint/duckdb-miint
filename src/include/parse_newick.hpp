#pragma once

#include "NewickTree.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

class ParseNewickTableFunction {
public:
	// Row of output data for a single node
	struct NodeRow {
		int64_t node_index;
		std::string name;
		double branch_length;
		std::optional<int64_t> edge_id;
		std::optional<int64_t> parent_index;
		bool is_tip;
	};

	struct Data : public TableFunctionData {
		std::vector<std::string> file_paths;
		bool include_filepath;
		bool uses_stdin;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(const std::vector<std::string> &paths, bool include_fp, bool stdin_used)
		    : file_paths(paths), include_filepath(include_fp), uses_stdin(stdin_used),
		      names({"node_index", "name", "branch_length", "edge_id", "parent_index", "is_tip"}),
		      types({LogicalType::BIGINT, LogicalType::VARCHAR, LogicalType::DOUBLE, LogicalType::BIGINT,
		             LogicalType::BIGINT, LogicalType::BOOLEAN}) {
			if (include_filepath) {
				names.emplace_back("filepath");
				types.emplace_back(LogicalType::VARCHAR);
			}
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::vector<std::string> file_paths;
		size_t next_file_idx;
		bool uses_stdin;

		idx_t MaxThreads() const override {
			if (uses_stdin) {
				return 1;
			}
			return std::min<idx_t>(file_paths.size(), std::thread::hardware_concurrency());
		}

		GlobalState(const std::vector<std::string> &paths, bool stdin_used)
		    : file_paths(paths), next_file_idx(0), uses_stdin(stdin_used) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		size_t current_file_idx;
		bool has_file;
		std::vector<NodeRow> current_rows;
		size_t current_row_idx;
		std::string current_filepath;

		LocalState() : current_file_idx(0), has_file(false), current_row_idx(0) {
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	// Helper to read and parse a newick file (handles gzip)
	static std::string ReadNewickFile(const std::string &path);

	// Convert parsed tree to row format
	static std::vector<NodeRow> TreeToRows(const miint::NewickTree &tree);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
