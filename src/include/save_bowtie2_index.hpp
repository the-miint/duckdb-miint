#pragma once
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <string>
#include <vector>

namespace duckdb {

// save_bowtie2_index(subject_table, output_path, [threads]) builds a bowtie2
// index from a single-end subject table via the gpl-boundary daemon's
// bowtie2-build tool and persists it at `output_path` (a basename prefix:
// writes <output_path>.1.bt2 … .rev.2.bt2). The bowtie2 analogue of
// save_minimap2_index — except the index files are written by the daemon, so
// this requires gpl-boundary (it advertises bowtie2-build). No >= 0.4.2 gate:
// building never uses the align-time `memory_mapped` knob.
class SaveBowtie2IndexTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string subject_table;
		std::string output_path;
		int64_t threads = 1;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names({"success", "index_path", "num_subjects"}),
		      types({LogicalType::BOOLEAN, LogicalType::VARCHAR, LogicalType::BIGINT}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		bool done = false;
		int64_t num_subjects = 0;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	struct LocalState : public LocalTableFunctionState {};

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
