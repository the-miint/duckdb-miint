#pragma once

#include "VsearchChimeraWrapper.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <atomic>
#include <optional>
#include <string>
#include <vector>

namespace duckdb {

class UchimeRefTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string ref_table;
		miint::UchimeParams params;

		SequenceTableSchema query_schema;
		SequenceTableSchema ref_schema;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	// Each thread claims BATCH_SIZE query IDs atomically, then fetches sequences
	// and runs chimera detection independently with no shared mutable state.
	static constexpr idx_t BATCH_SIZE = 64;
	static constexpr idx_t MAX_THREADS = 8;

	struct GlobalState : public GlobalTableFunctionState {
		// Reference DB loaded at init, shared read-only across threads.
		miint::VsearchChimeraWrapper wrapper;

		// Pre-materialized query IDs (lightweight — just strings, no sequences).
		std::vector<std::string> all_query_ids;
		std::atomic<idx_t> next_batch_offset {0};

		// Query table name and schema for per-thread ReadBatchByIds calls.
		std::string query_table;
		SequenceTableSchema query_schema;

		idx_t MaxThreads() const override {
			return MAX_THREADS;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		// Per-thread vsearch detection handle (wraps chimera_info_s).
		// Initialized in InitLocal after the global DB is ready.
		std::optional<miint::VsearchChimeraWrapper::DetectHandle> handle;

		std::vector<miint::UchimeResult> result_buffer;
		idx_t buffer_offset = 0;
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
