#pragma once

#include "VsearchChimeraWrapper.hpp"
#include "per_sample_table_function.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <memory>
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

		bool has_sample_id = false;
		PerSampleBindInfo sample_info;
	};

	struct GlobalState : public PerSampleGlobalState {
		// Reference table is loaded and set on the wrapper exactly once; the wrapper
		// is shared read-only during detect_batch calls.
		miint::VsearchChimeraWrapper wrapper;

		// Non-sample path only: streams the full query table on a dedicated connection.
		std::unique_ptr<QuerySequenceStream> query_stream;

		// Non-sample path only: buffered results from the last detect_batch call.
		std::vector<miint::UchimeResult> result_buffer;
		idx_t result_offset = 0;
	};

	struct LocalState : public LocalTableFunctionState {
		unique_ptr<Connection> conn; // sample_id path only

		// sample_id path only: a fresh QuerySequenceStream per claimed sample,
		// bound to `conn` and reading from a per-sample TEMP VIEW created on
		// the same connection.
		std::unique_ptr<QuerySequenceStream> query_stream;
		std::vector<miint::UchimeResult> result_buffer;
		idx_t result_offset = 0;
		Value sample_value;
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
