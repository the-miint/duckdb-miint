#pragma once

#include "VsearchSearchWrapper.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <atomic>
#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace duckdb {

class SearchSequencesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string ref_table;
		miint::SearchParams params;

		SequenceTableSchema query_schema;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	static constexpr idx_t BATCH_SIZE = 64;
	static constexpr idx_t MAX_THREADS = 8;

	struct GlobalState : public GlobalTableFunctionState {
		miint::VsearchSearchWrapper wrapper;

		// Lazy streaming of query sequences — thread-safe, no upfront materialization.
		std::unique_ptr<QuerySequenceStream> query_stream;

		idx_t MaxThreads() const override {
			return MAX_THREADS;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		std::optional<miint::VsearchSearchWrapper::SearchHandle> handle;

		std::vector<miint::SearchResult> result_buffer;
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
