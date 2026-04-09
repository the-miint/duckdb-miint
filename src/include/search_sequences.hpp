#pragma once

#include "VsearchSearchWrapper.hpp"
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

	struct GlobalState : public GlobalTableFunctionState {
		miint::VsearchSearchWrapper wrapper;

		// Lazy streaming of query sequences.
		std::unique_ptr<QuerySequenceStream> query_stream;

		// Buffered results from the last search_batch call.
		std::vector<miint::SearchResult> result_buffer;
		idx_t result_offset = 0;

		idx_t MaxThreads() const override {
			return 1;
		}
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
