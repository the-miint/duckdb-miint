#pragma once

#include "ena_client.hpp"
#include "ena_parser.hpp"

#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace duckdb {

// `ena_searchable_fields(result_type)` — table function that hits ENA's
// `/returnFields?result=<result_type>&format=tsv` endpoint and returns the
// rows verbatim. Used as a discoverability tool: given a read_run / sample /
// study / experiment result type, list the field names, types, and
// descriptions that ENA's structured search accepts.
//
// Exactly one HTTP call is made on the first Execute invocation (response is
// ~60 rows, so no streaming / batching needed). Subsequent calls drain the
// materialized rows.
//
// Output schema:
//   field_name  VARCHAR
//   type        VARCHAR
//   description VARCHAR (nullable — ENA sometimes omits descriptions)
class ReadENASearchableFieldsTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string result_type;
		explicit Data(std::string result_type);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		std::unique_ptr<miint::ENAClient> client;
		bool fetched = false;
		// Rows of (field_name, type, description) after parsing the TSV.
		std::vector<std::array<std::string, 3>> rows;
		size_t offset = 0;

		explicit GlobalState(DatabaseInstance &db);

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
