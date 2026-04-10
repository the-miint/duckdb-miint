#pragma once

#include "ena_client.hpp"
#include "ena_parser.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <mutex>
#include <vector>

namespace duckdb {

class ReadENATableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions;
		std::string result_type;
		std::string fields;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions, const std::string &result_type, const std::string &fields);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		unique_ptr<miint::ENAClient> client;
		// Rows reordered to match the schema column order
		std::vector<std::vector<std::string>> rows;
		size_t next_accession_idx;
		size_t row_offset;

		GlobalState(DatabaseInstance &db, const std::vector<std::string> &accessions, const std::string &result_type,
		            const std::string &fields, const std::vector<std::string> &schema_names);

		bool FetchNextAccession();

		idx_t MaxThreads() const override {
			return 1;
		}

	private:
		std::vector<std::string> accessions;
		std::string result_type;
		std::string fields;
		std::vector<std::string> schema_names; // Column order from user's fields param
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
