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

class ReadENAAttributesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::vector<std::string> accessions;
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data(std::vector<std::string> accessions);
	};

	struct GlobalState : public GlobalTableFunctionState {
		mutex lock;
		unique_ptr<miint::ENAClient> client;
		std::vector<miint::SampleAttribute> attributes;
		size_t row_offset;
		bool fetched;

		GlobalState(DatabaseInstance &db);

		void FetchAttributes(const std::vector<std::string> &accessions);

		idx_t MaxThreads() const override {
			return 1;
		}

	private:
		std::vector<std::string> ResolveSampleAccessions(const std::vector<std::string> &accessions);
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
