#pragma once

#include "VsearchClusterWrapper.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <memory>
#include <string>
#include <vector>

namespace duckdb {

class ClusterSequencesTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string input_table;
		miint::ClusterParams params;

		// Storage type of the input read_id column (VARCHAR, BIGINT, or UUID),
		// captured at bind. Both output id columns (read_id and centroid_id, which
		// is itself one of the input read_ids) mirror this type. Defaults to
		// INVALID so any path that reaches the egress dispatcher without Bind
		// having set it fails loud rather than misreading the carrier.
		LogicalType id_type = LogicalType(LogicalTypeId::INVALID);

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	struct GlobalState : public GlobalTableFunctionState {
		// All cluster results, computed in InitGlobal.
		std::vector<miint::ClusterResult> results;
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
