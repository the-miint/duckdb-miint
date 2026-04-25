#pragma once

#include "VsearchChimeraWrapper.hpp"
#include "per_sample_table_function.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>
#include <vector>

namespace duckdb {

class UchimeDenovoTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string input_table;
		miint::UchimeParams params;

		// Column-name overrides for the input table. Defaults preserve historical
		// names so chains like deblur(count_col := 'abundance') →
		// detect_chimera_uchime_denovo(count_col := 'abundance') work without rename.
		std::string id_col = "read_id";
		std::string sequence_col = "sequence1";
		std::string count_col = "size";

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		bool has_sample_id = false;
		PerSampleBindInfo sample_info;
	};

	struct GlobalState : public PerSampleGlobalState {
		// Non-sample path: precomputed sequence set + incremental processing state.
		// Single-thread drain (max_threads=1); wrapper holds k-mer index as it grows.
		miint::VsearchChimeraWrapper wrapper;

		std::vector<std::string> labels;
		std::vector<std::string> sequences;
		std::vector<int64_t> sizes;

		idx_t input_offset = 0;
		idx_t indexed_count = 0;
		std::vector<miint::UchimeResult> results;
		idx_t result_offset = 0;
	};

	struct LocalState : public LocalTableFunctionState {
		unique_ptr<Connection> conn; // sample_id path only

		// Current sample's pre-computed result buffer (drained over multiple Execute calls).
		std::vector<miint::UchimeResult> results;
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
