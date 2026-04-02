#pragma once

#include "VsearchChimeraWrapper.hpp"
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

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	struct GlobalState : public GlobalTableFunctionState {
		// Vsearch wrapper with all sequences loaded and DUST-masked.
		// K-mer index starts empty — non-chimeras are indexed incrementally.
		miint::VsearchChimeraWrapper wrapper;

		std::vector<std::string> labels;
		std::vector<std::string> sequences;
		std::vector<int64_t> sizes;

		// Incremental processing state. Each Execute() call processes a batch of
		// input sequences and buffers the results, allowing DuckDB to check for
		// cancellation between calls.
		idx_t input_offset = 0;
		idx_t indexed_count = 0; // Number of sequences indexed so far
		std::vector<miint::UchimeResult> results;
		idx_t result_offset = 0;

		idx_t MaxThreads() const override {
			return 1; // De novo mode is inherently sequential
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
