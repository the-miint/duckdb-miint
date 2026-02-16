#pragma once

#ifdef RYPE_ARROW

#include "rype.h"

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/function/table/arrow.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <string>
#include <vector>

namespace duckdb {

// ============================================================================
// Shared bind data for both extraction functions
// ============================================================================
struct RypeExtractData : public TableFunctionData {
	std::string sequence_table;
	std::string id_column;
	size_t k;
	size_t w;
	uint64_t salt;

	std::vector<std::string> names;
	std::vector<LogicalType> types;

	RypeExtractData() : id_column("read_id"), k(0), w(0), salt(0) {
	}
};

// ============================================================================
// Shared global state for both extraction functions
// ============================================================================
struct RypeExtractGlobalState : public GlobalTableFunctionState {
	// Original read_ids indexed by row number (0-based).
	std::vector<std::string> read_ids;

	// Arrow output stream from RYpe extraction.
	// OWNERSHIP HIERARCHY (destruction must be in reverse order):
	// 1. current_chunk (shared_ptr — may outlive gstate via Vector ArrowAuxiliaryData)
	// 2. arrow_table (holds pointers into output_schema)
	// 3. output_schema
	// 4. output_stream
	ArrowArrayStream output_stream;
	ArrowSchema output_schema;
	ArrowTableSchema arrow_table;

	// Current batch wrapped in shared_ptr for zero-copy lifetime management.
	// ArrowArrayScanState::owned_data references this so DuckDB Vectors can
	// point directly into Arrow buffers without copying.
	shared_ptr<ArrowArrayWrapper> current_chunk;

	idx_t batch_offset = 0;
	bool schema_initialized = false;
	bool done = false;

	idx_t MaxThreads() const override {
		return 1;
	}

	~RypeExtractGlobalState();
};

// ============================================================================
// Shared local state
// ============================================================================
struct RypeExtractLocalState : public LocalTableFunctionState {
	std::unordered_map<idx_t, unique_ptr<ArrowArrayScanState>> array_states;
	ClientContext &context;

	explicit RypeExtractLocalState(ClientContext &ctx) : context(ctx) {
	}
	~RypeExtractLocalState();

	ArrowArrayScanState &GetState(idx_t col_idx);
	void ResetStates();
};

// ============================================================================
// rype_extract_minimizer_set table function
// ============================================================================
class RypeExtractMinimizerSetTableFunction {
public:
	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

// ============================================================================
// rype_extract_strand_minimizers table function
// ============================================================================
class RypeExtractStrandMinimizersTableFunction {
public:
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

#endif // RYPE_ARROW
