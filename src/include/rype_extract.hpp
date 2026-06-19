#pragma once

#include "rype.h"
#include "rype_common.hpp"

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
	LogicalType id_type = LogicalType(LogicalTypeId::INVALID); // Storage type of id_column; drives read_id output

	std::vector<std::string> names;
	std::vector<LogicalType> types;

	RypeExtractData() : id_column("read_id"), k(0), w(0), salt(6148914691236517205ULL) {
	}
};

// ============================================================================
// Shared global state for both extraction functions
// ============================================================================
struct RypeExtractGlobalState : public GlobalTableFunctionState {
	// Original read_ids indexed by row number (0-based).
	std::vector<std::string> read_ids;

	// Sub-connection used to query the sequence table. Must outlive the RYpe
	// output_stream because RYpe lazily consumes the input Arrow stream (backed
	// by a ResultArrowArrayStreamWrapper whose QueryResult holds a non-owning
	// optional_ptr<ClientContext> into this connection). Destroyed explicitly
	// in the destructor AFTER releasing RYpe streams.
	unique_ptr<Connection> input_connection;

	// Name of the per-call TEMP table that materializes (id, read_id, sequence1)
	// once on input_connection. Populated by BuildExtractionInputStream; the
	// destructor drops the table before tearing down input_connection.
	std::string tmp_table_name;

	// Arrow output stream from RYpe extraction.
	// OWNERSHIP HIERARCHY (destruction must be in reverse order):
	// 1. current_chunk (shared_ptr — may outlive gstate via Vector ArrowAuxiliaryData)
	// 2. arrow_table (holds pointers into output_schema)
	// 3. output_schema
	// 4. output_stream - until released, RYpe may still pull from the input Arrow stream
	//                    wrapper, which holds a non-owning ClientContext pointer into
	//                    input_connection — so this must be released before the DROP and
	//                    before input_connection.reset().
	// 5. tmp_table - DROPPED on input_connection (which owns the per-call TEMP). Must run
	//                while input_connection is alive AND after step 4 returns, since
	//                step 4 can re-enter the connection through RYpe's stream callbacks.
	// 6. input_connection - reset last; outlives steps 4 and 5.
	ArrowArrayStream output_stream;
	ArrowSchema output_schema;
	ArrowTableSchema arrow_table;

	// Current batch wrapped in shared_ptr for zero-copy lifetime management.
	// ArrowArrayScanState::owned_data references this so DuckDB Vectors can
	// point directly into Arrow buffers without copying.
	shared_ptr<ArrowArrayWrapper> current_chunk;

	idx_t batch_offset = 0;
	bool done = false;

	idx_t MaxThreads() const override {
		return 1;
	}

	~RypeExtractGlobalState();
};

// Use shared RypeArrowLocalState — identical across all RYpe table functions
using RypeExtractLocalState = RypeArrowLocalState;

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
