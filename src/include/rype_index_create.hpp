#pragma once

#include "rype.h"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {

// Build a RYpe .ryxdi index from an at-rest chunked sequence table (Qiita
// reference_sequence_chunks schema: feature_idx BIGINT, chunk_index INTEGER,
// chunk_data VARCHAR) plus an optional feature->bucket mapping. Wraps the rype
// FFI rype_index_build_from_arrow, feeding it two DuckDB query results as Arrow
// streams. Returns a single status row.
class RypeIndexCreateTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string chunk_table;
		std::string output_path;
		// Empty => synthesize a single 'unnamed-bucket' over all features.
		std::string mapping_table;
		int32_t k;
		int32_t w;
		uint64_t salt;
		bool orient;
		int64_t max_memory; // bytes; 0 = auto-detect

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : k(64), w(50), salt(0x5555555555555555ULL), orient(true), max_memory(0),
		      names({"output_path", "k", "w", "status"}),
		      types({LogicalType::VARCHAR, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::VARCHAR}) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		// The build is performed in InitGlobal. Execute only emits the one status
		// row; `done` guards against a second emit.
		bool done = false;

		// Sub-connection feeding both Arrow streams. The (small) bucket mapping is
		// materialized — RYpe reads it eagerly into memory anyway — while the (large)
		// sequence/chunk data is streamed via a SendQuery cursor and fetched lazily.
		// A materialized result does not occupy the connection's single
		// StreamQueryResult slot, so one connection serves both. It is held here so
		// it outlives RYpe's lazy consumption of the chunk cursor during the build.
		unique_ptr<Connection> input_connection;

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
