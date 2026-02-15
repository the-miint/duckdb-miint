#pragma once

#ifdef RYPE_ARROW

#include "rype.h"

#include "duckdb/common/arrow/arrow.hpp"
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

class RypeClassifyTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string index_path;
		std::string sequence_table;
		std::string id_column;
		double threshold;
		std::string negative_index_path;
		bool has_sequence2 = false; // Cached from ValidateSequenceTable

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : id_column("read_id"), threshold(0.1), names({"read_id", "bucket_id", "bucket_name", "score"}),
		      types({LogicalType::VARCHAR,  // read_id (original identifier)
		             LogicalType::UINTEGER, // bucket_id (UInt32)
		             LogicalType::VARCHAR,  // bucket_name
		             LogicalType::DOUBLE})  // score (Float64)
		{
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		// RYpe index (opaque pointer from C API)
		RypeIndex *index = nullptr;

		// Optional negative set for filtering
		RypeNegativeSet *negative_set = nullptr;

		// Original read_ids indexed by row number (0-based).
		// RYpe receives row indices as query_id, so read_ids[query_id] gives the original identifier.
		std::vector<std::string> read_ids;

		// Arrow output stream from RYpe.
		// OWNERSHIP HIERARCHY (destruction must be in reverse order):
		// 1. current_batch - obtained via get_next(), release before getting next or on destruction
		// 2. arrow_table - holds pointers INTO output_schema, clear before releasing schema
		// 3. output_schema - obtained via get_schema(), separately owned copy, release on destruction
		// 4. output_stream - returned by rype_classify_arrow(), release on destruction
		ArrowArrayStream output_stream;
		ArrowSchema output_schema;
		ArrowTableSchema arrow_table;
		ArrowArray current_batch;

		idx_t batch_offset = 0;
		bool schema_initialized = false;
		bool done = false;

		// No mutex needed - MaxThreads() returns 1, enforcing single-threaded execution.

		idx_t MaxThreads() const override {
			return 1;
		}

		~GlobalState();
	};

	struct LocalState : public LocalTableFunctionState {
		std::unordered_map<idx_t, unique_ptr<ArrowArrayScanState>> array_states;
		ClientContext &context;

		explicit LocalState(ClientContext &ctx) : context(ctx) {
		}
		~LocalState();

		ArrowArrayScanState &GetState(idx_t col_idx);
		void ResetStates();
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

#endif // RYPE_ARROW
