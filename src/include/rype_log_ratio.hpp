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

class RypeLogRatioTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string numerator_path;
		std::string denominator_path;
		std::string sequence_table;
		std::string id_column;
		double skip_threshold;
		bool has_sequence2 = false;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : id_column("read_id"),
		      // Default 0.5: reads with numerator score >= 50% skip denominator classification
		      // and get +inf log_ratio immediately (fast-path). Set to 0 or negative to disable.
		      skip_threshold(0.5), names({"read_id", "log_ratio", "fast_path"}),
		      types({LogicalType::VARCHAR,  // read_id (original identifier)
		             LogicalType::DOUBLE,   // log_ratio (Float64)
		             LogicalType::INTEGER}) // fast_path (Int32)
		{
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		// RYpe indices (opaque pointers from C API)
		RypeIndex *numerator_index = nullptr;
		RypeIndex *denominator_index = nullptr;

		// Original read_ids indexed by row number (0-based).
		std::vector<std::string> read_ids;

		// Sub-connection used to query the sequence table. Must outlive the RYpe
		// output_stream because RYpe lazily consumes the input Arrow stream (backed
		// by a ResultArrowArrayStreamWrapper whose QueryResult holds a non-owning
		// optional_ptr<ClientContext> into this connection). Destroyed explicitly
		// in the destructor AFTER releasing RYpe streams.
		unique_ptr<Connection> input_connection;

		// Name of the per-call TEMP table that materializes (id, read_id, sequence1,
		// sequence2?) once on input_connection. Populated by InitGlobal; the
		// destructor drops the table before tearing down input_connection.
		std::string tmp_table_name;

		// Arrow output stream from RYpe.
		// OWNERSHIP HIERARCHY (destruction must be in reverse order):
		// 1. current_chunk (shared_ptr — may outlive gstate via Vector ArrowAuxiliaryData)
		// 2. arrow_table - holds pointers INTO output_schema, clear before releasing schema
		// 3. output_schema - obtained via get_schema(), separately owned copy, release on destruction
		// 4. output_stream - returned by rype_classify_arrow_log_ratio(), release on destruction.
		//                    Until released, RYpe may still pull from the input Arrow stream
		//                    wrapper, which holds a non-owning ClientContext pointer into
		//                    input_connection — so this must be released before the DROP and
		//                    before input_connection.reset().
		// 4a. numerator_index / denominator_index - rype_index_free; safe any time after
		//                    output_stream.release returns (RYpe is done with them).
		// 5. tmp_table - DROPPED on input_connection (which owns the per-call TEMP). Must run
		//                while input_connection is alive AND after step 4 returns, since
		//                step 4 can re-enter the connection through RYpe's stream callbacks.
		// 6. input_connection - reset last; outlives steps 4 and 5.
		ArrowArrayStream output_stream;
		ArrowSchema output_schema;
		ArrowTableSchema arrow_table;
		shared_ptr<ArrowArrayWrapper> current_chunk;

		idx_t batch_offset = 0;
		bool done = false;

		idx_t MaxThreads() const override {
			return 1;
		}

		~GlobalState();
	};

	// Use shared RypeArrowLocalState — identical across all RYpe table functions
	using LocalState = RypeArrowLocalState;

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
