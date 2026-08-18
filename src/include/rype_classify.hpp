#pragma once

#include "rype.h"
#include "rype_common.hpp"
#include "rype_id_map.hpp"
#include "rype_input_stream.hpp"

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

class RypeClassifyTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string index_path;
		std::string sequence_table;
		std::string id_column;
		double threshold;
		std::string negative_index_path;
		// Byte budget handed to RYpe's batch sizing; 0 = derive one (ResolveRypeMemoryBudget).
		int64_t max_memory = 0;
		// Report the chosen batch size and memory estimates via miint_warnings() (#204).
		bool debug = false;
		bool has_sequence2 = false;                                // Cached from ValidateSequenceTable
		LogicalType id_type = LogicalType(LogicalTypeId::INVALID); // Storage type of id_column; drives read_id output

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

		// Caller identifiers indexed by the surrogate id miint hands RYpe, which
		// RYpe echoes back as query_id. Filled incrementally by RypeInputStream as
		// RYpe pulls input batches, so it must outlive the input stream — see the
		// ownership hierarchy below. Constructed in InitGlobal once the id type is
		// known.
		unique_ptr<RypeIdMap> id_map;

		// Sub-connection used to query the sequence table. Must outlive the RYpe
		// output_stream because RYpe lazily consumes the input Arrow stream (whose
		// streaming QueryResult holds a non-owning optional_ptr<ClientContext> into
		// this connection). Destroyed explicitly in the destructor AFTER releasing
		// RYpe streams.
		unique_ptr<Connection> input_connection;

		// The Arrow input stream handed to RYpe. Owned here, not by RYpe — see
		// the ownership note in rype_input_stream.hpp for why the usual
		// transfer-to-callee arrangement is unsafe against a -1 return.
		unique_ptr<RypeInputStream> input_stream;

		// Arrow output stream from RYpe.
		// OWNERSHIP HIERARCHY (destruction must be in reverse order):
		// 1. current_chunk (shared_ptr — may outlive gstate via Vector ArrowAuxiliaryData)
		// 2. arrow_table - holds pointers INTO output_schema, clear before releasing schema
		// 3. output_schema - obtained via get_schema(), separately owned copy, release on destruction
		// 4. output_stream - returned by rype_classify_arrow_ex(), release on destruction.
		//                    Releasing it makes RYpe drop its reader, which calls the input
		//                    stream's release callback and stops any further writes to id_map.
		// 4a. negative_set / index - rype_negative_set_free / rype_index_free; safe any time
		//                    after output_stream.release returns (RYpe is done with them).
		// 5. input_stream - after step 4, because RYpe consumes it lazily and may still be
		//                    pulling from it until then.
		// 6. input_connection / id_map - torn down last. input_stream holds a streaming
		//                    QueryResult carrying a non-owning ClientContext pointer into
		//                    input_connection, so the connection must outlive step 5.
		ArrowArrayStream output_stream;
		ArrowSchema output_schema;
		ArrowTableSchema arrow_table;
		shared_ptr<ArrowArrayWrapper> current_chunk;

		idx_t batch_offset = 0;
		bool done = false;

		// No mutex needed - MaxThreads() returns 1, enforcing single-threaded execution.

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
