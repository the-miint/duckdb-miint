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

		// Sub-connection used to query the sequence table. Must outlive the RYpe
		// output_stream because RYpe lazily consumes the input Arrow stream (backed
		// by a ResultArrowArrayStreamWrapper whose QueryResult holds a non-owning
		// optional_ptr<ClientContext> into this connection). Destroyed explicitly
		// in the destructor AFTER releasing RYpe streams.
		unique_ptr<Connection> input_connection;

		// Arrow output stream from RYpe.
		// OWNERSHIP HIERARCHY (destruction must be in reverse order):
		// 1. current_chunk (shared_ptr — may outlive gstate via Vector ArrowAuxiliaryData)
		// 2. arrow_table - holds pointers INTO output_schema, clear before releasing schema
		// 3. output_schema - obtained via get_schema(), separately owned copy, release on destruction
		// 4. output_stream - returned by rype_classify_arrow(), release on destruction
		// 5. input_connection - must outlive output_stream (RYpe holds ref to input Arrow stream)
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
