#pragma once

#include "ChimeraDetector.hpp"
#include "WFA2Aligner.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <atomic>
#include <mutex>
#include <string>
#include <vector>

namespace duckdb {

class UchimeRefTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string ref_table;
		miint::UchimeParams params;

		// Validated schemas from Bind() (confirms read_id + sequence1 exist)
		SequenceTableSchema query_schema;
		SequenceTableSchema ref_schema;

		std::vector<std::string> names;
		std::vector<LogicalType> types;
	};

	struct GlobalState : public GlobalTableFunctionState {
		// Reference DB loaded at init, shared read-only across threads.
		miint::ChimeraDetector detector;

		// Lazy streaming reader for query sequences (thread-safe).
		std::unique_ptr<QuerySequenceStream> query_stream;

		idx_t num_threads = 1;

		idx_t MaxThreads() const override {
			return num_threads;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		// Per-thread WFA2 aligner with UCHIME-equivalent penalties.
		miint::WFA2Aligner aligner;

		// Buffered results from processing a sub-batch of queries.
		std::vector<miint::UchimeResult> result_buffer;
		idx_t buffer_offset = 0;

		LocalState()
		    : aligner(miint::UCHIME_WFA2_MISMATCH, miint::UCHIME_WFA2_GAP_OPEN, miint::UCHIME_WFA2_GAP_EXTEND) {
		}
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
