#pragma once

#include "SortMeRNAAligner.hpp"
#include "align_common.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <memory>
#include <mutex>
#include <vector>

namespace duckdb {

// SAM-schema counterpart of align_sortmerna_rrna. The streaming skeleton is
// identical; only the output projection changes — see OutputSortMeRNASamBatch
// in sortmerna_result_utils.hpp.
class AlignSortMeRNATableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::vector<std::string> ref_paths;
		SequenceTableSchema query_schema;
		miint::SortMeRNAConfig config;

		std::vector<std::string> names;
		std::vector<LogicalType> types;

		// SortMeRNA is VARCHAR-only for read_id / reference / mate_reference
		// in PR 1. BIGINT support requires updating the sortmerna emit helper
		// (out of scope here); see findings doc.
		Data()
		    : names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::unique_ptr<miint::SortMeRNAAligner> aligner;
		std::unique_ptr<QuerySequenceStream> query_stream;
		std::mutex lock;
		miint::SortMeRNAResultBatch result_buffer;
		idx_t buffer_offset = 0;

		// Sortmerna's g_run_mutex serializes calls process-wide; running with
		// MaxThreads() > 1 would just queue calls behind the mutex. Run in a
		// single DuckDB thread and let sortmerna's internal pool parallelize.
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
