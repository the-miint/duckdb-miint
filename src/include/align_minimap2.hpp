#pragma once
#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <mutex>
#include <vector>

namespace duckdb {

class AlignMinimap2TableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string subject_table;           // OPTIONAL (either this or index_path required)
		std::string index_path;              // OPTIONAL: path to .mmi file
		bool per_subject_database;
		miint::Minimap2Config config;
		SequenceTableSchema query_schema;
		std::vector<miint::AlignmentSubject> subjects; // Pre-loaded at bind time (empty if using index_path)

		// Helper to check if using pre-built index
		bool using_prebuilt_index() const {
			return !index_path.empty();
		}

		// Output schema (shared with align_minimap2_sharded)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : per_subject_database(false),
		      names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes()) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		idx_t current_query_offset;
		idx_t current_subject_idx;    // For per_subject mode
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset;
		bool done;

		// For per-subject mode: store all queries in memory to avoid re-reading
		miint::SequenceRecordBatch all_queries;
		bool queries_loaded;

		idx_t MaxThreads() const override {
			// Single-threaded for now - minimap2 has internal threading
			return 1;
		}

		GlobalState() : current_query_offset(0), current_subject_idx(0), buffer_offset(0), done(false),
		                queries_loaded(false) {
		}
	};

	struct LocalState : public LocalTableFunctionState {
		LocalState() {
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
