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
#include <atomic>
#include <chrono>
#include <mutex>
#include <vector>

namespace duckdb {

class AlignMinimap2TableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string subject_table; // OPTIONAL (either this or index_path required)
		std::string index_path;    // OPTIONAL: path to .mmi file
		bool per_subject_database;
		miint::Minimap2Config config;
		SequenceTableSchema query_schema;
		std::vector<miint::AlignmentSubject> subjects; // Pre-loaded at bind time (empty if using index_path)

		// Subject-side id type. Drives output `reference` + `mate_reference`.
		// Defaults to INVALID so a Bind path that forgets to set it fails loud
		// at the helper dispatch in id_column_utils.hpp rather than silently
		// producing wrong-typed output. Bind sets it explicitly: subject_table
		// mode pulls from the subject schema's id_type; index_path mode sets
		// VARCHAR (the .mmi file stores subject names as opaque bytes).
		LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);

		// Helper to check if using pre-built index
		bool using_prebuilt_index() const {
			return !index_path.empty();
		}

		// Output schema. Names are constant; types are mutated by Bind once
		// the query and subject id types are known.
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		bool debug = false;

		// types is rebuilt by Bind once the actual query/subject id types are
		// known; the placeholder VARCHAR/VARCHAR here is never observed.
		Data()
		    : per_subject_database(false), names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	// Standard mode state: multi-threaded, shared index, lazy sub-batch streaming
	struct StandardModeState {
		std::shared_ptr<miint::SharedMinimap2Index> shared_index;
		std::unique_ptr<QuerySequenceStream> query_stream; // Lazy streaming reader
	};

	// Per-subject mode state: single-threaded, builds index per subject
	struct PerSubjectModeState {
		std::mutex lock;
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		idx_t current_subject_idx = 0;
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
		bool done = false;
		miint::SequenceRecordBatch all_queries;
		bool queries_loaded = false;
	};

	struct GlobalState : public GlobalTableFunctionState {
		bool per_subject_mode = false;
		idx_t num_threads = 1;
		bool debug = false;
		std::chrono::steady_clock::time_point start_time;
		std::atomic<idx_t> init_local_count {0};

		// Exactly one of these is populated based on per_subject_mode
		std::unique_ptr<StandardModeState> standard;
		std::unique_ptr<PerSubjectModeState> per_subject;

		idx_t MaxThreads() const override {
			return num_threads;
		}
	};

	struct LocalState : public LocalTableFunctionState {
		// Per-thread aligner (standard mode)
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		// Per-thread output buffer
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
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
