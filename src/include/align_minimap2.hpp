#pragma once
#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_common.hpp"
#include "catalog_utils.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/typedefs.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/function.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "minimap2_part_cursor.hpp"
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

	// Standard mode state: multi-threaded, shared index, lazy sub-batch streaming.
	//
	// `parts` is null except when using_prebuilt_index() opened a multi-part
	// .mmi. In that case the cursor owns the CURRENT part and coordinates worker
	// threads through part transitions (memory bound, lazy attach and the
	// leader/waiter protocol are documented on Minimap2PartCursor), shared_index
	// is unused, and query_stream is the replay stream over the query snapshot
	// for that part — swapped by the leader in the cursor's `publish` step.
	// query_stream is a shared_ptr (not unique_ptr) so a thread that captured it
	// just before a swap can keep draining it safely instead of racing the
	// object's destruction.
	struct StandardModeState {
		std::shared_ptr<miint::SharedMinimap2Index> shared_index;
		std::shared_ptr<QuerySequenceStream> query_stream;
		std::unique_ptr<miint::Minimap2PartCursor> parts; // multi-part prebuilt index only
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

		// Multi-part prebuilt index only: TEMP snapshot of the query relation so
		// it can be replayed once per part (#229 — see
		// docs/internals/reading-tables-views.md § "Read the relation ONCE").
		// Empty/null for single-part indexes and subject_table mode, which keep
		// the original single streaming pass with no snapshot at all.
		//
		// `standard->query_stream` can hold a QuerySequenceStream built with the
		// Connection& overload, which points at this connection without owning
		// it, so it must die first — ~GlobalState resets `standard` explicitly
		// (it has to, for the DROP); the declaration order here only mirrors that.
		std::unique_ptr<Connection> snapshot_conn;
		std::string query_snapshot; // unquoted; empty => no snapshot to drop

		// Opens a replay stream over the snapshot. Every part's stream must be
		// built identically — one built differently mid-scan would change the
		// projection or sub-batching partway through — so both the first part
		// (InitGlobal) and every later one (the cursor's `prepare`) come through
		// here rather than repeating the construction.
		std::shared_ptr<QuerySequenceStream> OpenSnapshotStream(const SequenceTableSchema &schema) {
			return std::make_shared<QuerySequenceStream>(*snapshot_conn, query_snapshot, schema);
		}

		// Exactly one of these is populated based on per_subject_mode
		std::unique_ptr<StandardModeState> standard;
		std::unique_ptr<PerSubjectModeState> per_subject;

		idx_t MaxThreads() const override {
			return num_threads;
		}

		~GlobalState() override {
			// Release state that may hold a live stream over the snapshot table
			// BEFORE dropping the table itself — the destructor body runs before
			// member destruction, so without this explicit reset the DROP below
			// would run while `standard->query_stream` (an early-terminated query,
			// e.g. LIMIT, never reads it to exhaustion) still has an open
			// StreamQueryResult over that same table.
			standard.reset();
			if (snapshot_conn) {
				DropHelperTempRelation(*snapshot_conn, KeywordHelper::WriteOptionallyQuoted(query_snapshot));
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		// Per-thread aligner (standard mode)
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		// Per-thread output buffer
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
		// Multi-part prebuilt index only: which part this thread's aligner is
		// attached to. InitLocal deliberately does not attach in this mode — see
		// Minimap2PartCursor on lazy attach. Unused otherwise.
		miint::Minimap2PartCursor::Attachment part;
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
