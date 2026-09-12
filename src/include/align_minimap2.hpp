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
#include <atomic>
#include <chrono>
#include <condition_variable>
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
	// index_reader is null except when using_prebuilt_index() opened a
	// multi-part .mmi. In that case shared_index / query_stream hold GlobalState's
	// reference to the CURRENT part — the leader in AdvanceMinimap2Part drops it
	// before loading the next one. That bounds peak memory to roughly the
	// largest single part plus fixed overhead, NOT one part at a time exactly: a
	// worker still inside align() on the outgoing part keeps its own shared_ptr
	// copy alive while the next part loads, so the true transient peak is
	// (outgoing part still in flight + incoming part), and steady-state peak
	// tracks whichever part is biggest, not a fixed "2-3GB"-style constant —
	// measured at 5.28GB peak RSS for an evenly-split 4-part ~9GB human genome
	// index. part_lock/part_cv/part_generation/advancing coordinate worker
	// threads through the transition from one part to the next; see
	// AdvanceMinimap2Part in align_minimap2.cpp. query_stream is a shared_ptr
	// (not unique_ptr) so a thread that captured it just before a swap can keep
	// draining it safely instead of racing the object's destruction.
	struct StandardModeState {
		std::shared_ptr<miint::SharedMinimap2Index> shared_index;
		std::shared_ptr<QuerySequenceStream> query_stream;

		std::unique_ptr<miint::Minimap2IndexReader> index_reader; // multi-part prebuilt index only
		std::mutex part_lock;
		std::condition_variable part_cv;
		idx_t part_generation = 0; // bumped every time shared_index/query_stream are swapped to a new part
		bool advancing = false;    // true while one thread is loading the next part
		bool parts_exhausted = false;
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
		// Declared BEFORE `standard` so members destruct in the right order:
		// `standard->query_stream` can hold a QuerySequenceStream built with the
		// Connection& overload, which stores a raw `Connection*` into
		// snapshot_conn rather than owning it. Members destruct in reverse
		// declaration order, so snapshot_conn here means `standard` (and its
		// query_stream) is torn down FIRST, before the connection it points into
		// goes away.
		std::unique_ptr<Connection> snapshot_conn;
		std::string query_snapshot; // unquoted; empty => no snapshot to drop

		// Exactly one of these is populated based on per_subject_mode
		std::unique_ptr<StandardModeState> standard;
		std::unique_ptr<PerSubjectModeState> per_subject;

		idx_t MaxThreads() const override {
			return num_threads;
		}

		GlobalState() = default;
		~GlobalState() override {
			// Release state that may hold a live stream over the snapshot table
			// BEFORE dropping the table itself — the destructor body runs before
			// member destruction, so without this explicit reset the DROP below
			// would run while `standard->query_stream` (an early-terminated query,
			// e.g. LIMIT, never reads it to exhaustion) still has an open
			// StreamQueryResult over that same table.
			standard.reset();
			per_subject.reset();
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
		// Multi-part prebuilt index only: which StandardModeState::part_generation
		// this thread's aligner is currently attached to, and whether it has
		// attached at all yet. InitLocal deliberately does NOT attach eagerly in
		// this mode (unlike the single-part/subject_table paths): DuckDB may
		// initialize more LocalStates than ever receive real work (observed with
		// a handful of query rows and MaxThreads()=12 — most threads got exactly
		// one empty fetch and never ran again), and an eager attach would leave
		// each such idle thread holding a live shared_ptr to whatever part was
		// current at InitLocal time for the rest of the query — keeping that part
		// resident alongside every later one and defeating the one-part-at-a-time
		// memory bound streaming exists for. Attaching lazily on first real use
		// means a thread that never does real work never holds a part reference
		// at all. Always false/0 (and never consulted) outside multi-part mode.
		bool attached_to_part = false;
		idx_t part_generation = 0;
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
