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
#include "duckdb/parallel/task_scheduler.hpp"
#include "minimap2_part_cursor.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

// Information about a single shard
struct ShardInfo {
	std::string name;       // e.g., "shard_001"
	std::string index_path; // e.g., "/path/shards/shard_001.mmi"
	idx_t read_count;       // Number of reads for this shard (for priority ordering)
};

// The TEMP snapshots every shard stream reads from: the query relation and the
// routing relation, each read exactly once (#229) and then replayed once per
// shard and index part.
//
// Held by shared_ptr, and pinned by every QuerySequenceStream opened over them,
// because the DROP must wait for the last open stream. DuckDB does not guarantee
// that local states die before the global one, and a LocalState keeps its
// ActiveShard — with that shard's live StreamQueryResult — whenever the scan ends
// without draining (a LIMIT above the function is the ordinary case). Dropping a
// table out from under a live stream is not something the reader survives, so
// ownership decides when it happens: the stream that reads these tables holds
// them (see QuerySequenceStream's `source_keepalive`), rather than any holder
// maintaining a destruction order it cannot enforce.
//
// The connection is part of the handle for the same reason it was part of the
// state before: the tables were created on a connection that inherits the
// caller's TEMP catalog, so they live in the user's session, and a missed drop
// leaves them visible in SHOW TABLES rather than dying with us.
struct QuerySnapshots {
	std::unique_ptr<Connection> conn;
	std::string query_reads;   // unquoted; empty => nothing to drop
	std::string read_to_shard; // unquoted; empty => nothing to drop

	~QuerySnapshots() {
		if (!conn) {
			return;
		}
		// DropHelperTempRelation is no-op on an empty name and never propagates.
		DropHelperTempRelation(*conn, KeywordHelper::WriteOptionallyQuoted(read_to_shard));
		DropHelperTempRelation(*conn, KeywordHelper::WriteOptionallyQuoted(query_reads));
	}
};

// A shard that is currently being processed by one or more threads.
// `parts` owns the shard's .mmi and, for a multi-part shard, walks its workers
// through the parts one at a time (see Minimap2PartCursor). The shard's reads
// are never held whole: `stream` streams them from the query snapshot in
// sub-batches, and every part gets a fresh stream, opened by the leader in the
// cursor's `prepare` step and swapped in under the cursor lock in `publish`.
// `stream` is a shared_ptr so a thread that picked it up just before a swap can
// keep draining it safely. Worker tracking uses atomics and never holds the
// global lock.
struct ActiveShard {
	idx_t shard_idx;                                  // Index into Data::shards
	std::shared_ptr<QuerySequenceStream> stream;      // Current part's reads; read via parts->WithCurrentPart
	std::unique_ptr<miint::Minimap2PartCursor> parts; // Index parts; set once ready
	std::atomic<idx_t> active_workers {0};            // Threads currently on this shard
	std::atomic<bool> exhausted {false};              // Set when no more batches to read
	std::atomic<bool> ready {false};                  // Set when index is loaded and stream opened
	std::atomic<bool> progress_reconciled {false};    // Estimated read count corrected to the streamed one
	// Progress-only (read/written only when GlobalState::progress is true).
	std::atomic<idx_t> alignments_emitted {0};        // Mapped alignments produced for this shard
	std::chrono::steady_clock::time_point start_time; // Stamped when the shard becomes ready
};

class AlignMinimap2ShardedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string shard_directory;
		std::string read_to_shard_table;
		SequenceTableSchema query_schema;
		miint::Minimap2Config config;
		std::vector<ShardInfo> shards; // Sorted by read_count DESC (largest first)
		idx_t max_threads_per_shard = 4;
		bool debug = false;
		bool progress = false;
		bool include_shard_name = false;

		// Subject-side id type. Sharded mode always loads prebuilt .mmi indexes
		// whose subject names are opaque bytes — same contract as align_minimap2
		// `index_path` mode — so this defaults to VARCHAR once Bind runs. The
		// INVALID sentinel here mirrors align_minimap2.hpp's fail-loud default:
		// any path that forgets to set this triggers a clear error at the helper
		// dispatch in id_column_utils.hpp.
		LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);

		// Output schema (shared with align_minimap2). `names` are constant;
		// `types` is rebuilt by Bind once query_schema.id_type is known.
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		// types is rebuilt by Bind once the actual query/subject id types are
		// known; the placeholder VARCHAR/VARCHAR here is never observed.
		Data()
		    : names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes(LogicalType::VARCHAR, LogicalType::VARCHAR)) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		std::condition_variable cv;
		idx_t next_shard_idx = 0;
		idx_t shard_count = 0;
		idx_t max_threads_per_shard = 4;
		idx_t max_active_shards = 1; // ceil(db_threads / max_threads_per_shard)
		bool debug = false;
		bool progress = false;
		std::chrono::steady_clock::time_point start_time;
		std::vector<std::shared_ptr<ActiveShard>> active_shards;
		std::atomic<idx_t> total_associations {0};
		std::atomic<idx_t> associations_processed {0};

		// TEMP snapshots of the query and routing relations, which every shard's
		// stream reads from (BuildShardReadsSelect). Taken even for a single shard:
		// a multi-part shard replays its reads once per part, and re-reading the
		// user's relations instead would silently drop rows for any relation not
		// stable across re-evaluation (#229 — see
		// docs/internals/reading-tables-views.md § "Read the relation ONCE").
		std::shared_ptr<QuerySnapshots> snapshots;

		idx_t MaxThreads() const override {
			return max_active_shards * max_threads_per_shard;
		}

		GlobalState() = default;

		// No destructor: the snapshots are dropped by ~QuerySnapshots once the last
		// shard still streaming them is gone, which may be after this state dies.
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		std::shared_ptr<ActiveShard> current_active_shard;
		bool has_shard = false;
		miint::Minimap2PartCursor::Attachment part; // which part of current_active_shard aligner is on
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;
		std::string current_shard_name;

		LocalState() = default;
	};

	static unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input,
	                                     vector<LogicalType> &return_types, vector<std::string> &names);

	static unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input);

	static unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
	                                                     GlobalTableFunctionState *global_state);

	static void Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output);

	static double Progress(ClientContext &context, const FunctionData *bind_data,
	                       const GlobalTableFunctionState *global_state);

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);

private:
	// Claim work: join an existing active shard or start a new one.
	// Returns the ActiveShard to work on, or nullptr if no more work.
	// Index loading happens outside the lock.
	static std::shared_ptr<ActiveShard> ClaimWork(ClientContext &context, GlobalState &gstate, const Data &bind_data,
	                                              LocalState &lstate);

	// Release work: detach from current shard, clean up if last worker on exhausted shard.
	static void ReleaseWork(GlobalState &gstate, LocalState &lstate);
};

} // namespace duckdb
