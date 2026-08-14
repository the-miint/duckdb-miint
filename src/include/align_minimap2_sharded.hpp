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

// A shard that is currently being processed by one or more threads.
// The shared index is immutable after construction; atomic counters
// coordinate batch claiming and worker tracking without holding the global lock.
struct ActiveShard {
	idx_t shard_idx;                                   // Index into Data::shards
	idx_t batch_size;                                  // Per-shard batch size
	miint::SequenceRecordBatch shard_sequences;        // Pre-fetched sequences for this shard
	std::shared_ptr<miint::SharedMinimap2Index> index; // Shared index, immutable after construction
	std::atomic<idx_t> next_batch_offset {0};          // Threads atomically claim ranges into shard_sequences
	std::atomic<idx_t> active_workers {0};             // Threads currently on this shard
	std::atomic<bool> exhausted {false};               // Set when no more batches to read
	std::atomic<bool> ready {false};                   // Set when index is loaded and IDs materialized
	// Progress-only (read/written only when GlobalState::progress is true).
	std::atomic<idx_t> alignments_emitted {0};        // Mapped alignments produced for this shard
	idx_t total_reads = 0;                            // Reads pre-fetched for this shard
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

		// What ClaimWork reads a shard's sequences FROM: either a one-pass TEMP
		// snapshot of the shard-assigned reads (multi-shard, #229) or an inline
		// subquery over the query relation (single shard, where nothing is re-read so
		// a snapshot would be pure overhead). Set once in InitGlobal; keeping the
		// choice here is what lets ClaimWork stay branch-free.
		std::string shard_read_source;

		// Set only when shard_read_source is a snapshot. The connection is held for
		// the life of the state so the destructor can drop the TEMP table: it was
		// created on an inheriting connection, so it lives in the caller's catalog,
		// and a missed drop leaks a relation into the user's session (visible in
		// SHOW TABLES) instead of dying with us.
		std::unique_ptr<Connection> snapshot_conn;
		std::string query_snapshot; // unquoted; empty => no snapshot to drop

		idx_t MaxThreads() const override {
			return max_active_shards * max_threads_per_shard;
		}

		GlobalState() = default;

		~GlobalState() override {
			if (snapshot_conn) {
				DropHelperTempRelation(*snapshot_conn, KeywordHelper::WriteOptionallyQuoted(query_snapshot));
			}
		}
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		std::shared_ptr<ActiveShard> current_active_shard;
		bool has_shard = false;
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
