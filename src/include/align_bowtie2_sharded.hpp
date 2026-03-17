#pragma once
#include "Bowtie2Aligner.hpp"
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
#include "duckdb/parallel/task_scheduler.hpp"
#include <algorithm>
#include <atomic>
#include <chrono>
#include <condition_variable>
#include <mutex>
#include <thread>
#include <vector>

namespace duckdb {

// Information about a single bowtie2 shard
struct Bowtie2ShardInfo {
	std::string name;         // e.g., "shard_a"
	std::string index_prefix; // e.g., "/path/shards/shard_a/index"
	idx_t read_count;         // Number of reads for this shard (for priority ordering)
};

// A shard that is currently being processed by one or more threads.
// Each thread spawns its own bowtie2 subprocess (no shared index object).
// Atomic counters coordinate batch claiming and worker tracking without holding the global lock.
struct Bowtie2ActiveShard {
	idx_t shard_idx;                            // Index into Data::shards
	miint::SequenceRecordBatch shard_sequences; // Pre-fetched sequences for this shard
	std::atomic<idx_t> next_batch_offset {0};   // Threads atomically claim ranges into shard_sequences
	std::atomic<idx_t> active_workers {0};      // Threads currently on this shard
	std::atomic<bool> exhausted {false};        // No more batches to claim
	std::atomic<bool> ready {false};            // IDs materialized, ready for workers to join
};

class AlignBowtie2ShardedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string shard_directory;
		std::string read_to_shard_table;
		SequenceTableSchema query_schema;
		miint::Bowtie2Config config;
		std::vector<Bowtie2ShardInfo> shards; // Sorted by read_count DESC (largest first)
		idx_t max_threads_per_shard = 4;
		bool debug = false;
		bool include_shard_name = false;

		// Output schema (shared with align_bowtie2)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data() : names(GetAlignmentOutputNames()), types(GetAlignmentOutputTypes()) {
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
		std::chrono::steady_clock::time_point start_time;
		std::vector<std::shared_ptr<Bowtie2ActiveShard>> active_shards;
		std::atomic<idx_t> total_associations {0};
		std::atomic<idx_t> associations_processed {0};

		idx_t MaxThreads() const override {
			return max_active_shards * max_threads_per_shard;
		}

		GlobalState() = default;
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Bowtie2Aligner> aligner;
		std::shared_ptr<Bowtie2ActiveShard> current_active_shard;
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
	// Returns the Bowtie2ActiveShard to work on, or nullptr if no more work.
	// ID materialization happens outside the lock.
	static std::shared_ptr<Bowtie2ActiveShard> ClaimWork(ClientContext &context, GlobalState &gstate,
	                                                     const Data &bind_data, LocalState &lstate);

	// Release work: detach from current shard, clean up if last worker on exhausted shard.
	static void ReleaseWork(GlobalState &gstate, LocalState &lstate);
};

} // namespace duckdb
