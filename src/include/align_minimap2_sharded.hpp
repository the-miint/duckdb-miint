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
#include <algorithm>
#include <atomic>
#include <mutex>
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
	std::shared_ptr<miint::SharedMinimap2Index> index; // Shared index, immutable after construction
	std::atomic<idx_t> next_batch_offset {0};          // Threads atomically claim ranges
	std::atomic<idx_t> active_workers {0};             // Threads currently on this shard
	std::atomic<bool> exhausted {false};               // Set when no more batches to read
	std::atomic<bool> ready {false};                   // Set when index is loaded and published
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

		// Output schema (shared with align_minimap2)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data() : names(GetAlignmentOutputNames()), types(GetAlignmentOutputTypes()) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		idx_t next_shard_idx = 0;
		idx_t shard_count = 0;
		idx_t max_threads_per_shard = 4;
		std::vector<std::shared_ptr<ActiveShard>> active_shards;
		idx_t total_associations = 0;
		std::atomic<idx_t> associations_processed {0};

		idx_t MaxThreads() const override {
			// Cap to avoid oversubscription: at most a few shards are active
			// simultaneously, so requesting shard_count * max_threads_per_shard
			// would be misleading. We cap at a reasonable multiple.
			constexpr idx_t MAX_CONCURRENT_SHARDS = 4;
			idx_t concurrent = std::min(shard_count, MAX_CONCURRENT_SHARDS);
			return concurrent * max_threads_per_shard;
		}

		GlobalState() = default;
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		std::shared_ptr<ActiveShard> current_active_shard;
		bool has_shard = false;
		miint::SAMRecordBatch result_buffer;
		idx_t buffer_offset = 0;

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
	static std::shared_ptr<ActiveShard> ClaimWork(GlobalState &gstate, const Data &bind_data, LocalState &lstate);

	// Release work: detach from current shard, clean up if last worker on exhausted shard.
	static void ReleaseWork(GlobalState &gstate, LocalState &lstate);
};

} // namespace duckdb
