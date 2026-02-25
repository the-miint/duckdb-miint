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
#include <atomic>
#include <mutex>
#include <vector>

namespace duckdb {

// Information about a single bowtie2 shard
struct Bowtie2ShardInfo {
	std::string name;         // e.g., "shard_a"
	std::string index_prefix; // e.g., "/path/shards/shard_a/index"
	idx_t read_count;         // Number of reads for this shard (for priority ordering)
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

		// Output schema (shared with align_bowtie2)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data() : names(GetAlignmentOutputNames()), types(GetAlignmentOutputTypes()) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		idx_t next_shard_idx = 0;
		idx_t shard_count = 1;
		idx_t total_associations = 0;
		std::atomic<idx_t> associations_processed {0};

		idx_t MaxThreads() const override {
			// One DuckDB thread per shard, each running one single-threaded bowtie2 process.
			// The 'threads' parameter is ignored in sharded mode - parallelism comes from
			// running multiple shards concurrently, not from bowtie2's internal threading.
			return shard_count;
		}

		GlobalState() = default;
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Bowtie2Aligner> aligner;
		idx_t current_shard_idx = DConstants::INVALID_INDEX;
		bool has_shard = false;
		bool finished_aligning = false; // True after calling finish() for current shard
		idx_t current_read_offset = 0;  // For streaming queries
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
};

} // namespace duckdb
