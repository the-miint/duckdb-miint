#pragma once
#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "align_minimap2_common.hpp"
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

// Information about a single shard
struct ShardInfo {
	std::string name;        // e.g., "shard_001"
	std::string index_path;  // e.g., "/path/shards/shard_001.mmi"
	idx_t read_count;        // Number of reads for this shard (for priority ordering)
};

class AlignMinimap2ShardedTableFunction {
public:
	struct Data : public TableFunctionData {
		std::string query_table;
		std::string shard_directory;
		std::string read_to_shard_table;
		SequenceTableSchema query_schema;
		miint::Minimap2Config config;
		std::vector<ShardInfo> shards;  // Sorted by read_count DESC (largest first)

		// Output schema (shared with align_minimap2)
		std::vector<std::string> names;
		std::vector<LogicalType> types;

		Data()
		    : names(GetAlignmentOutputNames()),
		      types(GetAlignmentOutputTypes()) {
		}
	};

	struct GlobalState : public GlobalTableFunctionState {
		std::mutex lock;
		idx_t next_shard_idx = 0;
		idx_t shard_count = 1;

		idx_t MaxThreads() const override {
			// No cap - one thread per shard, let DuckDB scheduler manage
			return shard_count;
		}

		GlobalState() = default;
	};

	struct LocalState : public LocalTableFunctionState {
		std::unique_ptr<miint::Minimap2Aligner> aligner;
		idx_t current_shard_idx = DConstants::INVALID_INDEX;
		bool has_shard = false;
		idx_t current_read_offset = 0;  // For LIMIT/OFFSET streaming per shard
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

	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
