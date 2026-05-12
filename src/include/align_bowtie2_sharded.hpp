#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// align_bowtie2_sharded routes per-shard through the gpl-boundary daemon's
// bowtie2-align tool. Each shard is a pre-built bowtie2 index on disk; the
// table function visits shards sequentially, submitting that shard's
// matched reads with `config.index_path = <shard prefix>`. Cross-shard
// parallelism (Phase 6 follow-up) comes for free from the daemon's
// per-(tool, config_fingerprint) worker pool.
class AlignBowtie2ShardedTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
