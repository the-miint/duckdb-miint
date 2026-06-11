#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// Effective bowtie2 `-p` (nthreads) for one shard's batch submit. As worker
// threads finish the tail of small shards and exit, the survivors should grow
// their thread count to fill the cores those workers freed — otherwise a lone
// big shard runs at the base `-p` with the rest of the box idle (measured: a
// dominant shard pinned at -p4 leaves half an 8-core box unused for most of a
// run). bowtie2's `-p` only partitions reads across threads, so the alignments
// produced are independent of its value — this is a pure throughput knob.
//
// Gated on `all_shards_claimed`: while shards are still being handed out we keep
// every worker at the base so nobody oversubscribes the workers about to claim
// the rest. Once the last shard is claimed `active_workers` only decreases, so
// the result is monotonic non-decreasing per surviving worker (bounded daemon
// index reloads, no fingerprint flapping). Never drops below
// `base_threads_per_shard` (the user's floor) and never exceeds `db_threads`
// (the cores actually available).
inline idx_t EffectiveShardThreads(idx_t base_threads_per_shard, idx_t db_threads, idx_t active_workers,
                                   bool all_shards_claimed) {
	if (!all_shards_claimed || active_workers == 0) {
		return base_threads_per_shard;
	}
	const idx_t fair_share = db_threads / active_workers; // <= db_threads (active_workers >= 1)
	return fair_share > base_threads_per_shard ? fair_share : base_threads_per_shard;
}

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
