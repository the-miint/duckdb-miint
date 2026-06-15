#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

#include <cstdint>
#include <filesystem>
#include <iomanip>
#include <sstream>
#include <string>
#include <vector>

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

// Inject align_bowtie2_sharded's mm-off default: insert `memory_mapped=false`
// into the named-parameter map iff the user did not supply it. HPC telemetry
// (WOL3, 1000 shards, cold Lustre) measured sequential-fread index loads at 3.7x
// the throughput of bowtie2's `--mm` default (worker_majflt 13 vs 47.5M), so
// sharded mode defaults the knob OFF — but a user-supplied value (true or false)
// is preserved, since this is a perf default, not a hard policy. Output-invariant
// either way (mm vs fread produce identical alignments).
//
// Templated over the map type and a bool→value factory so the logic is unit-
// testable in the standalone Catch2 binary, which links no libduckdb (hence no
// `Value` / `named_parameter_map_t` symbols). Production instantiates it with
// `named_parameter_map_t` + `Value::BOOLEAN`; the test uses a plain
// `std::map<std::string,bool>`. Same rationale as keeping EffectiveShardThreads /
// ShardIndexFiles inline and duckdb-free below.
template <class Map, class MakeBool>
void InjectMemoryMappedDefault(Map &params, MakeBool make_bool) {
	if (params.find("memory_mapped") == params.end()) {
		params.emplace("memory_mapped", make_bool(false));
	}
}

// The bowtie2 index files that actually exist for a shard prefix — the subset
// of the 8 candidates (`<prefix>.1.bt2` … `.rev.2.bt2`, and the large-index
// `.bt2l` variants) present on disk. Used to warm them into the OS page cache
// (POSIX_FADV_WILLNEED) ahead of a worker claiming the shard, so the daemon's
// index load doesn't stall on a cold network-FS fault. Under the sharded mm-off
// default (see InjectMemoryMappedDefault) the worker freads the index
// sequentially from offset 0, which is exactly the access pattern WILLNEED
// readahead serves — a better match than the lazy random faults of `--mm`.
// Returns only files that exist: an incomplete index yields a partial list
// (warming a subset is harmless), an absent one yields empty. Inline (like
// EffectiveShardThreads) so the C++ unit test links it without pulling in the
// .cpp's daemon dependencies.
inline std::vector<std::string> ShardIndexFiles(const std::string &prefix) {
	// The 8 possible bowtie2 index files: 4 mandatory in small (.bt2) format,
	// or 4 in large (.bt2l) format for big references. We list whichever exist.
	static const char *const kSuffixes[] = {".1.bt2",  ".2.bt2",  ".rev.1.bt2",  ".rev.2.bt2",
	                                        ".1.bt2l", ".2.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"};
	std::vector<std::string> out;
	for (const char *ext : kSuffixes) {
		std::string p = prefix + ext;
		if (std::filesystem::exists(p)) {
			out.push_back(std::move(p));
		}
	}
	return out;
}

// =============================================================================
// Per-batch telemetry (env-gated by MIINT_BT2_TELEMETRY; see align_bowtie2_
// sharded.cpp). One record per daemon Submit, emitted as a TSV line, splitting
// the per-shard wall time into miint-side phases (fetch / encode / decode) vs.
// the opaque daemon round-trip (`t_submit_ms`) — so we measure where the
// per-shard seconds go instead of inferring it from straces/seff. Output-
// invariant: pure diagnostics, never affects alignment results.
//
// The record struct + formatter live in the header (like EffectiveShardThreads /
// ShardIndexFiles) so the C++ unit test can exercise the formatting — the one
// part that's deterministically testable — without the .cpp's daemon deps.
// =============================================================================
struct BatchTelemetry {
	double wall_ms = 0.0;          // steady-clock ms since scan start, stamped at emit
	int64_t worker_id = 0;         // daemon pid: tags which worker thread (and daemon)
	std::string shard;             // shard name this batch was aligned against
	idx_t batch_seq = 0;           // 0-based batch counter within this worker
	idx_t n_reads = 0;             // reads submitted in this batch
	idx_t input_bytes = 0;         // Arrow IPC bytes sent to the daemon
	idx_t n_alignments = 0;        // alignment rows decoded from the response
	idx_t output_bytes = 0;        // Arrow IPC bytes received from the daemon
	double t_open_stream_ms = 0.0; // cursor open (first batch of a shard; 0 after)
	double t_fetch_ms = 0.0;       // FetchShardBatch: read⋈read_to_shard pull
	double t_encode_ms = 0.0;      // BuildQueryIpc
	double t_submit_ms = 0.0;      // Session::Submit round-trip (daemon black box)
	double t_decode_ms = 0.0;      // IPC decode loop
	std::string metrics;           // raw daemon `metrics` JSON; empty if absent
};

// Tab-separated column names, in the exact order FormatBatchTelemetry writes the
// values. `metrics` is LAST so the daemon's raw JSON (braces/commas) can never be
// mistaken for additional columns. Keep this and FormatBatchTelemetry in lockstep
// — the unit test asserts the header and a data line have the same column count.
inline const char *BatchTelemetryColumns() {
	return "wall_ms\tworker_id\tshard\tbatch_seq\tn_reads\tinput_bytes\t"
	       "n_alignments\toutput_bytes\tt_open_stream_ms\tt_fetch_ms\t"
	       "t_encode_ms\tt_submit_ms\tt_decode_ms\tmetrics";
}

inline std::string BatchTelemetryHeader() {
	return std::string(BatchTelemetryColumns()) + "\n";
}

// Replace TSV-breaking control characters (tab, CR, LF) in a free-form string
// field with spaces, so an embedded one can't split a column or truncate the
// line. Numeric fields never need this; only the two std::string fields (`shard`,
// a user-supplied read_to_shard value, and `metrics`, daemon-provided) do. The
// daemon JSON is re-serialized by yyjson (which escapes control chars), so this
// is belt-and-braces there, but `shard` names are arbitrary user strings.
inline std::string TsvSanitize(std::string s) {
	for (char &c : s) {
		if (c == '\t' || c == '\n' || c == '\r') {
			c = ' ';
		}
	}
	return s;
}

// Render one record as a newline-terminated TSV line. Milliseconds get 3 decimal
// places (microsecond resolution); counts/sizes/ids print as integers. Pure (no
// timing, no I/O) so it is unit-testable in isolation.
inline std::string FormatBatchTelemetry(const BatchTelemetry &r) {
	std::ostringstream os;
	os << std::fixed << std::setprecision(3);
	os << r.wall_ms << '\t' << r.worker_id << '\t' << TsvSanitize(r.shard) << '\t' << r.batch_seq << '\t' << r.n_reads
	   << '\t' << r.input_bytes << '\t' << r.n_alignments << '\t' << r.output_bytes << '\t' << r.t_open_stream_ms
	   << '\t' << r.t_fetch_ms << '\t' << r.t_encode_ms << '\t' << r.t_submit_ms << '\t' << r.t_decode_ms << '\t'
	   << TsvSanitize(r.metrics) << '\n';
	return os.str();
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
