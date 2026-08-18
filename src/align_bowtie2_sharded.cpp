#include "align_bowtie2_sharded.hpp"
#include "align_bowtie2_daemon_common.hpp"
#include "align_common.hpp"
#include "catalog_utils.hpp"
#include "miint_log.hpp"
#include "sequence_table_reader.hpp"
#include "shard_progress.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parallel/task_scheduler.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <atomic>

#include "gpl_boundary/arrow_ipc.hpp"
#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"

#include <chrono>
#include <cstdlib>
#include <cstring>
#include <fcntl.h>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

// Telemetry timing uses a monotonic clock (wall-clock-independent, never goes
// backwards across NTP steps). All timing is gated behind a single fd>=0 check
// in Execute, so when telemetry is off these are never called.
using TelClock = std::chrono::steady_clock;
inline double TelMsSince(TelClock::time_point t0) {
	return std::chrono::duration<double, std::milli>(TelClock::now() - t0).count();
}

// Append one telemetry line to the (already-open) fd from a worker thread. A
// SINGLE write() per line is deliberate, and is what keeps this lock-free:
//   - regular file (a path in MIINT_BT2_TELEMETRY): O_APPEND makes the seek-to-
//     EOF + write atomic w.r.t. other writers on the same inode, so concurrent
//     workers' whole lines never interleave (Linux holds the inode lock across
//     the write — not subject to PIPE_BUF for regular files).
//   - stderr ("stderr"/"1"): if stderr is a pipe, atomicity holds only up to
//     PIPE_BUF (4096 B). Our lines are a few hundred bytes (the daemon `metrics`
//     object is a small rusage blob), so they stay well under that.
// A partial short-write would tear a line, but looping would split into multiple
// write()s and re-introduce interleaving — for a best-effort diagnostic the
// single write() is the right trade. The result is ignored: a failed telemetry
// write must never abort a real query.
inline void EmitBatchTelemetry(int fd, const BatchTelemetry &rec) {
	const std::string line = FormatBatchTelemetry(rec);
	const ssize_t n = ::write(fd, line.data(), line.size());
	(void)n;
}

// =============================================================================
// Per-caller known-parameter set. Common bowtie2-align knobs live in
// `bt2_daemon::kCommonAlignParams`; we union those with the sharded-specific
// miint-side knobs (`shard_directory`, `read_to_shard`, `threads`,
// `max_threads_per_shard`, `include_shard_name`). Anything outside this set
// is rejected at bind time so typos surface at SQL-compile rather than
// running silently with default semantics.
// =============================================================================

std::unordered_set<std::string> MakeKnownShardedParams() {
	auto s = bt2_daemon::kCommonAlignParams;
	s.insert("shard_directory");
	s.insert("read_to_shard");
	s.insert("threads"); // sharded-mode warning; not forwarded to daemon
	s.insert("max_threads_per_shard");
	s.insert("include_shard_name");
	s.insert("submit_batch_reads");
	s.insert("prefetch_ahead");
	s.insert("progress");
	return s;
}

// =============================================================================
// Shard index discovery. Bowtie2 emits 4 mandatory files per index in either
// small (.bt2) or large (.bt2l) format.
// =============================================================================

bool HasShardIndex(const std::string &prefix) {
	namespace fs = std::filesystem;
	// A valid index is the complete small (.bt2) set OR the complete large
	// (.bt2l) set — the same suffix lists ShardIndexFiles enumerates for prefetch
	// (single source of truth in align_bowtie2_sharded.hpp).
	auto all_exist = [&](const auto &suffixes) {
		for (const char *ext : suffixes) {
			if (!fs::exists(prefix + ext)) {
				return false;
			}
		}
		return true;
	};
	return all_exist(kBowtie2IndexSuffixesSmall) || all_exist(kBowtie2IndexSuffixesLarge);
}

// Best-effort: warm a shard's index files into the OS page cache ahead of the
// daemon's mmap'd (`--mm`) load, so the worker that claims this shard next skips
// the cold network-FS fault that otherwise dominates the long tail (5–99s for a
// few-read shard). POSIX_FADV_WILLNEED is non-blocking and consumes no alignment
// thread — the kernel does the readahead — so it adds load concurrency without
// touching the per-shard `-p` core budget (the only in-allocation way to overlap
// more than `max_active_shards` cold loads at once). The WILLNEED hint itself is
// non-blocking, but the per-file `open()` it requires IS a blocking cold-FS
// metadata round-trip on the claiming worker's critical path; left in place
// deliberately, pending an HPC measurement showing it actually stalls the claim
// frontier before adding a background-prefetch mechanism to hide it. A pure cache
// hint: errors are ignored, alignment results are unaffected. Compiles to a no-op
// where POSIX_FADV_WILLNEED is unavailable (macOS/WASM, which never run the daemon).
//
// Takes the precomputed file list (cached in InitGlobal) rather than re-deriving
// it: a fresh ShardIndexFiles() here would issue up to 8 stat()s per call on the
// very filesystem whose slowness we are working around.
void PrefetchShardIndexFiles(const std::vector<std::string> &index_files) {
#ifdef POSIX_FADV_WILLNEED
	for (const auto &path : index_files) {
		const int fd = ::open(path.c_str(), O_RDONLY | O_CLOEXEC);
		if (fd < 0) {
			continue;
		}
		// offset 0, length 0 = advise the whole file.
		(void)::posix_fadvise(fd, 0, 0, POSIX_FADV_WILLNEED);
		::close(fd);
	}
#else
	(void)index_files;
#endif
}

// read_to_shard schema validation is handled by align_common.hpp's catalog-
// backed `ValidateReadToShardSchema` (case-insensitive column names, shared
// with align_minimap2_sharded). Direct call below in Bind().

struct ShardInfo {
	std::string name;
	std::string index_prefix; // <shard_directory>/<shard_name>/index
	idx_t read_count;
};

// Enrich align_common.hpp's name+count rows with the bowtie2 index path. The
// underlying GROUP BY / NULL-guard / largest-first ordering is shared with
// align_minimap2_sharded via `ReadShardNameCounts`; only the
// `<dir>/<name>/index` path layout is bowtie2-specific (minimap2 uses a
// single per-shard file, not a 4-file prefix). Per-shard filesystem checks
// (`HasShardIndex`) are intentionally deferred to InitGlobal — hitting the
// filesystem from Bind scales poorly over NFS / hundreds of shards, and the
// check is racy regardless (a shard can be deleted between bind and execute).
std::vector<ShardInfo> EnumerateShards(ClientContext &context, const std::string &read_to_shard_table,
                                       const std::string &shard_directory) {
	const auto counts = ReadShardNameCounts(context, read_to_shard_table);
	std::vector<ShardInfo> out;
	out.reserve(counts.size());
	for (const auto &c : counts) {
		ShardInfo info;
		info.name = c.name;
		info.read_count = c.count;
		info.index_prefix = shard_directory;
		if (!info.index_prefix.empty() && info.index_prefix.back() != '/') {
			info.index_prefix += '/';
		}
		info.index_prefix += info.name + "/index";
		out.push_back(std::move(info));
	}
	return out;
}

// =============================================================================
// BindData / GlobalState
// =============================================================================

struct AlignBowtie2ShardedBindData : public TableFunctionData {
	std::string query_table;
	std::string shard_directory;
	std::string read_to_shard_table;
	named_parameter_map_t named_params;

	// Detected at bind time; affects the per-batch Arrow IPC encoding.
	bool query_has_sequence2 = false;
	bool query_has_qual1 = false;
	bool query_has_qual2 = false;

	// Id-column types captured by ValidateSequenceTableSchema. The query side
	// may be VARCHAR or BIGINT; the read_to_shard table's read_id must match
	// (enforced by ValidateReadToShardSchema with expected_read_id_type). The
	// subject side is always VARCHAR because sharded mode loads prebuilt
	// bowtie2 indexes whose reference names are opaque bytes — same contract
	// as align_minimap2_sharded.
	LogicalType query_id_type = LogicalType(LogicalTypeId::INVALID);
	LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);

	// Ordered shards (largest first) with index path resolved + verified.
	std::vector<ShardInfo> shards;

	// Default 1: one shard per lane — maximize shard-level concurrency with a
	// single-threaded bowtie2 each. For many-shard workloads on a cold network FS
	// (the dominant sharded case) per-shard index load dominates and is hidden best
	// by prefetch + mm-off across many parallel lanes; this measured faster than
	// fewer fat lanes (e.g. one shard per core beat ceil(cores/4) shards at `-p4`).
	// Raise for few-shard / few-core workloads where bowtie2's internal `-p`
	// threading matters more than the number of concurrent shards.
	idx_t max_threads_per_shard = 1;
	bool include_shard_name = false;
	// Opt-in per-shard progress to stderr (default false; see shard_progress.hpp).
	bool progress = false;

	// Lower-bound threshold on reads accumulated into one daemon Submit. Larger
	// batches amortize bowtie2's per-batch FM-index reload and keep `bowtie2 -p N`
	// fed; throughput knob only — alignment results are identical to any other
	// value because each read is aligned independently. Accumulation is
	// chunk-granular (one DuckDB chunk <= STANDARD_VECTOR_SIZE at a time in
	// FetchShardBatch), so a value below the chunk size still submits a full
	// chunk per round-trip.
	idx_t submit_batch_reads = 16384;

	// How many upcoming shards' index files to warm into page cache (via
	// POSIX_FADV_WILLNEED) ahead of the claim frontier — see
	// PrefetchShardIndexFiles. -1 is an INTERNAL "auto" sentinel (the field
	// default); users get auto by omitting the parameter, and the auto value is
	// resolved to max_active_shards in InitGlobal (GlobalState::prefetch_ahead).
	// User-supplied values are bind-validated to 0 (disabled) .. 4096. Throughput
	// knob only — alignment output is identical for any value (a cache hint, not
	// a semantic change).
	int64_t prefetch_ahead = -1;
};

struct AlignBowtie2ShardedGlobalState : public GlobalTableFunctionState {
	// Shared shard-claim counter. Each Execute thread atomically increments
	// to grab the next unclaimed shard; when this exceeds shards.size() the
	// scan is done. No per-shard mutex needed — the counter is the only
	// shared state across worker threads.
	std::atomic<idx_t> next_shard_idx {0};

	// Per-shard worker concurrency. Derived from DuckDB's thread pool and
	// `max_threads_per_shard`: ceil(db_threads / max_threads_per_shard),
	// clamped to the number of shards. Each miint worker thread orchestrates
	// one shard at a time; the daemon's per-fingerprint worker pool fans out
	// internally on `nthreads` = max_threads_per_shard, so total CPU usage
	// is roughly max_active_shards × max_threads_per_shard ≈ db_threads.
	idx_t max_active_shards = 1;

	// Live count of worker threads currently holding a shard (++ on claim, -- on
	// release in Execute). Read at submit time to redistribute idle cores to
	// survivors via `EffectiveShardThreads` once `next_shard_idx` shows all
	// shards are claimed. Advisory only (never affects alignment results), so
	// relaxed ordering is sufficient.
	std::atomic<idx_t> active_workers {0};

	// DuckDB thread-pool size, captured in InitGlobal. The ceiling for a
	// survivor's bowtie2 `-p` when it grows to fill freed cores.
	idx_t db_threads = 1;

	// Resolved prefetch lookahead (BindData.prefetch_ahead, with its -1 "auto"
	// sentinel collapsed to max_active_shards here in InitGlobal where db_threads
	// is known). 0 = disabled. Read directly in Execute — no per-claim ternary.
	idx_t prefetch_ahead = 0;

	// Per-shard index file paths, parallel to bd.shards, computed once in
	// InitGlobal (the filesystem-metadata cache is already warm there from the
	// HasShardIndex pass). Claim-time prefetch indexes into this so it issues NO
	// stats on the (slow, the whole reason we prefetch) shard filesystem.
	std::vector<std::vector<std::string>> shard_index_files;

	// Per-batch telemetry sink (env-gated by MIINT_BT2_TELEMETRY; resolved once
	// in InitGlobal). -1 = disabled (the only cost when off is the fd>=0 check at
	// each batch). When enabled it's either an O_APPEND file fd we own (and must
	// close) or STDERR_FILENO (which we must NOT close). `telemetry_start` anchors
	// the per-line `wall_ms`. See InitGlobal for the gate and EmitBatchTelemetry
	// for the write. Like the other fields above (db_threads, max_active_shards,
	// prefetch_ahead), these are written once in InitGlobal — before any Execute
	// runs, under the scheduler's thread-launch barrier — then read-only from the
	// worker threads, so no atomics are needed.
	int telemetry_fd = -1;
	bool telemetry_owns_fd = false;
	TelClock::time_point telemetry_start;

	// Opt-in per-shard progress (mirrors BindData::progress; set in InitGlobal,
	// read-only from workers). When false, no progress lines are emitted.
	bool progress = false;

	// How many reads each shard's cursor should deliver, keyed by shard name,
	// computed once in InitGlobal (#229). Empty for a single shard, where nothing is
	// re-read and the check cannot detect anything a lone cursor would have missed.
	//
	// Deliberately NOT ShardInfo::read_count, which is COUNT(*) over read_to_shard
	// alone: a mapping that legitimately lists reads absent from the query relation
	// would then look like data loss. These counts come from the same
	// `query JOIN read_to_shard` the cursors use, so a mapping superset reduces
	// expectation and actual equally, and only a genuine short read shows up.
	//
	// One extra scan of the query relation, projecting read_id only, and no
	// materialization — so the cursors keep their bounded memory profile.
	std::unordered_map<std::string, idx_t> expected_reads_per_shard;

	idx_t MaxThreads() const override {
		// DuckDB clamps to its own scheduler concurrency anyway, but
		// returning the per-shard ceiling here keeps the planner honest
		// when the thread pool is larger than the work supports.
		return max_active_shards;
	}

	~AlignBowtie2ShardedGlobalState() override {
		if (telemetry_owns_fd && telemetry_fd >= 0) {
			::close(telemetry_fd);
		}
	}
};

// Per-thread state — each miint worker thread owns its own gpl-boundary
// daemon process (spawned eagerly in InitLocal), its own DuckDB connection
// for the streaming input cursor, and its own decode buffers. This keeps
// daemon I/O lock-free: Session::Submit is request-response synchronous on
// a single pipe pair, so two threads sharing one Session would race on
// stdin writes and stdout reads. One Session per thread sidesteps that
// entirely.
//
// Member-declaration order is load-bearing for clean teardown. C++ destroys
// members in reverse declaration order, so:
//   1. `current_batches` (Arrow arrays referencing SubmitResult SHM regions)
//      destroys first.
//   2. `current_decoder` (which can still write to its source bytes) next.
//   3. `current_result` (owner of the SHM regions) last among the decode set.
//   4. `input_stream` (holds an open cursor on `input_conn`) destroys before
//      `input_conn`. Reordering these two would leave a dangling cursor on
//      a dead connection during teardown.
//   5. `session` destroys last so any I/O races during shutdown are bounded
//      to a per-thread daemon that's already in its own process.
// DO NOT reorder these fields without re-reading the above; the silence on
// failure makes the bug class very expensive to find.
struct AlignBowtie2ShardedLocalState : public LocalTableFunctionState {
	std::unique_ptr<gb::Session> session;
	std::unique_ptr<Connection> input_conn;

	// Current shard claim. Sentinel value DConstants::INVALID_INDEX means
	// "no shard claimed yet"; Execute will claim next on the next iteration.
	idx_t current_shard_idx = DConstants::INVALID_INDEX;
	std::unique_ptr<QueryResult> input_stream;
	std::string current_shard_name; // copied for `include_shard_name`
	// Set once `input_stream->Fetch()` returns EOF. A streaming QueryResult
	// throws "closed pending query result" if Fetch() is called again after it
	// has returned its final (empty) chunk, so once a FetchShardBatch loop
	// consumes EOF we must not Fetch this stream again — the next call returns
	// false immediately. Reset when a new shard stream is opened. Plain bool,
	// order-independent — safe to declare outside the SHM/Arrow ordering block.
	bool stream_exhausted = false;

	// Latest Submit response + its decoder. The ArrowArrayWrapper batches
	// inside `current_batches` hold pointers into SubmitResult's SHM regions,
	// so we keep both alive via unique_ptr ordering — see member-order note
	// above the struct.
	std::unique_ptr<gb::SubmitResult> current_result;
	std::unique_ptr<gb::IpcStreamDecoder> current_decoder;
	ArrowSchemaWrapper current_schema;
	std::vector<ArrowArrayWrapper> current_batches;
	idx_t batch_index = 0;
	idx_t row_in_batch = 0;
	bool schema_validated = false;

	// Per-claim config-JSON cache. BuildAlignConfigJson is a pure function of
	// (shard.index_prefix, nthreads): the param map is constant for the whole
	// query and index_prefix is fixed by current_shard_idx, so the config string
	// only changes when this worker moves to a new shard or the EffectiveShardThreads
	// ramp steps nthreads up. We key on (shard_idx, nthreads) and rebuild only on a
	// miss — avoiding a full param-map re-serialization on every batch (the common
	// case is one shard, many batches, constant nthreads). Output-identical to
	// rebuilding each batch. Keying on shard_idx (not just nthreads) is load-bearing
	// for correctness: two consecutive shards a worker processes can share the same
	// nthreads, and reusing the prior shard's config would send reads to the wrong
	// index. Plain values, order-independent — outside the SHM/Arrow teardown block.
	std::string cached_config_json;
	idx_t cached_config_shard_idx = DConstants::INVALID_INDEX;
	idx_t cached_config_nthreads = 0;

	// Telemetry-only (read/written solely on the MIINT_BT2_TELEMETRY path).
	// `telemetry_batch_seq` orders this worker's batches in the TSV;
	// `telemetry_open_stream_ms` carries the cursor-open cost from a shard claim
	// forward onto that shard's first batch line (0 on subsequent batches).
	// Progress-only per-shard accumulators (used when GlobalState::progress is
	// true). Reset when a shard stream opens; read when it exhausts.
	// `shard_reads` is accumulated unconditionally, not just when progress is on:
	// it is also compared against the shard's expected count when the cursor
	// exhausts (#229). One add per batch, so it is free.
	idx_t shard_reads = 0;
	idx_t shard_alignments = 0;
	TelClock::time_point shard_start;

	idx_t telemetry_batch_seq = 0;
	double telemetry_open_stream_ms = 0.0;

	~AlignBowtie2ShardedLocalState() override {
		if (session) {
			try {
				session->Shutdown();
			} catch (...) {
				// Destructor must not throw; child reaped by ~ChildProcess.
			}
		}
	}
};

// Detect optional sequence2/qual1/qual2 columns on the query table. Validates
// types here so a wrong-type column fails at bind time rather than later via a
// cryptic worker crash (see ../docs/internals/embedded-tools.md: the daemon
// receives qual1/qual2 as Phred+33 Utf8 strings; the canonical miint
// representation produced by `read_fastx` is `UTINYINT[]` raw Phred values,
// and `FetchShardBatch` converts on the fly. Any other column shape is an
// error — we never want users to be in the business of pre-encoding ASCII
// quality strings themselves).
void DetectQueryColumns(ClientContext &context, AlignBowtie2ShardedBindData &bd) {
	auto conn = MakeReadOnlyHelperConnection(context);
	const std::string sql = "SELECT column_name, column_type FROM (DESCRIBE " +
	                        KeywordHelper::WriteOptionallyQuoted(bd.query_table) +
	                        ") WHERE column_name IN ('sequence2','qual1','qual2')";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("align_bowtie2_sharded: failed to introspect query table '%s': %s", bd.query_table,
		                            result->GetError());
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); ++i) {
			const auto col = chunk->GetValue(0, i).ToString();
			const auto typ = chunk->GetValue(1, i).ToString();
			if (col == "sequence2") {
				bd.query_has_sequence2 = true;
			} else if (col == "qual1") {
				bd.query_has_qual1 = true;
				if (typ != "UTINYINT[]") {
					throw BinderException("align_bowtie2_sharded: column 'qual1' in '%s' must be UTINYINT[] (raw Phred "
					                      "values, as produced by read_fastx); got %s",
					                      bd.query_table, typ);
				}
			} else if (col == "qual2") {
				bd.query_has_qual2 = true;
				if (typ != "UTINYINT[]") {
					throw BinderException("align_bowtie2_sharded: column 'qual2' in '%s' must be UTINYINT[] (raw Phred "
					                      "values, as produced by read_fastx); got %s",
					                      bd.query_table, typ);
				}
			}
		}
	}
}

// Build the bowtie2-align config_json for one shard. Each shard has a
// distinct `index_path` → distinct daemon worker fingerprint, so concurrent
// Submits from different miint worker threads fan out to independent
// per-fingerprint workers in the daemon — that's how cross-shard
// parallelism falls out for free.
std::string BuildAlignConfigJson(const named_parameter_map_t &named_params, const std::string &index_prefix,
                                 idx_t max_threads_per_shard) {
	bt2_daemon::ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_prefix);
	cfg.append_int("nthreads", static_cast<int64_t>(max_threads_per_shard));

	bt2_daemon::AppendBowtie2AlignParams(cfg, named_params, "align_bowtie2_sharded");

	// Sharded contract: emit only mapped reads. The pre-migration direct-subprocess
	// path called FilterMappedOnly on every batch; on the daemon path we delegate
	// to bowtie2's own --no-unal so unaligned records never cross the boundary.
	// `no_unal` is intentionally excluded from kCommonAlignParams so users can't
	// override this here.
	cfg.append_bool("no_unal", true);

	return cfg.build();
}

// =============================================================================
// Bind
// =============================================================================

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<std::string> &names) {
	if (input.inputs.size() < 1) {
		throw BinderException("align_bowtie2_sharded requires query_table parameter");
	}
	auto bd = make_uniq<AlignBowtie2ShardedBindData>();
	bd->query_table = input.inputs[0].ToString();

	auto shard_dir_param = input.named_parameters.find("shard_directory");
	if (shard_dir_param == input.named_parameters.end() || shard_dir_param->second.IsNull()) {
		throw BinderException("align_bowtie2_sharded requires shard_directory parameter");
	}
	bd->shard_directory = shard_dir_param->second.ToString();

	auto read_to_shard_param = input.named_parameters.find("read_to_shard");
	if (read_to_shard_param == input.named_parameters.end() || read_to_shard_param->second.IsNull()) {
		throw BinderException("align_bowtie2_sharded requires read_to_shard parameter");
	}
	bd->read_to_shard_table = read_to_shard_param->second.ToString();

	auto &fs = FileSystem::GetFileSystem(context);
	if (!fs.DirectoryExists(bd->shard_directory)) {
		throw BinderException("Shard directory does not exist: %s", bd->shard_directory);
	}

	// Validate query_table (existence + read_id type). Captures VARCHAR or
	// BIGINT into bd->query_id_type for the output schema and the egress
	// codec dispatch in EmitChunkRows. Replaces the earlier DESCRIBE-only
	// existence probe.
	auto query_schema = ValidateSequenceTableSchema(context, bd->query_table, /*allow_bigint=*/true);
	bd->query_id_type = query_schema.id_type;

	// Subject side is always VARCHAR — sharded mode loads prebuilt bowtie2
	// indexes whose reference names are opaque bytes (same contract as
	// align_minimap2_sharded). Captured here so the emit path and output
	// schema use the same value the bind committed to.
	bd->subject_id_type = LogicalType::VARCHAR;

	// read_to_shard.read_id must match the query table's read_id type. The
	// strict equality check prevents the downstream JOIN inside
	// OpenCurrentShardStream from relying on implicit casts.
	ValidateReadToShardSchema(context, bd->read_to_shard_table, bd->query_id_type);

	// Reject unknown params at bind time.
	for (const auto &kv : input.named_parameters) {
		static const auto kKnown = MakeKnownShardedParams();
		if (kKnown.find(kv.first) == kKnown.end()) {
			throw InvalidInputException("align_bowtie2_sharded: unknown named parameter '%s'", kv.first);
		}
	}
	bd->named_params = input.named_parameters;

	// Sharded default: mm-off (sequential fread) unless the user set memory_mapped.
	// Injected after the unknown-param check and the input copy, so it reaches the
	// daemon via AppendBowtie2AlignParams (which consumes bd->named_params at submit
	// time) without tripping the typo guard above. See InjectMemoryMappedDefault.
	InjectMemoryMappedDefault(bd->named_params, [](bool b) { return Value::BOOLEAN(b); });

	// `threads` is ignored at the table-function level: with the Phase 6
	// fan-out, cross-shard parallelism is driven by DuckDB's own scheduler
	// (`SET threads=N`), and per-shard bowtie2 internal threading is driven
	// by `max_threads_per_shard` (the daemon's `nthreads`).
	// Preserved the "Parameter 'threads' is ignored in sharded mode" prefix
	// because miint_warnings_bowtie2.test asserts on that substring.
	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		const int64_t threads_val = bt2_daemon::ValueAsInt("align_bowtie2_sharded", "threads", threads_param->second);
		if (threads_val != 1) {
			::miint::EmitWarning(
			    context, "WARNING: Parameter 'threads' is ignored in sharded mode. "
			             "Use `SET threads=N` to control cross-shard parallelism (one worker thread per shard); "
			             "use `max_threads_per_shard` to control bowtie2's internal threads per shard.");
		}
	}

	auto max_tps_param = input.named_parameters.find("max_threads_per_shard");
	if (max_tps_param != input.named_parameters.end() && !max_tps_param->second.IsNull()) {
		const int64_t val =
		    bt2_daemon::ValueAsInt("align_bowtie2_sharded", "max_threads_per_shard", max_tps_param->second);
		if (val < 1 || val > 64) {
			throw BinderException("max_threads_per_shard must be between 1 and 64 (got %lld)",
			                      static_cast<long long>(val));
		}
		bd->max_threads_per_shard = static_cast<idx_t>(val);
	}

	auto include_shard_param = input.named_parameters.find("include_shard_name");
	if (include_shard_param != input.named_parameters.end() && !include_shard_param->second.IsNull()) {
		bd->include_shard_name =
		    bt2_daemon::ValueAsBool("align_bowtie2_sharded", "include_shard_name", include_shard_param->second);
	}

	// Opt-in per-shard progress to stderr (default false): a programmatic caller
	// that never passes progress:=true emits nothing.
	auto progress_param = input.named_parameters.find("progress");
	if (progress_param != input.named_parameters.end() && !progress_param->second.IsNull()) {
		bd->progress = bt2_daemon::ValueAsBool("align_bowtie2_sharded", "progress", progress_param->second);
	}

	auto batch_param = input.named_parameters.find("submit_batch_reads");
	if (batch_param != input.named_parameters.end() && !batch_param->second.IsNull()) {
		const int64_t val = bt2_daemon::ValueAsInt("align_bowtie2_sharded", "submit_batch_reads", batch_param->second);
		if (val < 1 || val > 1000000) {
			throw BinderException("submit_batch_reads must be between 1 and 1000000 (got %lld)",
			                      static_cast<long long>(val));
		}
		bd->submit_batch_reads = static_cast<idx_t>(val);
	}

	auto prefetch_param = input.named_parameters.find("prefetch_ahead");
	if (prefetch_param != input.named_parameters.end() && !prefetch_param->second.IsNull()) {
		const int64_t val = bt2_daemon::ValueAsInt("align_bowtie2_sharded", "prefetch_ahead", prefetch_param->second);
		if (val < 0 || val > 4096) {
			throw BinderException("prefetch_ahead must be between 0 (disabled) and 4096 (got %lld)",
			                      static_cast<long long>(val));
		}
		bd->prefetch_ahead = val;
	}

	DetectQueryColumns(context, *bd);

	bd->shards = EnumerateShards(context, bd->read_to_shard_table, bd->shard_directory);

	// Output schema reflects captured id types: read_id mirrors the query
	// side; reference and mate_reference are always VARCHAR (subject side
	// is opaque bytes in the prebuilt bowtie2 index).
	bt2_daemon::PopulateOutputSchema(names, return_types, bd->query_id_type, bd->subject_id_type);
	if (bd->include_shard_name) {
		names.emplace_back("shard_name");
		return_types.emplace_back(LogicalType::VARCHAR);
	}

	return std::move(bd);
}

// =============================================================================
// InitGlobal
// =============================================================================

std::unique_ptr<gb::Session> SpawnAndCheckSession() {
	const std::string gpl_path = gb::FindGplBoundary();
	if (gpl_path.empty()) {
		throw IOException("align_bowtie2_sharded: gpl-boundary binary not found. To install:\n"
		                  "  Easiest:  SELECT install_gpl_boundary();\n"
		                  "  Manual:   curl -fsSL "
		                  "https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh | sh\n"
		                  "If gpl-boundary is installed at a non-standard location, set "
		                  "MIINT_GPL_BOUNDARY_PATH=<absolute path>.");
	}
	// Fail loud if the resolved daemon predates bowtie2 `memory_mapped` support
	// (>= 0.4.2): the sharded mm-off default depends on it, and older daemons
	// silently ignore the field — reintroducing the cold-FS regression.
	bt2_daemon::RequireGplBoundaryVersion(gpl_path, "align_bowtie2_sharded");
	std::vector<std::string> argv = {gpl_path};
	gb::ChildProcess child(argv);
	auto session = std::make_unique<gb::Session>(std::move(child));
	session->Initialize();
	if (!session->has_tool("bowtie2-align")) {
		throw IOException("align_bowtie2_sharded: gpl-boundary daemon does not advertise bowtie2-align. "
		                  "Upgrade to v0.2.0 or later (`SELECT install_gpl_boundary()`).");
	}
	const uint32_t got = session->tool_schema_version("bowtie2-align");
	if (got != bt2_daemon::kAlignSchemaVersion) {
		throw IOException("align_bowtie2_sharded: daemon reports bowtie2-align schema_version=%u, miint expects %u.",
		                  got, bt2_daemon::kAlignSchemaVersion);
	}
	return session;
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bd = input.bind_data->Cast<AlignBowtie2ShardedBindData>();

	auto gs = make_uniq<AlignBowtie2ShardedGlobalState>();

	// Resolve concurrency + prefetch depth FIRST, before the validation loop.
	// Each miint worker thread orchestrates one shard at a time; the daemon's
	// per-fingerprint pool fans out internally on `nthreads`
	// (= max_threads_per_shard). So our slice of the DuckDB thread pool is
	// `ceil(db_threads / max_threads_per_shard)`, clamped by the shard count
	// (no point asking for more workers than shards). Mirrors the formula
	// `align_minimap2_sharded` uses to honor `SET threads` without a separate
	// knob. The prefetch_ahead -1 "auto" sentinel resolves to max_active_shards
	// (one full wave of lead); we resolve it here so the loop below can skip the
	// per-shard index-file enumeration entirely when prefetch is disabled.
	const idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	const idx_t derived = (db_threads + bd.max_threads_per_shard - 1) / bd.max_threads_per_shard;
	gs->max_active_shards = std::max<idx_t>(1, std::min<idx_t>(derived, bd.shards.size()));
	gs->db_threads = db_threads;
	gs->progress = bd.progress;
	gs->prefetch_ahead = bd.prefetch_ahead < 0 ? gs->max_active_shards : static_cast<idx_t>(bd.prefetch_ahead);
	const bool prefetch_enabled = gs->prefetch_ahead > 0;

	// Verify each shard's bowtie2 index files exist on disk. Done here rather
	// than in Bind because filesystem stats from the planner cost real time
	// on NFS / very wide shard sets, and the check is best-effort anyway
	// (the daemon will fail at Submit if a shard disappears mid-query). When
	// prefetch is enabled, cache the existing index-file paths now — HasShardIndex
	// has just warmed the metadata cache for them — so claim-time prefetch never
	// re-stats the slow shard filesystem. When prefetch is disabled, skip the
	// enumeration entirely: shard_index_files stays empty and the (ahead>0)-guarded
	// prefetch in Execute never reads it.
	if (prefetch_enabled) {
		gs->shard_index_files.reserve(bd.shards.size());
	}
	for (const auto &shard : bd.shards) {
		if (!HasShardIndex(shard.index_prefix)) {
			throw IOException("No valid bowtie2 index found at prefix: %s. "
			                  "Expected files like %s.1.bt2, %s.rev.1.bt2, etc.",
			                  shard.index_prefix, shard.index_prefix, shard.index_prefix);
		}
		if (prefetch_enabled) {
			gs->shard_index_files.push_back(ShardIndexFiles(shard.index_prefix));
		}
	}

	// Telemetry gate (resolved once, here). MIINT_BT2_TELEMETRY is a file path,
	// or the literal "stderr"/"1" for fd 2. Purely diagnostic and output-
	// invariant, so a bad path warns and disables rather than failing the query.
	if (const char *tel_env = std::getenv("MIINT_BT2_TELEMETRY")) {
		if (tel_env[0] != '\0') {
			const std::string spec(tel_env);
			if (spec == "stderr" || spec == "1") {
				gs->telemetry_fd = STDERR_FILENO;
				gs->telemetry_owns_fd = false;
			} else {
				const int fd = ::open(spec.c_str(), O_WRONLY | O_CREAT | O_APPEND, 0644);
				if (fd < 0) {
					::miint::EmitWarning(context,
					                     "align_bowtie2_sharded: MIINT_BT2_TELEMETRY set but could not open '" + spec +
					                         "' for append; telemetry disabled.");
				} else {
					gs->telemetry_fd = fd;
					gs->telemetry_owns_fd = true;
					// Header only for a freshly-created (empty) file, so a file
					// appended across runs stays a single parseable TSV. (stderr
					// gets no header — it interleaves with daemon output.)
					struct stat st;
					if (::fstat(fd, &st) == 0 && st.st_size == 0) {
						const std::string hdr = BatchTelemetryHeader();
						const ssize_t n = ::write(fd, hdr.data(), hdr.size());
						(void)n;
					}
				}
			}
			if (gs->telemetry_fd >= 0) {
				gs->telemetry_start = TelClock::now();
			}
		}
	}

	// Per-shard expected read counts, from the same JOIN the cursors use (#229), so
	// the comparison in Execute turns what used to be silent truncation into a hard
	// error. Skipped for a single shard: its lone cursor reads the relation once, so
	// there is nothing to detect and no reason to pay an extra scan — matching the
	// single-shard skip in align_minimap2_sharded.
	if (bd.shards.size() > 1) {
		for (auto &sc : ReadShardNameCounts(context, bd.read_to_shard_table, bd.query_table)) {
			gs->expected_reads_per_shard[sc.name] = sc.count;
		}
	}

	return std::move(gs);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &,
                                              GlobalTableFunctionState *) {
	auto ls = make_uniq<AlignBowtie2ShardedLocalState>();
	// One DuckDB connection per worker thread — `QueryResult::Fetch` isn't
	// thread-safe across a shared connection, and each shard needs its own
	// streaming cursor.
	ls->input_conn = std::make_unique<Connection>(DatabaseInstance::GetDatabase(context.client));
	InheritTempObjects(context.client, *ls->input_conn);
	// Spawn the per-thread daemon eagerly: if gpl-boundary is missing /
	// wrong-version, we want that error to surface before any row flows.
	// Lazy-spawning inside Execute() would let some threads stream rows
	// while another thread fails its first Spawn, which would leave the
	// user with partial Parquet output plus a confusing mid-scan abort.
	// The ~50ms fork+Initialize cost is paid up front per worker; this is
	// dominated by per-shard index-load latency in any real workload.
	ls->session = SpawnAndCheckSession();
	return std::move(ls);
}

// =============================================================================
// Execute — per-shard sequential processing
// =============================================================================

// Open a streaming cursor on the reads matched to the current shard.
void OpenCurrentShardStream(AlignBowtie2ShardedLocalState &local, const AlignBowtie2ShardedBindData &bd) {
	const auto &shard = bd.shards[local.current_shard_idx];
	std::string select = "SELECT q.read_id, q.sequence1";
	if (bd.query_has_sequence2) {
		select += ", q.sequence2";
	}
	if (bd.query_has_qual1) {
		select += ", q.qual1";
	}
	if (bd.query_has_qual2) {
		select += ", q.qual2";
	}
	select += " FROM " + KeywordHelper::WriteOptionallyQuoted(bd.query_table) + " q";
	select += " JOIN " + KeywordHelper::WriteOptionallyQuoted(bd.read_to_shard_table) + " rts";
	select += " ON q.read_id = rts.read_id";
	// shard.name is a user-supplied row value from `read_to_shard`. Direct
	// concatenation into the WHERE clause would be a SQL-injection vector
	// (e.g. shard names containing a single quote could rewrite the predicate
	// and silently send reads to the wrong index). WriteQuoted wraps the
	// value in single quotes and doubles any embedded single quote, matching
	// the convention already used in sequence_table_reader.cpp:377.
	select += " WHERE rts.shard_name = " + KeywordHelper::WriteQuoted(shard.name, '\'');

	local.input_stream = local.input_conn->SendQuery(select);
	local.stream_exhausted = false;
	if (local.input_stream->HasError()) {
		throw InvalidInputException("align_bowtie2_sharded: failed to open cursor for shard '%s': %s", shard.name,
		                            local.input_stream->GetError());
	}
	local.current_shard_name = shard.name;
}

// Quality decoding lives in align_bowtie2_daemon_common — see
// bt2_daemon::DecodeListQualToPhred33. Shared with align_bowtie2.

// Accumulate reads from the current shard's stream into one daemon Submit:
// pull whole DuckDB chunks (each <= STANDARD_VECTOR_SIZE) until reaching
// `submit_batch_reads` or the stream exhausts. Larger batches amortize
// bowtie2's per-batch FM-index reload and keep `bowtie2 -p N` fed; results are
// identical to one-chunk-per-Submit because each read is aligned independently.
bool FetchShardBatch(AlignBowtie2ShardedLocalState &local, const AlignBowtie2ShardedBindData &bd,
                     bt2_daemon::QueryBatch &out) {
	if (local.stream_exhausted) {
		return false; // never Fetch() past EOF — see stream_exhausted note
	}
	out.read_ids.clear();
	out.sequence1.clear();
	out.sequence2.clear();
	out.sequence2_valid.clear();
	out.qual1.clear();
	out.qual1_valid.clear();
	out.qual2.clear();
	out.qual2_valid.clear();
	int col = 0;
	const int col_read_id = col++;
	const int col_sequence1 = col++;
	const int col_sequence2 = bd.query_has_sequence2 ? col++ : -1;
	const int col_qual1 = bd.query_has_qual1 ? col++ : -1;
	const int col_qual2 = bd.query_has_qual2 ? col++ : -1;
	std::vector<uint8_t> qual_scratch; // reused across rows to amortize allocs
	std::string qual_encoded;
	idx_t accumulated = 0;
	while (accumulated < bd.submit_batch_reads) {
		auto chunk = local.input_stream->Fetch();
		if (!chunk || chunk->size() == 0) {
			local.stream_exhausted = true;
			break;
		}
		const idx_t n = chunk->size();
		out.read_ids.reserve(out.read_ids.size() + n);
		out.sequence1.reserve(out.sequence1.size() + n);
		// Reserve the optional columns too: with batching these vectors span
		// many chunks, so growing them by push_back alone would repeatedly
		// reallocate (and copy the std::strings) as the batch fills.
		if (col_sequence2 >= 0) {
			out.sequence2.reserve(out.sequence2.size() + n);
			out.sequence2_valid.reserve(out.sequence2_valid.size() + n);
		}
		if (col_qual1 >= 0) {
			out.qual1.reserve(out.qual1.size() + n);
			out.qual1_valid.reserve(out.qual1_valid.size() + n);
		}
		if (col_qual2 >= 0) {
			out.qual2.reserve(out.qual2.size() + n);
			out.qual2_valid.reserve(out.qual2_valid.size() + n);
		}
		for (idx_t i = 0; i < n; ++i) {
			auto rid = chunk->GetValue(col_read_id, i);
			auto s1 = chunk->GetValue(col_sequence1, i);
			if (rid.IsNull() || s1.IsNull()) {
				throw InvalidInputException("align_bowtie2_sharded: NULL read_id or sequence1 in query table '%s'",
				                            bd.query_table);
			}
			out.read_ids.push_back(rid.GetValue<std::string>());
			out.sequence1.push_back(s1.GetValue<std::string>());
			if (col_sequence2 >= 0) {
				auto v = chunk->GetValue(col_sequence2, i);
				if (v.IsNull()) {
					out.sequence2.emplace_back();
					out.sequence2_valid.push_back(0);
				} else {
					out.sequence2.push_back(v.GetValue<std::string>());
					out.sequence2_valid.push_back(1);
				}
			}
			if (col_qual1 >= 0) {
				auto v = chunk->GetValue(col_qual1, i);
				if (v.IsNull()) {
					out.qual1.emplace_back();
					out.qual1_valid.push_back(0);
				} else {
					bt2_daemon::DecodeListQualToPhred33(v, "qual1", bd.query_table, qual_encoded, qual_scratch);
					out.qual1.push_back(qual_encoded);
					out.qual1_valid.push_back(1);
				}
			}
			if (col_qual2 >= 0) {
				auto v = chunk->GetValue(col_qual2, i);
				if (v.IsNull()) {
					out.qual2.emplace_back();
					out.qual2_valid.push_back(0);
				} else {
					bt2_daemon::DecodeListQualToPhred33(v, "qual2", bd.query_table, qual_encoded, qual_scratch);
					out.qual2.push_back(qual_encoded);
					out.qual2_valid.push_back(1);
				}
			}
		}
		accumulated += n;
	}
	return accumulated > 0;
}

// `tel` is non-null only when telemetry is enabled; when set, this stamps the
// encode / submit-round-trip / decode phase timings, the input/output byte
// sizes, the decoded alignment count, and the raw daemon `metrics` JSON into it.
void SubmitAndDecode(AlignBowtie2ShardedLocalState &local, const AlignBowtie2ShardedBindData &bd,
                     const AlignBowtie2ShardedGlobalState &gs, const bt2_daemon::QueryBatch &qb, BatchTelemetry *tel) {
	const auto &shard = bd.shards[local.current_shard_idx];
	// Once every shard has been handed out, surviving workers grow their `-p` to
	// reclaim the cores that finished workers freed (see EffectiveShardThreads).
	// Changing nthreads re-fingerprints this shard in the daemon, spawning a new
	// bowtie2 worker that reloads the FM-index. Under the sharded mm-off default
	// (memory_mapped=false) that reload is a full sequential fread of the (multi-GB)
	// index, NOT the cheap warm-mmap minor-fault reload — far from free on a cold
	// network FS. It stays bounded, though: nthreads only ever steps UP, and only
	// when db_threads/active_workers crosses an integer, so a surviving worker pays
	// it a handful of times at most, amortized over that shard's remaining batches
	// at the higher `-p`. Kept per-batch deliberately: claim-locking nthreads would
	// pin a big head shard (claimed first under largest-first ordering, while shards
	// are still being handed out) at the base `-p`, surrendering the very tail-core
	// reclaim this ramp exists for — claim-lock was evaluated for this workload and
	// rejected.
	const bool all_shards_claimed = gs.next_shard_idx.load(std::memory_order_relaxed) >= bd.shards.size();
	const idx_t nthreads = EffectiveShardThreads(bd.max_threads_per_shard, gs.db_threads,
	                                             gs.active_workers.load(std::memory_order_relaxed), all_shards_claimed);
	// Reuse the cached config_json unless this worker moved to a new shard or the
	// ramp stepped nthreads (see the cache fields on LocalState). config_json is a
	// pure function of (shard.index_prefix, nthreads); rebuilding only on a miss
	// avoids re-serializing the full param map every batch and is output-identical.
	if (local.current_shard_idx != local.cached_config_shard_idx || nthreads != local.cached_config_nthreads) {
		local.cached_config_json = BuildAlignConfigJson(bd.named_params, shard.index_prefix, nthreads);
		local.cached_config_shard_idx = local.current_shard_idx;
		local.cached_config_nthreads = nthreads;
	}
	const std::string &config_json = local.cached_config_json;
	bt2_daemon::QueryArrowSchema schema_flags {bd.query_has_sequence2, bd.query_has_qual1, bd.query_has_qual2};
	const auto t_encode0 = tel ? TelClock::now() : TelClock::time_point {};
	const auto ipc = bt2_daemon::BuildQueryIpc(qb, schema_flags);
	if (tel) {
		tel->t_encode_ms = TelMsSince(t_encode0);
		tel->input_bytes = static_cast<idx_t>(ipc.size());
	}
	const auto t_submit0 = tel ? TelClock::now() : TelClock::time_point {};
	// Opt into the daemon's per-batch worker metrics (getrusage: ru_majflt, CPU
	// vs wall, RSS, reused-flag) only when telemetry is on — they land in the
	// telemetry line's `metrics` column. Off ⇒ no flag, no daemon-side cost.
	auto submit_result =
	    local.session->Submit("bowtie2-align", config_json, ipc.data(), ipc.size(), /*request_metrics=*/tel != nullptr);
	if (tel) {
		tel->t_submit_ms = TelMsSince(t_submit0);
	}
	if (submit_result.outputs.empty()) {
		throw IOException("align_bowtie2_sharded: daemon returned zero shm_outputs for shard '%s'", shard.name);
	}
	if (submit_result.schema_version != bt2_daemon::kAlignSchemaVersion) {
		throw IOException("align_bowtie2_sharded: daemon returned schema_version=%u, expected %u",
		                  submit_result.schema_version, bt2_daemon::kAlignSchemaVersion);
	}

	local.current_batches.clear();
	if (local.current_schema.arrow_schema.release) {
		local.current_schema.arrow_schema.release(&local.current_schema.arrow_schema);
	}
	local.current_decoder.reset();
	local.current_result.reset();
	local.batch_index = 0;
	local.row_in_batch = 0;

	local.current_result = std::make_unique<gb::SubmitResult>(std::move(submit_result));
	const auto &out0 = local.current_result->outputs[0];
	if (tel) {
		tel->output_bytes = static_cast<idx_t>(out0.size_bytes());
		tel->metrics = local.current_result->metrics_json; // empty on daemons without metrics
	}
	local.current_decoder = std::make_unique<gb::IpcStreamDecoder>(out0.bytes(), out0.size_bytes());
	local.current_decoder->GetSchema(&local.current_schema.arrow_schema);

	if (!local.schema_validated) {
		bt2_daemon::ValidateOutputSchema(local.current_schema.arrow_schema);
		local.schema_validated = true;
	}

	const auto t_decode0 = tel ? TelClock::now() : TelClock::time_point {};
	for (;;) {
		ArrowArrayWrapper w;
		if (!local.current_decoder->NextBatch(&w.arrow_array)) {
			break;
		}
		local.current_batches.push_back(std::move(w));
	}
	if (tel) {
		tel->t_decode_ms = TelMsSince(t_decode0);
		idx_t rows = 0;
		for (const auto &w : local.current_batches) {
			rows += static_cast<idx_t>(w.arrow_array.length);
		}
		tel->n_alignments = rows;
	}
}

// Synthesize `shard_name` into the last output column when include_shard_name=true.
void FillShardNameColumn(DataChunk &output, idx_t to_emit, const std::string &shard_name) {
	auto &v = output.data[bt2_daemon::kNumOutputColumns]; // 21st (0-indexed) column
	auto *out_data = FlatVector::GetData<string_t>(v);
	auto &validity = FlatVector::Validity(v);
	for (idx_t i = 0; i < to_emit; ++i) {
		out_data[i] = StringVector::AddString(v, shard_name);
		validity.SetValid(i);
	}
}

void Execute(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	(void)context;
	auto &bd = data.bind_data->Cast<AlignBowtie2ShardedBindData>();
	auto &gs = data.global_state->Cast<AlignBowtie2ShardedGlobalState>();
	auto &local = data.local_state->Cast<AlignBowtie2ShardedLocalState>();

	while (true) {
		// 1. Drain any decoded rows from the current Submit response.
		while (local.batch_index < local.current_batches.size()) {
			auto &batch = local.current_batches[local.batch_index].arrow_array;
			if (batch.offset != 0) {
				throw IOException(
				    "align_bowtie2_sharded: decoded batch has non-zero parent offset (%lld); not yet handled",
				    static_cast<long long>(batch.offset));
			}
			const idx_t total = static_cast<idx_t>(batch.length);
			const idx_t remaining = total - local.row_in_batch;
			if (remaining == 0) {
				local.batch_index += 1;
				local.row_in_batch = 0;
				continue;
			}
			const idx_t to_emit = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
			output.SetCardinality(to_emit);
			bt2_daemon::EmitChunkRows(output, to_emit, local.row_in_batch, batch, bd.query_id_type, bd.subject_id_type);
			if (bd.include_shard_name) {
				FillShardNameColumn(output, to_emit, local.current_shard_name);
			}
			local.row_in_batch += to_emit;
			return;
		}

		// 2. Current Submit is drained. Pull next batch from the shard this
		//    thread is currently working on (if any). `FetchShardBatch`
		//    returns false iff the stream is fully exhausted; a successful
		//    fetch always populates at least one row (the NULL guard inside
		//    the loop throws rather than returning an empty batch), so there
		//    is no need for an "empty but non-exhausted" code path.
		if (local.current_shard_idx != DConstants::INVALID_INDEX) {
			bt2_daemon::QueryBatch qb;
			BatchTelemetry rec;
			BatchTelemetry *tel = gs.telemetry_fd >= 0 ? &rec : nullptr;
			const auto t_fetch0 = tel ? TelClock::now() : TelClock::time_point {};
			const bool got = FetchShardBatch(local, bd, qb);
			if (tel) {
				rec.t_fetch_ms = TelMsSince(t_fetch0);
			}
			if (got) {
				if (tel) {
					rec.worker_id = static_cast<int64_t>(local.session->daemon_pid());
					rec.shard = local.current_shard_name;
					rec.batch_seq = local.telemetry_batch_seq++;
					rec.n_reads = static_cast<idx_t>(qb.read_ids.size());
					// Attribute the shard's one-time cursor-open cost to its first
					// batch, then clear so later batches of the same shard read 0.
					rec.t_open_stream_ms = local.telemetry_open_stream_ms;
					local.telemetry_open_stream_ms = 0.0;
				}
				SubmitAndDecode(local, bd, gs, qb, tel);
				if (tel) {
					rec.wall_ms = TelMsSince(gs.telemetry_start);
					EmitBatchTelemetry(gs.telemetry_fd, rec);
				}
				local.shard_reads += static_cast<idx_t>(qb.read_ids.size());
				if (gs.progress) {
					for (const auto &w : local.current_batches) {
						local.shard_alignments += static_cast<idx_t>(w.arrow_array.length);
					}
				}
				continue; // loop back to drain decoded rows
			}
			// Shard exhausted. Before releasing it, verify the cursor actually
			// delivered every read mapped to this shard (#229).
			//
			// Each shard opens its own cursor over the query relation, so a relation
			// that is not stable across re-evaluation — a volatile view, a view over
			// a changing table, a registered single-pass Arrow stream — returns a
			// different row set to a later cursor, its JOIN matches less (often
			// nothing), and those reads were previously dropped with no error.
			// Failing here is the whole point: a short read is a wrong answer, and a
			// wrong answer that looks right is worse than a stopped query.
			//
			// Two honest limits. (1) This compares VOLUME, not identity: a relation
			// that re-reads to a *different* set of the same size passes. It catches
			// the realistic failures (an exhausted stream yields nothing; a shifted
			// key yields fewer join matches) and costs one add per batch, where a
			// per-read digest would put hashing on the alignment path. (2) It fires
			// mid-scan, so shards already completed have emitted rows — the partial
			// output that InitLocal's eager daemon spawn (see its comment) exists to
			// avoid. Accepted deliberately: the alternative is knowing at bind
			// whether a relation is replayable, which is not decidable here.
			{
				auto expected = gs.expected_reads_per_shard.find(local.current_shard_name);
				if (expected != gs.expected_reads_per_shard.end() && local.shard_reads < expected->second) {
					throw InvalidInputException(
					    "align_bowtie2_sharded: shard '%s' delivered %llu of %llu mapped reads. The query relation "
					    "returned different rows when re-read for this shard, so reads would be silently dropped. "
					    "This happens when '%s' is not stable across repeated reads (a volatile view, a view over a "
					    "changing table, or a registered single-pass Arrow stream). Materialize it first, e.g. "
					    "CREATE TEMP TABLE q AS SELECT * FROM %s, and pass that instead.",
					    local.current_shard_name, static_cast<uint64_t>(local.shard_reads),
					    static_cast<uint64_t>(expected->second), bd.query_table, bd.query_table);
				}
			}
			// Release for re-claim attempt. Drop this worker from the active count
			// so a surviving worker can grow its `-p`.
			if (gs.progress) {
				const double elapsed_s = std::chrono::duration<double>(TelClock::now() - local.shard_start).count();
				shard_progress::Emit("bowtie2", shard_progress::FormatShardDone(
				                                    static_cast<uint64_t>(local.current_shard_idx) + 1,
				                                    static_cast<uint64_t>(bd.shards.size()), local.current_shard_name,
				                                    local.shard_reads, local.shard_alignments, elapsed_s));
			}
			local.input_stream.reset();
			local.current_shard_idx = DConstants::INVALID_INDEX;
			gs.active_workers.fetch_sub(1, std::memory_order_relaxed);
			// A shard that matched zero reads emits no batch line, so its
			// cursor-open cost has nothing to attach to; clear it so it can't be
			// mistaken for the next shard's open. (The next OpenCurrentShardStream
			// overwrites it regardless, so this is belt-and-braces, but it keeps
			// the per-batch open-cost attribution correct under any reordering.)
			local.telemetry_open_stream_ms = 0.0;
		}

		// 3. Claim the next shard atomically. fetch_add is the entire
		//    coordination — no mutex needed. When the counter exceeds the
		//    shard count this thread is done.
		//    `memory_order_relaxed` is sufficient because the only shared
		//    state we read after claiming is `bd.shards[shard_idx]`, and
		//    `bd.shards` is built once in Bind() and never mutated
		//    afterwards — there is no store that any thread needs to
		//    observe via this atomic. Don't "upgrade" this to acq_rel
		//    without re-checking that invariant.
		const idx_t shard_idx = gs.next_shard_idx.fetch_add(1, std::memory_order_relaxed);
		if (shard_idx >= bd.shards.size()) {
			output.SetCardinality(0);
			return;
		}
		local.current_shard_idx = shard_idx;
		// Count this worker as active for the duration of the shard (balanced by
		// the fetch_sub on exhaustion above). Done only after a valid claim, so
		// the failed-claim path that ends the worker never touches the count.
		gs.active_workers.fetch_add(1, std::memory_order_relaxed);
		// Warm the index of a shard `ahead` positions later into page cache while
		// this worker loads+aligns its own, so by the time some worker claims that
		// shard its mmap'd index isn't a cold fault. `shard_idx` is unique per
		// atomic claim and `ahead` is constant, so each target is warmed at most
		// once (no redundant opens, no shared state). The first `ahead` shards
		// have no earlier claim to warm them, so they load cold — unavoidable at
		// the start. `gs.prefetch_ahead` is the resolved lookahead (0 disables);
		// pure CPU-free kernel readahead, output-invariant.
		const idx_t ahead = gs.prefetch_ahead;
		if (ahead > 0) {
			// shard_idx < shards.size() (claim semantics) + ahead <= 4096
			// (bind-validated) => no meaningful wrap; the guard below is what
			// bounds the access regardless.
			const idx_t target = shard_idx + ahead;
			if (target < bd.shards.size()) {
				PrefetchShardIndexFiles(gs.shard_index_files[target]);
			}
		}
		// Time the cursor open so its cost shows up (on this shard's first
		// batch line) distinctly from the daemon round-trip. Stamped into
		// LocalState; consumed by the next FetchShardBatch's telemetry above.
		const bool tel_on = gs.telemetry_fd >= 0;
		const auto t_open0 = tel_on ? TelClock::now() : TelClock::time_point {};
		OpenCurrentShardStream(local, bd);
		if (tel_on) {
			local.telemetry_open_stream_ms = TelMsSince(t_open0);
		}
		// Per-shard progress accounting: reset accumulators and stamp the start
		// of this shard, then announce it. Read count is unknown at open (reads
		// stream in batches), so the "started" line omits it; the "done" line
		// below reports the totals.
		// Must reset per shard: shard_reads is compared against a per-shard expected
		// count when the cursor exhausts, so carrying it across shards would make
		// that check pass vacuously.
		local.shard_reads = 0;
		local.shard_alignments = 0;
		local.shard_start = TelClock::now();
		if (gs.progress) {
			shard_progress::Emit("bowtie2", shard_progress::FormatShardStart(
			                                    static_cast<uint64_t>(local.current_shard_idx) + 1,
			                                    static_cast<uint64_t>(bd.shards.size()), local.current_shard_name,
			                                    /*n_reads=*/-1));
		}
	}
}

} // namespace

// =============================================================================
// Public registration
// =============================================================================

TableFunction AlignBowtie2ShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	bt2_daemon::RegisterBowtie2AlignNamedParameterTypes(tf);
	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["threads"] = LogicalType::INTEGER; // ignored in sharded mode; warning at bind
	tf.named_parameters["max_threads_per_shard"] = LogicalType::INTEGER;
	tf.named_parameters["include_shard_name"] = LogicalType::BOOLEAN;
	tf.named_parameters["submit_batch_reads"] = LogicalType::INTEGER;
	tf.named_parameters["prefetch_ahead"] = LogicalType::INTEGER;
	tf.named_parameters["progress"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignBowtie2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
