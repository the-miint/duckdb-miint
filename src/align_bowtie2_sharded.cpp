#include "align_bowtie2_sharded.hpp"
#include "align_bowtie2_daemon_common.hpp"
#include "align_common.hpp"
#include "miint_log.hpp"
#include "sequence_table_reader.hpp"

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
// `oversubscribe`, `big_shard_min_batches`, `include_shard_name`). Anything
// outside this set is rejected at bind time so typos surface at SQL-compile
// rather than running silently with default semantics.
// =============================================================================

std::unordered_set<std::string> MakeKnownShardedParams() {
	auto s = bt2_daemon::kCommonAlignParams;
	s.insert("shard_directory");
	s.insert("read_to_shard");
	s.insert("threads"); // sharded-mode warning; not forwarded to daemon
	s.insert("oversubscribe");
	s.insert("big_shard_min_batches");
	s.insert("include_shard_name");
	s.insert("submit_batch_reads");
	s.insert("prefetch_ahead");
	return s;
}

// =============================================================================
// Shard index discovery. Bowtie2 emits 4 mandatory files per index in either
// small (.bt2) or large (.bt2l) format.
// =============================================================================

bool HasShardIndex(const std::string &prefix) {
	namespace fs = std::filesystem;
	const std::vector<std::string> small = {".1.bt2", ".2.bt2", ".rev.1.bt2", ".rev.2.bt2"};
	bool small_ok = true;
	for (const auto &ext : small) {
		if (!fs::exists(prefix + ext)) {
			small_ok = false;
			break;
		}
	}
	if (small_ok) {
		return true;
	}
	const std::vector<std::string> large = {".1.bt2l", ".2.bt2l", ".rev.1.bt2l", ".rev.2.bt2l"};
	for (const auto &ext : large) {
		if (!fs::exists(prefix + ext)) {
			return false;
		}
	}
	return true;
}

// Best-effort: warm a shard's index files into the OS page cache ahead of the
// daemon's sequential fread load (mm-off default; see InjectMemoryMappedDefault),
// so the worker that claims this shard next skips the cold network-FS fault that
// otherwise dominates the long tail (5–99s for a few-read shard). WILLNEED
// readahead reads from offset 0 forward, which exactly matches the mm-off
// fread access pattern (the mismatch that made it useless under `--mm`'s random
// faults is gone). POSIX_FADV_WILLNEED is non-blocking and consumes no alignment
// thread — the kernel does the readahead — so it adds load concurrency without
// touching a lane's `-p` core budget, overlapping more cold tail loads than the
// number of concurrent tail lanes alone would. A pure cache hint: errors are
// ignored, alignment results are unaffected. Compiles to a no-op where
// POSIX_FADV_WILLNEED is unavailable (macOS/WASM, which never run the daemon).
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

	// EXTRA bowtie2 OS threads to run beyond the core budget (`SET threads`), as
	// single-thread (`-p1`) TAIL lanes. The big lane takes ~75% of the budget; the
	// budget remainder plus `oversubscribe` become tail lanes, so total bowtie2
	// threads = db_threads + oversubscribe (see ResolveLaneSplit). Default 2 — a
	// light oversubscription that overlaps a few cold tail loads without starving
	// the big lane (the first HPC run's 1×-p8 + 7×-p1 = 15 threads on 8 cores ran
	// the big lane at -p4 speed). Bind-validated to 0..256 (0 = no oversubscription:
	// total threads = db_threads). Throughput knob only — alignment output invariant.
	int64_t oversubscribe = 2;

	// Breakpoint between the big and tail pools: a shard is "big" iff it needs
	// >= big_shard_min_batches daemon Submits (ceil(read_count / submit_batch_reads);
	// see IsBigShard). Default 2 ⇒ every single-batch shard is tail. Bind-validated
	// to >= 1 (1 ⇒ every non-empty shard is big). Throughput knob only.
	idx_t big_shard_min_batches = 2;

	bool include_shard_name = false;

	// Lower-bound threshold on reads accumulated into one daemon Submit. Larger
	// batches amortize bowtie2's per-batch FM-index reload and keep `bowtie2 -p N`
	// fed; throughput knob only — alignment results are identical to any other
	// value because each read is aligned independently. Accumulation is
	// chunk-granular (one DuckDB chunk <= STANDARD_VECTOR_SIZE at a time in
	// FetchShardBatch), so a value below the chunk size still submits a full
	// chunk per round-trip.
	idx_t submit_batch_reads = 16384;

	// How many upcoming shards' index files to warm into page cache (via
	// POSIX_FADV_WILLNEED) ahead of the per-pool claim frontier — see
	// PrefetchShardIndexFiles. -1 is an INTERNAL "auto" sentinel (the field
	// default); users get auto by omitting the parameter. In InitGlobal this
	// resolves per pool (GlobalState::prefetch_ahead_big / _tail): big = 1, tail =
	// tail_lanes under auto, or the explicit value for both pools' tail frontier.
	// 0 disables prefetch entirely. User-supplied values are bind-validated to
	// 0 .. 4096. Throughput knob only — alignment output is identical for any
	// value (a cache hint, not a semantic change).
	int64_t prefetch_ahead = -1;
};

struct AlignBowtie2ShardedGlobalState : public GlobalTableFunctionState {
	// Two-pool scheduler. Shards are partitioned in InitGlobal into BIG
	// (compute-bound; run one-at-a-time on the full-thread `-p db_threads` lane)
	// and TAIL (load/FS-bound; run single-threaded `-p1`, many concurrently so
	// their cold index loads overlap). Both vectors hold ORIGINAL indices into
	// bd.shards, each preserving bd.shards' largest-first order. Built once in
	// InitGlobal before any Execute runs; read-only thereafter.
	std::vector<idx_t> big_shards;
	std::vector<idx_t> tail_shards;

	// Lock-free claim frontiers, one per pool (same idiom as the old single
	// next_shard_idx). A worker grabs the next unclaimed position via fetch_add;
	// when the returned position is >= the pool size the pool is drained.
	//   - next_big_idx is touched ONLY by the big worker (ordinal 0), so it needs
	//     no atomicity for correctness, but stays atomic to keep one idiom.
	//   - next_tail_idx is shared by the tail workers AND the big worker (which
	//     falls through to help drain the tail once bigs are done); fetch_add makes
	//     that race-free.
	// relaxed ordering suffices: uniqueness of the value each thread gets is
	// guaranteed by the atomicity of fetch_add itself (an indivisible RMW — two
	// threads can never read the same position), NOT by the memory order. The
	// order only governs visibility of surrounding stores, and the only shared
	// state read after a claim is big_shards/tail_shards/bd.shards, all built
	// before Execute under the scheduler's thread-launch barrier (which provides
	// that release/acquire) — so no store needs to be observed via these atomics.
	std::atomic<idx_t> next_big_idx {0};
	std::atomic<idx_t> next_tail_idx {0};

	// Hands each InitLocal call a unique 0-based ordinal. Ordinal 0 is the single
	// BIG worker; all others are TAIL workers. fetch_add gives the uniqueness.
	std::atomic<idx_t> next_worker_ordinal {0};

	// DuckDB thread-pool size (`SET threads`) = the core budget, captured in InitGlobal.
	idx_t db_threads = 1;

	// Resolved thread split (ResolveLaneSplit, in InitGlobal once db_threads is known).
	idx_t big_lane_threads = 1; // big lane's bowtie2 `-p` (~75% of db_threads)
	idx_t tail_lanes = 0;       // count of concurrent `-p1` tail lanes (budget remainder + oversubscribe)

	// Resolved prefetch lookahead, per pool. 0 disables that pool's prefetch.
	// big = 1 (one shard ahead; the big lane is sequential, single-consumer) and
	// is NOT user-tunable (YAGNI); tail = oversubscribe under auto (warm one full
	// wave ahead) or the explicit prefetch_ahead value. Read directly in Execute.
	idx_t prefetch_ahead_big = 0;
	idx_t prefetch_ahead_tail = 0;

	// Per-shard index file paths, parallel to bd.shards (indexed by ORIGINAL idx),
	// computed once in InitGlobal (the filesystem-metadata cache is already warm
	// there from the HasShardIndex pass). Claim-time prefetch indexes into this so
	// it issues NO stats on the (slow, the whole reason we prefetch) shard FS.
	std::vector<std::vector<std::string>> shard_index_files;

	// Per-batch telemetry sink (env-gated by MIINT_BT2_TELEMETRY; resolved once
	// in InitGlobal). -1 = disabled (the only cost when off is the fd>=0 check at
	// each batch). When enabled it's either an O_APPEND file fd we own (and must
	// close) or STDERR_FILENO (which we must NOT close). `telemetry_start` anchors
	// the per-line `wall_ms`. See InitGlobal for the gate and EmitBatchTelemetry
	// for the write. Like the resolved knobs above (db_threads, big_lane_threads,
	// tail_lanes, prefetch_ahead_*), these are written once in InitGlobal — before
	// any Execute runs, under the scheduler's thread-launch barrier — then read-only
	// from the worker threads, so no atomics are needed.
	int telemetry_fd = -1;
	bool telemetry_owns_fd = false;
	TelClock::time_point telemetry_start;

	idx_t MaxThreads() const override {
		// One big lane + `tail_lanes` tail lanes, clamped to the shard count (no
		// point asking for more workers than shards). DuckDB further clamps to its
		// own `SET threads`, so the effective concurrent lanes are
		// min(1 + tail_lanes, SET threads, n_shards). max(1, ...) guards the empty
		// shard set so we never return 0 (which the planner treats oddly).
		const idx_t n_shards = big_shards.size() + tail_shards.size();
		return std::max<idx_t>(1, std::min<idx_t>(1 + tail_lanes, n_shards));
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

	// This worker's role in the two-pool scheduler, assigned once in InitLocal:
	// ordinal 0 is the single BIG worker (drains big_shards on the `-p db_threads`
	// lane, then helps drain the tail), every other ordinal is a TAIL worker
	// (tail_shards only, `-p1`). Set once, then read-only. See Execute's claim.
	idx_t worker_ordinal = 0;
	// Whether the CURRENTLY-claimed shard is a big-pool shard (⇒ `-p db_threads`)
	// or a tail-pool shard (⇒ `-p1`). Set by each claim in Execute; read by
	// SubmitAndDecode to pick `nthreads`. The big worker flips this to false when
	// it falls through to tail shards.
	bool current_is_big = false;
	// Set once (big worker only) when next_big_idx first runs past big_shards: after
	// that the big worker claims only from the tail pool, so it must stop bumping
	// next_big_idx — otherwise every later claim does a useless atomic RMW on it.
	// Single-threaded (this worker writes and reads it), so a plain bool is enough.
	bool big_pool_exhausted = false;

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

	// Telemetry-only (read/written solely on the MIINT_BT2_TELEMETRY path).
	// `telemetry_batch_seq` orders this worker's batches in the TSV;
	// `telemetry_open_stream_ms` carries the cursor-open cost from a shard claim
	// forward onto that shard's first batch line (0 on subsequent batches).
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
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
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
                                 idx_t nthreads) {
	bt2_daemon::ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_prefix);
	cfg.append_int("nthreads", static_cast<int64_t>(nthreads));

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

	// `threads` is ignored at the table-function level: the two-pool scheduler
	// derives the big lane's bowtie2 `-p` from DuckDB's own scheduler (`SET
	// threads=N`), and the tail lanes always run `-p1`. `oversubscribe` and
	// `big_shard_min_batches` are the sharded-mode tuning knobs.
	// Preserved the "Parameter 'threads' is ignored in sharded mode" prefix
	// because miint_warnings_bowtie2.test asserts on that substring.
	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		const int64_t threads_val = bt2_daemon::ValueAsInt("align_bowtie2_sharded", "threads", threads_param->second);
		if (threads_val != 1) {
			::miint::EmitWarning(context,
			                     "WARNING: Parameter 'threads' is ignored in sharded mode. "
			                     "Use `SET threads=N` to set the core budget (the big lane runs ~75% of it); "
			                     "use `oversubscribe` for extra single-thread tail threads beyond the budget, and "
			                     "`big_shard_min_batches` to set the big/tail breakpoint.");
		}
	}

	// Extra bowtie2 threads beyond the core budget, run as `-p1` tail lanes
	// (see ResolveLaneSplit). 0..256 (0 = no oversubscription). Default 2.
	auto oversub_param = input.named_parameters.find("oversubscribe");
	if (oversub_param != input.named_parameters.end() && !oversub_param->second.IsNull()) {
		const int64_t val = bt2_daemon::ValueAsInt("align_bowtie2_sharded", "oversubscribe", oversub_param->second);
		if (val < 0 || val > 256) {
			throw BinderException("oversubscribe must be between 0 and 256 (got %lld)", static_cast<long long>(val));
		}
		bd->oversubscribe = val;
	}

	// Big/tail breakpoint: a shard is "big" iff it needs >= this many daemon
	// Submits (see IsBigShard). >= 1 (1 ⇒ every non-empty shard is big).
	auto min_batches_param = input.named_parameters.find("big_shard_min_batches");
	if (min_batches_param != input.named_parameters.end() && !min_batches_param->second.IsNull()) {
		const int64_t val =
		    bt2_daemon::ValueAsInt("align_bowtie2_sharded", "big_shard_min_batches", min_batches_param->second);
		if (val < 1) {
			throw BinderException("big_shard_min_batches must be >= 1 (got %lld)", static_cast<long long>(val));
		}
		bd->big_shard_min_batches = static_cast<idx_t>(val);
	}

	auto include_shard_param = input.named_parameters.find("include_shard_name");
	if (include_shard_param != input.named_parameters.end() && !include_shard_param->second.IsNull()) {
		bd->include_shard_name =
		    bt2_daemon::ValueAsBool("align_bowtie2_sharded", "include_shard_name", include_shard_param->second);
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

	// Verify each shard's bowtie2 index files exist on disk. Done here rather
	// than in Bind because filesystem stats from the planner cost real time
	// on NFS / very wide shard sets, and the check is best-effort anyway
	// (the daemon will fail at Submit if a shard disappears mid-query).
	auto gs = make_uniq<AlignBowtie2ShardedGlobalState>();
	gs->shard_index_files.reserve(bd.shards.size());
	for (const auto &shard : bd.shards) {
		if (!HasShardIndex(shard.index_prefix)) {
			throw IOException("No valid bowtie2 index found at prefix: %s. "
			                  "Expected files like %s.1.bt2, %s.rev.1.bt2, etc.",
			                  shard.index_prefix, shard.index_prefix, shard.index_prefix);
		}
		// Cache the file list for prefetch now, while HasShardIndex has just warmed
		// the metadata cache for these paths — so the claim-time prefetch never
		// re-stats the slow shard filesystem.
		gs->shard_index_files.push_back(ShardIndexFiles(shard.index_prefix));
	}

	// Partition shards into the two pools, preserving bd.shards' largest-first
	// order within each (so the big lane processes the heaviest shard first, and
	// tail lanes drain largest-first too). Both vectors hold ORIGINAL indices into
	// bd.shards so shard_index_files[idx] / bd.shards[idx] stay valid lookups.
	for (idx_t i = 0; i < bd.shards.size(); ++i) {
		if (IsBigShard(bd.shards[i].read_count, bd.submit_batch_reads, bd.big_shard_min_batches)) {
			gs->big_shards.push_back(i);
		} else {
			gs->tail_shards.push_back(i);
		}
	}

	// `SET threads` is the core budget. The big lane takes ~75% of it; the budget
	// remainder plus `oversubscribe` extra threads become `-p1` tail lanes, so total
	// bowtie2 threads = db_threads + oversubscribe (see ResolveLaneSplit). max(1, ..):
	// a thread pool always has >=1 worker, but guard a degenerate NumberOfThreads()==0
	// so the big lane never sends nthreads=0 to the daemon.
	const idx_t db_threads =
	    std::max<idx_t>(1, NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads()));
	gs->db_threads = db_threads;
	const auto split = ResolveLaneSplit(db_threads, static_cast<idx_t>(bd.oversubscribe));
	gs->big_lane_threads = split.big_lane_threads;
	gs->tail_lanes = split.tail_lanes;

	// Resolve per-pool prefetch lookahead. prefetch_ahead==0 disables both pools.
	// Otherwise the big lane warms 1 ahead (hardcoded; it is sequential and not a
	// knob), and the tail warms `tail_lanes` ahead under auto (-1) — one full wave
	// of lead — or the explicit value the user gave.
	if (bd.prefetch_ahead == 0) {
		gs->prefetch_ahead_big = 0;
		gs->prefetch_ahead_tail = 0;
	} else {
		gs->prefetch_ahead_big = 1;
		gs->prefetch_ahead_tail = bd.prefetch_ahead < 0 ? gs->tail_lanes : static_cast<idx_t>(bd.prefetch_ahead);
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

	return std::move(gs);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &,
                                              GlobalTableFunctionState *gstate) {
	auto ls = make_uniq<AlignBowtie2ShardedLocalState>();
	// Assign this worker's two-pool role. fetch_add hands out unique 0-based
	// ordinals across the (sequential) InitLocal calls; ordinal 0 is the single
	// BIG worker, all others are TAIL workers. relaxed is fine — InitLocal runs
	// before Execute under the scheduler's launch barrier, and the value is only
	// read by this worker's own Execute thereafter.
	auto &gs = gstate->Cast<AlignBowtie2ShardedGlobalState>();
	ls->worker_ordinal = gs.next_worker_ordinal.fetch_add(1, std::memory_order_relaxed);
	// One DuckDB connection per worker thread — `QueryResult::Fetch` isn't
	// thread-safe across a shared connection, and each shard needs its own
	// streaming cursor.
	ls->input_conn = std::make_unique<Connection>(DatabaseInstance::GetDatabase(context.client));
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
	// Per-role bowtie2 `-p`, fixed for the life of this shard's claim (no mid-run
	// ramp, so no daemon re-fingerprint churn): the big lane gets ~75% of the cores;
	// tail lanes run single-threaded and overlap each other on I/O.
	const idx_t nthreads = local.current_is_big ? gs.big_lane_threads : 1;
	if (tel) {
		tel->is_big = local.current_is_big; // tag the lane + its resolved -p so offline
		tel->nthreads = nthreads;           // analysis attributes time to big vs tail
	}
	const std::string config_json = BuildAlignConfigJson(bd.named_params, shard.index_prefix, nthreads);
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
				continue; // loop back to drain decoded rows
			}
			// Shard exhausted — release so the loop claims this worker's next shard.
			local.input_stream.reset();
			local.current_shard_idx = DConstants::INVALID_INDEX;
			// A shard that matched zero reads emits no batch line, so its
			// cursor-open cost has nothing to attach to; clear it so it can't be
			// mistaken for the next shard's open. (The next OpenCurrentShardStream
			// overwrites it regardless, so this is belt-and-braces, but it keeps
			// the per-batch open-cost attribution correct under any reordering.)
			local.telemetry_open_stream_ms = 0.0;
		}

		// 3. Claim the next shard from this worker's reachable pool(s). fetch_add
		//    is the entire coordination — no mutex needed. relaxed ordering is
		//    sufficient because the only shared state read after a claim is
		//    big_shards/tail_shards/bd.shards, all built before Execute (under the
		//    scheduler's thread-launch barrier) and never mutated — there is no
		//    store any thread needs to observe via these atomics. Don't "upgrade"
		//    to acq_rel without re-checking that invariant.
		//
		// The BIG worker (ordinal 0) drains big_shards first, then FALLS THROUGH to
		// help drain the tail; TAIL workers (ordinal != 0) only ever claim tail
		// shards. A worker is done when its reachable pool(s) are exhausted.
		const std::vector<idx_t> *pool = nullptr; // chosen pool list (for prefetch)
		idx_t pool_pos = 0;                       // claimed position within `pool`
		idx_t ahead = 0;                          // resolved lookahead for `pool`
		bool claimed = false;
		if (local.worker_ordinal == 0 && !local.big_pool_exhausted) {
			const idx_t p = gs.next_big_idx.fetch_add(1, std::memory_order_relaxed);
			if (p < gs.big_shards.size()) {
				local.current_shard_idx = gs.big_shards[p];
				local.current_is_big = true;
				pool = &gs.big_shards;
				pool_pos = p;
				ahead = gs.prefetch_ahead_big;
				claimed = true;
			} else {
				// Big pool drained. Latch so future iterations skip the (now always
				// failing) next_big_idx bump and go straight to the tail claim below.
				local.big_pool_exhausted = true;
			}
		}
		if (!claimed) {
			const idx_t q = gs.next_tail_idx.fetch_add(1, std::memory_order_relaxed);
			if (q < gs.tail_shards.size()) {
				local.current_shard_idx = gs.tail_shards[q];
				local.current_is_big = false;
				pool = &gs.tail_shards;
				pool_pos = q;
				ahead = gs.prefetch_ahead_tail;
				claimed = true;
			}
		}
		if (!claimed) {
			output.SetCardinality(0);
			return;
		}
		// Warm the index of the shard `ahead` positions later IN THE SAME POOL into
		// page cache while this worker loads+aligns its own, so when a worker claims
		// that position its index load isn't a cold fault. `pool_pos` is unique per
		// atomic claim within a pool and `ahead` is constant, so each target is
		// warmed at most once. The first `ahead` positions of each pool have no
		// earlier claim to warm them, so they load cold — unavoidable at the start.
		// 0 disables; pure CPU-free kernel readahead, output-invariant.
		if (ahead > 0) {
			const idx_t target = pool_pos + ahead;
			if (target < pool->size()) {
				PrefetchShardIndexFiles(gs.shard_index_files[(*pool)[target]]);
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
	tf.named_parameters["oversubscribe"] = LogicalType::INTEGER;
	tf.named_parameters["big_shard_min_batches"] = LogicalType::INTEGER;
	tf.named_parameters["include_shard_name"] = LogicalType::BOOLEAN;
	tf.named_parameters["submit_batch_reads"] = LogicalType::INTEGER;
	tf.named_parameters["prefetch_ahead"] = LogicalType::INTEGER;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignBowtie2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
