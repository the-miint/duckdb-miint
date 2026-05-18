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

#include <cstring>
#include <filesystem>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

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

	idx_t max_threads_per_shard = 4;
	bool include_shard_name = false;
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

	idx_t MaxThreads() const override {
		// DuckDB clamps to its own scheduler concurrency anyway, but
		// returning the per-shard ceiling here keeps the planner honest
		// when the thread pool is larger than the work supports.
		return max_active_shards;
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
	for (const auto &shard : bd.shards) {
		if (!HasShardIndex(shard.index_prefix)) {
			throw IOException("No valid bowtie2 index found at prefix: %s. "
			                  "Expected files like %s.1.bt2, %s.rev.1.bt2, etc.",
			                  shard.index_prefix, shard.index_prefix, shard.index_prefix);
		}
	}

	auto gs = make_uniq<AlignBowtie2ShardedGlobalState>();
	// Each miint worker thread orchestrates one shard at a time; the daemon's
	// per-fingerprint pool fans out internally on `nthreads`
	// (= max_threads_per_shard). So our slice of the DuckDB thread pool is
	// `ceil(db_threads / max_threads_per_shard)`, clamped by the shard count
	// (no point asking for more workers than shards). Mirrors the formula
	// `align_minimap2_sharded` uses to honor `SET threads` without a separate
	// knob.
	const idx_t db_threads = NumericCast<idx_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	const idx_t derived = (db_threads + bd.max_threads_per_shard - 1) / bd.max_threads_per_shard;
	gs->max_active_shards = std::max<idx_t>(1, std::min<idx_t>(derived, bd.shards.size()));
	return std::move(gs);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &context, TableFunctionInitInput &,
                                              GlobalTableFunctionState *) {
	auto ls = make_uniq<AlignBowtie2ShardedLocalState>();
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
	if (local.input_stream->HasError()) {
		throw InvalidInputException("align_bowtie2_sharded: failed to open cursor for shard '%s': %s", shard.name,
		                            local.input_stream->GetError());
	}
	local.current_shard_name = shard.name;
}

// Quality decoding lives in align_bowtie2_daemon_common — see
// bt2_daemon::DecodeListQualToPhred33. Shared with align_bowtie2.

// Pull one DuckDB chunk worth of reads from the current shard's stream.
bool FetchShardBatch(AlignBowtie2ShardedLocalState &local, const AlignBowtie2ShardedBindData &bd,
                     bt2_daemon::QueryBatch &out) {
	auto chunk = local.input_stream->Fetch();
	if (!chunk || chunk->size() == 0) {
		return false;
	}
	const idx_t n = chunk->size();
	out.read_ids.clear();
	out.sequence1.clear();
	out.sequence2.clear();
	out.sequence2_valid.clear();
	out.qual1.clear();
	out.qual1_valid.clear();
	out.qual2.clear();
	out.qual2_valid.clear();
	out.read_ids.reserve(n);
	out.sequence1.reserve(n);
	int col = 0;
	const int col_read_id = col++;
	const int col_sequence1 = col++;
	const int col_sequence2 = bd.query_has_sequence2 ? col++ : -1;
	const int col_qual1 = bd.query_has_qual1 ? col++ : -1;
	const int col_qual2 = bd.query_has_qual2 ? col++ : -1;
	std::vector<uint8_t> qual_scratch; // reused across rows to amortize allocs
	std::string qual_encoded;
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
	return true;
}

void SubmitAndDecode(AlignBowtie2ShardedLocalState &local, const AlignBowtie2ShardedBindData &bd,
                     const bt2_daemon::QueryBatch &qb) {
	const auto &shard = bd.shards[local.current_shard_idx];
	const std::string config_json = BuildAlignConfigJson(bd.named_params, shard.index_prefix, bd.max_threads_per_shard);
	bt2_daemon::QueryArrowSchema schema_flags {bd.query_has_sequence2, bd.query_has_qual1, bd.query_has_qual2};
	const auto ipc = bt2_daemon::BuildQueryIpc(qb, schema_flags);
	auto submit_result = local.session->Submit("bowtie2-align", config_json, ipc.data(), ipc.size());
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
	local.current_decoder = std::make_unique<gb::IpcStreamDecoder>(out0.bytes(), out0.size_bytes());
	local.current_decoder->GetSchema(&local.current_schema.arrow_schema);

	if (!local.schema_validated) {
		bt2_daemon::ValidateOutputSchema(local.current_schema.arrow_schema);
		local.schema_validated = true;
	}

	for (;;) {
		ArrowArrayWrapper w;
		if (!local.current_decoder->NextBatch(&w.arrow_array)) {
			break;
		}
		local.current_batches.push_back(std::move(w));
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
			if (FetchShardBatch(local, bd, qb)) {
				SubmitAndDecode(local, bd, qb);
				continue; // loop back to drain decoded rows
			}
			// Shard exhausted — release for re-claim attempt.
			local.input_stream.reset();
			local.current_shard_idx = DConstants::INVALID_INDEX;
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
		OpenCurrentShardStream(local, bd);
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
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignBowtie2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
