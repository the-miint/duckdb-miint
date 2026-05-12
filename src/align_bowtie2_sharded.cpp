#include "align_bowtie2_sharded.hpp"
#include "align_bowtie2_daemon_common.hpp"
#include "miint_log.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/materialized_query_result.hpp"
#include "duckdb/main/query_result.hpp"
#include "duckdb/parser/keyword_helper.hpp"

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
// Config builders. Duplicated from align_bowtie2.cpp's analogous helpers
// because the validators throw with table-function-specific error messages
// ("align_bowtie2_sharded:" vs "align_bowtie2:"). The amount of duplication
// is small enough that a shared parameter-mapper hasn't earned its keep
// yet (KISS — extract on the third caller).
// =============================================================================

struct ConfigJsonBuilder {
	std::ostringstream os;
	bool first = true;
	void append_raw(const std::string &key, const std::string &raw_value) {
		if (!first) {
			os << ",";
		}
		first = false;
		os << "\"" << key << "\":" << raw_value;
	}
	void append_str(const std::string &key, const std::string &value) {
		append_raw(key, "\"" + value + "\"");
	}
	void append_int(const std::string &key, int64_t value) {
		append_raw(key, std::to_string(value));
	}
	void append_bool(const std::string &key, bool value) {
		append_raw(key, value ? "true" : "false");
	}
	std::string build() const {
		return "{" + os.str() + "}";
	}
};

const std::unordered_set<std::string> kKnownShardedParams = {
    "shard_directory", "read_to_shard",         "preset", "local", "threads",
    "max_secondary",   "max_threads_per_shard", "quiet",  "debug", "include_shard_name"};

const std::unordered_set<std::string> kKnownPresets = {"very-fast", "fast", "sensitive", "very-sensitive"};

bool ValueAsBool(const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::BOOLEAN) {
		throw InvalidInputException("align_bowtie2_sharded: parameter '%s' expects a boolean, got %s", name,
		                            v.type().ToString());
	}
	return v.GetValue<bool>();
}

int64_t ValueAsInt(const std::string &name, const Value &v) {
	auto t = v.type().id();
	if (t == LogicalTypeId::INTEGER || t == LogicalTypeId::BIGINT || t == LogicalTypeId::SMALLINT ||
	    t == LogicalTypeId::TINYINT || t == LogicalTypeId::HUGEINT || t == LogicalTypeId::UINTEGER ||
	    t == LogicalTypeId::UBIGINT || t == LogicalTypeId::USMALLINT || t == LogicalTypeId::UTINYINT) {
		return v.GetValue<int64_t>();
	}
	throw InvalidInputException("align_bowtie2_sharded: parameter '%s' expects an integer, got %s", name,
	                            v.type().ToString());
}

std::string ValueAsStr(const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("align_bowtie2_sharded: parameter '%s' expects a string, got %s", name,
		                            v.type().ToString());
	}
	return v.GetValue<std::string>();
}

// =============================================================================
// Shard index discovery — replaces the dependency on
// Bowtie2Aligner::is_index_prefix. Bowtie2 emits 4 mandatory files per index
// in either small (.bt2) or large (.bt2l) format.
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

// =============================================================================
// read_to_shard schema validation — preserves the existing test contract
// ("read_to_shard table missing required column 'read_id'" / 'shard_name').
// We re-implement here rather than pull align_common.hpp's helper because
// the latter is wedded to the direct-subprocess Bowtie2Aligner namespace.
// =============================================================================

struct ShardInfo {
	std::string name;
	std::string index_prefix; // <shard_directory>/<shard_name>/index
	idx_t read_count;
};

void ValidateReadToShardSchema(ClientContext &context, const std::string &table_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const std::string sql = "DESCRIBE " + KeywordHelper::WriteOptionallyQuoted(table_name);
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw BinderException("Table or view '%s' does not exist", table_name);
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	bool has_read_id = false;
	bool has_shard_name = false;
	bool read_id_is_varchar = false;
	bool shard_name_is_varchar = false;
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); ++i) {
			const auto col = chunk->GetValue(0, i).ToString();
			const auto typ = chunk->GetValue(1, i).ToString();
			if (col == "read_id") {
				has_read_id = true;
				read_id_is_varchar = typ == "VARCHAR";
			} else if (col == "shard_name") {
				has_shard_name = true;
				shard_name_is_varchar = typ == "VARCHAR";
			}
		}
	}
	if (!has_read_id) {
		throw BinderException("read_to_shard table '%s' missing required column 'read_id'", table_name);
	}
	if (!read_id_is_varchar) {
		throw BinderException("Column 'read_id' in read_to_shard table '%s' must be VARCHAR", table_name);
	}
	if (!has_shard_name) {
		throw BinderException("read_to_shard table '%s' missing required column 'shard_name'", table_name);
	}
	if (!shard_name_is_varchar) {
		throw BinderException("Column 'shard_name' in read_to_shard table '%s' must be VARCHAR", table_name);
	}
}

// Per-shard read counts, ordered largest-first. Drives the per-shard
// processing order and lets us skip empty shards quickly.
std::vector<ShardInfo> EnumerateShards(ClientContext &context, const std::string &read_to_shard_table,
                                       const std::string &shard_directory) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const std::string sql = "SELECT shard_name, COUNT(*) AS cnt FROM " +
	                        KeywordHelper::WriteOptionallyQuoted(read_to_shard_table) +
	                        " GROUP BY shard_name ORDER BY cnt DESC";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("align_bowtie2_sharded: failed to read shard counts from '%s': %s",
		                            read_to_shard_table, result->GetError());
	}
	auto &materialized = result->Cast<MaterializedQueryResult>();
	std::vector<ShardInfo> out;
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); ++i) {
			auto name_val = chunk->GetValue(0, i);
			auto count_val = chunk->GetValue(1, i);
			if (name_val.IsNull()) {
				throw InvalidInputException("align_bowtie2_sharded: read_to_shard '%s' contains NULL shard_name",
				                            read_to_shard_table);
			}
			ShardInfo info;
			info.name = name_val.GetValue<std::string>();
			info.read_count = static_cast<idx_t>(count_val.GetValue<int64_t>());
			info.index_prefix = shard_directory;
			if (!info.index_prefix.empty() && info.index_prefix.back() != '/') {
				info.index_prefix += '/';
			}
			info.index_prefix += info.name + "/index";
			if (!HasShardIndex(info.index_prefix)) {
				throw BinderException("No valid bowtie2 index found at prefix: %s. "
				                      "Expected files like %s.1.bt2, %s.rev.1.bt2, etc.",
				                      info.index_prefix, info.index_prefix, info.index_prefix);
			}
			out.push_back(std::move(info));
		}
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

	// Ordered shards (largest first) with index path resolved + verified.
	std::vector<ShardInfo> shards;

	idx_t max_threads_per_shard = 4;
	bool include_shard_name = false;
};

struct AlignBowtie2ShardedGlobalState : public GlobalTableFunctionState {
	std::unique_ptr<gb::Session> session;

	// Iteration cursor over shards. We process one shard at a time
	// (Phase 5); Phase 6 may add cross-shard parallelism.
	idx_t shard_index = 0;

	// Streaming cursor over the current shard's matched reads.
	std::unique_ptr<Connection> input_conn;
	std::unique_ptr<QueryResult> input_stream;
	bool current_shard_open = false; // true while input_stream is live for shard_index
	std::string current_shard_name;  // copied for `include_shard_name`

	// Latest Submit response. Same shape as AlignBowtie2GlobalState's;
	// owned via unique_ptr so the ArrowArrayWrapper batches stay valid.
	std::unique_ptr<gb::SubmitResult> current_result;
	std::unique_ptr<gb::IpcStreamDecoder> current_decoder;
	ArrowSchemaWrapper current_schema;
	std::vector<ArrowArrayWrapper> current_batches;
	idx_t batch_index = 0;
	idx_t row_in_batch = 0;
	bool schema_validated = false;

	idx_t MaxThreads() const override {
		return 1; // sequential per shard for Phase 5
	}
};

struct AlignBowtie2ShardedLocalState : public LocalTableFunctionState {};

// Detect optional sequence2/qual1/qual2 columns on the query table.
void DetectQueryColumns(ClientContext &context, AlignBowtie2ShardedBindData &bd) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const std::string sql = "SELECT column_name FROM (DESCRIBE " +
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
			if (col == "sequence2") {
				bd.query_has_sequence2 = true;
			} else if (col == "qual1") {
				bd.query_has_qual1 = true;
			} else if (col == "qual2") {
				bd.query_has_qual2 = true;
			}
		}
	}
}

// Build the bowtie2-align config_json for one shard. Each shard has a
// distinct `index_path` → distinct daemon worker fingerprint, which is
// how Phase 6 will get cross-shard parallelism for free.
std::string BuildAlignConfigJson(const named_parameter_map_t &named_params, const std::string &index_prefix,
                                 idx_t max_threads_per_shard) {
	ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_prefix);
	cfg.append_int("nthreads", static_cast<int64_t>(max_threads_per_shard));

	auto get = [&](const std::string &k) -> const Value * {
		auto it = named_params.find(k);
		return (it == named_params.end() || it->second.IsNull()) ? nullptr : &it->second;
	};

	bool local = false;
	if (auto *v = get("local")) {
		local = ValueAsBool("local", *v);
		cfg.append_bool("local_align", local);
	}
	if (auto *v = get("preset")) {
		const std::string p = ValueAsStr("preset", *v);
		if (kKnownPresets.find(p) == kKnownPresets.end()) {
			throw InvalidInputException("align_bowtie2_sharded: invalid preset '%s' "
			                            "(expected one of very-fast/fast/sensitive/very-sensitive)",
			                            p);
		}
		const std::string preset_value = local ? (p + "-local") : p;
		cfg.append_str("preset", preset_value);
	}
	if (auto *v = get("max_secondary")) {
		const int64_t n = ValueAsInt("max_secondary", *v);
		if (n < 0) {
			throw InvalidInputException("align_bowtie2_sharded: max_secondary must be >= 0 (got %lld)",
			                            static_cast<long long>(n));
		}
		cfg.append_int("k", n);
	}
	bool quiet = true;
	if (auto *v = get("quiet")) {
		quiet = ValueAsBool("quiet", *v);
	}
	cfg.append_bool("verbose", !quiet);

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

	// Validate query_table exists.
	{
		auto &db = DatabaseInstance::GetDatabase(context);
		Connection conn(db);
		auto probe = conn.Query("DESCRIBE " + KeywordHelper::WriteOptionallyQuoted(bd->query_table));
		if (probe->HasError()) {
			throw InvalidInputException("align_bowtie2_sharded: query table '%s' does not exist or is not readable: %s",
			                            bd->query_table, probe->GetError());
		}
	}

	ValidateReadToShardSchema(context, bd->read_to_shard_table);

	// Reject unknown params at bind time.
	for (const auto &kv : input.named_parameters) {
		if (kKnownShardedParams.find(kv.first) == kKnownShardedParams.end()) {
			throw InvalidInputException("align_bowtie2_sharded: unknown named parameter '%s'", kv.first);
		}
	}
	bd->named_params = input.named_parameters;

	// `threads` is ignored at the table-function level — each daemon
	// worker uses bowtie2's internal threading via `nthreads`
	// (max_threads_per_shard). Preserve the existing warning.
	auto threads_param = input.named_parameters.find("threads");
	if (threads_param != input.named_parameters.end() && !threads_param->second.IsNull()) {
		const int64_t threads_val = ValueAsInt("threads", threads_param->second);
		if (threads_val != 1) {
			::miint::EmitWarning(context, "WARNING: Parameter 'threads' is ignored in sharded mode. "
			                              "In sharded mode, each bowtie2 process uses a single thread. "
			                              "Parallelism comes from running multiple processes per shard "
			                              "(max_threads_per_shard) across multiple shards.");
		}
	}

	auto max_tps_param = input.named_parameters.find("max_threads_per_shard");
	if (max_tps_param != input.named_parameters.end() && !max_tps_param->second.IsNull()) {
		const int64_t val = ValueAsInt("max_threads_per_shard", max_tps_param->second);
		if (val < 1 || val > 64) {
			throw BinderException("max_threads_per_shard must be between 1 and 64 (got %lld)",
			                      static_cast<long long>(val));
		}
		bd->max_threads_per_shard = static_cast<idx_t>(val);
	}

	auto include_shard_param = input.named_parameters.find("include_shard_name");
	if (include_shard_param != input.named_parameters.end() && !include_shard_param->second.IsNull()) {
		bd->include_shard_name = ValueAsBool("include_shard_name", include_shard_param->second);
	}

	DetectQueryColumns(context, *bd);

	bd->shards = EnumerateShards(context, bd->read_to_shard_table, bd->shard_directory);

	bt2_daemon::PopulateOutputSchema(names, return_types);
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
	auto gs = make_uniq<AlignBowtie2ShardedGlobalState>();
	gs->session = SpawnAndCheckSession();
	gs->input_conn = std::make_unique<Connection>(DatabaseInstance::GetDatabase(context));
	(void)bd;
	return std::move(gs);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &, TableFunctionInitInput &,
                                              GlobalTableFunctionState *) {
	return make_uniq<AlignBowtie2ShardedLocalState>();
}

// =============================================================================
// Execute — per-shard sequential processing
// =============================================================================

// Open a streaming cursor on the reads matched to the current shard.
void OpenCurrentShardStream(AlignBowtie2ShardedGlobalState &gs, const AlignBowtie2ShardedBindData &bd) {
	const auto &shard = bd.shards[gs.shard_index];
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
	select += " WHERE rts.shard_name = '" + shard.name + "'";

	gs.input_stream = gs.input_conn->SendQuery(select);
	if (gs.input_stream->HasError()) {
		throw InvalidInputException("align_bowtie2_sharded: failed to open cursor for shard '%s': %s", shard.name,
		                            gs.input_stream->GetError());
	}
	gs.current_shard_name = shard.name;
	gs.current_shard_open = true;
}

// Pull one DuckDB chunk worth of reads from the current shard's stream.
bool FetchShardBatch(AlignBowtie2ShardedGlobalState &gs, const AlignBowtie2ShardedBindData &bd,
                     bt2_daemon::QueryBatch &out) {
	auto chunk = gs.input_stream->Fetch();
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
				out.qual1.push_back(v.GetValue<std::string>());
				out.qual1_valid.push_back(1);
			}
		}
		if (col_qual2 >= 0) {
			auto v = chunk->GetValue(col_qual2, i);
			if (v.IsNull()) {
				out.qual2.emplace_back();
				out.qual2_valid.push_back(0);
			} else {
				out.qual2.push_back(v.GetValue<std::string>());
				out.qual2_valid.push_back(1);
			}
		}
	}
	return true;
}

void SubmitAndDecode(AlignBowtie2ShardedGlobalState &gs, const AlignBowtie2ShardedBindData &bd,
                     const bt2_daemon::QueryBatch &qb) {
	const auto &shard = bd.shards[gs.shard_index];
	const std::string config_json = BuildAlignConfigJson(bd.named_params, shard.index_prefix, bd.max_threads_per_shard);
	bt2_daemon::QueryArrowSchema schema_flags {bd.query_has_sequence2, bd.query_has_qual1, bd.query_has_qual2};
	const auto ipc = bt2_daemon::BuildQueryIpc(qb, schema_flags);
	auto submit_result = gs.session->Submit("bowtie2-align", config_json, ipc.data(), ipc.size());
	if (submit_result.outputs.empty()) {
		throw IOException("align_bowtie2_sharded: daemon returned zero shm_outputs for shard '%s'", shard.name);
	}
	if (submit_result.schema_version != bt2_daemon::kAlignSchemaVersion) {
		throw IOException("align_bowtie2_sharded: daemon returned schema_version=%u, expected %u",
		                  submit_result.schema_version, bt2_daemon::kAlignSchemaVersion);
	}

	gs.current_batches.clear();
	if (gs.current_schema.arrow_schema.release) {
		gs.current_schema.arrow_schema.release(&gs.current_schema.arrow_schema);
	}
	gs.current_decoder.reset();
	gs.current_result.reset();
	gs.batch_index = 0;
	gs.row_in_batch = 0;

	gs.current_result = std::make_unique<gb::SubmitResult>(std::move(submit_result));
	const auto &out0 = gs.current_result->outputs[0];
	gs.current_decoder = std::make_unique<gb::IpcStreamDecoder>(out0.bytes(), out0.size_bytes());
	gs.current_decoder->GetSchema(&gs.current_schema.arrow_schema);

	if (!gs.schema_validated) {
		bt2_daemon::ValidateOutputSchema(gs.current_schema.arrow_schema);
		gs.schema_validated = true;
	}

	for (;;) {
		ArrowArrayWrapper w;
		if (!gs.current_decoder->NextBatch(&w.arrow_array)) {
			break;
		}
		gs.current_batches.push_back(std::move(w));
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

	while (true) {
		// 1. Drain any decoded rows from the current Submit response.
		while (gs.batch_index < gs.current_batches.size()) {
			auto &batch = gs.current_batches[gs.batch_index].arrow_array;
			if (batch.offset != 0) {
				throw IOException(
				    "align_bowtie2_sharded: decoded batch has non-zero parent offset (%lld); not yet handled",
				    static_cast<long long>(batch.offset));
			}
			const idx_t total = static_cast<idx_t>(batch.length);
			const idx_t remaining = total - gs.row_in_batch;
			if (remaining == 0) {
				gs.batch_index += 1;
				gs.row_in_batch = 0;
				continue;
			}
			const idx_t to_emit = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
			output.SetCardinality(to_emit);
			bt2_daemon::EmitChunkRows(output, to_emit, gs.row_in_batch, batch);
			if (bd.include_shard_name) {
				FillShardNameColumn(output, to_emit, gs.current_shard_name);
			}
			gs.row_in_batch += to_emit;
			return;
		}

		// 2. Current Submit is drained. Either pull next chunk for current
		//    shard, or advance to next shard.
		if (!gs.current_shard_open) {
			if (gs.shard_index >= bd.shards.size()) {
				output.SetCardinality(0);
				return;
			}
			OpenCurrentShardStream(gs, bd);
		}

		bt2_daemon::QueryBatch qb;
		if (FetchShardBatch(gs, bd, qb)) {
			if (!qb.read_ids.empty()) {
				SubmitAndDecode(gs, bd, qb);
				continue; // loop back to drain
			}
			continue;
		}
		// Current shard exhausted. Advance.
		gs.input_stream.reset();
		gs.current_shard_open = false;
		gs.shard_index += 1;
	}
}

} // namespace

// =============================================================================
// Public registration
// =============================================================================

TableFunction AlignBowtie2ShardedTableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2_sharded", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["shard_directory"] = LogicalType::VARCHAR;
	tf.named_parameters["read_to_shard"] = LogicalType::VARCHAR;
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["local"] = LogicalType::BOOLEAN;
	tf.named_parameters["threads"] = LogicalType::INTEGER;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["max_threads_per_shard"] = LogicalType::INTEGER;
	tf.named_parameters["quiet"] = LogicalType::BOOLEAN;
	tf.named_parameters["debug"] = LogicalType::BOOLEAN;
	tf.named_parameters["include_shard_name"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignBowtie2ShardedTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
