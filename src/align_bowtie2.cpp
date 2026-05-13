#include "align_bowtie2.hpp"
#include "align_bowtie2_daemon_common.hpp"

#include "duckdb/common/arrow/arrow.hpp"
#include "duckdb/common/arrow/arrow_wrapper.hpp"
#include "duckdb/common/exception.hpp"
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

#include "nanoarrow/nanoarrow.h"
#include "yyjson.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <filesystem>
#include <memory>
#include <mutex>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <unordered_set>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;
namespace yj = duckdb_yyjson;

// Schema/output helpers live in align_bowtie2_daemon_common.{hpp,cpp} so
// align_bowtie2_sharded reuses the same wire-schema-validation,
// EmitChunkRows decoder, and BuildQueryIpc encoder.

// -----------------------------------------------------------------------------
// JSON config builder — mirrors phylogeny_fasttree's ConfigJsonBuilder. Shared
// pattern intentional: a third tool eventually justifies extracting into a
// gpl_boundary helper, but two does not (KISS).
// -----------------------------------------------------------------------------

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
		// Route through gb::JsonEscape: today the only string values are
		// validated presets and a mkdtemp-built index path, but the builder
		// is shared infrastructure and a future user-controlled string is
		// one parameter away. Cheap defense over an audit each time.
		append_raw(key, "\"" + gb::JsonEscape(value) + "\"");
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

// User-facing parameter set for align_bowtie2. Anything outside this set is
// rejected at bind time so typos surface at SQL-compile rather than running
// silently with default semantics.
const std::unordered_set<std::string> kKnownAlignParams = {"preset", "local", "threads", "max_secondary", "quiet"};

// Allowed values for the `preset` parameter; preserved from the direct-
// subprocess interface. The daemon accepts the same set (plus *-local
// variants we trigger via local=true).
const std::unordered_set<std::string> kKnownPresets = {"very-fast", "fast", "sensitive", "very-sensitive"};

bool ValueAsBool(const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::BOOLEAN) {
		throw InvalidInputException("align_bowtie2: parameter '%s' expects a boolean, got %s", name,
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
	throw InvalidInputException("align_bowtie2: parameter '%s' expects an integer, got %s", name, v.type().ToString());
}

std::string ValueAsStr(const std::string &name, const Value &v) {
	if (v.type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("align_bowtie2: parameter '%s' expects a string, got %s", name,
		                            v.type().ToString());
	}
	return v.GetValue<std::string>();
}

// Build the bowtie2-align config_json. The daemon's typed schema is
// comprehensive (~36 knobs); we expose the small subset miint's SQL surface
// used in the direct-subprocess era. `extra_args` is intentionally dropped
// (no test or doc used it).
std::string BuildAlignConfigJson(const named_parameter_map_t &named_params, const std::string &index_basename) {
	ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_basename);

	for (const auto &kv : named_params) {
		if (kKnownAlignParams.find(kv.first) == kKnownAlignParams.end()) {
			throw InvalidInputException("align_bowtie2: unknown named parameter '%s'. "
			                            "Supported: preset, local, threads, max_secondary, quiet.",
			                            kv.first);
		}
	}

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
			throw InvalidInputException("align_bowtie2: invalid preset '%s' "
			                            "(expected one of very-fast/fast/sensitive/very-sensitive)",
			                            p);
		}
		// Daemon understands `*-local` preset variants; compose if user
		// asked for both. If only `local=true` (no preset), the daemon's
		// local_align field is enough.
		const std::string preset_value = local ? (p + "-local") : p;
		cfg.append_str("preset", preset_value);
	}
	if (auto *v = get("threads")) {
		const int64_t n = ValueAsInt("threads", *v);
		if (n < 1) {
			throw InvalidInputException("align_bowtie2: threads must be >= 1 (got %lld)", static_cast<long long>(n));
		}
		cfg.append_int("nthreads", n);
	}
	if (auto *v = get("max_secondary")) {
		const int64_t n = ValueAsInt("max_secondary", *v);
		if (n < 0) {
			throw InvalidInputException("align_bowtie2: max_secondary must be >= 0 (got %lld)",
			                            static_cast<long long>(n));
		}
		// Daemon's `k` field: number of alignments to report.
		// Direct-subprocess parity: max_secondary=0 ⇒ default single
		// alignment per read; max_secondary=N ⇒ report up to N.
		cfg.append_int("k", n);
	}
	// quiet=true (default) ⇒ verbose=false on the daemon. quiet=false flips.
	bool quiet = true;
	if (auto *v = get("quiet")) {
		quiet = ValueAsBool("quiet", *v);
	}
	cfg.append_bool("verbose", !quiet);

	return cfg.build();
}

std::string BuildBuildConfigJson(const std::string &index_basename, int64_t nthreads) {
	ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_basename);
	if (nthreads > 1) {
		cfg.append_int("nthreads", nthreads);
	}
	return cfg.build();
}

// -----------------------------------------------------------------------------
// Temp-dir for bowtie2-build outputs. Respects $TMPDIR (default /tmp). One
// dir per query; cleaned up in GlobalState destructor.
// -----------------------------------------------------------------------------

std::string MakeTempIndexDir() {
	const char *tmp = std::getenv("TMPDIR");
	if (!tmp || !*tmp) {
		tmp = "/tmp";
	}
	std::string tmpl = std::string(tmp) + "/miint-bt2-XXXXXX";
	std::vector<char> buf(tmpl.begin(), tmpl.end());
	buf.push_back('\0');
	if (::mkdtemp(buf.data()) == nullptr) {
		throw IOException("align_bowtie2: failed to create temp dir for bowtie2-build under '%s' (errno=%d)",
		                  std::string(tmp), errno);
	}
	return std::string(buf.data());
}

// -----------------------------------------------------------------------------
// Read subjects via separate connection. Reuses the existing miint contract:
// subjects table has (read_id, sequence1); sequence2 if present must be
// all-NULL. The daemon's bowtie2-build expects (name, sequence), so we
// rename when materializing the Arrow batch.
// -----------------------------------------------------------------------------

struct LoadedSubjects {
	std::vector<std::string> names;
	std::vector<std::string> sequences;
};

LoadedSubjects LoadSubjects(ClientContext &context, const std::string &table_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	// Pull sequence2 too — we have to validate it's all-NULL to match the
	// direct-subprocess contract (test/sql/align_bowtie2.test:191-200 expects
	// "non-NULL sequence2" error for paired subjects). Use COLUMNS to
	// gracefully handle missing sequence2 column.
	const std::string sql =
	    "SELECT read_id, sequence1, "
	    "       COALESCE(CAST(NULL AS VARCHAR), NULL) AS sequence2_unused FROM " // placeholder, replaced below
	    + KeywordHelper::WriteOptionallyQuoted(table_name);
	(void)sql; // not used — we use the two-query approach below

	// Two-query approach: first probe column presence with information_schema,
	// then issue the actual SELECT. Simpler than COLUMNS hacks.
	const std::string columns_sql = "SELECT column_name FROM (DESCRIBE " +
	                                KeywordHelper::WriteOptionallyQuoted(table_name) +
	                                ") WHERE column_name IN ('sequence2')";
	auto columns_res = conn.Query(columns_sql);
	if (columns_res->HasError()) {
		throw InvalidInputException("align_bowtie2: failed to introspect subject table '%s': %s", table_name,
		                            columns_res->GetError());
	}
	const bool has_sequence2 = columns_res->RowCount() > 0;

	std::string select_sql = "SELECT read_id, sequence1";
	if (has_sequence2) {
		select_sql += ", sequence2";
	}
	select_sql += " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);

	auto result = conn.Query(select_sql);
	if (result->HasError()) {
		throw InvalidInputException("align_bowtie2: failed to read subject table '%s' "
		                            "(must have columns read_id VARCHAR, sequence1 VARCHAR): %s",
		                            table_name, result->GetError());
	}

	LoadedSubjects out;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	out.names.reserve(materialized.RowCount());
	out.sequences.reserve(materialized.RowCount());
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); ++i) {
			auto name_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			if (name_val.IsNull() || seq_val.IsNull()) {
				throw InvalidInputException("align_bowtie2: NULL read_id or sequence1 in subject table '%s'",
				                            table_name);
			}
			if (has_sequence2) {
				auto s2_val = chunk->GetValue(2, i);
				if (!s2_val.IsNull()) {
					throw InvalidInputException("align_bowtie2: subject table '%s' has non-NULL sequence2 — "
					                            "subjects must be single-end (sequence2 must be NULL)",
					                            table_name);
				}
			}
			out.names.push_back(name_val.GetValue<std::string>());
			out.sequences.push_back(seq_val.GetValue<std::string>());
		}
	}
	if (out.names.empty()) {
		throw InvalidInputException("align_bowtie2: subject table '%s' is empty", table_name);
	}
	return out;
}

// -----------------------------------------------------------------------------
// Arrow IPC encoders. Subjects ⇒ {name, sequence}; queries ⇒ {read_id,
// sequence1, sequence2?, qual1?, qual2?}. Mirrors phylogeny_fasttree's
// BuildInputIpcStream (nanoarrow's Init/StartAppending/AppendString/
// FinishElement/FinishBuildingDefault + EncodeIpcStream).
// -----------------------------------------------------------------------------

std::vector<uint8_t> BuildSubjectsIpc(const LoadedSubjects &subjects) {
	ArrowSchema schema {};
	auto rc = ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT);
	if (rc != NANOARROW_OK) {
		throw InternalException("align_bowtie2: ArrowSchemaInit failed");
	}
	rc = ArrowSchemaAllocateChildren(&schema, 2);
	if (rc != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowSchemaAllocateChildren failed");
	}
	ArrowSchemaInitFromType(schema.children[0], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[0], "name");
	ArrowSchemaInitFromType(schema.children[1], NANOARROW_TYPE_STRING);
	ArrowSchemaSetName(schema.children[1], "sequence");

	ArrowArray array {};
	ArrowError err {};
	if (ArrowArrayInitFromSchema(&array, &schema, &err) != NANOARROW_OK) {
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayInit failed: %s", err.message);
	}
	if (ArrowArrayStartAppending(&array) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayStartAppending failed");
	}
	for (size_t i = 0; i < subjects.names.size(); ++i) {
		ArrowStringView nv {subjects.names[i].data(), static_cast<int64_t>(subjects.names[i].size())};
		ArrowStringView sv {subjects.sequences[i].data(), static_cast<int64_t>(subjects.sequences[i].size())};
		if (ArrowArrayAppendString(array.children[0], nv) != NANOARROW_OK ||
		    ArrowArrayAppendString(array.children[1], sv) != NANOARROW_OK ||
		    ArrowArrayFinishElement(&array) != NANOARROW_OK) {
			array.release(&array);
			schema.release(&schema);
			throw InternalException("align_bowtie2: subjects append failed at row %d", static_cast<int>(i));
		}
	}
	if (ArrowArrayFinishBuildingDefault(&array, &err) != NANOARROW_OK) {
		array.release(&array);
		schema.release(&schema);
		throw InternalException("align_bowtie2: ArrowArrayFinishBuilding failed: %s", err.message);
	}

	std::vector<uint8_t> bytes;
	try {
		bytes = gb::EncodeIpcStream(&schema, &array, 1);
	} catch (...) {
		array.release(&array);
		schema.release(&schema);
		throw;
	}
	if (array.release) {
		array.release(&array);
	}
	if (schema.release) {
		schema.release(&schema);
	}
	return bytes;
}

// -----------------------------------------------------------------------------
// Parse bowtie2-build's `result.index_files` array into a vector of paths.
// Returns paths in registry order; we use the first one to derive the index
// basename. Daemon contract: 6 absolute paths.
// -----------------------------------------------------------------------------

std::vector<std::string> ParseIndexFilesFromResultJson(const std::string &result_json) {
	using YyjsonDocPtr = std::unique_ptr<yj::yyjson_doc, decltype(&yj::yyjson_doc_free)>;
	YyjsonDocPtr doc(yj::yyjson_read(result_json.data(), result_json.size(), 0), &yj::yyjson_doc_free);
	if (!doc) {
		throw IOException("align_bowtie2: bowtie2-build result_json was not valid JSON: %s", result_json);
	}
	yj::yyjson_val *root = yj::yyjson_doc_get_root(doc.get());
	if (!yj::yyjson_is_obj(root)) {
		throw IOException("align_bowtie2: bowtie2-build result was not a JSON object: %s", result_json);
	}
	yj::yyjson_val *arr = yj::yyjson_obj_get(root, "index_files");
	if (!arr || !yj::yyjson_is_arr(arr)) {
		throw IOException("align_bowtie2: bowtie2-build result missing 'index_files' array: %s", result_json);
	}
	std::vector<std::string> out;
	const size_t n = yj::yyjson_arr_size(arr);
	out.reserve(n);
	for (size_t i = 0; i < n; ++i) {
		yj::yyjson_val *item = yj::yyjson_arr_get(arr, i);
		if (!yj::yyjson_is_str(item)) {
			throw IOException("align_bowtie2: bowtie2-build index_files contained non-string entry");
		}
		out.emplace_back(yj::yyjson_get_str(item));
	}
	return out;
}

// -----------------------------------------------------------------------------
// Bind / GlobalState / Execute
// -----------------------------------------------------------------------------

struct AlignBowtie2BindData : public TableFunctionData {
	std::string query_table;
	std::string subject_table;
	// Carry the user's named_params forward so InitGlobal can build the
	// align config_json once it knows the index basename.
	named_parameter_map_t named_params;

	// Detected at bind time, used by Execute when building the query Arrow
	// batches. Auto-detection of paired-end is per-batch on the daemon
	// side, but only columns that exist in the source can be emitted.
	bool query_has_sequence2 = false;
	bool query_has_qual1 = false;
	bool query_has_qual2 = false;
};

struct AlignBowtie2GlobalState : public GlobalTableFunctionState {
	std::unique_ptr<gb::Session> session;
	std::string config_json_align;

	// Index lifetime — populated after bowtie2-build returns, unlinked in
	// the destructor (below Session::Shutdown).
	std::string temp_dir;
	std::vector<std::string> index_files;

	// Streaming input cursor.
	std::unique_ptr<Connection> input_conn;
	std::unique_ptr<QueryResult> input_stream;
	bool input_exhausted = false;
	std::string query_select_sql;

	// Latest Submit response + its decoder. Held as unique_ptrs because
	// SubmitResult / IpcStreamDecoder own resources whose lifetimes must
	// outlive the ArrowArrayWrapper batches inside `current_batches`.
	std::unique_ptr<gb::SubmitResult> current_result;
	std::unique_ptr<gb::IpcStreamDecoder> current_decoder;
	ArrowSchemaWrapper current_schema;
	std::vector<ArrowArrayWrapper> current_batches;
	idx_t batch_index = 0;
	idx_t row_in_batch = 0;
	bool schema_validated = false;

	idx_t MaxThreads() const override {
		return 1;
	}

	~AlignBowtie2GlobalState() override {
		// Shut down the daemon first so it releases its index handles
		// before we unlink the files. POSIX allows unlink-while-open
		// (inode persists until last close), but ordering is cleaner.
		if (session) {
			try {
				session->Shutdown();
			} catch (...) {
				// Destructor must not throw.
			}
			session.reset();
		}
		if (!temp_dir.empty()) {
			// Sweep every entry under the temp dir rather than relying on the
			// `index_files` list. Covers the daemon-crash-mid-build case
			// where partial `.bt2` files were written before the result JSON
			// came back: without this, `rmdir` would silently fail and leak
			// the directory + files.
			try {
				std::error_code ec;
				std::filesystem::remove_all(temp_dir, ec);
				// `ec` is ignored on purpose — destructor can't propagate it
				// and there's no useful action on cleanup failure beyond
				// letting the OS reclaim on reboot.
			} catch (...) {
				// Destructor must not throw.
			}
		}
	}
};

struct AlignBowtie2LocalState : public LocalTableFunctionState {};

// Detect query-table column shape via DESCRIBE. Cheaper than
// information_schema, and works against views.
void DetectQueryColumns(ClientContext &context, AlignBowtie2BindData &bd) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const std::string sql = "SELECT column_name FROM (DESCRIBE " +
	                        KeywordHelper::WriteOptionallyQuoted(bd.query_table) +
	                        ") WHERE column_name IN ('sequence2','qual1','qual2')";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("align_bowtie2: failed to introspect query table '%s': %s", bd.query_table,
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

unique_ptr<FunctionData> Bind(ClientContext &context, TableFunctionBindInput &input, vector<LogicalType> &return_types,
                              vector<std::string> &names) {
	if (input.inputs.size() < 2) {
		throw BinderException("align_bowtie2 requires query_table and subject_table parameters");
	}
	auto bd = make_uniq<AlignBowtie2BindData>();
	bd->query_table = input.inputs[0].ToString();
	bd->subject_table = input.inputs[1].ToString();

	// Validate table existence at bind time by checking columns are
	// describable. Defers heavy reads to InitGlobal.
	auto &db = DatabaseInstance::GetDatabase(context);
	{
		Connection conn(db);
		auto probe = conn.Query("DESCRIBE " + KeywordHelper::WriteOptionallyQuoted(bd->query_table));
		if (probe->HasError()) {
			throw InvalidInputException("align_bowtie2: query table '%s' does not exist or is not readable: %s",
			                            bd->query_table, probe->GetError());
		}
	}
	{
		Connection conn(db);
		auto probe = conn.Query("DESCRIBE " + KeywordHelper::WriteOptionallyQuoted(bd->subject_table));
		if (probe->HasError()) {
			throw InvalidInputException("align_bowtie2: subject table '%s' does not exist or is not readable: %s",
			                            bd->subject_table, probe->GetError());
		}
	}

	DetectQueryColumns(context, *bd);
	bd->named_params = input.named_parameters;

	// Reject unknown named parameters at bind time (BuildAlignConfigJson
	// duplicates this check, but doing it here too keeps the failure mode
	// at SQL-compile rather than at execution).
	for (const auto &kv : bd->named_params) {
		if (kKnownAlignParams.find(kv.first) == kKnownAlignParams.end()) {
			throw InvalidInputException("align_bowtie2: unknown named parameter '%s'. "
			                            "Supported: preset, local, threads, max_secondary, quiet.",
			                            kv.first);
		}
	}

	bt2_daemon::PopulateOutputSchema(names, return_types);
	return std::move(bd);
}

// Open the daemon, validate tool support + schema_version. Hard fail (loud)
// rather than producing garbled output if the daemon doesn't have bowtie2
// or has drifted to a different schema_version.
std::unique_ptr<gb::Session> SpawnAndCheckSession() {
	const std::string gpl_path = gb::FindGplBoundary();
	if (gpl_path.empty()) {
		throw IOException("align_bowtie2: gpl-boundary binary not found. To install:\n"
		                  "  Easiest:  SELECT install_gpl_boundary();   "
		                  "-- downloads a prebuilt binary into miint's cache dir\n"
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
		throw IOException("align_bowtie2: gpl-boundary daemon does not advertise bowtie2-align. "
		                  "Upgrade to v0.2.0 or later (`SELECT install_gpl_boundary()`).");
	}
	if (!session->has_tool("bowtie2-build")) {
		throw IOException("align_bowtie2: gpl-boundary daemon does not advertise bowtie2-build. "
		                  "Upgrade to v0.2.0 or later (`SELECT install_gpl_boundary()`).");
	}
	const uint32_t got = session->tool_schema_version("bowtie2-align");
	if (got != bt2_daemon::kAlignSchemaVersion) {
		throw IOException("align_bowtie2: daemon reports bowtie2-align schema_version=%u, miint expects %u. "
		                  "Mismatched gpl-boundary release.",
		                  got, bt2_daemon::kAlignSchemaVersion);
	}
	return session;
}

unique_ptr<GlobalTableFunctionState> InitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &bd = input.bind_data->Cast<AlignBowtie2BindData>();
	auto gs = make_uniq<AlignBowtie2GlobalState>();

	// 1. Spawn the daemon + handshake. Fail loud if not present or stale.
	gs->session = SpawnAndCheckSession();

	// 2. Load subjects + materialize as Arrow IPC.
	const auto subjects = LoadSubjects(context, bd.subject_table);
	const auto subjects_ipc = BuildSubjectsIpc(subjects);

	// 3. Build the bowtie2 index under a fresh temp dir.
	gs->temp_dir = MakeTempIndexDir();
	const std::string index_basename = gs->temp_dir + "/idx";

	// nthreads for build comes from the user's `threads` param (same
	// rationale as the old direct-subprocess path: one knob covers both).
	int64_t threads = 1;
	{
		auto it = bd.named_params.find("threads");
		if (it != bd.named_params.end() && !it->second.IsNull()) {
			threads = ValueAsInt("threads", it->second);
			if (threads < 1) {
				throw InvalidInputException("align_bowtie2: threads must be >= 1 (got %lld)",
				                            static_cast<long long>(threads));
			}
		}
	}
	const std::string build_config = BuildBuildConfigJson(index_basename, threads);
	auto build_result = gs->session->Submit("bowtie2-build", build_config, subjects_ipc.data(), subjects_ipc.size());
	gs->index_files = ParseIndexFilesFromResultJson(build_result.result_json);
	if (gs->index_files.empty()) {
		throw IOException("align_bowtie2: bowtie2-build returned zero index_files");
	}

	// 4. Now that we have the index basename, finish building the align
	//    config JSON.
	gs->config_json_align = BuildAlignConfigJson(bd.named_params, index_basename);

	// 5. Open a streaming cursor on the query table. SendQuery returns a
	//    StreamQueryResult that fetches chunks lazily.
	gs->input_conn = std::make_unique<Connection>(DatabaseInstance::GetDatabase(context));
	std::string select = "SELECT read_id, sequence1";
	if (bd.query_has_sequence2) {
		select += ", sequence2";
	}
	if (bd.query_has_qual1) {
		select += ", qual1";
	}
	if (bd.query_has_qual2) {
		select += ", qual2";
	}
	select += " FROM " + KeywordHelper::WriteOptionallyQuoted(bd.query_table);
	gs->query_select_sql = select;
	gs->input_stream = gs->input_conn->SendQuery(select);
	if (gs->input_stream->HasError()) {
		throw InvalidInputException("align_bowtie2: failed to open streaming cursor on query table '%s': %s",
		                            bd.query_table, gs->input_stream->GetError());
	}

	return std::move(gs);
}

unique_ptr<LocalTableFunctionState> InitLocal(ExecutionContext &, TableFunctionInitInput &,
                                              GlobalTableFunctionState *) {
	return make_uniq<AlignBowtie2LocalState>();
}

// Pull one DuckDB chunk worth of queries from the streaming cursor. Returns
// true if rows were produced (caller should submit); false if the cursor is
// exhausted (caller emits done).
bool FetchNextQueryBatch(AlignBowtie2GlobalState &gs, const AlignBowtie2BindData &bd, bt2_daemon::QueryBatch &out) {
	if (gs.input_exhausted) {
		return false;
	}
	auto chunk = gs.input_stream->Fetch();
	if (!chunk || chunk->size() == 0) {
		gs.input_exhausted = true;
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
	if (bd.query_has_sequence2) {
		out.sequence2.reserve(n);
		out.sequence2_valid.reserve(n);
	}
	if (bd.query_has_qual1) {
		out.qual1.reserve(n);
		out.qual1_valid.reserve(n);
	}
	if (bd.query_has_qual2) {
		out.qual2.reserve(n);
		out.qual2_valid.reserve(n);
	}
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
			throw InvalidInputException(
			    "align_bowtie2: NULL read_id or sequence1 in query table '%s' (daemon requires both columns non-null)",
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

// Submit a query batch through the daemon and decode the response into
// gs.current_*. Replaces any prior decoded batches.
void SubmitAndDecode(AlignBowtie2GlobalState &gs, const AlignBowtie2BindData &bd, const bt2_daemon::QueryBatch &qb) {
	bt2_daemon::QueryArrowSchema schema_flags {bd.query_has_sequence2, bd.query_has_qual1, bd.query_has_qual2};
	const auto ipc = bt2_daemon::BuildQueryIpc(qb, schema_flags);
	auto submit_result = gs.session->Submit("bowtie2-align", gs.config_json_align, ipc.data(), ipc.size());
	if (submit_result.outputs.empty()) {
		throw IOException("align_bowtie2: daemon returned zero shm_outputs for bowtie2-align batch");
	}
	if (submit_result.schema_version != bt2_daemon::kAlignSchemaVersion) {
		throw IOException("align_bowtie2: daemon returned schema_version=%u, expected %u", submit_result.schema_version,
		                  bt2_daemon::kAlignSchemaVersion);
	}

	// Tear down previous batch state. ArrowArrayWrapper releases reach into
	// the prior SubmitResult's mmap, so reset them before the SubmitResult.
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

	// Validate schema on the first batch only — subsequent batches under
	// the same Submit cycle share the daemon's output_schema.
	if (!gs.schema_validated) {
		bt2_daemon::ValidateOutputSchema(gs.current_schema.arrow_schema);
		gs.schema_validated = true;
	}

	// Drain all batches for this Submit eagerly. bowtie2-align emits one
	// per call in v0.2.0, but the API allows multiple.
	for (;;) {
		ArrowArrayWrapper w;
		if (!gs.current_decoder->NextBatch(&w.arrow_array)) {
			break;
		}
		gs.current_batches.push_back(std::move(w));
	}
}

void Execute(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	(void)context;
	auto &bd = data.bind_data->Cast<AlignBowtie2BindData>();
	auto &gs = data.global_state->Cast<AlignBowtie2GlobalState>();

	// Drain any rows still pending in the current decoded batch first.
	while (gs.batch_index < gs.current_batches.size()) {
		auto &batch = gs.current_batches[gs.batch_index].arrow_array;
		if (batch.offset != 0) {
			throw IOException("align_bowtie2: decoded batch has non-zero parent offset (%lld); not yet handled",
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
		gs.row_in_batch += to_emit;
		return;
	}

	// Current Submit response is fully drained. Pull another query chunk and
	// submit. Loop until we produce rows or input is exhausted (bowtie2 may
	// drop reads that don't align under `no_unal`, though we don't set it;
	// looping defensively is cheap).
	bt2_daemon::QueryBatch qb;
	while (true) {
		if (!FetchNextQueryBatch(gs, bd, qb)) {
			output.SetCardinality(0);
			return;
		}
		if (qb.read_ids.empty()) {
			continue;
		}
		SubmitAndDecode(gs, bd, qb);
		if (!gs.current_batches.empty()) {
			break;
		}
		// Daemon returned zero output rows for this input chunk. Re-loop.
	}

	auto &batch = gs.current_batches[gs.batch_index].arrow_array;
	if (batch.offset != 0) {
		throw IOException("align_bowtie2: decoded batch has non-zero parent offset (%lld); not yet handled",
		                  static_cast<long long>(batch.offset));
	}
	const idx_t total = static_cast<idx_t>(batch.length);
	const idx_t to_emit = MinValue<idx_t>(total, STANDARD_VECTOR_SIZE);
	output.SetCardinality(to_emit);
	bt2_daemon::EmitChunkRows(output, to_emit, 0, batch);
	gs.row_in_batch = to_emit;
}

} // namespace

// =============================================================================
// Public registration
// =============================================================================

TableFunction AlignBowtie2TableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                        InitLocal);
	tf.named_parameters["preset"] = LogicalType::VARCHAR;
	tf.named_parameters["local"] = LogicalType::BOOLEAN;
	tf.named_parameters["threads"] = LogicalType::INTEGER;
	tf.named_parameters["max_secondary"] = LogicalType::INTEGER;
	tf.named_parameters["quiet"] = LogicalType::BOOLEAN;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void AlignBowtie2TableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

// =============================================================================
// bowtie2_available() scalar — probes the daemon's tool registry via
// `gpl-boundary --list-tools` (mirrors PhylogenyFastTreeAvailableImpl).
// =============================================================================

static void Bowtie2AvailableFunction(DataChunk &args, ExpressionState &state, Vector &result) {
	(void)args;
	(void)state;
	static std::once_flag flag;
	static bool available = false;

	std::call_once(flag, []() {
		const std::string path = gb::FindGplBoundary();
		if (path.empty()) {
			available = false;
			return;
		}
		try {
			std::vector<std::string> argv = {path, "--list-tools"};
			gb::ChildProcess child(argv);
			// Drain stdout fully. `--list-tools` writes a small JSON array
			// (a few hundred bytes for the v0.2.0 registry) and exits — well
			// below the pipe-buffer limit, so no risk of the child blocking
			// here. Then drain stderr to keep the child from blocking on a
			// stderr write before exit. `gpl-boundary --list-tools` doesn't
			// write stderr in the success path (see GPL-boundary
			// `src/main.rs:88-94`), but error paths do; the second drain is
			// defensive and bounded.
			auto drain = [](int fd, std::string &out) {
				if (fd < 0) {
					return;
				}
				char buf[256];
				ssize_t n;
				while ((n = ::read(fd, buf, sizeof(buf))) > 0) {
					out.append(buf, static_cast<size_t>(n));
				}
			};
			std::string out;
			std::string err;
			drain(child.stdout_fd(), out);
			drain(child.stderr_fd(), err);
			const int status = child.Wait();
			available = WIFEXITED(status) && WEXITSTATUS(status) == 0 &&
			            out.find("bowtie2-align") != std::string::npos &&
			            out.find("bowtie2-build") != std::string::npos;
		} catch (...) {
			available = false;
		}
	});

	result.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto &constant_validity = ConstantVector::Validity(result);
	constant_validity.SetAllValid(1);
	*ConstantVector::GetData<bool>(result) = available;
}

void RegisterBowtie2AvailableFunction(ExtensionLoader &loader) {
	ScalarFunction bowtie2_available("bowtie2_available", {}, LogicalType::BOOLEAN, Bowtie2AvailableFunction);
	loader.RegisterFunction(bowtie2_available);
}

} // namespace duckdb
