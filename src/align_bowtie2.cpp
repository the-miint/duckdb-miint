#include "align_bowtie2.hpp"
#include "align_bowtie2_daemon_common.hpp"
#include "sequence_table_reader.hpp"

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

// Schema/output helpers, ConfigJsonBuilder, ValueAs* coercers, and the
// shared bowtie2-align param mapper all live in
// align_bowtie2_daemon_common.{hpp,cpp}. This file only keeps the
// align_bowtie2-specific glue: subject loading, bowtie2-build invocation,
// query streaming, and the wrapper that decides which params are valid for
// the non-sharded SQL surface.

// User-facing parameter set for align_bowtie2: the common bowtie2-align
// knobs (preset, local, max_secondary, quiet, plus the typed daemon knobs
// added via the migration — see bt2_daemon::kCommonAlignParams) plus the
// miint-side `threads` knob (mapped to nthreads). Anything outside this
// set is rejected at bind time so typos surface at SQL-compile rather than
// running silently with default semantics.
std::unordered_set<std::string> MakeKnownAlignParams() {
	auto s = bt2_daemon::kCommonAlignParams;
	s.insert("threads"); // miint-side; mapped to daemon nthreads below.
	return s;
}

// Build the bowtie2-align config_json. Common bowtie2 knobs flow through
// `bt2_daemon::AppendBowtie2AlignParams`; the non-sharded path additionally
// handles `threads` (its own knob, since the sharded variant warns about it
// instead).
std::string BuildAlignConfigJson(const named_parameter_map_t &named_params, const std::string &index_basename) {
	static const auto kKnown = MakeKnownAlignParams();
	for (const auto &kv : named_params) {
		if (kKnown.find(kv.first) == kKnown.end()) {
			throw InvalidInputException("align_bowtie2: unknown named parameter '%s'.", kv.first);
		}
	}

	bt2_daemon::ConfigJsonBuilder cfg;
	cfg.append_str("index_path", index_basename);

	auto threads_it = named_params.find("threads");
	if (threads_it != named_params.end() && !threads_it->second.IsNull()) {
		const int64_t n = bt2_daemon::ValueAsInt("align_bowtie2", "threads", threads_it->second);
		if (n < 1) {
			throw InvalidInputException("align_bowtie2: threads must be >= 1 (got %lld)", static_cast<long long>(n));
		}
		cfg.append_int("nthreads", n);
	}

	bt2_daemon::AppendBowtie2AlignParams(cfg, named_params, "align_bowtie2");
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

	// Id-column types captured by ValidateSequenceTableSchema in Bind. The
	// daemon's wire schema is always strings — these types drive only the
	// C++/DuckDB output schema and the egress codec dispatch in EmitChunkRows.
	// Default INVALID per the fail-loud convention: any Execute path that
	// runs before Bind populated these throws at the codec dispatch.
	LogicalType query_id_type = LogicalType(LogicalTypeId::INVALID);
	LogicalType subject_id_type = LogicalType(LogicalTypeId::INVALID);
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
// information_schema, and works against views. Validates qual1/qual2 types
// here so a wrong-type column fails fast at bind rather than causing a worker
// crash downstream — the daemon's bowtie2-align input schema wants Phred+33
// Utf8 strings, and miint's canonical raw representation is UTINYINT[] per
// read_fastx. FetchNextQueryBatch converts on the fly.
void DetectQueryColumns(ClientContext &context, AlignBowtie2BindData &bd) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const std::string sql = "SELECT column_name, column_type FROM (DESCRIBE " +
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
			const auto typ = chunk->GetValue(1, i).ToString();
			if (col == "sequence2") {
				bd.query_has_sequence2 = true;
			} else if (col == "qual1") {
				bd.query_has_qual1 = true;
				if (typ != "UTINYINT[]") {
					throw BinderException("align_bowtie2: column 'qual1' in '%s' must be UTINYINT[] (raw Phred values, "
					                      "as produced by read_fastx); got %s",
					                      bd.query_table, typ);
				}
			} else if (col == "qual2") {
				bd.query_has_qual2 = true;
				if (typ != "UTINYINT[]") {
					throw BinderException("align_bowtie2: column 'qual2' in '%s' must be UTINYINT[] (raw Phred values, "
					                      "as produced by read_fastx); got %s",
					                      bd.query_table, typ);
				}
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

	// Validate query + subject tables (existence + read_id type). Both may
	// be VARCHAR or BIGINT; the captured types drive the output schema and
	// the egress codec dispatch in EmitChunkRows. Replaces the earlier
	// DESCRIBE-only existence probes — ValidateSequenceTableSchema does
	// existence + column-type validation in one pass.
	auto query_schema = ValidateSequenceTableSchema(context, bd->query_table, /*allow_bigint=*/true);
	auto subject_schema = ValidateSequenceTableSchema(context, bd->subject_table, /*allow_bigint=*/true);
	bd->query_id_type = query_schema.id_type;
	bd->subject_id_type = subject_schema.id_type;

	DetectQueryColumns(context, *bd);
	bd->named_params = input.named_parameters;

	// Reject unknown named parameters at bind time (BuildAlignConfigJson
	// duplicates this check, but doing it here too keeps the failure mode
	// at SQL-compile rather than at execution).
	static const auto kKnown = MakeKnownAlignParams();
	for (const auto &kv : bd->named_params) {
		if (kKnown.find(kv.first) == kKnown.end()) {
			throw InvalidInputException("align_bowtie2: unknown named parameter '%s'.", kv.first);
		}
	}

	bt2_daemon::PopulateOutputSchema(names, return_types, bd->query_id_type, bd->subject_id_type);
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
	// Fail loud if the resolved daemon predates bowtie2 `memory_mapped` support
	// (>= 0.4.2), since older daemons silently ignore the field.
	bt2_daemon::RequireGplBoundaryVersion(gpl_path, "align_bowtie2");
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

	// 2. Load subjects + materialize as Arrow IPC (shared bowtie2-build glue).
	const auto subjects = bt2_daemon::LoadSingleEndSubjects(context, bd.subject_table, "align_bowtie2");
	const auto subjects_ipc = bt2_daemon::BuildSubjectsIpc(subjects, "align_bowtie2");

	// 3. Build the bowtie2 index under a fresh temp dir.
	gs->temp_dir = MakeTempIndexDir();
	const std::string index_basename = gs->temp_dir + "/idx";

	// nthreads for build comes from the user's `threads` param (same
	// rationale as the old direct-subprocess path: one knob covers both).
	int64_t threads = 1;
	{
		auto it = bd.named_params.find("threads");
		if (it != bd.named_params.end() && !it->second.IsNull()) {
			threads = bt2_daemon::ValueAsInt("align_bowtie2", "threads", it->second);
			if (threads < 1) {
				throw InvalidInputException("align_bowtie2: threads must be >= 1 (got %lld)",
				                            static_cast<long long>(threads));
			}
		}
	}
	const std::string build_config = bt2_daemon::BuildBowtie2BuildConfigJson(index_basename, threads);
	auto build_result = gs->session->Submit("bowtie2-build", build_config, subjects_ipc.data(), subjects_ipc.size());
	gs->index_files = bt2_daemon::ParseBowtie2BuildIndexFiles(build_result.result_json, "align_bowtie2");
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
	std::vector<uint8_t> qual_scratch; // reused across rows to amortize per-row allocations
	std::string qual_encoded;
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
		bt2_daemon::EmitChunkRows(output, to_emit, gs.row_in_batch, batch, bd.query_id_type, bd.subject_id_type);
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
	bt2_daemon::EmitChunkRows(output, to_emit, 0, batch, bd.query_id_type, bd.subject_id_type);
	gs.row_in_batch = to_emit;
}

} // namespace

// =============================================================================
// Public registration
// =============================================================================

TableFunction AlignBowtie2TableFunction::GetFunction() {
	auto tf = TableFunction("align_bowtie2", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind, InitGlobal,
	                        InitLocal);
	bt2_daemon::RegisterBowtie2AlignNamedParameterTypes(tf);
	tf.named_parameters["threads"] = LogicalType::INTEGER; // miint-side; maps to daemon nthreads
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
