#include "save_bowtie2_index.hpp"

#include "align_bowtie2_daemon_common.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"

#include <filesystem>
#include <memory>
#include <string>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

// Spawn the daemon and require the bowtie2-build tool. Deliberately does NOT
// call bt2_daemon::RequireGplBoundaryVersion: building an index never uses the
// align-time `memory_mapped` knob (the reason align_bowtie2* gate on >= 0.4.2),
// so index-building stays usable on any daemon that advertises bowtie2-build.
std::unique_ptr<gb::Session> SpawnBuildSession() {
	const std::string gpl_path = gb::FindGplBoundary();
	if (gpl_path.empty()) {
		throw IOException("save_bowtie2_index: gpl-boundary binary not found. To install:\n"
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
	if (!session->has_tool("bowtie2-build")) {
		throw IOException("save_bowtie2_index: gpl-boundary daemon does not advertise bowtie2-build. "
		                  "Upgrade to v0.2.0 or later (`SELECT install_gpl_boundary()`).");
	}
	return session;
}

} // namespace

unique_ptr<FunctionData> SaveBowtie2IndexTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	auto data = make_uniq<Data>();

	if (input.inputs.size() < 2) {
		throw BinderException("save_bowtie2_index requires subject_table and output_path parameters");
	}
	data->subject_table = input.inputs[0].ToString();
	data->output_path = input.inputs[1].ToString();

	// Validate the subject table exists and has the right column types (read_id
	// VARCHAR/BIGINT, sequence1 VARCHAR). Throws BinderException on a missing
	// table or wrong schema — the daemon wire schema is strings, so the returned
	// schema isn't otherwise needed here.
	ValidateSequenceTableSchema(context, data->subject_table, /*allow_bigint=*/true);

	// Carry the named parameters to InitGlobal, where nthreads is resolved against
	// the live DuckDB thread budget (matching align_bowtie2). Resolving there, not
	// here, means a prepared statement re-executed after `SET threads=N` tracks the
	// new count instead of the value captured at prepare time.
	data->named_params = input.named_parameters;

	for (const auto &name : data->names) {
		names.emplace_back(name);
	}
	for (const auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> SaveBowtie2IndexTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Resolve nthreads before doing any work (spawn/load/build): an explicit
	// `threads` wins — validated >= 1, so a bad value still throws before the
	// daemon starts — else default to DuckDB's configured thread budget so the
	// index builds in parallel. Read here (not at bind) so a re-executed prepared
	// statement tracks the current `SET threads`.
	const int64_t db_threads = static_cast<int64_t>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	const int64_t nthreads = bt2_daemon::ResolveNthreadsFromParams(data.named_params, db_threads, "save_bowtie2_index");

	// 1. Spawn the daemon + require bowtie2-build (no >= 0.4.2 gate).
	auto session = SpawnBuildSession();

	// 2. Load subjects (rejects empty / paired) + encode as Arrow IPC.
	const auto subjects = bt2_daemon::LoadSingleEndSubjects(context, data.subject_table, "save_bowtie2_index");
	const auto subjects_ipc = bt2_daemon::BuildSubjectsIpc(subjects, "save_bowtie2_index");

	// 3. Ensure the parent directory of the index basename exists, then build
	//    the index directly at the user's path (kept, unlike align_bowtie2 which
	//    builds into a temp dir and unlinks it).
	std::filesystem::path index_path(data.output_path);
	if (index_path.has_parent_path()) {
		std::error_code ec;
		std::filesystem::create_directories(index_path.parent_path(), ec);
		if (ec) {
			throw IOException("save_bowtie2_index: failed to create directory '%s': %s",
			                  index_path.parent_path().string(), ec.message());
		}
	}

	const std::string build_config = bt2_daemon::BuildBowtie2BuildConfigJson(data.output_path, nthreads);
	auto build_result = session->Submit("bowtie2-build", build_config, subjects_ipc.data(), subjects_ipc.size());
	const auto index_files = bt2_daemon::ParseBowtie2BuildIndexFiles(build_result.result_json, "save_bowtie2_index");
	if (index_files.empty()) {
		throw IOException("save_bowtie2_index: bowtie2-build returned zero index_files");
	}

	gstate->num_subjects = static_cast<int64_t>(subjects.names.size());

	// The daemon has finished writing the .bt2 files; release it.
	session->Shutdown();

	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState>
SaveBowtie2IndexTableFunction::InitLocal(ExecutionContext &, TableFunctionInitInput &, GlobalTableFunctionState *) {
	return make_uniq<LocalState>();
}

void SaveBowtie2IndexTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.done) {
		output.SetCardinality(0);
		return;
	}

	output.data[0].SetValue(0, Value::BOOLEAN(true));               // success
	output.data[1].SetValue(0, Value(bind_data.output_path));       // index_path
	output.data[2].SetValue(0, Value::BIGINT(gstate.num_subjects)); // num_subjects
	output.SetCardinality(1);
	gstate.done = true;
}

TableFunction SaveBowtie2IndexTableFunction::GetFunction() {
	auto tf = TableFunction("save_bowtie2_index", {LogicalType::VARCHAR, LogicalType::VARCHAR}, Execute, Bind,
	                        InitGlobal, InitLocal);
	tf.named_parameters["threads"] = LogicalType::INTEGER;
	return tf;
}

void SaveBowtie2IndexTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
