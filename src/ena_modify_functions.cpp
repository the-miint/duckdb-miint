// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 MODIFY actions. Currently scoped to
// projects (L4a); samples/experiments/runs follow in subsequent phases.
// See ena_modify_functions.hpp for the contract.

#include "ena_modify_functions.hpp"

#include "ena_client.hpp"
#include "ena_envelope_builder.hpp"
#include "ena_insert_common.hpp"
#include "ena_post_fn.hpp"
#include "ena_projects_insert.hpp"
#include "ena_storage.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace miint {

namespace {

using duckdb::BinderException;
using duckdb::ClientContext;
using duckdb::DataChunk;
using duckdb::ExtensionLoader;
using duckdb::FunctionData;
using duckdb::GlobalTableFunctionState;
using duckdb::InvalidInputException;
using duckdb::LogicalType;
using duckdb::make_uniq;
using duckdb::string;
using duckdb::TableFunction;
using duckdb::TableFunctionBindInput;
using duckdb::TableFunctionData;
using duckdb::TableFunctionInitInput;
using duckdb::TableFunctionInput;
using duckdb::unique_ptr;
using duckdb::Value;
using duckdb::vector;

struct ModifyProjectBindData : public TableFunctionData {
	string secret_name;
	string accession;
	string alias;
	string title;
	string description;
	string project_type;
	bool is_umbrella = false;
	string catalog_name;
	bool catalog_explicit = false;
};

unique_ptr<FunctionData> BindModifyProject(ClientContext &, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<string> &names) {
	auto bd = make_uniq<ModifyProjectBindData>();

	auto get_str = [&](const char *key, string &out, bool required, bool *was_set = nullptr) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("ena_modify_project: required named parameter '%s' is missing", key);
			}
			return;
		}
		out = it->second.ToString();
		if (was_set) {
			*was_set = true;
		}
	};
	auto get_bool = [&](const char *key, bool &out) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			return;
		}
		out = it->second.GetValue<bool>();
	};

	get_str("secret", bd->secret_name, /*required=*/true);
	get_str("accession", bd->accession, /*required=*/true);
	get_str("alias", bd->alias, /*required=*/true);
	// `title` is required on MODIFY: Webin V2 expects "re-submit the full
	// updated XML" and `AppendProject` emits `"title":` unconditionally, so
	// an unset/empty title would silently overwrite the existing title with
	// "" on the server. Catching the omission at bind time is the only
	// place where the "user forgot title vs. user wants empty title"
	// distinction is preserved (both look like `""` to AppendProject).
	get_str("title", bd->title, /*required=*/true);
	get_str("description", bd->description, /*required=*/false);
	get_str("project_type", bd->project_type, /*required=*/false);
	get_bool("is_umbrella", bd->is_umbrella);
	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}

	// Bind-time guards mirror the lifecycle table functions: catch the
	// obvious user errors here so the diagnostic doesn't surface mid-Execute
	// after the secret has been resolved and the envelope half-built.
	if (duckdb::IsENAStringWhitespaceOnly(bd->accession)) {
		throw BinderException("ena_modify_project: 'accession' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->alias)) {
		throw BinderException("ena_modify_project: 'alias' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->title)) {
		throw BinderException("ena_modify_project: 'title' must not be empty or whitespace-only "
		                      "(MODIFY re-submits the full object; an empty title would overwrite "
		                      "the existing one)");
	}

	names = {"action", "target", "success", "alias", "era_accession", "error_messages", "duration_ms"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT};
	return std::move(bd);
}

struct ModifyProjectGlobalState : public GlobalTableFunctionState {
	bool emitted = false;
	miint::ENASubmissionOutcome outcome;
};

// Init only allocates state. The HTTP POST is intentionally deferred to
// Execute: MODIFY has side effects on already-registered records, so a
// scan that produces zero rows (LIMIT 0, WHERE false) MUST not POST. Init
// runs even when the planner knows the scan is empty; Execute runs only
// when a row is actually consumed.
unique_ptr<GlobalTableFunctionState> InitModifyProjectGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ModifyProjectGlobalState>();
}

void ExecuteModifyProject(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gs = data.global_state->Cast<ModifyProjectGlobalState>();
	if (gs.emitted) {
		output.SetCardinality(0);
		return;
	}
	auto &bd = data.bind_data->Cast<ModifyProjectBindData>();

	auto creds = duckdb::ResolveENACredentialsByName(context, "ena_modify_project", bd.secret_name);
	if (creds.endpoint_url.empty()) {
		creds.endpoint_url = duckdb::ResolveDefaultENAEndpointURL(creds.endpoint);
	}

	miint::ProjectSpec spec;
	spec.alias = bd.alias;
	spec.accession = bd.accession;
	spec.title = bd.title;
	spec.description = bd.description;
	spec.project_type = bd.project_type;
	spec.is_umbrella = bd.is_umbrella;

	miint::ENAProjectInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.action = miint::ENAAction::MODIFY;

	miint::ENAClient client(*context.db);
	miint::ENAPostFn post_fn = [&client](const string &url, const string &body, const string &user,
	                                     const string &password, const string &content_type) {
		// Projects use JSON (V2's JSON dispatcher works for project + sample;
		// experiment/run/analysis must be XML, deferred to L4b/c/d). The
		// content-type dispatch matches the INSERT-path post functor so any
		// future routing change lands in one place.
		if (content_type == "application/xml") {
			return client.PostXML(url, body, user, password);
		}
		return client.PostJSONReceiveXML(url, body, user, password);
	};

	std::vector<miint::ProjectSpec> specs = {spec};

	bool transport_failure = false;
	string transport_error;
	try {
		gs.outcome = miint::SubmitProjectInsertOutcome(specs, opts, post_fn);
	} catch (const std::exception &e) {
		transport_failure = true;
		transport_error = e.what();
	}

	const bool success = !transport_failure && gs.outcome.success;
	auto error_messages = gs.outcome.error_messages;
	if (transport_failure && error_messages.empty()) {
		error_messages.push_back(transport_error);
	}
	gs.outcome.error_messages = error_messages;

	if (auto *catalog =
	        duckdb::FindAttachedENACatalog(context, "ena_modify_project", bd.catalog_name, bd.catalog_explicit)) {
		duckdb::SubmissionLogPayload payload;
		payload.object_type = "projects";
		payload.action = miint::ActionName(miint::ENAAction::MODIFY);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = gs.outcome.duration_ms;
		payload.envelope_payload = gs.outcome.envelope_payload;
		payload.raw_receipt = gs.outcome.raw_receipt;
		payload.era_accession = gs.outcome.era_accession;
		payload.error_messages = error_messages;
		// `target` for MODIFY is the accession the user is updating —
		// matches the convention from the lifecycle ops, lets the audit
		// reader filter by accession across action types.
		payload.target = bd.accession;
		// On a successful MODIFY the receipt carries the per-object child;
		// populate from outcome.rows. On failure outcome.rows is empty and
		// the lists stay empty (matches the failed-ADD audit shape).
		payload.object_aliases.reserve(gs.outcome.rows.size());
		payload.object_accessions.reserve(gs.outcome.rows.size());
		for (const auto &row : gs.outcome.rows) {
			payload.object_aliases.push_back(row.alias);
			payload.object_accessions.push_back(row.prjeb_accession);
		}
		duckdb::RecordSubmissionLog(*catalog, creds, payload);
	}

	if (transport_failure) {
		throw InvalidInputException("ena_modify_project: %s", transport_error);
	}
	if (!success) {
		string detail;
		for (const auto &m : error_messages) {
			if (!detail.empty()) {
				detail += "; ";
			}
			detail += m;
		}
		if (detail.empty()) {
			detail = "no error detail";
		}
		throw InvalidInputException("ena_modify_project: %s", detail);
	}

	output.SetCardinality(1);
	output.data[0].SetValue(0, Value(string(miint::ActionName(miint::ENAAction::MODIFY))));
	output.data[1].SetValue(0, Value(bd.accession));
	output.data[2].SetValue(0, Value::BOOLEAN(success));
	output.data[3].SetValue(0, Value(bd.alias));
	output.data[4].SetValue(0, Value(gs.outcome.era_accession));

	vector<Value> errs;
	errs.reserve(error_messages.size());
	for (const auto &m : error_messages) {
		errs.emplace_back(m);
	}
	output.data[5].SetValue(0, Value::LIST(LogicalType::VARCHAR, errs));
	output.data[6].SetValue(0, Value::BIGINT(gs.outcome.duration_ms));

	gs.emitted = true;
}

void AddModifyProjectNamedParameters(TableFunction &tf) {
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["accession"] = LogicalType::VARCHAR;
	tf.named_parameters["alias"] = LogicalType::VARCHAR;
	tf.named_parameters["title"] = LogicalType::VARCHAR;
	tf.named_parameters["description"] = LogicalType::VARCHAR;
	tf.named_parameters["project_type"] = LogicalType::VARCHAR;
	tf.named_parameters["is_umbrella"] = LogicalType::BOOLEAN;
	tf.named_parameters["catalog"] = LogicalType::VARCHAR;
}

} // namespace

void RegisterENAModifyTableFunctions(ExtensionLoader &loader) {
	TableFunction tf("ena_modify_project", {}, ExecuteModifyProject, BindModifyProject, InitModifyProjectGlobal);
	AddModifyProjectNamedParameters(tf);
	loader.RegisterFunction(tf);
}

} // namespace miint
