// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 MODIFY actions. Wires the four MODIFY family
// members: ena_modify_project (L4a) + ena_modify_sample (L4b) +
// ena_modify_experiment (L4c) + ena_modify_run (L4d). See
// ena_modify_functions.hpp for the contract.

#include "ena_modify_functions.hpp"

#include "ena_client.hpp"
#include "ena_envelope_builder.hpp"
#include "ena_experiments_insert.hpp"
#include "ena_insert_common.hpp"
#include "ena_post_fn.hpp"
#include "ena_projects_insert.hpp"
#include "ena_runs_insert.hpp"
#include "ena_samples_insert.hpp"
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
		// Projects use JSON (V2's JSON dispatcher handles projects correctly).
		// Samples (post L4b-fix), experiments, and runs go through XML —
		// samples because V2's JSON-to-XML dispatcher misorders DESCRIPTION
		// vs. SAMPLE_NAME (XSD violation), the SRA-side objects because the
		// JSON dispatcher NPEs on them. The content-type dispatch matches
		// the INSERT-path post functor so any future routing change lands
		// in one place.
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

// ===========================================================================
// ena_modify_sample
// ===========================================================================

struct ModifySampleBindData : public TableFunctionData {
	string secret_name;
	string accession;
	string alias;
	int64_t taxon_id = 0;
	string scientific_name;
	string title;
	string description;
	string checklist;
	std::vector<std::pair<std::string, std::string>> attributes;
	std::vector<std::pair<std::string, std::string>> attribute_units;
	string catalog_name;
	bool catalog_explicit = false;
};

unique_ptr<FunctionData> BindModifySample(ClientContext &, TableFunctionBindInput &input,
                                          vector<LogicalType> &return_types, vector<string> &names) {
	auto bd = make_uniq<ModifySampleBindData>();

	auto get_str = [&](const char *key, string &out, bool required, bool *was_set = nullptr) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("ena_modify_sample: required named parameter '%s' is missing", key);
			}
			return;
		}
		out = it->second.ToString();
		if (was_set) {
			*was_set = true;
		}
	};

	get_str("secret", bd->secret_name, /*required=*/true);
	get_str("accession", bd->accession, /*required=*/true);
	get_str("alias", bd->alias, /*required=*/true);

	// taxon_id is required: AppendSample emits "taxonId":"<value>"
	// unconditionally and ValidateSampleSpec rejects taxon_id <= 0. Both an
	// unset (default 0) and an explicit non-positive value would surface as
	// the same envelope-build error mid-Execute; bind-time rejection is
	// faster and clearer.
	{
		auto it = input.named_parameters.find("taxon_id");
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			throw BinderException("ena_modify_sample: required named parameter 'taxon_id' is missing");
		}
		bd->taxon_id = it->second.GetValue<int64_t>();
		if (bd->taxon_id <= 0) {
			// Format the int64_t into a temporary so the message is stable
			// across LP64 / LLP64 platforms (no `%lld` + `long long` cast
			// dance, no `<cinttypes>` PRId64 macro).
			throw BinderException("ena_modify_sample: 'taxon_id' must be > 0 (got " + std::to_string(bd->taxon_id) +
			                      ")");
		}
	}

	get_str("scientific_name", bd->scientific_name, /*required=*/false);
	get_str("title", bd->title, /*required=*/false);
	get_str("description", bd->description, /*required=*/false);
	bool checklist_was_set = false;
	get_str("checklist", bd->checklist, /*required=*/false, &checklist_was_set);
	// Extract MAPs at bind so a malformed entry surfaces as a `BinderException`
	// before secret resolution / network. The named-parameter LogicalType is
	// MAP(VARCHAR, VARCHAR), so DuckDB has already validated the outer shape;
	// `ExtractENAKeyValueMap` catches the per-entry guarantees (struct shape,
	// non-empty key).
	{
		auto it = input.named_parameters.find("attributes");
		if (it != input.named_parameters.end()) {
			bd->attributes = duckdb::ExtractENAKeyValueMap(it->second, "ena_modify_sample", "attributes");
		}
	}
	{
		auto it = input.named_parameters.find("attribute_units");
		if (it != input.named_parameters.end()) {
			bd->attribute_units = duckdb::ExtractENAKeyValueMap(it->second, "ena_modify_sample", "attribute_units");
		}
	}

	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}

	if (duckdb::IsENAStringWhitespaceOnly(bd->accession)) {
		throw BinderException("ena_modify_sample: 'accession' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->alias)) {
		throw BinderException("ena_modify_sample: 'alias' must not be whitespace-only");
	}
	// `checklist` only when the user explicitly set it (passing nothing keeps
	// `bd->checklist` empty, which `AppendSample` correctly omits). Whitespace-
	// only is a different failure mode — it WOULD be emitted as the
	// ENA-CHECKLIST attribute value verbatim, and the server's confusing
	// "checklist '   ' not found" surfaces only after a round-trip.
	if (checklist_was_set && duckdb::IsENAStringWhitespaceOnly(bd->checklist)) {
		throw BinderException("ena_modify_sample: 'checklist' must not be whitespace-only");
	}

	names = {"action", "target", "success", "alias", "era_accession", "error_messages", "duration_ms"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT};
	return std::move(bd);
}

struct ModifySampleGlobalState : public GlobalTableFunctionState {
	bool emitted = false;
	miint::ENASamplesSubmissionOutcome outcome;
};

unique_ptr<GlobalTableFunctionState> InitModifySampleGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ModifySampleGlobalState>();
}

void ExecuteModifySample(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gs = data.global_state->Cast<ModifySampleGlobalState>();
	if (gs.emitted) {
		output.SetCardinality(0);
		return;
	}
	auto &bd = data.bind_data->Cast<ModifySampleBindData>();

	auto creds = duckdb::ResolveENACredentialsByName(context, "ena_modify_sample", bd.secret_name);
	if (creds.endpoint_url.empty()) {
		creds.endpoint_url = duckdb::ResolveDefaultENAEndpointURL(creds.endpoint);
	}

	miint::SampleSpec spec;
	spec.alias = bd.alias;
	spec.accession = bd.accession;
	spec.taxon_id = bd.taxon_id;
	spec.scientific_name = bd.scientific_name;
	spec.title = bd.title;
	spec.description = bd.description;
	spec.checklist = bd.checklist;
	spec.attributes = bd.attributes;
	spec.attribute_units = bd.attribute_units;

	miint::ENASampleInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.action = miint::ENAAction::MODIFY;

	miint::ENAClient client(*context.db);
	miint::ENAPostFn post_fn = [&client](const string &url, const string &body, const string &user,
	                                     const string &password, const string &content_type) {
		if (content_type == "application/xml") {
			return client.PostXML(url, body, user, password);
		}
		return client.PostJSONReceiveXML(url, body, user, password);
	};

	std::vector<miint::SampleSpec> specs = {spec};

	bool transport_failure = false;
	string transport_error;
	try {
		gs.outcome = miint::SubmitSampleInsertOutcome(specs, opts, post_fn);
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
	        duckdb::FindAttachedENACatalog(context, "ena_modify_sample", bd.catalog_name, bd.catalog_explicit)) {
		duckdb::SubmissionLogPayload payload;
		payload.object_type = "samples";
		payload.action = miint::ActionName(miint::ENAAction::MODIFY);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = gs.outcome.duration_ms;
		payload.envelope_payload = gs.outcome.envelope_payload;
		payload.raw_receipt = gs.outcome.raw_receipt;
		payload.era_accession = gs.outcome.era_accession;
		payload.error_messages = error_messages;
		// `target` is the user-supplied accession (audit ground truth),
		// matching the L4a ena_modify_project convention (decision #23).
		payload.target = bd.accession;
		payload.object_aliases.reserve(gs.outcome.rows.size());
		payload.object_accessions.reserve(gs.outcome.rows.size());
		for (const auto &row : gs.outcome.rows) {
			payload.object_aliases.push_back(row.alias);
			payload.object_accessions.push_back(row.ers_accession);
		}
		duckdb::RecordSubmissionLog(*catalog, creds, payload);
	}

	if (transport_failure) {
		throw InvalidInputException("ena_modify_sample: %s", transport_error);
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
		throw InvalidInputException("ena_modify_sample: %s", detail);
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

void AddModifySampleNamedParameters(TableFunction &tf) {
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["accession"] = LogicalType::VARCHAR;
	tf.named_parameters["alias"] = LogicalType::VARCHAR;
	tf.named_parameters["taxon_id"] = LogicalType::BIGINT;
	tf.named_parameters["scientific_name"] = LogicalType::VARCHAR;
	tf.named_parameters["title"] = LogicalType::VARCHAR;
	tf.named_parameters["description"] = LogicalType::VARCHAR;
	tf.named_parameters["checklist"] = LogicalType::VARCHAR;
	tf.named_parameters["attributes"] = LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR);
	tf.named_parameters["attribute_units"] = LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR);
	tf.named_parameters["catalog"] = LogicalType::VARCHAR;
}

// ===========================================================================
// ena_modify_experiment
// ===========================================================================
//
// MODIFY semantics for experiments: re-submit the full updated experiment XML
// identified by ERX accession. The L4a/L4b template applies — bind-time
// require everything `ValidateExperimentSpec` and `AppendXmlExperiment` would
// otherwise surface mid-Execute. Cross-references (`study_ref`,
// `sample_descriptor`) are XSD-mandatory; each is a single VARCHAR that
// accepts either an ENA accession or a parent's alias. The naming and
// disambiguation matches the INSERT path (`INSERT INTO ena.experiments`):
// strings shaped like `<PREFIX><NUMERIC>` (e.g. PRJEB123, SAMEA456) route to
// `RefDescriptor::accession`; everything else is treated as a refname.
// Sharing the disambiguator (`duckdb::ResolveENARefDescriptor`) keeps the
// INSERT/MODIFY user surface symmetric so a user pulling refs from
// `ena.submission_log` doesn't have to think about which form to pass.

struct ModifyExperimentBindData : public TableFunctionData {
	string secret_name;
	string accession;
	string alias;
	string title;
	// Cross-references — single VARCHAR each, matching the INSERT-path column
	// names. Disambiguated to RefDescriptor at Execute time via
	// duckdb::ResolveENARefDescriptor.
	string study_ref;
	string sample_descriptor;
	string design_description;
	string library_name;
	string library_strategy;
	string library_source;
	string library_selection;
	bool library_layout_paired = false;
	bool library_layout_was_set = false;
	string platform;
	string instrument_model;
	string catalog_name;
	bool catalog_explicit = false;
};

unique_ptr<FunctionData> BindModifyExperiment(ClientContext &, TableFunctionBindInput &input,
                                              vector<LogicalType> &return_types, vector<string> &names) {
	auto bd = make_uniq<ModifyExperimentBindData>();

	auto get_str = [&](const char *key, string &out, bool required, bool *was_set = nullptr) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("ena_modify_experiment: required named parameter '%s' is missing", key);
			}
			return;
		}
		out = it->second.ToString();
		if (was_set) {
			*was_set = true;
		}
	};

	get_str("secret", bd->secret_name, /*required=*/true);
	get_str("accession", bd->accession, /*required=*/true);
	get_str("alias", bd->alias, /*required=*/true);
	get_str("title", bd->title, /*required=*/false);
	get_str("study_ref", bd->study_ref, /*required=*/true);
	get_str("sample_descriptor", bd->sample_descriptor, /*required=*/true);
	get_str("design_description", bd->design_description, /*required=*/false);
	get_str("library_name", bd->library_name, /*required=*/false);
	get_str("library_strategy", bd->library_strategy, /*required=*/true);
	get_str("library_source", bd->library_source, /*required=*/true);
	get_str("library_selection", bd->library_selection, /*required=*/true);
	get_str("platform", bd->platform, /*required=*/true);
	get_str("instrument_model", bd->instrument_model, /*required=*/true);

	// library_layout: optional, defaults to SINGLE. Accept "SINGLE" or "PAIRED"
	// (case-insensitive) — that's the SRA.experiment.xsd vocabulary.
	{
		auto it = input.named_parameters.find("library_layout");
		if (it != input.named_parameters.end() && !it->second.IsNull()) {
			bd->library_layout_was_set = true;
			string layout = it->second.ToString();
			std::string upper;
			upper.reserve(layout.size());
			for (char c : layout) {
				upper.push_back(static_cast<char>(c >= 'a' && c <= 'z' ? c - ('a' - 'A') : c));
			}
			if (upper == "PAIRED") {
				bd->library_layout_paired = true;
			} else if (upper == "SINGLE") {
				bd->library_layout_paired = false;
			} else {
				throw BinderException("ena_modify_experiment: 'library_layout' must be 'SINGLE' or 'PAIRED' (got '%s')",
				                      layout.c_str());
			}
		}
	}

	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}

	// Bind-time guards — same pattern as ena_modify_sample. Catch the obvious
	// errors here so the diagnostic doesn't surface mid-Execute after the
	// secret has been resolved and the envelope half-built.
	if (duckdb::IsENAStringWhitespaceOnly(bd->accession)) {
		throw BinderException("ena_modify_experiment: 'accession' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->alias)) {
		throw BinderException("ena_modify_experiment: 'alias' must not be whitespace-only");
	}

	// Whitespace-only ref values would route through ResolveENARefDescriptor
	// as a refname (since whitespace-only doesn't match any accession prefix)
	// and emit `<STUDY_REF refname="   "/>` to the server — a round-trip
	// rejection at best, garbage stored at worst.
	if (duckdb::IsENAStringWhitespaceOnly(bd->study_ref)) {
		throw BinderException("ena_modify_experiment: 'study_ref' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->sample_descriptor)) {
		throw BinderException("ena_modify_experiment: 'sample_descriptor' must not be whitespace-only");
	}

	// Required free-form fields cannot be whitespace-only — AppendXmlElement
	// would emit them verbatim and the server would reject after a round-trip
	// (or worse, accept them and store garbage).
	if (duckdb::IsENAStringWhitespaceOnly(bd->library_strategy)) {
		throw BinderException("ena_modify_experiment: 'library_strategy' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->library_source)) {
		throw BinderException("ena_modify_experiment: 'library_source' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->library_selection)) {
		throw BinderException("ena_modify_experiment: 'library_selection' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->platform)) {
		throw BinderException("ena_modify_experiment: 'platform' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->instrument_model)) {
		throw BinderException("ena_modify_experiment: 'instrument_model' must not be whitespace-only");
	}

	names = {"action", "target", "success", "alias", "era_accession", "error_messages", "duration_ms"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT};
	return std::move(bd);
}

struct ModifyExperimentGlobalState : public GlobalTableFunctionState {
	bool emitted = false;
	miint::ENAExperimentSubmissionOutcome outcome;
};

unique_ptr<GlobalTableFunctionState> InitModifyExperimentGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ModifyExperimentGlobalState>();
}

void ExecuteModifyExperiment(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gs = data.global_state->Cast<ModifyExperimentGlobalState>();
	if (gs.emitted) {
		output.SetCardinality(0);
		return;
	}
	auto &bd = data.bind_data->Cast<ModifyExperimentBindData>();

	auto creds = duckdb::ResolveENACredentialsByName(context, "ena_modify_experiment", bd.secret_name);
	if (creds.endpoint_url.empty()) {
		creds.endpoint_url = duckdb::ResolveDefaultENAEndpointURL(creds.endpoint);
	}

	miint::ExperimentSpec spec;
	spec.alias = bd.alias;
	spec.accession = bd.accession;
	spec.title = bd.title;
	// Disambiguate `study_ref` / `sample_descriptor` strings into
	// accession-vs-refname via the per-kind wrappers in ena_envelope_builder.
	// Same logic the INSERT path uses, so the prefix lists live in one place.
	spec.study_ref = miint::ResolveENAStudyRef(bd.study_ref);
	spec.sample_ref = miint::ResolveENASampleRef(bd.sample_descriptor);
	spec.design_description = bd.design_description;
	spec.library_name = bd.library_name;
	spec.library_strategy = bd.library_strategy;
	spec.library_source = bd.library_source;
	spec.library_selection = bd.library_selection;
	spec.library_layout = bd.library_layout_paired ? miint::ENALibraryLayout::PAIRED : miint::ENALibraryLayout::SINGLE;
	spec.platform = bd.platform;
	spec.instrument_model = bd.instrument_model;

	miint::ENAExperimentInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.action = miint::ENAAction::MODIFY;

	miint::ENAClient client(*context.db);
	miint::ENAPostFn post_fn = [&client](const string &url, const string &body, const string &user,
	                                     const string &password, const string &content_type) {
		// Experiments are XML-only — V2's JSON dispatcher NPEs for SRA-side
		// objects. `ExperimentSubmitTraits::ContentType()` always returns
		// "application/xml", so a non-XML content_type here is a programmer
		// error (a future trait/dispatch refactor that broke the
		// round-trip). Throw an InternalException rather than fall through
		// to a JSON path that would itself fail server-side, so the
		// diagnostic surfaces at the layer where the contract was violated.
		if (content_type != "application/xml") {
			throw duckdb::InternalException(
			    "ena_modify_experiment: post_fn received content_type='%s'; experiments are XML-only", content_type);
		}
		return client.PostXML(url, body, user, password);
	};

	std::vector<miint::ExperimentSpec> specs = {spec};

	bool transport_failure = false;
	string transport_error;
	try {
		gs.outcome = miint::SubmitExperimentInsertOutcome(specs, opts, post_fn);
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
	        duckdb::FindAttachedENACatalog(context, "ena_modify_experiment", bd.catalog_name, bd.catalog_explicit)) {
		duckdb::SubmissionLogPayload payload;
		payload.object_type = "experiments";
		payload.action = miint::ActionName(miint::ENAAction::MODIFY);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = gs.outcome.duration_ms;
		payload.envelope_payload = gs.outcome.envelope_payload;
		payload.raw_receipt = gs.outcome.raw_receipt;
		payload.era_accession = gs.outcome.era_accession;
		payload.error_messages = error_messages;
		// `target` is the user-supplied accession (audit ground truth,
		// decision #23 — same convention as L4a/L4b).
		payload.target = bd.accession;
		payload.object_aliases.reserve(gs.outcome.rows.size());
		payload.object_accessions.reserve(gs.outcome.rows.size());
		for (const auto &row : gs.outcome.rows) {
			payload.object_aliases.push_back(row.alias);
			payload.object_accessions.push_back(row.erx_accession);
		}
		duckdb::RecordSubmissionLog(*catalog, creds, payload);
	}

	if (transport_failure) {
		throw InvalidInputException("ena_modify_experiment: %s", transport_error);
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
		throw InvalidInputException("ena_modify_experiment: %s", detail);
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

void AddModifyExperimentNamedParameters(TableFunction &tf) {
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["accession"] = LogicalType::VARCHAR;
	tf.named_parameters["alias"] = LogicalType::VARCHAR;
	tf.named_parameters["title"] = LogicalType::VARCHAR;
	// Cross-references: matches the `INSERT INTO ena.experiments` column
	// names. Single VARCHAR each; accession-vs-refname disambiguated at
	// Execute by `ResolveENARefDescriptor`.
	tf.named_parameters["study_ref"] = LogicalType::VARCHAR;
	tf.named_parameters["sample_descriptor"] = LogicalType::VARCHAR;
	tf.named_parameters["design_description"] = LogicalType::VARCHAR;
	tf.named_parameters["library_name"] = LogicalType::VARCHAR;
	tf.named_parameters["library_strategy"] = LogicalType::VARCHAR;
	tf.named_parameters["library_source"] = LogicalType::VARCHAR;
	tf.named_parameters["library_selection"] = LogicalType::VARCHAR;
	tf.named_parameters["library_layout"] = LogicalType::VARCHAR;
	tf.named_parameters["platform"] = LogicalType::VARCHAR;
	tf.named_parameters["instrument_model"] = LogicalType::VARCHAR;
	tf.named_parameters["catalog"] = LogicalType::VARCHAR;
}

// ===========================================================================
// ena_modify_run
// ===========================================================================
//
// MODIFY semantics for runs: re-submit the full updated run XML identified by
// ERR accession. The L4a/L4b/L4c template applies. `files` is
// `LIST<STRUCT<filename VARCHAR, filetype VARCHAR, md5 VARCHAR>>`,
// extracted at bind via `duckdb::ExtractENARunFilesList` so the per-entry
// invariants (non-empty filename/filetype/checksum, no NULL entries) surface
// as `BinderException`s before secret resolution. `experiment_ref` matches the
// INSERT-path column name; single VARCHAR each, accession-vs-refname
// disambiguated at Execute via `miint::ResolveENAExperimentRef` (decision #36).

struct ModifyRunBindData : public TableFunctionData {
	string secret_name;
	string accession;
	string alias;
	string title;
	string experiment_ref;
	std::vector<miint::RunFile> files;
	string catalog_name;
	bool catalog_explicit = false;
};

unique_ptr<FunctionData> BindModifyRun(ClientContext &, TableFunctionBindInput &input,
                                       vector<LogicalType> &return_types, vector<string> &names) {
	auto bd = make_uniq<ModifyRunBindData>();

	auto get_str = [&](const char *key, string &out, bool required, bool *was_set = nullptr) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("ena_modify_run: required named parameter '%s' is missing", key);
			}
			return;
		}
		out = it->second.ToString();
		if (was_set) {
			*was_set = true;
		}
	};

	get_str("secret", bd->secret_name, /*required=*/true);
	get_str("accession", bd->accession, /*required=*/true);
	get_str("alias", bd->alias, /*required=*/true);
	get_str("title", bd->title, /*required=*/false);
	get_str("experiment_ref", bd->experiment_ref, /*required=*/true);

	// `files` is required: AppendXmlRun emits at least one <FILE/> child and
	// ValidateRunSpec rejects an empty list. Catch the omission at bind so
	// the diagnostic doesn't surface mid-Execute after secret resolution.
	{
		auto it = input.named_parameters.find("files");
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			throw BinderException("ena_modify_run: required named parameter 'files' is missing");
		}
		bd->files = duckdb::ExtractENARunFilesList(it->second, "ena_modify_run");
		if (bd->files.empty()) {
			throw BinderException("ena_modify_run: 'files' must contain at least one entry");
		}
	}

	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}

	if (duckdb::IsENAStringWhitespaceOnly(bd->accession)) {
		throw BinderException("ena_modify_run: 'accession' must not be whitespace-only");
	}
	if (duckdb::IsENAStringWhitespaceOnly(bd->alias)) {
		throw BinderException("ena_modify_run: 'alias' must not be whitespace-only");
	}
	// Whitespace-only experiment_ref would route through ResolveENAExperimentRef
	// as a refname (whitespace doesn't match `ERX<digits>`) and emit
	// `<EXPERIMENT_REF refname="   "/>` to the server.
	if (duckdb::IsENAStringWhitespaceOnly(bd->experiment_ref)) {
		throw BinderException("ena_modify_run: 'experiment_ref' must not be whitespace-only");
	}

	names = {"action", "target", "success", "alias", "era_accession", "error_messages", "duration_ms"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT};
	return std::move(bd);
}

struct ModifyRunGlobalState : public GlobalTableFunctionState {
	bool emitted = false;
	miint::ENARunSubmissionOutcome outcome;
};

unique_ptr<GlobalTableFunctionState> InitModifyRunGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<ModifyRunGlobalState>();
}

void ExecuteModifyRun(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gs = data.global_state->Cast<ModifyRunGlobalState>();
	if (gs.emitted) {
		output.SetCardinality(0);
		return;
	}
	auto &bd = data.bind_data->Cast<ModifyRunBindData>();

	auto creds = duckdb::ResolveENACredentialsByName(context, "ena_modify_run", bd.secret_name);
	if (creds.endpoint_url.empty()) {
		creds.endpoint_url = duckdb::ResolveDefaultENAEndpointURL(creds.endpoint);
	}

	miint::RunSpec spec;
	spec.alias = bd.alias;
	spec.accession = bd.accession;
	spec.title = bd.title;
	spec.experiment_ref = miint::ResolveENAExperimentRef(bd.experiment_ref);
	spec.files = bd.files;

	miint::ENARunInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.action = miint::ENAAction::MODIFY;

	miint::ENAClient client(*context.db);
	miint::ENAPostFn post_fn = [&client](const string &url, const string &body, const string &user,
	                                     const string &password, const string &content_type) {
		// Runs are XML-only — V2's JSON dispatcher NPEs for SRA-side objects.
		// `RunSubmitTraits::ContentType()` always returns "application/xml",
		// so a non-XML content_type here indicates a future trait/dispatch
		// refactor that broke the contract; throw at the violation point
		// rather than route to a JSON path that itself fails server-side.
		// Mirrors decision #37 (ena_modify_experiment).
		if (content_type != "application/xml") {
			throw duckdb::InternalException("ena_modify_run: post_fn received content_type='%s'; runs are XML-only",
			                                content_type);
		}
		return client.PostXML(url, body, user, password);
	};

	std::vector<miint::RunSpec> specs = {spec};

	bool transport_failure = false;
	string transport_error;
	try {
		gs.outcome = miint::SubmitRunInsertOutcome(specs, opts, post_fn);
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
	        duckdb::FindAttachedENACatalog(context, "ena_modify_run", bd.catalog_name, bd.catalog_explicit)) {
		duckdb::SubmissionLogPayload payload;
		payload.object_type = "runs";
		payload.action = miint::ActionName(miint::ENAAction::MODIFY);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = gs.outcome.duration_ms;
		payload.envelope_payload = gs.outcome.envelope_payload;
		payload.raw_receipt = gs.outcome.raw_receipt;
		payload.era_accession = gs.outcome.era_accession;
		payload.error_messages = error_messages;
		// `target` is the user-supplied accession (audit ground truth, decision #23).
		payload.target = bd.accession;
		payload.object_aliases.reserve(gs.outcome.rows.size());
		payload.object_accessions.reserve(gs.outcome.rows.size());
		for (const auto &row : gs.outcome.rows) {
			payload.object_aliases.push_back(row.alias);
			payload.object_accessions.push_back(row.err_accession);
		}
		duckdb::RecordSubmissionLog(*catalog, creds, payload);
	}

	if (transport_failure) {
		throw InvalidInputException("ena_modify_run: %s", transport_error);
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
		throw InvalidInputException("ena_modify_run: %s", detail);
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

void AddModifyRunNamedParameters(TableFunction &tf) {
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["accession"] = LogicalType::VARCHAR;
	tf.named_parameters["alias"] = LogicalType::VARCHAR;
	tf.named_parameters["title"] = LogicalType::VARCHAR;
	// Cross-reference: matches the `INSERT INTO ena.runs` column name. Single
	// VARCHAR; accession-vs-refname disambiguated at Execute by
	// `ResolveENAExperimentRef`.
	tf.named_parameters["experiment_ref"] = LogicalType::VARCHAR;
	// `files` is LIST<STRUCT<filename, filetype, md5>>. Field names match the
	// `INSERT INTO ena.runs.files` column shape so users can reuse the
	// `ena_upload_reads` RETURNING projection (or anything else built from
	// the INSERT-path schema) in `ena_modify_run` calls verbatim. DuckDB
	// validates the outer LogicalType at bind; per-entry invariants
	// (non-empty fields, no NULL entries) caught by ExtractENARunFilesList.
	tf.named_parameters["files"] = LogicalType::LIST(LogicalType::STRUCT(
	    {{"filename", LogicalType::VARCHAR}, {"filetype", LogicalType::VARCHAR}, {"md5", LogicalType::VARCHAR}}));
	tf.named_parameters["catalog"] = LogicalType::VARCHAR;
}

} // namespace

void RegisterENAModifyTableFunctions(ExtensionLoader &loader) {
	{
		TableFunction tf("ena_modify_project", {}, ExecuteModifyProject, BindModifyProject, InitModifyProjectGlobal);
		AddModifyProjectNamedParameters(tf);
		loader.RegisterFunction(tf);
	}
	{
		TableFunction tf("ena_modify_sample", {}, ExecuteModifySample, BindModifySample, InitModifySampleGlobal);
		AddModifySampleNamedParameters(tf);
		loader.RegisterFunction(tf);
	}
	{
		TableFunction tf("ena_modify_experiment", {}, ExecuteModifyExperiment, BindModifyExperiment,
		                 InitModifyExperimentGlobal);
		AddModifyExperimentNamedParameters(tf);
		loader.RegisterFunction(tf);
	}
	{
		TableFunction tf("ena_modify_run", {}, ExecuteModifyRun, BindModifyRun, InitModifyRunGlobal);
		AddModifyRunNamedParameters(tf);
		loader.RegisterFunction(tf);
	}
}

} // namespace miint
