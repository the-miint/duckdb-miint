// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 targeted lifecycle actions: ena_cancel,
// ena_release, ena_hold. See ena_lifecycle_functions.hpp.

#include "ena_lifecycle_functions.hpp"

#include "ena_client.hpp"
#include "ena_envelope_builder.hpp"
#include "ena_insert_common.hpp"
#include "ena_lifecycle_submit.hpp"
#include "ena_storage.hpp"

#include "duckdb/catalog/catalog_transaction.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/attached_database.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/database_manager.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/main/secret/secret.hpp"
#include "duckdb/main/secret/secret_manager.hpp"

namespace miint {

namespace {

using duckdb::AttachedDatabase;
using duckdb::BinderException;
using duckdb::CatalogTransaction;
using duckdb::ClientContext;
using duckdb::DatabaseManager;
using duckdb::DataChunk;
using duckdb::ExtensionLoader;
using duckdb::FlatVector;
using duckdb::FunctionData;
using duckdb::GlobalTableFunctionState;
using duckdb::InvalidInputException;
using duckdb::KeyValueSecret;
using duckdb::ListVector;
using duckdb::LogicalType;
using duckdb::make_uniq;
using duckdb::OperatorResultType;
using duckdb::optional_ptr;
using duckdb::SecretManager;
using duckdb::SourceResultType;
using duckdb::string;
using duckdb::string_t;
using duckdb::StringVector;
using duckdb::TableFunction;
using duckdb::TableFunctionBindInput;
using duckdb::TableFunctionData;
using duckdb::TableFunctionInitInput;
using duckdb::TableFunctionInput;
using duckdb::unique_ptr;
using duckdb::Value;
using duckdb::vector;

// Look up the named ATTACHed database and return it iff its catalog is an
// ENACatalog. Behaviour depends on whether the catalog name was the
// (silent) default "ena" or explicitly requested by the user:
//   - default + missing       → nullptr (one-shot CANCEL without ATTACH)
//   - default + wrong type    → nullptr (something else named "ena" — silent)
//   - explicit + missing      → throw (user expected logging; missing audit gap is invisible if silent)
//   - explicit + wrong type   → throw (same reason)
duckdb::ENACatalog *FindAttachedENACatalog(ClientContext &context, const string &caller, const string &catalog_name,
                                           bool explicit_name) {
	if (catalog_name.empty()) {
		return nullptr;
	}
	auto &db_manager = DatabaseManager::Get(context);
	auto db = db_manager.GetDatabase(context, catalog_name);
	if (!db) {
		if (explicit_name) {
			throw InvalidInputException("%s: catalog '%s' is not attached", caller, catalog_name);
		}
		return nullptr;
	}
	auto &catalog = db->GetCatalog();
	if (catalog.GetCatalogType() != "ena") {
		if (explicit_name) {
			throw InvalidInputException("%s: catalog '%s' is type '%s', not 'ena'", caller, catalog_name,
			                            catalog.GetCatalogType());
		}
		return nullptr;
	}
	return &catalog.Cast<duckdb::ENACatalog>();
}

bool IsWhitespaceOnly(const string &s) {
	for (char c : s) {
		if (c != ' ' && c != '\t' && c != '\r' && c != '\n') {
			return false;
		}
	}
	return true;
}

struct LifecycleBindData : public TableFunctionData {
	ENAAction action;
	string fn_name; // for error messages — owned string so plan-cache serialization is safe
	string secret_name;
	string accession;
	string until_date;             // HOLD only; ignored for CANCEL/RELEASE
	string catalog_name;           // resolved name; "ena" by default
	bool catalog_explicit = false; // user passed `catalog =>` vs defaulted
};

// One bind helper for all three lifecycle functions — they differ only in
// `action` and the user-facing function name. DuckDB's TableFunction
// constructor takes a raw function pointer for the bind, so we can't use
// a captured lambda; instead each registered function points at a thin
// thunk (BindCancel/BindRelease/BindHold below) that calls BindLifecycle
// with the right action+name.
unique_ptr<FunctionData> BindLifecycle(ClientContext &, TableFunctionBindInput &input,
                                       vector<LogicalType> &return_types, vector<string> &names, ENAAction action,
                                       const char *fn_name) {
	auto bd = make_uniq<LifecycleBindData>();
	bd->action = action;
	bd->fn_name = fn_name;

	auto get_str = [&](const char *key, string &out, bool required, bool *was_set = nullptr) {
		auto it = input.named_parameters.find(key);
		if (it == input.named_parameters.end() || it->second.IsNull()) {
			if (required) {
				throw BinderException("%s: required named parameter '%s' is missing", fn_name, key);
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
	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}
	if (action == ENAAction::HOLD) {
		get_str("until", bd->until_date, /*required=*/true);
	}

	// Reject whitespace-only accession at bind time. The envelope builder
	// also rejects, but its error surfaces at execute time after we've
	// already resolved the secret and built up Init state — failing in Bind
	// gives the cleanest user feedback. Empty (`""`) is already caught by
	// the required-parameter check above.
	if (IsWhitespaceOnly(bd->accession)) {
		throw BinderException("%s: 'accession' must not be whitespace-only", fn_name);
	}

	names = {"action", "target", "success", "era_accession", "hold_until_date", "error_messages", "duration_ms"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BOOLEAN,
	                LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::LIST(LogicalType::VARCHAR),
	                LogicalType::BIGINT};
	return std::move(bd);
}

struct LifecycleGlobalState : public GlobalTableFunctionState {
	bool emitted = false;
	LifecycleOutcome outcome;
};

// Init only allocates state. The HTTP POST is intentionally deferred to
// Execute: lifecycle actions have side effects (CANCEL is irreversible),
// so a query that filters or limits the function's output to zero rows
// MUST not perform the POST. Init runs even when the planner knows the
// scan is empty; Execute runs only when a row is actually consumed.
unique_ptr<GlobalTableFunctionState> InitLifecycleGlobal(ClientContext &, TableFunctionInitInput &) {
	return make_uniq<LifecycleGlobalState>();
}

void ExecuteLifecycle(ClientContext &context, TableFunctionInput &data, DataChunk &output) {
	auto &gs = data.global_state->Cast<LifecycleGlobalState>();
	if (gs.emitted) {
		output.SetCardinality(0);
		return;
	}
	auto &bd = data.bind_data->Cast<LifecycleBindData>();

	auto creds = duckdb::ResolveENACredentialsByName(context, bd.fn_name, bd.secret_name);
	// The shared resolver returns endpoint_url empty when the secret had no
	// override; derive the default Webin URL from the endpoint label.
	if (creds.endpoint_url.empty()) {
		creds.endpoint_url = duckdb::ResolveDefaultENAEndpointURL(creds.endpoint);
	}

	LifecycleSubmitOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.target.accession = bd.accession;
	opts.hold_until_date = bd.until_date;

	ENAClient client(*context.db);
	ENAPostFn post_fn = [&client](const string &url, const string &body, const string &user, const string &password,
	                              const string &content_type) {
		// SubmitLifecycle always uses application/xml; ENAClient::PostXML
		// matches that. Keep the dispatch tolerant in case a future call
		// path goes through here with a different content type.
		if (content_type == "application/xml") {
			return client.PostXML(url, body, user, password);
		}
		return client.PostJSONReceiveXML(url, body, user, password);
	};

	bool transport_failure = false;
	string transport_error;
	try {
		gs.outcome = SubmitLifecycle(bd.action, opts, post_fn);
	} catch (const std::exception &e) {
		transport_failure = true;
		transport_error = e.what();
		// Pre-POST exceptions (envelope-validation or transport setup) leave
		// `gs.outcome` default-initialized — fill in the fields we know the
		// log row needs so the audit trail still identifies the action+target.
		gs.outcome.action = bd.action;
		gs.outcome.target = bd.accession;
		gs.outcome.hold_until_date = bd.until_date;
	}

	const bool success = !transport_failure && gs.outcome.success;
	auto error_messages = gs.outcome.error_messages;
	if (transport_failure && error_messages.empty()) {
		error_messages.push_back(transport_error);
	}
	// Make the emitted row's `error_messages` reflect what we actually
	// computed (pre-existing logic was a no-op, but the assignment matters
	// once we add transport_error to the list above).
	gs.outcome.error_messages = error_messages;

	// submission_log integration. Default catalog "ena" is silently skipped
	// when not attached (one-shot CANCEL from a no-ATTACH session). An
	// explicitly-named catalog that's missing or wrong-typed throws so the
	// audit gap can't go unnoticed.
	if (auto *catalog = FindAttachedENACatalog(context, bd.fn_name, bd.catalog_name, bd.catalog_explicit)) {
		duckdb::ResolvedENACredentials log_creds = creds;
		duckdb::SubmissionLogPayload payload;
		payload.object_type = ""; // lifecycle ops aren't bound to a single object_type
		payload.action = ActionName(gs.outcome.action);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = gs.outcome.duration_ms;
		payload.envelope_payload = gs.outcome.envelope_payload;
		payload.raw_receipt = gs.outcome.raw_receipt;
		payload.era_accession = gs.outcome.era_accession;
		payload.error_messages = error_messages;
		payload.target = gs.outcome.target;
		duckdb::RecordSubmissionLog(*catalog, log_creds, payload);
	}

	if (transport_failure) {
		throw InvalidInputException("%s: %s", bd.fn_name, transport_error);
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
		throw InvalidInputException("%s: %s", bd.fn_name, detail);
	}

	output.SetCardinality(1);

	const auto &o = gs.outcome;
	output.data[0].SetValue(0, Value(string(ActionName(o.action))));
	output.data[1].SetValue(0, Value(o.target));
	output.data[2].SetValue(0, Value::BOOLEAN(o.success));
	output.data[3].SetValue(0, Value(o.era_accession));
	output.data[4].SetValue(0, Value(o.hold_until_date));

	// error_messages: build a list value
	vector<Value> errs;
	errs.reserve(o.error_messages.size());
	for (const auto &m : o.error_messages) {
		errs.emplace_back(m);
	}
	output.data[5].SetValue(0, Value::LIST(LogicalType::VARCHAR, errs));
	output.data[6].SetValue(0, Value::BIGINT(o.duration_ms));

	gs.emitted = true;
}

unique_ptr<FunctionData> BindCancel(ClientContext &ctx, TableFunctionBindInput &input, vector<LogicalType> &rt,
                                    vector<string> &n) {
	return BindLifecycle(ctx, input, rt, n, ENAAction::CANCEL, "ena_cancel");
}
unique_ptr<FunctionData> BindRelease(ClientContext &ctx, TableFunctionBindInput &input, vector<LogicalType> &rt,
                                     vector<string> &n) {
	return BindLifecycle(ctx, input, rt, n, ENAAction::RELEASE, "ena_release");
}
unique_ptr<FunctionData> BindHold(ClientContext &ctx, TableFunctionBindInput &input, vector<LogicalType> &rt,
                                  vector<string> &n) {
	return BindLifecycle(ctx, input, rt, n, ENAAction::HOLD, "ena_hold");
}

void AddLifecycleNamedParameters(TableFunction &tf, bool include_until) {
	tf.named_parameters["secret"] = LogicalType::VARCHAR;
	tf.named_parameters["accession"] = LogicalType::VARCHAR;
	tf.named_parameters["catalog"] = LogicalType::VARCHAR;
	if (include_until) {
		tf.named_parameters["until"] = LogicalType::VARCHAR;
	}
}

} // namespace

void RegisterENALifecycleTableFunctions(ExtensionLoader &loader) {
	{
		TableFunction tf("ena_cancel", {}, ExecuteLifecycle, BindCancel, InitLifecycleGlobal);
		AddLifecycleNamedParameters(tf, /*include_until=*/false);
		loader.RegisterFunction(tf);
	}
	{
		TableFunction tf("ena_release", {}, ExecuteLifecycle, BindRelease, InitLifecycleGlobal);
		AddLifecycleNamedParameters(tf, /*include_until=*/false);
		loader.RegisterFunction(tf);
	}
	{
		TableFunction tf("ena_hold", {}, ExecuteLifecycle, BindHold, InitLifecycleGlobal);
		AddLifecycleNamedParameters(tf, /*include_until=*/true);
		loader.RegisterFunction(tf);
	}
}

} // namespace miint
