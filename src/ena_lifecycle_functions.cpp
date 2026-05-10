// SPDX-License-Identifier: MIT
//
// SQL surface for ENA Webin V2 targeted lifecycle actions: ena_cancel,
// ena_release, ena_hold. See ena_lifecycle_functions.hpp.

#include "ena_lifecycle_functions.hpp"

#include "ena_alias_check.hpp" // AliasObjectKindFromTableName
#include "ena_client.hpp"
#include "ena_envelope_builder.hpp"
#include "ena_insert_common.hpp"
#include "ena_lifecycle_submit.hpp"
#include "ena_reports_client.hpp"
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

struct LifecycleBindData : public TableFunctionData {
	ENAAction action;
	string fn_name; // for error messages — owned string so plan-cache serialization is safe
	string secret_name;
	string accession;              // server-assigned ID — wins over refname when both set
	string refname;                // user-supplied alias; translated → accession at execute via Reports API (L5)
	string refname_kind;           // required when refname is set: 'projects' / 'samples' / 'experiments' / 'runs'
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
	get_str("accession", bd->accession, /*required=*/false);
	get_str("refname", bd->refname, /*required=*/false);
	get_str("kind", bd->refname_kind, /*required=*/false);
	get_str("catalog", bd->catalog_name, /*required=*/false, &bd->catalog_explicit);
	if (bd->catalog_name.empty()) {
		bd->catalog_name = "ena";
	}
	if (action == ENAAction::HOLD) {
		get_str("until", bd->until_date, /*required=*/true);
	}

	// Accession xor refname: exactly one identifies the target. RefDescriptor's
	// "accession wins" precedence applies if both are set, but we reject the
	// double-set case at bind to flag user confusion early (e.g. accidentally
	// passing both from a templated query). At least one is required.
	if (bd->accession.empty() && bd->refname.empty()) {
		throw BinderException("%s: required named parameter 'accession' or 'refname' is missing", fn_name);
	}
	if (!bd->accession.empty() && !bd->refname.empty()) {
		throw BinderException("%s: pass exactly one of 'accession' or 'refname' (got both)", fn_name);
	}

	// Reject whitespace-only accession / refname at bind time. The envelope
	// builder rejects whitespace-only target too, but that surfaces at execute
	// time after secret resolution; failing in Bind gives the cleanest user
	// feedback. Empty (`""`) is already caught by the xor check above.
	if (!bd->accession.empty() && duckdb::IsENAStringWhitespaceOnly(bd->accession)) {
		throw BinderException("%s: 'accession' must not be whitespace-only", fn_name);
	}
	if (!bd->refname.empty() && duckdb::IsENAStringWhitespaceOnly(bd->refname)) {
		throw BinderException("%s: 'refname' must not be whitespace-only", fn_name);
	}

	// `kind` is required iff `refname` is set — the Reports API URL is
	// kind-tagged (/projects vs /samples vs /experiments vs /runs) and the
	// alias-uniqueness scope is per-account-per-kind, so probing all four
	// would be ambiguous when the same alias is reused across kinds.
	// Accession-targeted lifecycle is kind-agnostic (the prefix encodes it).
	if (!bd->refname.empty() && bd->refname_kind.empty()) {
		throw BinderException("%s: 'kind' is required when 'refname' is set "
		                      "(one of 'projects' / 'samples' / 'experiments' / 'runs')",
		                      fn_name);
	}
	if (bd->refname.empty() && !bd->refname_kind.empty()) {
		throw BinderException("%s: 'kind' is only meaningful with 'refname' (got 'kind' without 'refname')", fn_name);
	}
	if (!bd->refname_kind.empty()) {
		// Validate the kind value lands on one of the four submittable
		// table names. AliasObjectKindFromTableName throws std::invalid_argument
		// on an unknown name; rewrap as BinderException for the user-facing
		// surface.
		try {
			(void)miint::AliasObjectKindFromTableName(bd->refname_kind);
		} catch (const std::invalid_argument &) {
			throw BinderException("%s: invalid 'kind' '%s' "
			                      "(must be one of 'projects' / 'samples' / 'experiments' / 'runs')",
			                      fn_name, bd->refname_kind);
		}
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

	ENAClient client(*context.db);

	// L5: when the user supplied refname (alias) instead of accession, hit
	// the Webin Reports API to translate before envelope build. The envelope
	// builder rejects refname on lifecycle actions (decision #12) — that's
	// defense-in-depth; by the time we hand `LifecycleSubmitOptions` over,
	// `target.accession` is set and `target.refname` is empty.
	string resolved_accession = bd.accession;
	if (!bd.refname.empty()) {
		const auto kind = miint::AliasObjectKindFromTableName(bd.refname_kind);
		// Reports base must match the endpoint the lifecycle POST will hit
		// (each server only sees its own account's records); see decision
		// notes in ena_reports_client.hpp.
		const auto reports_base = miint::ResolveReportsBaseForEndpoint(creds.endpoint);
		// Lifetime contract for the fetcher closure: `client` lives on this
		// stack frame and the fetcher is invoked synchronously by
		// LookupAccessionByAlias before this scope exits — the URLFetcher
		// type-erases through std::function but we never store / copy /
		// defer the closure beyond this call. user+password copied so the
		// closure body owns its credentials independent of `creds`.
		ENAClient *client_ptr = &client;
		const auto user = creds.user;
		const auto password = creds.password;
		miint::URLFetcher fetcher = [client_ptr, user, password](const string &url) {
			return client_ptr->AuthenticatedGet(url, user, password);
		};
		try {
			resolved_accession = miint::LookupAccessionByAlias(reports_base, kind, bd.refname, fetcher);
		} catch (const std::exception &e) {
			throw InvalidInputException("%s: Reports API lookup for refname '%s' (kind=%s) failed: %s", bd.fn_name,
			                            bd.refname, bd.refname_kind, e.what());
		}
		if (resolved_accession.empty()) {
			throw InvalidInputException("%s: refname '%s' not found in this submission account "
			                            "(kind=%s); aliases are unique per kind — verify the kind is correct, "
			                            "or supply 'accession' directly",
			                            bd.fn_name, bd.refname, bd.refname_kind);
		}
	}

	LifecycleSubmitOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.target.accession = resolved_accession;
	opts.hold_until_date = bd.until_date;
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
		// `resolved_accession` carries the post-Reports-translation value so
		// refname-targeted calls log the actual accession we attempted to act on.
		gs.outcome.action = bd.action;
		gs.outcome.target = resolved_accession;
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
	if (auto *catalog = duckdb::FindAttachedENACatalog(context, bd.fn_name, bd.catalog_name, bd.catalog_explicit)) {
		duckdb::ResolvedENACredentials log_creds = creds;
		duckdb::SubmissionLogPayload payload;
		// Audit-log object_type: when the user supplied `kind => '…'` (the L5
		// alias-targeted path), use the explicit kind so the audit row matches
		// the per-table convention used by the DELETE-by-alias path. The
		// accession-targeted path leaves it empty — the lifecycle table fns
		// don't know the kind from the accession alone (the prefix encodes it,
		// but we'd be re-implementing prefix→kind here for an audit-only field
		// when the user could have used DELETE FROM ena.X for the same effect).
		payload.object_type = bd.refname_kind;
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
	// L5: alias-based UX. `refname` is translated to accession at execute via
	// the Webin Reports API; `kind` is required alongside (one of 'projects' /
	// 'samples' / 'experiments' / 'runs'). Aliases are unique per-account-per-
	// kind, so the kind disambiguates a reused alias.
	tf.named_parameters["refname"] = LogicalType::VARCHAR;
	tf.named_parameters["kind"] = LogicalType::VARCHAR;
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
