// SPDX-License-Identifier: MIT
//
// Implementation of `DELETE FROM ena.X WHERE col=val` → SubmitLifecycle CANCEL.
// See ena_lifecycle_delete.hpp.

#include "ena_lifecycle_delete.hpp"

#include "ena_alias_check.hpp" // AliasObjectKindFromTableName
#include "ena_client.hpp"
#include "ena_envelope_builder.hpp"
#include "ena_insert_common.hpp"
#include "ena_lifecycle_submit.hpp"
#include "ena_reports_client.hpp"
#include "ena_storage.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/execution/physical_operator.hpp"
#include "duckdb/execution/physical_plan_generator.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/operator/logical_delete.hpp"
#include "duckdb/planner/operator/logical_filter.hpp"
#include "duckdb/planner/operator/logical_get.hpp"

namespace duckdb {

namespace {

// Friendly name for diagnostics ("DELETE FROM ena.<table>").
string DeleteCallerName(ENATableEntry &table) {
	return "DELETE FROM ena." + table.name;
}

// Map the WHERE column name to a target on RefDescriptor.
// - `<accession_col> = '...'` → target.accession (no Reports lookup needed).
// - `alias = '...'`           → target.refname  (translated to accession at
//                                                execute via the Webin Reports
//                                                API; the kind is implicit
//                                                from the LogicalDelete table).
//
// Kept in lock-step with the per-table accession columns declared in
// BuildENATableInfo (ena_storage.cpp); when adding new accession
// columns there, mirror them here.
//
// Returns true and fills `target` if the column maps to a valid CANCEL key.
bool ResolveDeleteTarget(ENATableKind kind, const string &column_name, const string &value,
                         miint::RefDescriptor &target) {
	if (column_name == "alias") {
		target.refname = value;
		return true;
	}
	switch (kind) {
	case ENATableKind::PROJECTS:
		if (column_name == "prjeb_accession" || column_name == "erp_accession") {
			target.accession = value;
			return true;
		}
		break;
	case ENATableKind::SAMPLES:
		if (column_name == "ers_accession" || column_name == "samea_accession") {
			target.accession = value;
			return true;
		}
		break;
	case ENATableKind::EXPERIMENTS:
		if (column_name == "erx_accession") {
			target.accession = value;
			return true;
		}
		break;
	case ENATableKind::RUNS:
		if (column_name == "err_accession") {
			target.accession = value;
			return true;
		}
		break;
	default:
		break;
	}
	return false;
}

// `payload.object_type` mirrors the singular form recorded by INSERT (matches
// `ENASamplesInsert::ObjectName()` etc.). Reusing those identifiers keeps
// downstream consumers — the `submission_log` view, audit dashboards — from
// having to special-case lifecycle vs add rows.
const char *ObjectTypeForKind(ENATableKind kind) {
	switch (kind) {
	case ENATableKind::PROJECTS:
		return "projects";
	case ENATableKind::SAMPLES:
		return "samples";
	case ENATableKind::EXPERIMENTS:
		return "experiments";
	case ENATableKind::RUNS:
		return "runs";
	default:
		return "";
	}
}

// `WHERE col = 'value'` → (column_name, value). All other shapes throw.
struct DeletePredicate {
	string column_name;
	string value;
};

[[noreturn]] void ThrowBadPredicate(const string &caller) {
	throw BinderException("%s: WHERE clause must be equality on alias or accession", caller);
}

// `BoundColumnRefExpression` (set by the binder) is rewritten to
// `BoundReferenceExpression` by the column-binding-resolver pass that runs
// before PhysicalPlanGenerator dispatches PlanDelete. The alias field —
// "samea_accession", "alias", etc. — survives both forms, so column-name
// matching works either way; just accept both classes.
bool IsBoundColumn(const Expression &expr) {
	const auto cls = expr.GetExpressionClass();
	return cls == ExpressionClass::BOUND_COLUMN_REF || cls == ExpressionClass::BOUND_REF;
}

DeletePredicate ExtractDeletePredicate(LogicalDelete &op, const string &caller) {
	// LogicalDelete's child is either:
	//   - LogicalFilter(predicate) → LogicalGet  (the standard WHERE case)
	//   - LogicalGet                              (no WHERE — reject)
	// The optimizer might in theory push the filter into LogicalGet's
	// `table_filters`, but ena_virtual_scan declares no filter pushdown
	// support, so the LogicalFilter survives optimization.
	if (op.children.size() != 1) {
		throw InternalException("%s: LogicalDelete unexpectedly has %llu children", caller,
		                        static_cast<unsigned long long>(op.children.size()));
	}
	auto &child = *op.children[0];
	if (child.type == LogicalOperatorType::LOGICAL_GET) {
		throw BinderException("%s: must filter by alias or accession (a WHERE clause is required)", caller);
	}
	if (child.type != LogicalOperatorType::LOGICAL_FILTER) {
		ThrowBadPredicate(caller);
	}
	auto &filter = child.Cast<LogicalFilter>();
	// `SplitPredicates` runs during binding for AND-combined clauses. Multiple
	// expressions therefore mean a compound WHERE — reject; users reach for
	// `ena_cancel` for those.
	if (filter.expressions.size() != 1) {
		ThrowBadPredicate(caller);
	}
	auto &expr = *filter.expressions[0];
	if (expr.GetExpressionClass() != ExpressionClass::BOUND_COMPARISON) {
		ThrowBadPredicate(caller);
	}
	auto &cmp = expr.Cast<BoundComparisonExpression>();
	if (cmp.GetExpressionType() != ExpressionType::COMPARE_EQUAL) {
		ThrowBadPredicate(caller);
	}
	// Allow `col = const` and `const = col` — neither the binder nor the
	// optimizer canonicalises operand order, so we have to.
	Expression *col_expr = nullptr;
	BoundConstantExpression *const_expr = nullptr;
	if (IsBoundColumn(*cmp.left) && cmp.right->GetExpressionClass() == ExpressionClass::BOUND_CONSTANT) {
		col_expr = cmp.left.get();
		const_expr = &cmp.right->Cast<BoundConstantExpression>();
	} else if (cmp.left->GetExpressionClass() == ExpressionClass::BOUND_CONSTANT && IsBoundColumn(*cmp.right)) {
		col_expr = cmp.right.get();
		const_expr = &cmp.left->Cast<BoundConstantExpression>();
	} else {
		ThrowBadPredicate(caller);
	}
	if (const_expr->value.IsNull()) {
		throw BinderException("%s: WHERE value must be non-NULL", caller);
	}
	// Type-guard against integer/date/timestamp constants that survived the
	// binder's coercion (e.g. `WHERE accession = 12345`). `ToString()` on
	// those would silently produce a numeric string and we'd ship it as a
	// CANCEL target — better to fail explicitly with the actual type.
	if (const_expr->value.type().id() != LogicalTypeId::VARCHAR) {
		throw BinderException("%s: WHERE value must be a string literal (got %s)", caller,
		                      const_expr->value.type().ToString());
	}
	DeletePredicate result;
	// `alias` is set to the unqualified column name by the binder
	// (`ColumnRefExpression::GetName`) and preserved through the
	// column-binding-resolver rewrite to BoundReferenceExpression.
	result.column_name = col_expr->alias;
	// If a future DuckDB optimizer ever drops the alias on the rewrite to
	// BoundReferenceExpression we'd silently see an empty column name, fall
	// through every branch in `ResolveDeleteTarget`, and emit a misleading
	// "cannot DELETE on column ''" error. Catch it loudly here instead so
	// the failure points at the resolver, not the user's predicate.
	if (result.column_name.empty()) {
		throw InternalException("%s: column reference in WHERE has no alias (expression class=%d) — please file a bug",
		                        caller, static_cast<int>(col_expr->GetExpressionClass()));
	}
	result.value = const_expr->value.ToString();
	if (result.value.empty()) {
		throw BinderException("%s: WHERE value must not be empty", caller);
	}
	return result;
}

// Source-only physical operator. Defers the side-effectful HTTP POST and
// submission_log write to first GetData call, so EXPLAIN never triggers them
// (verified via the existing LIMIT-0 test: planner skips source pipelines
// when output is unused, so the POST stays unfired).
class ENALifecycleDelete : public PhysicalOperator {
public:
	ENALifecycleDelete(PhysicalPlan &physical_plan, ENACatalog &catalog_p, string object_type_p,
	                   miint::RefDescriptor target_p)
	    : PhysicalOperator(physical_plan, PhysicalOperatorType::EXTENSION, {LogicalType::BIGINT}, 1),
	      catalog(catalog_p), object_type(std::move(object_type_p)), target(std::move(target_p)) {
	}

	ENACatalog &catalog;
	string object_type;
	miint::RefDescriptor target;

public:
	string GetName() const override {
		return "ENA_LIFECYCLE_DELETE";
	}

	bool IsSource() const override {
		return true;
	}

	class GlobalState : public GlobalSourceState {
	public:
		bool emitted = false;
	};

	unique_ptr<GlobalSourceState> GetGlobalSourceState(ClientContext &) const override {
		return make_uniq<GlobalState>();
	}

protected:
	SourceResultType GetDataInternal(ExecutionContext &context, DataChunk &chunk,
	                                 OperatorSourceInput &input) const override {
		auto &gs = input.global_state.Cast<GlobalState>();
		if (gs.emitted) {
			chunk.SetCardinality(0);
			return SourceResultType::FINISHED;
		}

		auto &client_context = context.client;
		auto creds = ResolveENACredentials(client_context, catalog);

		miint::ENAClient client(*client_context.db);

		// L5: alias-targeted DELETE — translate refname → accession via the
		// Webin Reports API before envelope build. The envelope-builder
		// rejection of `target_refname` (decision #12) stays as defense-in-
		// depth; by the time we hand the spec over, target.accession is set
		// and target.refname is empty. The kind is implicit from the table
		// being DELETEd from (vs. the lifecycle table fns where the kind is
		// an explicit named param because those fns are kind-agnostic).
		miint::RefDescriptor resolved_target = target;
		if (!resolved_target.refname.empty() && resolved_target.accession.empty()) {
			const auto kind = miint::AliasObjectKindFromTableName(object_type);
			// Reports base must match the endpoint the lifecycle POST will hit
			// (each server only sees its own account's records); see decision
			// notes in ena_reports_client.hpp.
			const auto reports_base = miint::ResolveReportsBaseForEndpoint(creds.endpoint);
			// Lifetime contract for the fetcher closure: `client` lives on
			// this stack frame and the fetcher is invoked synchronously by
			// LookupAccessionByAlias before this scope exits — the URLFetcher
			// type-erases through std::function but we never store / copy /
			// defer the closure beyond this call. Mirrors the comment in
			// ena_lifecycle_functions.cpp::ExecuteLifecycle.
			miint::ENAClient *client_ptr = &client;
			const auto user = creds.user;
			const auto password = creds.password;
			miint::URLFetcher fetcher = [client_ptr, user, password](const string &url) {
				return client_ptr->AuthenticatedGet(url, user, password);
			};
			std::string accession;
			try {
				accession = miint::LookupAccessionByAlias(reports_base, kind, resolved_target.refname, fetcher);
			} catch (const std::exception &e) {
				throw InvalidInputException("DELETE FROM ena.%s: Reports API lookup for alias '%s' failed: %s",
				                            object_type, resolved_target.refname, e.what());
			}
			if (accession.empty()) {
				throw InvalidInputException("DELETE FROM ena.%s: alias '%s' not found in this submission account",
				                            object_type, resolved_target.refname);
			}
			resolved_target.accession = accession;
			resolved_target.refname.clear();
		}

		miint::LifecycleSubmitOptions opts;
		opts.endpoint_url = creds.endpoint_url + "/submit";
		opts.user = creds.user;
		opts.password = creds.password;
		opts.target = resolved_target;
		// CANCEL takes no date — leave hold_until_date empty.
		miint::ENAPostFn post_fn = [&client](const std::string &url, const std::string &body, const std::string &user,
		                                     const std::string &password, const std::string &content_type) {
			if (content_type == "application/xml") {
				return client.PostXML(url, body, user, password);
			}
			return client.PostJSONReceiveXML(url, body, user, password);
		};

		miint::LifecycleOutcome outcome;
		bool transport_failure = false;
		string transport_error;
		try {
			outcome = miint::SubmitLifecycle(miint::ENAAction::CANCEL, opts, post_fn);
		} catch (const std::exception &e) {
			transport_failure = true;
			transport_error = e.what();
			outcome.action = miint::ENAAction::CANCEL;
			// `resolved_target.accession` is the post-Reports-translation
			// value on the alias path and the user-supplied accession on the
			// direct-accession path. Either way, by this point .refname is
			// empty (cleared after Reports translation; never set on the
			// direct path) and .accession is non-empty — the Reports lookup
			// throws on miss before SubmitLifecycle is reached, so we only
			// land here with a known accession in hand.
			outcome.target = resolved_target.accession;
		}

		const bool success = !transport_failure && outcome.success;
		auto error_messages = outcome.error_messages;
		if (transport_failure && error_messages.empty()) {
			error_messages.push_back(transport_error);
		}

		// Always log — even on failure — so the audit trail captures every
		// attempted CANCEL. The submit + log + error-aggregation block below
		// duplicates the equivalent block in
		// `ena_lifecycle_functions.cpp::ExecuteLifecycle`. We don't share a
		// helper because the two call sites disagree on catalog plumbing —
		// the table-fn version resolves the catalog by user-supplied name
		// (and silently skips logging when not attached); here the catalog
		// is baked into the operator at plan time. A shared helper would
		// require parameterising over both shapes and isn't worth it for
		// a ~30-line block. See L1b reviewer-fixes in
		// localdocs/ena-lifecycle-plan.md.
		SubmissionLogPayload payload;
		payload.object_type = object_type;
		payload.action = miint::ActionName(outcome.action);
		payload.n_objects = 1;
		payload.success = success;
		payload.duration_ms = outcome.duration_ms;
		payload.envelope_payload = outcome.envelope_payload;
		payload.raw_receipt = outcome.raw_receipt;
		payload.era_accession = outcome.era_accession;
		payload.error_messages = error_messages;
		payload.target = outcome.target;
		RecordSubmissionLog(catalog, creds, payload);

		gs.emitted = true;
		if (transport_failure) {
			throw InvalidInputException("DELETE FROM ena.%s: %s", object_type, transport_error);
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
			throw InvalidInputException("DELETE FROM ena.%s: %s", object_type, detail);
		}

		// LogicalDelete's "Count" return type is BIGINT. We treat one
		// successful CANCEL as a single row deleted — matches the user's
		// mental model of `DELETE … WHERE accession='X'` even though our
		// catalog is virtual.
		chunk.SetCardinality(1);
		chunk.data[0].SetValue(0, Value::BIGINT(1));
		return SourceResultType::FINISHED;
	}
};

} // namespace

} // namespace duckdb

namespace miint {

duckdb::PhysicalOperator &PlanENALifecycleDelete(duckdb::ClientContext &, duckdb::PhysicalPlanGenerator &planner,
                                                 duckdb::ENACatalog &catalog, duckdb::LogicalDelete &op) {
	auto &table_entry = op.table.Cast<duckdb::ENATableEntry>();
	const auto kind = table_entry.GetKind();
	const std::string caller = duckdb::DeleteCallerName(table_entry);

	if (op.return_chunk) {
		throw duckdb::BinderException("%s: RETURNING is not supported", caller);
	}
	const char *object_type = duckdb::ObjectTypeForKind(kind);
	if (object_type[0] == '\0') {
		throw duckdb::BinderException(
		    "%s: DELETE is not implemented for this table (CANCEL only applies to projects/samples/experiments/runs)",
		    caller);
	}
	// Surface the "ATTACH didn't supply a SECRET" case at plan time as a
	// BinderException — the right exception type for the planner phase.
	// The execute-time `ResolveENACredentials` call will repeat the check
	// (in case the secret is dropped between plan and execute), so this is
	// purely about giving the user a clean diagnostic on the common case.
	if (catalog.GetSecretName().empty()) {
		throw duckdb::BinderException("%s: ENA catalog has no SECRET — re-attach with (TYPE ENA, SECRET 'name')",
		                              caller);
	}

	auto pred = duckdb::ExtractDeletePredicate(op, caller);

	RefDescriptor target;
	if (!duckdb::ResolveDeleteTarget(kind, pred.column_name, pred.value, target)) {
		throw duckdb::BinderException("%s: cannot DELETE on column '%s' (use 'alias' or a per-table accession column)",
		                              caller, pred.column_name);
	}

	return planner.Make<duckdb::ENALifecycleDelete>(catalog, std::string(object_type), std::move(target));
}

} // namespace miint
