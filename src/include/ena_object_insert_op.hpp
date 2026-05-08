// SPDX-License-Identifier: MIT
//
// CRTP base for ENA virtual-catalog INSERT operators (ENAProjectsInsert,
// ENASamplesInsert, future ENAExperimentsInsert / ENARunsInsert).
//
// The Sink (chunk buffering) and Source (RETURNING / count) machinery is
// identical across object types. The Finalize body is also nearly identical:
// resolve credentials, build specs from the buffered ColumnDataCollection,
// fill the per-table OptsT, POST via the shared ENAClient, log to
// ena.submission_log, and append RETURNING rows.
//
// Derived types provide:
//   - static const char *ObjectName()         // "projects" / "samples"
//   - static const char *ThrowPrefix()        // "INSERT INTO ena.projects"
//   - static string OperatorName()            // "ENA_PROJECTS_INSERT"
//   - static BuiltSpecs<SpecT> BuildFromBuffer(ColumnDataCollection &,
//                                              const physical_index_vector_t<idx_t> &)
//   - static OutcomeT Submit(const vector<SpecT> &, const OptsT &, const miint::ENAPostFn &)
//   - static void AppendReturningRows(ColumnDataCollection &, const vector<LogicalType> &,
//                                     const vector<SpecT> &, const OutcomeT &)

#pragma once

#include "ena_alias_check.hpp"
#include "ena_client.hpp"
#include "ena_envelope_builder.hpp" // miint::ENAAction, miint::ActionName
#include "ena_insert_common.hpp"
#include "ena_post_fn.hpp"
#include "ena_storage.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/index_vector.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/types/column/column_data_scan_states.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/execution/physical_operator.hpp"
#include "duckdb/execution/progress_data.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/event.hpp"
#include "duckdb/planner/operator/logical_insert.hpp"

#include <exception>
#include <string>
#include <vector>

namespace duckdb {

class ENATableEntry;

// Build-from-buffer return shape. `hold_until_date` is empty for object types
// that don't expose a per-row HOLD column (e.g. samples).
template <class SpecT>
struct ENABuiltSpecs {
	std::vector<SpecT> specs;
	std::string hold_until_date;
};

// Common Sink state. Buffers input chunks; Finalize drains them, builds the
// envelope, POSTs it, parses the receipt, populates a return chunk (when
// RETURNING is present), and writes a submission_log row.
class ENAObjectInsertGlobalSinkState : public GlobalSinkState {
public:
	ENAObjectInsertGlobalSinkState(ClientContext &context, vector<LogicalType> input_types,
	                               vector<LogicalType> return_types_p)
	    : buffered(context, std::move(input_types)), return_types(std::move(return_types_p)),
	      return_collection(context, return_types) {
	}

	mutex lock;
	ColumnDataCollection buffered;
	vector<LogicalType> return_types;
	ColumnDataCollection return_collection;
	idx_t insert_count = 0;
};

class ENAObjectInsertGlobalSourceState : public GlobalSourceState {
public:
	ColumnDataScanState scan_state;
	bool finished = false;
};

// Fill the user/password/endpoint/hold/action fields shared by every per-table
// OptsT (`ENAProjectInsertOptions`, `ENASampleInsertOptions`, ...). Centralised
// so a future field addition lands in one place. `action` is required: the
// CRTP base picks ADD or VALIDATE based on the session pragma, and a future
// MODIFY phase will pass MODIFY here too.
template <class OptsT>
void FillCommonENAInsertOptions(OptsT &opts, const ResolvedENACredentials &creds, const std::string &hold,
                                miint::ENAAction action) {
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.hold_until_date = hold;
	opts.action = action;
}

// Read the session-scoped `miint_ena_validate_only` flag. Registered as a
// BOOLEAN extension option in `LoadInternal` (see src/miint_extension.cpp).
// True → flip the next ena.* INSERT to VALIDATE (server-side dry-run, no
// accessions returned). Default false. Reads via TryGetCurrentSetting so a
// missing-option fallback would surface as a programmer error rather than
// silently disabling the feature.
bool IsENAValidateOnlyEnabled(ClientContext &context);

template <class Derived, class SpecT, class OptsT, class OutcomeT>
class ENAObjectInsertOperator : public PhysicalOperator {
public:
	static constexpr const PhysicalOperatorType TYPE = PhysicalOperatorType::EXTENSION;

public:
	ENAObjectInsertOperator(PhysicalPlan &physical_plan, LogicalOperator &op, ENATableEntry &table_p,
	                        physical_index_vector_t<idx_t> column_index_map_p, bool return_chunk_p)
	    : PhysicalOperator(physical_plan, PhysicalOperatorType::EXTENSION, op.types, 1), table(table_p),
	      column_index_map(std::move(column_index_map_p)), return_chunk(return_chunk_p) {
	}

	ENATableEntry &table;
	physical_index_vector_t<idx_t> column_index_map;
	bool return_chunk;

	bool IsSource() const override {
		return true;
	}
	bool IsSink() const override {
		return true;
	}
	bool ParallelSink() const override {
		return false;
	}
	string GetName() const override {
		return Derived::OperatorName();
	}

	// Mark the sink progress invalid. Our Finalize blocks for several
	// seconds across multiple synchronous HTTPS round-trips (checklist
	// fetch on samples, alias collision check, envelope POST, receipt
	// parse). DuckDB's progress-bar UKF
	// (duckdb/src/common/progress_bar/unscented_kalman_filter.cpp) becomes
	// unstable under sparse measurements like that — its innovation
	// covariance `S[0][0]` can converge to 0, and `S_inv = 1.0 / S[0][0]`
	// then produces NaN/Inf that surfaces as SIGFPE in the duckdb shell
	// (observed during live wwwdev smoke). Reporting invalid progress
	// short-circuits ProgressBar::Update before it feeds the UKF.
	ProgressData GetSinkProgress(ClientContext & /*context*/, GlobalSinkState & /*gstate*/,
	                             const ProgressData source_progress) const override {
		ProgressData pd = source_progress;
		pd.SetInvalid();
		return pd;
	}

	unique_ptr<GlobalSinkState> GetGlobalSinkState(ClientContext &context) const override {
		D_ASSERT(!children.empty());
		auto input_types = children[0].get().GetTypes();
		return make_uniq<ENAObjectInsertGlobalSinkState>(context, std::move(input_types), types);
	}

	SinkResultType Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const override {
		auto &state = input.global_state.Cast<ENAObjectInsertGlobalSinkState>();
		lock_guard<mutex> guard(state.lock);
		state.buffered.Append(chunk);
		state.insert_count += chunk.size();
		return SinkResultType::NEED_MORE_INPUT;
	}

	SinkFinalizeType Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
	                          OperatorSinkFinalizeInput &input) const override {
		auto &state = input.global_state.Cast<ENAObjectInsertGlobalSinkState>();
		if (state.buffered.Count() == 0) {
			return SinkFinalizeType::READY;
		}

		auto &catalog = table.catalog.Cast<ENACatalog>();
		auto creds = ResolveENACredentials(context, catalog);

		auto built = Derived::BuildFromBuffer(state.buffered, column_index_map);

		const bool validate_only = IsENAValidateOnlyEnabled(context);
		const auto action = validate_only ? miint::ENAAction::VALIDATE : miint::ENAAction::ADD;

		OptsT opts;
		FillCommonENAInsertOptions(opts, creds, built.hold_until_date, action);

		miint::ENAClient client(*context.db);

		// Per-table client-side validation. ENASamplesInsert overrides this
		// to validate user-supplied attributes against the chosen ENA
		// checklist (mandatory fields, units, CV). Default no-op for other
		// object types. Runs before the alias check so the user gets the
		// most informative error first (validation issues say "fix your
		// attributes"; alias collision says "this name is taken").
		Derived::ValidateBuiltSpecs(built.specs, client);

		// Pre-INSERT: ask the ENA portal API whether any of these aliases
		// already exist in the submission account. Aliases are unique per
		// (account, object_type) on the server side; reuse is a hard error
		// in the receipt. Catching it here lets us surface the offending
		// aliases by name and saves the envelope POST + receipt parse.
		//
		// Throws InvalidInputException early on collision — intentionally
		// before RecordSubmissionLog below, since this is a client-side
		// validation failure (user error) and not a recorded submission
		// attempt. Pinned by test/sql/ena_alias_collision_mock.test which
		// asserts no submission_log row appears for a blocked alias.
		//
		// Skipped under VALIDATE: dry-run validation against an alias the
		// user has already registered (e.g. before a future MODIFY) is a
		// legitimate workflow — the collision check would otherwise block
		// the very thing the user is trying to dry-run. The server-side
		// receipt still surfaces any genuine alias issue if the validation
		// itself rejects it.
		if (!validate_only) {
			RunAliasCollisionCheck(built.specs, client, opts);
		}
		// Dispatch by Content-Type. Post L4b-fix (2026-05-07), V2 JSON is
		// used only for projects; samples + experiment + run + analysis go
		// through XML. Samples were on JSON until live wwwdev surfaced an
		// XSD-ordering bug in V2's JSON-to-XML dispatcher (DESCRIPTION
		// emitted before SAMPLE_NAME); experiment/run/analysis already
		// required XML because the JSON dispatcher NPEs for SRA-side
		// objects. Both paths request an XML receipt — our parser only
		// consumes the canonical XSD shape.
		auto post_fn = [&client](const std::string &url, const std::string &body, const std::string &user,
		                         const std::string &password, const std::string &content_type) {
			if (content_type == "application/xml") {
				return client.PostXML(url, body, user, password);
			}
			return client.PostJSONReceiveXML(url, body, user, password);
		};

		// Logical failures (server success=false, parse errors, missing aliases)
		// surface as `outcome.success=false` with `error_messages` populated and
		// `raw_receipt` preserved — that lets the submission_log row keep the
		// server response even when parsing failed. Transport failures (network,
		// 5xx after exhausted retries, malformed POST) come back as exceptions.
		OutcomeT outcome;
		bool transport_failure = false;
		std::string transport_error;
		try {
			outcome = Derived::Submit(built.specs, opts, post_fn);
		} catch (const std::exception &e) {
			transport_failure = true;
			transport_error = e.what();
		}

		const bool success = !transport_failure && outcome.success;
		auto error_messages = outcome.error_messages;
		if (transport_failure && error_messages.empty()) {
			error_messages.push_back(transport_error);
		}

		SubmissionLogPayload log_payload;
		log_payload.object_type = Derived::ObjectName();
		log_payload.action = miint::ActionName(opts.action);
		log_payload.n_objects = static_cast<int32_t>(built.specs.size());
		log_payload.success = success;
		log_payload.duration_ms = outcome.duration_ms;
		log_payload.envelope_payload = outcome.envelope_payload;
		log_payload.raw_receipt = outcome.raw_receipt;
		log_payload.era_accession = outcome.era_accession;
		log_payload.error_messages = error_messages;
		// Capture per-object alias / primary accession into the audit log so
		// users can recover an accession from an alias within the session
		// (cross-submission lifecycle ops on Webin V2 require accession
		// targeting — see docs/ena.md). On a parse failure or partial
		// receipt `outcome.rows` is empty, so the lists stay empty too.
		log_payload.object_aliases.reserve(outcome.rows.size());
		log_payload.object_accessions.reserve(outcome.rows.size());
		for (const auto &row : outcome.rows) {
			log_payload.object_aliases.push_back(row.alias);
			log_payload.object_accessions.push_back(Derived::PrimaryAccession(row));
		}
		RecordSubmissionLog(catalog, creds, log_payload);

		if (transport_failure) {
			throw InvalidInputException("%s failed: %s", Derived::ThrowPrefix(), transport_error);
		}
		if (!outcome.success) {
			std::string detail;
			if (error_messages.empty()) {
				detail = "no error detail";
			} else {
				for (const auto &m : error_messages) {
					if (!detail.empty()) {
						detail += "; ";
					}
					detail += m;
				}
			}
			throw InvalidInputException("%s failed: %s", Derived::ThrowPrefix(), detail);
		}

		// VALIDATE leaves outcome.rows empty (no per-object accessions in the
		// receipt). AppendReturningRows asserts rows.size() == specs.size() in
		// debug builds, so guard the call rather than weaken that invariant —
		// it's load-bearing for the ADD path.
		if (return_chunk && !outcome.rows.empty()) {
			Derived::AppendReturningRows(state.return_collection, state.return_types, built.specs, outcome);
		}
		return SinkFinalizeType::READY;
	}

	unique_ptr<GlobalSourceState> GetGlobalSourceState(ClientContext &context) const override {
		auto state = make_uniq<ENAObjectInsertGlobalSourceState>();
		if (return_chunk && sink_state) {
			auto &gstate = sink_state->Cast<ENAObjectInsertGlobalSinkState>();
			gstate.return_collection.InitializeScan(state->scan_state);
		}
		return std::move(state);
	}

	SourceResultType GetDataInternal(ExecutionContext &context, DataChunk &chunk,
	                                 OperatorSourceInput &input) const override {
		auto &source_state = input.global_state.Cast<ENAObjectInsertGlobalSourceState>();
		D_ASSERT(sink_state);
		auto &gstate = sink_state->Cast<ENAObjectInsertGlobalSinkState>();

		if (!return_chunk) {
			if (source_state.finished) {
				chunk.SetCardinality(0);
				return SourceResultType::FINISHED;
			}
			chunk.SetCardinality(1);
			chunk.SetValue(0, 0, Value::BIGINT(NumericCast<int64_t>(gstate.insert_count)));
			source_state.finished = true;
			return SourceResultType::FINISHED;
		}

		gstate.return_collection.Scan(source_state.scan_state, chunk);
		return chunk.size() == 0 ? SourceResultType::FINISHED : SourceResultType::HAVE_MORE_OUTPUT;
	}

public:
	// Default no-op validation hook. Derived classes may override (static
	// dispatch via Derived::ValidateBuiltSpecs in Finalize). Called between
	// BuildFromBuffer and the alias collision check.
	static void ValidateBuiltSpecs(const std::vector<SpecT> & /*specs*/, miint::ENAClient & /*client*/) {
	}

private:
	// Pre-INSERT alias collision check. Static so unit tests at the per-table
	// layer can also reach it without instantiating a PhysicalOperator.
	//
	// Anonymous portal-API GET — the ENA portal only indexes public records,
	// and adding Authorization: Basic to a portal query returns HTTP 500
	// (the endpoint has no authenticated mode). Same approach as
	// ena-upload-cli's `check_remote.py`. Caveat: collisions against the
	// user's OWN HOLD/private records that haven't yet been published won't
	// be caught here; those surface at envelope-POST time as a server error.
	// The Reports API (`/ena/submit/report/...`) would close that gap; not
	// wired yet.
	static void RunAliasCollisionCheck(const std::vector<SpecT> &specs, miint::ENAClient &client,
	                                   const OptsT & /*opts*/) {
		if (specs.empty()) {
			return;
		}
		std::vector<std::string> aliases;
		aliases.reserve(specs.size());
		for (const auto &s : specs) {
			aliases.push_back(s.alias);
		}
		const auto kind = miint::AliasObjectKindFromTableName(Derived::ObjectName());
		const auto portal_base = miint::ResolvePortalBaseFromEnv();
		miint::URLFetcher fetcher = [&client](const std::string &url) {
			return client.FetchURL(url);
		};
		const auto hits = miint::CheckAliasCollisions(portal_base, kind, aliases, fetcher);
		if (hits.empty()) {
			return;
		}
		std::string detail;
		for (size_t i = 0; i < hits.size(); i++) {
			if (i > 0) {
				detail += ", ";
			}
			detail += "'";
			detail += hits[i];
			detail += "'";
		}
		const char *noun = hits.size() == 1 ? "alias" : "aliases";
		const char *verb = hits.size() == 1 ? "exists" : "exist";
		throw InvalidInputException("%s: %s already %s in submission account: %s", Derived::ThrowPrefix(), noun, verb,
		                            detail);
	}
};

} // namespace duckdb
