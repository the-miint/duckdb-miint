// SPDX-License-Identifier: MIT
//
// Implementation of ENAProjectsInsert. See ena_projects_insert_op.hpp.

#include "ena_projects_insert_op.hpp"

#include "ena_client.hpp"
#include "ena_projects_insert.hpp"
#include "ena_storage.hpp"

#include "duckdb/catalog/catalog_transaction.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/types/column/column_data_scan_states.hpp"
#include "duckdb/common/types/date.hpp"
#include "duckdb/common/types/timestamp.hpp"
#include "duckdb/common/types/uuid.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/main/attached_database.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/main/secret/secret.hpp"
#include "duckdb/main/secret/secret_manager.hpp"
#include "duckdb/parallel/event.hpp"
#include "duckdb/planner/operator/logical_insert.hpp"

#include <chrono>

namespace duckdb {

namespace {

// Per-statement Sink state. Buffers input chunks; Finalize drains them, builds
// the envelope, POSTs it, parses the receipt, and assembles a return chunk
// (when RETURNING is present) plus a submission_log row.
class ENAProjectsInsertGlobalSinkState : public GlobalSinkState {
public:
	ENAProjectsInsertGlobalSinkState(ClientContext &context, vector<LogicalType> input_types,
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

class ENAProjectsInsertGlobalSourceState : public GlobalSourceState {
public:
	ColumnDataScanState scan_state;
	bool finished = false;
};

idx_t ResolveInputColumn(const physical_index_vector_t<idx_t> &column_index_map, idx_t input_chunk_columns,
                         idx_t table_column_index) {
	if (column_index_map.empty()) {
		// Input chunk has the full table layout in declaration order.
		return table_column_index < input_chunk_columns ? table_column_index : DConstants::INVALID_INDEX;
	}
	// `physical_index_vector_t` is sparse; out-of-range reads silently return
	// `INVALID_INDEX` (the "not provided" sentinel). Treat an undersized map
	// the same way to keep behaviour deterministic if a future refactor of
	// LogicalInsert ever truncates the map.
	if (table_column_index >= column_index_map.size()) {
		return DConstants::INVALID_INDEX;
	}
	const auto mapped = column_index_map[PhysicalIndex(table_column_index)];
	if (mapped == DConstants::INVALID_INDEX) {
		return DConstants::INVALID_INDEX;
	}
	return mapped;
}

string ValueToVarchar(const Value &v) {
	if (v.IsNull()) {
		return string();
	}
	return v.ToString();
}

string ValueToDateString(const Value &v) {
	if (v.IsNull()) {
		return string();
	}
	if (v.type().id() != LogicalTypeId::DATE) {
		throw InvalidInputException("ENA projects: hold_until_date must be a DATE");
	}
	return Date::ToString(v.GetValue<date_t>());
}

// Fetch the user/password/endpoint_url to use for this submission. The
// ENACatalog stores the secret name supplied at ATTACH time; we look it up
// just-in-time here so secret edits between ATTACH and INSERT take effect.
struct ResolvedENACredentials {
	string user;
	string password;
	string endpoint_url;
	string endpoint;
	string secret_name;
};

unique_ptr<SecretEntry> LookupSecret(ClientContext &context, const string &secret_name) {
	auto &secret_manager = SecretManager::Get(context);
	auto transaction = CatalogTransaction::GetSystemCatalogTransaction(context);
	// Empty storage string asks the manager to search every registered backend
	// (memory, local_file, plus any extension-provided backends). The narrow
	// memory+local_file path used by duckdb-postgres misses secrets stored
	// elsewhere — confusing for users whose `duckdb_secrets()` shows the secret
	// but `INSERT` complains it's missing.
	return secret_manager.GetSecretByName(transaction, secret_name);
}

ResolvedENACredentials ResolveCredentials(ClientContext &context, ENACatalog &catalog) {
	ResolvedENACredentials creds;
	creds.endpoint = catalog.GetEndpoint();
	creds.endpoint_url = catalog.GetEndpointURL();
	creds.secret_name = catalog.GetSecretName();
	if (creds.secret_name.empty()) {
		throw BinderException("ENA INSERT requires a SECRET — re-attach with (TYPE ENA, SECRET 'name')");
	}
	auto entry = LookupSecret(context, creds.secret_name);
	if (!entry) {
		throw BinderException("ENA secret '%s' not found", creds.secret_name);
	}
	auto kv = dynamic_cast<const KeyValueSecret *>(entry->secret.get());
	if (!kv) {
		throw InvalidInputException("ENA secret '%s' is not a KeyValueSecret — was it created with TYPE ENA?",
		                            creds.secret_name);
	}
	auto user_val = kv->TryGetValue("user");
	auto password_val = kv->TryGetValue("password");
	if (user_val.IsNull() || password_val.IsNull()) {
		throw BinderException("ENA secret '%s' is missing user or password", creds.secret_name);
	}
	creds.user = user_val.ToString();
	creds.password = password_val.ToString();
	auto endpoint_url_val = kv->TryGetValue("endpoint_url");
	if (!endpoint_url_val.IsNull()) {
		creds.endpoint_url = endpoint_url_val.ToString();
	}
	return creds;
}

// Build ProjectSpecs from the buffered ColumnDataCollection and surface a
// uniform hold_until_date if present in any row. The ENA Webin V2 envelope
// applies HOLD at the submission level, not per-object, so for now we require
// agreement when more than one row has it set.
struct BuiltProjects {
	vector<miint::ProjectSpec> projects;
	string hold_until_date;
};

constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_TITLE = 1;
constexpr idx_t COL_DESCRIPTION = 2;
constexpr idx_t COL_PROJECT_TYPE = 3;
constexpr idx_t COL_UMBRELLA_CHILDREN = 4;
constexpr idx_t COL_ATTRIBUTES = 5;
constexpr idx_t COL_PRJEB_ACCESSION = 6;
constexpr idx_t COL_ERP_ACCESSION = 7;
constexpr idx_t COL_HOLD_UNTIL_DATE = 8;
constexpr idx_t TABLE_COLUMN_COUNT = 9;

BuiltProjects BuildProjectsFromBuffer(ColumnDataCollection &buffer,
                                      const physical_index_vector_t<idx_t> &column_index_map) {
	BuiltProjects out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.projects requires the 'alias' column");
	}

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::ProjectSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.projects: 'alias' must be non-empty");
			}
			const auto title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			const auto description_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESCRIPTION);
			if (description_idx != DConstants::INVALID_INDEX) {
				spec.description = ValueToVarchar(chunk.data[description_idx].GetValue(row));
			}
			const auto project_type_idx = ResolveInputColumn(column_index_map, input_columns, COL_PROJECT_TYPE);
			if (project_type_idx != DConstants::INVALID_INDEX) {
				spec.project_type = ValueToVarchar(chunk.data[project_type_idx].GetValue(row));
			}
			const auto hold_idx = ResolveInputColumn(column_index_map, input_columns, COL_HOLD_UNTIL_DATE);
			if (hold_idx != DConstants::INVALID_INDEX) {
				const auto row_hold = ValueToDateString(chunk.data[hold_idx].GetValue(row));
				if (!row_hold.empty()) {
					if (out.hold_until_date.empty()) {
						out.hold_until_date = row_hold;
					} else if (out.hold_until_date != row_hold) {
						throw InvalidInputException(
						    "INSERT INTO ena.projects: per-row hold_until_date values must agree across the "
						    "statement (got '%s' and '%s')",
						    out.hold_until_date, row_hold);
					}
				}
			}
			// umbrella_children / attributes / prjeb_accession / erp_accession: not
			// yet wired in Phase 4. Reject explicit non-null umbrella_children so
			// users don't think it's silently honoured.
			const auto umbrella_idx = ResolveInputColumn(column_index_map, input_columns, COL_UMBRELLA_CHILDREN);
			if (umbrella_idx != DConstants::INVALID_INDEX) {
				const auto v = chunk.data[umbrella_idx].GetValue(row);
				if (!v.IsNull() && !ListValue::GetChildren(v).empty()) {
					throw NotImplementedException(
					    "INSERT INTO ena.projects: 'umbrella_children' is not supported in this build");
				}
			}
			const auto attrs_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTES);
			if (attrs_idx != DConstants::INVALID_INDEX) {
				const auto v = chunk.data[attrs_idx].GetValue(row);
				if (!v.IsNull() && !ListValue::GetChildren(v).empty()) {
					throw NotImplementedException(
					    "INSERT INTO ena.projects: 'attributes' is not supported in this build");
				}
			}
			out.projects.push_back(std::move(spec));
		}
	}
	return out;
}

void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
                         const vector<miint::ProjectSpec> &projects, const miint::ENASubmissionOutcome &outcome) {
	// `LogicalInsert::ResolveTypes` sets `op.types = table.GetTypes()` when
	// return_chunk is true, so we always receive the full 9-column projects
	// schema. The DuckDB executor inserts a PhysicalProjection above us to
	// trim down to the user-requested RETURNING columns. Hard-asserting keeps
	// us honest if that contract ever changes upstream.
	D_ASSERT(return_types.size() == TABLE_COLUMN_COUNT);
	// `outcome.rows` is built in alias-order from `projects`, so positional
	// indexing is sufficient and avoids an O(n) per-row lookup.
	D_ASSERT(outcome.rows.size() == projects.size());

	DataChunk chunk;
	chunk.Initialize(Allocator::DefaultAllocator(), return_types);
	for (idx_t i = 0; i < outcome.rows.size(); i++) {
		const auto &row = outcome.rows[i];
		const auto &spec = projects[i];
		const auto idx = chunk.size();
		chunk.data[COL_ALIAS].SetValue(idx, Value(row.alias));
		chunk.data[COL_TITLE].SetValue(idx, Value(spec.title));
		chunk.data[COL_DESCRIPTION].SetValue(idx, Value(spec.description));
		chunk.data[COL_PROJECT_TYPE].SetValue(idx, Value(spec.project_type));
		chunk.data[COL_UMBRELLA_CHILDREN].SetValue(idx, Value(LogicalType::LIST(LogicalType::VARCHAR)));
		chunk.data[COL_ATTRIBUTES].SetValue(idx, Value(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)));
		chunk.data[COL_PRJEB_ACCESSION].SetValue(idx, Value(row.prjeb_accession));
		chunk.data[COL_ERP_ACCESSION].SetValue(idx, row.erp_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                      : Value(row.erp_accession));
		if (row.hold_until_date.empty()) {
			chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value(LogicalType::DATE));
		} else {
			date_t d;
			idx_t pos = 0;
			bool special = false;
			if (Date::TryConvertDate(row.hold_until_date.c_str(), row.hold_until_date.size(), pos, d, special) ==
			    DateCastResult::SUCCESS) {
				chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value::DATE(d));
			} else {
				chunk.data[COL_HOLD_UNTIL_DATE].SetValue(idx, Value(LogicalType::DATE));
			}
		}
		chunk.SetCardinality(idx + 1);
		if (chunk.size() == STANDARD_VECTOR_SIZE) {
			return_collection.Append(chunk);
			chunk.Reset();
		}
	}
	if (chunk.size() > 0) {
		return_collection.Append(chunk);
	}
}

string GenerateUUIDString() {
	auto value = UUID::GenerateRandomUUID();
	return UUID::ToString(value);
}

void RecordSubmissionLog(ENACatalog &catalog, const ResolvedENACredentials &creds, idx_t n_objects,
                         const miint::ENASubmissionOutcome &outcome, bool success,
                         const std::vector<std::string> &error_messages) {
	ENASubmissionLogRow row;
	row.submission_id = GenerateUUIDString();
	row.submitted_at = timestamp_tz_t(Timestamp::GetCurrentTimestamp());
	row.endpoint = creds.endpoint;
	row.secret_name = creds.secret_name;
	row.action = "ADD";
	row.object_type = "projects";
	row.n_objects = static_cast<int32_t>(n_objects);
	row.success = success;
	row.era_accession = outcome.era_accession;
	row.request_payload = outcome.envelope_payload;
	row.receipt = outcome.raw_receipt;
	row.error_messages.clear();
	row.error_messages.reserve(error_messages.size());
	for (const auto &m : error_messages) {
		row.error_messages.push_back(m);
	}
	row.duration_ms = static_cast<int32_t>(outcome.duration_ms);
	catalog.GetSubmissionLog().Append(row);
}

} // namespace

ENAProjectsInsert::ENAProjectsInsert(PhysicalPlan &physical_plan, LogicalOperator &op, ENATableEntry &table_p,
                                     physical_index_vector_t<idx_t> column_index_map_p, bool return_chunk_p)
    : PhysicalOperator(physical_plan, PhysicalOperatorType::EXTENSION, op.types, 1), table(table_p),
      column_index_map(std::move(column_index_map_p)), return_chunk(return_chunk_p) {
}

unique_ptr<GlobalSinkState> ENAProjectsInsert::GetGlobalSinkState(ClientContext &context) const {
	// The input chunk's column count and types must come from the child plan,
	// not the operator's `types` (which describes the RETURNING projection).
	D_ASSERT(!children.empty());
	auto input_types = children[0].get().GetTypes();
	return make_uniq<ENAProjectsInsertGlobalSinkState>(context, std::move(input_types), types);
}

SinkResultType ENAProjectsInsert::Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const {
	auto &state = input.global_state.Cast<ENAProjectsInsertGlobalSinkState>();
	lock_guard<mutex> guard(state.lock);
	state.buffered.Append(chunk);
	state.insert_count += chunk.size();
	return SinkResultType::NEED_MORE_INPUT;
}

SinkFinalizeType ENAProjectsInsert::Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
                                             OperatorSinkFinalizeInput &input) const {
	auto &state = input.global_state.Cast<ENAProjectsInsertGlobalSinkState>();
	if (state.buffered.Count() == 0) {
		return SinkFinalizeType::READY;
	}

	auto &catalog = table.catalog.Cast<ENACatalog>();
	auto creds = ResolveCredentials(context, catalog);

	auto built = BuildProjectsFromBuffer(state.buffered, column_index_map);

	miint::ENAProjectInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.hold_until_date = built.hold_until_date;

	miint::ENAClient client(*context.db);
	auto post_fn = [&client](const std::string &url, const std::string &body, const std::string &user,
	                         const std::string &password, const std::string &content_type) {
		(void)content_type; // PostJSON always sends application/json
		return client.PostJSON(url, body, user, password);
	};

	// `SubmitProjectInsertOutcome` only throws on transport failures (network,
	// 5xx after exhausted retries, malformed POST). Logical failures (server
	// returned success=false, parse errors, missing receipt entries) come back
	// as `outcome.success == false` with `error_messages` populated and
	// `raw_receipt` preserved — that lets the submission_log row keep the
	// server response even when parsing failed.
	miint::ENASubmissionOutcome outcome;
	bool transport_failure = false;
	string transport_error;
	try {
		outcome = miint::SubmitProjectInsertOutcome(built.projects, opts, post_fn);
	} catch (const std::exception &e) {
		transport_failure = true;
		transport_error = e.what();
	}

	const bool success = !transport_failure && outcome.success;
	auto error_messages = outcome.error_messages;
	if (transport_failure && error_messages.empty()) {
		error_messages.push_back(transport_error);
	}
	RecordSubmissionLog(catalog, creds, built.projects.size(), outcome, success, error_messages);

	if (transport_failure) {
		throw InvalidInputException("INSERT INTO ena.projects failed: %s", transport_error);
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
		throw InvalidInputException("INSERT INTO ena.projects failed: %s", detail);
	}

	if (return_chunk) {
		AppendReturningRows(state.return_collection, state.return_types, built.projects, outcome);
	}
	return SinkFinalizeType::READY;
}

unique_ptr<GlobalSourceState> ENAProjectsInsert::GetGlobalSourceState(ClientContext &context) const {
	auto state = make_uniq<ENAProjectsInsertGlobalSourceState>();
	if (return_chunk && sink_state) {
		auto &gstate = sink_state->Cast<ENAProjectsInsertGlobalSinkState>();
		gstate.return_collection.InitializeScan(state->scan_state);
	}
	return std::move(state);
}

SourceResultType ENAProjectsInsert::GetDataInternal(ExecutionContext &context, DataChunk &chunk,
                                                    OperatorSourceInput &input) const {
	auto &source_state = input.global_state.Cast<ENAProjectsInsertGlobalSourceState>();
	// The executor schedules the sink pipeline before invoking source-side reads
	// for sink-source operators, so `sink_state` is always populated here.
	D_ASSERT(sink_state);
	auto &gstate = sink_state->Cast<ENAProjectsInsertGlobalSinkState>();

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

} // namespace duckdb
