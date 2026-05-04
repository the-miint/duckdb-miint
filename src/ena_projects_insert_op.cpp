// SPDX-License-Identifier: MIT
//
// Implementation of ENAProjectsInsert. See ena_projects_insert_op.hpp.
//
// The chunk-buffering Sink and the RETURNING/count Source path are shared
// with ENASamplesInsert (and future per-table operators) via the helpers in
// ena_insert_common.{hpp,cpp}. Keep this file focused on the project-specific
// data mapping (chunk row → ProjectSpec, RETURNING row population).

#include "ena_projects_insert_op.hpp"

#include "ena_client.hpp"
#include "ena_insert_common.hpp"
#include "ena_projects_insert.hpp"
#include "ena_storage.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/mutex.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/common/types/column/column_data_scan_states.hpp"
#include "duckdb/common/types/date.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"
#include "duckdb/parallel/event.hpp"
#include "duckdb/planner/operator/logical_insert.hpp"

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
	// `column_index_map` is per-statement, not per-row; resolve every optional
	// column index up front so the row loop doesn't repeat the work.
	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
	const idx_t description_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESCRIPTION);
	const idx_t project_type_idx = ResolveInputColumn(column_index_map, input_columns, COL_PROJECT_TYPE);
	const idx_t hold_idx = ResolveInputColumn(column_index_map, input_columns, COL_HOLD_UNTIL_DATE);
	const idx_t umbrella_idx = ResolveInputColumn(column_index_map, input_columns, COL_UMBRELLA_CHILDREN);
	const idx_t attrs_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTES);

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
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			if (description_idx != DConstants::INVALID_INDEX) {
				spec.description = ValueToVarchar(chunk.data[description_idx].GetValue(row));
			}
			if (project_type_idx != DConstants::INVALID_INDEX) {
				spec.project_type = ValueToVarchar(chunk.data[project_type_idx].GetValue(row));
			}
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
			if (umbrella_idx != DConstants::INVALID_INDEX) {
				const auto v = chunk.data[umbrella_idx].GetValue(row);
				if (!v.IsNull() && !ListValue::GetChildren(v).empty()) {
					throw NotImplementedException(
					    "INSERT INTO ena.projects: 'umbrella_children' is not supported in this build");
				}
			}
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

} // namespace

ENAProjectsInsert::ENAProjectsInsert(PhysicalPlan &physical_plan, LogicalOperator &op, ENATableEntry &table_p,
                                     physical_index_vector_t<idx_t> column_index_map_p, bool return_chunk_p)
    : PhysicalOperator(physical_plan, PhysicalOperatorType::EXTENSION, op.types, 1), table(table_p),
      column_index_map(std::move(column_index_map_p)), return_chunk(return_chunk_p) {
}

unique_ptr<GlobalSinkState> ENAProjectsInsert::GetGlobalSinkState(ClientContext &context) const {
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
	auto creds = ResolveENACredentials(context, catalog);

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

	SubmissionLogPayload log_payload;
	log_payload.object_type = "projects";
	log_payload.action = "ADD";
	log_payload.n_objects = static_cast<int32_t>(built.projects.size());
	log_payload.success = success;
	log_payload.duration_ms = outcome.duration_ms;
	log_payload.envelope_payload = outcome.envelope_payload;
	log_payload.raw_receipt = outcome.raw_receipt;
	log_payload.era_accession = outcome.era_accession;
	log_payload.error_messages = error_messages;
	RecordSubmissionLog(catalog, creds, log_payload);

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
