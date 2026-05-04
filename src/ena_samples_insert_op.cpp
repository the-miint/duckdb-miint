// SPDX-License-Identifier: MIT
//
// Implementation of ENASamplesInsert. See ena_samples_insert_op.hpp.
//
// Mirrors ENAProjectsInsert. Shared chunk-buffering Sink, RETURNING/count
// Source machinery, credential resolution, and submission_log scaffolding
// live in ena_insert_common.{hpp,cpp}; this file focuses on the
// sample-specific data mapping.

#include "ena_samples_insert_op.hpp"

#include "ena_client.hpp"
#include "ena_insert_common.hpp"
#include "ena_samples_insert.hpp"
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

class ENASamplesInsertGlobalSinkState : public GlobalSinkState {
public:
	ENASamplesInsertGlobalSinkState(ClientContext &context, vector<LogicalType> input_types,
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

class ENASamplesInsertGlobalSourceState : public GlobalSourceState {
public:
	ColumnDataScanState scan_state;
	bool finished = false;
};

// Build SampleSpecs from the buffered ColumnDataCollection. The HOLD-action
// agreement check matches ENAProjectsInsert: V2 envelopes apply HOLD at the
// submission level, so all rows must agree if it's set on more than one.
struct BuiltSamples {
	vector<miint::SampleSpec> samples;
	string hold_until_date;
};

constexpr idx_t COL_ALIAS = 0;
constexpr idx_t COL_TITLE = 1;
constexpr idx_t COL_DESCRIPTION = 2;
constexpr idx_t COL_TAXON_ID = 3;
constexpr idx_t COL_SCIENTIFIC_NAME = 4;
constexpr idx_t COL_CHECKLIST = 5;
constexpr idx_t COL_ATTRIBUTES = 6;
constexpr idx_t COL_ERS_ACCESSION = 7;
constexpr idx_t COL_SAMEA_ACCESSION = 8;
constexpr idx_t TABLE_COLUMN_COUNT = 9;

// Pull a MAP(VARCHAR,VARCHAR) value into the (tag,value) attribute list that
// SampleSpec expects. NULL maps and empty maps both yield no attributes.
// Iteration order follows DuckDB's map-storage order, which matches the
// caller's MAP-literal order in practice; this is the same guarantee the
// envelope builder already documents for project attributes.
std::vector<std::pair<std::string, std::string>> ExtractAttributesMap(const Value &v) {
	std::vector<std::pair<std::string, std::string>> out;
	if (v.IsNull()) {
		return out;
	}
	const auto &entries = ListValue::GetChildren(v);
	out.reserve(entries.size());
	for (const auto &entry : entries) {
		const auto &kv_children = StructValue::GetChildren(entry);
		if (kv_children.size() != 2) {
			throw InvalidInputException("INSERT INTO ena.samples: 'attributes' MAP entries must have key+value pairs");
		}
		const auto key = ValueToVarchar(kv_children[0]);
		const auto value = ValueToVarchar(kv_children[1]);
		if (key.empty()) {
			throw InvalidInputException("INSERT INTO ena.samples: 'attributes' map keys must be non-empty");
		}
		out.emplace_back(key, value);
	}
	return out;
}

BuiltSamples BuildSamplesFromBuffer(ColumnDataCollection &buffer,
                                    const physical_index_vector_t<idx_t> &column_index_map) {
	BuiltSamples out;
	const idx_t input_columns = buffer.ColumnCount();
	const idx_t alias_idx = ResolveInputColumn(column_index_map, input_columns, COL_ALIAS);
	if (alias_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.samples requires the 'alias' column");
	}
	const idx_t taxon_idx = ResolveInputColumn(column_index_map, input_columns, COL_TAXON_ID);
	if (taxon_idx == DConstants::INVALID_INDEX) {
		throw InvalidInputException("INSERT INTO ena.samples requires the 'taxon_id' column");
	}
	// `column_index_map` is per-statement, not per-row; resolve every optional
	// column index up front so the row loop doesn't repeat the work.
	const idx_t title_idx = ResolveInputColumn(column_index_map, input_columns, COL_TITLE);
	const idx_t description_idx = ResolveInputColumn(column_index_map, input_columns, COL_DESCRIPTION);
	const idx_t scientific_idx = ResolveInputColumn(column_index_map, input_columns, COL_SCIENTIFIC_NAME);
	const idx_t checklist_idx = ResolveInputColumn(column_index_map, input_columns, COL_CHECKLIST);
	const idx_t attrs_idx = ResolveInputColumn(column_index_map, input_columns, COL_ATTRIBUTES);

	ColumnDataScanState scan_state;
	buffer.InitializeScan(scan_state);
	DataChunk chunk;
	buffer.InitializeScanChunk(chunk);

	while (buffer.Scan(scan_state, chunk)) {
		for (idx_t row = 0; row < chunk.size(); row++) {
			miint::SampleSpec spec;
			spec.alias = ValueToVarchar(chunk.data[alias_idx].GetValue(row));
			if (spec.alias.empty()) {
				throw InvalidInputException("INSERT INTO ena.samples: 'alias' must be non-empty");
			}
			const auto taxon_val = chunk.data[taxon_idx].GetValue(row);
			if (taxon_val.IsNull()) {
				throw InvalidInputException("INSERT INTO ena.samples: 'taxon_id' must be non-null for alias '%s'",
				                            spec.alias);
			}
			spec.taxon_id = taxon_val.GetValue<int64_t>();
			if (spec.taxon_id <= 0) {
				throw InvalidInputException("INSERT INTO ena.samples: 'taxon_id' must be > 0 for alias '%s' (got %lld)",
				                            spec.alias, static_cast<long long>(spec.taxon_id));
			}
			if (title_idx != DConstants::INVALID_INDEX) {
				spec.title = ValueToVarchar(chunk.data[title_idx].GetValue(row));
			}
			if (description_idx != DConstants::INVALID_INDEX) {
				spec.description = ValueToVarchar(chunk.data[description_idx].GetValue(row));
			}
			if (scientific_idx != DConstants::INVALID_INDEX) {
				spec.scientific_name = ValueToVarchar(chunk.data[scientific_idx].GetValue(row));
			}
			if (checklist_idx != DConstants::INVALID_INDEX) {
				spec.checklist = ValueToVarchar(chunk.data[checklist_idx].GetValue(row));
			}
			if (attrs_idx != DConstants::INVALID_INDEX) {
				spec.attributes = ExtractAttributesMap(chunk.data[attrs_idx].GetValue(row));
			}
			// HOLD on samples is not exposed in Phase 5 (no hold_until_date
			// column on ena.samples).
			out.samples.push_back(std::move(spec));
		}
	}
	return out;
}

void AppendReturningRows(ColumnDataCollection &return_collection, const vector<LogicalType> &return_types,
                         const vector<miint::SampleSpec> &samples, const miint::ENASamplesSubmissionOutcome &outcome) {
	D_ASSERT(return_types.size() == TABLE_COLUMN_COUNT);
	D_ASSERT(outcome.rows.size() == samples.size());

	DataChunk chunk;
	chunk.Initialize(Allocator::DefaultAllocator(), return_types);
	for (idx_t i = 0; i < outcome.rows.size(); i++) {
		const auto &row = outcome.rows[i];
		const auto &spec = samples[i];
		const auto idx = chunk.size();
		chunk.data[COL_ALIAS].SetValue(idx, Value(row.alias));
		chunk.data[COL_TITLE].SetValue(idx, Value(spec.title));
		chunk.data[COL_DESCRIPTION].SetValue(idx, Value(spec.description));
		chunk.data[COL_TAXON_ID].SetValue(idx, Value::INTEGER(NumericCast<int32_t>(spec.taxon_id)));
		chunk.data[COL_SCIENTIFIC_NAME].SetValue(idx, Value(spec.scientific_name));
		chunk.data[COL_CHECKLIST].SetValue(idx, Value(spec.checklist));
		// attributes column emits NULL on RETURNING in Phase 5 — the user-
		// supplied attribute list is preserved verbatim in
		// `ena.submission_log.request_payload`, and the server's echo (which
		// may include system-injected attributes the user did not write) is
		// in `submission_log.receipt`. Same trade-off as ena.projects in
		// Phase 4. Real MAP-value emission is a future task.
		chunk.data[COL_ATTRIBUTES].SetValue(idx, Value(LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)));
		chunk.data[COL_ERS_ACCESSION].SetValue(idx, Value(row.ers_accession));
		chunk.data[COL_SAMEA_ACCESSION].SetValue(idx, row.samea_accession.empty() ? Value(LogicalType::VARCHAR)
		                                                                          : Value(row.samea_accession));
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

ENASamplesInsert::ENASamplesInsert(PhysicalPlan &physical_plan, LogicalOperator &op, ENATableEntry &table_p,
                                   physical_index_vector_t<idx_t> column_index_map_p, bool return_chunk_p)
    : PhysicalOperator(physical_plan, PhysicalOperatorType::EXTENSION, op.types, 1), table(table_p),
      column_index_map(std::move(column_index_map_p)), return_chunk(return_chunk_p) {
}

unique_ptr<GlobalSinkState> ENASamplesInsert::GetGlobalSinkState(ClientContext &context) const {
	D_ASSERT(!children.empty());
	auto input_types = children[0].get().GetTypes();
	return make_uniq<ENASamplesInsertGlobalSinkState>(context, std::move(input_types), types);
}

SinkResultType ENASamplesInsert::Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const {
	auto &state = input.global_state.Cast<ENASamplesInsertGlobalSinkState>();
	lock_guard<mutex> guard(state.lock);
	state.buffered.Append(chunk);
	state.insert_count += chunk.size();
	return SinkResultType::NEED_MORE_INPUT;
}

SinkFinalizeType ENASamplesInsert::Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
                                            OperatorSinkFinalizeInput &input) const {
	auto &state = input.global_state.Cast<ENASamplesInsertGlobalSinkState>();
	if (state.buffered.Count() == 0) {
		return SinkFinalizeType::READY;
	}

	auto &catalog = table.catalog.Cast<ENACatalog>();
	auto creds = ResolveENACredentials(context, catalog);

	auto built = BuildSamplesFromBuffer(state.buffered, column_index_map);

	miint::ENASampleInsertOptions opts;
	opts.endpoint_url = creds.endpoint_url + "/submit";
	opts.user = creds.user;
	opts.password = creds.password;
	opts.hold_until_date = built.hold_until_date;

	miint::ENAClient client(*context.db);
	auto post_fn = [&client](const std::string &url, const std::string &body, const std::string &user,
	                         const std::string &password, const std::string &content_type) {
		(void)content_type;
		return client.PostJSON(url, body, user, password);
	};

	miint::ENASamplesSubmissionOutcome outcome;
	bool transport_failure = false;
	string transport_error;
	try {
		outcome = miint::SubmitSampleInsertOutcome(built.samples, opts, post_fn);
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
	log_payload.object_type = "samples";
	log_payload.action = "ADD";
	log_payload.n_objects = static_cast<int32_t>(built.samples.size());
	log_payload.success = success;
	log_payload.duration_ms = outcome.duration_ms;
	log_payload.envelope_payload = outcome.envelope_payload;
	log_payload.raw_receipt = outcome.raw_receipt;
	log_payload.era_accession = outcome.era_accession;
	log_payload.error_messages = error_messages;
	RecordSubmissionLog(catalog, creds, log_payload);

	if (transport_failure) {
		throw InvalidInputException("INSERT INTO ena.samples failed: %s", transport_error);
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
		throw InvalidInputException("INSERT INTO ena.samples failed: %s", detail);
	}

	if (return_chunk) {
		AppendReturningRows(state.return_collection, state.return_types, built.samples, outcome);
	}
	return SinkFinalizeType::READY;
}

unique_ptr<GlobalSourceState> ENASamplesInsert::GetGlobalSourceState(ClientContext &context) const {
	auto state = make_uniq<ENASamplesInsertGlobalSourceState>();
	if (return_chunk && sink_state) {
		auto &gstate = sink_state->Cast<ENASamplesInsertGlobalSinkState>();
		gstate.return_collection.InitializeScan(state->scan_state);
	}
	return std::move(state);
}

SourceResultType ENASamplesInsert::GetDataInternal(ExecutionContext &context, DataChunk &chunk,
                                                   OperatorSourceInput &input) const {
	auto &source_state = input.global_state.Cast<ENASamplesInsertGlobalSourceState>();
	D_ASSERT(sink_state);
	auto &gstate = sink_state->Cast<ENASamplesInsertGlobalSinkState>();

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
