#include "uchime_ref.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>

namespace duckdb {

static constexpr idx_t UCHIME_MAX_THREADS = 8;

// Output column names and types for the 18-column uchimeout format.
static std::vector<std::string> GetUchimeOutputNames() {
	return {"score",        "query",      "parent_a", "parent_b",      "closest_parent", "id_query_model",
	        "id_query_a",   "id_query_b", "id_a_b",   "id_query_top",  "left_yes",       "left_no",
	        "left_abstain", "right_yes",  "right_no", "right_abstain", "divergence",     "flag"};
}

static std::vector<LogicalType> GetUchimeOutputTypes() {
	return {LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	        LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
	        LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::DOUBLE,  LogicalType::VARCHAR};
}

// Write a batch of UchimeResults into a DataChunk.
static idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset,
                                 idx_t count) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = 0;

	// score
	auto score_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		score_data[i] = results[offset + i].score;
	}

	// query, parent_a, parent_b, closest_parent (VARCHAR columns)
	auto &query_vec = output.data[col++];
	auto &parent_a_vec = output.data[col++];
	auto &parent_b_vec = output.data[col++];
	auto &closest_parent_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		FlatVector::GetData<string_t>(query_vec)[i] = StringVector::AddString(query_vec, r.query_label);
		FlatVector::GetData<string_t>(parent_a_vec)[i] = StringVector::AddString(parent_a_vec, r.parent_a_label);
		FlatVector::GetData<string_t>(parent_b_vec)[i] = StringVector::AddString(parent_b_vec, r.parent_b_label);
		FlatVector::GetData<string_t>(closest_parent_vec)[i] =
		    StringVector::AddString(closest_parent_vec, r.closest_parent_label);
	}

	// id_query_model, id_query_a, id_query_b, id_a_b, id_query_top (DOUBLE columns)
	auto id_qm = FlatVector::GetData<double>(output.data[col++]);
	auto id_qa = FlatVector::GetData<double>(output.data[col++]);
	auto id_qb = FlatVector::GetData<double>(output.data[col++]);
	auto id_ab = FlatVector::GetData<double>(output.data[col++]);
	auto id_qt = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		id_qm[i] = r.id_query_model;
		id_qa[i] = r.id_query_a;
		id_qb[i] = r.id_query_b;
		id_ab[i] = r.id_a_b;
		id_qt[i] = r.id_query_top;
	}

	// left_yes, left_no, left_abstain, right_yes, right_no, right_abstain (INTEGER columns)
	auto ly = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ln = FlatVector::GetData<int32_t>(output.data[col++]);
	auto la = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ry = FlatVector::GetData<int32_t>(output.data[col++]);
	auto rn = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ra = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		ly[i] = r.left_yes;
		ln[i] = r.left_no;
		la[i] = r.left_abstain;
		ry[i] = r.right_yes;
		rn[i] = r.right_no;
		ra[i] = r.right_abstain;
	}

	// divergence
	auto div_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		div_data[i] = results[offset + i].divergence;
	}

	// flag
	auto &flag_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(flag_vec)[i] = StringVector::AddString(flag_vec, results[offset + i].flag);
	}

	D_ASSERT(col == output.ColumnCount());
	output.SetCardinality(actual);
	return actual;
}

unique_ptr<FunctionData> UchimeRefTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<Data>();

	// Positional argument: query table name
	data->query_table = input.inputs[0].GetValue<std::string>();

	// Required named parameter: db (reference table name)
	auto db_it = input.named_parameters.find("db");
	if (db_it == input.named_parameters.end()) {
		throw BinderException("uchime_ref requires 'db' parameter (reference table name)");
	}
	data->ref_table = db_it->second.GetValue<std::string>();

	// Validate both tables exist and have required columns (read_id, sequence1).
	// This catches schema errors at bind time, not execution time.
	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);
	data->ref_schema = ValidateSequenceTableSchema(context, data->ref_table);

	// Optional scoring parameters with range validation
	auto get_double = [&](const std::string &name, double &out, double min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<double>();
			if (out < min_val) {
				throw InvalidInputException("uchime_ref: %s must be %s (got %g)", name, constraint, out);
			}
		}
	};
	auto get_int = [&](const std::string &name, int &out, int min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val) {
				throw InvalidInputException("uchime_ref: %s must be %s (got %d)", name, constraint, out);
			}
		}
	};

	get_double("minh", data->params.minh, 0.0, ">= 0");
	get_double("xn", data->params.xn, 1.0, ">= 1.0");
	get_double("dn", data->params.dn, 0.0, ">= 0");
	get_double("mindiv", data->params.mindiv, 0.0, ">= 0");
	get_int("mindiffs", data->params.mindiffs, 1, ">= 1");

	// Set output schema
	data->names = GetUchimeOutputNames();
	data->types = GetUchimeOutputTypes();
	for (auto &n : data->names) {
		names.push_back(n);
	}
	for (auto &t : data->types) {
		return_types.push_back(t);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> UchimeRefTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	gstate->detector = miint::ChimeraDetector(data.params);

	// Load reference sequences via a separate connection (avoids deadlocking the current context).
	// NOTE: The entire reference database is held in memory. For large reference
	// databases (e.g., full SILVA at 2GB+), ensure sufficient process memory.
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection ref_conn(db);
	auto ref_result =
	    ref_conn.Query("SELECT read_id, sequence1 FROM " + KeywordHelper::WriteOptionallyQuoted(data.ref_table));
	if (ref_result->HasError()) {
		throw InvalidInputException("Failed to read reference table '%s': %s", data.ref_table, ref_result->GetError());
	}

	std::vector<std::string> ref_labels;
	std::vector<std::string> ref_sequences;
	auto &materialized = ref_result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto read_id_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			if (read_id_val.IsNull() || seq_val.IsNull()) {
				continue; // Skip rows with NULL read_id or sequence1
			}
			auto seq_str = seq_val.GetValue<std::string>();
			if (seq_str.empty()) {
				continue; // Skip empty sequences
			}
			ref_labels.push_back(read_id_val.GetValue<std::string>());
			ref_sequences.push_back(std::move(seq_str));
		}
	}

	if (ref_labels.empty()) {
		throw InvalidInputException("Reference table '%s' is empty (or contains only NULL/empty sequences)",
		                            data.ref_table);
	}

	gstate->detector.set_reference(ref_labels, ref_sequences);

	// Set up lazy streaming reader for query sequences using validated schema.
	gstate->query_stream = std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema);

	gstate->num_threads = std::min(static_cast<idx_t>(std::thread::hardware_concurrency()), UCHIME_MAX_THREADS);
	if (gstate->num_threads == 0) {
		gstate->num_threads = 4;
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> UchimeRefTableFunction::InitLocal(ExecutionContext &context,
                                                                      TableFunctionInitInput &input,
                                                                      GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void UchimeRefTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();

	while (true) {
		// 1. Drain buffered results
		if (lstate.buffer_offset < lstate.result_buffer.size()) {
			idx_t remaining = lstate.result_buffer.size() - lstate.buffer_offset;
			idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputUchimeResults(output, lstate.result_buffer, lstate.buffer_offset, count);
			lstate.buffer_offset += count;
			return;
		}

		// 2. Buffer exhausted — fetch next sub-batch from stream (thread-safe)
		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Run chimera detection on each query in the sub-batch
		lstate.result_buffer.reserve(query_batch.size());
		for (idx_t i = 0; i < query_batch.size(); i++) {
			auto result = gstate.detector.detect(query_batch.read_ids[i], query_batch.sequences1[i], lstate.aligner);
			lstate.result_buffer.push_back(std::move(result));
		}
	}
}

TableFunction UchimeRefTableFunction::GetFunction() {
	auto tf = TableFunction("uchime_ref", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	tf.named_parameters["db"] = LogicalType::VARCHAR;
	tf.named_parameters["minh"] = LogicalType::DOUBLE;
	tf.named_parameters["xn"] = LogicalType::DOUBLE;
	tf.named_parameters["dn"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiv"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiffs"] = LogicalType::INTEGER;

	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void UchimeRefTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
