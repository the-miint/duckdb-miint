#include "uchime_ref.hpp"
#include "table_function_common.hpp"
#include "uchime_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>

namespace duckdb {

unique_ptr<FunctionData> UchimeRefTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->query_table = input.inputs[0].GetValue<std::string>();

	auto db_it = input.named_parameters.find("db");
	if (db_it == input.named_parameters.end()) {
		throw BinderException("uchime_ref requires 'db' parameter (reference table name)");
	}
	data->ref_table = db_it->second.GetValue<std::string>();

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);
	data->ref_schema = ValidateSequenceTableSchema(context, data->ref_table);

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

	gstate->wrapper = miint::VsearchChimeraWrapper(data.params);
	gstate->query_table = data.query_table;
	gstate->query_schema = data.query_schema;

	// Load reference sequences via a separate connection.
	// NOTE: The entire reference database is held in memory.
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
				continue;
			}
			auto seq_str = seq_val.GetValue<std::string>();
			if (seq_str.empty()) {
				continue;
			}
			ref_labels.push_back(read_id_val.GetValue<std::string>());
			ref_sequences.push_back(std::move(seq_str));
		}
	}

	if (ref_labels.empty()) {
		throw InvalidInputException("Reference table '%s' is empty (or contains only NULL/empty sequences)",
		                            data.ref_table);
	}

	gstate->wrapper.set_reference(ref_labels, ref_sequences);

	// Materialize all query IDs (lightweight — just strings, no sequences).
	// This allows lock-free atomic batch claiming in Execute().
	// NULL read_ids are rejected — every query must have an identifier.
	Connection id_conn(db);
	auto id_result = id_conn.Query("SELECT read_id FROM " + KeywordHelper::WriteOptionallyQuoted(data.query_table));
	if (id_result->HasError()) {
		throw InvalidInputException("Failed to read query table '%s': %s", data.query_table, id_result->GetError());
	}
	auto &id_mat = id_result->Cast<MaterializedQueryResult>();
	while (auto chunk = id_mat.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto val = chunk->GetValue(0, i);
			if (val.IsNull()) {
				throw InvalidInputException("Query table '%s' contains NULL read_id values", data.query_table);
			}
			gstate->all_query_ids.push_back(val.GetValue<std::string>());
		}
	}

	return gstate;
}

unique_ptr<LocalTableFunctionState> UchimeRefTableFunction::InitLocal(ExecutionContext &context,
                                                                      TableFunctionInitInput &input,
                                                                      GlobalTableFunctionState *global_state) {
	auto &gstate = global_state->Cast<GlobalState>();
	auto lstate = make_uniq<LocalState>();
	lstate->handle.emplace(gstate.wrapper.create_detect_handle());
	return lstate;
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

		// 2. Buffer exhausted — claim next batch of query IDs atomically (no mutex).
		lstate.result_buffer.clear();
		lstate.buffer_offset = 0;

		idx_t offset = gstate.next_batch_offset.fetch_add(BATCH_SIZE);
		if (offset >= gstate.all_query_ids.size()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Fetch sequences for claimed IDs via ReadBatchByIds (temp table JOIN).
		miint::SequenceRecordBatch query_batch;
		ReadBatchByIds(context, gstate.query_table, gstate.query_schema, gstate.all_query_ids, offset, BATCH_SIZE,
		               query_batch);

		if (query_batch.empty()) {
			// Batch claimed valid IDs but join returned nothing — skip and try next batch.
			// This can happen if the table was modified between InitGlobal and Execute.
			continue;
		}

		// 4. Run chimera detection on each query in the batch
		lstate.result_buffer.reserve(query_batch.size());
		for (idx_t i = 0; i < query_batch.size(); i++) {
			auto result = lstate.handle->detect(query_batch.read_ids[i], query_batch.sequences1[i]);
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
