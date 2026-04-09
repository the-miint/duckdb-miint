#include "uchime_ref.hpp"
#include "table_function_common.hpp"
#include "uchime_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"

#include <algorithm>

namespace duckdb {

unique_ptr<FunctionData> UchimeRefTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->query_table = input.inputs[0].GetValue<std::string>();

	auto db_it = input.named_parameters.find("db");
	if (db_it == input.named_parameters.end()) {
		throw BinderException("detect_chimera_uchime requires 'db' parameter (reference table name)");
	}
	data->ref_table = db_it->second.GetValue<std::string>();

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);
	data->ref_schema = ValidateSequenceTableSchema(context, data->ref_table);

	auto get_double = [&](const std::string &name, double &out, double min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<double>();
			if (out < min_val) {
				throw InvalidInputException("detect_chimera_uchime: %s must be %s (got %g)", name, constraint, out);
			}
		}
	};
	auto get_int = [&](const std::string &name, int &out, int min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val) {
				throw InvalidInputException("detect_chimera_uchime: %s must be %s (got %d)", name, constraint, out);
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

	auto ref = LoadSingleEndSequences(context, data.ref_table, "detect_chimera_uchime");
	gstate->wrapper.set_reference(ref.labels, ref.sequences);

	// Lazy streaming reader — one STANDARD_VECTOR_SIZE chunk per Execute call.
	gstate->query_stream =
	    std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema, STANDARD_VECTOR_SIZE);

	return gstate;
}

void UchimeRefTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	while (true) {
		// 1. Drain buffered results from last batch
		if (gstate.result_offset < gstate.result_buffer.size()) {
			idx_t remaining = gstate.result_buffer.size() - gstate.result_offset;
			idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputUchimeResults(output, gstate.result_buffer, gstate.result_offset, count);
			gstate.result_offset += count;
			return;
		}

		// 2. Buffer exhausted — fetch next chunk and detect via batch API.
		gstate.result_buffer.clear();
		gstate.result_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Detect chimeras for entire chunk via vsearch's internal thread pool.
		gstate.wrapper.detect_batch(query_batch.read_ids, query_batch.sequences1, gstate.result_buffer);
	}
}

TableFunction UchimeRefTableFunction::GetFunction() {
	auto tf = TableFunction("detect_chimera_uchime", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);

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
	auto tf = GetFunction();
	loader.RegisterFunction(tf);
}

} // namespace duckdb
