#include "search_sequences.hpp"
#include "table_function_common.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>

namespace duckdb {

// Output column names and types for search results.
static std::vector<std::string> GetSearchOutputNames() {
	return {"read_id", "target_id",        "identity",     "matches",       "mismatches",
	        "gaps",    "alignment_length", "query_length", "target_length", "accepted"};
}

static std::vector<LogicalType> GetSearchOutputTypes() {
	// accepted: true for strong hits passing the identity threshold,
	// false for "weak" hits (vsearch returns both accepted and weak hits
	// from search_single; weak hits are near-misses that didn't pass).
	return {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::BOOLEAN};
}

// Write a batch of SearchResults into a DataChunk.
static idx_t OutputSearchResults(DataChunk &output, const std::vector<miint::SearchResult> &results, idx_t offset,
                                 idx_t count) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = 0;

	auto &read_id_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(read_id_vec)[i] =
		    StringVector::AddString(read_id_vec, results[offset + i].query_id);
	}

	auto &target_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(target_vec)[i] =
		    StringVector::AddString(target_vec, results[offset + i].target_id);
	}

	auto identity_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		identity_data[i] = results[offset + i].identity;
	}

	auto matches_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		matches_data[i] = results[offset + i].matches;
	}

	auto mismatches_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		mismatches_data[i] = results[offset + i].mismatches;
	}

	auto gaps_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		gaps_data[i] = results[offset + i].gaps;
	}

	auto alnlen_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		alnlen_data[i] = results[offset + i].alignment_length;
	}

	auto qlen_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		qlen_data[i] = results[offset + i].query_length;
	}

	auto tlen_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		tlen_data[i] = results[offset + i].target_length;
	}

	auto accepted_data = FlatVector::GetData<bool>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		accepted_data[i] = results[offset + i].accepted;
	}

	D_ASSERT(col == output.ColumnCount());
	output.SetCardinality(actual);
	return actual;
}

unique_ptr<FunctionData> SearchSequencesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                            vector<LogicalType> &return_types,
                                                            vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->query_table = input.inputs[0].GetValue<std::string>();

	auto db_it = input.named_parameters.find("db");
	if (db_it == input.named_parameters.end()) {
		throw BinderException("search_sequences_vsearch requires 'db' parameter (reference table name)");
	}
	data->ref_table = db_it->second.GetValue<std::string>();

	// id is required — no silent 0.0 default that accepts everything.
	auto id_it = input.named_parameters.find("id");
	if (id_it == input.named_parameters.end()) {
		throw BinderException("search_sequences_vsearch requires 'id' parameter (identity threshold, e.g. id:=0.97)");
	}
	data->params.id = id_it->second.GetValue<double>();
	if (data->params.id < 0.0 || data->params.id > 1.0) {
		throw InvalidInputException("search_sequences_vsearch: id must be between 0.0 and 1.0 (got %g)",
		                            data->params.id);
	}

	data->query_schema = ValidateSequenceTableSchema(context, data->query_table);

	// Validate ref table exists and has required columns
	ValidateSequenceTableSchema(context, data->ref_table);

	auto get_int = [&](const std::string &name, int &out, int min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val) {
				throw InvalidInputException("search_sequences_vsearch: %s must be %s (got %d)", name, constraint, out);
			}
		}
	};

	get_int("maxaccepts", data->params.maxaccepts, 1, ">= 1");
	get_int("maxrejects", data->params.maxrejects, 1, ">= 1");

	data->names = GetSearchOutputNames();
	data->types = GetSearchOutputTypes();
	for (auto &n : data->names) {
		names.push_back(n);
	}
	for (auto &t : data->types) {
		return_types.push_back(t);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> SearchSequencesTableFunction::InitGlobal(ClientContext &context,
                                                                              TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	gstate->wrapper = miint::VsearchSearchWrapper(data.params);

	auto ref = LoadSingleEndSequences(context, data.ref_table, "search_sequences_vsearch");
	gstate->wrapper.set_database(ref.labels, ref.sequences);

	// Lazy streaming reader — one STANDARD_VECTOR_SIZE chunk per Execute call.
	gstate->query_stream =
	    std::make_unique<QuerySequenceStream>(context, data.query_table, data.query_schema, STANDARD_VECTOR_SIZE);

	return gstate;
}

void SearchSequencesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	while (true) {
		// 1. Drain buffered results from last batch
		if (gstate.result_offset < gstate.result_buffer.size()) {
			idx_t remaining = gstate.result_buffer.size() - gstate.result_offset;
			idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			OutputSearchResults(output, gstate.result_buffer, gstate.result_offset, count);
			gstate.result_offset += count;
			return;
		}

		// 2. Buffer exhausted — fetch next chunk and search via batch API.
		gstate.result_buffer.clear();
		gstate.result_offset = 0;

		auto query_batch = gstate.query_stream->FetchSubBatch();
		if (query_batch.empty()) {
			output.SetCardinality(0);
			return;
		}

		// 3. Search entire chunk via vsearch's internal thread pool.
		gstate.wrapper.search_batch(query_batch.read_ids, query_batch.sequences1, gstate.result_buffer);
	}
}

TableFunction SearchSequencesTableFunction::GetFunction() {
	auto tf = TableFunction("search_sequences_vsearch", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);

	tf.named_parameters["db"] = LogicalType::VARCHAR;
	tf.named_parameters["id"] = LogicalType::DOUBLE;
	tf.named_parameters["maxaccepts"] = LogicalType::INTEGER;
	tf.named_parameters["maxrejects"] = LogicalType::INTEGER;

	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void SearchSequencesTableFunction::Register(ExtensionLoader &loader) {
	auto tf = GetFunction();
	loader.RegisterFunction(tf);
}

} // namespace duckdb
