#include "cluster_sequences.hpp"
#include "id_column_utils.hpp"
#include "sequence_table_reader.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <algorithm>

namespace duckdb {

static std::vector<std::string> GetClusterOutputNames() {
	return {"read_id", "is_centroid", "cluster_id", "centroid_id", "identity", "cigar", "cigar_truncated"};
}

// read_id (col 0) and centroid_id (col 3) mirror the input read_id type; the
// centroid_id values are themselves drawn from the input read_ids, so the two
// id columns always share one type.
static std::vector<LogicalType> GetClusterOutputTypes(const LogicalType &id_type) {
	return {id_type,
	        LogicalType::BOOLEAN,
	        LogicalType::INTEGER,
	        id_type,
	        LogicalType::DOUBLE,
	        LogicalType::VARCHAR,
	        LogicalType::BOOLEAN};
}

static idx_t OutputClusterResults(DataChunk &output, const std::vector<miint::ClusterResult> &results, idx_t offset,
                                  idx_t count, const LogicalType &id_type) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = 0;

	auto &read_id_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		EmitIdCell(read_id_vec, i, results[offset + i].read_id, id_type);
	}

	auto is_centroid_data = FlatVector::GetData<bool>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		is_centroid_data[i] = results[offset + i].is_centroid;
	}

	auto cluster_id_data = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		cluster_id_data[i] = results[offset + i].cluster_id;
	}

	auto &centroid_id_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		EmitIdCell(centroid_id_vec, i, results[offset + i].centroid_id, id_type);
	}

	auto identity_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		identity_data[i] = results[offset + i].identity;
	}

	auto &cigar_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(cigar_vec)[i] = StringVector::AddString(cigar_vec, results[offset + i].cigar);
	}

	auto cigar_trunc_data = FlatVector::GetData<bool>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		cigar_trunc_data[i] = results[offset + i].cigar_truncated;
	}

	D_ASSERT(col == output.ColumnCount());
	output.SetCardinality(actual);
	return actual;
}

unique_ptr<FunctionData> ClusterSequencesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<LogicalType> &return_types,
                                                             vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->input_table = input.inputs[0].GetValue<std::string>();

	// Validate table has read_id (VARCHAR, BIGINT, or UUID) and sequence1
	// (VARCHAR). Capture the read_id type so the output id columns mirror it.
	auto schema = ValidateSequenceTableSchema(context, data->input_table, /*allow_bigint=*/true);
	data->id_type = schema.id_type;

	// id is required — no silent default.
	auto id_it = input.named_parameters.find("id");
	if (id_it == input.named_parameters.end()) {
		throw BinderException("cluster_sequences_vsearch requires 'id' parameter (identity threshold, e.g. id:=0.97)");
	}
	data->params.id = id_it->second.GetValue<double>();
	if (data->params.id < 0.0 || data->params.id > 1.0) {
		throw InvalidInputException("cluster_sequences_vsearch: id must be between 0.0 and 1.0 (got %g)",
		                            data->params.id);
	}

	// strand: 'plus' (default) or 'both'
	auto strand_it = input.named_parameters.find("strand");
	if (strand_it != input.named_parameters.end()) {
		auto strand_str = strand_it->second.GetValue<std::string>();
		if (strand_str == "plus") {
			data->params.strand = 1;
		} else if (strand_str == "both") {
			data->params.strand = 2;
		} else {
			throw InvalidInputException("cluster_sequences_vsearch: strand must be 'plus' or 'both' (got '%s')",
			                            strand_str);
		}
	}

	auto get_int = [&](const std::string &name, int &out, int min_val, int max_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val || out > max_val) {
				throw InvalidInputException("cluster_sequences_vsearch: %s must be %s (got %d)", name, constraint, out);
			}
		}
	};

	// Default to DuckDB's configured thread count so `SET threads=N` is honored.
	// Explicit `threads:=N` overrides. Capped at vsearch's CLI ceiling (n_threads_max
	// = 1024 in ext/vsearch/src/vsearch.cc) so we fail with a clear message rather
	// than blowing up inside pthread_create.
	data->params.threads = NumericCast<int>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	get_int("threads", data->params.threads, 1, 1024, ">= 1 and <= 1024");

	data->names = GetClusterOutputNames();
	data->types = GetClusterOutputTypes(data->id_type);
	for (auto &n : data->names) {
		names.push_back(n);
	}
	for (auto &t : data->types) {
		return_types.push_back(t);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ClusterSequencesTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	// Clustering is materialized upfront: vsearch's clustering algorithm is
	// inherently sequential (each centroid must be indexed before the next
	// sequence is processed), and the entire session must complete before
	// releasing the vsearch mutex. Execute() paginates the results.
	miint::VsearchClusterWrapper wrapper(data.params);

	auto loaded = LoadSingleEndSequences(context, data.input_table, "cluster_sequences_vsearch", /*strict=*/true);
	wrapper.set_sequences(loaded.labels, loaded.sequences);
	gstate->results = wrapper.cluster_all();

	return gstate;
}

void ClusterSequencesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	if (gstate.result_offset >= gstate.results.size()) {
		output.SetCardinality(0);
		return;
	}

	idx_t remaining = gstate.results.size() - gstate.result_offset;
	idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	OutputClusterResults(output, gstate.results, gstate.result_offset, count, data.id_type);
	gstate.result_offset += count;
}

TableFunction ClusterSequencesTableFunction::GetFunction() {
	auto tf = TableFunction("cluster_sequences_vsearch", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);

	tf.named_parameters["id"] = LogicalType::DOUBLE;
	tf.named_parameters["strand"] = LogicalType::VARCHAR;
	tf.named_parameters["threads"] = LogicalType::INTEGER;

	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void ClusterSequencesTableFunction::Register(ExtensionLoader &loader) {
	auto tf = GetFunction();
	loader.RegisterFunction(tf);
}

} // namespace duckdb
