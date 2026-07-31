#include "deblur_table_function.hpp"
#include "catalog_utils.hpp"
#include "deblur.hpp"
#include "per_sample_table_function.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {

struct DeblurData : public TableFunctionData {
	std::string input_table;
	miint::DeblurParams params;

	// Column-name overrides for the input table. Defaults preserve the historical
	// names so existing call sites work unchanged. Allow chaining (e.g. align_mafft
	// → deblur with sequence_col := 'aligned_sequence') without an intermediate rename.
	std::string id_col = "read_id";
	std::string sequence_col = "sequence1";
	std::string count_col = "abundance";

	bool has_sample_id = false;
	PerSampleBindInfo sample_info;
};

struct DeblurGlobalState : public PerSampleGlobalState {
	// Non-sample path only: pre-computed at InitGlobal, immutable thereafter.
	// The single Execute thread (enforced by max_threads=1) reads it via its lstate cursor.
	std::vector<miint::DeblurResult> shared_results;
};

struct DeblurLocalState : public LocalTableFunctionState {
	unique_ptr<Connection> conn; // sample_id path only

	// Per-thread state. Non-sample path: cursor into gstate.shared_results.
	// Sample path: own buffer for the currently-claimed sample; cursor is reset per sample.
	std::vector<miint::DeblurResult> sample_results;
	idx_t current_row = 0;
	Value sample_value;
};

// Validate that the table/view has the resolved id (VARCHAR), sequence (VARCHAR),
// and count (any integer type) columns. Names are resolved from named parameters
// at bind time.
static void ValidateDeblurTableSchema(ClientContext &context, const std::string &table_name, const std::string &id_col,
                                      const std::string &sequence_col, const std::string &count_col) {
	auto cols = GetTableOrViewColumns(context, table_name, "Deblur input");

	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < cols.names.size(); i++) {
		name_to_idx[StringUtil::Lower(cols.names[i])] = i;
	}

	auto it = name_to_idx.find(StringUtil::Lower(id_col));
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column '%s' (VARCHAR)", table_name, id_col);
	}
	if (cols.types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("deblur: column '%s' in table '%s' must be VARCHAR", id_col, table_name);
	}

	it = name_to_idx.find(StringUtil::Lower(sequence_col));
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column '%s' (VARCHAR)", table_name, sequence_col);
	}
	if (cols.types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("deblur: column '%s' in table '%s' must be VARCHAR", sequence_col, table_name);
	}

	it = name_to_idx.find(StringUtil::Lower(count_col));
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column '%s' (integer type)", table_name, count_col);
	}
	auto atype = cols.types[it->second].id();
	if (atype != LogicalTypeId::INTEGER && atype != LogicalTypeId::BIGINT && atype != LogicalTypeId::SMALLINT &&
	    atype != LogicalTypeId::TINYINT && atype != LogicalTypeId::HUGEINT && atype != LogicalTypeId::UINTEGER &&
	    atype != LogicalTypeId::UBIGINT && atype != LogicalTypeId::USMALLINT && atype != LogicalTypeId::UTINYINT &&
	    atype != LogicalTypeId::UHUGEINT) {
		throw BinderException("deblur: column '%s' in table '%s' must be an integer type (got %s)", count_col,
		                      table_name, cols.types[it->second].ToString());
	}
}

static unique_ptr<FunctionData> DeblurBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<DeblurData>();
	data->input_table = input.inputs[0].GetValue<std::string>();

	auto get_col_override = [&](const std::string &param_name, std::string &out) {
		auto it = input.named_parameters.find(param_name);
		if (it != input.named_parameters.end()) {
			auto val = it->second.GetValue<std::string>();
			if (val.empty()) {
				throw InvalidInputException("deblur: %s must not be empty", param_name);
			}
			out = std::move(val);
		}
	};
	get_col_override("id_col", data->id_col);
	get_col_override("sequence_col", data->sequence_col);
	get_col_override("count_col", data->count_col);

	ValidateDeblurTableSchema(context, data->input_table, data->id_col, data->sequence_col, data->count_col);

	// Extract named parameters
	{
		auto it = input.named_parameters.find("mean_error");
		if (it != input.named_parameters.end()) {
			data->params.mean_error = it->second.GetValue<double>();
			if (data->params.mean_error <= 0.0 || data->params.mean_error >= 1.0) {
				throw InvalidInputException("deblur: mean_error must be > 0 and < 1 (got %g)", data->params.mean_error);
			}
		}
	}

	{
		auto it = input.named_parameters.find("indel_prob");
		if (it != input.named_parameters.end()) {
			data->params.indel_prob = it->second.GetValue<double>();
			if (data->params.indel_prob < 0.0 || data->params.indel_prob > 1.0) {
				throw InvalidInputException("deblur: indel_prob must be >= 0 and <= 1 (got %g)",
				                            data->params.indel_prob);
			}
		}
	}

	{
		auto it = input.named_parameters.find("indel_max");
		if (it != input.named_parameters.end()) {
			data->params.indel_max = it->second.GetValue<int>();
			if (data->params.indel_max < 0) {
				throw InvalidInputException("deblur: indel_max must be >= 0 (got %d)", data->params.indel_max);
			}
		}
	}

	// error_profile as LIST(DOUBLE)
	{
		auto it = input.named_parameters.find("error_profile");
		if (it != input.named_parameters.end()) {
			auto &list_val = it->second;
			auto children = ListValue::GetChildren(list_val);
			if (children.empty()) {
				throw InvalidInputException("deblur: error_profile must be a non-empty list");
			}
			data->params.error_dist.reserve(children.size());
			for (idx_t i = 0; i < children.size(); i++) {
				double v = children[i].GetValue<double>();
				if (v < 0.0) {
					throw InvalidInputException(
					    "deblur: error_profile values must be non-negative (got %g at index %d)", v,
					    static_cast<int>(i));
				}
				data->params.error_dist.push_back(v);
			}
		}
	}

	auto sample_it = input.named_parameters.find("sample_id");
	if (sample_it != input.named_parameters.end()) {
		data->has_sample_id = true;
		data->sample_info.sample_id_col = sample_it->second.GetValue<string>();
	}

	if (data->has_sample_id) {
		auto conn = MakeReadOnlyHelperConnection(context);
		// Reserved output-column names the sample_id column must not collide with.
		DiscoverSamples(conn, data->input_table, data->sample_info.sample_id_col, {"read_id", "sequence", "abundance"},
		                "deblur", data->sample_info);

		names.push_back(data->sample_info.sample_id_col);
		return_types.push_back(data->sample_info.sample_id_type);
	}

	names.emplace_back("read_id");
	names.emplace_back("sequence");
	names.emplace_back("abundance");
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::BIGINT);

	return data;
}

static unique_ptr<GlobalTableFunctionState> DeblurInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<DeblurData>();
	auto gstate = make_uniq<DeblurGlobalState>();

	if (data.has_sample_id) {
		// Re-entrant algorithm: default hint lets DuckDB threads work on distinct samples.
		InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size());
		return gstate;
	}

	gstate->max_threads = 1;

	// Non-sample path: load whole table once, deblur, hold for single-threaded drain.
	auto conn = MakeReadOnlyHelperConnection(context);
	auto q_id = KeywordHelper::WriteOptionallyQuoted(data.id_col);
	auto q_seq = KeywordHelper::WriteOptionallyQuoted(data.sequence_col);
	auto q_count = KeywordHelper::WriteOptionallyQuoted(data.count_col);
	auto result = conn.Query("SELECT " + q_id + ", " + q_seq + ", CAST(" + q_count + " AS BIGINT) FROM " +
	                         KeywordHelper::WriteOptionallyQuoted(data.input_table) + " ORDER BY " + q_count + " DESC");
	if (result->HasError()) {
		throw InvalidInputException("deblur: failed to read table '%s': %s", data.input_table, result->GetError());
	}

	std::vector<miint::DeblurSequence> sequences;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto read_id_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			auto abundance_val = chunk->GetValue(2, i);
			if (read_id_val.IsNull() || seq_val.IsNull() || abundance_val.IsNull()) {
				continue;
			}
			auto seq_str = seq_val.GetValue<std::string>();
			if (seq_str.empty()) {
				continue;
			}
			sequences.push_back({read_id_val.GetValue<std::string>(), std::move(seq_str),
			                     static_cast<double>(abundance_val.GetValue<int64_t>())});
		}
	}

	if (!sequences.empty()) {
		try {
			gstate->shared_results = miint::deblur(std::move(sequences), data.params);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("deblur: %s", e.what());
		}
	}

	return gstate;
}

static unique_ptr<LocalTableFunctionState> DeblurInitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                           GlobalTableFunctionState * /*global_state*/) {
	auto &data = input.bind_data->Cast<DeblurData>();
	auto lstate = make_uniq<DeblurLocalState>();
	if (data.has_sample_id) {
		lstate->conn = make_uniq<Connection>(DatabaseInstance::GetDatabase(context.client));
		InheritTempObjects(context.client, *lstate->conn);
	}
	return lstate;
}

// Run deblur for a single sample on the per-thread connection, returning the result set.
// Per-sample partitioning uses `CAST(col AS VARCHAR) = CAST(literal AS VARCHAR)`, mirroring
// massql and woltka_ogu. This is correct for VARCHAR, integer, and date types — the common
// sample_id choices. DECIMAL sample_id columns are not recommended: DuckDB's VARCHAR cast
// preserves trailing zeros (e.g. `3.500`) while `ToSQLString()` strips them (`3.5`), so
// equality can silently miss rows. If DECIMAL sample IDs become a requirement we should
// switch this (and the other per-sample call sites) to type-aware literals.
static std::vector<miint::DeblurResult> RunDeblurForSample(Connection &conn, const DeblurData &data,
                                                           const Value &sample_value) {
	auto q_src = KeywordHelper::WriteOptionallyQuoted(data.input_table);
	auto q_col = KeywordHelper::WriteOptionallyQuoted(data.sample_info.sample_id_col);
	auto q_id = KeywordHelper::WriteOptionallyQuoted(data.id_col);
	auto q_seq = KeywordHelper::WriteOptionallyQuoted(data.sequence_col);
	auto q_count = KeywordHelper::WriteOptionallyQuoted(data.count_col);
	auto sample_literal = sample_value.ToSQLString();

	auto sql = "SELECT " + q_id + ", " + q_seq + ", CAST(" + q_count + " AS BIGINT) FROM " + q_src + " WHERE CAST(" +
	           q_col + " AS VARCHAR) = CAST(" + sample_literal + " AS VARCHAR) ORDER BY " + q_count + " DESC";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("deblur: failed to read table '%s' for sample %s: %s", data.input_table,
		                            sample_literal, result->GetError());
	}

	std::vector<miint::DeblurSequence> sequences;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto read_id_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			auto abundance_val = chunk->GetValue(2, i);
			if (read_id_val.IsNull() || seq_val.IsNull() || abundance_val.IsNull()) {
				continue;
			}
			auto seq_str = seq_val.GetValue<std::string>();
			if (seq_str.empty()) {
				continue;
			}
			sequences.push_back({read_id_val.GetValue<std::string>(), std::move(seq_str),
			                     static_cast<double>(abundance_val.GetValue<int64_t>())});
		}
	}

	if (sequences.empty()) {
		return {};
	}
	try {
		return miint::deblur(std::move(sequences), data.params);
	} catch (const std::invalid_argument &e) {
		throw InvalidInputException("deblur: %s for sample %s", e.what(), sample_literal);
	}
}

// Emit up to STANDARD_VECTOR_SIZE rows from `source` starting at lstate.current_row.
// When has_sample_id, col 0 is the sample column and is filled via a ConstantVector
// reference to lstate.sample_value (constant across the chunk, zero-copy).
static void EmitRows(const DeblurData &data, DeblurLocalState &lstate, const std::vector<miint::DeblurResult> &source,
                     DataChunk &output) {
	idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, source.size() - lstate.current_row);

	idx_t col = 0;
	if (data.has_sample_id) {
		output.data[col++].Reference(lstate.sample_value);
	}
	auto &read_id_vec = output.data[col++];
	auto &seq_vec = output.data[col++];
	auto &abundance_vec = output.data[col++];
	auto read_id_data = FlatVector::GetData<string_t>(read_id_vec);
	auto seq_data = FlatVector::GetData<string_t>(seq_vec);
	auto abundance_data = FlatVector::GetData<int64_t>(abundance_vec);
	for (idx_t i = 0; i < count; i++) {
		auto &r = source[lstate.current_row + i];
		read_id_data[i] = StringVector::AddString(read_id_vec, r.label);
		seq_data[i] = StringVector::AddString(seq_vec, r.sequence);
		abundance_data[i] = r.abundance;
	}
	lstate.current_row += count;
	output.SetCardinality(count);
}

static void DeblurExecute(ClientContext & /*context*/, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<DeblurData>();
	auto &gstate = data_p.global_state->Cast<DeblurGlobalState>();
	auto &lstate = data_p.local_state->Cast<DeblurLocalState>();

	if (!data.has_sample_id) {
		// Single thread (MaxThreads()=1) drains gstate.shared_results via lstate's cursor.
		if (lstate.current_row >= gstate.shared_results.size()) {
			output.SetCardinality(0);
			return;
		}
		EmitRows(data, lstate, gstate.shared_results, output);
		return;
	}

	while (true) {
		if (lstate.current_row < lstate.sample_results.size()) {
			EmitRows(data, lstate, lstate.sample_results, output);
			return;
		}
		// Buffer drained; claim the next sample. Samples that reduce to zero rows
		// fall through to claim again.
		idx_t sample_idx;
		if (!ClaimNextSample(gstate, data.sample_info.sample_values.size(), sample_idx)) {
			output.SetCardinality(0);
			return;
		}
		lstate.sample_value = data.sample_info.sample_values[sample_idx];
		lstate.sample_results = RunDeblurForSample(*lstate.conn, data, lstate.sample_value);
		lstate.current_row = 0;
	}
}

TableFunction DeblurTableFunction::GetFunction() {
	auto tf =
	    TableFunction("deblur", {LogicalType::VARCHAR}, DeblurExecute, DeblurBind, DeblurInitGlobal, DeblurInitLocal);
	tf.named_parameters["mean_error"] = LogicalType::DOUBLE;
	tf.named_parameters["error_profile"] = LogicalType::LIST(LogicalType::DOUBLE);
	tf.named_parameters["indel_prob"] = LogicalType::DOUBLE;
	tf.named_parameters["indel_max"] = LogicalType::INTEGER;
	tf.named_parameters["sample_id"] = LogicalType::VARCHAR;
	tf.named_parameters["id_col"] = LogicalType::VARCHAR;
	tf.named_parameters["sequence_col"] = LogicalType::VARCHAR;
	tf.named_parameters["count_col"] = LogicalType::VARCHAR;
	tf.order_preservation_type = OrderPreservationType::NO_ORDER;
	return tf;
}

void DeblurTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
