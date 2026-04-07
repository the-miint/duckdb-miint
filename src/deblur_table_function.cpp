#include "deblur_table_function.hpp"
#include "catalog_utils.hpp"
#include "deblur.hpp"

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
};

struct DeblurGlobalState : public GlobalTableFunctionState {
	std::vector<miint::DeblurResult> results;
	idx_t current_row = 0;

	idx_t MaxThreads() const override {
		return 1;
	}
};

// Validate that the table/view has read_id (VARCHAR), sequence1 (VARCHAR),
// and abundance (any integer type).
static void ValidateDeblurTableSchema(ClientContext &context, const std::string &table_name) {
	auto cols = GetTableOrViewColumns(context, table_name, "Deblur input");

	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < cols.names.size(); i++) {
		name_to_idx[StringUtil::Lower(cols.names[i])] = i;
	}

	auto it = name_to_idx.find("read_id");
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column 'read_id' (VARCHAR)", table_name);
	}
	if (cols.types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("deblur: column 'read_id' in table '%s' must be VARCHAR", table_name);
	}

	it = name_to_idx.find("sequence1");
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column 'sequence1' (VARCHAR)", table_name);
	}
	if (cols.types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("deblur: column 'sequence1' in table '%s' must be VARCHAR", table_name);
	}

	it = name_to_idx.find("abundance");
	if (it == name_to_idx.end()) {
		throw BinderException("deblur: table '%s' missing required column 'abundance' (integer type)", table_name);
	}
	auto atype = cols.types[it->second].id();
	if (atype != LogicalTypeId::INTEGER && atype != LogicalTypeId::BIGINT && atype != LogicalTypeId::SMALLINT &&
	    atype != LogicalTypeId::TINYINT && atype != LogicalTypeId::HUGEINT && atype != LogicalTypeId::UINTEGER &&
	    atype != LogicalTypeId::UBIGINT && atype != LogicalTypeId::USMALLINT && atype != LogicalTypeId::UTINYINT &&
	    atype != LogicalTypeId::UHUGEINT) {
		throw BinderException("deblur: column 'abundance' in table '%s' must be an integer type (got %s)", table_name,
		                      cols.types[it->second].ToString());
	}
}

static unique_ptr<FunctionData> DeblurBind(ClientContext &context, TableFunctionBindInput &input,
                                           vector<LogicalType> &return_types, vector<std::string> &names) {
	auto data = make_uniq<DeblurData>();
	data->input_table = input.inputs[0].GetValue<std::string>();

	ValidateDeblurTableSchema(context, data->input_table);

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

	names = {"read_id", "sequence", "abundance"};
	return_types = {LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::BIGINT};

	return data;
}

static unique_ptr<GlobalTableFunctionState> DeblurInitGlobal(ClientContext &context, TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<DeblurData>();
	auto gstate = make_uniq<DeblurGlobalState>();

	// Load all sequences via a separate connection to avoid deadlocking
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.Query("SELECT \"read_id\", \"sequence1\", CAST(\"abundance\" AS BIGINT) FROM " +
	                         KeywordHelper::WriteOptionallyQuoted(data.input_table) + " ORDER BY \"abundance\" DESC");
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
			gstate->results = miint::deblur(std::move(sequences), data.params);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("deblur: %s", e.what());
		}
	}

	return gstate;
}

static void DeblurExecute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<DeblurGlobalState>();

	idx_t total = gstate.results.size();
	if (gstate.current_row >= total) {
		output.SetCardinality(0);
		return;
	}

	idx_t count = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.current_row);

	for (idx_t i = 0; i < count; i++) {
		idx_t row = gstate.current_row + i;
		auto &r = gstate.results[row];

		FlatVector::GetData<string_t>(output.data[0])[i] = StringVector::AddString(output.data[0], r.label);
		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], r.sequence);
		FlatVector::GetData<int64_t>(output.data[2])[i] = r.abundance;
	}

	gstate.current_row += count;
	output.SetCardinality(count);
}

TableFunction DeblurTableFunction::GetFunction() {
	auto tf = TableFunction("deblur", {LogicalType::VARCHAR}, DeblurExecute, DeblurBind, DeblurInitGlobal);
	tf.named_parameters["mean_error"] = LogicalType::DOUBLE;
	tf.named_parameters["error_profile"] = LogicalType::LIST(LogicalType::DOUBLE);
	tf.named_parameters["indel_prob"] = LogicalType::DOUBLE;
	tf.named_parameters["indel_max"] = LogicalType::INTEGER;
	return tf;
}

void DeblurTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
