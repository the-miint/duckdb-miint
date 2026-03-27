#include "uchime_denovo.hpp"
#include "table_function_common.hpp"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>
#include <numeric>

namespace duckdb {

// Reuse the same output schema as uchime_ref (defined in uchime_ref.cpp).
// Duplicated here to avoid cross-file dependency — both produce identical 18-column output.
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

// Same OutputUchimeResults as uchime_ref — duplicated to avoid cross-file linkage.
static idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset,
                                 idx_t count) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = 0;

	auto score_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		score_data[i] = results[offset + i].score;
	}

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

	auto div_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		div_data[i] = results[offset + i].divergence;
	}

	auto &flag_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(flag_vec)[i] = StringVector::AddString(flag_vec, results[offset + i].flag);
	}

	D_ASSERT(col == output.ColumnCount());
	output.SetCardinality(actual);
	return actual;
}

// Validate that a table has read_id (VARCHAR), sequence1 (VARCHAR), and size (INTEGER/BIGINT).
static void ValidateDenovoTableSchema(ClientContext &context, const std::string &table_name) {
	// First validate read_id + sequence1 via existing infrastructure
	ValidateSequenceTableSchema(context, table_name);

	// Additionally check for size column
	EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, table_name, QueryErrorContext());
	auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);
	if (!entry) {
		throw BinderException("Table or view '%s' does not exist", table_name);
	}

	vector<string> col_names;
	vector<LogicalType> col_types;
	if (entry->type == CatalogType::TABLE_ENTRY) {
		auto &table = entry->Cast<TableCatalogEntry>();
		auto &columns = table.GetColumns();
		for (idx_t i = 0; i < columns.LogicalColumnCount(); i++) {
			auto &c = columns.GetColumn(LogicalIndex(i));
			col_names.push_back(c.Name());
			col_types.push_back(c.Type());
		}
	} else if (entry->type == CatalogType::VIEW_ENTRY) {
		auto &view = entry->Cast<ViewCatalogEntry>();
		view.BindView(context);
		auto col_info = view.GetColumnInfo();
		col_names = col_info->names;
		col_types = col_info->types;
	}

	bool found_size = false;
	for (idx_t i = 0; i < col_names.size(); i++) {
		if (StringUtil::Lower(col_names[i]) == "size") {
			auto tid = col_types[i].id();
			if (tid != LogicalTypeId::INTEGER && tid != LogicalTypeId::BIGINT && tid != LogicalTypeId::SMALLINT &&
			    tid != LogicalTypeId::TINYINT && tid != LogicalTypeId::HUGEINT) {
				throw BinderException("Column 'size' in table '%s' must be an integer type", table_name);
			}
			found_size = true;
			break;
		}
	}
	if (!found_size) {
		throw BinderException("Table '%s' missing required column 'size' (INTEGER)", table_name);
	}
}

unique_ptr<FunctionData> UchimeDenovoTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->input_table = input.inputs[0].GetValue<std::string>();

	// Validate table schema: read_id (VARCHAR), sequence1 (VARCHAR), size (INTEGER)
	ValidateDenovoTableSchema(context, data->input_table);

	// Optional scoring parameters with range validation
	auto get_double = [&](const std::string &name, double &out, double min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<double>();
			if (out < min_val) {
				throw InvalidInputException("uchime_denovo: %s must be %s (got %g)", name, constraint, out);
			}
		}
	};
	auto get_int = [&](const std::string &name, int &out, int min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val) {
				throw InvalidInputException("uchime_denovo: %s must be %s (got %d)", name, constraint, out);
			}
		}
	};

	get_double("minh", data->params.minh, 0.0, ">= 0");
	get_double("xn", data->params.xn, 1.0, ">= 1.0");
	get_double("dn", data->params.dn, 0.0, ">= 0");
	get_double("mindiv", data->params.mindiv, 0.0, ">= 0");
	get_double("abskew", data->params.abskew, 1.0, ">= 1.0");
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

unique_ptr<GlobalTableFunctionState> UchimeDenovoTableFunction::InitGlobal(ClientContext &context,
                                                                           TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = make_uniq<GlobalState>();

	gstate->detector = miint::ChimeraDetector(data.params);

	// Load all sequences with their abundance via a separate connection.
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	auto result = conn.Query("SELECT read_id, sequence1, size FROM " +
	                         KeywordHelper::WriteOptionallyQuoted(data.input_table) + " ORDER BY size DESC");
	if (result->HasError()) {
		throw InvalidInputException("Failed to read table '%s': %s", data.input_table, result->GetError());
	}

	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		for (idx_t i = 0; i < chunk->size(); i++) {
			auto read_id_val = chunk->GetValue(0, i);
			auto seq_val = chunk->GetValue(1, i);
			auto size_val = chunk->GetValue(2, i);
			if (read_id_val.IsNull() || seq_val.IsNull() || size_val.IsNull()) {
				continue;
			}
			auto seq_str = seq_val.GetValue<std::string>();
			if (seq_str.empty()) {
				continue;
			}
			gstate->labels.push_back(read_id_val.GetValue<std::string>());
			gstate->sequences.push_back(std::move(seq_str));
			gstate->sizes.push_back(size_val.GetValue<int64_t>());
		}
	}

	if (gstate->labels.empty()) {
		throw InvalidInputException("Table '%s' is empty (or contains only NULL/empty sequences)", data.input_table);
	}

	return gstate;
}

void UchimeDenovoTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<GlobalState>();

	// Process all sequences on first Execute call (single-threaded, sequential).
	// De novo mode is inherently sequential: each non-chimera is added to the
	// reference DB before processing the next query (lower abundance).
	if (!gstate.processed) {
		gstate.processed = true;
		gstate.results.reserve(gstate.labels.size());

		for (idx_t i = 0; i < gstate.labels.size(); i++) {
			miint::UchimeResult result;

			if (gstate.detector.ref_sequences().size() < 2) {
				// Not enough references yet — cannot detect chimeras.
				// Add this sequence to the reference and report as non-chimeric.
				result.query_label = gstate.labels[i];
				result.parent_a_label = "*";
				result.parent_b_label = "*";
				result.closest_parent_label = "*";
				result.flag = "N";
			} else {
				result = gstate.detector.detect_denovo(gstate.labels[i], gstate.sequences[i], gstate.sizes[i],
				                                       gstate.aligner);
			}

			// De novo: non-chimeras (and borderline) are added to the reference DB
			if (result.flag != "Y") {
				gstate.detector.add_to_reference(gstate.labels[i], gstate.sequences[i], gstate.sizes[i]);
			}

			gstate.results.push_back(std::move(result));
		}
	}

	// Drain results in STANDARD_VECTOR_SIZE chunks
	if (gstate.result_offset < gstate.results.size()) {
		idx_t remaining = gstate.results.size() - gstate.result_offset;
		idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
		OutputUchimeResults(output, gstate.results, gstate.result_offset, count);
		gstate.result_offset += count;
	} else {
		output.SetCardinality(0);
	}
}

TableFunction UchimeDenovoTableFunction::GetFunction() {
	auto tf = TableFunction("uchime_denovo", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);

	tf.named_parameters["minh"] = LogicalType::DOUBLE;
	tf.named_parameters["xn"] = LogicalType::DOUBLE;
	tf.named_parameters["dn"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiv"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiffs"] = LogicalType::INTEGER;
	tf.named_parameters["abskew"] = LogicalType::DOUBLE;

	return tf;
}

void UchimeDenovoTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
