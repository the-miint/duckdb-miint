#include "uchime_denovo.hpp"
#include "table_function_common.hpp"
#include "uchime_common.hpp"

#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <algorithm>

namespace duckdb {

// Validate that a table has read_id (VARCHAR), sequence1 (VARCHAR), and size (integer type).
// Single catalog lookup — checks all three columns in one pass.
static void ValidateDenovoTableSchema(ClientContext &context, const std::string &table_name) {
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
	} else {
		throw BinderException("'%s' is not a table or view", table_name);
	}

	// Build case-insensitive lookup
	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	// Check read_id
	auto it = name_to_idx.find("read_id");
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column 'read_id' (VARCHAR)", table_name);
	}
	if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'read_id' in table '%s' must be VARCHAR", table_name);
	}

	// Check sequence1
	it = name_to_idx.find("sequence1");
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column 'sequence1' (VARCHAR)", table_name);
	}
	if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column 'sequence1' in table '%s' must be VARCHAR", table_name);
	}

	// Check size (any integer type)
	it = name_to_idx.find("size");
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column 'size' (integer type)", table_name);
	}
	auto size_type = col_types[it->second].id();
	if (size_type != LogicalTypeId::INTEGER && size_type != LogicalTypeId::BIGINT &&
	    size_type != LogicalTypeId::SMALLINT && size_type != LogicalTypeId::TINYINT &&
	    size_type != LogicalTypeId::HUGEINT) {
		throw BinderException("Column 'size' in table '%s' must be an integer type (got %s)", table_name,
		                      col_types[it->second].ToString());
	}
}

unique_ptr<FunctionData> UchimeDenovoTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->input_table = input.inputs[0].GetValue<std::string>();

	// Validate table schema: read_id (VARCHAR), sequence1 (VARCHAR), size (integer)
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
	// Sorted by size DESC so the most abundant sequences are processed first.
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

	// De novo mode is inherently sequential: each non-chimera must be added to the
	// reference DB before processing the next (lower abundance) query. We process
	// a small batch per Execute() call to allow DuckDB's cancellation mechanism to
	// work between calls, rather than processing everything in a single blocking call.

	// 1. Drain any buffered results first
	if (gstate.result_offset < gstate.results.size()) {
		idx_t remaining = gstate.results.size() - gstate.result_offset;
		idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
		OutputUchimeResults(output, gstate.results, gstate.result_offset, count);
		gstate.result_offset += count;
		return;
	}

	// 2. Process next batch of input sequences
	gstate.results.clear();
	gstate.result_offset = 0;

	idx_t batch_size = STANDARD_VECTOR_SIZE;
	idx_t processed = 0;

	while (gstate.input_offset < gstate.labels.size() && processed < batch_size) {
		idx_t i = gstate.input_offset++;
		miint::UchimeResult result;

		if (gstate.detector.ref_sequences().size() < 2) {
			// Bootstrap: the first two sequences (highest abundance) cannot be chimeras
			// because no more-abundant parents exist. They are unconditionally added to
			// the reference DB to seed the k-mer index for subsequent queries.
			result.query_label = gstate.labels[i];
			result.parent_a_label = "*";
			result.parent_b_label = "*";
			result.closest_parent_label = "*";
			result.flag = "N";
		} else {
			result =
			    gstate.detector.detect_denovo(gstate.labels[i], gstate.sequences[i], gstate.sizes[i], gstate.aligner);
		}

		// Non-chimeras (and borderline ?) are added to the reference DB
		if (result.flag != "Y") {
			gstate.detector.add_to_reference(gstate.labels[i], gstate.sequences[i], gstate.sizes[i]);
		}

		gstate.results.push_back(std::move(result));
		processed++;
	}

	if (gstate.results.empty()) {
		output.SetCardinality(0);
		return;
	}

	// Output the first chunk of results
	idx_t count = std::min(static_cast<idx_t>(gstate.results.size()), static_cast<idx_t>(STANDARD_VECTOR_SIZE));
	OutputUchimeResults(output, gstate.results, 0, count);
	gstate.result_offset = count;
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
