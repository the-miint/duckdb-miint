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

// Validate that a table has the resolved id (VARCHAR), sequence (VARCHAR), and count
// (integer type) columns. Names are resolved from named parameters at bind time.
static void ValidateDenovoTableSchema(ClientContext &context, const std::string &table_name, const std::string &id_col,
                                      const std::string &sequence_col, const std::string &count_col) {
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

	std::unordered_map<string, idx_t> name_to_idx;
	for (idx_t i = 0; i < col_names.size(); i++) {
		name_to_idx[StringUtil::Lower(col_names[i])] = i;
	}

	auto it = name_to_idx.find(StringUtil::Lower(id_col));
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column '%s' (VARCHAR)", table_name, id_col);
	}
	if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column '%s' in table '%s' must be VARCHAR", id_col, table_name);
	}

	it = name_to_idx.find(StringUtil::Lower(sequence_col));
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column '%s' (VARCHAR)", table_name, sequence_col);
	}
	if (col_types[it->second].id() != LogicalTypeId::VARCHAR) {
		throw BinderException("Column '%s' in table '%s' must be VARCHAR", sequence_col, table_name);
	}

	it = name_to_idx.find(StringUtil::Lower(count_col));
	if (it == name_to_idx.end()) {
		throw BinderException("Table '%s' missing required column '%s' (integer type)", table_name, count_col);
	}
	auto size_type = col_types[it->second].id();
	if (size_type != LogicalTypeId::INTEGER && size_type != LogicalTypeId::BIGINT &&
	    size_type != LogicalTypeId::SMALLINT && size_type != LogicalTypeId::TINYINT &&
	    size_type != LogicalTypeId::HUGEINT && size_type != LogicalTypeId::UINTEGER &&
	    size_type != LogicalTypeId::UBIGINT && size_type != LogicalTypeId::USMALLINT &&
	    size_type != LogicalTypeId::UTINYINT && size_type != LogicalTypeId::UHUGEINT) {
		throw BinderException("Column '%s' in table '%s' must be an integer type (got %s)", count_col, table_name,
		                      col_types[it->second].ToString());
	}
}

// Pull (id, sequence, count) rows for a single sample (or the whole table when
// where_sql is empty) via the given connection, ordered by count DESC, dropping
// NULL/empty rows. Column names are caller-supplied to honor user overrides.
static void LoadDenovoSequences(Connection &conn, const std::string &table_name, const std::string &id_col,
                                const std::string &sequence_col, const std::string &count_col,
                                const std::string &where_sql, std::vector<std::string> &out_labels,
                                std::vector<std::string> &out_sequences, std::vector<int64_t> &out_sizes) {
	auto q_id = KeywordHelper::WriteOptionallyQuoted(id_col);
	auto q_seq = KeywordHelper::WriteOptionallyQuoted(sequence_col);
	auto q_count = KeywordHelper::WriteOptionallyQuoted(count_col);
	auto sql =
	    "SELECT " + q_id + ", " + q_seq + ", " + q_count + " FROM " + KeywordHelper::WriteOptionallyQuoted(table_name);
	if (!where_sql.empty()) {
		sql += " WHERE " + where_sql;
	}
	sql += " ORDER BY " + q_count + " DESC";
	auto result = conn.Query(sql);
	if (result->HasError()) {
		throw InvalidInputException("Failed to read table '%s': %s", table_name, result->GetError());
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
			out_labels.push_back(read_id_val.GetValue<std::string>());
			out_sequences.push_back(std::move(seq_str));
			out_sizes.push_back(size_val.GetValue<int64_t>());
		}
	}
}

// Run the denovo bootstrap + incremental k-mer indexing over a pre-loaded set of
// size-sorted sequences, returning the full result vector. Used by the per-sample
// path to materialize one sample's results before emission.
static std::vector<miint::UchimeResult> RunDenovoForSet(const miint::UchimeParams &params,
                                                        const std::vector<std::string> &labels,
                                                        const std::vector<std::string> &sequences,
                                                        const std::vector<int64_t> &sizes) {
	miint::VsearchChimeraWrapper wrapper(params);
	wrapper.prepare_denovo(labels, sequences, sizes);

	std::vector<miint::UchimeResult> results;
	results.reserve(labels.size());

	idx_t indexed_count = 0;
	for (idx_t i = 0; i < labels.size(); i++) {
		miint::UchimeResult r;
		if (indexed_count < 2) {
			// Bootstrap: seed the k-mer index with the two most abundant sequences.
			r.query_label = labels[i];
			r.flag = "N";
		} else {
			r = wrapper.detect_denovo(labels[i], sequences[i], sizes[i]);
		}
		if (r.flag != "Y") {
			wrapper.index_sequence(i);
			indexed_count++;
		}
		results.push_back(std::move(r));
	}
	return results;
}

unique_ptr<FunctionData> UchimeDenovoTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                         vector<LogicalType> &return_types,
                                                         vector<std::string> &names) {
	auto data = make_uniq<Data>();

	data->input_table = input.inputs[0].GetValue<std::string>();

	auto get_col_override = [&](const std::string &param_name, std::string &out) {
		auto it = input.named_parameters.find(param_name);
		if (it != input.named_parameters.end()) {
			auto val = it->second.GetValue<std::string>();
			if (val.empty()) {
				throw InvalidInputException("detect_chimera_uchime_denovo: %s must not be empty", param_name);
			}
			out = std::move(val);
		}
	};
	get_col_override("id_col", data->id_col);
	get_col_override("sequence_col", data->sequence_col);
	get_col_override("count_col", data->count_col);

	ValidateDenovoTableSchema(context, data->input_table, data->id_col, data->sequence_col, data->count_col);

	auto get_double = [&](const std::string &name, double &out, double min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<double>();
			if (out < min_val) {
				throw InvalidInputException("detect_chimera_uchime_denovo: %s must be %s (got %g)", name, constraint,
				                            out);
			}
		}
	};
	auto get_int = [&](const std::string &name, int &out, int min_val, const char *constraint) {
		auto it = input.named_parameters.find(name);
		if (it != input.named_parameters.end()) {
			out = it->second.GetValue<int>();
			if (out < min_val) {
				throw InvalidInputException("detect_chimera_uchime_denovo: %s must be %s (got %d)", name, constraint,
				                            out);
			}
		}
	};

	get_double("minh", data->params.minh, 0.0, ">= 0");
	get_double("xn", data->params.xn, 1.0, ">= 1.0");
	get_double("dn", data->params.dn, 0.0, ">= 0");
	get_double("mindiv", data->params.mindiv, 0.0, ">= 0");
	get_double("abskew", data->params.abskew, 1.0, ">= 1.0");
	get_int("mindiffs", data->params.mindiffs, 1, ">= 1");

	auto sample_it = input.named_parameters.find("sample_id");
	if (sample_it != input.named_parameters.end()) {
		data->has_sample_id = true;
		data->sample_info.sample_id_col = sample_it->second.GetValue<string>();
	}

	data->names = GetUchimeOutputNames();
	data->types = GetUchimeOutputTypes();

	if (data->has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context);
		Connection conn(db);
		// Reserved = 18 uchime output columns (case-insensitive, lowercase already).
		DiscoverSamples(conn, data->input_table, data->sample_info.sample_id_col, data->names,
		                "detect_chimera_uchime_denovo", data->sample_info);

		names.push_back(data->sample_info.sample_id_col);
		return_types.push_back(data->sample_info.sample_id_type);
	}

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

	if (data.has_sample_id) {
		// Wrapper is explicitly not thread-safe per its header; clamp to 1.
		InitPerSampleGlobal(context, *gstate, data.sample_info.sample_values.size(), /*max_threads_hint=*/1);
		return gstate;
	}

	gstate->max_threads = 1;
	gstate->wrapper = miint::VsearchChimeraWrapper(data.params);

	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	LoadDenovoSequences(conn, data.input_table, data.id_col, data.sequence_col, data.count_col, /*where_sql=*/"",
	                    gstate->labels, gstate->sequences, gstate->sizes);

	if (gstate->labels.empty()) {
		throw InvalidInputException("Table '%s' is empty (or contains only NULL/empty sequences)", data.input_table);
	}

	// Load all sequences into vsearch's DB with DUST masking. K-mer index is
	// allocated but empty — non-chimeras are indexed incrementally in Execute().
	gstate->wrapper.prepare_denovo(gstate->labels, gstate->sequences, gstate->sizes);

	return gstate;
}

unique_ptr<LocalTableFunctionState> UchimeDenovoTableFunction::InitLocal(ExecutionContext &context,
                                                                         TableFunctionInitInput &input,
                                                                         GlobalTableFunctionState * /*gstate*/) {
	auto &data = input.bind_data->Cast<Data>();
	auto lstate = make_uniq<LocalState>();
	if (data.has_sample_id) {
		auto &db = DatabaseInstance::GetDatabase(context.client);
		lstate->conn = make_uniq<Connection>(db);
	}
	return lstate;
}

void UchimeDenovoTableFunction::Execute(ClientContext & /*context*/, TableFunctionInput &data_p, DataChunk &output) {
	auto &data = data_p.bind_data->Cast<Data>();
	auto &gstate = data_p.global_state->Cast<GlobalState>();
	auto &lstate = data_p.local_state->Cast<LocalState>();

	if (!data.has_sample_id) {
		// Non-sample path: incremental processing preserves the existing cancellation
		// window (one STANDARD_VECTOR_SIZE batch of queries per Execute call).

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

			if (gstate.indexed_count < 2) {
				result.query_label = gstate.labels[i];
				result.flag = "N";
			} else {
				result = gstate.wrapper.detect_denovo(gstate.labels[i], gstate.sequences[i], gstate.sizes[i]);
			}

			if (result.flag != "Y") {
				gstate.wrapper.index_sequence(i);
				gstate.indexed_count++;
			}

			gstate.results.push_back(std::move(result));
			processed++;
		}

		if (gstate.results.empty()) {
			output.SetCardinality(0);
			return;
		}

		idx_t count = std::min(static_cast<idx_t>(gstate.results.size()), static_cast<idx_t>(STANDARD_VECTOR_SIZE));
		OutputUchimeResults(output, gstate.results, 0, count);
		gstate.result_offset = count;
		return;
	}

	// Sample path: run the full denovo compute for one sample, emit it in chunks,
	// then claim the next. Samples that reduce to zero rows fall through the loop.
	while (true) {
		if (lstate.result_offset < lstate.results.size()) {
			idx_t remaining = lstate.results.size() - lstate.result_offset;
			idx_t count = std::min(remaining, static_cast<idx_t>(STANDARD_VECTOR_SIZE));
			output.data[0].Reference(lstate.sample_value);
			OutputUchimeResults(output, lstate.results, lstate.result_offset, count, /*start_col=*/1);
			lstate.result_offset += count;
			return;
		}

		idx_t sample_idx;
		if (!ClaimNextSample(gstate, data.sample_info.sample_values.size(), sample_idx)) {
			output.SetCardinality(0);
			return;
		}
		lstate.sample_value = data.sample_info.sample_values[sample_idx];
		auto sample_literal = lstate.sample_value.ToSQLString();
		auto q_col = KeywordHelper::WriteOptionallyQuoted(data.sample_info.sample_id_col);
		// CAST-as-VARCHAR equality: see note in deblur_table_function.cpp re: DECIMAL.
		auto where_sql = "CAST(" + q_col + " AS VARCHAR) = CAST(" + sample_literal + " AS VARCHAR)";

		std::vector<std::string> labels, sequences;
		std::vector<int64_t> sizes;
		LoadDenovoSequences(*lstate.conn, data.input_table, data.id_col, data.sequence_col, data.count_col, where_sql,
		                    labels, sequences, sizes);
		lstate.results.clear();
		lstate.result_offset = 0;
		if (!labels.empty()) {
			lstate.results = RunDenovoForSet(data.params, labels, sequences, sizes);
		}
	}
}

TableFunction UchimeDenovoTableFunction::GetFunction() {
	auto tf =
	    TableFunction("detect_chimera_uchime_denovo", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);

	tf.named_parameters["minh"] = LogicalType::DOUBLE;
	tf.named_parameters["xn"] = LogicalType::DOUBLE;
	tf.named_parameters["dn"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiv"] = LogicalType::DOUBLE;
	tf.named_parameters["mindiffs"] = LogicalType::INTEGER;
	tf.named_parameters["abskew"] = LogicalType::DOUBLE;
	tf.named_parameters["sample_id"] = LogicalType::VARCHAR;
	tf.named_parameters["id_col"] = LogicalType::VARCHAR;
	tf.named_parameters["sequence_col"] = LogicalType::VARCHAR;
	tf.named_parameters["count_col"] = LogicalType::VARCHAR;

	tf.order_preservation_type = OrderPreservationType::NO_ORDER;

	return tf;
}

void UchimeDenovoTableFunction::Register(ExtensionLoader &loader) {
	auto tf = GetFunction();
	loader.RegisterFunction(tf);
}

} // namespace duckdb
