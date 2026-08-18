#include "read_ena.hpp"
#include "duckdb/common/vector_size.hpp"
#include <algorithm>
#include <sstream>

namespace duckdb {

// ---- Helpers ----

static std::string TrimWhitespace(const std::string &s) {
	auto start = s.find_first_not_of(" \t\n\r");
	if (start == std::string::npos) {
		return "";
	}
	auto end = s.find_last_not_of(" \t\n\r");
	return s.substr(start, end - start + 1);
}

// Find a column value by name in a parsed TSV result.
static std::string FindColumnValue(const miint::ENATSVResult &parsed, const std::string &column_name) {
	for (size_t i = 0; i < parsed.column_names.size(); i++) {
		if (parsed.column_names[i] == column_name && !parsed.rows.empty() && i < parsed.rows[0].size()) {
			return parsed.rows[0][i];
		}
	}
	return "";
}

// Resolve an accession to the correct query accession for a given result type.
static std::string ResolveAccession(miint::ENAClient &client, const std::string &accession,
                                    const std::string &result_type) {
	auto acc_type = miint::ENAParser::DetectAccessionType(accession);

	if (result_type == "read_run") {
		return accession;
	}

	if (result_type == "sample") {
		if (acc_type == miint::ENAAccessionType::SAMPLE || acc_type == miint::ENAAccessionType::STUDY) {
			return accession;
		}
		auto tsv = client.Search(accession, "read_run", "sample_accession");
		auto parsed = miint::ENAParser::ParseTSV(tsv);
		auto resolved = FindColumnValue(parsed, "sample_accession");
		if (!resolved.empty()) {
			return resolved;
		}
		throw IOException("read_ena: could not resolve '%s' to a sample accession", accession);
	}

	if (result_type == "study") {
		if (acc_type == miint::ENAAccessionType::STUDY) {
			return accession;
		}
		auto tsv = client.Search(accession, "read_run", "study_accession");
		auto parsed = miint::ENAParser::ParseTSV(tsv);
		auto resolved = FindColumnValue(parsed, "study_accession");
		if (!resolved.empty()) {
			return resolved;
		}
		throw IOException("read_ena: could not resolve '%s' to a study accession", accession);
	}

	return accession;
}

// Reorder a TSV row to match the schema column order.
// ENA returns columns in its own order; we need them in the user-requested order.
static std::vector<std::string> ReorderRow(const std::vector<std::string> &row,
                                           const std::vector<std::string> &tsv_column_names,
                                           const std::vector<int> &column_map) {
	std::vector<std::string> reordered(column_map.size());
	for (size_t i = 0; i < column_map.size(); i++) {
		int src_idx = column_map[i];
		if (src_idx >= 0 && src_idx < static_cast<int>(row.size())) {
			reordered[i] = row[src_idx];
		}
	}
	return reordered;
}

// Build a mapping from schema column positions to TSV column positions.
static std::vector<int> BuildColumnMap(const std::vector<std::string> &schema_names,
                                       const std::vector<std::string> &tsv_column_names) {
	std::vector<int> map;
	for (const auto &name : schema_names) {
		int found = -1;
		for (size_t i = 0; i < tsv_column_names.size(); i++) {
			if (tsv_column_names[i] == name) {
				found = static_cast<int>(i);
				break;
			}
		}
		map.push_back(found);
	}
	return map;
}

// Map a field name to its LogicalType. Fields not listed here default to VARCHAR.
static LogicalType FieldTypeForColumn(const std::string &field_name) {
	if (field_name == "read_count" || field_name == "base_count" || field_name == "tax_id") {
		return LogicalType::BIGINT;
	}
	if (field_name == "fastq_ftp" || field_name == "fastq_aspera" || field_name == "fastq_md5") {
		return LogicalType::LIST(LogicalType::VARCHAR);
	}
	if (field_name == "fastq_bytes") {
		return LogicalType::LIST(LogicalType::BIGINT);
	}
	return LogicalType::VARCHAR;
}

// Split a ;-delimited string into parts (mirrors ENAParser::SplitSemicolon semantics).
static std::vector<std::string> SplitSemicolon(const std::string &field) {
	std::vector<std::string> parts;
	if (field.empty()) {
		return parts;
	}
	std::string::size_type start = 0;
	while (true) {
		auto pos = field.find(';', start);
		if (pos == std::string::npos) {
			parts.push_back(field.substr(start));
			break;
		}
		parts.push_back(field.substr(start, pos - start));
		start = pos + 1;
	}
	return parts;
}

// ---- Data ----

ReadENATableFunction::Data::Data(std::vector<std::string> accessions, const std::string &result_type,
                                 const std::string &fields)
    : accessions(std::move(accessions)), result_type(result_type), fields(fields) {
	// Column names and types from the user-specified fields string (comma-separated)
	std::istringstream stream(fields);
	std::string field;
	while (std::getline(stream, field, ',')) {
		names.push_back(field);
		types.push_back(FieldTypeForColumn(field));
	}
}

// ---- GlobalState ----

ReadENATableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::vector<std::string> &accessions,
                                               const std::string &result_type, const std::string &fields,
                                               const std::vector<std::string> &schema_names)
    : client(make_uniq<miint::ENAClient>(db)), next_accession_idx(0), row_offset(0), accessions(accessions),
      result_type(result_type), fields(fields), schema_names(schema_names) {
}

bool ReadENATableFunction::GlobalState::FetchNextAccession() {
	if (next_accession_idx >= accessions.size()) {
		return false;
	}

	const auto &accession = accessions[next_accession_idx];
	next_accession_idx++;

	auto resolved = ResolveAccession(*client, accession, result_type);
	auto tsv_body = client->Search(resolved, result_type, fields);
	auto parsed = miint::ENAParser::ParseTSV(tsv_body);

	// Clear previously consumed rows to avoid unbounded memory growth
	rows.clear();
	row_offset = 0;

	// Reorder each row to match the schema column order
	auto column_map = BuildColumnMap(schema_names, parsed.column_names);
	for (auto &row : parsed.rows) {
		rows.push_back(ReorderRow(row, parsed.column_names, column_map));
	}

	return true;
}

// ---- Bind ----

unique_ptr<FunctionData> ReadENATableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                    vector<LogicalType> &return_types, vector<std::string> &names) {
	std::vector<std::string> accessions;
	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_ena: accession cannot be NULL");
	}
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_ena: accession list cannot contain NULL");
			}
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ena: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena: at least one accession must be provided");
	}
	for (auto &acc : accessions) {
		acc = TrimWhitespace(acc);
		if (acc.empty()) {
			throw InvalidInputException("read_ena: accession cannot be empty");
		}
	}

	// Parse result type (default: read_run)
	std::string result_type = "read_run";
	auto result_param = input.named_parameters.find("result");
	if (result_param != input.named_parameters.end() && !result_param->second.IsNull()) {
		result_type = result_param->second.ToString();
		if (result_type != "read_run" && result_type != "sample" && result_type != "study") {
			throw InvalidInputException("read_ena: result must be one of: read_run, sample, study");
		}
	}

	// Parse fields (default: from ENAParser::DefaultFields)
	std::string fields;
	auto fields_param = input.named_parameters.find("fields");
	if (fields_param != input.named_parameters.end() && !fields_param->second.IsNull()) {
		fields = fields_param->second.ToString();
	} else {
		fields = miint::ENAParser::DefaultFields(result_type);
	}

	if (fields.empty()) {
		throw InvalidInputException("read_ena: fields parameter cannot be empty");
	}

	auto data = make_uniq<Data>(std::move(accessions), result_type, fields);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENATableFunction::InitGlobal(ClientContext &context,
                                                                      TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.accessions, data.result_type, data.fields, data.names);
}

unique_ptr<LocalTableFunctionState> ReadENATableFunction::InitLocal(ExecutionContext &context,
                                                                    TableFunctionInitInput &input,
                                                                    GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENATableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	while (global_state.row_offset >= global_state.rows.size()) {
		if (!global_state.FetchNextAccession()) {
			output.SetCardinality(0);
			return;
		}
	}

	idx_t remaining = global_state.rows.size() - global_state.row_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	auto &rows = global_state.rows;
	size_t offset = global_state.row_offset;
	idx_t num_cols = bind_data.names.size();

	for (idx_t col = 0; col < num_cols; col++) {
		const auto &col_type = bind_data.types[col];
		auto &validity = FlatVector::Validity(output.data[col]);

		if (col_type.id() == LogicalTypeId::BIGINT) {
			auto col_data = FlatVector::GetData<int64_t>(output.data[col]);
			for (idx_t i = 0; i < count; i++) {
				auto &row = rows[offset + i];
				if (col < row.size()) {
					const auto &val = row[col];
					if (val.empty()) {
						validity.SetInvalid(i);
					} else {
						try {
							col_data[i] = std::stoll(val);
						} catch (const std::exception &) {
							throw IOException("read_ena: cannot parse '%s' as BIGINT for column '%s'", val,
							                  bind_data.names[col]);
						}
					}
				} else {
					validity.SetInvalid(i);
				}
			}
		} else if (col_type.id() == LogicalTypeId::LIST) {
			// List columns: LIST(VARCHAR) or LIST(BIGINT)
			auto &child_type = ListType::GetChildType(col_type);
			auto list_entries = ListVector::GetData(output.data[col]);
			auto &child_vec = ListVector::GetEntry(output.data[col]);

			// Reserve capacity for the worst case (every row contributes at least one element)
			idx_t worst_case = 0;
			for (idx_t i = 0; i < count; i++) {
				auto &row = rows[offset + i];
				if (col < row.size() && !row[col].empty()) {
					auto parts = SplitSemicolon(row[col]);
					worst_case += parts.size();
				}
			}
			ListVector::Reserve(output.data[col], worst_case);

			idx_t child_offset = 0;
			for (idx_t i = 0; i < count; i++) {
				auto &row = rows[offset + i];
				if (col < row.size() && !row[col].empty()) {
					auto parts = SplitSemicolon(row[col]);
					list_entries[i] = list_entry_t(child_offset, parts.size());

					if (child_type.id() == LogicalTypeId::BIGINT) {
						auto child_data = FlatVector::GetData<int64_t>(child_vec);
						for (size_t j = 0; j < parts.size(); j++) {
							try {
								child_data[child_offset + j] = std::stoll(parts[j]);
							} catch (const std::exception &) {
								throw IOException("read_ena: cannot parse '%s' as BIGINT for column '%s'", parts[j],
								                  bind_data.names[col]);
							}
						}
					} else {
						auto child_data = FlatVector::GetData<string_t>(child_vec);
						for (size_t j = 0; j < parts.size(); j++) {
							child_data[child_offset + j] = StringVector::AddString(child_vec, parts[j]);
						}
					}
					child_offset += parts.size();
				} else {
					// Empty/missing → NULL list
					validity.SetInvalid(i);
					list_entries[i] = list_entry_t(child_offset, 0);
				}
			}
			ListVector::SetListSize(output.data[col], child_offset);
		} else {
			// Default: VARCHAR
			auto col_data = FlatVector::GetData<string_t>(output.data[col]);
			for (idx_t i = 0; i < count; i++) {
				auto &row = rows[offset + i];
				if (col < row.size()) {
					const auto &val = row[col];
					if (val.empty()) {
						validity.SetInvalid(i);
					} else {
						col_data[i] = StringVector::AddString(output.data[col], val);
					}
				} else {
					validity.SetInvalid(i);
				}
			}
		}
	}

	global_state.row_offset += count;
	output.SetCardinality(count);
}

// ---- Registration ----

TableFunction ReadENATableFunction::GetFunction() {
	auto tf = TableFunction("read_ena", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["result"] = LogicalType::VARCHAR;
	tf.named_parameters["fields"] = LogicalType::VARCHAR;
	return tf;
}

void ReadENATableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
