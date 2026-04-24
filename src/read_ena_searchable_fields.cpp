#include "read_ena_searchable_fields.hpp"

#include "duckdb/common/vector_size.hpp"

#include <sstream>

namespace duckdb {

// ---- Data / GlobalState ----

ReadENASearchableFieldsTableFunction::Data::Data(std::string result_type) : result_type(std::move(result_type)) {
}

ReadENASearchableFieldsTableFunction::GlobalState::GlobalState(DatabaseInstance &db)
    : client(make_uniq<miint::ENAClient>(db)) {
}

// ---- Bind ----

unique_ptr<FunctionData> ReadENASearchableFieldsTableFunction::Bind(ClientContext &context,
                                                                    TableFunctionBindInput &input,
                                                                    vector<LogicalType> &return_types,
                                                                    vector<std::string> &names) {
	if (input.inputs.empty() || input.inputs[0].IsNull()) {
		throw InvalidInputException(
		    "ena_searchable_fields: result_type is required (e.g., 'sample', 'read_run', 'study', 'experiment')");
	}
	if (input.inputs[0].type().id() != LogicalTypeId::VARCHAR) {
		throw InvalidInputException("ena_searchable_fields: result_type must be VARCHAR");
	}
	std::string result_type = input.inputs[0].ToString();

	// Trim whitespace and reject empty.
	auto first = result_type.find_first_not_of(" \t\n\r");
	if (first == std::string::npos) {
		throw InvalidInputException("ena_searchable_fields: result_type cannot be empty");
	}
	auto last = result_type.find_last_not_of(" \t\n\r");
	result_type = result_type.substr(first, last - first + 1);

	// Validate: ENA's /returnFields endpoint accepts a result-type identifier.
	// Apply the same charset policy as `ENAParser::ValidateAccession` (alnum +
	// _-.) to prevent URL injection.
	for (char c : result_type) {
		if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_' && c != '-' && c != '.') {
			throw InvalidInputException(
			    "ena_searchable_fields: result_type contains invalid characters (allowed: alphanumeric, _ - .)");
		}
	}

	names.emplace_back("field_name");
	names.emplace_back("type");
	names.emplace_back("description");
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::VARCHAR);
	return_types.emplace_back(LogicalType::VARCHAR);

	return make_uniq<Data>(std::move(result_type));
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENASearchableFieldsTableFunction::InitGlobal(ClientContext &context,
                                                                                      TableFunctionInitInput &input) {
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db);
}

unique_ptr<LocalTableFunctionState>
ReadENASearchableFieldsTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                                GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENASearchableFieldsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                   DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global = data_p.global_state->Cast<GlobalState>();
	lock_guard<mutex> guard(global.lock);

	// One-shot fetch on first call. Set `fetched` before making the request so
	// that a network error doesn't trap us in a per-chunk retry loop.
	if (!global.fetched) {
		global.fetched = true;
		std::ostringstream url;
		url << miint::ENAParser::PORTAL_BASE << "/returnFields?result=" << bind_data.result_type << "&format=tsv";
		auto tsv = global.client->FetchURL(url.str());
		auto parsed = miint::ENAParser::ParseTSV(tsv);

		// Locate the columns by header name so we're robust to ENA reordering.
		int col_id = -1;
		int col_type = -1;
		int col_desc = -1;
		for (size_t i = 0; i < parsed.column_names.size(); i++) {
			const auto &name = parsed.column_names[i];
			if (name == "columnId" || name == "field_name") {
				col_id = static_cast<int>(i);
			} else if (name == "type") {
				col_type = static_cast<int>(i);
			} else if (name == "description") {
				col_desc = static_cast<int>(i);
			}
		}
		if (col_id < 0) {
			throw IOException(
			    "ena_searchable_fields: /returnFields response missing expected 'columnId' header (result_type='%s')",
			    bind_data.result_type);
		}
		for (const auto &row : parsed.rows) {
			std::array<std::string, 3> out {};
			if (col_id >= 0 && static_cast<size_t>(col_id) < row.size()) {
				out[0] = row[col_id];
			}
			if (col_type >= 0 && static_cast<size_t>(col_type) < row.size()) {
				out[1] = row[col_type];
			}
			if (col_desc >= 0 && static_cast<size_t>(col_desc) < row.size()) {
				out[2] = row[col_desc];
			}
			if (out[0].empty()) {
				// Skip rows with no field name — unlikely but a corrupt TSV
				// row shouldn't show up as a blank in the output.
				continue;
			}
			global.rows.push_back(std::move(out));
		}
	}

	idx_t remaining = global.rows.size() - global.offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);
	if (count == 0) {
		output.SetCardinality(0);
		return;
	}

	auto name_data = FlatVector::GetData<string_t>(output.data[0]);
	auto type_data = FlatVector::GetData<string_t>(output.data[1]);
	auto desc_data = FlatVector::GetData<string_t>(output.data[2]);
	auto &desc_validity = FlatVector::Validity(output.data[2]);

	for (idx_t i = 0; i < count; i++) {
		const auto &row = global.rows[global.offset + i];
		name_data[i] = StringVector::AddString(output.data[0], row[0]);
		type_data[i] = StringVector::AddString(output.data[1], row[1]);
		if (row[2].empty()) {
			desc_validity.SetInvalid(i);
		} else {
			desc_data[i] = StringVector::AddString(output.data[2], row[2]);
		}
	}

	global.offset += count;
	output.SetCardinality(count);
}

// ---- Registration ----

TableFunction ReadENASearchableFieldsTableFunction::GetFunction() {
	return TableFunction("ena_searchable_fields", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal, InitLocal);
}

void ReadENASearchableFieldsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
