#include "read_ncbi.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

// Data constructor - sets up schema
ReadNCBITableFunction::Data::Data(std::vector<std::string> accessions, const std::string &api_key)
    : accessions(std::move(accessions)), api_key(api_key) {

	// Schema for GenBank metadata
	names = {"accession", "version", "description", "organism", "taxonomy_id", "length", "molecule_type", "update_date"};
	types = {LogicalType::VARCHAR,  // accession
	         LogicalType::INTEGER,  // version
	         LogicalType::VARCHAR,  // description
	         LogicalType::VARCHAR,  // organism
	         LogicalType::BIGINT,   // taxonomy_id
	         LogicalType::BIGINT,   // length
	         LogicalType::VARCHAR,  // molecule_type
	         LogicalType::DATE};    // update_date
}

// GlobalState constructor
ReadNCBITableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                  const std::vector<std::string> &accessions)
    : client(make_uniq<miint::NCBIClient>(db, api_key)), next_accession_idx(0), result_offset(0),
      accessions(accessions) {
}

bool ReadNCBITableFunction::GlobalState::FetchNextAccession() {
	if (next_accession_idx >= accessions.size()) {
		return false;
	}

	const auto &accession = accessions[next_accession_idx];
	next_accession_idx++;

	// Fetch GenBank XML from NCBI
	std::string xml = client->FetchGenBankXML(accession);

	// Check for empty response (indicates invalid accession or no data)
	if (xml.empty()) {
		throw IOException("read_ncbi: no data returned for accession '%s' - verify accession is valid", accession);
	}

	// Parse into metadata
	auto metadata = miint::NCBIParser::ParseGenBankXML(xml);

	// Only add if we got valid data
	if (!metadata.accession.empty()) {
		metadata_results.push_back(std::move(metadata));
	}

	return true;
}

unique_ptr<FunctionData> ReadNCBITableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<LogicalType> &return_types,
                                                      vector<std::string> &names) {
	// Parse accession(s) - can be VARCHAR or VARCHAR[]
	std::vector<std::string> accessions;

	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ncbi: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ncbi: at least one accession must be provided");
	}

	// Validate that no accession is empty
	for (const auto &acc : accessions) {
		if (acc.empty()) {
			throw InvalidInputException("read_ncbi: accession cannot be empty");
		}
	}

	// Parse api_key parameter (optional)
	std::string api_key;
	auto api_key_param = input.named_parameters.find("api_key");
	if (api_key_param != input.named_parameters.end() && !api_key_param->second.IsNull()) {
		api_key = api_key_param->second.ToString();
	}

	auto data = make_uniq<Data>(std::move(accessions), api_key);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadNCBITableFunction::InitGlobal(ClientContext &context,
                                                                         TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.accessions);
}

unique_ptr<LocalTableFunctionState> ReadNCBITableFunction::InitLocal(ExecutionContext &context,
                                                                       TableFunctionInitInput &input,
                                                                       GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadNCBITableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// Keep fetching until we have data or run out of accessions
	while (global_state.result_offset >= global_state.metadata_results.size()) {
		if (!global_state.FetchNextAccession()) {
			// No more accessions
			output.SetCardinality(0);
			return;
		}
	}

	// Determine how many records to output
	idx_t remaining = global_state.metadata_results.size() - global_state.result_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	// Fill output vectors
	size_t offset = global_state.result_offset;

	for (idx_t i = 0; i < count; i++) {
		const auto &meta = global_state.metadata_results[offset + i];

		// accession (column 0)
		FlatVector::GetData<string_t>(output.data[0])[i] =
		    StringVector::AddString(output.data[0], meta.accession);

		// version (column 1)
		FlatVector::GetData<int32_t>(output.data[1])[i] = meta.version;

		// description (column 2)
		FlatVector::GetData<string_t>(output.data[2])[i] =
		    StringVector::AddString(output.data[2], meta.description);

		// organism (column 3)
		FlatVector::GetData<string_t>(output.data[3])[i] =
		    StringVector::AddString(output.data[3], meta.organism);

		// taxonomy_id (column 4)
		FlatVector::GetData<int64_t>(output.data[4])[i] = meta.taxonomy_id;

		// length (column 5)
		FlatVector::GetData<int64_t>(output.data[5])[i] = meta.length;

		// molecule_type (column 6)
		FlatVector::GetData<string_t>(output.data[6])[i] =
		    StringVector::AddString(output.data[6], meta.molecule_type);

		// update_date (column 7) - parse YYYY-MM-DD to DATE
		if (!meta.update_date.empty()) {
			date_t date;
			idx_t pos;
			bool special;
			auto result = Date::TryConvertDate(meta.update_date.c_str(), meta.update_date.size(), pos, date, special);
			if (result == DateCastResult::SUCCESS) {
				FlatVector::GetData<date_t>(output.data[7])[i] = date;
			} else {
				FlatVector::Validity(output.data[7]).SetInvalid(i);
			}
		} else {
			FlatVector::Validity(output.data[7]).SetInvalid(i);
		}
	}

	global_state.result_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBITableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	return tf;
}

void ReadNCBITableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
