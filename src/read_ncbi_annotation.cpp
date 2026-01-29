#include "read_ncbi_annotation.hpp"
#include "duckdb/common/vector_size.hpp"
#include <sstream>

namespace duckdb {

// Data constructor - sets up schema matching read_gff
ReadNCBIAnnotationTableFunction::Data::Data(std::vector<std::string> accessions, const std::string &api_key,
                                            bool include_filepath)
    : accessions(std::move(accessions)), api_key(api_key), include_filepath(include_filepath) {

	// Schema matches read_gff macro output
	names = {"seqid", "source", "type", "position", "stop_position", "score", "strand", "phase", "attributes"};
	types = {LogicalType::VARCHAR,                                          // seqid
	         LogicalType::VARCHAR,                                          // source
	         LogicalType::VARCHAR,                                          // type
	         LogicalType::INTEGER,                                          // position
	         LogicalType::INTEGER,                                          // stop_position
	         LogicalType::DOUBLE,                                           // score
	         LogicalType::VARCHAR,                                          // strand
	         LogicalType::INTEGER,                                          // phase
	         LogicalType::MAP(LogicalType::VARCHAR, LogicalType::VARCHAR)}; // attributes

	if (include_filepath) {
		names.push_back("filepath");
		types.push_back(LogicalType::VARCHAR);
	}
}

// GlobalState constructor
ReadNCBIAnnotationTableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                          const std::vector<std::string> &accessions)
    : client(make_uniq<miint::NCBIClient>(db, api_key)), next_accession_idx(0), batch_offset(0),
      accessions(accessions) {
}

bool ReadNCBIAnnotationTableFunction::GlobalState::FetchNextAccession() {
	if (next_accession_idx >= accessions.size()) {
		return false;
	}

	current_accession = accessions[next_accession_idx];
	next_accession_idx++;

	// Fetch feature table from NCBI
	std::string feature_table = client->FetchFeatureTable(current_accession);

	// Check for empty response (indicates invalid accession or no data)
	if (feature_table.empty()) {
		throw IOException("read_ncbi_annotation: no data returned for accession '%s' - verify accession is valid",
		                  current_accession);
	}

	// Parse into annotation batch
	current_batch = miint::NCBIParser::ParseFeatureTable(feature_table);
	batch_offset = 0;

	// Empty batch after parsing is valid (accession exists but has no features)
	return !current_batch.empty();
}

unique_ptr<FunctionData> ReadNCBIAnnotationTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
		throw InvalidInputException("read_ncbi_annotation: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ncbi_annotation: at least one accession must be provided");
	}

	// Validate that no accession is empty
	for (const auto &acc : accessions) {
		if (acc.empty()) {
			throw InvalidInputException("read_ncbi_annotation: accession cannot be empty");
		}
	}

	// Parse api_key parameter (optional)
	std::string api_key;
	auto api_key_param = input.named_parameters.find("api_key");
	if (api_key_param != input.named_parameters.end() && !api_key_param->second.IsNull()) {
		api_key = api_key_param->second.ToString();
	}

	// Parse include_filepath parameter (optional)
	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	auto data = make_uniq<Data>(std::move(accessions), api_key, include_filepath);

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

unique_ptr<GlobalTableFunctionState> ReadNCBIAnnotationTableFunction::InitGlobal(ClientContext &context,
                                                                                 TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.accessions);
}

unique_ptr<LocalTableFunctionState> ReadNCBIAnnotationTableFunction::InitLocal(ExecutionContext &context,
                                                                               TableFunctionInitInput &input,
                                                                               GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadNCBIAnnotationTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// If current batch is exhausted, fetch next accession
	while (global_state.current_batch.empty() || global_state.batch_offset >= global_state.current_batch.size()) {
		if (!global_state.FetchNextAccession()) {
			// No more accessions
			output.SetCardinality(0);
			return;
		}
	}

	// Determine how many records to output
	idx_t remaining = global_state.current_batch.size() - global_state.batch_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	// Fill output vectors
	auto &features = global_state.current_batch.features;
	size_t offset = global_state.batch_offset;

	// Pre-calculate total attribute count for MapVector reservation (#12)
	size_t total_attrs = 0;
	for (idx_t i = 0; i < count; i++) {
		total_attrs += features[offset + i].attrs.size();
	}

	// Get map vectors and pre-reserve capacity
	auto &map_vec = output.data[8];
	auto &map_keys = MapVector::GetKeys(map_vec);
	auto &map_values = MapVector::GetValues(map_vec);
	auto initial_offset = ListVector::GetListSize(map_vec);
	ListVector::Reserve(map_vec, initial_offset + total_attrs);

	size_t attr_offset = initial_offset;

	for (idx_t i = 0; i < count; i++) {
		const auto &feat = features[offset + i];

		// seqid (column 0)
		FlatVector::GetData<string_t>(output.data[0])[i] = StringVector::AddString(output.data[0], feat.seqid);

		// source (column 1)
		FlatVector::GetData<string_t>(output.data[1])[i] = StringVector::AddString(output.data[1], feat.source);

		// type (column 2)
		FlatVector::GetData<string_t>(output.data[2])[i] = StringVector::AddString(output.data[2], feat.type);

		// position (column 3) - cast from int64_t to int32_t for schema compatibility
		FlatVector::GetData<int32_t>(output.data[3])[i] = static_cast<int32_t>(feat.position);

		// stop_position (column 4) - cast from int64_t to int32_t for schema compatibility
		FlatVector::GetData<int32_t>(output.data[4])[i] = static_cast<int32_t>(feat.stop_position);

		// score (column 5) - may be NULL
		if (feat.has_score) {
			FlatVector::GetData<double>(output.data[5])[i] = feat.score;
		} else {
			FlatVector::Validity(output.data[5]).SetInvalid(i);
		}

		// strand (column 6)
		FlatVector::GetData<string_t>(output.data[6])[i] = StringVector::AddString(output.data[6], feat.strand);

		// phase (column 7) - may be NULL
		if (feat.phase >= 0) {
			FlatVector::GetData<int32_t>(output.data[7])[i] = feat.phase;
		} else {
			FlatVector::Validity(output.data[7]).SetInvalid(i);
		}

		// attributes (column 8) - MAP(VARCHAR, VARCHAR)
		auto list_size = feat.attrs.size();

		// Set up list entry for this row
		auto entry_data = FlatVector::GetData<list_entry_t>(map_vec);
		entry_data[i].offset = attr_offset;
		entry_data[i].length = list_size;

		// Add key-value pairs
		for (size_t j = 0; j < list_size; j++) {
			FlatVector::GetData<string_t>(map_keys)[attr_offset + j] =
			    StringVector::AddString(map_keys, feat.attrs[j].first);
			FlatVector::GetData<string_t>(map_values)[attr_offset + j] =
			    StringVector::AddString(map_values, feat.attrs[j].second);
		}

		attr_offset += list_size;
	}

	// Update total size once at the end
	ListVector::SetListSize(map_vec, attr_offset);

	// filepath (column 9) - optional
	if (bind_data.include_filepath) {
		// Use NCBI URL as filepath (use stored current_accession to avoid off-by-one)
		std::ostringstream url;
		url << miint::NCBIClient::EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << global_state.current_accession
		    << "&rettype=ft";
		output.data[9].SetVectorType(VectorType::CONSTANT_VECTOR);
		auto filepath_data = ConstantVector::GetData<string_t>(output.data[9]);
		*filepath_data = StringVector::AddString(output.data[9], url.str());
	}

	global_state.batch_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBIAnnotationTableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi_annotation", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadNCBIAnnotationTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
