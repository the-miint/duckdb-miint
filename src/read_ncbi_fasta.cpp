#include "read_ncbi_fasta.hpp"
#include "duckdb/common/vector_size.hpp"
#include <sstream>

namespace duckdb {

// Data constructor - sets up schema matching read_fastx
ReadNCBIFastaTableFunction::Data::Data(std::vector<std::string> accessions, const std::string &api_key,
                                        bool include_filepath)
    : accessions(std::move(accessions)), api_key(api_key), include_filepath(include_filepath) {

	// Schema matches read_fastx exactly
	names = {"sequence_index", "read_id", "comment", "sequence1", "sequence2", "qual1", "qual2"};
	types = {LogicalType::BIGINT,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::VARCHAR,
	         LogicalType::LIST(LogicalType::UTINYINT),
	         LogicalType::LIST(LogicalType::UTINYINT)};

	if (include_filepath) {
		names.push_back("filepath");
		types.push_back(LogicalType::VARCHAR);
	}
}

// GlobalState constructor
ReadNCBIFastaTableFunction::GlobalState::GlobalState(DatabaseInstance &db, const std::string &api_key,
                                                      const std::vector<std::string> &accessions)
    : client(make_uniq<miint::NCBIClient>(db, api_key)), next_accession_idx(0), batch_offset(0), sequence_index(0),
      accessions(accessions) {
}

bool ReadNCBIFastaTableFunction::GlobalState::FetchNextAccession() {
	if (next_accession_idx >= accessions.size()) {
		return false;
	}

	current_accession = accessions[next_accession_idx];
	next_accession_idx++;

	std::string fasta_text;

	// Check accession type
	auto acc_type = miint::NCBIParser::DetectAccessionType(current_accession);

	if (acc_type == miint::AccessionType::ASSEMBLY) {
		// Fetch from Datasets API (returns ZIP with FASTA)
		fasta_text = client->FetchAssemblyFasta(current_accession);
	} else {
		// Fetch from E-utilities (existing behavior)
		fasta_text = client->FetchFasta(current_accession);
	}

	// Check for empty response (indicates invalid accession or no data)
	if (fasta_text.empty()) {
		throw IOException("read_ncbi_fasta: no data returned for accession '%s' - verify accession is valid",
		                  current_accession);
	}

	// Parse into batch
	current_batch = miint::NCBIParser::ParseFasta(fasta_text);
	batch_offset = 0;

	// Empty batch after parsing means no valid sequences found
	return !current_batch.empty();
}

unique_ptr<FunctionData> ReadNCBIFastaTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
		throw InvalidInputException("read_ncbi_fasta: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ncbi_fasta: at least one accession must be provided");
	}

	// Validate that no accession is empty
	for (const auto &acc : accessions) {
		if (acc.empty()) {
			throw InvalidInputException("read_ncbi_fasta: accession cannot be empty");
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

unique_ptr<GlobalTableFunctionState> ReadNCBIFastaTableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db, data.api_key, data.accessions);
}

unique_ptr<LocalTableFunctionState> ReadNCBIFastaTableFunction::InitLocal(ExecutionContext &context,
                                                                           TableFunctionInitInput &input,
                                                                           GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

void ReadNCBIFastaTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// If current batch is exhausted, fetch next accession
	while (global_state.current_batch.empty() ||
	       global_state.batch_offset >= global_state.current_batch.size()) {
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
	auto &batch = global_state.current_batch;
	size_t offset = global_state.batch_offset;

	// sequence_index (column 0)
	auto seq_idx_data = FlatVector::GetData<int64_t>(output.data[0]);
	for (idx_t i = 0; i < count; i++) {
		seq_idx_data[i] = global_state.sequence_index++;
	}

	// read_id (column 1)
	auto read_id_data = FlatVector::GetData<string_t>(output.data[1]);
	for (idx_t i = 0; i < count; i++) {
		read_id_data[i] = StringVector::AddString(output.data[1], batch.read_ids[offset + i]);
	}

	// comment (column 2) - nullable
	auto comment_data = FlatVector::GetData<string_t>(output.data[2]);
	auto &comment_validity = FlatVector::Validity(output.data[2]);
	for (idx_t i = 0; i < count; i++) {
		const auto &comment = batch.comments[offset + i];
		comment_data[i] = StringVector::AddString(output.data[2], comment);
		if (comment.empty()) {
			comment_validity.SetInvalid(i);
		}
	}

	// sequence1 (column 3)
	auto seq1_data = FlatVector::GetData<string_t>(output.data[3]);
	for (idx_t i = 0; i < count; i++) {
		seq1_data[i] = StringVector::AddString(output.data[3], batch.sequences1[offset + i]);
	}

	// sequence2 (column 4) - always NULL for NCBI FASTA
	output.data[4].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[4], true);

	// qual1 (column 5) - always NULL for NCBI FASTA (no quality scores)
	output.data[5].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[5], true);

	// qual2 (column 6) - always NULL
	output.data[6].SetVectorType(VectorType::CONSTANT_VECTOR);
	ConstantVector::SetNull(output.data[6], true);

	// filepath (column 7) - optional
	if (bind_data.include_filepath) {
		// Use NCBI URL as filepath (use stored current_accession to avoid off-by-one)
		std::ostringstream url;
		url << miint::NCBIClient::EUTILS_BASE << "/efetch.fcgi?db=nuccore&id=" << global_state.current_accession
		    << "&rettype=fasta";
		output.data[7].SetVectorType(VectorType::CONSTANT_VECTOR);
		auto filepath_data = ConstantVector::GetData<string_t>(output.data[7]);
		*filepath_data = StringVector::AddString(output.data[7], url.str());
	}

	global_state.batch_offset += count;
	output.SetCardinality(count);
}

TableFunction ReadNCBIFastaTableFunction::GetFunction() {
	auto tf = TableFunction("read_ncbi_fasta", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["api_key"] = LogicalType::VARCHAR;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadNCBIFastaTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
