#include "read_ena_attributes.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

// ---- Data ----

ReadENAAttributesTableFunction::Data::Data(std::vector<std::string> accessions)
    : accessions(std::move(accessions)), names({"sample_accession", "tag", "value"}),
      types({LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}) {
}

// ---- GlobalState ----

ReadENAAttributesTableFunction::GlobalState::GlobalState(DatabaseInstance &db)
    : client(make_uniq<miint::ENAClient>(db)), row_offset(0), fetched(false) {
}

std::vector<std::string>
ReadENAAttributesTableFunction::GlobalState::ResolveSampleAccessions(const std::vector<std::string> &accessions) {
	std::vector<std::string> sample_accs;

	for (const auto &acc : accessions) {
		auto acc_type = miint::ENAParser::DetectAccessionType(acc);

		if (acc_type == miint::ENAAccessionType::SAMPLE) {
			sample_accs.push_back(acc);
		} else if (acc_type == miint::ENAAccessionType::STUDY) {
			// Get all sample accessions for the study
			auto tsv = client->Search(acc, "sample", "sample_accession");
			auto parsed = miint::ENAParser::ParseTSV(tsv);
			for (size_t i = 0; i < parsed.column_names.size(); i++) {
				if (parsed.column_names[i] == "sample_accession") {
					for (const auto &row : parsed.rows) {
						if (i < row.size() && !row[i].empty()) {
							sample_accs.push_back(row[i]);
						}
					}
					break;
				}
			}
		} else {
			// Run or experiment: resolve to sample_accession via read_run
			auto tsv = client->Search(acc, "read_run", "sample_accession");
			auto parsed = miint::ENAParser::ParseTSV(tsv);
			for (size_t i = 0; i < parsed.column_names.size(); i++) {
				if (parsed.column_names[i] == "sample_accession") {
					for (const auto &row : parsed.rows) {
						if (i < row.size() && !row[i].empty()) {
							sample_accs.push_back(row[i]);
						}
					}
					break;
				}
			}
		}
	}

	return sample_accs;
}

void ReadENAAttributesTableFunction::GlobalState::FetchAttributes(const std::vector<std::string> &accessions) {
	// Mark as fetched immediately to avoid infinite retry loops on network errors
	fetched = true;

	auto sample_accs = ResolveSampleAccessions(accessions);

	if (sample_accs.empty()) {
		return;
	}

	// NOTE: This materializes all attributes before returning any rows.
	// For large studies (thousands of samples), this may use significant memory
	// and block until all XML batches are fetched. A streaming approach would
	// be better for very large studies.
	static constexpr size_t BATCH_SIZE = 50;
	for (size_t start = 0; start < sample_accs.size(); start += BATCH_SIZE) {
		size_t end = std::min(start + BATCH_SIZE, sample_accs.size());
		std::vector<std::string> batch(sample_accs.begin() + start, sample_accs.begin() + end);

		auto xml = client->FetchXML(batch);
		auto batch_attrs = miint::ENAParser::ParseSampleAttributesXML(xml);

		for (auto &attr : batch_attrs) {
			attributes.push_back(std::move(attr));
		}
	}
}

// ---- Bind ----

unique_ptr<FunctionData> ReadENAAttributesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                              vector<LogicalType> &return_types,
                                                              vector<std::string> &names) {
	std::vector<std::string> accessions;
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ena_attributes: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena_attributes: at least one accession must be provided");
	}
	for (auto &acc : accessions) {
		// Trim whitespace
		auto start = acc.find_first_not_of(" \t\n\r");
		if (start == std::string::npos) {
			acc.clear();
		} else {
			acc = acc.substr(start, acc.find_last_not_of(" \t\n\r") - start + 1);
		}
		if (acc.empty()) {
			throw InvalidInputException("read_ena_attributes: accession cannot be empty");
		}
	}

	auto data = make_uniq<Data>(std::move(accessions));

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENAAttributesTableFunction::InitGlobal(ClientContext &context,
                                                                                TableFunctionInitInput &input) {
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db);
}

unique_ptr<LocalTableFunctionState> ReadENAAttributesTableFunction::InitLocal(ExecutionContext &context,
                                                                              TableFunctionInitInput &input,
                                                                              GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENAAttributesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// Fetch all attributes on first call
	if (!global_state.fetched) {
		global_state.FetchAttributes(bind_data.accessions);
	}

	if (global_state.row_offset >= global_state.attributes.size()) {
		output.SetCardinality(0);
		return;
	}

	idx_t remaining = global_state.attributes.size() - global_state.row_offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	auto &attrs = global_state.attributes;
	size_t offset = global_state.row_offset;

	// sample_accession (column 0)
	auto acc_data = FlatVector::GetData<string_t>(output.data[0]);
	for (idx_t i = 0; i < count; i++) {
		acc_data[i] = StringVector::AddString(output.data[0], attrs[offset + i].sample_accession);
	}

	// tag (column 1)
	auto tag_data = FlatVector::GetData<string_t>(output.data[1]);
	for (idx_t i = 0; i < count; i++) {
		tag_data[i] = StringVector::AddString(output.data[1], attrs[offset + i].tag);
	}

	// value (column 2)
	auto val_data = FlatVector::GetData<string_t>(output.data[2]);
	auto &val_validity = FlatVector::Validity(output.data[2]);
	for (idx_t i = 0; i < count; i++) {
		const auto &val = attrs[offset + i].value;
		if (val.empty()) {
			val_validity.SetInvalid(i);
		} else {
			val_data[i] = StringVector::AddString(output.data[2], val);
		}
	}

	global_state.row_offset += count;
	output.SetCardinality(count);
}

// ---- Registration ----

TableFunction ReadENAAttributesTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_attributes", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	return tf;
}

void ReadENAAttributesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
