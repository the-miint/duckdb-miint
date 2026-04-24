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
    : client(make_uniq<miint::ENAClient>(db)) {
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

void ReadENAAttributesTableFunction::GlobalState::ResolveAccessions(const std::vector<std::string> &accessions) {
	// Mark as resolved immediately so that a network error during the
	// Search call below doesn't cause successive Execute callbacks to keep
	// retrying the same failing request in a loop.
	resolved = true;
	pending_sample_accs = ResolveSampleAccessions(accessions);
	total_samples_expected.store(pending_sample_accs.size(), std::memory_order_relaxed);
	samples_fetched.store(0, std::memory_order_relaxed);
}

bool ReadENAAttributesTableFunction::GlobalState::FetchNextBatch() {
	if (pending_sample_accs.empty()) {
		return false;
	}
	size_t n = std::min<size_t>(BATCH_SIZE, pending_sample_accs.size());
	std::vector<std::string> batch(pending_sample_accs.begin(), pending_sample_accs.begin() + n);
	pending_sample_accs.erase(pending_sample_accs.begin(), pending_sample_accs.begin() + n);

	auto xml = client->FetchXML(batch);
	current_batch = miint::ENAParser::ParseSampleAttributesXML(xml);
	current_batch_offset = 0;
	samples_fetched.fetch_add(n, std::memory_order_relaxed);
	return true;
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

	// First call: resolve input accessions to a list of sample accessions.
	// Subsequent fetches happen lazily one batch at a time below.
	if (!global_state.resolved) {
		global_state.ResolveAccessions(bind_data.accessions);
	}

	// Advance to the next fetched batch if we've drained the current one.
	// Loop to skip empty batches (a 50-sample batch whose XML has no sample
	// with any attributes would otherwise end the scan prematurely).
	while (global_state.current_batch_offset >= global_state.current_batch.size()) {
		if (!global_state.FetchNextBatch()) {
			output.SetCardinality(0);
			return;
		}
	}

	auto &attrs = global_state.current_batch;
	size_t offset = global_state.current_batch_offset;
	idx_t remaining = attrs.size() - offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

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

	global_state.current_batch_offset += count;
	output.SetCardinality(count);
}

// ---- Progress ----

double ReadENAAttributesTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                                const GlobalTableFunctionState *global_state) {
	if (!global_state) {
		return -1.0;
	}
	auto &state = global_state->Cast<GlobalState>();
	// Before ResolveAccessions runs we don't know the scale of the work. Once
	// resolved, report samples-fetched / total. This gives users visible
	// progress for large studies like PRJEB11419 (33k+ samples) where a full
	// scan is long-running by nature.
	size_t total = state.total_samples_expected.load(std::memory_order_relaxed);
	if (total == 0) {
		return -1.0;
	}
	size_t fetched = state.samples_fetched.load(std::memory_order_relaxed);
	double pct = 100.0 * static_cast<double>(fetched) / static_cast<double>(total);
	return std::min(100.0, pct);
}

// ---- Registration ----

TableFunction ReadENAAttributesTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_attributes", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.table_scan_progress = Progress;
	return tf;
}

void ReadENAAttributesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
