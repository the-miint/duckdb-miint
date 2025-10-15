#include "SequenceReader.hpp"
#include "SequenceRecord.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include <read_fastx.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadFastxTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<duckdb::LogicalType> &return_types,
                                                      vector<std::string> &names) {

	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> sequence1_paths;

	// Handle VARCHAR or VARCHAR[] input for sequence1
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		sequence1_paths.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			sequence1_paths.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_fastx: first argument must be VARCHAR or VARCHAR[]");
	}

	// Detect stdin usage
	bool uses_stdin = false;
	auto is_stdin = [](const std::string &path) { return path == "-" || path == "/dev/stdin"; };

	if (sequence1_paths.size() == 1 && is_stdin(sequence1_paths[0])) {
		uses_stdin = true;
		// Normalize stdin to /dev/stdin for kseq++
		if (sequence1_paths[0] == "-") {
			sequence1_paths[0] = "/dev/stdin";
		}
	} else if (sequence1_paths.size() > 1) {
		// Check if stdin is in an array (not allowed)
		for (const auto &path : sequence1_paths) {
			if (is_stdin(path)) {
				throw InvalidInputException("stdin ('-' or '/dev/stdin') must be a single file path, not part of an array");
			}
		}
	}

	// Validate all sequence1 files exist (skip stdin)
	for (const auto &path : sequence1_paths) {
		if (!is_stdin(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	// Parse sequence2 parameter (optional)
	std::optional<std::vector<std::string>> sequence2_paths;
	auto path2_param = input.named_parameters.find("sequence2");
	if (path2_param != input.named_parameters.end() && !path2_param->second.IsNull()) {
		// stdin cannot be used with sequence2
		if (uses_stdin) {
			throw InvalidInputException("stdin cannot be used with sequence2 parameter (paired-end not supported)");
		}

		const auto &path2_value = path2_param->second;

		std::vector<std::string> seq2_paths;
		if (path2_value.type().id() == LogicalTypeId::VARCHAR) {
			seq2_paths.push_back(path2_value.ToString());
		} else if (path2_value.type().id() == LogicalTypeId::LIST) {
			auto &list_children = ListValue::GetChildren(path2_value);
			for (const auto &child : list_children) {
				seq2_paths.push_back(child.ToString());
			}
		} else {
			throw InvalidInputException("sequence2 parameter must be VARCHAR or VARCHAR[]");
		}

		// Validate array length consistency
		if (seq2_paths.size() != sequence1_paths.size()) {
			throw InvalidInputException("Mismatched array lengths: sequence1 has " +
			                            std::to_string(sequence1_paths.size()) + " files, sequence2 has " +
			                            std::to_string(seq2_paths.size()) + " files");
		}

		// Validate all sequence2 files exist
		for (const auto &path : seq2_paths) {
			if (!fs.FileExists(path)) {
				throw IOException("File not found: " + path);
			}
		}

		sequence2_paths = seq2_paths;
	}

	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	uint8_t qual_offset = 33;
	auto offset_param = input.named_parameters.find("qual_offset");
	if (offset_param != input.named_parameters.end() && !offset_param->second.IsNull()) {
		int64_t offset_value = offset_param->second.GetValue<int64_t>();
		if (offset_value != 33 && offset_value != 64) {
			throw InvalidInputException("qual_offset must be 33 or 64");
		}
		qual_offset = static_cast<uint8_t>(offset_value);
	}

	auto data = duckdb::make_uniq<Data>(sequence1_paths, sequence2_paths, include_filepath, uses_stdin, qual_offset);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadFastxTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto gstate = duckdb::make_uniq<GlobalState>(data.sequence1_paths, data.sequence2_paths, data.uses_stdin);

	return std::move(gstate);
}

void ReadFastxTableFunction::SetResultVector(Vector &result_vector, const miint::SequenceRecordField &field,
                                             const std::vector<miint::SequenceRecord> &records, uint8_t qual_offset) {
	auto is_paired = records[0].is_paired;
	if (!is_paired && (field == miint::SequenceRecordField::SEQUENCE2 || field == miint::SequenceRecordField::QUAL2)) {
		SetResultVectorNull(result_vector);
	} else if (field == miint::SequenceRecordField::READ_ID || field == miint::SequenceRecordField::SEQUENCE1 ||
	           field == miint::SequenceRecordField::SEQUENCE2) {
		SetResultVectorString(result_vector, field, records);
	} else if (field == miint::SequenceRecordField::COMMENT) {
		SetResultVectorStringNullable(result_vector, field, records);
	} else if (field == miint::SequenceRecordField::QUAL1 || field == miint::SequenceRecordField::QUAL2) {
		SetResultVectorListUInt8(result_vector, field, records, qual_offset);
	} else {
		throw NotImplementedException("field not supported");
	}
}

void ReadFastxTableFunction::SetResultVectorNull(Vector &result_vector) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	result_vector.SetValue(0, Value());
}

void ReadFastxTableFunction::SetResultVectorString(Vector &result_vector, const miint::SequenceRecordField &field,
                                                   const std::vector<miint::SequenceRecord> &records) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, records[j].GetString(field));
	}
}

void ReadFastxTableFunction::SetResultVectorStringNullable(Vector &result_vector,
                                                           const miint::SequenceRecordField &field,
                                                           const std::vector<miint::SequenceRecord> &records) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(records.size());

	for (idx_t j = 0; j < records.size(); j++) {
		const auto &value = records[j].GetString(field);
		result_data[j] = StringVector::AddString(result_vector, value);

		if (!value.empty()) {
			validity.SetValid(j);
		}
	}
}

void ReadFastxTableFunction::SetResultVectorListUInt8(Vector &result_vector, const miint::SequenceRecordField &field,
                                                      const std::vector<miint::SequenceRecord> &records, uint8_t qual_offset) {
	// Check if this is FASTA (quality scores are empty)
	bool is_fasta = records[0].qual1.as_string().empty();
	bool is_qual2 = (field == miint::SequenceRecordField::QUAL2);

	if (is_fasta || (is_qual2 && !records[0].is_paired)) {
		SetResultVectorNull(result_vector);
		return;
	}

	// compute total number of elements
	idx_t total_child_elements = 0;
	for (auto &rec : records) {
		total_child_elements += rec.GetLength(field);
	}

	ListVector::Reserve(result_vector, total_child_elements);
	ListVector::SetListSize(result_vector, total_child_elements);

	auto &child_vector = ListVector::GetEntry(result_vector);
	auto child_data = FlatVector::GetData<uint8_t>(child_vector);
	auto list_entries = FlatVector::GetData<list_entry_t>(result_vector);

	const auto output_count = records.size();
	idx_t value_offset = 0;
	idx_t row_offset = 0;
	for (idx_t row_offset = 0; row_offset < output_count; row_offset++) {
		const auto qual_data = records[row_offset].GetQual(field).as_vec(qual_offset);
		list_entries[row_offset].offset = value_offset;
		list_entries[row_offset].length = qual_data.size();

		for (auto &value : qual_data) {
			child_data[value_offset++] = value;
		}
	}

	// all entries are not null
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllValid(output_count);

	auto &child_validity = FlatVector::Validity(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

void ReadFastxTableFunction::SetResultVectorFilepath(Vector &result_vector, const std::string &filepath,
                                                     size_t num_records) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void ReadFastxTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	std::vector<miint::SequenceRecord> records;
	std::string current_filepath;
	uint64_t start_sequence_index;

	{
		lock_guard<mutex> read_lock(global_state.lock);

		if (global_state.finished) {
			output.SetCardinality(0);
			return;
		}

		// Try to read from current file
		while (global_state.current_file_idx < global_state.readers.size()) {
			records = global_state.readers[global_state.current_file_idx]->read(STANDARD_VECTOR_SIZE);

			if (records.size() > 0) {
				current_filepath = global_state.sequence1_filepaths[global_state.current_file_idx];
				// Capture starting sequence_index for this chunk and increment counter
				start_sequence_index = global_state.sequence_index_counter;
				global_state.sequence_index_counter += records.size();
				break;
			}

			// Current file exhausted, move to next
			global_state.current_file_idx++;
			if (global_state.current_file_idx < global_state.readers.size()) {
				global_state.current_filepath = global_state.sequence1_filepaths[global_state.current_file_idx];
			}
		}

		// All files exhausted
		if (records.empty()) {
			global_state.finished = true;
			output.SetCardinality(0);
			return;
		}
	}

	// Set sequence_index column (first column, index 0)
	auto &sequence_index_vector = output.data[0];
	auto sequence_index_data = FlatVector::GetData<int64_t>(sequence_index_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		sequence_index_data[j] = static_cast<int64_t>(start_sequence_index + j);
	}

	// Set result vectors for sequence fields (now starting at index 1)
	for (idx_t i = 0; i < bind_data.fields.size(); i++) {
		auto &result_vector = output.data[i + 1];
		auto &field = bind_data.fields[i];
		SetResultVector(result_vector, field, records, bind_data.qual_offset);
	}

	if (bind_data.include_filepath) {
		auto &result_vector = output.data[1 + bind_data.fields.size()];
		SetResultVectorFilepath(result_vector, current_filepath, records.size());
	}

	output.SetCardinality(records.size());
}

TableFunction ReadFastxTableFunction::GetFunction() {
	// IMPORTANT: "set preserve_insertion_order=false" to run in parallel
	// otherwise it will default to serial for many (not all) operations
	auto tf = TableFunction("read_fastx", {LogicalType::ANY}, Execute, Bind, InitGlobal);
	tf.named_parameters["sequence2"] = LogicalType::ANY;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	return tf;
}
}; // namespace duckdb
