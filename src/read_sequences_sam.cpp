#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_sequences_sam.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadSequencesSamTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<duckdb::LogicalType> &return_types,
                                                             vector<std::string> &names) {

	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> file_paths;

	// Handle VARCHAR (single path, potentially a glob) or VARCHAR[] (array of literal paths)
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		file_paths = std::move(result.paths);
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			file_paths.push_back(child.ToString());
		}
		if (file_paths.empty()) {
			throw InvalidInputException("read_sequences_sam: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_sequences_sam: first argument must be VARCHAR or VARCHAR[]");
	}

	// Detect stdin usage
	bool uses_stdin = false;
	if (file_paths.size() == 1 && IsStdinPath(file_paths[0])) {
		uses_stdin = true;
		if (file_paths[0] == "-") {
			file_paths[0] = "/dev/stdin";
		}
	} else if (file_paths.size() > 1) {
		for (const auto &path : file_paths) {
			if (IsStdinPath(path)) {
				throw InvalidInputException(
				    "stdin ('-' or '/dev/stdin') must be a single file path, not part of an array");
			}
		}
	}

	// Validate all files exist (skip stdin)
	for (const auto &path : file_paths) {
		if (!IsStdinPath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	bool include_filepath = ParseIncludeFilepathParameter(input.named_parameters);

	auto data = duckdb::make_uniq<Data>(file_paths, include_filepath, uses_stdin);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return data;
}

unique_ptr<GlobalTableFunctionState> ReadSequencesSamTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	return duckdb::make_uniq<GlobalState>(data.file_paths, data.uses_stdin);
}

unique_ptr<LocalTableFunctionState> ReadSequencesSamTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadSequencesSamTableFunction::SetResultVectorNull(Vector &result_vector) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	result_vector.SetValue(0, Value());
}

void ReadSequencesSamTableFunction::SetResultVectorString(Vector &result_vector,
                                                          const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
	}
}

void ReadSequencesSamTableFunction::SetResultVectorStringNullable(Vector &result_vector,
                                                                  const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
		if (!values[j].empty()) {
			validity.SetValid(j);
		}
	}
}

void ReadSequencesSamTableFunction::SetResultVectorListUInt8(Vector &result_vector,
                                                             const std::vector<miint::QualScore> &values) {
	// SAMReader stores quality as Phred33 ASCII; as_vec converts to raw Phred scores (0-93)
	constexpr uint8_t PHRED33_OFFSET = 33;

	// Per-record nullability: BAM records may individually have or lack quality scores
	auto &validity = FlatVector::Validity(result_vector);
	auto list_entries = FlatVector::GetData<list_entry_t>(result_vector);

	// Compute total child elements, accounting for records without quality
	idx_t total_child_elements = 0;
	for (auto &qual : values) {
		total_child_elements += qual.as_string().length();
	}

	ListVector::Reserve(result_vector, total_child_elements);
	ListVector::SetListSize(result_vector, total_child_elements);

	auto &child_vector = ListVector::GetEntry(result_vector);
	auto child_data = FlatVector::GetData<uint8_t>(child_vector);

	const auto output_count = values.size();
	validity.SetAllValid(output_count);
	idx_t value_offset = 0;
	for (idx_t row_offset = 0; row_offset < output_count; row_offset++) {
		if (values[row_offset].as_string().empty()) {
			list_entries[row_offset].offset = value_offset;
			list_entries[row_offset].length = 0;
			validity.SetInvalid(row_offset);
		} else {
			const auto qual_data = values[row_offset].as_vec(PHRED33_OFFSET);
			list_entries[row_offset].offset = value_offset;
			list_entries[row_offset].length = qual_data.size();

			for (auto &value : qual_data) {
				child_data[value_offset++] = value;
			}
		}
	}

	auto &child_validity = FlatVector::Validity(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

void ReadSequencesSamTableFunction::SetResultVectorFilepath(Vector &result_vector, const std::string &filepath) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void ReadSequencesSamTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SAMRecordBatch batch;
	std::string current_filepath;

	// Loop until we get data or run out of files
	while (true) {
		// If this thread doesn't have a file, claim one
		if (!local_state.has_file) {
			lock_guard<mutex> read_lock(global_state.lock);

			if (global_state.next_file_idx >= global_state.readers.size()) {
				output.SetCardinality(0);
				return;
			}

			local_state.current_file_idx = global_state.next_file_idx;
			global_state.next_file_idx++;
			local_state.has_file = true;
		}

		batch = global_state.readers[local_state.current_file_idx]->read(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		break;
	}

	// Get sequence indices for this chunk
	uint64_t start_sequence_index = global_state.file_sequence_counters[local_state.current_file_idx];
	global_state.file_sequence_counters[local_state.current_file_idx] += batch.size();

	// Set sequence_index column
	auto &sequence_index_vector = output.data[0];
	auto sequence_index_data = FlatVector::GetData<int64_t>(sequence_index_vector);
	for (idx_t j = 0; j < batch.size(); j++) {
		sequence_index_data[j] = static_cast<int64_t>(start_sequence_index + j);
	}

	size_t field_idx = 1;

	// read_id from BAM QNAME
	SetResultVectorString(output.data[field_idx++], batch.read_ids);

	// comment: always NULL (BAM has no comment field)
	SetResultVectorNull(output.data[field_idx++]);

	// sequence1 from BAM SEQ
	SetResultVectorString(output.data[field_idx++], batch.sequences);

	// sequence2: always NULL (single-end)
	SetResultVectorNull(output.data[field_idx++]);

	// qual1 from BAM QUAL (per-record nullability)
	SetResultVectorListUInt8(output.data[field_idx++], batch.quals);

	// qual2: always NULL
	SetResultVectorNull(output.data[field_idx++]);

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[field_idx++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

TableFunction ReadSequencesSamTableFunction::GetFunction() {
	auto tf = TableFunction("read_sequences_sam", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadSequencesSamTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
} // namespace duckdb
