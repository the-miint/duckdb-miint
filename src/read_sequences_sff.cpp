#include "SFFReader.hpp"
#include "SequenceRecord.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_sequences_sff.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadSequencesSFFTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
			throw InvalidInputException("read_sequences_sff: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_sequences_sff: first argument must be VARCHAR or VARCHAR[]");
	}

	// SFF is a binary format that requires seeking - stdin is not supported
	for (const auto &path : file_paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_sequences_sff: stdin is not supported (SFF is a binary format "
			                            "that requires file seeking)");
		}
	}

	// Validate all files exist
	for (const auto &path : file_paths) {
		if (!fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	bool include_filepath = ParseIncludeFilepathParameter(input.named_parameters);

	// Parse trim parameter (default: true)
	bool trim = true;
	auto trim_param = input.named_parameters.find("trim");
	if (trim_param != input.named_parameters.end() && !trim_param->second.IsNull()) {
		trim = trim_param->second.GetValue<bool>();
	}

	auto data = duckdb::make_uniq<Data>(file_paths, include_filepath, trim);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return data;
}

unique_ptr<GlobalTableFunctionState> ReadSequencesSFFTableFunction::InitGlobal(ClientContext &context,
                                                                               TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	return duckdb::make_uniq<GlobalState>(data.file_paths, data.trim);
}

unique_ptr<LocalTableFunctionState> ReadSequencesSFFTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadSequencesSFFTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SequenceRecordBatch batch;
	std::string current_filepath;

	// Loop until we get data or run out of files
	while (true) {
		// If this thread doesn't have a file, claim one
		if (!local_state.has_file) {
			lock_guard<mutex> read_lock(global_state.lock);

			if (global_state.next_file_idx >= global_state.filepaths.size()) {
				output.SetCardinality(0);
				return;
			}

			local_state.current_file_idx = global_state.next_file_idx;
			global_state.next_file_idx++;
			local_state.has_file = true;
		}

		// Lazily open the file (safe without lock - this thread has exclusive access to this file index)
		if (!global_state.readers[local_state.current_file_idx]) {
			global_state.readers[local_state.current_file_idx] = std::make_unique<miint::SFFReader>(
			    global_state.filepaths[local_state.current_file_idx], global_state.trim);
		}

		batch = global_state.readers[local_state.current_file_idx]->read(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		break;
	}

	// No lock needed - this thread has exclusive access to this file index
	uint64_t start_sequence_index = global_state.file_sequence_counters[local_state.current_file_idx];
	global_state.file_sequence_counters[local_state.current_file_idx] += batch.size();

	// Set sequence_index column
	auto &sequence_index_vector = output.data[0];
	auto sequence_index_data = FlatVector::GetData<int64_t>(sequence_index_vector);
	for (idx_t j = 0; j < batch.size(); j++) {
		sequence_index_data[j] = static_cast<int64_t>(start_sequence_index + j);
	}

	size_t field_idx = 1;

	// read_id
	SetResultVectorString(output.data[field_idx++], batch.read_ids);

	// comment: always NULL (SFF has no comment field)
	SetResultVectorNull(output.data[field_idx++]);

	// sequence1
	SetResultVectorString(output.data[field_idx++], batch.sequences1);

	// sequence2: always NULL (SFF is single-end)
	SetResultVectorNull(output.data[field_idx++]);

	// qual1: SFF always has quality scores, use offset 33 (stored internally as Phred+33)
	SetResultVectorListUInt8(output.data[field_idx++], batch.quals1, 33);

	// qual2: always NULL
	SetResultVectorNull(output.data[field_idx++]);

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[field_idx++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

TableFunction ReadSequencesSFFTableFunction::GetFunction() {
	auto tf = TableFunction("read_sequences_sff", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["trim"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadSequencesSFFTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
} // namespace duckdb
