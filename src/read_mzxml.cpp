#include "MzXMLReader.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_mzxml.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadMzXMLTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<duckdb::LogicalType> &return_types,
                                                      vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> file_paths;

	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		file_paths = std::move(result.paths);
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			file_paths.push_back(child.ToString());
		}
		if (file_paths.empty()) {
			throw InvalidInputException("read_mzxml: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_mzxml: first argument must be VARCHAR or VARCHAR[]");
	}

	for (const auto &path : file_paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_mzxml: stdin is not supported");
		}
	}

	for (const auto &path : file_paths) {
		if (!miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	bool include_filepath = ParseIncludeFilepathParameter(input.named_parameters);

	auto data = duckdb::make_uniq<Data>(file_paths, include_filepath);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return data;
}

unique_ptr<GlobalTableFunctionState> ReadMzXMLTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &fs = FileSystem::GetFileSystem(context);
	return duckdb::make_uniq<GlobalState>(data.file_paths, fs);
}

unique_ptr<LocalTableFunctionState> ReadMzXMLTableFunction::InitLocal(ExecutionContext &context,
                                                                      TableFunctionInitInput &input,
                                                                      GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadMzXMLTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::MzMLSpectrumBatch batch;
	std::string current_filepath;

	while (true) {
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

		if (!global_state.readers[local_state.current_file_idx]) {
			global_state.readers[local_state.current_file_idx] = std::make_unique<miint::MzXMLReader>(
			    global_state.fs, global_state.filepaths[local_state.current_file_idx]);
		}

		batch = global_state.readers[local_state.current_file_idx]->read_spectra(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		break;
	}

	PopulateSpectrumBatchOutput(output, batch, bind_data.include_filepath, current_filepath);
}

TableFunction ReadMzXMLTableFunction::GetFunction() {
	auto tf = TableFunction("read_mzxml", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadMzXMLTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
} // namespace duckdb
