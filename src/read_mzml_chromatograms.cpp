#include "MzMLReader.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_mzml_chromatograms.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadMzMLChromatogramsTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
			throw InvalidInputException("read_mzml_chromatograms: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_mzml_chromatograms: first argument must be VARCHAR or VARCHAR[]");
	}

	for (const auto &path : file_paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_mzml_chromatograms: stdin is not supported (mzML requires file seeking)");
		}
	}

	for (const auto &path : file_paths) {
		if (!fs.FileExists(path)) {
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

unique_ptr<GlobalTableFunctionState> ReadMzMLChromatogramsTableFunction::InitGlobal(ClientContext &context,
                                                                                    TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	return duckdb::make_uniq<GlobalState>(data.file_paths);
}

unique_ptr<LocalTableFunctionState>
ReadMzMLChromatogramsTableFunction::InitLocal(ExecutionContext &context, TableFunctionInitInput &input,
                                              GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadMzMLChromatogramsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p,
                                                 DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::MzMLChromatogramBatch batch;
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

		// Safe without lock: each thread claims an exclusive file index via next_file_idx++
		// under the lock above, so no two threads access the same reader slot.
		if (!global_state.readers[local_state.current_file_idx]) {
			global_state.readers[local_state.current_file_idx] =
			    std::make_unique<miint::MzMLReader>(global_state.filepaths[local_state.current_file_idx]);
		}

		batch = global_state.readers[local_state.current_file_idx]->read_chromatograms(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		break;
	}

	size_t col = 0;

	// chromatogram_index (INTEGER, NOT NULL)
	SetResultVectorInt32(output.data[col++], batch.chromatogram_index);

	// chromatogram_id (VARCHAR, NOT NULL)
	SetResultVectorString(output.data[col++], batch.chromatogram_id);

	// chromatogram_type (VARCHAR, nullable — empty string maps to NULL)
	SetResultVectorStringNullable(output.data[col++], batch.chromatogram_type);

	// precursor_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.precursor_mz, batch.precursor_mz_valid);

	// product_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.product_mz, batch.product_mz_valid);

	// time_array (LIST(DOUBLE))
	SetResultVectorListDouble(output.data[col++], batch.time_array);

	// intensity_array (LIST(DOUBLE))
	SetResultVectorListDouble(output.data[col++], batch.intensity_array);

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[col++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

TableFunction ReadMzMLChromatogramsTableFunction::GetFunction() {
	auto tf = TableFunction("read_mzml_chromatograms", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadMzMLChromatogramsTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
} // namespace duckdb
