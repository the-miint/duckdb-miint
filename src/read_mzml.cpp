#include "MzMLReader.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_mzml.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadMzMLTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
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
			throw InvalidInputException("read_mzml: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_mzml: first argument must be VARCHAR or VARCHAR[]");
	}

	for (const auto &path : file_paths) {
		if (IsStdinPath(path)) {
			throw InvalidInputException("read_mzml: stdin is not supported (mzML requires file seeking)");
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

unique_ptr<GlobalTableFunctionState> ReadMzMLTableFunction::InitGlobal(ClientContext &context,
                                                                       TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	return duckdb::make_uniq<GlobalState>(data.file_paths);
}

unique_ptr<LocalTableFunctionState> ReadMzMLTableFunction::InitLocal(ExecutionContext &context,
                                                                     TableFunctionInitInput &input,
                                                                     GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadMzMLTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
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

		// Safe without lock: each thread claims an exclusive file index via next_file_idx++
		// under the lock above, so no two threads access the same reader slot.
		if (!global_state.readers[local_state.current_file_idx]) {
			global_state.readers[local_state.current_file_idx] =
			    std::make_unique<miint::MzMLReader>(global_state.filepaths[local_state.current_file_idx]);
		}

		batch = global_state.readers[local_state.current_file_idx]->read_spectra(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		break;
	}

	size_t col = 0;

	// spectrum_index (INTEGER, NOT NULL)
	SetResultVectorInt32(output.data[col++], batch.spectrum_index);

	// spectrum_id (VARCHAR, NOT NULL)
	SetResultVectorString(output.data[col++], batch.spectrum_id);

	// ms_level (INTEGER, NOT NULL)
	SetResultVectorInt32(output.data[col++], batch.ms_level);

	// retention_time (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.retention_time, batch.retention_time_valid);

	// spectrum_type (VARCHAR, nullable — empty string maps to NULL)
	SetResultVectorStringNullable(output.data[col++], batch.spectrum_type);

	// polarity (VARCHAR, nullable — empty string maps to NULL)
	SetResultVectorStringNullable(output.data[col++], batch.polarity);

	// base_peak_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.base_peak_mz, batch.base_peak_mz_valid);

	// base_peak_intensity (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.base_peak_intensity, batch.base_peak_intensity_valid);

	// total_ion_current (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.total_ion_current, batch.total_ion_current_valid);

	// lowest_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.lowest_mz, batch.lowest_mz_valid);

	// highest_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.highest_mz, batch.highest_mz_valid);

	// default_array_length (INTEGER, NOT NULL)
	SetResultVectorInt32(output.data[col++], batch.default_array_length);

	// precursor_mz (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.precursor_mz, batch.precursor_mz_valid);

	// precursor_charge (INTEGER, nullable)
	SetResultVectorInt32Nullable(output.data[col++], batch.precursor_charge, batch.precursor_charge_valid);

	// precursor_intensity (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.precursor_intensity, batch.precursor_intensity_valid);

	// isolation_window_target (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_target,
	                              batch.isolation_window_target_valid);

	// isolation_window_lower (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_lower, batch.isolation_window_lower_valid);

	// isolation_window_upper (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.isolation_window_upper, batch.isolation_window_upper_valid);

	// activation_method (VARCHAR, nullable — empty string maps to NULL)
	SetResultVectorStringNullable(output.data[col++], batch.activation_method);

	// collision_energy (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.collision_energy, batch.collision_energy_valid);

	// mz_array (LIST(DOUBLE))
	SetResultVectorListDouble(output.data[col++], batch.mz_array);

	// intensity_array (LIST(DOUBLE))
	SetResultVectorListDouble(output.data[col++], batch.intensity_array);

	// filter_string (VARCHAR, nullable — empty string maps to NULL)
	SetResultVectorStringNullable(output.data[col++], batch.filter_string);

	// scan_window_lower (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.scan_window_lower, batch.scan_window_lower_valid);

	// scan_window_upper (DOUBLE, nullable)
	SetResultVectorDoubleNullable(output.data[col++], batch.scan_window_upper, batch.scan_window_upper_valid);

	// ms1_scan_index (INTEGER, nullable)
	SetResultVectorInt32Nullable(output.data[col++], batch.ms1_scan_index, batch.ms1_scan_index_valid);

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[col++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

TableFunction ReadMzMLTableFunction::GetFunction() {
	auto tf = TableFunction("read_mzml", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadMzMLTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
} // namespace duckdb
