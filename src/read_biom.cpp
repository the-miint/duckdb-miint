#include "read_biom.hpp"
#include "BIOMReader.hpp"
#include "duckdb.h"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

unique_ptr<FunctionData> ReadBIOMTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                     vector<LogicalType> &return_types, vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> biom_paths;

	// Handle VARCHAR or VARCHAR[] input for file paths
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		biom_paths.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			biom_paths.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_biom: first argument must be VARCHAR or VARCHAR[]");
	}

	// Validate that at least one file was provided
	if (biom_paths.empty()) {
		throw InvalidInputException("read_biom: first argument must be VARCHAR or VARCHAR[]");
	}

	// Validate all files exist and are BIOM
	for (const auto &path : biom_paths) {
		if (!fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}

		if (!miint::BIOMReader::IsBIOM(path)) {
			throw IOException("File is not a BIOM file: " + path);
		}
	}

	// Parse include_filepath parameter (optional BOOLEAN, default false)
	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	auto data = duckdb::make_uniq<Data>(biom_paths, include_filepath);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadBIOMTableFunction::InitGlobal(ClientContext &context,
                                                                       TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();

	auto gstate = duckdb::make_uniq<GlobalState>(data.biom_paths);

	return std::move(gstate);
}

unique_ptr<LocalTableFunctionState> ReadBIOMTableFunction::InitLocal(ExecutionContext &context,
                                                                     TableFunctionInitInput &input,
                                                                     GlobalTableFunctionState *global_state) {
	auto local_state = make_uniq<LocalState>();
	auto &gstate = global_state->Cast<GlobalState>();

	// Each thread grabs its first file
	local_state->GetNextFile(gstate);

	return local_state;
}

void ReadBIOMTableFunction::SetResultVector(Vector &result_vector, const miint::BIOMTableField &field,
                                            const miint::BIOMTable &record, const size_t &current_row,
                                            const size_t &n_rows) {
	switch (field) {
	case miint::BIOMTableField::SAMPLE_ID:
	case miint::BIOMTableField::FEATURE_ID:
		SetResultVectorString(result_vector, field, record, current_row, n_rows);
		break;
	case miint::BIOMTableField::VALUE:
		SetResultVectorDouble(result_vector, field, record, current_row, n_rows);
		break;
	default:
		throw NotImplementedException("field not supported");
	}
}

void ReadBIOMTableFunction::SetResultVectorString(Vector &result_vector, const miint::BIOMTableField &field,
                                                  const miint::BIOMTable &record, const size_t &current_row,
                                                  const size_t &n_rows) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);

	// Use const reference to avoid copying the entire string vector
	const auto &current_names = (field == miint::BIOMTableField::SAMPLE_ID)
	                            ? record.COOSamples()
	                            : record.COOFeatures();

	for (size_t i = 0; i < n_rows; i++) {
		result_data[i] = StringVector::AddString(result_vector, current_names[current_row + i]);
	}
}

void ReadBIOMTableFunction::SetResultVectorDouble(Vector &result_vector, const miint::BIOMTableField &field,
                                                  const miint::BIOMTable &record, const size_t &current_row,
                                                  const size_t &n_rows) {
	auto result_data = FlatVector::GetData<double>(result_vector);
	auto &data = record.COOValues();

	for (size_t i = 0; i < n_rows; i++) {
		result_data[i] = data[current_row + i];
	}
}

void ReadBIOMTableFunction::SetResultVectorFilepath(Vector &result_vector, const std::string &filepath,
                                                    size_t num_records) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void ReadBIOMTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	// Check if we need to load the next file
	while (local_state.current_row >= local_state.total_rows) {
		// Current file is exhausted, try to get next file
		bool got_file = local_state.GetNextFile(global_state);
		if (!got_file) {
			// No more files to process
			output.SetCardinality(0);
			return;
		}
		// GetNextFile resets current_row to 0 and loads new table
	}

	// Calculate how many rows we can emit in this chunk
	size_t n_rows =
	    std::min(local_state.current_row + STANDARD_VECTOR_SIZE, local_state.total_rows) - local_state.current_row;

	// Populate output vectors
	for (size_t i = 0; i < bind_data.fields.size(); i++) {
		auto &result_vector = output.data[i];
		auto &field = bind_data.fields[i];
		SetResultVector(result_vector, field, local_state.table, local_state.current_row, n_rows);
	}

	if (bind_data.include_filepath) {
		auto &result_vector = output.data[bind_data.fields.size()];
		SetResultVectorFilepath(result_vector, local_state.path, n_rows);
	}

	output.SetCardinality(n_rows);
	local_state.current_row += n_rows;
}

TableFunction ReadBIOMTableFunction::GetFunction() {
	auto tf = TableFunction("read_biom", {LogicalType::ANY}, Execute, Bind, InitGlobal);
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.init_local = InitLocal;
	return tf;
}

void ReadBIOMTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
}; // namespace duckdb
