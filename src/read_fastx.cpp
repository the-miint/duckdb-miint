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

	auto path1 = input.inputs[0].ToString();
	if (!fs.FileExists(path1)) {
		throw IOException("File not found: " + path1);
	}

	auto path2 = input.named_parameters["reverse"];
	auto path2_value = path2.IsNull() ? "" : path2.GetValue<std::string>();
	if (!path2.IsNull()) {
		if (!fs.FileExists(path2_value)) {
			throw IOException("File not found: " + path2_value);
		}
	}

	auto data = duckdb::make_uniq<Data>(path1, path2_value);
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
	auto gstate = duckdb::make_uniq<GlobalState>(data.read1_path, data.read2_path);

	return std::move(gstate);
}

void ReadFastxTableFunction::SetResultVector(Vector &result_vector, const miint::SequenceRecordField &field,
                                             const std::vector<miint::SequenceRecord> &records) {
	auto is_paired = records[0].is_paired;
	if (!is_paired && (field == miint::SequenceRecordField::READ2 || field == miint::SequenceRecordField::QUAL2)) {
		SetResultVectorNull(result_vector);
	} else if (field == miint::SequenceRecordField::READ_ID || field == miint::SequenceRecordField::READ1 ||
	           field == miint::SequenceRecordField::READ2) {
		SetResultVectorString(result_vector, field, records);
	} else if (field == miint::SequenceRecordField::COMMENT) {
		SetResultVectorStringNullable(result_vector, field, records);
	} else if (field == miint::SequenceRecordField::QUAL1 || field == miint::SequenceRecordField::QUAL2) {
		SetResultVectorListUInt8(result_vector, field, records);
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
                                                      const std::vector<miint::SequenceRecord> &records) {
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
		const auto qual_data = records[row_offset].GetQual(field).as_vec();
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

void ReadFastxTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	auto records = std::vector<miint::SequenceRecord>();
	{
		lock_guard<mutex> read_lock(global_state.lock);
		records = global_state.reader.read(STANDARD_VECTOR_SIZE);
	}

	if (!records.size()) {
		output.SetCardinality(0);
		return; // No more work
	}

	for (idx_t i = 0; i < bind_data.fields.size(); i++) {
		auto &result_vector = output.data[i];
		auto &field = bind_data.fields[i];
		SetResultVector(result_vector, field, records);
	}

	output.SetCardinality(records.size());
}

TableFunction ReadFastxTableFunction::GetFunction() {
	// IMPORTANT: "set preserve_insertion_order=false" to run in parallel
	// otherwise it will default to serial for many (not all) operations
	auto tf = TableFunction("read_fastx", {LogicalType::VARCHAR}, Execute, Bind, InitGlobal);
	tf.named_parameters["reverse"] = LogicalType::VARCHAR;
	return tf;
}
}; // namespace duckdb
