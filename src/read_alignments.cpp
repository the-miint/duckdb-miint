#include "read_alignments.hpp"
#include "reference_table_reader.hpp"
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"

namespace duckdb {

unique_ptr<FunctionData> ReadAlignmentsTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                           vector<LogicalType> &return_types,
                                                           vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> sam_paths;

	// Handle VARCHAR or VARCHAR[] input for file paths
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		sam_paths.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			sam_paths.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_alignments: first argument must be VARCHAR or VARCHAR[]");
	}

	// Validate all files exist
	for (const auto &path : sam_paths) {
		if (!fs.FileExists(path)) {
			throw IOException("File not found: " + path);
		}
	}

	// Parse reference_lengths parameter (optional VARCHAR - table name)
	std::optional<std::string> reference_lengths_table;
	auto ref_param = input.named_parameters.find("reference_lengths");
	if (ref_param != input.named_parameters.end() && !ref_param->second.IsNull()) {
		const auto &table_value = ref_param->second;

		if (table_value.type().id() != LogicalTypeId::VARCHAR) {
			throw InvalidInputException("reference_lengths must be a VARCHAR (table name)");
		}

		reference_lengths_table = table_value.ToString();

		// Validate table exists
		auto catalog_entry = Catalog::GetEntry<TableCatalogEntry>(
		    context, INVALID_CATALOG, INVALID_SCHEMA, reference_lengths_table.value(), OnEntryNotFound::RETURN_NULL);
		if (!catalog_entry) {
			throw InvalidInputException("Table '%s' does not exist", reference_lengths_table.value());
		}
	}

	// Parse include_filepath parameter (optional BOOLEAN, default false)
	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	// Check header consistency across all files
	bool first_file_has_header = false;
	for (size_t i = 0; i < sam_paths.size(); i++) {
		const auto &path = sam_paths[i];

		// Check if file has header by attempting to read it
		miint::SAMFilePtr test_fp(sam_open(path.c_str(), "r"));
		if (!test_fp) {
			throw IOException("Failed to open SAM file: " + path);
		}
		miint::SAMHeaderPtr test_hdr(sam_hdr_read(test_fp.get()));
		bool has_header = (test_hdr && test_hdr->n_targets > 0);

		if (i == 0) {
			first_file_has_header = has_header;

			// Validate first file: if no header, require reference_lengths
			if (!has_header && !reference_lengths_table.has_value()) {
				throw IOException("File lacks a header, and no reference information provided");
			}

			// Validate first file: if has header, reject reference_lengths
			if (has_header && reference_lengths_table.has_value()) {
				// Get actual format type from HTSlib
				const htsFormat *fmt = hts_get_format(test_fp.get());
				const char *format_name = (fmt && fmt->format == bam) ? "BAM" : "SAM";
				throw IOException(std::string(format_name) +
				                  " file has header, but reference_lengths parameter was provided");
			}
		} else {
			// Validate subsequent files have same header status
			if (has_header != first_file_has_header) {
				if (first_file_has_header) {
					throw IOException("Inconsistent headers across files: '" + sam_paths[0] + "' has header, '" + path +
					                  "' does not");
				} else {
					throw IOException("Inconsistent headers across files: '" + sam_paths[0] + "' lacks header, '" +
					                  path + "' has header");
				}
			}
		}
	}

	// If headerless, validate that all references in SAM are provided
	// This is deferred to InitGlobal since we need to read records to see which references are used

	auto data = duckdb::make_uniq<Data>(sam_paths, reference_lengths_table, include_filepath);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ReadAlignmentsTableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();

	// Read reference table if provided
	std::optional<std::unordered_map<std::string, uint64_t>> reference_lengths;
	if (data.reference_lengths_table.has_value()) {
		reference_lengths = ReadReferenceTable(context, data.reference_lengths_table.value());
	}

	auto gstate = duckdb::make_uniq<GlobalState>(data.sam_paths, reference_lengths);

	return std::move(gstate);
}

void ReadAlignmentsTableFunction::SetResultVector(Vector &result_vector, const miint::SAMRecordField &field,
                                                  const std::vector<miint::SAMRecord> &records) {
	switch (field) {
	case miint::SAMRecordField::READ_ID:
	case miint::SAMRecordField::REFERENCE:
	case miint::SAMRecordField::CIGAR:
	case miint::SAMRecordField::MATE_REFERENCE:
		SetResultVectorString(result_vector, field, records);
		break;
	case miint::SAMRecordField::TAG_YT:
	case miint::SAMRecordField::TAG_MD:
	case miint::SAMRecordField::TAG_SA:
		SetResultVectorStringNullable(result_vector, field, records);
		break;
	case miint::SAMRecordField::MAPQ:
		SetResultVectorUInt8(result_vector, field, records);
		break;
	case miint::SAMRecordField::FLAGS:
		SetResultVectorUInt16(result_vector, field, records);
		break;
	case miint::SAMRecordField::POSITION:
	case miint::SAMRecordField::STOP_POSITION:
	case miint::SAMRecordField::MATE_POSITION:
	case miint::SAMRecordField::TEMPLATE_LENGTH:
		SetResultVectorInt64(result_vector, field, records);
		break;
	case miint::SAMRecordField::TAG_AS:
	case miint::SAMRecordField::TAG_XS:
	case miint::SAMRecordField::TAG_YS:
	case miint::SAMRecordField::TAG_XN:
	case miint::SAMRecordField::TAG_XM:
	case miint::SAMRecordField::TAG_XO:
	case miint::SAMRecordField::TAG_XG:
	case miint::SAMRecordField::TAG_NM:
		SetResultVectorInt64Nullable(result_vector, field, records);
		break;
	default:
		throw NotImplementedException("field not supported");
	}
}

void ReadAlignmentsTableFunction::SetResultVectorString(Vector &result_vector, const miint::SAMRecordField &field,
                                                        const std::vector<miint::SAMRecord> &records) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, records[j].GetString(field));
	}
}

void ReadAlignmentsTableFunction::SetResultVectorStringNullable(Vector &result_vector,
                                                                const miint::SAMRecordField &field,
                                                                const std::vector<miint::SAMRecord> &records) {
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

void ReadAlignmentsTableFunction::SetResultVectorUInt8(Vector &result_vector, const miint::SAMRecordField &field,
                                                       const std::vector<miint::SAMRecord> &records) {
	auto result_data = FlatVector::GetData<uint8_t>(result_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		result_data[j] = records[j].GetUInt8(field);
	}
}

void ReadAlignmentsTableFunction::SetResultVectorUInt16(Vector &result_vector, const miint::SAMRecordField &field,
                                                        const std::vector<miint::SAMRecord> &records) {
	auto result_data = FlatVector::GetData<uint16_t>(result_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		result_data[j] = records[j].GetUInt16(field);
	}
}

void ReadAlignmentsTableFunction::SetResultVectorInt64(Vector &result_vector, const miint::SAMRecordField &field,
                                                       const std::vector<miint::SAMRecord> &records) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	for (idx_t j = 0; j < records.size(); j++) {
		result_data[j] = records[j].GetInt64(field);
	}
}

void ReadAlignmentsTableFunction::SetResultVectorInt64Nullable(Vector &result_vector,
                                                               const miint::SAMRecordField &field,
                                                               const std::vector<miint::SAMRecord> &records) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(records.size());

	for (idx_t j = 0; j < records.size(); j++) {
		const auto value = records[j].GetInt64(field);
		result_data[j] = value;

		// Tags return -1 when not present
		if (value != -1) {
			validity.SetValid(j);
		}
	}
}

void ReadAlignmentsTableFunction::SetResultVectorFilepath(Vector &result_vector, const std::string &filepath,
                                                          size_t num_records) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void ReadAlignmentsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	std::vector<miint::SAMRecord> records;
	std::string current_filepath;

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
				current_filepath = global_state.filepaths[global_state.current_file_idx];
				break;
			}

			// Current file exhausted, move to next
			global_state.current_file_idx++;
		}

		// All files exhausted
		if (records.empty()) {
			global_state.finished = true;
			output.SetCardinality(0);
			return;
		}
	}

	// Set result vectors (outside lock - this is CPU work)
	for (idx_t i = 0; i < bind_data.fields.size(); i++) {
		auto &result_vector = output.data[i];
		auto &field = bind_data.fields[i];
		SetResultVector(result_vector, field, records);
	}

	// Set filepath column if requested
	if (bind_data.include_filepath) {
		auto &result_vector = output.data[bind_data.fields.size()];
		SetResultVectorFilepath(result_vector, current_filepath, records.size());
	}

	output.SetCardinality(records.size());
}

TableFunction ReadAlignmentsTableFunction::GetFunction() {
	auto tf = TableFunction("read_alignments", {LogicalType::ANY}, Execute, Bind, InitGlobal);
	tf.named_parameters["reference_lengths"] = LogicalType::ANY;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadAlignmentsTableFunction::Register(ExtensionLoader &loader) {
	// Register the main function
	loader.RegisterFunction(GetFunction());

	// Register backward compatibility alias
	auto read_sam_alias = TableFunction("read_sam", {LogicalType::ANY}, Execute, Bind, InitGlobal);
	read_sam_alias.named_parameters["reference_lengths"] = LogicalType::ANY;
	read_sam_alias.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(read_sam_alias);
}

}; // namespace duckdb
