#include "read_alignments.hpp"
#include "reference_table_reader.hpp"
#include "table_function_common.hpp"
#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "duckdb/catalog/catalog.hpp"
#include "duckdb/catalog/catalog_entry/table_catalog_entry.hpp"
#include "duckdb/catalog/catalog_entry/view_catalog_entry.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"

// Helper function to detect stdin paths
static bool IsStdinPath(const std::string &path) {
	return path == "-" || path == "/dev/stdin";
}

namespace duckdb {

unique_ptr<FunctionData> ReadAlignmentsTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                           vector<LogicalType> &return_types,
                                                           vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> sam_paths;

	// Handle VARCHAR (single path, potentially a glob) or VARCHAR[] (array of literal paths)
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		// Single string - could be a glob pattern
		sam_paths = ExpandGlobPattern(fs, context, input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		// Array of strings - literal paths only (no glob expansion)
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			sam_paths.push_back(child.ToString());
		}
		if (sam_paths.empty()) {
			throw InvalidInputException("read_alignments: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_alignments: first argument must be VARCHAR or VARCHAR[]");
	}

	// Check if any path is stdin
	bool has_stdin = false;
	for (const auto &path : sam_paths) {
		if (IsStdinPath(path)) {
			has_stdin = true;
			break;
		}
	}

	// Error if stdin is used with multiple files
	if (has_stdin && sam_paths.size() > 1) {
		throw InvalidInputException("Cannot use stdin (/dev/stdin or -) with multiple files");
	}

	// Validate all files exist (skip stdin)
	for (const auto &path : sam_paths) {
		if (!IsStdinPath(path) && !fs.FileExists(path)) {
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

		// Validate table or view exists (use TABLE_ENTRY lookup which returns either)
		EntryLookupInfo lookup_info(CatalogType::TABLE_ENTRY, reference_lengths_table.value(), QueryErrorContext());
		auto entry = Catalog::GetEntry(context, INVALID_CATALOG, INVALID_SCHEMA, lookup_info, OnEntryNotFound::RETURN_NULL);
		if (!entry) {
			throw InvalidInputException("Table or view '%s' does not exist", reference_lengths_table.value());
		}
		if (entry->type != CatalogType::TABLE_ENTRY && entry->type != CatalogType::VIEW_ENTRY) {
			throw InvalidInputException("'%s' is not a table or view", reference_lengths_table.value());
		}
	}

	// Parse include_filepath parameter (optional BOOLEAN, default false)
	bool include_filepath = false;
	auto fp_param = input.named_parameters.find("include_filepath");
	if (fp_param != input.named_parameters.end() && !fp_param->second.IsNull()) {
		include_filepath = fp_param->second.GetValue<bool>();
	}

	// Parse include_seq_qual parameter (optional BOOLEAN, default false)
	bool include_seq_qual = false;
	auto seq_param = input.named_parameters.find("include_seq_qual");
	if (seq_param != input.named_parameters.end() && !seq_param->second.IsNull()) {
		include_seq_qual = seq_param->second.GetValue<bool>();
	}

	auto data = duckdb::make_uniq<Data>(sam_paths, reference_lengths_table, include_filepath, include_seq_qual);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return data;
}

unique_ptr<GlobalTableFunctionState> ReadAlignmentsTableFunction::InitGlobal(ClientContext &context,
                                                                             TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();

	// Read reference table if provided
	std::optional<std::unordered_map<std::string, uint64_t>> reference_lengths;
	if (data.reference_lengths_table.has_value()) {
		reference_lengths = ReadReferenceTable(context, data.reference_lengths_table.value());
	}

	// Check if reading from stdin (single file that is a stdin path)
	bool reading_from_stdin = (data.sam_paths.size() == 1 && IsStdinPath(data.sam_paths[0]));

	// For stdin, skip header checking to avoid consuming data
	// Validation will occur in SAMReader constructor when GlobalState is created
	// - If stdin has header without reference_lengths: SAMReader(path) succeeds
	// - If stdin is headerless without reference_lengths: SAMReader(path) throws "SAM file missing required header"
	// - If stdin is headerless with reference_lengths: SAMReader(path, refs) succeeds
	// - If stdin has header with reference_lengths: User error, cannot detect early (documented limitation)
	if (!reading_from_stdin) {
		// Check header consistency across all files
		bool first_file_has_header = false;
		for (size_t i = 0; i < data.sam_paths.size(); i++) {
			const auto &path = data.sam_paths[i];

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
				if (!has_header && !reference_lengths.has_value()) {
					throw IOException("File lacks a header, and no reference information provided");
				}

				// Validate first file: if has header, reject reference_lengths
				if (has_header && reference_lengths.has_value()) {
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
						throw IOException("Inconsistent headers across files: '" + data.sam_paths[0] + "' has header, '" + path +
						                  "' does not");
					} else {
						throw IOException("Inconsistent headers across files: '" + data.sam_paths[0] + "' lacks header, '" +
						                  path + "' has header");
					}
				}
			}
		}
	}

	auto gstate = duckdb::make_uniq<GlobalState>(data.sam_paths, reference_lengths, data.include_seq_qual);

	return gstate;
}

unique_ptr<LocalTableFunctionState> ReadAlignmentsTableFunction::InitLocal(ExecutionContext &context,
                                                                            TableFunctionInitInput &input,
                                                                            GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadAlignmentsTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SAMRecordBatch batch;
	std::string current_filepath;

	// Loop until we get data or run out of files
	while (true) {
		// If this thread doesn't have a file, claim one
		if (!local_state.has_file) {
			lock_guard<mutex> lock(global_state.lock);

			// Check if all files exhausted
			if (global_state.next_file_idx >= global_state.readers.size()) {
				output.SetCardinality(0);
				return;
			}

			// Claim next available file
			local_state.current_file_idx = global_state.next_file_idx;
			global_state.next_file_idx++;
			local_state.has_file = true;
		}

		// Read from claimed file (no lock needed - exclusive access)
		batch = global_state.readers[local_state.current_file_idx]->read(STANDARD_VECTOR_SIZE);
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		// If this file is exhausted, release it and try to claim another
		if (batch.empty()) {
			local_state.has_file = false;
			continue;
		}

		// Got data, break out of loop
		break;
	}

	// Set result vectors by directly copying from SOA batch
	size_t field_idx = 0;
	SetResultVectorString(output.data[field_idx++], batch.read_ids);
	SetResultVectorUInt16(output.data[field_idx++], batch.flags);
	SetResultVectorString(output.data[field_idx++], batch.references);
	SetResultVectorInt64(output.data[field_idx++], batch.positions);
	SetResultVectorInt64(output.data[field_idx++], batch.stop_positions);
	SetResultVectorUInt8(output.data[field_idx++], batch.mapqs);
	SetResultVectorString(output.data[field_idx++], batch.cigars);
	SetResultVectorString(output.data[field_idx++], batch.mate_references);
	SetResultVectorInt64(output.data[field_idx++], batch.mate_positions);
	SetResultVectorInt64(output.data[field_idx++], batch.template_lengths);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_as_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xs_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_ys_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xn_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xm_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xo_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_xg_values);
	SetResultVectorInt64Nullable(output.data[field_idx++], batch.tag_nm_values);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_yt_values);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_md_values);
	SetResultVectorStringNullable(output.data[field_idx++], batch.tag_sa_values);

	// Set SEQUENCE and QUAL columns if requested
	if (bind_data.include_seq_qual) {
		SetResultVectorString(output.data[field_idx++], batch.sequences);
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals);
	}

	// Set filepath column if requested
	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[field_idx++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

void ReadAlignmentsTableFunction::SetResultVectorString(Vector &result_vector,
                                                         const std::vector<std::string> &values) {
	auto result_data = FlatVector::GetData<string_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = StringVector::AddString(result_vector, values[j]);
	}
}

void ReadAlignmentsTableFunction::SetResultVectorStringNullable(Vector &result_vector,
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

void ReadAlignmentsTableFunction::SetResultVectorUInt8(Vector &result_vector, const std::vector<uint8_t> &values) {
	auto result_data = FlatVector::GetData<uint8_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void ReadAlignmentsTableFunction::SetResultVectorUInt16(Vector &result_vector, const std::vector<uint16_t> &values) {
	auto result_data = FlatVector::GetData<uint16_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void ReadAlignmentsTableFunction::SetResultVectorInt64(Vector &result_vector, const std::vector<int64_t> &values) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
	}
}

void ReadAlignmentsTableFunction::SetResultVectorInt64Nullable(Vector &result_vector,
                                                                const std::vector<int64_t> &values) {
	auto result_data = FlatVector::GetData<int64_t>(result_vector);
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllInvalid(values.size());

	for (idx_t j = 0; j < values.size(); j++) {
		result_data[j] = values[j];
		// Tags return -1 when not present
		if (values[j] != -1) {
			validity.SetValid(j);
		}
	}
}

void ReadAlignmentsTableFunction::SetResultVectorFilepath(Vector &result_vector, const std::string &filepath) {
	result_vector.SetVectorType(VectorType::CONSTANT_VECTOR);
	auto result_data = ConstantVector::GetData<string_t>(result_vector);
	*result_data = StringVector::AddString(result_vector, filepath);
}

void ReadAlignmentsTableFunction::SetResultVectorListUInt8(Vector &result_vector,
                                                             const std::vector<miint::QualScore> &values) {
	// Compute total number of elements
	idx_t total_child_elements = 0;
	for (auto &qual : values) {
		total_child_elements += qual.as_string().length();
	}

	ListVector::Reserve(result_vector, total_child_elements);
	ListVector::SetListSize(result_vector, total_child_elements);

	auto &child_vector = ListVector::GetEntry(result_vector);
	auto child_data = FlatVector::GetData<uint8_t>(child_vector);
	auto list_entries = FlatVector::GetData<list_entry_t>(result_vector);

	const auto output_count = values.size();
	idx_t value_offset = 0;
	for (idx_t row_offset = 0; row_offset < output_count; row_offset++) {
		// Always use offset=33 since we store internally as Phred33 ASCII
		const auto qual_data = values[row_offset].as_vec(33);
		list_entries[row_offset].offset = value_offset;
		list_entries[row_offset].length = qual_data.size();

		for (auto &value : qual_data) {
			child_data[value_offset++] = value;
		}
	}

	// All entries are not null
	auto &validity = FlatVector::Validity(result_vector);
	validity.SetAllValid(output_count);

	auto &child_validity = FlatVector::Validity(child_vector);
	child_validity.SetAllValid(total_child_elements);
}

TableFunction ReadAlignmentsTableFunction::GetFunction() {
	auto tf = TableFunction("read_alignments", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["reference_lengths"] = LogicalType::ANY;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["include_seq_qual"] = LogicalType::BOOLEAN;
	return tf;
}

void ReadAlignmentsTableFunction::Register(ExtensionLoader &loader) {
	// Register the main function
	loader.RegisterFunction(GetFunction());

	// Register backward compatibility alias
	auto read_sam_alias = TableFunction("read_sam", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	read_sam_alias.named_parameters["reference_lengths"] = LogicalType::ANY;
	read_sam_alias.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	read_sam_alias.named_parameters["include_seq_qual"] = LogicalType::BOOLEAN;
	loader.RegisterFunction(read_sam_alias);
}

}; // namespace duckdb
