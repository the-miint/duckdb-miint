#include "SequenceReader.hpp"
#include "SequenceRecord.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/storage/statistics/node_statistics.hpp"
#include <filesystem>
#include <read_fastx.hpp>

namespace duckdb {

unique_ptr<FunctionData> ReadFastxTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                      vector<duckdb::LogicalType> &return_types,
                                                      vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> sequence1_paths;
	bool sequence1_is_glob = false;

	// Helper to detect stdin paths
	auto is_stdin = [](const std::string &path) {
		return path == "-" || path == "/dev/stdin" || path == "/dev/fd/0" || path == "/proc/self/fd/0";
	};

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_fastx: first argument cannot be NULL");
	}
	// Handle VARCHAR (single path, potentially a glob) or VARCHAR[] (array of literal paths)
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		// Single string - could be a glob pattern
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		sequence1_paths = std::move(result.paths);
		sequence1_is_glob = result.is_glob;
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		// Array of strings - literal paths only (no glob expansion)
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_fastx: file path list cannot contain NULL");
			}
			sequence1_paths.push_back(child.ToString());
		}
		if (sequence1_paths.empty()) {
			throw InvalidInputException("read_fastx: at least one file path must be provided");
		}
	} else {
		throw InvalidInputException("read_fastx: first argument must be VARCHAR or VARCHAR[]");
	}

	// Detect stdin usage
	// Note: We detect common stdin paths, but cannot catch all cases (e.g., symlinks to stdin,
	// /proc/self/fd/0 on Linux). If parallel threads attempt to read from an undetected stdin,
	// the behavior is undefined and may result in data corruption or crashes.
	bool uses_stdin = false;

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
				throw InvalidInputException(
				    "stdin ('-' or '/dev/stdin') must be a single file path, not part of an array");
			}
		}
	}

	// Validate all sequence1 files exist (skip stdin and remote paths)
	for (const auto &path : sequence1_paths) {
		if (!is_stdin(path) && !miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
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
		bool sequence2_is_glob = false;

		if (path2_value.type().id() == LogicalTypeId::VARCHAR) {
			// Single string - could be a glob pattern
			auto result = ExpandGlobPatternWithInfo(fs, context, path2_value.ToString());
			seq2_paths = std::move(result.paths);
			sequence2_is_glob = result.is_glob;
		} else if (path2_value.type().id() == LogicalTypeId::LIST) {
			// Array of strings - literal paths only
			auto &list_children = ListValue::GetChildren(path2_value);
			for (const auto &child : list_children) {
				if (child.IsNull()) {
					throw InvalidInputException("read_fastx: sequence2 list cannot contain NULL");
				}
				seq2_paths.push_back(child.ToString());
			}
		} else {
			throw InvalidInputException("sequence2 parameter must be VARCHAR or VARCHAR[]");
		}

		// For paired-end with globs: both must be globs or both must be literals
		if (sequence1_is_glob != sequence2_is_glob) {
			throw InvalidInputException(
			    "When using glob patterns, both sequence1 and sequence2 must use glob patterns");
		}

		// Validate file counts match
		if (seq2_paths.size() != sequence1_paths.size()) {
			if (sequence1_is_glob && sequence2_is_glob) {
				throw InvalidInputException("Glob patterns matched different number of files: sequence1 matched " +
				                            std::to_string(sequence1_paths.size()) + " files, sequence2 matched " +
				                            std::to_string(seq2_paths.size()) + " files");
			} else {
				throw InvalidInputException("Mismatched array lengths: sequence1 has " +
				                            std::to_string(sequence1_paths.size()) + " files, sequence2 has " +
				                            std::to_string(seq2_paths.size()) + " files");
			}
		}

		// Validate all sequence2 files exist (skip remote paths)
		for (const auto &path : seq2_paths) {
			if (!miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
				throw IOException("File not found: " + path);
			}
		}

		sequence2_paths = seq2_paths;
	}

	// Interleaved paired-end from a single stream (record 2k-1 = R1, 2k = R2). Mutually
	// exclusive with the two-file sequence2 mode.
	bool interleaved = false;
	auto il_param = input.named_parameters.find("interleaved");
	if (il_param != input.named_parameters.end() && !il_param->second.IsNull()) {
		interleaved = il_param->second.GetValue<bool>();
	}
	if (interleaved && sequence2_paths.has_value()) {
		throw InvalidInputException("read_fastx: interleaved cannot be combined with the sequence2 parameter");
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

	// Default byte budget: 512 MiB. Large enough that typical short-read workflows
	// always fill STANDARD_VECTOR_SIZE (~600 KB per chunk) without being capped, but
	// small enough to prevent multi-GB chunks when reading large sequences
	// (e.g. bacterial genomes at ~5 MB/record would otherwise blow past RAM at 2048x).
	uint64_t max_batch_bytes = 512ULL * 1024 * 1024;
	auto mbb_param = input.named_parameters.find("max_batch_bytes");
	if (mbb_param != input.named_parameters.end() && !mbb_param->second.IsNull()) {
		std::string raw = mbb_param->second.ToString();
		idx_t parsed = 0;
		std::string err = StringUtil::TryParseFormattedBytes(raw, parsed);
		if (!err.empty()) {
			throw InvalidInputException("max_batch_bytes: %s", err);
		}
		if (parsed == 0) {
			throw InvalidInputException("max_batch_bytes must be greater than 0");
		}
		max_batch_bytes = static_cast<uint64_t>(parsed);
	}

	auto data = duckdb::make_uniq<Data>(sequence1_paths, sequence2_paths, include_filepath, uses_stdin, interleaved,
	                                    qual_offset, max_batch_bytes);
	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}
	return data;
}

unique_ptr<GlobalTableFunctionState> ReadFastxTableFunction::InitGlobal(ClientContext &context,
                                                                        TableFunctionInitInput &input) {
	auto &data = input.bind_data->Cast<Data>();
	auto &fs = FileSystem::GetFileSystem(context);
	return duckdb::make_uniq<GlobalState>(fs, data.sequence1_paths, data.sequence2_paths, data.uses_stdin,
	                                      data.interleaved);
}

unique_ptr<LocalTableFunctionState> ReadFastxTableFunction::InitLocal(ExecutionContext &context,
                                                                      TableFunctionInitInput &input,
                                                                      GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadFastxTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	miint::SequenceRecordBatch batch;
	std::string current_filepath;
	uint64_t start_sequence_index;

	// Loop until we get data or run out of files
	while (true) {
		// If this thread doesn't have a file, claim one
		if (!local_state.has_file) {
			{
				lock_guard<mutex> read_lock(global_state.lock);

				// Check if all files exhausted
				if (global_state.next_file_idx >= global_state.sequence1_filepaths.size()) {
					output.SetCardinality(0);
					return;
				}

				// Claim next available file index
				local_state.current_file_idx = global_state.next_file_idx;
				global_state.next_file_idx++;
				local_state.has_file = true;
			} // Lock released before I/O

			// Open reader without holding the lock — safe because this thread
			// has exclusive ownership of the claimed file index
			global_state.OpenReader(local_state.current_file_idx);
		}

		// Read from claimed file (no lock needed)
		batch =
		    global_state.readers[local_state.current_file_idx]->read(STANDARD_VECTOR_SIZE, bind_data.max_batch_bytes);
		current_filepath = global_state.sequence1_filepaths[local_state.current_file_idx];

		// If this file is exhausted, release it and try to claim another
		if (batch.empty()) {
			local_state.has_file = false;
			continue; // Loop to claim next file
		}

		// Got data, break out of loop
		break;
	}

	// Get sequence indices for this chunk from the current file's counter
	// No atomic operation needed - this thread has exclusive access to this file
	start_sequence_index = global_state.file_sequence_counters[local_state.current_file_idx];
	global_state.file_sequence_counters[local_state.current_file_idx] += batch.size();

	// Set sequence_index column (first column, index 0)
	auto &sequence_index_vector = output.data[0];
	auto sequence_index_data = FlatVector::GetData<int64_t>(sequence_index_vector);
	for (idx_t j = 0; j < batch.size(); j++) {
		sequence_index_data[j] = static_cast<int64_t>(start_sequence_index + j);
	}

	// Set result vectors for sequence fields (now starting at index 1)
	size_t field_idx = 1;
	SetResultVectorString(output.data[field_idx++], batch.read_ids);
	SetResultVectorStringNullable(output.data[field_idx++], batch.comments);
	SetResultVectorString(output.data[field_idx++], batch.sequences1);

	// SEQUENCE2 - nullable for unpaired
	if (!batch.is_paired) {
		SetResultVectorNull(output.data[field_idx++]);
	} else {
		SetResultVectorString(output.data[field_idx++], batch.sequences2);
	}

	// QUAL1 - null for FASTA (all records in a file are guaranteed same format)
	if (batch.quals1[0].as_string().empty()) {
		SetResultVectorNull(output.data[field_idx++]);
	} else {
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals1, bind_data.qual_offset);
	}

	// QUAL2 - null for unpaired or FASTA
	if (!batch.is_paired || batch.quals2.empty() || batch.quals2[0].as_string().empty()) {
		SetResultVectorNull(output.data[field_idx++]);
	} else {
		SetResultVectorListUInt8(output.data[field_idx++], batch.quals2, bind_data.qual_offset);
	}

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[field_idx++], current_filepath);
	}

	output.SetCardinality(batch.size());
}

// Cardinality estimate so the query optimizer knows roughly how big the scan is. Without
// this, DuckDB defaults unknown table-function scans to a tiny cardinality, which causes
// joins to pick read_fastx as the BUILD side -- catastrophic when the file is multi-GB of
// genomic sequence (e.g. bacterial genome catalogs): DuckDB ends up hashing ~50 GB of
// sequence into memory instead of streaming.
//
// We deliberately err on the side of OVER-estimating record count:
//   - if the estimate is too high, the optimizer treats the scan as "large" and probes it
//     -- always the safe choice for a streaming file source.
//   - if it's too low (the status quo), the optimizer builds a hash table on it -> OOM.
//
// Heuristic: bytes-on-disk / AVG_RECORD_BYTES. 50 bytes/record is intentionally small so
// that bacterial-genome FASTA (~5 MB/record) and short-read FASTQ (~300 B/record) both
// resolve to "big" when the file is big. Gzipped files use a 4x expansion factor.
// Remote files and stdin fall back to a large constant (exact record count is unknown but
// streaming is always the right plan).
static unique_ptr<NodeStatistics> ReadFastxCardinality(ClientContext &context, const FunctionData *bind_data_p) {
	constexpr idx_t AVG_BYTES_PER_RECORD = 50;
	constexpr idx_t GZIP_EXPANSION_FACTOR = 4;
	constexpr idx_t REMOTE_OR_STDIN_ESTIMATE = 10'000'000; // err large

	auto &bind_data = bind_data_p->Cast<ReadFastxTableFunction::Data>();

	if (bind_data.uses_stdin) {
		return make_uniq<NodeStatistics>(REMOTE_OR_STDIN_ESTIMATE);
	}

	idx_t total_bytes = 0;
	auto accumulate = [&](const std::string &path) -> bool {
		if (miint::RemoteFileHelper::IsRemotePath(path)) {
			return false; // signal "unknown, fall back"
		}
		std::error_code ec;
		auto sz = std::filesystem::file_size(path, ec);
		if (ec) {
			return false;
		}
		if (IsGzipped(path)) {
			sz *= GZIP_EXPANSION_FACTOR;
		}
		total_bytes += sz;
		return true;
	};

	for (const auto &path : bind_data.sequence1_paths) {
		if (!accumulate(path)) {
			return make_uniq<NodeStatistics>(REMOTE_OR_STDIN_ESTIMATE);
		}
	}
	if (bind_data.sequence2_paths.has_value()) {
		for (const auto &path : bind_data.sequence2_paths.value()) {
			if (!accumulate(path)) {
				return make_uniq<NodeStatistics>(REMOTE_OR_STDIN_ESTIMATE);
			}
		}
	}

	idx_t estimated = total_bytes / AVG_BYTES_PER_RECORD;
	if (estimated < 1) {
		estimated = 1;
	}
	return make_uniq<NodeStatistics>(estimated);
}

TableFunction ReadFastxTableFunction::GetFunction() {
	auto tf = TableFunction("read_fastx", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.named_parameters["sequence2"] = LogicalType::ANY;
	tf.named_parameters["interleaved"] = LogicalType::BOOLEAN;
	tf.named_parameters["include_filepath"] = LogicalType::BOOLEAN;
	tf.named_parameters["qual_offset"] = LogicalType::BIGINT;
	tf.named_parameters["max_batch_bytes"] = LogicalType::VARCHAR;
	tf.cardinality = ReadFastxCardinality;
	return tf;
}

void ReadFastxTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
}; // namespace duckdb
