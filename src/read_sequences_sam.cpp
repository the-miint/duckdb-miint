#include "SAMReader.hpp"
#include "SAMRecord.hpp"
#include "remote_file_helper.hpp"
#include "table_function_common.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include <read_sequences_sam.hpp>
#include <cstring>
#include <stdexcept>

namespace duckdb {

unique_ptr<FunctionData> ReadSequencesSamTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                             vector<duckdb::LogicalType> &return_types,
                                                             vector<std::string> &names) {
	FileSystem &fs = FileSystem::GetFileSystem(context);

	std::vector<std::string> file_paths;

	if (input.inputs[0].IsNull()) {
		throw InvalidInputException("read_sequences_sam: first argument cannot be NULL");
	}
	// Handle VARCHAR (single path, potentially a glob) or VARCHAR[] (array of literal paths)
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		auto result = ExpandGlobPatternWithInfo(fs, context, input.inputs[0].ToString());
		file_paths = std::move(result.paths);
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			if (child.IsNull()) {
				throw InvalidInputException("read_sequences_sam: file path list cannot contain NULL");
			}
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

	// Validate all files exist (skip stdin and remote paths)
	for (const auto &path : file_paths) {
		if (!IsStdinPath(path) && !miint::RemoteFileHelper::IsRemotePath(path) && !fs.FileExists(path)) {
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
	auto &fs = FileSystem::GetFileSystem(context);

	return duckdb::make_uniq<GlobalState>(data.file_paths, fs, data.uses_stdin);
}

unique_ptr<LocalTableFunctionState> ReadSequencesSamTableFunction::InitLocal(ExecutionContext &context,
                                                                             TableFunctionInitInput &input,
                                                                             GlobalTableFunctionState *global_state) {
	return duckdb::make_uniq<LocalState>();
}

void ReadSequencesSamTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();
	auto &local_state = data_p.local_state->Cast<LocalState>();

	// Output schema (fixed): 0=sequence_index 1=read_id 2=comment 3=sequence1 4=sequence2
	// 5=qual1 6=qual2 [7=filepath]. read_id / sequence1 / qual1 are decoded straight from
	// the raw bam1_t into these vectors -- no intermediate std::string, no SAMRecordBatch
	// SOA, and none of the alignment metadata read_sequences_sam discards.
	auto read_id_data = FlatVector::GetData<string_t>(output.data[1]);
	auto sequence1_data = FlatVector::GetData<string_t>(output.data[3]);
	auto &qual1_vector = output.data[5];
	auto qual1_entries = FlatVector::GetData<list_entry_t>(qual1_vector);
	auto &qual_scratch = local_state.qual_scratch;

	std::string current_filepath;
	idx_t count = 0;

	// Loop until we decode at least one record or run out of files.
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

		auto &reader = *global_state.readers[local_state.current_file_idx];
		current_filepath = global_state.filepaths[local_state.current_file_idx];

		qual_scratch.clear();
		count = 0;
		while (count < STANDARD_VECTOR_SIZE) {
			const bam1_t *aln = reader.read_raw();
			if (!aln) {
				break; // end of this file
			}

			const int seq_len = aln->core.l_qseq;

			// Primary/unmapped reads must carry sequence (matches parse_record_to_batch).
			const bool is_unmapped = (aln->core.flag & 0x4) != 0;
			const bool is_secondary = (aln->core.flag & 0x100) != 0;
			const bool is_supplementary = (aln->core.flag & 0x800) != 0;
			if (seq_len == 0 && (is_unmapped || (!is_secondary && !is_supplementary))) {
				throw std::runtime_error("Primary/unmapped read missing sequence (SEQ='*'): " +
				                         std::string(reinterpret_cast<const char *>(aln->data)));
			}

			// read_id: copy QNAME once into the vector's string heap.
			read_id_data[count] = StringVector::AddString(output.data[1], reinterpret_cast<const char *>(aln->data));

			// sequence1: decode 4-bit nucleotides straight into the vector buffer (one pass).
			string_t seq_str = StringVector::EmptyString(output.data[3], seq_len);
			char *seq_dest = seq_str.GetDataWriteable();
			const uint8_t *seq = bam_get_seq(aln);
			for (int i = 0; i < seq_len; i++) {
				seq_dest[i] = seq_nt16_str[bam_seqi(seq, i)];
			}
			seq_str.Finalize();
			sequence1_data[count] = seq_str;

			// qual1: raw Phred scores copied straight across (no +33/-33 round trip).
			// Absent quality (first byte 0xFF) or empty sequence -> NULL list.
			const uint8_t *qual = bam_get_qual(aln);
			qual1_entries[count].offset = qual_scratch.size();
			if (seq_len > 0 && qual[0] != 0xFF) {
				qual_scratch.insert(qual_scratch.end(), qual, qual + seq_len);
				qual1_entries[count].length = static_cast<uint64_t>(seq_len);
			} else {
				qual1_entries[count].length = 0;
			}

			count++;
		}

		if (count == 0) {
			// File exhausted before producing any record; release and try the next.
			local_state.has_file = false;
			continue;
		}
		break;
	}

	// Finalize qual1 child with a single bulk copy of the accumulated raw Phred bytes.
	ListVector::Reserve(qual1_vector, qual_scratch.size());
	ListVector::SetListSize(qual1_vector, qual_scratch.size());
	if (!qual_scratch.empty()) {
		auto &qual_child = ListVector::GetEntry(qual1_vector);
		auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
		std::memcpy(qual_child_data, qual_scratch.data(), qual_scratch.size());
		FlatVector::Validity(qual_child).SetAllValid(qual_scratch.size());
	}
	// Per-row qual1 nullability: zero-length entries (absent quality) are NULL.
	auto &qual1_validity = FlatVector::Validity(qual1_vector);
	qual1_validity.SetAllValid(count);
	for (idx_t j = 0; j < count; j++) {
		if (qual1_entries[j].length == 0) {
			qual1_validity.SetInvalid(j);
		}
	}

	// sequence_index: per-file monotonic counter.
	uint64_t start_sequence_index = global_state.file_sequence_counters[local_state.current_file_idx];
	global_state.file_sequence_counters[local_state.current_file_idx] += count;
	auto sequence_index_data = FlatVector::GetData<int64_t>(output.data[0]);
	for (idx_t j = 0; j < count; j++) {
		sequence_index_data[j] = static_cast<int64_t>(start_sequence_index + j);
	}

	// Constant-NULL columns: comment, sequence2, qual2 (single-end BAM, no comment field).
	SetResultVectorNull(output.data[2]);
	SetResultVectorNull(output.data[4]);
	SetResultVectorNull(output.data[6]);

	if (bind_data.include_filepath) {
		SetResultVectorFilepath(output.data[7], current_filepath);
	}

	output.SetCardinality(count);
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
