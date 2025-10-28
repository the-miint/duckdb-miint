#include "copy_fastq.hpp"
#include "copy_format_common.hpp"
#include "QualScore.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include <sstream>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct FastqCopyBindData : public FunctionData {
	bool interleave = false;
	bool id_as_sequence_index = false;
	bool include_comment = false;
	uint8_t qual_offset = 33;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	idx_t flush_size = 1024 * 1024; // NOLINT
	string file_path;
	bool is_paired = false;
	vector<string> names;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<FastqCopyBindData>();
		result->interleave = interleave;
		result->id_as_sequence_index = id_as_sequence_index;
		result->include_comment = include_comment;
		result->qual_offset = qual_offset;
		result->compression = compression;
		result->flush_size = flush_size;
		result->file_path = file_path;
		result->is_paired = is_paired;
		result->names = names;
		return std::move(result);
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastqCopyBindData>();
		return interleave == other.interleave && id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment && qual_offset == other.qual_offset &&
		       compression == other.compression && file_path == other.file_path && is_paired == other.is_paired;
	}
};

//===--------------------------------------------------------------------===//
// Helper Functions
//===--------------------------------------------------------------------===//
static string EncodeQuality(const vector<uint8_t> &qual, uint8_t offset) {
	string result;
	result.reserve(qual.size());
	for (auto q : qual) {
		result.push_back(static_cast<char>(q + offset));
	}
	return result;
}

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> FastqCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                              const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<FastqCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;

	// Detect columns
	ColumnIndices indices;
	indices.FindIndices(names);

	bool has_sequence1 = indices.sequence1_idx != DConstants::INVALID_INDEX;
	bool has_sequence2 = indices.sequence2_idx != DConstants::INVALID_INDEX;
	bool has_read_id = indices.read_id_idx != DConstants::INVALID_INDEX;
	bool has_qual1 = indices.qual1_idx != DConstants::INVALID_INDEX;
	bool has_qual2 = indices.qual2_idx != DConstants::INVALID_INDEX;
	bool has_sequence_index = indices.sequence_index_idx != DConstants::INVALID_INDEX;

	// Validate required columns
	ValidateRequiredColumns(has_read_id, has_sequence1, "FASTQ");
	if (!has_qual1) {
		throw BinderException("COPY FORMAT FASTQ requires 'qual1' column");
	}

	result->is_paired = has_sequence2 && has_qual2;

	// Parse common parameters
	CommonCopyParameters common_params;

	Value qual_offset_param;
	bool has_interleave_param = false;

	for (auto &option : input.info.options) {
		if (StringUtil::CIEquals(option.first, "interleave")) {
			has_interleave_param = true;
		} else if (StringUtil::CIEquals(option.first, "qual_offset")) {
			qual_offset_param = option.second[0];
		} else if (!StringUtil::CIEquals(option.first, "id_as_sequence_index") &&
		           !StringUtil::CIEquals(option.first, "include_comment") &&
		           !StringUtil::CIEquals(option.first, "compression")) {
			throw BinderException("Unknown option for COPY FORMAT FASTQ: %s", option.first);
		}
	}

	common_params.ParseFromOptions(input.info.options, result->file_path);

	result->interleave = common_params.interleave;
	result->id_as_sequence_index = common_params.id_as_sequence_index;
	result->include_comment = common_params.include_comment;
	result->compression = common_params.compression;
	result->flush_size = common_params.flush_size;

	// Validate paired-end parameters
	ValidatePairedEndParameters(result->is_paired, has_interleave_param, result->interleave, result->file_path);

	// Validate sequence_index parameter
	ValidateSequenceIndexParameter(result->id_as_sequence_index, has_sequence_index);

	// Parse FASTQ-specific qual_offset parameter
	if (!qual_offset_param.IsNull()) {
		int64_t offset = qual_offset_param.GetValue<int64_t>();
		if (offset != 33 && offset != 64) {
			throw BinderException("QUAL_OFFSET must be 33 or 64");
		}
		result->qual_offset = static_cast<uint8_t>(offset);
	}

	return std::move(result);
}

//===--------------------------------------------------------------------===//
// Global State
//===--------------------------------------------------------------------===//
struct FastqCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	unique_ptr<CopyFileHandle> file_r1;
	unique_ptr<CopyFileHandle> file_r2;
	bool is_paired = false;
	bool interleave = false;
};

static unique_ptr<GlobalFunctionData> FastqCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                const string &file_path) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	auto &fs = FileSystem::GetFileSystem(context);

	auto gstate = make_uniq<FastqCopyGlobalState>();
	gstate->is_paired = fdata.is_paired;
	gstate->interleave = fdata.interleave;

	if (fdata.is_paired && !fdata.interleave) {
		// Split mode: open two files
		string path_r1 = SubstituteOrientation(file_path, "R1");
		string path_r2 = SubstituteOrientation(file_path, "R2");

		gstate->file_r1 = make_uniq<CopyFileHandle>(fs, path_r1, fdata.compression);
		gstate->file_r2 = make_uniq<CopyFileHandle>(fs, path_r2, fdata.compression);
	} else {
		// Interleaved or single-end: open one file
		gstate->file_r1 = make_uniq<CopyFileHandle>(fs, file_path, fdata.compression);
	}

	return std::move(gstate);
}

//===--------------------------------------------------------------------===//
// Local State
//===--------------------------------------------------------------------===//
struct FastqCopyLocalState : public LocalFunctionData {
	unique_ptr<FormatWriterState> writer_state_r1;
	unique_ptr<FormatWriterState> writer_state_r2; // For split paired-end
};

static unique_ptr<LocalFunctionData> FastqCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	auto lstate = make_uniq<FastqCopyLocalState>();

	lstate->writer_state_r1 = make_uniq<FormatWriterState>(context.client, fdata.flush_size);

	if (fdata.is_paired && !fdata.interleave) {
		lstate->writer_state_r2 = make_uniq<FormatWriterState>(context.client, fdata.flush_size);
	}

	return std::move(lstate);
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void WriteFastqRecordToBuffer(MemoryStream &stream, const string &id, const string &seq,
                                     const vector<uint8_t> &qual, uint8_t qual_offset, const string &comment) {
	// Build record string
	string record = "@" + id;
	if (!comment.empty()) {
		record += " " + comment;
	}
	record += "\n" + seq + "\n+\n" + EncodeQuality(qual, qual_offset) + "\n";

	// Write to stream
	stream.WriteData(const_data_ptr_cast(record.c_str()), record.size());
}

static void FastqCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                          LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	auto &gstate = gstate_p.Cast<FastqCopyGlobalState>();
	auto &lstate = lstate_p.Cast<FastqCopyLocalState>();

	// Get column indices (no lock needed yet!)
	ColumnIndices indices;
	indices.FindIndices(fdata.names);

	// Process each row
	UnifiedVectorFormat read_id_data, sequence_index_data, comment_data;
	UnifiedVectorFormat seq1_data, seq2_data, qual1_data, qual2_data;

	input.data[indices.read_id_idx].ToUnifiedFormat(input.size(), read_id_data);
	auto read_ids = UnifiedVectorFormat::GetData<string_t>(read_id_data);

	if (fdata.id_as_sequence_index) {
		input.data[indices.sequence_index_idx].ToUnifiedFormat(input.size(), sequence_index_data);
	}

	if (fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX) {
		input.data[indices.comment_idx].ToUnifiedFormat(input.size(), comment_data);
	}

	input.data[indices.sequence1_idx].ToUnifiedFormat(input.size(), seq1_data);
	auto seq1_strings = UnifiedVectorFormat::GetData<string_t>(seq1_data);

	input.data[indices.qual1_idx].ToUnifiedFormat(input.size(), qual1_data);

	bool is_paired = fdata.is_paired;
	if (is_paired) {
		input.data[indices.sequence2_idx].ToUnifiedFormat(input.size(), seq2_data);
		input.data[indices.qual2_idx].ToUnifiedFormat(input.size(), qual2_data);
	}

	// Get references to local buffers
	auto &stream_r1 = *lstate.writer_state_r1->stream;
	MemoryStream *stream_r2 = nullptr;
	if (is_paired && !fdata.interleave) {
		stream_r2 = lstate.writer_state_r2->stream.get();
	}

	// Build all records into local buffer(s) - NO LOCK
	for (idx_t row = 0; row < input.size(); row++) {
		auto row_idx = read_id_data.sel->get_index(row);

		// Get ID
		string id;
		if (fdata.id_as_sequence_index) {
			auto seq_idx_data = UnifiedVectorFormat::GetData<int64_t>(sequence_index_data);
			auto seq_idx_row = sequence_index_data.sel->get_index(row);
			id = to_string(seq_idx_data[seq_idx_row]);
		} else {
			id = read_ids[row_idx].GetString();
		}

		// Get comment
		string comment;
		if (fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX) {
			auto comment_strings = UnifiedVectorFormat::GetData<string_t>(comment_data);
			auto comment_row = comment_data.sel->get_index(row);
			if (comment_data.validity.RowIsValid(comment_row)) {
				comment = comment_strings[comment_row].GetString();
			}
		}

		// Get R1 data
		auto seq1_row = seq1_data.sel->get_index(row);
		string seq1 = seq1_strings[seq1_row].GetString();

		auto qual1_row = qual1_data.sel->get_index(row);
		auto qual1_list = ListVector::GetEntry(input.data[indices.qual1_idx]);
		auto qual1_list_data = FlatVector::GetData<uint8_t>(qual1_list);
		auto qual1_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual1_data);
		vector<uint8_t> qual1_vec;
		for (idx_t i = 0; i < qual1_entries[qual1_row].length; i++) {
			qual1_vec.push_back(qual1_list_data[qual1_entries[qual1_row].offset + i]);
		}

		// Write R1 record to local buffer
		WriteFastqRecordToBuffer(stream_r1, id, seq1, qual1_vec, fdata.qual_offset, comment);
		lstate.writer_state_r1->written_anything = true;

		// Handle R2 for paired-end
		if (is_paired) {
			auto seq2_strings = UnifiedVectorFormat::GetData<string_t>(seq2_data);
			auto seq2_row = seq2_data.sel->get_index(row);
			string seq2 = seq2_strings[seq2_row].GetString();

			auto qual2_row = qual2_data.sel->get_index(row);
			auto qual2_list = ListVector::GetEntry(input.data[indices.qual2_idx]);
			auto qual2_list_data = FlatVector::GetData<uint8_t>(qual2_list);
			auto qual2_entries = UnifiedVectorFormat::GetData<list_entry_t>(qual2_data);
			vector<uint8_t> qual2_vec;
			for (idx_t i = 0; i < qual2_entries[qual2_row].length; i++) {
				qual2_vec.push_back(qual2_list_data[qual2_entries[qual2_row].offset + i]);
			}

			if (fdata.interleave) {
				// Write R2 to same buffer
				WriteFastqRecordToBuffer(stream_r1, id, seq2, qual2_vec, fdata.qual_offset, comment);
			} else {
				// Write R2 to separate buffer
				WriteFastqRecordToBuffer(*stream_r2, id, seq2, qual2_vec, fdata.qual_offset, comment);
				lstate.writer_state_r2->written_anything = true;
			}
		}
	}

	// Check if we need to flush (buffer exceeded threshold)
	if (lstate.writer_state_r1->stream->GetPosition() >= lstate.writer_state_r1->flush_size) {
		FlushFormatBuffer(*lstate.writer_state_r1, *gstate.file_r1, gstate.lock);
	}

	if (stream_r2 && lstate.writer_state_r2->stream->GetPosition() >= lstate.writer_state_r2->flush_size) {
		FlushFormatBuffer(*lstate.writer_state_r2, *gstate.file_r2, gstate.lock);
	}
}

//===--------------------------------------------------------------------===//
// Combine
//===--------------------------------------------------------------------===//
static void FastqCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                             LocalFunctionData &lstate_p) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	auto &gstate = gstate_p.Cast<FastqCopyGlobalState>();
	auto &lstate = lstate_p.Cast<FastqCopyLocalState>();

	// Flush any remaining data in local buffers
	FlushFormatBuffer(*lstate.writer_state_r1, *gstate.file_r1, gstate.lock);

	if (fdata.is_paired && !fdata.interleave) {
		FlushFormatBuffer(*lstate.writer_state_r2, *gstate.file_r2, gstate.lock);
	}
}

//===--------------------------------------------------------------------===//
// Finalize
//===--------------------------------------------------------------------===//
static void FastqCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<FastqCopyGlobalState>();
	lock_guard<mutex> glock(gstate.lock);

	if (gstate.file_r1) {
		gstate.file_r1->Close();
	}
	if (gstate.file_r2) {
		gstate.file_r2->Close();
	}
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
CopyFunction CopyFastqFunction::GetFunction() {
	CopyFunction func("fastq");
	func.copy_to_bind = FastqCopyBind;
	func.copy_to_initialize_local = FastqCopyInitializeLocal;
	func.copy_to_initialize_global = FastqCopyInitializeGlobal;
	func.copy_to_sink = FastqCopySink;
	func.copy_to_combine = FastqCopyCombine;
	func.copy_to_finalize = FastqCopyFinalize;
	return func;
}

void CopyFastqFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
