#include "copy_fasta.hpp"
#include "copy_format_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include <sstream>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data (inherits from shared base - no FASTA-specific fields needed)
//===--------------------------------------------------------------------===//
struct FastaCopyBindData : public SequenceCopyBindData {
	// No FASTA-specific fields (unlike FASTQ which has qual_offset)

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<FastaCopyBindData>();
		// Copy all base class fields
		result->interleave = interleave;
		result->id_as_sequence_index = id_as_sequence_index;
		result->include_comment = include_comment;
		result->compression = compression;
		result->flush_size = flush_size;
		result->file_path = file_path;
		result->is_paired = is_paired;
		result->names = names;
		result->indices = indices;
		// No FASTA-specific fields to copy
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastaCopyBindData>();
		return interleave == other.interleave && id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment && compression == other.compression &&
		       file_path == other.file_path && is_paired == other.is_paired && flush_size == other.flush_size &&
		       names == other.names;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> FastaCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                              const vector<string> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<FastaCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = names;

	// Detect and store column indices (computed once at bind time)
	result->indices.FindIndices(names);

	bool has_sequence1 = result->indices.sequence1_idx != DConstants::INVALID_INDEX;
	bool has_sequence2 = result->indices.sequence2_idx != DConstants::INVALID_INDEX;
	bool has_read_id = result->indices.read_id_idx != DConstants::INVALID_INDEX;
	bool has_sequence_index = result->indices.sequence_index_idx != DConstants::INVALID_INDEX;

	// Validate required columns
	ValidateRequiredColumns(has_read_id, has_sequence1, "FASTA");

	result->is_paired = has_sequence2;

	// Parse common parameters
	CommonCopyParameters common_params;
	bool has_interleave_param = false;

	for (auto &option : input.info.options) {
		if (StringUtil::CIEquals(option.first, "interleave")) {
			has_interleave_param = true;
		} else if (!StringUtil::CIEquals(option.first, "id_as_sequence_index") &&
		           !StringUtil::CIEquals(option.first, "include_comment") &&
		           !StringUtil::CIEquals(option.first, "compression")) {
			throw BinderException("Unknown option for COPY FORMAT FASTA: %s", option.first);
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

	return result;
}

//===--------------------------------------------------------------------===//
// Global State (use shared implementation)
//===--------------------------------------------------------------------===//
using FastaCopyGlobalState = SequenceCopyGlobalState;

static unique_ptr<GlobalFunctionData> FastaCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                const string &file_path) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	return SequenceCopyInitializeGlobal(context, fdata, file_path);
}

//===--------------------------------------------------------------------===//
// Local State (use shared implementation)
//===--------------------------------------------------------------------===//
using FastaCopyLocalState = SequenceCopyLocalState;

static unique_ptr<LocalFunctionData> FastaCopyInitializeLocal(ExecutionContext &context, FunctionData &bind_data) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	return SequenceCopyInitializeLocal(context, fdata);
}

//===--------------------------------------------------------------------===//
// Sink
//===--------------------------------------------------------------------===//
static void WriteFastaRecordToBuffer(MemoryStream &stream, const string &id, const string &seq, const string &comment) {
	// Pre-calculate total size to avoid reallocations
	idx_t size = 1 + id.size() + 1 + seq.size() + 1; // > + id + \n + seq + \n
	if (!comment.empty()) {
		size += 1 + comment.size(); // space + comment
	}

	// Build record string with pre-reserved capacity
	string record;
	record.reserve(size);
	record += '>';
	record += id;
	if (!comment.empty()) {
		record += ' ';
		record += comment;
	}
	record += '\n';
	record += seq;
	record += '\n';

	// Write to stream
	stream.WriteData(const_data_ptr_cast(record.c_str()), record.size());
}

static void FastaCopySink(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                          LocalFunctionData &lstate_p, DataChunk &input) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	auto &lstate = lstate_p.Cast<FastaCopyLocalState>();

	// Use pre-computed column indices from bind data
	auto &indices = fdata.indices;

	// Process each row
	UnifiedVectorFormat read_id_data, sequence_index_data, comment_data;
	UnifiedVectorFormat seq1_data, seq2_data;

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

	bool is_paired = fdata.is_paired;
	if (is_paired) {
		input.data[indices.sequence2_idx].ToUnifiedFormat(input.size(), seq2_data);
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

		// Validate required fields are not NULL
		if (!read_id_data.validity.RowIsValid(row_idx)) {
			throw InvalidInputException("NULL value in read_id column (row %llu)", row);
		}

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
		if (!seq1_data.validity.RowIsValid(seq1_row)) {
			throw InvalidInputException("NULL value in sequence1 column (row %llu)", row);
		}
		string seq1 = seq1_strings[seq1_row].GetString();

		// Write R1 record to local buffer
		WriteFastaRecordToBuffer(stream_r1, id, seq1, comment);
		lstate.writer_state_r1->written_anything = true;

		// Handle R2 for paired-end
		if (is_paired) {
			auto seq2_strings = UnifiedVectorFormat::GetData<string_t>(seq2_data);
			auto seq2_row = seq2_data.sel->get_index(row);
			if (!seq2_data.validity.RowIsValid(seq2_row)) {
				throw InvalidInputException("NULL value in sequence2 column (row %llu)", row);
			}
			string seq2 = seq2_strings[seq2_row].GetString();

			if (fdata.interleave) {
				// Write R2 to same buffer
				WriteFastaRecordToBuffer(stream_r1, id, seq2, comment);
			} else {
				// Write R2 to separate buffer
				WriteFastaRecordToBuffer(*stream_r2, id, seq2, comment);
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
// Combine (use shared implementation)
//===--------------------------------------------------------------------===//
static void FastaCopyCombine(ExecutionContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p,
                             LocalFunctionData &lstate_p) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	auto &lstate = lstate_p.Cast<FastaCopyLocalState>();
	SequenceCopyCombine(fdata, gstate, lstate);
}

//===--------------------------------------------------------------------===//
// Finalize (use shared implementation)
//===--------------------------------------------------------------------===//
static void FastaCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	SequenceCopyFinalize(gstate);
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
CopyFunction CopyFastaFunction::GetFunction() {
	CopyFunction func("fasta");
	func.copy_to_bind = FastaCopyBind;
	func.copy_to_initialize_local = FastaCopyInitializeLocal;
	func.copy_to_initialize_global = FastaCopyInitializeGlobal;
	func.copy_to_sink = FastaCopySink;
	func.copy_to_combine = FastaCopyCombine;
	func.copy_to_finalize = FastaCopyFinalize;
	return func;
}

void CopyFastaFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}
}; // namespace duckdb
