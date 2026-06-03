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
		result->has_r2_columns = has_r2_columns;
		result->split_output = split_output;
		result->names = names;
		result->indices = indices;
		// No FASTA-specific fields to copy
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastaCopyBindData>();
		return interleave == other.interleave && id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment && compression == other.compression &&
		       file_path == other.file_path && has_r2_columns == other.has_r2_columns &&
		       split_output == other.split_output && flush_size == other.flush_size && names == other.names;
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

	// Column presence only; whether each record is actually paired is decided per-row at write
	// time from sequence2 NULL-ness (read_fastx always emits sequence2, NULL for single-end).
	result->has_r2_columns = has_sequence2;

	// Parse common parameters
	CommonCopyParameters common_params;

	for (auto &option : input.info.options) {
		if (!StringUtil::CIEquals(option.first, "interleave") &&
		    !StringUtil::CIEquals(option.first, "id_as_sequence_index") &&
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

	// {ORIENTATION} in the path => split R1/R2 files; otherwise a single (optionally interleaved) file.
	ResolveOutputMode(result->file_path, result->interleave, result->split_output);

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
// Write a FASTA record directly to the stream without building an intermediate string.
// Avoids O(seq.size()) extra copies per record -- the previous implementation copied the
// sequence bytes twice (into an intermediate `record` string, then into the stream). For
// a FASTA of bacterial genomes (~5 MB/record) that's ~10 MB of extra memory bandwidth per
// record; direct writes let the sequence flow kseq->MemoryStream with a single memcpy.
static void WriteFastaRecordToBuffer(MemoryStream &stream, const char *id, idx_t id_size, const char *seq,
                                     idx_t seq_size, const char *comment, idx_t comment_size) {
	stream.WriteData(const_data_ptr_cast(">"), 1);
	stream.WriteData(const_data_ptr_cast(id), id_size);
	if (comment_size > 0) {
		stream.WriteData(const_data_ptr_cast(" "), 1);
		stream.WriteData(const_data_ptr_cast(comment), comment_size);
	}
	stream.WriteData(const_data_ptr_cast("\n"), 1);
	stream.WriteData(const_data_ptr_cast(seq), seq_size);
	stream.WriteData(const_data_ptr_cast("\n"), 1);
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

	// The schema may carry a sequence2 column even for single-end data (read_fastx always emits it
	// as NULL). Whether a given row is paired is decided per-row below from sequence2 NULL-ness.
	bool has_r2 = fdata.has_r2_columns;
	if (has_r2) {
		input.data[indices.sequence2_idx].ToUnifiedFormat(input.size(), seq2_data);
	}

	// Get references to local buffers
	auto &stream_r1 = *lstate.writer_state_r1->stream;
	MemoryStream *stream_r2 = fdata.split_output ? lstate.writer_state_r2->stream.get() : nullptr;

	// Build all records into local buffer(s) - NO LOCK.
	// Hot-path rule: never call string_t::GetString() on the sequence -- it makes a heap
	// copy of the (potentially multi-MB) sequence per record. Use GetData()/GetSize() so
	// the sequence bytes flow directly from the DuckDB string heap into the MemoryStream.
	string id_buf; // only used when id_as_sequence_index is true
	for (idx_t row = 0; row < input.size(); row++) {
		auto row_idx = read_id_data.sel->get_index(row);

		// Validate required fields are not NULL
		if (!read_id_data.validity.RowIsValid(row_idx)) {
			throw InvalidInputException("NULL value in read_id column (row %llu)", row);
		}

		// Get ID
		const char *id_ptr;
		idx_t id_size;
		if (fdata.id_as_sequence_index) {
			auto seq_idx_data = UnifiedVectorFormat::GetData<int64_t>(sequence_index_data);
			auto seq_idx_row = sequence_index_data.sel->get_index(row);
			id_buf = to_string(seq_idx_data[seq_idx_row]);
			id_ptr = id_buf.data();
			id_size = id_buf.size();
		} else {
			id_ptr = read_ids[row_idx].GetData();
			id_size = read_ids[row_idx].GetSize();
		}

		// Get comment (may be absent or NULL)
		const char *comment_ptr = nullptr;
		idx_t comment_size = 0;
		if (fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX) {
			auto comment_strings = UnifiedVectorFormat::GetData<string_t>(comment_data);
			auto comment_row = comment_data.sel->get_index(row);
			if (comment_data.validity.RowIsValid(comment_row)) {
				comment_ptr = comment_strings[comment_row].GetData();
				comment_size = comment_strings[comment_row].GetSize();
			}
		}

		// Get R1 data
		auto seq1_row = seq1_data.sel->get_index(row);
		if (!seq1_data.validity.RowIsValid(seq1_row)) {
			throw InvalidInputException("NULL value in sequence1 column (row %llu)", row);
		}
		const char *seq1_ptr = seq1_strings[seq1_row].GetData();
		idx_t seq1_size = seq1_strings[seq1_row].GetSize();

		// Write R1 record to local buffer
		WriteFastaRecordToBuffer(stream_r1, id_ptr, id_size, seq1_ptr, seq1_size, comment_ptr, comment_size);
		lstate.writer_state_r1->written_anything = true;

		// Handle R2. A record is paired iff sequence2 is non-NULL; a NULL sequence2 means a
		// single-end record.
		if (has_r2) {
			auto seq2_row = seq2_data.sel->get_index(row);
			if (seq2_data.validity.RowIsValid(seq2_row)) {
				lstate.saw_paired = true;

				auto seq2_strings = UnifiedVectorFormat::GetData<string_t>(seq2_data);
				const char *seq2_ptr = seq2_strings[seq2_row].GetData();
				idx_t seq2_size = seq2_strings[seq2_row].GetSize();

				if (fdata.split_output) {
					// Write R2 to the separate R2 buffer (file created lazily on first flush).
					WriteFastaRecordToBuffer(*stream_r2, id_ptr, id_size, seq2_ptr, seq2_size, comment_ptr,
					                         comment_size);
					lstate.writer_state_r2->written_anything = true;
				} else if (fdata.interleave) {
					// Write R2 interleaved into the same buffer as R1.
					WriteFastaRecordToBuffer(stream_r1, id_ptr, id_size, seq2_ptr, seq2_size, comment_ptr,
					                         comment_size);
				} else {
					// Single-file output but this record is paired and there is nowhere to put R2.
					throw InvalidInputException(
					    "Row %llu has paired-end data (sequence2 is set) but the output is a single file. Use the "
					    "{ORIENTATION} placeholder in the output path to write split R1/R2 files, or set "
					    "INTERLEAVE=true to interleave R1/R2 into one file.",
					    row);
				}
			} else {
				lstate.saw_single = true;
			}
		}

		// Flush per-record once the threshold is crossed, not per-chunk. Large records
		// (bacterial genomes, plant chromosomes) can each be multiple MB; batching them
		// until the end of a DataChunk would let the buffer grow to many GB before the
		// first flush, which blows past gzip's 4 GiB write limit and chews memory.
		if (lstate.writer_state_r1->stream->GetPosition() >= lstate.writer_state_r1->flush_size) {
			FlushFormatBuffer(*lstate.writer_state_r1, *gstate.file_r1, gstate.lock);
		}
		if (stream_r2 && lstate.writer_state_r2->stream->GetPosition() >= lstate.writer_state_r2->flush_size) {
			FlushR2Buffer(*lstate.writer_state_r2, gstate);
		}
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
