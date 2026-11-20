#include "copy_fasta.hpp"
#include "copy_format_common.hpp"
#include "id_column_utils.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/planner/binder.hpp"
#include <sstream>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct FastaCopyBindData : public SequenceCopyBindData {
	// LINE_WIDTH: wrap each sequence at this many characters; 0 (no LINE_WIDTH) writes it on one line.
	idx_t line_width = 0;

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<FastaCopyBindData>();
		// Copy all base class fields
		result->interleave = interleave;
		result->id_as_sequence_index = id_as_sequence_index;
		result->include_comment = include_comment;
		result->compression = compression;
		result->compression_level = compression_level;
		result->file_path = file_path;
		result->has_r2_columns = has_r2_columns;
		result->split_output = split_output;
		result->names = names;
		result->indices = indices;
		result->read_id_type = read_id_type;
		result->line_width = line_width;
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastaCopyBindData>();
		return interleave == other.interleave && id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment && compression == other.compression &&
		       compression_level == other.compression_level && file_path == other.file_path &&
		       has_r2_columns == other.has_r2_columns && split_output == other.split_output && names == other.names &&
		       read_id_type == other.read_id_type && line_width == other.line_width;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> FastaCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                              const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<FastaCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = IdentifiersToStrings(names);

	// Detect and store column indices (computed once at bind time)
	result->indices.FindIndices(IdentifiersToStrings(names));

	bool has_sequence1 = result->indices.sequence1_idx != DConstants::INVALID_INDEX;
	bool has_sequence2 = result->indices.sequence2_idx != DConstants::INVALID_INDEX;
	bool has_read_id = result->indices.read_id_idx != DConstants::INVALID_INDEX;
	bool has_sequence_index = result->indices.sequence_index_idx != DConstants::INVALID_INDEX;

	// Validate required columns
	ValidateRequiredColumns(has_read_id, has_sequence1, "FASTA");

	// read_id may be VARCHAR, BIGINT, or UUID (mirrors COPY FORMAT SAM). Validate at
	// bind so a bad type is a clean BinderError, not a fatal internal error when writing.
	if (!IsAllowedIdType(sql_types[result->indices.read_id_idx])) {
		throw BinderException("Column 'read_id' must be %s", AllowedIdTypeList());
	}
	result->read_id_type = sql_types[result->indices.read_id_idx];

	// Column presence only; whether each record is actually paired is decided per-row at write
	// time from sequence2 NULL-ness (read_fastx always emits sequence2, NULL for single-end).
	result->has_r2_columns = has_sequence2;

	// Parse common parameters
	CommonCopyParameters common_params;
	common_params.ParseFromOptions(input.info.options, result->file_path);

	result->interleave = common_params.interleave;
	result->id_as_sequence_index = common_params.id_as_sequence_index;
	result->include_comment = common_params.include_comment;
	result->compression = common_params.compression;
	AddGzipFileExtension(input, result->compression);
	result->compression_level = common_params.compression_level;

	for (auto &option : input.info.options) {
		if (option.first == "line_width") {
			// Declared BIGINT, so DuckDB hands over exactly one non-NULL BIGINT.
			const auto width = option.second[0].GetValue<int64_t>();
			if (width < 1) {
				throw BinderException("LINE_WIDTH must be at least 1, got %lld", static_cast<long long>(width));
			}
			result->line_width = static_cast<idx_t>(width);
		}
	}

	// {ORIENTATION} in the path => split R1/R2 files; otherwise a single (optionally interleaved) file.
	ResolveOutputMode(result->file_path, result->interleave, result->split_output);

	// Validate sequence_index parameter
	ValidateSequenceIndexParameter(result->id_as_sequence_index, has_sequence_index);

	// Same rationale as the read_id check above, for every other column the writer reads with
	// a typed accessor: a non-VARCHAR sequence1 otherwise aborts with a fatal internal error
	// that invalidates the database. Quality columns are not validated: the FASTA writer
	// ignores any qual1/qual2 present in the relation rather than reading them, so rejecting
	// them would fail queries that are actually fine.
	ValidateSequenceColumnTypes(result->indices, sql_types, /*validate_quals=*/false,
	                            /*validate_comment=*/result->include_comment,
	                            /*validate_sequence_index=*/result->id_as_sequence_index);

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
// Record text
//===--------------------------------------------------------------------===//
// Writes a FASTA record straight to the stream, with no intermediate string: the sequence is copied once, which for
// genome-length records (~5 MB each) saves one copy of every record.
// With line_width > 0 the sequence is written in lines of line_width, the last one shorter, and no empty line
// after a sequence whose length is a multiple of it (as Prodigal writes -a and -d); an empty sequence is an empty
// line either way.
static void WriteFastaRecordToBuffer(MemoryStream &stream, const char *id, idx_t id_size, const char *seq,
                                     idx_t seq_size, const char *comment, idx_t comment_size, idx_t line_width) {
	stream.WriteData(const_data_ptr_cast(">"), 1);
	stream.WriteData(const_data_ptr_cast(id), id_size);
	if (comment_size > 0) {
		stream.WriteData(const_data_ptr_cast(" "), 1);
		stream.WriteData(const_data_ptr_cast(comment), comment_size);
	}
	stream.WriteData(const_data_ptr_cast("\n"), 1);
	if (line_width == 0 || seq_size <= line_width) {
		stream.WriteData(const_data_ptr_cast(seq), seq_size);
		stream.WriteData(const_data_ptr_cast("\n"), 1);
		return;
	}
	for (idx_t offset = 0; offset < seq_size; offset += line_width) {
		stream.WriteData(const_data_ptr_cast(seq + offset), MinValue(line_width, seq_size - offset));
		stream.WriteData(const_data_ptr_cast("\n"), 1);
	}
}

//===--------------------------------------------------------------------===//
// Finalize (use shared implementation)
//===--------------------------------------------------------------------===//
static void FastaCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<FastaCopyGlobalState>();
	SequenceCopyFinalize(gstate);
}

//===--------------------------------------------------------------------===//
// Batch API (SequencePrepareBatch / SequenceFlushBatch)
//===--------------------------------------------------------------------===//
struct FastaRecord {
	const char *id;
	idx_t id_size;
	const char *comment;
	idx_t comment_size;
	const char *seq;
	idx_t seq_size;
};

// The FORMAT of SequencePrepareBatch for FASTA: a row's R1 record, then its R2 record when sequence2 is set.
class FastaBatchFormat {
public:
	using Record = FastaRecord;

	explicit FastaBatchFormat(const FastaCopyBindData &fdata_p) : fdata(fdata_p) {
	}

	template <class EMIT>
	void Visit(DataChunk &input, SequenceBatchData &batch, EMIT &&emit) {
		SequenceChunkReader reader(fdata, input);
		for (idx_t row = 0; row < input.size(); row++) {
			FastaRecord record;
			reader.Header(row, record.id, record.id_size, record.comment, record.comment_size);
			const auto &seq1 = reader.Sequence1(row);
			record.seq = seq1.GetData();
			record.seq_size = seq1.GetSize();
			emit(SequenceStream::R1, record);

			// A record is paired iff sequence2 is non-NULL.
			SequenceStream r2_stream;
			if (!reader.RouteR2(row, reader.HasSequence2(row), batch, r2_stream)) {
				continue;
			}
			const auto &seq2 = reader.Sequence2(row);
			record.seq = seq2.GetData();
			record.seq_size = seq2.GetSize();
			emit(r2_stream, record);
		}
	}

	// What WriteFastaRecordToBuffer writes: ">" id [" " comment] "\n", then the sequence in lines, each ending "\n".
	idx_t Size(const FastaRecord &record) const {
		const auto width = fdata.line_width;
		const idx_t lines = (width == 0 || record.seq_size <= width) ? 1 : (record.seq_size + width - 1) / width;
		return 1 + record.id_size + (record.comment_size > 0 ? 1 + record.comment_size : 0) + 1 + record.seq_size +
		       lines;
	}

	void Write(MemoryStream &stream, const FastaRecord &record) const {
		WriteFastaRecordToBuffer(stream, record.id, record.id_size, record.seq, record.seq_size, record.comment,
		                         record.comment_size, fdata.line_width);
	}

	// Write raises no data error, so there is nothing to check ahead of it.
	void Check(const FastaRecord &) const {
	}

private:
	const FastaCopyBindData &fdata;
};

// Reads only the bind data: see SequencePrepareBatch.
static unique_ptr<PreparedBatchData> FastaPrepareBatch(ClientContext &context, FunctionData &bind_data,
                                                       GlobalFunctionData &, unique_ptr<ColumnDataCollection> rows) {
	FastaBatchFormat format(bind_data.Cast<FastaCopyBindData>());
	return SequencePrepareBatch(context, format, std::move(rows));
}

static void FastaFlushBatch(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate,
                            PreparedBatchData &batch) {
	auto &fdata = bind_data.Cast<FastaCopyBindData>();
	FastaBatchFormat format(fdata);
	SequenceFlushBatch(context, gstate.Cast<FastaCopyGlobalState>(), batch.Cast<SequenceBatchData>(), format);
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
static constexpr const char *FASTA_WRITER_FORMAT = "_miint_fasta";

static void FastaCopyOptions(ClientContext &context, CopyOptionsInput &input) {
	DeclareSequenceCopyOptions(input);
	DeclareCopyOption(input, "line_width", LogicalType::BIGINT);
}

static BoundStatement FastaCopyPlan(Binder &binder, CopyStatement &stmt) {
	return SequenceCopyPlan(binder, stmt, FASTA_WRITER_FORMAT);
}

// FORMAT FASTA: the options and the plan callback, which rebinds as the writer below (see SequenceCopyPlan).
CopyFunction CopyFastaFunction::GetFunction() {
	CopyFunction func("fasta");
	func.copy_options = FastaCopyOptions;
	func.plan = FastaCopyPlan;
	func.extension = "fasta";
	return func;
}

// FORMAT _miint_fasta: the writer. Named directly, it writes without the default batch cap.
static CopyFunction GetWriterFunction() {
	CopyFunction func(FASTA_WRITER_FORMAT);
	func.copy_options = FastaCopyOptions;
	func.copy_to_bind = FastaCopyBind;
	func.copy_to_initialize_global = FastaCopyInitializeGlobal;
	func.copy_to_finalize = FastaCopyFinalize;
	func.execution_mode = SequenceCopyExecutionMode;
	func.prepare_batch = FastaPrepareBatch;
	func.flush_batch = FastaFlushBatch;
	func.extension = "fasta";
	return func;
}

void CopyFastaFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
	loader.RegisterFunction(GetWriterFunction());
}
}; // namespace duckdb
