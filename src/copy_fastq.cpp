#include "copy_fastq.hpp"
#include "copy_format_common.hpp"
#include "id_column_utils.hpp"
#include "fastq_encoder.hpp"
#include "QualScore.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/vector/list_vector.hpp"
#include "duckdb/planner/binder.hpp"
#include <sstream>

namespace duckdb {

//===--------------------------------------------------------------------===//
// Bind Data
//===--------------------------------------------------------------------===//
struct FastqCopyBindData : public SequenceCopyBindData {
	uint8_t qual_offset = 33; // FASTQ-specific: quality score offset

	unique_ptr<FunctionData> Copy() const override {
		auto result = make_uniq<FastqCopyBindData>();
		// Copy base class fields
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
		// Copy FASTQ-specific field
		result->qual_offset = qual_offset;
		return result;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<FastqCopyBindData>();
		return interleave == other.interleave && id_as_sequence_index == other.id_as_sequence_index &&
		       include_comment == other.include_comment && qual_offset == other.qual_offset &&
		       compression == other.compression && compression_level == other.compression_level &&
		       file_path == other.file_path && has_r2_columns == other.has_r2_columns &&
		       split_output == other.split_output && names == other.names && read_id_type == other.read_id_type;
	}
};

//===--------------------------------------------------------------------===//
// Bind
//===--------------------------------------------------------------------===//
static unique_ptr<FunctionData> FastqCopyBind(ClientContext &context, CopyFunctionBindInput &input,
                                              const vector<Identifier> &names, const vector<LogicalType> &sql_types) {
	auto result = make_uniq<FastqCopyBindData>();
	result->file_path = input.info.file_path;
	result->names = IdentifiersToStrings(names);

	// Detect and store column indices (computed once at bind time)
	result->indices.FindIndices(IdentifiersToStrings(names));

	bool has_sequence1 = result->indices.sequence1_idx != DConstants::INVALID_INDEX;
	bool has_sequence2 = result->indices.sequence2_idx != DConstants::INVALID_INDEX;
	bool has_read_id = result->indices.read_id_idx != DConstants::INVALID_INDEX;
	bool has_qual1 = result->indices.qual1_idx != DConstants::INVALID_INDEX;
	bool has_qual2 = result->indices.qual2_idx != DConstants::INVALID_INDEX;
	bool has_sequence_index = result->indices.sequence_index_idx != DConstants::INVALID_INDEX;

	// Validate required columns
	ValidateRequiredColumns(has_read_id, has_sequence1, "FASTQ");
	if (!has_qual1) {
		throw BinderException("COPY FORMAT FASTQ requires 'qual1' column");
	}

	// read_id may be VARCHAR, BIGINT, or UUID (mirrors COPY FORMAT SAM). Validate at
	// bind so a bad type is a clean BinderError, not a fatal internal error when writing.
	if (!IsAllowedIdType(sql_types[result->indices.read_id_idx])) {
		throw BinderException("Column 'read_id' must be %s", AllowedIdTypeList());
	}
	result->read_id_type = sql_types[result->indices.read_id_idx];

	// Column presence only; whether each record is actually paired is decided per-row at write
	// time from R2 NULL-ness (read_fastx always emits sequence2/qual2, NULL for single-end).
	result->has_r2_columns = has_sequence2 && has_qual2;

	// Parse common parameters
	CommonCopyParameters common_params;

	Value qual_offset_param;

	for (auto &option : input.info.options) {
		if (option.first == "qual_offset") {
			qual_offset_param = option.second[0];
		}
	}

	common_params.ParseFromOptions(input.info.options, result->file_path);

	result->interleave = common_params.interleave;
	result->id_as_sequence_index = common_params.id_as_sequence_index;
	result->include_comment = common_params.include_comment;
	result->compression = common_params.compression;
	AddGzipFileExtension(input, result->compression);
	result->compression_level = common_params.compression_level;

	// {ORIENTATION} in the path => split R1/R2 files; otherwise a single (optionally interleaved) file.
	ResolveOutputMode(result->file_path, result->interleave, result->split_output);

	// Validate sequence_index parameter
	ValidateSequenceIndexParameter(result->id_as_sequence_index, has_sequence_index);

	// Same rationale as the read_id check above, for every other column the writer reads with
	// a typed accessor. Without this, an INTEGER[] qual1 (what [40 for _ in range(n)]
	// yields) or a non-VARCHAR sequence1 reaches the writer and aborts with a fatal internal
	// error that invalidates the database. Must run after the options are
	// parsed: comment and sequence_index are only read when their option is on, and a
	// present-but-unread column of the wrong type is harmless.
	ValidateSequenceColumnTypes(result->indices, sql_types, /*validate_quals=*/true,
	                            /*validate_comment=*/result->include_comment,
	                            /*validate_sequence_index=*/result->id_as_sequence_index);

	// Parse FASTQ-specific qual_offset parameter
	if (!qual_offset_param.IsNull()) {
		int64_t offset = qual_offset_param.GetValue<int64_t>();
		if (offset != 33 && offset != 64) {
			throw BinderException("QUAL_OFFSET must be 33 or 64");
		}
		result->qual_offset = static_cast<uint8_t>(offset);
	}

	return result;
}

//===--------------------------------------------------------------------===//
// Global State (use shared implementation)
//===--------------------------------------------------------------------===//
using FastqCopyGlobalState = SequenceCopyGlobalState;

static unique_ptr<GlobalFunctionData> FastqCopyInitializeGlobal(ClientContext &context, FunctionData &bind_data,
                                                                const string &file_path) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	return SequenceCopyInitializeGlobal(context, fdata, file_path);
}

//===--------------------------------------------------------------------===//
// Finalize (use shared implementation)
//===--------------------------------------------------------------------===//
static void FastqCopyFinalize(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate_p) {
	auto &gstate = gstate_p.Cast<FastqCopyGlobalState>();
	SequenceCopyFinalize(gstate);
}

//===--------------------------------------------------------------------===//
// Batch API (SequencePrepareBatch / SequenceFlushBatch)
//===--------------------------------------------------------------------===//
struct FastqRecord {
	const char *id;
	idx_t id_size;
	const char *comment;
	idx_t comment_size;
	const char *seq;
	idx_t seq_size;
	const uint8_t *qual;
	idx_t qual_size;
	// For a quality error, which Write and Check find: the record's row and whether it is the mate.
	const SequenceChunkReader *reader;
	idx_t row;
	bool r2;
};

// The FORMAT of SequencePrepareBatch for FASTQ: a row's R1 record, then its R2 record when sequence2 and qual2 are set.
class FastqBatchFormat {
public:
	using Record = FastqRecord;

	explicit FastqBatchFormat(const FastqCopyBindData &fdata_p) : fdata(fdata_p), encoder(fdata_p.qual_offset) {
	}

	template <class EMIT>
	void Visit(DataChunk &input, SequenceBatchData &batch, EMIT &&emit) {
		auto &indices = fdata.indices;
		SequenceChunkReader reader(fdata, input);
		QualColumn qual1, qual2;
		qual1.Init(input.data[indices.qual1_idx]);
		if (fdata.has_r2_columns) {
			qual2.Init(input.data[indices.qual2_idx]);
		}

		for (idx_t row = 0; row < input.size(); row++) {
			FastqRecord record;
			reader.Header(row, record.id, record.id_size, record.comment, record.comment_size);
			const auto &seq1 = reader.Sequence1(row);
			auto qual1_row = qual1.data.sel->get_index(row);
			if (!qual1.data.validity.RowIsValid(qual1_row)) {
				reader.ThrowNull("qual1", row);
			}
			SetRead(record, seq1, qual1, qual1_row, reader, row, false);
			emit(SequenceStream::R1, record);

			// A record is paired iff sequence2 and qual2 are both non-NULL; exactly one NULL is malformed.
			const bool paired = reader.HasSequence2(row);
			if (fdata.has_r2_columns) {
				auto qual2_row = qual2.data.sel->get_index(row);
				if (paired != qual2.data.validity.RowIsValid(qual2_row)) {
					throw InvalidInputException(
					    "read_id '%s': sequence2 and qual2 must both be set (paired-end) or both be NULL (single-end)",
					    reader.ReadId(row));
				}
				if (paired) {
					SetRead(record, reader.Sequence2(row), qual2, qual2_row, reader, row, true);
				}
			}
			SequenceStream r2_stream;
			if (!reader.RouteR2(row, paired, batch, r2_stream)) {
				continue;
			}
			emit(r2_stream, record);
		}
	}

	// What FastqEncoder::Encode writes: "@" id [" " comment] "\n" sequence "\n+\n" qualities "\n".
	idx_t Size(const FastqRecord &record) const {
		return 1 + record.id_size + (record.comment_size > 0 ? 1 + record.comment_size : 0) + 1 + record.seq_size + 3 +
		       record.qual_size + 1;
	}

	void Write(MemoryStream &stream, const FastqRecord &record) {
		// The sink wraps WriteData so the encoder remains DuckDB-free (its unit tests link against the header alone).
		Encode(record,
		       [&stream](const char *data, std::size_t size) { stream.WriteData(const_data_ptr_cast(data), size); });
	}

	// For a batch formatted only when it is written: raises now, while the batch is prepared, the error Write would.
	void Check(const FastqRecord &record) {
		Encode(record, [](const char *, std::size_t) {});
	}

private:
	// A UTINYINT[] quality column of one chunk.
	struct QualColumn {
		UnifiedVectorFormat data;
		const uint8_t *scores = nullptr;
		const ValidityMask *score_validity = nullptr;
		// Whether any score in the chunk is NULL, found 64 scores at a time; only then are a read's scores checked.
		bool has_null_score = false;

		void Init(Vector &quals) {
			quals.ToUnifiedFormat(data);
			auto &child = ListVector::GetChildMutable(quals);
			scores = FlatVector::GetData<uint8_t>(child);
			score_validity = &FlatVector::Validity(child);
			has_null_score = !score_validity->CheckAllValid(ListVector::GetListSize(quals));
		}
	};

	// The encoder refuses a score that would encode past '~'; the error then names the record.
	void Encode(const FastqRecord &record, const FastqEncoder::Sink &sink) {
		try {
			encoder.Encode(sink, record.id, record.id_size, record.comment, record.comment_size, record.seq,
			               record.seq_size, record.qual, record.qual_size);
		} catch (const Exception &) {
			throw; // DuckDB's own, such as a full stream: not about the data
		} catch (const std::runtime_error &e) {
			throw InvalidInputException("%s for read_id '%s'%s", e.what(), record.reader->ReadId(record.row),
			                            record.r2 ? " (R2)" : "");
		}
	}

	// Points the record at a read's sequence (a string_t in the vector, which must outlive the record) and its
	// qualities (one row of a UTINYINT[] column), and checks them: one score per base and no NULL score. Errors name
	// the row's read_id, and R2 for the mate.
	void SetRead(FastqRecord &record, const string_t &seq, const QualColumn &quals, idx_t qual_row,
	             const SequenceChunkReader &reader, idx_t row, bool r2) const {
		auto &entry = UnifiedVectorFormat::GetData<list_entry_t>(quals.data)[qual_row];
		record.seq = seq.GetData();
		record.seq_size = seq.GetSize();
		record.qual = quals.scores + entry.offset;
		record.qual_size = entry.length;
		record.reader = &reader;
		record.row = row;
		record.r2 = r2;
		if (record.qual_size != record.seq_size) {
			throw InvalidInputException(
			    "Quality score length (%llu) does not match sequence length (%llu) for read_id '%s'%s",
			    record.qual_size, record.seq_size, reader.ReadId(row), r2 ? " (R2)" : "");
		}
		if (quals.has_null_score) {
			for (idx_t k = 0; k < entry.length; k++) {
				if (!quals.score_validity->RowIsValid(entry.offset + k)) {
					throw InvalidInputException("NULL quality score in %s (read_id '%s')", r2 ? "qual2" : "qual1",
					                            reader.ReadId(row));
				}
			}
		}
	}

	const FastqCopyBindData &fdata;
	FastqEncoder encoder; // reuses its quality buffer across records
};

// Reads only the bind data: see SequencePrepareBatch.
static unique_ptr<PreparedBatchData> FastqPrepareBatch(ClientContext &context, FunctionData &bind_data,
                                                       GlobalFunctionData &, unique_ptr<ColumnDataCollection> rows) {
	FastqBatchFormat format(bind_data.Cast<FastqCopyBindData>());
	return SequencePrepareBatch(context, format, std::move(rows));
}

static void FastqFlushBatch(ClientContext &context, FunctionData &bind_data, GlobalFunctionData &gstate,
                            PreparedBatchData &batch) {
	auto &fdata = bind_data.Cast<FastqCopyBindData>();
	FastqBatchFormat format(fdata);
	SequenceFlushBatch(context, gstate.Cast<FastqCopyGlobalState>(), batch.Cast<SequenceBatchData>(), format);
}

//===--------------------------------------------------------------------===//
// Register Function
//===--------------------------------------------------------------------===//
static constexpr const char *FASTQ_WRITER_FORMAT = "_miint_fastq";

static void FastqCopyOptions(ClientContext &context, CopyOptionsInput &input) {
	DeclareSequenceCopyOptions(input);
	DeclareCopyOption(input, "qual_offset", LogicalType::BIGINT);
}

static BoundStatement FastqCopyPlan(Binder &binder, CopyStatement &stmt) {
	return SequenceCopyPlan(binder, stmt, FASTQ_WRITER_FORMAT);
}

// FORMAT FASTQ: the options and the plan callback, which rebinds as the writer below (see SequenceCopyPlan).
CopyFunction CopyFastqFunction::GetFunction() {
	CopyFunction func("fastq");
	func.copy_options = FastqCopyOptions;
	func.plan = FastqCopyPlan;
	func.extension = "fastq";
	return func;
}

// FORMAT _miint_fastq: the writer. Named directly, it writes without the default batch cap.
static CopyFunction GetWriterFunction() {
	CopyFunction func(FASTQ_WRITER_FORMAT);
	func.copy_options = FastqCopyOptions;
	func.copy_to_bind = FastqCopyBind;
	func.copy_to_initialize_global = FastqCopyInitializeGlobal;
	func.copy_to_finalize = FastqCopyFinalize;
	func.execution_mode = SequenceCopyExecutionMode;
	func.prepare_batch = FastqPrepareBatch;
	func.flush_batch = FastqFlushBatch;
	func.extension = "fastq";
	return func;
}

void CopyFastqFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
	loader.RegisterFunction(GetWriterFunction());
}

} // namespace duckdb
