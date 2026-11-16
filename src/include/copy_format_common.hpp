#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/serializer/buffered_file_writer.hpp"
#include "duckdb/common/serializer/memory_stream.hpp"
#include "duckdb/common/enums/file_compression_type.hpp"
#include "duckdb/common/vector/flat_vector.hpp"
#include "duckdb/common/types/column/column_data_collection.hpp"
#include "duckdb/planner/bound_statement.hpp"
// ResolveSequenceRecordId is header-inline (hot path) so it needs the id codec,
// UUID::ToString, AllowedIdTypeList, exception types, and UnifiedVectorFormat --
// id_column_utils.hpp transitively provides all of them.
#include "id_column_utils.hpp"
#include "copy_option_parsing.hpp"

// htslib's BGZF handle; defined in <htslib/bgzf.h>, only included in the .cpp.
struct BGZF;

namespace duckdb {

class Binder;
class CopyStatement;

//===--------------------------------------------------------------------===//
// Constants
//===--------------------------------------------------------------------===//
constexpr idx_t DEFAULT_COPY_FLUSH_SIZE = 1024 * 1024; // 1MB default buffer size
// The batch cap FASTA/FASTQ add when the statement gives none (ROW_GROUP_SIZE, ROW_GROUP_SIZE_BYTES). DuckDB flushes a
// batch once it holds at least this many rows, or once it reaches this many bytes; only the byte cap cuts an input
// chunk. The bytes bound the memory a batch of long records takes; the rows gather small input chunks. 8 MB rather
// than more: with order kept one thread copies, formats and writes every batch, and FORMAT FASTA of a genome set from
// 2-record Parquet files ran 7% slower with 64 MB batches than with 8 MB ones (8 threads, three runs each, spread 3%);
// Revio HiFi reads agreed. Smaller batches also hold less memory.
constexpr idx_t SEQUENCE_DEFAULT_BATCH_ROWS = 64;
constexpr const char *SEQUENCE_DEFAULT_BATCH_BYTES = "8MB";
// A batch whose output is larger than this is formatted in flush_batch and written out each time
// DEFAULT_COPY_FLUSH_SIZE has gathered, rather than formatted into one buffer when it is prepared: that buffer would
// sit outside DuckDB's memory_limit. A record larger than DEFAULT_COPY_FLUSH_SIZE is still held whole.
constexpr idx_t SEQUENCE_DEFER_FORMAT_BYTES = 64ULL * 1024 * 1024;
// FASTQ/FASTA gzip level when COMPRESSION_LEVEL is not given. COPY usually writes whole read sets,
// such as every read that survives a host filter; there, level 1 compressed about three times as fast
// as level 6 (wall time of a COPY from read_fastx; 10 M pairs of 150 bp reads, split R1/R2 output,
// 8 threads, libdeflate) and wrote files about 9% larger.
constexpr int SEQUENCE_DEFAULT_GZIP_LEVEL = 1;

//===--------------------------------------------------------------------===//
// File Handle Wrapper (handles both compressed and uncompressed)
//===--------------------------------------------------------------------===//
class CopyFileHandle {
public:
	// compression_threads > 1 enables htslib bgzf multithreaded compression for
	// LOCAL gzip output (bgzf is gzip-compatible and bgzf_mt preserves block
	// order). Uncompressed output and remote (scheme://) gzip targets keep the
	// BufferedFileWriter path. compression_level applies to local gzip output only; it has no
	// default, so no caller silently gets a level it did not choose.
	CopyFileHandle(FileSystem &fs, const string &path, FileCompressionType compression, int compression_threads,
	               int compression_level);
	~CopyFileHandle();

	void Write(const_data_ptr_t data, idx_t size);
	void WriteString(const string &data);
	void Close();

private:
	// Exactly one is active for the handle's lifetime: bgzf_file for local gzip (htslib bgzf),
	// file_writer for uncompressed output and the remote-gzip fallback.
	unique_ptr<BufferedFileWriter> file_writer;
	struct ::BGZF *bgzf_file = nullptr;
};

//===--------------------------------------------------------------------===//
// Format Writer State (a buffer written out every flush_size bytes)
//===--------------------------------------------------------------------===//
struct FormatWriterState {
	FormatWriterState(ClientContext &context, idx_t flush_size);
	~FormatWriterState();

	void Reset();

	idx_t flush_size;
	unique_ptr<MemoryStream> stream;
};

//===--------------------------------------------------------------------===//
// Common Helper Functions
//===--------------------------------------------------------------------===//
// COMPRESSION for the formats that write text: the value given, or, with none (a NULL compression_param),
// gzip for a path ending in .gz. A format whose COMPRESSION takes other values checks its own.
FileCompressionType DetectCompressionType(const string &file_path, const Value &compression_param);
// Whether a COMPRESSION value names gzip: 'gzip' or 'gz', in any case.
bool IsGzipCompressionName(const string &name);
// A COMPRESSION value of a format that writes gzip or nothing: true for 'gzip' or 'gz', false for 'none' or
// 'uncompressed', in any case; any other value is refused.
bool CompressionValueIsGzip(const Value &compression);
string SubstituteOrientation(const string &path, const string &orientation);
bool HasOrientationPlaceholder(const string &path);

//===--------------------------------------------------------------------===//
// Common Column Index Finder
//===--------------------------------------------------------------------===//
struct ColumnIndices {
	idx_t read_id_idx = DConstants::INVALID_INDEX;
	idx_t sequence_index_idx = DConstants::INVALID_INDEX;
	idx_t comment_idx = DConstants::INVALID_INDEX;
	idx_t sequence1_idx = DConstants::INVALID_INDEX;
	idx_t sequence2_idx = DConstants::INVALID_INDEX;
	idx_t qual1_idx = DConstants::INVALID_INDEX;
	idx_t qual2_idx = DConstants::INVALID_INDEX;

	void FindIndices(const vector<string> &names);
};

//===--------------------------------------------------------------------===//
// Common Parameter Parsing
//===--------------------------------------------------------------------===//
struct CommonCopyParameters {
	bool interleave = false;
	bool id_as_sequence_index = false;
	bool include_comment = false;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	int compression_level = SEQUENCE_DEFAULT_GZIP_LEVEL;

	void ParseFromOptions(const identifier_map_t<vector<Value>> &options, const string &file_path);
};

//! The options COPY FORMAT FASTQ and FASTA share, for their copy_options.
void DeclareSequenceCopyOptions(CopyOptionsInput &input);

//! Files DuckDB names itself (PARTITION_BY, PER_THREAD_OUTPUT) end in .gz when the output is gzip, as DuckDB's CSV
//! writer names its own; an extension that already ends in .gz is kept.
void AddGzipFileExtension(CopyFunctionBindInput &input, FileCompressionType compression);

//! COMPRESSION_LEVEL for output written through CopyFileHandle. Refused where it would be ignored:
//! uncompressed output, and remote targets, which DuckDB's gzip writer compresses at its own level.
int GzipCompressionLevel(const vector<Value> &values, FileCompressionType compression, const string &file_path);

//===--------------------------------------------------------------------===//
// Common Validation Functions
//===--------------------------------------------------------------------===//
void ValidateRequiredColumns(bool has_read_id, bool has_sequence1, const string &format_name);

// Validate the storage type of every column the FASTA/FASTQ writers read with a typed
// accessor: sequence1/sequence2/comment via GetData<string_t>, sequence_index via
// GetData<int64_t>, and (FASTQ only) qual1/qual2 via FlatVector::GetData<uint8_t> on the
// list child. Those reads perform no type dispatch, so a mismatched column type is not a
// wrong answer -- it trips a DuckDB assertion, which is an internal error that
// invalidates the database for the rest of the session, burying the real cause under a
// cascade of "database has been invalidated" failures.
//
// Validating at bind turns that into a clean BinderException naming the column and the
// expected type. This mirrors the read_id check each writer already performs and the
// checks COPY FORMAT UBAM makes in its bind (UBAMCopyBind); the message wording is kept
// identical to UBAM's so the two never drift.
//
// Only columns the given writer actually reads for this query are validated, because a
// column that is present but never read is harmless and rejecting it would break queries
// that work today (e.g. a stray non-BIGINT sequence_index surviving a SELECT *):
//   - validate_quals          false for FASTA, which ignores qual columns entirely
//   - validate_comment        pass include_comment; comment is only read when it is set
//   - validate_sequence_index pass id_as_sequence_index; likewise
// sequence1/sequence2 are always read when present, so they are always validated.
// Consequently this must be called after the COPY options have been parsed.
//
// read_id is not checked here -- callers validate it via IsAllowedIdType because they also
// need to capture the resolved type.
void ValidateSequenceColumnTypes(const ColumnIndices &indices, const vector<LogicalType> &sql_types,
                                 bool validate_quals, bool validate_comment, bool validate_sequence_index);
// Resolves the output structure from the path + INTERLEAVE option. The presence of the
// {ORIENTATION} placeholder is the split-vs-single switch: when present, R1/R2 are written to
// separate files; when absent, output is a single file (interleaved iff INTERLEAVE=true).
// Sets out_split_output. Throws on the one contradictory combination ({ORIENTATION}+INTERLEAVE=true).
void ResolveOutputMode(const string &file_path, bool interleave, bool &out_split_output);
void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index);

// Resolve a read_id cell to (out_ptr, out_size), dispatching on the SQL type
// captured at bind (VARCHAR / BIGINT / UUID). VARCHAR is zero-copy into the DuckDB
// string heap; BIGINT/UUID are stringified into `buf` (a caller-owned buffer reused
// across rows). The caller is responsible for NULL handling before calling. Shared
// by the FASTA and FASTQ writers so id-type support lives in one place. Inline: this
// is on the per-record write path, so the VARCHAR fast path stays call-free.
inline void ResolveSequenceRecordId(const UnifiedVectorFormat &read_id_data, idx_t row, const LogicalType &read_id_type,
                                    string &buf, const char *&out_ptr, idx_t &out_size) {
	auto idx = read_id_data.sel->get_index(row);
	switch (read_id_type.id()) {
	case LogicalTypeId::VARCHAR: {
		// Zero-copy: point directly into the DuckDB string heap (hot path).
		auto data = UnifiedVectorFormat::GetData<string_t>(read_id_data);
		out_ptr = data[idx].GetData();
		out_size = data[idx].GetSize();
		return;
	}
	case LogicalTypeId::BIGINT: {
		auto data = UnifiedVectorFormat::GetData<int64_t>(read_id_data);
		buf = ::miint::FormatIdFromInt64(data[idx]);
		out_ptr = buf.data();
		out_size = buf.size();
		return;
	}
	case LogicalTypeId::UUID: {
		auto data = UnifiedVectorFormat::GetData<hugeint_t>(read_id_data);
		buf = UUID::ToString(data[idx]);
		out_ptr = buf.data();
		out_size = buf.size();
		return;
	}
	default:
		// Bind validates the type via IsAllowedIdType, so this is unreachable.
		throw InternalException("ResolveSequenceRecordId: unsupported read_id type '%s' (must be %s)",
		                        read_id_type.ToString(), AllowedIdTypeList());
	}
}

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Structures (FASTA/FASTQ)
//===--------------------------------------------------------------------===//

// Base bind data for sequence formats (FASTA/FASTQ)
struct SequenceCopyBindData : public FunctionData {
	bool interleave = false; // INTERLEAVE=true: write R1/R2 records alternating into one file
	bool id_as_sequence_index = false;
	bool include_comment = false;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	int compression_level = SEQUENCE_DEFAULT_GZIP_LEVEL; // COMPRESSION_LEVEL, for local gzip output
	string file_path;
	// Whether the input *schema* carries R2 columns (sequence2 [+ qual2 for FASTQ]). This does
	// NOT mean the data is paired -- read_fastx always emits these columns (NULL for single-end).
	// Whether a given record is paired is decided per-row at write time from R2 NULL-ness.
	bool has_r2_columns = false;
	// Whether the output path contains the {ORIENTATION} placeholder, i.e. split R1/R2 files.
	bool split_output = false;
	vector<string> names;
	ColumnIndices indices; // Pre-computed column indices
	// SQL type of the read_id column, captured at bind (VARCHAR / BIGINT / UUID).
	// The writers dispatch on this so BIGINT/UUID ids are stringified rather than
	// (incorrectly) read as string_t -- see ResolveSequenceRecordId.
	LogicalType read_id_type;

	// Subclasses must implement Copy() and Equals()
};

// Shared global state for sequence formats (100% identical for FASTA/FASTQ)
struct SequenceCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	unique_ptr<CopyFileHandle> file_r1;
	unique_ptr<CopyFileHandle> file_r2; // split mode only; created lazily on first R2 write
	// Deferred R2-file creation parameters (split mode). The R2 file is only created once a
	// non-NULL R2 record is actually written, so an all-single-end dataset never produces an
	// empty R2 file even when {ORIENTATION} is used.
	FileSystem *fs = nullptr;
	string r2_path;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	// Worker threads for bgzf gzip compression (follows DuckDB's thread count) and the gzip level;
	// captured at global-init so the lazily-created R2 file uses the same values.
	int compression_threads = 1;
	int compression_level = SEQUENCE_DEFAULT_GZIP_LEVEL;
	// Consistency tracking across all threads: a single COPY must be either all single-end or
	// all paired-end. Both set true -> inconsistent input (checked at finalize).
	bool saw_paired = false;
	bool saw_single = false;
};

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Functions
//===--------------------------------------------------------------------===//

// Initialize global state for sequence formats
unique_ptr<GlobalFunctionData> SequenceCopyInitializeGlobal(ClientContext &context, const SequenceCopyBindData &fdata,
                                                            const string &file_path);

// Finalize sequence copy
void SequenceCopyFinalize(SequenceCopyGlobalState &gstate);

//===--------------------------------------------------------------------===//
// Sequence Copy on DuckDB's Batch API (FASTA/FASTQ)
//===--------------------------------------------------------------------===//
// FORMAT FASTA/FASTQ are two registrations each: the public one ("fasta") declares the options, the extension and a
// plan callback, SequenceCopyPlan, which rebinds the statement as the writer ("_miint_fasta"). See
// docs/internals/architecture.md, "Record writers on DuckDB's batch API".

// REGULAR when insertion order is kept, PARALLEL when it is not; never BATCH (see the architecture notes).
CopyFunctionExecutionMode SequenceCopyExecutionMode(bool preserve_insertion_order, bool supports_batch_index);

// The public registration's plan callback: adds ROW_GROUP_SIZE and ROW_GROUP_SIZE_BYTES where absent, then binds the
// statement as `writer_format`.
BoundStatement SequenceCopyPlan(Binder &binder, CopyStatement &stmt, const string &writer_format);

// Where a record's text goes, named for the {ORIENTATION} files: R1 (which also takes the R2 records of interleaved
// output) or R2 (split output only).
enum class SequenceStream : uint8_t { R1, R2 };

// A batch as prepare_batch leaves it for flush_batch: its text for each file, formatted, or (when that would exceed
// SEQUENCE_DEFER_FORMAT_BYTES) the rows themselves. The paired/single flags travel with the batch to the global state.
struct SequenceBatchData : public PreparedBatchData {
	// Each file's text: a buffer of exactly its size, written through a stream that does not own it (an owning
	// MemoryStream needs a power-of-two capacity).
	AllocatedData r1_data;
	AllocatedData r2_data;
	unique_ptr<MemoryStream> r1;
	unique_ptr<MemoryStream> r2;
	unique_ptr<ColumnDataCollection> deferred;
	bool saw_paired = false;
	bool saw_single = false;
};

// Writes `data` to the R1 or R2 file, creating the R2 file on its first write. Caller holds gstate.lock.
void WriteSequenceStream(SequenceCopyGlobalState &gstate, SequenceStream stream, MemoryStream &data);

// The columns FASTA and FASTQ read alike, for one chunk: the id (read_id, or sequence_index with ID_AS_SEQUENCE_INDEX),
// comment, sequence1 and sequence2, with the checks and the R2 routing the two writers share. Each writer's Visit
// adds its own columns (FASTQ's qualities).
class SequenceChunkReader {
public:
	SequenceChunkReader(const SequenceCopyBindData &fdata_p, DataChunk &input) : fdata(fdata_p) {
		auto &indices = fdata.indices;
		input.data[indices.read_id_idx].ToUnifiedFormat(read_id_data);
		if (fdata.id_as_sequence_index) {
			input.data[indices.sequence_index_idx].ToUnifiedFormat(sequence_index_data);
		}
		has_comment = fdata.include_comment && indices.comment_idx != DConstants::INVALID_INDEX;
		if (has_comment) {
			input.data[indices.comment_idx].ToUnifiedFormat(comment_data);
		}
		input.data[indices.sequence1_idx].ToUnifiedFormat(seq1_data);
		if (fdata.has_r2_columns) {
			input.data[indices.sequence2_idx].ToUnifiedFormat(seq2_data);
		}
	}

	// The row's header: the id written for it and its comment (none when there is no comment, or it is NULL). The id
	// stays valid until the next call. A NULL read_id, or a NULL sequence_index when that is the id, is an error.
	void Header(idx_t row, const char *&id, idx_t &id_size, const char *&comment, idx_t &comment_size) {
		if (!read_id_data.validity.RowIsValid(read_id_data.sel->get_index(row))) {
			throw InvalidInputException("NULL value in read_id column");
		}
		if (fdata.id_as_sequence_index) {
			auto seq_idx_row = sequence_index_data.sel->get_index(row);
			if (!sequence_index_data.validity.RowIsValid(seq_idx_row)) {
				ThrowNull("sequence_index", row);
			}
			id_buf = to_string(UnifiedVectorFormat::GetData<int64_t>(sequence_index_data)[seq_idx_row]);
			id = id_buf.data();
			id_size = id_buf.size();
		} else {
			ResolveSequenceRecordId(read_id_data, row, fdata.read_id_type, id_buf, id, id_size);
		}
		comment = nullptr;
		comment_size = 0;
		if (has_comment) {
			auto comment_row = comment_data.sel->get_index(row);
			if (comment_data.validity.RowIsValid(comment_row)) {
				auto comment_strings = UnifiedVectorFormat::GetData<string_t>(comment_data);
				comment = comment_strings[comment_row].GetData();
				comment_size = comment_strings[comment_row].GetSize();
			}
		}
	}

	// sequence1; NULL is an error. A reference into the vector: a string of up to 12 bytes is stored inside the
	// string_t, so a copy's GetData() would not outlive the copy.
	const string_t &Sequence1(idx_t row) const {
		auto seq1_row = seq1_data.sel->get_index(row);
		if (!seq1_data.validity.RowIsValid(seq1_row)) {
			ThrowNull("sequence1", row);
		}
		return UnifiedVectorFormat::GetData<string_t>(seq1_data)[seq1_row];
	}

	// Whether the row has a sequence2 (false without R2 columns).
	bool HasSequence2(idx_t row) const {
		return fdata.has_r2_columns && seq2_data.validity.RowIsValid(seq2_data.sel->get_index(row));
	}

	// sequence2, for a row HasSequence2 accepted. A reference into the vector, as Sequence1.
	const string_t &Sequence2(idx_t row) const {
		auto seq2_row = seq2_data.sel->get_index(row);
		D_ASSERT(fdata.has_r2_columns && seq2_data.validity.RowIsValid(seq2_row));
		return UnifiedVectorFormat::GetData<string_t>(seq2_data)[seq2_row];
	}

	// Records whether the row is paired (`paired`, as the writer decides it) and returns in `stream` where its R2
	// record goes. Returns false for a single-end row. A paired row with nowhere to put R2 (single-file output) is an
	// error.
	bool RouteR2(idx_t row, bool paired, SequenceBatchData &batch, SequenceStream &stream) const {
		if (!fdata.has_r2_columns) {
			return false;
		}
		if (!paired) {
			batch.saw_single = true;
			return false;
		}
		batch.saw_paired = true;
		if (!fdata.split_output && !fdata.interleave) {
			throw InvalidInputException(
			    "Record read_id '%s' has paired-end data (sequence2 is set) but the output is a single file. Use the "
			    "{ORIENTATION} placeholder in the output path to write split R1/R2 files, or set INTERLEAVE=true to "
			    "interleave R1/R2 into one file.",
			    ReadId(row));
		}
		stream = fdata.split_output ? SequenceStream::R2 : SequenceStream::R1;
		return true;
	}

	// The row's read_id as text, for a data error to name.
	string ReadId(idx_t row) const;

	// The data error for a NULL in `column`, naming the row's read_id.
	[[noreturn]] void ThrowNull(const char *column, idx_t row) const {
		throw InvalidInputException("NULL value in %s column (read_id '%s')", column, ReadId(row));
	}

private:
	const SequenceCopyBindData &fdata;
	UnifiedVectorFormat read_id_data, sequence_index_data, comment_data, seq1_data, seq2_data;
	bool has_comment;
	string id_buf; // a BIGINT/UUID read_id or a sequence_index as text
};

// FASTA and FASTQ differ only in how rows become records and records become text, which FORMAT supplies:
//   FORMAT::Record              one record's fields
//   Visit(chunk, batch, emit)   checks each row, calls emit(stream, record) for its records in output order, and
//                               sets batch.saw_paired / saw_single
//   Size(record)                the record's exact length as text
//   Write(stream, record)       appends the record's text; may raise a data error (FASTQ: a score past '~')
//   Check(record)               raises the data error Write would, without writing
// prepare_batch runs concurrently and may be handed another output file's global state (PARTITION_BY), so it reads only
// the bind data.
template <class FORMAT>
unique_ptr<PreparedBatchData> SequencePrepareBatch(ClientContext &context, FORMAT &format,
                                                   unique_ptr<ColumnDataCollection> collection) {
	auto batch = make_uniq<SequenceBatchData>();
	// Pass 1 checks every row and sizes each file's text. Data errors are raised in prepare_batch, on the preparing
	// thread: here, or by Write in pass 2.
	idx_t r1_size = 0;
	idx_t r2_size = 0;
	for (auto &chunk : collection->Chunks()) {
		format.Visit(chunk, *batch, [&](SequenceStream stream, const typename FORMAT::Record &record) {
			(stream == SequenceStream::R1 ? r1_size : r2_size) += format.Size(record);
		});
	}
	if (r1_size + r2_size > SEQUENCE_DEFER_FORMAT_BYTES) {
		// Write runs in flush_batch for this batch, so the errors it would raise are raised now.
		for (auto &chunk : collection->Chunks()) {
			format.Visit(chunk, *batch,
			             [&](SequenceStream, const typename FORMAT::Record &record) { format.Check(record); });
		}
		batch->deferred = std::move(collection);
		return std::move(batch);
	}
	// Pass 2 formats into buffers of exactly that size.
	auto &allocator = Allocator::Get(context);
	if (r1_size > 0) {
		batch->r1_data = allocator.Allocate(r1_size);
		batch->r1 = make_uniq<MemoryStream>(batch->r1_data.get(), r1_size);
	}
	if (r2_size > 0) {
		batch->r2_data = allocator.Allocate(r2_size);
		batch->r2 = make_uniq<MemoryStream>(batch->r2_data.get(), r2_size);
	}
	for (auto &chunk : collection->Chunks()) {
		format.Visit(chunk, *batch, [&](SequenceStream stream, const typename FORMAT::Record &record) {
			auto &out = stream == SequenceStream::R1 ? batch->r1 : batch->r2;
			D_ASSERT(out); // pass 1 gave this file's text a size, so it has a buffer
			format.Write(*out, record);
		});
	}
	// Size() and Write() must agree. Text longer than its size fails in WriteData; this catches text that is shorter.
	if ((batch->r1 && batch->r1->GetPosition() != r1_size) || (batch->r2 && batch->r2->GetPosition() != r2_size)) {
		throw InternalException("FASTA/FASTQ record text is shorter than its computed size");
	}
	return std::move(batch);
}

// flush_batch: writes the batch's text to the R1 and R2 files under the global lock. DuckDB calls it for one batch at a
// time per file, in input order when order is kept.
template <class FORMAT>
void SequenceFlushBatch(ClientContext &context, SequenceCopyGlobalState &gstate, SequenceBatchData &batch,
                        FORMAT &format) {
	lock_guard<mutex> glock(gstate.lock);
	if (batch.deferred) {
		FormatWriterState r1(context, DEFAULT_COPY_FLUSH_SIZE);
		FormatWriterState r2(context, DEFAULT_COPY_FLUSH_SIZE);
		for (auto &chunk : batch.deferred->Chunks()) {
			format.Visit(chunk, batch, [&](SequenceStream stream, const typename FORMAT::Record &record) {
				auto &writer = stream == SequenceStream::R1 ? r1 : r2;
				format.Write(*writer.stream, record);
				if (writer.stream->GetPosition() >= writer.flush_size) {
					WriteSequenceStream(gstate, stream, *writer.stream);
					writer.Reset();
				}
			});
		}
		WriteSequenceStream(gstate, SequenceStream::R1, *r1.stream);
		WriteSequenceStream(gstate, SequenceStream::R2, *r2.stream);
	} else {
		if (batch.r1) {
			WriteSequenceStream(gstate, SequenceStream::R1, *batch.r1);
		}
		if (batch.r2) {
			WriteSequenceStream(gstate, SequenceStream::R2, *batch.r2);
		}
	}
	gstate.saw_paired |= batch.saw_paired;
	gstate.saw_single |= batch.saw_single;
}

} // namespace duckdb
