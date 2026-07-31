#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/serializer/buffered_file_writer.hpp"
#include "duckdb/common/serializer/memory_stream.hpp"
#include "duckdb/common/enums/file_compression_type.hpp"
// ResolveSequenceRecordId is header-inline (hot path) so it needs the id codec,
// UUID::ToString, AllowedIdTypeList, exception types, and UnifiedVectorFormat --
// id_column_utils.hpp transitively provides all of them.
#include "id_column_utils.hpp"

// htslib's BGZF handle; defined in <htslib/bgzf.h>, only included in the .cpp.
struct BGZF;

namespace duckdb {

//===--------------------------------------------------------------------===//
// Constants
//===--------------------------------------------------------------------===//
constexpr idx_t DEFAULT_COPY_FLUSH_SIZE = 1024 * 1024; // 1MB default buffer size

//===--------------------------------------------------------------------===//
// File Handle Wrapper (handles both compressed and uncompressed)
//===--------------------------------------------------------------------===//
class CopyFileHandle {
public:
	// compression_threads > 1 enables htslib bgzf multithreaded compression for
	// LOCAL gzip output (bgzf is gzip-compatible and bgzf_mt preserves block
	// order). Uncompressed output and remote (scheme://) gzip targets keep the
	// BufferedFileWriter path.
	CopyFileHandle(FileSystem &fs, const string &path, FileCompressionType compression, int compression_threads = 1);
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
// Format Writer State (Per-Thread Buffer)
//===--------------------------------------------------------------------===//
struct FormatWriterState {
	FormatWriterState(ClientContext &context, idx_t flush_size);
	~FormatWriterState();

	void Reset();

	idx_t flush_size;
	unique_ptr<MemoryStream> stream;
	bool written_anything = false;
};

//===--------------------------------------------------------------------===//
// Format Writer Helper Functions
//===--------------------------------------------------------------------===//
void FlushFormatBuffer(FormatWriterState &local_state, CopyFileHandle &file, mutex &lock);

//===--------------------------------------------------------------------===//
// Common Helper Functions
//===--------------------------------------------------------------------===//
FileCompressionType DetectCompressionType(const string &file_path, const Value &compression_param);
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
	idx_t flush_size = DEFAULT_COPY_FLUSH_SIZE;

	void ParseFromOptions(const case_insensitive_map_t<vector<Value>> &options, const string &file_path);
};

//===--------------------------------------------------------------------===//
// Common Validation Functions
//===--------------------------------------------------------------------===//
void ValidateRequiredColumns(bool has_read_id, bool has_sequence1, const string &format_name);

// Validate the storage type of every column the FASTA/FASTQ sinks read with a TYPED
// accessor: sequence1/sequence2/comment via GetData<string_t>, sequence_index via
// GetData<int64_t>, and (FASTQ only) qual1/qual2 via FlatVector::GetData<uint8_t> on the
// list child. Those reads perform no type dispatch, so a mismatched column type is not a
// wrong answer -- it trips a DuckDB assertion, which is an INTERNAL error that
// INVALIDATES THE DATABASE for the rest of the session, burying the real cause under a
// cascade of "database has been invalidated" failures (issue #191).
//
// Validating at bind turns that into a clean BinderException naming the column and the
// expected type. This mirrors the read_id check each writer already performs and the
// identical checks COPY FORMAT UBAM does (src/copy_ubam.cpp:203-213) -- FASTA/FASTQ were
// the outliers. Message wording is kept identical to UBAM's so the two never drift.
//
// Only columns the given writer actually reads for THIS query are validated, because a
// column that is present but never read is harmless and rejecting it would break queries
// that work today (e.g. a stray non-BIGINT sequence_index surviving a SELECT *):
//   - validate_quals          false for FASTA, which ignores qual columns entirely
//   - validate_comment        pass include_comment; comment is only read when it is set
//   - validate_sequence_index pass id_as_sequence_index; likewise
// sequence1/sequence2 are always read when present, so they are always validated.
// Consequently this must be called AFTER the COPY options have been parsed.
//
// read_id is NOT checked here -- callers validate it via IsAllowedIdType because they also
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
// by the FASTA and FASTQ sinks so id-type support lives in one place. Inline: this
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
	idx_t flush_size = DEFAULT_COPY_FLUSH_SIZE;
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
	// The sink dispatches on this so BIGINT/UUID ids are stringified rather than
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
	// Worker threads for bgzf gzip compression (follows DuckDB's thread count);
	// captured at global-init so the lazily-created R2 file uses the same value.
	int compression_threads = 1;
	// Consistency tracking across all threads: a single COPY must be either all single-end or
	// all paired-end. Both set true -> inconsistent input (checked at finalize).
	bool saw_paired = false;
	bool saw_single = false;
};

// Shared local state for sequence formats (100% identical for FASTA/FASTQ)
struct SequenceCopyLocalState : public LocalFunctionData {
	unique_ptr<FormatWriterState> writer_state_r1;
	unique_ptr<FormatWriterState> writer_state_r2; // For split paired-end
	bool saw_paired = false;                       // this thread wrote >=1 paired record
	bool saw_single = false;                       // this thread wrote >=1 single-end record
};

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Functions
//===--------------------------------------------------------------------===//

// Initialize global state for sequence formats
unique_ptr<GlobalFunctionData> SequenceCopyInitializeGlobal(ClientContext &context, const SequenceCopyBindData &fdata,
                                                            const string &file_path);

// Initialize local state for sequence formats
unique_ptr<LocalFunctionData> SequenceCopyInitializeLocal(ExecutionContext &context, const SequenceCopyBindData &fdata);

// Flush a local R2 buffer in split mode, lazily creating the R2 file on first non-empty flush.
void FlushR2Buffer(FormatWriterState &local_state, SequenceCopyGlobalState &gstate);

// Combine (flush) local buffers for sequence formats
void SequenceCopyCombine(const SequenceCopyBindData &fdata, SequenceCopyGlobalState &gstate,
                         SequenceCopyLocalState &lstate);

// Finalize sequence copy
void SequenceCopyFinalize(SequenceCopyGlobalState &gstate);

} // namespace duckdb
