#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/serializer/buffered_file_writer.hpp"
#include "duckdb/common/serializer/memory_stream.hpp"
#include "duckdb/common/enums/file_compression_type.hpp"

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
// Resolves the output structure from the path + INTERLEAVE option. The presence of the
// {ORIENTATION} placeholder is the split-vs-single switch: when present, R1/R2 are written to
// separate files; when absent, output is a single file (interleaved iff INTERLEAVE=true).
// Sets out_split_output. Throws on the one contradictory combination ({ORIENTATION}+INTERLEAVE=true).
void ResolveOutputMode(const string &file_path, bool interleave, bool &out_split_output);
void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index);

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
