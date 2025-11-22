#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include "duckdb/common/serializer/buffered_file_writer.hpp"
#include "duckdb/common/serializer/memory_stream.hpp"
#include "duckdb/common/enums/file_compression_type.hpp"

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
	CopyFileHandle(FileSystem &fs, const string &path, FileCompressionType compression);
	~CopyFileHandle();

	void Write(const_data_ptr_t data, idx_t size);
	void WriteString(const string &data);
	void Close();

private:
	unique_ptr<BufferedFileWriter> file_writer;
	FileCompressionType compression;
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
void ValidatePairedEndParameters(bool is_paired, bool has_interleave_param, bool interleave, const string &file_path);
void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index);

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Structures (FASTA/FASTQ)
//===--------------------------------------------------------------------===//

// Base bind data for sequence formats (FASTA/FASTQ)
struct SequenceCopyBindData : public FunctionData {
	bool interleave = false;
	bool id_as_sequence_index = false;
	bool include_comment = false;
	FileCompressionType compression = FileCompressionType::UNCOMPRESSED;
	idx_t flush_size = DEFAULT_COPY_FLUSH_SIZE;
	string file_path;
	bool is_paired = false;
	vector<string> names;
	ColumnIndices indices;  // Pre-computed column indices

	// Subclasses must implement Copy() and Equals()
};

// Shared global state for sequence formats (100% identical for FASTA/FASTQ)
struct SequenceCopyGlobalState : public GlobalFunctionData {
	mutex lock;
	unique_ptr<CopyFileHandle> file_r1;
	unique_ptr<CopyFileHandle> file_r2;
	bool is_paired = false;
	bool interleave = false;
};

// Shared local state for sequence formats (100% identical for FASTA/FASTQ)
struct SequenceCopyLocalState : public LocalFunctionData {
	unique_ptr<FormatWriterState> writer_state_r1;
	unique_ptr<FormatWriterState> writer_state_r2; // For split paired-end
};

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Functions
//===--------------------------------------------------------------------===//

// Initialize global state for sequence formats
unique_ptr<GlobalFunctionData> SequenceCopyInitializeGlobal(ClientContext &context, const SequenceCopyBindData &fdata,
                                                             const string &file_path);

// Initialize local state for sequence formats
unique_ptr<LocalFunctionData> SequenceCopyInitializeLocal(ExecutionContext &context, const SequenceCopyBindData &fdata);

// Combine (flush) local buffers for sequence formats
void SequenceCopyCombine(const SequenceCopyBindData &fdata, SequenceCopyGlobalState &gstate,
                         SequenceCopyLocalState &lstate);

// Finalize sequence copy
void SequenceCopyFinalize(SequenceCopyGlobalState &gstate);

} // namespace duckdb
