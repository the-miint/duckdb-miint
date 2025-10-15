#pragma once

#include "duckdb/function/copy_function.hpp"
#include "duckdb/common/file_system.hpp"
#include <zlib.h>

namespace duckdb {

//===--------------------------------------------------------------------===//
// File Handle Wrapper (handles both compressed and uncompressed)
//===--------------------------------------------------------------------===//
class CopyFileHandle {
public:
	CopyFileHandle(FileSystem &fs, const string &path, bool use_compression);
	~CopyFileHandle();
	
	void Write(const string &data);
	void Close();
	
private:
	unique_ptr<FileHandle> file_handle;
	gzFile gz_file = nullptr;
	bool use_compression;
};

//===--------------------------------------------------------------------===//
// Common Helper Functions
//===--------------------------------------------------------------------===//
bool DetectCompression(const string &file_path, const Value &compression_param);
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
	bool use_compression = false;
	
	void ParseFromOptions(const case_insensitive_map_t<vector<Value>> &options,
	                      const string &file_path);
};

//===--------------------------------------------------------------------===//
// Common Validation Functions
//===--------------------------------------------------------------------===//
void ValidateRequiredColumns(bool has_read_id, bool has_sequence1, const string &format_name);
void ValidatePairedEndParameters(bool is_paired, bool has_interleave_param, bool interleave,
                                 const string &file_path);
void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index);

} // namespace duckdb
