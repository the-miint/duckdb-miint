#include "copy_format_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"

namespace duckdb {

//===--------------------------------------------------------------------===//
// FormatWriterState Implementation
//===--------------------------------------------------------------------===//
FormatWriterState::FormatWriterState(ClientContext &context, idx_t flush_size_p)
    : flush_size(flush_size_p), written_anything(false) {
	stream = make_uniq<MemoryStream>();
}

FormatWriterState::~FormatWriterState() {
	// MemoryStream auto-cleans
}

void FormatWriterState::Reset() {
	stream->Rewind();
	written_anything = false;
}

//===--------------------------------------------------------------------===//
// Format Writer Helper Functions
//===--------------------------------------------------------------------===//
void FlushFormatBuffer(FormatWriterState &local_state, CopyFileHandle &file, mutex &lock) {
	if (!local_state.written_anything) {
		return;
	}
	
	lock_guard<mutex> glock(lock);
	
	// Write accumulated buffer to file
	file.Write(local_state.stream->GetData(), local_state.stream->GetPosition());
	
	// Reset local buffer
	local_state.Reset();
}

//===--------------------------------------------------------------------===//
// CopyFileHandle Implementation
//===--------------------------------------------------------------------===//
CopyFileHandle::CopyFileHandle(FileSystem &fs, const string &path, FileCompressionType compression_p)
    : compression(compression_p) {
	auto flags = FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW | 
	             FileLockType::WRITE_LOCK | compression;
	
	// BufferedFileWriter handles both file opening and buffering
	file_writer = make_uniq<BufferedFileWriter>(fs, path, flags);
}

CopyFileHandle::~CopyFileHandle() {
	Close();
}

void CopyFileHandle::Write(const_data_ptr_t data, idx_t size) {
	if (file_writer) {
		file_writer->WriteData(data, size);
	}
}

void CopyFileHandle::WriteString(const string &data) {
	Write(const_data_ptr_cast(data.c_str()), data.size());
}

void CopyFileHandle::Close() {
	if (file_writer) {
		file_writer->Close();
		file_writer.reset();
	}
}

//===--------------------------------------------------------------------===//
// Common Helper Functions
//===--------------------------------------------------------------------===//
FileCompressionType DetectCompressionType(const string &file_path, const Value &compression_param) {
	// Explicit parameter takes precedence
	if (!compression_param.IsNull()) {
		string comp_str = StringUtil::Lower(compression_param.ToString());
		if (comp_str == "gzip" || comp_str == "gz") {
			return FileCompressionType::GZIP;
		} else if (comp_str == "zstd" || comp_str == "zst") {
			return FileCompressionType::ZSTD;
		} else if (comp_str == "none") {
			return FileCompressionType::UNCOMPRESSED;
		} else {
			throw InvalidInputException("compression must be 'gzip', 'gz', 'zstd', 'zst', or 'none'");
		}
	}
	
	// Auto-detect from extension
	if (file_path.size() >= 3 && file_path.substr(file_path.size() - 3) == ".gz") {
		return FileCompressionType::GZIP;
	}
	if (file_path.size() >= 4 && file_path.substr(file_path.size() - 4) == ".zst") {
		return FileCompressionType::ZSTD;
	}
	
	return FileCompressionType::UNCOMPRESSED;
}

string SubstituteOrientation(const string &path, const string &orientation) {
	size_t pos = path.find("{ORIENTATION}");
	if (pos == string::npos) {
		return path;
	}
	string result = path;
	result.replace(pos, 13, orientation);
	return result;
}

bool HasOrientationPlaceholder(const string &path) {
	return path.find("{ORIENTATION}") != string::npos;
}

//===--------------------------------------------------------------------===//
// ColumnIndices Implementation
//===--------------------------------------------------------------------===//
void ColumnIndices::FindIndices(const vector<string> &names) {
	for (idx_t i = 0; i < names.size(); i++) {
		auto &name = names[i];
		if (name == "read_id") read_id_idx = i;
		else if (name == "sequence_index") sequence_index_idx = i;
		else if (name == "comment") comment_idx = i;
		else if (name == "sequence1") sequence1_idx = i;
		else if (name == "sequence2") sequence2_idx = i;
		else if (name == "qual1") qual1_idx = i;
		else if (name == "qual2") qual2_idx = i;
	}
}

//===--------------------------------------------------------------------===//
// CommonCopyParameters Implementation
//===--------------------------------------------------------------------===//
void CommonCopyParameters::ParseFromOptions(const case_insensitive_map_t<vector<Value>> &options,
                                            const string &file_path) {
	Value interleave_param;
	Value id_as_sequence_index_param;
	Value include_comment_param;
	Value compression_param;
	
	for (auto &option : options) {
		if (StringUtil::CIEquals(option.first, "interleave")) {
			interleave_param = option.second[0];
		} else if (StringUtil::CIEquals(option.first, "id_as_sequence_index")) {
			id_as_sequence_index_param = option.second[0];
		} else if (StringUtil::CIEquals(option.first, "include_comment")) {
			include_comment_param = option.second[0];
		} else if (StringUtil::CIEquals(option.first, "compression")) {
			compression_param = option.second[0];
		}
	}
	
	if (!interleave_param.IsNull()) {
		interleave = interleave_param.GetValue<bool>();
	}
	
	if (!id_as_sequence_index_param.IsNull()) {
		id_as_sequence_index = id_as_sequence_index_param.GetValue<bool>();
	}
	
	if (!include_comment_param.IsNull()) {
		include_comment = include_comment_param.GetValue<bool>();
	}
	
	compression = DetectCompressionType(file_path, compression_param);
}

//===--------------------------------------------------------------------===//
// Common Validation Functions
//===--------------------------------------------------------------------===//
void ValidateRequiredColumns(bool has_read_id, bool has_sequence1, const string &format_name) {
	if (!has_read_id) {
		throw BinderException("COPY FORMAT %s requires 'read_id' column", format_name);
	}
	if (!has_sequence1) {
		throw BinderException("COPY FORMAT %s requires 'sequence1' column", format_name);
	}
}

void ValidatePairedEndParameters(bool is_paired, bool has_interleave_param, bool interleave,
                                 const string &file_path) {
	// INTERLEAVE parameter required for paired-end data
	if (is_paired && !has_interleave_param) {
		throw BinderException("INTERLEAVE parameter required for paired-end data");
	}
	
	// Validate {ORIENTATION} usage
	bool has_orientation = HasOrientationPlaceholder(file_path);
	if (is_paired && !interleave && !has_orientation) {
		throw BinderException("Paired-end data with INTERLEAVE=false requires {ORIENTATION} in file path");
	}
	if (!is_paired && has_orientation) {
		throw BinderException("Single-end data cannot use {ORIENTATION} in file path");
	}
}

void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index) {
	if (id_as_sequence_index && !has_sequence_index) {
		throw BinderException("ID_AS_SEQUENCE_INDEX=true requires 'sequence_index' column");
	}
}

} // namespace duckdb
