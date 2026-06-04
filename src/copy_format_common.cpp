#include "copy_format_common.hpp"
#include "remote_file_helper.hpp"
#include "bgzf_write_mode.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/numeric_utils.hpp"
#include "duckdb/parallel/task_scheduler.hpp"

#include <htslib-1.22.1/htslib/bgzf.h>
#include <cerrno>

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
// Write the accumulated local buffer to the file in <= 1 GiB slices, then reset the buffer.
// Caller must hold the global lock. This guards every CopyFileHandle::Write backend, not just one:
// the remote-gzip fallback goes through DuckDB's MiniZStreamWrapper, which casts the write size to
// unsigned int and throws "Information loss on integer cast" on a single >4 GiB write; chunking
// keeps that (and any 32-bit limit in the bgzf/BufferedFileWriter backends) safe regardless of how
// much the local MemoryStream accumulated.
static void WriteBufferToFile(CopyFileHandle &file, FormatWriterState &local_state) {
	constexpr idx_t MAX_WRITE_CHUNK = 1ULL * 1024 * 1024 * 1024; // 1 GiB
	const_data_ptr_t data = local_state.stream->GetData();
	idx_t remaining = local_state.stream->GetPosition();
	idx_t offset = 0;
	while (remaining > 0) {
		idx_t to_write = remaining < MAX_WRITE_CHUNK ? remaining : MAX_WRITE_CHUNK;
		file.Write(data + offset, to_write);
		offset += to_write;
		remaining -= to_write;
	}

	// Reset local buffer
	local_state.Reset();
}

void FlushFormatBuffer(FormatWriterState &local_state, CopyFileHandle &file, mutex &lock) {
	if (!local_state.written_anything) {
		return;
	}

	lock_guard<mutex> glock(lock);
	WriteBufferToFile(file, local_state);
}

void FlushR2Buffer(FormatWriterState &local_state, SequenceCopyGlobalState &gstate) {
	if (!local_state.written_anything) {
		return;
	}

	lock_guard<mutex> glock(gstate.lock);
	if (!gstate.file_r2) {
		// First non-NULL R2 record across all threads -> create the R2 file now. This is what
		// keeps an all-single-end dataset from producing an empty R2 file in split mode.
		gstate.file_r2 =
		    make_uniq<CopyFileHandle>(*gstate.fs, gstate.r2_path, gstate.compression, gstate.compression_threads);
	}
	WriteBufferToFile(*gstate.file_r2, local_state);
}

//===--------------------------------------------------------------------===//
// CopyFileHandle Implementation
//===--------------------------------------------------------------------===//
CopyFileHandle::CopyFileHandle(FileSystem &fs, const string &path, FileCompressionType compression_p,
                               int compression_threads) {
	// Local gzip output goes through htslib bgzf: it is gzip-compatible (reads back through any gzip
	// reader, including read_fastx) and bgzf_mt parallelizes deflate while emitting blocks in order,
	// so a single-threaded in-order feed still produces a deterministically-ordered file. DuckDB's
	// own gzip writer is single-threaded, which is the bottleneck this avoids. Remote targets keep
	// the BufferedFileWriter path: htslib's hFILE cannot open DuckDB's virtual filesystem, and we use
	// the same RemoteFileHelper::IsRemotePath test the reader (read_fastx) uses so writer and reader
	// agree on what counts as local.
	if (compression_p == FileCompressionType::GZIP && !miint::RemoteFileHelper::IsRemotePath(path)) {
		// BgzfWriteMode picks "wx6" for regular files (atomic no-clobber/no-TOCTOU, the
		// FILE_CREATE_NEW contract) and "w6" for stdout / pipes / devices, which always exist and
		// would otherwise fail O_EXCL with EEXIST. bgzf_open maps "-" (and /dev/stdout) to stdout.
		errno = 0;
		bgzf_file = bgzf_open(path.c_str(), miint::BgzfWriteMode(path));
		if (!bgzf_file) {
			if (errno == EEXIST) {
				throw IOException("Failed to open \"%s\" for writing: file already exists", path);
			}
			throw IOException("Failed to open \"%s\" for bgzf writing", path);
		}
		if (compression_threads > 1) {
			// 256 sub-blocks per thread batch is htslib's typical default. Best-effort: on failure
			// bgzf transparently falls back to single-threaded compression.
			bgzf_mt(bgzf_file, compression_threads, 256);
		}
		return;
	}

	auto flags =
	    FileFlags::FILE_FLAGS_WRITE | FileFlags::FILE_FLAGS_FILE_CREATE_NEW | FileLockType::WRITE_LOCK | compression_p;

	// BufferedFileWriter handles both file opening and buffering
	file_writer = make_uniq<BufferedFileWriter>(fs, path, flags);
}

CopyFileHandle::~CopyFileHandle() {
	// Never let a flush error escape a destructor; Finalize calls Close() explicitly
	// for the error-reporting path.
	try {
		Close();
	} catch (...) { // NOLINT: destructor must not throw
	}
}

void CopyFileHandle::Write(const_data_ptr_t data, idx_t size) {
	if (bgzf_file) {
		auto written = bgzf_write(bgzf_file, data, size);
		if (written < 0 || static_cast<idx_t>(written) != size) {
			throw IOException("bgzf_write failed (wrote %lld of %llu bytes)", static_cast<long long>(written),
			                  static_cast<unsigned long long>(size));
		}
		return;
	}
	if (file_writer) {
		file_writer->WriteData(data, size);
	}
}

void CopyFileHandle::WriteString(const string &data) {
	Write(const_data_ptr_cast(data.c_str()), data.size());
}

void CopyFileHandle::Close() {
	if (bgzf_file) {
		auto ret = bgzf_close(bgzf_file); // flushes pending blocks + writes BGZF EOF
		bgzf_file = nullptr;
		if (ret < 0) {
			throw IOException("bgzf_close failed");
		}
		return;
	}
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
			// TODO: zstd support temporarily disabled due to unicode encoding issues during read-back
			throw NotImplementedException("zstd compression is temporarily disabled for FASTX/SAM formats. "
			                              "Please use gzip compression instead. "
			                              "zstd support will be re-enabled in a future release.");
		} else if (comp_str == "none") {
			return FileCompressionType::UNCOMPRESSED;
		} else {
			throw InvalidInputException("compression must be 'gzip', 'gz', or 'none' (zstd temporarily disabled)");
		}
	}

	// Auto-detect from extension
	if (file_path.size() >= 3 && file_path.substr(file_path.size() - 3) == ".gz") {
		return FileCompressionType::GZIP;
	}
	if (file_path.size() >= 4 && file_path.substr(file_path.size() - 4) == ".zst") {
		// TODO: zstd support temporarily disabled due to unicode encoding issues during read-back
		throw NotImplementedException(
		    "zstd compression (.zst extension) is temporarily disabled for FASTX/SAM formats. "
		    "Please use .gz extension or COMPRESSION='gzip' instead. "
		    "zstd support will be re-enabled in a future release.");
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
		if (name == "read_id")
			read_id_idx = i;
		else if (name == "sequence_index")
			sequence_index_idx = i;
		else if (name == "comment")
			comment_idx = i;
		else if (name == "sequence1")
			sequence1_idx = i;
		else if (name == "sequence2")
			sequence2_idx = i;
		else if (name == "qual1")
			qual1_idx = i;
		else if (name == "qual2")
			qual2_idx = i;
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

void ResolveOutputMode(const string &file_path, bool interleave, bool &out_split_output) {
	// The {ORIENTATION} placeholder is the split-vs-single switch. INTERLEAVE=true means "one file
	// with R1/R2 alternating", which is the opposite of split output -- so the two are contradictory.
	bool has_orientation = HasOrientationPlaceholder(file_path);
	if (has_orientation && interleave) {
		throw BinderException(
		    "Cannot combine the {ORIENTATION} placeholder with INTERLEAVE=true: use {ORIENTATION} in the path "
		    "to write split R1/R2 files, or INTERLEAVE=true (without {ORIENTATION}) to interleave R1/R2 into a "
		    "single file");
	}
	out_split_output = has_orientation;
}

void ValidateSequenceIndexParameter(bool id_as_sequence_index, bool has_sequence_index) {
	if (id_as_sequence_index && !has_sequence_index) {
		throw BinderException("ID_AS_SEQUENCE_INDEX=true requires 'sequence_index' column");
	}
}

//===--------------------------------------------------------------------===//
// Shared Sequence Copy Functions
//===--------------------------------------------------------------------===//

unique_ptr<GlobalFunctionData> SequenceCopyInitializeGlobal(ClientContext &context, const SequenceCopyBindData &fdata,
                                                            const string &file_path) {
	auto &fs = FileSystem::GetFileSystem(context);

	auto gstate = make_uniq<SequenceCopyGlobalState>();
	gstate->fs = &fs;
	gstate->compression = fdata.compression;
	// Follow DuckDB's configured thread count for bgzf compression workers (only matters for the
	// gzip path; ignored for uncompressed output). Split output opens two bgzf pools (R1 + R2), so
	// halve the budget to avoid running ~2N deflate workers on an N-core host.
	int db_threads = NumericCast<int>(TaskScheduler::GetScheduler(context).NumberOfThreads());
	gstate->compression_threads = (fdata.split_output && db_threads > 1) ? db_threads / 2 : db_threads;

	// R1 is always written. SubstituteOrientation is a no-op when {ORIENTATION} is absent, so this
	// resolves to the raw path for single-file/interleaved output and to the R1 file in split mode.
	string path_r1 = SubstituteOrientation(file_path, "R1");
	gstate->file_r1 = make_uniq<CopyFileHandle>(fs, path_r1, fdata.compression, gstate->compression_threads);

	if (fdata.split_output) {
		// R2 file is created lazily on the first non-NULL R2 record (see FlushR2Buffer), so a
		// single-end dataset written with {ORIENTATION} never produces an empty R2 file.
		gstate->r2_path = SubstituteOrientation(file_path, "R2");
	}

	return gstate;
}

unique_ptr<LocalFunctionData> SequenceCopyInitializeLocal(ExecutionContext &context,
                                                          const SequenceCopyBindData &fdata) {
	auto lstate = make_uniq<SequenceCopyLocalState>();

	lstate->writer_state_r1 = make_uniq<FormatWriterState>(context.client, fdata.flush_size);

	if (fdata.split_output) {
		lstate->writer_state_r2 = make_uniq<FormatWriterState>(context.client, fdata.flush_size);
	}

	return lstate;
}

void SequenceCopyCombine(const SequenceCopyBindData &fdata, SequenceCopyGlobalState &gstate,
                         SequenceCopyLocalState &lstate) {
	// Flush any remaining data in local buffers
	FlushFormatBuffer(*lstate.writer_state_r1, *gstate.file_r1, gstate.lock);

	if (fdata.split_output) {
		FlushR2Buffer(*lstate.writer_state_r2, gstate);
	}

	// Publish this thread's single/paired observations for the cross-thread consistency check.
	lock_guard<mutex> glock(gstate.lock);
	gstate.saw_paired |= lstate.saw_paired;
	gstate.saw_single |= lstate.saw_single;
}

void SequenceCopyFinalize(SequenceCopyGlobalState &gstate) {
	lock_guard<mutex> glock(gstate.lock);

	// A single COPY must be either all single-end or all paired-end. Seeing both means the input
	// mixed paired and unpaired records, which would silently misalign split/interleaved output.
	if (gstate.saw_paired && gstate.saw_single) {
		throw InvalidInputException(
		    "Inconsistent paired-end data: some records have sequence2 set and others do not. A single COPY must "
		    "be either all single-end or all paired-end.");
	}

	if (gstate.file_r1) {
		gstate.file_r1->Close();
	}
	if (gstate.file_r2) {
		gstate.file_r2->Close();
	}
}

} // namespace duckdb
