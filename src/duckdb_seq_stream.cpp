#include "duckdb_seq_stream.hpp"
#include "duckdb/common/exception.hpp"

#ifdef MIINT_STATIC_BUILD

namespace miint {

DuckDBSeqStream::DuckDBSeqStream()
    : is_gzipped(false), zs({}), zs_initialized(false), compressed_avail(0), compressed_next(nullptr), input_eof(false),
      stream_end(false) {
}

DuckDBSeqStream::~DuckDBSeqStream() {
	if (zs_initialized) {
		inflateEnd(&zs);
	}
}

int duckdb_seq_read(DuckDBSeqStream *stream, void *dst, unsigned int len) {
	if (!stream->is_gzipped) {
		auto n = stream->handle->Read(dst, len);
		if (n < 0) {
			// Upstream HTTP/file read failed. Throwing here surfaces as a
			// mid-stream error in read_ena_sequences (loud warning + run
			// added to skipped_runs), rather than the silent-empty-batch
			// path that previously hid truncated downloads.
			throw duckdb::IOException("duckdb_seq_read: upstream read failed");
		}
		return static_cast<int>(n);
	}

	auto read_raw = [stream](void *buf, size_t sz) -> int {
		// Propagate negative results so InflateFromSource can return -1 and
		// the gzip-truncation detection below fires. The previous mask
		// `(n <= 0) ? 0 : n` actively hid HTTP read errors.
		auto n = stream->handle->Read(buf, sz);
		return static_cast<int>(n);
	};

	int result = InflateFromSource(stream->zs, stream->compressed_buf, DuckDBSeqStream::COMPRESSED_BUF_SIZE,
	                               stream->compressed_avail, stream->compressed_next, stream->input_eof,
	                               stream->stream_end, read_raw, dst, len);
	if (result < 0) {
		// stream_end is set only when Z_STREAM_END was observed. If it's still
		// false here, the gzip stream ended without a valid trailer — the
		// signature of a half-downloaded file. If it's true, decompression
		// itself failed (corrupt data). Either way: throw so the caller treats
		// it as a real failure, not a clean EOF.
		throw duckdb::IOException(stream->stream_end ? "duckdb_seq_read: gzip decompression error"
		                                             : "duckdb_seq_read: gzip stream truncated before end marker "
		                                               "(download likely incomplete)");
	}
	return result;
}

int duckdb_seq_close(DuckDBSeqStream *stream) {
	delete stream;
	return 0;
}

} // namespace miint

#include "table_function_common.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/file_open_flags.hpp"
#include "duckdb/common/file_system.hpp"
#include <memory>

namespace duckdb {

miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path) {
	return CreateDuckDBSeqStream(fs, path, IsGzipped(path));
}

miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path, bool is_gzipped) {
	auto stream = std::make_unique<miint::DuckDBSeqStream>();
	auto handle = fs.OpenFile(path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));
	stream->handle = std::shared_ptr<FileHandle>(handle.release());
	stream->is_gzipped = is_gzipped;
	if (stream->is_gzipped) {
		if (inflateInit2(&stream->zs, 16 + MAX_WBITS) != Z_OK) {
			throw IOException("Failed to initialize zlib for: " + path);
		}
		stream->zs_initialized = true;
	}
	return stream.release();
}

} // namespace duckdb

#endif // MIINT_STATIC_BUILD
