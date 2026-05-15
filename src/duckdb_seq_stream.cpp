#include "duckdb_seq_stream.hpp"

#ifdef MIINT_STATIC_BUILD

namespace miint {

thread_local std::string g_seq_read_error;

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
	// On any error: write a message to the thread-local error channel and
	// return -1. SequenceReader::read_stream picks it up after kseq++ returns
	// and raises IOException from a frame we control. See g_seq_read_error
	// docs in the header for the rationale.
	try {
		if (!stream->is_gzipped) {
			auto n = stream->handle->Read(dst, len);
			if (n < 0) {
				g_seq_read_error = "duckdb_seq_read: upstream read failed";
				return -1;
			}
			return static_cast<int>(n);
		}

		auto read_raw = [stream](void *buf, size_t sz) -> int {
			auto n = stream->handle->Read(buf, sz);
			return static_cast<int>(n);
		};

		int result = InflateFromSource(stream->zs, stream->compressed_buf, DuckDBSeqStream::COMPRESSED_BUF_SIZE,
		                               stream->compressed_avail, stream->compressed_next, stream->input_eof,
		                               stream->stream_end, read_raw, dst, len);
		if (result < 0) {
			g_seq_read_error = stream->stream_end ? "duckdb_seq_read: gzip decompression error"
			                                      : "duckdb_seq_read: gzip stream truncated before end marker "
			                                        "(download likely incomplete)";
			return -1;
		}
		return result;
	} catch (const std::exception &e) {
		// DuckDB's FileHandle::Read can throw on HTTP errors (libcurl
		// CURLE_RECV_ERROR, etc.). Catch here so the throw never crosses
		// kseq++'s template machinery.
		g_seq_read_error = std::string("duckdb_seq_read: ") + e.what();
		return -1;
	} catch (...) {
		g_seq_read_error = "duckdb_seq_read: unknown error";
		return -1;
	}
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
