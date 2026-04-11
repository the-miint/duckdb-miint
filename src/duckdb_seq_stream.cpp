#include "duckdb_seq_stream.hpp"

#ifdef MIINT_STATIC_BUILD

namespace miint {

DuckDBSeqStream::DuckDBSeqStream()
    : is_gzipped(false), zs({}), zs_initialized(false), compressed_avail(0), compressed_next(nullptr),
      input_eof(false) {
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
			return -1;
		}
		return static_cast<int>(n);
	}

	auto read_raw = [stream](void *buf, size_t sz) -> int {
		auto n = stream->handle->Read(buf, sz);
		return (n <= 0) ? 0 : static_cast<int>(n);
	};

	return InflateFromSource(stream->zs, stream->compressed_buf, DuckDBSeqStream::COMPRESSED_BUF_SIZE,
	                         stream->compressed_avail, stream->compressed_next, stream->input_eof, read_raw, dst, len);
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
	auto stream = std::make_unique<miint::DuckDBSeqStream>();
	auto handle = fs.OpenFile(path, FileOpenFlags(FileOpenFlags::FILE_FLAGS_READ));
	stream->handle = std::shared_ptr<FileHandle>(handle.release());
	stream->is_gzipped = IsGzipped(path);
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
