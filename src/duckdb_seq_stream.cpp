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
		// Uncompressed: delegate directly to FileHandle
		auto n = stream->handle->Read(dst, len);
		if (n < 0) {
			return -1;
		}
		return static_cast<int>(n);
	}

	// Gzipped: decompress on the fly
	stream->zs.avail_out = len;
	stream->zs.next_out = reinterpret_cast<Bytef *>(dst);

	while (stream->zs.avail_out > 0) {
		// Refill compressed buffer if needed
		if (stream->compressed_avail == 0 && !stream->input_eof) {
			auto n = stream->handle->Read(stream->compressed_buf, DuckDBSeqStream::COMPRESSED_BUF_SIZE);
			if (n <= 0) {
				stream->input_eof = true;
			} else {
				stream->compressed_avail = static_cast<int>(n);
				stream->compressed_next = stream->compressed_buf;
			}
		}

		stream->zs.avail_in = static_cast<uInt>(stream->compressed_avail);
		stream->zs.next_in = reinterpret_cast<Bytef *>(stream->compressed_next);

		int ret = inflate(&stream->zs, Z_NO_FLUSH);

		// Update consumed compressed bytes
		int consumed = stream->compressed_avail - static_cast<int>(stream->zs.avail_in);
		stream->compressed_avail -= consumed;
		stream->compressed_next += consumed;

		if (ret == Z_STREAM_END) {
			break;
		}
		if (ret != Z_OK) {
			return -1;
		}
		if (stream->input_eof && stream->compressed_avail == 0) {
			break;
		}
	}

	return static_cast<int>(len - stream->zs.avail_out);
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
