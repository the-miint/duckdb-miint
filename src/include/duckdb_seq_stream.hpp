#pragma once
#include <zlib.h>
#include <kseq++/kseq++.hpp>

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/file_system.hpp"
#endif

namespace miint {

// Shared zlib inflate helper used by both DuckDBSeqStream and AsperaSeqStream.
// ReadRawFn signature: int(void *dst, size_t len) — returns bytes read (>0), 0=EOF, <0=error.
// Returns: decompressed bytes written to dst.
template <typename ReadRawFn>
int InflateFromSource(z_stream &zs, char *compressed_buf, size_t buf_size, int &compressed_avail,
                      char *&compressed_next, bool &input_eof, ReadRawFn read_raw, void *dst, unsigned int len) {
	zs.avail_out = len;
	zs.next_out = reinterpret_cast<Bytef *>(dst);

	while (zs.avail_out > 0) {
		if (compressed_avail == 0 && !input_eof) {
			auto n = read_raw(compressed_buf, buf_size);
			if (n <= 0) {
				input_eof = true;
			} else {
				compressed_avail = static_cast<int>(n);
				compressed_next = compressed_buf;
			}
		}

		zs.avail_in = static_cast<uInt>(compressed_avail);
		zs.next_in = reinterpret_cast<Bytef *>(compressed_next);

		int ret = inflate(&zs, Z_NO_FLUSH);

		int consumed = compressed_avail - static_cast<int>(zs.avail_in);
		compressed_avail -= consumed;
		compressed_next += consumed;

		if (ret == Z_STREAM_END) {
			break;
		}
		if (ret != Z_OK) {
			return -1;
		}
		if (input_eof && compressed_avail == 0) {
			break;
		}
	}

	return static_cast<int>(len - zs.avail_out);
}

#ifdef MIINT_STATIC_BUILD

// Stream adapter that reads from DuckDB FileHandle with optional gzip decompression.
// Designed to be used as kseq++'s TFile parameter for streaming remote files.
// Ownership: allocated with `new`, freed by duckdb_seq_close (called from kseq++ destructor).
struct DuckDBSeqStream {
	std::shared_ptr<duckdb::FileHandle> handle;
	bool is_gzipped;

	// Gzip decompression state
	z_stream zs;
	bool zs_initialized;
	static constexpr size_t COMPRESSED_BUF_SIZE = 65536;
	char compressed_buf[COMPRESSED_BUF_SIZE];
	int compressed_avail;
	char *compressed_next;
	bool input_eof;

	DuckDBSeqStream();
	~DuckDBSeqStream();

	DuckDBSeqStream(const DuckDBSeqStream &) = delete;
	DuckDBSeqStream &operator=(const DuckDBSeqStream &) = delete;
};

// kseq++ read callback: >0 bytes read, 0 = EOF, <0 = error
int duckdb_seq_read(DuckDBSeqStream *stream, void *dst, unsigned int len);

// kseq++ close callback: cleans up zlib state and deletes the stream
int duckdb_seq_close(DuckDBSeqStream *stream);

// kseq++ KStreamIn instantiation for DuckDB streams
using DuckDBSeqStreamIn = klibpp::KStreamIn<DuckDBSeqStream *, int (*)(DuckDBSeqStream *, void *, unsigned int)>;

#endif // MIINT_STATIC_BUILD

} // namespace miint

#ifdef MIINT_STATIC_BUILD
namespace duckdb {
// Create a DuckDBSeqStream for reading a remote (or local) file through DuckDB's FileSystem.
// Handles gzip detection and decompression initialization.
// Caller takes ownership of the returned pointer (kseq++ close callback deletes it).
miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path);

// Overload that accepts an explicit gzip flag instead of inferring from path extension.
// Use when the file path does not reflect the actual compression (e.g., temp files).
miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path, bool is_gzipped);
} // namespace duckdb
#endif
