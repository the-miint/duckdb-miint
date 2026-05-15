#pragma once
#include <zlib.h>
#include <memory>
#include <kseq++/kseq++.hpp>

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/file_system.hpp"
#endif

namespace miint {

// Shared zlib inflate helper used by both DuckDBSeqStream and AsperaSeqStream.
// ReadRawFn signature: int(void *dst, size_t len) — returns bytes read (>0), 0=EOF, <0=error.
// Returns: decompressed bytes written to dst, or -1 on error (decompression error, upstream
// read error, or — critically — EOF reached before Z_STREAM_END, which means the gzip
// stream was truncated).
//
// `stream_end` tracks whether Z_STREAM_END has been observed across calls. Once set, the
// function returns 0 cleanly on subsequent invocations rather than treating the trailing
// "no input, no output" state as truncation. Caller initializes to false.
template <typename ReadRawFn>
int InflateFromSource(z_stream &zs, char *compressed_buf, size_t buf_size, int &compressed_avail,
                      char *&compressed_next, bool &input_eof, bool &stream_end, ReadRawFn read_raw, void *dst,
                      unsigned int len) {
	if (stream_end) {
		return 0;
	}

	zs.avail_out = len;
	zs.next_out = reinterpret_cast<Bytef *>(dst);

	while (zs.avail_out > 0) {
		if (compressed_avail == 0 && !input_eof) {
			auto n = read_raw(compressed_buf, buf_size);
			if (n < 0) {
				// Upstream read failure — propagate. Do NOT set input_eof; a
				// retried call would otherwise look like a clean EOF and
				// silently truncate (the bug this whole helper now guards
				// against).
				return -1;
			}
			if (n == 0) {
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
			stream_end = true;
			break;
		}
		if (ret != Z_OK) {
			return -1;
		}
		if (input_eof && compressed_avail == 0) {
			// EOF reached before the gzip trailer (CRC32 + ISIZE) signaled
			// Z_STREAM_END. Stream was truncated mid-block. Return -1 so the
			// caller throws — silent acceptance here is exactly how
			// half-downloaded ENA runs used to slip through.
			return -1;
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
	bool stream_end; // Z_STREAM_END observed → subsequent reads are legitimately EOF

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

// Owning handle for a raw DuckDBSeqStream until it is transferred into a kseq++
// KStreamIn (i.e., into a SequenceReader). The deleter is duckdb_seq_close,
// which does `delete stream` — the same operation kseq++'s close callback
// performs, so handing the handle off to the SequenceReader by value (which
// .release()s once the kstream wrapper has taken ownership) prevents the
// double-free that arose from raw-pointer + try/catch ownership dancing.
using DuckDBSeqStreamHandle = std::unique_ptr<DuckDBSeqStream, int (*)(DuckDBSeqStream *)>;

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
