#pragma once
#include <zlib.h>
#include <memory>
#include <string>
#include <kseq++/kseq++.hpp>

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/file_system.hpp"
#include "stream_md5.hpp"
#endif

namespace miint {

// Per-thread error channel for the kseq++ read callbacks (`duckdb_seq_read`,
// `aspera_seq_read`). We deliberately avoid throwing from inside those
// callbacks: they're invoked deep inside kseq++'s template machinery, and
// throwing across that boundary risked a recursive-terminate abort observed in
// the field. Instead, on error the callback writes a human-readable message
// here and returns -1 (which kseq++ tolerates cleanly via its err() flag), and
// `SequenceReader::read_stream` checks this string after kseq++ returns and
// raises `duckdb::IOException` from a frame we control.
//
// Lifetime: thread_local, so concurrent reads on different threads don't
// stomp each other. `read_stream` clears it before each kseq call and consumes
// it afterward.
extern thread_local std::string g_seq_read_error;

// Shared zlib inflate helper used by both DuckDBSeqStream and AsperaSeqStream.
// ReadRawFn signature: int(void *dst, size_t len) — returns bytes read (>0), 0=EOF, <0=error.
// Returns: decompressed bytes written to dst, or -1 on error (decompression error, upstream
// read error, or — critically — EOF reached before Z_STREAM_END, which means the gzip
// stream was truncated).
//
// `stream_end` tracks whether the FINAL gzip member has ended. gzip is a concatenation of
// members (BGZF/block-gzip emits one member per ~64 KiB block, and tools like `cat a.gz b.gz`
// produce multi-member files), so a single Z_STREAM_END is a member boundary, not necessarily
// end-of-file: this helper resets and continues into the next member, matching zlib's gzread.
// Only once no further gzip member follows is `stream_end` set; subsequent calls then return 0
// cleanly rather than treating the trailing "no input, no output" state as truncation. Caller
// initializes to false.
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
			// Member boundary, not necessarily EOF. gzip is concatenated members and BGZF emits
			// one per block, so continue into the next member instead of stopping at the first
			// boundary (the bug that silently truncated bgzf .gz to its first ~64 KiB block when
			// read through this path). A new member starts with the gzip magic 0x1f 0x8b; anything
			// else (or no more input) is the true end -- and trailing non-gzip padding is tolerated,
			// not treated as corruption, matching zlib's gzread.
			// The 2 magic bytes can straddle a read (and a read may return fewer bytes than asked),
			// so keep compacting the leftover to the front and topping up until we have both bytes
			// or hit EOF.
			while (compressed_avail < 2 && !input_eof) {
				for (int k = 0; k < compressed_avail; k++) {
					compressed_buf[k] = compressed_next[k];
				}
				compressed_next = compressed_buf;
				auto n = read_raw(compressed_buf + compressed_avail, buf_size - static_cast<size_t>(compressed_avail));
				if (n < 0) {
					return -1;
				}
				if (n == 0) {
					input_eof = true;
				} else {
					compressed_avail += static_cast<int>(n);
				}
			}
			bool next_is_gzip_member = compressed_avail >= 2 &&
			                           static_cast<unsigned char>(compressed_next[0]) == 0x1f &&
			                           static_cast<unsigned char>(compressed_next[1]) == 0x8b;
			if (next_is_gzip_member) {
				if (inflateReset(&zs) != Z_OK) {
					return -1;
				}
				continue; // decode the next member into the same output buffer
			}
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

	// Optional md5 verification tap over the raw (pre-decompression) bytes as
	// they arrive from `handle`. Null when verification isn't requested for
	// this stream (the common case — read_fastx callers never set it).
	// Shared (not owned) so the constructing caller (PerRunReader) can retain
	// its own reference and call VerifyOrThrow once the read loop reaches
	// true EOF; a fresh DuckDBSeqStream (e.g. the empty-R2 reopen) gets a
	// fresh StreamMd5 passed in, so there's nothing to "reset" here.
	std::shared_ptr<miint::StreamMd5> md5_tap;

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
// `md5_tap`, when non-null, receives the raw (pre-decompression) bytes read from `handle` — see
// duckdb_seq_read. Optional; read_fastx and other callers that don't need verification omit it.
miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path,
                                              std::shared_ptr<miint::StreamMd5> md5_tap = nullptr);

// Overload that accepts an explicit gzip flag instead of inferring from path extension.
// Use when the file path does not reflect the actual compression (e.g., temp files).
miint::DuckDBSeqStream *CreateDuckDBSeqStream(FileSystem &fs, const std::string &path, bool is_gzipped,
                                              std::shared_ptr<miint::StreamMd5> md5_tap = nullptr);
} // namespace duckdb
#endif
