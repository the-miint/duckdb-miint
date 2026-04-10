#pragma once
#include <zlib.h>
#include <kseq++/kseq++.hpp>

#ifdef MIINT_STATIC_BUILD
#include "duckdb/common/file_system.hpp"
#endif

namespace miint {

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
} // namespace duckdb
#endif
