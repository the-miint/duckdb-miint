// Unit tests for InflateFromSource (src/include/duckdb_seq_stream.hpp), the streaming gzip
// decompressor used by the REMOTE read path (read_fastx over http/s3, read_ena_sequences).
//
// Regression focus: gzip is a concatenation of members and BGZF (block-gzip, what the FASTQ/FASTA
// COPY writer now produces) emits one member per ~64 KiB block. The decoder must continue across
// member boundaries; an earlier version stopped at the first Z_STREAM_END and silently truncated a
// remotely-read bgzf file to its first block. The local read path uses zlib gzread (multi-member
// capable), so only this path was affected and only over the network -- exactly the gap SQL
// round-trip tests on local files could not see.

#include <catch2/catch_test_macros.hpp>

#include <algorithm>
#include <cstring>
#include <string>
#include <vector>
#include <zlib.h>

#include <duckdb_seq_stream.hpp>

namespace {

// Compress `data` into one standalone gzip member (RFC1952) -- the unit BGZF and `cat a.gz b.gz`
// concatenate.
std::string GzipMember(const std::string &data) {
	z_stream zs;
	std::memset(&zs, 0, sizeof(zs));
	REQUIRE(deflateInit2(&zs, 6, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) == Z_OK);
	std::string out;
	out.resize(deflateBound(&zs, data.size()) + 64);
	zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(data.data()));
	zs.avail_in = static_cast<uInt>(data.size());
	zs.next_out = reinterpret_cast<Bytef *>(&out[0]);
	zs.avail_out = static_cast<uInt>(out.size());
	REQUIRE(deflate(&zs, Z_FINISH) == Z_STREAM_END);
	out.resize(out.size() - zs.avail_out);
	deflateEnd(&zs);
	return out;
}

// Drive InflateFromSource over `compressed`, serving it in `chunk`-byte reads so member boundaries
// straddle the read buffer (stressing the magic-byte lookahead/compaction). Returns the fully
// decompressed bytes; sets `ok` false on any reported error.
std::string Decode(const std::string &compressed, size_t chunk, bool &ok) {
	size_t off = 0;
	auto read_raw = [&](void *buf, size_t sz) -> int {
		size_t n = std::min({sz, chunk, compressed.size() - off});
		std::memcpy(buf, compressed.data() + off, n);
		off += n;
		return static_cast<int>(n);
	};

	z_stream zs;
	std::memset(&zs, 0, sizeof(zs));
	if (inflateInit2(&zs, 16 + MAX_WBITS) != Z_OK) {
		ok = false;
		return {};
	}
	std::vector<char> cbuf(256);
	int cavail = 0;
	char *cnext = cbuf.data();
	bool input_eof = false, stream_end = false;
	std::string decoded;
	char dst[1024];
	ok = true;
	for (;;) {
		int n = miint::InflateFromSource(zs, cbuf.data(), cbuf.size(), cavail, cnext, input_eof, stream_end, read_raw,
		                                 dst, sizeof(dst));
		if (n < 0) {
			ok = false;
			break;
		}
		if (n == 0) {
			break;
		}
		decoded.append(dst, static_cast<size_t>(n));
	}
	inflateEnd(&zs);
	return decoded;
}

// Non-trivially-compressible content (varies per byte) so the gzip members are real, multi-block.
std::string VariedBlock(char base, size_t len) {
	std::string s(len, base);
	for (size_t i = 0; i < len; i++) {
		s[i] = static_cast<char>(base + static_cast<char>(i % 23));
	}
	return s;
}

} // namespace

TEST_CASE("InflateFromSource decodes a single-member gzip stream", "[seqstream]") {
	std::string a = VariedBlock('A', 40000);
	bool ok = false;
	std::string decoded = Decode(GzipMember(a), 4096, ok);
	REQUIRE(ok);
	CHECK(decoded == a);
}

TEST_CASE("InflateFromSource decodes multi-member (bgzf-style concatenated) gzip", "[seqstream]") {
	// Two members back-to-back, exactly the shape BGZF produces (one member per block).
	std::string a = VariedBlock('A', 20000);
	std::string b = VariedBlock('a', 20000);
	std::string stream = GzipMember(a) + GzipMember(b);

	// chunk=1 guarantees the 2-byte gzip magic of member 2 arrives across two separate reads,
	// exercising the boundary compaction; the others cover the common cases.
	for (size_t chunk : {static_cast<size_t>(1), static_cast<size_t>(7), static_cast<size_t>(4096)}) {
		bool ok = false;
		std::string decoded = Decode(stream, chunk, ok);
		REQUIRE(ok);
		// Pre-fix this returned just `a` (first Z_STREAM_END treated as end-of-stream).
		CHECK(decoded == a + b);
	}
}

TEST_CASE("InflateFromSource consumes a trailing empty member (BGZF EOF marker)", "[seqstream]") {
	// BGZF terminates with a 28-byte empty block: a valid gzip member that decompresses to nothing.
	std::string a = VariedBlock('A', 30000);
	std::string stream = GzipMember(a) + GzipMember("");
	bool ok = false;
	std::string decoded = Decode(stream, 8, ok);
	REQUIRE(ok);
	CHECK(decoded == a);
}

TEST_CASE("InflateFromSource still reports truncation rather than silent short read", "[seqstream]") {
	std::string a = VariedBlock('A', 30000);
	std::string member = GzipMember(a);
	std::string truncated = member.substr(0, member.size() - 6); // drop the gzip trailer
	bool ok = true;
	Decode(truncated, 4096, ok);
	CHECK_FALSE(ok); // EOF before Z_STREAM_END must surface as an error
}
