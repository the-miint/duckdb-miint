// Tests for the local gzip writer's open-mode selection (src/include/bgzf_write_mode.hpp), the fix
// that lets COPY ... (COMPRESSION gzip) stream to stdout / named pipes.
//
// Regression: the writer opened bgzf with mode "wx6" (O_EXCL). O_EXCL rejects any target that
// already exists -- stdout, /dev/stdout, pre-created pipes -- with EEXIST, so
//   duckdb -c "COPY (...) TO '/dev/stdout' (FORMAT FASTQ, COMPRESSION gzip)" | downstream_tool
// failed with "file already exists" even though the uncompressed path streams to them fine. The fix
// keeps O_EXCL only for regular-file targets. This verifies the mode decision and that gzip data
// actually round-trips through a real FIFO opened with the chosen mode.
//
// POSIX-only (mkfifo / O_NONBLOCK): not built on Windows or WASM.

#if !defined(_WIN32) && !defined(__EMSCRIPTEN__)

#include <catch2/catch_test_macros.hpp>

#include <bgzf_write_mode.hpp>
#include <htslib-1.22.1/htslib/bgzf.h>
#include <zlib.h>

#include <cerrno>
#include <cstring>
#include <fcntl.h>
#include <string>
#include <sys/stat.h>
#include <unistd.h>

namespace {

// Inflate the first gzip member (the lone bgzf data block for a small payload) back to plaintext.
std::string GunzipFirstMember(const std::string &in) {
	z_stream zs;
	std::memset(&zs, 0, sizeof(zs));
	REQUIRE(inflateInit2(&zs, 16 + MAX_WBITS) == Z_OK);
	zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(in.data()));
	zs.avail_in = static_cast<uInt>(in.size());
	std::string out;
	char buf[8192];
	int ret;
	do {
		zs.next_out = reinterpret_cast<Bytef *>(buf);
		zs.avail_out = sizeof(buf);
		ret = inflate(&zs, Z_NO_FLUSH);
		REQUIRE((ret == Z_OK || ret == Z_STREAM_END));
		out.append(buf, sizeof(buf) - zs.avail_out);
	} while (ret != Z_STREAM_END);
	inflateEnd(&zs);
	return out;
}

} // namespace

TEST_CASE("BgzfWriteMode keeps O_EXCL only for regular files", "[copy][bgzf][pipe]") {
	// stdout / stderr / fd aliases stream -> no O_EXCL.
	REQUIRE(std::string(miint::BgzfWriteMode("-")) == "w6");
	REQUIRE(std::string(miint::BgzfWriteMode("/dev/stdout")) == "w6");
	REQUIRE(std::string(miint::BgzfWriteMode("/dev/stderr")) == "w6");
	REQUIRE(std::string(miint::BgzfWriteMode("/dev/fd/1")) == "w6");
	// /dev/null is an existing non-regular (char) device -> no O_EXCL.
	REQUIRE(std::string(miint::BgzfWriteMode("/dev/null")) == "w6");
	// A not-yet-existing regular path keeps the atomic no-clobber guard.
	REQUIRE(std::string(miint::BgzfWriteMode("/tmp/miint_nonexistent_xyz.fq.gz")) == "wx6");
}

TEST_CASE("gzip COPY streams to a named pipe without O_EXCL EEXIST", "[copy][bgzf][pipe]") {
	char dir_tmpl[] = "/tmp/miint_pipe_XXXXXX";
	const char *dir = mkdtemp(dir_tmpl);
	REQUIRE(dir != nullptr);
	const std::string fifo = std::string(dir) + "/out.fq.gz";
	REQUIRE(mkfifo(fifo.c_str(), 0600) == 0);

	// The fix: an existing FIFO must select "w6" (no O_EXCL). Pre-fix this was "wx6".
	REQUIRE(std::string(miint::BgzfWriteMode(fifo)) == "w6");

	const std::string payload = "@r1\nACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIII\n";

	// Open the read end non-blocking FIRST so the writer's blocking O_WRONLY open succeeds
	// immediately and a regression (open failing) FAILS the test instead of hanging on a reader with
	// no writer. The payload is far below the pipe buffer, so the writer never blocks on a full pipe.
	const int rfd = open(fifo.c_str(), O_RDONLY | O_NONBLOCK);
	REQUIRE(rfd >= 0);

	// Open + write through the real bgzf path using the production-selected mode. Pre-fix, "wx6" on
	// the existing FIFO returned a null handle (EEXIST); the REQUIRE below would have failed.
	BGZF *bg = bgzf_open(fifo.c_str(), miint::BgzfWriteMode(fifo));
	REQUIRE(bg != nullptr);
	REQUIRE(bgzf_write(bg, payload.data(), payload.size()) == static_cast<ssize_t>(payload.size()));
	REQUIRE(bgzf_close(bg) == 0);

	std::string captured;
	char buf[8192];
	for (;;) {
		ssize_t n = read(rfd, buf, sizeof(buf));
		if (n > 0) {
			captured.append(buf, static_cast<size_t>(n));
			continue;
		}
		if (n == 0) {
			break; // EOF: writer closed
		}
		if (errno == EAGAIN || errno == EWOULDBLOCK) {
			break;
		}
		break;
	}
	close(rfd);

	REQUIRE(!captured.empty());
	REQUIRE(GunzipFirstMember(captured) == payload);

	unlink(fifo.c_str());
	rmdir(dir);
}

#endif // !_WIN32 && !__EMSCRIPTEN__
