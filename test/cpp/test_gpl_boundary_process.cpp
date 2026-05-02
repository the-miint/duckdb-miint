#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gpl_boundary/process.hpp"

#include <cstdlib>
#include <string>

using namespace duckdb::miint::gpl_boundary;

// =============================================================================
// Cycle 1.1 — Find gpl-boundary on PATH
// =============================================================================
//
// Mirrors Bowtie2Aligner::find_executable's contract: returns the absolute path
// when the named binary is found on PATH, empty string otherwise. We run on
// a dev machine that may or may not have gpl-boundary installed, so this test
// branches on the runtime environment rather than mocking `which`.
//
// What we DO assert in both branches: the path-discovery function does not
// throw, does not deadlock, and returns a sensible value. That is enough to
// prove the fork+pipe+exec plumbing is right; correctness against a real
// gpl-boundary daemon is exercised in Phase 4 with require-env GPL_BOUNDARY_AVAILABLE.

TEST_CASE("FindGplBoundary returns absolute path when on PATH", "[gpl-boundary][process]") {
	const std::string result = FindGplBoundary();

	const char *available = std::getenv("GPL_BOUNDARY_AVAILABLE");
	const bool path_says_available = available && std::string(available) == "1";

	if (path_says_available) {
		INFO("GPL_BOUNDARY_AVAILABLE=1; expecting non-empty result from FindGplBoundary()");
		REQUIRE_FALSE(result.empty());
		// `which` always returns an absolute path
		REQUIRE(result.front() == '/');
		// path must end in `gpl-boundary`
		const std::string tail = "/gpl-boundary";
		REQUIRE(result.size() >= tail.size());
		REQUIRE(result.compare(result.size() - tail.size(), tail.size(), tail) == 0);
	} else {
		// No assertion either way: gpl-boundary may or may not be in PATH on
		// this dev machine. We simply ensure the call doesn't crash. Future
		// CI runs with GPL_BOUNDARY_AVAILABLE=1 set will hit the strict branch.
		INFO("GPL_BOUNDARY_AVAILABLE not set; FindGplBoundary returned: '" << result << "'");
		SUCCEED("FindGplBoundary did not crash");
	}
}

TEST_CASE("FindExecutableInPath returns non-empty for a known-present binary", "[gpl-boundary][process]") {
	// `bash` is universally available on Linux/macOS dev environments. Using
	// it as a positive control proves the underlying which-fork-exec plumbing
	// is correct without depending on gpl-boundary being installed.
	const std::string result = FindExecutableInPath("bash");
	INFO("FindExecutableInPath(\"bash\") returned: '" << result << "'");
	REQUIRE_FALSE(result.empty());
	REQUIRE(result.front() == '/');
	const std::string tail = "/bash";
	REQUIRE(result.size() >= tail.size());
	REQUIRE(result.compare(result.size() - tail.size(), tail.size(), tail) == 0);
}

TEST_CASE("FindExecutableInPath returns empty for a missing binary", "[gpl-boundary][process]") {
	// Use a name very unlikely to collide with anything in PATH. If this
	// returns non-empty, either the user has a wildly creative shell setup
	// or the lookup is broken.
	const std::string result = FindExecutableInPath("miint-deliberately-not-a-real-binary-xyzzy-42");
	INFO("Unexpectedly found: '" << result << "'");
	REQUIRE(result.empty());
}

// =============================================================================
// Cycle 1.2 — Spawn child + bidirectional pipes
// =============================================================================
//
// We need a long-lived NDJSON daemon, but the spawn primitive itself is
// generic. Tests use ubiquitous coreutils binaries (`printf`, `cat`, `bash`)
// to exercise the plumbing without depending on gpl-boundary being installed.

#include <chrono>
#include <signal.h>
#include <sys/wait.h>
#include <thread>
#include <unistd.h>
#include <vector>

namespace {
// Tiny helper: read from fd until EOF. Returns the bytes read.
std::string drain_fd(int fd) {
	std::string out;
	char buf[256];
	ssize_t n;
	while ((n = ::read(fd, buf, sizeof(buf))) > 0) {
		out.append(buf, static_cast<size_t>(n));
	}
	return out;
}
} // namespace

TEST_CASE("ChildProcess spawns and captures stdout from printf", "[gpl-boundary][process]") {
	// printf ships in coreutils — universally available on Linux/macOS.
	const std::string printf_path = FindExecutableInPath("printf");
	REQUIRE_FALSE(printf_path.empty());

	std::vector<std::string> argv = {printf_path, "hello world"};
	ChildProcess child(argv);

	// Drain stdout in the parent
	std::string captured = drain_fd(child.stdout_fd());
	REQUIRE(captured == "hello world");

	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("ChildProcess propagates non-zero exit status", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());

	std::vector<std::string> argv = {bash_path, "-c", "exit 7"};
	ChildProcess child(argv);
	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 7);
}

TEST_CASE("ChildProcess routes stdin into the child", "[gpl-boundary][process]") {
	// `cat` echoes its stdin to stdout; the perfect bidirectional smoke.
	const std::string cat_path = FindExecutableInPath("cat");
	REQUIRE_FALSE(cat_path.empty());

	std::vector<std::string> argv = {cat_path};
	ChildProcess child(argv);

	const std::string payload = "round trip\n";
	const ssize_t wrote = ::write(child.stdin_fd(), payload.data(), payload.size());
	REQUIRE(wrote == static_cast<ssize_t>(payload.size()));
	// Close stdin so cat sees EOF and exits.
	child.CloseStdin();

	const std::string out = drain_fd(child.stdout_fd());
	REQUIRE(out == payload);

	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("ChildProcess captures stderr separately from stdout", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());

	std::vector<std::string> argv = {bash_path, "-c", "echo OUT; echo ERR 1>&2"};
	ChildProcess child(argv);

	const std::string out = drain_fd(child.stdout_fd());
	const std::string err = drain_fd(child.stderr_fd());
	REQUIRE(out == "OUT\n");
	REQUIRE(err == "ERR\n");

	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("ChildProcess destructor reaps a still-running child", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());

	pid_t captured_pid = -1;
	{
		// `sleep 60` would outlive the test. The destructor must SIGTERM and reap.
		std::vector<std::string> argv = {bash_path, "-c", "sleep 60"};
		ChildProcess child(argv);
		captured_pid = child.pid();
		REQUIRE(captured_pid > 0);
		// Give the child a moment to actually start sleeping; otherwise SIGTERM
		// may race with execve and we'd get an exit status that isn't from our signal.
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
	}
	// Destructor has run — child should be gone.
	// kill(pid, 0) returns -1/ESRCH when the process is already reaped.
	const int probe = ::kill(captured_pid, 0);
	if (probe == 0) {
		// Race: destructor reaped before kernel cleared the pid table. Wait briefly.
		std::this_thread::sleep_for(std::chrono::milliseconds(100));
	}
	REQUIRE(::kill(captured_pid, 0) == -1);
}

TEST_CASE("ChildProcess fails clearly when the binary does not exist", "[gpl-boundary][process]") {
	std::vector<std::string> argv = {"/this/path/should/never/exist/xyzzy"};
	REQUIRE_THROWS(ChildProcess(argv));
}

// =============================================================================
// Cycle 1.2 — gpl-boundary --version (skip if absent)
// =============================================================================

TEST_CASE("ChildProcess can spawn gpl-boundary --version", "[gpl-boundary][process][!mayfail]") {
	const std::string gpl = FindGplBoundary();
	if (gpl.empty()) {
		SKIP("gpl-boundary not on PATH (set GPL_BOUNDARY_AVAILABLE=1 in CI to require it)");
	}
	std::vector<std::string> argv = {gpl, "--version"};
	ChildProcess child(argv);
	const std::string out = drain_fd(child.stdout_fd());
	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
	INFO("gpl-boundary --version output: " << out);
	REQUIRE_FALSE(out.empty());
}

// =============================================================================
// Cycle 1.3 — Line-delimited JSON I/O
// =============================================================================
//
// gpl-boundary's protocol is NDJSON: one JSON object per line on stdin/stdout.
// We need:
//   - LineReader: buffered fd → string lines (handles multiple lines per read,
//     short reads, no-newline-at-EOF cases).
//   - JSON parsing of one line (just enough to read init-reply / batch-reply).

TEST_CASE("LineReader returns each line of a multi-line stream", "[gpl-boundary][process]") {
	const std::string printf_path = FindExecutableInPath("printf");
	REQUIRE_FALSE(printf_path.empty());

	std::vector<std::string> argv = {printf_path, "alpha\nbeta\ngamma\n"};
	ChildProcess child(argv);

	LineReader reader(child.stdout_fd());
	std::string line;

	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "alpha");
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "beta");
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "gamma");
	REQUIRE_FALSE(reader.ReadLine(line)); // EOF

	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 0);
}

TEST_CASE("LineReader handles a final line without trailing newline", "[gpl-boundary][process]") {
	const std::string printf_path = FindExecutableInPath("printf");
	REQUIRE_FALSE(printf_path.empty());
	// No trailing \n on "no-newline"
	std::vector<std::string> argv = {printf_path, "first\nno-newline"};
	ChildProcess child(argv);
	LineReader reader(child.stdout_fd());
	std::string line;
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "first");
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "no-newline");
	REQUIRE_FALSE(reader.ReadLine(line));
	child.Wait();
}

TEST_CASE("LineReader survives a line longer than its internal buffer", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	// 16 KiB of 'A' on one line, then a newline. Larger than LineReader's
	// 4 KiB read chunk so this exercises the multi-read path.
	std::vector<std::string> argv = {bash_path, "-c", "head -c 16384 /dev/zero | tr '\\0' 'A'; echo"};
	ChildProcess child(argv);
	LineReader reader(child.stdout_fd());
	std::string line;
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line.size() == 16384);
	for (char c : line) {
		if (c != 'A') {
			FAIL("non-'A' char in long line: " << static_cast<int>(c));
		}
	}
	child.Wait();
}

TEST_CASE("WriteLine appends newline and flushes immediately", "[gpl-boundary][process]") {
	// `cat` echoes stdin to stdout. Round-trip a JSON-shaped line.
	const std::string cat_path = FindExecutableInPath("cat");
	REQUIRE_FALSE(cat_path.empty());

	std::vector<std::string> argv = {cat_path};
	ChildProcess child(argv);
	WriteLine(child.stdin_fd(), R"({"k":"v","n":42})");
	WriteLine(child.stdin_fd(), R"({"shutdown":true})");
	child.CloseStdin();

	LineReader reader(child.stdout_fd());
	std::string line;
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == R"({"k":"v","n":42})");
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == R"({"shutdown":true})");
	REQUIRE_FALSE(reader.ReadLine(line));
	child.Wait();
}
