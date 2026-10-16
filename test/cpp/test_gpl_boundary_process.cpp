#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gpl_boundary/process.hpp"

#include <cstdlib>
#include <string>
#include <unistd.h> // mkdtemp lives in <unistd.h> on macOS; <cstdlib>/<stdlib.h> is enough on glibc

using namespace duckdb::miint::gpl_boundary;

// =============================================================================
// Cycle 1.1 — Find gpl-boundary on PATH
// =============================================================================
//
// Contract: returns the absolute path when the named binary is found on PATH,
// empty string otherwise. We run on
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
// Cycle 1.1b — Cache-dir lookup + lookup priority for FindGplBoundary
// =============================================================================
//
// `install_gpl_boundary()` deposits binaries into miint's cache dir
// (`$XDG_CACHE_HOME/miint/bin` or `$HOME/.cache/miint/bin`). FindGplBoundary
// has to find them there even when the user hasn't put the cache dir on PATH,
// otherwise the install function is useless. Priority order:
//   1. MIINT_GPL_BOUNDARY_PATH override (no validation)
//   2. cache dir
//   3. PATH (`which`)

#include <sys/stat.h>

namespace {
// Save+restore an env var across a scope. `unsetenv("X")` if value is nullptr.
struct ScopedEnv {
	std::string name;
	bool had_prior;
	std::string prior;
	ScopedEnv(const char *n, const char *value) : name(n) {
		const char *p = ::getenv(n);
		had_prior = p != nullptr;
		if (had_prior) {
			prior = p;
		}
		if (value) {
			::setenv(n, value, 1);
		} else {
			::unsetenv(n);
		}
	}
	~ScopedEnv() {
		if (had_prior) {
			::setenv(name.c_str(), prior.c_str(), 1);
		} else {
			::unsetenv(name.c_str());
		}
	}
};

std::string make_test_tmpdir() {
	char tmpl[] = "/tmp/miint_test_gpl_XXXXXX";
	const char *dir = ::mkdtemp(tmpl);
	REQUIRE(dir != nullptr);
	return std::string(dir);
}

void rm_rf(const std::string &path) {
	// Cheap recursive remove via shell; tests own the path so the injection
	// surface is nil. Failure is non-fatal — we just leak the tmpdir.
	const std::string cmd = "rm -rf '" + path + "'";
	(void)::system(cmd.c_str());
}
} // namespace

TEST_CASE("MiintGplBoundaryCacheDir prefers XDG_CACHE_HOME", "[gpl-boundary][process]") {
	ScopedEnv xdg("XDG_CACHE_HOME", "/tmp/xdg-test");
	ScopedEnv home("HOME", "/tmp/home-test");
	REQUIRE(MiintGplBoundaryCacheDir() == "/tmp/xdg-test/miint/bin");
}

TEST_CASE("MiintGplBoundaryCacheDir falls back to HOME when XDG unset", "[gpl-boundary][process]") {
	ScopedEnv xdg("XDG_CACHE_HOME", nullptr);
	ScopedEnv home("HOME", "/tmp/home-test");
	REQUIRE(MiintGplBoundaryCacheDir() == "/tmp/home-test/.cache/miint/bin");
}

TEST_CASE("MiintGplBoundaryCacheDir returns empty when neither HOME nor XDG set", "[gpl-boundary][process]") {
	ScopedEnv xdg("XDG_CACHE_HOME", nullptr);
	ScopedEnv home("HOME", nullptr);
	REQUIRE(MiintGplBoundaryCacheDir().empty());
}

TEST_CASE("MiintGplBoundaryCacheBinary appends /gpl-boundary", "[gpl-boundary][process]") {
	ScopedEnv xdg("XDG_CACHE_HOME", "/tmp/xdg-test");
	ScopedEnv home("HOME", "/tmp/home-test");
	REQUIRE(MiintGplBoundaryCacheBinary() == "/tmp/xdg-test/miint/bin/gpl-boundary");
}

TEST_CASE("FindGplBoundary honors MIINT_GPL_BOUNDARY_PATH override", "[gpl-boundary][process]") {
	// Override is intentionally NOT validated to exist — the caller surfaces
	// a clearer error than "not found in any location" when it points at a
	// dead path. Test verifies the literal value comes back regardless of
	// what's in cache or PATH.
	ScopedEnv override_env("MIINT_GPL_BOUNDARY_PATH", "/this/does/not/exist/gpl-boundary");
	REQUIRE(FindGplBoundary() == "/this/does/not/exist/gpl-boundary");
}

TEST_CASE("FindGplBoundary finds binary in cache dir when not on PATH", "[gpl-boundary][process]") {
	const std::string tmpdir = make_test_tmpdir();
	const std::string bin_dir = tmpdir + "/miint/bin";
	const std::string bin_path = bin_dir + "/gpl-boundary";

	// Drop a fake +x file at the expected cache location.
	const std::string mkdir_cmd = "mkdir -p '" + bin_dir + "'";
	REQUIRE(::system(mkdir_cmd.c_str()) == 0);
	const std::string write_cmd =
	    "printf '#!/bin/sh\\necho fake\\n' > '" + bin_path + "' && chmod +x '" + bin_path + "'";
	REQUIRE(::system(write_cmd.c_str()) == 0);

	// Point cache-dir resolution at the tmpdir; clear the override path so
	// rule (1) doesn't short-circuit. PATH stays whatever the host has —
	// FindGplBoundary should return the cache hit BEFORE consulting PATH.
	ScopedEnv override_env("MIINT_GPL_BOUNDARY_PATH", nullptr);
	ScopedEnv xdg("XDG_CACHE_HOME", tmpdir.c_str());

	const std::string result = FindGplBoundary();
	INFO("Expected cache hit at " << bin_path << "; got '" << result << "'");
	REQUIRE(result == bin_path);

	rm_rf(tmpdir);
}

TEST_CASE("FindGplBoundary falls back to PATH when cache empty", "[gpl-boundary][process]") {
	// Point cache resolution at a dir that doesn't contain gpl-boundary.
	const std::string tmpdir = make_test_tmpdir();
	ScopedEnv override_env("MIINT_GPL_BOUNDARY_PATH", nullptr);
	ScopedEnv xdg("XDG_CACHE_HOME", tmpdir.c_str());

	// Result depends on whether the host has gpl-boundary on PATH. We're
	// only verifying the lookup falls through cleanly (doesn't crash, returns
	// the same answer FindExecutableInPath("gpl-boundary") would).
	const std::string result = FindGplBoundary();
	const std::string path_only = FindExecutableInPath("gpl-boundary");
	REQUIRE(result == path_only);

	rm_rf(tmpdir);
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
// Cycle 1.2b — ReapAndDescribe: signal-decoded termination diagnostics
// =============================================================================
//
// Session::Submit's EPIPE/EOF paths use ReapAndDescribe to turn a bare "Broken
// pipe" into "the daemon was killed by signal 15". These pin the decoding so a
// SIGTERM (the suspected PR_SET_PDEATHSIG cause) reports distinctly from a
// SIGKILL (OOM) or a normal exit — that distinction is the whole point.

TEST_CASE("ReapAndDescribe decodes a SIGTERM kill", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	std::vector<std::string> argv = {bash_path, "-c", "sleep 60"};
	ChildProcess child(argv);
	REQUIRE(child.pid() > 0);
	// Let the child reach `sleep` so the signal lands on a running process, not
	// a racing execve (which could surface as a different status).
	std::this_thread::sleep_for(std::chrono::milliseconds(50));
	REQUIRE(::kill(child.pid(), SIGTERM) == 0);
	const std::string desc = child.ReapAndDescribe();
	INFO("ReapAndDescribe: " << desc);
	REQUIRE_THAT(desc, Catch::Matchers::ContainsSubstring("signal 15"));
}

TEST_CASE("ReapAndDescribe decodes a normal exit code", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	std::vector<std::string> argv = {bash_path, "-c", "exit 7"};
	ChildProcess child(argv);
	const std::string desc = child.ReapAndDescribe(); // grace window covers the quick exit
	INFO("ReapAndDescribe: " << desc);
	REQUIRE_THAT(desc, Catch::Matchers::ContainsSubstring("exited with code 7"));
}

TEST_CASE("ReapAndDescribe is idempotent and leaves Wait consistent", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	std::vector<std::string> argv = {bash_path, "-c", "exit 3"};
	ChildProcess child(argv);
	const std::string first = child.ReapAndDescribe();
	REQUIRE_THAT(first, Catch::Matchers::ContainsSubstring("exited with code 3"));
	// Second call returns the cached description; Wait() must not block or differ.
	REQUIRE(child.ReapAndDescribe() == first);
	const int status = child.Wait();
	REQUIRE(WIFEXITED(status));
	REQUIRE(WEXITSTATUS(status) == 3);
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

// =============================================================================
// Cycle 1.4 — Hosts that live beside other hosts and clean up after themselves
// =============================================================================
//
// skb-arrow's host shuts down on stdin EOF, notices a dead caller by its
// stdout losing its reader, and must be killable with everything it started.
// Several can run at once, from DuckDB worker threads. These pin the
// ChildProcess behavior that makes that safe.

#include <fcntl.h>
#include <cerrno>
#include <memory>
#include <stdexcept>

#include "process_test_helpers.hpp"

using miint_test::EofWithin;
using miint_test::Witness;

TEST_CASE("ChildProcess pipe ends do not leak into a later child", "[gpl-boundary][process]") {
	// Each child's parent-side pipe ends must be close-on-exec. Otherwise a
	// child spawned later inherits them, and the first child never sees stdin
	// EOF while the later one lives: skb-arrow's close-stdin shutdown would
	// hang, and its dead-caller check would never fire.
	const std::string cat_path = FindExecutableInPath("cat");
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(cat_path.empty());
	REQUIRE_FALSE(bash_path.empty());
	ChildProcess first({cat_path});
	ChildProcess later({bash_path, "-c", "sleep 30"});
	first.CloseStdin();
	const std::string desc = first.ReapAndDescribe(2000);
	INFO("first child: " << desc);
	REQUIRE_THAT(desc, Catch::Matchers::ContainsSubstring("exited with code 0"));
}

TEST_CASE("ChildProcess wires the child's stdin when the parent's fd 0 is closed", "[gpl-boundary][process]") {
	// A DuckDB embedded in a daemon may run with fd 0 closed. pipe() then hands
	// out fd 0 itself, and the child must still find its stdin pipe on fd 0.
	const std::string cat_path = FindExecutableInPath("cat");
	REQUIRE_FALSE(cat_path.empty());
	std::string out;
	{
		const int saved = ::dup(STDIN_FILENO);
		REQUIRE(saved >= 0);
		struct Restore {
			int saved;
			~Restore() {
				::dup2(saved, STDIN_FILENO);
				::close(saved);
			}
		} restore {saved};
		::close(STDIN_FILENO);
		ChildProcess child({cat_path});
		WriteLine(child.stdin_fd(), "round trip");
		child.CloseStdin();
		out = drain_fd(child.stdout_fd());
		child.Wait();
	}
	REQUIRE(out == "round trip\n");
}

TEST_CASE("ChildProcess new_session makes the child lead its own process group", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ChildProcessOptions opts;
	opts.new_session = true;
	ChildProcess child({bash_path, "-c", "sleep 30"}, opts);
	REQUIRE(::getpgid(child.pid()) == child.pid());
}

TEST_CASE("ChildProcess by default stays in the caller's process group", "[gpl-boundary][process]") {
	// The control for the test above: without new_session the group is ours.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ChildProcess child({bash_path, "-c", "sleep 30"});
	REQUIRE(::getpgid(child.pid()) == ::getpgid(0));
}

TEST_CASE("KillGroup kills the child and everything it started", "[gpl-boundary][process]") {
	// Cancelling a query must stop the host and anything it started: a
	// grandchild left behind keeps the host's pipes (and its memory) alive.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	Witness witness;
	ChildProcessOptions opts;
	opts.new_session = true;
	ChildProcess child({bash_path, "-c", "sleep 30 & echo started; wait"}, opts);
	witness.CloseWriteEnd();
	LineReader reader(child.stdout_fd());
	std::string line;
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "started");
	child.KillGroup();
	// The witness reaches EOF only once bash AND its `sleep` have exited.
	const bool all_gone = EofWithin(witness.read_fd, 2000);
	REQUIRE(all_gone);
	const int status = child.Wait();
	REQUIRE(WIFSIGNALED(status));
	REQUIRE(WTERMSIG(status) == SIGKILL);
}

TEST_CASE("KillGroup requires new_session", "[gpl-boundary][process]") {
	// Without its own group, killpg would hit the caller's group: DuckDB itself.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ChildProcess child({bash_path, "-c", "sleep 30"});
	REQUIRE_THROWS_AS(child.KillGroup(), std::logic_error);
}

TEST_CASE("Destroying a new_session ChildProcess kills its whole group", "[gpl-boundary][process]") {
	// An exception unwinding past a host must not leave its descendants running.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	Witness witness;
	{
		ChildProcessOptions opts;
		opts.new_session = true;
		ChildProcess child({bash_path, "-c", "sleep 30 & echo started; wait"}, opts);
		witness.CloseWriteEnd();
		LineReader reader(child.stdout_fd());
		std::string line;
		REQUIRE(reader.ReadLine(line));
		REQUIRE(line == "started");
	}
	REQUIRE(EofWithin(witness.read_fd, 2000));
}

TEST_CASE("ChildProcessOptions set_env and unset_env shape the child's environment", "[gpl-boundary][process]") {
	// skb-arrow's thread caps travel in the host's environment, and a
	// Python-hosted DuckDB's PYTHONHOME/PYTHONPATH must not reach the host.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ScopedEnv keep("MIINT_TEST_KEEP", "kept");
	ScopedEnv replace("MIINT_TEST_REPLACE", "old");
	ScopedEnv drop("MIINT_TEST_DROP", "present");
	ChildProcessOptions opts;
	opts.set_env = {{"MIINT_TEST_REPLACE", "new"}, {"MIINT_TEST_ADD", "3"}};
	opts.unset_env = {"MIINT_TEST_DROP"};
	const RunResult r =
	    RunCapture({bash_path, "-c",
	                "echo \"$MIINT_TEST_KEEP|$MIINT_TEST_REPLACE|$MIINT_TEST_ADD|${MIINT_TEST_DROP-unset}\"; "
	                "env | grep -c '^MIINT_TEST_REPLACE='"},
	               opts);
	INFO("stderr: " << r.stderr_text);
	// The count proves the replaced variable appears once, not twice.
	REQUIRE(r.stdout_text == "kept|new|3|unset\n1\n");
}

TEST_CASE("RunCapture drains stderr while the child still writes", "[gpl-boundary][process]") {
	// 200 KB on stderr before anything on stdout. Reading stdout to EOF first
	// leaves the child blocked on a full stderr pipe, and both sides wait
	// forever; uv's installer is chatty enough on stderr to hit this.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	const RunResult r = RunCapture({bash_path, "-c", "head -c 200000 /dev/zero >&2; echo ok"});
	REQUIRE(r.spawned);
	REQUIRE(r.exited_normally);
	REQUIRE(r.exit_code == 0);
	REQUIRE(r.stdout_text == "ok\n");
	REQUIRE(r.stderr_text.size() == 200000);
}

#if defined(__linux__)
namespace {
// Spawn on a thread that exits right after, as a DuckDB worker thread may.
std::unique_ptr<ChildProcess> SpawnOnExitingThread(const std::vector<std::string> &argv,
                                                   const ChildProcessOptions &opts) {
	std::unique_ptr<ChildProcess> child;
	std::thread spawner([&]() { child = std::make_unique<ChildProcess>(argv, opts); });
	spawner.join();
	return child;
}
} // namespace

TEST_CASE("Without parent_death_signal a child outlives the thread that spawned it", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ChildProcessOptions opts;
	opts.parent_death_signal = false;
	auto child = SpawnOnExitingThread({bash_path, "-c", "sleep 30"}, opts);
	REQUIRE_THAT(child->ReapAndDescribe(300), Catch::Matchers::ContainsSubstring("still running"));
}

TEST_CASE("By default a child dies with the thread that spawned it", "[gpl-boundary][process]") {
	// The control for the test above, and the hazard it avoids: prctl(2)'s
	// PR_SET_PDEATHSIG fires when the forking *thread* exits.
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	auto child = SpawnOnExitingThread({bash_path, "-c", "sleep 30"}, {});
	REQUIRE_THAT(child->ReapAndDescribe(2000), Catch::Matchers::ContainsSubstring("signal 15"));
}
#endif

// =============================================================================
// Cycle 1.5 — LineReader for a caller that polls between reads
// =============================================================================
//
// skb-arrow calls can run for minutes. The caller polls the host's stdout so
// it can check for a cancelled query between waits, and must never block in
// ReadLine on half a line.

namespace {
// FillOnce until a line (or EOF) is buffered, at most `max_reads` times.
bool FillUntilLine(LineReader &reader, int max_reads) {
	for (int i = 0; i < max_reads && !reader.HasLine(); ++i) {
		reader.FillOnce();
	}
	return reader.HasLine();
}
} // namespace

TEST_CASE("LineReader::FillOnce buffers part of a line without blocking for the rest", "[gpl-boundary][process]") {
	const std::string bash_path = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash_path.empty());
	ChildProcess child({bash_path, "-c", "printf 'par'; sleep 0.3; printf 'tial\\n'; printf 'next\\n'"});
	LineReader reader(child.stdout_fd());
	REQUIRE_FALSE(reader.HasLine());
	reader.FillOnce(); // "par": returns without waiting for the newline
	REQUIRE_FALSE(reader.HasLine());
	REQUIRE(FillUntilLine(reader, 100));
	std::string line;
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "partial");
	REQUIRE(FillUntilLine(reader, 100));
	REQUIRE(reader.ReadLine(line));
	REQUIRE(line == "next");
	// At EOF, HasLine is true so the caller's ReadLine reports EOF at once.
	REQUIRE(FillUntilLine(reader, 100));
	REQUIRE_FALSE(reader.ReadLine(line));
	child.Wait();
}
