#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"

#include <atomic>
#include <chrono>
#include <csignal>
#include <cstring>
#include <fcntl.h>
#include <fstream>
#include <stdexcept>
#include <string>
#include <sys/mman.h>
#include <sys/stat.h>
#include <thread>
#include <unistd.h>
#include <vector>

using namespace duckdb::miint::gpl_boundary;
using Catch::Matchers::ContainsSubstring;

// =============================================================================
// Cycle 1.4 — Init handshake + shutdown (Session)
// =============================================================================
//
// Session orchestrates the gpl-boundary daemon's session lifecycle:
//
//   parent ---> {init: {...}}                   stdin
//   parent <--- {success:true, protocol_version:3, tools:[...]}   stdout
//   parent ---> {shutdown: true}                stdin
//   child  exits 0
//
// We don't have gpl-boundary on every dev box, so we test against tiny bash
// shims that mimic the protocol exactly.
//
// Invariants:
//  - Initialize() succeeds when child responds with protocol_version=3.
//  - Initialize() throws clearly when protocol_version != 3.
//  - Initialize() throws clearly when child responds with success=false.
//  - Initialize() throws clearly when child closes stdin without responding.
//  - Initialize() populates `tools()` from the init reply's `tools` array.
//  - Shutdown() sends {shutdown:true} and waits for the child to exit cleanly.

namespace {
// Wrap a bash -c invocation as a ChildProcess.
ChildProcess spawn_shim(const std::string &script) {
	const std::string bash = FindExecutableInPath("bash");
	REQUIRE_FALSE(bash.empty());
	std::vector<std::string> argv = {bash, "-c", script};
	return ChildProcess(argv);
}
} // namespace

TEST_CASE("Session::Initialize succeeds on protocol_version 3", "[gpl-boundary][session]") {
	// Read one line (the init), echo a successful reply, then read another
	// (the shutdown) and exit.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
read -r shutdown_line
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE_NOTHROW(session.Shutdown());
}

TEST_CASE("Session::Initialize throws on protocol_version != 3", "[gpl-boundary][session]") {
	// Use v2 (the previous required version) so the rejection message
	// makes sense to a reader debugging a stale-daemon scenario.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":2}'
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_THROWS_WITH(session.Initialize(), ContainsSubstring("protocol_version") && ContainsSubstring("3"));
}

TEST_CASE("Session::Initialize throws on explicit error response", "[gpl-boundary][session]") {
	const std::string script =
	    R"(read -r init_line
echo '{"success":false,"error":"daemon was rude"}'
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_THROWS_WITH(session.Initialize(), ContainsSubstring("daemon was rude"));
}

TEST_CASE("Session::Initialize throws when child closes stdin without responding", "[gpl-boundary][session]") {
	// Child reads the init line and then exits without writing anything.
	const std::string script = R"(read -r init_line
exit 1)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_THROWS(session.Initialize());
}

TEST_CASE("Session::Initialize is idempotent (second call is a no-op)", "[gpl-boundary][session]") {
	// Shim must NOT respond to a second init line — if Initialize is buggy
	// and re-sends, the test will hang on the second `read -r`. The shim
	// only handshakes once, then waits for shutdown.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
read -r shutdown_line
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE_NOTHROW(session.Initialize()); // must not double-send
	REQUIRE(session.initialized());
	REQUIRE_NOTHROW(session.Shutdown());
}

TEST_CASE("Session::tools is populated from init reply", "[gpl-boundary][session]") {
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
read -r shutdown_line
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());

	// All three tools the bowtie2 migration depends on land in the registry.
	REQUIRE(session.has_tool("bowtie2-align"));
	REQUIRE(session.has_tool("bowtie2-build"));
	REQUIRE(session.has_tool("fasttree"));
	// Schema version round-trips so Bind validators can compare against
	// the version their code was written for.
	REQUIRE(session.tool_schema_version("bowtie2-align") == 2);
	REQUIRE(session.tool_schema_version("fasttree") == 2);
	// Absent tools and absent versions both map to false/0, not throws.
	REQUIRE_FALSE(session.has_tool("nonexistent-tool"));
	REQUIRE(session.tool_schema_version("nonexistent-tool") == 0);
	REQUIRE(session.tools().size() == 3);

	REQUIRE_NOTHROW(session.Shutdown());
}

TEST_CASE("Session::tools is empty when init reply omits the field", "[gpl-boundary][session]") {
	// Defensive: a v3 daemon ought to send `tools`, but a malformed reply
	// should not crash Initialize() — empty registry is the right
	// behavior, and downstream `has_tool()` returns false.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3}'
read -r shutdown_line
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE(session.tools().empty());
	REQUIRE_FALSE(session.has_tool("bowtie2-align"));
	REQUIRE_NOTHROW(session.Shutdown());
}

TEST_CASE("Session::Shutdown is idempotent and survives early child exit", "[gpl-boundary][session]") {
	// Mock a daemon that handshakes correctly, then exits before we send
	// shutdown. Shutdown() must not throw on a missing pipe.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	// Give the bash script time to exit so that Shutdown() faces a stdin pipe
	// whose other end is already gone (write would EPIPE).
	std::this_thread::sleep_for(std::chrono::milliseconds(50));
	REQUIRE_NOTHROW(session.Shutdown());
	REQUIRE_NOTHROW(session.Shutdown()); // idempotent
}

TEST_CASE("Session against real gpl-boundary --version round-trip", "[gpl-boundary][session][!mayfail]") {
	const std::string gpl = FindGplBoundary();
	if (gpl.empty()) {
		SKIP("gpl-boundary not on PATH");
	}
	std::vector<std::string> argv = {gpl}; // no args → daemon mode
	ChildProcess child(argv);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE_NOTHROW(session.Shutdown());
}

// =============================================================================
// Cycle 4.1 — Submit-batch round trip (mock daemon)
// =============================================================================
//
// Strategy: pre-create a shm output segment with known bytes (using
// `fake_daemon_segment`-style helper, mirroring Phase 2). Spawn a bash shim
// that handshakes init, captures the batch JSON to a side file for inspection,
// and replies with a fixed batch response pointing to the pre-created segment.
//
// The test verifies:
//   1. The batch JSON we sent contains all required fields (incl. shm_input_size).
//   2. Submit() returns a SubmitResult with the expected output bytes.
//   3. Output shm is unlinked when SubmitResult is destroyed.

namespace {

// Same name-uniqueness machinery as test_gpl_boundary_shm.cpp uses, so
// concurrent CI shards / repeat runs don't collide.
std::string make_test_segment_name(const char *label) {
	static std::atomic<unsigned> counter {0};
	const unsigned n = counter.fetch_add(1, std::memory_order_relaxed);
	std::string name = "/miint-test-session-";
	name += std::to_string(static_cast<unsigned long>(::getpid()));
	name += "-";
	name += std::to_string(n);
	name += "-";
	name += label;
	return name;
}

// Pre-create an output-style shm segment containing `payload`.
struct DaemonStubSegment {
	std::string name;
	size_t size;
};
DaemonStubSegment fake_daemon_segment(const std::string &name, const std::string &payload) {
	const int fd = ::shm_open(name.c_str(), O_CREAT | O_EXCL | O_RDWR, 0600);
	if (fd < 0) {
		throw std::runtime_error("fake_daemon_segment shm_open failed for " + name + ": " +
		                         std::string(::strerror(errno)));
	}
	if (::ftruncate(fd, static_cast<off_t>(payload.size())) != 0) {
		::close(fd);
		::shm_unlink(name.c_str());
		throw std::runtime_error("fake_daemon_segment ftruncate failed");
	}
	void *p = ::mmap(nullptr, payload.size(), PROT_READ | PROT_WRITE, MAP_SHARED, fd, 0);
	::close(fd);
	if (p == MAP_FAILED) {
		::shm_unlink(name.c_str());
		throw std::runtime_error("fake_daemon_segment mmap failed");
	}
	std::memcpy(p, payload.data(), payload.size());
	::munmap(p, payload.size());
	return {name, payload.size()};
}

// Read entire contents of a file as a string.
std::string slurp(const std::string &path) {
	std::ifstream f(path);
	std::string out((std::istreambuf_iterator<char>(f)), std::istreambuf_iterator<char>());
	return out;
}

} // namespace

TEST_CASE("Session::Submit happy path: batch JSON and shm round trip", "[gpl-boundary][session]") {
	// 1. Pre-create the output segment that the mock daemon will point miint at.
	const std::string out_name = make_test_segment_name("submit-basic");
	const std::string payload = "this is a fake fasttree result blob for the test";
	auto stub = fake_daemon_segment(out_name, payload);

	// 2. Build the bash shim. It captures the batch line for assertions.
	const std::string capture_path =
	    "/tmp/miint-test-batch-line-" + std::to_string(::getpid()) + "-" + out_name.substr(1);
	const std::string response_template = R"({"success":true,"schema_version":2,"batch_id":1,)"
	                                      R"("shm_outputs":[{"name":")" +
	                                      out_name + R"(","label":"tree","size":)" + std::to_string(stub.size) +
	                                      R"(}],"result":{"n_nodes":7}})";
	std::string script;
	script += "read -r init_line\n";
	script +=
	    R"(echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}')";
	script += "\n";
	script += "read -r batch_line\n";
	script += "printf '%s\\n' \"$batch_line\" > '" + capture_path + "'\n";
	script += "echo '" + response_template + "'\n";
	script += "read -r shutdown_line\n";
	script += "exit 0\n";

	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());

	// 3. Submit. Input bytes are arbitrary (mock daemon doesn't read shm input).
	const std::string input_bytes = "ARROW_IPC_PLACEHOLDER_BYTES";
	SubmitResult result;
	REQUIRE_NOTHROW(result = session.Submit("fasttree", R"({"seed":42})", input_bytes.data(), input_bytes.size()));

	// 4. Validate the response we got.
	REQUIRE(result.batch_id == 1);
	REQUIRE(result.schema_version == 2);
	REQUIRE(result.outputs.size() == 1);
	REQUIRE(result.outputs[0].label == "tree");
	REQUIRE(result.outputs[0].size_bytes() == stub.size);
	REQUIRE(std::memcmp(result.outputs[0].bytes(), payload.data(), payload.size()) == 0);
	REQUIRE_FALSE(result.result_json.empty());
	REQUIRE(result.result_json.find("\"n_nodes\"") != std::string::npos);

	// 5. Inspect the captured batch JSON. All required fields must be present.
	const std::string captured = slurp(capture_path);
	INFO("captured batch line: " << captured);
	REQUIRE(captured.find("\"tool\":\"fasttree\"") != std::string::npos);
	REQUIRE(captured.find("\"config\"") != std::string::npos);
	REQUIRE(captured.find("\"shm_input\"") != std::string::npos);
	// The critical post-19306f6 field — without it, the daemon errors out.
	REQUIRE(captured.find("\"shm_input_size\"") != std::string::npos);
	REQUIRE(captured.find("\"batch_id\"") != std::string::npos);
	REQUIRE(captured.find("\"shm_input_size\":" + std::to_string(input_bytes.size())) != std::string::npos);

	REQUIRE_NOTHROW(session.Shutdown());

	// 6. After SubmitResult is destroyed (it goes out of scope at end of test),
	//    the output segment must be unlinked. Verify here while it's alive
	//    that the segment is there; the unlink-on-destroy property is already
	//    covered by `OutputShmRegion is auto-unlinked on destruction`.
	::unlink(capture_path.c_str());
}

TEST_CASE("Session::Submit opts into per-batch metrics only when requested", "[gpl-boundary][session]") {
	// gpl-boundary v0.4.2 added an opt-in `metrics` object (bowtie2 worker
	// getrusage: ru_majflt, CPU/faults/RSS, worker_reused). It is gated on a
	// top-level `"metrics":true` in the BatchRequest — non-opted batches pay
	// nothing (no extra syscalls daemon-side). miint must send the flag ONLY when
	// its telemetry is enabled, so an ordinary run is unaffected. This pins the
	// wire contract: request_metrics=true emits `"metrics":true`; the default
	// omits the key ENTIRELY (not `"metrics":false`), so a pre-0.4.2 daemon — which
	// ignores unknown request fields — sees exactly the request shape it always did.
	const std::string out_a = make_test_segment_name("metrics-on");
	const std::string out_b = make_test_segment_name("metrics-off");
	const std::string payload = "x"; // tiny output; the assertion is on the REQUEST line
	auto stub_a = fake_daemon_segment(out_a, payload);
	auto stub_b = fake_daemon_segment(out_b, payload);

	const std::string capture_path =
	    "/tmp/miint-test-metrics-line-" + std::to_string(::getpid()) + "-" + out_a.substr(1);
	::unlink(capture_path.c_str()); // start clean — the shim appends two lines

	auto resp = [](const std::string &name, size_t sz, int64_t bid) {
		return R"({"success":true,"schema_version":2,"batch_id":)" + std::to_string(bid) +
		       R"(,"shm_outputs":[{"name":")" + name + R"(","label":"sam","size":)" + std::to_string(sz) + R"(}]})";
	};
	std::string script;
	script += "read -r init_line\n";
	script += R"(echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2}]}')";
	script += "\n";
	script += "read -r batch_line\n";
	script += "printf '%s\\n' \"$batch_line\" >> '" + capture_path + "'\n";
	script += "echo '" + resp(out_a, stub_a.size, 1) + "'\n";
	script += "read -r batch_line\n";
	script += "printf '%s\\n' \"$batch_line\" >> '" + capture_path + "'\n";
	script += "echo '" + resp(out_b, stub_b.size, 2) + "'\n";
	script += "read -r shutdown_line\n";
	script += "exit 0\n";

	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());

	const std::string input_bytes = "ARROW_IPC_PLACEHOLDER";
	// Hold both results alive so neither output segment is unlinked early.
	SubmitResult r_on = session.Submit("bowtie2-align", R"({"index_path":"x"})", input_bytes.data(), input_bytes.size(),
	                                   /*request_metrics=*/true);
	SubmitResult r_off =
	    session.Submit("bowtie2-align", R"({"index_path":"x"})", input_bytes.data(), input_bytes.size());
	REQUIRE_NOTHROW(session.Shutdown());

	const std::string captured = slurp(capture_path);
	const auto nl = captured.find('\n');
	REQUIRE(nl != std::string::npos);
	const std::string line_on = captured.substr(0, nl);
	const std::string line_off = captured.substr(nl + 1);
	INFO("opted-in line: " << line_on);
	INFO("default line:  " << line_off);
	REQUIRE(line_on.find("\"metrics\":true") != std::string::npos);
	REQUIRE(line_off.find("\"metrics\"") == std::string::npos);

	::unlink(capture_path.c_str());
}

TEST_CASE("Session::Submit propagates daemon error responses", "[gpl-boundary][session]") {
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
read -r batch_line
echo '{"success":false,"error":"tool blew up","batch_id":1}'
read -r shutdown_line
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE_THROWS_WITH(session.Submit("fasttree", "{}", "x", 1), Catch::Matchers::ContainsSubstring("tool blew up"));
}

TEST_CASE("Session::Submit reports the daemon's death signal when it is killed mid-batch", "[gpl-boundary][session]") {
	// Field repro: the daemon handshakes, consumes our batch line, then dies by
	// signal (mimicking PR_SET_PDEATHSIG's SIGTERM firing when the DuckDB worker
	// thread that forked it unwinds) WITHOUT replying. The pre-fix error was a
	// bare "daemon closed stdout while waiting for batch response", which says
	// nothing about *why* the daemon vanished. Submit must now splice in the
	// signal so the cause is visible in the surfaced SQL error.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2}]}'
read -r batch_line
kill -TERM $$
sleep 5)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	REQUIRE_THROWS_WITH(session.Submit("bowtie2-align", "{}", "x", 1), ContainsSubstring("signal 15"));
}

TEST_CASE("Session::Submit drains daemon stderr each batch so a verbose daemon can't deadlock",
          "[gpl-boundary][session]") {
	// Field repro for the align_bowtie2_sharded hang: with quiet=false, bowtie2
	// prints a ~20-line alignment summary to stderr on every batch. The daemon's
	// stderr is a pipe back to miint that was only ever read on the failure
	// paths — never during normal Submit traffic. After enough batches the
	// daemon's writes crossed the ~64 KB OS pipe buffer, the daemon blocked in
	// write(2) before it could answer, and Submit hung forever in ReadLine (CPU
	// flatlined; no crash). Submit now drains stderr after every response so the
	// pipe can never back up.
	//
	// This shim writes 48 KB to stderr per batch across 8 batches (384 KB, far
	// past 64 KB). Each single batch's 48 KB fits an empty pipe, so the fix
	// (drain-per-batch) keeps every write unblocked; without it the pipe fills
	// by batch ~2 and the shim wedges. A watchdog bounds the wait and SIGKILLs
	// the shim on timeout so a regression fails loudly instead of hanging CI.
	const int kBatches = 8;
	std::string script =
	    R"SHIM(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2}]}'
bid=1
while [ $bid -le )SHIM" +
	    std::to_string(kBatches) + R"SHIM( ]; do
  read -r batch_line || break
  head -c 49152 /dev/zero | tr '\0' 'x' >&2
  echo "{\"success\":true,\"schema_version\":2,\"batch_id\":$bid}"
  bid=$((bid+1))
done
read -r shutdown_line
exit 0
)SHIM";

	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	const pid_t shim_pid = session.daemon_pid();

	std::atomic<bool> done {false};
	std::atomic<int> completed {0};
	std::thread worker([&]() {
		try {
			for (int i = 0; i < kBatches; ++i) {
				session.Submit("bowtie2-align", "{}", "x", 1);
				completed.fetch_add(1, std::memory_order_relaxed);
			}
		} catch (...) {
			// Watchdog SIGKILL (or any error) surfaces here as a thrown
			// exception once ReadLine sees the closed pipe. `completed` already
			// records how far we got; let the assertions below judge.
		}
		done.store(true, std::memory_order_release);
	});

	const auto deadline = std::chrono::steady_clock::now() + std::chrono::seconds(15);
	while (!done.load(std::memory_order_acquire) && std::chrono::steady_clock::now() < deadline) {
		std::this_thread::sleep_for(std::chrono::milliseconds(20));
	}
	const bool finished_in_time = done.load(std::memory_order_acquire);
	if (!finished_in_time) {
		// Pre-fix: the worker is blocked in ReadLine because the shim is blocked
		// in write(2,stderr). Killing the shim closes its stdout so ReadLine
		// returns and the worker can be joined cleanly.
		::kill(shim_pid, SIGKILL);
	}
	worker.join();

	REQUIRE(finished_in_time); // false on the pre-fix deadlock
	REQUIRE(completed.load() == kBatches);
	REQUIRE_NOTHROW(session.Shutdown());
}

TEST_CASE("Session::Submit rejects calls before Initialize", "[gpl-boundary][session]") {
	const std::string script = R"(sleep 5
exit 0)"; // shim that never responds
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_THROWS(session.Submit("fasttree", "{}", "x", 1));
}

TEST_CASE("Session::Submit rejects malformed config_json", "[gpl-boundary][session]") {
	// Mock that handshakes correctly but is never expected to receive a batch.
	const std::string script =
	    R"(read -r init_line
echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}'
sleep 5
exit 0)";
	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	// Caller-side validation must fail BEFORE any bytes hit the daemon.
	REQUIRE_THROWS(session.Submit("fasttree", "not json at all", "x", 1));
	// Object-injection attempt: a stray closing brace turns this into two
	// top-level values, which would let a malicious caller smuggle fields
	// into the batch envelope.
	REQUIRE_THROWS(session.Submit("fasttree", R"({"seed":42},"injected":"value")", "x", 1));
	// JSON arrays are valid JSON but not valid configs (we expect an object).
	REQUIRE_THROWS(session.Submit("fasttree", R"([1,2,3])", "x", 1));
}

TEST_CASE("Session::Submit increments batch_id across consecutive calls "
          "and reuses the same daemon process",
          "[gpl-boundary][session]") {
	// Cycle 4.2: same Session, two Submit calls — daemon PID must be stable
	// (same child process), and our batch_id must increment.
	const std::string out_name1 = make_test_segment_name("multi-1");
	const std::string out_name2 = make_test_segment_name("multi-2");
	auto stub1 = fake_daemon_segment(out_name1, "first response");
	auto stub2 = fake_daemon_segment(out_name2, "second response");

	const std::string resp1 = R"({"success":true,"schema_version":2,"batch_id":1,"shm_outputs":[{"name":")" +
	                          out_name1 + R"(","label":"tree","size":)" + std::to_string(stub1.size) + R"(}]})";
	const std::string resp2 = R"({"success":true,"schema_version":2,"batch_id":2,"shm_outputs":[{"name":")" +
	                          out_name2 + R"(","label":"tree","size":)" + std::to_string(stub2.size) + R"(}]})";
	std::string script;
	script += "read -r init_line\n";
	script +=
	    R"(echo '{"success":true,"protocol_version":3,"tools":[{"name":"bowtie2-align","schema_version":2},{"name":"bowtie2-build","schema_version":1},{"name":"fasttree","schema_version":2}]}')";
	script += "\n";
	script += "read -r batch1\n";
	script += "echo '" + resp1 + "'\n";
	script += "read -r batch2\n";
	script += "echo '" + resp2 + "'\n";
	script += "read -r shutdown_line\n";
	script += "exit 0\n";

	auto child = spawn_shim(script);
	Session session(std::move(child));
	REQUIRE_NOTHROW(session.Initialize());
	const pid_t pid_before = session.daemon_pid();

	auto r1 = session.Submit("fasttree", "{}", "x", 1);
	REQUIRE(r1.batch_id == 1);
	REQUIRE(r1.outputs.size() == 1);

	const pid_t pid_after_first = session.daemon_pid();
	REQUIRE(pid_after_first == pid_before);

	auto r2 = session.Submit("fasttree", "{}", "y", 1);
	REQUIRE(r2.batch_id == 2);

	REQUIRE(session.daemon_pid() == pid_before);
	REQUIRE_NOTHROW(session.Shutdown());
}

// =============================================================================
// ParseGplBoundaryVersion — the bowtie2 version gate parses `--version` JSON to
// detect a daemon older than the bowtie2 `memory_mapped` minimum (0.4.2), since
// the IPC handshake can't report the release version.
// =============================================================================

TEST_CASE("ParseGplBoundaryVersion: extracts the semver from --version JSON", "[gpl-boundary][session]") {
	int major = -1, minor = -1, patch = -1;

	// Real shape: {"gpl_boundary": "0.4.2", "tools": [...]}.
	REQUIRE(ParseGplBoundaryVersion(R"({"gpl_boundary": "0.4.2", "tools": []})", major, minor, patch));
	REQUIRE(major == 0);
	REQUIRE(minor == 4);
	REQUIRE(patch == 2);

	// A higher release parses (no whitespace around the colon either).
	REQUIRE(ParseGplBoundaryVersion(R"({"gpl_boundary":"1.2.3"})", major, minor, patch));
	REQUIRE((major == 1 && minor == 2 && patch == 3));

	// A pre-release/build suffix after patch is tolerated.
	REQUIRE(ParseGplBoundaryVersion(R"({"gpl_boundary":"0.4.2-rc1"})", major, minor, patch));
	REQUIRE((major == 0 && minor == 4 && patch == 2));
}

TEST_CASE("ParseGplBoundaryVersion: unparseable output is rejected (caller fails loud)", "[gpl-boundary][session]") {
	int major = 9, minor = 9, patch = 9;
	// Field absent, plain-text (older/foreign --version), garbage value, and empty
	// all return false → the caller treats them as "can't confirm >= 0.4.2" and
	// throws rather than proceeding against an unknown daemon.
	REQUIRE_FALSE(ParseGplBoundaryVersion(R"({"tools": []})", major, minor, patch));
	REQUIRE_FALSE(ParseGplBoundaryVersion("gpl-boundary 0.4.1\n", major, minor, patch));
	REQUIRE_FALSE(ParseGplBoundaryVersion(R"({"gpl_boundary": "not-a-version"})", major, minor, patch));
	REQUIRE_FALSE(ParseGplBoundaryVersion("", major, minor, patch));
}
