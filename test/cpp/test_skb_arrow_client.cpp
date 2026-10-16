#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "skb_arrow/client.hpp"
#include "skb_arrow/install.hpp"

#include "env_test_helpers.hpp"
#include "gpl_boundary/arrow_ipc.hpp"
#include "nanoarrow/nanoarrow.h"
#include "process_test_helpers.hpp"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <dirent.h>
#include <fcntl.h>
#include <fstream>
#include <limits>
#include <sstream>
#include <string>
#include <sys/stat.h>
#include <unistd.h>
#include <vector>

using namespace duckdb::miint::skb_arrow;
namespace gb = duckdb::miint::gpl_boundary;
using Catch::Matchers::ContainsSubstring;
using miint_test::EnvGuard;
using miint_test::EofWithin;
using miint_test::Witness;

// =============================================================================
// Fixtures
// =============================================================================

namespace {

const InterruptCheck kNever = [] {
};

bool Exists(const std::string &path) {
	struct stat st;
	return ::lstat(path.c_str(), &st) == 0;
}

std::vector<std::string> Entries(const std::string &dir) {
	std::vector<std::string> out;
	DIR *d = ::opendir(dir.c_str());
	if (!d) {
		return out;
	}
	while (const dirent *e = ::readdir(d)) {
		const std::string name = e->d_name;
		if (name != "." && name != "..") {
			out.push_back(name);
		}
	}
	::closedir(d);
	return out;
}

std::string Slurp(const std::string &path) {
	std::ifstream in(path, std::ios::binary);
	std::stringstream ss;
	ss << in.rdbuf();
	return ss.str();
}

// A scratch directory under /tmp, removed with its contents at scope exit.
struct TempDir {
	std::string path;
	TempDir() {
		char tmpl[] = "/tmp/miint_test_skb_XXXXXX";
		const char *p = ::mkdtemp(tmpl);
		REQUIRE(p != nullptr);
		path = p;
	}
	~TempDir() {
		const std::string cmd = "rm -rf '" + path + "'";
		(void)::system(cmd.c_str());
	}
};

// `s` as one single-quoted shell word.
std::string ShQuote(const std::string &s) {
	std::string out = "'";
	for (const char c : s) {
		out += c == '\'' ? std::string("'\\''") : std::string(1, c);
	}
	return out + "'";
}

// The real host's ready line, as a fake host sends it.
constexpr const char *kReady = R"({"type":"ready","protocol_version":1,"host_version":"9.9.9",)"
                               R"("capabilities":{"echo":1,"mantel":1},"id":null,"warnings":[]})";

// Writes fake skb-arrow hosts: bash scripts that, like the real host, exit 2
// unless started as `--segment-dir DIR` with DIR an empty directory, and that
// log every request line to `log`.
struct FakeHosts {
	TempDir tmp;
	std::string log = tmp.path + "/requests.log";
	std::string sessions = tmp.path + "/sessions";
	FakeHosts() {
		REQUIRE(::mkdir(sessions.c_str(), 0700) == 0);
	}

	// `body` runs after the handshake (`ready` is the reply to init). $DIR is
	// the session directory; `next` reads and logs one request into $line.
	std::string Write(const std::string &body, const std::string &ready = kReady) {
		const std::string path = tmp.path + "/fake-skb-arrow";
		std::ofstream out(path);
		out << "#!/bin/bash\n"
		    << "[ \"$1\" = --segment-dir ] && [ -d \"$2\" ] && [ -z \"$(ls -A \"$2\")\" ] || exit 2\n"
		    << "DIR=\"$2\"\n"
		    << "next() { IFS= read -r line || exit 0; printf '%s\\n' \"$line\" >> '" << log << "'; }\n"
		    << "next\n"
		    << "printf '%s\\n' " << ShQuote(ready) << "\n"
		    << body << "\n";
		out.close();
		REQUIRE(::chmod(path.c_str(), 0755) == 0);
		return path;
	}

	HostOptions Options(const std::string &binary) const {
		HostOptions o;
		o.binary = binary;
		o.segment_parent = sessions;
		return o;
	}

	std::vector<std::string> Requests() const {
		std::vector<std::string> lines;
		std::istringstream in(Slurp(log));
		for (std::string line; std::getline(in, line);) {
			lines.push_back(line);
		}
		return lines;
	}
};

// Reply `response`, then idle until stdin closes and exit 0, as the real host does.
std::string Respond(const std::string &response) {
	return "printf '%s\\n' " + ShQuote(response) + "\nwhile IFS= read -r line; do :; done\nexit 0";
}

// After the handshake: answer one call with `response`.
std::string AnswerOnce(const std::string &response) {
	return "next\n" + Respond(response);
}

// nanoarrow fixtures. `int64 a, utf8 s` tables, built in batches.
struct OwnedSchema {
	ArrowSchema schema {};
	OwnedSchema() {
		REQUIRE(ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT) == NANOARROW_OK);
		REQUIRE(ArrowSchemaAllocateChildren(&schema, 2) == NANOARROW_OK);
		REQUIRE(ArrowSchemaInitFromType(schema.children[0], NANOARROW_TYPE_INT64) == NANOARROW_OK);
		REQUIRE(ArrowSchemaSetName(schema.children[0], "a") == NANOARROW_OK);
		REQUIRE(ArrowSchemaInitFromType(schema.children[1], NANOARROW_TYPE_STRING) == NANOARROW_OK);
		REQUIRE(ArrowSchemaSetName(schema.children[1], "s") == NANOARROW_OK);
	}
	~OwnedSchema() {
		if (schema.release) {
			schema.release(&schema);
		}
	}
};

struct OwnedArray {
	ArrowArray array {};
	OwnedArray(const ArrowSchema &schema, int64_t first, int64_t count) {
		ArrowError err {};
		REQUIRE(ArrowArrayInitFromSchema(&array, &schema, &err) == NANOARROW_OK);
		REQUIRE(ArrowArrayStartAppending(&array) == NANOARROW_OK);
		for (int64_t v = first; v < first + count; ++v) {
			const std::string s = "row" + std::to_string(v);
			REQUIRE(ArrowArrayAppendInt(array.children[0], v) == NANOARROW_OK);
			REQUIRE(ArrowArrayAppendString(array.children[1], {s.data(), static_cast<int64_t>(s.size())}) ==
			        NANOARROW_OK);
			REQUIRE(ArrowArrayFinishElement(&array) == NANOARROW_OK);
		}
		REQUIRE(ArrowArrayFinishBuildingDefault(&array, &err) == NANOARROW_OK);
	}
	~OwnedArray() {
		if (array.release) {
			array.release(&array);
		}
	}
};

// Decode one IPC stream: its column names and the `a` values of every batch.
struct Decoded {
	std::vector<std::string> names;
	std::vector<int64_t> a;
	int batches = 0;
};
Decoded DecodeIntTable(const std::string &bytes) {
	Decoded out;
	gb::IpcStreamDecoder decoder(bytes.data(), bytes.size());
	ArrowSchema schema {};
	decoder.GetSchema(&schema);
	for (int64_t i = 0; i < schema.n_children; ++i) {
		out.names.push_back(schema.children[i]->name);
	}
	ArrowArray batch {};
	while (decoder.NextBatch(&batch)) {
		++out.batches;
		const ArrowArray *col = batch.children[0];
		const auto *values = static_cast<const int64_t *>(col->buffers[1]);
		for (int64_t r = 0; r < batch.length; ++r) {
			out.a.push_back(values[col->offset + batch.offset + r]);
		}
		batch.release(&batch);
	}
	schema.release(&schema);
	return out;
}

} // namespace

// =============================================================================
// Discovery and placement
// =============================================================================

TEST_CASE("FindSkbArrow honors MIINT_SKB_ARROW_PATH without validating it", "[skb-arrow]") {
	// An explicit path wins, so a caller can say "that binary doesn't exist"
	// instead of "not found anywhere".
	EnvGuard override_path("MIINT_SKB_ARROW_PATH", "/this/does/not/exist/skb-arrow");
	REQUIRE(FindSkbArrow() == "/this/does/not/exist/skb-arrow");
}

TEST_CASE("FindSkbArrow finds miint's cached binary before PATH", "[skb-arrow]") {
	TempDir tmp;
	const std::string bin_dir = tmp.path + "/miint/bin";
	REQUIRE(::system(("mkdir -p '" + bin_dir + "'").c_str()) == 0);
	const std::string bin = bin_dir + "/skb-arrow";
	REQUIRE(::system(("printf '#!/bin/sh\\n' > '" + bin + "' && chmod +x '" + bin + "'").c_str()) == 0);
	EnvGuard override_path("MIINT_SKB_ARROW_PATH", nullptr);
	EnvGuard xdg("XDG_CACHE_HOME", tmp.path.c_str());
	REQUIRE(MiintSkbArrowCacheBinary() == bin);
	REQUIRE(FindSkbArrow() == bin);
}

TEST_CASE("FindSkbArrow falls back to PATH when the cache has no binary", "[skb-arrow]") {
	TempDir tmp;
	EnvGuard override_path("MIINT_SKB_ARROW_PATH", nullptr);
	EnvGuard xdg("XDG_CACHE_HOME", tmp.path.c_str());
	REQUIRE(FindSkbArrow() == gb::FindExecutableInPath("skb-arrow"));
}

TEST_CASE("SegmentDirParent honors MIINT_SKB_ARROW_SEGMENT_DIR", "[skb-arrow]") {
	// The escape from Docker's 64 MiB /dev/shm.
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", "/somewhere/roomy");
	REQUIRE(SegmentDirParent() == "/somewhere/roomy");
}

TEST_CASE("SegmentDirParent defaults to tmpfs on Linux, else TMPDIR", "[skb-arrow]") {
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", nullptr);
	EnvGuard tmpdir("TMPDIR", "/my/tmp");
	struct stat st;
	const bool shm = ::stat("/dev/shm", &st) == 0 && S_ISDIR(st.st_mode) && ::access("/dev/shm", W_OK | X_OK) == 0;
#if defined(__linux__)
	REQUIRE(SegmentDirParent() == (shm ? "/dev/shm" : "/my/tmp"));
#else
	(void)shm;
	REQUIRE(SegmentDirParent() == "/my/tmp");
#endif
}

// =============================================================================
// Session directory
// =============================================================================

TEST_CASE("SessionDir is a fresh private directory, removed with its contents", "[skb-arrow]") {
	TempDir parent;
	std::string path;
	{
		SessionDir dir(parent.path);
		path = dir.path();
		REQUIRE_THAT(path, ContainsSubstring(parent.path + "/miint-skb-"));
		struct stat st;
		REQUIRE(::stat(path.c_str(), &st) == 0);
		REQUIRE(S_ISDIR(st.st_mode));
		REQUIRE((st.st_mode & 0777) == 0700);
		REQUIRE(Entries(path).empty());
		// A segment the host never consumed is left behind on failure paths.
		std::ofstream(path + "/leftover") << "x";
		dir.Remove();
		REQUIRE_FALSE(Exists(path));
		dir.Remove(); // idempotent
	}
	REQUIRE_FALSE(Exists(path));
}

TEST_CASE("SessionDir under a missing parent is a resource error naming the override", "[skb-arrow]") {
	try {
		SessionDir dir("/this/does/not/exist");
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "resource");
		REQUIRE_THAT(e.what(), ContainsSubstring("/this/does/not/exist"));
		REQUIRE_THAT(e.what(), ContainsSubstring("MIINT_SKB_ARROW_SEGMENT_DIR"));
	}
}

// =============================================================================
// Segments
// =============================================================================

TEST_CASE("SegmentWriter rolls over to a new segment at the cap", "[skb-arrow]") {
	TempDir dir;
	OwnedSchema schema;
	OwnedArray b0(schema.schema, 0, 3), b1(schema.schema, 3, 2), b2(schema.schema, 5, 4);
	std::vector<std::string> names;
	{
		// A 1-byte cap: every batch starts a new segment.
		SegmentWriter writer(dir.path, "t", schema.schema, 1);
		writer.Append(b0.array);
		writer.Append(b1.array);
		writer.Append(b2.array);
		names = writer.Finish();
	}
	REQUIRE(names == std::vector<std::string> {"t-0", "t-1", "t-2"});
	std::vector<int64_t> all;
	for (const auto &name : names) {
		struct stat st;
		REQUIRE(::stat((dir.path + "/" + name).c_str(), &st) == 0);
		REQUIRE((st.st_mode & 0777) == 0600);
		const Decoded d = DecodeIntTable(Slurp(dir.path + "/" + name));
		REQUIRE(d.names == std::vector<std::string> {"a", "s"});
		REQUIRE(d.batches == 1);
		all.insert(all.end(), d.a.begin(), d.a.end());
	}
	REQUIRE(all == std::vector<int64_t> {0, 1, 2, 3, 4, 5, 6, 7, 8});
}

TEST_CASE("SegmentWriter keeps batches together under a roomy cap", "[skb-arrow]") {
	TempDir dir;
	OwnedSchema schema;
	OwnedArray b0(schema.schema, 0, 3), b1(schema.schema, 3, 2);
	SegmentWriter writer(dir.path, "t", schema.schema);
	writer.Append(b0.array);
	writer.Append(b1.array);
	const auto names = writer.Finish();
	REQUIRE(names == std::vector<std::string> {"t-0"});
	const Decoded d = DecodeIntTable(Slurp(dir.path + "/t-0"));
	REQUIRE(d.batches == 2);
	REQUIRE(d.a == std::vector<int64_t> {0, 1, 2, 3, 4});
}

TEST_CASE("SegmentWriter writes an empty table as one schema-only segment", "[skb-arrow]") {
	TempDir dir;
	OwnedSchema schema;
	SegmentWriter writer(dir.path, "empty", schema.schema);
	const auto names = writer.Finish();
	REQUIRE(names == std::vector<std::string> {"empty-0"});
	const Decoded d = DecodeIntTable(Slurp(dir.path + "/empty-0"));
	REQUIRE(d.names == std::vector<std::string> {"a", "s"});
	REQUIRE(d.batches == 0);
}

TEST_CASE("An unfinished SegmentWriter unlinks what it wrote", "[skb-arrow]") {
	// A failure between writing and calling must not strand segments.
	TempDir dir;
	OwnedSchema schema;
	OwnedArray b0(schema.schema, 0, 3);
	{
		SegmentWriter writer(dir.path, "t", schema.schema, 1);
		writer.Append(b0.array);
		writer.Append(b0.array);
		REQUIRE(Entries(dir.path).size() == 2);
	}
	REQUIRE(Entries(dir.path).empty());
}

TEST_CASE("SegmentWriter never overwrites, or removes, a file it did not create", "[skb-arrow]") {
	TempDir dir;
	std::ofstream(dir.path + "/t-0") << "someone else's";
	OwnedSchema schema;
	OwnedArray b0(schema.schema, 0, 3);
	{
		SegmentWriter writer(dir.path, "t", schema.schema);
		REQUIRE_THROWS_AS(writer.Append(b0.array), Error);
	}
	REQUIRE(Slurp(dir.path + "/t-0") == "someone else's");
}

TEST_CASE("SegmentWriter refuses names the protocol forbids", "[skb-arrow]") {
	// Names match [A-Za-z0-9_-][A-Za-z0-9._-]{0,127}: no path separators, no
	// leading dot, at most 128 characters.
	TempDir dir;
	OwnedSchema schema;
	for (const std::string &prefix : std::vector<std::string> {"a/b", ".hidden", std::string(130, 'x')}) {
		INFO("prefix " << prefix);
		REQUIRE_THROWS_AS(SegmentWriter(dir.path, prefix, schema.schema).Finish(), Error);
	}
	REQUIRE(Entries(dir.path).empty());
}

// =============================================================================
// Params
// =============================================================================

TEST_CASE("Params refuses a number that is not finite", "[skb-arrow]") {
	// NaN and Infinity aren't JSON; a writer that turned them into null would
	// silently mean "use the default".
	Params params;
	for (const double bad :
	     {std::nan(""), std::numeric_limits<double>::infinity(), -std::numeric_limits<double>::infinity()}) {
		try {
			params.SetNumber("tol", bad);
			FAIL("expected an Error");
		} catch (const Error &e) {
			REQUIRE(e.kind() == "invalid_param");
			REQUIRE_THAT(e.what(), ContainsSubstring("tol"));
		}
	}
	REQUIRE(params.entries().empty());
}

// =============================================================================
// Host: handshake
// =============================================================================

TEST_CASE("Host completes the handshake and reports the host version", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce("unused"));
	Host host(fakes.Options(bin), kNever);
	REQUIRE(host.host_version() == "9.9.9");
	REQUIRE(Entries(host.dir()).empty());
	REQUIRE(fakes.Requests() == std::vector<std::string> {R"({"type":"init","protocol_version":1})"});
	host.Shutdown();
}

TEST_CASE("Host reports a refused handshake with the host's kind and message", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(
	    "",
	    R"({"type":"error","kind":"host_incompatible","message":"protocol_version 1 is not 2","id":null,"warnings":[]})");
	try {
		Host host(fakes.Options(bin), kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "host_incompatible");
		REQUIRE_THAT(e.what(), ContainsSubstring("protocol_version 1 is not 2"));
	}
	REQUIRE(Entries(fakes.sessions).empty());
}

TEST_CASE("Host refuses a host without a required capability", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce("unused"));
	HostOptions options = fakes.Options(bin);
	options.required = {{"ancombc", 1}};
	try {
		Host host(options, kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "host_incompatible");
		REQUIRE_THAT(e.what(), ContainsSubstring("ancombc"));
		REQUIRE_THAT(e.what(), ContainsSubstring("9.9.9"));
		REQUIRE_THAT(e.what(), ContainsSubstring("install_skb_arrow()"));
	}
	REQUIRE(Entries(fakes.sessions).empty());
}

TEST_CASE("Host refuses a capability at another schema_version", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce("unused"));
	HostOptions options = fakes.Options(bin);
	options.required = {{"mantel", 2}};
	try {
		Host host(options, kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "host_incompatible");
		REQUIRE_THAT(e.what(), ContainsSubstring("mantel"));
		REQUIRE_THAT(e.what(), ContainsSubstring("schema_version 1"));
		REQUIRE_THAT(e.what(), ContainsSubstring("2"));
	}
}

TEST_CASE("Host reports a host that exits at startup, with its exit code and stderr", "[skb-arrow]") {
	// The real host exits 2 when it can't use its session directory; the
	// operator needs that reason, not a bare broken pipe.
	FakeHosts fakes;
	const std::string bin = fakes.tmp.path + "/exits-at-start";
	std::ofstream(bin) << "#!/bin/bash\necho 'skb-arrow: --segment-dir is not usable' >&2\nexit 2\n";
	REQUIRE(::chmod(bin.c_str(), 0755) == 0);
	try {
		Host host(fakes.Options(bin), kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "internal");
		REQUIRE_THAT(e.what(), ContainsSubstring("exited with code 2"));
		REQUIRE_THAT(e.what(), ContainsSubstring("--segment-dir is not usable"));
	}
	REQUIRE(Entries(fakes.sessions).empty());
}

TEST_CASE("Host reports a host it can't start as a broken install, not a miint bug", "[skb-arrow]") {
	// A launcher whose interpreter is gone (a removed uv-managed Python):
	// execve fails. Reinstalling fixes it; reporting it as a bug doesn't.
	FakeHosts fakes;
	const std::string bin = fakes.tmp.path + "/orphaned-launcher";
	std::ofstream(bin) << "#!/nonexistent/python3\n";
	REQUIRE(::chmod(bin.c_str(), 0755) == 0);
	try {
		Host host(fakes.Options(bin), kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "host_incompatible");
		REQUIRE_THAT(e.what(), ContainsSubstring(bin));
		REQUIRE_THAT(e.what(), ContainsSubstring("SELECT install_skb_arrow(true)"));
	}
	REQUIRE(Entries(fakes.sessions).empty());
}

TEST_CASE("Host reports a host killed mid-response as a dead host, not a malformed response", "[skb-arrow]") {
	// The half line a host leaves when killed while writing is not a
	// response: the user needs how it died, not a JSON parse error.
	FakeHosts fakes;
	const std::string bin = fakes.Write("next\nprintf '%s' '{\"type\":\"res'\nkill -9 $$");
	Host host(fakes.Options(bin), kNever);
	const std::string dir = host.dir();
	try {
		host.Call("echo", Params {}, {{"table", {"t-0"}}}, kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "internal");
		REQUIRE_THAT(e.what(), ContainsSubstring("killed by signal 9"));
	}
	REQUIRE_FALSE(Exists(dir));
}

TEST_CASE("Host passes thread caps and strips Python's environment", "[skb-arrow]") {
	// A Python-hosted DuckDB must not leak its own interpreter's settings into
	// the host, and the host's thread pools follow the caller's budget.
	EnvGuard pythonpath("PYTHONPATH", "/somewhere/else");
	EnvGuard pythonhome("PYTHONHOME", "/another/python");
	FakeHosts fakes;
	const std::string env_log = fakes.tmp.path + "/env.log";
	const std::string bin =
	    fakes.Write("echo \"$NUMBA_NUM_THREADS $OMP_NUM_THREADS $OPENBLAS_NUM_THREADS $VECLIB_MAXIMUM_THREADS "
	                "${PYTHONPATH-unset} ${PYTHONHOME-unset} $PYTHONNOUSERSITE\" > '" +
	                env_log + "'\n" + AnswerOnce("unused"));
	HostOptions options = fakes.Options(bin);
	options.threads = 3;
	Host host(options, kNever);
	host.Shutdown();
	REQUIRE(Slurp(env_log) == "3 3 3 3 unset unset 1\n");
}

// =============================================================================
// Host: calls
// =============================================================================

TEST_CASE("Host sends a call naming its params and input segments", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce(R"({"type":"result","output":["skbout-0"],"id":1,"warnings":[]})"));
	Host host(fakes.Options(bin), kNever);
	Params params;
	params.SetString("method", "spearman");
	params.SetInt("permutations", 99);
	params.SetNumber("alpha", 0.25);
	params.SetBool("flag", true);
	params.SetStrings("posthoc", {"global", "pairwise"});
	const CallResult result = host.Call("mantel", params, {{"x", {"x-0", "x-1"}}, {"y", {"y-0"}}}, kNever);
	REQUIRE(result.output == std::vector<std::string> {"skbout-0"});
	REQUIRE(result.warnings.empty());
	host.Shutdown();
	const auto requests = fakes.Requests();
	REQUIRE(requests.size() == 2);
	REQUIRE(requests[1] == R"({"type":"call","id":1,"capability":"mantel","params":{"method":"spearman",)"
	                       R"("permutations":99,"alpha":0.25,"flag":true,"posthoc":["global","pairwise"]},)"
	                       R"("input":{"x":["x-0","x-1"],"y":["y-0"]}})");
}

TEST_CASE("Host omits params when none are set", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce(R"({"type":"result","output":["skbout-0"],"id":1,"warnings":[]})"));
	Host host(fakes.Options(bin), kNever);
	host.Call("echo", Params {}, {{"table", {"t-0"}}}, kNever);
	host.Shutdown();
	REQUIRE(fakes.Requests().at(1) == R"({"type":"call","id":1,"capability":"echo","input":{"table":["t-0"]}})");
}

TEST_CASE("Host returns a result's warnings", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write(AnswerOnce(
	    R"({"type":"result","output":["skbout-0","skbout-1"],"id":1,"warnings":[{"category":"scipy.stats._warnings_errors.ConstantInputWarning","message":"All values in data are exactly equal; results may not be reliable.","count":2}],"later":"ignored"})"));
	Host host(fakes.Options(bin), kNever);
	const CallResult result = host.Call("mantel", Params {}, {{"x", {"x-0"}}, {"y", {"y-0"}}}, kNever);
	REQUIRE(result.output == std::vector<std::string> {"skbout-0", "skbout-1"});
	REQUIRE(result.warnings.size() == 1);
	REQUIRE(result.warnings[0].category == "scipy.stats._warnings_errors.ConstantInputWarning");
	REQUIRE(result.warnings[0].message == "All values in data are exactly equal; results may not be reliable.");
	REQUIRE(result.warnings[0].count == 2);
	host.Shutdown();
}

TEST_CASE("Host raises an error response with its kind; only internal keeps a traceback", "[skb-arrow]") {
	struct Case {
		std::string response;
		std::string kind;
		std::string message;
		std::string traceback;
	};
	const std::vector<Case> cases = {
	    {R"({"type":"error","kind":"invalid_param","message":"mantel param 'method' must be one of pearson, spearman, kendalltau","id":1,"warnings":[]})",
	     "invalid_param", "mantel param 'method' must be one of pearson, spearman, kendalltau", ""},
	    {R"({"type":"error","kind":"invalid_input","message":"x: 1 of 6 pairs are missing","traceback":"not for this kind","id":1,"warnings":[]})",
	     "invalid_input", "x: 1 of 6 pairs are missing", ""},
	    {R"({"type":"error","kind":"internal","message":"ZeroDivisionError: division by zero","traceback":"Traceback (most recent call last): ...","id":1,"warnings":[]})",
	     "internal", "ZeroDivisionError: division by zero", "Traceback (most recent call last): ..."},
	    // A kind this client doesn't know is internal, keeping what the host said.
	    {R"({"type":"error","kind":"frobnicated","message":"a newer kind","id":1,"warnings":[]})", "internal",
	     "frobnicated", ""},
	};
	for (const auto &c : cases) {
		INFO(c.response);
		FakeHosts fakes;
		const std::string bin = fakes.Write(AnswerOnce(c.response));
		Host host(fakes.Options(bin), kNever);
		try {
			host.Call("mantel", Params {}, {{"x", {"x-0"}}, {"y", {"y-0"}}}, kNever);
			FAIL("expected an Error");
		} catch (const Error &e) {
			REQUIRE(e.kind() == c.kind);
			REQUIRE_THAT(e.what(), ContainsSubstring(c.message));
			REQUIRE(e.traceback() == c.traceback);
		}
		// An error response leaves the host serving; it still shuts down cleanly.
		host.Shutdown();
	}
}

TEST_CASE("Host treats a response to another request as host_incompatible", "[skb-arrow]") {
	for (const std::string response :
	     {std::string(R"({"type":"result","output":["skbout-0"],"id":7,"warnings":[]})"),
	      std::string("this is not json"), std::string(R"({"type":"surprise","id":1,"warnings":[]})")}) {
		INFO(response);
		FakeHosts fakes;
		const std::string bin = fakes.Write(AnswerOnce(response));
		Host host(fakes.Options(bin), kNever);
		try {
			host.Call("echo", Params {}, {{"table", {"t-0"}}}, kNever);
			FAIL("expected an Error");
		} catch (const Error &e) {
			REQUIRE(e.kind() == "host_incompatible");
		}
	}
}

TEST_CASE("Host reports a host that dies mid-call, with its last words, and cleans up", "[skb-arrow]") {
	// The OOM killer, or a crash in compiled code: EOF on stdout mid-call.
	FakeHosts fakes;
	const std::string bin = fakes.Write("next\necho 'last words before dying' >&2\nkill -9 $$");
	Host host(fakes.Options(bin), kNever);
	const std::string dir = host.dir();
	try {
		host.Call("echo", Params {}, {{"table", {"t-0"}}}, kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "internal");
		REQUIRE_THAT(e.what(), ContainsSubstring("killed by signal 9"));
		REQUIRE_THAT(e.what(), ContainsSubstring("last words before dying"));
	}
	REQUIRE_FALSE(Exists(dir));
}

TEST_CASE("Host kills the host and everything it started when the call is interrupted", "[skb-arrow]") {
	// No cancel message exists: cancelling a query kills the process group.
	// The witness pipe reaches EOF only once the host and its `sleep` are gone.
	Witness witness;
	FakeHosts fakes;
	const std::string bin = fakes.Write("next\nsleep 30 &\nwait");
	std::string dir;
	std::chrono::steady_clock::duration took {};
	{
		Host host(fakes.Options(bin), kNever);
		witness.CloseWriteEnd();
		dir = host.dir();
		const auto start = std::chrono::steady_clock::now();
		const InterruptCheck cancel_soon = [start] {
			if (std::chrono::steady_clock::now() - start > std::chrono::milliseconds(200)) {
				throw std::runtime_error("query interrupted");
			}
		};
		REQUIRE_THROWS_WITH(host.Call("echo", Params {}, {{"table", {"t-0"}}}, cancel_soon), "query interrupted");
		took = std::chrono::steady_clock::now() - start;
		// Already gone when the exception arrives, not only at destruction.
		REQUIRE_FALSE(Exists(dir));
		REQUIRE(EofWithin(witness.read_fd, 2000));
	}
	REQUIRE(took < std::chrono::seconds(2));
}

TEST_CASE("Host maps an output segment and unlinks it", "[skb-arrow]") {
	// The test stages a real IPC stream; the fake host moves it into DIR as
	// its output, as the real host creates skbout-<n>.
	FakeHosts fakes;
	OwnedSchema schema;
	OwnedArray batch(schema.schema, 10, 3);
	const auto bytes = gb::EncodeIpcStream(&schema.schema, &batch.array, 1);
	const std::string staged = fakes.tmp.path + "/staged";
	std::ofstream(staged, std::ios::binary).write(reinterpret_cast<const char *>(bytes.data()), bytes.size());
	const std::string bin = fakes.Write("next\nmv '" + staged + "' \"$DIR/skbout-0\"\n" +
	                                    Respond(R"({"type":"result","output":["skbout-0"],"id":1,"warnings":[]})"));
	Host host(fakes.Options(bin), kNever);
	const CallResult result = host.Call("echo", Params {}, {{"table", {"t-0"}}}, kNever);
	MappedSegment segment = host.MapOutput(result.output.at(0));
	REQUIRE_FALSE(Exists(host.dir() + "/skbout-0"));
	REQUIRE(segment.size() == bytes.size());
	const Decoded d = DecodeIntTable(std::string(reinterpret_cast<const char *>(segment.data()), segment.size()));
	REQUIRE(d.a == std::vector<int64_t> {10, 11, 12});
	host.Shutdown();
}

TEST_CASE("Host shutdown lets the host exit, then removes the session directory", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string exited = fakes.tmp.path + "/exited";
	const std::string bin = fakes.Write("while IFS= read -r line; do :; done\ntouch '" + exited + "'\nexit 0");
	Host host(fakes.Options(bin), kNever);
	const std::string dir = host.dir();
	host.Shutdown();
	REQUIRE(Exists(exited)); // it saw EOF and finished on its own
	REQUIRE_FALSE(Exists(dir));
}

TEST_CASE("Host shutdown kills a host that won't exit", "[skb-arrow]") {
	FakeHosts fakes;
	const std::string bin = fakes.Write("exec 0</dev/null\nsleep 30");
	HostOptions options = fakes.Options(bin);
	options.shutdown_grace_ms = 200;
	Host host(options, kNever);
	const std::string dir = host.dir();
	const auto start = std::chrono::steady_clock::now();
	host.Shutdown();
	REQUIRE(std::chrono::steady_clock::now() - start < std::chrono::seconds(2));
	REQUIRE_FALSE(Exists(dir));
}

// =============================================================================
// Checking and installing
// =============================================================================

namespace {
// A ready line offering every capability miint calls.
constexpr const char *kReadyAll = R"({"type":"ready","protocol_version":1,"host_version":"9.9.9",)"
                                  R"("capabilities":{"ancombc":1,"echo":1,"mantel":1},"id":null,"warnings":[]})";

// After the handshake: idle until stdin closes, then exit 0.
constexpr const char *kIdle = "while IFS= read -r line; do :; done\nexit 0";

void WriteExecutable(const std::string &path, const std::string &content) {
	std::ofstream out(path);
	out << content;
	out.close();
	REQUIRE(::chmod(path.c_str(), 0755) == 0);
}

std::vector<std::string> Lines(const std::string &path) {
	std::vector<std::string> lines;
	std::istringstream in(Slurp(path));
	for (std::string line; std::getline(in, line);) {
		lines.push_back(line);
	}
	return lines;
}

// Writes a fake uv: it records its arguments, a line each, then the
// variables that place an install, then runs `body`.
struct FakeUv {
	TempDir tmp;
	std::string log = tmp.path + "/uv.log";
	std::string Write(const std::string &body) {
		const std::string path = tmp.path + "/uv";
		WriteExecutable(path, "#!/bin/bash\n{ printf '%s\\n' \"$@\"; printf '%s\\n' "
		                      "\"UV_TOOL_BIN_DIR=$UV_TOOL_BIN_DIR\" \"UV_TOOL_DIR=$UV_TOOL_DIR\" "
		                      "\"UV_PYTHON_INSTALL_DIR=$UV_PYTHON_INSTALL_DIR\"; } > " +
		                          ShQuote(log) + "\n" + body + "\n");
		return path;
	}
	bool Ran() const {
		return Exists(log);
	}
};

// A fake uv body that installs `host` as the tool's binary.
std::string Installs(const std::string &host) {
	return "mkdir -p \"$UV_TOOL_BIN_DIR\" && cp " + ShQuote(host) + " \"$UV_TOOL_BIN_DIR/skb-arrow\"";
}
} // namespace

TEST_CASE("Probe returns the version of a host offering every capability miint calls", "[skb-arrow]") {
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	const std::string bin = fake.Write(kIdle, kReadyAll);
	REQUIRE(Probe(bin, kNever) == "9.9.9");
	// A real handshake, not a version string: the contract is what's checked.
	REQUIRE(fake.Requests() == std::vector<std::string> {R"({"type":"init","protocol_version":1})"});
	REQUIRE(Entries(fake.sessions).empty());
}

TEST_CASE("Probe refuses a host without a capability miint calls", "[skb-arrow]") {
	// skb_arrow_available() must not vouch for a host that can't serve
	// skb_ancombc, even though it serves skb_mantel.
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	const std::string bin = fake.Write(kIdle); // kReady: echo and mantel only
	try {
		Probe(bin, kNever);
		FAIL("expected an Error");
	} catch (const Error &e) {
		REQUIRE(e.kind() == "host_incompatible");
		REQUIRE_THAT(e.what(), ContainsSubstring("'ancombc'"));
	}
	REQUIRE(Entries(fake.sessions).empty());
}

TEST_CASE("FindUv looks on PATH, then where uv's installer puts it", "[skb-arrow]") {
	TempDir tmp;
	const std::string which = gb::FindExecutableInPath("which");
	REQUIRE_FALSE(which.empty());
	// PATH holds only `which`, which FindExecutableInPath runs, so no uv on
	// this machine can leak in.
	const std::string path_dir = tmp.path + "/path";
	const std::string home = tmp.path + "/home";
	REQUIRE(::mkdir(path_dir.c_str(), 0700) == 0);
	REQUIRE(::symlink(which.c_str(), (path_dir + "/which").c_str()) == 0);
	REQUIRE(::system(("mkdir -p " + ShQuote(home + "/.local/bin")).c_str()) == 0);
	EnvGuard path("PATH", path_dir.c_str());
	EnvGuard home_env("HOME", home.c_str());

	REQUIRE(FindUv().empty());
	WriteExecutable(home + "/.local/bin/uv", "#!/bin/sh\n");
	REQUIRE(FindUv() == home + "/.local/bin/uv");
	WriteExecutable(path_dir + "/uv", "#!/bin/sh\n");
	REQUIRE(FindUv() == path_dir + "/uv");
}

TEST_CASE("Install returns at once when the existing binary passes the probe", "[skb-arrow]") {
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.existing = fake.Write(kIdle, kReadyAll);
	o.uv = uv.Write("exit 0");
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE(r.installed);
	REQUIRE(r.path == o.existing);
	REQUIRE(r.version == "9.9.9");
	REQUIRE_THAT(r.message, ContainsSubstring("already available"));
	REQUIRE_FALSE(uv.Ran());
}

TEST_CASE("Install runs uv on the pinned release, keeping the tool, its Python and its binary in miint's cache",
          "[skb-arrow]") {
	// --managed-python keeps conda's or the system's interpreter out;
	// UV_PYTHON_INSTALL_DIR keeps uv's own out of ~/.local/share/uv; the pin
	// is the release miint is tested with. uv runs only when there is no
	// working skb-arrow, so it always reinstalls: see the next tests.
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.uv = uv.Write(Installs(fake.Write(kIdle, kReadyAll)));
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE(Lines(uv.log) == std::vector<std::string> {"tool", "install", "--managed-python", "--python", "3.14",
	                                                   "skb-arrow==0.1.0", "--reinstall", "--force",
	                                                   "UV_TOOL_BIN_DIR=" + cache.path + "/bin",
	                                                   "UV_TOOL_DIR=" + cache.path + "/uv-tools",
	                                                   "UV_PYTHON_INSTALL_DIR=" + cache.path + "/uv-python"});
	REQUIRE(r.installed);
	REQUIRE(r.path == cache.path + "/bin/skb-arrow");
	REQUIRE(r.version == "9.9.9");
	REQUIRE_THAT(r.message, ContainsSubstring("Installed skb-arrow 9.9.9 to " + r.path));
}

TEST_CASE("Install installs over an existing binary that fails the probe", "[skb-arrow]") {
	FakeHosts stale, fresh;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fresh.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.existing = stale.Write(kIdle); // no ancombc
	o.uv = uv.Write(Installs(fresh.Write(kIdle, kReadyAll)));
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE(uv.Ran());
	REQUIRE(r.installed);
	REQUIRE(r.path == cache.path + "/bin/skb-arrow");
}

TEST_CASE("Install reinstalls a tool uv still has after its binary is gone", "[skb-arrow]") {
	// uv skips a tool it already has, exiting 0 without restoring a deleted
	// or broken binary ("`skb-arrow==0.1.0` is already installed"); only
	// --reinstall replaces it. This fake uv behaves that way.
	FakeHosts fresh;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fresh.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o; // no existing binary: FindSkbArrow found none
	o.uv =
	    uv.Write("if printf '%s\\n' \"$@\" | grep -qx -- --reinstall; then " + Installs(fresh.Write(kIdle, kReadyAll)) +
	             "; else echo '`skb-arrow==0.1.0` is already installed' >&2; fi");
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE(r.installed);
	REQUIRE(r.path == cache.path + "/bin/skb-arrow");
}

TEST_CASE("Install with force skips the probe and has uv reinstall and replace the binary", "[skb-arrow]") {
	// force is the escape from a binary that passes the probe but is broken
	// some other way, or one uv didn't put there (`--force` overwrites it).
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.existing = fake.Write(kIdle, kReadyAll);
	o.uv = uv.Write(Installs(o.existing));
	o.cache_dir = cache.path;
	o.force = true;
	const InstallReport r = Install(o, kNever);
	const std::vector<std::string> args = Lines(uv.log);
	REQUIRE(args.size() == 11);
	REQUIRE(std::vector<std::string>(args.begin() + 6, args.begin() + 8) ==
	        std::vector<std::string> {"--reinstall", "--force"});
	REQUIRE(r.installed);
	REQUIRE(r.path == cache.path + "/bin/skb-arrow");
}

TEST_CASE("Install without uv says how to get it", "[skb-arrow]") {
	TempDir cache;
	InstallOptions o;
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE_FALSE(r.installed);
	REQUIRE_THAT(r.message, ContainsSubstring("uv not found"));
	REQUIRE_THAT(r.message, ContainsSubstring("https://astral.sh/uv/install.sh"));
	REQUIRE_THAT(r.message, ContainsSubstring("SELECT install_skb_arrow()"));
}

TEST_CASE("Install without a cache location says so", "[skb-arrow]") {
	FakeUv uv;
	InstallOptions o;
	o.uv = uv.Write("exit 0");
	const InstallReport r = Install(o, kNever);
	REQUIRE_FALSE(r.installed);
	REQUIRE_THAT(r.message, ContainsSubstring("XDG_CACHE_HOME"));
	REQUIRE_FALSE(uv.Ran());
}

TEST_CASE("Install reports uv's failure with what uv said", "[skb-arrow]") {
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.uv = uv.Write("echo 'error: Failed to fetch: no network' >&2\nexit 2");
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE_FALSE(r.installed);
	REQUIRE_THAT(r.message, ContainsSubstring("exit code 2"));
	REQUIRE_THAT(r.message, ContainsSubstring("error: Failed to fetch: no network"));
}

TEST_CASE("Install reports a uv success that left no binary", "[skb-arrow]") {
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.uv = uv.Write("exit 0");
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE_FALSE(r.installed);
	REQUIRE_THAT(r.message, ContainsSubstring("no skb-arrow at " + cache.path + "/bin/skb-arrow"));
}

TEST_CASE("Install reports an installed binary that fails the probe", "[skb-arrow]") {
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	FakeUv uv;
	TempDir cache;
	InstallOptions o;
	o.uv = uv.Write(Installs(fake.Write(kIdle))); // no ancombc
	o.cache_dir = cache.path;
	const InstallReport r = Install(o, kNever);
	REQUIRE_FALSE(r.installed);
	REQUIRE_THAT(r.message, ContainsSubstring("'ancombc'"));
}

TEST_CASE("Install puts the binary where FindSkbArrow looks", "[skb-arrow]") {
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	FakeUv uv;
	TempDir xdg;
	EnvGuard xdg_env("XDG_CACHE_HOME", xdg.path.c_str());
	EnvGuard override_path("MIINT_SKB_ARROW_PATH", nullptr);
	REQUIRE(MiintCacheDir() == xdg.path + "/miint");
	InstallOptions o;
	o.uv = uv.Write(Installs(fake.Write(kIdle, kReadyAll)));
	o.cache_dir = MiintCacheDir();
	const InstallReport r = Install(o, kNever);
	REQUIRE(r.installed);
	REQUIRE(r.path == MiintSkbArrowCacheBinary());
	REQUIRE(FindSkbArrow() == r.path);
}

TEST_CASE("Install lets an interrupt through rather than reporting it", "[skb-arrow]") {
	// A swallowed interrupt would leave Ctrl-C waiting on a hung binary.
	FakeHosts fake;
	EnvGuard segdir("MIINT_SKB_ARROW_SEGMENT_DIR", fake.sessions.c_str());
	TempDir tmp;
	InstallOptions o;
	o.existing = tmp.path + "/silent-skb-arrow";
	WriteExecutable(o.existing, "#!/bin/bash\nexec sleep 30\n");
	o.cache_dir = tmp.path;
	const auto interrupt = [] {
		throw std::runtime_error("interrupted");
	};
	const auto start = std::chrono::steady_clock::now();
	REQUIRE_THROWS_WITH(Install(o, interrupt), "interrupted");
	REQUIRE(std::chrono::steady_clock::now() - start < std::chrono::seconds(5));
	REQUIRE(Entries(fake.sessions).empty());
}

// =============================================================================
// Against the real host (skipped without one)
// =============================================================================

namespace {
// struct<id_a utf8, id_b utf8, distance double> over the pairs of a..d.
struct DistanceTable {
	ArrowSchema schema {};
	ArrowArray array {};
	explicit DistanceTable(const std::vector<double> &distances) {
		REQUIRE(ArrowSchemaInitFromType(&schema, NANOARROW_TYPE_STRUCT) == NANOARROW_OK);
		REQUIRE(ArrowSchemaAllocateChildren(&schema, 3) == NANOARROW_OK);
		const char *names[] = {"id_a", "id_b", "distance"};
		for (int i = 0; i < 3; ++i) {
			REQUIRE(ArrowSchemaInitFromType(schema.children[i],
			                                i < 2 ? NANOARROW_TYPE_STRING : NANOARROW_TYPE_DOUBLE) == NANOARROW_OK);
			REQUIRE(ArrowSchemaSetName(schema.children[i], names[i]) == NANOARROW_OK);
		}
		ArrowError err {};
		REQUIRE(ArrowArrayInitFromSchema(&array, &schema, &err) == NANOARROW_OK);
		REQUIRE(ArrowArrayStartAppending(&array) == NANOARROW_OK);
		const char *pairs[][2] = {{"a", "b"}, {"a", "c"}, {"a", "d"}, {"b", "c"}, {"b", "d"}, {"c", "d"}};
		for (std::size_t i = 0; i < distances.size(); ++i) {
			REQUIRE(ArrowArrayAppendString(array.children[0], {pairs[i][0], 1}) == NANOARROW_OK);
			REQUIRE(ArrowArrayAppendString(array.children[1], {pairs[i][1], 1}) == NANOARROW_OK);
			REQUIRE(ArrowArrayAppendDouble(array.children[2], distances[i]) == NANOARROW_OK);
			REQUIRE(ArrowArrayFinishElement(&array) == NANOARROW_OK);
		}
		REQUIRE(ArrowArrayFinishBuildingDefault(&array, &err) == NANOARROW_OK);
	}
	~DistanceTable() {
		array.release(&array);
		schema.release(&schema);
	}
};

std::vector<std::string> WriteTable(const Host &host, const std::string &prefix, const ArrowSchema &schema,
                                    const std::vector<const ArrowArray *> &batches, std::size_t cap = kSegmentBytes) {
	SegmentWriter writer(host.dir(), prefix, schema, cap);
	for (const ArrowArray *batch : batches) {
		writer.Append(*batch);
	}
	return writer.Finish();
}

HostOptions RealHostOptions(const std::string &binary) {
	HostOptions o;
	o.binary = binary;
	o.segment_parent = SegmentDirParent();
	o.threads = 2;
	o.required = {{"echo", 1}, {"mantel", 1}};
	return o;
}
} // namespace

TEST_CASE("The real host echoes a two-segment table", "[skb-arrow]") {
	const std::string bin = FindSkbArrow();
	if (bin.empty() || !Exists(bin)) {
		SKIP("skb-arrow not found (SELECT install_skb_arrow(), or set MIINT_SKB_ARROW_PATH)");
	}
	Host host(RealHostOptions(bin), kNever);
	OwnedSchema schema;
	OwnedArray b0(schema.schema, 0, 3), b1(schema.schema, 3, 2);
	const auto segments = WriteTable(host, "table", schema.schema, {&b0.array, &b1.array}, 1);
	REQUIRE(segments.size() == 2);
	const CallResult result = host.Call("echo", Params {}, {{"table", segments}}, kNever);
	REQUIRE(Entries(host.dir()).size() == result.output.size()); // inputs consumed, outputs waiting
	std::vector<int64_t> all;
	for (const auto &name : result.output) {
		MappedSegment segment = host.MapOutput(name);
		const Decoded d = DecodeIntTable(std::string(reinterpret_cast<const char *>(segment.data()), segment.size()));
		REQUIRE(d.names == std::vector<std::string> {"a", "s"});
		all.insert(all.end(), d.a.begin(), d.a.end());
	}
	REQUIRE(all == std::vector<int64_t> {0, 1, 2, 3, 4});
	const std::string dir = host.dir();
	host.Shutdown();
	REQUIRE_FALSE(Exists(dir));
}

TEST_CASE("The real host's mantel matches the handoff's known answer exactly", "[skb-arrow]") {
	// skb-arrow's handoff transcript, call 2: recorded against the published
	// 0.1.0. Mantel's statistic and permutation p-value don't depend on thread
	// counts, so they must match bit for bit.
	const std::string bin = FindSkbArrow();
	if (bin.empty() || !Exists(bin)) {
		SKIP("skb-arrow not found (SELECT install_skb_arrow(), or set MIINT_SKB_ARROW_PATH)");
	}
	Host host(RealHostOptions(bin), kNever);
	DistanceTable x({1, 2, 3, 4, 5, 6});
	DistanceTable y({1.5, 2, 3.5, 4, 6, 6.5});
	const auto xs = WriteTable(host, "x", x.schema, {&x.array});
	const auto ys = WriteTable(host, "y", y.schema, {&y.array});
	const CallResult result = host.Call("mantel", Params {}, {{"x", xs}, {"y", ys}}, kNever);
	REQUIRE(result.output.size() == 1);
	REQUIRE(result.warnings.empty());
	MappedSegment segment = host.MapOutput(result.output[0]);
	gb::IpcStreamDecoder decoder(segment.data(), segment.size());
	ArrowSchema schema {};
	decoder.GetSchema(&schema);
	REQUIRE(schema.n_children == 3);
	REQUIRE(std::string(schema.children[0]->name) == "statistic");
	REQUIRE(std::string(schema.children[1]->name) == "pvalue");
	REQUIRE(std::string(schema.children[2]->name) == "n");
	ArrowArray batch {};
	REQUIRE(decoder.NextBatch(&batch));
	REQUIRE(batch.length == 1);
	const double statistic = static_cast<const double *>(batch.children[0]->buffers[1])[batch.children[0]->offset];
	const double pvalue = static_cast<const double *>(batch.children[1]->buffers[1])[batch.children[1]->offset];
	const int64_t n = static_cast<const int64_t *>(batch.children[2]->buffers[1])[batch.children[2]->offset];
	batch.release(&batch);
	schema.release(&schema);
	REQUIRE(statistic == 0.9849397709505465);
	REQUIRE(pvalue == 0.036);
	REQUIRE(n == 4);
	host.Shutdown();
}

TEST_CASE("The real host passes the probe at the release miint pins", "[skb-arrow]") {
	// The pinned release must offer every capability miint calls, or
	// install_skb_arrow() would install a host skb_arrow_available() rejects.
	const std::string bin = FindSkbArrow();
	if (bin.empty() || !Exists(bin)) {
		SKIP("skb-arrow not found (SELECT install_skb_arrow(), or set MIINT_SKB_ARROW_PATH)");
	}
	REQUIRE("skb-arrow==" + Probe(bin, kNever) == kPinnedRelease);
}
