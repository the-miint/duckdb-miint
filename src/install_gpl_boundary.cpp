#include "install_gpl_boundary.hpp"

#include "gpl_boundary/process.hpp"
#include "yyjson.hpp"

#include <array>
#include <cerrno>
#include <cstring>
#include <mutex>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>
#include <vector>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

constexpr const char *kInstallScriptUrl =
    "https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh";

// --max-time on OUR fetch of install.sh. Once we hand off to /bin/sh running
// install.sh, we're at the mercy of upstream's internal curl call for the
// tarball (which has no timeout we can set from outside). The user-facing
// docs at docs/utilities.md call this out as a known limitation.
constexpr int kInstallScriptFetchTimeoutSec = 60;

// =============================================================================
// Result struct
// =============================================================================
struct InstallReport {
	bool installed = false;
	std::string path;
	std::string version;
	std::string message;
};

// =============================================================================
// Helpers
// =============================================================================

// Trim trailing whitespace/newlines from a captured command output.
std::string strip_trailing_ws(std::string s) {
	while (!s.empty() && (s.back() == '\n' || s.back() == '\r' || s.back() == ' ' || s.back() == '\t')) {
		s.pop_back();
	}
	return s;
}

// POSIX shell single-quote escaping. For any byte sequence `s`, the result is
// inert input to /bin/sh: it expands to `s` verbatim with no metacharacter
// interpretation. Each `'` in the input is encoded as `'\''` (close-quote,
// escaped-quote, re-open-quote).
//
// This is the right answer for splicing user-controlled paths (from $HOME /
// $XDG_CACHE_HOME / $MIINT_GPL_BOUNDARY_PATH) into a shell -c string.
// Without it, a value like `/tmp/x" ; rm -rf $HOME ; echo "` would execute.
std::string sh_quote(const std::string &s) {
	std::string out;
	out.reserve(s.size() + 2);
	out += '\'';
	for (char c : s) {
		if (c == '\'') {
			out += "'\\''";
		} else {
			out += c;
		}
	}
	out += '\'';
	return out;
}

// Drain everything from `fd` until EOF. Used to pull stdout / stderr off a
// `gb::ChildProcess`. Robust to short reads and EINTR.
std::string drain_fd(int fd) {
	std::string out;
	std::array<char, 4096> buf {};
	for (;;) {
		ssize_t n = ::read(fd, buf.data(), buf.size());
		if (n > 0) {
			out.append(buf.data(), static_cast<size_t>(n));
			continue;
		}
		if (n == 0) {
			break;
		}
		if (errno == EINTR) {
			continue;
		}
		break; // any other error: bail with whatever we got
	}
	return out;
}

// Result of running a child process with separate stdout/stderr capture.
struct ProcResult {
	bool spawned = false;         // true iff fork+exec succeeded
	bool exited_normally = false; // true iff child returned via exit()
	int exit_code = -1;           // valid iff exited_normally
	int term_signal = 0;          // non-zero iff WIFSIGNALED — child was killed
	std::string stdout_text;
	std::string stderr_text;
	std::string spawn_error; // populated iff !spawned
};

// Fork+exec the binary at argv[0] with argv[1..], capture stdout and stderr
// separately, wait, and return a structured result. No shell — no quoting
// required. This is the injection-safe version-probe path.
ProcResult run_argv(const std::vector<std::string> &argv) {
	ProcResult r;
	try {
		gb::ChildProcess child(argv);
		// Close stdin so the child immediately sees EOF and won't block
		// waiting for input on a one-shot probe.
		child.CloseStdin();
		r.stdout_text = drain_fd(child.stdout_fd());
		r.stderr_text = drain_fd(child.stderr_fd());
		const int status = child.Wait();
		r.spawned = true;
		if (WIFEXITED(status)) {
			r.exited_normally = true;
			r.exit_code = WEXITSTATUS(status);
		} else if (WIFSIGNALED(status)) {
			r.term_signal = WTERMSIG(status);
		}
	} catch (const std::exception &e) {
		r.spawn_error = e.what();
	}
	return r;
}

// gpl-boundary --version emits a JSON document like
//   { "gpl_boundary": "0.1.0", "tools": [ ... ] }
// Pull just the top-level "gpl_boundary" string. Uses yyjson (already linked
// via duckdb's bundled third_party) so substring-confusion ("non_gpl_boundary"
// matching as a suffix) and escape-handling edge cases are handled correctly.
// Falls back to "(unparsed)" if the JSON shape surprises us — the full payload
// is still exposed via `version` so callers can json_extract it themselves.
std::string extract_short_version(const std::string &full) {
	const auto fallback = "(unparsed)";
	auto *doc = duckdb_yyjson::yyjson_read(full.data(), full.size(), 0);
	if (!doc) {
		return fallback;
	}
	auto *root = duckdb_yyjson::yyjson_doc_get_root(doc);
	if (!root || !duckdb_yyjson::yyjson_is_obj(root)) {
		duckdb_yyjson::yyjson_doc_free(doc);
		return fallback;
	}
	auto *v = duckdb_yyjson::yyjson_obj_get(root, "gpl_boundary");
	if (!v || !duckdb_yyjson::yyjson_is_str(v)) {
		duckdb_yyjson::yyjson_doc_free(doc);
		return fallback;
	}
	std::string out(duckdb_yyjson::yyjson_get_str(v), duckdb_yyjson::yyjson_get_len(v));
	duckdb_yyjson::yyjson_doc_free(doc);
	return out;
}

// =============================================================================
// Probe + install
// =============================================================================

// Probe whatever FindGplBoundary returns (override path, cache dir, or PATH);
// if a working binary is there, populate `r` and return true. A "working
// binary" responds to `--version` with exit 0. Mirrors the stricter gating
// in `phylogeny_fasttree_available()`.
//
// Uses gb::ChildProcess directly (no shell) so the path can contain any byte
// sequence including spaces, quotes, or shell metacharacters without becoming
// an injection vector.
bool probe_existing(InstallReport &r) {
	const std::string p = gb::FindGplBoundary();
	if (p.empty() || ::access(p.c_str(), X_OK) != 0) {
		return false;
	}
	const ProcResult result = run_argv({p, "--version"});
	if (!result.spawned || !result.exited_normally || result.exit_code != 0) {
		return false; // present but unusable; treat as missing so install can fix it
	}
	r.installed = true;
	r.path = p;
	r.version = strip_trailing_ws(result.stdout_text);
	r.message =
	    "gpl-boundary " + extract_short_version(r.version) + " already available at " + p + "; no install performed";
	return true;
}

// Process-local mutex over `install_impl`. Note: this serializes within ONE
// DuckDB process. Two separate DuckDB processes calling install_gpl_boundary()
// concurrently both run install.sh against the same cache dir; the upstream
// installer's final `mv` is atomic on the same filesystem so the on-disk
// binary is one of the two valid downloads (both are bit-identical for the
// same `latest` release), but the tmpdirs and download bandwidth are wasted.
// Multi-process coordination would need a file lock on the cache dir;
// deferred until anyone reports a real problem.
std::mutex &install_mutex() {
	static std::mutex m;
	return m;
}

InstallReport install_impl(bool force) {
	InstallReport r;

	// Double-checked-locking: probe once outside the lock to short-circuit
	// the common case where gpl-boundary is already on PATH. Concurrent
	// callers may all run the version probe in parallel, but the probe is
	// idempotent and cheap. We re-probe inside the lock to close the race
	// where an in-flight install completed between our outer check and lock
	// acquisition.
	//
	// `force=true` skips the probe and goes straight to download. Useful
	// when a stale binary on PATH (e.g. an old ~/.cargo/bin/gpl-boundary
	// from `cargo install`) is shadowing the latest release.
	if (!force && probe_existing(r)) {
		return r;
	}

	std::lock_guard<std::mutex> lock(install_mutex());

	// 1. Re-probe inside the lock; another thread may have just installed.
	if (!force && probe_existing(r)) {
		return r;
	}

	// 2. Resolve install location.
	const std::string cache_dir = gb::MiintGplBoundaryCacheDir();
	const std::string cache_bin = gb::MiintGplBoundaryCacheBinary();
	if (cache_dir.empty() || cache_bin.empty()) {
		r.message = "install_gpl_boundary: cannot resolve install location "
		            "($XDG_CACHE_HOME and $HOME are both unset)";
		return r;
	}

	// 3. curl is required by upstream install.sh. Check it's on PATH so we
	//    can give a clear error before fork+exec spits out something cryptic.
	if (gb::FindExecutableInPath("curl").empty()) {
		r.message = "install_gpl_boundary: 'curl' not found on PATH. install.sh requires curl; "
		            "either install curl or fetch the gpl-boundary release tarball manually from "
		            "https://github.com/the-miint/GPL-boundary/releases";
		return r;
	}

	// 4. Make a tmpdir to hold install.sh.
	char tmpl[] = "/tmp/miint_install_gpl_XXXXXX";
	if (::mkdtemp(tmpl) == nullptr) {
		r.message = "install_gpl_boundary: mkdtemp failed: ";
		r.message += std::strerror(errno);
		return r;
	}
	const std::string tmpdir = tmpl;

	// 5. Run the installer. Path values are sh_quote'd inline rather than
	//    setenv'd, so:
	//      - we don't pollute the process's environment for unrelated code
	//      - we don't race with concurrent getenv() calls in other threads
	//        (POSIX setenv/getenv are not thread-safe wrt each other)
	//      - any byte in $HOME / $XDG_CACHE_HOME / mkdtemp output is inert
	//        to /bin/sh, including quotes, spaces, and shell metacharacters
	//
	//    Upstream install.sh: detects platform, downloads the tarball +
	//    SHA256SUMS, verifies SHA256, extracts to $INSTALL_DIR/gpl-boundary.
	//    Our role is fetching install.sh itself + setting INSTALL_DIR.
	const std::string install_sh = tmpdir + "/install.sh";
	const std::string sh_cmd = std::string("set -e; ") + "mkdir -p " + sh_quote(cache_dir) + " && " +
	                           "curl -fsSL --proto '=https' --tlsv1.2 --max-time " +
	                           std::to_string(kInstallScriptFetchTimeoutSec) + " " + sh_quote(kInstallScriptUrl) +
	                           " -o " + sh_quote(install_sh) + " && " + "INSTALL_DIR=" + sh_quote(cache_dir) + " sh " +
	                           sh_quote(install_sh);

	const ProcResult install = run_argv({"/bin/sh", "-c", sh_cmd});

	// Best-effort cleanup of the tmpdir.
	(void)run_argv({"/bin/rm", "-rf", tmpdir});

	if (!install.spawned) {
		r.message = "install_gpl_boundary: failed to spawn /bin/sh: " + install.spawn_error;
		return r;
	}
	if (!install.exited_normally) {
		r.message = "install_gpl_boundary: installer killed by signal " + std::to_string(install.term_signal) +
		            "\nstdout:\n" + install.stdout_text + "\nstderr:\n" + install.stderr_text;
		return r;
	}
	if (install.exit_code != 0) {
		r.message = "install_gpl_boundary: installer failed (exit=" + std::to_string(install.exit_code) +
		            "). install.sh prints platform-detection / SHA-mismatch diagnostics on stderr:\n" + "stdout:\n" +
		            install.stdout_text + "\nstderr:\n" + install.stderr_text;
		return r;
	}

	// 6. Verify the binary now exists where we expected.
	if (::access(cache_bin.c_str(), X_OK) != 0) {
		r.message = "install_gpl_boundary: installer reported success but no executable at " + cache_bin +
		            "\nstdout:\n" + install.stdout_text + "\nstderr:\n" + install.stderr_text;
		return r;
	}

	// 7. Probe version on the fresh binary (argv-based — no shell, no
	//    injection risk regardless of what's in cache_bin).
	const ProcResult ver = run_argv({cache_bin, "--version"});

	r.installed = true;
	r.path = cache_bin;
	r.version =
	    (ver.spawned && ver.exited_normally && ver.exit_code == 0) ? strip_trailing_ws(ver.stdout_text) : "(unknown)";
	r.message = "Installed gpl-boundary " + extract_short_version(r.version) + " to " + cache_bin;
	return r;
}

// =============================================================================
// DuckDB scalar plumbing
// =============================================================================

LogicalType InstallReturnType() {
	return LogicalType::STRUCT({{"installed", LogicalType::BOOLEAN},
	                            {"path", LogicalType::VARCHAR},
	                            {"version", LogicalType::VARCHAR},
	                            {"message", LogicalType::VARCHAR}});
}

// Resolve the `force` argument from the call site. The 0-arg overload defaults
// to false; the 1-arg overload reads a BOOLEAN. A NULL input is treated as
// false rather than throwing — same convention as `mask_dust(seq, NULL)`.
bool ResolveForceArg(DataChunk &args) {
	if (args.ColumnCount() == 0) {
		return false;
	}
	UnifiedVectorFormat fmt;
	args.data[0].ToUnifiedFormat(args.size(), fmt);
	const auto idx = fmt.sel->get_index(0);
	if (!fmt.validity.RowIsValid(idx)) {
		return false;
	}
	return UnifiedVectorFormat::GetData<bool>(fmt)[idx];
}

void InstallGplBoundaryExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	(void)state;
	// Run the install once and emit a constant vector; even if the query
	// surface is `SELECT install_gpl_boundary() FROM range(N)`, we don't
	// re-probe N times. We only have to write index 0 because subsequent
	// rows are aliased through CONSTANT_VECTOR semantics.
	const bool force = ResolveForceArg(args);
	InstallReport report = install_impl(force);

	auto &entries = StructVector::GetEntries(result);
	FlatVector::GetData<bool>(*entries[0])[0] = report.installed;
	FlatVector::GetData<string_t>(*entries[1])[0] = StringVector::AddString(*entries[1], report.path);
	FlatVector::GetData<string_t>(*entries[2])[0] = StringVector::AddString(*entries[2], report.version);
	FlatVector::GetData<string_t>(*entries[3])[0] = StringVector::AddString(*entries[3], report.message);
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
}

} // namespace

ScalarFunction InstallGplBoundaryScalar::GetFunction() {
	return ScalarFunction("install_gpl_boundary", {}, InstallReturnType(), InstallGplBoundaryExecute);
}

void InstallGplBoundaryScalar::Register(ExtensionLoader &loader) {
	ScalarFunctionSet set("install_gpl_boundary");
	// 0-arg: idempotent find-or-install. Short-circuits on any binary found
	// at the override path, miint cache, or PATH.
	set.AddFunction(ScalarFunction({}, InstallReturnType(), InstallGplBoundaryExecute));
	// 1-arg: force=true bypasses the probe and re-downloads latest into the
	// miint cache. Used to escape a stale binary on PATH that the probe
	// considers "good enough".
	set.AddFunction(ScalarFunction({LogicalType::BOOLEAN}, InstallReturnType(), InstallGplBoundaryExecute));
	loader.RegisterFunction(set);
}

} // namespace duckdb
