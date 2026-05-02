#include "install_gpl_boundary.hpp"

#include "gpl_boundary/process.hpp"

#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <mutex>
#include <string>
#include <sys/stat.h>
#include <sys/wait.h>
#include <unistd.h>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

constexpr const char *kInstallScriptUrl =
    "https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh";

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

// gpl-boundary --version emits a JSON document like
//   { "gpl_boundary": "0.1.0", "tools": [ ... ] }
// Pull just the top-level "gpl_boundary" string for compact display in
// `version` / `message`. Falls back to "(unparsed)" if the JSON shape
// surprises us — we still expose the full payload via the `version` field
// so callers can parse it with json_extract themselves.
std::string extract_short_version(const std::string &full) {
	const std::string key = "\"gpl_boundary\"";
	auto pos = full.find(key);
	if (pos == std::string::npos) {
		return "(unparsed)";
	}
	pos = full.find(':', pos + key.size());
	if (pos == std::string::npos) {
		return "(unparsed)";
	}
	auto open = full.find('"', pos);
	if (open == std::string::npos) {
		return "(unparsed)";
	}
	auto close = full.find('"', open + 1);
	if (close == std::string::npos) {
		return "(unparsed)";
	}
	return full.substr(open + 1, close - open - 1);
}

// Run `cmd` via /bin/sh -c, capturing combined stdout+stderr. Returns the
// captured output; writes the exit status (or -1 on spawn failure) into
// `*exit_status`. Used for shell-driven steps that don't need bidirectional
// I/O — `gb::ChildProcess` is the right tool for long-lived daemons but
// overkill for a one-shot installer invocation.
std::string popen_capture(const std::string &cmd, int *exit_status) {
	FILE *p = ::popen(cmd.c_str(), "r");
	if (!p) {
		*exit_status = -1;
		return {};
	}
	std::string out;
	char buf[256];
	for (;;) {
		std::size_t n = ::fread(buf, 1, sizeof(buf), p);
		if (n == 0) {
			break;
		}
		out.append(buf, n);
	}
	int rc = ::pclose(p);
	*exit_status = WIFEXITED(rc) ? WEXITSTATUS(rc) : -1;
	return out;
}

// Probe whatever FindGplBoundary returns (override path, cache dir, or PATH);
// if a working binary is there, populate `r` and return true. A "working
// binary" is one that responds to `--version` with exit 0. Mirrors the
// stricter gating in `phylogeny_fasttree_available()`.
bool probe_existing(InstallReport &r) {
	const std::string p = gb::FindGplBoundary();
	if (p.empty() || ::access(p.c_str(), X_OK) != 0) {
		return false;
	}
	int rc = -1;
	const std::string cmd = "\"" + p + "\" --version 2>/dev/null";
	std::string out = popen_capture(cmd, &rc);
	if (rc != 0) {
		return false; // present but unusable; treat as missing so install can fix it
	}
	r.installed = true;
	r.path = p;
	r.version = strip_trailing_ws(out);
	r.message =
	    "gpl-boundary " + extract_short_version(r.version) + " already available at " + p + "; no install performed";
	return true;
}

// =============================================================================
// Install logic
// =============================================================================
//
// Mutex serializes concurrent installs (two parallel queries calling
// install_gpl_boundary() can't both write to the same cache dir).
std::mutex &install_mutex() {
	static std::mutex m;
	return m;
}

InstallReport install_impl() {
	std::lock_guard<std::mutex> lock(install_mutex());

	InstallReport r;

	// 1. Already installed? Short-circuit without touching the network.
	if (probe_existing(r)) {
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

	// 5. Run the installer. We pass the cache dir + tmpdir via env vars
	//    rather than splicing them into the shell command, so paths
	//    containing quotes/spaces don't break the shell parsing. The
	//    upstream install.sh:
	//      - downloads the platform-specific tarball,
	//      - downloads SHA256SUMS,
	//      - verifies the SHA256,
	//      - extracts to $INSTALL_DIR/gpl-boundary.
	//    Our sole responsibility is to fetch install.sh itself and run it
	//    with INSTALL_DIR pointed at miint's cache.
	::setenv("MIINT_INSTALL_TMP", tmpdir.c_str(), 1);
	::setenv("MIINT_INSTALL_DIR_INTERNAL", cache_dir.c_str(), 1);

	const std::string cmd = std::string("set -e; ") + "mkdir -p \"$MIINT_INSTALL_DIR_INTERNAL\" && " +
	                        "curl -fsSL --proto '=https' --tlsv1.2 --max-time 60 " + "\"" + kInstallScriptUrl +
	                        "\" -o \"$MIINT_INSTALL_TMP/install.sh\" && " +
	                        "INSTALL_DIR=\"$MIINT_INSTALL_DIR_INTERNAL\" sh \"$MIINT_INSTALL_TMP/install.sh\" 2>&1";

	int install_rc = -1;
	std::string install_out = popen_capture(cmd, &install_rc);

	// Best-effort cleanup of the tmpdir; we don't care if it fails.
	{
		const std::string clean = "rm -rf \"$MIINT_INSTALL_TMP\"";
		(void)::system(clean.c_str());
	}

	if (install_rc != 0) {
		r.message = "install_gpl_boundary: installer failed (exit=" + std::to_string(install_rc) +
		            "). The upstream install.sh prints a platform-detection or SHA-mismatch "
		            "diagnostic on stderr; captured output below:\n" +
		            install_out;
		return r;
	}

	// 6. Verify the binary now exists where we expected.
	if (::access(cache_bin.c_str(), X_OK) != 0) {
		r.message = "install_gpl_boundary: installer reported success but no executable at " + cache_bin +
		            ". Captured output:\n" + install_out;
		return r;
	}

	// 7. Probe version on the fresh binary.
	int v_rc = -1;
	const std::string v_cmd = "\"" + cache_bin + "\" --version 2>/dev/null";
	std::string v_out = popen_capture(v_cmd, &v_rc);

	r.installed = true;
	r.path = cache_bin;
	r.version = (v_rc == 0) ? strip_trailing_ws(v_out) : "(unknown)";
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

void InstallGplBoundaryExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	(void)args;
	(void)state;
	// Run the install once for this Execute call; replicate the result across
	// all output rows. Calling install_gpl_boundary() inside a multi-row query
	// (e.g., `SELECT install_gpl_boundary() FROM range(N)`) would otherwise
	// re-probe (and possibly re-install) N times.
	InstallReport report = install_impl();

	auto &entries = StructVector::GetEntries(result);
	auto installed_data = FlatVector::GetData<bool>(*entries[0]);
	auto &path_vec = *entries[1];
	auto &version_vec = *entries[2];
	auto &message_vec = *entries[3];

	const idx_t n = args.size();
	for (idx_t i = 0; i < n; i++) {
		installed_data[i] = report.installed;
		FlatVector::GetData<string_t>(path_vec)[i] = StringVector::AddString(path_vec, report.path);
		FlatVector::GetData<string_t>(version_vec)[i] = StringVector::AddString(version_vec, report.version);
		FlatVector::GetData<string_t>(message_vec)[i] = StringVector::AddString(message_vec, report.message);
	}
	result.SetVectorType(VectorType::CONSTANT_VECTOR);
}

} // namespace

ScalarFunction InstallGplBoundaryScalar::GetFunction() {
	return ScalarFunction("install_gpl_boundary", {}, InstallReturnType(), InstallGplBoundaryExecute);
}

void InstallGplBoundaryScalar::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
