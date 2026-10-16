#include "skb_arrow/install.hpp"

#include <cstdlib>
#include <mutex>
#include <unistd.h>
#include <vector>

#include "gpl_boundary/process.hpp"

namespace duckdb {
namespace miint {
namespace skb_arrow {

namespace gb = ::duckdb::miint::gpl_boundary;

namespace {

// uv's own output is short, but a resolver failure can list many packages.
constexpr std::size_t kUvOutputCap = 8192;

std::mutex &InstallMutex() {
	static std::mutex m;
	return m;
}

} // namespace

std::string Probe(const std::string &binary, const InterruptCheck &check) {
	HostOptions options;
	options.binary = binary;
	options.segment_parent = SegmentDirParent();
	options.required = MiintCapabilities();
	Host host(options, check);
	host.Shutdown();
	return host.host_version();
}

std::string FindUv() {
	std::string on_path = gb::FindExecutableInPath("uv");
	if (!on_path.empty()) {
		return on_path;
	}
	if (const char *home = ::getenv("HOME"); home && home[0] != '\0') {
		const std::string installed = std::string(home) + "/.local/bin/uv";
		if (::access(installed.c_str(), X_OK) == 0) {
			return installed;
		}
	}
	return {};
}

std::string MiintCacheDir() {
	const std::string bin = gb::MiintGplBoundaryCacheDir();
	return bin.empty() ? bin : bin.substr(0, bin.rfind('/'));
}

InstallReport Install(const InstallOptions &options, const InterruptCheck &check) {
	InstallReport r;
	if (!options.force && !options.existing.empty()) {
		try {
			r.version = Probe(options.existing, check);
			r.installed = true;
			r.path = options.existing;
			r.message = "skb-arrow " + r.version + " already available at " + r.path + "; no install performed";
			return r;
		} catch (const Error &) {
			// Present but unusable: install over it.
		}
	}

	std::lock_guard<std::mutex> lock(InstallMutex());
	if (options.cache_dir.empty()) {
		r.message = "install_skb_arrow: cannot resolve install location ($XDG_CACHE_HOME and $HOME are both unset)";
		return r;
	}
	if (options.uv.empty()) {
		r.message = "install_skb_arrow: uv not found (on PATH or at ~/.local/bin/uv). skb-arrow is installed with "
		            "uv; get it with `curl -LsSf https://astral.sh/uv/install.sh | sh` "
		            "(https://docs.astral.sh/uv/), then run SELECT install_skb_arrow(); again";
		return r;
	}

	// uv runs only when there is no working skb-arrow, and it skips a tool it
	// already has even when that tool's binary is gone or broken. So always
	// reinstall, replacing a binary uv didn't put there (--force).
	const std::vector<std::string> argv = {options.uv,    "tool",         "install",     "--managed-python", "--python",
	                                       kPinnedPython, kPinnedRelease, "--reinstall", "--force"};
	const std::string bin_dir = options.cache_dir + "/bin";
	gb::ChildProcessOptions spawn;
	spawn.set_env = {{"UV_TOOL_BIN_DIR", bin_dir},
	                 {"UV_TOOL_DIR", options.cache_dir + "/uv-tools"},
	                 {"UV_PYTHON_INSTALL_DIR", options.cache_dir + "/uv-python"}};
	gb::RunResult uv = gb::RunCapture(argv, spawn);
	gb::TruncateHead(uv.stderr_text, kUvOutputCap);
	if (!uv.spawned) {
		r.message = "install_skb_arrow: failed to start uv at " + options.uv + ": " + uv.spawn_error;
		return r;
	}
	if (!uv.exited_normally) {
		r.message =
		    "install_skb_arrow: uv was killed by signal " + std::to_string(uv.term_signal) + ":\n" + uv.stderr_text;
		return r;
	}
	if (uv.exit_code != 0) {
		r.message = "install_skb_arrow: uv failed (exit code " + std::to_string(uv.exit_code) + "):\n" + uv.stderr_text;
		return r;
	}

	const std::string binary = bin_dir + "/skb-arrow";
	if (::access(binary.c_str(), X_OK) != 0) {
		r.message =
		    "install_skb_arrow: uv reported success, but there is no skb-arrow at " + binary + ":\n" + uv.stderr_text;
		return r;
	}
	r.path = binary;
	try {
		r.version = Probe(binary, check);
	} catch (const Error &e) {
		r.message = std::string("install_skb_arrow: installed skb-arrow at ") + binary + ", but it fails its " +
		            "handshake: " + e.what();
		return r;
	}
	r.installed = true;
	r.message = "Installed skb-arrow " + r.version + " to " + binary;
	return r;
}

} // namespace skb_arrow
} // namespace miint
} // namespace duckdb
