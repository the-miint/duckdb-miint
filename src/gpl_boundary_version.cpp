#include "gpl_boundary_version.hpp"

#include "gpl_boundary/process.hpp"
#include "gpl_boundary/session.hpp"

#include "duckdb/common/exception.hpp"

#include <map>
#include <mutex>
#include <tuple>
#include <unistd.h>

namespace duckdb {

namespace {

namespace gb = ::duckdb::miint::gpl_boundary;

// Run `<binary> --version` and parse its {"gpl_boundary":"X.Y.Z", ...}. False when it cannot be run or parsed.
bool ProbeRelease(const std::string &binary_path, GplBoundaryRelease &out) {
	try {
		std::vector<std::string> argv = {binary_path, "--version"};
		gb::ChildProcess child(argv);
		auto drain = [](int fd, std::string &dst) {
			if (fd < 0) {
				return;
			}
			char buf[256];
			ssize_t n;
			while ((n = ::read(fd, buf, sizeof(buf))) > 0) {
				dst.append(buf, static_cast<size_t>(n));
			}
		};
		std::string out_text;
		std::string err_text;
		drain(child.stdout_fd(), out_text);
		drain(child.stderr_fd(), err_text);
		child.Wait();
		return gb::ParseGplBoundaryVersion(out_text, out.major, out.minor, out.patch);
	} catch (...) {
		// Any spawn/read failure is "release undetermined", which the caller refuses.
		return false;
	}
}

std::string ToString(const GplBoundaryRelease &r) {
	return std::to_string(r.major) + "." + std::to_string(r.minor) + "." + std::to_string(r.patch);
}

} // namespace

void RequireGplBoundaryAtLeast(const std::string &binary_path, const char *caller, GplBoundaryRelease required,
                               const char *needed_for) {
	static std::mutex lock;
	static std::map<std::string, GplBoundaryRelease> known;

	const auto satisfies = [&](const GplBoundaryRelease &r) {
		return std::tie(r.major, r.minor, r.patch) >= std::tie(required.major, required.minor, required.patch);
	};
	{
		std::lock_guard<std::mutex> guard(lock);
		auto it = known.find(binary_path);
		if (it != known.end() && satisfies(it->second)) {
			return;
		}
	}
	// Not known to be new enough: ask the binary again, so one upgraded in this process (install_gpl_boundary(true)
	// replaces it at the same path) is seen without a restart.
	GplBoundaryRelease found {0, 0, 0};
	const bool parsed = ProbeRelease(binary_path, found);
	if (parsed) {
		std::lock_guard<std::mutex> guard(lock);
		known[binary_path] = found;
	}
	if (parsed && satisfies(found)) {
		return;
	}
	const std::string found_text =
	    parsed ? ToString(found) : std::string("unknown (could not parse `gpl-boundary --version`)");
	throw IOException("%s requires gpl-boundary >= %s (for %s), but the binary at %s reports version %s. Upgrade with "
	                  "`SELECT install_gpl_boundary(true);`, or point MIINT_GPL_BOUNDARY_PATH at a >= %s "
	                  "binary.",
	                  caller, ToString(required), needed_for, binary_path, found_text, ToString(required));
}

} // namespace duckdb
