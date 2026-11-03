#pragma once

#include <string>

namespace duckdb {

// A gpl-boundary release, as `gpl-boundary --version` reports it.
struct GplBoundaryRelease {
	int major;
	int minor;
	int patch;
};

// Throw IOException (naming `caller`) unless the gpl-boundary at `binary_path` reports a release of at least
// `required`, or when its release cannot be determined. `needed_for` says what needs it, for the message.
//
// The IPC handshake carries only protocol_version and each tool's schema_version, not the daemon's release, so a
// feature added without a schema bump (bowtie2's memory_mapped in 0.4.2) is detectable only by running
// `--version`. A binary's release is cached once it is new enough for a caller; otherwise (too old, or unreadable)
// it is probed again on the next call, so upgrading it in place takes effect without a restart.
void RequireGplBoundaryAtLeast(const std::string &binary_path, const char *caller, GplBoundaryRelease required,
                               const char *needed_for);

// The gpl-boundary binary FindGplBoundary finds, or an IOException (naming `caller`) that says how to install one.
std::string FindGplBoundaryOrThrow(const char *caller);

} // namespace duckdb
