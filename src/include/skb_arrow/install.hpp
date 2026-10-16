#pragma once

// Checking and installing skb-arrow, for `skb_arrow_available()` and
// `install_skb_arrow()`. DuckDB-free, like client.hpp, so the Catch2 `tests`
// binary drives it with fake uv and fake hosts.

#include <string>

#include "skb_arrow/client.hpp"

namespace duckdb {
namespace miint {
namespace skb_arrow {

/// The release `Install` pins (the one miint is tested with), and its Python.
constexpr const char *kPinnedRelease = "skb-arrow==0.1.0";
constexpr const char *kPinnedPython = "3.14";

/// Handshake with the skb-arrow at `binary`, then shut it down. Returns its
/// `host_version`. Throws `Error` unless it speaks `kProtocolVersion` and
/// offers every capability in `MiintCapabilities()`. Its session directory
/// goes under `SegmentDirParent()`.
std::string Probe(const std::string &binary, const InterruptCheck &check);

/// Locate uv: PATH, then `$HOME/.local/bin/uv`, where uv's installer puts it.
/// Empty when neither has it.
std::string FindUv();

/// miint's cache directory, `<cache>/miint`: the parent of
/// `MiintGplBoundaryCacheDir()`. Empty when there is no cache dir.
std::string MiintCacheDir();

struct InstallOptions {
	/// A binary to try before installing, e.g. `FindSkbArrow()`. Empty: none.
	std::string existing;
	/// uv's path, e.g. `FindUv()`. Empty: uv is missing.
	std::string uv;
	/// e.g. `MiintCacheDir()`. The binary goes in its `bin`, beside
	/// gpl-boundary; uv's tool environment in `uv-tools`; the interpreter in
	/// `uv-python`. Empty: there is nowhere to install.
	std::string cache_dir;
	/// Skip `existing`: install even over one that passes `Probe`.
	bool force = false;
};

/// What `install_skb_arrow()` returns.
struct InstallReport {
	/// The skb-arrow at `path` passed `Probe`.
	bool installed = false;
	std::string path;
	/// Its `host_version`.
	std::string version;
	std::string message;
};

/// Unless `force`, return at once when `existing` passes `Probe`. Otherwise
/// run `uv tool install --managed-python --python 3.14 skb-arrow==0.1.0
/// --reinstall --force` with everything under `cache_dir`, then `Probe` what
/// it installed.
/// Serialized within the process. Failures come back in `message`; only
/// `check` throws.
InstallReport Install(const InstallOptions &options, const InterruptCheck &check);

} // namespace skb_arrow
} // namespace miint
} // namespace duckdb
