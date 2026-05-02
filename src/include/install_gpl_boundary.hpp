#pragma once

#include "duckdb.hpp"

namespace duckdb {

/// Scalar function `install_gpl_boundary()` — downloads the upstream
/// install.sh from
/// https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh
/// and runs it with INSTALL_DIR pointed at miint's cache (so the result is
/// findable by `FindGplBoundary` without the user editing PATH).
///
/// Returns STRUCT(installed BOOLEAN, path VARCHAR, version VARCHAR, message VARCHAR).
/// `installed=true` means the binary is now available at `path`. `installed=false`
/// means the install was not performed (or failed); `message` carries the reason.
///
/// Idempotent: when `gpl-boundary` is already discoverable (override path,
/// cache dir, or PATH), short-circuits with `installed=true` and a message
/// noting no changes were made.
class InstallGplBoundaryScalar {
public:
	static ScalarFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
