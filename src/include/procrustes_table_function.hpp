#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// Register the `procrustes(reference, other)` table function. Aligns two
// long-format ordination tables (sample_id, axis, coordinate) into a shared
// standardized frame and emits both transformed clouds plus the fit disparity
// M^2 and (full-overlap case) a Monte Carlo PROTEST p-value. Depends only on the
// Eigen-backed procrustes core, so it is registered unconditionally (no UniFrac
// gate). See src/procrustes_core.cpp for the numeric core.
void RegisterProcrustes(ExtensionLoader &loader);

} // namespace duckdb
