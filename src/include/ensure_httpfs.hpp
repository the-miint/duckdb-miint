#pragma once

namespace duckdb {
class ClientContext;
}

namespace miint {

// Ensure the httpfs extension is loaded. https:// access via both the DuckDB
// FileSystem and the core HTTPUtil client depends on httpfs, which the miint
// extension does not auto-load. Loads it via a separate connection (works in both
// static and loadable extension builds and does not depend on the
// autoload_known_extensions setting). Throws IOException if it cannot be loaded.
void EnsureHttpfsLoaded(duckdb::ClientContext &context);

} // namespace miint
