#pragma once

#include "duckdb.hpp"

namespace duckdb {

/// Registers:
/// - `install_skb_arrow()` and `install_skb_arrow(force BOOLEAN)`: find an
///   skb-arrow that passes the handshake, or install the pinned release with
///   uv into miint's cache (`skb_arrow::Install`). Returns
///   STRUCT(installed BOOLEAN, path VARCHAR, version VARCHAR, message VARCHAR);
///   failures come back in `message`. Volatile: it has side effects.
/// - `skb_arrow_available()`: whether miint finds an skb-arrow that passes
///   the handshake (`skb_arrow::Probe`).
void RegisterSkbArrowInstallFunctions(ExtensionLoader &loader);

} // namespace duckdb
