#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// sourcetracker(feature_table, sample_metadata, ...): SourceTracker source
// attribution over a long-form feature table, computed by the embedded st3
// library. Registered only when the extension is built with st3 (MIINT_HAS_ST3).
void RegisterSourcetracker(ExtensionLoader &loader);

} // namespace duckdb
