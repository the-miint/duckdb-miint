#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// lineage_rank(ncbi_rank VARCHAR) -> VARCHAR: the shared lineage schema's column that an
// NCBI rank collapses into (miint::FindRankColumn), or NULL for a rank the schema drops.
class LineageRankFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
