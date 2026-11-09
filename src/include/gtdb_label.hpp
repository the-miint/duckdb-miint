#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// gtdb_label(name VARCHAR) -> STRUCT(base, suffix), for GTDB names above species, and
// gtdb_species_label(name VARCHAR) -> STRUCT(genus_base, genus_suffix, epithet,
// epithet_suffix), for species. Each splits a name into its base and the suffix GTDB
// marks a split taxon with (Bacillota_A: Bacillota and A); a name that does not fit gives
// NULL. They only split: grouping keys are the caller's to build from the parts.
class GtdbLabelFunctions {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
