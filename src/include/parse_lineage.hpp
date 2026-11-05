#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// parse_lineage(lineage VARCHAR, domain_letter := 'd') -> the shared lineage schema's
// eight rank columns as a STRUCT. The parsing rules are miint::ParseLineageString.
class ParseLineageFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
