#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// taxonomy_violations(taxonomy)
//
// Every problem in a taxonomy in the shared lineage schema (feature_id and the eight rank
// columns), one row each, never raised: an empty or NULL id, a value that is not a bare name, a
// repeated id, a gap, and a label under more than one parent. The checks are TaxonomyStore's
// (miint::FindRowViolations and miint::FindMultipleParents), over the names as given. feature_id
// may be VARCHAR, BIGINT or UUID and comes back in that type. See docs/insdc_ncbi.md.
class TaxonomyViolationsTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
