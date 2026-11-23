#pragma once

#include "duckdb/function/table_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

// taxonomy_violations(taxonomy)
//
// Every problem in a taxonomy (feature_id and the shared lineage schema's eight rank columns), one
// row each, never raised; docs/insdc_ncbi.md lists the kinds. The checks are TaxonomyStore's
// (miint::FindRowViolations and miint::FindMultipleParents), over the names as given. feature_id
// may be VARCHAR, BIGINT or UUID and comes back in that type.
class TaxonomyViolationsTableFunction {
public:
	static TableFunction GetFunction();
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
