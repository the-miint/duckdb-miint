#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_fit_classifier(data, metadata, ...) / sc_fit_regressor(data, metadata, ...)
//!
//! Table functions that fit an sc RandomForest over a long-format count table
//! and a long-format metadata table, returning the trained model as a BLOB.
//!
//!   data     -- (sample_id, feature_id, value), e.g. read_biom(...) or
//!               woltka_ogu_per_sample(...). sc is coupled to this triple, not
//!               to any file format.
//!   metadata -- (sample_id, <target_column>), one row per sample. Long-format
//!               metadata should be narrowed to a single variable first, e.g.
//!               `WHERE variable = 'month'`.
//!
//! Both relations must cover exactly the same sample set. Mismatches, duplicate
//! cells, duplicate labels and NULLs are errors rather than silent repairs: in
//! SQL the caller can express whichever join or aggregation they meant in one
//! line, so anything reaching here unbalanced is a mistake worth surfacing.
class ScFitFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
