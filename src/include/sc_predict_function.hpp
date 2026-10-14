#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_predict(data, model) -> (sample_id, prediction, sample_coverage)
//!
//! One function for both tasks, mirroring sc's own single `sc_predict`, which
//! dispatches on the model's recorded task. The prediction column is VARCHAR
//! for a classifier and DOUBLE for a regressor, resolved at bind time from the
//! model relation's `task` column -- the same way `read_csv` decides its schema
//! by sniffing the file. Asking the caller to pick would only create a way to
//! pick wrong.
//!
//! No metadata argument: a fitted model already carries the trees, the task and
//! the feature vocabulary. Ground truth joins on afterwards in SQL, which keeps
//! scoring open to any metric.
//!
//! Features are encoded against the model's vocabulary, never re-derived from
//! the prediction data: two datasets differing by one feature in each direction
//! have the same width, and sc validates only the width, so re-deriving would
//! return confident nonsense with no error anywhere. `sample_coverage` reports
//! matched/observed per sample so a caller can see how much of their data the
//! model could actually use.
class ScPredictFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
