#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_model_features(model[, name]) -> (feature_id, column_index)
//!
//! The vocabulary a model was trained on, in its own column order. This is what
//! `sc_predict` encodes incoming data against, so it is the reference for
//! answering "why is sample_coverage low" -- the usual cause is a naming
//! mismatch (case, prefix, a different reference database) rather than missing
//! biology.
//!
//!   -- in the model but absent from new data
//!   SELECT feature_id FROM sc_model_features('models', name := 'rf')
//!   EXCEPT SELECT DISTINCT feature_id FROM new_data;
//!
//!   -- in new data but unknown to the model, and therefore dropped
//!   SELECT DISTINCT feature_id FROM new_data
//!   EXCEPT SELECT feature_id FROM sc_model_features('models', name := 'rf');
//!
//! `column_index` is the position the id occupies in the model's matrix. It is
//! the thing that has to line up between fit and predict, and sc validates only
//! the feature *count*, so seeing the mapping is the only way to check it.
//!
//! The same list is reachable via sc_feature_importances, but that walks every
//! tree in the forest to compute MDI -- wasted work when the question is which
//! features exist, and not free at 200k features across 500 trees.
class ScModelFeaturesFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
