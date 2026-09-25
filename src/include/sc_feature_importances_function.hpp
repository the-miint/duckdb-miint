#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_feature_importances(model) -> (feature_id, importance)
//!
//! Mean-decrease-in-impurity importances for a fitted model, labelled with the
//! feature ids the model was trained on. sc returns the two as separate arrays
//! -- a bare `Float64` vector and the training vocabulary -- and neither is
//! useful alone: the numbers have no names, and the names have no weights.
//!
//! MDI is what the forest *used*, measured on training data, and it carries the
//! usual biases: it favours features with more distinct values, and correlated
//! features split their importance between them. Microbiome data hits both --
//! abundant taxa have more split points than rare ones, and compositional data
//! is heavily co-varying -- so a low score is not evidence a taxon is
//! irrelevant.
class ScFeatureImportancesFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
