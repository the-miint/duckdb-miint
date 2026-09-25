#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_predict_proba(data, model[, name]) -> (sample_id, class, probability, sample_coverage)
//!
//! Class probabilities from a classifier, one row per (sample, class) rather
//! than a wide row per sample. Long format because the class set is a property
//! of the model, not of the query: a wide shape would need its column names
//! decided at bind time from the model's classes, and every downstream join
//! would have to know them. Long rows pivot on demand and join on `class`.
//!
//! Classifiers only. `sc_predict` is the regressor's answer, and a regressor
//! has no classes to give probabilities over.
//!
//! Probabilities for a sample sum to 1. A sample whose features the model does
//! not recognise still gets a full set of rows -- predicted from an all-zero
//! matrix -- with `sample_coverage` reporting how much of its data was usable.
class ScPredictProbaFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
