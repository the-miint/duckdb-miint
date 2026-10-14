#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_shap(data, model[, name, predicted_class_only, top_k, max_attributions, batch_size, n_threads])
//!   -> (sample_id, class, feature_id, shap_value, base_value, sample_coverage)
//!
//! Path-dependent TreeSHAP: how much each feature pushed one sample's prediction
//! away from the model's baseline. Additive by construction --
//! `base_value + sum(shap_value)` over a sample's features equals the
//! prediction: the class probability for a classifier, the predicted value for a
//! regressor.
//!
//! A classifier explains every class, but only the PREDICTED class is returned
//! by default. With two classes the other one is an exact mirror (every value
//! negated), so nothing is lost. `predicted_class_only := false` returns them
//! all, which matters with three or more classes, where "why not the runner-up"
//! is not recoverable from the winner alone. `class` is NULL for a regressor.
//!
//! `top_k := k` keeps k features per sample (per class when all classes are
//! returned): the strongest pushes toward the prediction and the strongest
//! against it, split evenly, with the positive side taking any odd one out. If a
//! sample runs short on one side the other fills the gap, so k rows come back
//! whenever the sample has k non-zero attributions. Off by default, because
//! truncating silently hides data and breaks the additivity above.
//!
//! Features the model knows but the sample lacks are explained too: absent is
//! zero, and a zero can push a prediction as hard as a count. Features the model
//! never saw are dropped, exactly as in sc_predict -- no tree splits on them, so
//! their attribution would be 0 anyway. `sample_coverage` is the share of the
//! sample's features the model knows, as in sc_predict; at 0 the attributions
//! still add up, but they explain an all-zero row.
//!
//! Rows come out in waterfall order,
//!   ORDER BY sample_id, class, shap_value DESC
//! with ties in model column order. That is the order rows are produced in, not
//! a property DuckDB tracks: a plain SELECT or CREATE TABLE AS keeps it, joins
//! and aggregates may not, and an outer ORDER BY simply re-sorts.
//!
//! SHAP is computed densely -- every class of every sample in a call, before any
//! `top_k` filtering -- so samples are explained in batches and rows stream out
//! one batch at a time; memory is one batch however many samples there are.
//! `max_attributions` (default 10M, ~80 MB) is how many attributions a batch may
//! hold, and the batch size follows: max_attributions / (classes x features)
//! samples. `batch_size := n` fixes it instead; passing both is an error. A
//! sample that alone needs more than max_attributions still runs, by itself,
//! with a warning. Batching changes no value and no row order, because each
//! sample is explained independently. `top_k` does not reduce what is computed;
//! it keeps the output small, which matters once the result is stored.
class ScShapFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
