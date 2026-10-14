#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

//! sc_cross_validate_classifier / sc_cross_validate_regressor(data, metadata[,
//!   target_column, n_folds, parameter_tuning, random_state, ... forest params])
//!   -> (sample_id, prediction, actual, [probabilities, classes,] metric,
//!       mean_score, std_score, fold_scores)
//!
//! How well does this model generalise to samples it has never seen?
//!
//! The samples are shuffled once and cut into `cv` folds. Each fold is
//! held out in turn while a forest is fit on the rest and used to predict it, so
//! every sample ends up with exactly one prediction from a model that never saw
//! it -- an out-of-fold (OOF) prediction. Folds are stratified for a classifier
//! (each fold keeps the class proportions, so no fold can miss a class) and
//! plain for a regressor, matching q2 and sklearn.
//!
//! `actual` is the target the fold was scored against, carried through so the
//! result scores itself: a confusion matrix, a residual, or an accuracy needs no
//! join back to the metadata relation.
//!
//! The per-fold forests are scored and discarded: this measures the RECIPE, not
//! a model. Fit the model you will actually use with sc_fit on all the data.
//!
//! `probabilities` is the OOF probability of each class, and `classes` names
//! them: both are lists in the same order on every row, so `probabilities[i]` is
//! the probability of `classes[i]`, and the class with the highest probability
//! is the one in `prediction`. Positional lists rather than a map keep the class
//! order sc reported, and keep one row per sample -- a (sample, class) grain
//! would multiply every score column too. A regressor has no classes to be
//! probable, and the task is fixed by which function was called, so the two
//! columns are absent from its schema rather than present and always NULL.
//!
//! `metric` is 'accuracy' (classifier) or 'mse' (regressor). `mean_score` and
//! `std_score` are over the per-fold scores -- the same mean and population
//! standard deviation q2 prints -- and `fold_scores` holds each fold's score, so
//! `unnest(fold_scores)` shows the spread. All four are properties of the whole
//! run, so they repeat on every row; `SELECT DISTINCT` collapses them. Two grains
//! in one result because cross-validation fits `cv` forests: asking for the
//! scores and the predictions separately would pay that cost twice.
//!
//! Scoring is the mean of the per-fold scores, not a pooled score over all
//! predictions. They differ whenever folds are unequal in size, and the
//! fold-mean is what sklearn's cross_val_score reports.
//!
//! `parameter_tuning := true` searches the built-in hyperparameter grid inside
//! each training fold and refits with the winner (q2's `parameter_tuning`). The
//! search never sees the held-out fold, which is what keeps the score honest;
//! it is also cv x grid times the work.
class ScCrossValidateFunction {
public:
	static void Register(ExtensionLoader &loader);
};

} // namespace duckdb
