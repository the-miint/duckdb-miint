# Sample classification and regression

Predict a sample's metadata from its feature table: which body site a sample came from, how old the infant was, whether a treatment took. A random forest is trained on a `(sample_id, feature_id, value)` count table plus one label per sample, stored in a DuckDB table, and then used to predict, to score itself by cross-validation, and to explain individual predictions.

These functions are powered by [sc](https://github.com/the-miint/sc), a dependency-free Rust reimplementation of the random forest, TreeSHAP and nested-CV parts of [q2-sample-classifier](https://github.com/qiime2/q2-sample-classifier). Defaults match scikit-learn's `RandomForestClassifier` / `RandomForestRegressor`, and results are checked against scikit-learn 1.4.2 fixtures.

## Get started

```sql
CREATE TABLE counts AS SELECT * FROM read_biom('train.biom');
CREATE TABLE meta AS SELECT sample_id, body_site FROM read_csv('metadata.tsv') WHERE body_site IS NOT NULL;

-- Training needs the same samples on both sides: keep only those in both.
CREATE VIEW train_counts AS SELECT * FROM counts SEMI JOIN meta USING (sample_id);
CREATE VIEW train_meta   AS SELECT * FROM meta SEMI JOIN counts USING (sample_id);

CREATE TABLE models AS SELECT * FROM fit_classifier('train_counts', 'train_meta', name := 'site');
CREATE TABLE new_counts AS SELECT * FROM read_biom('new.biom');
SELECT * FROM predict('new_counts', 'models', name := 'site');
```

## Table of Contents

- [Input tables](#input-tables)
- [The model table](#the-model-table)
- [Fitting](#fitting) - `fit_classifier` / `fit_regressor`
- [Predicting](#predicting) - `predict`, `predict_proba`
- [Cross-validation](#cross-validation) - `cross_validate_classifier` / `cross_validate_regressor`
- [Explaining a prediction](#explaining-a-prediction) - `shap`
- [Feature importances](#feature-importances) - `feature_importances`, `model_features`
- [Unknown and absent features](#unknown-and-absent-features) - what `sample_coverage` means
- [Reproducibility](#reproducibility)
- [Scale and memory](#scale-and-memory)
- [Citations](#citations)

### Input tables

The **count table** is the long-form `(sample_id, feature_id, value)` relation that [`read_biom`](reading.md#biom) and [`woltka_ogu`](profiling.md) produce. Absent cells are zero. NULLs are an error rather than a zero, because a NULL means a broken join upstream. Duplicate `(sample_id, feature_id)` cells are also an error rather than a silent sum; the message shows their values, since identical values suggest a join fanout and differing ones repeat measurements.

`sample_id` and `feature_id` may be VARCHAR, BIGINT or UUID; any other type is an error naming the cast to write. Ids come back as the type they went in as, so results join to typed tables without a cast.

The **metadata table** has one row per sample: `sample_id` plus the target. With one non-`sample_id` column it is used automatically; with more, name it with `target_column :=`. Both tables must describe exactly the same samples, and a mismatch is an error listing the offenders.

### The model table

A fit returns one row. Several models can share a table and are selected by `name`:

```sql
CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', target_column := 'body_site', name := 'site');
INSERT INTO models SELECT * FROM fit_regressor('counts', 'meta', target_column := 'age_months', name := 'age');
```

| column | |
|---|---|
| `name` | required at fit; how every later function selects the model |
| `model` | the algorithm, currently always `random_forest` |
| `model_blob` | the serialized forest |
| `n_samples`, `n_features`, `n_trees` | what it was trained on |
| `task` | `classification` or `regression`; decides the type of the prediction column |
| `random_state` | the seed, so a fit can be reproduced |
| `feature_id_type` | `VARCHAR`, `BIGINT` or `UUID`: the `feature_id` type data must have to use this model |

`predict`, `predict_proba` and `shap` require the data's `feature_id` to have exactly the model's `feature_id_type`. Features are matched by their text, so a BIGINT `7` would match a VARCHAR model's `'7'` but not its `'007'`; a mismatch is an error that names the cast to apply.

### Fitting

```sql
SELECT * FROM fit_classifier('counts', 'meta', name := 'site', n_estimators := 500, random_state := 42);
```

| parameter | default | |
|---|---|---|
| `name` | *required* | the model's handle |
| `target_column` | auto | required when the metadata has more than one candidate column |
| `model` | `random_forest` | the only algorithm so far |
| `n_estimators` | 100 | trees |
| `random_state` | 0 | seed |
| `n_threads` | 0 | 0 = one thread per core |
| `max_depth` | unbounded | `<= 0` is sklearn's `None` |
| `criterion` | `gini` / `squared_error` | classifier / regressor |
| `max_features` | `sqrt` / `all` | features each split may consider: `'sqrt'`, `'log2'`, `'all'`, a fraction in (0, 1], or an integer count |
| `min_samples_split` | 2 | integer count or fraction |
| `min_samples_leaf` | 1 | integer count or fraction |
| `min_weight_fraction_leaf`, `min_impurity_decrease` | 0.0 | |
| `bootstrap` | true | |
| `max_samples` | all | only with `bootstrap := true` |

As in sklearn, integer versus float matters: `max_features := 10` means ten features, `0.3` means thirty percent.

`max_features` is not a filter on the table. At every split of every tree, a fresh random set of that many candidate features is drawn and the best split among them is taken. That randomness decorrelates the trees, which is where a forest's advantage comes from; `'all'` removes it and is slower without usually scoring better.

### Predicting

```sql
SELECT * FROM predict('new_counts', 'models', name := 'site');
-- sample_id, prediction, sample_coverage

SELECT * FROM predict_proba('new_counts', 'models', name := 'site');
-- sample_id, class, probability, sample_coverage   (one row per sample and class)
```

`prediction` is VARCHAR for a classifier and DOUBLE for a regressor, read from the model's `task` when the query is planned. A wrong model name, or `predict_proba` on a regressor, therefore fails before any data is read. A forest's probabilities are vote fractions, not calibrated posteriors.

### Cross-validation

How well does this recipe generalize to unseen samples? The samples are cut into `n_folds` folds (stratified for a classifier). Each fold is held out in turn while a forest fit on the rest predicts it, so every sample gets one prediction from a model that never saw it.

```sql
CREATE TABLE cv AS SELECT * FROM cross_validate_classifier('counts', 'meta', target_column := 'body_site', n_folds := 5, random_state := 42);
SELECT DISTINCT metric, mean_score, std_score, fold_scores FROM cv;
```

| column | |
|---|---|
| `sample_id`, `prediction`, `actual` | the out-of-fold prediction and the target it was scored against |
| `probabilities`, `classes` | classifier only; lists aligned position by position |
| `metric` | `accuracy` (classifier) or `mse` (regressor; lower is better) |
| `mean_score`, `std_score` | over the per-fold scores |
| `fold_scores` | each fold's score |

**`mean_score` is the generalization estimate**: cross-validated accuracy, or MSE, where `sqrt(mean_score)` is an RMSE in the target's units. It averages folds, not samples, as sklearn does.

- `std_score` is computed with `ddof=0` (population standard deviation) to match `q2-sample-classifier`.
- **The per-fold forests are discarded.** This measures the recipe, not a model; fit the model you will use with `fit_classifier` or `fit_regressor` on all the data.
- **The four summary columns repeat on every row**, because a table function returns one relation and a separate summary function would have to rerun every fold. `SELECT DISTINCT` recovers the single summary row.
- **`n_folds` is also accepted as `cv`**, q2's and sklearn's name for it; passing both is an error.
- **`parameter_tuning := true`** runs a randomized search inside each training fold, so the held-out fold stays unseen. It costs about `n_folds × 20` extra fits, and the chosen parameters are not yet reported.

```sql
-- Confusion matrix, no join required.
SELECT actual, prediction, count(*) FROM cv GROUP BY ALL ORDER BY 1, 2;
-- The least confident calls.
SELECT sample_id, prediction, actual, list_max(probabilities) AS confidence FROM cv ORDER BY confidence LIMIT 10;
```

### Explaining a prediction

TreeSHAP attributes one prediction to individual features: how far each pushed this sample from the model's baseline. It is additive: `base_value + sum(shap_value)` over a sample's features equals the predicted value (the continuous target for regression, or the class vote fraction for classification).

```sql
SELECT * FROM shap('counts', 'models', name := 'site', top_k := 10);
-- sample_id, class, feature_id, shap_value, base_value, sample_coverage
```

| parameter | default | |
|---|---|---|
| `predicted_class_only` | true | with two classes the other is an exact mirror; with three or more, `false` also explains the other classes |
| `top_k` | off | keep the k strongest features per sample, split evenly between those pushing toward the prediction and against it |
| `max_attributions` | 10,000,000 | attributions one batch may hold (samples × classes × features, 8 bytes each, so ~80 MB) |
| `batch_size` | derived | samples per batch, overriding `max_attributions`; passing both is an error |

- **`base_value`** is the model's average output over the training data, the same on every row.
- **Absent features get attributions too.** A feature the model knows but the sample lacks is a zero, and a zero can push a prediction as hard as a large count; the absence of adult taxa is evidence of an infant.
- **`top_k` shrinks the output, not the computation**, since every attribution is computed before the strongest are picked.
- **Rows come out in waterfall order**: `sample_id`, `class`, then `shap_value` descending.
- **Memory is bounded by the batch**, so explaining 20,000 samples costs no more memory than 2,000.

### Feature importances

```sql
SELECT * FROM feature_importances('models', name := 'site') ORDER BY importance DESC LIMIT 20;
-- feature_id, importance
SELECT * FROM model_features('models', name := 'site');
-- feature_id, column_index
```

`feature_importances` is mean decrease in impurity, summing to 1. It describes the model, not any prediction, and is biased toward features with many distinct values; the mean absolute SHAP value per feature is a less biased global ranking. `model_features` lists the features the model expects, in its own column order.

### Unknown and absent features

- **A feature the model never saw is dropped**, with a warning; no tree splits on it.
- **A feature the model knows but the sample lacks is zero**, which is what absence means in a sparse table.

`sample_coverage` is the share of a sample's observed cells whose feature the model knows. A low value means the model is being applied to different data (another reference database, pipeline or 16S region). A 0.0 means the sample is predicted from an all-zero row: a well-formed but uninformed prediction.

### Reproducibility

`random_state` fixes everything sc does: bootstrap draws, feature subsampling, cross-validation folds and the tuning search. Output is bit-identical at any `n_threads`. Sample and feature ids are sorted before encoding, so the order DuckDB scans rows in cannot change the result, and prediction always uses the model's own feature order.

```sql
-- Same seed, different thread counts: identical models.
CREATE TABLE a AS SELECT * FROM fit_regressor('counts', 'meta', name := 'a', random_state := 7);
CREATE TABLE b AS SELECT * FROM fit_regressor('counts', 'meta', name := 'b', random_state := 7, n_threads := 8);
SELECT (SELECT md5(model_blob) FROM a) = (SELECT md5(model_blob) FROM b) AS identical;
```

### Scale and memory

Cost scales with the feature count, and microbiome tables are long-tailed, so a prevalence filter is the biggest lever:

```sql
CREATE VIEW common AS
  SELECT c.* FROM counts c
  SEMI JOIN (SELECT feature_id FROM counts GROUP BY 1
             HAVING count(DISTINCT sample_id) >= 0.01 * (SELECT count(DISTINCT sample_id) FROM counts)) USING (feature_id);
```

- **Fitting** spends most of its time loading the table into sc, so the cell count matters more than `n_estimators`. Memory grows with threads; `n_threads := 1` trades time for memory.
- **`shap`** memory is capped by `max_attributions`, whatever the sample count; fewer features is the way to go further.
- **wasm** shares one 4 GiB memory with DuckDB and runs single-threaded, so expect a much lower ceiling than native.

### Citations

- Breiman, L. (2001). Random Forests. *Machine Learning* 45(1), 5–32.
- Lundberg, S. M., Erion, G. G., & Lee, S.-I. (2018). Consistent Individualized Feature Attribution for Tree Ensembles. [arXiv:1802.03888](https://arxiv.org/abs/1802.03888) — the path-dependent TreeSHAP algorithm `shap` implements. Cite this if you publish SHAP values.
- Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley, G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning tools for microbiome classification and regression. *Journal of Open Source Software* 3(30), 934.
- Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR* 12, 2825–2830 — the reference implementation sc's defaults and oracle fixtures follow.
