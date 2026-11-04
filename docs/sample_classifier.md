# Sample classification and regression

Predict a sample's metadata from its feature table: which body site a sample came from, how old the infant was, whether a treatment took. A random forest is trained on a `(sample_id, feature_id, value)` count table plus one label per sample, stored in a DuckDB table, and then used to predict, to score itself by cross-validation, and to explain individual predictions.

These functions are powered by [sc](https://github.com/the-miint/sc), a Rust reimplementation of the random forest and nested cross-validation in [q2-sample-classifier](https://github.com/qiime2/q2-sample-classifier) (QIIME 2's sample classifier, "q2" below), plus SHAP (SHapley Additive exPlanations) feature attributions, computed exactly for tree ensembles with path-dependent TreeSHAP, to explain predictions. Defaults match scikit-learn's `RandomForestClassifier` / `RandomForestRegressor`, except that `random_state` defaults to 0 rather than unseeded, so a fit is reproducible.

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
- [Unknown and absent features](#unknown-and-absent-features) - what `known_feature_fraction` means
- [Reproducibility](#reproducibility) - `random_state`, threads, and the sc commit
- [Scale and memory](#scale-and-memory)
- [Citations](#citations)

### Input tables

The **count table** is the long-form `(sample_id, feature_id, value)` relation that [`read_biom`](reading.md#biom) and [`woltka_ogu`](profiling.md) produce. Absent cells are zero. NULLs are an error rather than a zero, because a NULL means a broken join upstream. Duplicate `(sample_id, feature_id)` cells are also an error rather than a silent sum; the message shows their values, since identical values suggest a join that duplicated rows and differing ones repeat measurements.

`sample_id` and `feature_id` may be VARCHAR, BIGINT or UUID; any other type is an error naming the cast to write. Ids come back as the type they went in as, so results join to typed tables without a cast. They are ordered by their text, not their value: a BIGINT `10` sorts before `2`, and rows come back in that order. An outer `ORDER BY` restores numeric order.

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
| `model` | the algorithm: `random_forest` |
| `model_blob` | the serialized forest |
| `n_samples`, `n_features`, `n_estimators` | what it was trained on |
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
| `model` | `random_forest` | the algorithm; `random_forest` is the one accepted |
| `n_estimators` | 100 | trees, at most 100,000 |
| `random_state` | 0 | seed |
| `threads` | 0 | 0 = DuckDB's `threads` setting, otherwise a count; ignored in DuckDB-Wasm, which has one |
| `max_depth` | unbounded | `<= 0` is scikit-learn's `None` |
| `criterion` | `gini` / `squared_error` | classifier: `'gini'` or `'entropy'`; regressor: `'squared_error'` |
| `max_features` | `sqrt` / `all` | features each split may consider: `'sqrt'`, `'log2'`, `'all'` (or `'none'`), a fraction in (0, 1], or an integer count |
| `min_samples_split` | 2 | an integer count >= 2, or a fraction in (0, 1] |
| `min_samples_leaf` | 1 | an integer count >= 1, or a fraction in (0, 1] |
| `min_weight_fraction_leaf`, `min_impurity_decrease` | 0.0 | in [0, 0.5], and finite and >= 0 |
| `bootstrap` | true | |
| `max_samples` | all | each tree's bootstrap sample: `'all'` (or `'none'`), a fraction in (0, 1], or a positive count; only with `bootstrap := true` |

As in scikit-learn, integer versus float matters: `max_features := 10` means ten features, `0.3` means thirty percent.

`max_features` is not a filter on the table. At every split of every tree, a fresh random set of that many candidate features is drawn and the best split among them is taken. That randomness decorrelates the trees, which is where a forest's advantage comes from; `'all'` removes it.

### Predicting

```sql
SELECT * FROM predict('new_counts', 'models', name := 'site');
-- sample_id, prediction, known_feature_fraction

SELECT * FROM predict_proba('new_counts', 'models', name := 'site');
-- sample_id, class, probability, known_feature_fraction   (one row per sample and class)
```

| parameter | default | |
|---|---|---|
| `name` | | the model to use; required when the model relation holds more than one |

`prediction` is VARCHAR for a classifier and DOUBLE for a regressor, read from the model's `task` when the query is planned. A wrong model name, or `predict_proba` on a regressor, therefore fails before any data is read. A forest's probability is the mean of its trees' leaf class proportions, not a calibrated posterior.

### Cross-validation

How well does this recipe generalize to unseen samples? The samples are cut into `n_folds` folds (stratified for a classifier). Each fold is held out in turn while a forest fit on the rest predicts it, so every sample gets one prediction from a model that never saw it.

```sql
CREATE TABLE cv AS SELECT * FROM cross_validate_classifier('counts', 'meta', target_column := 'body_site', n_folds := 5, random_state := 42);
SELECT DISTINCT metric, mean_score, std_score, fold_scores FROM cv;
```

| parameter | default | |
|---|---|---|
| `target_column` | auto | as at fit |
| `n_folds` | 5 | at least 2; also accepted as `cv`, q2's and scikit-learn's name for it; passing both is an error |
| `parameter_tuning` | false | see below |
| `n_estimators`, `random_state`, `threads` and every forest parameter from [Fitting](#fitting) | as at fit | the recipe each fold's forest is fit with |

| column | |
|---|---|
| `sample_id`, `prediction`, `actual` | the out-of-fold prediction and the target it was scored against |
| `probabilities`, `classes` | classifier only; lists aligned position by position |
| `metric` | `accuracy` (classifier) or `mse` (regressor; lower is better) |
| `mean_score`, `std_score` | over the per-fold scores |
| `fold_scores` | each fold's score |

**`mean_score` is the generalization estimate**: cross-validated accuracy, or MSE, where `sqrt(mean_score)` is an RMSE in the target's units. It averages folds, not samples, as scikit-learn does.

- `std_score` is computed with `ddof=0` (population standard deviation) to match `q2-sample-classifier`.
- **The per-fold forests are discarded.** This measures the recipe, not a model; fit the model you will use with `fit_classifier` or `fit_regressor` on all the data.
- **The four summary columns repeat on every row.** `SELECT DISTINCT` recovers the single summary row.
- **`parameter_tuning := true`** runs a randomized search inside each training fold, so the held-out fold stays unseen. Each training fold scores 20 candidate settings with an inner `n_folds`-fold cross-validation, so it costs about 20 × `n_folds`² extra fits (500 at five folds). The chosen parameters are not returned.

```sql
-- Confusion matrix, no join required.
SELECT actual, prediction, count(*) FROM cv GROUP BY ALL ORDER BY 1, 2;
-- The least confident calls.
SELECT sample_id, prediction, actual, list_max(probabilities) AS confidence FROM cv ORDER BY confidence LIMIT 10;
```

### Explaining a prediction

TreeSHAP attributes one prediction to individual features: how far each pushed this sample from the model's baseline. It is additive: `base_value + sum(shap_value)` over a sample's features equals the predicted value (the continuous target for regression, or the class probability for classification).

```sql
SELECT * FROM shap('counts', 'models', name := 'site', top_k := 10);
-- sample_id, class, feature_id, shap_value, base_value, known_feature_fraction
```

| parameter | default | |
|---|---|---|
| `name` | | as at prediction |
| `predicted_class_only` | true | with two classes the other is an exact mirror; with three or more, `false` also explains the other classes |
| `top_k` | off | keep the k strongest features per sample, or per sample and class with `predicted_class_only := false`, split evenly between those pushing toward the prediction and against it |
| `max_attributions` | 10,000,000 | attributions one batch may hold (samples × classes × features, 8 bytes each, so ~80 MB) |
| `batch_size` | derived | samples per batch, overriding `max_attributions`; passing both is an error |
| `threads` | 0 | as at fit |

- **`base_value`** is the model's average output over the training data for the row's class (a regressor has one), so it repeats on every row of that class.
- **Absent features get attributions too.** A feature the model knows but the sample lacks is a zero, and a zero can push a prediction as hard as a large count; the absence of adult taxa is evidence of an infant.
- **`top_k` shrinks the output, not the computation**, since every attribution is computed before the strongest are picked.
- **Rows come out in waterfall order**: `sample_id` (as text; see [Input tables](#input-tables)), `class`, then `shap_value` descending.
- **Attributions are computed one batch at a time**, sized by `max_attributions`, so attribution memory stays under the cap however many samples there are. The exception is a model whose classes × features exceed the cap on their own: every sample then runs alone, above the cap, with a warning.

### Feature importances

```sql
SELECT * FROM feature_importances('models', name := 'site') ORDER BY importance DESC LIMIT 20;
-- feature_id, importance
SELECT * FROM model_features('models', name := 'site');
-- feature_id, column_index
```

Both take `name`, as `predict` does. `feature_importances` is mean decrease in impurity, summing to 1. It describes the model, not any prediction, and is biased toward features with many distinct values; the mean absolute SHAP value per feature is a less biased global ranking. `model_features` lists the features the model expects, in its own column order.

### Unknown and absent features

- **A feature the model never saw is dropped**, with a warning; no tree splits on it.
- **A feature the model knows but the sample lacks is zero**, which is what absence means in a sparse table.

`known_feature_fraction` is the share of a sample's observed cells whose feature the model knows. A low value means the model is being applied to different data (another reference database, pipeline or 16S region). A 0.0 means the sample is predicted from an all-zero row: a well-formed but uninformed prediction.

### Reproducibility

`random_state` fixes everything sc does: bootstrap draws, feature subsampling, cross-validation folds and the tuning search. Output is bit-identical at any `threads`. Sample and feature ids are sorted before encoding, so the order DuckDB scans rows in cannot change the result, and prediction always uses the model's own feature order.

```sql
-- Same seed, different thread counts: identical models.
CREATE TABLE a AS SELECT * FROM fit_regressor('counts', 'meta', name := 'a', random_state := 7);
CREATE TABLE b AS SELECT * FROM fit_regressor('counts', 'meta', name := 'b', random_state := 7, threads := 8);
SELECT (SELECT md5(model_blob) FROM a) = (SELECT md5(model_blob) FROM b) AS identical;
```

Reproducing a fit also needs the same sc. `SELECT version FROM miint_versions() WHERE library = 'sc'` reports the sc commit this build embeds.

### Scale and memory

Cost scales with the feature count, and microbiome tables are long-tailed, so a prevalence filter is the biggest lever:

```sql
CREATE VIEW common AS
  SELECT c.* FROM counts c
  SEMI JOIN (SELECT feature_id FROM counts GROUP BY 1
             HAVING count(DISTINCT sample_id) >= 0.01 * (SELECT count(DISTINCT sample_id) FROM counts)) USING (feature_id);
```

- **Cancelling:** Ctrl-C and a `max_execution_time` deadline stop a function while it reads its data relation, and `shap` between batches too; `fit_*` and `cross_validate_*` cannot be stopped once sc is running.
- **Very large values are refused.** A regression target above about 1e154 in magnitude overflows sc's node statistics, and a feature value beyond the float32 range (about ±3.4e38) cannot be split on, so the fit errors and asks you to rescale. scikit-learn fits the first; it refuses the second too.
- **`shap`** holds the input table in full, as prediction does, plus one batch of attributions (see `max_attributions` under [Explaining a prediction](#explaining-a-prediction)); fewer features is the way to go further.
- **DuckDB-Wasm** gives DuckDB and the extension one 4 GiB memory and one thread, so expect a much lower ceiling than native. `threads` has no effect there; results are the same as native. `cross_validate_*` needs the `wasm_eh` build, the one miint publishes, and refuses to run in a build without wasm-native exceptions.

### Citations

- Breiman, L. (2001). Random Forests. *Machine Learning* 45(1), 5–32.
- Lundberg, S. M., Erion, G. G., & Lee, S.-I. (2018). Consistent Individualized Feature Attribution for Tree Ensembles. [arXiv:1802.03888](https://arxiv.org/abs/1802.03888) — the path-dependent TreeSHAP algorithm `shap` implements. Cite this if you publish SHAP values.
- Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley, G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning tools for microbiome classification and regression. *Journal of Open Source Software* 3(30), 934.
- Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR* 12, 2825–2830 — the reference implementation sc's defaults and oracle fixtures follow.
