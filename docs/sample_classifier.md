# Sample classification and regression

Predict a sample's metadata from its feature table: which body site a sample came from, how old the infant was, whether a treatment took. A random forest is trained on a `(sample_id, feature_id, value)` count table plus one label per sample, stored in a DuckDB table as a blob, and then used to predict, to score itself honestly by cross-validation, and to explain individual predictions.

These functions are powered by [sc](https://github.com/the-miint/sc) (sample classifier), a dependency-free Rust reimplementation of the random forest, TreeSHAP and nested-CV parts of [q2-sample-classifier](https://github.com/qiime2/q2-sample-classifier). Defaults match scikit-learn's `RandomForestClassifier` / `RandomForestRegressor`, and the results are checked against scikit-learn 1.4.2 fixtures in sc's own oracle tests.

## Table of Contents

- [Input tables](#input-tables) - the count table and the metadata table both functions expect
- [The model table](#the-model-table) - what a fit produces and how models are selected later
- [Fitting](#fitting) - `sc_fit_classifier` / `sc_fit_regressor`
- [Predicting](#predicting) - `sc_predict`, `sc_predict_proba`
- [Cross-validation](#cross-validation) - `sc_cross_validate_classifier` / `sc_cross_validate_regressor`
- [Explaining a prediction](#explaining-a-prediction) - `sc_shap`
- [Feature importances](#feature-importances) - `sc_feature_importances`, `sc_model_features`
- [Identifier and label types](#identifier-and-label-types) - what type each returned column takes, and why
- [Unknown and absent features](#unknown-and-absent-features) - what `sample_coverage` means
- [Reproducibility](#reproducibility) - what `random_state` does and does not fix
- [Scale and memory](#scale-and-memory) - what costs what, and the prevalence filter that matters most
- [Citations](#citations)

### Input tables

The **count table** is the same long-form `(sample_id, feature_id, value)` relation [`read_biom`](reading.md#biom) and [`woltka_ogu`](profiling.md) produce. `sample_id` and `feature_id` may be **VARCHAR, BIGINT or UUID** — anything else is a bind error — and `value` any numeric type. Absent cells are zero; NULLs are rejected rather than treated as zero, because a NULL means a broken join upstream.

```sql
CREATE TABLE counts AS SELECT * FROM read_biom('table.biom');
```

Duplicate `(sample_id, feature_id)` cells are an error, not a silent sum. Identical values usually mean a join fanout, differing values usually mean repeat measurements, and the two need opposite repairs — so the message shows the values and leaves the decision to you.

The **metadata table** is one row per sample: a `sample_id` column and the target column. With exactly one non-`sample_id` column it is used automatically; with more, name it with `target_column :=` rather than letting the wrong variable be picked silently.

```sql
CREATE TABLE meta AS SELECT sample_id, body_site, age_months FROM read_csv('metadata.tsv');
```

Both relations must describe exactly the same samples. A mismatch in either direction is an error listing the offenders, because a truncated metadata export otherwise trains a model on fewer samples than you think.

### The model table

A fit returns one row. Store it — several models can live in one table, selected by name:

```sql
CREATE TABLE models AS SELECT * FROM sc_fit_classifier('counts', 'meta', target_column := 'body_site', name := 'site');
INSERT INTO models SELECT * FROM sc_fit_regressor('counts', 'meta', target_column := 'age_months', name := 'age');
```

| column | |
|---|---|
| `name` | required at fit; how every later function selects the model |
| `model` | the algorithm, currently always `random_forest` |
| `model_blob` | the serialized forest |
| `n_samples`, `n_features`, `n_trees` | what it was trained on |
| `task` | `classification` or `regression`; decides the type of a prediction column |
| `random_state` | echoed so a fit can be reproduced without remembering the seed |
| `feature_id_type`, `target_type` | the types this model was fit from, so later calls can hand ids and labels back unchanged |

`name` is mandatory because it removes a state from the system: a model is never unnamed, so `name :=` always works and nobody has to reason about selecting a model by position.

### Fitting

```sql
SELECT * FROM sc_fit_classifier('counts', 'meta', name := 'site', n_estimators := 500, random_state := 42);
```

| parameter | default | |
|---|---|---|
| `name` | *required* | the model's handle |
| `target_column` | auto | required when the metadata has more than one candidate column |
| `model` | `random_forest` | the algorithm; any other value is an error, and the value is echoed into the model table's `model` column |
| `n_estimators` | 100 | trees |
| `random_state` | 0 | seed |
| `n_threads` | 0 | 0 = sc's global pool, one thread per core (see future work 4) |
| `max_depth` | unbounded | `<= 0` is sklearn's `None` |
| `criterion` | `gini` / `squared_error` | classifier / regressor |
| `max_features` | `sqrt` / `all` | how many features each split may consider: `'sqrt'`, `'log2'`, `'all'`, a fraction in (0, 1], or an integer count |
| `min_samples_split` | 2 | integer count or fraction |
| `min_samples_leaf` | 1 | integer count or fraction |
| `min_weight_fraction_leaf`, `min_impurity_decrease` | 0.0 | |
| `bootstrap` | true | |
| `max_samples` | all | only with `bootstrap := true` |

The integer-versus-float distinction is load-bearing and follows sklearn: `max_features := 10` means ten features, `max_features := 0.3` means thirty percent of them.

`max_features` deserves care, because it is not a filter on the table. At **every split of every tree**, the tree draws that many *candidate* features uniformly at random — a partial Fisher–Yates shuffle over all feature indices — and picks the best split among only those. It is not the most abundant or most prevalent features, and it is not one subset chosen once: a fresh draw happens at every node, so over a forest every feature gets many chances. That randomness is what decorrelates the trees, which is where a forest's advantage over a single tree comes from. `max_features := 'all'` evaluates every feature at every split, consumes no randomness at all, and yields correlated trees that fit more slowly and usually score no better.

### Predicting

One function for both tasks; the `prediction` column takes the model's own label type for a classifier and DOUBLE for a regressor, resolved from the model's `task` when the query is bound.

```sql
SELECT * FROM sc_predict('new_counts', 'models', name := 'site');
-- sample_id, prediction, sample_coverage
```

For class probabilities, one row per sample and class:

```sql
SELECT * FROM sc_predict_proba('new_counts', 'models', name := 'site');
-- sample_id, class, probability, sample_coverage
```

`sc_predict_proba` rejects a regressor at **bind time** — while DuckDB plans the query, before a single row is read — rather than partway through execution. Binding is the same phase that fixes each column's type, which is how `sc_predict`'s `prediction` column can be VARCHAR for one model and DOUBLE for another: both the schema and this error come from reading the model's `task` during planning. In practice it means a wrong model name or a wrong task fails instantly, without scanning anything. Note that a forest's probabilities are vote fractions, not calibrated posteriors: a 0.9 does not mean nine of ten such samples belong to that class.

### Cross-validation

How well does this recipe generalize to samples it has never seen? The samples are shuffled once and cut into `n_folds` folds — stratified for a classifier, so no fold can miss a class, plain for a regressor. Each fold is held out in turn while a forest is fit on the rest and used to predict it, so every sample ends up with exactly one prediction from a model that never saw it.

```sql
CREATE TABLE cv AS SELECT * FROM sc_cross_validate_classifier('counts', 'meta', target_column := 'body_site', n_folds := 5, random_state := 42);
SELECT DISTINCT metric, mean_score, std_score, fold_scores FROM cv;
```

| column | |
|---|---|
| `sample_id`, `prediction`, `actual` | the out-of-fold prediction and the target it was scored against |
| `probabilities`, `classes` | classifier only; lists aligned position by position |
| `metric` | `accuracy` (classifier) or `mse` (regressor) |
| `mean_score`, `std_score` | over the per-fold scores |
| `fold_scores` | each fold's score |

**`mean_score` is the generalization estimate** — cross-validated accuracy for a classifier, cross-validated MSE for a regressor, where `sqrt(mean_score)` is an RMSE in the target's own units. `std_score` says how much that estimate moved between folds, which is the part a single number hides. The training error of a model fit on all the data is not comparable and will look far better; that gap is the reason for the held-out folds.

Some details worth knowing:

- **The per-fold forests are scored and discarded.** This measures the recipe, not a model. Fit the model you will actually use with `sc_fit` on all the data.
- **`mean_score` averages folds, not samples.** When the sample count does not divide evenly, folds carry unequal weight and the number differs slightly from scoring the whole `prediction` column yourself. This is sklearn's convention.
- **`std_score` is the population standard deviation** (ddof=0) of the fold scores, matching `np.std` and q2's printed value.
- **MSE is an error, so lower is better** — the opposite direction from the classifier's accuracy in the same column, which is why `metric` is part of the output.
- **The last four columns describe the run, not the sample**, so the same four values repeat on every row. There are two grains in one result — per-sample predictions, and one experiment-level summary — and a table function can return only one relation. The alternative, a second function returning just the summary, would have to run the entire cross-validation again: k more forest fits to recompute numbers already in hand. So the summary rides along, and `SELECT DISTINCT metric, mean_score, std_score, fold_scores` recovers the single row it really is. If sc ever exports which fold each sample was in (future work 3), the summary becomes a `GROUP BY fold` over the per-sample rows and these columns can go.
- **`actual` is carried through** so the result scores itself — a confusion matrix or a residual needs no join back to the metadata.
- **`n_folds` is spelled `cv`** by q2, scikit-learn's `cross_val_score` and sc's own C struct. Both names work; passing both is an error.
- **`parameter_tuning := true`** runs a randomized search inside each training fold and refits with the winner. The search never sees the held-out fold, which is what keeps the score honest. It costs roughly `n_folds × 20` extra fits, and see future work 1 before trusting the result.  The parameters that give the best score are not currently reported.

```sql
-- Confusion matrix, no join required.
SELECT actual, prediction, count(*) FROM cv GROUP BY ALL ORDER BY 1, 2;
-- The least confident calls.
SELECT sample_id, prediction, actual, list_max(probabilities) AS confidence FROM cv ORDER BY confidence LIMIT 10;
```

### Explaining a prediction

TreeSHAP attributes a single prediction to individual features: how much each one pushed this sample away from the model's baseline. It is additive by construction — `base_value + sum(shap_value)` over a sample's features equals the prediction, which is the class probability for a classifier and the predicted value for a regressor.

```sql
SELECT * FROM sc_shap('counts', 'models', name := 'site', top_k := 10);
-- sample_id, class, feature_id, shap_value, base_value, sample_coverage
```

| parameter | default | |
|---|---|---|
| `predicted_class_only` | true | with two classes the other is an exact mirror, so nothing is lost; with three or more, `false` shows why not the runner-up for the given sample |
| `top_k` | off | keep k features per sample: the strongest pushes toward the prediction and against it, split evenly, the positive side taking any odd one out, either side filling in when the other runs short |
| `max_attributions` | 10,000,000 | how many attributions one batch may hold (~80 MB). 8 bytes per attribution.  attributions = samples * classes * features |
| `batch_size` | derived | samples per batch, overriding the budget; passing both is an error. batch_size = max_attributions / (classes * features) |

- **`base_value` is the model's average output over the training data**, before any feature is known. The same on every row, since it depends on the trees alone.
- **Absent features get attributions too.** A feature the model knows but the sample lacks is fed to the trees as zero, and a zero can push a prediction as hard as a large count. On ECAM, one sample's explanation was 70% carried by features it did not contain, which in microbiome terms is real signal: the absence of adult taxa is evidence of an infant.
- **`top_k` does not reduce computation.** Everything must be computed before the strongest can be picked; `top_k` keeps the *output* small, which is what matters when the result is stored (`CREATE TABLE AS` costs roughly 72 bytes a row before compression).
- **Rows come out in waterfall order** — `sample_id`, then `class`, then `shap_value` descending — which a plain `SELECT` or `CREATE TABLE AS` preserves and an outer `ORDER BY` simply re-sorts.
- **Memory is bounded by the batch, not the run.** Samples are explained in batches and rows stream out one batch at a time, so explaining 20,000 samples batched costs the same as explaining 2,000. Every attribution is one `double`, so the 10M-attribution default budget is 10M × 8 bytes = 80 MB of attributions per batch.

### Feature importances

```sql
SELECT * FROM sc_feature_importances('models', name := 'site') ORDER BY importance DESC LIMIT 20;
-- feature_id, importance
```

Mean decrease in impurity, averaged over the trees and summing to 1. It is a property of the *model*, not of any prediction, and it is biased toward features with many distinct values. For "why this sample", use `sc_shap`; for a global feature ranking that does not share that bias, the mean absolute SHAP value over all samples for a given feature is the better measure.

```sql
SELECT * FROM sc_model_features('models', name := 'site');
-- feature_id, column_index
```

`sc_model_features` lists the training vocabulary in the model's own column order — useful when a prediction table shares few features with the model and you want to see exactly what the model expects.

### Identifier and label types

Ids are matched **as text** — that is how a BIGINT `42` and a VARCHAR `'42'` name one feature — but they are returned as the type they arrived as, so a result joins and sorts like its source:

| column | takes its type from |
|---|---|
| `sample_id` | the data relation of that call |
| `feature_id` | the relation the model was fit from |
| `prediction`, `class`, `classes`, `actual` (classifier) | the metadata column the model was fit from |
| `prediction`, `actual` (regressor) | DOUBLE — a forest predicts a continuous value whatever the column was |

`feature_id` follows the model rather than the call because a sample is explained against features it does not contain, so those ids are the model's, not the data's. The fit records both types on the model row (`feature_id_type`, `target_type`); a model row written without them reports VARCHAR, which is what these functions did before.

Two consequences worth knowing:

- **The type must be an id type.** `sample_id` and `feature_id` must be VARCHAR, BIGINT or UUID. An INTEGER or DATE id is rejected at bind with the cast to write; targets carry no such restriction, since any type that renders to text can label a sample.
- **Matching is still textual, so spelling matters.** `'042'`, `'42.0'`, a trailing space, or an uppercase UUID held in a VARCHAR column are all different features from `42`. When nothing matches, the error retries the lookup under trimming, lowercasing and canonical digits, and names the cast that would fix it.

### Unknown and absent features

Prediction data rarely matches training data exactly, and the two mismatches are handled differently:

- **A feature the model never saw is dropped**, with a warning. No tree splits on it, so it could not affect the prediction anyway.
- **A feature the model knows but the sample lacks is zero**, which is what absence means in a sparse table, and the trees split on it normally.

`sample_coverage` reports the share of a sample's own observed cells whose feature the model knows. A low value means the model is being applied to different data — another reference database, another pipeline, another 16S region — and a 0.0 means the sample kept nothing at all and is being predicted from an all-zero row. The predictions in that case are perfectly well-formed and completely uninformed, which is exactly why the column exists.

### Reproducibility

`random_state` fixes everything sc does: the bootstrap draws, the feature subsampling at each split, the cross-validation folds, and the tuning search. Results are also invariant to `n_threads` — the same seed gives bit-identical output at one thread or eight, which sc tests directly.

One thing outside the seed's control matters more than it looks: both dictionaries are **sorted** before encoding. Feature ids are sorted so that column *i* means the same feature in every table, and sample ids are sorted so that a parallel scan delivering chunks in a different order cannot change which samples the bootstrap draws. At prediction time the model's own vocabulary order is used and never re-sorted, since re-deriving an encoding from prediction data would silently point every learned split at the wrong feature.

To check reproducibility on your own data — two fits with one seed produce the identical serialized forest, whatever the thread count:

```sql
CREATE TABLE a AS SELECT * FROM sc_fit_regressor('counts', 'meta', name := 'a', random_state := 7);
CREATE TABLE b AS SELECT * FROM sc_fit_regressor('counts', 'meta', name := 'b', random_state := 7, n_threads := 8);
SELECT (SELECT md5(model_blob) FROM a) = (SELECT md5(model_blob) FROM b) AS identical;
```

### Scale and memory

The lever that matters most is not a parameter here but the shape of the table. Microbiome feature tables are long-tailed, with most features appearing in a handful of samples and contributing little, so a prevalence filter cuts the feature count hard for little cost — and everything below scales with that count:

```sql
CREATE VIEW common AS
  SELECT c.* FROM counts c
  SEMI JOIN (SELECT feature_id FROM counts GROUP BY 1
             HAVING count(DISTINCT sample_id) >= 0.01 * (SELECT count(DISTINCT sample_id) FROM counts)) USING (feature_id);
```

What to expect:

- **Fitting** spends most of its time moving the table into sc rather than growing trees, so the cell count matters more than `n_estimators`. Memory grows with the thread count, since each thread builds a tree at once; `n_threads := 1` is the way to trade time for memory.
- **`sc_shap`** holds one `double` per (sample, class, feature) in the batch it is working on, so memory stays flat as the sample count grows — `max_attributions` (80 MB by default) is the ceiling, not the size of the run. Narrow the feature count first if that is too much: attributions scale with it directly.
- **wasm** is the tightest environment: DuckDB and this extension share a single 4 GiB memory and run single-threaded, so expect a far lower ceiling on table size than the same query takes natively.

### Citations

- Breiman, L. (2001). Random Forests. *Machine Learning* 45(1), 5–32.
- Lundberg, S. M., Erion, G. G., & Lee, S.-I. (2018). Consistent Individualized Feature Attribution for Tree Ensembles. [arXiv:1802.03888](https://arxiv.org/abs/1802.03888) — the path-dependent TreeSHAP algorithm `sc_shap` implements. Cite this if you publish SHAP values.
- Bokulich, N. A., Dillon, M. R., Bolyen, E., Kaehler, B. D., Huttley, G. A., & Caporaso, J. G. (2018). q2-sample-classifier: machine-learning tools for microbiome classification and regression. *Journal of Open Source Software* 3(30), 934.
- Pedregosa, F. et al. (2011). Scikit-learn: Machine Learning in Python. *JMLR* 12, 2825–2830 — the reference implementation sc's defaults and oracle fixtures follow.
