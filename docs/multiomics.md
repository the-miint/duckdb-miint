# Multi-omics integration (MMvec)

Learn which **metabolites co-occur with which microbes**, from paired count tables of the same
samples.

MMvec (Morton et al. 2019) fits a neural-network-free probabilistic model of the conditional
probability of each metabolite given each microbe, and reports those conditionals as
co-occurrence scores. Unlike correlation, it is not confounded by the compositional nature of
sequencing data — the reason correlating two count tables directly produces so many spurious
associations.

> The examples on this page are executed as tests. Every statement and every number below is
> mirrored in `test/sql/mmvec_docs.test`; if you change one, change the other.

## Table of Contents

- [Feature tables](#feature-tables) - input expectations for both modalities.
- [The model](#the-model) - what is actually being fitted, in one paragraph.
- [`mmvec_fit`](#mmvec_fit) - fit a model.
- [The model relation](#the-model-relation) - what `mmvec_fit` returns and how to read it.
- [`mmvec_ranks`](#mmvec_ranks) - the headline output: co-occurrence scores.
- [`mmvec_predict`](#mmvec_predict) - predicted metabolite profiles for new samples.
- [`mmvec_score`](#mmvec_score) - Q² against a held-out set.
- [`mmvec_train_test_split`](#mmvec_train_test_split) - reproducible sample-wise splits.
- [Worked example: cross-validation](#worked-example-cross-validation) - split → fit → filter → predict → score.
- [Things that surprise people](#things-that-surprise-people) - read this before filing a bug.
- [Reproducibility](#reproducibility) - what is guaranteed, and what is not.
- [Differences from scikit-bio](#differences-from-scikit-bio) - a complete list.
- [Citations](#citations)

## Feature tables

Two long-form `(sample_id, feature_id, value)` relations — the same schema as
[`read_biom`](reading.md#biom) and everything in [diversity](diversity.md). One holds the
**X** modality (canonically microbes), the other **Y** (canonically metabolites).

Requirements, all enforced with errors that name the offending id:

- Both tables must describe **exactly the same samples**, matched by `sample_id`. Not by row
  order — long-form tables have none.
- Y needs **at least 2 features** (see [the reference category](#the-reference-category)).
- Values must be non-negative and finite. MMvec is a multinomial model on counts; a
  variance-stabilized or CLR-transformed table is a modelling error here, not an input to
  silently fit.
- No duplicate `(sample_id, feature_id)` cells. They are rejected rather than summed, because
  summing them would change the fit.
- No all-zero rows or columns. An all-zero feature contributes nothing to the likelihood, so
  its parameters would be driven entirely by the prior and the fit would "succeed" while that
  feature's output was meaningless.

`feature_id` and `sample_id` may be `VARCHAR`, `BIGINT` or `UUID`; the type is mirrored onto the
output so results join back to your metadata without a cast. The two modalities may use
different id types.

```sql
CREATE TABLE microbes AS SELECT * FROM (VALUES
  ('s1','m1',12),('s1','m2', 3),('s1','m3', 1),('s1','m4', 4),
  ('s2','m1', 9),('s2','m2', 7),('s2','m3', 2),('s2','m4', 2),
  ('s3','m1', 2),('s3','m2',11),('s3','m3', 6),('s3','m4', 1),
  ('s4','m1', 1),('s4','m2', 8),('s4','m3',10),('s4','m4', 1),('s4','m5', 5),
  ('s5','m1', 6),('s5','m2', 2),('s5','m3', 3),('s5','m4', 9),('s5','m5', 2),
  ('s6','m1', 3),('s6','m2', 4),('s6','m3', 8),('s6','m4', 5)
) AS t(sample_id, feature_id, value);

CREATE TABLE metabolites AS SELECT * FROM (VALUES
  ('s1','b1',300),('s1','b2', 40),('s1','b3', 10),
  ('s2','b1',250),('s2','b2', 80),('s2','b3', 20),
  ('s3','b1', 60),('s3','b2',260),('s3','b3', 30),
  ('s4','b1', 30),('s4','b2',300),('s4','b3', 70),
  ('s5','b1',180),('s5','b2', 50),('s5','b3',170),
  ('s6','b1', 90),('s6','b2',120),('s6','b3',190)
) AS t(sample_id, feature_id, value);
```

## The model

For each microbe *i* and metabolite *j*, MMvec models

```
P(Y_j | X_i) = softmax_j( x_main[i] · y_main[:,j] + x_bias[i] + y_bias[j] )
```

where `x_main` and `y_main` are latent embeddings of dimension `dimensions`. Fitting maximizes
the posterior under Gaussian priors on every parameter. **The sample axis is summed away before
the first iteration** — the entire contribution of the data is the cross-product `Xᵀ·Y` and its
row sums — so cost per iteration depends on the number of *features*, not the number of samples.

## `mmvec_fit`

```sql
CREATE TABLE fit AS
SELECT * FROM mmvec_fit('microbes', 'metabolites', dimensions := 1, max_iter := 200, seed := 42);
```

**Parameters:**

- `x_table` (VARCHAR): name of the X feature table
- `y_table` (VARCHAR): name of the Y feature table
- `dimensions` (INTEGER, default 3): latent dimension
- `optimizer` (VARCHAR, default `'lbfgs'`): `'lbfgs'` or `'adam'`
- `max_iter` (BIGINT, default 1000): L-BFGS **iterations**, or Adam **epochs** — one parameter
  for both, matching scikit-bio, which is why there is no separate `epochs`
- `seed` (BIGINT, default 0): seeds the parameter initialization
- `x_prior_mean`, `x_prior_scale`, `y_prior_mean`, `y_prior_scale` (DOUBLE, defaults 0 / 1 / 0 / 1):
  the Gaussian priors
- Adam only, ignored by L-BFGS: `learning_rate` (default 0.001), `batch_size` (default 50),
  `beta_1` (default 0.9), `beta_2` (default **0.95**, not the 0.999 usual elsewhere),
  `clipnorm` (default 10.0), `batch_norm` (`'unbiased'` default, or `'legacy'`)

There is deliberately **no `threads` parameter**: the fit is single-threaded and pins Eigen to
one thread, which is what makes a seeded fit bit-reproducible.

**`batch_norm` is not a step size.** `'legacy'` reproduces upstream mmvec by scaling the
likelihood by `n_samples / sum(X)` while leaving the priors untouched — a differently
regularized model, typically far more heavily regularized, not a rescaled version of the same
one. `'unbiased'` makes the minibatch loss an unbiased estimator of the full-batch objective, so
Adam and L-BFGS target the same posterior. Use `'legacy'` only to reproduce published mmvec
numbers.

## The model relation

`mmvec_fit` returns one long-form relation carrying the parameters *and* the fit diagnostics:

| column | type | meaning |
|---|---|---|
| `modality` | VARCHAR | `'x'`, `'y'`, or `'loss'` |
| `x_feature_id` | mirrors X's id type | set when `modality = 'x'`, else NULL |
| `y_feature_id` | mirrors Y's id type | set when `modality = 'y'`, else NULL |
| `axis` | INTEGER | `0` is the bias, `1..dimensions` are embedding coordinates; for `'loss'`, the 1-based evaluation ordinal |
| `value` | DOUBLE | the coordinate, bias, or loss |
| `converged` | BOOLEAN | broadcast onto every row |
| `n_iter` | BIGINT | broadcast |
| `final_loss` | DOUBLE | broadcast |
| `max_abs_grad` | DOUBLE | broadcast |
| `message` | VARCHAR | broadcast; human-readable outcome |

`axis = 0` being the bias matches scikit-bio's `x_embeddings = hstack([x_main, x_bias])` column
order.

```sql
SELECT modality, count(*) AS n FROM fit GROUP BY modality ORDER BY modality;
```

```
loss    62
x       10
y        6
```

Five microbes × (1 bias + 1 coordinate) = 10; three metabolites × 2 = 6; and 62 loss rows, one
per objective *evaluation*.

The five diagnostic columns are broadcast onto every row rather than stored as extra rows, so
each keeps its natural type and needs no pivot:

```sql
SELECT DISTINCT converged, n_iter, round(final_loss, 3) AS final_loss FROM fit;
```

```
true    57    47967.431
```

## `mmvec_ranks`

The headline output — one row per (microbe, metabolite) pair.

```sql
SELECT x_feature_id, y_feature_id, round(rank, 4) AS rank, round(prob, 4) AS prob
FROM mmvec_ranks('fit')
ORDER BY x_feature_id, y_feature_id;
```

**Output schema:** `(x_feature_id, y_feature_id, rank DOUBLE, prob DOUBLE)`, `d1 × d2` rows.

- **`rank`** is the log-conditional-probability, **row-centred**: each microbe's ranks sum to
  zero. Centring is what makes one microbe's scores comparable to another's, and what removes
  the arbitrariness of which metabolite was the reference category. Positive means "more
  associated with this microbe than that microbe's average".
- **`prob`** is `P(Y_j | X_i)` itself; each microbe's probabilities sum to 1.

Ranks are the quantity to interpret and to compare across runs. The **embeddings are not** —
they are identified only up to an orthogonal rotation, so two runs that agree perfectly on every
rank can have completely unrelated `x_main` / `y_main`. That is why the model relation exposes
them but nothing downstream compares them.

The usual question — *what does each microbe co-occur with most?*

```sql
SELECT x_feature_id, y_feature_id, round(rank, 4) AS rank
FROM (SELECT *, row_number() OVER (PARTITION BY x_feature_id ORDER BY rank DESC) AS rn
      FROM mmvec_ranks('fit'))
WHERE rn = 1
ORDER BY x_feature_id;
```

```
m1  b1  0.7666
m2  b2  0.5093
m3  b2  0.4845
m4  b1  0.3551
m5  b2  0.6610
```

Which is the structure built into the fixture: `m1` is abundant in the samples richest in `b1`,
`m2` and `m3` in those richest in `b2`.

## `mmvec_predict`

Expected metabolite proportions for samples the model has not seen.

```sql
SELECT sample_id, y_feature_id, round(proportion, 4) AS proportion
FROM mmvec_predict('model', 'x_eval')
ORDER BY sample_id, y_feature_id;
```

**Output schema:** `(sample_id, y_feature_id, proportion DOUBLE)`. Each sample's proportions sum
to 1.

The prediction is the sample's microbe **proportions** times the conditional probabilities, so
it is invariant to sequencing depth: multiplying a sample's counts by any constant leaves its
prediction bit-for-bit unchanged. Samples in the X table need not appear in the model — the
model constrains which *features* may appear, never which samples.

## `mmvec_score`

Q², the fraction of variance in the held-out metabolite proportions that the prediction
explains.

```sql
SELECT round(q_squared, 6) AS q_squared FROM mmvec_score('model', 'x_eval', 'y_eval');
```

**Output schema:** `(q_squared DOUBLE)`, a single row.

```
SS_res = Σ (y_props − predicted)²
SS_tot = Σ (y_props − colmean(y_props))²
Q²     = 1 − SS_res / SS_tot
```

`colmean` is per-metabolite across samples, so the baseline being beaten is "the mean community",
not a global constant. **Q² is unbounded below** — unlike R², it is evaluated on data the model
did not fit, so a badly fitting model scores arbitrarily negative. A value near 0 means "no
better than predicting the average community".

## `mmvec_train_test_split`

```sql
CREATE TABLE split AS SELECT * FROM mmvec_train_test_split('microbes', 0.34, 7);
```

**Parameters:** `(relation, test_fraction, seed)`. **Output schema:** `(sample_id, split)`, one
row per distinct sample, `split` being `'train'` or `'test'`. `sample_id` passes through with its
own type so it joins back without a cast.

One relation, not two: `mmvec_fit` requires X and Y to describe exactly the same samples and
validates that itself, so a second relation would carry no information — split either one and
filter both by the result.

Exactly `round(n × test_fraction)` samples are assigned to test. The assignment is a
deterministic function of `(sample_id, seed)` alone — samples are ordered by
`md5(seed || ':' || sample_id)` — so it does not depend on row order or on how the scan was
parallelized, and re-running it in another session, or on a permuted table, gives the same split.
md5 rather than DuckDB's `hash()` deliberately: `hash()` is an implementation detail free to
change between versions, whereas md5 is specified, so a split recorded in a paper still
reproduces after an upgrade.

`test_fraction` outside `[0, 1]`, or a NULL `test_fraction` or `seed`, is an error rather than a
silent all-train or all-test split.

## Worked example: cross-validation

The whole point of `predict` and `score` — and the one place people get stuck.

```sql
CREATE TABLE split AS SELECT * FROM mmvec_train_test_split('microbes', 0.34, 7);

CREATE TABLE x_train AS SELECT m.* FROM microbes m JOIN split s USING (sample_id) WHERE s.split = 'train';
CREATE TABLE y_train AS SELECT m.* FROM metabolites m JOIN split s USING (sample_id) WHERE s.split = 'train';
CREATE TABLE x_test  AS SELECT m.* FROM microbes m JOIN split s USING (sample_id) WHERE s.split = 'test';
CREATE TABLE y_test  AS SELECT m.* FROM metabolites m JOIN split s USING (sample_id) WHERE s.split = 'test';

CREATE TABLE model AS
SELECT * FROM mmvec_fit('x_train', 'y_train', dimensions := 1, max_iter := 200, seed := 42);
```

**Now the step that is not optional.** A sample-wise split routinely leaves *test-only features*
behind — here `m5` appears only in `s4` and `s5`, both of which landed in the test set, so the
model was never fitted on it:

```sql
SELECT (SELECT count(DISTINCT x_feature_id) FROM model WHERE modality = 'x') AS model_microbes,
       (SELECT count(DISTINCT feature_id) FROM x_test) AS test_microbes;
```

```
4    5
```

Predicting without restricting the test table is an error, by design — there is no conditional
probability for a feature the model has never seen, and silently dropping it would quietly change
the sample's proportions:

```
Invalid Input Error: mmvec: the X table has 1 feature(s) the model was never fitted on (m5);
there is no conditional probability for a feature the model has not seen, so restrict the table
to the model's own features first -- e.g. WHERE feature_id IN (SELECT DISTINCT x_feature_id
FROM <model> WHERE modality = 'x') -- or refit including them
```

So filter first:

```sql
CREATE TABLE x_eval AS SELECT * FROM x_test
  WHERE feature_id IN (SELECT DISTINCT x_feature_id FROM model WHERE modality = 'x');
CREATE TABLE y_eval AS SELECT * FROM y_test
  WHERE feature_id IN (SELECT DISTINCT y_feature_id FROM model WHERE modality = 'y');

SELECT round(q_squared, 6) AS q_squared FROM mmvec_score('model', 'x_eval', 'y_eval');
```

```
0.037772
```

The reverse case is fine and needs no handling: a model feature *missing* from the test table is
accepted and simply contributes no weight.

## Things that surprise people

### The reference category

The **lexicographically first Y feature has all-zero parameters, by construction.** Its row in
the model relation is not missing data:

```sql
SELECT y_feature_id, axis, value FROM fit WHERE modality = 'y' AND y_feature_id = 'b1' ORDER BY axis;
```

```
b1  0  0.0
b1  1  0.0
```

A softmax over *d* categories has only *d−1* free sets of parameters, so one category is pinned
at zero to identify the rest. `mmvec_ranks` still reports a real, non-zero rank for it — the
row-centring puts it back on the same footing as the others.

**This choice is not statistically inert.** The Gaussian priors are not shift-invariant, so which
feature holds the reference slot genuinely changes the fitted optimum — measured to move ranks by
up to 31% of their magnitude. It is deterministic and documented rather than arbitrary, but if
you rename your features you will get a different (though equally valid) fit.

### BIGINT and UUID ids sort by their text form

Feature order is lexicographic on the id **as text**, whatever the column type. So a table keyed
`9, 10` puts `10` first, and the reference category is `10`, not `9`:

```sql
CREATE TABLE bx AS SELECT sample_id, CAST(f AS BIGINT) AS feature_id, value FROM (VALUES
  ('s1',9,5),('s1',10,3),('s2',9,2),('s2',10,7),('s3',9,4),('s3',10,6)) t(sample_id,f,value);
CREATE TABLE bmet AS SELECT sample_id, CAST(f AS BIGINT) AS feature_id, value FROM (VALUES
  ('s1',9,50),('s1',10,20),('s2',9,10),('s2',10,60),('s3',9,30),('s3',10,40)) t(sample_id,f,value);

CREATE TABLE bfit AS SELECT * FROM mmvec_fit('bx', 'bmet', dimensions := 1, max_iter := 50);
SELECT y_feature_id, axis, value FROM bfit WHERE modality = 'y' AND value = 0.0 ORDER BY axis;
```

```
10  0  0.0
10  1  0.0
```

Zero-pad numeric ids if you want numeric order to be the reference order.

### Two of the id columns are always NULL

Every row of the model relation has exactly one of `x_feature_id` / `y_feature_id` populated, and
`'loss'` rows have neither. That is the price of keeping both id types exact in one relation; the
alternative, a single VARCHAR id column, would lose the type the moment you inspected the model.

The practical consequence: **`JOIN ... USING (x_feature_id, y_feature_id)` against the model
relation silently returns nothing**, because NULL never matches NULL. Use `IS NOT DISTINCT FROM`,
or filter by `modality` first. (Output of `mmvec_ranks` has both columns populated on every row,
so ordinary joins work there.)

### `converged = false` is the normal outcome on real data

At the `max_iter` scikit-bio uses for them, real datasets routinely stop on the iteration limit
rather than on a convergence test. The published soils dataset needs roughly 2430 iterations to
actually reach a stationary point, against the 500 scikit-bio runs it for. That is not a failure,
and neither scikit-bio nor MIINT raises — but MIINT surfaces it explicitly rather than leaving you
to assume otherwise.

**`converged` means "the optimizer stopped on its own terms", not "the gradient is zero".** The
stopping rule tests the relative change in the objective, not the gradient, and SciPy reports
success on that basis too. The two genuinely disagree: a fit can report `converged = true` with
`max_abs_grad` in the thousands. If you need a stationary point, read `max_abs_grad`, which is a
separate column for exactly this reason.

Adam always reports `converged = false`: it has no convergence test at all and simply runs its
epoch schedule to the end. `message` says so.

### `n_iter` is not scikit-bio's `n_iter_`

`n_iter` here is genuine optimizer iterations. scikit-bio's `n_iter_` is the length of its loss
curve, which counts objective **evaluations** including line-search probes — a larger number
(130 for `max_iter=200` on their fixtures). For the scikit-bio quantity:

```sql
SELECT count(*) FROM fit WHERE modality = 'loss';
```

After a line-search failure `n_iter` is `0`, and that zero means *unknown*, not *none*: the
solver reports its count only on a normal return. `message` says so explicitly on that path.

### Counts, not proportions — but depth does not matter

Pass raw counts. Both tables are internally read as proportions, so sequencing depth cancels:
multiplying one sample's counts by ten leaves its prediction bit-for-bit identical. What *does*
matter is the relative depth *between* features within a sample, which is the signal.

## Reproducibility

**On a given build, a given seed reproduces bit-for-bit.** The fit is single-threaded, Eigen is
pinned to one thread, and the RNG is a hand-written transform over `std::mt19937_64` rather than
`<random>`'s distributions — which are implementation-defined and would otherwise differ between
standard libraries.

**Across build targets, agreement is about 1e-9, not bit-identity.** The RNG stream *is*
bit-identical everywhere, and a fixed-parameter forward pass (`ranks`, `probs`, `predict`,
`score` on an existing model) is bit-identical too. But an *iterative fit* diverges, because the
WebAssembly build has no `-msimd128` and so uses scalar Eigen where a native build uses 2-wide
SSE2 — a different reduction order, amplified over hundreds of iterations. Compare fits across
targets with a tolerance, not with equality.

Since the parameters are only identified up to a rotation, compare **ranks**, not embeddings:

```sql
WITH a AS (SELECT * FROM mmvec_ranks('fit_a')), b AS (SELECT * FROM mmvec_ranks('fit_b'))
SELECT max(abs(a.rank - b.rank))
FROM a JOIN b USING (x_feature_id, y_feature_id);
```

## Differences from scikit-bio

MIINT reimplements scikit-bio 0.7.3's `mmvec` in C++ to reproduce its results. No scikit-bio code
is vendored, translated or executed. The deliberate differences:

- **Samples are joined by `sample_id`**, where scikit-bio checks only that the two tables have the
  same number of rows and pairs them positionally. Long-form SQL tables have no row order to pair
  on, and a mismatch here is a loud error naming the offending ids.
- **Feature order is lexicographic on the id**, which fixes the reference category deterministically.
- **`seed` defaults to 0**, not to system entropy: mmvec is an estimator whose reproducibility is
  the point.
- **Validation is stricter**: duplicate cells, all-zero rows and columns, negative and non-finite
  values, and a single-feature Y are all rejected rather than fitted.
- **Unknown features in a test table are an error**, not a silent drop.
- **`converged` and `max_abs_grad` are reported separately**, and `n_iter` counts iterations
  rather than evaluations, with the evaluation trace available as the `'loss'` rows.
- **Q²'s zero-variance guard is exact equality**, matching scikit-bio, so it fires only for a
  single sample. Several near-identical samples leave a tiny non-zero `SS_tot` and score hugely
  negative rather than returning 0.

## Citations

**MMvec**
Morton JT, Aksenov AA, Nothias LF, Foulds JR, Quinn RA, Badri MH, Swenson TL, Van Goethem MW,
Northen TR, Vazquez-Baeza Y, Wang M, Bokulich NA, Watters A, Song SJ, Bonneau R, Dorrestein PC,
Knight R. "Learning representations of microbe-metabolite interactions." *Nature Methods*
16(12):1306-1314, 2019. [doi:10.1038/s41592-019-0616-3](https://doi.org/10.1038/s41592-019-0616-3)

**scikit-bio**, whose `skbio.stats.ordination.mmvec` this reimplements and whose results are the
parity target. See `THIRD_PARTY_LICENSES.md`.

**L-BFGS**
Liu DC, Nocedal J. "On the limited memory BFGS method for large scale optimization."
*Mathematical Programming* 45:503-528, 1989. The solver is the vendored
[LBFGS++](https://github.com/yixuan/LBFGSpp).

**Adam**
Kingma DP, Ba J. "Adam: A Method for Stochastic Optimization." *ICLR*, 2015.
[arXiv:1412.6980](https://arxiv.org/abs/1412.6980)
