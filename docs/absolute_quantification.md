# Absolute quantification (synDNA spike-ins)

Shotgun metagenomic read counts are *compositional*: they tell you the relative make-up of
a sample, not how much of anything was actually there. Two samples can produce identical
profiles while differing a hundred-fold in absolute microbial load, and a taxon that holds
a constant absolute abundance will appear to fall whenever anything else rises.

The synDNA spike-in method of Zaramela et al. 2022 breaks that tie. A known mass of a
pool of synthetic DNA constructs is added to each sample before extraction. Because the
mass of each construct going in is known and its read count coming out is measured, the
constructs form a **per-sample standard curve** mapping read count to input mass — and
that curve then converts any feature's read count into an absolute mass.

This page covers fitting those standard curves.

## Table of Contents

- [Input relations](#input-relations) - the three relations every function here takes.
- [Sample identifier types](#sample-identifier-types) - how VARCHAR/BIGINT/UUID sample ids are handled.
- [`absquant_fit_models`](#absquant_fit_models) - fit the per-sample standard curve.
- [Which samples get a model](#which-samples-get-a-model) - what is dropped, and how to find out.
- [Applying a model](#applying-a-model) - turning a fitted curve into masses.
- [Differences from pysyndna](#differences-from-pysyndna) - for readers coming from the reference implementation.
- [Citations](#citations)

## Input relations

All three are passed **by name**, as quoted string literals — see
[Passing relations by name](table_of_contents.md#passing-relations-by-name). Tables,
views, `TEMP` tables and `TEMP` views all work.

### synDNA counts

A long-form `(sample_id, feature_id VARCHAR, value DOUBLE)` relation holding reads
assigned to each synDNA construct in each sample — the same schema as a
[feature table](diversity.md#feature-table), and the same schema
[`read_biom`](reading.md#biom) produces.

```sql
CREATE TABLE syndna_counts AS SELECT * FROM read_csv('data/syndna/fit_counts.csv');
-- sample_id,feature_id,value
-- A,p126,93135.0
-- A,p136,15190.0
-- ...
```

Rows with a NULL `sample_id`, `feature_id` or `value` are dropped, as are zero counts —
the sparse-storage invariant shared with the diversity functions. A zero read count
carries no information for a log-scale fit, so dropping it changes nothing.

Counts must be **finite and not negative**, and each `(sample_id, feature_id)` may appear
at most once. Both are errors rather than filters; see
[Differences from pysyndna](#differences-from-pysyndna).

### synDNA pool composition

`(feature_id VARCHAR, syndna_indiv_ng_ul DOUBLE)` — one row per construct in the pool,
giving its concentration. Concentrations set the *relative* mass of each construct, so
the units cancel: any consistent unit works, and scaling every value by a constant leaves
the fitted models unchanged to within floating-point rounding.

```sql
CREATE TABLE syndna_pool AS SELECT * FROM read_csv('data/syndna/fit_concentrations.csv');
-- feature_id,syndna_indiv_ng_ul
-- p126,1.0
-- p136,0.1
-- ...
```

Every value must be present, finite and not negative, `feature_id` must be unique, and
the values must sum to something positive and finite. This relation is configuration
rather than measurement, so a hole in it is reported, never propagated.

### Sample parameters

`(sample_id, mass_syndna_input_ng DOUBLE)` — the total mass of synDNA pool added to each
sample.

```sql
CREATE TABLE syndna_params AS SELECT * FROM read_csv('data/syndna/fit_params.csv');
-- sample_id,mass_syndna_input_ng
-- A,0.25
-- B,0.2
```

`sample_id` must be unique. A NULL, negative or infinite `mass_syndna_input_ng` marks
that sample as un-quantifiable and drops it, with a warning — it is a legitimate "no
measurement for this one", not malformed input.

## Sample identifier types

The counts relation's `sample_id` may be `VARCHAR`, `BIGINT` or `UUID`; any other type is
accepted and read as text. `absquant_fit_models` **mirrors that type onto its output**
(`BIGINT → BIGINT`, `UUID → UUID`, otherwise `VARCHAR`), so models join back to typed
metadata without a cast — the same convention as
[`unifrac_distances`](diversity.md#sample-identifier-types) and `align_minimap2`.

The reference relations are read as text and need not agree on a declared type: a
`BIGINT` `sample_id` in the counts relation matches a `VARCHAR` `'1'` in the parameters
relation, and the output is still `BIGINT`.

---

## `absquant_fit_models`

Fit, per sample, the linear model

```
log10(synDNA mass in ng)  ~  log10(read count)
```

by ordinary least squares, and return its coefficients. One row per sample that produced
a usable model.

```sql
CREATE TABLE models AS
SELECT * FROM absquant_fit_models(
    'syndna_counts', 'syndna_pool', 'syndna_params',
    1.0,
    min_syndna_counts := 50);
```

```
┌───────────┬────────────────────┬────────────────────┬────────────────────┬────────────────────────┬─────────────────────┬─────────────────────┐
│ sample_id │       slope        │     intercept      │       rvalue       │         pvalue         │       stderr        │  intercept_stderr   │
│  varchar  │       double       │       double       │       double       │         double         │       double        │       double        │
├───────────┼────────────────────┼────────────────────┼────────────────────┼────────────────────────┼─────────────────────┼─────────────────────┤
│ A         │ 1.2448765237913189 │ -7.355939160548428 │ 0.9865030975156573 │ 1.4284435606598037e-07 │ 0.07305408550335017 │   0.271274537363401 │
│ B         │ 1.2467591360440702 │ -7.450040830377357 │ 0.9863241797356327 │ 1.5053811468097085e-07 │ 0.07365795255302396 │ 0.27294119993260013 │
└───────────┴────────────────────┴────────────────────┴────────────────────┴────────────────────────┴─────────────────────┴─────────────────────┘
```

**Parameters:**
- `syndna_counts` (VARCHAR): name of the [synDNA counts](#syndna-counts) relation
- `syndna_pool` (VARCHAR): name of the [pool composition](#syndna-pool-composition) relation
- `sample_params` (VARCHAR): name of the [sample parameters](#sample-parameters) relation
- `syndna_contributing_fraction` (DOUBLE): fraction of the input pool mass that actually contributed reads — the mass entering the model is `mass_syndna_input_ng × syndna_contributing_fraction`. Must be `> 0` and `<= 1`. Use `1.0` when the whole spike-in is assumed to be sequenced.
- `min_syndna_counts` (BIGINT, default 1): drop any construct whose read count **summed over every sample** falls below this. Must be `>= 1`. The comparison is strict, so a construct landing exactly on the threshold is kept.

**Output schema:** `(sample_id, slope DOUBLE, intercept DOUBLE, rvalue DOUBLE, pvalue DOUBLE, stderr DOUBLE, intercept_stderr DOUBLE)` — the six fields `scipy.stats.linregress` returns, in its order. `sample_id` mirrors the input type (see [above](#sample-identifier-types)). Rows are **not** returned in any particular order; add `ORDER BY` if you need one.

**Reading the fit.** `slope` and `intercept` are the model; the other four describe how much to trust it. `rvalue` is Pearson's *r* — its square is the fraction of variance explained, and `r² < 0.8` is the conventional cut for discarding a sample's curve. `stderr` and `intercept_stderr` are the standard errors of the two coefficients, and `pvalue` is the two-sided test of a zero slope.

**How the mass of each construct is derived.** Each construct's share of the pool is its
concentration over the sum of **all** concentrations in the pool relation:

```
mass_pool(s)     = mass_syndna_input_ng(s) × syndna_contributing_fraction
fraction(i)      = syndna_indiv_ng_ul(i) / Σ syndna_indiv_ng_ul
syndna_ng(i, s)  = mass_pool(s) × fraction(i)
```

That denominator is the **full** pool relation and is never rescaled, including when
`min_syndna_counts` drops a construct: dropping a construct removes its data points from
the regression but does not change how much mass the remaining ones carry. This matters
— rescaling to the survivors shifts every fitted intercept.

Note the consequence for `syndna_contributing_fraction`: because it multiplies the mass
of every construct equally, it shifts the intercept by `log10` of itself and leaves the
slope untouched. Halving it is not the same as halving your answer.

## Which samples get a model

Samples disappear from the output for four different reasons, and **every one of them is
reported through [`miint_warnings()`](utilities.md)**. A short result set is never
silent:

| Reason | Warning |
|---|---|
| `mass_syndna_input_ng` is NULL, negative or not finite | `dropped N sample(s) whose mass_syndna_input_ng is NULL, negative or not finite: …` |
| No non-zero synDNA counts (the sample is absent from, or all-zero in, the counts relation) | `N sample(s) in '…' contributed no nonzero synDNA counts from '…': …` |
| Too few usable points, or the fit is degenerate (constant read counts, constant masses, fewer than two constructs left) | `no model could be fit for N sample(s): …` |
| — | constructs dropped by `min_syndna_counts` are reported too: `dropped N synDNA feature(s) with fewer than K total reads across all samples: …` |

```sql
CREATE TABLE models AS SELECT * FROM absquant_fit_models(…);
SELECT message FROM miint_warnings();
```

A sample present in the counts relation but **missing** from the parameters relation is
an error, not a warning — reads you cannot convert to mass are unusable. The reverse
(parameters for a sample you did not sequence) is merely untidy and is warned about.

## Applying a model

The fitted curve inverts to give the mass corresponding to any read count in that sample.
Below, `feature_counts` is **your own** feature table — the biological features you want
quantified, in the same long-form `(sample_id, feature_id, value)` shape as the synDNA
counts but a different relation:

```sql
-- ng of DNA attributable to each feature, per sample.
SELECT c.sample_id,
       c.feature_id,
       pow(10, m.slope * log10(c.value) + m.intercept) AS ng
FROM feature_counts c
JOIN models m USING (sample_id)
WHERE c.value > 0;
```

Two cautions. The model is fitted on `log10` of both axes, so it must be applied on that
scale and only within the range of read counts the constructs actually spanned —
extrapolating a power law beyond its calibration is how a spike-in method turns into
noise. And it is worth filtering on fit quality first, since a sample whose curve is poor
will still produce numbers:

```sql
SELECT * FROM models WHERE rvalue * rvalue >= 0.8;
```

## Differences from pysyndna

miint's models reproduce
[pysyndna](https://github.com/biocore/pysyndna)'s `fit_linear_regression_models` on its
own published fixtures to the full precision pysyndna emits — it serializes coefficients
truncated to 12 decimal places, and the parity tests hold to `1e-11` relative.

The behavioral differences below are deliberate. The first five change only *which* input
is accepted and how a rejection is reported, never a coefficient. The last one is a
numerical difference, and it is the one to know about if you diff the two
implementations.

- **Structurally impossible input is refused rather than propagated.** A negative or
  non-finite read count, a NULL/negative/non-finite concentration, a duplicate key in any
  of the three relations, or a pool summing to zero all raise. pysyndna lets several of
  these become NaN and silently void a whole sample.
- **An infinite `mass_syndna_input_ng` is treated as a bad parameter**, not as a value to
  compute with. pysyndna accepts it and then discards the sample anyway when the
  arithmetic produces NaN — the same outcome, reported as the wrong kind of problem.
- **A construct named by the pool but absent from the counts relation is accepted**, and
  is dropped by `min_syndna_counts` as the zero-count construct it is. pysyndna requires
  the two id sets to match exactly, which it can afford because its input is a dense
  table where a construct that sequenced nowhere is still present. In a long-form sparse
  relation "zero in every sample" and "absent" are the same input, so requiring an exact
  match would reject legitimate data. Because `min_syndna_counts` is at least 1, pysyndna
  drops such a construct too: the models are identical, only the error behavior differs.
  The other direction — a construct with reads but no concentration — is still an error,
  since its mass cannot be computed.
- **Samples with no usable model are omitted from the output**, with a warning. pysyndna
  represents this two different ways depending on how the fit failed.
- **`stderr` is computed from the residual sum of squares** rather than from scipy's
  `sqrt((1 − r²)·ssym/ssxm/df)`. The two are algebraically identical, but `1 − r²` is
  pure cancellation: it loses a digit for every decade the residuals shrink, and on a
  very clean standard curve it degrades badly (returning `0.0` for data with visible
  residuals at `r = 1`). If you compare against pysyndna on near-perfect data and see
  `stderr` disagree, this is why, and miint is the accurate side.

## Citations

Zaramela, L.S.; Tjuanta, M.; Moyne, O.; Neal, M.; and Zengler, K. (2022) "synDNA—a
Synthetic DNA Spike-in Method for Absolute Quantification of Shotgun Metagenomic
Sequencing", mSystems, 7(6), e00447-22. doi: 10.1128/msystems.00447-22

pysyndna — the reference implementation, by Amanda Birmingham (biocore). See
`THIRD_PARTY_LICENSES.md`.
