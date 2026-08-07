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

This page covers fitting those standard curves and applying them to get absolute cell
counts. It also covers a third function, [`absquant_orf_copies`](#absquant_orf_copies),
which asks the same kind of question of a **metatranscriptome** — how many copies of each
gene's RNA were present per gram of sample — and gets there without a spike-in at all.

## Table of Contents

- [Input relations](#input-relations) - the three relations the fit takes.
- [Sample identifier types](#sample-identifier-types) - how VARCHAR/BIGINT/UUID sample ids are handled.
- [`absquant_fit_models`](#absquant_fit_models) - fit the per-sample standard curve.
- [Which samples get a model](#which-samples-get-a-model) - what is dropped, and how to find out.
- [`absquant_cell_counts`](#absquant_cell_counts) - apply a curve to get absolute cell counts.
- [Choosing a metric](#choosing-a-metric) - per gram of gDNA, or per gram, microlitre or cm² of sample.
- [Which cells you get back](#which-cells-you-get-back) - the three filters, and how to see what they removed.
- [Composing the two functions](#composing-the-two-functions) - end to end, from spike-ins to cell counts.
- [Applying a model by hand](#applying-a-model-by-hand) - if you only want the mass.
- [`absquant_orf_copies`](#absquant_orf_copies) - copies of each ORF's ssRNA, per gram of sample.
- [The three RNA relations](#the-three-rna-relations) - what `absquant_orf_copies` reads.
- [ORF coordinates from `read_gff`](#orf-coordinates-from-read_gff) - and the off-by-one it will cost you.
- [Which ORF cells you get back](#which-orf-cells-you-get-back) - one filter, one report, and what is refused.
- [Differences from pysyndna](#differences-from-pysyndna) - for readers coming from the reference implementation.
- [Citations](#citations)

## Input relations

These three are what `absquant_fit_models` takes; `absquant_cell_counts` takes
[five of its own](#absquant_cell_counts). All are passed **by name**, as quoted string
literals — see
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

---

## `absquant_cell_counts`

Turn each feature's read count into **absolute cell counts**, by applying a fitted
standard curve. This is the payoff of the whole method: the output is an ordinary feature
table whose values are absolute rather than compositional, so it joins, pivots and writes
to BIOM like any other.

```sql
absquant_cell_counts(counts, models, coverage, lengths, params,
                     metric, min_coverage, min_rsquared := 0.8)
```

The chain, per surviving `(sample, feature)`:

```
gdna_mass_ng   = 10 ^ (slope * log10(read_count) + intercept)
genomes        = gdna_mass_ng * 6.02214076e23 / (ogu_len_in_bp * 650 * 1e9)
cells_per_g    = genomes / sequenced_sample_gdna_mass_ng * 1e9
```

`650` is the average molar mass of a base pair of double-stranded DNA, in g/mole. The
last step assumes **one genome per cell**, which is the published method's own
simplifying assumption and is wrong for dividing and polyploid organisms — treat the
output as genome equivalents if that matters to you.

### Choosing a metric

`cells_per_g` above is per gram of the gDNA that went into the **sequencer** — an
instrument quantity. The other three metrics restate it per unit of the material you
actually **collected**, which is what makes numbers comparable between samples:

```
extracted_gdna_mass_g = extracted_gdna_concentration_ng_ul * vol_extracted_elution_ul / 1e9
value                 = cells_per_g * extracted_gdna_mass_g / <the metric's denominator>
```

`extracted_gdna_mass_g` is the DNA recovered from the whole extraction, not the portion
sequenced — that is exactly what converts an instrument-side quantity into a sample-side
one.

| `metric` | Per unit of | Denominator column required in `params` |
|---|---|---|
| `cells_per_g_of_gdna` | gram of sequenced gDNA | — |
| `cells_per_g_of_sample` | gram of sample | `calc_mass_sample_aliquot_input_g` |
| `cells_per_ul_of_sample` | microlitre of sample | `sample_volume_ul` |
| `cells_per_cm2_of_sample` | square centimetre of sample | `sample_surface_area_cm2` |

Metric names are case-insensitive. Only the requested metric's denominator has to exist:
a study that recorded volumes but never surface areas can ask for
`cells_per_ul_of_sample` without carrying a `sample_surface_area_cm2` column.

Pick the metric that matches how the sample was taken — a swab has an area, a stool
aliquot a mass, a saliva sample a volume. `cells_per_g_of_gdna` is the right answer only
when you genuinely mean *per unit of DNA sequenced*, which is rarely the biological
question.

### The five relations

| Argument | Columns read | Notes |
|---|---|---|
| `counts` | `(sample_id, feature_id, value)` | your biological feature table, not the synDNA one |
| `models` | `(sample_id, slope, intercept, rvalue)` | `absquant_fit_models`' output, unprojected |
| `coverage` | `(sample_id, feature_id, coverage)` | a **fraction** in `[0, 1]`, never a percent |
| `lengths` | `(feature_id, ogu_len_in_bp)` | genome length in base pairs; must be positive |
| `params` | `(sample_id, sequenced_sample_gdna_mass_ng, extracted_gdna_concentration_ng_ul, vol_extracted_elution_ul)`, plus the metric's denominator column | see [Choosing a metric](#choosing-a-metric) |

Extra columns are ignored everywhere, so relations produced for other purposes can be
passed straight in. Each relation's key must be unique.

For `cells_per_g_of_gdna`, only `sequenced_sample_gdna_mass_ng` enters the arithmetic.
The other two parameter columns are required and screened anyway, matching pysyndna,
which filters on them regardless of which metric you ask for — so a sample with a NULL
elution volume is dropped from a calculation that never touches it. Relaxing that would
change which samples appear in your output while every value stayed identical, which is
why it is kept.

The metric's own denominator is screened the same way, and **only for the metric that
divides by it**. A sample missing `calc_mass_sample_aliquot_input_g` disappears from
`cells_per_g_of_sample` and stays put in the other three. That is pysyndna's rule, and
it means changing the metric can change which samples you get back, not only their
units — check the warnings when you switch.

**Where `sequenced_sample_gdna_mass_ng` comes from.** This is the mass of *sample* gDNA
that went into the library prep. It is not the gDNA recovered from the extraction, and it
does not include the mass of the spike-in pool added alongside it. This is the easiest
column on the page to fill with the wrong quantity, because the obvious candidate —
`extracted_gdna_concentration_ng_ul * vol_extracted_elution_ul / 1e9`, the very expression
[Choosing a metric](#choosing-a-metric) uses — is a *different* number, and the arithmetic
accepts either without complaint. Getting it wrong scales every value in the sample.

If you measured the mass, pass the measurement. If you did not, it follows from the
spike-in design, since the pool is added at a known fraction of the sample gDNA mass:

```
sequenced_sample_gdna_mass_ng = mass_syndna_input_ng / syndna_mass_fraction_of_sample
```

With the customary `syndna_mass_fraction_of_sample` of `0.05` that is 20× the pool mass
you already supplied to the fit. pysyndna's Qiita entry points derive the column this way
(`calc_cell_counts.py:414`); miint takes it as an input column instead, matching
`calc_ogu_cell_counts_biom`, the pysyndna function `absquant_cell_counts` mirrors.

### Output

`(sample_id, feature_id, value DOUBLE)`, with both identifier types
[mirrored](#sample-identifier-types) from the `counts` relation.

```sql
SELECT * FROM absquant_cell_counts(
    'ogu_counts', 'models', 'ogu_coverage', 'ogu_lengths', 'sample_params',
    'cells_per_g_of_gdna', 0.1)
ORDER BY sample_id, feature_id LIMIT 5;
```

```
┌───────────┬─────────────────────────────┬────────────────────┐
│ sample_id │         feature_id          │       value        │
│  varchar  │           varchar           │       double       │
├───────────┼─────────────────────────────┼────────────────────┤
│ example1  │ Escherichia coli            │  2260023103719.978 │
│ example1  │ Fusobacterium periodonticum │   732784863204.921 │
│ example1  │ Lactobacillus gasseri       │  5438576851832.264 │
│ example1  │ Leptolyngbya valderiana     │ 1158313590928.3315 │
│ example1  │ Neisseria flavescens        │  958792227127.9537 │
└───────────┴─────────────────────────────┴────────────────────┘
```

### Where coverage comes from

`coverage` is the fraction of each genome that aligned reads actually covered — a guard
against calling a genome present on the strength of a handful of reads landing in one
conserved region. [`genome_coverage_per_sample`](alignment_analysis.md#per-sample-genome-coverage)
produces it directly:

```sql
CREATE TABLE ogu_coverage AS
SELECT sample_id, genome_id AS feature_id, proportion_covered AS coverage
FROM genome_coverage_per_sample(alignments, genome_lengths, contig_to_genome);
```

Per-sample is the only supported shape. If you have one coverage vector for the whole
study, broadcast it explicitly so that the choice is visible in your query rather than
implied:

```sql
CREATE TABLE ogu_coverage AS
SELECT s.sample_id, c.feature_id, c.coverage
FROM (SELECT DISTINCT sample_id FROM ogu_counts) s
CROSS JOIN study_wide_coverage c;
```

**Coverage is a fraction, always.** pysyndna accepts fractions or percents and leaves it
to you to pass a matching `min_coverage`; get that pairing wrong and the filter silently
keeps everything. miint rejects any coverage outside `[0, 1]`, so `92.597` is an error
rather than a disabled filter.

## Which cells you get back

Three filters run, in this order, and **the order matters for the warnings, not the
values**:

| # | Filter | Warning |
|---|---|---|
| 1 | Sample has a NULL, negative or non-finite required parameter — including the requested metric's denominator | `dropped N sample(s) with a NULL, negative or non-finite required parameter: …` |
| 2 | Cell's `coverage < min_coverage` (strict, so a cell exactly on the threshold is kept) | `dropped N feature(s) whose coverage fell below X in at least one sample: …` |
| 3 | Sample has no usable model, or its `rvalue² < min_rsquared` | `no usable model in '…' for N sample(s): …` / `dropped N sample(s) whose model r^2 fell below X: …` |

A sample whose every cell fails the coverage filter produces nothing at all, and is
reported separately (`dropped N sample(s) with no feature at or above X coverage: …`) so
it does not vanish without explanation. All of these go to
[`miint_warnings()`](utilities.md).

One more way to get nothing back: if a sample's every cell computes to exactly zero, the
sparse output omits all of them and the sample disappears. In practice this means a zero
ratio on one of the sample-level metrics — `extracted_gdna_concentration_ng_ul` or
`vol_extracted_elution_ul` is **zero**. Zero is a legitimate measurement, a blank really
can extract no DNA, and pysyndna's screen tests `< 0` so zero passes it, which is why
this is a warning rather than an error:
`N sample(s) produced only zero-valued cells for metric '…' and appear in no output row: …`.

Only the *whole-sample* case is reported. A single zero cell among others is simply
omitted, because a missing row and a dense `0.0` say the same thing (see below) and the
sample is still in the output to be read. It is the all-zero sample that the sparse form
cannot otherwise express.

A **zero denominator** is different and is an error, because dividing by it is not a
measurement problem but an impossible one. pysyndna divides anyway and emits `inf` for
every feature in the sample without saying so. miint refuses, naming the column the
metric you asked for actually divides by:

```
absquant_cell_counts: calc_mass_sample_aliquot_input_g must not be zero; offending sample_id 's1'
```

The ordering is inherited from pysyndna and is deliberate: a sample dropped on its
parameters never contributes its poorly covered features to filter 2's list, because no
surviving sample was going to use them.

A sample with no model is a **warning, not an error** — `absquant_fit_models` routinely
returns fewer models than it was given samples, so the two functions have to compose
without that being fatal.

### Zero-valued cells are omitted

Like every long-form relation in miint, the output is sparse: a cell whose value is zero
has no row. pysyndna pivots and then fills with `0`, so its equivalent output carries an
explicit `0.0` for a feature that survived coverage in one sample but not another. The
two are the same feature table — `test/sql/absquant_cell_counts_biom.test` proves it by
writing both forms to BIOM and showing they produce identical files.

If you need the dense form, ask for it:

```sql
SELECT s.sample_id, f.feature_id, COALESCE(c.value, 0.0) AS value
FROM (SELECT DISTINCT sample_id FROM cells) s
CROSS JOIN (SELECT DISTINCT feature_id FROM cells) f
LEFT JOIN cells c ON c.sample_id = s.sample_id AND c.feature_id = f.feature_id;
```

## Composing the two functions

`absquant_fit_models`' output *is* the `models` argument, with no projection in between:

```sql
CREATE TABLE models AS
SELECT * FROM absquant_fit_models('syndna_counts', 'syndna_pool', 'syndna_params', 1.0,
                                  min_syndna_counts := 50);

CREATE TABLE cells AS
SELECT * FROM absquant_cell_counts('ogu_counts', 'models', 'ogu_coverage', 'ogu_lengths',
                                   'sample_params', 'cells_per_g_of_gdna', 0.1);
```

The two relations must use **the same `sample_id` values**. That sounds obvious and is
the easiest thing to get wrong: the synDNA counts and the OGU counts usually come from
different pipelines, and if their sample naming differs by so much as a prefix, every
sample lands in the "no usable model" warning and the result is empty rather than wrong.
Check the warnings before you check the numbers.

Note also that the two functions are gated independently. `absquant_fit_models` returns
every model it could fit, however poor; `min_rsquared` here is what actually excludes a
weak curve from being used.

## Applying a model by hand

If you want the DNA mass rather than a cell count — the intermediate step of the chain
above — the fitted curve inverts directly:

```sql
-- ng of DNA attributable to each feature, per sample.
SELECT c.sample_id,
       c.feature_id,
       pow(10, m.slope * log10(c.value) + m.intercept) AS ng
FROM feature_counts c
JOIN models m USING (sample_id)
WHERE c.value > 0;
```

Two cautions, which apply to `absquant_cell_counts` equally. The model is fitted on
`log10` of both axes, so it must be applied on that scale and only within the range of
read counts the constructs actually spanned — extrapolating a power law beyond its
calibration is how a spike-in method turns into noise. And it is worth filtering on fit
quality, since a sample whose curve is poor will still produce numbers:

```sql
SELECT * FROM models WHERE rvalue * rvalue >= 0.8;
```

## `absquant_orf_copies`

The third workflow, and the only one that is not about DNA. Given metatranscriptomic read
counts against OGU+ORF references — woltka's output, typically — it reports **copies of
each ORF's ssRNA per gram of sample**.

```sql
absquant_orf_copies(counts, coords, params)
```

There is no standard curve here, no spike-in, and nothing to configure: no `metric`, no
coverage filter, no r² gate, no named parameters. The accounting is direct. What fraction
of the sample's biological reads landed on this ORF, times how much total ssRNA came out
of the extraction, times how many copies a gram of that ORF's ssRNA represents, divided by
how much sample went in:

```
ogu_orf_len        = |ogu_orf_end - ogu_orf_start| + 1
copies_per_g_ssrna = 6.02214076e23 / (ogu_orf_len * 340)
value              = read_count / total_biological_reads_r1r2
                     * total_rna_concentration_ng_ul * vol_extracted_elution_ul / 1e9
                     * copies_per_g_ssrna
                     / calc_mass_sample_aliquot_input_g
```

`340` is the average molar mass of one base of **single-stranded** RNA, in g/mole — the
single-stranded counterpart of the `650` per base *pair* above. Which constant applies is
decided by the molecule, not by the function.

**The quantity is copies of the ORF's own ssRNA, not of the transcript containing it.**
A transcript may carry several ORFs and so be heavier; this counts the ORF. If you want
per-transcript abundance, this is not it, and the difference is not a scaling factor you
can apply afterwards.

`ogu_orf_start > ogu_orf_end` is legal and marks a **reverse-strand** ORF — woltka writes
them that way, and the length is unsigned, so both orientations give the same answer.

### The three RNA relations

All three are passed **by name**, as quoted string literals, like every other relation
argument on this page.

`counts` — long-form OGU+ORF read counts. The same shape `absquant_cell_counts` takes.

| column | type | notes |
|---|---|---|
| `sample_id` | VARCHAR / BIGINT / UUID | mirrored onto the output |
| `feature_id` | VARCHAR / BIGINT / UUID | the OGU+ORF id; mirrored onto the output |
| `value` | DOUBLE | reads assigned to this ORF in this sample |

`coords` — where each ORF sits on its genome. One row per ORF.

| column | type | notes |
|---|---|---|
| `feature_id` | joins to `counts.feature_id` | |
| `ogu_orf_start` | DOUBLE | 1-based, **inclusive** |
| `ogu_orf_end` | DOUBLE | 1-based, **inclusive**; may be less than the start |

`params` — one row per sample. All four parameter columns are required and all four are
used.

| column | type | notes |
|---|---|---|
| `sample_id` | joins to `counts.sample_id` | |
| `calc_mass_sample_aliquot_input_g` | DOUBLE | grams of sample that went into the extraction |
| `total_rna_concentration_ng_ul` | DOUBLE | ng/µL of total ssRNA in the eluate |
| `vol_extracted_elution_ul` | DOUBLE | µL of eluate |
| `total_biological_reads_r1r2` | DOUBLE | total biological reads for the sample, R1 + R2 |

`total_biological_reads_r1r2` is the whole sequencing effort for the sample, not the reads
that mapped to ORFs. Passing the mapped total instead inflates every value by the ratio of
the two, and nothing about the output will look wrong.

Output is `(sample_id, feature_id, value DOUBLE)`, with both id types mirrored from
`counts` — see [Sample identifier types](#sample-identifier-types).

### ORF coordinates from `read_gff`

miint reads no bespoke coordinate format: `coords` is an ordinary relation, so an
annotation in GFF3 becomes one with a projection. **It needs a conversion, and leaving it
out is silent.**

`read_gff` and `read_ncbi_annotation` emit `stop_position` as 1-based **half-open** —
normalized `+1` from GFF3's closed `end`, which is the project-wide convention (issue
#196). `absquant_orf_copies` takes 1-based **closed** coordinates, the convention woltka's
`coords.txt` uses. So:

```sql
CREATE TABLE coords AS
SELECT attributes['ID'] AS feature_id,
       position           AS ogu_orf_start,
       stop_position - 1  AS ogu_orf_end   -- half-open -> closed
FROM read_gff('annotations.gff')
WHERE type = 'CDS';
```

Without the `- 1` every ORF is one base too long, so every copy count comes out low by
`len / (len + 1)` — about 0.3% on a 300 bp ORF. It type-checks, it runs, it raises
nothing, and the error is small enough to survive a sanity check. If your coordinates come
from woltka's `coords.txt` they are already closed and need no adjustment.

### Which ORF cells you get back

A sample is **filtered, with a warning**, when any of the four parameter columns is NULL,
negative or non-finite. This matches pysyndna, which screens the same four.

```sql
SELECT * FROM miint_warnings()
WHERE message LIKE 'absquant_orf_copies:%';
```

A sample is **reported, with a warning**, when it passes that filter and still produces no
row, because every one of its cells came out zero. In practice that means a zero
`total_rna_concentration_ng_ul` or `vol_extracted_elution_ul` — a blank extraction really
can yield no ssRNA, and zero passes the NULL/negative screen legitimately. The output is
sparse, so without the warning that sample would be indistinguishable from one that was
never in the input.

The following are **errors**:

- a zero `total_biological_reads_r1r2` or `calc_mass_sample_aliquot_input_g` — both are
  divided by, and zero is exactly the case the NULL/negative screen cannot catch;
- a read count that is negative or non-finite (zero is fine, and simply means no copies —
  there is no `log10` in this chain);
- an `ogu_orf_start` or `ogu_orf_end` that is not a finite **whole number**. Negative
  coordinates are accepted: the length is well defined either way and woltka's format says
  nothing about the origin;
- duplicate keys in any of the three relations;
- a counted ORF with no `coords` row, or a counted sample with no `params` row. The
  reverse is fine — a `coords` relation is a whole annotation and a `params` relation a
  whole study, so both may describe far more than the counts do, and their unused rows are
  not screened at all;
- a computed value that overflows to `inf`.

Zero-valued cells are omitted from the output, as everywhere else in miint — see
[Zero-valued cells are omitted](#zero-valued-cells-are-omitted).

## Differences from pysyndna

miint's models reproduce
[pysyndna](https://github.com/biocore/pysyndna)'s `fit_linear_regression_models` on its
own published fixtures to the full precision pysyndna emits — it serializes coefficients
truncated to 12 decimal places, and the parity tests hold to `1e-11` relative. The cell
counts reproduce `calc_ogu_cell_counts_biom` on the same fixtures to `1e-9` relative, and
the ORF copies reproduce `calc_copies_of_ogu_orf_ssrna_per_g_sample` to `1e-9` as well —
bit for bit, in fact, on the committed fixtures, because that chain has no `pow` or
`log10` in it, only multiplies and divides, and miint applies them in pysyndna's own
association order. The asserted bound stays `1e-9` regardless: whether the compiler
contracts `a * b / c` into a fused multiply-add is a property of the build, not of the
arithmetic, so bit equality is not something to rely on across platforms.

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

For `absquant_cell_counts` specifically:

- **Coverage must be a fraction in `[0, 1]`.** pysyndna takes fractions or percents and
  documents that matching `min_coverage` to them is the caller's problem. Passing percents
  against a fractional threshold is then a filter that keeps everything, silently, and
  looks like a successful run.
- **A counted cell with no coverage row is an error.** pysyndna validates the sample id
  set and the feature id set *separately*, so a sample and a feature can each be present
  while their cell is not. Its join then yields NaN; `NaN >= min_coverage` is false so the
  cell is dropped, and `NaN < min_coverage` is false too so it never reaches the
  low-coverage log. The cell disappears with no diagnostic of any kind.
- **Every denominator must be positive.** A zero `ogu_len_in_bp`,
  `sequenced_sample_gdna_mass_ng`, or sample-level denominator column is refused.
  pysyndna's parameter screen catches NULL and negative but structurally cannot catch
  zero, which then divides through to `inf` cell counts for every affected feature, with
  no warning.
- **The truncated `6.022e23` Avogadro is not available.** pysyndna carries it behind an
  `is_test` flag to reproduce a rounding in the original published notebook. Offering a
  knob that makes results *less* accurate invites leaving it on.
- **Output is sparse**; see [Zero-valued cells are omitted](#zero-valued-cells-are-omitted).
- **A sample whose every feature fails the coverage filter is reported.** pysyndna emits
  nothing for it — such a sample never reaches its per-sample loop — so it simply vanishes
  from the output. The values agree; miint says so.
- **A sample whose every cell comes out zero is reported.** Reached in practice through a
  sample-level metric with a zero extraction concentration or elution volume. pysyndna's
  output is dense, so it writes the zeros out and the question never arises; a sparse
  output has to say in words what it cannot say in rows.
- **A cell count that overflows to `inf` is an error.** Nothing bounds a model's
  coefficients, so a large enough intercept overflows `10^x`, and an extraction
  concentration times an elution volume can overflow before the divide. pysyndna emits
  the `inf`. An infinite cell count is not a measurement, and the same reasoning applies
  as for the zero denominators above.

For `absquant_orf_copies` specifically:

- **A zero `total_biological_reads_r1r2` or `calc_mass_sample_aliquot_input_g` is an
  error.** Both are divided by, and zero passes pysyndna's NULL/negative screen — it then
  divides through to `inf` copies for every ORF in the sample, silently. A NULL or
  negative value is *not* an error and drops the sample with a warning, which is what
  pysyndna does too: `REQUIRED_PARAM_KEYS` is the union of its sample-info and RNA-prep
  key lists and goes to `filter_data_by_sample_info` whole.
- **A coordinate must be a finite whole number of bases.** pysyndna's `cast_cols` accepts
  `100.5` and produces a fractional ORF length, which divides into Avogadro's number and
  yields a copy count nothing downstream can distinguish from a real one. Negative
  coordinates are accepted, since the length is well defined for any pair of whole
  numbers.
- **A counted ORF with no coordinates is an error.** pysyndna does not validate this
  direction: it reaches the coordinates through a bare `.at[]` inside a `biom.transform`
  lambda, so the same input surfaces as a pandas `KeyError` from inside a callback. The
  sample-id direction it does check, and asymmetrically — extra parameter rows are fine,
  which miint matches.
- **A zero read count is accepted**, unlike in `absquant_cell_counts`. The difference is
  the `log10`: there is none in this chain, so a zero count simply means no copies of that
  ORF. Carrying the stricter rule across would be a rule invented for no reason.
- **A sample whose every cell comes out zero is reported**, for the same reason as above —
  a sparse output has to say in words what it cannot say in rows. Reached through a zero
  `total_rna_concentration_ng_ul` or `vol_extracted_elution_ul`.
- **A copy count that overflows to `inf` is an error.** The extracted ssRNA mass is a
  product formed before the ng→g divide, and the copies-per-gram multiply comes after
  that, so the ceiling sits well below where the two input columns alone would reach.
- **Duplicate keys in any relation are an error.** pysyndna cannot express one: its input
  is a `biom.Table`, unique by construction. A long-form relation can, and a repeat means
  the caller's join fanned out.

## Citations

Zaramela, L.S.; Tjuanta, M.; Moyne, O.; Neal, M.; and Zengler, K. (2022) "synDNA—a
Synthetic DNA Spike-in Method for Absolute Quantification of Shotgun Metagenomic
Sequencing", mSystems, 7(6), e00447-22. doi: 10.1128/msystems.00447-22

pysyndna — the reference implementation, by Amanda Birmingham (biocore). See
`THIRD_PARTY_LICENSES.md`.
