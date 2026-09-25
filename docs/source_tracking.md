# Source tracking (SourceTracker)

Estimate **what fraction of each sink community came from each source environment**, from a
feature table and a sample-metadata relation that says which samples are sources of which
environment and which are sinks.

SourceTracker (Knights et al. 2011) treats every sink as a mixture of the known source
environments plus one **Unknown** source, and estimates the mixing proportions by collapsed
Gibbs sampling: each sequence in the sink is withdrawn and re-assigned to the environment most
likely to have produced its feature, given every other assignment, over and over, and the
proportions are averaged over the retained draws. The implementation is SourceTracker3
([st3](https://github.com/the-miint/st3)), a Rust reimplementation of the SourceTracker2 model;
every default below is SourceTracker2's.

> The examples on this page are executed as tests. Every statement and every number below is
> mirrored in `test/sql/sourcetracker_docs.test`; if you change one, change the other. The
> proportions are Monte-Carlo means, so the test holds them to within 0.05 of what is printed
> rather than to the last digit.

## Table of Contents

- [Inputs](#inputs) - the feature table and the sample metadata.
- [`sourcetracker`](#sourcetracker) - parameters, and why the default call refuses shallow data.
- [The result relation](#the-result-relation) - one row per sink and source.
- [Per-feature assignments](#per-feature-assignments) - which features carried each source's share.
- [Leave-one-out](#leave-one-out) - validate the sources against each other.
- [Rarefaction](#rarefaction) - what the two depths apply to.
- [Things that surprise people](#things-that-surprise-people) - read this before filing a bug.
- [Parallelism](#parallelism)
- [Reproducibility](#reproducibility)
- [Differences from SourceTracker2](#differences-from-sourcetracker2) - a complete list.
- [Citations](#citations)

## Inputs

Two relations, both passed **by name** (see [passing relations by name](table_of_contents.md#passing-relations-by-name)).

**The feature table** is the long-form `(sample_id, feature_id, value)` relation used by
[`read_biom`](reading.md#biom) and everything in [diversity](diversity.md). Values are sequence
counts: SourceTracker withdraws and re-assigns individual sequences, so a cell must be a whole
number, finite and non-negative. A table of relative abundances is refused rather than
truncated to zeros. NULL and zero cells are ignored; a duplicate `(sample_id, feature_id)`
cell is an error, because summing it would change the answer. `sample_id` and `feature_id`
may be `VARCHAR`, `BIGINT` or `UUID`; the types are mirrored onto `sink_id` and onto the keys
of the assignments map, so results join back to typed metadata without a cast.

**The sample metadata** must expose three columns, matched by name without regard to case;
any other column is ignored:

| column | meaning |
|---|---|
| `sample_id` | the sample, matching the feature table's `sample_id` |
| `source_sink` | `source` or `sink` (case-insensitive) |
| `env` | the source's environment; required for every source, ignored for sinks |

Requirements, all enforced with errors that name the offending sample:

- The feature table and the metadata must describe **exactly the same samples**. A sample on
  one side only is an error, in both directions. SourceTracker2 silently intersects the two;
  that is the difference most likely to turn a typo into a wrong answer, so here it is loud.
- At least one source, and (unless `loo := true`) at least one sink.
- No sample with only zero cells: there would be nothing to attribute.
- Every sample listed once in the metadata.

The examples use two soil sources, two gut sources, one water source and two sinks, six
features, about 200 to 270 sequences per sample:

```sql
CREATE TABLE counts AS SELECT * FROM (VALUES
  ('soil1','f1',120),('soil1','f2',100),('soil1','f3',5),('soil1','f4',3),('soil1','f5',2),
  ('soil2','f1',110),('soil2','f2',130),('soil2','f3',4),('soil2','f4',2),('soil2','f6',4),
  ('gut1','f1',3),('gut1','f2',2),('gut1','f3',130),('gut1','f4',110),('gut1','f5',5),
  ('gut2','f2',4),('gut2','f3',120),('gut2','f4',140),('gut2','f5',2),('gut2','f6',4),
  ('water1','f1',2),('water1','f3',3),('water1','f4',5),('water1','f5',140),('water1','f6',120),
  ('sink_a','f1',90),('sink_a','f2',80),('sink_a','f3',20),('sink_a','f4',10),
  ('sink_b','f1',5),('sink_b','f2',5),('sink_b','f3',60),('sink_b','f4',50),('sink_b','f5',40),('sink_b','f6',40)
) AS t(sample_id, feature_id, value);

CREATE TABLE samples AS SELECT * FROM (VALUES
  ('soil1','source','soil'),('soil2','source','soil'),
  ('gut1','source','gut'),('gut2','source','gut'),
  ('water1','source','water'),
  ('sink_a','sink',NULL),('sink_b','sink',NULL)
) AS t(sample_id, source_sink, env);
```

By construction `sink_a` is mostly soil with some gut, and `sink_b` is roughly half gut and
a third water.

## `sourcetracker`

```sql
SELECT * FROM sourcetracker('counts', 'samples',
    loo := false, assignments := false,
    alpha1 := 0.001, alpha2 := 0.1, beta := 10.0,
    restarts := 10, draws_per_restart := 1, burnin := 100, delay := 1,
    source_rarefaction_depth := 1000, sink_rarefaction_depth := 1000,
    with_replacement := false, collapse := 'mean',
    seed := -1, threads := 0);
```

Every named parameter is optional; the values shown are the defaults, which are
SourceTracker2's.

| parameter | default | meaning |
|---|---|---|
| `loo` | `false` | [Leave-one-out](#leave-one-out): predict each source sample from the other sources instead of predicting the sinks. |
| `assignments` | `false` | Fill the [`assignments`](#per-feature-assignments) map. Not available with `loo`. |
| `alpha1` | `0.001` | Prior count of each feature in each source environment. Higher values trust the source data less and smooth its feature distributions; SourceTracker2 calls 0.01 a more conservative choice. |
| `alpha2` | `0.1` | Prior count of each feature in the Unknown environment, as a fraction of the sink's depth. Higher values make Unknown smoother and less prone to absorbing whatever the sources do not explain. |
| `beta` | `10` | Count added to each environment, Unknown included, when the mixing proportions are drawn. |
| `restarts` | `10` | Independent Markov chains per sink. `restarts * draws_per_restart` is the number of retained draws the mean is taken over. |
| `draws_per_restart` | `1` | Draws retained from each chain. |
| `burnin` | `100` | Passes over the sink's sequences before the first draw is retained. |
| `delay` | `1` | Passes between retained draws (thinning). |
| `source_rarefaction_depth` | `1000` | Subsample each collapsed source environment to this many sequences; `0` disables. See [Rarefaction](#rarefaction). |
| `sink_rarefaction_depth` | `1000` | Subsample each sink to this many sequences; `0` disables. Ignored under `loo`. |
| `with_replacement` | `false` | Rarefy with replacement. |
| `collapse` | `'mean'` | How the source samples of one environment are combined before rarefaction: `'mean'` averages their counts (SourceTracker2's behaviour), `'sum'` adds them. |
| `seed` | `-1` | A value `>= 0` fixes the random stream; `-1` draws a fresh seed on every execution. See [Reproducibility](#reproducibility). |
| `threads` | `0` | Sinks sampled concurrently; `0` follows DuckDB's `threads` setting. See [Parallelism](#parallelism). |

Bad parameter values are refused at bind: a `collapse` other than `mean` or `sum`, a `seed`
below `-1`, a prior that is negative or not finite, a chain parameter below 1, a negative
depth or thread count, or `assignments := true` together with `loo := true`.

### The default call refuses shallow data

The example samples are about 200 sequences deep, and the default depths are 1000:

```sql
SELECT * FROM sourcetracker('counts', 'samples');
-- Invalid Input Error: sourcetracker: You requested rarefaction of sink samples at 1000,
-- but 2 sink samples have fewer sequences than that. The shallowest of these is 'sink_a'
-- with 200 sequences.
```

That is SourceTracker2's message, and the same decision: a sink that cannot be rarefied to
the requested depth stops the run rather than being dropped or analysed at its own depth. Set
the depths to something the data supports, or to `0` to analyse every sample as is:

```sql
SELECT sink_id, source, round(proportion, 2) AS proportion, round(proportion_std, 2) AS proportion_std
FROM sourcetracker('counts', 'samples', seed := 42,
                   source_rarefaction_depth := 0, sink_rarefaction_depth := 0)
ORDER BY sink_id, source;
```

| sink_id | source | proportion | proportion_std |
|---|---|---|---|
| sink_a | Unknown | 0.06 | 0.02 |
| sink_a | gut | 0.11 | 0.01 |
| sink_a | soil | 0.83 | 0.02 |
| sink_a | water | 0.00 | 0.00 |
| sink_b | Unknown | 0.14 | 0.05 |
| sink_b | gut | 0.49 | 0.02 |
| sink_b | soil | 0.03 | 0.01 |
| sink_b | water | 0.35 | 0.03 |

## The result relation

One row per **sink** and **source environment**, including `Unknown`; a sink's proportions
sum to 1.

| column | type | meaning |
|---|---|---|
| `sink_id` | mirrors the input `sample_id` type | the sink; under `loo`, the held-out source sample |
| `source` | `VARCHAR` | a source environment, or `Unknown` |
| `proportion` | `DOUBLE` | the mean, over the retained draws, of the fraction of the sink's sequences assigned to this source |
| `proportion_std` | `DOUBLE` | the standard deviation of that fraction over the same draws: how much the chains disagree, not a confidence interval |
| `assignments` | `MAP(feature_id type, DOUBLE)` | `NULL` unless `assignments := true`; see below |

`Unknown` is always present. It collects the sequences that none of the source environments
explains well, and its share is the first thing to look at: a large `Unknown` means the sinks
were drawn from somewhere the sources do not cover.

## Per-feature assignments

`assignments := true` records, for every (sink, source) row, the mean number of the sink's
sequences of each feature that were assigned to that source. This is what SourceTracker2
writes as one feature table per sink; here it is a `MAP` in the same row, with only the
features that received any mass:

```sql
SELECT sink_id, source, assignments
FROM sourcetracker('counts', 'samples', assignments := true, seed := 42,
                   source_rarefaction_depth := 0, sink_rarefaction_depth := 0)
WHERE sink_id = 'sink_a' ORDER BY source;
```

| sink_id | source | assignments |
|---|---|---|
| sink_a | Unknown | `{f1=3.9, f2=3.6, f3=3.1, f4=1.6}` |
| sink_a | gut | `{f1=0.2, f2=0.2, f3=14.3, f4=7.9}` |
| sink_a | soil | `{f1=85.9, f2=76.2, f3=2.6, f4=0.4}` |
| sink_a | water | `{f4=0.1}` |

Every one of the sink's sequences is attributed somewhere, so a sink's maps together sum to its
depth (200 here, or the sink rarefaction depth when one is set), and one row's map sums to
`proportion * depth`. Unnest it to work with the cells as rows:

```sql
SELECT sink_id, source,
       unnest(map_keys(assignments)) AS feature_id,
       round(unnest(map_values(assignments)), 1) AS mean_count
FROM sourcetracker('counts', 'samples', assignments := true, seed := 42,
                   source_rarefaction_depth := 0, sink_rarefaction_depth := 0)
WHERE sink_id = 'sink_a' AND source = 'soil' ORDER BY feature_id;
```

| sink_id | source | feature_id | mean_count |
|---|---|---|---|
| sink_a | soil | f1 | 85.9 |
| sink_a | soil | f2 | 76.2 |
| sink_a | soil | f3 | 2.6 |
| sink_a | soil | f4 | 0.4 |

The map is `NULL` when assignments were not requested and an **empty map** when they were
requested but the source received none of the sink's sequences, so the two cases stay
distinguishable. The tally is collected from the same draws as the proportions, so requesting
it changes nothing else in the result.

## Leave-one-out

`loo := true` holds each **source** sample out in turn and predicts it from the remaining
sources, which is the standard check that the environments are distinguishable at all: a
source sample that does not come back as its own environment says the sources overlap, or the
sample is mislabelled. The rows are source samples, named in `sink_id`; the sinks in the
metadata are ignored, as is `sink_rarefaction_depth`.

```sql
SELECT sink_id AS source_sample, source, round(proportion, 2) AS proportion
FROM sourcetracker('counts', 'samples', loo := true, seed := 42, source_rarefaction_depth := 0)
ORDER BY source_sample, source;
```

| source_sample | source | proportion |
|---|---|---|
| gut1 | Unknown | 0.04 |
| gut1 | gut | 0.93 |
| gut1 | soil | 0.01 |
| gut1 | water | 0.01 |
| gut2 | Unknown | 0.05 |
| gut2 | gut | 0.93 |
| gut2 | soil | 0.01 |
| gut2 | water | 0.01 |
| soil1 | Unknown | 0.04 |
| soil1 | gut | 0.02 |
| soil1 | soil | 0.94 |
| soil1 | water | 0.01 |
| soil2 | Unknown | 0.04 |
| soil2 | gut | 0.01 |
| soil2 | soil | 0.94 |
| soil2 | water | 0.01 |
| water1 | Unknown | 0.99 |
| water1 | gut | 0.01 |
| water1 | soil | 0.00 |
| water1 | water | 0.00 |

`water1` is the only water source. Holding it out leaves no water to attribute to, so its
`water` cell is exactly zero, not merely small, and nearly all of it lands in `Unknown`. That
is the correct answer for a sole-member environment, and the reason a leave-one-out run needs
at least two samples per environment to say anything about that environment.

There is no per-feature tally in leave-one-out; `loo := true` with `assignments := true` is
refused at bind.

## Rarefaction

The two depths apply to different things, in SourceTracker2's order:

- **Sinks** are each subsampled to `sink_rarefaction_depth`. A sink shallower than that stops
  the run, naming the shallowest sink.
- **Sources** are first **collapsed** by environment (`collapse := 'mean'` averages the
  environment's samples feature by feature, `'sum'` adds them), and each collapsed environment
  is then subsampled to `source_rarefaction_depth`. The check is therefore on the collapsed
  totals: the two soil samples above (230 and 250 sequences) sum to 480 but average to 239,
  because the mean is taken feature by feature and rounded down to whole sequences, as
  SourceTracker2 does. A source depth of 250 therefore passes with `'sum'` and is refused
  with `'mean'`.

Either depth may be `0` to skip that subsampling. Rarefaction is a source of randomness in
its own right, so a seeded run with rarefaction on is reproducible but two seeds will differ a
little more than with it off:

```sql
SELECT sink_id, source, round(proportion, 2) AS proportion
FROM sourcetracker('counts', 'samples', seed := 42, collapse := 'sum',
                   source_rarefaction_depth := 200, sink_rarefaction_depth := 200)
ORDER BY sink_id, source;
```

| sink_id | source | proportion |
|---|---|---|
| sink_a | Unknown | 0.06 |
| sink_a | gut | 0.12 |
| sink_a | soil | 0.82 |
| sink_a | water | 0.00 |
| sink_b | Unknown | 0.13 |
| sink_b | gut | 0.49 |
| sink_b | soil | 0.02 |
| sink_b | water | 0.35 |

## Things that surprise people

### The default call fails on most amplicon tables

The default depths are 1000 because SourceTracker2's are, and a table with a sink below
1000 sequences produces the error shown [above](#the-default-call-refuses-shallow-data). This
is deliberate: choosing a depth is part of the analysis, and quietly analysing a shallow sink
at whatever depth it has would hide that the sinks are not comparable.

### Counts must be whole numbers

The sampler moves individual sequences between environments, so a cell of `2.5` has no
meaning and is refused by name. A relative-abundance or normalised table has to be scaled back
to counts first, and doing so with `round()` invents sequences: prefer the original counts.

### The sample sets must match

A sample in the feature table without a metadata row, or a metadata row without any non-zero
cells, is an error naming the sample. SourceTracker2 intersects the two silently. If a run
that used to work now fails here, a filter upstream dropped a sample from one side.

### The sampler cannot be interrupted

The Gibbs sampler runs as one call into st3, and Ctrl-C is honoured before it starts and
after it returns, not inside it. On a large table with many sinks and high `restarts`,
a cancelled query keeps its cores busy until the sampler finishes.

## Parallelism

Sinks (or, under `loo`, held-out source samples) are sampled concurrently. `threads := 0`
follows DuckDB's `threads` setting; an explicit value overrides it for this call. The result
does not depend on the thread count: each sink's chains draw from their own seeded stream, so
`threads := 1` and `threads := 8` give identical rows. The DuckDB-Wasm build compiles the
sampler in but cannot call into it yet; see the
[installation notes](installation.md#optional-feature-flags).

## Reproducibility

With `seed` set to a value `>= 0`, the same inputs, parameters, and st3 version give the same
rows, to the bit, on any thread count. The default `seed := -1` draws a fresh seed for every
execution, including every `EXECUTE` of a prepared statement, which is what SourceTracker2
does (it has no seed parameter).

The proportions are Monte-Carlo estimates. On the example above, a different seed moves them
in the second decimal place; on real data with overlapping environments the spread is larger.
`proportion_std` reports the spread across the retained draws of a single run, which is the
first thing to raise `restarts` and `draws_per_restart` against. st3's own gate against the
SourceTracker2 reference is an absolute 0.02 per cell at 1000 retained draws.

## Differences from SourceTracker2

Everything not listed here follows SourceTracker2: the model, the defaults, collapse-then-
rarefy, the `Unknown` source, and the shallow-sample error.

- **Sample-set mismatches are errors**, in both directions. SourceTracker2 intersects the
  feature table and the metadata silently.
- **`proportion_std` is the standard deviation over draws** of each source's share.
  SourceTracker2's `proportions_std` divides each draw's counts by the total over all draws
  rather than by the sink's depth, so it understates the spread by a factor of the number of
  draws; the means are unaffected and agree.
- **`collapse := 'sum'`** is reachable. SourceTracker2's collapse helper knows a sum as well, but
  its `gibbs` entry point and command line always average.
- **`seed`** exists. SourceTracker2 has none.
- **The random stream differs.** st3 uses its own generator, so results are statistically
  equivalent to SourceTracker2's rather than bit-identical to them.
- **Per-feature assignments are a column**, not one file per sink, and list only the features
  that received mass.
- **Counts must be whole numbers** here; SourceTracker2 requires integer tables too, but
  refuses them at load rather than by cell.

## Citations

- Knights D, Kuczynski J, Charlson ES, Zaneveld J, Mozer MC, Collman RG, Bushman FD, Knight R,
  Kelley ST. Bayesian community-wide culture-independent microbial source tracking.
  *Nature Methods* 8, 761–763 (2011). https://doi.org/10.1038/nmeth.1650
- SourceTracker2: https://github.com/biota/sourcetracker2
- SourceTracker3 (st3), the embedded implementation: https://github.com/the-miint/st3
