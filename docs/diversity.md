# Diversity calculations and statistics

Methods to estimate alpha and beta diversity, and supporting statistics.

## Table of Contents

- [Feature table](#feature-table) - Feature table expectations.
- [Tree](#tree) - Tree expectations.
- [Metadata](#metadata) - Metadata expectations.
- [Sample identifier types](#sample-identifier-types) - How VARCHAR/BIGINT/UUID sample ids are handled.
- [UniFrac algorithm variants](#unifrac-algorithm-variants) - Detail on different UniFrac algorithms and how to specify them.
- [Rarefaction](#rarefaction) - Rarefaction detail with UniFrac and Faith PD.
- [UniFrac distances](#unifrac-distances) - Condensed (pairwise) UniFrac distances in long form
- [Community distances (non-phylogenetic)](#community-distances-non-phylogenetic) - taxon-based β-diversity (Bray-Curtis, Jaccard, Morisita-Horn, χ², Gower, …) from a feature table
- [Beta-distance macros](#beta-distance-macros) - within/between-group distributions and k-nearest-neighbors over a distance table
- [PCoA (from a distance table)](#pcoa-from-a-distance-table) - metric-agnostic PCoA over any condensed distance table
- [PERMANOVA (from a distance table)](#permanova-from-a-distance-table) - metric-agnostic PERMANOVA over any condensed distance table
- [Two-sample Kolmogorov-Smirnov test](#two-sample-kolmogorov-smirnov-test) - `ks_2samp` over two numeric lists, with SciPy-exact p-values
- [Sample clustering (k-means and UPGMA)](#sample-clustering-k-means-and-upgma) - group samples from ordination coordinates or a distance table
- [UniFrac PCoA](#unifrac-pcoa) - UniFrac distance + Principal Coordinates Analysis
- [UniFrac PERMANOVA](#unifrac-permanova) - UniFrac distance + PERMANOVA pseudo-F + p-value
- [Faith PD](#faith-pd) - Faith's phylogenetic diversity per sample
- [Citations](#citations) - primary sources for the metrics, clustering, and ordination methods

These methods are powered by the embedded [`unifrac-binaries`](https://github.com/biocore/unifrac-binaries) and [`scikit-bio-binaries`](https://github.com/scikit-bio/scikit-bio-binaries) libraries (see `docs/internals/embedded-tools.md` for the build details).

> Relating **two** feature tables of the same samples — microbes against metabolites, say — is a different problem from the ordination and distance methods here, because correlating compositional tables directly produces spurious associations. See [multi-omics integration](multiomics.md) for MMvec.

### Feature table

A long-form `(sample_id VARCHAR, feature_id VARCHAR, value DOUBLE)` relation — the same schema produced by [`read_biom`](reading.md#biom). Duplicate `(sample_id, feature_id)` rows are summed; rows with NULL or zero values are silently dropped (sparse-storage invariant).

```sql
CREATE TABLE observations AS SELECT * FROM read_biom('data/biom/test.biom');
```

### Tree

A relation matching the [`read_newick`](reading.md#newick) schema: `(node_index BIGINT, parent_index BIGINT, name VARCHAR, branch_length DOUBLE, edge_id BIGINT)`. The tree's tip names must include every `feature_id` present in the feature table — otherwise binding fails with an error naming the missing tip.

```sql
CREATE TABLE tree AS SELECT * FROM read_newick('data/unifrac/gg_otu_tree.nwk');
```

### Metadata 

Wide-form: one row per sample, with a `sample_id` column (case-insensitive) plus one column per variable. Values are cast to VARCHAR. The `variables` parameter selects which columns to test; default is every non-`sample_id` column.

```sql
CREATE TABLE metadata AS SELECT * FROM read_csv('data/unifrac/small_metadata.csv');
-- sample_id,body_site,treatment
-- Sample1,gut,A
-- Sample2,gut,B
-- ...
```

Samples appearing in the metadata but not in the feature-table are silently ignored. Samples missing from the metadata fail loudly with an error naming both the missing sample and the variable.

### Sample identifier types

The feature table's `sample_id` column may be `VARCHAR`, `BIGINT`, or `UUID` (any other type is accepted and read as text). `unifrac_distances`, `unifrac_pcoa`, and `unifrac_faith_pd` **mirror the input `sample_id` type onto their output identifier columns** — `BIGINT → BIGINT`, `UUID → UUID`, otherwise `VARCHAR` — so results join back to typed metadata without a cast (parity with `align_minimap2`). `unifrac_permanova` emits no per-sample identifier column, so this does not apply to it.

### UniFrac algorithm variants 

`unifrac_distances`, `unifrac_pcoa`, and `unifrac_permanova` accept the following `variant` values (passed unquoted to libssu, with `_fp32` appended internally):

| Variant | Description |
|---|---|
| `weighted_normalized` (default) | Weighted UniFrac normalized to [0, 1] |
| `weighted_unnormalized` | Weighted UniFrac without normalization |
| `unweighted` | Classic unweighted UniFrac |
| `unweighted_unnormalized` | Unweighted without normalization |
| `generalized` | GUniFrac (controlled by `alpha`) |

### Rarefaction

`unifrac_distances`, `unifrac_pcoa`, `unifrac_permanova`, and `unifrac_faith_pd` all accept `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, and `seed` parameters.

- `subsample_depth := 0` (default): no subsampling; iterations would be identical so `n_subsamples > 1` is rejected at bind.
- `subsample_depth > 0`: libssu rarefies every sample to exactly `subsample_depth` total counts. **Samples whose total count is below the depth are dropped** from the result; depending on the function's invariants this can also cause a bind-time error (e.g., `unifrac_pcoa` requires `n_dims ≤ n_samples - 1` after subsampling).
- `subsample_with_replacement := true` uses multinomial sampling; the default `false` uses permutation without replacement.
- `seed := -1` (default) uses system entropy; any `seed >= 0` is deterministic. Per-iteration seed is `seed + iteration_index`; bind rejects `seed + n_subsamples - 1 > INT_MAX`.

---

### UniFrac distances

Compute the UniFrac distance matrix and return it in **condensed long form** — one row per unordered sample pair, `(iteration, sample_a, sample_b, distance)`. This is the primitive underneath [PCoA](#unifrac-pcoa) and [PERMANOVA](#unifrac-permanova); exposing it directly lets you build within/between-group distributions, k-nearest-neighbors, and correlations between distance matrices with ordinary SQL (see [Beta-distance macros](#beta-distance-macros)).

```sql
SELECT * FROM unifrac_distances('observations', 'tree',
    variant := 'weighted_normalized',
    variance_adjust := false,
    alpha := 1.0,
    bypass_tips := false,
    normalize_sample_counts := true,
    subsample_depth := 0,
    subsample_with_replacement := false,
    n_subsamples := 1,
    seed := -1);
```

**Parameters:**
- `observations` (VARCHAR): name of the feature-table relation
- `tree` (VARCHAR): name of the tree relation
- `variant` (VARCHAR, default `'weighted_normalized'`): one of the [variants](#unifrac-algorithm-variants) above
- `variance_adjust` (BOOLEAN, default false): apply Chang & Liu's variance adjustment
- `alpha` (DOUBLE, default 1.0): GUniFrac alpha, only relevant when `variant := 'generalized'`
- `bypass_tips` (BOOLEAN, default false): skip tip-level contributions (≈50% faster, mildly different result)
- `normalize_sample_counts` (BOOLEAN, default true): normalize each sample to relative abundances before computing distances
- `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, `seed`: see [Rarefaction](#rarefaction)

**Output schema:**
- `iteration` (INTEGER): 0-indexed iteration; constant when `n_subsamples = 1`
- `sample_a`, `sample_b` (mirror input type — see [Sample identifier types](#sample-identifier-types)): the two samples of the pair
- `distance` (DOUBLE): the UniFrac distance between them (fp32-computed, widened to DOUBLE)

**Behavior:**
- **Condensed / upper triangle:** each unordered pair is emitted exactly once, no self-pairs. Pairs are ordered by the canonical *lexical* sample-id form, so `sample_a < sample_b` holds for VARCHAR; for BIGINT/UUID the direction follows the text form (e.g. `10` may sort before `2`) — uniqueness holds regardless.
- **Streaming:** the full N×N matrix is held in memory once (O(N²) — intrinsic to the distance computation), but the pairs are streamed out, so the O(N²) *result rows* are never all materialized at once.
- **Fewer than two samples:** an empty result (the condensed form of a <2-sample matrix has no pairs). This differs from `unifrac_pcoa`, which errors — an ordination is undefined for <2 samples, whereas "no pairs" is a valid distance answer.
- **Rarefaction:** with `subsample_depth > 0`, samples whose total count falls below the depth are dropped from that iteration; each iteration emits its own triangle, tagged by `iteration`.
- **Reproducibility:** same `seed` → same distances at fp32 precision (compare with a `1e-5` tolerance; see the [PCoA reproducibility note](#reproducibility)).

> **Output size:** the result has `N·(N-1)/2` rows per iteration — ~50M rows at 10k samples, growing quadratically. Aggregate or filter it (the macros below, a `GROUP BY`, a `WHERE`) rather than `SELECT *`-ing it into a client.

**Examples:**

```sql
-- Condensed distances for the whole table
SELECT * FROM unifrac_distances('observations', 'tree', seed := 42)
ORDER BY sample_a, sample_b;

-- Within/between-group distribution via a plain SQL join to metadata
SELECT CASE WHEN a.body_site = b.body_site THEN 'within' ELSE 'between' END AS comparison,
       count(*) AS n, avg(distance) AS mean,
       quantile_cont(distance, [0.25, 0.5, 0.75]) AS quartiles
FROM unifrac_distances('observations', 'tree') d
JOIN metadata a ON d.sample_a = a.sample_id
JOIN metadata b ON d.sample_b = b.sample_id
GROUP BY comparison;
```

---

### Community distances (non-phylogenetic)

`community_distances(feature_table, metric, ...)` computes classic **taxon-based** (aphylogenetic) sample×sample distances from a feature table — the counterpart to [`unifrac_distances`](#unifrac-distances) for when you have no tree, or deliberately want a metric that ignores one. Its output is the same condensed `(sample_a, sample_b, distance)` triple, so it feeds [`pcoa`](#pcoa-from-a-distance-table), [`permanova`](#permanova-from-a-distance-table) and the [beta-distance macros](#beta-distance-macros) unchanged.

It takes a relation *name*, so materialize the feature table first:

```sql
CREATE TABLE ft AS SELECT * FROM read_biom('table.biom');
CREATE TABLE dm AS SELECT * FROM community_distances('ft', 'bray_curtis');

SELECT * FROM pcoa('dm', n_dims := 3, seed := 42);
```

**Parameters:**
- `feature_table` (VARCHAR): name of a long-form `(sample_id, feature_id, value)` relation — see [Feature table](#feature-table)
- `metric` (VARCHAR): one of the metrics below (case-insensitive)
- `threads` (INTEGER, default 0): worker threads for the pairwise loop; `0` follows DuckDB's thread count, `1` runs serial

**Output schema:** `(sample_a, sample_b, distance DOUBLE)` — the upper triangle only, one row per unordered pair, no self-pairs. `sample_a`/`sample_b` mirror the input `sample_id` type — see [Sample identifier types](#sample-identifier-types).

**Metrics.** With `x`, `y` two sample rows, sums over features, `X = Σx`, `Y = Σy`:

| metric | formula | range | notes | primary source |
|---|---|---|---|---|
| `bray_curtis` | `Σ\|xₖ−yₖ\| / Σ(xₖ+yₖ)` | [0,1] | empty pair → 0 | Bray & Curtis 1957 |
| `euclidean` | `sqrt(Σ(xₖ−yₖ)²)` | [0,∞) | | |
| `jaccard` | binary presence/absence `(b+c)/(a+b+c)` | [0,1] | **presence/absence**, not abundance; empty pair → 0 | Jaccard 1912 |
| `soergel` | `Σ\|xₖ−yₖ\| / Σ max(xₖ,yₖ)` | [0,1] | empty pair → 0 | |
| `morisita_horn` | `1 − 2Σ(xₖyₖ) / ((Σxₖ²/X² + Σyₖ²/Y²)·X·Y)` | [0,1] | Horn's Cλ on relative abundances; both-empty → 0, one-empty → 1 | Morisita 1959; Horn 1966; Magurran 2004 p.246 |
| `pearson` | `1 − r` over features | [0,2] | constant row → 0 vs another constant row, 1 vs a non-constant one | |
| `chisq` | `sqrt(Σₖ (GT/colₖ)(xₖ/X − yₖ/Y)²)` | [0,∞) | correspondence-analysis χ²; `GT` = grand total; zero-sum row → 0 vs another empty row, 1 otherwise | Faith, Minchin & Belbin 1987 |
| `gower` | `Σₖ \|xₖ−yₖ\| / rangeₖ` | [0,∞) | un-normalized; `rangeₖ` over all samples | Gower 1971; Faith, Minchin & Belbin 1987 |

The zero-variance and zero-row-sum conventions (the `pearson` and `chisq` notes
above) follow `cogent3.maths.distance_transform`, the reference implementation of
the metrics used by Kuczynski et al. 2010 — see [Citations](#citations).

**Behavior:**
- **Raw values, no pre-normalization.** Abundances are used exactly as given; each metric applies whatever internal normalization its own definition requires. Feeding counts versus relative abundance is your modeling choice. (At equal per-sample depth the two differ only by a constant factor for most metrics, which does not change ordination geometry.)
- **Matrix-wide metrics.** `chisq` and `gower` depend on *global* column statistics (column sums and column ranges across all samples), so a pair's distance is a function of the whole matrix, not just that pair. Subsetting samples changes them.
- **Sparse contract:** a sample whose values are all zero or NULL is dropped before computation, matching `unifrac_distances`. Distances are therefore emitted only among samples that actually carry signal.
- **Threading is exact:** results are **bit-identical for any thread count** (each pair writes to a fixed condensed slot), so `threads` is purely a performance knob.
- Fewer than two samples yields no rows.

**Example — comparing what different metrics see:**

```sql
-- Presence/absence vs abundance overlap can disagree sharply; that disagreement
-- is often the biologically interesting part.
CREATE TABLE dj AS SELECT * FROM community_distances('ft', 'jaccard');
CREATE TABLE dmh AS SELECT * FROM community_distances('ft', 'morisita_horn');

SELECT j.sample_a, j.sample_b, j.distance AS jaccard, m.distance AS morisita_horn
FROM dj j JOIN dmh m USING (sample_a, sample_b)
ORDER BY abs(j.distance - m.distance) DESC
LIMIT 10;
```

> Available in builds with the UniFrac feature enabled (the default) — it shares that feature's table readers.

---

### Beta-distance macros

Three SQL macros operate over any **condensed distance table** with columns `(sample_a, sample_b, distance)` — the [`unifrac_distances`](#unifrac-distances) output, or any distance relation you build yourself. They take relation *names* (a table or view), following the `query_table` convention, so materialize the distances first:

```sql
-- Macros take a relation name, not a table function call
CREATE TABLE dm AS SELECT * FROM unifrac_distances('observations', 'tree', seed := 42);
```

> **One iteration at a time:** these macros are `iteration`-blind. If the distance table carries multiple iterations (`n_subsamples > 1`), filter to a single iteration first (e.g. build `dm` with `n_subsamples := 1`, or `... WHERE iteration = 0`) — otherwise distributions and neighbor rankings pool across bootstrap replicates. Pass the distances **as a materialized table**, not a view over `unifrac_distances`, since it is referenced more than once.

#### `beta_group_distances(distances, groups)`

Labels each pair within- or between-group. `groups` is a `(sample_id, grouping)` relation — reduce your metadata to the one grouping column of interest:

```sql
CREATE VIEW groups AS SELECT sample_id, body_site AS grouping FROM metadata;

-- Aggregated within/between distribution (the group-cohesion question)
SELECT comparison, count(*) AS n,
       quantile_cont(distance, [0.25, 0.5, 0.75]) AS quartiles
FROM beta_group_distances(dm, groups)
GROUP BY comparison ORDER BY comparison;
```

Returns `(sample_a, sample_b, distance, group_a, group_b, comparison)` where `comparison` is `'within'` when `group_a = group_b` else `'between'`. Notes: `sample_id` must be unique in `groups` (a duplicate fans the joins out and inflates counts); a pair whose either endpoint is absent from `groups` is dropped; two NULL-group samples are labeled `'between'` (SQL `NULL = NULL` is NULL).

**Per-group-vs-rest:** the condensed table holds each pair once, so a naive `GROUP BY group_a` attributes a between-pair only to `sample_a`'s group. To get, for each group, its within-group distances and its distances to every other group, attribute each *between* pair to **both** its groups but each *within* pair to its group only once (both endpoints share it):

```sql
WITH bg AS (SELECT * FROM beta_group_distances(dm, groups))
SELECT grp, comparison, count(*) AS n, avg(distance) AS mean FROM (
    -- within pairs: count once for their (single) group
    SELECT group_a AS grp, comparison, distance FROM bg WHERE comparison = 'within'
    UNION ALL
    -- between pairs: count once for each of the two groups
    SELECT group_a AS grp, comparison, distance FROM bg WHERE comparison = 'between'
    UNION ALL
    SELECT group_b AS grp, comparison, distance FROM bg WHERE comparison = 'between'
) GROUP BY grp, comparison ORDER BY grp, comparison;
```

#### `beta_knn(distances, k)`

The `k` nearest neighbors of every sample (both orientations of the condensed table are considered). Returns `(sample_id, neighbor, distance, rank)`, `rank` in `1..k`, ties broken by neighbor id:

```sql
SELECT * FROM beta_knn(dm, 5) ORDER BY sample_id, rank;
```

#### `beta_knn_from_sample(distances, k, source)`

The `k` samples nearest one `source` sample. Returns `(neighbor, distance)`, nearest first:

```sql
SELECT * FROM beta_knn_from_sample(dm, 10, 'Sample1');
```

#### Correlating two distance matrices (Mantel)

The Mantel **r-statistic** — the correlation between two distance matrices over the same samples — is a one-liner over two condensed tables (both must be over the same sample set, so their `(sample_a, sample_b)` pairs align):

```sql
SELECT corr(x.distance, y.distance) AS mantel_r
FROM dm_weighted x JOIN dm_unweighted y USING (sample_a, sample_b);
```

The permutation-based Mantel **p-value** is not yet implemented; it is tracked in [issue #160](https://github.com/the-miint/duckdb-miint/issues/160).

---

### PCoA (from a distance table)

`pcoa(distances, ...)` runs Principal Coordinates Analysis over **any** condensed distance table with columns `(sample_a, sample_b, distance)` — the [`unifrac_distances`](#unifrac-distances) output, a [beta-distance](#beta-distance-macros) relation, or a precomputed Bray-Curtis / Jaccard / Euclidean table you built yourself. It is **metric-agnostic**: the ordination is decoupled from UniFrac. [`unifrac_pcoa`](#unifrac-pcoa) remains the end-to-end convenience wrapper (feature table + tree → coordinates); `pcoa` is the same ordination over a distance matrix you already have.

It takes a relation *name* (a table or view), so materialize the distances first:

```sql
CREATE TABLE dm AS SELECT sample_a, sample_b, distance
    FROM unifrac_distances('observations', 'tree', seed := 42);

SELECT * FROM pcoa('dm', n_dims := 3, seed := 42);
```

**Parameters:**
- `distances` (VARCHAR): name of the condensed distance relation exposing `(sample_a, sample_b, distance)`
- `n_dims` (INTEGER, default 3): number of PCoA axes to compute; must be `≤ n_samples - 1`
- `seed` (INTEGER, default -1): FSVD randomization seed; `-1` = unseeded
- `threads` (INTEGER, default 0): OpenMP threads; `0` follows DuckDB's thread count

**Output schema:** identical to [`unifrac_pcoa`](#unifrac-pcoa) — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained)`. `iteration` is always `0` (kept for schema parity; a distance table is a single fixed matrix, so there is no subsampling). `sample_id` mirrors the input `sample_a` type — see [Sample identifier types](#sample-identifier-types).

**Behavior:**
- **Input shape:** distinct ids are collected from *both* `sample_a` and `sample_b`, sorted lexicographically, and assembled into a dense symmetric matrix with a zero diagonal. `sample_a` and `sample_b` must resolve to the same output type (a BIGINT/VARCHAR mix is rejected at bind).
- <a name="distance-table-completeness"></a>**Completeness (fail loud):** the matrix must be complete — every unordered pair present exactly once. A missing pair, a negative or non-finite distance, a nonzero self-distance, or a pair given two conflicting values is an error that names the offending pair. Rows with a NULL `sample_a`/`sample_b` are skipped; a NULL/NaN `distance` is treated as "not provided" (its ids are still recorded, so a sample whose every distance is NULL/NaN surfaces as an incompleteness error rather than silently vanishing). `unifrac_distances` always emits the full triangle, so its output is complete by construction; a hand-built table you must complete yourself.
- **Equivalence:** `pcoa` over `unifrac_distances(obs, tree, variant := v, seed := s)` reproduces `unifrac_pcoa(obs, tree, variant := v, n_dims := d, seed := s)` — same matrix, same FSVD call, same seed. Coordinates match up to axis sign (compare `abs(coordinate)`) within the fp32 tolerance noted under [Reproducibility](#reproducibility).
- **Fewer than two samples**, or **`n_dims > n_samples - 1`**: an error (an ordination is undefined).

**Example:**

```sql
-- Bray-Curtis PCoA with no tree involved: build the distance table however you
-- like (here a placeholder), then ordinate it.
CREATE TABLE bc AS SELECT sample_a, sample_b, distance FROM my_bray_curtis_distances;
SELECT sample_id, axis, coordinate
FROM pcoa('bc', n_dims := 2, seed := 42)
ORDER BY axis, sample_id;
```

---

### PERMANOVA (from a distance table)

`permanova(distances, metadata, ...)` runs PERMANOVA (pseudo-F + a permutation p-value) over any condensed distance table and a wide-form metadata relation. Like [`pcoa`](#pcoa-from-a-distance-table) it is **metric-agnostic** — the omnibus test is decoupled from UniFrac. [`unifrac_permanova`](#unifrac-permanova) remains the end-to-end wrapper.

```sql
CREATE TABLE dm AS SELECT sample_a, sample_b, distance
    FROM unifrac_distances('observations', 'tree', seed := 42);

SELECT * FROM permanova('dm', 'metadata',
    variables := ['body_site'], n_permutations := 999, seed := 42);
```

**Parameters:**
- `distances` (VARCHAR): name of the condensed distance relation exposing `(sample_a, sample_b, distance)`
- `metadata` (VARCHAR): name of the wide-form metadata relation — a `sample_id` column plus one column per variable (see [Metadata](#metadata))
- `n_permutations` (INTEGER, default 999): permutations for the p-value
- `variables` (VARCHAR[], default all non-`sample_id` columns): metadata columns to test
- `seed` (INTEGER, default -1): permutation seed; `-1` = unseeded
- `threads` (INTEGER, default 0): OpenMP threads; `0` follows DuckDB's thread count

**Output schema:** identical to [`unifrac_permanova`](#unifrac-permanova) — `(iteration, variable, n_groups, f_stat, p_value, n_permutations)`. `iteration` is always `0`. There is no `sample_id` output, so no id-type mirroring.

**Behavior:**
- **Input shape / completeness:** the same dense-matrix construction and [completeness requirement](#distance-table-completeness) as `pcoa`.
- **Metadata alignment:** every sample in the distance matrix must have a row for each tested variable, else an error naming the missing sample and variable; metadata samples not present in the distance matrix are ignored. A NULL metadata value is treated as an (empty-string) group value — filter the metadata first if you want NULL samples dropped.
- **Equivalence:** `permanova` over `unifrac_distances(obs, tree, variant := v, seed := s)` reproduces `unifrac_permanova(obs, tree, metadata, variant := v, seed := s)` under the same `seed` and `n_permutations`: `p_value` is byte-identical (same permutations) and `f_stat` matches within `1e-5` (fp32; see the [reproducibility note](#reproducibility)).

---

### Two-sample Kolmogorov-Smirnov test

`ks_2samp(a, b [, method])` is the two-sided two-sample Kolmogorov–Smirnov test: a scalar over two numeric lists, returning the statistic and an **exact** p-value.

It is the last step of the micov coverage-curve comparison — with [`cumulative_coverage_curve`](alignment_analysis.md#cumulative-coverage) producing curves in-database, comparing two metadata groups runs entirely in SQL.

```sql
SELECT ks_2samp([1.0, 2.0, 3.0], [1.0, 2.0, 3.0]);
-- {'statistic': 0.0, 'pvalue': 1.0}

SELECT (ks_2samp([1.0, 2.0], [10.0, 20.0])).statistic;
-- 1.0
```

**Parameters:**
- `a`, `b` (`LIST(DOUBLE)`): the two samples. They need not be the same length. `LIST(INTEGER)`/`LIST(BIGINT)`/`LIST(DECIMAL)` are accepted by implicit cast, so the coverage macros' output feeds in directly.
- `method` (VARCHAR, optional, default `'auto'`): `'auto'` or `'exact'` (case-insensitive). See [Method and the size ceiling](#method-and-the-size-ceiling).

**Returns** a `STRUCT`:

| Field | Type | Description |
|---|---|---|
| `statistic` | DOUBLE | the two-sided statistic `D = sup\|F_a(x) − F_b(x)\|`, in `[0, 1]` |
| `pvalue` | DOUBLE | `P(D ≥ statistic)` under H₀ that both samples are drawn from the same continuous distribution |

**Behavior:**
- **NULL `a`, `b`, or `method` → NULL.**
- **Either list empty → NULL**, not `0.0`. `0.0` would read as "the two samples are identical", the opposite of what an empty group means.
- **NULL elements *inside* a list are dropped** and the test runs on the remainder, matching [`compress_intervals`](alignment_analysis.md) and the QC aggregates. The sample size shrinks with them, so the p-value reflects the surviving observations. A list whose every element is NULL reduces to the empty case, and therefore to NULL.
- **NaN is an error**, not a dropped value: it has no position in the sort order, so the ECDF would be undefined, and silently discarding it would change the sample size the p-value is computed against. `±Infinity` is fine — it sorts.
- **Ties are handled correctly in the statistic**, which matters here more than usual: coverage curves are mostly ties. The ECDFs are compared only after *every* observation equal to the current pooled value is consumed, so a sample compared against itself gives `D = 0` even when it contains repeats.
- **But the p-value is conservative under ties.** The exact null distribution assumes all `C(n₁+n₂, n₁)` interleavings are equally likely, which requires a continuous distribution; ties break that, so the reported p-value is somewhat *larger* than the tie-aware permutation p-value — the test is under-powered, never anti-conservative. This matters because the headline use case is tie-heavy by construction. SciPy has the identical limitation, so parity is unaffected.
- **The p-value can lose accuracy, and eventually underflow, when `n` is large *and* `D` is large.** The value itself stays good far further down than you might expect: it is accurate to ~1e-15 relative across the normal range and remains **bit-identical to SciPy at `6.6e-318`** for equal-sized disjoint samples, well inside the subnormal range. What breaks it is not the magnitude of the *answer* but underflow of the sweep's *intermediate* probabilities, which needs a large `n` together with a large `D`. Measured: 9000 vs 5000 at `D = 0.411` gives `1.86e-321` where SciPy gives `8.90e-321` — a factor of ~4.8, because both accumulate in probability space but in different orders. Past that, the answer underflows to exactly `0.0` (two disjoint samples of 10000 have a true p-value near `1e-6019`), which means "below double precision", not "impossible under H₀", and propagates through a Bonferroni or FDR correction as `0.0`. Treat a p-value that small as "vanishingly small" rather than as a number.
- **Symmetric to the bit:** `ks_2samp(a, b)` and `ks_2samp(b, a)` return identical values, not merely values within 1 ULP.
- **Order-independent:** the lists are sorted internally, so `list(v)` without an `ORDER BY` is fine.

#### Method and the size ceiling

The p-value is **exact** — the null distribution over all `C(n₁+n₂, n₁)` interleavings, not an asymptotic approximation. `'auto'` and `'exact'` therefore behave identically today; both names are accepted for signature parity with SciPy.

**`ks_2samp` raises when `max(n₁, n₂) > 10000`** rather than returning an approximation. Two distinct things meet at that number and only one of them is SciPy's:

- For **`'auto'`**, 10000 is SciPy's own `MAX_AUTO_N`: above it SciPy switches to `'asymp'`, which miint does not implement — see below.
- For **`'exact'`**, SciPy has *no* size cap; it computes exactly at any size, bailing only when `lcm(n₁,n₂) ≥ 2³¹`. There the ceiling is **miint's own cost bound**, not SciPy behaviour, and a conservative one — the lcm at 10000×9999 is `1.0e8`, three orders below SciPy's guard. Raising it needs only a cost decision, not new mathematics.

The `'asymp'` reason is worth stating in full, because the obvious substitute is wrong:

SciPy's `ks_2samp` switches `method='auto'` to `'asymp'` above that same size, and its asymptotic branch is *not* the textbook Kolmogorov series `2·Σ(−1)^(k−1)·exp(−2k²λ²)`. It is `kstwo.sf(d, round(n₁n₂/(n₁+n₂)))` — the exact *one-sample* KS distribution at an effective sample size — which SciPy evaluates by selecting per `(n, d)` region among the Durbin-matrix, Pomeranz, and Pelz-Good expansions. Measured against SciPy 1.18.0, the plain Kolmogorov series is 0.7% off in the most favourable reachable case and 47% off at small effective `n`, which `'auto'` does reach when the samples are lopsided. Since micov's published statistics came from SciPy, shipping an approximation would silently break reproducibility of already-published results — so miint declines to answer instead. [Issue #218](https://github.com/the-miint/duckdb-miint/issues/218) records what implementing it would take.

This ceiling is not a practical constraint for the curve-comparison use case: `n` there is the number of **samples in a group** (one curve rank per sample), so 10–100 against a ceiling of 10000.

Cost grows with both the sample sizes and the observed statistic — the exact computation sweeps a band of an `n₁ × n₂` lattice whose width is proportional to `D`, so a *small* `D` at a large `n` is cheap and only a large `D` at a large `n` is not. Measured: 100 vs 100 is 6 ms, 1000 vs 1000 is 4 ms, 10000 vs 10000 with a small `D` is 2 ms, and 10000 vs 10000 with `D = 1` — the widest band, and the worst case at the ceiling — is 513 ms. 2000 separate 100-vs-100 comparisons in one query total 88 ms. A lopsided pair is cheap in memory as well as time: the working state is `O(min(n₁, n₂))`, so 10000 vs 5 needs six doubles, not ten thousand.

#### Examples

Long-form input, one list per group:

```sql
SELECT ks_2samp(list(v) FILTER (WHERE grp = 'A'),
                list(v) FILTER (WHERE grp = 'B')) AS ks
FROM observations;
```

Pairwise comparison of cumulative coverage curves between metadata groups, with a Bonferroni correction done in SQL. No multiple-comparison function is provided — the correction is one expression, and which one to apply is the analyst's call:

```sql
WITH curves AS (
    SELECT group_id, list(proportion_covered ORDER BY rank) AS curve
    FROM cumulative_coverage_curve(positions, roster, 4719737)
    GROUP BY group_id
),
pairs AS (
    SELECT a.group_id AS group_a, b.group_id AS group_b,
           ks_2samp(a.curve, b.curve) AS ks
    FROM curves a JOIN curves b ON a.group_id < b.group_id
)
SELECT group_a, group_b,
       ks.statistic,
       ks.pvalue,
       least(1.0, ks.pvalue * COUNT(*) OVER ()) AS pvalue_bonferroni
FROM pairs
ORDER BY ks.pvalue;
```

#### SciPy parity

`statistic` and `pvalue` match `scipy.stats.ks_2samp` to double precision wherever `ks_2samp` answers, except where the sweep's intermediate probabilities underflow (see the accuracy note above). Parity is pinned by a committed golden of SciPy 1.18.0 outputs (`test/sql/ks_2samp_parity.test`, fixtures under `data/ks2samp/`) in two grids:

- **32 small cases** covering unequal and coprime sizes, lopsided pairs in both orientations, heavy ties, curve-shaped input, `n = 1`, identical samples, and disjoint samples down to `p = 2.2e-59`.
- **5 large cases** at `n₁` of 5000–10000 (with `n₂` down to 3000), including one at `p = 1.7e-98`, so the claim is evidenced across the size range rather than only at the small end. Their inputs are arithmetic ranges the test rebuilds in SQL, which is why reaching `n = 10000` costs ten oracle rows instead of forty thousand fixture rows.

The generator is committed at `test/scripts/generate_ks2samp_oracle.py` and reproduces both goldens byte for byte, so a future SciPy bump can be audited rather than guessed at. The test itself needs no SciPy installed.

One-sided variants (`less` / `greater`) are not implemented; only the two-sided form is.

---

### Sample clustering (k-means and UPGMA)

Two clustering functions group **samples** — as opposed to [sequence clustering](clustering.md), which groups reads into OTUs. Both are metric-agnostic and sit downstream of an ordination or a distance table, which makes them useful for scoring how well a distance choice recovers known structure.

#### `cluster_kmeans(coords_table, k, ...)`

Lloyd's k-means (Lloyd 1982) with k-means++ seeding (Arthur & Vassilvitskii 2007) over ordination coordinates — the `(sample_id, axis, coordinate)` long form emitted by [`pcoa`](#pcoa-from-a-distance-table) / [`unifrac_pcoa`](#unifrac-pcoa).

```sql
CREATE TABLE dm AS SELECT * FROM community_distances('ft', 'jaccard');
CREATE TABLE co AS SELECT * FROM pcoa('dm', n_dims := 3, seed := 42);

SELECT * FROM cluster_kmeans('co', k := 3, seed := 42);
```

**Parameters:**
- `coords_table` (VARCHAR): name of a relation exposing `(sample_id, axis, coordinate)`
- `k` (INTEGER, required): number of clusters
- `seed` (BIGINT, default 0): seeds all randomness
- `max_iter` (INTEGER, default 100): Lloyd iterations per restart
- `n_init` (INTEGER, default 10): k-means++ restarts; the lowest-inertia one wins
- `n_dims` (INTEGER, default 0): use only the first `n_dims` axes; `0` uses every axis present

**Output schema:** `(sample_id, cluster INTEGER)`, one row per sample. `sample_id` mirrors the input type. Cluster ids are canonicalized in order of first-seen sample, so equivalent solutions get a stable labelling.

**Behavior:**
- **Cluster count matches scikit-learn.** When the points span at least `k` distinct locations you get exactly `k` non-empty clusters. When there are fewer distinct locations than `k`, you get fewer — the same thing `sklearn.KMeans` does (it warns), since `k` non-empty clusters are then impossible.
- **Determinism:** fixed `seed` reproduces on a given build. The `std::uniform_*` distribution mappings are not standardized across C++ standard libraries, so exact draws may differ between e.g. libstdc++ and libc++.
- Empty clusters are re-seeded to the farthest still-unclaimed point.
- Errors on a duplicate `(sample_id, axis)` pair, a missing cell in the sample×axis grid, `k < 1`, or `k >` the number of samples.

#### `cluster_upgma(distances)`

Average-linkage hierarchical clustering (UPGMA; Sokal & Michener 1958) over any condensed distance table.

```sql
SELECT * FROM cluster_upgma('dm');
```

**Parameters:**
- `distances` (VARCHAR): name of a condensed `(sample_a, sample_b, distance)` relation

**Output schema:** the [`read_newick`](reading.md#newick) tree-table schema — `(node_index BIGINT, name VARCHAR, branch_length DOUBLE, edge_id BIGINT, parent_index BIGINT, is_tip BOOLEAN)` — deliberately, so every existing tree utility works on the result without conversion. Leaf `name` is the sample id as text; `parent_index` is NULL at the root. Tips occupy `node_index` `0 .. n−1`, internal nodes `n .. 2n−2`, root last.

**Behavior:**
- The tree is **ultrametric**: a merge at distance `d` places its node at height `d/2`, and branch lengths are height differences.
- Input must satisfy the same [completeness requirement](#distance-table-completeness) as `pcoa`.
- Dense O(n²) memory and O(n³) time, with no size cap — the same characteristics as `pcoa`/`permanova`.

**Example — is a known grouping recovered as clades?**

```sql
-- Walk each leaf's ancestors, then ask which nodes have a descendant leaf set
-- exactly equal to one of the true groups.
WITH RECURSIVE t AS (SELECT * FROM cluster_upgma('dm')),
anc(leaf, node) AS (
    SELECT name, node_index FROM t WHERE is_tip
    UNION ALL
    SELECT a.leaf, p.node_index FROM anc a JOIN t c ON c.node_index = a.node
                                           JOIN t p ON p.node_index = c.parent_index
)
SELECT node, count(*) AS n_leaves, list_sort(list(leaf)) AS leaves
FROM anc GROUP BY node ORDER BY n_leaves;
```

---

### UniFrac PCoA 

Compute a UniFrac distance matrix and reduce it to PCoA coordinates via [scikit-bio's randomized FSVD](https://arxiv.org/abs/1007.5510).

```sql
SELECT * FROM unifrac_pcoa('observations', 'tree',
    variant := 'weighted_normalized',
    n_dims := 3,
    variance_adjust := false,
    alpha := 1.0,
    bypass_tips := false,
    normalize_sample_counts := true,
    subsample_depth := 0,
    subsample_with_replacement := false,
    n_subsamples := 1,
    seed := -1);
```

**Parameters:**
- `observations` (VARCHAR): name of the feature-table relation
- `tree` (VARCHAR): name of the tree relation
- `variant` (VARCHAR, default `'weighted_normalized'`): one of the [variants](#unifrac-algorithm-variants) above
- `n_dims` (INTEGER, default 3): number of PCoA axes to compute; must be `≤ n_samples - 1`
- `variance_adjust` (BOOLEAN, default false): apply Chang & Liu's variance adjustment
- `alpha` (DOUBLE, default 1.0): GUniFrac alpha, only relevant when `variant := 'generalized'`
- `bypass_tips` (BOOLEAN, default false): skip tip-level contributions (≈50% faster, mildly different result)
- `normalize_sample_counts` (BOOLEAN, default true): normalize each sample to relative abundances before computing distances; set false for absolute-quants pipelines
- `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, `seed`: see [Subsampling](#rarefaction)

**Output schema:**
- `iteration` (INTEGER): 0-indexed iteration; constant when `n_subsamples = 1`
- `sample_id` (mirrors input type — see [Sample identifier types](#sample-identifier-types)): sample identifier (post-subsample ordering)
- `axis` (INTEGER): 0-indexed PCoA axis; `axis = 0` is the largest eigenvalue
- `coordinate` (DOUBLE): the sample's coordinate on this axis
- `eigenvalue` (DOUBLE): eigenvalue for this axis — **replicated across all samples within `(iteration, axis)`**
- `proportion_explained` (DOUBLE): fraction of total variance — **replicated across all samples within `(iteration, axis)`**

#### Replication note

`eigenvalue` and `proportion_explained` are properties of the axis, not the sample. They are emitted once per row for SQL ergonomics (a single long-form relation joins cleanly with other per-sample data), but their values are constant within each `(iteration, axis)` pair. To get the per-axis summary use `DISTINCT`:

```sql
-- Per-axis variance summary (one row per axis per iteration)
SELECT DISTINCT iteration, axis, eigenvalue, proportion_explained
FROM unifrac_pcoa('observations', 'tree', n_dims := 3, seed := 42)
ORDER BY iteration, axis;
```

#### Reproducibility

Same seed → same output **at fp32 precision**. The randomized FSVD plus cblas/LAPACKE reduction ordering produces ~1e-7 bit-level noise across runs even with `OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1`. For seed-reproducibility checks, compare with a tolerance:

```sql
WITH a AS (SELECT * FROM unifrac_pcoa('observations', 'tree', seed := 42)),
     b AS (SELECT * FROM unifrac_pcoa('observations', 'tree', seed := 42))
SELECT bool_and(abs(a.coordinate - b.coordinate) < 1e-5)
FROM a JOIN b USING (iteration, sample_id, axis);
```

Axis **signs** are also not stable across calls (a known property of any randomized SVD); compare coordinates as `abs(coordinate)` or perform Procrustes alignment for cross-iteration aggregation.

#### Examples

```sql
-- Single-iteration default
SELECT * FROM unifrac_pcoa('observations', 'tree', seed := 42)
ORDER BY iteration, sample_id, axis;

-- Multi-iteration subsampled PCoA
SELECT * FROM unifrac_pcoa('observations', 'tree',
    n_dims := 3, subsample_depth := 3, n_subsamples := 100, seed := 42);

-- 2D scatter of the first two axes
SELECT sample_id,
       MAX(CASE WHEN axis = 0 THEN coordinate END) AS pc1,
       MAX(CASE WHEN axis = 1 THEN coordinate END) AS pc2
FROM unifrac_pcoa('observations', 'tree', n_dims := 2, seed := 42)
GROUP BY sample_id;
```

#### UniFrac PERMANOVA

Test whether group assignments explain variation in UniFrac distances, using the standard PERMANOVA pseudo-F statistic.

```sql
SELECT * FROM unifrac_permanova('observations', 'tree', 'metadata',
    variant := 'weighted_normalized',
    n_permutations := 999,
    variables := ['body_site', 'treatment'],
    variance_adjust := false,
    alpha := 1.0,
    bypass_tips := false,
    normalize_sample_counts := true,
    subsample_depth := 0,
    subsample_with_replacement := false,
    n_subsamples := 1,
    seed := -1);
```

**Parameters:**
- `observations`, `tree` (VARCHAR): feature-table and tree relation names
- `metadata` (VARCHAR): wide-form metadata relation name (see [Metadata](#metadata))
- `variant`, `variance_adjust`, `alpha`, `bypass_tips`, `normalize_sample_counts`: same as `unifrac_pcoa`
- `n_permutations` (INTEGER, default 999): number of permutations for the p-value
- `variables` (VARCHAR[], default all non-`sample_id` columns): metadata columns to test
- `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, `seed`: see [Subsampling](#rarefaction)

**Output schema:**
- `iteration` (INTEGER): 0-indexed iteration
- `variable` (VARCHAR): canonical metadata column name (preserves case from the table)
- `n_groups` (INTEGER): number of distinct values for the variable in the canonical sample order
- `f_stat` (DOUBLE): pseudo-F statistic
- `p_value` (DOUBLE): one-sided p-value, fraction of permuted F-stats ≥ observed
- `n_permutations` (INTEGER): copy of the input parameter for downstream record-keeping

#### Behavior

- **Variable resolution:** `variables` is case-insensitive but the emitted `variable` column uses the canonical column name from the metadata table, not the user-supplied spelling.
- **Identical groupings produce identical f_stat:** factorization is value-driven and order-stable (encoded by first-appearance in canonical sample order). Two variables that partition samples identically under different labels — e.g. `treatment = {A, B, A, B, ...}` and `condition = {0, 1, 0, 1, ...}` — get the same `labels[]` array and therefore the same `f_stat` for a given seed.
- **Single-group variables** (`n_groups = 1`) produce a row with non-finite `f_stat` (libskbb returns `-inf` for the degenerate case) and `p_value = 1.0`. Filter those out at the SQL layer if undesired.
- **Missing metadata fails loud:** a sample present in the feature-table but absent for a requested variable throws a bind-time error naming both the sample and the variable.
- **Seed reproducibility:** `skbb_permanova_fp32` reconstructs its RNG from `seed` at every call, so the permutation sequence is deterministic. `p_value`, `n_groups`, and `n_permutations` are byte-identical across invocations under the same seed. **`f_stat` is not** — it's computed on the same fp32 UniFrac distance matrix that PCoA uses, so the cblas/LAPACKE reduction-ordering noise (~1e-7) carries through. Use a `abs(f_stat_a - f_stat_b) < 1e-5` tolerance when comparing across invocations.

#### Examples

```sql
-- Default: test every metadata column
SELECT * FROM unifrac_permanova('observations', 'tree', 'metadata', seed := 42);

-- Only specific variables, more permutations
SELECT * FROM unifrac_permanova('observations', 'tree', 'metadata',
    variables := ['body_site'], n_permutations := 9999, seed := 42);

-- Multi-iteration over a rarefaction: average pseudo-F across subsamples
SELECT variable, avg(f_stat) AS mean_f, median(p_value) AS median_p
FROM unifrac_permanova('observations', 'tree', 'metadata',
    subsample_depth := 3, n_subsamples := 50, seed := 42)
GROUP BY variable;
```

#### Faith PD 

Compute Faith's phylogenetic diversity (sum of branch lengths on the minimal subtree spanning a sample's tips) for each sample.

```sql
SELECT * FROM unifrac_faith_pd('observations', 'tree',
    subsample_depth := 0,
    subsample_with_replacement := false,
    n_subsamples := 1,
    seed := -1);
```

**Parameters:**
- `observations`, `tree` (VARCHAR): feature-table and tree relation names
- `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, `seed`: see [Subsampling](#rarefaction)

**Output schema:**
- `iteration` (INTEGER): 0-indexed iteration; always 0 when `subsample_depth = 0`
- `sample_id` (mirrors input type — see [Sample identifier types](#sample-identifier-types)): sample identifier
- `faith_pd` (DOUBLE): sum of branch lengths on the spanning subtree (non-negative)

#### Behavior

- **No subsampling** (`subsample_depth = 0`): one row per input sample; output is deterministic regardless of `seed`. `n_subsamples > 1` is rejected because iterations would be identical.
- **Subsampling** (`subsample_depth > 0`): per iteration, `subsample_table_inmem_seeded(seed + i)` is called and Faith PD is computed on the bridged result. Same seed → byte-identical reconstruction (no fp32 tolerance needed).
- **Subsampling drops low-count samples:** if a sample's total count falls below `subsample_depth`, it does not appear in that iteration's output.

#### Examples

```sql
-- Per-sample Faith PD on the full table
SELECT * FROM unifrac_faith_pd('observations', 'tree') ORDER BY sample_id;

-- Rarefied Faith PD: mean across 100 iterations at depth 3.
-- Samples whose total count falls below subsample_depth are silently
-- dropped from individual iterations, so `count(*)` over the grouping
-- below tells you how many iterations actually contributed per sample.
-- Add `HAVING count(*) = 100` (or similar) to filter out samples that
-- weren't deep enough to be rarefied to depth 3.
SELECT sample_id, count(*) AS n_iter, avg(faith_pd) AS mean_pd,
       stddev_samp(faith_pd) AS sd_pd
FROM unifrac_faith_pd('observations', 'tree',
    subsample_depth := 3, n_subsamples := 100, seed := 42)
GROUP BY sample_id;
```

---

### Citations

If you use these methods, please cite the primary sources.

**Community distances.** Bray, J.R. and Curtis, J.T. (1957) "An ordination of the
upland forest communities of southern Wisconsin", *Ecological Monographs* 27(4),
325-349. · Jaccard, P. (1912) "The distribution of the flora in the alpine zone",
*New Phytologist* 11(2), 37-50. · Gower, J.C. (1971) "A general coefficient of
similarity and some of its properties", *Biometrics* 27(4), 857-871. · Faith,
D.P., Minchin, P.R. and Belbin, L. (1987) "Compositional dissimilarity as a robust
measure of ecological distance", *Vegetatio* 69, 57-68. · Morisita, M. (1959)
"Measuring of interspecific association and similarity between communities",
*Memoirs of the Faculty of Science, Kyushu University, Series E* 3, 65-80. ·
Horn, H.S. (1966) "Measurement of 'overlap' in comparative ecological studies",
*The American Naturalist* 100(914), 419-424. · Magurran, A.E. (2004) *Measuring
Biological Diversity*, Blackwell (p.246 for the Morisita-Horn form used here).

**Sample clustering.** Lloyd, S.P. (1982) "Least squares quantization in PCM",
*IEEE Transactions on Information Theory* 28(2), 129-137. · Arthur, D. and
Vassilvitskii, S. (2007) "k-means++: the advantages of careful seeding",
*Proceedings of the 18th Annual ACM-SIAM Symposium on Discrete Algorithms*,
1027-1035. · Sokal, R.R. and Michener, C.D. (1958) "A statistical method for
evaluating systematic relationships", *University of Kansas Science Bulletin* 38,
1409-1438. · Rand, W.M. (1971) "Objective criteria for the evaluation of
clustering methods", *Journal of the American Statistical Association* 66(336),
846-850 (the Rand index used to score cluster recovery).

**Ordination.** Torgerson, W.S. (1952) "Multidimensional scaling: I. Theory and
method", *Psychometrika* 17, 401-419. · Gower, J.C. (1966) "Some distance
properties of latent root and vector methods used in multivariate analysis",
*Biometrika* 53(3-4), 325-338. · Anderson, M.J. (2001) "A new method for
non-parametric multivariate analysis of variance", *Austral Ecology* 26(1), 32-46
(PERMANOVA). · Lozupone, C. and Knight, R. (2005) "UniFrac: a new phylogenetic
method for comparing microbial communities", *Applied and Environmental
Microbiology* 71(12), 8228-8235. · Faith, D.P. (1992) "Conservation evaluation and
phylogenetic diversity", *Biological Conservation* 61(1), 1-10.

**Reference implementations.** The `cogent3` / PyCogent, NumPy, SciPy,
scikit-learn and scikit-bio projects — consulted for metric conventions and used
as parity oracles — are credited with their licenses and citations in
[`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md). Note that `ks_2samp` is
checked against `scipy.stats.ks_2samp` but is **not** derived from it: the exact
p-value comes from the lattice-path formulation in Hodges (1958) and is computed
as escaping probability mass rather than as `1 - P(stay inside)`. See the SciPy
entry there, and `src/include/KsTwoSample.hpp`.
