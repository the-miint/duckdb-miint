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
- [Beta-distance macros](#beta-distance-macros) - within/between-group distributions and k-nearest-neighbors over a distance table
- [PCoA (from a distance table)](#pcoa-from-a-distance-table) - metric-agnostic PCoA over any condensed distance table
- [PERMANOVA (from a distance table)](#permanova-from-a-distance-table) - metric-agnostic PERMANOVA over any condensed distance table
- [UniFrac PCoA](#unifrac-pcoa) - UniFrac distance + Principal Coordinates Analysis
- [UniFrac PERMANOVA](#unifrac-permanova) - UniFrac distance + PERMANOVA pseudo-F + p-value
- [Faith PD](#faith-pd) - Faith's phylogenetic diversity per sample

These methods are powered by the embedded [`unifrac-binaries`](https://github.com/biocore/unifrac-binaries) and [`scikit-bio-binaries`](https://github.com/scikit-bio/scikit-bio-binaries) libraries (see `docs/internals/embedded-tools.md` for the build details).

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
