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
- [Progressive PCoA (from a distance table)](#progressive-pcoa-from-a-distance-table) - scalable reference-anchored PCoA without a dense N×N decomposition
- [Progressive PCoA (from UniFrac)](#progressive-pcoa-from-unifrac) - scalable UniFrac PCoA computing distances on the fly per batch (true-10M path)
- [PERMANOVA (from a distance table)](#permanova-from-a-distance-table) - metric-agnostic PERMANOVA over any condensed distance table
- [UniFrac PCoA](#unifrac-pcoa) - UniFrac distance + Principal Coordinates Analysis
- [UniFrac PERMANOVA](#unifrac-permanova) - UniFrac distance + PERMANOVA pseudo-F + p-value
- [Faith PD](#faith-pd) - Faith's phylogenetic diversity per sample
- [Procrustes](#procrustes-align-two-ordinations) - align two ordinations into a common frame (disparity M² + PROTEST p-value)

Most of these methods are powered by the embedded [`unifrac-binaries`](https://github.com/biocore/unifrac-binaries) and [`scikit-bio-binaries`](https://github.com/scikit-bio/scikit-bio-binaries) libraries (see `docs/internals/embedded-tools.md` for the build details); [`procrustes`](#procrustes-align-two-ordinations) is a self-contained Eigen-backed port of SciPy and works on any ordination table.

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

<a name="subsample-thread-count"></a>**A seeded subsample reproduces per thread count, not across thread counts.** libssu distributes the draw across the OpenMP team — one generator per thread, seeded in turn — so both how many generators exist and which observation consumes which one depend on the team size. The OpenMP width comes from DuckDB's `SET threads = N` (or the per-call `threads` parameter), so the *same* query with the *same* `seed` and `subsample_depth > 0` can produce a different rarefied table, and therefore different distances, on a differently configured server.

Widths that happen to divide the work identically can coincidentally agree — in one measured case widths 1, 2, and 4 agreed while 8 differed — so matching results at two thread counts is not a guarantee of matching at a third. If you need a rarefied result to be reproducible across machines, fix the width explicitly (`threads := N`, or `SET threads = N`) alongside the seed, and record it with the result. This is a property of the draw itself and applies to every function in this section; it is upstream of, and independent of, the fp32 reduction-ordering noise described under [Reproducibility](#reproducibility).

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
- <a name="pcoa-memory"></a>**Memory (and the size guard):** the ordination itself is dense, needing roughly `5 × n_samples²` bytes — about 3 GB at 25,000 samples, 18 GB at 60,000 — independent of how many pairs the input has. The distances are *streamed* into that matrix (the relation is scanned twice: once to enumerate the sample ids, once to fill), so nothing per-row is held. Because the matrix is allocated outside DuckDB's buffer manager, `memory_limit` cannot bound it on its own; instead `pcoa` estimates the requirement up front and **refuses with a clear error** when it exceeds what `memory_limit` leaves unused, naming both levers (raise `memory_limit`, or switch to [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table), which never forms the dense matrix). Raising `memory_limit` above physical RAM re-enables the OOM this guard exists to prevent.
- **The relation must be stable across scans.** Since it is read twice, a relation that returns different rows each time — a view over `random()` or `nextval()`, or a file rewritten mid-query — is rejected with an error naming the unexpected id, rather than silently producing a wrong or incomplete matrix. Materialize such input into a table first.

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

### Progressive PCoA (from a distance table)

`progressive_pcoa_from_distances(distances, ...)` runs a **scalable, reference-anchored** PCoA over a condensed `(sample_a, sample_b, distance)` relation for sample counts where the dense N×N eigendecomposition [`pcoa`](#pcoa-from-a-distance-table) performs is infeasible, but the condensed distances still fit on disk.

It works by ordinating a small set of shared **anchor** samples once (this fixes a common reference frame), then streaming the remaining samples in batches: each batch is ordinated *together with the anchors* and aligned back onto the reference frame by a [partial procrustes](#procrustes-align-two-ordinations) fit on the anchor overlap. No decomposition larger than `(n_anchors + batch_size)²` is ever computed, and the dense N×N matrix is never materialized.

```sql
CREATE TABLE dm AS SELECT sample_a, sample_b, distance
    FROM unifrac_distances('observations', 'tree', seed := 42);

SELECT * FROM progressive_pcoa_from_distances('dm',
    n_dims := 3, n_anchors := 100, batch_size := 1000, seed := 42);
```

**Parameters:**
- `distances` (VARCHAR): name of the condensed distance relation exposing `(sample_a, sample_b, distance)`
- `n_dims` (INTEGER, default 3): number of PCoA axes; must be `≤ n_samples - 1` and `≤ n_anchors - 1`
- `n_anchors` (INTEGER, default 100): number of anchor samples defining the reference frame; must be `≥ n_dims + 1` and `≤ n_samples`. Anchors are chosen at random (seeded)
- `batch_size` (INTEGER, default 1000): non-anchor samples ordinated per batch (`≥ 1`)
- `seed` (INTEGER, default -1): seeds both the anchor draw and the FSVD randomization; `-1` = unseeded (nondeterministic)

**Output schema:** [`pcoa`](#pcoa-from-a-distance-table)'s six columns — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained)` — plus two appended diagnostic columns, `(batch, batch_anchor_m2)`. `iteration` is always `0`. `sample_id` mirrors the input `sample_a` type — see [Sample identifier types](#sample-identifier-types). **Caveat:** `eigenvalue` and `proportion_explained` are the *anchor* reference ordination's (they describe the anchor subspace, not the full sample set); the per-sample `coordinate`s span all samples.

<a name="progressive-batch-diagnostics"></a>**Reading `batch` / `batch_anchor_m2` (per-run quality evidence).** Each batch is placed into the shared frame by a procrustes fit on its anchor overlap; `batch_anchor_m2` is that fit's disparity — how well this batch's own view of the anchors agreed with the reference view — and `batch` is the 0-based batch that placed the sample. Anchor rows report **NULL** for both: they *define* the frame rather than being fitted into it, so a `0.0` there would be a fabricated perfect fit. This matters because at the scale these functions exist for you cannot check the result against a full `pcoa` — the diagnostic is the accuracy signal you *do* have, and it costs nothing (the fit is computed anyway).

```sql
-- audit a run: worst-fitting batches first
SELECT DISTINCT batch, batch_anchor_m2
FROM progressive_pcoa_from_distances('dm', n_anchors := 1000, batch_size := 1000)
WHERE batch IS NOT NULL ORDER BY batch_anchor_m2 DESC LIMIT 10;
```

Values near 0 mean batches slotted cleanly into the frame. A large value means those samples are poorly determined by the anchor set — typically too few anchors, or anchors that don't span the region that batch occupies. Raise `n_anchors` (accuracy improves monotonically with anchor count) or supply a better `anchors` set.

**Read it as a relative signal, not an error bar.** It measures *frame consistency* and is deliberately conservative: it includes each batch's own ordination noise alongside any genuine disagreement. On a 25,145-sample real dataset the median batch reported `0.101` while the run's actual disparity against a full `pcoa` was `0.0075` — roughly 10× pessimistic. Use it to rank batches within a run and find the badly-placed ones, not to estimate total error. It also cannot detect an anchor set that is self-consistent but collectively unrepresentative of the sample space: every batch can agree with a frame that is itself skewed.

**Behavior:**
- **Accuracy:** the result reproduces a full [`pcoa`](#pcoa-from-a-distance-table) up to a similarity transform — exactly (to numerical precision) for Euclidean-embeddable distances, and closely for others. Each batch is aligned to the reference *independently*, so alignment error does not compound across batches. Validate on your own data by aligning against a full `pcoa` with [`procrustes`](#procrustes-align-two-ordinations) and checking the disparity `m2`.
- **Anchor coordinates are batch-invariant:** for a fixed anchor set and seed the anchor coordinates (and eigenvalues/proportions) do not depend on `batch_size`.
- **Completeness (fail loud):** each batch needs its full `(anchors + batch)²` block present in the distance relation; a missing pair within a block is an error naming the offending pair. NULL sample ids and NULL/NaN distances are skipped. `sample_a`/`sample_b` must resolve to the same output type (a BIGINT/VARCHAR mix is rejected at bind).
- **Fewer than two samples**, `n_anchors` outside `[n_dims + 1, n_samples]`, `n_dims > n_samples - 1`, or `batch_size < 1`: an error.

**Example:**

```sql
-- Ordinate a large precomputed distance table without a dense N×N decomposition,
-- then confirm it matches a full pcoa (small procrustes disparity) on a subset.
SELECT sample_id, axis, coordinate
FROM progressive_pcoa_from_distances('dm', n_dims := 3, n_anchors := 100, batch_size := 1000, seed := 42)
ORDER BY axis, sample_id;
```

> **Note:** the current implementation sources each batch's distance block with its own query against the relation (bounded memory — one block at a time), which re-reads the anchor rows per batch. A future revision caches the anchor-touching rows once and reads batch rows via contiguous ranges.

---

### Progressive PCoA (from UniFrac)

`progressive_pcoa_from_unifrac(feature_table, tree, ...)` is the end-to-end scalable path: it runs the same reference-anchored progressive PCoA as [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table), but computes UniFrac **on the fly, one batch at a time**, directly from a feature table and tree — so the full N×N UniFrac matrix is never formed. It is to [`unifrac_pcoa`](#unifrac-pcoa) what `progressive_pcoa_from_distances` is to [`pcoa`](#pcoa-from-a-distance-table): the memory-bounded ordination for sample counts where the dense decomposition is infeasible.

This is correct because **UniFrac is pairwise-local** — the distance between two samples depends only on their own abundance vectors and the tree, never on which other samples share the table — so a UniFrac block computed over `(anchors + batch)` is identical to the corresponding slice of the full matrix. (Empirically, its output matches `progressive_pcoa_from_distances` over `unifrac_distances(...)` to within a procrustes disparity of ~1e-9.)

```sql
CREATE TABLE observations AS SELECT * FROM read_biom('data/biom/test.biom');
CREATE TABLE tree AS SELECT * FROM read_newick('data/unifrac/gg_otu_tree.nwk');

SELECT * FROM progressive_pcoa_from_unifrac('observations', 'tree',
    variant := 'weighted_normalized', n_dims := 3, n_anchors := 100, batch_size := 1000, seed := 42);
```

**Parameters:**
- `feature_table` (VARCHAR): name of the feature relation exposing `(sample_id, feature_id, value)` (see [Feature table](#feature-table))
- `tree` (VARCHAR): name of the tree relation (see [Tree](#tree))
- `n_dims` (3), `n_anchors` (100), `batch_size` (1000), `seed` (-1), `threads` (0): as in [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table) — `n_anchors` must be in `[n_dims + 1, n_samples]`
- `variant`, `variance_adjust`, `alpha`, `bypass_tips`, `normalize_sample_counts`: the UniFrac controls, identical to [`unifrac_pcoa`](#unifrac-pcoa)

There is deliberately **no `subsample_depth`**: rarefaction and progressive alignment do not compose cleanly (each batch would rarefy independently against a different RNG draw). Rarefy upstream if needed.

**Output schema:** [`unifrac_pcoa`](#unifrac-pcoa)'s six columns — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained)`, `iteration` always `0` — plus the appended `(batch, batch_anchor_m2)` diagnostics, identical in meaning to `progressive_pcoa_from_distances`' (see [Reading `batch` / `batch_anchor_m2`](#progressive-batch-diagnostics)). As with `progressive_pcoa_from_distances`, `eigenvalue`/`proportion_explained` are the *anchor* reference ordination's (a documented caveat). `sample_id` mirrors the input type — see [Sample identifier types](#sample-identifier-types).

**Behavior:**
- **Accuracy / batch-invariance:** same guarantees as [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table) — reproduces a full [`unifrac_pcoa`](#unifrac-pcoa) up to a similarity transform, with alignment error that does not compound across batches.
- <a name="progressive-unifrac-threads"></a>**`threads` controls block CONCURRENCY here, and it is the main performance knob.** Unlike `unifrac_pcoa`, a single block cannot use more than one core no matter how high you set it: libssu's UniFrac parallel degree is `ceil(n_samples / 2048)` over the *block*, and a block is only `n_anchors + batch_size` samples, so at any ordinary setting that is one stripe. So `threads` is spent running that many *blocks* at once instead, each pinned narrow. On 17,483 samples (`n_dims := 5`, `n_anchors := 200`, `batch_size := 500`) this took roughly **250 s at `threads := 1` versus ~32 s at `threads := 14`** — about **8×**, moving from ~1.1 to ~9.8 cores busy. (Timings here and below are indicative: they come from a shared machine, so treat the ratios and the ordering as the result, not the exact seconds.)
- <a name="progressive-unifrac-batch-size"></a>**Tune `m = n_anchors + batch_size`, and keep it small — the opposite of what a serial run wants.** A block's cost grows with `m` faster than linearly, and once blocks run concurrently it is total CPU, not per-block latency, that sets wall time; smaller blocks also mean *more* of them, which is what there is to parallelize. Serially the ordering reverses (big blocks amortize a fixed per-block cost), so a `batch_size` tuned on a one-core run is the wrong one here. The optimum tracks `m`, not `batch_size` alone — at `n_anchors := 200`, `batch_size` 244 → ~31 s, 500 → ~34 s, 1000 → ~56 s, 1849 → ~73 s, 3000 → ~103 s (17,483 samples, `threads := 14`), with total CPU rising ~3.5× across that range. Total CPU is the more trustworthy signal here, since it is far less sensitive to other load on the machine than wall clock. At `n_anchors := 1000` the same table's optimum moves out to `batch_size` 500–1000, because `n_anchors` then dominates `m` and shrinking `batch_size` only adds blocks without shrinking the block.
- **`batch_size` also trades accuracy, and which way depends on how many axes you actually use.** On the real EMP 90 bp deblur table (23,814 samples, `unweighted`, `n_dims := 10`, `n_anchors := 1000`, procrustes vs a full `unifrac_pcoa`): at **d = 3** smaller is better — M² 0.0087 / 0.0106 / 0.0188 for `batch_size` 500 / 1000 / 3000 — while at **d = 10** it inverts, 0.197 / 0.195 / 0.135. Bigger batches align the low-variance axes better; smaller batches align the leading ones better, and per-batch `batch_anchor_m2` rises with `batch_size` (0.081 / 0.117 / 0.160). Since the usual practice is *compute wide, compare narrow* (see the note under `n_dims` below), a smaller `batch_size` is usually the right call — but if you intend to interpret all ten axes, raise it. `batch_size` changes coordinates either way, so pin it alongside the seed when comparing runs.
- **Coordinates are not bit-comparable across `threads` values, but the ordination is preserved.** Two effects, neither a correctness bug. (1) `threads` also sets the OpenMP width of the *anchor* block's reference PCoA, and skbb's centering is an OpenMP `reduction(+:)` whose summation order depends on that width — so the reference frame shifts in its last bits and every coordinate inherits the shift. (2) Concurrent *first* calls race on libssu/skbb's non-atomic CPU/accelerator detection caches, so the **first** progressive run in a session differs slightly from later ones and nothing after it does. `threads := 1` is bit-exact run-to-run.

  **Judge agreement on the ordination, not on coordinates.** The drift is amplified by a randomized SVD and a procrustes fit, so it is larger than one rounding step and it grows with the table: 3.1e-06 max abs on a 40-sample fixture, but 2.4e-05 on 23,814 samples — about 5e-4 relative to the largest coordinate there. What held exactly is the thing that matters: procrustes M² against a full `unifrac_pcoa` was **identical to five decimals (0.19484) at `threads := 1` and `threads := 14`**. So do not diff coordinates across thread counts or use them as a cache key; compare ordinations with [`procrustes`](#procrustes). If you do need identical coordinates, fix `threads := N` alongside the seed and re-run in one session, or use `threads := 1`. See [Reproducibility](#reproducibility) for the same fp32 caveat on the non-progressive paths.
- **Samples & features:** samples with no nonzero feature are excluded (they cannot be ordinated). The tree must cover every feature in the table (validated once at bind); each batch's features are a subset, so the check makes them all safe.
- **Fewer than two samples**, `n_anchors` outside `[n_dims + 1, n_samples]`, `n_dims > n_samples - 1`, `batch_size < 1`, an unknown `variant`, or a tree missing a feature: an error.

> **Note:** like `progressive_pcoa_from_distances`, this computes each batch's block with its own feature-table slice query, so memory stays bounded by one `(anchors + batch)` block. The anchor samples' feature rows are read once and cached for the whole run, so each batch's query fetches only its own samples.

<a name="feature-table-sort-order"></a>**Performance: store the feature table sorted by `sample_id`.** Samples are batched in sorted id order, so a batch is a contiguous id range. If the table's physical layout is also sorted, DuckDB/Parquet prune by each row group's min/max `sample_id` and a batch's slice query reads only that batch's own rows; if it is not, every slice query scans the whole `sample_id` column, so the run's slicing cost grows as `n_batches × table_rows` instead of one pass.

Measured on the 17,483-sample EMP deblur table remapped to short (WoL-style) feature ids, and an 8× replication of it, `batch_size := 1000`, one batch's slice query:

| feature table | sorted by `sample_id` | unsorted |
|---|---|---|
| 17,483 samples / 13.1 M rows | 6 ms | 7 ms |
| 139,864 samples / 104.4 M rows | 9 ms | 33 ms |

Sorted stays flat as the table grows; unsorted scales with it. The gap widens with the width of `feature_id`: repeating the 13.1 M-row measurement with the original 150 bp ASV sequences as ids gave 8 ms sorted vs **140 ms** unsorted, because an unpruned scan pays for the wide string column too. For ASV-keyed tables at scale this is the difference between minutes and hours.

`read_biom` emits each sample's rows contiguously, but in the BIOM file's own sample order — **grouped is not sorted**, and grouping alone prunes nothing, since each row group then spans nearly the whole id range. Add the `ORDER BY` when you materialize:

```sql
CREATE TABLE observations AS
    SELECT * FROM read_biom('table.biom') ORDER BY sample_id;
-- or, writing Parquet
COPY (SELECT * FROM read_biom('table.biom') ORDER BY sample_id)
    TO 'observations.parquet' (FORMAT PARQUET);
```

**With a BIGINT `sample_id`, sort by the text form.** Samples are enumerated and batched by `sample_id::VARCHAR`, so batches are *lexical* ranges (`1, 10, 100, 1000, …, 2, 20`) and a numerically sorted table will not line up with them — use `ORDER BY sample_id::VARCHAR`. VARCHAR and UUID ids sort identically either way (see [Sample identifier types](#sample-identifier-types)).

Sort order affects only how much a run reads — never the coordinates it produces.

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
- **Subsampling** (`subsample_depth > 0`): per iteration, `subsample_table_inmem_seeded(seed + i)` is called and Faith PD is computed on the bridged result. Same seed **at the same thread count** → byte-identical reconstruction (no fp32 tolerance needed); across thread counts the draw itself changes, see [seeded subsampling and thread count](#subsample-thread-count).
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

### Procrustes (align two ordinations)

`procrustes(reference, other, ...)` superimposes one ordination onto another with the optimal similarity transform (translation, uniform scaling, and an orthogonal rotation/reflection), so two ordinations of the same (or overlapping) samples can be compared or plotted in a common frame. Unlike the UniFrac-powered methods above, `procrustes` is self-contained — an Eigen-backed port of [`scipy.spatial.procrustes`](https://docs.scipy.org/doc/scipy/reference/generated/scipy.spatial.procrustes.html) / `scipy.linalg.orthogonal_procrustes` (BSD-3) — so it works on **any** long-form ordination table, not just PCoA output.

Both inputs are relation *names* exposing `(sample_id, axis, coordinate)` — the same long-form shape [`pcoa`](#pcoa-from-a-distance-table) and [`unifrac_pcoa`](#unifrac-pcoa) emit (`axis` 0-indexed).

```sql
-- Two PCoA ordinations of the same samples (e.g. two metrics), aligned:
CREATE TABLE ord_a AS SELECT sample_id, axis, coordinate FROM pcoa('dm_bray',    n_dims := 3, seed := 42);
CREATE TABLE ord_b AS SELECT sample_id, axis, coordinate FROM pcoa('dm_unifrac', n_dims := 3, seed := 42);

SELECT * FROM procrustes('ord_a', 'ord_b', permutations := 999, seed := 42);
```

**Parameters:**
- `reference` (VARCHAR): name of the reference ordination relation `(sample_id, axis, coordinate)`; its frame is the target
- `other` (VARCHAR): name of the ordination to transform onto the reference
- `pairing` (VARCHAR, optional): name of a `(reference_id, other_id)` relation. Absent → **full** mode (both ordinations must describe the same samples). Present → **partial** mode: fit the transform on just the paired anchor rows, then apply it to *every* row of both ordinations (the q2-diversity `partial_procrustes` technique, qiime2/q2-diversity#338)
- `n_dims` (INTEGER, default: all available axes): number of leading axes to use; must be `≤` the axes present in each input
- `permutations` (INTEGER, default 999): Monte Carlo permutations for the PROTEST p-value (full mode only); `0` disables the test (p-value is NULL)
- `seed` (INTEGER, default -1): permutation seed; `-1` = unseeded

**Output schema:** mirrors the [`unifrac_pcoa`](#unifrac-pcoa) shape, with a leading `matrix` discriminator in place of `iteration` and the two fit-level scalars in the trailing slots:
- `matrix` (VARCHAR): `'reference'` (the standardized reference) or `'other'` (the transformed other)
- `sample_id` (mirrors input type — see [Sample identifier types](#sample-identifier-types)): sample identifier. The output type is the shared BIGINT/UUID type when `reference` and `other` agree, otherwise VARCHAR
- `axis` (INTEGER): 0-indexed axis
- `coordinate` (DOUBLE): the sample's coordinate on this axis, in the shared standardized frame
- `m2` (DOUBLE): the Procrustes disparity M² — **replicated on every row** (it is a property of the whole fit). `M² = 1 − (Σσ)² ∈ [0, 1]`; `0` = a perfect superimposition
- `pvalue` (DOUBLE): PROTEST Monte Carlo p-value — **replicated on every row**; `NULL` in partial mode and when `permutations := 0`

**Behavior:**
- **Full mode:** `reference` and `other` must describe the **same** sample set and carry the **same** number of axes (mirrors `scipy.spatial.procrustes`, which requires identical shapes). The transform is fit on all shared samples; `m2` is the disparity; `pvalue` is the PROTEST test.
- **Partial mode (`pairing`):** the pairing must be **1:1** — a repeated `reference_id` or `other_id` is rejected (it would drop a row or feed one physical sample into the fit as several anchors). The fit uses the matched anchor rows (at least `n_dims + 1` usable pairs are required); it is then applied to every row of both inputs. `m2` is the disparity over the anchors; `pvalue` is `NULL` (matching q2's `partial_procrustes`, which defines no Monte Carlo test).
- **Sample ids:** VARCHAR/BIGINT/UUID are all accepted (see [Sample identifier types](#sample-identifier-types)); full-mode matching and pairing lookups are by the id's string form, so a BIGINT `reference` and a VARCHAR `other` with the same ids align and emit under VARCHAR.
- **Fail loud:** a missing `(sample_id, axis, coordinate)` column, a ragged ordination (a sample missing an axis), a duplicate `(sample_id, axis)`, an out-of-range `axis`, fewer than `n_dims + 1` points/anchors, or (full mode) differing sample sets or axis counts each raise a named error.
- **P-value reproducibility:** the PROTEST test is reproducible under a fixed `seed` within one build, but it is **not** bit-for-bit comparable to q2's `procrustes_analysis` — q2 uses an *unseeded* RNG, and the C++ PRNG differs from NumPy's, so agreement is statistical (within Monte Carlo error), not exact. Disparity and coordinates, by contrast, match SciPy to machine precision.

**Examples:**

```sql
-- Per-fit disparity + p-value (one row via DISTINCT, since both are replicated):
SELECT DISTINCT m2, pvalue
FROM procrustes('ord_a', 'ord_b', permutations := 999, seed := 42);

-- Partial (anchored) alignment: project a new batch `ord_b` into the reference
-- frame using a set of shared anchor samples, carrying the batch's extra samples
-- through the same transform.
CREATE TABLE anchors AS SELECT ref_id AS reference_id, batch_id AS other_id FROM anchor_map;
SELECT sample_id, axis, coordinate
FROM procrustes('ord_a', 'ord_b', pairing := 'anchors')
WHERE matrix = 'other'
ORDER BY sample_id, axis;
```
