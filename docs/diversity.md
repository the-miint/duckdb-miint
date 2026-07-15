# Diversity calculations and statistics

Methods to estimate alpha and beta diversity, and supporting statistics.

## Table of Contents

- [Feature table](#feature-table) - Feature table expectations.
- [Tree](#tree) - Tree expectations.
- [Metadata](#metadata) - Metadata expectations.
- [UniFrac algorithm variants](#unifrac-algorithm-variants) - Detail on different UniFrac algorithms and how to specify them.
- [Rarefaction](#rarefaction) - Rarefaction detail with UniFrac and Faith PD.
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

### UniFrac algorithm variants 

`unifrac_pcoa` and `unifrac_permanova` accept the following `variant` values (passed unquoted to libssu, with `_fp32` appended internally):

| Variant | Description |
|---|---|
| `weighted_normalized` (default) | Weighted UniFrac normalized to [0, 1] |
| `weighted_unnormalized` | Weighted UniFrac without normalization |
| `unweighted` | Classic unweighted UniFrac |
| `unweighted_unnormalized` | Unweighted without normalization |
| `generalized` | GUniFrac (controlled by `alpha`) |

### Rarefaction

All three functions accept `subsample_depth`, `subsample_with_replacement`, `n_subsamples`, and `seed` parameters.

- `subsample_depth := 0` (default): no subsampling; iterations would be identical so `n_subsamples > 1` is rejected at bind.
- `subsample_depth > 0`: libssu rarefies every sample to exactly `subsample_depth` total counts. **Samples whose total count is below the depth are dropped** from the result; depending on the function's invariants this can also cause a bind-time error (e.g., `unifrac_pcoa` requires `n_dims ≤ n_samples - 1` after subsampling).
- `subsample_with_replacement := true` uses multinomial sampling; the default `false` uses permutation without replacement.
- `seed := -1` (default) uses system entropy; any `seed >= 0` is deterministic. Per-iteration seed is `seed + iteration_index`; bind rejects `seed + n_subsamples - 1 > INT_MAX`.

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
- `sample_id` (VARCHAR): sample identifier (post-subsample ordering)
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
- `sample_id` (VARCHAR): sample identifier
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
