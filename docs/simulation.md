# Community simulation

Generate synthetic OTU/feature tables with known ground truth, for benchmarking resemblance, ordination, and hypothesis-testing methods (β-diversity, PCoA, PERMANOVA). The two simulators reproduce the non-phylogenetic models of Kuczynski et al. 2010 (*Nature Methods* 7:813–819): an **environmental gradient** (species with unimodal response curves along a 1-D gradient) and **discrete clusters** (samples drawn around perturbed cluster centroids). Both take an empirical (or synthetic) relative-abundance vector as input and emit a long/COO table of nonzero counts labelled with the generating ground truth.

## Table of Contents

- [Gradient communities](#gradient-communities) - Samples along a 1-D environmental gradient.
- [Clustered communities](#clustered-communities) - Samples grouped into discrete clusters.
- [Supplying the abundance vector](#supplying-the-abundance-vector) - Passing a `LIST(DOUBLE)` argument.

### Gradient communities

Each species is assigned a random optimum position drawn uniformly on `[0, 1]` and a Gaussian response curve of fixed width `sp_width`, scaled by its input relative abundance. The gradient is sampled at `num_samples` evenly spaced positions on `[range_lo, range_hi]`; each sample's per-species profile is optionally perturbed by noise, then `seqs_per_sample` reads are drawn multinomially. Species absent from every sample are dropped (the output is nonzero-only).

**Function signature**:

`simulate_gradient_otus(abundances, num_samples, seqs_per_sample, [options])`

**Parameters:**
- `abundances` (LIST(DOUBLE), required, positional): Per-species relative abundances (peak heights). Must be non-empty, finite, non-negative, and not all zero. See [Supplying the abundance vector](#supplying-the-abundance-vector).
- `num_samples` (INTEGER, required, positional): Number of samples along the gradient (≥ 1).
- `seqs_per_sample` (BIGINT, required, positional): Sequencing depth — reads drawn per sample (≥ 1).
- `sp_width` (DOUBLE, default `0.1`): Gaussian width (σ) of every species' response curve; must be finite and positive.
- `noise` (DOUBLE, default `0.0`): Noise magnitude. `0.0` = unperturbed. The perturbation standard deviation is `noise ×` the width source selected by `noise_type`.
- `noise_type` (VARCHAR, default `'+species'`): selects where each perturbation's standard deviation comes from — a species' own abundance, or the sample total shared across all species.
  - `'+species'` (default): adds `Normal(0, noise·x)` to each species' abundance `x`. The width is per species, so the perturbation is a relative jitter of the same proportion everywhere, and a species with zero abundance takes zero width — the perturbation step leaves it at zero. Note this does not by itself guarantee such a species stays absent from the output: the floor-shift below lifts the whole profile whenever some other species is driven negative, which becomes likely as `noise` approaches (or exceeds) 1.
  - `'*sample'`: multiplies each species by `Normal(1, noise·Σ)`, where `Σ` is the sample's abundance sum.
  - `'+sample'`: adds `Normal(0, noise·Σ)` to each species.

  After perturbation the profile is shifted so its minimum is 0 (a uniform per-sample shift — not a per-species clip; there is no clip option for gradient noise) and then rescaled to preserve the sample's mean abundance. Only validated when `noise ≠ 0`.
- `range_lo` (DOUBLE, default `0.1`), `range_hi` (DOUBLE, default `0.9`): Gradient sampling interval; both finite with `range_lo ≤ range_hi`.
- `seed` (BIGINT, default `-1`): RNG seed. `-1` = nondeterministic; `≥ 0` = reproducible.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | INTEGER | 0-based sample index (ordered along the gradient) |
| `otu_id` | INTEGER | 0-based species index into `abundances` |
| `count` | BIGINT | Reads assigned to this (sample, OTU); always ≥ 1 (nonzero-only) |
| `gradient_position` | DOUBLE | The sample's position on the gradient (`range_lo`…`range_hi`) |

```sql
-- A prominent gradient (low noise): PCoA on Bray-Curtis should recover position.
SET VARIABLE ab = (SELECT list(1.0/(i+1) ORDER BY i) FROM range(150) t(i));
CREATE TABLE grad AS
  SELECT * FROM simulate_gradient_otus(getvariable('ab'), 40, 1000, noise := 0.5, seed := 42);

-- Per-sample depth is exactly seqs_per_sample.
SELECT sample_id, sum("count") FROM grad GROUP BY sample_id;
```

### Clustered communities

Starting from the input vector (renormalized to sum 1), one **centroid** is built per cluster by perturbing the base by magnitude `cluster_spacing`; each within-cluster sample is then drawn by perturbing its centroid by magnitude `sample_spacing`. Each perturbation multiplies (or adds) `Normal` noise, clips/rescales to non-negative, and renormalizes to sum 1; then `seqs_per_sample` reads are drawn multinomially. Species absent from every sample are dropped.

**Function signature**:

`simulate_cluster_otus(abundances, seqs_per_sample, [options])`

**Parameters:**
- `abundances` (LIST(DOUBLE), required, positional): Base per-species relative abundances (as above).
- `seqs_per_sample` (BIGINT, required, positional): Sequencing depth per sample (≥ 1).
- `cluster_sizes` (LIST(INTEGER), default `[30, 30, 30]`): Number of samples in each cluster; non-empty, each ≥ 1. The number of clusters is the list length; total samples is the sum.
- `cluster_spacing` (DOUBLE, default `1.0`): Perturbation magnitude from base → cluster centroid (larger = clusters further apart). Finite, ≥ 0.
- `sample_spacing` (DOUBLE, default `0.5`): Perturbation magnitude from centroid → sample (larger = looser clusters). Finite, ≥ 0.
- `noise_type` (VARCHAR, default `'*sample'`): `'*sample'` (multiplicative) or `'+sample'` (additive).
- `normalization` (VARCHAR, default `'clip'`): `'clip'` floors negatives at 0; `'rescale'` shifts the minimum to 0. Both then renormalize to sum 1.
- `seed` (BIGINT, default `-1`): RNG seed (`-1` = nondeterministic).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | INTEGER | 0-based sample index (grouped by cluster, in `cluster_sizes` order) |
| `otu_id` | INTEGER | 0-based species index into `abundances` |
| `count` | BIGINT | Reads assigned to this (sample, OTU); always ≥ 1 (nonzero-only) |
| `cluster_id` | INTEGER | 0-based generating cluster of the sample |

```sql
-- Prominent (1.0/0.5) vs subtle (0.1/0.1) clustering, K=3 clusters of 30.
SET VARIABLE ab = (SELECT list(1.0/(i+1) ORDER BY i) FROM range(150) t(i));
CREATE TABLE clust AS
  SELECT * FROM simulate_cluster_otus(getvariable('ab'), 1000,
                                      cluster_spacing := 1.0, sample_spacing := 0.5, seed := 42);

-- Recover the 3 clusters and their sizes.
SELECT cluster_id, count(DISTINCT sample_id) AS n_samples FROM clust GROUP BY cluster_id;
```

### Supplying the abundance vector

Both functions take `abundances` as a `LIST(DOUBLE)` positional argument. DuckDB table-function arguments must be **constant-foldable** (a subquery raises `Table function cannot contain subqueries`), so build the list into a session variable first and pass it via `getvariable`:

```sql
-- From a table of (idx, abund) rows:
SET VARIABLE ab = (SELECT list(abund ORDER BY idx) FROM my_abundances);
SELECT * FROM simulate_gradient_otus(getvariable('ab'), 20, 1000, seed := 42);
```

A small vector can also be passed as a literal, e.g. `[0.4, 0.25, 0.2, 0.1, 0.05]::DOUBLE[]`.

**Reproducibility.** With `seed ≥ 0` the output is deterministic for a given build. Sample statistics (richness, distance-decay, cluster separation) are stable across builds, but exact per-cell counts are **not** portable across C++ standard-library implementations (the standard does not fix the output of `std::normal_distribution`/`std::discrete_distribution`); compare simulations by statistical summary, not byte-for-byte.
