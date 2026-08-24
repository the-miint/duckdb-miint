# Diversity calculations and statistics

Methods to estimate alpha and beta diversity, and supporting statistics.

## Table of Contents

- [Feature table](#feature-table) - Feature table expectations.
- [Tree](#tree) - Tree expectations.
- [Metadata](#metadata) - Metadata expectations.
- [Sample identifier types](#sample-identifier-types) - How VARCHAR/BIGINT/UUID sample ids are handled.
- [UniFrac algorithm variants](#unifrac-algorithm-variants) - Detail on different UniFrac algorithms and how to specify them.
- [Rarefaction](#rarefaction) - Rarefaction detail with UniFrac and Faith PD.
- [Rarefy a feature table](#rarefy-a-feature-table) - produce an even-depth feature table as a standalone result
- [Pick anchors](#pick-anchors) - subset selection from ordination coordinates: proportional stratified sampling (the measured best rule for progressive-PCoA anchors) or greedy farthest-point
- [UniFrac distances](#unifrac-distances) - Condensed (pairwise) UniFrac distances in long form
- [Community distances (non-phylogenetic)](#community-distances-non-phylogenetic) - taxon-based β-diversity (Bray-Curtis, Jaccard, Morisita-Horn, χ², Gower, …) from a feature table
- [Beta-distance macros](#beta-distance-macros) - within/between-group distributions and k-nearest-neighbors over a distance table
- [PCoA (from a distance table)](#pcoa-from-a-distance-table) - metric-agnostic PCoA over any condensed distance table
- [Progressive PCoA (from a distance table)](#progressive-pcoa-from-a-distance-table) - scalable reference-anchored PCoA without a dense N×N decomposition
- [Progressive PCoA (from UniFrac)](#progressive-pcoa-from-unifrac) - scalable UniFrac PCoA computing distances on the fly per batch (true-10M path)
- [Progressive PCoA (from community distances)](#progressive-pcoa-from-community-distances) - the same, tree-free: block-wise Bray-Curtis/Jaccard/… straight from a feature table
- [PERMANOVA (from a distance table)](#permanova-from-a-distance-table) - metric-agnostic PERMANOVA over any condensed distance table
- [Two-sample Kolmogorov-Smirnov test](#two-sample-kolmogorov-smirnov-test) - `ks_2samp` over two numeric lists, with SciPy-exact p-values
- [Sample clustering (k-means and UPGMA)](#sample-clustering-k-means-and-upgma) - group samples from ordination coordinates or a distance table
- [UniFrac PCoA](#unifrac-pcoa) - UniFrac distance + Principal Coordinates Analysis
- [UniFrac PERMANOVA](#unifrac-permanova) - UniFrac distance + PERMANOVA pseudo-F + p-value
- [Faith PD](#faith-pd) - Faith's phylogenetic diversity per sample
- [Procrustes](#procrustes-align-two-ordinations) - align two ordinations into a common frame (disparity M² + PROTEST p-value)
- [Citations](#citations) - primary sources for the metrics, clustering, and ordination methods

Most of these methods are powered by the embedded [`unifrac-binaries`](https://github.com/biocore/unifrac-binaries) and [`scikit-bio-binaries`](https://github.com/scikit-bio/scikit-bio-binaries) libraries (see `docs/internals/embedded-tools.md` for the build details); [`procrustes`](#procrustes-align-two-ordinations) is a self-contained Eigen-backed port of SciPy and works on any ordination table.

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

<a name="subsample-thread-count"></a>**A seeded subsample reproduces per thread count, not across thread counts.** libssu distributes the draw across the OpenMP team — one generator per thread, seeded in turn — so both how many generators exist and which observation consumes which one depend on the team size. The OpenMP width comes from DuckDB's `SET threads = N` (or the per-call `threads` parameter), so the *same* query with the *same* `seed` and `subsample_depth > 0` can produce a different rarefied table, and therefore different distances, on a differently configured server.

Widths that happen to divide the work identically can coincidentally agree — in one measured case widths 1, 2, and 4 agreed while 8 differed — so matching results at two thread counts is not a guarantee of matching at a third. If you need a rarefied result to be reproducible across machines, fix the width explicitly (`threads := N`, or `SET threads = N`) alongside the seed, and record it with the result. This is a property of the draw itself and applies to every function in this section; it is upstream of, and independent of, the fp32 reduction-ordering noise described under [Reproducibility](#reproducibility).

---

### Pick anchors

Select `n_anchors` of the samples in an **ordination coordinate table** — the `(sample_id, axis, coordinate)` long form emitted by [`pcoa`](#pcoa-from-a-distance-table) and [`unifrac_pcoa`](#unifrac-pcoa), the same contract [`cluster_kmeans`](#sample-clustering-k-means-and-upgma) consumes — returning `(anchor_rank, sample_id)` in selection order.

```sql
-- ordinate a candidate subset, pick anchors from those coordinates, anchor the full run
CREATE TABLE cand AS SELECT sample_id, axis, coordinate FROM pcoa('sub_dm', n_dims := 10, seed := 42);
SET VARIABLE anch = (SELECT list(sample_id) FROM pick_anchors('cand', n_anchors := 1000, seed := 42));
SELECT * FROM progressive_pcoa_from_unifrac('ft', 'tree', variant := 'unweighted',
    n_dims := 10, seed := 42, anchors := getvariable('anch'));
```

Working from coordinates rather than from a distance matrix is what makes this usable at scale: both rules are **linear in the number of samples**, so a 10M-sample selection needs an N×d table, not an N×N one.

**Parameters:**
- `coordinate_table` (VARCHAR): name of the relation exposing `(sample_id, axis INTEGER, coordinate DOUBLE)`
- `n_anchors` (INTEGER, **required**): how many samples to select. No default — there is no defensible universal anchor count.
- `method` (VARCHAR, default `'stratified'`): `'stratified'` or `'farthest_point'`; see below.
- `seed` (BIGINT, default 0): salts the within-stratum draw. `'stratified'` only — `'farthest_point'` takes no seed.
- `n_dims` (INTEGER, default 3): use the leading `n_dims` axes by ascending axis label; `0` uses every axis present. For `'stratified'` these are the axes the strata grid is built on; for `'farthest_point'` they are the axes distances are measured over.
- `n_bins` (INTEGER, default 4): equal-frequency bins per axis, so the grid holds up to `n_bins ^ n_dims` strata (64 at the defaults). `'stratified'` only.

**Output:** `(anchor_rank INTEGER, sample_id)`; `sample_id` mirrors the input id type. Selection is **prefix-stable** under both methods: the first *m* of *k* anchors are exactly the *k = m* result, so you can pick a generous `k` once and trim.

#### `method := 'stratified'` (default) — proportional stratified sampling

Bin each of the leading `n_dims` axes into `n_bins` equal-frequency bins, take the resulting grid cell as the stratum, and draw from every stratum in proportion to its size. Within a stratum the order is a salted hash of the sample id, so the draw carries no geometric preference; across strata, samples are ordered by `(within-stratum rank) / (stratum size)`, so taking the first `k` draws proportionally and totals exactly `k`. Classic survey sampling — Neyman 1934, Cochran 1977 ch. 5.

**This is the measured best rule for progressive-PCoA anchors.** Several literature-backed rules were compared on the rarefied EMP 90 bp table (23,814 samples, `unweighted`, `n_anchors := 1000`, `batch_size := 1000`, `n_dims := 10`), scored by procrustes M² against a full `unifrac_pcoa` at d=3, all drawing from an identical 4,000-sample candidate pool so the rule is the only difference. Lower is better:

| rule | literature | M² vs the full ordination (d=3) |
|---|---|---|
| **proportional stratified** | survey sampling (Neyman 1934; Cochran 1977) | **0.0051** (mean of 15 draws, 0.0038–0.0079) |
| plain seeded random (the built-in default) | — | 0.0113 (mean of 10 draws, 0.0042–0.0220) |
| stratum medoids | k-medoids / PAM (Kaufman & Rousseeuw 1990) | 0.0176 |
| leverage-proportional | CUR / column subset selection (Drineas & Mahoney) | 0.0195 |
| farthest-point | k-center (Gonzalez 1985) | 0.11 (over coordinates) / 0.17 (over full distances) |

Stratified is **2.19× better on the mean**, and — the more useful property — its spread is **5.3× tighter** (sd 0.0011 vs 0.0061). That is what you are buying: random anchoring is not reliably bad, it is *unreliable*. Its best draw of ten (0.0042) beat stratified's median, but its worst (0.0220) was nearly 3× worse than stratified's worst, and **no stratified draw in fifteen was worse than random's median** while only one random draw in ten beat stratified's median. Over all 150 pairings a stratified draw beats a random one 84% of the time. The ranges do overlap, so a single lucky random draw can match it; the point is that you cannot tell in advance which draw you got, and stratified removes most of that exposure.

The mechanism is worth stating, because it is what picks this rule out of the five: medoids, leverage and farthest-point each systematically prefer a *kind* of point — central, high-influence, extreme — and **all three lose to an unbiased draw**. [`progressive_pcoa_from_unifrac`](#progressive-pcoa-from-unifrac) builds its reference frame from a PCoA *of the anchors*, so the anchors' leading axes have to match the full ordination's, and any rule that biases *which* points are chosen rotates that frame away from the truth. Stratified sampling is the only one of the four that stays unbiased *within* a stratum and merely equalizes coverage *across* strata. Measured directly on the same pool: the fraction drawn per grid cell has sd 0.058 under stratified against 0.092 under random, around a 25% target. So: don't bias which points are chosen, only reduce the lumpiness of where they land.

Keep the number of non-empty strata well below `n_anchors` — that is what makes the allocation proportional. With far more strata than anchors the draw degenerates gracefully toward simple random sampling over the most-populated cells, which is the baseline, not a failure.

#### `method := 'farthest_point'` — greedy k-center

Repeatedly take the sample whose nearest already-selected neighbour is farthest away (Gonzalez 1985), a 2-approximation to minimizing the maximum distance from any sample to its nearest selected one. Deterministic and seedless: rank 0 is the most peripheral sample and every tie breaks to the lowest sorted id, so a selection is a reproducible property of the data. "Most peripheral" is the sample farthest from the centroid, which in Euclidean space is exactly the sample with the largest total squared distance to all others — the same rule a distance-matrix implementation reaches through row sums, at O(N·d) instead of O(N²·d).

> **Do not use this for progressive-PCoA anchors** — that is what it was built for, and it measured **an order of magnitude worse than either the seeded random default or `'stratified'`**: M² of 0.11 selecting over the leading 3 coordinate axes, and 0.17 in an earlier formulation selecting over the full distance matrix, against a random band of 0.004–0.022 on the same pool. It is not simply outlier-picking: mean distance-to-all was 0.869 for these anchors against 0.861 for random and a 0.862 pool mean. A reference frame must be **representative**, not merely **covering**.

It is kept because farthest-point selection is a legitimate, well-studied rule for **diverse subset selection** generally — a maximally spread reference panel, a coverage-based subsample, a diverse review set. Use it for those.

**Cost:** `'stratified'` is O(N·d·log N); `'farthest_point'` is O(N·d) then O(N·k·d). Neither materializes a distance matrix. To pass the result into a function that takes `anchors :=`, hand it through a variable — a table function's named parameter cannot contain a subquery.

---

### Rarefy a feature table

Rarefy a feature table to an even per-sample depth and return **the rarefied table itself** — `(sample_id, feature_id, value)`, the same shape every function on this page consumes.

```sql
-- rarefy once, then reuse that one table for every downstream analysis
CREATE TABLE ft_even AS
    SELECT * FROM rarefy_feature_table('observations', depth := 1000, seed := 42, threads := 8)
    ORDER BY sample_id;

SELECT * FROM progressive_pcoa_from_unifrac('ft_even', 'tree', variant := 'unweighted', n_dims := 10);
```

**Why this exists separately** from the `subsample_depth` parameter on the other functions: those rarefy *internally* and never hand back the table, so the drawn table cannot be inspected, shared between analyses, or fed to a function that has no `subsample_depth` of its own. [`progressive_pcoa_from_unifrac`](#progressive-pcoa-from-unifrac) is exactly that case — it deliberately omits `subsample_depth` because each batch would otherwise draw independently — so rarefying through this function is the only way to give it an even-depth table. It also lets you draw **once** and reuse it, so a PCoA, a PERMANOVA, and a Faith PD all describe the same rarefied data rather than three different draws.

**Parameters:**
- `feature_table` (VARCHAR): name of the feature relation exposing `(sample_id, feature_id, value)` (see [Feature table](#feature-table))
- `depth` (INTEGER, **required**): target per-sample total count. No default — there is no defensible universal depth, and a silent one would change every downstream result.
- `with_replacement` (BOOLEAN, default `false`): `true` draws multinomially (a cell may exceed its input count); `false` permutes without replacement.
- `seed` (INTEGER, default `-1`): `-1` uses system entropy; any `seed >= 0` is deterministic — subject to the thread-count caveat below.
- `threads` (INTEGER, default 0): OpenMP threads; `0` follows DuckDB's thread count.

The parameter is `depth` rather than the `subsample_depth` used elsewhere on this page: the function name already supplies that context, so `subsample_depth` would be redundant here.

**Output schema:** `(sample_id, feature_id, value DOUBLE)`. **Both** id columns mirror their input types (BIGINT/UUID preserved, everything else VARCHAR — see [Sample identifier types](#sample-identifier-types)), so a rarefied table still joins to typed metadata and taxonomy exactly as its input did.

**Behavior:**
- **Every surviving sample comes out at exactly `depth`.** That is the point: without it, sequencing depth acts as a covariate, and unweighted UniFrac in particular is strongly depth-sensitive.
- **Samples whose total count is below `depth` are DROPPED** (never padded or partially filled), matching `subsample_depth` on the other functions. Features left with no nonzero count anywhere are dropped too, preserving the sparse-storage invariant. Check what you lost before trusting the result: compare `count(DISTINCT sample_id)` against the input.
- **If no sample reaches `depth`, that is an error**, not an empty table — an empty result reads too easily as a successful analysis of nothing.
- **Reproducibility:** a seeded draw reproduces per thread count, not across thread counts — see [the caveat above](#subsample-thread-count). Fix `threads := N` alongside the seed if the rarefied table has to be reproducible elsewhere.
- **Rarefy once, upstream.** Storing the result (as in the example) is preferable to calling this inline in several queries: each call is a fresh draw unless you pin both `seed` and `threads`.

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

| metric | formula | range | block-wise? | notes | primary source |
|---|---|---|---|---|---|
| `bray_curtis` | `Σ\|xₖ−yₖ\| / Σ(xₖ+yₖ)` | [0,1] | ✅ | empty pair → 0 | Bray & Curtis 1957 |
| `euclidean` | `sqrt(Σ(xₖ−yₖ)²)` | [0,∞) | ✅ | | |
| `jaccard` | binary presence/absence `(b+c)/(a+b+c)` | [0,1] | ✅ | **presence/absence**, not abundance; empty pair → 0 | Jaccard 1912 |
| `soergel` | `Σ\|xₖ−yₖ\| / Σ max(xₖ,yₖ)` | [0,1] | ✅ | empty pair → 0 | |
| `morisita_horn` | `1 − 2Σ(xₖyₖ) / ((Σxₖ²/X² + Σyₖ²/Y²)·X·Y)` | [0,1] | ✅ | Horn's Cλ on relative abundances; both-empty → 0, one-empty → 1 | Morisita 1959; Horn 1966; Magurran 2004 p.246 |
| `pearson` | `1 − r` over features | [0,2] | ❌ | the per-sample mean is `Σx / n_features`, and features zero in *both* samples still enter the covariance — so the value moves with the feature space | |
| `chisq` | `sqrt(Σₖ (GT/colₖ)(xₖ/X − yₖ/Y)²)` | [0,∞) | ❌ | needs column sums and the grand total | Faith, Minchin & Belbin 1987 |
| `gower` | `Σₖ \|xₖ−yₖ\| / rangeₖ` | [0,∞) | ❌ | needs each feature's range over all samples | Gower 1971; Faith, Minchin & Belbin 1987 |

**"Block-wise?"** marks the metrics that [`progressive_pcoa_from_features`](#progressive-pcoa-from-community-distances) accepts. A ✅ metric is *pairwise-local*: `d(i,j)` depends only on samples `i` and `j`, reads no statistic taken over the other rows, and is unchanged by dropping features that are zero in both — so it can be computed one block of samples at a time and get exactly the answer the full matrix would give. The three ❌ metrics read the whole matrix, so a per-block value would silently be a *different distance in every block*; they are refused at bind rather than computed. All eight remain available from `community_distances` itself, which always sees the whole table.

The zero-variance and zero-row-sum conventions (the `pearson` and `chisq` notes
above) follow `cogent3.maths.distance_transform`, the reference implementation of
the metrics used by Kuczynski et al. 2010 — see [Citations](#citations).

**Behavior:**
- **Raw values, no pre-normalization.** Abundances are used exactly as given; each metric applies whatever internal normalization its own definition requires. Feeding counts versus relative abundance is your modeling choice. (At equal per-sample depth the two differ only by a constant factor for most metrics, which does not change ordination geometry.)
- **Matrix-wide metrics.** `chisq` and `gower` depend on *global* column statistics (column sums and column ranges across all samples), so a pair's distance is a function of the whole matrix, not just that pair. Subsetting samples changes them.
- **Sparse contract:** a sample whose values are all zero or NULL is dropped before computation, matching `unifrac_distances`. Distances are therefore emitted only among samples that actually carry signal.
- **Threading is exact:** results are **bit-identical for any thread count** (each pair writes to a fixed condensed slot), so `threads` is purely a performance knob.
- <a name="community-distances-dense-kernel"></a>**Dense inputs switch kernels automatically.** The default kernel merges each pair's *nonzeros*, which is the right shape for a microbiome table (typically far under 1% of cells occupied). On a genuinely dense matrix that is the wrong trade, so `euclidean`, `jaccard` and `morisita_horn` switch to a Gram-matrix formulation above ~1.5% density — each is a function of one inner product plus per-sample scalars, so `n²` distances come out of a single matrix product. Measured on 20,000 samples × 2,000 features, one thread: **2.1× at 2% density, 5.1× at 5%, 13× at 15%, 25× at 40%**, and 9.7× on a 38.8%-dense image matrix with only 673 features. `bray_curtis` and `soergel` never switch — `Σ|x−y|` and `Σ max(x,y)` are not inner products, so no such formulation exists for them — and neither do matrices small enough (under 64 samples) or wide enough that the dense operand would exceed 512 MB. The switch does cost memory — the dense operand, plus a sample×sample Gram matrix that is twice the condensed result you are already asking for (measured: +29 MB on a 2,000-sample block, against the 32 MB that predicts) — so it is a bounded factor on an allocation you have already chosen, not a new risk. It changes speed, not results: distances agree with the merge kernel to within floating-point last bits (exactly, for `jaccard`), and both entry points dispatch on the same rule, so the sparse and dense APIs stay bit-identical to each other.
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
- `global_rotation` (BOOLEAN, default `true`): express the finished configuration in its own principal axes rather than the anchors'. See [Global rotation](#progressive-global-rotation) — it cannot change accuracy, but it costs streaming and the ability to extend the ordination later

**Output schema:** [`pcoa`](#pcoa-from-a-distance-table)'s six columns — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained)` — plus two appended diagnostic columns, `(batch, batch_anchor_m2)`. `iteration` is always `0`. `sample_id` mirrors the input `sample_a` type — see [Sample identifier types](#sample-identifier-types). **Caveat:** `eigenvalue` and `proportion_explained` are the *anchor* reference ordination's (they describe the anchor subspace, not the full sample set); the per-sample `coordinate`s span all samples.

<a name="progressive-global-rotation"></a>**`global_rotation` — whose axes are these?** Without it the reference frame is the *anchors'* principal axes and every sample is expressed in them. That is a consistent frame, but not the one a PCoA promises: the anchors are a sample of the data, so their leading axis only *estimates* the full configuration's and is off by however much the draw was unlucky. The symptom is that axis 0 of a progressive run and axis 0 of a full [`pcoa`](#pcoa-from-a-distance-table) are close but not the same direction — so per-axis work (plotting PC1 against PC2, correlating an axis with metadata, comparing two runs axis by axis) disagrees by more than the ordinations themselves do.

With `global_rotation := true` (the default) the run finishes by centring the assembled configuration and rotating it onto its own covariance eigenvectors. The axes then mean what they mean in `pcoa`: axis 0 is the direction of greatest variance *in the returned coordinates*, axis variances are strictly descending, and every axis has mean 0.

**Measured against a full `pcoa`.** Three 5,000-sample tables (a rarefied microbiome table under `bray_curtis` and `jaccard`, and a dense image matrix under `euclidean`), scored as the mean |correlation| between each returned axis and the *same-numbered* axis of a full `pcoa` — **without** procrustes alignment, because aligning first would rotate the two into agreement and hide the very thing being measured:

| table / metric | `n_anchors` | mean \|r\|, axes 0–2 (on → off) |
|---|---|---|
| image / `euclidean` | 100 | **0.998** → 0.668 |
| image / `euclidean` | 200 | **0.998** → 0.652 |
| image / `euclidean` | 1000 | 0.999 → 0.987 |
| microbiome / `bray_curtis` | 100 | **0.919** → 0.496 |
| microbiome / `bray_curtis` | 200 | 0.751 → **0.927** |
| microbiome / `bray_curtis` | 1000 | **0.969** → 0.694 |
| microbiome / `jaccard` | 100 | **0.970** → 0.885 |
| microbiome / `jaccard` | 200 | **0.983** → 0.891 |
| microbiome / `jaccard` | 1000 | 0.995 → 0.989 |

Better in eight of the nine, and **the gain is largest exactly where anchors are scarce** — the regime this function exists for. With 1,000 anchors on 5,000 samples the anchor frame is already close to the right frame, and rotating buys almost nothing.

**It is not a per-axis guarantee, and the one inversion above shows why.** `bray_curtis` on that table has adjacent axis variances in the ratio 0.94 (axes 1–2) and 0.99 (axes 8–9). When two axes carry nearly the same variance, *which* of them is "axis 1" is not determined — not for this function and not for a full `pcoa` either, since any rotation inside a degenerate eigenspace is equally principal. Where the spectrum is flat, compare subspaces (or use [`procrustes`](#procrustes-align-two-ordinations)), not individual axes.

**It cannot change accuracy — only which coordinate carries which structure.** A rotation about the centroid is a rigid motion and procrustes M² is invariant to rigid motions, so the disparity against a full `pcoa` is identical either way. Measured on the same run: M² between the rotated and unrotated results `4.8e-30`, and all pairwise distances agreeing to `5.6e-17`. And in all nine configurations above, the disparity against the full `pcoa` is the same with the flag on and off to 8–16 significant digits — the axes move, the ordination does not. If a comparison you care about is frame-invariant, this setting is a no-op for it.

**Two costs, both structural rather than numerical:**
- **The run stops streaming.** The rotation is a property of the whole configuration, so every coordinate is staged before any row is emitted, and the first row arrives at the *end* of the run rather than after the first batch. Concretely, `SELECT * FROM progressive_pcoa_from_features('t', 'bray_curtis') LIMIT 10` used to return as soon as the anchor frame existed; with rotation on it pays for the whole run first. If that matters — interactive exploration, inspecting the first batch before committing to the rest — pass `global_rotation := false`. Staging is buffer-managed, so `memory_limit` bounds it and it spills to `temp_directory` instead of growing unbounded; the volume is small anyway (about 100 bytes per sample plus the id, at `n_dims := 10`). **The wall-clock cost is not the issue** — measured by running the same input, seed and binary with the flag on and off, interleaved: **1.2% at 200,000 samples**, and 0–6% elsewhere with peak RSS unchanged. What changes is *when* the first row arrives, not how long the run takes.
- **A finished ordination can no longer be extended.** Adding samples later changes the configuration, hence its principal axes, hence the frame — so coordinates returned earlier would no longer be comparable to the new ones. If you intend to place new samples into an existing ordination, use `global_rotation := false`: the anchor frame is fixed by the anchors alone and is therefore stable across runs, which is exactly the property that makes it extensible.

`eigenvalue` and `proportion_explained` remain the *anchor* ordination's in both modes (the caveat above) — they are not re-derived from the rotated configuration, so do not read them as the variance of the axes actually returned.

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
- **Reproducibility:** a seeded run is byte-reproducible only at `threads := 1`; above that, two runs of the same query can differ in every coordinate by a last-bit frame shift. Compare ordinations with [`procrustes`](#procrustes-align-two-ordinations), not coordinates — the full explanation is under [`progressive_pcoa_from_unifrac`](#progressive-repro) and applies unchanged here, since the effect is in the ordination rather than in how blocks are obtained.
- **Completeness (fail loud):** each batch needs its full `(anchors + batch)²` block present in the distance relation; a missing pair within a block is an error naming the offending pair. NULL sample ids and NULL/NaN distances are skipped. `sample_a`/`sample_b` must resolve to the same output type (a BIGINT/VARCHAR mix is rejected at bind).
- <a name="progressive-cancellation"></a>**Cancellable:** cancellation is checked before each wave's scan and before every batch, so Ctrl-C (or a client cancel) stops the run within about one block — or, if it lands while a wave's scan is running, once that scan finishes. Rows already returned are discarded with the rest of the query, as for any cancelled statement.
- **Fewer than two samples**, `n_anchors` outside `[n_dims + 1, n_samples]`, `n_dims > n_samples - 1`, or `batch_size < 1`: an error. Parameter and schema problems are reported when the query is bound; problems in the *data* (an incomplete block, conflicting distances) surface once the run reaches the batch that hits them.

**Example:**

```sql
-- Ordinate a large precomputed distance table without a dense N×N decomposition,
-- then confirm it matches a full pcoa (small procrustes disparity) on a subset.
SELECT sample_id, axis, coordinate
FROM progressive_pcoa_from_distances('dm', n_dims := 3, n_anchors := 100, batch_size := 1000, seed := 42)
ORDER BY axis, sample_id;
```

> **Note:** blocks are filled a *wave* at a time — one pass over the relation serves many blocks, instead of one scan per batch — and the wave's width is chosen from what `memory_limit` currently leaves free. The run is driven from the scan rather than computed up front, so each batch's rows are handed to the query as that batch finishes and live memory is one wave of blocks plus that wave's rows, never the whole result.

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
- `n_dims` (3), `n_anchors` (100), `batch_size` (1000), `seed` (-1), `threads` (0), `global_rotation` (`true`), `anchors`: as in [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table) — `n_anchors` must be in `[n_dims + 1, n_samples]`. `anchors := [...]` supplies the anchor set explicitly and takes precedence over `n_anchors`/`seed`; note the seed still seeds each block's ordination, so a seeded run is required for reproducibility even with explicit anchors. The seeded random default is a sound choice; the one rule measured to beat it is proportional stratified sampling via [`pick_anchors`](#pick-anchors), and that section also records which rules are *worse* than random and why.
- `variant`, `variance_adjust`, `alpha`, `bypass_tips`, `normalize_sample_counts`: the UniFrac controls, identical to [`unifrac_pcoa`](#unifrac-pcoa)

There is deliberately **no `subsample_depth`**: rarefaction and progressive alignment do not compose cleanly (each batch would rarefy independently against a different RNG draw). Rarefy upstream instead with [`rarefy_feature_table`](#rarefy-a-feature-table) and pass the resulting table here — worth doing, since unweighted UniFrac is strongly depth-sensitive and an uneven-depth table makes sequencing depth a covariate of the ordination.

**Output schema:** [`unifrac_pcoa`](#unifrac-pcoa)'s six columns — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained)`, `iteration` always `0` — plus the appended `(batch, batch_anchor_m2)` diagnostics, identical in meaning to `progressive_pcoa_from_distances`' (see [Reading `batch` / `batch_anchor_m2`](#progressive-batch-diagnostics)). As with `progressive_pcoa_from_distances`, `eigenvalue`/`proportion_explained` are the *anchor* reference ordination's (a documented caveat). `sample_id` mirrors the input type — see [Sample identifier types](#sample-identifier-types).

**Behavior:**
- **Accuracy / batch-invariance:** same guarantees as [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table) — reproduces a full [`unifrac_pcoa`](#unifrac-pcoa) up to a similarity transform, with alignment error that does not compound across batches.
- <a name="progressive-unifrac-threads"></a>**`threads` controls block CONCURRENCY here, and it is the main performance knob.** Unlike `unifrac_pcoa`, a single block cannot use more than one core no matter how high you set it: libssu's UniFrac parallel degree is `ceil(n_samples / 2048)` over the *block*, and a block is only `n_anchors + batch_size` samples, so at any ordinary setting that is one stripe. So `threads` is spent running that many *blocks* at once instead, each pinned narrow. On 17,483 samples (`n_dims := 5`, `n_anchors := 200`, `batch_size := 500`) this took roughly **250 s at `threads := 1` versus ~32 s at `threads := 14`** — about **8×**, moving from ~1.1 to ~9.8 cores busy. (Timings here and below are indicative: they come from a shared machine, so treat the ratios and the ordering as the result, not the exact seconds.)

  **Expect the returns to stop well before your core count.** Measured 1 → 14 threads on 23,814 EMP samples (`n_dims := 10`, `n_anchors := 1000`, `batch_size := 1000`) on a 14-core host: **142 s / 74 s / 41 s / 25 s / 24 s** at 1 / 2 / 4 / 8 / 14 threads — that is 1.9× / 3.5× / 5.8× / 6.0×, so **96%, 88%, 72% and 43%** of linear. The last step buys **1.04× for 1.75× the threads**. Blocks are memory-bandwidth-bound, so past a point the extra workers contend rather than add: the same block took 8.45 s when fourteen ran at once and 6.23 s when nine did. Set `threads` to roughly the point where your own curve flattens (~8 on that host) rather than to the core count, and leave the rest for whoever else is on the machine.
- <a name="progressive-unifrac-batch-count"></a>**With few batches, `threads` above your core count can still help — nothing else will.** Blocks are handed to workers from one pool, so a worker takes its next block the moment it frees; blocks are uneven (each shears the tree to its own features — 4.4 s to 9.5 s on EMP), and pulling from a pool is what keeps that unevenness from idling anyone. What no scheduler can remove is that *k* indivisible blocks over *t* workers take `ceil(k/t)` rounds: on 23,814 EMP samples `batch_size := 1000` is 23 batches, so at `threads := 14` the second round runs nine blocks on fourteen workers. Since a block is one core here, raising `threads` past the core count is what closes that — measured ~19.4 s at `threads := 14` versus ~17.2 s at `threads := 28`, same result — at the cost of holding one block per concurrent batch, so only do it on a machine you are not sharing and with the memory to spare. Adjusting `batch_size` to make the batch count a multiple of `threads` is **not** worth it: `batch_size := 815` (28 batches, two exact rounds) measured no faster than 1000 (~18.6 s vs ~18.9 s). All of this fades as the run grows; by a few hundred batches the last round is under 1% of the work.
- <a name="progressive-unifrac-batch-size"></a>**Tune `m = n_anchors + batch_size`, and keep it small — the opposite of what a serial run wants.** A block's cost grows with `m` faster than linearly, and once blocks run concurrently it is total CPU, not per-block latency, that sets wall time; smaller blocks also mean *more* of them, which is what there is to parallelize. Serially the ordering reverses (big blocks amortize a fixed per-block cost), so a `batch_size` tuned on a one-core run is the wrong one here. The optimum tracks `m`, not `batch_size` alone — at `n_anchors := 200`, `batch_size` 244 → ~31 s, 500 → ~34 s, 1000 → ~56 s, 1849 → ~73 s, 3000 → ~103 s (17,483 samples, `threads := 14`), with total CPU rising ~3.5× across that range. Total CPU is the more trustworthy signal here, since it is far less sensitive to other load on the machine than wall clock. At `n_anchors := 1000` the same table's optimum moves out to `batch_size` 500–1000, because `n_anchors` then dominates `m` and shrinking `batch_size` only adds blocks without shrinking the block.
- **`batch_size` also trades accuracy, and which way depends on how many axes you actually use.** On the real EMP 90 bp deblur table (23,814 samples, `unweighted`, `n_dims := 10`, `n_anchors := 1000`, procrustes vs a full `unifrac_pcoa`): at **d = 3** smaller is better — M² 0.0087 / 0.0106 / 0.0188 for `batch_size` 500 / 1000 / 3000 — while at **d = 10** it inverts, 0.197 / 0.195 / 0.135. Bigger batches align the low-variance axes better; smaller batches align the leading ones better, and per-batch `batch_anchor_m2` rises with `batch_size` (0.081 / 0.117 / 0.160). Since the usual practice is *compute wide, compare narrow* (see the note under `n_dims` below), a smaller `batch_size` is usually the right call — but if you intend to interpret all ten axes, raise it. `batch_size` changes coordinates either way, so pin it alongside the seed when comparing runs.
- <a name="progressive-repro"></a>**Coordinates are reproducible only at `threads := 1`; the ordination is preserved everywhere.** Three effects, none of them a correctness bug, all of them reasons not to diff coordinates. (1) `threads` also sets the OpenMP width of the *anchor* block's reference PCoA, and skbb's centering is an OpenMP `reduction(+:)` whose summation order depends on that width — so the reference frame shifts in its last bits and every coordinate inherits the shift. (2) Concurrent *first* calls race on libssu/skbb's non-atomic CPU/accelerator detection caches, so the **first** progressive run in a session differs slightly from later ones and nothing after it does. (3) **At any width above 1, two runs of the same seeded query need not agree bit-for-bit even on the same machine** — the summation order is not pinned by the seed, only the draw is. Observed on `linux_amd64`: the same query, same `seed`, run twice in one session, differing in *every* coordinate by the usual last-bit shift. Not observed on `macos/arm64`, so treat it as platform-dependent and not as something you can test your way out of.

  **So do not expect a seeded run to be byte-reproducible unless you also pin `threads := 1`** — and note that even then the general fp32 caveat under [Reproducibility](#reproducibility) applies, since the FSVD and its BLAS reductions carry ~1e-7 noise of their own. If you need to compare two runs — different storage layout, different tree, different build — compare them as ordinations with [`procrustes`](#procrustes), whose M² is frame-invariant and therefore immune to all three effects. That is the only comparison this project is willing to promise. The same applies to `progressive_pcoa_from_distances`, which exhibits it without any UniFrac being involved at all.

  **Judge agreement on the ordination, not on coordinates.** The drift is amplified by a randomized SVD and a procrustes fit, so it is larger than one rounding step and it grows with the table: 3.1e-06 max abs on a 40-sample fixture, but 2.4e-05 on 23,814 samples — about 5e-4 relative to the largest coordinate there. What held exactly is the thing that matters: procrustes M² against a full `unifrac_pcoa` was **identical to five decimals (0.19484) at `threads := 1` and `threads := 14`**. So do not diff coordinates across thread counts or use them as a cache key; compare ordinations with [`procrustes`](#procrustes). Fixing `threads := N` and re-running in one session is **not** enough to get identical coordinates back at N > 1 — see the point above — so pin `threads := 1` if you need to compare coordinates at all. See [Reproducibility](#reproducibility) for the same fp32 caveat on the non-progressive paths.
- **Samples & features:** samples with no nonzero feature are excluded (they cannot be ordinated). The tree must cover every feature in the table (validated once at bind); each batch's features are a subset, so the check makes them all safe.
- **Cancellable:** as in [`progressive_pcoa_from_distances`](#progressive-cancellation) — cancellation is checked before every batch, so a run that would take hours stops within about one block.
- **Fewer than two samples**, `n_anchors` outside `[n_dims + 1, n_samples]`, `n_dims > n_samples - 1`, `batch_size < 1`, an unknown `variant`, or a tree missing a feature: an error. Parameter, tree and schema problems are reported when the query is bound; a failure inside a block surfaces when the run reaches it.

> **Note:** like `progressive_pcoa_from_distances`, this computes each batch's block with its own feature-table slice query, so live blocks are bounded by one `(anchors + batch)` block per worker. Rows reach the query a group of batches at a time — at least `threads` of them, plus as many more as a 64 MB coordinate budget allows — so the coordinates are never all held at once. The anchor samples' feature rows are read once and cached for the whole run, so each batch's query fetches only its own samples.

<a name="feature-table-sort-order"></a>**Performance: store the feature table sorted by `sample_id`.** Samples are batched in sorted id order, so a batch is a contiguous id range. If the table's physical layout is also sorted, DuckDB/Parquet prune by each row group's min/max `sample_id` and a batch's slice query reads only that batch's own rows; if it is not, every slice query scans the whole `sample_id` column, so the run's slicing cost grows as `n_batches × table_rows` instead of one pass.

Measured on the 17,483-sample EMP deblur table remapped to short (WoL-style) feature ids, and an 8× replication of it, `batch_size := 1000`, one batch's slice query:

| feature table | sorted by `sample_id` | unsorted |
|---|---|---|
| 17,483 samples / 13.1 M rows | 6 ms | 7 ms |
| 139,864 samples / 104.4 M rows | 9 ms | 33 ms |

Sorted stays flat as the table grows; unsorted scales with it. The gap widens with the width of `feature_id`: repeating the 13.1 M-row measurement with the original 150 bp ASV sequences as ids gave 8 ms sorted vs **140 ms** unsorted, because an unpruned scan pays for the wide string column too. For ASV-keyed tables at scale this is the difference between minutes and hours.

**You will be told.** Because nothing in the output reveals which case you are in, a multi-batch run checks the stored order up front and emits a warning naming the table when it is not sorted — readable with [`miint_warnings()`](utilities.md#miint_warnings), and on stderr. The check counts places where `sample_id` steps backwards between physically adjacent rows, which is exact (a sequence is sorted precisely when it never steps back) and costs one pass over the `sample_id` column alone: 0.09 s over 4.7 M rows, 19.6 s over 1.85 B. A single-batch run is never warned about — it reads the table once whatever the order. The warning is about cost only; ignoring it changes nothing about the coordinates.

`read_biom` emits each sample's rows contiguously, but in the BIOM file's own sample order — **grouped is not sorted**, and grouping alone prunes nothing, since each row group then spans nearly the whole id range. Add the `ORDER BY` when you materialize:

```sql
CREATE TABLE observations AS
    SELECT * FROM read_biom('table.biom') ORDER BY sample_id;
-- or, writing Parquet
COPY (SELECT * FROM read_biom('table.biom') ORDER BY sample_id)
    TO 'observations.parquet' (FORMAT PARQUET);
```

**`ORDER BY sample_id` is right for every id type**, including `BIGINT`. Samples are enumerated in the id column's own order, so batches are ranges in that same order and a table stored by `sample_id` lines up with them. The unsorted-table warning judges the column the same way, so an integer key is compared numerically rather than as text.

Sort order affects only how much a run reads — never the coordinates it produces.

A non-`VARCHAR` `sample_id` is also *compared* in its own type, so a batch's slice query is a filter DuckDB can push into the scan and prune row groups with, rather than an expression it must evaluate per row. There is nothing to configure; it follows from the column's type.

---

### Progressive PCoA (from community distances)

`progressive_pcoa_from_features(feature_table, metric, ...)` is the tree-free counterpart to [`progressive_pcoa_from_unifrac`](#progressive-pcoa-from-unifrac): the same reference-anchored progressive PCoA, but each batch's block is a [`community_distances`](#community-distances-non-phylogenetic) matrix computed on the fly over `(anchors + batch)`, so the full N×N matrix is never formed. Before this existed, a non-phylogenetic analysis could only reach [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table), which needs every one of the N²/2 pairs materialized first — the very cost progressive PCoA exists to avoid.

```sql
CREATE TABLE observations AS
    SELECT * FROM read_biom('table.biom') ORDER BY sample_id;

SELECT * FROM progressive_pcoa_from_features('observations', 'bray_curtis',
    n_dims := 10, n_anchors := 1000, batch_size := 1000, seed := 42, threads := 8);
```

**Parameters:**
- `feature_table` (VARCHAR): name of the feature relation exposing `(sample_id, feature_id, value)` (see [Feature table](#feature-table))
- `metric` (VARCHAR): one of `bray_curtis`, `euclidean`, `jaccard`, `soergel`, `morisita_horn` (case-insensitive). See below for why the other three `community_distances` metrics are not accepted.
- `n_dims` (3), `n_anchors` (100), `batch_size` (1000), `seed` (-1), `threads` (0), `global_rotation` (`true`), `anchors`: exactly as in [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table). [`pick_anchors`](#pick-anchors) applies here too.

**Output schema:** identical to the other two progressive functions — `(iteration, sample_id, axis, coordinate, eigenvalue, proportion_explained, batch, batch_anchor_m2)`, with `iteration` always `0` and anchor rows reporting `NULL` for both diagnostics (see [Reading `batch` / `batch_anchor_m2`](#progressive-batch-diagnostics)). `eigenvalue`/`proportion_explained` describe the *anchor* reference ordination, the same documented caveat as the sibling functions.

<a name="progressive-features-admissible"></a>**Only pairwise-local metrics are accepted, and the rest are refused rather than approximated.** A block carries an arbitrary subset of the samples and only the features that subset happens to use. A metric can be computed that way if and only if `d(i,j)` depends on samples `i` and `j` alone — no statistic over the other rows, and no sensitivity to features that are zero in both. `pearson`, `chisq` and `gower` all fail that test (see the block-wise column under [Community distances](#community-distances-non-phylogenetic)), and the failure is silent: each block would quietly measure a *different* distance and the ordination would come out wrong with no error raised. So they are rejected when the query is bound, with a message naming the metric and pointing at `progressive_pcoa_from_distances` over `community_distances(...)`, which forms the whole matrix and is therefore limited to sample counts that fit in memory.

That admissibility is asserted, not assumed: for every accepted metric the test suite checks that `progressive_pcoa_from_features` reproduces `progressive_pcoa_from_distances` over `community_distances(...)` *coordinate for coordinate*, at two different `batch_size` values and with both seeded and explicit anchors. On **dense** inputs read that as agreement to floating-point last bits rather than to the bit, for `euclidean` and `morisita_horn` only: a block and a whole table have different densities, so one can qualify for the [Gram kernel](#community-distances-dense-kernel) while the other does not, and the two formulations associate the same arithmetic differently. `jaccard` is exact either way (its inner products are sums of 0/1 products), and on any table sparse enough to be worth this function the question does not arise.

**Behavior:**
- **Accuracy:** as [`progressive_pcoa_from_distances`](#progressive-pcoa-from-a-distance-table) — reproduces the full ordination up to a similarity transform, with alignment error that does not compound across batches.
- **`threads` buys concurrent blocks**, as on the UniFrac path: each worker computes one block's pair loop, and the per-block width is divided out of `threads` so the two levels cannot oversubscribe.
- **Reproducibility:** seeded runs are byte-reproducible only at `threads := 1`; above that the ordination is preserved but coordinates shift in their last bits. The full explanation under [`progressive_pcoa_from_unifrac`](#progressive-repro) applies here unchanged — compare ordinations with [`procrustes`](#procrustes-align-two-ordinations), not coordinates.
- **Store the feature table sorted by `sample_id`** — batches are contiguous id ranges, and an unsorted table makes every batch scan the whole relation. Same rule, same warning, and the same measurements as [under `progressive_pcoa_from_unifrac`](#feature-table-sort-order).
- **Samples with no nonzero feature are dropped**, matching `community_distances` and `unifrac_distances`; the sparse contract cannot distinguish an all-zero sample from an absent one.
- **Cancellable** — checked before every batch, so a long run stops within about one block.
- **Errors:** an unknown or non-pairwise-local metric, fewer than two samples, `n_anchors` outside `[n_dims + 1, n_samples]`, `n_dims > n_samples - 1`, or `batch_size < 1`. Parameter and schema problems surface at bind; a failure inside a block surfaces when the run reaches it.

<a name="progressive-features-perf"></a>**Measured scaling.** On a 1,212,988-sample table (3,251,297 features, 128M nonzero cells, rarefied to even 1k), `n_dims := 10`, `n_anchors := 1000`, `batch_size := 1000`, `bray_curtis`, on a 32-core host with `memory_limit=16GB`:

| samples | 1 thread | 2 | 4 | 8 | speedup at 8 | peak RSS at 8 |
|---|---|---|---|---|---|---|
| 10,000 | 7.6 s | 4.3 s | 2.7 s | 1.9 s | 4.0× | 1.0 GB |
| 50,000 | 41.3 s | 21.6 s | 11.4 s | 6.7 s | 6.2× | 1.4 GB |
| 200,000 | 175.6 s | 90.7 s | 47.1 s | 26.6 s | 6.6× | 2.3 GB |
| **1,212,988** | **1093.9 s** | **569.0 s** | **296.9 s** | **165.1 s** | **6.6×** | **5.3 GB** |

Cost is linear in N (5.00× the samples cost 3.51× the time from 10k to 50k — sublinear while the fixed 1000-anchor cost is still being amortized — then 4.00× → 3.99× and 6.06× → 6.20×). Memory is set by `threads`, not by N or by the feature space: live memory is one block per worker. `jaccard` tracks `bray_curtis` within 2% throughout.

For scale, a dense `community_distances` matrix over that same table would be roughly 4 TB, so it is not a slow query but an impossible one.

**Blocks are sparse, which is why the feature count does not matter.** A block's pair loop merges each pair's nonzeros rather than walking the whole feature space. On this table a default block spans ~11k features but averages only ~89 nonzeros per sample, so the merge does roughly 62× less arithmetic than a dense loop would — and the block's memory is bounded by nonzeros, not by the 3.25M-feature dictionary.

**If your blocks are *not* sparse, the kernel changes under you — for the better.** Density, not feature count, is what decides: above ~1.5% of a block's cells occupied, `euclidean`, `jaccard` and `morisita_horn` switch to a [Gram-matrix kernel](#community-distances-dense-kernel) worth 2×–25× depending on how dense the block is. Real microbiome blocks sit two orders of magnitude below that line (0.19% and 0.21% on the 10k and 50k subsets here, 0.8% on the 1.2M table), so this path is not what makes those runs fast — the merge is. It matters when `progressive_pcoa_from_features` is pointed at something other than a microbiome table, where the same function would otherwise be paying microbiome-shaped costs on data that is nothing like it. The crossover was measured at **0.96%** and the switch deliberately waits until 1.5%, because right at the crossover the two kernels cost the same and only the Gram one allocates a dense block.

**Against the phylogenetic path.** Same 10,000-sample subset, same parameters, 8 threads, `unweighted` UniFrac over a 3.45M-tip tree: `progressive_pcoa_from_unifrac` took 14.5 s and 4.1 GB, `progressive_pcoa_from_features` with `bray_curtis` took 1.9 s and 1.0 GB. Choosing a tree-free metric is not a scale compromise here — it is the cheaper path.

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
- `pairing` (VARCHAR, optional): name of a `(reference_id, other_id)` relation. Absent → **full** mode (both ordinations must describe the same samples). Present → **partial** mode: fit the transform on just the paired anchor rows, then apply it to *every* row of both ordinations — a partial-procrustes technique (McDonald, qiime2/q2-diversity#338)
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
- **Partial mode (`pairing`):** the pairing must be **1:1** — a repeated `reference_id` or `other_id` is rejected (it would drop a row or feed one physical sample into the fit as several anchors). The fit uses the matched anchor rows (at least `n_dims + 1` usable pairs are required); it is then applied to every row of both inputs. `m2` is the disparity over the anchors; `pvalue` is `NULL` — a partial fit has no null model to permute against, and `q2-diversity`'s `partial_procrustes` likewise defines no Monte Carlo test.
- **Sample ids:** VARCHAR/BIGINT/UUID are all accepted (see [Sample identifier types](#sample-identifier-types)); full-mode matching and pairing lookups are by the id's string form, so a BIGINT `reference` and a VARCHAR `other` with the same ids align and emit under VARCHAR.
- **Fail loud:** a missing `(sample_id, axis, coordinate)` column, a ragged ordination (a sample missing an axis), a duplicate `(sample_id, axis)`, an out-of-range `axis`, fewer than `n_dims + 1` points/anchors, or (full mode) differing sample sets or axis counts each raise a named error.
- **P-value reproducibility:** the PROTEST test is reproducible under a fixed `seed` within one build, but it is **not** bit-for-bit comparable to `q2-diversity`'s `procrustes_analysis` — that implementation uses an *unseeded* RNG, and the C++ PRNG differs from NumPy's, so agreement is statistical (within Monte Carlo error), not exact. Disparity and coordinates, by contrast, match SciPy to machine precision.

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

**Subset selection.** Neyman, J. (1934) "On the two different aspects of the
representative method", *Journal of the Royal Statistical Society* 97(4), 558-625.
· Cochran, W.G. (1977) *Sampling Techniques*, 3rd ed., Wiley, ch. 5 (proportional
stratified allocation, the default rule in [`pick_anchors`](#pick-anchors)). ·
Gonzalez, T.F. (1985) "Clustering to minimize the maximum intercluster distance",
*Theoretical Computer Science* 38, 293-306 (greedy farthest-point / k-center).

**Ordination.** Torgerson, W.S. (1952) "Multidimensional scaling: I. Theory and
method", *Psychometrika* 17, 401-419. · Gower, J.C. (1966) "Some distance
properties of latent root and vector methods used in multivariate analysis",
*Biometrika* 53(3-4), 325-338. · Anderson, M.J. (2001) "A new method for
non-parametric multivariate analysis of variance", *Austral Ecology* 26(1), 32-46
(PERMANOVA). · Lozupone, C. and Knight, R. (2005) "UniFrac: a new phylogenetic
method for comparing microbial communities", *Applied and Environmental
Microbiology* 71(12), 8228-8235. · Faith, D.P. (1992) "Conservation evaluation and
phylogenetic diversity", *Biological Conservation* 61(1), 1-10.

**Ordination alignment.** Gower, J.C. (1975) "Generalized procrustes analysis",
*Psychometrika* 40(1), 33-51. · Jackson, D.A. (1995) "PROTEST: a PROcrustean
Randomization TEST of community environment concordance", *Écoscience* 2(3),
297-303 (the permutation test behind `pvalue`). · Peres-Neto, P.R. and Jackson,
D.A. (2001) "How well do multivariate data sets match? The advantages of a
Procrustean superimposition approach over the Mantel test", *Oecologia* 129(2),
169-178. · McDonald, D., qiime2/q2-diversity#338 (the partial-procrustes
technique — fit on a paired subset, apply to all rows — used by
[`procrustes`](#procrustes-align-two-ordinations) in partial mode and by
[`progressive_pcoa_from_unifrac`](#progressive-pcoa-from-unifrac) to place each
batch into the anchor frame).

**Reference implementations.** The `cogent3` / PyCogent, NumPy, SciPy,
scikit-learn and scikit-bio projects — consulted for metric conventions and used
as parity oracles — are credited with their licenses and citations in
[`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md). `q2-diversity` is used
the same way, as an optional test-time cross-check only (see
`test/scripts/gen_procrustes_oracle.py --check-q2`); no code from it is
incorporated. Note that `ks_2samp` is checked against `scipy.stats.ks_2samp` but
is **not** derived from it: the exact p-value comes from the lattice-path
formulation in Hodges (1958) and is computed as escaping probability mass rather
than as `1 - P(stay inside)`. See the SciPy entry there, and
`src/include/KsTwoSample.hpp`.
