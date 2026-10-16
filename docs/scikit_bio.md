# scikit-bio methods

Some statistics in MIINT are computed by [scikit-bio](https://scikit.bio) itself, so their answers are scikit-bio's, not a reimplementation's. They run in [skb-arrow](https://github.com/the-miint/skb-arrow), a separate program that MIINT starts for the query and stops when it ends. Your data reaches it as Arrow files and its results come back as rows; there is no Python inside DuckDB.

## Table of Contents

- [Installing](#installing) - `install_skb_arrow` / `skb_arrow_available`
- [How a call runs](#how-a-call-runs) - processes, temporary files, threads, interrupts, errors
- [Mantel test](#skb_mantel) - `skb_mantel`: the correlation of two distance matrices and its permutation p-value
- [ANCOM-BC](#skb_ancombc) - `skb_ancombc`: differential abundance with bias correction, and its post-hoc tests

## Installing

skb-arrow is installed with [uv](https://docs.astral.sh/uv/), which brings its own Python, so neither your Python nor conda's is touched. With uv installed (`curl -LsSf https://astral.sh/uv/install.sh | sh`):

```sql
SELECT (install_skb_arrow()).message;
SELECT skb_arrow_available();   -- true
```

### `install_skb_arrow([force])`

Find a working skb-arrow, or install the release MIINT is tested with (skb-arrow 0.1.0, on Python 3.14) into MIINT's cache (`$XDG_CACHE_HOME/miint` or `~/.cache/miint`). It runs `uv tool install --managed-python --python 3.14 skb-arrow==0.1.0` with the program, its environment and its interpreter all inside that cache, about 850 MB on Linux x86_64.

**Parameters:**
- `force` (BOOLEAN, optional, default `false`): skip the check for a working skb-arrow, and have uv reinstall it into the cache and replace the program there.

**Returns:** `STRUCT(installed BOOLEAN, path VARCHAR, version VARCHAR, message VARCHAR)`

| Field | Meaning |
|---|---|
| `installed` | `true` when the skb-arrow at `path` works: it answers MIINT's handshake with every capability MIINT uses. |
| `path` | The program. Empty when there is none. |
| `version` | skb-arrow's version, e.g. `0.1.0`. |
| `message` | What happened: "already available", "Installed skb-arrow 0.1.0 to ...", or why it failed. A failure is reported here, not raised. |

**Behavior:**
- **Idempotent.** A working skb-arrow (found as below) returns `installed = true` at once, offline. A broken or incompatible one is installed over.
- **Lookup order**, used by every function on this page: `$MIINT_SKB_ARROW_PATH`, then MIINT's cache, then `PATH`. An installed copy in the cache therefore shadows one on `PATH`; `$MIINT_SKB_ARROW_PATH` overrides both, even after an install.
- **Needs uv, and the network.** Without uv, `message` gives uv's install command. uv is looked for on `PATH`, then at `~/.local/bin/uv`, where its installer puts it.
- **Platforms:** skb-arrow supports Linux x86_64 and macOS arm64. On WASM and Windows builds, MIINT has no skb-arrow support: `install_skb_arrow()` returns `installed = false` with that reason, `skb_arrow_available()` returns `false`, and the functions below do not exist.
- **Cannot be cancelled.** Ctrl-C waits for uv to finish.
- **Provenance.** uv installs from PyPI, at the pinned release whose own dependencies are pinned exactly; MIINT does not verify file hashes. Where you need attested provenance, install skb-arrow your own way and set `MIINT_SKB_ARROW_PATH` to it.

### `skb_arrow_available()`

`true` when MIINT finds an skb-arrow that answers its handshake. A `true` is remembered for the life of the process; a `false` is checked again on the next call, so running `install_skb_arrow()` in between is seen.

## How a call runs

- **One skb-arrow per query.** A query calling `skb_mantel` or `skb_ancombc` starts its own skb-arrow and stops it when the query ends. Starting takes 0.5–1 s (Python and scikit-bio importing); `skb_mantel` then compiles its numba kernel, about 0.3 s more. `EXPLAIN` starts nothing.
- **Data travels as temporary files** in a private directory, created for the query and removed after it:
  - On Linux, the directory goes in `/dev/shm`, which is memory, so nothing touches disk. It needs room for the inputs and outputs.
  - On macOS, it goes in `$TMPDIR`.
  - To put it elsewhere, set `MIINT_SKB_ARROW_SEGMENT_DIR` to a directory with room.
  - **In Docker**, `/dev/shm` defaults to 64 MiB. Run the container with `--shm-size=1g` (or more), or set `MIINT_SKB_ARROW_SEGMENT_DIR`. When space runs out, the error says so and names the variable.
- **Threads.** `threads := n` caps the threads skb-arrow uses (numba's, OpenMP's and BLAS's). The default, 0, follows DuckDB's `threads` setting. It changes speed, not answers: `skb_mantel`'s don't depend on it at all, and `skb_ancombc`'s only to floating-point rounding (see [Reproducibility](#reproducibility)).
- **Interrupts.** Ctrl-C, or a `max_execution_time` deadline, stops the query, skb-arrow, and everything skb-arrow started, and removes the temporary directory.
- **Checks and errors.** skb-arrow checks every input and parameter before it computes, and MIINT reports its refusals as `Invalid Input Error: <function>: <skb-arrow's message>`, which names the problem, e.g. `x: 1 of 6 pairs are missing, e.g. ('c', 'd')`. A message beginning `skb-arrow internal error (please report)` is a bug; it carries skb-arrow's traceback. A failure never leaves the database unusable.
- **Warnings** from scikit-bio and its dependencies appear in `miint_warnings()`, as `<function>: <category>: <message>`.
- **Parameters** you leave unset take skb-arrow's defaults, listed below.

### Reproducibility

Every stochastic method takes an integer `seed` with a fixed default, 0, so a repeated query gives the same answer. For a given skb-arrow version, answers are bit-for-bit reproducible on the same machine and thread count, and agree to floating-point tolerance across machines: BLAS (which `skb_ancombc` uses) can move the last bits.

## `skb_mantel`

`skb_mantel(x, y, ...)` is the Mantel test: the correlation between two distance matrices over the same samples, and a permutation p-value for it. It runs scikit-bio's `mantel`.

```sql
CREATE TABLE observations AS SELECT * FROM read_biom('data/biom/test.biom');
CREATE TABLE dm_bc AS SELECT * FROM community_distances('observations', 'bray_curtis');
CREATE TABLE dm_jac AS SELECT * FROM community_distances('observations', 'jaccard');

SELECT * FROM skb_mantel('dm_bc', 'dm_jac');
-- statistic 0.7742362028809426, pvalue 0.002, n 6

SELECT * FROM skb_mantel('dm_bc', 'dm_jac', method := 'spearman', permutations := 9999, seed := 42);
```

**Parameters:**
- `x`, `y` (VARCHAR): names of two condensed distance relations, each exposing `(sample_a, sample_b, distance)`. The output of [`unifrac_distances`](diversity.md#unifrac-distances) and [`community_distances`](diversity.md#community-distances-non-phylogenetic) is one. `x` is the one permuted.
- `method` (VARCHAR, default `'pearson'`): `'pearson'`, `'spearman'`, or `'kendalltau'`.
- `permutations` (BIGINT, default 999): permutations for the p-value; 0 skips it.
- `alternative` (VARCHAR, default `'two-sided'`): `'two-sided'`, `'greater'` (a positive correlation), or `'less'` (a negative one).
- `seed` (BIGINT, default 0).
- `threads` (INTEGER, default 0): see [How a call runs](#how-a-call-runs).

**Output:** one row, `(statistic DOUBLE, pvalue DOUBLE, n BIGINT)`. `n` is the number of samples.

**Inputs:**
- **Condensed and complete.** Each unordered pair of distinct samples appears exactly once, in either orientation, and no row pairs a sample with itself. A square matrix (both orientations and the diagonal) is refused; pass its upper triangle, e.g. a view `SELECT * FROM dm WHERE sample_a < sample_b`. A missing pair is an error naming it.
- **Every distance present.** A NULL, NaN, infinite or negative distance is an error. (This is stricter than [`pcoa`](diversity.md#pcoa-from-a-distance-table), which skips NULLs.)
- **The same samples in both.** `x` and `y` must hold the same samples, as the same kind of ID: two integer-ID matrices, or two text-ID ones. Within each, `sample_a` and `sample_b` must be one kind too; an INTEGER `sample_a` beside a VARCHAR `sample_b` is refused. Integer IDs of up to 64 bits go to skb-arrow as integers, and every other type as text. A UBIGINT ID beyond BIGINT's range is an error; cast such IDs to VARCHAR.
- Distances of any numeric type are read as DOUBLE.

**Behavior:**
- **Constant input.** When either matrix is constant, the correlation is undefined: `statistic` and `pvalue` are NaN. For `pearson` and `spearman`, scipy's `ConstantInputWarning` appears in `miint_warnings()`.
- **Memory.** scikit-bio draws every permutation up front: (`permutations` + 1) × samples × 8 bytes. That is 16 MB for 999 permutations of 2,000 samples, and 16 GB for a million.
- **`kendalltau` is slow.** Its permutations run in Python, not numba: about a minute at 1,000 samples and 999 permutations.
- **It agrees with SQL.** With `method := 'pearson'`, `statistic` is the Pearson correlation of the paired distances, `corr(x.distance, y.distance)` over the two relations joined on `(sample_a, sample_b)`. See [Correlating two distance matrices](diversity.md#correlating-two-distance-matrices-mantel). Only the p-value needs `skb_mantel`.

## `skb_ancombc`

`skb_ancombc(table, metadata, formula, ...)` is ANCOM-BC (Lin & Peddada 2020): differential abundance testing on compositional counts, with a correction for each sample's sampling fraction. It runs scikit-bio's `ancombc`, with its global, pairwise and Dunnett post-hoc tests.

```sql
CREATE TABLE observations AS SELECT * FROM read_biom('data/biom/test.biom');
CREATE TABLE metadata AS SELECT * FROM read_csv('data/unifrac/small_metadata.csv');

-- treatment B against A, per feature. The table has zeros, so it needs a pseudocount.
SELECT feature_id, term, lfc, qvalue, signif
FROM skb_ancombc('observations', 'metadata', 'treatment', pseudocount := 1)
WHERE term = 'treatment[T.B]';

-- With a covariate of three or more groups (here a `region` column of a
-- metadata relation `visits`), the post-hoc tests compare them:
SELECT * FROM skb_ancombc('counts', 'visits', 'region + age',
    grouping := 'region', posthoc := ['global', 'pairwise'], pseudocount := 1);
```

**Parameters:**
- `table` (VARCHAR): name of a long feature table, `(sample_id, feature_id, value)`. See [Inputs](#ancom-bc-inputs).
- `metadata` (VARCHAR): name of a wide metadata relation: a `sample_id` column and a column per covariate.
- `formula` (VARCHAR): the model's right-hand side, in patsy's syntax, naming metadata columns: `'treatment'`, `'region + age'`, `'region * age'`. See [Formulas](#formulas).
- `grouping` (VARCHAR, default none): the covariate the post-hoc tests compare. It must be a text column, a term of its own in `formula`, and have at least three groups among the table's samples; `formula` must have an intercept.
- `posthoc` (VARCHAR[], default `[]`): any of `'global'`, `'pairwise'`, `'dunnett'`. Needs `grouping`.
- `pseudocount` (DOUBLE, default 0): added to every cell. A table with zeros needs one (R's ANCOMBC tutorial uses 1).
- `max_iter` (BIGINT, default 100), `tol` (DOUBLE, default 1e-5): the bias-correction EM's iteration limit and convergence tolerance.
- `alpha` (DOUBLE, default 0.05): the significance level: `signif` is `qvalue <= alpha`.
- `p_adjust` (VARCHAR, default `'holm'`): `'holm'`, `'bonferroni'`, `'bh'` or `'by'`.
- `bootstraps` (BIGINT, default 100), `seed` (BIGINT, default 0): for `'dunnett'` only, the one stochastic test.
- `threads` (INTEGER, default 0): see [How a call runs](#how-a-call-runs).

A DOUBLE parameter must be finite: NaN or infinity is an error.

**Output:** one long table, in skb-arrow's order: by `test`, then feature, then term.

| Column | Type | Meaning |
|---|---|---|
| `feature_id` | as the input's (below) | the feature |
| `test` | VARCHAR | `'main'`, then each requested post-hoc test in the order `'global'`, `'pairwise'`, `'dunnett'` |
| `term` | VARCHAR | the design column or comparison: `'Intercept'`, `'region[T.NE]'`, `'age'`, `'region[T.c]_region[T.b]'` (c minus b); NULL for `'global'` |
| `lfc`, `se` | DOUBLE | the natural-log fold change and its standard error; NULL for `'global'` |
| `w` | DOUBLE | the test statistic: a z-statistic for `'main'`, chi-square (or F) for `'global'` |
| `pvalue`, `qvalue` | DOUBLE | the p-value, and as adjusted by `p_adjust` |
| `signif` | BOOLEAN | `qvalue <= alpha` |

`'pairwise'` compares each group with the reference, then each pair of the others; `'dunnett'` compares each with the reference. Both correct only across the features a screen keeps, and give a dropped feature `pvalue` 1. NaN and ±infinity are scikit-bio's, passed through. Nothing flags a coefficient the data can't estimate: a collinear design gives finite numbers.

`feature_id` keeps the input's type when it is BIGINT or UUID, and is VARCHAR for anything else (INTEGER included), as with [other MIINT functions](diversity.md#sample-identifier-types).

### ANCOM-BC inputs

- **Feature table.** Each `(sample_id, feature_id)` at most once, and no NULL or negative `value`; an absent pair counts as 0, before the `pseudocount` is added. Unlike the UniFrac functions, a duplicate row is an error, not summed: aggregate first (`GROUP BY sample_id, feature_id`). `value` of any numeric type is read as DOUBLE.
- **Metadata.** Every sample in the table needs a row, and at most one; other rows are ignored. Its `sample_id` must be the same kind of ID as the table's: both integers or both text. A UBIGINT ID beyond BIGINT's range is an error; cast such IDs to VARCHAR. No covariate the formula uses may be NULL, NaN or infinite for a sample in the table.
- **Only the columns the formula or `grouping` names are read**, so other columns, of any type, don't matter.
- **Covariate types:**
  - BOOLEAN and text (VARCHAR) are categorical, with levels sorted; the first is the reference level.
  - Integers, FLOAT, DOUBLE, DECIMAL and HUGEINT are numeric.
  - An **ENUM is refused**. Its level order would set the reference level, but MIINT can't yet send one. Cast it to VARCHAR in a view (`CREATE VIEW m AS SELECT * REPLACE (region::VARCHAR AS region) FROM metadata`), and note that the reference level then becomes the alphabetically first level, not the ENUM's first. To choose another, rename the levels so the one you want sorts first.
  - Other types (DATE, lists, ...) are refused.

### Formulas

- A formula is checked before anything evaluates it. Every factor must be a metadata column's bare name, and names are case-sensitive.
- The operators `+ - * : / **` are allowed, and so are `0` and `1` for the intercept.
- Function calls are refused: `log(age)`, `C(region)`, `Q('a b')`. Transform or rename columns in a view instead.
- A column named `Intercept` collides with the one patsy adds; rename it.

## Citations

- Lin, H., & Peddada, S. D. (2020). Analysis of compositions of microbiomes with bias correction. *Nature Communications*, 11, 3514. https://doi.org/10.1038/s41467-020-17041-7
- Mantel, N. (1967). The detection of disease clustering and a generalized regression approach. *Cancer Research*, 27(2), 209–220.
