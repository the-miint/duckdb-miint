# Phylogeny

Methods for estimating and operating on phylogenies.

## Table of Contents

- [Shear (subset to tips)](#shear-subset-to-tips) - Prune a tree down to a set of tips.
- [Resolve multifurcations](#resolve-multifurcations) - Resolve polytomies into a strictly bifurcating tree.
- [Resolve placements](#resolve-placements) - Fully resolve sequence placements against a backbone.
- [FastTree](#fasttree) - Estimate a phylogeny from a MSA with FastTree.
- [Independent contrasts (PIC)](#independent-contrasts-pic) - Felsenstein (1985) phylogenetic independent contrasts.

### Shear (subset to tips)

Subset ("shear" / prune) a tree down to a specified set of tips. This is the operation you want when a large reference phylogeny is persisted (e.g. in Parquet via `COPY ... (FORMAT NEWICK)`) and you need to restrict it to the tips relevant to one analysis.

**Function signature**:

`shear_tree(tree_table, tips_table, collapse := true, ignore_missing := false)`

**Parameters:**
- `tree_table` (VARCHAR): Name of a table or view containing tree data in [`read_newick`](reading.md#newick) schema (`node_index`, `parent_index`, `name`, `branch_length`, `edge_id`).
- `tips_table` (VARCHAR): Name of a table or view with a `name` (VARCHAR) column listing the tip names to keep. Passing the tips as a table (rather than a list literal) keeps the call practical when the kept set is large. Duplicate and `NULL` names are ignored.
- `collapse` (BOOLEAN, default `true`): How to treat internal nodes left with a single child after pruning.
  - `true` — standard phylogenetic shear: remove single-descendant ancestors and **sum their branch lengths onto the surviving descendant edge**, so tip-to-tip distances are preserved. The lowest common ancestor (LCA) of the kept tips becomes the new root; nodes above it are dropped.
  - `false` — preserve every internal node that lies on a path to a kept tip, with original parent links and branch lengths unchanged (unifurcations retained).
- `ignore_missing` (BOOLEAN, default `false`): When `false`, a requested tip name not present as a tip in the tree is an error (the message lists the missing names). When `true`, absent names are skipped. Either way, if *no* requested tip matches, the call errors (a tree cannot be sheared to nothing).

**Output schema:** Same as [`read_newick`](reading.md#newick) (without filepath), reindexed 0-based:
- `node_index` (BIGINT), `name` (VARCHAR), `branch_length` (DOUBLE, nullable), `edge_id` (BIGINT, nullable), `parent_index` (BIGINT, nullable, NULL for root), `is_tip` (BOOLEAN).

The schema is `UNION ALL`-compatible with `read_newick` and can be written straight back out with `COPY ... (FORMAT NEWICK)`.

**Behavior / semantics:**
- Only **tip** names are matched — an internal-node label in `tips_table` is treated as missing.
- Branch-length summation (collapse mode) treats an unspecified length (NaN → `NULL`) as `0`, but a collapsed chain that is *entirely* unspecified stays `NULL`, so topology-only trees (cladograms) round-trip as cladograms.
- A collapsed edge keeps the surviving (lower) node's `edge_id`; the `edge_id`s of merged intermediate nodes are dropped.
- The new root has no incoming edge, so its `branch_length` is `NULL`.
- Reads both `tree_table` and `tips_table` via a fresh connection (works with tables and views); single-threaded build.

**Examples:**
```sql
-- Persisted reference tree in Parquet, and the tips this analysis cares about.
CREATE TABLE ref_tree AS SELECT * FROM 'reference_tree.parquet';
CREATE TABLE keep AS SELECT feature_id AS name FROM my_feature_table;

-- Standard shear: collapse single-child ancestors, preserve distances.
SELECT * FROM shear_tree('ref_tree', 'keep');

-- Keep every internal node on the retained paths instead.
SELECT * FROM shear_tree('ref_tree', 'keep', collapse := false);

-- Shear many trees against one shared tip list, tolerating tips absent from a
-- given tree.
SELECT * FROM shear_tree('ref_tree', 'keep', ignore_missing := true);

-- Write the sheared tree back out to Newick.
COPY (
    SELECT node_index, name, branch_length, edge_id, parent_index
    FROM shear_tree('ref_tree', 'keep')
) TO 'sheared.nwk' (FORMAT NEWICK);
```

**Error conditions:**
- `tree_table` or `tips_table` does not exist.
- `tree_table` missing required `node_index` / `parent_index` columns; `tips_table` missing required `name` column.
- A requested tip is not a tip in the tree and `ignore_missing := false`.
- No requested tip matches any tip in the tree.

**Performance note (large trees):** building and shearing a large tree is
allocation-bound — profiling a multi-million-node shear shows ~40% of the time
in the system allocator (many small allocations for node names, child lists, and
index maps). glibc's allocator is comparatively slow for this pattern. Because
miint runs inside a host process, the safe way to switch allocators is a
whole-process launch option rather than anything baked into the extension:
preload jemalloc so *every* allocation in the process uses it uniformly.

```bash
# Linux: ~30% faster on large tree build/shear in our benchmarks.
LD_PRELOAD=/usr/lib/x86_64-linux-gnu/libjemalloc.so.2 duckdb   # or python, etc.
```

This applies to any allocation-heavy miint workload (large `read_newick`
builds, `tree_resolve_placement`, QC), not just `shear_tree`. It does not reduce
peak memory materially — it is a speed optimization.

### Resolve multifurcations

Resolve every multifurcation (a node with more than two children — a "polytomy") in a tree into a series of bifurcations, producing a strictly bifurcating tree. This is the operation you want before running a method that requires a binary tree — notably [`phylo_independent_contrasts`](#independent-contrasts-pic) — on a tree that contains polytomies (e.g. an unrooted tree whose root is a trifurcation, as `phylogeny_fasttree` emits, or a tree flattened from a taxonomy).

**Function signature**:

`tree_resolve_multifurcations(tree_table)`

**Parameters:**
- `tree_table` (VARCHAR): Name of a table or view containing tree data in [`read_newick`](reading.md#newick) schema (requires `node_index` and `parent_index`; `name`, `branch_length`, `edge_id` optional).

**Output schema:** Same as [`read_newick`](reading.md#newick) (without filepath), reindexed 0-based: `node_index`, `name`, `branch_length`, `edge_id`, `parent_index`, `is_tip`. The schema is `UNION ALL`-compatible with `read_newick`, can be written back out with `COPY ... (FORMAT NEWICK)`, and feeds directly into `phylo_independent_contrasts`.

**Behavior / semantics:**
- A node with `m ≥ 3` children `c0, c1, …, c(m-1)` (in `node_index` order) is resolved into a **deterministic left-comb**: the node keeps `c0` and a new internal node `N1`; `N1` keeps `c1` and `N2`; … ; `N(m-2)` keeps `c(m-2)` and `c(m-1)`. This inserts `m-2` new internal nodes.
- New connector nodes are **unnamed**, have **branch length `0`**, and no `edge_id`. Because the connectors are zero-length, the original edge lengths and all root-to-tip distances are unchanged.
- Nodes with two or fewer children are left untouched. **Single-child unifurcations are not collapsed** — that is a different operation; use [`shear_tree`](#shear-subset-to-tips) with `collapse := true` for that.
- A multifurcating root is resolved as well.
- Child order follows `node_index` (the canonical order `read_newick` assigns, encoding the original Newick left-to-right order); the reader sorts input rows by `node_index`, so the result is reproducible regardless of how the tree relation is physically stored.

**Important — the resolution is statistically arbitrary.** A polytomy carries no information about the order in which its lineages diverged, so any bifurcating resolution is as valid as any other. When the resolved tree feeds a downstream analysis such as independent contrasts, each artificially-created dichotomy contributes a contrast that is not real independent information; the effective degrees of freedom should be reduced accordingly. Resolving is a pragmatic step to satisfy a binary-tree requirement, not a reconstruction of unknown branching order.

**Examples:**
```sql
-- A tree with a trifurcating (unrooted) root, e.g. from phylogeny_fasttree.
CREATE TABLE ref_tree AS SELECT * FROM read_newick('unrooted.nwk');

-- Resolve polytomies into bifurcations (zero-length connectors).
SELECT * FROM tree_resolve_multifurcations('ref_tree');

-- Prep a tree for PIC: resolve, persist, then run independent contrasts.
CREATE TABLE binary_tree AS SELECT * FROM tree_resolve_multifurcations('ref_tree');
SELECT * FROM phylo_independent_contrasts('binary_tree', 'traits');

-- Write the resolved tree back out.
COPY (
    SELECT node_index, name, branch_length, edge_id, parent_index
    FROM tree_resolve_multifurcations('ref_tree')
) TO 'binary.nwk' (FORMAT NEWICK);
```

**Error conditions:**
- `tree_table` does not exist.
- `tree_table` missing required `node_index` / `parent_index` columns.

### Resolve placements

Resolve phylogenetic placements into a reference tree, returning a fully resolved tree with placed fragments as new tips. This exposes the `insert_fully_resolved` algorithm as a SQL-accessible table function.

**Function signature**:

`tree_resolve_placement(tree_table, placements_table)`

**Parameters:**
- `tree_table` (VARCHAR): Name of a table or view containing tree data in `read_newick` schema (requires `node_index` and `parent_index` columns; `name`, `branch_length`, `edge_id` are optional)
- `placements_table` (VARCHAR): Name of a table or view containing placement data (requires `fragment_id`, `edge_id`, `like_weight_ratio`, `distal_length`, `pendant_length` columns)

**Output schema:** Same as [`read_newick`](reading.md#newick) (without filepath):
- `node_index` (BIGINT): 0-based index of node in resolved tree
- `name` (VARCHAR): Node label (placed fragments use their fragment_id as name)
- `branch_length` (DOUBLE, nullable): Branch length
- `edge_id` (BIGINT, nullable): Edge identifier (NULL for newly created nodes)
- `parent_index` (BIGINT, nullable): Parent node's node_index (NULL for root)
- `is_tip` (BOOLEAN): Whether node is a tip/leaf

**Behavior:**
- Reads tree data from the tree table and builds a `NewickTree`
- Reads placements and inserts each fragment as a new tip on the specified edge
- Each placement creates 2 new nodes: an internal node (splitting the edge) and a fragment tip
- Deduplicates placements by `fragment_id` (keeps highest `like_weight_ratio`, then lowest `pendant_length`)
- Multiple placements on the same edge are sorted by `distal_length` and inserted as a chain
- Preserves original tip-to-tip distances in the tree
- Works with both tables and views for either parameter
- Schema is UNION ALL-compatible with `read_newick`

**Examples:**
```sql
-- Basic workflow: load tree, load placements, resolve
CREATE TABLE ref_tree AS
SELECT * FROM read_newick('reference.nwk');

CREATE TABLE placements AS
SELECT * FROM (VALUES
    ('seq1', 0::BIGINT, 0.95::DOUBLE, 0.05::DOUBLE, 0.001::DOUBLE),
    ('seq2', 1::BIGINT, 0.80::DOUBLE, 0.10::DOUBLE, 0.002::DOUBLE)
) AS t(fragment_id, edge_id, like_weight_ratio, distal_length, pendant_length);

SELECT * FROM tree_resolve_placement('ref_tree', 'placements');

-- Full jplace workflow: extract tree and placements from jplace file
CREATE TABLE jplace_tree AS
SELECT * FROM read_jplace_newick('results.jplace');

CREATE TABLE jplace_placements AS
SELECT fragment AS fragment_id, edge_num::BIGINT AS edge_id,
       like_weight_ratio, distal_length, pendant_length
FROM read_jplace('results.jplace');

-- Resolve and inspect the result
SELECT name FROM tree_resolve_placement('jplace_tree', 'jplace_placements')
WHERE is_tip = true
ORDER BY name;

-- Write resolved tree to Newick format
COPY (
    SELECT node_index, name, branch_length, edge_id, parent_index
    FROM tree_resolve_placement('ref_tree', 'placements')
) TO 'resolved.nwk' (FORMAT NEWICK);

-- Use a view for filtered placements
CREATE VIEW confident_placements AS
SELECT * FROM jplace_placements WHERE like_weight_ratio > 0.8;

SELECT COUNT(*) FROM tree_resolve_placement('jplace_tree', 'confident_placements')
WHERE is_tip = true;
```

**Error conditions:**
- Tree table or placements table does not exist
- Tree table missing required `node_index` or `parent_index` columns
- Placements table missing required columns (`fragment_id`, `edge_id`, `like_weight_ratio`, `distal_length`, `pendant_length`)
- Placement references an `edge_id` not present in the tree
- `distal_length` exceeds the edge's branch length
- Negative `distal_length` or `pendant_length`

### FastTree

Build an approximately-maximum-likelihood phylogenetic tree from a multiple sequence alignment using [FastTree 2](http://www.microbesonline.org/fasttree/) (Price, Dehal, and Arkin 2010, Mol. Biol. Evol. 27:1641-1650). Returns one row per node (tips and internal nodes) with parent links suitable for joining back to the input or rendering to Newick via `COPY ... (FORMAT 'newick')`.

**Function signature**:
`phylogeny_fasttree(table_name, [options])`

**Architecture:** FastTree is statically linked into a separate process (`gpl-boundary`) that miint launches per query. The daemon is GPL-licensed; miint stays BSD-clean by communicating with it over POSIX shared memory and Arrow IPC. See [`docs/internals/embedded-tools.md`](internals/embedded-tools.md) for the protocol details and process-lifecycle invariants.

**Requirements:**
- The `gpl-boundary` binary must be on `PATH`. Use [`phylogeny_fasttree_available()`](utilities.md#phylogeny_fasttree_available) to check at runtime.

**Parameters:**

- `table_name` (VARCHAR): Name of a table or view containing the input alignment. Must have columns `name` (VARCHAR) and `sequence` (VARCHAR), one row per tip. Minimum 3 sequences (FastTree cannot build a meaningful tree below that). All sequences should be the same length (a true MSA — pad with `-` if needed; see `align_mafft`).

The remaining parameters mirror the FastTree CLI flags. Defaults are the daemon's FastTree defaults; omitted parameters are not emitted into the wire config.

*Sequence type and model:*

- `seq_type` (VARCHAR, default `'auto'`): One of `'auto'`, `'nucleotide'`, `'protein'`. Drives which substitution models are valid.
- `model` (VARCHAR, default `'auto'`): Substitution model. One of `'auto'`, `'jtt'`, `'lg'`, `'wag'` (protein), `'jc'`, `'gtr'` (nucleotide).
- `gtrrates` (DOUBLE[], optional): GTR rate parameters (6 entries: AC, AG, AT, CG, CT, GT, all ≥ 0). Only valid with `model='gtr'`.
- `gtrfreq` (DOUBLE[], optional): GTR base frequencies (4 entries: A, C, G, T, all ≥ 0). Only valid with `model='gtr'`.
- `gamma` (BOOLEAN, optional): Use a discretized gamma distribution for among-site rate variation.
- `cat` (BIGINT, optional, ≥ 1): Number of rate categories for CAT approximation.

*Tree search:*

- `nni` (BIGINT, optional, ≥ 0): Number of NNI rounds. `0` disables NNI.
- `spr` (BIGINT, optional, ≥ 0): Number of SPR rounds (subtree prune-regraft). `0` disables SPR.
- `mlnni` (BIGINT, optional, ≥ 0): Number of ML NNI rounds.
- `mlacc` (BIGINT, optional, ≥ 1): ML accuracy iterations per branch.
- `noml` (BOOLEAN, optional): Skip the ML refinement phase entirely. Mutually exclusive with `mlnni > 0`.
- `slow` (BOOLEAN, optional): Use slower / more accurate neighbor-joining (disables top hits).
- `bionj` (BOOLEAN, optional): Use BIONJ joins. Mutually exclusive with `nj`.
- `nj` (BOOLEAN, optional): Use vanilla NJ joins. Mutually exclusive with `bionj`.

*Top-hits heuristics:*

- `top` (BOOLEAN, optional): Force top-hits heuristics on.
- `notop` (BOOLEAN, optional): Force top-hits heuristics off.
- `topm` (DOUBLE, optional, > 0): Top-hits multiplier (typical: `1.0`-`2.0`). Only meaningful when top hits are enabled (rejected if `notop=true` or `slow=true`).

*Bootstrap / support:*

- `bootstrap` (BIGINT, optional, ≥ 0): Number of resamples for SH-test branch support. `0` disables. Mutually exclusive with `nosupport=true`.
- `nosupport` (BOOLEAN, optional): Suppress branch support output. Mutually exclusive with `bootstrap > 0`.

*Pseudocounts (low-confidence trees):*

- `pseudo` (BOOLEAN, optional): Enable pseudocounts.
- `pseudo_weight` (DOUBLE, optional, ≥ 0): Pseudocount weight. Requires `pseudo=true`.

*Misc:*

- `seed` (BIGINT, optional): RNG seed. Pin for reproducibility (combine with `threads=1`).
- `threads` (BIGINT, optional, ≥ 1): Thread count. **Output is non-deterministic across thread counts** — FastTree's parallel NJ/NNI/SPR sections produce floating-point variation. Force `threads=1` for reproducibility.
- `verbose` (BOOLEAN, optional): Forward FastTree's diagnostic output to the daemon log.
- `quote` (BOOLEAN, optional): Quote tip names in the rendered Newick (when going through `COPY FORMAT 'newick'`).

**Cross-parameter validation rules** (raised at SQL bind time, before any process spawn):

| Constraint | Error |
|---|---|
| `bootstrap > 0` + `nosupport=true` | "bootstrap=N cannot be combined with nosupport=true" |
| `pseudo_weight` without `pseudo=true` | "pseudo_weight requires pseudo=true" |
| `mlnni > 0` + `noml=true` | "mlnni=N cannot be combined with noml=true" |
| Protein model + `seq_type='nucleotide'` | "model='X' is a protein model and cannot be combined with seq_type='nucleotide'" |
| Nucleotide model + `seq_type='protein'` | "model='X' is a nucleotide model and cannot be combined with seq_type='protein'" |
| `gtrrates` / `gtrfreq` without `model='gtr'` | "gtrrates/gtrfreq are only valid with model='gtr'" |
| `bionj=true` + `nj=true` | "bionj=true and nj=true are mutually exclusive" |
| `top=true` + `notop=true` | "top=true and notop=true are mutually exclusive" |
| `slow=true` + `top=true` | "slow=true disables top hits; cannot also set top=true" |
| `topm` + (`notop=true` or `slow=true`) | "topm is only meaningful when top hits are enabled" |
| `fastest=true` | Currently rejected — the embedded library API is a strict subset of the CLI's `-fastest`. Will be supported once the upstream surface catches up. |

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `node_index` | BIGINT | 0-based row index, also the node identifier within this tree |
| `parent_index` | BIGINT (nullable) | `node_index` of this node's parent; `NULL` for the root |
| `edge_id` | BIGINT (nullable) | Edge identifier (matches the `{N}` placement notation, when used with `read_jplace`) |
| `branch_length` | DOUBLE (nullable) | Length of the edge from this node to its parent (NULL on the root) |
| `support` | DOUBLE (nullable) | SH-like branch support in `[0,1]` when `bootstrap > 0`; NULL otherwise |
| `is_tip` | BOOLEAN | True if this is a leaf (tip) node |
| `name` | VARCHAR (nullable) | Tip name (matches the input `name` column) for tips; NULL or empty for internal nodes |
| `n_children` | BIGINT | Number of direct children (0 for tips, 2 for binary internals, 3 for the unrooted root) |

The schema is `UNION ALL`-compatible with `read_newick` if you select `node_index, parent_index, name, branch_length, edge_id, is_tip` (i.e., drop `support` and `n_children`).

**Behavior:**

- One `gpl-boundary` daemon per query call (no cross-query connection caching yet).
- Sequences are read from `table_name` via a fresh DuckDB connection (avoids re-entrant `context.Query()` deadlock — see [`docs/internals/reading-tables-views.md`](internals/reading-tables-views.md)).
- The full alignment is materialized in shared memory and submitted as a single Arrow IPC batch (FastTree is not a streaming algorithm).
- The daemon is sent SIGTERM in the function's destructor; SIGKILL after a 30 × 10ms grace period if it doesn't exit cleanly.
- FastTree's text output uses `%.9f` for branch lengths; miint's intermediate Arrow representation is full IEEE-754 double precision. Round-trip parity to the CLI binary is verified in `test/sql/phylogeny_fasttree_parity.test`.

**Examples:**

```sql
-- Basic workflow: align with MAFFT, then build a tree
CREATE TABLE seqs AS SELECT read_id AS name, sequence1 AS sequence FROM read_fastx('16s.fasta');

-- Default tree (auto-detected seq_type, FastTree defaults for everything else)
SELECT * FROM phylogeny_fasttree('seqs');

-- Reproducible run: pin seed and force threads=1
SELECT * FROM phylogeny_fasttree('seqs', seed := 42, threads := 1);

-- Nucleotide GTR with branch support from 100 bootstraps
SELECT * FROM phylogeny_fasttree(
    'seqs',
    seq_type := 'nucleotide',
    model := 'gtr',
    bootstrap := 100,
    seed := 42
);

-- GTR with explicit rate parameters
SELECT * FROM phylogeny_fasttree(
    'seqs',
    seq_type := 'nucleotide',
    model := 'gtr',
    gtrrates := [1.0, 2.5, 1.0, 1.0, 2.5, 1.0],
    gtrfreq  := [0.25, 0.25, 0.25, 0.25]
);

-- Render the tree to a Newick file
COPY (
    SELECT node_index, parent_index, name, branch_length, edge_id
    FROM phylogeny_fasttree('seqs', seed := 42, threads := 1)
) TO 'tree.nwk' (FORMAT 'newick');

-- Tip-only summary
SELECT name, branch_length AS terminal_branch_length
FROM phylogeny_fasttree('seqs', seed := 42)
WHERE is_tip
ORDER BY branch_length DESC;

-- Pipeline: align → tree → resolve placements (round-trip via read_newick)
COPY (
    SELECT node_index, parent_index, name, branch_length, edge_id
    FROM phylogeny_fasttree('seqs', seed := 42, threads := 1)
) TO 'ref.nwk' (FORMAT 'newick');

CREATE TABLE ref_tree AS SELECT * FROM read_newick('ref.nwk');
SELECT * FROM tree_resolve_placement('ref_tree', 'placements');
```

**Error conditions:**

- `gpl-boundary` not on PATH (use `phylogeny_fasttree_available()` to gate)
- Input table does not exist or is missing `name` / `sequence` columns
- Fewer than 3 input sequences
- Any of the cross-parameter rules above
- Daemon process exits non-zero or returns a non-OK protocol response (exit message includes the daemon's stderr)

**Reproducibility:** With `threads=1` and a fixed `seed`, the tree is bit-deterministic across runs of the same `gpl-boundary` build. With `threads>1`, FastTree's parallel sections produce floating-point variation in branch lengths and (rarely) topology — pin `threads=1` for deterministic output.

### Independent contrasts (PIC)

Compute Felsenstein's (1985) **phylogenetic independent contrasts** for one or more numeric per-tip traits over a tree. Per-tip traits are not statistically independent — related tips share evolutionary history — so naively correlating two traits across tips is biased. Independent contrasts is the standard correction: it transforms the correlated tip values into `n-1` contrasts that are independent and identically distributed under a Brownian-motion model, which can then be correlated or regressed (through the origin) in downstream SQL.

Reference: Felsenstein, J. (1985) "Phylogenies and the Comparative Method." *The American Naturalist* 125(1):1–15.

**Function signature**:

`phylo_independent_contrasts(tree_table, traits_table)`

**Parameters:**
- `tree_table` (VARCHAR): Name of a table or view containing tree data in [`read_newick`](reading.md#newick) schema (requires `node_index` and `parent_index`; `name` and `branch_length` are needed for a meaningful contrast; `edge_id` optional). The tree must be **strictly bifurcating** with **finite, non-negative branch lengths** (see requirements below).
- `traits_table` (VARCHAR): Name of a table or view in **long form** with columns `name`, `trait`, `value`:
  - `name`: the tip identifier, matched to the tree's tip labels. Type-flexible — `VARCHAR`, `UUID`, and `BIGINT` are all accepted and matched to tip labels by their canonical text form (the same posture as `read_id` elsewhere), so a numeric or UUID feature key works directly.
  - `trait` (VARCHAR): the trait name; one trait is computed per distinct value.
  - `value` (DOUBLE): the numeric trait value for that tip.

  Every trait must cover **exactly** the tree's tip set — no missing tips, no extra names, no `NULL` values. To analyze a subset of tips, [`shear_tree`](#shear-subset-to-tips) the tree to those tips first: pruning re-sums branch lengths, which changes the contrasts, so it cannot be skipped by simply omitting tips from the traits.

**Output schema:** long form, one row per (internal node, trait):

| Column | Type | Description |
|--------|------|-------------|
| `node_index` | BIGINT | The input tree's `node_index` for the internal node at which the contrast is computed |
| `trait` | VARCHAR | The trait name |
| `contrast` | DOUBLE | The standardized independent contrast `(X_i − X_j) / sqrt(v_i + v_j)` |
| `ancestral_estimate` | DOUBLE | The reconstructed (inverse-variance-weighted) trait value at the node |
| `contrast_variance` | DOUBLE | `v_i + v_j`, the contrast's expected variance (a property of the tree, identical across traits) |

For a rooted bifurcating tree with `n` tips there are `n-1` internal nodes and thus `n-1` contrasts per trait.

**Behavior / semantics:**
- Computed by Felsenstein's post-order pruning. At each internal node with children `i`, `j` and variance-extended branch lengths `v_i`, `v_j`: `contrast = (X_i − X_j)/sqrt(v_i+v_j)`; the ancestral estimate is the inverse-variance-weighted mean `(X_i·v_j + X_j·v_i)/(v_i+v_j)`; and the node's branch length is extended by `(v_i·v_j)/(v_i+v_j)` for use higher in the tree.
- **The sign of a contrast is arbitrary** — it depends on which child is `i` vs `j` (here, `node_index` order). This is why the standard two-trait analysis is a regression *through the origin* over the paired contrasts (or a correlation after fixing one axis's sign). Results are reproducible: child order follows `node_index`.
- `contrast_variance` is identical across traits for a given node, which is what makes multi-trait computation cheap (the tree is validated and its variances precomputed once).
- The two-trait correlation / regression is intentionally **not** built in — layer it in SQL over the per-node output (example below).

**Requirements and error conditions:**
- **Strictly bifurcating.** Every internal node must have exactly two children. A polytomy (`> 2` children) errors and points to [`tree_resolve_multifurcations`](#resolve-multifurcations); a unifurcation (`1` child) errors and points to [`shear_tree`](#shear-subset-to-tips). Common gotcha: an *unrooted* Newick tree has a trifurcating root and must be resolved (or rooted) first.
- **Finite, non-negative branch lengths on every non-root edge.** PIC divides by `v_i + v_j`, so an unspecified (`NULL`/NaN) or infinite length errors, and a negative length errors. The root's own branch length is unused. Zero-length internal edges are allowed (e.g. the connectors from `tree_resolve_multifurcations`); only a contrast whose variance `v_i + v_j` is exactly zero (both children zero-length) errors.
- **Trait coverage.** A tip with no trait value, a trait value for a name that is not a tip, a duplicate `(name, trait)` pair, or a `NULL` `value` all error.
- Tip names must be unique (they key the traits); `tree_table` / `traits_table` must exist and have the required columns.

**Examples:**
```sql
-- Tree (must be bifurcating) and per-tip traits in long form.
CREATE TABLE tree AS SELECT * FROM read_newick('tree.nwk');
CREATE TABLE traits AS SELECT * FROM (VALUES
    ('sp1', 'body_mass', 3.2), ('sp2', 'body_mass', 1.1), ('sp3', 'body_mass', 2.7),
    ('sp1', 'metabolic_rate', 0.9), ('sp2', 'metabolic_rate', 0.4), ('sp3', 'metabolic_rate', 0.7)
) AS v(name, trait, value);

-- Contrasts for every trait.
SELECT * FROM phylo_independent_contrasts('tree', 'traits') ORDER BY trait, node_index;

-- Two-trait association via regression THROUGH THE ORIGIN over paired contrasts
-- (through the origin because each contrast's sign is arbitrary).
WITH c AS (SELECT * FROM phylo_independent_contrasts('tree', 'traits'))
SELECT sum(x.contrast * y.contrast) / sum(x.contrast * x.contrast) AS pic_slope
FROM c AS x JOIN c AS y USING (node_index)
WHERE x.trait = 'body_mass' AND y.trait = 'metabolic_rate';

-- If the tree has polytomies, resolve first (see the arbitrariness caveat there).
CREATE TABLE binary_tree AS SELECT * FROM tree_resolve_multifurcations('tree');
SELECT * FROM phylo_independent_contrasts('binary_tree', 'traits');
```

**Cleanroom note:** this is a from-scratch implementation of the algorithm as published in Felsenstein (1985); it does not derive from any GPL comparative-methods package.
