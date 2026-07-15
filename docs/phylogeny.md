# Phylogeny

Methods for estimating and operating on phylogenies.

## Table of Contents

- [Resolve placements](#resolve-placements) - Fully resolve sequence placements against a backbone.
- [FastTree](#fasttree) - Estimate a phylogeny from a MSA with FastTree.

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
