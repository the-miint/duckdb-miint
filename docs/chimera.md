# Chimera checking

Detect chimeric sequences — PCR artifacts formed by joining fragments of two different parent sequences — so you can flag or remove them before downstream analysis. Use reference-based detection when you have a trusted chimera-free database to compare against, or de novo detection when you only have your own sequences and their abundances.

These methods are powered by the embedded [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584), which implements the UCHIME algorithm (Edgar et al. 2011, Bioinformatics 27:2194-2200).

## Table of Contents

- [Reference-based (UCHIME)](#reference-based-uchime) - Detect chimeras by comparing queries against a trusted reference database.
- [De novo (UCHIME)](#de-novo-uchime) - Detect chimeras without a reference, using abundance information.

### Reference-based (UCHIME)

Reference-based chimera detection using the UCHIME algorithm, detecting chimeric sequences by comparing queries against a trusted chimera-free reference database.

**Function signature**:

`detect_chimera_uchime(query_table, db='refs_table', [sample_id='col'], [options])`

**Parameters:**
- `query_table` (VARCHAR): Name of a table or view containing query sequences. Must have `read_id` (VARCHAR, BIGINT, or UUID — see *Identifier-column types* below) and `sequence1` (VARCHAR) columns.
- `db` (VARCHAR, required): Name of a table or view containing reference sequences. Same schema requirements as `query_table` (its `read_id` type is independent of the query's).
- `sample_id` (VARCHAR, optional): Name of a column in `query_table` to partition by. When provided, queries are scored per-sample against the (shared, load-once) reference database, and the sample column is prepended to the output. Execution is serialized (the vsearch wrapper is not thread-safe across concurrent calls). A TEMP table or view cannot be used as the source in this mode: the per-sample path builds a fixed-name temporary view on a private per-thread connection, which must stay isolated so parallel workers do not collide on that name. Materialize the source as an ordinary table first. Tracked as #207.
- `minh` (DOUBLE, default 0.28): Minimum h-score to flag as chimeric. Range [0, 1].
- `xn` (DOUBLE, default 8.0): Weight of "no" votes in h-score computation. Must be >= 1.0.
- `dn` (DOUBLE, default 1.4): Pseudo-count prior on "no" votes. Must be >= 0.
- `mindiv` (DOUBLE, default 0.8): Minimum divergence (percentage points) from closest parent. Must be >= 0.
- `mindiffs` (INTEGER, default 3): Minimum number of diffs in each segment (left and right). Must be >= 1.
- `threads` (INTEGER, optional): Number of threads vsearch uses for the internal `chimera_detect_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema (18 columns, compatible with vsearch `--uchimeout`):**

| Column | Type | Description |
|--------|------|-------------|
| `score` | DOUBLE | Chimera h-score (higher = more likely chimeric) |
| `read_id` | VARCHAR, BIGINT, or UUID | Query sequence identifier (mirrors `query_table.read_id`) |
| `parent_a_id` | VARCHAR, BIGINT, or UUID | Parent A identifier (mirrors `db.read_id`; NULL if non-chimeric) |
| `parent_b_id` | VARCHAR, BIGINT, or UUID | Parent B identifier (mirrors `db.read_id`; NULL if non-chimeric) |
| `closest_parent_id` | VARCHAR, BIGINT, or UUID | Closest parent to query (mirrors `db.read_id`; NULL if non-chimeric) |
| `id_query_model` | DOUBLE | Query-to-chimeric-model identity % |
| `id_query_a` | DOUBLE | Query-to-parent-A identity % |
| `id_query_b` | DOUBLE | Query-to-parent-B identity % |
| `id_a_b` | DOUBLE | Parent-A-to-parent-B identity % |
| `id_query_top` | DOUBLE | Query-to-closest-parent identity % |
| `left_yes` | INTEGER | Left segment yes votes |
| `left_no` | INTEGER | Left segment no votes |
| `left_abstain` | INTEGER | Left segment abstain votes |
| `right_yes` | INTEGER | Right segment yes votes |
| `right_no` | INTEGER | Right segment no votes |
| `right_abstain` | INTEGER | Right segment abstain votes |
| `divergence` | DOUBLE | Model divergence (id_query_model - id_query_top) |
| `flag` | VARCHAR | Classification: `Y` (chimera), `N` (non-chimera), `?` (borderline) |

**Identifier-column types (`read_id`, `parent_a_id`, `parent_b_id`, `closest_parent_id`):**
- The query and reference `read_id` columns may each be `VARCHAR`, `BIGINT`, or `UUID`, independently. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` mirrors the query table's type; the three parent columns (reference labels) mirror the reference table's type, and they can differ (e.g. a `UUID` query against a `BIGINT` reference). `BIGINT`/`UUID` ids round-trip through their decimal / canonical-lowercase forms. The parent columns stay SQL NULL for non-chimeric rows regardless of id type.

**Behavior:**
- One output row per query sequence (all queries are reported, not just chimeras)
- Non-chimeric sequences (flag=`N`) have NULL for parent and identity columns, and 0 for vote columns
- `id_a_b` (parent A vs parent B identity) is only computed for chimeric (`Y`) and borderline (`?`) results. Non-chimeric results report `id_a_b=0.0`. This avoids an extra pairwise alignment per query that is not needed for classification. Note: vsearch computes `id_a_b` for all queries unconditionally.
- Multi-threaded: vsearch's internal `chimera_detect_batch` parallelizes across the `threads` parameter (defaults to DuckDB's configured thread count; override per-call with `threads:=N`)
- The reference database is fully materialized in memory at init time
- Tables and views are both supported for query and reference inputs

**Algorithm:**
1. For each query, partition into 4 chunks and search the reference DB using an 8-mer index for candidate parents (up to 16)
2. Align query to each candidate using vsearch's SIMD-optimized Needleman-Wunsch global alignment
3. Select the 2 best parents via smoothed identity (32bp sliding window)
4. Build a 3-way star alignment and classify each column (match-A, match-B, no-vote, abstain)
5. Sweep all breakpoints left-to-right, computing h-score at each position
6. Classify based on h-score, divergence, and minimum diff thresholds

**Examples:**
```sql
-- Load sequences from FASTA files into tables
CREATE TABLE refs AS SELECT read_id, sequence1 FROM read_fastx('gold.fasta');
CREATE TABLE queries AS SELECT read_id, sequence1 FROM read_fastx('amplicons.fasta');

-- Detect chimeras
SELECT * FROM detect_chimera_uchime('queries', db:='refs');

-- Filter chimeric sequences
CREATE TABLE clean_seqs AS
SELECT q.* FROM queries q
JOIN detect_chimera_uchime('queries', db:='refs') u ON q.read_id = u.read_id
WHERE u.flag = 'N';

-- Count chimeras
SELECT flag, count(*) FROM detect_chimera_uchime('queries', db:='refs') GROUP BY flag;
```

Sequences are loaded here with [`read_fastx`](reading.md); for finding best-matching references more generally, see [sequence search](search.md).

**Error conditions:**
- Error if `db` parameter is missing
- Error if query or reference table does not exist or lacks `read_id`/`sequence1` columns
- Error if reference table is empty
- Error if query table contains NULL `read_id` values
- Error if scoring parameters are out of valid range

### De novo (UCHIME)

De novo chimera detection using the UCHIME algorithm, detecting chimeric sequences without a reference database by using abundance information: more abundant sequences are assumed to be non-chimeric and serve as parents for less abundant sequences.

**Function signature**:

`detect_chimera_uchime_denovo(input_table, [sample_id='col'], [options])`

**Parameters:**
- `input_table` (VARCHAR): Name of a table or view containing sequences with abundance. By default must have `read_id` (VARCHAR, BIGINT, or UUID — see *Identifier-column types* below), `sequence1` (VARCHAR), and `size` (integer type) columns; use `id_col`/`sequence_col`/`count_col` to override.
- `sample_id` (VARCHAR, optional): Name of a column in `input_table` to partition by. Each sample gets its own k-mer index and bootstrap; a read_id that appears in multiple samples is therefore scored independently. The sample column is prepended to the output. Execution is serialized per the vsearch wrapper's thread-safety constraints. A TEMP table or view cannot be used as the source in this mode: the per-sample path builds a fixed-name temporary view on a private per-thread connection, which must stay isolated so parallel workers do not collide on that name. Materialize the source as an ordinary table first. Tracked as #207.
- `id_col` (VARCHAR, default `'read_id'`): Name of the read identifier column in `input_table`.
- `sequence_col` (VARCHAR, default `'sequence1'`): Name of the sequence column.
- `count_col` (VARCHAR, default `'size'`): Name of the per-sequence count column. Set to `'abundance'` to chain `deblur(...)` directly into this function.
- `abskew` (DOUBLE, default 2.0): Abundance skew. Candidate parents must have abundance >= abskew * query abundance. Must be >= 1.0.
- `minh`, `xn`, `dn`, `mindiv`, `mindiffs`: Same as [reference-based UCHIME](#reference-based-uchime). (No `threads` parameter — de novo detection is sequential by construction; vsearch is run with `opt_threads=1`.)

**Output schema:** Same 18 columns as [reference-based UCHIME](#reference-based-uchime).

**Identifier-column types:** The `id_col` column may be `VARCHAR`, `BIGINT`, or `UUID` (other numeric types — INTEGER, UBIGINT, HUGEINT, DOUBLE — are rejected at bind time with `Column '<id_col>' in table '<table>' must be VARCHAR, BIGINT, or UUID`). De novo detection has a single id source, so `read_id` and all three parent columns share that one type; the parent columns stay NULL for non-chimeric rows. `BIGINT`/`UUID` ids round-trip through their decimal / canonical-lowercase forms. The `;size=` abundance annotation is never embedded in the id — the count comes from `count_col`, so a numeric id never collides with abundance.

**Behavior:**
- Sequences are processed in decreasing abundance order (highest first)
- The first two sequences (most abundant) are unconditionally treated as non-chimeric to seed the reference database
- Non-chimeric and borderline sequences are added to the reference DB incrementally
- Chimeric sequences are NOT added to the reference DB
- Single-threaded (inherently sequential — each result depends on previous classifications)
- One output row per input sequence

**Examples:**
```sql
-- Load sequences and compute abundance
CREATE TABLE seqs AS
SELECT read_id, sequence1, count(*) AS size
FROM read_fastx('amplicons.fasta')
GROUP BY read_id, sequence1;

-- De novo chimera detection
SELECT * FROM detect_chimera_uchime_denovo('seqs');

-- Filter out chimeras
SELECT s.* FROM seqs s
JOIN detect_chimera_uchime_denovo('seqs') u ON s.read_id = u.read_id
WHERE u.flag != 'Y';
```

The `count_col := 'abundance'` option lets you chain the output of [denoising](denoising.md) directly into de novo chimera detection.

**Error conditions:**
- Error if table does not exist or lacks the resolved id/sequence/count columns
- Error if the count column is not an integer type
- Error if table is empty
- Error if scoring parameters are out of valid range
- Error if any of `id_col`/`sequence_col`/`count_col` is the empty string
