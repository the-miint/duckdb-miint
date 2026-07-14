# Sequence search

Find, for each of your query sequences, the best-matching reference sequences in a database by global-alignment percent identity. Use this when you have a set of unknown sequences (e.g. reads or amplicons) and want to assign each one to the closest entry in a curated reference set at a chosen identity cutoff.

## Table of Contents

- [Global sequence search](#global-sequence-search) - Search query sequences against a reference database by global-alignment identity.

### Global sequence search

Global pairwise sequence search, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Finds the best matching sequences in a reference database for each query sequence using SIMD-optimized Needleman-Wunsch alignment with k-mer candidate filtering. Requires the extension to be built with vsearch support.

**Function signature**:

`search_sequences_vsearch(query_table, db='ref_table', id=threshold, [options])`

**Parameters:**
- `query_table` (VARCHAR): Name of a table or view containing query sequences. Must have `read_id` (VARCHAR, BIGINT, or UUID — see *Identifier-column types* below) and `sequence1` (VARCHAR) columns.
- `db` (VARCHAR, required): Name of a table or view containing reference sequences. Same schema requirements as `query_table` (its `read_id` type is independent of the query's).
- `id` (DOUBLE, required): Minimum identity threshold (0.0-1.0). No silent default — must be specified explicitly.
- `maxaccepts` (INTEGER, default 1): Maximum number of accepted hits per query. Must be >= 1.
- `maxrejects` (INTEGER, default 32): Maximum rejected targets before stopping search. Must be >= 1.
- `threads` (INTEGER, optional): Number of threads vsearch uses for its internal `search_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR, BIGINT, or UUID | Query sequence identifier (mirrors `query_table.read_id`) |
| `target_id` | VARCHAR, BIGINT, or UUID | Reference sequence identifier (mirrors `db.read_id`) |
| `identity` | DOUBLE | Percent identity (0-100) |
| `matches` | INTEGER | Number of matching columns |
| `mismatches` | INTEGER | Number of mismatching columns |
| `gaps` | INTEGER | Number of gap columns |
| `alignment_length` | INTEGER | Total alignment length |
| `query_length` | INTEGER | Query sequence length |
| `target_length` | INTEGER | Target sequence length |
| `accepted` | BOOLEAN | True if hit passes identity threshold |

**Identifier-column types (`read_id`, `target_id`):**
- The query and reference `read_id` columns may each be `VARCHAR`, `BIGINT`, or `UUID`, independently. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- The output `read_id` mirrors the query table's type; `target_id` (a reference label) mirrors the reference table's type. They can differ (e.g. a `UUID` query against a `BIGINT` reference). `BIGINT`/`UUID` ids round-trip through their decimal / canonical-lowercase forms. Search emits one row per hit, so neither id is ever NULL.

**Examples:**
```sql
CREATE TABLE refs AS SELECT read_id, sequence1 FROM read_fastx('database.fasta');
CREATE TABLE queries AS SELECT read_id, sequence1 FROM read_fastx('queries.fasta');

-- Search at 97% identity
SELECT * FROM search_sequences_vsearch('queries', db:='refs', id:=0.97);

-- Top 3 hits per query at 90% identity
SELECT * FROM search_sequences_vsearch('queries', db:='refs', id:=0.90, maxaccepts:=3);

-- Count queries with hits
SELECT count(DISTINCT read_id) FROM search_sequences_vsearch('queries', db:='refs', id:=0.97);
```

**Behavior:**
- Each query produces 0 to `maxaccepts` output rows
- Results include both accepted (above threshold) and weak (near-miss) hits; use the `accepted` column to filter
- Plus-strand only (no reverse complement search)
- Multi-threaded: vsearch's internal `search_batch` parallelizes across the `threads` parameter (defaults to DuckDB's configured thread count; override per-call with `threads:=N`)
- Reference database is fully materialized in memory at init time
- RNA sequences (U) are automatically converted to DNA (T)

**Error conditions:**
- Error if `db` or `id` parameter is missing
- Error if `id` is not between 0.0 and 1.0
- Error if query or reference table does not exist or lacks required columns
- Error if reference table is empty

## See also

- [Reading sequences](reading.md) — load FASTA/FASTQ into the `read_id` / `sequence1` schema these tables expect (`read_fastx`).
- [Sequence clustering](clustering.md) — group sequences into clusters by identity (`cluster_sequences_vsearch`, also requires vsearch).
