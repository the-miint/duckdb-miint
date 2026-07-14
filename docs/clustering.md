# Sequence clustering

Group similar sequences into clusters (OTUs) at a chosen identity threshold, picking one representative "centroid" per cluster. Use this when you want to collapse a set of sequences (e.g. amplicons or reads) into operational taxonomic units, dereplicate near-duplicates, or reduce a large sequence set to a smaller set of representatives before downstream analysis.

## Table of Contents

- [Greedy sequence clustering](#greedy-sequence-clustering) - Cluster sequences by identity, assigning each to a centroid.

### Greedy sequence clustering

Greedy sequence clustering, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Clusters sequences by iterating in input order: each sequence is compared against existing centroids, and either joins the best matching cluster (if above the identity threshold) or becomes a new centroid. Requires the extension to be built with vsearch support.

**Function signature**:

`cluster_sequences_vsearch(input_table, id=threshold, [options])`

**Parameters:**
- `input_table` (VARCHAR): Name of a table or view containing sequences. Must have `read_id` (VARCHAR, BIGINT, or UUID — see *Identifier-column types* below) and `sequence1` (VARCHAR) columns.
- `id` (DOUBLE, required): Minimum identity threshold (0.0-1.0). No silent default — must be specified explicitly.
- `strand` (VARCHAR, default `'plus'`): `'plus'` for plus-strand only, `'both'` to also search reverse complements.
- `threads` (INTEGER, optional): Number of threads vsearch uses for its internal `cluster_assign_batch` parallel scan. Defaults to DuckDB's configured thread count (`SET threads=N`) at bind time; pass an explicit value to override. Must be 1–1024 (matching vsearch's CLI ceiling).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | VARCHAR, BIGINT, or UUID | Input sequence identifier (mirrors the input `read_id` type) |
| `is_centroid` | BOOLEAN | True if this sequence started a new cluster |
| `cluster_id` | INTEGER | Cluster number (0-based) |
| `centroid_id` | VARCHAR, BIGINT, or UUID | Identifier of the cluster's centroid (mirrors the input `read_id` type) |
| `identity` | DOUBLE | Percent identity to centroid (100.0 if centroid) |
| `cigar` | VARCHAR | CIGAR alignment to centroid (empty if centroid) |
| `cigar_truncated` | BOOLEAN | True if CIGAR was truncated (>4096 chars) |

**Identifier-column types (`read_id`, `centroid_id`):**
- The input `read_id` column may be `VARCHAR`, `BIGINT`, or `UUID`. Other numeric types (INTEGER, UBIGINT, HUGEINT, DOUBLE) are rejected at bind time with the message `Column 'read_id' in table '<table>' must be VARCHAR, BIGINT, or UUID`.
- Both output id columns mirror the input type. `centroid_id` is itself one of the input `read_id`s, so it always shares `read_id`'s type. `BIGINT` ids round-trip through their decimal form; `UUID` ids round-trip through their canonical 36-char lowercase form (uppercase input comes back lowercase).

**Sort order is the caller's responsibility.** The function clusters sequences in the order they appear in the input table. For `cluster_fast`-equivalent behavior (longest first), sort by length descending. For `cluster_size`-equivalent behavior (most abundant first), sort by abundance descending:

```sql
-- Load and sort by length descending (like vsearch --cluster_fast)
CREATE TABLE sorted_seqs AS
  SELECT * FROM read_fastx('sequences.fasta')
  ORDER BY length(sequence1) DESC;
SELECT * FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97);

-- Cluster by abundance (like vsearch --cluster_size)
CREATE TABLE by_abundance AS
  SELECT read_id, sequence1, count(*) AS size
  FROM read_fastx('amplicons.fasta')
  GROUP BY read_id, sequence1
  ORDER BY size DESC;
SELECT * FROM cluster_sequences_vsearch('by_abundance', id:=0.97);

-- Count clusters
SELECT count(*) FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97) WHERE is_centroid;

-- Get cluster sizes
SELECT centroid_id, count(*) AS size
FROM cluster_sequences_vsearch('sorted_seqs', id:=0.97)
GROUP BY centroid_id ORDER BY size DESC;
```

**Behavior:**
- One output row per input sequence
- Single-threaded (inherently sequential — each centroid must be indexed before the next sequence is processed)
- RNA sequences (U) are automatically converted to DNA (T)
- All results are materialized before returning (vsearch session mutex held for the duration)

**Error conditions:**
- Error if `id` parameter is missing
- Error if `id` is not between 0.0 and 1.0
- Error if `strand` is not `'plus'` or `'both'`
- Error if table does not exist or lacks required columns
- Error if table is empty or contains NULL read_ids, NULL sequences, or empty sequences
- Error if any read_id exceeds 1023 characters

## See also

- [Reading sequences](reading.md) — load FASTA/FASTQ into the `read_id` / `sequence1` schema this function expects (`read_fastx`).
- [Sequence search](search.md) — find the best-matching reference sequences per query by global-alignment identity (`search_sequences_vsearch`, also requires vsearch).
- [Chimera checking](chimera.md) — detect chimeric sequences before or after clustering (also requires vsearch).
