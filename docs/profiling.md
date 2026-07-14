# Profiling and feature tables

Estimate *who* is in a community and *how much* of each is present, and turn those estimates into a feature table (a `sample × feature` abundance matrix). Two complementary methods: `woltka_ogu` builds an OGU feature table from reads you have already aligned to a reference, while `sylph_profile` profiles raw shotgun reads directly against a pre-built reference database. Both emit long-form output that feeds straight into [diversity analysis](diversity.md) and can be written to BIOM with [`COPY ... (FORMAT BIOM)`](writing.md).

## Table of Contents

- [`woltka_ogu`](#woltka_ogu) - OGU (Operational Genomic Unit) feature table from alignments, global or per-sample.
- [`sylph_profile`](#sylph_profile) - FracMinHash relative-abundance profiling of shotgun reads against a `.syldb` database.

### `woltka_ogu`

Compute [Woltka](https://github.com/qiyunzhu/woltka) OGU (Operational Genomic Unit) counts over SAM-like alignment data. Implements Woltka's classification algorithm, assigning reads to taxonomic units while fractionally distributing multi-mapped reads. When the optional `sample_id` named parameter is supplied, the aggregation runs in parallel across distinct sample values — one DuckDB query per sample on a dedicated per-thread connection — which bounds memory to a single sample's footprint.

**Function signature**:

`woltka_ogu(relation, sequence_id_field [, sample_id])`

**Parameters:**
- `relation` (VARCHAR): Name of a table or view containing SAM-like alignment data (must be a resolvable catalog name; pass as a string literal)
- `sequence_id_field` (VARCHAR): Name of the column holding sequence identifiers — typically `read_id`, or a numeric index column for better hash performance
- `sample_id` (VARCHAR, named, optional): Name of the column holding sample identifiers. When supplied, adds that column to the output and runs the aggregation per distinct sample value in parallel.

**Required columns in relation:**
- Column named by `sequence_id_field`: read/sequence identifier
- `reference` (VARCHAR): reference sequence name (becomes `feature_id`)
- `flags` (USMALLINT): SAM alignment flags
- When `sample_id` is supplied: the named column (any comparable type) — NULLs are rejected at bind time

**Output schema:**
- When `sample_id` is omitted: `(feature_id VARCHAR, value DOUBLE)`
- When `sample_id` is supplied: `(<sample_id_column> <its_type>, feature_id VARCHAR, value DOUBLE)` — the first column's name matches the value you passed to `sample_id`.

This long-form `(sample_id, feature_id, value)` shape is exactly what the [diversity functions](diversity.md) consume and what [`COPY ... (FORMAT BIOM)`](writing.md) writes.

**Algorithm:**
1. Orients reads using alignment flags (forward/reverse via [`alignment_is_read1`](alignment_analysis.md#sam-flag-functions)).
2. For each read orientation, divides 1 by the number of unique features aligned to.
3. Aggregates fractional counts per feature (and per sample, when `sample_id` is used).

**Correctness assumption for per-sample mode:** read IDs are unique across samples. The per-sample subset then yields the same distribution as the global aggregation.

**Examples:**
```sql
-- Global aggregation across the whole relation
SELECT * FROM woltka_ogu('my_alignments', 'read_id');

-- Per-sample aggregation — one aggregation per distinct sample value, in parallel
SELECT * FROM woltka_ogu(
    'my_alignments',
    'read_id',
    sample_id := 'sample_id'
);

-- Filter high-quality alignments via a view, then classify
CREATE OR REPLACE VIEW primary_alignments AS
    SELECT * FROM read_alignments('alignments.bam')
    WHERE alignment_is_primary(flags) AND mapq >= 20;

SELECT * FROM woltka_ogu('primary_alignments', 'read_id');

-- Multi-sample: union sources under a sample_id column, then classify per sample
CREATE OR REPLACE VIEW all_samples AS
    SELECT *, 'sample1' AS sample_id FROM read_alignments('sample1.bam')
    UNION ALL
    SELECT *, 'sample2' AS sample_id FROM read_alignments('sample2.bam')
    UNION ALL
    SELECT *, 'sample3' AS sample_id FROM read_alignments('sample3.bam');

SELECT * FROM woltka_ogu('all_samples', 'read_id', sample_id := 'sample_id')
ORDER BY sample_id, feature_id;

-- Export per-sample results to BIOM for downstream analysis
COPY (
    SELECT * FROM woltka_ogu('my_alignments', 'read_id', sample_id := 'sample_id')
) TO 'ogu_table.biom' (FORMAT BIOM, COMPRESSION 'gzip');

-- Export global results to BIOM (add a sample_id column)
COPY (
    SELECT feature_id, 'MySample' AS sample_id, value
    FROM woltka_ogu('my_alignments', 'read_id')
) TO 'ogu_single.biom' (FORMAT BIOM);
```

**Behavior:**
- Arguments `relation` and `sequence_id_field` (and `sample_id` if used) must be quoted string literals; they are catalog/column names resolved at bind time.
- Multi-mapped reads are fractionally assigned: each mapping receives weight 1/N where N is the number of unique references the read maps to, within the same orientation.
- Read orientation (forward/reverse) is considered separately via the [`alignment_is_read1()`](alignment_analysis.md#sam-flag-functions) flag — paired-end reads are handled automatically.
- For better performance on large datasets, add a numeric index column and pass it as `sequence_id_field` instead of `read_id`.
- Output row order is non-deterministic when `sample_id` is used (parallel per-sample execution). Use an explicit `ORDER BY` if stable ordering is required.

### `sylph_profile`

FracMinHash-based relative-abundance profiling of paired-end shotgun metagenomic reads against a pre-built `.syldb` reference database, using [sylph](https://github.com/bluenote-1577/sylph) (Shaw & Yu 2024, *Nature Biotechnology*). Sequences come from a DuckDB table or view — there is no FASTQ path argument — and the database is loaded once per call, mmap-backed. The result is streamed back via the Arrow C Data Interface (zero-copy).

Embedded as a Rust static library (sylph 0.9.0-miint fork; MIT). Linux and macOS only; the function is not registered on WASM or Windows builds.

**Function signature**:

`sylph_profile(source_table, syldb_path, [sample_id='col'], [options])`

**Parameters:**

- `source_table` (VARCHAR, positional): Name of a table or view with columns `read_id` (VARCHAR), `sequence1` (VARCHAR), and optionally `sequence2` (VARCHAR). When `sequence2` is present and non-empty per row, the read pair is processed paired-end; otherwise the call is single-end.
- `syldb_path` (VARCHAR, positional): Path to a sylph `.syldb` reference database, built offline via the upstream `sylph sketch` CLI. The file is opened read-only and shared mmap-style across the call.
- `sample_id` (VARCHAR, optional): Name of a column on `source_table` to partition by. When set, the function fans out per-sample (parallelized via the per-sample helper used by [`deblur`](denoising.md) / [`align_mafft`](alignment_multiple.md) / [`detect_chimera_uchime_denovo`](chimera.md)) and prepends a `sample_id` column to the output. Without this option, the entire table is processed as a single sample.
- `min_ani` (DOUBLE, default = sylph default): Minimum adjusted ANI (percent, 0..100) for a genome to be reported. Negative or unset = use sylph's built-in default.
- `min_number_kmers` (UINTEGER, default = 50): Minimum number of matching k-mers required to report a genome.
- `min_count_correct` (DOUBLE, default = 3.0): Minimum corrected k-mer count for genome detection. Lower values increase sensitivity at the cost of false positives.
- `estimate_unknown` (BOOLEAN, default false): Renormalize `taxonomic_abundance` to fraction-of-reads-explained, accounting for unknown / unmatched fraction.
- `dedup_paired_reads` (BOOLEAN, default true): Enable approximate paired-read deduplication during sketching. Matches sylph CLI behavior.
- `dedup_fpr` (DOUBLE, default 0.0001): Bloom-filter false-positive rate used during sketch-time deduplication. Matches `--fpr` on the upstream CLI.
- `threads` (UINTEGER, default 0 = auto): Rayon thread count for the inner sylph compute. In per-sample mode, the inner count auto-balances against DuckDB's worker count so that `outer × inner ≈ db_threads` — same total core utilization as `sylph profile -t db_threads`.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | (matches source column type) | Only present when `sample_id` is set; the partition value for the row |
| `genome_index` | UINTEGER | 0-based row offset of the genome within the `.syldb` |
| `genome_name` | VARCHAR | Genome identifier from the sketch (typically the reference FASTA path or accession) |
| `contig_name` | VARCHAR | First contig name from the genome's sketch — disambiguator when multiple genomes share `genome_name` |
| `sequence_abundance` | DOUBLE | Fraction of the sample's reads explained by this genome (0..1) |
| `taxonomic_abundance` | DOUBLE | Coverage-corrected relative abundance (sums to ~1.0 across all reported genomes, or to "explained fraction" when `estimate_unknown := true`) |
| `adjusted_ani` | DOUBLE | Lambda-adjusted average nucleotide identity (percent) between sample reads and the genome sketch |
| `eff_cov` | DOUBLE | Effective sequencing coverage of the genome by the sample |
| `naive_ani` | DOUBLE | Unadjusted ANI (percent) — useful for comparison against `adjusted_ani` to spot low-coverage / high-correction cases |
| `kmers_reassigned` | UBIGINT | Number of k-mers reassigned away from this genome by the winner-take-all step |

Adding columns at the end is non-breaking; reordering is.

**Behavior:**

- **Single-vs-paired:** auto-detected per-row. Rows whose `sequence2` is NULL or empty are sketched single-end; rows with both sequences set are sketched paired. Mixing within one sample is allowed.
- **Database lifecycle:** the `.syldb` is loaded once per `sylph_profile` call into the GlobalState and reused across all execute invocations on that call. Back-to-back calls against the same DB pay the load cost each time (a future session-scoped cache is on the roadmap).
- **Ordering:** `order_preservation_type = NO_ORDER` — rows may interleave across genomes when running per-sample. Sort downstream if needed.
- **Numerical parity:** bit-identical to upstream `sylph profile --reads` on the same data, modulo a small documented set of divergences in the embedded library.

**Examples:**

```sql
-- Profile a metagenome against the GTDB sylph DB (single sample).
CREATE TABLE reads AS
  SELECT read_id, sequence1, sequence2
  FROM read_fastx('sample_R1.fq.gz', sequence2 := 'sample_R2.fq.gz');

SELECT genome_name, taxonomic_abundance, adjusted_ani, eff_cov
FROM sylph_profile('reads', '/refs/gtdb-r220-c200-dbv1.syldb')
ORDER BY taxonomic_abundance DESC
LIMIT 20;
```

```sql
-- Per-sample profiling: one row stream per sample_accession, parallelized.
CREATE TABLE reads AS
  SELECT sample_accession, read_id, sequence1, sequence2
  FROM 'PRJNA000000.parquet';

SELECT *
FROM sylph_profile('reads', '/refs/gtdb-r220-c200-dbv1.syldb',
                   sample_id := 'sample_accession',
                   estimate_unknown := true);
```

**Error conditions:**

- Error if `source_table` does not exist or is missing required columns (`read_id`, `sequence1`).
- Error if `syldb_path` cannot be opened, is corrupt, or is not a sylph `.syldb` file.
- Error if `sample_id` is set but the named column does not exist on `source_table` or collides with an output column name.
- IO error surfaces the underlying sylph diagnostic string (e.g., truncated database, unsupported version).

## See also

- [Reading files](reading.md) — `read_fastx` (reads) and `read_alignments` (alignments) that feed these functions.
- [Writing files](writing.md) — exporting a feature table to BIOM.
- [Alpha and beta diversity](diversity.md) — consuming the long-form feature table produced here.
- [Taxonomic classification](classification.md) — per-sequence k-mer classification with RYpe.
