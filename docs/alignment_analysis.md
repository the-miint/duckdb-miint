# Alignment analysis

A considerable amount of analysis on alignment data can be performed with native SQL. MIINT further provides a set of optimized methods to support specific complex operations.

## Table of contents

- [SAM flag functions](#sam-flag-functions) - Test individual SAM flag bits
- [Per-base pileup](#per-base-pileup) - A per-base CIGAR walker.
- [Alignment slicing](#alignment-slicing) - Slice an alignment based on start and stop coordinates with implicit CIGAR update.
- [Sequence identity](#sequence-identity) - Sequence identity calculation (multi-mode, requires NM/MD for legacy CIGAR)
- [Sequence identity from CIGAR](#sequence-identity-from-cigar) - Sequence identity from extended CIGAR alone
- [Query length](#cigar-query-length) - Query length from CIGAR
- [Query coverage](#cigar-query-coverage) - Query coverage from CIGAR
- [Merge overlapping intervals](#merge-overlapping-intervals) - Merge overlapping genomic intervals (aggregate)
- [Coverage depth](#coverage-depth) - Per-position depth of coverage (aggregate)
- [Genome coverage](#genome-coverage) - Proportion of each genome covered by alignments
- [Barcode matching](#barcode-matching) - Hamming-distance matcher for short fixed-length barcodes
- [MSA column consensus](#msa-column-consensus) - Quality-aware consensus from a multiple alignment

### SAM flag functions

Test individual SAM flag bits. Each function takes a `USMALLINT` (the flags column from `read_alignments`) and returns a `BOOLEAN`.

**Primary names:**
- `alignment_is_paired(flags)` - Read is paired (0x1)
- `alignment_is_proper_pair(flags)` - Read is properly paired (0x2)
- `alignment_is_unmapped(flags)` - Read is unmapped (0x4)
- `alignment_is_mate_unmapped(flags)` - Mate is unmapped (0x8)
- `alignment_is_reverse(flags)` - Read is reverse strand (0x10)
- `alignment_is_mate_reverse(flags)` - Mate is reverse strand (0x20)
- `alignment_is_read1(flags)` - Read is first in pair (0x40)
- `alignment_is_read2(flags)` - Read is second in pair (0x80)
- `alignment_is_secondary(flags)` - Secondary alignment (0x100)
- `alignment_is_primary(flags)` - Primary alignment (neither secondary nor supplementary)
- `alignment_is_qc_failed(flags)` - QC failure (0x200)
- `alignment_is_duplicate(flags)` - PCR/optical duplicate (0x400)
- `alignment_is_supplementary(flags)` - Supplementary alignment (0x800)

**HTSlib-compatible aliases:**
`is_paired`, `is_proper_pair`, `is_unmapped`, `is_munmap`, `is_reverse`, `is_mreverse`, `is_read1`, `is_read2`, `is_secondary`, `is_qcfail`, `is_dup`, `is_supplementary`

**Example:**
```sql
SELECT read_id, flags
FROM read_alignments('alignments.sam')
WHERE alignment_is_paired(flags)
  AND NOT alignment_is_unmapped(flags);
```

### Per-base pileup

Per-base CIGAR walker. Expands each alignment into one row per reference position covered, with the reference base, the query base, and the per-base query quality. Replaces `samtools view | bcftools mpileup` for single-sample variant-call positions inside SQL pipelines.

**Function signature**:

`compute_pileup(alignments_table, reference_table)`

CIGAR op handling:
- `M` / `=` / `X` — emit one row per ref position; advance both cursors.
- `D` / `N` — emit one row per ref position with `query_base = NULL` and
  `query_qual = NULL`; advance ref only.
- `I` — emit one row per inserted base with `ref_base = NULL`,
  `ref_pos` = preceding reference position, `insert_pos` = 1..N.
- `S` — soft-clipped bases are dropped (advance query only).
- `H` / `P` — no-op.

Reuses miint's existing CIGAR parser, not htslib. The reference base at
each position is looked up from the supplied `reference_table` — no MD-tag
input is required.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `alignments_table` | VARCHAR | Table with `(read_id, reference, position, cigar, sequence, qual)` columns. Compatible with the schema produced by `read_alignments(include_seq_qual := true)`. |
| `reference_table` | VARCHAR | Table with `(ref_id, sequence)` columns. Provides the reference sequence per contig. |

`position` is 1-based (matches SAM POS). `reference` is the contig name
that should appear as `ref_id` in the reference table.

**Returns:** TABLE with columns:

| Column | Type | Description |
|--------|------|-------------|
| `ref_id` | VARCHAR | Reference contig name |
| `ref_pos` | BIGINT | 1-based reference position (for insertion rows: preceding ref position) |
| `read_id` | VARCHAR | The originating read |
| `ref_base` | VARCHAR | Reference base at this position (NULL on insertion rows) |
| `query_base` | VARCHAR | Read base at this position (NULL on D/N) |
| `query_qual` | UTINYINT | Read qual at this position (NULL on D/N or when input qual is NULL) |
| `insert_pos` | INTEGER | 0 for reference-aligned rows; 1-based position within an insertion event |

**Behavior:**
- Throws if an alignment references a contig missing from
  `reference_table`, or if a CIGAR's reference span runs past the end of
  the reference.
- Throws on CIGAR / seq length mismatches.
- Empty / `*` CIGAR rows produce zero output rows (silently skipped).
- Insertion rows use `ref_pos` = the preceding reference-consuming position
  (SAM convention). For leading insertions (before any ref-consuming op),
  `ref_pos` = `position - 1`. Use `WHERE insert_pos = 0` to filter out
  insertion rows and recover reference-aligned-only output.
- Alignments are streamed (not fully materialized). Memory usage scales with
  the pileup rows per alignment chunk, not the total number of alignments.

**Example:**
```sql
-- Per-position pileup 
SELECT ref_pos, query_base, query_qual
FROM compute_pileup('bin_a_alignments', 'bin_a_reference')
WHERE ref_id = 'binA_cons'
ORDER BY ref_pos, read_id;

-- Allele counts per position (filter out NULL query bases from deletions)
SELECT ref_pos, query_base, COUNT(*) AS n
FROM compute_pileup('alignments', 'reference')
WHERE query_base IS NOT NULL
GROUP BY ref_pos, query_base
ORDER BY ref_pos, query_base;
```


### Alignment slicing

Slice alignment data from a table or view to a genomic region. Each alignment is trimmed to the specified region `[start, stop)`, with trimmed portions represented as hard clips (H) in the output CIGAR.

**Function signature:**

`alignment_slice(table_name, start, stop, [include_deletions=false])`

**Parameters:**
- `table_name` (VARCHAR): Name of a table or view containing alignment data
- `start` (BIGINT): Region start position (1-based, inclusive)
- `stop` (BIGINT): Region stop position (1-based, exclusive)
- `include_deletions` (BOOLEAN, default `false`): Whether deletion (D) operations count as overlap evidence. When `false`, reads whose only overlap with the region is through deletions are excluded.

**Coordinates:** 1-based half-open `[start, stop)`, consistent with `position` and `stop_position` from `read_alignments`.

**Required input columns:** `cigar` (VARCHAR), `position` (BIGINT), `stop_position` (BIGINT)

**Output schema:** Same columns as found in the input table (from the recognized alignment column set), with adjusted values for overlapping reads.

**Behavior:**
- Reads that don't overlap the region are excluded
- Soft clips (S) do not count as overlap evidence
- Trimmed portions of the CIGAR become hard clips (H)
- Tags (`tag_as` through `tag_sa`) are set to NULL when trimming occurs
- `template_length` is set to NULL when trimming occurs
- `mapq` and mate fields are preserved
- If the input table has a `reference` column, all rows must have the same reference (single-region slicing)
- Rows with NULL `cigar`, `position`, or `stop_position` are skipped

**Examples:**

```sql
-- Load alignments and filter to a single reference
CREATE VIEW chr1_alns AS
  SELECT * FROM read_alignments('sample.bam')
  WHERE reference = 'chr1';

-- Slice to a region of interest
SELECT * FROM alignment_slice('chr1_alns', 1000, 2000);

-- Include reads that overlap only via deletions
SELECT * FROM alignment_slice('chr1_alns', 1000, 2000, include_deletions := true);

-- Composable: compute coverage depth on sliced alignments
SELECT compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')
FROM alignment_slice('chr1_alns', 1000, 2000);
```

**Notes:**
- The function reads from the input table via a separate connection, so uncommitted changes in the current transaction may not be visible
- Multi-region slicing (different regions per reference) is not yet supported; use separate queries per region
- `alignment_seq_identity` with the `'cigar'` method works on sliced output because it reads identity directly from `=`/`X` CIGAR ops without needing tags. Other methods (`gap_compressed`, `blast`, `gap_excluded`) require NM or MD tags which are NULLed after trimming.

### Sequence identity 

Calculate sequence identity between read and reference using three different methods. They are derived from Heng Li's [blog post](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

**Function signature:**
`alignment_seq_identity(cigar, nm, md, type)`

*IMPORTANT*: these calculations assume use of '='/'X' (e.g., `--xeq` with bowtie2, `-eqx` with minimap2, etc).

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment.
- `nm` (BIGINT): NM tag value (edit distance)
- `md` (VARCHAR): MD tag value (mismatch positions)
- `type` (VARCHAR, default='gap_compressed'): Identity calculation method

**Identity types:**

1. **`'gap_excluded'`**: Ignore gaps, only consider match/mismatch positions
   - Formula: `#matches / (#matches + #mismatches)`
   - Requires: CIGAR + MD tag
   - Use case: Genetic divergence between species

2. **`'blast'`**: Traditional BLAST-style identity
   - Formula: `#matches / alignment_columns`
   - Requires: CIGAR + NM tag
   - Use case: General similarity measurement
   - Note: Large indels significantly lower identity

3. **`'gap_compressed'`** (default): Count consecutive gaps as single events
   - Formula: `(m - n + g) / (m + o)` where m=M_columns, n=NM, g=gap_bases, o=gap_opens
   - Equivalent to: `1 - (n - g + o) / (m + o)` from the blog post
   - Requires: CIGAR + NM tag
   - Use case: Filtering alignments (recommended)
   - Note: More robust to structural variations

4. **`'cigar'`**: Identity from extended CIGAR operations only (no tags needed)
   - Formula: `match_ops / alignment_columns` where match_ops = count of `=` ops
   - Requires: CIGAR with `=`/`X` ops (not `M`). NM and MD tags are ignored.
   - Returns NULL if CIGAR uses `M` (ambiguous) or mixes `M` with `=`/`X` (inconsistent)
   - Use case: Computing identity on trimmed alignments from `alignment_slice`, where tags have been invalidated

**Returns:** DOUBLE between 0.0 and 1.0, or NULL for unmapped reads or missing required tags

**Example:**
```sql
-- Calculate gap-compressed identity (default)
SELECT read_id, alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS identity
FROM read_alignments('alignments.sam')
WHERE tag_nm IS NOT NULL;

-- Filter high-quality alignments
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE alignment_seq_identity(cigar, tag_nm, tag_md) > 0.95;

-- Compare different identity methods
SELECT read_id,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_excluded') AS gap_excl,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'blast') AS blast,
  alignment_seq_identity(cigar, tag_nm, tag_md, 'gap_compressed') AS gap_comp
FROM read_alignments('alignments.sam')
WHERE tag_nm IS NOT NULL AND tag_md IS NOT NULL;

-- Compute identity on sliced alignments (tags NULLed, cigar method still works)
SELECT read_id, alignment_seq_identity(cigar, tag_nm, tag_md, 'cigar') AS identity
FROM alignment_slice('my_alignments', 1000, 2000);
```

**Reference:** [On the definition of sequence identity](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity) by Heng Li

### Sequence identity from CIGAR

One-arg convenience wrapper over `alignment_seq_identity(cigar, NULL, NULL, 'cigar')`. Use this when your CIGAR uses extended `=`/`X` ops and you don't need NM/MD-based identity flavors.

**Function signature:**

`cigar_sequence_identity(cigar)`

**Parameters:**
- `cigar` (VARCHAR): CIGAR string with extended ops (`=` for match, `X` for mismatch). Aligners must be invoked with `--xeq` (bowtie2) or `-eqx` (minimap2) for this to work; legacy `M`-only CIGAR cannot be used.

**Formula:** `match_ops / alignment_columns` where `match_ops` is the count of `=` operations and `alignment_columns` is `M + I + D` (with `=`/`X` counting as `M`).

**Returns:** DOUBLE in [0.0, 1.0], or NULL when identity can't be determined from CIGAR alone:
- CIGAR is NULL, empty, or `*` (unmapped)
- CIGAR uses only `M` (legacy — can't distinguish matches from mismatches; use `alignment_seq_identity` with NM or MD instead)
- CIGAR mixes `M` with `=`/`X` (inconsistent encoding)
- CIGAR has no `=`/`X` ops at all (e.g. pure `I`, `D`, `S`, `H`, `N`, `P`)

**Example:**
```sql
-- Modern CIGAR (=/X), no NM/MD needed
SELECT read_id, cigar_sequence_identity(cigar) AS identity
FROM read_alignments('alignments.sam');

-- Filter high-identity alignments
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE cigar_sequence_identity(cigar) > 0.95;

-- Same as alignment_seq_identity(cigar, NULL, NULL, 'cigar') but shorter:
SELECT read_id, cigar_sequence_identity(cigar) AS identity
FROM alignment_slice('my_alignments', 1000, 2000);
```

**See also:** `alignment_seq_identity` for legacy `M`-only CIGAR or BLAST/gap-compressed/gap-excluded identity flavors.

## `cigar_query_length(cigar, [include_hard_clips=true])`

Calculate the total query length from a CIGAR string. This is useful for understanding read lengths and query coverage.

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment
- `include_hard_clips` (BOOLEAN, default=true): Whether to include hard-clipped bases in the length

**Returns:** BIGINT - Total query length, or 0 for NULL/unmapped reads

**Behavior:**
- When `include_hard_clips=true`: Returns M + I + S + = + X + H (all query-consuming operations)
- When `include_hard_clips=false`: Returns M + I + S + = + X (matches HTSlib's `bam_cigar2qlen`)
- Soft clips (S) are always included (they're present in the sequence field)
- Hard clips (H) are only included when parameter is true
- Deletions (D) and reference skips (N) don't consume query, so they're never counted

**Examples:**
```sql
-- Get query length including hard clips (default)
SELECT read_id, cigar_query_length(cigar) AS query_len
FROM read_alignments('alignments.sam');

-- Get query length excluding hard clips (matches bam_cigar2qlen)
SELECT read_id, cigar_query_length(cigar, false) AS query_len
FROM read_alignments('alignments.sam');

-- Compare lengths with and without hard clips
SELECT read_id, cigar,
  cigar_query_length(cigar, true) AS len_with_hard,
  cigar_query_length(cigar, false) AS len_without_hard
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%H%';

-- Calculate average query length per reference
SELECT reference, AVG(cigar_query_length(cigar)) AS avg_query_len
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Note:** When `include_hard_clips=false`, this function's output matches HTSlib's `bam_cigar2qlen` behavior, which counts M, I, S, =, and X operations.

### CIGAR query length

Calculate the total query length from a CIGAR string. This is useful for understanding read lengths and query coverage.

**Function signature:**

`cigar_query_length(cigar, [include_hard_clips=true])`

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment
- `include_hard_clips` (BOOLEAN, default=true): Whether to include hard-clipped bases in the length

**Returns:** BIGINT - Total query length, or 0 for NULL/unmapped reads

**Behavior:**
- When `include_hard_clips=true`: Returns M + I + S + = + X + H (all query-consuming operations)
- When `include_hard_clips=false`: Returns M + I + S + = + X (matches HTSlib's `bam_cigar2qlen`)
- Soft clips (S) are always included (they're present in the sequence field)
- Hard clips (H) are only included when parameter is true
- Deletions (D) and reference skips (N) don't consume query, so they're never counted

**Examples:**
```sql
-- Get query length including hard clips (default)
SELECT read_id, cigar_query_length(cigar) AS query_len
FROM read_alignments('alignments.sam');

-- Get query length excluding hard clips (matches bam_cigar2qlen)
SELECT read_id, cigar_query_length(cigar, false) AS query_len
FROM read_alignments('alignments.sam');

-- Compare lengths with and without hard clips
SELECT read_id, cigar,
  cigar_query_length(cigar, true) AS len_with_hard,
  cigar_query_length(cigar, false) AS len_without_hard
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%H%';

-- Calculate average query length per reference
SELECT reference, AVG(cigar_query_length(cigar)) AS avg_query_len
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Note:** When `include_hard_clips=false`, this function's output matches HTSlib's `bam_cigar2qlen` behavior, which counts M, I, S, =, and X operations.

### CIGAR query coverage

Calculate the proportion of query bases covered by the reference alignment. This helps assess how much of a read actually aligns versus being clipped.

**Function signature:**

`cigar_query_coverage(cigar, [type='aligned'])`

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from alignment
- `type` (VARCHAR, default='aligned'): Coverage calculation method

**Coverage types:**

1. **`'aligned'`** (default): Bases that align to the reference
   - Formula: `(M + = + X) / (M + I + S + = + X + H)`
   - Only counts bases that match/mismatch the reference
   - Insertions and clips reduce coverage

2. **`'mapped'`**: Bases that are mapped (not clipped)
   - Formula: `(M + I + = + X) / (M + I + S + = + X + H)`
   - Counts insertions as "mapped" even though they don't align
   - Only clips reduce coverage

**Returns:** DOUBLE between 0.0 and 1.0, or NULL for NULL CIGAR

**Behavior:**
- Query length denominator always includes hard clips
- Returns 0.0 for reads with only clipping operations (no alignment)
- Deletions (D) and reference skips (N) don't affect coverage (they don't consume query)

**Examples:**
```sql
-- Get aligned coverage (default)
SELECT read_id, cigar_query_coverage(cigar) AS aligned_cov
FROM read_alignments('alignments.sam');

-- Get mapped coverage (includes insertions)
SELECT read_id, cigar_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam');

-- Compare aligned vs mapped coverage
SELECT read_id, cigar,
  cigar_query_coverage(cigar, 'aligned') AS aligned_cov,
  cigar_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%I%';  -- Reads with insertions show the difference

-- Filter reads with high query coverage
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE cigar_query_coverage(cigar, 'aligned') > 0.9;

-- Find heavily clipped reads
SELECT read_id, cigar, cigar_query_coverage(cigar) AS coverage
FROM read_alignments('alignments.sam')
WHERE cigar_query_coverage(cigar) < 0.5
ORDER BY coverage;

-- Calculate average coverage per reference
SELECT reference,
  AVG(cigar_query_coverage(cigar, 'aligned')) AS avg_aligned_cov,
  AVG(cigar_query_coverage(cigar, 'mapped')) AS avg_mapped_cov
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Use cases:**
- **Aligned coverage**: Assess quality of alignment (how much of read actually matches reference)
- **Mapped coverage**: Identify heavily clipped reads (adapters, chimeras, low-quality ends)
- Filter reads based on alignment quality thresholds
- QC metrics for sequencing runs

### Merge overlapping intervals

`compress_intervals(start, stop)`

Aggregate function that merges overlapping genomic intervals into a minimal set of non-overlapping intervals. Useful for computing coverage regions, reducing redundant intervals, and analyzing read depth.

**Parameters:**
- `start` (BIGINT): Start position of interval
- `stop` (BIGINT): Stop position of interval

**Returns:** `LIST<STRUCT(start BIGINT, stop BIGINT)>` — array of merged intervals, sorted by start position

**Behavior:**
- Intervals that overlap or touch (stop1 = start2) are merged into a single interval
- Automatically compresses state when accumulating >1M intervals (prevents memory issues)
- Thread-safe: works correctly with parallel GROUP BY operations
- Returns NULL for empty groups

**Examples:**
```sql
-- Calculate coverage regions per reference from SAM/BAM alignments
SELECT reference,
       compress_intervals(position, stop_position) AS coverage
FROM read_alignments('alignments.bam')
GROUP BY reference;

-- Calculate total covered bases per reference
SELECT reference,
       SUM(c.stop - c.start) AS total_coverage
FROM (
  SELECT reference,
         UNNEST(compress_intervals(position, stop_position)) AS c
  FROM read_alignments('alignments.bam')
  GROUP BY reference
);

-- Find gaps in coverage (regions with no reads)
WITH coverage AS (
  SELECT reference,
         UNNEST(compress_intervals(position, stop_position)) AS interval
  FROM read_alignments('alignments.bam')
  GROUP BY reference
)
SELECT reference,
       interval.stop AS gap_start,
       LEAD(interval.start) OVER (PARTITION BY reference ORDER BY interval.start) AS gap_end
FROM coverage
WHERE LEAD(interval.start) OVER (PARTITION BY reference ORDER BY interval.start) IS NOT NULL;
```

**Performance notes:**
- Automatic periodic compression at 1M intervals prevents memory bloat with large datasets
- Multi-threaded aggregation: each thread maintains its own state, merged at finalization
- Algorithm: sorts intervals by start position, then single-pass merge (O(n log n))

### Coverage depth

`compute_coverage_depth(position, stop_position, cigar, reference_length, mode)`

Aggregate function that computes per-position depth of coverage across a reference sequence. Follows `samtools depth` semantics: two modes for handling deletion (D) operations, and N (ref-skip) operations are always excluded.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `position` | BIGINT | 1-based inclusive start position |
| `stop_position` | BIGINT | 1-based exclusive end position (half-open, matches `read_alignments` output) |
| `cigar` | VARCHAR | CIGAR string |
| `reference_length` | BIGINT | Determines output list size (must be positive, max 2,000,000,000) |
| `mode` | VARCHAR | `'include_deletions'` or `'exclude_deletions'` (must be a constant) |

**Returns:** `LIST(UINTEGER)` — element `i` is the depth at 1-based position `i+1`.

**Mode semantics:**
- `'include_deletions'`: M/=/X/D count as coverage, N excluded (equivalent to `samtools depth -J`)
- `'exclude_deletions'`: only M/=/X count, D and N excluded (equivalent to `samtools depth` default)

**Behavior:**
- NULL rows are ignored
- Empty groups (all NULLs) return NULL
- All rows in a group must have the same `reference_length`
- Mode parameter must be a bind-time constant string (not a column reference)

**Examples:**
```sql
-- Per-position depth for a single reference
SELECT compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')
FROM read_alignments('alignments.bam')
WHERE reference = 'chr1';

-- Mean depth via UNNEST
SELECT AVG(depth) AS mean_depth
FROM (
  SELECT UNNEST(compute_coverage_depth(position, stop_position, cigar, 1000, 'exclude_deletions')) AS depth
  FROM read_alignments('alignments.bam')
  WHERE reference = 'chr1'
);
```

**Performance notes:**
- Memory: allocates `reference_length * 4` bytes per group (e.g., ~1GB for a 250M-position human chromosome)
- Multi-threaded aggregation: each thread maintains independent state, merged at finalization
- Fast path: reads with only M/=/X operations (no N or excluded D) skip CIGAR walking
- Maximum reference_length is 2,000,000,000; use [`compress_intervals`](#merge-overlapping-intervals) for larger references

### Genome coverage

`genome_coverage(alignments, subject_total_length, subject_genome_id)`

Table macro that computes genome coverage from alignment data. It compresses overlapping alignment intervals per reference contig using [`compress_intervals`](#merge-overlapping-intervals), maps contigs to genomes, sums covered bases, and joins with total genome lengths to compute the proportion covered.

**Parameters:**

All three parameters are unquoted table/view names (not string literals):

- `alignments`: A relation with columns `reference` (VARCHAR), `position` (BIGINT), `stop_position` (BIGINT)
- `subject_total_length`: A relation with columns `genome_id` (VARCHAR), `total_length` (BIGINT)
- `subject_genome_id`: A relation with columns `contig_id` (VARCHAR), `genome_id` (VARCHAR)

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `genome_id` | VARCHAR | Genome identifier |
| `covered` | BIGINT | Total number of bases covered |
| `proportion_covered` | DOUBLE | Fraction of genome covered (`covered / total_length`) |

**Behavior:**
- Overlapping alignments on the same contig are merged before counting (via `compress_intervals`)
- Multiple contigs mapping to the same genome have their coverage summed
- Contigs in `alignments` that are not present in `subject_genome_id` are excluded from output
- Uses half-open coordinates consistent with `read_alignments` output

**Examples:**
```sql
CREATE TABLE alignments AS
  SELECT reference, position, stop_position
  FROM read_alignments('alignments.bam');

CREATE TABLE genome_lengths (genome_id VARCHAR, total_length BIGINT);
INSERT INTO genome_lengths VALUES ('genomeA', 100000), ('genomeB', 200000);

CREATE TABLE contig_to_genome (contig_id VARCHAR, genome_id VARCHAR);
INSERT INTO contig_to_genome VALUES
  ('contig1', 'genomeA'), ('contig2', 'genomeA'),
  ('contig3', 'genomeB');

-- Compute genome coverage, filter to genomes with >50% coverage
SELECT * FROM genome_coverage(alignments, genome_lengths, contig_to_genome)
WHERE proportion_covered > 0.5;
```

### Barcode matching

`match_short_barcodes(query_table, ref_table, max_nm:=N, [report_all:=true])`

Hamming-distance matcher for short fixed-length barcodes. Each query sequence is compared to every reference sequence at the same length; pairs within `max_nm` Hamming distance are emitted. The implementation packs each barcode into one or two `uint64_t` (covers up to 32 bp) and uses `__builtin_popcountll` over even/odd nibble masks to count per-nucleotide differences in 1–2 word ops; brute-force O(N×M) is fast enough for the typical input scale (thousands of queries × thousands of references at 8–36 bp).

Intended for UMI binning, sample-index demultiplexing, and any length-uniform barcode lookup. Not appropriate for variable-length adapter matching (use [`extract_linked_amplicon`](utilities.md#extracting-linked-amplicons)) or for similarity search over longer sequences (use [`search_sequences_vsearch`](search.md)).

**Positional parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `query_table` | VARCHAR | Table or view with `(id VARCHAR, sequence VARCHAR)` columns |
| `ref_table` | VARCHAR | Table or view with `(id VARCHAR, sequence VARCHAR)` columns |

**Named parameters:**

| Parameter | Type | Required | Description |
|-----------|------|----------|-------------|
| `max_nm` | INTEGER | yes | Max allowed Hamming distance (≥ 0); pairs above this are not emitted |
| `report_all` | BOOLEAN | no (default `true`) | If `false`, keep only the best hit per query (lowest `nm`; tie-broken by `ref_id` ASC) |

Both input tables must use the column names `id` and `sequence` exactly. `max_nm` is mandatory and must be passed by name (`max_nm := N`); there is no positional form.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `query_id` | VARCHAR | The query row's `id` |
| `ref_id` | VARCHAR | The reference row's `id` |
| `nm` | INTEGER | Hamming distance between the two sequences |

**Behavior:**
- Throws if a query and a candidate reference have different lengths (sequences must be length-uniform within a comparison).
- Empty query or ref tables produce zero rows.
- Empty sequences are rejected at the row level.
- `report_all := false` collapses ties deterministically by `ref_id`.

**Example:**
```sql
-- Bin reads by their 8-bp UMI against a set of canonical UMIs
CREATE TABLE read_umis AS
SELECT read_id AS id, substr(extracted_umi, 1, 8) AS sequence
FROM extracted;

CREATE TABLE umi_refs AS SELECT * FROM (VALUES
    ('binA', 'AACCAACC'),
    ('binB', 'TTGGTTGG')
) t(id, sequence);

-- One best hit per read, allowing 1 mismatch
SELECT query_id AS read_id, ref_id AS bin_id, nm
FROM match_short_barcodes('read_umis', 'umi_refs',
                          max_nm := 1, report_all := false);
```

### MSA column consensus

`compute_msa_consensus(aligned_seq, qual)`

Aggregate function that computes a quality-aware column consensus from a multiple sequence alignment (MSA). Designed as the consensus-building step for long-read amplicon UMI pipelines (Karst-protocol consensus on Revio HiFi), replacing Racon polishing for HiFi-grade per-base quality. General enough to be useful anywhere you have a per-bin MSA + per-read quality.

Algorithm:

1. **Per-column 5-state log-likelihood vote** over `{A, C, G, T, '-'}`. Each base observation uses `p_err = 10^(-q/10)`; gap observations use a fixed `p_err = 0.05`. The argmax base is emitted; columns where gap wins are suppressed from the output.
2. **Posterior Q** derived from `log(sum_{k≠best} exp(ll[k]) / sum_all)`, clamped to `[0, 60]` and emitted as UTINYINT alongside the base.
3. **HP-aware post-correction**: each homopolymer run (length ≥ 2) in the column-voted consensus is replaced by `median(per-read HP length at the corresponding ungapped locus)`. Critical: MAFFT places HP gaps inconsistently and naive column voting biases HP length short.
4. **Bin of size 1** bypasses the MSA pipeline entirely and emits the lone read's gap-stripped sequence + unchanged qual (the MSA is degenerate at n=1).

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `aligned_seq` | VARCHAR | MSA row for one read (gaps encoded as `-`), as produced by [`align_mafft`](alignment_multiple.md) |
| `qual` | LIST(UTINYINT) | The read's original *ungapped* per-base Phred qual (no ASCII offset; matches `read_fastx`) |

The deliberate type asymmetry means callers join the MAFFT row back to the original quality list by `read_id`; Q never enters MAFFT, so it never has to be reassembled across gap insertions.

**Returns:** STRUCT with fields:

| Field | Type | Description |
|-------|------|-------------|
| `seq` | VARCHAR | Consensus sequence (no gaps; HP-corrected) |
| `qual` | LIST(UTINYINT) | Posterior per-base Q list (parallel to `seq`, capped at 60) |

**Behavior:**
- NULL rows are ignored (`IgnoreNull() = true`). Groups where every row is NULL return NULL.
- All rows within a group must share the same MSA width; throws on a width mismatch.
- Each row's `qual` length must equal its sequence's ungapped count; throws on mismatch.
- Tie-break is deterministic: log-likelihood difference < 1e-9 falls back to per-candidate Q-sum, then alphabetical (`A < C < G < T < -`).
- Lower-median is used for even HP-length samples (a length some read actually observed).
- Eager materialization: all `(aligned_seq, qual)` rows are held in the aggregate state. Memory is `O(rows × MSA_width)` per group, which is comfortable for typical UMI bin sizes (5–30 reads × ~5 kb amplicons).

**Example:**
```sql
-- Per-bin consensus via the standard MAFFT → consensus chain.
-- 1. Per-bin MSA
CREATE TABLE mafft_out AS
SELECT * FROM align_mafft('amplicons_with_bin', sample_id := 'bin_id');

-- 2. Join MAFFT's gapped seq back to the original qual list
CREATE VIEW consensus_input AS
SELECT m.bin_id, m.read_id, m.aligned_sequence, a.qual
FROM mafft_out m
JOIN read_amplicons a USING (read_id);

-- 3. Aggregate per bin
SELECT bin_id, (compute_msa_consensus(aligned_sequence, qual)).seq AS consensus
FROM consensus_input
GROUP BY bin_id;
```

**See also:** [`match_short_barcodes`](#barcode-matching), [`compute_pileup`](#per-base-pileup), [`extract_linked_amplicon`](utilities.md#extracting-linked-amplicons), and [`align_mafft`](alignment_multiple.md) — together these primitives compose into a long-read UMI consensus pipeline.
