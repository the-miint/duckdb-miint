# Alignment analysis

A considerable amount of analysis on alignment data can be performed with native SQL. MIINT further provides a set of optimized methods to support specific complex operations.

## Table of contents

- [SAM flag functions](#sam-flag-functions) - Test individual SAM flag bits
- [Per-base pileup](#per-base-pileup) - A per-base CIGAR walker.
- [Alignment slicing](#alignment-slicing) - Slice an alignment based on start and stop coordinates with implicit CIGAR update.
- [Sequence identity](#sequence-identity) - Sequence identity calculation (multi-mode, requires NM/MD for legacy CIGAR)
- [Sequence identity from CIGAR](#sequence-identity-from-cigar) - Sequence identity from extended CIGAR alone
- [Pooled sequence identity](#pooled-sequence-identity) - One identity for a read the aligner split into several records (aggregate)
- [Query length](#cigar-query-length) - Query length from CIGAR
- [Query coverage](#cigar-query-coverage) - Query coverage from CIGAR
- [Query intervals](#cigar-query-intervals) - Which read positions an alignment covers, on the read's own axis
- [Merge overlapping intervals](#merge-overlapping-intervals) - Merge overlapping genomic intervals (aggregate)
- [Genome bins](#genome-bins) - Divide a genome into equal-width bins; map positions and intervals to them
- [Coverage depth](#coverage-depth) - Per-position depth of coverage (aggregate)
- [Genome coverage](#genome-coverage) - Proportion of each genome covered by alignments
- [Per-sample genome coverage](#per-sample-genome-coverage) - The same, reported separately for each sample
- [Region presence](#region-presence) - Three-state (present / absent / not applicable) region calls per sample
- [Region coverage](#region-coverage) - Breadth of coverage of a sub-genome region, with the region's own length as the denominator
- [Circular query coverage](#circular-query-coverage) - Query coverage pooled across the fragments of one read, for reads spanning a circular reference's origin
- [Barcode matching](#barcode-matching) - Hamming-distance matcher for short fixed-length barcodes
- [MSA column consensus](#msa-column-consensus) - Quality-aware consensus from a multiple alignment
- [Concordant-pair identity filtering](#worked-example-joint-identity-filtering-of-concordant-read-pairs) - Filter paired-end alignments as a unit

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
- `alignment_is_primary(flags)` - Primary *line* — neither secondary (0x100) nor supplementary (0x800) set, i.e. `FLAG & 0x900 == 0`. ⚠️ **TRUE for an unmapped read** (see note below).
- `alignment_is_mapped_primary(flags)` - Primary line **and** mapped — `FLAG & 0x904 == 0`; equivalent to `alignment_is_primary(flags) AND NOT alignment_is_unmapped(flags)`. Use this when you mean "a real, mapped primary alignment".
- `alignment_is_qc_failed(flags)` - QC failure (0x200)
- `alignment_is_duplicate(flags)` - PCR/optical duplicate (0x400)
- `alignment_is_supplementary(flags)` - Supplementary alignment (0x800)

**HTSlib-compatible aliases:**
`is_paired`, `is_proper_pair`, `is_unmapped`, `is_munmap`, `is_reverse`, `is_mreverse`, `is_read1`, `is_read2`, `is_secondary`, `is_qcfail`, `is_dup`, `is_supplementary`

> **Note on `alignment_is_primary` and unmapped reads.** `alignment_is_primary` intentionally matches the SAM spec's definition of the *primary line*: the one line per read with neither the SECONDARY (0x100) nor SUPPLEMENTARY (0x800) bit set (`FLAG & 0x900 == 0`). The spec guarantees exactly one such line per read, and that line exists **whether or not the read is mapped** — so `alignment_is_primary(flags)` is `true` for an unmapped read (`flags = 0x4`). This is consistent with the SAM specification and with HTSlib/samtools. If you actually want "a mapped primary alignment", use `alignment_is_mapped_primary(flags)` (or spell out `alignment_is_primary(flags) AND NOT alignment_is_unmapped(flags)`).

**Example:**
```sql
SELECT read_id, flags
FROM read_alignments('alignments.sam')
WHERE alignment_is_mapped_primary(flags);
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

**Output schema:** Same columns as found in the input table (from the recognized alignment column set), with adjusted values for overlapping reads. The identifier columns `read_id`, `reference`, and `mate_reference` keep their input storage type — `VARCHAR`, `BIGINT`, or `UUID` — rather than being coerced to `VARCHAR`, so a `BIGINT`/`UUID` id round-trips through slicing unchanged (consistent with `align_minimap2`).

**Behavior:**
- Reads that don't overlap the region are excluded
- Soft clips (S) do not count as overlap evidence
- Trimmed portions of the CIGAR become hard clips (H)
- Tags (`tag_as` through `tag_sa`) are set to NULL when trimming occurs
- `template_length` is set to NULL when trimming occurs
- `mapq` and mate fields are preserved
- `read_id`, `reference`, and `mate_reference` must each be `VARCHAR`, `BIGINT`, or `UUID` if present; another type is rejected at bind time
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

**See also:** `alignment_seq_identity` for legacy `M`-only CIGAR or BLAST/gap-compressed/gap-excluded identity flavors, and [`cigar_pooled_identity`](#pooled-sequence-identity) when one read is described by several alignment records.

### Pooled sequence identity

Sequence identity for a read the aligner reported as **several alignment records** — a
supplementary alignment, a read spanning a circular reference's origin, any split alignment.
`cigar_sequence_identity` is correct per record but cannot see the others, so no record's
figure describes the molecule. This aggregate pools them.

**Function signature:**

`cigar_pooled_identity(cigar)` — aggregate

**Parameters:**
- `cigar` (VARCHAR): CIGAR strings with extended ops (`=`/`X`), one row per alignment record. Group by whatever identifies one read: `read_id`, or `(read_id, alignment_is_read1(flags))` for paired data.

**Formula:** `sum(match_ops) / sum(alignment_columns)` over the group — the same ratio as
[`cigar_sequence_identity`](#sequence-identity-from-cigar), with both terms summed across the
records first. Clipping contributes to neither, so the split itself costs nothing: a perfect
read reported as two records is 1.0, not 1.0 minus a fabricated gap.

**Returns:** DOUBLE in [0.0, 1.0], or NULL under the same conditions as
`cigar_sequence_identity` applied to the totals:
- every CIGAR in the group is NULL, empty, or `*`
- the group uses only `M` (legacy — can't distinguish matches from mismatches). Re-aligning with extended CIGARs is the only way to recover a *pooled* figure; [`alignment_seq_identity`](#sequence-identity) with NM or MD reconstructs what `M` omits but only per record, and those per-record figures cannot be averaged into the read's identity (see below)
- the group mixes `M` with `=`/`X`, **including across records** — one record in legacy `M` and another extended means the `=`/`X` counts describe only part of the read, so no identity is reported for any of it
- the group has no `=`/`X` ops at all

NULL CIGARs are skipped rather than poisoning the group. On a group of one record the result
is exactly `cigar_sequence_identity` of that record, including every case where that is NULL,
so a query can move from one to the other without changing what identity means.

**Example:**
```sql
-- One identity per read, however many records the aligner emitted
SELECT read_id, cigar_pooled_identity(cigar) AS identity
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags) AND NOT alignment_is_secondary(flags)
GROUP BY read_id;

-- Paired data: keep the mates apart, they are different molecules
SELECT read_id, alignment_is_read1(flags) AS is_read1,
       cigar_pooled_identity(cigar) AS identity
FROM read_alignments('alignments.bam')
GROUP BY read_id, is_read1;
```

A split read's records mislead individually, and averaging them is not the answer either. Take
a read split into a 1000-column record at 100% and a 100-column record with one mismatch
(`1000=` and `99=1X`): the records report 1.0 and 0.99, their mean is 0.995, and the read is
**0.99909** — one mismatch in 1100 aligned columns. The mean overstates the error more than
fivefold, because it gives a record ten times shorter equal weight. Pooling weights each
record by the columns it actually aligned, which is the only reading under which the figure
means "identity of this read against this reference".

> Do **not** approximate this with `cigar_sequence_identity(string_agg(cigar, ''))`.
> Concatenating CIGARs builds a string the SAM spec forbids, because clipping is only legal at
> an end, and it fails outright on a group containing an unmapped record. It also does not
> generalise: anything adjacency-sensitive — gap-compressed identity above all — is *wrong* on
> concatenated input, because joining two records manufactures a gap-open event that never
> happened.

**See also:** [`circular_query_coverage`](#circular-query-coverage), which reports this as its
`identity` column alongside the coverage and topology evidence for the same read.

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

> ⚠️ **This is a per-row measure and cannot see a read's other alignment records.** When a
> read is split across several records — it spans the origin of a circular reference, or has
> supplementary alignments — each record covers only part of the read, so **no row passes a
> coverage floor like `> 0.9` even when the read is fully explained**. That is a silent loss,
> and for assembled circular genomes it is concentrated at the origin. Use
> [`circular_query_coverage`](#circular-query-coverage) to pool a read's records; it returns
> the same value as this function for reads with a single alignment.

### CIGAR query intervals

Return which positions of the read an alignment covers, as half-open `[start, stop)`
intervals on the **read's own axis**.

`cigar_query_coverage` answers "how much of the read does this one alignment record
explain" and collapses the answer to a single number. When a read is split across several
records — an origin-spanning read on a circular reference, a supplementary alignment, a
chimera — you need to know *which* part of the read each record explains before you can
combine them. That is what this returns.

**Function signature:**

`cigar_query_intervals(cigar, flags, [type='aligned'])`

**Parameters:**
- `cigar` (VARCHAR): CIGAR string from the alignment
- `flags` (USMALLINT): SAM flags. Only the reverse-strand bit (`0x10`) is consulted
- `type` (VARCHAR, default `'aligned'`): `'aligned'` counts `M`/`=`/`X`; `'mapped'` also counts `I`. Same vocabulary as [`cigar_query_coverage`](#cigar-query-coverage)

**Returns:** `LIST<STRUCT(start BIGINT, stop BIGINT)>` — 0-based, half-open, ascending and
non-overlapping. Empty list for an unmapped or empty CIGAR. NULL if `cigar` or `flags` is
NULL.

**Behavior:**
- **Intervals are in the original read's orientation, not the reference's.** SAM writes
  CIGARs in reference orientation, so on a reverse-strand record the leading clip sits at
  the read's 3′ end. The reverse bit mirrors each interval onto the read axis, which is what
  makes intervals from different records of the same read directly comparable.
- `S` and `H` advance the cursor but are never covered, so the denominator is the full read
  length — the same quantity `cigar_query_length(cigar, true)` reports.
- `D`, `N` and `P` consume no query, so the runs they separate are contiguous on the read
  and are returned merged.
- Under `'aligned'`, an insertion is a gap in coverage and splits the interval; under
  `'mapped'` it does not.
- The field names `start`/`stop` match [`compress_intervals`](#merge-overlapping-intervals),
  so the two compose directly.

**Examples:**
```sql
-- A trailing soft clip is not covered
SELECT cigar_query_intervals('3000=3000S', 0);
-- [{'start': 0, 'stop': 3000}]

-- Same alignment on the reverse strand: the covered block is at the other end of the read
SELECT cigar_query_intervals('3000=3000S', 16);
-- [{'start': 3000, 'stop': 6000}]

-- An insertion splits the covered region under 'aligned' but not under 'mapped'
SELECT cigar_query_intervals('100=10I100=', 0);            -- two intervals
SELECT cigar_query_intervals('100=10I100=', 0, 'mapped');  -- one interval

-- Pool the fragments of each read onto its own axis
SELECT read_id, compress_intervals(iv.start, iv.stop) AS covered
FROM (
  SELECT read_id, UNNEST(cigar_query_intervals(cigar, flags)) AS iv
  FROM read_alignments('alignments.bam')
  WHERE NOT alignment_is_unmapped(flags) AND NOT alignment_is_secondary(flags)
)
GROUP BY read_id;
```

For the common case of pooling coverage across a read's fragments, prefer
[`circular_query_coverage`](#circular-query-coverage), which does this and reports the
evidence needed to tell a genuine wrap from a chimera.

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

### Genome bins

`bin_of(pos, n_bins, genome_length)`
`bin_start(b, n_bins, genome_length)`
`interval_bins(start, stop, n_bins, genome_length)`

Scalar macros that divide a genome into `n_bins` near-equal-width bins. `bin_of` maps a position to its bin; `bin_start` maps a bin back to its first position; `interval_bins` reports which bins a half-open interval touches, designed to be `UNNEST`ed so a read spanning a boundary is counted in every bin it overlaps.

Together these support finding genomic regions whose coverage differs between sample groups: bin the alignments, count samples or reads per bin per group, and rank bins by the standard deviation of the per-group counts. Only the fan-out needed a primitive — the ranking half is plain SQL.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `pos` | BIGINT | A 1-based genomic position, in `[1, genome_length]` |
| `start` | BIGINT | 1-based inclusive interval start, matching `read_alignments` |
| `stop` | BIGINT | Exclusive interval end (half-open) |
| `b` | BIGINT | 0-based bin index, in `[0, n_bins]` |
| `n_bins` | BIGINT | Number of bins (must be ≥ 1); widths differ by at most one base |
| `genome_length` | BIGINT | Total length being binned (must be ≥ 1) |

**Returns:** `interval_bins` → `LIST(BIGINT)`, the bin indices touched, ascending and deduplicated. `bin_of` and `bin_start` → `BIGINT`. Every operand is normalized to `BIGINT` internally — in the bounds checks as well as the arithmetic — so `INTEGER` inputs, which is what `read_gff` and `read_ncbi_annotation` emit, cannot overflow. A `genome_length` or `n_bins` arriving as `DOUBLE`/`DECIMAL` (an inferred CSV or Parquet column) is rounded once and behaves exactly like that rounded integer.

**Bin index convention:** bins are **0-based**. They are array positions rather than genomic coordinates, so this deliberately cuts against the 1-based coordinate convention used elsewhere. The bin index is the join key you will carry into downstream queries, so it is worth stating explicitly in your own code too.

`bin_of` is the definition and `bin_start` is its left inverse:

```
bin_of(p)    = ((p - 1) * n_bins) // genome_length
bin_start(b) = ceil(b * genome_length / n_bins) + 1
```

The guarantee is **containment**, `bin_start(bin_of(p)) <= p < bin_start(bin_of(p) + 1)` for every `p` in `[1, genome_length]`. It is not a round trip in both directions: when `n_bins > genome_length`, zero-width bins make `bin_start` non-injective, so `bin_of(bin_start(b))` need not be `b` (with `n_bins=5, genome_length=3`, `bin_start(2) = 3` but `bin_of(3) = 3`).

**Bin edges are half-open, like every other interval in miint.** Bin `b` spans `[bin_start(b), bin_start(b + 1))` and its width is the plain difference — no `+1`, matching the `stop_position` rule in [`docs/internals/architecture.md`](internals/architecture.md). That makes `bin_start(0) = 1` and `bin_start(n_bins) = genome_length + 1`, so a bin drops straight into anything taking a half-open `(start, stop)` pair — `compress_intervals`, or your own region table — with no adjustment. Bin widths differ by at most one base, and every base belongs to exactly one bin.

> **If your source uses an inclusive end**, add one before calling: `interval_bins(start, end + 1, ...)`. Getting this wrong is silent — every interval's last base lands one bin early and the ranking still looks plausible. `read_alignments`, `compress_intervals`, `read_gff` and `read_ncbi_annotation` are all already half-open, so no adjustment is needed for those.

> **Do not pass `stop_position` to `bin_of`.** `bin_of` takes a position; `stop_position` is an exclusive end, so a read reaching the reference end carries `genome_length + 1` and `bin_of` raises. Use `bin_of(stop_position - 1, ...)` for the last covered base, or `interval_bins(position, stop_position, ...)` for the whole read. The error message says as much, but it aborts the query rather than one row — and it only triggers on reads at the very end of a reference, so it will pass on a small fixture and fail on a real BAM.

**Behavior:**
- NULL in any argument → NULL
- `stop <= start` (empty or inverted interval) → empty list, so `UNNEST` drops the row
- Intervals are clamped to `[1, genome_length]`; bins outside the genome are never emitted
- When `n_bins` exceeds `genome_length` some bins are zero-width. This is legitimate when a fixed bin count is applied across genomes of very different sizes. Such bins contain no bases and are **not** emitted, even for an interval that spans across them.
- `n_bins < 1` or `genome_length < 1` **raises**, rather than returning an empty list, so a bin count or length that computed to zero surfaces instead of quietly emptying the result.
- A NULL `genome_length` yields NULL, **not** an error — and `UNNEST(NULL)` emits no rows, so those reads disappear with no diagnostic. This is the shape a missed join takes (a miss produces NULL, not 0), so the positivity check above does not cover it. If a genome may be absent from your lengths table, guard the join yourself rather than relying on an error.
- `bin_of` raises on a position outside `[1, genome_length]`; `bin_start` raises on a bin index outside `[0, n_bins]`. Both would otherwise return a plausible-looking but fictional coordinate.

**Examples:**
```sql
-- Rank bins by how variable their per-group sample counts are, then turn the
-- top-ranked bins into half-open regions. One statement, so top_bins stays in scope.
CREATE TABLE differential_regions AS
WITH binned AS (
    SELECT p.genome_id, md.country, p.sample_id,
           UNNEST(interval_bins(p.start, p.stop, 500, g.length)) AS bin_index
    FROM positions p
    JOIN sample_metadata md USING (sample_id)
    JOIN genome_lengths  g  USING (genome_id)
    WHERE p.genome_id = 'G000436435'
), per_group AS (
    SELECT genome_id, bin_index, country, COUNT(DISTINCT sample_id) AS sample_hits
    FROM binned GROUP BY genome_id, bin_index, country
), top_bins AS (
    SELECT genome_id, bin_index, stddev_samp(sample_hits) AS sample_hits_std
    FROM per_group GROUP BY genome_id, bin_index
    ORDER BY sample_hits_std DESC NULLS LAST, bin_index
    LIMIT 20
)
SELECT genome_id,
       bin_start(bin_index, 500, 4719737)     AS region_start,
       bin_start(bin_index + 1, 500, 4719737) AS region_stop,   -- exclusive
       'bin_' || bin_index                    AS region_id,
       sample_hits_std
FROM top_bins;

-- Which bin does a single position fall in?
SELECT bin_of(481323, 500, 4719737);
```

> The `ORDER BY` above carries `bin_index` as a tiebreak. Without it, ties in `sample_hits_std` — common, since many bins share a std dev of 0 — make the reported top-20 vary between runs and thread counts.

> **A group with no sample in a bin contributes no row**, so `stddev_samp` sees one value and returns NULL rather than treating the absent group as zero. If absent-means-zero is what you want, densify against your sample roster first — that is a metadata decision, not a binning one.

**Performance notes:**
- When `n_bins <= genome_length` (the normal case) every bin is at least one base wide, so the result is a contiguous bin range and the cost is the number of bins the interval touches
- When `n_bins > genome_length` the macro enumerates positions rather than bins, so the cost is the interval's length rather than `n_bins`. Binning a 3 bp plasmid into 10⁶ bins costs three positions, not a million

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
- `subject_total_length`: A relation with columns `genome_id`, `total_length` (BIGINT)
- `subject_genome_id`: A relation with columns `contig_id` (VARCHAR), `genome_id`

`genome_id` is passed through without casting, so VARCHAR, BIGINT and UUID identifiers all work, and the output column keeps the type given in `subject_genome_id`.

The two reference relations do not have to declare `genome_id` with the *same* type — it is a join key, so DuckDB applies its usual implicit casts, and e.g. a BIGINT `77` matches a VARCHAR `'77'`. Prefer matching types anyway: when a value cannot be converted you get a cast error naming the column (`Could not convert string 'genome_A' to INT64 …`) rather than a clear diagnosis of the mismatch.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `genome_id` | as supplied | Genome identifier, in the type it was given |
| `covered` | BIGINT | Total number of bases covered |
| `proportion_covered` | DOUBLE | Fraction of genome covered (`covered / total_length`) |

**Behavior:**
- Overlapping alignments on the same contig are merged before counting (via `compress_intervals`)
- Multiple contigs mapping to the same genome have their coverage summed
- Contigs in `alignments` that are not present in `subject_genome_id` are excluded from output
- `subject_genome_id` is deduplicated on `(contig_id, genome_id)` before use, so a repeated mapping row cannot double-count coverage. A contig mapped to two *different* genomes is not a duplicate and still counts toward both
- Uses half-open coordinates consistent with `read_alignments` output
- **Multi-sample input is pooled.** This macro has no sample dimension: given alignments from several samples it returns one row per genome, computed as though every read came from a single sample. Use [`genome_coverage_per_sample`](#per-sample-genome-coverage) instead

**Errors:**

A genome that has coverage but no usable denominator raises, rather than being dropped or divided into a meaningless number:

- `no total_length entry for genome_id '<X>'` — a genome mapped in `subject_genome_id` and covered by alignments has no row in `subject_total_length`. Previously such genomes were silently omitted along with their covered bases
- `total_length must be positive for genome_id '<X>' (got <v>)` — the length is zero, negative, or NULL. A zero length previously yielded `proportion_covered = inf`

A genome listed in `subject_total_length` that has *no* coverage is not an error; it simply produces no row. Passing a complete reference catalogue is therefore fine.

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

### Per-sample genome coverage

`genome_coverage_per_sample(alignments, subject_total_length, subject_genome_id)`

Table macro that computes genome coverage **separately for each sample**. It is the sibling of [`genome_coverage`](#genome-coverage) and performs identical arithmetic, but partitions every step by a `sample_id` column that `alignments` must additionally carry: intervals are compressed per `(sample_id, reference)` and covered bases summed per `(sample_id, genome_id)`.

Use this whenever your alignment relation holds more than one sample. `genome_coverage` has no sample dimension and will pool them into a single row per genome — two samples covering 22 and 11 bases of a 100 bp genome return one row at `0.22`, over-reporting the second sample by 2×.

**Parameters:**

All three parameters are unquoted table/view names (not string literals):

- `alignments`: A relation with columns `sample_id`, `reference` (VARCHAR), `position` (BIGINT), `stop_position` (BIGINT)
- `subject_total_length`: A relation with columns `genome_id`, `total_length` (BIGINT)
- `subject_genome_id`: A relation with columns `contig_id` (VARCHAR), `genome_id`

`sample_id` and `genome_id` are passed through without casting, so VARCHAR, BIGINT and UUID identifiers all work.

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `sample_id` | as supplied | Sample identifier, in the type it was given |
| `genome_id` | as supplied | Genome identifier, in the type it was given |
| `covered` | BIGINT | Total number of bases covered in that sample |
| `proportion_covered` | DOUBLE | Fraction of genome covered (`covered / total_length`) |

**Behavior:**
- Intervals are merged *within* a sample, never *across* samples
- A `(sample, genome)` pair with no qualifying alignments produces no row, rather than a zero-coverage row
- On single-sample input the result matches `genome_coverage` exactly
- All other behavior — contig merging, deduplication of `subject_genome_id`, half-open coordinates — matches `genome_coverage`

**Errors:**

Raises on the same unusable-denominator cases as `genome_coverage`, plus:

- `NULL values in sample_id column 'sample_id'` — a NULL sample cannot form a meaningful partition, and is rejected rather than reported as its own group. This matches every other per-sample function in the extension

**Examples:**
```sql
CREATE TABLE alignments AS
  SELECT 'sampleA' AS sample_id, reference, position, stop_position
  FROM read_alignments('sampleA.bam')
  UNION ALL
  SELECT 'sampleB', reference, position, stop_position
  FROM read_alignments('sampleB.bam');

SELECT * FROM genome_coverage_per_sample(alignments, genome_lengths, contig_to_genome)
ORDER BY sample_id, genome_id;
```

Feeding per-sample coverage into [absolute quantification](absolute_quantification.md#absquant_cell_counts), which expects a `(sample_id, feature_id, coverage)` shape:

```sql
SELECT sample_id, genome_id AS feature_id, proportion_covered AS coverage
FROM genome_coverage_per_sample(alignments, genome_lengths, contig_to_genome);
```

**If you need the pooled number per sample.** Some workflows want each sample compared against coverage computed across the whole study rather than within the sample. There is no parameter for this — but it is one join, because `genome_coverage` already produces exactly that number:

```sql
-- Broadcast the pooled (all-sample) coverage to every sample.
SELECT s.sample_id, g.genome_id, g.covered, g.proportion_covered
FROM genome_coverage(alignments, genome_lengths, contig_to_genome) g
CROSS JOIN (SELECT DISTINCT sample_id FROM alignments) s;
```

Be deliberate about which one you want: the pooled value is identical for every sample and is *not* a property of any single sample.
### Region presence

`region_presence(positions, regions, samples)`

Table macro answering, per sample, whether a genomic region is present. The answer is **three-state**, and the third state is what makes downstream statistics honest:

| State | Meaning |
|---|---|
| `present` | The sample has at least one covered interval overlapping the region. |
| `absent` | The sample has coverage of the genome, but none of it overlaps the region. |
| `not applicable` | The sample has no coverage of the genome at all. |

Collapsing `not applicable` into `absent` conflates "the organism is here and this region is missing from it" — a strain-content claim — with "the organism is not here", a detection failure. Those must not be pooled, and the distinction is easy to lose when hand-rolling the SQL.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `positions` | relation | Columns `sample_id` (any type), `genome_id` (any type), `start` (BIGINT), `stop` (BIGINT). Intervals 1-based half-open, matching [`compress_intervals`](#merge-overlapping-intervals) output. |
| `regions` | relation | Columns `genome_id` (any type), `region_start` (BIGINT), `region_stop` (BIGINT), `region_id` (any type). Half-open on the same convention. |
| `samples` | relation | A `sample_id` column: the full cohort roster. |

> **These are not `read_alignments`' column names.** `positions` is a *covered-interval* relation, so rename on the way in — `reference` → `genome_id`, `position` → `start`, `stop_position` → `stop`. The first example below shows the rename. (`genome_coverage` takes the raw alignment names instead, because it consumes alignments directly.)

**Returns:** one row per region row per sample.

| Column | Source |
|---|---|
| `sample_id` | from `samples`, uncast |
| `genome_id` | from **`regions`**, uncast |
| `region_id` | from `regions`, uncast |
| `region_start`, `region_stop` | from `regions` |
| `state` | VARCHAR — `present` / `absent` / `not applicable` |

`sample_id`, `genome_id` and `region_id` are carried through uncast, so VARCHAR, BIGINT and UUID identifiers all survive end to end. Coordinates are normalized to `BIGINT`, so they come back as `BIGINT` whatever they went in as — including VARCHAR, which is what a BED file or an inferred CSV column gives you.

> **Match your identifier types across relations.** `genome_id` joins `positions` to `regions`, and `sample_id` joins `positions` to `samples`. DuckDB implicit-casts, so a BIGINT `77` does match a VARCHAR `'77'` — but two things follow. The output `genome_id` takes the **`regions`** type, not the `positions` type, so a carefully typed coverage table and a loosely typed region table give you the loose one back. And a single unconvertible value aborts the whole query with a bare `Conversion Error: Could not convert string 'abc' to INT64 … source column sample_id` that names neither `region_presence` nor which relation the bad value came from. One malformed roster id is enough.

The coordinates are emitted, not just `region_id`, because the presence join is keyed on them. Two rows sharing a `region_id` but covering different spans are scored independently, and without the coordinates in the output those rows would be indistinguishable — a `PIVOT` would silently collapse them and a join back to metadata would duplicate samples. When `region_id` is unique, the normal case, `(sample_id, region_id)` still keys the output.

**Why `samples` is required.** `not applicable` cannot be derived from `positions` alone: a sample with no coverage contributes no rows, so it is indistinguishable from a sample that is not in the study. Without the roster this function could only emit `present` and `absent` — the very conflation it exists to prevent. Keeping the roster explicit also means the caller decides what "the cohort" is.

**Behavior:**
- The overlap test is `start < region_stop AND stop > region_start`. A read starting exactly where a region ends shares no base with it and is `absent`.
- `positions` does **not** need to be pre-compressed. Presence asks only whether *any* interval overlaps, so overlapping input intervals give the same answer as their union — no `compress_intervals` pass needed.
- Samples appearing in `positions` but not in the roster are ignored — the roster defines the cohort.
- Duplicate rows in `positions`, `regions` or `samples` do not multiply the output — which is what makes passing uncompressed alignments safe.
- **NULL and empty intervals in `positions`.** A NULL `start` or `stop` has a defined meaning — a sample that exists with zero coverage — so the row is dropped and does *not* count as covering the genome. An interval with `start >= stop` covers no bases and is dropped for the same reason. Neither promotes a sample from `not applicable` to `absent`. Note the macro takes no `genome_length`, so it can recognise an *empty* interval but not an out-of-range one: `[-5, -1)` is well formed and does count as coverage. Clamp upstream if that matters.
- **A NULL `sample_id` or `genome_id` in `positions` raises.** Unlike the interval columns, these have no defined NULL meaning; silently dropping an unattributable interval would downgrade real coverage to a non-detection.
- Raises on a NULL `sample_id` in the roster, on any NULL column in `regions`, and on an empty or inverted region (`region_stop <= region_start`) — each would otherwise report every sample `absent`, indistinguishable from a real negative result. Both `regions` errors name the offending row.
- An empty roster, an empty `regions`, or empty `positions` are all legitimate: the first two produce no rows, and the third makes every pair `not applicable`.

> **Validation covers the rows the query reads.** The guards are ordinary SQL expressions, so a filter such as `WHERE region_id = 'PC351'` may be pushed beneath them and a malformed row *outside* that filter will not raise. This never changes the rows you get back — each output row depends only on its own region and sample — but do not treat a clean run of a filtered query as a validation pass over the whole relation.

**Examples:**
```sql
-- Build `positions` from alignments: compress per (sample, contig), then rename
-- into the column contract above.
CREATE TABLE positions AS
SELECT sample_id, reference AS genome_id, ci.start, ci.stop
FROM (
    SELECT sample_id, reference, UNNEST(compress_intervals(position, stop_position)) AS ci
    FROM alignments
    GROUP BY sample_id, reference
);

CREATE TABLE regions AS SELECT * FROM (VALUES
    ('G000436435', 481323, 486671, 'PC351')
) t(genome_id, region_start, region_stop, region_id);

CREATE VIEW roster AS SELECT sample_id FROM sample_metadata;

-- Three-state calls per sample
SELECT * FROM region_presence(positions, regions, roster);

-- How many samples carry the region, keeping non-detections separate
SELECT state, COUNT(*) FROM region_presence(positions, regions, roster)
GROUP BY state;

-- Feed presence into PERMANOVA as a metadata variable. Non-detections are
-- DROPPED rather than pooled with 'absent'.
CREATE TABLE region_md AS
    SELECT sample_id, state AS pc351
    FROM region_presence(positions, regions, roster)
    WHERE region_id = 'PC351' AND state IN ('present', 'absent');

SELECT * FROM permanova('dm', 'region_md', variables := ['pc351'],
                        n_permutations := 999999, seed := 42);
```

> **Filter positively — `state IN ('present', 'absent')`, not `state <> 'not applicable'`.** The state labels are plain strings and DuckDB will happily compare against one that does not exist, so a typo cannot be caught for you. It fails differently in the two forms: a typo in `<>` matches *everything*, silently pooling non-detections back into the cohort, whereas a typo in `IN` drops the mistyped level — an empty or visibly halved cohort rather than a quietly wrong one. Neither is safe if you don't check the row count, but only one of them fails loudly enough to notice.

Long form is deliberate — it is what `read_biom`, `woltka_ogu` and the diversity functions already consume. For a sample x region matrix, `PIVOT`:

```sql
PIVOT (SELECT sample_id, region_id, state FROM region_presence(positions, regions, roster))
ON region_id USING first(state);
```

> `first(state)` picks one row per `(sample_id, region_id)` cell, so pivot only when `region_id` is unique across your `regions` relation. If it is not, pivot on the coordinates too.

Regions can come straight from [`bin_start`](#genome-bins) — `(bin_start(b), bin_start(b + 1))` is already half-open on the same convention:

```sql
CREATE TABLE bin_regions AS
SELECT genome_id,
       bin_start(bin_index, 500, length)     AS region_start,
       bin_start(bin_index + 1, 500, length) AS region_stop,
       genome_id || ':bin_' || bin_index     AS region_id
FROM top_bins JOIN genome_lengths USING (genome_id)
WHERE bin_start(bin_index + 1, 500, length) > bin_start(bin_index, 500, length);
```

> `region_id` is qualified with `genome_id` on purpose. A bare `'bin_' || bin_index` repeats across genomes — every genome has a bin 3 — which breaks the uniqueness the `PIVOT` recipe above depends on.

> That `WHERE` matters only when `n_bins > genome_length`, where `bin_start` legitimately produces zero-width bins. A zero-width region is unanswerable — every sample would score `absent` — so `region_presence` rejects it rather than inventing a negative result, and the error names this filter as the remedy.

### Region coverage

`region_coverage(positions, regions)`

Table macro giving, per sample, how much of a sub-genome region is covered — with the **region's** length as the denominator. [`genome_coverage`](#genome-coverage) divides by the genome, so a fully covered 5 kb region inside a 5 Mb genome reports 0.1% there and 100% here. That denominator is the point of the function.

Two things have to be right for the number to mean anything, and the macro does both:

- Alignments are **clipped** to the region before they are merged, so one hanging over the boundary contributes only its in-region portion and one strictly containing the region contributes exactly `region_length`, not its own length.
- The clipped intervals are then **merged**, so a base under two alignments is counted once. An unmerged sum inflates the proportion — past 1.0 on a short region, invisibly on a long one.

**Parameters:**

| Parameter | Type | Description |
|-----------|------|-------------|
| `positions` | relation | Columns `sample_id` (any type), `genome_id` (any type), `start`, `stop` (any integer type). 1-based half-open. |
| `regions` | relation | Columns `genome_id` (any type), `region_start`, `region_stop` (any integer type), `region_id` (any type). Half-open on the same convention. |

Coordinates are normalized to `BIGINT` internally and come back as `BIGINT`, so `INTEGER` inputs — what `read_gff` and `read_ncbi_annotation` emit — need no cast, and are widened before they can overflow. Numeric coordinates arriving as `VARCHAR`, which is what an unconfigured BED or CSV read gives you, are also normalized rather than compared as text; that is defensive, not an invitation — type your coordinate columns.

> **Same column contract as [`region_presence`](#region-presence)**, deliberately — the two are called on the same relation in the same pipeline, and a rename between two consecutive lines of a query would be a needless trap. It is *not* `genome_coverage`'s contract: coming from `read_alignments` you rename once, in a view (`reference` → `genome_id`, `position` → `start`, `stop_position` → `stop`), and both region functions then take it. `genome_coverage` keeps the alignment names because it predates the sample dimension.

`positions` does not need to be pre-compressed — `compress_intervals` runs inside — so raw alignment rows are safe input.

**Returns:** one row per `(sample, region)` pair **that has overlap**. See [No zero rows](#no-zero-rows) below.

| Column | Type | Description |
|---|---|---|
| `sample_id` | from `positions`, uncast | |
| `genome_id` | from **`regions`**, uncast | |
| `region_id` | from `regions`, uncast | |
| `region_start`, `region_stop` | BIGINT | The region's bounds, as supplied |
| `covered` | BIGINT | Covered bases **within the region** |
| `region_length` | BIGINT | `region_stop - region_start` (half-open, no `+1`) |
| `proportion_covered` | DOUBLE | `covered / region_length` |

Identifiers are carried through uncast, so VARCHAR, BIGINT and UUID all survive end to end. The same [type-matching caveat as `region_presence`](#region-presence) applies — `genome_id` takes the **`regions`** type in the output, and one unconvertible identifier aborts the query with an error that names neither the function nor the relation.

The coordinates are emitted alongside `region_id` because the join is keyed on them: two rows sharing a `region_id` but covering different spans are scored independently, and without the coordinates those rows are indistinguishable downstream.

<a name="no-zero-rows"></a>
**No zero rows, on purpose.** A `(sample, region)` pair with no overlap gets no row — not a row with `covered = 0`. A fabricated zero asserts the region *was measured and found empty*, which is true of a sample that covers the genome elsewhere and false of a sample where the organism simply is not present. That is the conflation [`region_presence`](#region-presence) exists to prevent, and it cannot be undone downstream. To fill the grid honestly, compose the two — `region_presence` supplies the roster and the three states, `region_coverage` supplies the number where there is one:

```sql
SELECT p.sample_id, p.state,
       COALESCE(c.covered, 0)             AS covered,
       COALESCE(c.proportion_covered, 0)  AS proportion_covered
FROM region_presence(positions, regions, roster) p
LEFT JOIN region_coverage(positions, regions) c
       ON c.sample_id = p.sample_id AND c.genome_id = p.genome_id
      AND c.region_id = p.region_id
      AND c.region_start = p.region_start AND c.region_stop = p.region_stop
WHERE p.state IN ('present', 'absent');   -- drop non-detections; do not pool them
```

> **Join on the coordinates too, not just `region_id`.** Both functions emit `region_start`/`region_stop` because that is the tuple identifying a region, and two spans sharing a `region_id` are scored independently by design. Keyed on `region_id` alone this join cross-assigns: with regions `('g1',50,60,'RD')` and `('g1',700,750,'RD')` it reports the `[700,750)` row as `absent` while handing it `covered = 5` — the value belonging to `[50,60)` — and duplicates every row of a sample present in both spans. An `absent` row carrying coverage is a contradiction with no downstream tell.

**Behavior:**
- The overlap test is `start < region_stop AND stop > region_start`. An alignment starting exactly where a region ends shares no base with it and contributes nothing.
- **Overlapping regions are each scored against their own bounds**, so a base can count in two regions. Well defined, and required for sliding-window regions — but it means `SUM(covered)` across regions is not a partition of the genome.
- **Multi-contig genomes.** `regions.genome_id` is matched against `positions.genome_id` directly, so regions are effectively contig-level. Whole-genome coordinates for a multi-contig assembly are ambiguous; pick the contig yourself.
- Duplicate rows in either relation do not inflate the result — the merge and the final grouping absorb them.
- **NULL and degenerate intervals in `positions` drop; a NULL `sample_id` or `genome_id` raises.** The interval columns have a defined NULL meaning (a sample that exists with no coverage); the identity columns do not, and silently dropping an unattributable interval would subtract real coverage from whichever sample it belonged to.
- **An unusable region raises** — any NULL column, a zero-width region (no denominator), or inverted coordinates. Each error names the offending `region_id` and its bounds, and zero-width and inverted are reported separately because the remedies differ: filtering is right for a bin expansion and wrong for transposed coordinates.
- Empty `positions` or empty `regions` are legitimate and produce no rows.

> **Validation covers the rows the query reads.** As with `region_presence`, the guards are ordinary SQL expressions, so a filter such as `WHERE region_id = 'PC351'` may be pushed beneath them and a malformed row *outside* that filter will not raise. This never changes the rows you get back.

**Examples:**
```sql
-- Rename read_alignments' columns into the contract once, in a view.
CREATE VIEW positions AS
SELECT sample_id, reference AS genome_id, position AS start, stop_position AS stop
FROM alignments;

CREATE TABLE regions AS SELECT * FROM (VALUES
    ('G000436435', 481323, 486671, 'PC351')
) t(genome_id, region_start, region_stop, region_id);

-- Per-sample breadth of one differential region, region-relative denominator
SELECT * FROM region_coverage(positions, regions)
ORDER BY proportion_covered DESC;

-- Rank regions by how consistently the cohort covers them
SELECT region_id, AVG(proportion_covered) AS mean_prop, COUNT(*) AS n_samples
FROM region_coverage(positions, regions)
GROUP BY region_id
ORDER BY mean_prop DESC;
```

**When each genome is a single contig**, a whole-genome "region" reproduces `genome_coverage` with a sample dimension, which makes a useful equivalence check:

```sql
CREATE TABLE whole AS
SELECT genome_id, 1 AS region_start, total_length + 1 AS region_stop, genome_id AS region_id
FROM genome_lengths;

SELECT * FROM region_coverage(positions, whole);
```

> **That equivalence does not survive a multi-contig genome**, and neither does the obvious workaround. `genome_coverage` maps contigs to genomes through `subject_genome_id` and sums across them; `region_coverage` matches `regions.genome_id` against `positions.genome_id` directly, so its regions are contig-level. Hand it whole-genome coordinates for a 3-contig assembly and it matches nothing — **zero rows**, where `genome_coverage` returns the real proportion.
>
> Relabelling the contigs to the genome name is worse, not better: each contig has its own 1-based frame, so three `[10, 60)` intervals stack onto one another, the merge collapses them, and you get `covered = 50` where the truth is `150`. Zero rows are at least visibly wrong; 50 looks like an answer.
>
> Keep regions **contig-level** — one region per contig, in that contig's own frame — and sum afterwards:
>
> ```sql
> SELECT SUM(covered) AS covered,
>        SUM(covered)::DOUBLE / SUM(region_length) AS proportion_covered
> FROM region_coverage(positions, contig_regions)
> GROUP BY sample_id;
> ```

Regions can come straight from [`bin_start`](#genome-bins) — `(bin_start(b), bin_start(b + 1))` is already half-open on the same convention. Use the [`bin_regions` recipe](#region-presence) shown for `region_presence`; both functions take it unchanged.

### Circular query coverage

`circular_query_coverage(alignments, reference_lengths, [coverage_type='aligned'])`

> `reference_lengths` needs an `is_circular` column. Nothing in an alignment records whether
> a reference is circular, and the wrap test is only meaningful if it is — so the macro asks
> rather than assumes.

How much of each read is explained by each reference, **pooling every alignment record the
aligner split the read into**.

#### Why this exists alongside `cigar_query_coverage`

[`cigar_query_coverage`](#cigar-query-coverage) answers a per-row question — *how much of the
read does **this record** explain* — and cannot see the read's other records. A read that
spans the **origin of a circular reference** held as a linearised contig is emitted as two or
more records, so no single row reports much beyond half the read:

| record | position–stop | cigar | `cigar_query_coverage` |
|---|---|---|---|
| primary | 27001–30001 | `3000=3000S` | 0.5 |
| supplementary | 1–3001 | `3000H3000=` | 0.5 |

The read is a perfect, unambiguous match to the reference it came from, yet a query-coverage
floor of 0.90 discards it. For an assembled circular genome that is a systematic, silent loss
localised at the origin — and it is worst for exactly the elements most often recovered as
complete circles, small plasmids and phages.

This macro answers the question that actually matters — *how much of the read does **this
reference** explain* — and returns 1.0 for the read above.

**On a read with a single alignment the two agree exactly**, so you can use this
unconditionally rather than branching on whether a reference is circular.

Note that coverage is the **union** of the fragments' footprints on the read, not their sum.
Summing overshoots: a junction read with a couple of bases deleted across the origin produces
fragments that overlap by a base or two, and their coverages sum to more than 1.0. The union
is bounded by 1.0, and is unaffected by a read wrapping the reference more than once.

**Parameters:**

The first two are unquoted table/view names, not string literals:

- `alignments`: a relation with columns `read_id`, `flags`, `reference`, `position` (BIGINT), `stop_position` (BIGINT), `cigar` (VARCHAR) — the schema produced by `read_alignments` and by the `align_*` functions. `flags` may be any integer type; it is cast internally
- `reference_lengths`: a relation with columns `reference`, `length` (BIGINT) and `is_circular` (BOOLEAN). [`read_alignment_header`](reading.md#read_alignment_headerpath) supplies the first two; `is_circular` you add, because circularity is not recorded anywhere in an alignment and cannot be inferred from one:
  ```sql
  -- whole assembly is circular
  SELECT reference, length, true AS is_circular FROM read_alignment_header('aln.bam')
  -- mixed assembly: circular chromosome, linear contigs
  SELECT reference, length, reference IN ('chr', 'plasmid_1') AS is_circular
  FROM read_alignment_header('aln.bam')
  ```
  A reference whose `is_circular` is NULL is rejected rather than assumed circular — see below
- `coverage_type` (VARCHAR, named argument, default `'aligned'`): what counts as covered, using the same vocabulary as [`cigar_query_coverage`](#cigar-query-coverage) and [`cigar_query_intervals`](#cigar-query-intervals) — `'aligned'` counts `M`/`=`/`X`, `'mapped'` additionally counts inserted bases. Pass it as `coverage_type := 'mapped'`. Note this knob also moves `identity` and `n_fragments`: a record covering no query positions under `'aligned'` is excluded from the read's group entirely, and `'mapped'` admits it, so its `=`/`X` counts join the pooled identity too

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `read_id` | (inherits) | grouping key |
| `is_read1` | BOOLEAN | the SAM read1 bit; keeps paired mates apart |
| `reference` | (inherits) | grouping key |
| `coverage` | DOUBLE | union of the fragments' read footprints, over the read length |
| `identity` | DOUBLE | [`cigar_pooled_identity`](#pooled-sequence-identity) over the fragments — **see the caveat below** |
| `query_length` | BIGINT | full read length, including hard-clipped bases |
| `n_fragments` | BIGINT | alignment records pooled |
| `mixed_strand` | BOOLEAN | some pair of query-adjacent fragments lie on opposite strands |
| `max_ref_gap` | BIGINT | largest reference gap, modulo the reference length, between query-adjacent same-strand fragments. NULL when there is no such pair to measure |

#### The evidence columns are reported, not enforced

Two fragments of one read can reach `coverage` 1.0 and `identity` 1.0 while **not being an
origin span at all** — an inverted repeat, a chimera, an inversion, a misassembly. Identity
does not catch those: a chimera the aligner splits scores 1.0 on every fragment. So the macro
returns what it observed and you set the policy:

```sql
SELECT * FROM circular_query_coverage(alignments, reference_lengths)
WHERE coverage >= 0.90
  AND NOT mixed_strand
  AND COALESCE(max_ref_gap, 0) <= 100;
```

`max_ref_gap` is the distance between one fragment's reference end and the next fragment's
reference start, taken in **read order** — modulo the reference length on a reference declared
circular, and as a plain distance on one declared linear. On a circular reference, linearising
the circle is the only reason the aligner could not chain those fragments itself, so a genuine
origin span closes to 0 — or to whatever indel sits at the junction. On synthetic reads this measures
**0–2 bases for real origin spans and 48999–55000 for split chimeras**, so any tolerance from
roughly 50 to 1000 behaves identically. It is not a knife-edge.

`mixed_strand` exists because a genuine origin-spanning read is one contiguous molecule;
linearisation merely cuts it, so both fragments must lie on the same strand. Fragments on
opposite strands are not a wrap and pooling them would manufacture coverage.

> Testing abutment against the contig ends instead — `position = 1 AND stop_position = L+1` —
> looks simpler and is wrong twice over. It rejects a junction read with a couple of bases
> deleted at the origin, whose supplementary then starts at position 3 rather than 1. And it
> degenerates on a read that wraps more than once, where every fragment starts at 1 *and* ends
> at the contig end simultaneously.

#### ⚠️ `identity` requires an extended CIGAR

`identity` is **NULL** unless the aligner wrote an extended CIGAR using `=` and `X`. bowtie2,
bwa-mem and minimap2 without `eqx` all emit legacy `M`, which records that a base aligned but
not whether it matched — so identity is not recoverable from it, and
[`cigar_sequence_identity`](#sequence-identity-from-cigar) returns NULL on the same input for
the same reason. It is also NULL when the fragments **disagree** about the encoding — one in
legacy `M`, another extended — because the `=`/`X` counts would then describe only part of the
read.

Because `NULL >= 0.99` is NULL rather than false, **adding `AND identity >= 0.99` to the gate
above does not merely lose the identity check — it rejects every row**, silently discarding the
origin-spanning reads this macro exists to rescue. To gate on identity, either align with
extended CIGARs:

```sql
SELECT * FROM align_minimap2('reads', subject_table := 'refs',
                             preset := 'map-hifi', eqx := true);
```

Re-aligning is the only route to a pooled identity for a legacy-`M` read.
[`alignment_seq_identity`](#sequence-identity) with the NM or MD tag reconstructs what `M`
omits, but only one record at a time, and this macro computes `identity` internally — so it
cannot be substituted here, and the per-record figures it gives cannot be averaged into the
read's identity (see [Pooled sequence identity](#pooled-sequence-identity) for why). Every
other column is computable from legacy `M` CIGARs, so coverage and the topology evidence work
regardless of aligner.

**Preconditions and exclusions:**
- **Secondary records are excluded.** They re-place query bases a sibling already covers, so
  they cannot raise a union but would inflate `n_fragments`.
- **Unmapped records are excluded**, as are references absent from `reference_lengths` — the
  latter mirrors [`genome_coverage`](#genome-coverage) excluding contigs absent from its
  mapping relation.
- **Records with a NULL `read_id` are excluded.** `read_alignments` maps the `*` QNAME sentinel
  to NULL, and such records cannot be grouped; pooling every unnamed record together would
  merge unrelated molecules.
- **`is_circular` is required per reference, and NULL is rejected.** Only the caller knows
  which references are circular. The identical pair of same-strand fragments — one ending
  where the contig ends, the next starting at its origin — is a wrap on a circular reference
  and an end-join on a linear one, and adapter chimeras and assembly-graph artifacts produce
  exactly that shape. Assuming circular would admit them as perfect origin spans, with the
  assumption invisible at the call site. A reference recorded as both circular and linear
  raises, like a reference with two different lengths.
- **Linear references keep their rows.** Only the *wrap* reading of the gap depends on
  circularity; pooling a read's records is worth doing either way, since a read split across
  a supplementary is one molecule whatever the reference's topology. So `coverage`,
  `identity` and the rest are computed identically, and only `max_ref_gap` changes: on a
  linear reference it is the plain distance between the fragments, so an end-join reports the
  whole contig rather than 0 and fails the gate.
- **A read explained by two references produces one row per reference.** Each says how much
  of the read that reference explains, which is the question being asked — but it means the
  rows do not partition the read. A recruitment gate looks at one reference at a time and is
  unaffected. **Do not sum `coverage` across references to get abundance**: a read explained
  by three references contributes 3.0. For abundance use
  [`woltka_ogu`](profiling.md#woltka_ogu), which distributes multi-mapped reads fractionally
  instead of double-counting them. To reduce to one reference per read first:
  ```sql
  -- best reference per read, ties broken deterministically
  SELECT * FROM circular_query_coverage(alignments, reference_lengths)
  QUALIFY ROW_NUMBER() OVER (PARTITION BY read_id, is_read1
                             ORDER BY coverage DESC, identity DESC NULLS LAST, reference) = 1;

  -- or keep only reads that one reference explains unambiguously
  SELECT * FROM circular_query_coverage(alignments, reference_lengths)
  QUALIFY COUNT(*) OVER (PARTITION BY read_id, is_read1) = 1;
  ```
  Both are policies, not facts, which is why the macro reports the rows and leaves the choice
  to you.
- **Records with a NULL `position` or `stop_position` are excluded.** They cannot be placed on
  the reference, so they cannot take part in the contiguity test. Keeping them would be
  actively unsafe: such a record counts toward `coverage` and `n_fragments` while its gap is
  unmeasurable, and since `max_ref_gap` ignores unmeasurable pairs, a chimera could reach
  coverage 1.0 with `max_ref_gap` NULL and pass `COALESCE(max_ref_gap, 0) <= 100` on a gap
  that was never checked. Dropping the record instead lowers the read's coverage by whatever
  it would have explained, so the read fails a coverage floor rather than slipping through.
  Well-formed mapped records always carry both columns.
- **Paired mates are never pooled.** R1 and R2 share a `read_id` but are different molecules,
  so `is_read1` is part of the grouping key and appears in the output.
- **`reference_lengths` may repeat a reference but not contradict itself.** Duplicate identical
  rows are collapsed; two different lengths for one reference raise. This matters because
  `read_alignment_header()` over a glob, or a `UNION ALL` of headers, yields one row per contig
  *per file*.
- **A missing or non-positive `length` raises.** It cannot be tolerated: modulo by NULL or zero
  yields NULL in DuckDB, `max_ref_gap` would become NULL, and `COALESCE(max_ref_gap, 0)` would
  then read "unmeasured" as "perfectly abutting" — turning the one column that distinguishes an
  origin span from a chimera into a false accept.
- **All records of one read must report the same query length.** Disagreement raises: it means
  the relation grouped rows from more than one read under a single `read_id`.

**Example:**
```sql
CREATE VIEW recruited AS
  SELECT * FROM align_minimap2('sample_reads', subject_table := 'sample_assembly',
                               preset := 'map-hifi', eqx := true);

-- Circularity is a claim about the assembly, not something the alignment records, so it is
-- stated here. Most assemblers flag circular contigs in the FASTA header or a companion
-- table; substitute your own predicate for the name match.
CREATE VIEW contig_lengths AS
  SELECT read_id AS reference, length(sequence1) AS length,
         read_id LIKE '%circular%' AS is_circular
  FROM sample_assembly;

-- Reads confidently recruited to their own sample's assembly, including those crossing
-- the origin of a circular contig.
SELECT read_id, reference, coverage, identity
FROM circular_query_coverage(recruited, contig_lengths)
WHERE coverage >= 0.90 AND identity >= 0.99
  AND NOT mixed_strand AND COALESCE(max_ref_gap, 0) <= 100;
```

> **Calling convention.** Like [`genome_coverage`](#genome-coverage), the relation arguments are
> resolved by name, so they must be tables or views — not inline subqueries. They also cannot
> appear in a lateral position: a scalar subquery calling this macro from inside a query whose
> `FROM` already names the same relation fails to bind (`query_table` accepts only literals).
> Use a CTE or a join instead.

**Performance notes:**
- Grouping is done by DuckDB's own window and hash-aggregate operators rather than by
  bespoke machinery, and it handles per-read grouping over large alignment sets well —
  roughly a second per two million records on one machine. It is **not** a streaming
  pipeline, though: ordering fragments along the read requires a window operator, which
  materialises and sorts its input before emitting anything, so size memory for the whole
  input rather than for a constant-space scan. Per-group state is small and bounded
  (interval endpoints and four counters per read).
- Reference-side questions are a different axis: for how many times a read covers a reference,
  or per-position depth, use [`compress_intervals`](#merge-overlapping-intervals) or
  [`compute_coverage_depth`](#coverage-depth). Query coverage is capped at 1.0 by definition and
  says nothing about multiplicity.

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

### Worked example: joint identity filtering of concordant read pairs

A per-row identity filter breaks paired-end data: `WHERE cigar_sequence_identity(cigar) >= 0.99` is evaluated one alignment at a time, so it can keep one mate and drop the other, leaving orphaned half-pairs. This recipe applies an identity threshold to the **pair as a unit** — both mates of a concordant alignment are kept or dropped together — using a windowed `QUALIFY` instead of a `WHERE`.

**Applies to** any concordant paired-end alignment source with extended (`=`/`X`) CIGAR: `read_alignments()` on a paired BAM, [`align_bowtie2`](alignment_reference.md#bowtie2), or `align_bowtie2_sharded`. The aligner must emit `=`/`X` ops (`--xeq` / `xeq := true`) so `cigar_sequence_identity` can work from the CIGAR alone; a legacy `M`-only CIGAR yields NULL. Run with `no_discordant`/`no_mixed` (or otherwise restrict to concordant pairs) so every emitted record is part of a mapped pair.

**The pair key.** Group the two mate rows of one concordant placement with:

```sql
PARTITION BY read_id, reference,
             LEAST(position, mate_position),
             GREATEST(position, mate_position)
```

- `read_id` is shared by both mates (QNAME); `reference` is the shared contig of a concordant pair.
- The two mate rows carry `position`/`mate_position` **swapped** (mate 1: `pos=P1, mate_pos=P2`; mate 2: `pos=P2, mate_pos=P1`). `LEAST`/`GREATEST` normalizes both to the same unordered coordinate pair `{P1, P2}` so they land in the same partition.
- With `report_all`, a read pair can align concordantly at several places. bowtie2 emits each placement as its own self-consistent 2-record pair (duplicating the shared mate, secondary bit set on extras), so each placement gets a distinct `{P1, P2}` and is judged independently.

**Choosing the pair-level metric.** Each partition holds exactly the two mates, so the aggregate reduces to a two-value combine:

| Intent | Aggregate in `QUALIFY` |
|---|---|
| Both mates must pass (strict) | `MIN(cigar_sequence_identity(cigar)) OVER pair >= t` |
| At least one mate passes | `MAX(cigar_sequence_identity(cigar)) OVER pair >= t` |
| Mean of the two mate identities (equal weight per mate) | `AVG(cigar_sequence_identity(cigar)) OVER pair >= t` |
| Fragment-pooled identity (weighted by aligned length) | [`cigar_pooled_identity(cigar)`](#pooled-sequence-identity)` OVER pair >= t` |

The pooled form sums both mates' counters and scores them once, giving exactly `(matches₁+matches₂)/(cols₁+cols₂)`. It differs from `AVG` only when the two mates align over different numbers of columns (the longer mate gets more weight) — the correct behavior when read/aligned lengths are unequal.

**Example** — keep only concordant pairs whose fragment-pooled identity is ≥ 0.99:

```sql
WITH aln AS (
  SELECT * FROM align_bowtie2_sharded('reads',
    shard_directory := 'shards', read_to_shard := 'read_to_shard',
    report_all := true, xeq := true,
    no_discordant := true, no_mixed := true /* , other bowtie2 params ... */)
)
SELECT *
FROM aln
QUALIFY cigar_pooled_identity(cigar) OVER (
          PARTITION BY read_id, reference,
                       LEAST(position, mate_position),
                       GREATEST(position, mate_position)) >= 0.99;
```

**Gotchas:**
- Clause order: `QUALIFY` must come **after** any named `WINDOW` clause. To avoid the ordering entirely, inline the window with `OVER ( ... )` as above.
- If one mate is unmapped the pooled figure degrades to the other mate's identity, which is safe under `no_mixed`/`no_discordant` (both mates always mapped) but worth a guard if you relax those.
- This filter is a partitioned window aggregate (roughly one pass over the data); it scales linearly and comfortably handles tens of millions of alignment rows.

