# Scalar Functions

Scalar functions for alignment analysis and sequence processing.

## Table of Contents

- [SAM Flag Functions](#sam-flag-functions) - Test individual SAM flag bits
- [`alignment_seq_identity`](#alignment_seq_identitycigar-nm-md-type) - Sequence identity calculation
- [`alignment_query_length`](#alignment_query_lengthcigar-include_hard_clipstrue) - Query length from CIGAR
- [`alignment_query_coverage`](#alignment_query_coveragecigar-typealigned) - Query coverage from CIGAR
- [`mask_dust`](#mask_dustsequence-hardmaskfalse) - DUST low-complexity masking
- [`merge_pairs_vsearch`](#merge_pairsfwd_seq-fwd_qual-rev_seq-rev_qual-options) - Paired-end read merging

## SAM Flag Functions

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

## `alignment_seq_identity(cigar, nm, md, type)`

Calculate sequence identity between read and reference using three different methods. They are derived from Heng Li's [blog post](https://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity).

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

## `alignment_query_length(cigar, [include_hard_clips=true])`

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
SELECT read_id, alignment_query_length(cigar) AS query_len
FROM read_alignments('alignments.sam');

-- Get query length excluding hard clips (matches bam_cigar2qlen)
SELECT read_id, alignment_query_length(cigar, false) AS query_len
FROM read_alignments('alignments.sam');

-- Compare lengths with and without hard clips
SELECT read_id, cigar,
  alignment_query_length(cigar, true) AS len_with_hard,
  alignment_query_length(cigar, false) AS len_without_hard
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%H%';

-- Calculate average query length per reference
SELECT reference, AVG(alignment_query_length(cigar)) AS avg_query_len
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Note:** When `include_hard_clips=false`, this function's output matches HTSlib's `bam_cigar2qlen` behavior, which counts M, I, S, =, and X operations.

## `alignment_query_coverage(cigar, [type='aligned'])`

Calculate the proportion of query bases covered by the reference alignment. This helps assess how much of a read actually aligns versus being clipped.

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
SELECT read_id, alignment_query_coverage(cigar) AS aligned_cov
FROM read_alignments('alignments.sam');

-- Get mapped coverage (includes insertions)
SELECT read_id, alignment_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam');

-- Compare aligned vs mapped coverage
SELECT read_id, cigar,
  alignment_query_coverage(cigar, 'aligned') AS aligned_cov,
  alignment_query_coverage(cigar, 'mapped') AS mapped_cov
FROM read_alignments('alignments.sam')
WHERE cigar LIKE '%I%';  -- Reads with insertions show the difference

-- Filter reads with high query coverage
SELECT COUNT(*)
FROM read_alignments('alignments.bam')
WHERE alignment_query_coverage(cigar, 'aligned') > 0.9;

-- Find heavily clipped reads
SELECT read_id, cigar, alignment_query_coverage(cigar) AS coverage
FROM read_alignments('alignments.sam')
WHERE alignment_query_coverage(cigar) < 0.5
ORDER BY coverage;

-- Calculate average coverage per reference
SELECT reference,
  AVG(alignment_query_coverage(cigar, 'aligned')) AS avg_aligned_cov,
  AVG(alignment_query_coverage(cigar, 'mapped')) AS avg_mapped_cov
FROM read_alignments('alignments.bam')
WHERE NOT alignment_is_unmapped(flags)
GROUP BY reference;
```

**Use cases:**
- **Aligned coverage**: Assess quality of alignment (how much of read actually matches reference)
- **Mapped coverage**: Identify heavily clipped reads (adapters, chimeras, low-quality ends)
- Filter reads based on alignment quality thresholds
- QC metrics for sequencing runs

---

## `mask_dust(sequence, [hardmask=false])`

DUST low-complexity masking, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Identifies low-complexity regions (e.g., homopolymers, dinucleotide repeats) and masks them.

**Parameters:**
- `sequence` (VARCHAR): DNA sequence to mask.
- `hardmask` (BOOLEAN, default false): If false, soft-mask by lowercasing. If true, replace low-complexity regions with `N`.

**Returns:** VARCHAR — the masked sequence.

**Example:**
```sql
-- Soft-mask (lowercase low-complexity regions)
SELECT mask_dust('AAAAAAAAAAAACCCCCCCC');

-- Hard-mask (replace with N)
SELECT mask_dust('AAAAAAAAAAAACCCCCCCC', true);

-- Mask sequences in a table
SELECT read_id, mask_dust(sequence1) AS masked_seq
FROM read_fastx('sequences.fasta');
```

**Behavior:**
- High-complexity sequences are returned unchanged
- NULL input returns NULL
- Empty string returns empty string

---

## `merge_pairs_vsearch(fwd_seq, fwd_qual, rev_seq, rev_qual, [options])`

Paired-end read merging, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Merges overlapping forward and reverse reads into a single consensus sequence with merged quality scores.

**Parameters (4-arg form, default tuning):**
- `fwd_seq` (VARCHAR): Forward read sequence.
- `fwd_qual` (LIST(UTINYINT)): Forward read quality scores (numeric Phred, as output by `read_fastx`).
- `rev_seq` (VARCHAR): Reverse read sequence.
- `rev_qual` (LIST(UTINYINT)): Reverse read quality scores.

**Parameters (10-arg form, with tuning):**
- Same first 4 arguments, plus:
- `minovlen` (INTEGER): Minimum overlap length.
- `maxdiffs` (INTEGER): Maximum differences in overlap region.
- `maxdiffpct` (DOUBLE): Maximum difference percentage in overlap.
- `maxee` (DOUBLE): Maximum expected errors in merged sequence.
- `minlen` (INTEGER): Minimum merged sequence length.
- `maxlen` (INTEGER): Maximum merged sequence length.

**Returns:** STRUCT with fields:

| Field | Type | Description |
|-------|------|-------------|
| `merged` | BOOLEAN | Whether merging succeeded |
| `sequence` | VARCHAR | Merged sequence (NULL if not merged) |
| `quality` | LIST(UTINYINT) | Merged quality scores (NULL if not merged) |
| `ee_merged` | DOUBLE | Expected errors of merged sequence |
| `ee_fwd` | DOUBLE | Expected errors of forward read |
| `ee_rev` | DOUBLE | Expected errors of reverse read |
| `fwd_errors` | INTEGER | Errors in forward read overlap region |
| `rev_errors` | INTEGER | Errors in reverse read overlap region |
| `overlap` | INTEGER | Length of overlap region |

**Example:**
```sql
-- Merge paired-end reads
SELECT read_id,
       (merge_pairs_vsearch(sequence1, qual1, sequence2, qual2)).merged AS merged,
       (merge_pairs_vsearch(sequence1, qual1, sequence2, qual2)).sequence AS merged_seq
FROM read_fastx('forward.fq', 'reverse.fq');

-- Filter to only successfully merged reads
SELECT read_id, m.*
FROM read_fastx('forward.fq', 'reverse.fq'),
     LATERAL (SELECT merge_pairs_vsearch(sequence1, qual1, sequence2, qual2)) AS m(result)
WHERE m.result.merged;
```

**Behavior:**
- NULL inputs return a STRUCT with `merged=false` and NULL fields
- Input length guard: throws if forward + reverse > 9,999 bases (vsearch fixed buffer)
- Quality inputs and outputs use numeric Phred (LIST(UTINYINT)), matching `read_fastx` output
