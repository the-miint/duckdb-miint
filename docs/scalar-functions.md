# Scalar Functions

Scalar functions for alignment analysis and sequence processing.

## Table of Contents

- [SAM Flag Functions](#sam-flag-functions) - Test individual SAM flag bits
- [`alignment_seq_identity`](#alignment_seq_identitycigar-nm-md-type) - Sequence identity calculation (multi-mode, requires NM/MD for legacy CIGAR)
- [`cigar_sequence_identity`](#cigar_sequence_identitycigar) - Sequence identity from extended CIGAR alone
- [`cigar_query_length`](#cigar_query_lengthcigar-include_hard_clipstrue) - Query length from CIGAR
- [`cigar_query_coverage`](#cigar_query_coveragecigar-typealigned) - Query coverage from CIGAR
- [`mask_dust`](#mask_dustsequence-hardmaskfalse) - DUST low-complexity masking
- [QC functions](qc.md) - Adapter / polyG / polyX / quality trimming and per-read filtering (fastp port)
- [`merge_pairs_vsearch`](#merge_pairsfwd_seq-fwd_qual-rev_seq-rev_qual-options) - Paired-end read merging
- [`extract_linked_amplicon`](#extract_linked_ampliconseq-qual-anchor5-anchor3-min_len-max_len-error_rate-min_overlap) - Cut out the interior between two flanking adapters (cutadapt `-g X...Y` equivalent; `^`/`$` sigils select anchored mode)
- [`phylogeny_fasttree_available`](#phylogeny_fasttree_available) - Probe for the gpl-boundary daemon at runtime
- [`install_gpl_boundary`](#install_gpl_boundary) - Download and install the gpl-boundary binary into miint's cache

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

## `cigar_sequence_identity(cigar)`

One-arg convenience wrapper over `alignment_seq_identity(cigar, NULL, NULL, 'cigar')`. Use this when your CIGAR uses extended `=`/`X` ops and you don't need NM/MD-based identity flavors.

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

## `cigar_query_coverage(cigar, [type='aligned'])`

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

## `extract_linked_amplicon(seq, qual, anchor5, anchor3, [min_len, max_len, error_rate, [min_overlap]])`

Cut out the interior of a sequence located between two flanking adapter
anchors. Equivalent to cutadapt's `-g ANCHOR5...ANCHOR3` mode but as a SQL
scalar: each row independently locates the 5' anchor, then the 3' anchor
in the remainder, and returns the slice between them along with its qual
substring. Powered by the WFA2 semi-global aligner (ends-free), so each
anchor allows up to `ceil(len(anchor) * error_rate)` errors (mismatch or
indel) — substitutions and 1-bp indels at the boundary are tolerated.

Anchors support **IUPAC degenerate bases** (N, R, Y, S, W, K, M, B, D, H, V).
A degenerate position matches its compatible bases at zero cost — it does not
consume from the error budget. For example, an `N` in the anchor matches any
base, and `R` matches A or G. This matches cutadapt's adapter-matching
behavior and is essential for real primer sets with degenerate positions.

Designed for amplicon prep workflows: UMI / index extraction, primer
trimming, locus extraction from long-read FASTQ.

### Anchored vs non-anchored mode (cutadapt sigils)

By default the function uses cutadapt's **non-anchored** semantics: each
anchor may match anywhere within its search window, and additionally the
5' anchor's 5' end may hang off the read's 5' edge (and the 3' anchor's 3'
end may hang off the read's 3' edge) — the partial overlap requires at
least `min_overlap` anchor bases to align. This is essential for long-read
amplicon protocols where edge clipping (e.g., PacBio CCS) leaves primers
truncated at the read termini.

Prefix `anchor5` with `^` and/or suffix `anchor3` with `$` to switch the
respective end to **anchored** mode, matching cutadapt's anchored-adapter
syntax. The sigil is stripped per-row before the anchor sequence is used,
so the anchor itself remains pure IUPAC bases.

| Sigil pattern | Cutadapt equivalent | Meaning |
|---|---|---|
| `extract_linked_amplicon(..., 'X', 'Y', ...)` | `-g X...Y` | Both non-anchored (default) — partial overlap allowed at the read's outer edges |
| `extract_linked_amplicon(..., '^X', 'Y', ...)` | `-g ^X...Y` | 5' anchored: `X` must match end-to-end inside the window |
| `extract_linked_amplicon(..., 'X', 'Y$', ...)` | `-g X...Y$` | 3' anchored: `Y` must match end-to-end inside the remainder |
| `extract_linked_amplicon(..., '^X', 'Y$', ...)` | `-g ^X...Y$` | Both anchored: no edge-overlap recovery |

Within each call, anchored matching is tried first; the partial-overlap
fallback runs only when no in-window match is found. This deterministically
prefers full anchor matches over partial-overlap matches (cutadapt's
default behavior) and avoids WFA2 tie-breaks that would otherwise pull
the alignment toward maximal free trim.

Only the first `^` on `anchor5` and the last `$` on `anchor3` are recognized
as sigils — they are stripped exactly once. Repeated sigils (e.g., `'^^X'`)
leave the second character in place and the result is searched literally.
Post-strip anchor characters must be IUPAC bases (`[ACGTURYSWKMBDHVN]` plus
DNA/RNA case-insensitive); other characters will simply never match and the
row will return NULL. If you are driving anchors from a column whose values
may include a literal `^` or `$` as part of the base sequence (non-standard
but possible in dirty data), strip such characters at the call site before
passing them in.

> **Migration note (changed in this release).** Before this release the
> function silently behaved like cutadapt's anchored-both `-g ^X...Y$` mode
> regardless of input — reads with primers partially clipped off the read
> edges returned NULL. With this release the default is cutadapt's
> non-anchored `-g X...Y` mode (matching the docs as written), so calls
> with no sigil will now match a class of reads they previously rejected.
> If you were relying on the old NULL-as-quality-signal behavior, restore
> it by switching to `'^X', 'Y$'` (a one-line per-call-site change).

**Parameters (4-arg form, default tuning):**

| Parameter | Type | Description |
|-----------|------|-------------|
| `seq` | VARCHAR | The full read sequence |
| `qual` | LIST(UTINYINT) | Per-base Phred qual (no ASCII offset; matches `read_fastx`) |
| `anchor5` | VARCHAR | 5' flanking adapter; prefix with `^` to anchor it to the read 5' (cutadapt `-g ^X...`) |
| `anchor3` | VARCHAR | 3' flanking adapter (found in the remainder *after* anchor5); suffix with `$` to anchor it to the read 3' (cutadapt `-g X...Y$`) |

**Parameters (7-arg form, with tuning):**

Same first 4 plus:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_len` | BIGINT | 0 | Minimum interior length (rejected as NULL if shorter) |
| `max_len` | BIGINT | 2³¹−1 | Maximum interior length (rejected as NULL if longer) |
| `error_rate` | DOUBLE | 0.10 | Per-anchor error budget — used as `ceil(len(anchor) * error_rate)` |

**Parameters (8-arg form, with min_overlap tuning):**

Same first 7 plus:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_overlap` | BIGINT | 3 | Minimum anchor bases required to align in non-anchored mode. Matches cutadapt's `-O`. Floors how short a partial-overlap-at-edge match may be; lower values increase recovery of edge-clipped reads at the cost of more false positives on long-read data. Must be `>= 1`; `0` is rejected at bind time because it would let the entire anchor be trimmed for free. Ignored on the side(s) anchored via `^` / `$`. If `min_overlap >= len(anchor)` the partial path silently degenerates to anchored on that side (e.g., a 4-bp degenerate anchor with `min_overlap=10` behaves as `^X`). |

**Returns:** STRUCT with fields:

| Field | Type | Description |
|-------|------|-------------|
| `sequence` | VARCHAR | The extracted interior (between the two anchors) |
| `qual` | LIST(UTINYINT) | Per-base qual slice over the same range |
| `start` | INTEGER | 0-based start position of the interior in `seq` |
| `stop` | INTEGER | 0-based exclusive end position of the interior in `seq` |

Returns `NULL` if either anchor cannot be located within its error budget,
or if the interior length falls outside `[min_len, max_len]`.

**Example:**
```sql
-- UMI extraction: pull out the 18-bp UMI between fixed flanking primers
SELECT read_id,
       (e).sequence AS umi,
       (e).qual AS umi_qual
FROM (
    SELECT read_id,
           extract_linked_amplicon(sequence1, qual1, 'GTAATACG', 'AGAGCACAC',
                                   min_len := 18, max_len := 18) AS e
    FROM read_fastx('reads.fq')
) WHERE umi IS NOT NULL;

-- Long-read amplicon primer trimming with edge-clip recovery: PacBio CCS
-- reads often have 4-15 bp clipped off the primer at the read termini.
-- The default non-anchored mode recovers these reads via partial overlap.
SELECT read_id,
       extract_linked_amplicon(sequence1, qual1,
                               'CCRAMCTGTCTCACGACG',
                               'CTGAGCCADRATCAAACYCT',
                               0, 1000, 0.10) AS amp
FROM read_fastx('reads.fq');

-- Strict mode: require the primer to be fully present at both ends
-- (replicates cutadapt's `-g ^X...Y$`).
SELECT read_id,
       extract_linked_amplicon(sequence1, qual1,
                               '^CCRAMCTGTCTCACGACG',
                               'CTGAGCCADRATCAAACYCT$',
                               0, 1000, 0.10) AS amp
FROM read_fastx('reads.fq');
```

**Behavior:**
- NULL `seq` or `qual` returns NULL.
- Anchors are searched in order: 5' first across the full sequence, then 3'
  in `seq[stop5:]`. The 3' search will not find anchors that overlap the
  5' anchor.
- Each row's qual list must have the same length as its sequence; otherwise
  the row throws.
- A per-thread `WFA2Aligner` is reused across rows via `FunctionLocalState`
  so the DP engine is allocated once per scan, not per row.
- In non-anchored mode the function makes up to two WFA2 calls per anchor
  per row (anchored first, partial fallback). Anchored mode (`^X` or `Y$`)
  and rows whose anchor matches internally use only one call.

## `phylogeny_fasttree_available()`

Returns `BOOLEAN` indicating whether the `gpl-boundary` binary (which embeds FastTree) is installed and reachable on `PATH`. Used to gate calls to [`phylogeny_fasttree`](table-functions.md#phylogeny_fasttreetable_name-options) when the daemon may not be present (e.g., distributed builds, environments without the optional binary).

**Behavior:**
- Cached after the first call (uses `std::call_once`); subsequent calls are O(1).
- The probe finds `gpl-boundary` on `PATH`, runs it with `--list-tools`, and checks whether the output advertises `fasttree`. Returns `false` if any of those steps fails.
- Returns `false` (compile-time constant) on builds where gpl-boundary support was disabled at CMake time (`MIINT_ENABLE_GPL_BOUNDARY=OFF`, including Emscripten and Windows).

**Example:**
```sql
-- Conditional fall-through: build a tree only when the daemon is available
SELECT
    CASE WHEN phylogeny_fasttree_available()
         THEN 'tree available'
         ELSE 'install gpl-boundary to enable tree-building'
    END AS status;

-- Use in a CHECK or guard
SELECT * FROM phylogeny_fasttree('seqs')
WHERE phylogeny_fasttree_available();   -- short-circuits to empty if not present
```

## `install_gpl_boundary([force])`

Download a prebuilt `gpl-boundary` binary from the upstream GitHub releases and install it into miint's cache so subsequent `phylogeny_fasttree(...)`, `align_bowtie2(...)`, and `align_bowtie2_sharded(...)` calls find it without the user editing `PATH`.

Drives the upstream `install.sh` (https://github.com/the-miint/GPL-boundary/releases/latest/download/install.sh), which detects the platform, downloads the matching tarball, verifies its SHA256 against the release's `SHA256SUMS`, and extracts the binary. Miint sets `INSTALL_DIR` to its cache (`$XDG_CACHE_HOME/miint/bin` or `$HOME/.cache/miint/bin`); `FindGplBoundary()` checks that location before falling through to `PATH`, so no user-side PATH editing is needed.

**Parameters:**
- `force` (BOOLEAN, optional, default `false`): when `true`, bypass the existing-binary probe and re-download latest into the miint cache unconditionally. Use this to escape a stale binary on `PATH` (for example `~/.cargo/bin/gpl-boundary` left behind by an old `cargo install`) that the probe considers good enough but predates a protocol bump. The cache install does NOT shadow `MIINT_GPL_BOUNDARY_PATH`; if that env var is set, unset it after a forced install or `FindGplBoundary()` will keep returning the override.

**Returns:** `STRUCT(installed BOOLEAN, path VARCHAR, version VARCHAR, message VARCHAR)`

| Field | Meaning |
|---|---|
| `installed` | `true` iff a working `gpl-boundary` is now reachable (either pre-existing or just installed). |
| `path` | Absolute path to the binary that satisfies the request, or empty string on failure. |
| `version` | Output of `gpl-boundary --version` — a JSON document containing `gpl_boundary` and the per-tool versions. Use `json_extract_string(...)` to pluck individual fields. |
| `message` | Human-readable description of what happened ("already available", "Installed gpl-boundary 0.1.0 to ...", or a diagnostic on failure). |

**Behavior:**
- **Idempotent (default):** if `gpl-boundary` is already discoverable (via `MIINT_GPL_BOUNDARY_PATH`, miint's cache, or `PATH`), the function probes it with `--version` and returns `installed=true` without touching the network. Pass `force=true` to skip this probe and re-download regardless.
- **Mutex-protected within a process:** concurrent calls within ONE DuckDB process serialize through a lock. Two SEPARATE DuckDB processes calling `install_gpl_boundary()` concurrently are not coordinated — they will both run install.sh against the same cache dir; install.sh's final `mv` is atomic on the same filesystem so the on-disk binary is always one of the two valid downloads (bit-identical for the same `latest` release), but tmpdirs and download bandwidth are wasted.
- **Network access required** for the download path (the idempotent fast path is offline-only).
- **Supported prebuilt platforms:** Linux x86_64, macOS arm64. macOS Intel and other targets must build from source — install.sh prints a "build from source" message and the function returns `installed=false` with that diagnostic in `message`.
- **Returns `installed=false`** (with a stub message) on builds where `MIINT_HAS_GPL_BOUNDARY` was off at compile time (Emscripten, Windows).

**Security model (read this if you're in a regulated environment):**

This function downloads and **executes an unverified shell script** (`install.sh`) from a GitHub releases URL, with the privileges of the DuckDB process. The shell script then downloads a tarball and verifies its SHA256 against a sibling `SHA256SUMS` file from the same release. Trust chain:

| Step | What's verified | What's trusted |
|---|---|---|
| 1. Fetch `install.sh` over HTTPS | TLS cert (OS CA store) | github.com release CDN integrity |
| 2. Fetch tarball + `SHA256SUMS` | TLS cert | install.sh's own logic |
| 3. Verify tarball SHA256 | sha256sum match | the SHA256SUMS file (same release as install.sh) |
| 4. Extract + place binary | nothing | tarball contents |

The script itself is **not** integrity-verified by miint. If you require attested provenance (signed releases, SLSA, etc.), do NOT use this function — fetch the binary out-of-band and either drop it on `PATH`, place it at `~/.cache/miint/bin/gpl-boundary`, or set `MIINT_GPL_BOUNDARY_PATH=<absolute path>` before launching DuckDB.

**Network timeouts:** miint applies `--max-time 60s` to its own `install.sh` fetch but cannot bound the script's internal curl call for the larger tarball (a few tens of MB). On a hung network, install.sh may take longer to fail than feels reasonable; if that bites you, kill the SQL session and retry.

**Cache-vs-PATH staleness:** once `install_gpl_boundary()` has populated miint's cache, [`FindGplBoundary`](#) checks the cache **before** consulting `PATH`. That means a system-level upgrade of `gpl-boundary` (via package manager, conda, etc.) will be silently shadowed by the cached version forever. To pick up a system upgrade you must either:
- delete the cached binary: `rm -rf ~/.cache/miint/bin/`
- override the lookup explicitly: `export MIINT_GPL_BOUNDARY_PATH=/usr/local/bin/gpl-boundary`

**Examples:**

```sql
-- Bootstrap before using phylogeny_fasttree
SELECT (install_gpl_boundary()).message;
-- → "gpl-boundary 0.1.0 already available at /home/x/.cargo/bin/gpl-boundary; no install performed"
-- or "Installed gpl-boundary 0.1.0 to /home/x/.cache/miint/bin/gpl-boundary"

-- Just check whether it succeeded
SELECT (install_gpl_boundary()).installed;

-- Pluck the gpl-boundary semver out of the JSON version string
SELECT json_extract_string((install_gpl_boundary()).version, '$.gpl_boundary') AS gpl_boundary_version;

-- Inspect the per-tool versions advertised by the daemon
SELECT json_extract((install_gpl_boundary()).version, '$.tools');

-- Force a re-download even though a (possibly stale) binary is already on PATH.
SELECT (install_gpl_boundary(true)).message;
-- → "Installed gpl-boundary 0.2.0 to /home/x/.cache/miint/bin/gpl-boundary"
```

**Override the install location:** set `MIINT_GPL_BOUNDARY_PATH=<absolute path>` to point at a binary already on disk (e.g., a system package install). `FindGplBoundary()` honors that override before consulting the cache or `PATH`.
