# Utilities

Odds and ends that support the main workflows: transforming sequence strings, a few amplicon read-prep primitives, reporting versions and warnings, and installing optional external tools. If you are looking for a specific task (reading, aligning, diversity, …) start from the [table of contents](table_of_contents.md); this page collects the smaller helpers that don't belong to one task.

## Table of Contents

- [Reverse complement](#reverse-complement) - `sequence_dna_reverse_complement` / `sequence_rna_reverse_complement`
- [IUPAC sequence to regex](#iupac-sequence-to-regex) - `sequence_dna_as_regexp` / `sequence_rna_as_regexp`
- [Splitting long sequences](#splitting-long-sequences) - `sequence_split`
- [Low-complexity masking](#low-complexity-masking) - `mask_dust`
- [Merging paired-end reads](#merging-paired-end-reads) - `merge_pairs_vsearch`
- [Extracting linked amplicons](#extracting-linked-amplicons) - `extract_linked_amplicon`
- [Version and diagnostics](#version-and-diagnostics) - `miint_version` / `miint_versions` / `miint_warnings`
- [Installing optional tools](#installing-optional-tools) - `phylogeny_fasttree_available` / `bowtie2_available` / `install_gpl_boundary`

## Reverse complement

`sequence_dna_reverse_complement(sequence)` and `sequence_rna_reverse_complement(sequence)`

Calculate the reverse complement of DNA or RNA sequences. Supports full IUPAC nucleotide ambiguity codes and preserves case.

**Parameters:**
- `sequence` (VARCHAR): DNA or RNA sequence string

**Returns:** VARCHAR — the reverse complement of the input sequence

**Behavior:**
- Reverses the sequence order (5' to 3' becomes 3' to 5')
- Complements each base according to Watson-Crick pairing rules
- Preserves uppercase/lowercase in the input
- Supports gap characters (`.` and `-`) which map to themselves
- Strict molecular type validation: DNA function rejects U bases, RNA function rejects T bases

**Supported bases:**
- **DNA**: A<->T, C<->G, plus IUPAC codes (R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N)
- **RNA**: A<->U, C<->G, plus IUPAC codes (R<->Y, S<->S, W<->W, K<->M, B<->V, D<->H, N<->N)

**Examples:**
```sql
-- Basic DNA reverse complement
SELECT sequence_dna_reverse_complement('ATCG');
-- Returns: CGAT

-- Basic RNA reverse complement
SELECT sequence_rna_reverse_complement('AUCG');
-- Returns: CGAU

-- Works with IUPAC ambiguity codes
SELECT sequence_dna_reverse_complement('ACGTMRWSYKVHDBN.-');
-- Returns: -.NVHDBMRSWYKACGT

-- Case is preserved
SELECT sequence_dna_reverse_complement('AcGt');
-- Returns: aCgT

-- Find palindromic sequences (equal to their reverse complement)
SELECT read_id, sequence1
FROM read_fastx('sequences.fastq')
WHERE sequence1 = sequence_dna_reverse_complement(sequence1);

-- Error: DNA function rejects U bases
SELECT sequence_dna_reverse_complement('AUCG');
-- Error: Invalid DNA base 'U'
```

**IUPAC ambiguity code reference:**
R = A/G · Y = C/T(U) · S = G/C · W = A/T(U) · K = G/T(U) · M = A/C · B = C/G/T(U) · D = A/G/T(U) · H = A/C/T(U) · V = A/C/G · N = any base

## IUPAC sequence to regex

`sequence_dna_as_regexp(sequence)` and `sequence_rna_as_regexp(sequence)`

Convert DNA or RNA sequences with IUPAC ambiguity codes to regular expression patterns. Useful for pattern matching with degenerate primers and probes.

**Parameters:**
- `sequence` (VARCHAR): DNA or RNA sequence string with IUPAC codes

**Returns:** VARCHAR — regular expression pattern

**Behavior:**
- Unambiguous bases (A, C, G, T/U) remain unchanged
- Ambiguous IUPAC codes expand to character classes (e.g., R -> `[AG]`, N -> `[ACGT]`)
- Preserves uppercase/lowercase in the output
- Gap characters (`-` and `.`) convert to `.` (regex wildcard matching any character)
- Strict molecular type validation: DNA function rejects U bases, RNA function rejects T bases

**Expansion rules:**
- **Unambiguous bases**: A, C, G, T (DNA) or U (RNA) -> no brackets
- **Two-base codes**: R -> `[AG]`, Y -> `[CT]` or `[CU]`, S -> `[CG]`, W -> `[AT]` or `[AU]`, K -> `[GT]` or `[GU]`, M -> `[AC]`
- **Three-base codes**: B -> `[CGT]` or `[CGU]`, D -> `[AGT]` or `[AGU]`, H -> `[ACT]` or `[ACU]`, V -> `[ACG]`
- **Any base**: N -> `[ACGT]` or `[ACGU]`
- **Gap characters**: `-` or `.` -> `.` (matches any character)

**Examples:**
```sql
-- Degenerate primer with ambiguous positions
SELECT sequence_dna_as_regexp('ATNGG');
-- Returns: AT[ACGT]GG

-- Multiple IUPAC codes
SELECT sequence_dna_as_regexp('RYMKSW');
-- Returns: [AG][CT][AC][GT][CG][AT]

-- Use with DuckDB's regexp_matches for pattern searching
SELECT read_id, sequence1
FROM read_fastx('sequences.fastq')
WHERE regexp_matches(sequence1, sequence_dna_as_regexp('ATNGG'));

-- Find sequences matching a degenerate probe
CREATE TABLE probes AS
  SELECT 'probe1' AS name, 'GCRAA' AS sequence
  UNION ALL SELECT 'probe2', 'ATNGG';

SELECT p.name, f.read_id, f.sequence1
FROM read_fastx('sequences.fastq') f
CROSS JOIN probes p
WHERE regexp_matches(f.sequence1, sequence_dna_as_regexp(p.sequence));

-- Error: DNA function rejects U bases
SELECT sequence_dna_as_regexp('AUNGG');
-- Error: Invalid DNA base 'U'
```

**Note:** Gap characters become `.` (regex wildcard), which matches any single character. This is useful for representing unknown or variable positions in alignments.

## Splitting long sequences

`sequence_split(sequence, chunk_size)`

Split a sequence into fixed-width chunks in a single linear pass. Returns a list of `{chunk_index, chunk_data}` structs; `UNNEST` it to get one row per chunk. Intended for emitting very large records (host reference genomes, hundreds of MB to multi-GB) as bounded-width rows for chunked storage, without the quadratic blow-up of the `list_transform` + `substring` idiom.

**Parameters:**
- `sequence` (VARCHAR): the sequence (or any string) to split
- `chunk_size` (INTEGER): chunk width in bytes; must be `> 0`

**Returns:** `LIST(STRUCT(chunk_index INTEGER, chunk_data VARCHAR))`

**Behavior:**
- Chunks are exactly `chunk_size` bytes; the last chunk is the remainder (no padding).
- `chunk_index` is dense, 0-based, and ascending within a sequence.
- Byte-exact: concatenating `chunk_data` in `chunk_index` order reproduces the input.
- Empty sequence -> empty list (zero chunks), so `UNNEST` yields no rows for it.
- `NULL` sequence or `NULL` `chunk_size` -> `NULL`.
- `chunk_size <= 0` raises an error.
- Linear time and bounded (~O(L) per record) memory: a 1 GB record chunks in well under a second, where the `substring`-in-`list_transform` macro is O(L²) and takes tens of minutes.

**Examples:**
```sql
-- Fixed-width chunks; last chunk is the remainder
SELECT c.chunk_index, c.chunk_data
FROM (SELECT UNNEST(sequence_split('ACGTACGTAC', 4)) AS c)
ORDER BY c.chunk_index;
-- 0 | ACGT
-- 1 | ACGT
-- 2 | AC

-- Emit a sequence table as 64 KB chunk rows for chunked-Parquet storage
SELECT read_id, c.chunk_index, c.chunk_data
FROM (SELECT read_id, UNNEST(sequence_split(sequence, 65536)) AS c FROM reads);

-- Byte-exact round-trip
SELECT string_agg(c.chunk_data, '' ORDER BY c.chunk_index)
FROM (SELECT UNNEST(sequence_split('ACGTACGTAC', 4)) AS c);
-- Returns: ACGTACGTAC
```

**Note:** `chunk_data` is `VARCHAR`. Chunking is by byte width, which is byte-exact for ASCII sequence data; splitting arbitrary multi-byte UTF-8 text at a byte boundary can produce invalid-UTF-8 chunks.

## Low-complexity masking

`mask_dust(sequence, [hardmask=false])`

DUST low-complexity masking, powered by the [vsearch](https://github.com/torognes/vsearch) library (Rognes et al. 2016, PeerJ 4:e2584). Identifies low-complexity regions (e.g., homopolymers, dinucleotide repeats) and masks them.

**Parameters:**
- `sequence` (VARCHAR): DNA sequence to mask.
- `hardmask` (BOOLEAN, default false): If false, soft-mask by lowercasing. If true, replace low-complexity regions with `N`.

**Returns:** VARCHAR — the masked sequence.

**Behavior:**
- High-complexity sequences are returned unchanged
- NULL input returns NULL
- Empty string returns empty string

**Examples:**
```sql
-- Soft-mask (lowercase low-complexity regions)
SELECT mask_dust('AAAAAAAAAAAACCCCCCCC');

-- Hard-mask (replace with N)
SELECT mask_dust('AAAAAAAAAAAACCCCCCCC', true);

-- Mask sequences in a table
SELECT read_id, mask_dust(sequence1) AS masked_seq
FROM read_fastx('sequences.fasta');
```

## Merging paired-end reads

`merge_pairs_vsearch(fwd_seq, fwd_qual, rev_seq, rev_qual, [options])`

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

**Behavior:**
- NULL inputs return a STRUCT with `merged=false` and NULL fields
- Input length guard: throws if forward + reverse > 9,999 bases (vsearch fixed buffer)
- Quality inputs and outputs use numeric Phred (LIST(UTINYINT)), matching `read_fastx` output

**Examples:**
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

## Extracting linked amplicons

`extract_linked_amplicon(seq, qual, anchor5, anchor3, [min_len, max_len, error_rate, [min_overlap]])`

Cut out the interior of a sequence located between two flanking adapter anchors. Equivalent to cutadapt's `-g ANCHOR5...ANCHOR3` mode but as a SQL scalar: each row independently locates the 5' anchor, then the 3' anchor in the remainder, and returns the slice between them along with its qual substring. Powered by the WFA2 semi-global aligner (ends-free), so each anchor allows up to `ceil(len(anchor) * error_rate)` errors (mismatch or indel) — substitutions and 1-bp indels at the boundary are tolerated.

Anchors support **IUPAC degenerate bases** (N, R, Y, S, W, K, M, B, D, H, V). A degenerate position matches its compatible bases at zero cost — it does not consume from the error budget. For example, an `N` in the anchor matches any base, and `R` matches A or G. This matches cutadapt's adapter-matching behavior and is essential for real primer sets with degenerate positions.

Designed for amplicon prep workflows: UMI / index extraction, primer trimming, locus extraction from long-read FASTQ. It is the WFA2-powered primitive referenced from the [quality control](qc.md) doc's long-read amplicon pipeline.

### Anchored vs non-anchored mode (cutadapt sigils)

By default the function uses cutadapt's **non-anchored** semantics: each anchor may match anywhere within its search window, and additionally the 5' anchor's 5' end may hang off the read's 5' edge (and the 3' anchor's 3' end may hang off the read's 3' edge) — the partial overlap requires at least `min_overlap` anchor bases to align. This is essential for long-read amplicon protocols where edge clipping (e.g., PacBio CCS) leaves primers truncated at the read termini.

Prefix `anchor5` with `^` and/or suffix `anchor3` with `$` to switch the respective end to **anchored** mode, matching cutadapt's anchored-adapter syntax. The sigil is stripped per-row before the anchor sequence is used, so the anchor itself remains pure IUPAC bases.

| Sigil pattern | Cutadapt equivalent | Meaning |
|---|---|---|
| `extract_linked_amplicon(..., 'X', 'Y', ...)` | `-g X...Y` | Both non-anchored (default) — partial overlap allowed at the read's outer edges |
| `extract_linked_amplicon(..., '^X', 'Y', ...)` | `-g ^X...Y` | 5' anchored: `X` must match end-to-end inside the window |
| `extract_linked_amplicon(..., 'X', 'Y$', ...)` | `-g X...Y$` | 3' anchored: `Y` must match end-to-end inside the remainder |
| `extract_linked_amplicon(..., '^X', 'Y$', ...)` | `-g ^X...Y$` | Both anchored: no edge-overlap recovery |

Within each call, anchored matching is tried first; the partial-overlap fallback runs only when no in-window match is found. This deterministically prefers full anchor matches over partial-overlap matches (cutadapt's default behavior) and avoids WFA2 tie-breaks that would otherwise pull the alignment toward maximal free trim.

Only the first `^` on `anchor5` and the last `$` on `anchor3` are recognized as sigils — they are stripped exactly once. Repeated sigils (e.g., `'^^X'`) leave the second character in place and the result is searched literally. Post-strip anchor characters must be IUPAC bases (`[ACGTURYSWKMBDHVN]` plus DNA/RNA case-insensitive); other characters will simply never match and the row will return NULL. If you are driving anchors from a column whose values may include a literal `^` or `$` as part of the base sequence (non-standard but possible in dirty data), strip such characters at the call site before passing them in.

> **Migration note (changed in this release).** Before this release the function silently behaved like cutadapt's anchored-both `-g ^X...Y$` mode regardless of input — reads with primers partially clipped off the read edges returned NULL. With this release the default is cutadapt's non-anchored `-g X...Y` mode (matching the docs as written), so calls with no sigil will now match a class of reads they previously rejected. If you were relying on the old NULL-as-quality-signal behavior, restore it by switching to `'^X', 'Y$'` (a one-line per-call-site change).

**Parameters (4-arg form, default tuning):**

| Parameter | Type | Description |
|-----------|------|-------------|
| `seq` | VARCHAR | The full read sequence |
| `qual` | LIST(UTINYINT) | Per-base Phred qual (no ASCII offset; matches `read_fastx`) |
| `anchor5` | VARCHAR | 5' flanking adapter; prefix with `^` to anchor it to the read 5' (cutadapt `-g ^X...`) |
| `anchor3` | VARCHAR | 3' flanking adapter (found in the remainder *after* anchor5); suffix with `$` to anchor it to the read 3' (cutadapt `-g X...Y$`) |

**Parameters (7-arg form, with tuning):** same first 4 plus:

| Parameter | Type | Default | Description |
|-----------|------|---------|-------------|
| `min_len` | BIGINT | 0 | Minimum interior length (rejected as NULL if shorter) |
| `max_len` | BIGINT | 2³¹−1 | Maximum interior length (rejected as NULL if longer) |
| `error_rate` | DOUBLE | 0.10 | Per-anchor error budget — used as `ceil(len(anchor) * error_rate)` |

**Parameters (8-arg form, with min_overlap tuning):** same first 7 plus:

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

Returns `NULL` if either anchor cannot be located within its error budget, or if the interior length falls outside `[min_len, max_len]`.

**Examples:**
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
- Anchors are searched in order: 5' first across the full sequence, then 3' in `seq[stop5:]`. The 3' search will not find anchors that overlap the 5' anchor.
- Each row's qual list must have the same length as its sequence; otherwise the row throws.
- A per-thread `WFA2Aligner` is reused across rows via `FunctionLocalState` so the DP engine is allocated once per scan, not per row.
- In non-anchored mode the function makes up to two WFA2 calls per anchor per row (anchored first, partial fallback). Anchored mode (`^X` or `Y$`) and rows whose anchor matches internally use only one call.

## Version and diagnostics

### `miint_version()`

Returns the MIINT extension version string.

**Returns:** VARCHAR — version string derived from the git tag at build time (e.g., `v0.1.0`, or `v0.1.0-3-gabcdef1` if built from a commit after a tag).

```sql
SELECT miint_version();
```

### `miint_versions()`

Returns a table of the pinned versions of miint and every embedded library, one row per component. Which rows appear depends on the build's compile-time features (HDF5, vsearch, abpoa, sylph, unifrac, libcurl, libdeflate are conditional).

**Output schema:**

| Column | Type | Description |
|--------|------|-------------|
| `library` | VARCHAR | Component name (`miint`, `htslib`, `minimap2`, `kseq++`, `LBFGS++`, `WFA2-lib`, `zlib`, `rype`, `mafft`, and conditionally `HDF5`, `libdeflate`, `libcurl`, `vsearch`, `abpoa`, `sylph`, `unifrac`, `scikit-bio-binaries`) |
| `version` | VARCHAR | Version / git description for that component |

```sql
-- All embedded component versions
SELECT * FROM miint_versions();

-- Audit specific libraries (e.g., before running a diversity workflow)
SELECT version FROM miint_versions()
WHERE library IN ('unifrac', 'scikit-bio-binaries');
```

### `miint_warnings()`

Returns miint's operational warnings as a queryable table. Populated by user-facing warning sites — skipped accessions in `read_ena_sequences`, the SFF `max_sequences` caveat, the `threads` parameter being ignored in `align_bowtie2_sharded`, a [feature table not stored in `sample_id` order](diversity.md#feature-table-sort-order) for progressive UniFrac, mid-stream download failures, etc. Every entry is also printed to stderr, so interactive users see no change; pipeline and `COPY TO` users now have a way to inspect warnings after the fact.

**Output schema:**

| Column | Type | Description |
|---|---|---|
| `timestamp` | `TIMESTAMP WITH TIME ZONE` | When the warning was emitted |
| `message` | `VARCHAR` | Human-readable warning text |

```sql
-- Scan a study, then see what got skipped
SELECT COUNT(*) FROM read_ena_sequences('PRJEB1234');
SELECT timestamp, message FROM miint_warnings();
```

**Behavior:**
- The log is process-scoped and in-memory by default; entries accumulate across queries within a single DuckDB session.
- Warnings are captured regardless of the `enable_logging` / `logging_level` settings — miint writes directly to DuckDB's global log sink so `miint_warnings()` works without any setup.
- Under the hood, entries have type `'MiintWarning'` and live alongside DuckDB's own logs in `duckdb_logs()`; the macro just filters on type.

**Interactions with DuckDB logging settings:**
- `SET logging_storage='stderr'` — entries go straight to stderr instead of the in-memory sink, so `miint_warnings()` stops returning rows in that mode (the stderr output is the workaround). Default storage is in-memory; leave it alone unless you have reason to change it.
- `SET warnings_as_errors=true` — miint suppresses the log-storage write for that query so a skip-warning doesn't abort the retry or the partial-data path it was designed to allow. The stderr message still fires, so nothing is silently dropped, but `miint_warnings()` will not see the row for that call.

## Installing optional tools

Some capabilities are provided by an external, GPL-licensed helper binary (`gpl-boundary`, which embeds FastTree and bowtie2) that miint keeps at arm's length to stay BSD-clean. These two functions let you probe for and install it at runtime.

### `phylogeny_fasttree_available()`

Returns `BOOLEAN` indicating whether the `gpl-boundary` binary (which embeds FastTree) is installed and reachable on `PATH`. Used to gate calls to [`phylogeny_fasttree`](phylogeny.md#fasttree) when the daemon may not be present (e.g., distributed builds, environments without the optional binary).

**Behavior:**
- Cached after the first call (uses `std::call_once`); subsequent calls are O(1).
- The probe finds `gpl-boundary` on `PATH`, runs it with `--list-tools`, and checks whether the output advertises `fasttree`. Returns `false` if any of those steps fails.
- Returns `false` (compile-time constant) on builds where gpl-boundary support was disabled at CMake time (`MIINT_ENABLE_GPL_BOUNDARY=OFF`, including Emscripten and Windows).

```sql
-- Conditional fall-through: build a tree only when the daemon is available
SELECT
    CASE WHEN phylogeny_fasttree_available()
         THEN 'tree available'
         ELSE 'install gpl-boundary to enable tree-building'
    END AS status;
```

### `bowtie2_available()`

Returns `BOOLEAN` indicating whether the `gpl-boundary` binary (which embeds bowtie2) is installed and reachable, so you can gate calls to [`align_bowtie2`](alignment_reference.md) / `align_bowtie2_sharded` when the daemon may not be present. Mirrors `phylogeny_fasttree_available()`.

**Behavior:**
- Cached after the first call (uses `std::call_once`); subsequent calls are O(1).
- The probe finds `gpl-boundary` (via `MIINT_GPL_BOUNDARY_PATH`, miint's cache, or `PATH`), runs it with `--list-tools`, and returns `true` only if the output advertises both `bowtie2-align` and `bowtie2-build`. Returns `false` if any step fails.
- Returns `false` (compile-time constant) on builds where gpl-boundary support was disabled at CMake time (Emscripten, Windows).

```sql
-- Align only when the daemon is present
SELECT CASE WHEN bowtie2_available()
            THEN 'ready'
            ELSE 'run install_gpl_boundary() first' END AS status;
```

If it returns `false`, install the daemon with [`install_gpl_boundary()`](#install_gpl_boundaryforce) below.

### `install_gpl_boundary([force])`

Download a prebuilt `gpl-boundary` binary from the upstream GitHub releases and install it into miint's cache so subsequent [`phylogeny_fasttree(...)`](phylogeny.md#fasttree), [`align_bowtie2(...)`](alignment_reference.md), and `align_bowtie2_sharded(...)` calls find it without the user editing `PATH`.

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

This function downloads and **executes an unverified shell script** (`install.sh`) from a GitHub releases URL, with the privileges of the DuckDB process. The shell script then downloads a tarball and verifies its SHA256 against a sibling `SHA256SUMS` file from the same release. The script itself is **not** integrity-verified by miint. If you require attested provenance (signed releases, SLSA, etc.), do NOT use this function — fetch the binary out-of-band and either drop it on `PATH`, place it at `~/.cache/miint/bin/gpl-boundary`, or set `MIINT_GPL_BOUNDARY_PATH=<absolute path>` before launching DuckDB.

**Cache-vs-PATH staleness:** once `install_gpl_boundary()` has populated miint's cache, `FindGplBoundary()` checks the cache **before** consulting `PATH`. That means a system-level upgrade of `gpl-boundary` (via package manager, conda, etc.) will be silently shadowed by the cached version. To pick up a system upgrade you must either delete the cached binary (`~/.cache/miint/bin/`) or override the lookup with `export MIINT_GPL_BOUNDARY_PATH=/usr/local/bin/gpl-boundary`.

**Examples:**
```sql
-- Bootstrap before using phylogeny_fasttree
SELECT (install_gpl_boundary()).message;

-- Just check whether it succeeded
SELECT (install_gpl_boundary()).installed;

-- Pluck the gpl-boundary semver out of the JSON version string
SELECT json_extract_string((install_gpl_boundary()).version, '$.gpl_boundary') AS gpl_boundary_version;

-- Force a re-download even though a (possibly stale) binary is already on PATH.
SELECT (install_gpl_boundary(true)).message;
```

**Override the install location:** set `MIINT_GPL_BOUNDARY_PATH=<absolute path>` to point at a binary already on disk (e.g., a system package install). `FindGplBoundary()` honors that override before consulting the cache or `PATH`.
