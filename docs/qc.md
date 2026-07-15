# Quality Control (QC) Scalar Functions

A family of scalar functions for FASTQ read quality control: adapter trimming,
polyG/polyX tail trimming, sliding-window quality trimming, and per-read
length/quality/N-base filtering. The algorithms are ported (not vendored)
from [fastp](https://github.com/OpenGene/fastp) — see
[`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) for attribution.

These functions compose in SQL to express the canonical fastp QC pipeline
over any source of sequence + quality data — typically [`read_fastx`](reading.md#fasta-and-fastq).

## Function catalog

| Function | Returns | Purpose |
|---|---|---|
| `trim_quality_5p(seq, qual [, window, mean_q])` | trim struct | fastp `cut_front` — drop low-quality 5' bases via sliding window |
| `trim_quality_3p(seq, qual [, window, mean_q])` | trim struct | fastp `cut_tail` — drop low-quality 3' bases via sliding window |
| `trim_quality_sliding(seq, qual [, window, mean_q])` | trim struct | fastp `cut_right` — drop window-and-everything-right at first bad window |
| `trim_polyg(seq, qual [, min_len, max_mm, max_window_mean_q])` | trim struct | 3' polyG run trim (NextSeq dark-cycle cleanup), optionally quality-gated |
| `trim_polyx(seq, qual [, min_len, max_mm])` | trim struct | 3' generic homopolymer trim (any most-frequent base) |
| `trim_adapters(seq, qual, adapter [, match_revcomp, min_match, allow_pre_start])` | trim struct | 3-phase adapter match (Hamming + 1-insert + 1-delete); `adapter` is `VARCHAR` or `LIST(VARCHAR)` |
| `trim_adapters_pe(seq1, qual1, seq2, qual2 [, adapters, overlap_require, overlap_diff_limit, overlap_diff_percent_limit, match_revcomp, min_match, allow_pre_start])` | PE trim struct | fastp paired-end overlap adapter trim: infer the insert from the R1/revcomp(R2) overlap; 11-arg form adds an adapter-by-sequence fallback |
| `filter_read(seq, qual [, min_length, max_length, qualified_q, max_unqualified_pct, max_n, min_avg_q])` | filter struct | Per-read pass/fail with metrics |
| `qc_version()` | `VARCHAR` | Version string |

All scalars accept `seq VARCHAR` + `qual LIST(UTINYINT)`. Quality is the
decoded Phred integer list emitted by `read_fastx` (not ASCII), so values
are 0..93 directly.

Beyond these per-read scalars, [`infer_trim`](#reconciling-external-qc-results--infer_trim)
is a **table macro** that reconciles reads trimmed by an *external* QC tool back
to 5'/3' trim coordinates — see its section below.

## Return structs

Every trim function returns a `STRUCT` with these fields:

| Field | Type | Meaning |
|---|---|---|
| `sequence` | `VARCHAR` | The trimmed sequence (a slice of the input) |
| `quality` | `LIST(UTINYINT)` | The trimmed quality, in lockstep with `sequence` |
| `trimmed_5p` | `UINTEGER` | Number of bases removed from the 5' end |
| `trimmed_3p` | `UINTEGER` | Number of bases removed from the 3' end |

`trimmed_5p` is meaningful only for `trim_quality_5p` (and rarely
`trim_quality_sliding` when invoked at a position that drops the leading
bases). All other trimmers operate at the 3' end only, so `trimmed_5p` is
always `0` for `trim_quality_3p`, `trim_polyg`, `trim_polyx`, and
`trim_adapters`. Don't `SUM((t).trimmed_5p)` across a 3'-only pipeline
expecting non-zero — it's structurally zero.

`trim_adapters_pe` is the exception to this layout: it operates on a read
*pair* and returns a wider paired struct (`sequence1`/`quality1`/`sequence2`/
`quality2` plus `overlap_len`, `adapter_trimmed`, `trimmed1_3p`, `trimmed2_3p`).
See [its reference](#trim_adapters_peseq1-qual1-seq2-qual2--adapters-overlap_require-overlap_diff_limit-overlap_diff_percent_limit-match_revcomp-min_match-allow_pre_start) below.

`filter_read` returns:

| Field | Type | Meaning |
|---|---|---|
| `passed` | `BOOLEAN` | Whether the read passed all enabled filters |
| `fail_reason` | `VARCHAR` | `length`, `quality`, `n_base`, `too_long`, or `NULL` when passed |
| `length` | `UINTEGER` | Input length (post-trim if called downstream of trimmers) |
| `n_bases` | `UINTEGER` | Count of `N`/`n` bases in `sequence` |
| `low_qual_bases` | `UINTEGER` | Count of quality bytes strictly below `qualified_q` |
| `mean_quality` | `FLOAT` | Mean Phred over `sequence` |

All metrics are populated even when the read passes, so you can audit-filter
in SQL by any combination of fields without re-walking the bases.

## Canonical pipeline

fastp applies its per-read steps in a fixed internal order: **sliding-window
quality cut → polyG → adapter → polyX → per-read filter**. Of these, only
adapter trimming and the per-read filter are **on by default**. polyG
auto-enables for two-color instruments (NextSeq / NovaSeq) and is otherwise
off; the sliding-window quality cuts (`cut_front` / `cut_tail` / `cut_right`)
and polyX are **opt-in**. So the default fastp pipeline effectively reduces to
**polyG (two-color only) → adapter → filter**.

In SQL, chain the scalars through aliases in a single `WITH` and access struct
fields via `(alias).field`. The example follows fastp's order and includes the
opt-in quality cut for illustration:

```sql
WITH pipeline AS (
    SELECT
        sequence_index,
        read_id,
        trim_quality_3p(sequence1, qual1, 4, 20)                    AS t1,  -- opt-in (cut_tail); off in fastp default
        trim_polyg((t1).sequence, (t1).quality)                     AS t2,  -- auto-on for two-color instruments only
        trim_adapters((t2).sequence, (t2).quality, 'AGATCGGAAGAGC')  AS t3,
        filter_read((t3).sequence, (t3).quality)                    AS f
    FROM read_fastx('reads.fastq.gz')
)
SELECT
    sequence_index, read_id,
    (t3).sequence                                                AS sequence,
    (t3).quality                                                 AS quality,
    (t1).trimmed_3p + (t2).trimmed_3p + (t3).trimmed_3p          AS bases_trimmed,
    (f).fail_reason
FROM pipeline
WHERE (f).passed;
```

For paired-end reads, prefer `trim_adapters_pe` (reference below) in place of a
per-mate `trim_adapters`: it infers the adapter boundary from the R1 / R2
overlap, which is fastp's primary and most accurate PE adapter mechanism.

Each scalar is invoked exactly once per row. DuckDB's CSE optimizer
deduplicates repeated struct extractions automatically, so accessing
`(t1).sequence` and `(t1).quality` and `(t1).trimmed_3p` does not
re-evaluate `trim_quality_3p`.

### Stats and summaries

The per-row metrics fall out of standard aggregates — no custom aggregate
needed:

```sql
SELECT
    f.fail_reason,
    COUNT(*) AS n_reads,
    SUM(adapter_bases) AS adapter_bases_trimmed,
    SUM(polyg_bases)   AS polyg_bases_trimmed
FROM ( ... pipeline ... )
GROUP BY f.fail_reason;
```

### Don't use `LATERAL` — it's ~50× slower

A natural-looking alternative is `LATERAL (SELECT scalar(...))` to materialize
the struct in its own subquery. Don't. On 10M rows, LATERAL takes ~240ms
where alias-in-SELECT takes ~5ms — about 50× slower — because LATERAL is
planned as a `LogicalDependentJoin` even when the subquery has no FROM
clause. The optimizer doesn't simplify that to a plain projection. The
`(alias).field` pattern above produces a flat projection and runs at full
vectorized speed.

## Function reference

### `trim_quality_3p(seq, qual [, window, mean_q])` / `trim_quality_5p` / `trim_quality_sliding`

Three flavors of sliding-window quality trimming, ported from fastp's
`cut_tail`, `cut_front`, and `cut_right` respectively.

- **`trim_quality_3p` (cut_tail):** slide a window from the 3' end toward
  the 5' end. While the window mean is below `mean_q`, recede; stop at the
  first window that passes. Result keeps `[0, end_of_last_passing_window)`.
- **`trim_quality_5p` (cut_front):** slide a window from the 5' end toward
  the 3' end. While the window mean is below `mean_q`, advance; stop at the
  first window that passes. Result keeps `[start_of_first_passing_window, n)`.
  Note: the bases inside the passing window are kept whole, including any
  low-quality anchor bases at the window's edge.
- **`trim_quality_sliding` (cut_right):** slide a window 5'→3'. At the
  first window whose mean is below `mean_q`, drop that window AND everything
  to its right. Result keeps `[0, first_bad_window_start)`.

Defaults: `window=4`, `mean_q=20`. `window` is clamped to `1..1000` (fastp's
documented range). Read shorter than the window is a no-op.

### `trim_polyg(seq, qual [, min_len, max_mm, max_window_mean_q])`

3' polyG-run trimming. Identifies a candidate polyG region at the 3' end
using fastp's mismatch budget (≤1 mismatch per 8 scanned bases, capped at
`max_mm`), then validates with a **quality gate**: if the mean Phred of the
candidate region exceeds `max_window_mean_q`, the trim is refused.

Defaults: `min_len=10`, `max_mm=5`, `max_window_mean_q=5`.

The quality gate is a deliberate divergence from fastp — it prevents
legitimate G-rich genomic regions (e.g. CpG islands, where Q≥30 is normal)
from being misclassified as NextSeq dark-cycle polyG (where Q≤2 is the
universal signature). To disable the gate and match fastp's pure sequence
behavior, pass `max_window_mean_q=93` (the max valid Phred).

### `trim_polyx(seq, qual [, min_len, max_mm])`

3' generic homopolymer trim. Identifies the dominant base across the
scanned 3' tail and trims from that base's leftmost position. Ties on
dominant base are broken in **ACGT order** (earliest wins) so results are
reproducible across runs and platforms — fastp's tie behavior depends on
input scan order.

Defaults: `min_len=10`, `max_mm=5`. No quality gate.

### `trim_adapters(seq, qual, adapter [, match_revcomp, min_match, allow_pre_start])`

3-phase adapter matching. Tries phases in order, returning the first phase
that finds a match (and within a phase, the leftmost match):

1. **Exact Hamming** with one mismatch per 8 compared bases.
2. **Hamming + 1 insertion** in seq (seq has one extra base; exhaustive
   search across all candidate insertion positions).
3. **Hamming + 1 deletion** in seq (seq missing one base; exhaustive).

`adapter` is either a single `VARCHAR` or a `LIST(VARCHAR)`. When a list,
all adapters (and their reverse complements if `match_revcomp=true`) are
tried; the **leftmost** trim point across all candidate matches wins —
i.e. the most aggressive trim. Exact-duplicate candidates (including RC pairs
that collapse to the same string) are deduplicated *after* `min_match` is
fixed — a matching-cost saving that never changes the result — and the search
stops early once a candidate matches at position 0 (the leftmost possible
trim, so no remaining candidate can beat it).

Defaults:
- `match_revcomp=false`. fastp expects `--adapter_sequence_r2` to be the
  pre-RC'd adapter; we offer the convenience flag instead.
- `min_match=0` → auto-scale (4 for 1 candidate, 5 for 2–3, 6 for ≥4).
  The count is taken **after** reverse-complement expansion but **before**
  dedup, so with `match_revcomp=true` the candidate set is doubled before the
  tier is chosen and enabling RC matching can bump `min_match` up a tier (e.g.
  2 adapters → 4 candidates → 6). Pass `min_match >= 1` to override. fastp
  auto-scales silently; we expose it.
- `allow_pre_start=false`. fastp's negative-offset start handling is on
  unconditionally to accommodate Illumina A-tailing; we made it opt-in
  because non-Illumina protocols don't need it and it can over-trim.
  When enabled and a pre-start match is found, `trim_start=0` is returned —
  the entire read should be dropped.

### `trim_adapters_pe(seq1, qual1, seq2, qual2 [, adapters, overlap_require, overlap_diff_limit, overlap_diff_percent_limit, match_revcomp, min_match, allow_pre_start])`

Paired-end adapter trimming by insert-size inference — fastp's primary PE
adapter mechanism (`OverlapAnalysis`). Reverse-complements R2, scans for the
offset that aligns it against R1, and — when the inferred insert is *shorter*
than the read length, so each mate reads through into adapter past the overlap
— trims both mates back to the insert.

Two overloads:

- **4-arg** `(seq1, qual1, seq2, qual2)` — overlap-only, fastp defaults
  (`overlap_len_require=30`, `overlap_diff_limit=5`,
  `overlap_diff_percent_limit=20`). No adapter-by-sequence fallback.
- **11-arg** `(…, adapters, overlap_require, overlap_diff_limit,
  overlap_diff_percent_limit, match_revcomp, min_match, allow_pre_start)` —
  overlap analysis first; if it does **not** trim and `adapters` (a
  `LIST(VARCHAR)`) is non-empty, fall back to per-mate adapter-by-sequence
  using the same matcher as `trim_adapters` (leftmost match wins, same
  `match_revcomp` / `min_match` / `allow_pre_start` semantics). This mirrors
  fastp's overlap step → by-sequence fallback. The adapter list and all tuning
  params must be **constant** expressions — they are resolved once at bind time;
  a column reference throws.

Returns a `STRUCT`:

| Field | Type | Meaning |
|---|---|---|
| `sequence1` | `VARCHAR` | Trimmed R1 sequence |
| `quality1` | `LIST(UTINYINT)` | Trimmed R1 quality, in lockstep with `sequence1` |
| `sequence2` | `VARCHAR` | Trimmed R2 sequence |
| `quality2` | `LIST(UTINYINT)` | Trimmed R2 quality, in lockstep with `sequence2` |
| `overlap_len` | `INTEGER` | Width of the detected overlap (`0` if none) |
| `adapter_trimmed` | `BOOLEAN` | Whether any adapter was trimmed (by overlap or by fallback) |
| `trimmed1_3p` | `UINTEGER` | Bases removed from the R1 3' end |
| `trimmed2_3p` | `UINTEGER` | Bases removed from the R2 3' end |

`overlap_len` is the overlap-analysis result and is **independent of the trim
source**: it is non-zero for a full-insert overlap at offset ≥ 0 (where
`adapter_trimmed` is `false`, because neither mate reads into adapter), and a
row trimmed by the 11-arg adapter-by-sequence fallback still reports whatever
overlap width analysis found. Overlap trimming is 3'-only on both mates, so
there is no `trimmed_5p`.

The overlap path is **byte-for-byte parity** with native fastp 1.3.3
(`test/sql/qc_trim_adapters_pe_parity.test`, a committed frozen golden — fastp
need not be installed to run it). See the divergence notes below for the
no-gap and `complete_compare_require=50` caveats, and for the fallback's
matcher-divergence on ambiguous cases.

### `filter_read(seq, qual [, min_length, max_length, qualified_q, max_unqualified_pct, max_n, min_avg_q])`

Per-read filter with metric reporting. Computes length, N count, low-qual
count, and mean quality in a single pass over `qual`, then applies thresholds
in fastp's documented precedence order (first failing check wins):

1. Empty seq → `length`
2. `low_qual_bases * 100 > max_unqualified_pct * length` → `quality`
3. `min_avg_q > 0 AND mean < min_avg_q` → `quality`
4. `n_bases > max_n` → `n_base`
5. `length < min_length` → `length`
6. `max_length > 0 AND length > max_length` → `too_long`

Defaults: `min_length=15`, `max_length=0` (off), `qualified_q=15`,
`max_unqualified_pct=40`, `max_n=5`, `min_avg_q=0` (off).

## Where we diverge from fastp

| Algorithm | Divergence | Why |
|---|---|---|
| `trim_polyg` | Optional quality gate (default Q5) | Prevents over-trimming legit GC-rich gene starts; NextSeq dark-cycle polyG is uniformly Q≤2 |
| `trim_polyx` | Deterministic ACGT tie-break | Reproducibility across runs and platforms |
| `trim_adapters` | Leftmost match wins (not highest-scoring) | The leftmost candidate is the most likely true adapter start; later high-scoring matches are typically genomic chatter that happens to align well by chance |
| `trim_adapters` | `allow_pre_start` opt-in default off | fastp's always-on negative-offset can over-trim non-Illumina reads |
| `trim_adapters` | `min_match` exposed | fastp's silent auto-scale is unintuitive when tuning |
| `trim_adapters` | Exhaustive indel-position search | fastp's greedy commit produces false negatives when a sequencing error precedes the true indel site |
| `trim_adapters_pe` | Overlap analysis is **no-gap only** | fastp's one-gap overlap path is deferred (backlogged); the no-gap path covers the overwhelming majority of real overlaps |
| `trim_adapters_pe` | 11-arg fallback uses miint's adapter matcher | The overlap path is byte-for-byte fastp parity; the by-sequence fallback matches fastp on clean exact-adapter cases but inherits `trim_adapters`' leftmost/exhaustive-indel behavior, so it can diverge on ambiguous mismatch/indel cases |
| `filter_read` | All metrics in return struct, single fail reason | Richer audit than fastp's enum return; `WHERE` clause can filter on individual metrics |
| All | Input validation at scalar boundary | fastp accepts whatever C-string lengths it's handed; we throw on `seq.length() != qual.length()` and on invalid IUPAC bases |
| All | `window_size` capped at 1..1000 | Matches fastp's documented range; prevents integer overflow in window-sum math |

## Known fastp behaviors we faithfully ported

These could be considered weaknesses in fastp but we ported them as-is to
maintain bit-for-bit parity in the common case:

- **`trim_polyg`'s mismatch counter is global to the entire 3'→5' scan**,
  not local to the polyG region. Non-G bases at the 5' end of the read
  (which are normal) consume the same budget as internal mismatches in the
  polyG run. In practice the `max_mm=5` cap is generous enough that this
  rarely produces visible defects, but it can cause spurious refusals to
  trim when the read prefix is itself mostly non-G. fastp v0.23+ behaves
  the same way.

- **`trim_adapters_pe` overlap acceptance is decided on only the first 50
  overlapping bases** (fastp's `complete_compare_require=50`), but for an
  accepted overlap longer than 50 the reported `diff` (mismatch count) is the
  *full* count over the whole overlap, which can exceed `overlap_diff_limit`.
  This is a fastp quirk, ported faithfully to keep byte-for-byte parity.

## Reconciling external QC results — `infer_trim`

Sometimes reads are quality-controlled by an **external** tool rather than the
scalars above. When that tool only *trims* read ends (no base edits) and may
*omit* reads entirely, `infer_trim` recovers the per-read 5'/3' trim coordinates
so you can persist two integers against the originals instead of storing a
second copy of every trimmed sequence.

`infer_trim(original_reads, qcd_reads)` is a **table macro** (not a scalar). It
LEFT JOINs the two relations on `sequence_index` and locates the contiguous
QC'd sequence inside the original. Both relations must expose:

| Column | Type | Role |
|---|---|---|
| `sequence_index` | `BIGINT` | Join key |
| `sequence` | `VARCHAR` | Read sequence |

It returns one row per original read:

| Column | Type | Meaning |
|---|---|---|
| `sequence_index` | `BIGINT` | Join key, passed through |
| `trimmed_5p` | `UINTEGER` | Bases removed from the 5' end; `NULL` if the read was omitted by QC |
| `trimmed_3p` | `UINTEGER` | Bases removed from the 3' end; `NULL` if the read was omitted by QC |

> **Getting the join key right.** `sequence_index` must be **globally unique**
> per read and must identify the *same* read on both sides — the external tool
> has to carry the key through and emit it alongside each surviving read. Two
> traps with `read_fastx` as the source: it assigns `sequence_index`
> **positionally and resets it to 1 per input file**, so it is a safe key only
> for a *single* file; and you **cannot** recover a matching key by re-running
> `read_fastx` on the QC'd output — QC drops reads, so a fresh positional
> numbering no longer lines up with the original. A duplicated `sequence_index`
> on the `qcd_reads` side fans the join out (more than one row per original);
> the macro treats uniqueness as a precondition and does not police it. If keys
> are misaligned `infer_trim` usually fails loud on the resulting non-substring,
> but a coincidental substring match would pass silently, so get the key right.

```sql
-- Original reads from a SINGLE file; sequence_index is read_fastx's positional id.
CREATE VIEW original AS
    SELECT sequence_index, sequence1 AS sequence FROM read_fastx('reads.fastq.gz');

-- qcd_reads: the external tool's surviving reads, each carrying the ORIGINAL
-- sequence_index (project/rename so both sides share `sequence_index` + `sequence`).
SELECT * FROM infer_trim(original, qcd_reads) ORDER BY sequence_index;
```

Semantics:

- **Omitted reads → `NULL` / `NULL`.** A read with no row in `qcd_reads`
  (dropped by QC) yields `NULL` coordinates, straight from the LEFT JOIN — the
  "any sequence omitted is null" case. A `qcd_reads` row whose `sequence_index`
  is absent from `original_reads` is dropped (one row per original, by design).
- **Fails loud, whatever you `SELECT`.** The integrity checks live in a
  row-level `WHERE` filter the optimizer cannot prune, so they fire even for
  `SELECT count(*)` or a single-column projection — not only for `SELECT *`.
  `infer_trim` throws, naming the offending `sequence_index`, when a QC'd
  sequence is **not a contiguous substring** of its original (the external tool
  edited bases, or the join key is mismatched), or when a *present* QC row has a
  **NULL or empty** sequence, or the **original** sequence is `NULL`. It is
  strictly for pure end-trimming; represent a dropped read by omitting its row
  (→ `NULL` / `NULL`), not by an empty sequence. Base-modifying QC (quality
  masking, error correction, reverse-complementing) needs alignment, not a
  substring search.
- **Leftmost match wins.** If the kept block occurs at more than one offset in a
  self-repetitive read, the leftmost occurrence is chosen — the true offset
  cannot be recovered from sequence alone, so the most conservative 5' trim is
  taken.
- **Persist coordinates, reconstruct on demand.** Store only `trimmed_5p` /
  `trimmed_3p` next to the originals; regenerate the trimmed read when needed
  with `substr(sequence, trimmed_5p + 1, length(sequence) - trimmed_5p - trimmed_3p)`.

Performance: a parallel hash join plus one `position()` substring search per
matched read, so cost is linear in total sequence bytes (≈ Σ read length) with
the join scaling in row count; both stages are fully parallel.

## Related — long-read amplicon / UMI primitives

The QC trimmers above target short-read fastp-style cleanup. For
long-read amplicon work (UMI binning, MSA consensus, per-base variant
positions) miint exposes a complementary set of primitives that compose
into a Karst-protocol UMI consensus pipeline:

- [`extract_linked_amplicon`](utilities.md#extracting-linked-amplicons) — cut out the interior between two
  flanking adapters (cutadapt `-g X...Y` equivalent), WFA2-powered.
- [`match_short_barcodes`](alignment_analysis.md#barcode-matching) — Hamming-distance matcher for
  fixed-length barcodes (UMIs, sample indices).
- [`compute_pileup`](alignment_analysis.md#per-base-pileup) — per-base CIGAR walker that emits
  `(read_id, ref_pos, ref_base, query_base, query_qual)` rows.
- [`compute_msa_consensus`](alignment_analysis.md#msa-column-consensus) — Q-aware MSA column consensus
  with HP post-correction (replaces Racon polishing for HiFi).

## Backlog (deferred from v1)

- DuckDB named-parameter overloads to avoid the positional `(seq, qual, X,
  Y, Z, ...)` calling style for the 6+-arg forms.
- Test coverage for the `(i+1)/8` mismatch-budget truncation edge cases.
- Test for `trim_polyx` with embedded N bases in the 3' tail.
- SIMD-vectorized `filter_read` metric pass (currently scalar; profile
  before optimizing).
- Aho-Corasick multi-adapter matcher for very large adapter lists
  (50+ adapters) — current implementation is O(candidates × read_len ×
  adapter_len) after exact-duplicate dedup and the position-0 early-exit.
  Collapsing prefix/substring-redundant candidates is also still backlogged:
  it is unsafe under the per-length mismatch budget (a shorter candidate can
  match where its longer superstring would not), so only *exact* duplicates
  are removed.
- Integration test on a real-world FASTQ at ≥1M reads to verify end-to-end
  throughput and numerical stability. The current `qc_pipeline.test`
  fixture (30 reads) gives bit-for-bit ground-truth verification at toy
  scale; scale-validation is a separate concern.
