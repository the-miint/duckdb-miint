# Quality Control (QC) Scalar Functions

A family of scalar functions for FASTQ read quality control: adapter trimming,
polyG/polyX tail trimming, sliding-window quality trimming, and per-read
length/quality/N-base filtering. The algorithms are ported (not vendored)
from [fastp](https://github.com/OpenGene/fastp) — see
[`THIRD_PARTY_LICENSES.md`](../THIRD_PARTY_LICENSES.md) for attribution.

These functions compose in SQL to express the canonical fastp QC pipeline
over any source of sequence + quality data — typically [`read_fastx`](table-functions.md#read_fastx).

## Function catalog

| Function | Returns | Purpose |
|---|---|---|
| `trim_quality_5p(seq, qual [, window, mean_q])` | trim struct | fastp `cut_front` — drop low-quality 5' bases via sliding window |
| `trim_quality_3p(seq, qual [, window, mean_q])` | trim struct | fastp `cut_tail` — drop low-quality 3' bases via sliding window |
| `trim_quality_sliding(seq, qual [, window, mean_q])` | trim struct | fastp `cut_right` — drop window-and-everything-right at first bad window |
| `trim_polyg(seq, qual [, min_len, max_mm, max_window_mean_q])` | trim struct | 3' polyG run trim (NextSeq dark-cycle cleanup), optionally quality-gated |
| `trim_polyx(seq, qual [, min_len, max_mm])` | trim struct | 3' generic homopolymer trim (any most-frequent base) |
| `trim_adapters(seq, qual, adapter [, match_revcomp, min_match, allow_pre_start])` | trim struct | 3-phase adapter match (Hamming + 1-insert + 1-delete); `adapter` is `VARCHAR` or `LIST(VARCHAR)` |
| `filter_read(seq, qual [, min_length, max_length, qualified_q, max_unqualified_pct, max_n, min_avg_q])` | filter struct | Per-read pass/fail with metrics |
| `qc_version()` | `VARCHAR` | Version string |

All scalars accept `seq VARCHAR` + `qual LIST(UTINYINT)`. Quality is the
decoded Phred integer list emitted by `read_fastx` (not ASCII), so values
are 0..93 directly.

## Return structs

Every trim function returns a `STRUCT` with these fields:

| Field | Type | Meaning |
|---|---|---|
| `sequence` | `VARCHAR` | The trimmed sequence (a slice of the input) |
| `quality` | `LIST(UTINYINT)` | The trimmed quality, in lockstep with `sequence` |
| `trimmed_5p` | `UINTEGER` | Number of bases removed from the 5' end |
| `trimmed_3p` | `UINTEGER` | Number of bases removed from the 3' end |

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

The fastp-recommended order is: adapter trim → polyG → quality trim → filter.
In SQL, chain the scalars through aliases in a single `WITH` and access
struct fields via `(alias).field`:

```sql
WITH pipeline AS (
    SELECT
        sequence_index,
        read_id,
        trim_adapters(sequence1, qual1, 'AGATCGGAAGAGC')        AS t1,
        trim_polyg((t1).sequence, (t1).quality)                 AS t2,
        trim_quality_3p((t2).sequence, (t2).quality, 4, 20)     AS t3,
        filter_read((t3).sequence, (t3).quality)                AS f
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

Each scalar is invoked exactly once per row. DuckDB's CSE optimizer
deduplicates repeated struct extractions automatically, so accessing
`(t1).sequence` and `(t1).quality` and `(t1).trimmed_3p` does not
re-evaluate `trim_adapters`.

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
i.e. the most aggressive trim.

Defaults:
- `match_revcomp=false`. fastp expects `--adapter_sequence_r2` to be the
  pre-RC'd adapter; we offer the convenience flag instead.
- `min_match=0` → auto-scale (4 for 1 candidate, 5 for 2–3, 6 for ≥4).
  Pass `min_match >= 1` to override. fastp auto-scales silently; we expose
  it.
- `allow_pre_start=false`. fastp's negative-offset start handling is on
  unconditionally to accommodate Illumina A-tailing; we made it opt-in
  because non-Illumina protocols don't need it and it can over-trim.
  When enabled and a pre-start match is found, `trim_start=0` is returned —
  the entire read should be dropped.

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
| `trim_adapters` | Best-match-wins reverted; leftmost-wins (fastp) is kept | Leftmost match is biologically the true adapter start; later matches are typically genomic chatter |
| `trim_adapters` | `allow_pre_start` opt-in default off | fastp's always-on negative-offset can over-trim non-Illumina reads |
| `trim_adapters` | `min_match` exposed | fastp's silent auto-scale is unintuitive when tuning |
| `trim_adapters` | Exhaustive indel-position search | fastp's greedy commit produces false negatives when a sequencing error precedes the true indel site |
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

## Backlog (deferred from v1)

- DuckDB named-parameter overloads to avoid the positional `(seq, qual, X,
  Y, Z, ...)` calling style for the 6+-arg forms.
- Test coverage for the `(i+1)/8` mismatch-budget truncation edge cases.
- Test for `trim_polyx` with embedded N bases in the 3' tail.
- SIMD-vectorized `filter_read` metric pass (currently scalar; profile
  before optimizing).
- Aho-Corasick multi-adapter matcher for very large adapter lists
  (50+ adapters) — current implementation is O(adapters × read_len × adapter_len).
