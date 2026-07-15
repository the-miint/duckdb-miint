# Pairwise alignment

Line up two sequences against each other, base by base, to measure how similar they are and exactly where they differ. Use this when you have a query sequence and a reference (or two reads) and want a similarity score, a compact description of the matches/mismatches/gaps (a CIGAR string), or the two sequences written out with gap characters so you can eyeball the alignment.

## Table of Contents

- [Choosing a family](#choosing-a-family) - Which algorithm to pick, and how their scores differ.
- [Detail levels: score, cigar, full](#detail-levels-score-cigar-full) - The three output shapes every family shares.
- [WFA2 (Wavefront)](#wfa2-wavefront) - `align_pairwise_wfa2_*`: extended CIGAR, penalty-style score.
- [KSW2 extz (standard affine)](#ksw2-extz-standard-affine) - `align_pairwise_ksw2_*`: general DNA alignment with optional banding / z-drop.
- [KSW2 extd (dual affine)](#ksw2-extd-dual-affine) - `align_pairwise_ksw2_dual_affine_*`: long-read alignment with a second gap model for long indels.
- [KSW2 exts (splice-aware)](#ksw2-exts-splice-aware) - `align_pairwise_ksw2_splice_*`: intron-skipping for RNA-seq (experimental).

Pairwise sequence alignment is exposed through four families, each backed by a different algorithm or KSW2 mode. All families share three detail levels with consistent return shapes (see [Detail levels](#detail-levels-score-cigar-full)).

NULL inputs produce NULL output. Alignment failure (e.g., z-drop early termination, excessive divergence) also produces NULL. Penalty parameters must be constant values, not column references.

## Choosing a family

| Family | Backend | Score semantic | CIGAR ops | Pick when |
|---|---|---|---|---|
| `align_pairwise_wfa2_*` | WFA2-lib (Wavefront) | Penalty: `0` = identical, larger = more divergent | Extended (`=` / `X`) | Short to long DNA; exact gap-affine global alignment |
| `align_pairwise_ksw2_*` | KSW2 `ksw_extz2_sse` (SIMD banded DP) | Native positive: identical = `qlen * match` | Extended (`=` / `X`, plus `I` / `D`) | General DNA alignment; optional bandwidth / z-drop tuning |
| `align_pairwise_ksw2_dual_affine_*` | KSW2 `ksw_extd2_sse` | Native positive | `=`, `X`, `I`, `D` | Long-read alignment; long indels amortize over a second affine pair |
| `align_pairwise_ksw2_splice_*` | KSW2 `ksw_exts2_sse` | Native positive | `=`, `X`, `I`, `D`, `N` (intron skip) | Splice-aware (RNA-seq); intron-open penalty + non-canonical-boundary penalty |

All families emit **extended** CIGAR (`=` for match, `X` for mismatch). KSW2 natively produces only `M`; the `_cigar` and `_full` outputs run an eqx post-pass that splits each `M` into `=` / `X` by comparing the aligned bases (case-insensitive, so soft-masked lowercase bases are not counted as mismatches). This means sequence identity can be read directly off any family's CIGAR with `cigar_sequence_identity`.

WFA2 scores and KSW2 scores are on different scales -- WFA2 is penalty-style (lower is better, identical = 0), KSW2 is additive (higher is better, positive contributions from matches). Do not compare scores across families.

## Detail levels: score, cigar, full

Every family comes in three variants, distinguished by a suffix on the function name:

- `*_score` -- alignment score only (fastest). Returns INTEGER.
- `*_cigar` -- score + CIGAR string, as `STRUCT(score INTEGER, cigar VARCHAR)`.
- `*_full` -- score + CIGAR + aligned sequences with `-` gap characters, as `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)`.

## WFA2 (Wavefront)

Gap-affine pairwise sequence alignment powered by [WFA2-lib](https://github.com/smarco/WFA2-lib). Scores are penalty-style: `0` means identical, and larger values mean more divergent. The CIGAR uses the extended alphabet (`=` for a match, `X` for a mismatch), so match and mismatch are distinguishable.

Functions: `align_pairwise_wfa2_score`, `align_pairwise_wfa2_cigar`, `align_pairwise_wfa2_full`.

**Function signatures:**

```sql
-- 2-arg: default penalties (mismatch=4, gap_open=6, gap_extend=2)
SELECT align_pairwise_wfa2_score(query, subject);

-- 5-arg: custom penalties
SELECT align_pairwise_wfa2_score(query, subject, 2, 6, 2);
```

**Parameters (5-arg form):**
- `query` (VARCHAR), `subject` (VARCHAR): the two sequences to align
- `mismatch` (INTEGER): mismatch penalty (must be > 0)
- `gap_open` (INTEGER): gap-opening penalty (must be >= 0)
- `gap_extend` (INTEGER): gap-extension penalty (must be > 0)

**Return type:**
- `align_pairwise_wfa2_score` returns INTEGER (the alignment penalty)
- `align_pairwise_wfa2_cigar` returns `STRUCT(score INTEGER, cigar VARCHAR)`
- `align_pairwise_wfa2_full` returns `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)`

**Behavior:**
- Score is penalty-style: `0` = identical, larger = more divergent.
- CIGAR uses the extended alphabet (`=` for match, `X` for mismatch).
- NULL inputs produce NULL output; alignment failure produces NULL.
- Penalty parameters must be constant values, not column references.

**Examples:**
```sql
SELECT align_pairwise_wfa2_score('ACGT', 'ACGT');           -- 0  (identical)
SELECT align_pairwise_wfa2_score('ACGT', 'ACAT');           -- 4  (one mismatch)
SELECT (align_pairwise_wfa2_cigar('ACGT', 'ACAT')).cigar;   -- 2=1X1=
SELECT (align_pairwise_wfa2_full('ACGT', 'AGT')).query_aligned,
       (align_pairwise_wfa2_full('ACGT', 'AGT')).subject_aligned;
-- aligned sequences with '-' for the indel
```

**Error conditions:**
- `mismatch` not > 0
- `gap_open` not >= 0
- `gap_extend` not > 0
- Penalty parameters passed as column references rather than constants

## KSW2 extz (standard affine)

Standard affine extension alignment via `ksw_extz2_sse` (bundled inside [minimap2](https://github.com/lh3/minimap2)). Scores are native positive: an identical alignment scores `qlen * match`, and mismatches/gaps subtract. The `_cigar` / `_full` outputs use the extended alphabet (`=` match, `X` mismatch): KSW2 natively emits `M`, and an eqx post-pass splits each `M` into `=` / `X`.

Functions: `align_pairwise_ksw2_score`, `align_pairwise_ksw2_cigar`, `align_pairwise_ksw2_full`.

**Function signatures:**

```sql
-- 2-arg: defaults (match=2, mismatch=4, gap_open=6, gap_extend=2; w=-1, zdrop=-1)
SELECT align_pairwise_ksw2_score(query, subject);

-- 6-arg: custom penalties (bandwidth + z-drop default to -1, disabled)
SELECT align_pairwise_ksw2_score(query, subject, 2, 4, 6, 2);

-- 8-arg: advanced (explicit bandwidth and z-drop)
SELECT align_pairwise_ksw2_score(query, subject, 2, 4, 6, 2, 100, 400);
```

**Parameters (6-/8-arg forms):**
- `query` (VARCHAR), `subject` (VARCHAR): the two sequences to align
- `match` (INTEGER): match score (must be > 0); positive contribution per matched base
- `mismatch` (INTEGER): mismatch penalty (must be > 0); subtracted per mismatched base
- `gap_open` (INTEGER): gap-opening penalty (must be >= 0)
- `gap_extend` (INTEGER): gap-extension penalty (must be > 0)
- `w` (INTEGER, 8-arg only): bandwidth; negative disables banding (full DP)
- `zdrop` (INTEGER, 8-arg only): z-drop threshold for early termination; negative disables

**Return type:**
- `align_pairwise_ksw2_score` returns INTEGER (the alignment score)
- `align_pairwise_ksw2_cigar` returns `STRUCT(score INTEGER, cigar VARCHAR)`
- `align_pairwise_ksw2_full` returns `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)`

**Behavior:**
- Score is native positive: an identical alignment scores `qlen * match`.
- CIGAR uses the extended alphabet (`=` / `X`); KSW2's native `M` is split into `=` / `X` by an eqx post-pass.
- NULL inputs produce NULL output; alignment failure (e.g., z-drop early termination) produces NULL.
- Penalty parameters must be constant values, not column references.

**Examples:**
```sql
SELECT align_pairwise_ksw2_score('ACGT', 'ACGT');           -- 8  (4 * match=2)
SELECT align_pairwise_ksw2_score('ACGT', 'ACAT');           -- 2  (3*2 - 4)
SELECT (align_pairwise_ksw2_cigar('ACGT', 'ACAT')).cigar;   -- 2=1X1= (eqx post-pass splits M into =/X)
```

**Error conditions:**
- `match` not > 0
- `mismatch` not > 0
- `gap_open` not >= 0
- `gap_extend` not > 0
- Penalty parameters passed as column references rather than constants

## KSW2 extd (dual affine)

Dual affine gap penalties via `ksw_extd2_sse`. For each gap of length `L`, KSW2 picks the cheaper of `gap_open + L*gap_extend` and `gap_open2 + L*gap_extend2`. Typical configuration: first pair has cheap open + moderate extend (short indels), second pair has expensive open + cheap extend (long indels / structural variants). This makes it well suited to long-read alignment, where long indels amortize over the second affine pair.

Functions: `align_pairwise_ksw2_dual_affine_score`, `align_pairwise_ksw2_dual_affine_cigar`, `align_pairwise_ksw2_dual_affine_full`.

**Function signatures:**

```sql
-- 2-arg: defaults (match=2, mismatch=4, gap_open=6, gap_extend=2, gap_open2=24, gap_extend2=1)
SELECT align_pairwise_ksw2_dual_affine_score(query, subject);

-- 8-arg: custom penalties
SELECT align_pairwise_ksw2_dual_affine_score(query, subject, 2, 4, 6, 2, 24, 1);

-- 10-arg: advanced (with bandwidth and z-drop)
SELECT align_pairwise_ksw2_dual_affine_score(query, subject, 2, 4, 6, 2, 24, 1, -1, -1);
```

**Parameters (8-/10-arg forms):** same first four as the extz family (`match`, `mismatch`, `gap_open`, `gap_extend`), plus:
- `gap_open2` (INTEGER): second-pair gap-opening penalty (must be >= 0); typically larger than `gap_open`
- `gap_extend2` (INTEGER): second-pair gap-extension penalty (must be > 0); typically smaller than `gap_extend`
- `w`, `zdrop` (10-arg only): bandwidth and z-drop as in the extz family

**Return type:**
- `align_pairwise_ksw2_dual_affine_score` returns INTEGER (the alignment score)
- `align_pairwise_ksw2_dual_affine_cigar` returns `STRUCT(score INTEGER, cigar VARCHAR)`
- `align_pairwise_ksw2_dual_affine_full` returns `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)`

**Behavior:**
- Score is native positive (same scale as the extz family).
- CIGAR ops are `=`, `X`, `I`, `D` (eqx post-pass splits KSW2's native `M`).
- For each gap of length `L`, the cheaper of `gap_open + L*gap_extend` and `gap_open2 + L*gap_extend2` is chosen.
- NULL inputs produce NULL output; alignment failure produces NULL.
- Penalty parameters must be constant values, not column references.

**Example -- long gap uses the second affine:**
```sql
-- 20 matched bases + 30-bp insertion
-- First-affine gap cost: 6 + 30*2 = 66.  Second-affine: 24 + 30*1 = 54.  Dual picks the cheaper.
SELECT align_pairwise_ksw2_dual_affine_score(repeat('A', 50), repeat('A', 20));
-- -14   (vs. -26 from align_pairwise_ksw2_score with the same penalties)
```

**Error conditions:**
- `match` not > 0
- `mismatch` not > 0
- `gap_open` not >= 0
- `gap_extend` not > 0
- `gap_open2` not >= 0
- `gap_extend2` not > 0
- Penalty parameters passed as column references rather than constants

## KSW2 exts (splice-aware)

Splice-aware extension alignment via `ksw_exts2_sse`. Intended for RNA-seq alignment where introns are large gaps and canonical (GT-AG) splice sites are preferred. Junction guidance (`junc` array) is not exposed in v1; alignment relies on the score model alone. Forward-strand splice flag (`KSW_EZ_SPLICE_FOR`) is fixed in v1.

Functions: `align_pairwise_ksw2_splice_score`, `align_pairwise_ksw2_splice_cigar`, `align_pairwise_ksw2_splice_full`.

> ⚠️ **v1 status: experimental.** Without junction guidance (`junc` array), splice-site detection relies entirely on the score model and is unreliable for real RNA-seq workloads — short test inputs will not exercise intron skipping at all, and long inputs may produce plausible-looking but incorrect intron boundaries. Use this family today for verifying the score model and integrating with custom splice-aware pipelines that supply external junction calls; **do not use it as a drop-in replacement for a real splice-aware aligner** until a future version exposes the `junc` array and per-strand control.

**Function signatures:**

```sql
-- 2-arg: defaults (minimap2 --splice preset shape: match=2, mismatch=4, gap_open=6,
--                  gap_extend=2, gap_open2=24, noncan=9)
SELECT align_pairwise_ksw2_splice_score(query, subject);

-- 8-arg: custom penalties
SELECT align_pairwise_ksw2_splice_score(query, subject, 2, 4, 6, 2, 24, 9);

-- 9-arg: advanced (with z-drop; no bandwidth parameter -- ksw_exts2_sse has none)
SELECT align_pairwise_ksw2_splice_score(query, subject, 2, 4, 6, 2, 24, 9, -1);
```

**Parameters (8-/9-arg forms):** same first four as the extz family (`match`, `mismatch`, `gap_open`, `gap_extend`), plus:
- `gap_open2` (INTEGER): intron-open penalty (must be >= 0); introns extend at the `gap_extend` rate
- `noncan` (INTEGER): penalty added when the chosen intron boundary is non-canonical (must be >= 0)
- `zdrop` (9-arg only): z-drop threshold

**Return type:**
- `align_pairwise_ksw2_splice_score` returns INTEGER (the alignment score)
- `align_pairwise_ksw2_splice_cigar` returns `STRUCT(score INTEGER, cigar VARCHAR)`
- `align_pairwise_ksw2_splice_full` returns `STRUCT(score INTEGER, cigar VARCHAR, query_aligned VARCHAR, subject_aligned VARCHAR)`

**Behavior:**
- Score is native positive (same scale as the extz family).
- CIGAR ops are `=`, `X`, `I`, `D`, and `N` (intron skip); the eqx post-pass splits `M` into `=` / `X`.
- The `_full` aligned-sequence output renders `N` (intron-skip) the same as `D` -- the intron appears as gap characters in `query_aligned`, with the corresponding subject bases in `subject_aligned`. The CIGAR string preserves the aligned-op (`=` / `X`) vs `N` distinction for downstream consumers.
- Forward-strand splice flag (`KSW_EZ_SPLICE_FOR`) is fixed in v1; no junction guidance is available.
- NULL inputs produce NULL output; alignment failure produces NULL.
- Penalty parameters must be constant values, not column references.

**Example -- splice mode behaves like extz for short input without intron-boundary patterns:**
```sql
SELECT align_pairwise_ksw2_splice_score('ACGT', 'ACAT');   -- 2  (3*match - mismatch)
```

**Error conditions:**
- `match` not > 0
- `mismatch` not > 0
- `gap_open` not >= 0
- `gap_extend` not > 0
- `gap_open2` not >= 0
- `noncan` not >= 0
- Penalty parameters passed as column references rather than constants

---

See also: [reading sequence files](reading.md) for loading the query/subject sequences, [multiple sequence alignment](alignment_multiple.md), and [alignment analysis](alignment_analysis.md).
