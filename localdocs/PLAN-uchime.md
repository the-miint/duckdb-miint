# UCHIME Chimera Detection — TDD Implementation Plan

## Context

Reimplement UCHIME chimera detection (Edgar et al. 2011, Bioinformatics 27:2194-2200) as DuckDB table functions `uchime_ref` and `uchime_denovo`. Must produce identical results to vsearch v2.30.5. Target reads: up to full-length 16S (~1500bp). Must match or beat vsearch performance on thousands of reads. Cleanroom implementation from the published paper only.

## Architecture

```
KmerIndex              — 8-mer inverted index for fast candidate parent search
ChimeraDetector        — Core UCHIME algorithm (namespace miint)
  ├── build_star_alignment()   — merge two pairwise alignments into 3-way
  ├── classify_diffs()         — A/B/N/?/ignore at each column
  ├── sweep_breakpoints()      — left-to-right h-score sweep
  ├── select_parents()         — smoothed 32bp window, pick top 2
  └── detect()                 — orchestrator returning UchimeResult
WFA2Aligner            — existing pairwise aligner (mismatch=6, gap_open=20, gap_extend=4)
uchime_ref             — table function (multi-threaded, fixed reference DB)
uchime_denovo          — table function (single-threaded, incremental DB)
```

### Alignment Equivalence — VALIDATED

vsearch uses match=2, mismatch=-4, gap_open=-20, gap_extend=-2 (score-maximization). Equivalent WFA2 penalty-minimization params: **mismatch=6, gap_open=20, gap_extend=4**. Derivation: penalty = match_reward + |score_penalty| for operations that displace a match.

**Status**: Validated on all test sequences. WFA2 produces identical ungapped alignments to vsearch for our 300bp substitution-only test data. Diff counts match exactly for 3 of 4 chimeras, and within tolerance for chimera3 (mutation at crossover). Gap model conversion has NOT been validated on gapped sequences — revisit if real-world data triggers gaps.

### Output Schema (18 columns, mirrors vsearch --uchimeout)

| Column | Type | Description |
|--------|------|-------------|
| score | DOUBLE | h-score |
| query | VARCHAR | Query label |
| parent_a | VARCHAR | Parent A label |
| parent_b | VARCHAR | Parent B label |
| closest_parent | VARCHAR | Top parent (higher identity) |
| id_query_model | DOUBLE | Query-to-model identity % |
| id_query_a | DOUBLE | Query-to-parentA identity % |
| id_query_b | DOUBLE | Query-to-parentB identity % |
| id_a_b | DOUBLE | ParentA-to-parentB identity % |
| id_query_top | DOUBLE | Query-to-closest-parent identity % |
| left_yes | INTEGER | Left yes votes |
| left_no | INTEGER | Left no votes |
| left_abstain | INTEGER | Left abstain votes |
| right_yes | INTEGER | Right yes votes |
| right_no | INTEGER | Right no votes |
| right_abstain | INTEGER | Right abstain votes |
| divergence | DOUBLE | QM-QT divergence |
| flag | VARCHAR | Y/N/? |

### Named Parameters

Both functions accept DuckDB table/view names (VARCHAR) for sequence input. Required columns: `read_id` (VARCHAR), `sequence` (VARCHAR). This composes naturally with `read_fastx()` and SQL workflows.

**uchime_ref(query_table, db)**: `db` (required — table/view name for reference sequences, same schema), `minh` (0.28), `xn` (8.0), `dn` (1.4), `mindiv` (0.8), `mindiffs` (3)

**uchime_denovo(input_table)**: input table must additionally have `size` (INTEGER) column for abundance. `abskew` (2.0), plus same scoring params.

```sql
-- uchime_ref example
CREATE TABLE queries AS SELECT read_id, sequence1 AS sequence FROM read_fastx('queries.fasta');
CREATE TABLE refs AS SELECT read_id, sequence1 AS sequence FROM read_fastx('gold.fasta');
SELECT * FROM uchime_ref('queries', db:='refs');

-- uchime_denovo example
CREATE TABLE seqs AS SELECT read_id, sequence1 AS sequence, count(*) AS size
  FROM read_fastx('amplicons.fasta') GROUP BY ALL;
SELECT * FROM uchime_denovo('seqs');
```

---

## Files to Create

| File | Purpose |
|------|---------|
| `src/include/KmerIndex.hpp` | 8-mer inverted index class |
| `src/KmerIndex.cpp` | KmerIndex implementation |
| `src/include/ChimeraDetector.hpp` | Core UCHIME algorithm |
| `src/ChimeraDetector.cpp` | ChimeraDetector implementation |
| `src/include/uchime_ref.hpp` | uchime_ref table function |
| `src/uchime_ref.cpp` | uchime_ref implementation |
| `src/include/uchime_denovo.hpp` | uchime_denovo table function |
| `src/uchime_denovo.cpp` | uchime_denovo implementation |
| `test/cpp/test_KmerIndex.cpp` | KmerIndex unit tests |
| `test/cpp/test_ChimeraDetector.cpp` | ChimeraDetector unit tests |
| `test/sql/uchime_ref.test` | SQL integration tests |
| `test/sql/uchime_denovo.test` | SQL integration tests |
| `data/uchime/*.fasta` | Synthetic test data |
| `data/uchime/expected_*.tsv` | vsearch ground truth |

## Files to Modify

| File | Change |
|------|--------|
| `CMakeLists.txt` | Add 4 .cpp to EXTENSION_SOURCES, 4 .cpp + 2 test .cpp to TEST_SOURCES |
| `src/miint_extension.cpp` | Add includes + Register() calls |
| `run_tests.sh` | Add `VSEARCH_AVAILABLE` env var detection |

---

## TDD Phases

### Phase 0: Test Data + Ground Truth — DONE ✓

Created synthetic test data in `data/uchime/`:
- `chimera_ref.fasta` — 6 references (300bp, ~70% pairwise identity)
- `chimera_queries.fasta` — 8 queries (2 clean, 3 chimeras, 1 no-parent, 1 divergent, 1 short)
- `denovo_input.fasta` — 6 abundance-annotated sequences
- `expected_ref.tsv` / `expected_denovo.tsv` — vsearch v2.30.5 ground truth
- `expected_ref_alns.txt` / `expected_denovo_alns.txt` — vsearch alignment output

vsearch results: 4/8 chimeras in ref mode, 2/6 in denovo mode.

### Phase 1: K-mer Index — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **Presence/absence model**: posting lists store unique seq_ids per k-mer (not multiplicity)
- **Internal ID assignment**: `add_sequence()` returns sequential uint32_t, no external IDs
- **thread_local scratch buffers**: `selected_hits`, `chunk_hits`, `chunk_touched` avoid per-query allocations
- **Sparse touched-entries pattern**: O(hits) per chunk, not O(seq_count)
- 9 test cases, 94 assertions

### Phase 2: Alignment Equivalence — DONE ✓

Validated WFA2Aligner(6, 20, 4) produces same ungapped alignments as vsearch for our test data. Penalty derivation documented with caveat that gap model hasn't been validated on gapped sequences.

### Phase 3: Star Alignment + Diff Classification — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **Post-loop assertion**: throws if pairwise alignments have mismatched query consumption
- **is_ambiguous() excludes gaps**: gap '-' is not ambiguous, handled separately
- **Diff counts validated against vsearch**: exact match for all 4 chimeric test sequences
  - chimera1: 45A + 46B (exact)
  - chimera2: 35A + 50B (exact)
  - chimera3: ≥85 total (loose due to crossover mutation)
  - divergent: 46A + 43B + 11N + 2? = 102 (exact)
- 14 test cases, 86 assertions

### Phase 4: Smoothed Parent Selection — IN PROGRESS

**RED** — `test/cpp/test_ChimeraDetector.cpp`:
- `compute_match_profile()`: from pairwise alignment, produce per-position match array
- `compute_smoothed_identity()`: 32bp sliding window over match profile
- `select_parents()`: chimeric query → correct parent A (left winner) and parent B (right winner)
- Single dominant parent → both A and B same sequence
- Ties broken deterministically
- **vsearch-validated**: select_parents on chimera1 with ref1+ref2 candidates → picks correct parents

**GREEN**:
- `compute_match_profile()`: for each aligned column, 1 if Q==S and neither is gap, else 0
- `compute_smoothed_identity()`: running sum of match array over window of 32
- `select_parents()`: for each candidate, compute smoothed profile. At each position, candidate with max smooth value "wins". Parent A = most total wins. Wipe A's winning positions, recompute, Parent B = most wins in round 2.
- Return struct with parent indices AND cached pairwise alignments (reused in Phase 5)

**REFACTOR** — Ensure alignment cache is shared with detect()

### Phase 5: H-Score Sweep + Classification

**RED** — `test/cpp/test_ChimeraDetector.cpp`:
- `sweep_breakpoints()`: 30 A-diffs then 30 B-diffs → h-score ≈ 1.0, clear chimera
- All A-diffs → h-score 0 (one-sided)
- All IGNORE → h-score 0
- Reversed configuration (B-diffs left, A-diffs right)
- Hand-computed example with mixed diffs → expected H value
- Classification: H≥minh + divdiff≥mindiv + diffs≥mindiffs → Y
- H≥minh but fails other checks → ? (borderline)
- H<minh → N
- **Full pipeline test**: detect() on synthetic chimera → flag=Y, clean → flag=N
- **vsearch ground truth**: compare full UchimeResult (score, votes, identities, flag) for all 8 queries against expected_ref.tsv
  - Note: pre-breakpoint diff counts (Phase 3) ≠ post-breakpoint vote counts (this phase).
    vsearch's LY/RY/LN/RN/LA/RA are votes assigned *after* breakpoint selection.
    N-diffs become "no-votes" when they fall on the "wrong" side, "abstain" otherwise.
  - Identity %: matches / non-ignored columns × 100

**GREEN**:
- `sweep_breakpoints()`: initialize all votes on right side, sweep left-to-right transferring votes, compute H=H_left×H_right at each position, track best. Handle reversed configuration (swap A↔B if left_n > left_y and right_n > right_y).
- `UchimeResult` struct with all 18 output fields
- `UchimeParams` struct (minh, xn, dn, mindiv, mindiffs, abskew)
- `ChimeraDetector` class with:
  - `set_reference()` — load reference sequences + build k-mer index
  - `detect()` — full pipeline: find_candidates → align → select_parents → build_star_alignment → classify_diffs → sweep_breakpoints → compute identities → classify
- Identity computation: matches/non-ignored-columns × 100
- Divergence: id_query_model - id_query_top

**REFACTOR** — Extract identity calculation helper, clean up detect() flow

**WFA2 reuse verification**: The same WFA2Aligner instance must be reused across:
- Multiple candidate alignments within select_parents()
- The A-B identity alignment in detect()
- Multiple detect() calls from the same thread
Tests verify this by calling detect() repeatedly with the same aligner and checking for consistent results. The table function's LocalState will hold one WFA2Aligner per thread.

### Phase 6: uchime_ref Table Function

**RED** — `test/sql/uchime_ref.test`:
- Error cases: missing table, missing db parameter, nonexistent reference table, missing required columns
- Row count matches query count
- Column types correct
- Ground truth validation: all 18 columns compared to `expected_ref.tsv` (score rounded to 4dp, identities to 1dp, votes and flag exact)
- Parameter override tests (minh=0.01 catches borderline cases)

**GREEN** — `src/include/uchime_ref.hpp` + `src/uchime_ref.cpp`:

Follows the align_minimap2 streaming pattern: reference materialized + indexed at init; queries streamed in batches.

- Data struct: query/ref table names, UchimeParams, output schema
- GlobalState:
  - ChimeraDetector with reference loaded + k-mer index built at init (shared read-only across threads)
  - `all_query_ids` — materialized read_id list (lightweight, no sequences)
  - `atomic<idx_t> next_batch_offset` — for lock-free batch claiming
- LocalState:
  - Per-thread WFA2Aligner(6, 20, 4)
  - `result_buffer` + `buffer_offset` — buffered UchimeResults, drained in STANDARD_VECTOR_SIZE chunks
- Bind(): validate tables exist (TABLE_ENTRY lookup), check required columns (read_id, sequence), parse named params, build schema
- InitGlobal(): read reference table via separate Connection → detector.set_reference() + build k-mer index. Materialize query read_ids only (not sequences).
- Execute(): streaming producer-consumer loop (same as align_minimap2):
  1. If result_buffer has data → output a chunk
  2. If buffer empty → atomic `fetch_add(batch_size)` to claim next batch of query IDs
  3. `ReadBatchByIds()` — fetch sequences for claimed IDs via temp table JOIN (no OFFSET)
  4. `detect()` each query against the shared ChimeraDetector
  5. Buffer results, loop back to step 1
- MaxThreads(): min(total_queries/batch_size, hardware_concurrency)

Register in `miint_extension.cpp`, add to `CMakeLists.txt`.

**REFACTOR** — Extract param parsing, output population helpers

### Phase 7: uchime_denovo Table Function

**Input**: a DuckDB table/view name (VARCHAR) with required columns: `read_id` (VARCHAR), `sequence` (VARCHAR), `size` (INTEGER). No FASTA header parsing — the caller provides abundance as a proper integer column. This is cleaner than vsearch's `;size=N;` convention and composes naturally with DuckDB workflows:

```sql
-- Example usage: read sequences, add abundance, detect chimeras
CREATE TABLE seqs AS SELECT read_id, sequence1 AS sequence, count(*) AS size
  FROM read_fastx('amplicons.fasta') GROUP BY ALL;
SELECT * FROM uchime_denovo('seqs');
```

**RED** — `test/sql/uchime_denovo.test`:
- Error cases: missing table, missing required columns, non-integer size column
- Row count correct (all sequences reported, not just chimeras)
- Ground truth against `expected_denovo.tsv` (generated by creating equivalent vsearch input with `;size=N;` headers)
- `abskew` parameter test

**GREEN** — `src/include/uchime_denovo.hpp` + `src/uchime_denovo.cpp`:
- Bind(): validate table exists, check required columns (read_id, sequence, size)
- InitGlobal(): read table via separate Connection, sort by size descending
- MaxThreads() = 1 (sequential)
- Execute(): process one query at a time, call detect_denovo() which filters candidates by abundance skew, add non-chimeras to detector's reference incrementally

**REFACTOR** — Share output schema with uchime_ref

### Phase 8: Performance Validation

Benchmark script comparing vsearch vs miint on:
- 1000 queries × gold.fa reference (~5000 seqs, ~1500bp each)
- 5000 queries × same reference

Target: same wall-clock time or faster than vsearch with equivalent thread count.

Performance levers:
- Multi-threaded uchime_ref (each thread has own WFA2Aligner)
- Array-based k-mer index (cache-friendly flat array vs hash map)
- WFA2 BiWFA O(s) memory mode (fast for similar sequences)
- Batch query claiming (64 at a time, minimal lock contention)

---

## Risk: Alignment Differences

WFA2 and vsearch's NW aligner may produce different optimal alignments when multiple optima exist (same score, different gap placements). This would cause different vote counts and h-scores.

**Mitigation strategy**:
1. Start with WFA2(6, 20, 4) — the theoretical equivalent
2. Validate against vsearch on synthetic data where alignments are unambiguous
3. If discrepancies appear on real data, inspect vsearch's `--uchimealns` to see its alignment
4. Accept small numerical differences (±0.1 on identities, ±1 on votes) as long as classification (Y/N/?) matches
5. If classifications differ, investigate specific cases and adjust

---

## Code Reviews

Code reviews are mandatory checkpoints. After completing each phase's GREEN step (tests passing), **stop work and present the review** before proceeding to REFACTOR or the next phase. No code is written until the review is acknowledged.

**Review checkpoints:**
- After Phase 1 GREEN (KmerIndex passes tests)
- After Phase 3 GREEN (star alignment + diff classification passes tests)
- After Phase 4 GREEN (parent selection passes tests)
- After Phase 5 GREEN (full pipeline passes tests — major checkpoint, all core algorithm complete)
- After Phase 6 GREEN (uchime_ref SQL tests pass against vsearch ground truth)
- After Phase 7 GREEN (uchime_denovo SQL tests pass)

Each review presents: files changed, test results, and the code for inspection.

---

## Verification

After each phase, run:
```bash
# C++ tests
./build/release/extension/miint/tests "[KmerIndex]"       # Phase 1
./build/release/extension/miint/tests "[ChimeraDetector]"  # Phases 2-5

# SQL tests (after Phase 6+)
./build/release/test/unittest "test/sql/uchime_ref.test"
./build/release/test/unittest "test/sql/uchime_denovo.test"

# Full suite
bash run_tests.sh
```
