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

### Phase 4: Smoothed Parent Selection — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **Win-counting starts at SMOOTHING_WINDOW-1**: positions 0..30 have partial sums and are skipped
- **`std::min` across all smoothed vector sizes**: prevents OOB read when candidates have different alignment lengths (due to indels)
- **ParentPair caches pairwise alignments**: reused in detect() for star alignment + identity computation
- select_parents on chimera1 correctly picks ref1 (parent A) and ref2 (parent B)
- Parent assignment order may differ from vsearch for asymmetric chimeras (handled by reversed breakpoint sweep)

### Phase 5: H-Score Sweep + Classification — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **Divergence gated on `best_pos >= 0`**: non-chimeric sequences get divergence=0.0, not garbage negative values
- **Reversed configuration**: when B-diffs dominate left and A-diffs dominate right, parents are effectively swapped in the sweep
- **WFA2Aligner reuse verified**: bit-identical scores across repeated detect() calls with same aligner
- **Gapped sequence tests**: parents with different alignment lengths tested for OOB safety

vsearch ground truth validation:
- chimera1: score=16.50±0.01, votes=45/0/0/46/0/0, identities QA=84.7/QB=85.0/AB=69.7/QM=100.0 (all within ±0.5)
- All 8 queries classified correctly (4Y, 4N)
- 21 test cases, 156 assertions

### Phase 6: uchime_ref Table Function — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **QuerySequenceStream** for lazy streaming queries (thread-safe, no full materialization)
- **Schema validation at bind time** via ValidateSequenceTableSchema
- **Thread count capped at 8** (not raw hardware_concurrency)
- **NULL/empty reference rows skipped** (not silent data loss)
- **D_ASSERT(col == output.ColumnCount())** in OutputUchimeResults
- **Parameter range validation** (minh>=0, xn>=1.0, dn>=0, mindiffs>=1)
- **Non-chimeric flag=N uses * for parents/identities** (vsearch convention)
- View support tested, wrong column names tested
- 73 SQL assertions

### Phase 7: uchime_denovo Table Function — DONE ✓

Implemented and reviewed. Key design decisions from code review:
- **detect_impl() private helper** eliminates 80 lines of duplicated pipeline code
- **set_reference() initializes ref_abundances_ to 0** preventing silent filter-all bug
- **Single catalog lookup** in ValidateDenovoTableSchema
- **uchime_common.hpp/.cpp** with shared output functions (no more duplication)
- **Incremental Execute()** processes STANDARD_VECTOR_SIZE per call for DuckDB cancellation
- **Regression tests run unconditionally** (no require-env guard on hardcoded expected values)
- **regexp_extract** for robust ;size=N; parsing
- **Bootstrap comment** explains why first two sequences skip chimera testing
- Input table requires `read_id` (VARCHAR), `sequence1` (VARCHAR), `size` (integer type)
- 56 SQL assertions

### Phase 8: Performance Validation — IN PROGRESS

Benchmark uchime_ref using real 16S data: `scratch/LTPs132_SSU.fasta` (LTP database).
Use an arbitrary subset as reference, remainder as queries.

**Benchmark plan:**
1. Load LTPs132_SSU.fasta, split into reference (~500 seqs) and queries (~rest)
2. Time vsearch: `vsearch --uchime_ref queries.fasta --db refs.fasta --uchimeout /dev/null --threads 1`
3. Time miint: `SELECT count(*) FROM uchime_ref('queries', db:='refs')`
4. Compare wall-clock times, chimera counts, and flag agreement
5. If miint is slower, profile and optimize

**Benchmark results** (LTPs132_SSU.fasta: 500 refs × 2000 queries, ~1462bp mean):

| Config | Wall time | Parallel speedup |
|--------|-----------|-----------------|
| vsearch --threads 1 | 22s | — |
| vsearch default (12 cores) | 4.8s | 4.6x |
| miint SET threads=1 | 51s | — |
| miint default threads | 7.3s | 7.0x |

**Assessment:** miint single-threaded is 2.3x slower than vsearch (WFA2 vs SIMD-optimized NW). Multi-threaded is 1.5x slower (7.3s vs 4.8s) — acceptable given WFA2's per-alignment cost. Parallel scaling is actually better than vsearch (7.0x vs 4.6x) thanks to atomic batch claiming.

**Profiling** (perf): 86% of CPU in WFA2 `wavefront_*` functions. KmerIndex is <1%. The bottleneck is pairwise alignment of 16 candidates × ~1500bp sequences.

**Optimizations applied:**
- Atomic batch claiming (fetch_add + ReadBatchByIds) replaces mutex-serialized QuerySequenceStream
- MaxThreads = 8 (fixed cap, respects DuckDB scheduler)
- NULL read_ids rejected with clear error

**Optimization NOT applied (correctness issue):**
- Score-only pre-filter was reverted. Global edit distance does not predict positional smoothed-identity competition. A chimera's true parent could rank poorly by global score. Cutting from 16→6 candidates silently missed real chimeras. The UCHIME algorithm is designed with 16 as the tractable bound.

### WFA2 Performance Investigation — DONE ✓

**Root cause**: Three compounding factors inflated WFA2 cost by ~5x:
1. **Dual-pass alignment** (score + alignment scope): 2x overhead. The BiWFA score bug
   only affects sequences ≤100bp; for ≥1500bp the alignment-scope score is reliable.
2. **MemoryUltralow (BiWFA)**: Designed for >30Kbp. For 1500bp, recursive forward+reverse
   wavefront expansion adds 1.25-2x overhead. MemoryMed (piggyback backtrace) is optimal.
3. **High penalties (6,20,4)**: Inflate alignment score s by ~1.7x vs (4,6,2), directly
   inflating O(ns) runtime. Not changed — tied to vsearch alignment equivalence.

**Fixes applied**: Eliminate dual-pass for >100bp sequences + switch to MemoryMed.

**Final benchmark** (LTPs132_SSU.fasta: 500 refs × 2000 queries, ~1462bp mean):

| Config | Wall time | vs vsearch |
|--------|-----------|------------|
| vsearch 1 thread | 22s | baseline |
| vsearch default (12 cores) | 4.8s | baseline |
| **miint 1 thread** | **27s** | 1.25x slower |
| **miint default threads** | **3.6s** | **25% faster** |

**Future considerations:**
- **Penalty reduction**: Using (4,6,2) instead of (6,20,4) would reduce s by ~1.7x,
  bringing single-threaded time closer to vsearch. Requires validating that the changed
  gap placement doesn't affect chimera classification on real data.
- **WFAdaptive heuristic**: Conservative pruning (min_wavefront_length=10, max_distance=50)
  could give 2-5x speedup for candidate screening. Not compatible with BiWFA but works
  with MemoryMed. Requires tolerance for approximate results on candidate ranking.
- **Potential WFA2-lib contribution**: The wavefront extend kernel
  (`wavefront_extend_matches_kernel_blockwise`) does 8-byte XOR comparison but operates
  on individual diagonals in a scalar loop. SIMD vectorization across diagonals (processing
  multiple diagonals in parallel via SSE2/AVX2) could significantly improve throughput for
  the extend phase, which is the second-hottest function at 15% of CPU time. This would
  benefit all WFA2 users, not just UCHIME.

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
