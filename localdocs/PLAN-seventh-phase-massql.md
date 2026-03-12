# MassQL Full Parity — TDD Plan

## Context

Phases 1-5 implemented 36+ DuckDB macros. Phase 6b added a `massql()` table function with parser and transpiler (25/26 cycles done). The transpiler now defaults to 0.1 Da tolerance (matching Python MassQL). Several MassQL features are still missing: INTENSITYTICPERCENT, INTENSITYMATCH/Y-variable system, ANY wildcard, MASSDEFECT qualifier, scanrangesum(), aminoaciddelta(), peptide(), OTHERSCAN, division operator, `<` qualifier operator, and MATCHCOUNT alias.

**Goal**: Full MassQL parity — every valid MassQL query produces identical results to the Python reference implementation. Priorities: (1) verifiable correctness, (2) usability for non-SQL users.

**Out of scope**: Mixed X-variable patterns (separate plan — parser must emit clear error). MOBILITY (no column in `read_mzml` — parser must emit clear error).

## Completed Work

| Step | Commit | Summary |
|------|--------|---------|
| Pre-work A | `1ee49d4` | Phase 6b transpiler: parser, transpiler, massql() + massql_to_sql() functions, 13 files |
| Pre-work B | `984e37b` | Ground truth baseline: 33 Python MassQL queries, JSON reference, SQL regression test |
| Cycle 1 | `3e0b57c` | Default tolerance 0.5→0.1 Da in 3 transpiler locations + 13 macros |
| Cycle 2 | (in progress) | MATCHCOUNT alias, `<` operator, `/` division. Python MassQL `<` operator is buggy (see `localdocs/massql-bug.md`) — we implement correct semantics |

---

## Pre-work A: Commit Phase 6b Transpiler

Commit the 10 new files + 4 modified files already implemented (parser, transpiler, function, formula_parser, tests). This is the foundation everything builds on.

**Files**: `CMakeLists.txt`, `src/miint_extension.cpp`, `src/formula_function.cpp`, + 10 new `src/massql_*.cpp`, `src/formula_parser.cpp`, `test/cpp/test_MassQLParser.cpp`, `test/sql/massql_transpiler.test`

---

## Pre-work B: Ground Truth Baseline

Before any feature work, establish a regression net against Python MassQL using queries that work today.

1. Create Python script `test/scripts/massql_ground_truth.py`:
   - Loads `basic_3spectra.mzML` into Python MassQL
   - Runs ~10 queries covering currently-supported features (simple conditions, AND, OR, X-offset, cross-level, EXCLUDED, FILTER)
   - Records result counts + key values to a reference file
2. Create SQL test file `test/sql/massql_ground_truth.test` gated by `require-env MASSQL_PYTHON_AVAILABLE`:
   - Assert `massql()` result counts match Python reference counts for the same ~10 queries
3. In `run_tests.sh`: auto-detect Python MassQL, export `MASSQL_PYTHON_AVAILABLE`

This baseline grows as new features land — each group adds its queries to the ground truth file.

---

## Pre-work C: Test Data Design for Group D (Y-variable)

Design concrete test mzML files before any Y-variable implementation begins.

### `data/mzml/intensity_match.mzML`

Y-variable testing. 4 MS2 spectra with controlled intensity ratios:

| Spectrum | Peaks (m/z → intensity) | Purpose |
|----------|------------------------|---------|
| 0 | 150→1000, 250→500, 350→100 | Y*0.5 match at 250 (ref=150), Y*0.1 match at 350 |
| 1 | 150→1000, 250→490, 350→100 | Y*0.5 within 10% tolerance at 250 |
| 2 | 150→1000, 250→200, 350→100 | Y*0.5 does NOT match at 250 (too low) |
| 3 | 150→1000, 250→500, 350→500 | Two peaks at same intensity — tests non-unique Y |

Precursor m/z: 400.0 for all spectra. Retention times: 1.0, 2.0, 3.0, 4.0 min.

### `data/mzml/isotope_pattern_xy.mzML`

X+Y isotope testing. 3 MS1 spectra:

| Spectrum | Peaks (m/z → intensity) | Purpose |
|----------|------------------------|---------|
| 0 | 500.0→10000, 501.003→3000, 502.006→500 | Typical isotope envelope (M, M+1, M+2). M+1 = Y*0.30, M+2 = Y*0.05 |
| 1 | 600.0→8000, 601.003→5000, 602.006→2000 | Different ratios: M+1 = Y*0.625, M+2 = Y*0.25 |
| 2 | 700.0→5000, 701.003→1500 | Only M and M+1 — M+2 missing, should NOT match 3-peak isotope pattern |

### `data/mzml/degenerate.mzML`

Edge cases:

| Spectrum | Peaks | Purpose |
|----------|-------|---------|
| 0 | (none) | Empty spectrum, 0 peaks |
| 1 | 100.0→500 | Single peak |
| 2 | 200.0→1000, 200.0→500 | Duplicate m/z |
| 3 | 15000.0→100 | Very large m/z |
| 4 (MS1) | precursor_mz=300 | Malformed: precursor on MS1 |

Generate all with `test/scripts/generate_mzml_test_data.py` (extend existing script). **This pre-work must complete before Cycle 10 begins.**

---

## Pre-work D: Out-of-Scope Error Messages

Add explicit parser errors for features that are out of scope, so users get useful feedback instead of garbage SQL.

**RED**: C++ parser tests:
- `MS1MZ=X AND MS2PROD=Y` (mixed X-variable) → error: `"Mixed X and Y variable patterns are not yet supported"`
- `MS2PROD=220:MOBILITY=5` → error: `"MOBILITY is not supported (read_mzml does not provide ion mobility data)"`

**GREEN**: Add checks in parser for these tokens, throw descriptive `InvalidInputException`.

---

## Group A: Foundation (Cycles 1-3)

### Cycle 1: Default tolerance 0.5 → 0.1 Da

Match Python MassQL's default.

**RED**:
- C++ test: `massql_to_sql("QUERY MS2DATA WHERE MS2PROD=220", "src")` contains `mz_within(mz, 220, 0.1)` not `0.5`
- SQL test: verify result counts unchanged with basic_3spectra.mzML (peaks at exact integers, so 0.1 still matches)

**GREEN** — 3 transpiler locations + 13 macro defaults:
- `src/massql_transpiler.cpp:107` — simple condition default
- `src/massql_transpiler.cpp:163` — X-offset pattern default
- `src/massql_transpiler.cpp:315` — cross-level pattern default
- `src/include/miint_macros.hpp` — 13 macros with `tolerance := 0.5` → `tolerance := 0.1`:
  - `mzml_ms1_where_ms2prod`, `mzml_ms2_where_ms1mz`, `mzml_ms1_where_ms2prec`
  - `mzml_ms2_where_ms2prod_and_ms1mz`
  - `mzml_x_offset_pair`, `mzml_x_offset_triplet`, `mzml_x_offset_ntuple`
  - `mzml_x_prec_prod`, `mzml_x_ms1_ms2_prec`, `mzml_x_offset_pair_range`
  - `mzml_excluded_ms2prod`, `mzml_excluded_ms1mz`, `mzml_excluded_ms2prec`

**REFACTOR**: Extract `DEFAULT_TOLERANCE_DA = 0.1` constant in transpiler. Fix any broken test assertions (tests using explicit `0.5` are fine; tests relying on the default need updating). Update design decision #3 in PLAN doc.

**Risk**: Tests in `massql_one_liners.test`, `massql_x_offset.test`, `massql_x_offset_ntuple.test`, `massql_x_prec_prod.test` use default tolerance — need to check if peak spacing requires wider tolerance. If so, add explicit `:TOLERANCEMZ=0.5` to those test queries.

---

### Cycle 2: Parser completeness (MATCHCOUNT + `<` operator + division)

Batch three trivial parser additions into a single cycle.

**RED**:
- C++ parser test: `MS2PROD=(150 OR 250):MATCHCOUNT=range(min=1,max=2)` parses identically to CARDINALITY
- C++ parser test: `MS2PROD=220:INTENSITYPERCENT<50` parses with `op = LESS_THAN`
- C++ parser test: `MS2PROD=400/2` folds to 200.0

**GREEN**:
- Add `"MATCHCOUNT"` → same handling as `"CARDINALITY"` in parser qualifier map
- Add `LESS_THAN` to `QualifierOp` enum in `massql_parser.hpp`
- Tokenizer: handle `<` and `/` tokens
- Parser: accept `<` in qualifier position
- Transpiler: emit `<` in intensity expressions (currently emits `>=` for `=` and `>` for `>`)
- Constant folding: handle division in `parse_term()`

**Transpiler snapshot**: Capture `massql_to_sql()` golden output for 5 representative Group A queries.

---

## Group B: Missing Qualifiers & Values (Cycles 3-5)

### Cycle 3: INTENSITYTICPERCENT

Commonly used qualifier. `mzml_peaks()` already provides `i_tic` column (intensity / total_ion_current).

**RED**:
- C++ parser test: `MS2PROD=220:INTENSITYTICPERCENT=5` parses as qualifier
- SQL test: query with INTENSITYTICPERCENT=5 returns only peaks where `i_tic >= 0.05`

**GREEN**:
- Parser: add `"INTENSITYTICPERCENT"` to qualifier map (numeric value)
- Transpiler: emit `AND i_tic >= {value}/100.0` (or `>` / `<` based on operator) in peak matching WHERE clause, analogous to existing INTENSITYPERCENT → `i_norm`

**Note**: Python MassQL caps intensity percentages at 99% to avoid matching nothing. Add same cap.

---

### Cycle 4: ANY wildcard + MASSDEFECT qualifier

`MS2PROD=ANY` matches any peak (no m/z filter). MASSDEFECT filters by fractional mass. These compose naturally: `MS2PROD=ANY:MASSDEFECT=massdefect(min=0.1,max=0.2)`.

**RED**:
- C++ parser test: `MS2PROD=ANY` parses with `has_any_wildcard = true`
- C++ parser test: `MS2PROD=ANY:MASSDEFECT=massdefect(min=0.1,max=0.2)` parses with massdefect min/max stored
- SQL test: `MS2PROD=ANY` returns all MS2 peaks (no m/z filter)
- SQL test: MASSDEFECT query returns only peaks where `mz - FLOOR(mz) BETWEEN 0.1 AND 0.2`

**GREEN**:
- Parser: recognize `ANY` token as special value, set flag on `ConditionValue`
- Transpiler: when ANY, omit the `mz_within()` clause — just use ms_level filter + any intensity/massdefect qualifiers
- Parser: handle `MASSDEFECT=massdefect(min=N,max=N)` qualifier syntax (similar to CARDINALITY range parsing)
- Transpiler: emit `AND (mz - FLOOR(mz)) > {min} AND (mz - FLOOR(mz)) < {max}` in peak matching WHERE clause
- For X-level mass defect (`X=massdefect(min,max)`): add `WHERE (mz - FLOOR(mz)) > min AND (mz - FLOOR(mz)) < max` to `__x_candidates` CTE

---

### Cycle 5: scanrangesum() aggregation

EBNF line 117. Bins peaks by m/z in tolerance-width bins and sums intensities. Takes an optional TOLERANCE parameter: `scanrangesum(MS2DATA, TOLERANCE=0.1)`.

**RED**:
- C++ parser test: `QUERY scanrangesum(MS2DATA) WHERE MS2PROD=220` parses with `agg = SCANRANGESUM`
- C++ parser test: `QUERY scanrangesum(MS2DATA, TOLERANCE=0.5)` parses with tolerance param
- SQL test: returns binned m/z + summed intensity per bin per spectrum

**GREEN**:
- Parser: add `SCANRANGESUM` to `AggFunction` enum; handle optional `,TOLERANCE=N` parameter in function call
- Transpiler: emit `SELECT spectrum_index, FLOOR(mz / {tol}) * {tol} AS mz_bin, SUM(intensity) AS total_intensity FROM __final GROUP BY spectrum_index, mz_bin`
- Default bin width = 0.1 Da (matching Python)

**Output schema**: `spectrum_index, mz_bin, total_intensity`

**Transpiler snapshot**: Capture golden output for 5 representative Group B queries. Add Group A+B queries to ground truth file.

---

## Group C: New Functions (Cycles 6-7)

### Cycle 6: aminoaciddelta() function

EBNF line 147. Computes sum of amino acid residue masses from single-letter codes. Example: `aminoaciddelta(AG)` = mass(Ala) + mass(Gly).

**RED**:
- C++ parser test: `MS2PROD=aminoaciddelta(AG)` folds to 128.05858 (71.03711 + 57.02147)

**GREEN**:
- Add amino acid monoisotopic residue mass table (20 standard amino acids) to a new `src/mass_tables.cpp` / `src/include/mass_tables.hpp`
- Implement `AminoacidDeltaMass(string)` → double
- Parser: recognize `aminoaciddelta(` token, parse uppercase sequence, call `AminoacidDeltaMass`, fold to constant
- Works in any expression position: `MS2NL=aminoaciddelta(AG)`, `MS2PROD=X+aminoaciddelta(K)`

**Mass table source**: Same as Python MassQL (pyteomics masses).

---

### Cycle 7: peptide() function

EBNF lines 150-160. Computes theoretical mass of a peptide fragment ion. `peptide(ACDE, charge=2, ion=b)`.

**RED**:
- C++ parser test: `MS2PROD=peptide(ACDE,charge=2,ion=b)` folds to expected mass

**GREEN**:
- Implement `PeptideFragmentMass(sequence, charge, ion_type)` → double in `src/mass_tables.cpp`
- Ion types: a, b, c (N-terminal); x, y, z (C-terminal)
- Formula: sum residue masses + ion type adjustment (b = +proton, y = +water+proton, etc.) / charge
- Parser: recognize `peptide(` token, parse `SEQUENCE,charge=N,ion=X)` syntax, fold to constant

**Transpiler snapshot**: Capture golden output for Group C queries. Add to ground truth file.

---

## Group D: Y-Variable / Intensity Matching (Cycles 8-12)

This is the most complex feature. Y is a per-scan intensity reference, not an m/z candidate like X. Three qualifiers work together: INTENSITYMATCHREFERENCE (marker), INTENSITYMATCH=expr (expected intensity), INTENSITYMATCHPERCENT=N (tolerance band).

**Reference implementation**: `mzml_isotope_pattern` macro (miint_macros.hpp:739-793) already demonstrates the X+Y SQL pattern with `ref_intensity`, `target_intensity`, `matched_targets` CTEs.

**Prerequisite**: Pre-work C (test data) must be complete before starting this group.

### Cycle 8: Y-variable parser

**RED**: C++ parser tests for:
1. `MS1MZ=X:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE` → qualifier with no value, marks reference
2. `MS1MZ=X+2:INTENSITYMATCH=Y*0.66:INTENSITYMATCHPERCENT=30` → Y expression parsed
3. `INTENSITYMATCH=Y*(0.0608+(0.000002*X))` → X-dependent Y expression parsed
4. `INTENSITYMATCH=Y` → coefficient 1.0

**GREEN**:
- Add fields to `Qualifier` struct:
  ```cpp
  double y_expr_constant = 1.0;   // Y * constant  (or Y * (constant + x_coeff*X))
  double y_expr_x_coeff = 0.0;    // coefficient of X in Y expression
  bool y_expr_has_x = false;
  ```
- Add `INTENSITYMATCHREFERENCE` as no-value qualifier (like EXCLUDED)
- Add `INTENSITYMATCH` as expression qualifier with dedicated `parse_intensity_match_expr()`:
  - Handles: `Y`, `Y*number`, `Y*(number)`, `Y*(number+number*X)`, `Y*(number+(number*X))`
- Add `INTENSITYMATCHPERCENT` as numeric qualifier

---

### Cycle 9: Y-only transpiler (no X variable)

Simple case: fixed m/z values with intensity relationship. Uses `data/mzml/intensity_match.mzML`.

**RED**: SQL test:
```sql
-- Spectrum 0: 150→1000, 250→500 → Y*0.5 matches
-- Spectrum 2: 150→1000, 250→200 → Y*0.5 does NOT match
SELECT COUNT(*) FROM massql(
  'QUERY scaninfo(MS2DATA) WHERE MS2PROD=150:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS2PROD=250:INTENSITYMATCH=Y*0.5:INTENSITYMATCHPERCENT=10',
  '__massql_test');
-- Expected: 2 (spectra 0 and 1, not 2)
```

**GREEN** — Transpiler Y-detection + CTE generation:
1. **Detect**: scan conditions for INTENSITYMATCHREFERENCE → `y_ref_idx`; scan for INTENSITYMATCH → `y_match_indices`
2. **Generate `__y_ref` CTE**: match reference m/z, `SUM(intensity) AS y_val` per spectrum_index
3. **Generate `__y_match_N` CTE per match condition**: compute actual intensity at target m/z, JOIN with `__y_ref`, apply tolerance band:
   ```sql
   WHERE actual_int > y_val * {y_expr_constant} * (1.0 - {pct}/100.0)
     AND actual_int < y_val * {y_expr_constant} * (1.0 + {pct}/100.0)
   ```
4. **Qualifying**: INTERSECT all `__y_match_N` spectrum_index sets

---

### Cycle 10: X+Y offset pattern transpiler

The common isotope detection case. Combines X-enumeration with Y-intensity matching. Uses `data/mzml/isotope_pattern_xy.mzML`.

**RED**: SQL test:
```sql
-- Spectrum 0: 500→10000, 501.003→3000, 502.006→500
-- M+1 = Y*0.30, so Y*0.30 ± 30% should match
SELECT COUNT(*) FROM massql(
  'QUERY scaninfo(MS1DATA) WHERE MS1MZ=X:TOLERANCEMZ=0.01:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE:INTENSITYPERCENT=25 AND MS1MZ=X+1.003:TOLERANCEMZ=0.01:INTENSITYMATCH=Y*0.30:INTENSITYMATCHPERCENT=30',
  '__massql_test');
```

**GREEN** — CTE structure (follows `mzml_isotope_pattern` macro pattern):
1. `__x_candidates` — recursive CTE enumerating distinct m/z values (existing pattern)
2. `__y_ref` — JOIN candidates with base peaks at `x_val + ref_offset`, compute `SUM(intensity) AS y_val` per `(x_val, spectrum_index)`, apply ref qualifiers (INTENSITYPERCENT, etc.)
3. `__y_match_N` — for each match condition: compute actual intensity at `x_val + match_offset`, JOIN with `__y_ref` on `(x_val, spectrum_index)`, apply Y tolerance band
4. `__x_qualifying` — JOIN all `__y_match_N` on `(x_val, spectrum_index)` to require ALL matches at same X
5. `__filtered` — existing pattern

**REFACTOR**: Factor shared CTE generation between Y-only and X+Y patterns into a `generate_y_ctes()` helper.

---

### Cycle 11: Y expressions with X dependency

Handle `INTENSITYMATCH=Y*(0.0608+(0.000002*X))` where the expected intensity depends on both Y and X.

**RED**: C++ parser test + SQL test for X-dependent Y expression.

**GREEN**:
- SQL expression helper:
  ```cpp
  // Y * (constant + x_coeff * X) → "r.y_val * (0.0608 + 0.000002 * xc.x_val)"
  static string y_expr_sql(const Qualifier &q, const string &y_col, const string &x_col);
  ```
- In `__y_match_N` CTE, the tolerance band uses `y_expr_sql()` output as the expected value

---

### Cycle 12: Integration of Y with non-Y conditions + cross-level

Ensure Y-conditions compose with:
- Non-Y WHERE conditions (metadata filters, plain peak conditions)
- Cross-level X patterns (`MS1MZ=X:INTENSITYMATCHREFERENCE AND MS2PREC=X`)
- FILTER clause

**RED**: SQL tests for:
1. Y-pattern + metadata filter (RTMIN)
2. Cross-level: `MS1MZ=X:INTENSITYMATCHREFERENCE AND MS1MZ=X+2:INTENSITYMATCH=Y*0.66:INTENSITYMATCHPERCENT=30 AND MS2PREC=X:TOLERANCEMZ=2`
3. Y-pattern with FILTER clause

**GREEN**: Ensure transpiler composes Y CTEs with existing metadata filtering and non-Y match CTEs via INTERSECT in `__filtered`.

**Transpiler snapshot**: Capture golden output for Group D queries. Add Y-variable queries to ground truth file.

---

## Group E: OTHERSCAN Qualifier (Cycle 13)

### Cycle 13: OTHERSCAN=rtrange(left,right)

Cross-scan matching by retention time proximity. `MS2PROD=226.18:OTHERSCAN=rtrange(left=0.5,right=0.5)` means "match if this OR a nearby scan (±0.5 min RT) has the peak".

**RED**:
- C++ parser test: qualifier parses with left/right values
- SQL test: query matches scan even when peak is in a neighboring scan within RT window

**GREEN**:
- Parser: handle `OTHERSCAN=rtrange(left=N,right=N)` syntax (similar to CARDINALITY range parsing)
- Transpiler: for OTHERSCAN conditions, expand the match CTE to include spectra whose retention_time is within `[rt - left, rt + right]` of a matching spectrum:
  ```sql
  __match_0 AS (
    SELECT DISTINCT b2.spectrum_index
    FROM __base b1
    JOIN __base b2 ON b2.retention_time >= b1.retention_time - {left}
                  AND b2.retention_time <= b1.retention_time + {right}
    WHERE mz_within(b1.mz, {target}, {tol})
  )
  ```

**Performance note**: DuckDB's IEJoin handles the 2-predicate range join (`>=` and `<=`) in O(n log n), not O(n²). No optimization needed.

**Transpiler snapshot**: Capture golden output for OTHERSCAN. Add to ground truth file.

---

## Group F: Quality & Verification (Cycles 14-17)

### Cycle 14: Error message quality

Improve parser error messages for non-SQL users.

**RED**: C++ tests asserting error messages contain position info and suggestions.

**GREEN**:
- **Position tracking**: Track character offset during tokenization. Errors report: `"Unexpected 'BLAH' at position 18 (after 'MS2DATA')"`
- **"Did you mean?" suggestions**: For common misspellings:
  - `MS2PRODUCT` → "Did you mean MS2PROD?"
  - `SCANINF` → "Did you mean scaninfo?"
  - `TOLERENCE` → "Did you mean TOLERANCEMZ?"
- **Unsupported feature messages**: `MOBILITY` → "MOBILITY is not supported (read_mzml does not provide ion mobility data)"

---

### Cycle 15: Integration ground truth vs Python MassQL (final expansion)

Expand the ground truth baseline (from Pre-work B) to cover all implemented features.

**Queries to verify** (at minimum):
- `QUERY MS2DATA WHERE MS2PROD=220:TOLERANCEMZ=0.1`
- `QUERY scaninfo(MS2DATA) WHERE MS2PREC=200:TOLERANCEMZ=0.1`
- `QUERY scannum(MS2DATA) WHERE MS2NL=80:TOLERANCEMZ=0.1`
- `QUERY MS2DATA WHERE MS2PROD=220:TOLERANCEMZ=0.1:INTENSITYPERCENT=25`
- `QUERY scaninfo(MS2DATA) WHERE MS2PROD=150:TOLERANCEMZ=0.1 AND MS2PROD=250:TOLERANCEMZ=0.1`
- `QUERY scaninfo(MS2DATA) WHERE MS2PROD=(150 OR 250):TOLERANCEMZ=0.1:CARDINALITY=range(min=1,max=2)`
- `QUERY scaninfo(MS2DATA) WHERE MS2PROD=220:TOLERANCEMZ=0.1:EXCLUDED`
- X-offset, cross-level, Y-variable isotope patterns
- FILTER clause queries
- scanrangesum queries
- OTHERSCAN queries
- ANY + MASSDEFECT queries
- aminoaciddelta / peptide constant-folding queries

---

### Cycle 16: Kitchen-sink composition tests

Test that features compose correctly in combinations the individual cycles don't cover.

**RED**: SQL tests for at least these queries:
1. `MS2PROD=ANY:MASSDEFECT=massdefect(min=0.1,max=0.2):INTENSITYMATCH=Y*0.5:INTENSITYMATCHPERCENT=20` — ANY + MASSDEFECT + Y-variable
2. `MS1MZ=X:INTENSITYMATCH=Y:INTENSITYMATCHREFERENCE AND MS1MZ=X+1.003:INTENSITYMATCH=Y*0.3:INTENSITYMATCHPERCENT=30:OTHERSCAN=rtrange(left=0.5,right=0.5)` — X+Y + OTHERSCAN
3. `MS2PROD=aminoaciddelta(AG):INTENSITYTICPERCENT=5:EXCLUDED` — function + INTENSITYTICPERCENT + EXCLUDED
4. `QUERY scanrangesum(MS2DATA) WHERE MS2PROD=ANY:MASSDEFECT=massdefect(min=0.9,max=0.99)` — scanrangesum + ANY + MASSDEFECT

**GREEN**: Fix any transpiler bugs revealed by composition failures.

---

### Cycle 17: Performance regression + degenerate edge cases

**Performance**:
- Create synthetic mzML with ~10,000 spectra, ~100 peaks each (`data/mzml/perf_10k.mzML`)
- Benchmark: simple query < 1s, X-offset < 5s, X+Y isotope < 10s, OTHERSCAN < 5s
- Run as separate perf test, not in main suite

**Degenerate edge cases** (uses `data/mzml/degenerate.mzML`):
- Empty spectrum (0 peaks)
- Duplicate precursor_mz across spectra
- Single-peak spectrum
- Duplicate peaks at same m/z
- Very large m/z (>10000)
- Precursor_mz with ms_level=1 (malformed)

---

## Files Modified (Summary)

| File | Changes |
|------|---------|
| `src/include/massql_parser.hpp` | Add LESS_THAN to QualifierOp; add y_expr fields to Qualifier; add ANY flag to ConditionValue; add SCANRANGESUM to AggFunction; add massdefect fields |
| `src/massql_parser.cpp` | New qualifiers, ANY/MASSDEFECT/division parsing, Y-expression parser, aminoaciddelta/peptide parsing, OTHERSCAN, error position tracking, out-of-scope error messages |
| `src/massql_transpiler.cpp` | Default 0.1; INTENSITYTICPERCENT/MASSDEFECT/ANY/OTHERSCAN SQL gen; Y-variable CTE generation; scanrangesum aggregation; `<` operator |
| `src/include/miint_macros.hpp` | 13 macros: `tolerance := 0.5` → `tolerance := 0.1` |
| `src/mass_tables.cpp` | NEW — amino acid residue mass table, `AminoacidDeltaMass()`, `PeptideFragmentMass()` |
| `src/include/mass_tables.hpp` | NEW — declarations for mass table functions |
| `test/cpp/test_MassQLParser.cpp` | New test cases for all parser additions |
| `test/sql/massql_transpiler.test` | New SQL assertions for all transpiler additions |
| Various `test/sql/massql_*.test` | Fix any tests broken by 0.5→0.1 default change |
| `test/scripts/generate_mzml_test_data.py` | New test data generation |
| `test/scripts/massql_ground_truth.py` | Python ground truth script |
| `run_tests.sh` | Auto-detect MASSQL_PYTHON_AVAILABLE |

---

## Implementation Order

```
Pre-work A:  Commit Phase 6b transpiler                              (no cycle)
Pre-work B:  Ground truth baseline                                   (no cycle)
Pre-work C:  Test data design for Group D                            (no cycle)
Pre-work D:  Out-of-scope error messages                             (no cycle)
Cycle 1:     Default tolerance 0.5 → 0.1 Da                         (Group A)
Cycle 2:     Parser completeness (MATCHCOUNT + < + division)         (Group A)
Cycle 3:     INTENSITYTICPERCENT                                     (Group B)
Cycle 4:     ANY wildcard + MASSDEFECT qualifier                     (Group B)
Cycle 5:     scanrangesum() aggregation                              (Group B)
Cycle 6:     aminoaciddelta() function                               (Group C)
Cycle 7:     peptide() function                                      (Group C)
Cycle 8:     Y-variable parser                                       (Group D)
Cycle 9:     Y-only transpiler                                       (Group D)
Cycle 10:    X+Y offset pattern transpiler                           (Group D)
Cycle 11:    Y expressions with X dependency                         (Group D)
Cycle 12:    Y + non-Y conditions + cross-level integration          (Group D)
Cycle 13:    OTHERSCAN qualifier                                     (Group E)
Cycle 14:    Error message quality                                   (Group F)
Cycle 15:    Integration ground truth (final expansion)              (Group F)
Cycle 16:    Kitchen-sink composition tests                          (Group F)
Cycle 17:    Performance regression + degenerate edge cases          (Group F)
```

4 pre-work steps + 17 RED/GREEN/REFACTOR cycles. Build + test after every cycle.
Transpiler snapshots captured after Groups A, B, C, D, and E.
Ground truth file grows incrementally after each group.

---

## Verification Strategy

1. **After every cycle**: `bash build.sh && bash run_tests.sh && ./build/release/extension/miint/tests`
2. **After Group A**: Verify all existing tests still pass with new 0.1 default + capture transpiler snapshot
3. **After Group B**: Capture transpiler snapshot + expand ground truth
4. **After Group C**: Capture transpiler snapshot + expand ground truth
5. **After Group D**: Run isotope pattern queries through both `massql()` and `mzml_isotope_pattern` macro, verify identical results. Capture transpiler snapshot + expand ground truth
6. **After Group E**: Capture transpiler snapshot + expand ground truth
7. **Cycle 15**: Full ground truth verification against Python MassQL
8. **Cycle 16**: Kitchen-sink composition queries — verify no transpiler bugs from feature interactions
9. **Final**: Run every example query from `scratch/MassQueryLanguage/tests/` through `massql()` and compare
