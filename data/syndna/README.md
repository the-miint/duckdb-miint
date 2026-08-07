# synDNA absolute-quantification fixtures and parity goldens

Committed inputs and expected outputs for the `absquant_*` functions, generated
**offline** by running [pysyndna](https://github.com/biocore/pysyndna). Nothing here is
produced at build or test time, and duckdb-miint never invokes pysyndna — these are
numbers only.

pysyndna implements the synDNA spike-in method of Zaramela et al. 2022 (mSystems 7(6)
e0044722, PMC9765022) for turning relative metagenomic read counts into absolute cell
counts. See `THIRD_PARTY_LICENSES.md` for attribution.

## Provenance

| | |
|---|---|
| pysyndna | `a64687d4fb37ef7939b1cef8406c0b9758ebb8d7` (BSD-3-Clause, © Amanda Birmingham) |
| generator | `gen_syndna_oracle.py` — dev-only, **not committed** (same convention as the scikit-learn / scipy / `ape` generators behind `data/simsurvey/` and `data/asr/`) |
| python | 3.11.15 |
| scipy | 1.17.1 — computes `linregress`, so it defines the fit goldens exactly |
| numpy | 2.4.6 |
| pandas | 3.0.5 |
| biom-format | 2.1.17 |
| pyyaml | 6.0.3 |

## Regenerating

```sh
conda create -n pysyndna-oracle -c conda-forge -y \
    python=3.11 biom-format numpy pandas pyyaml scipy scikit-learn pytest
git -C ~/dev/pysyndna checkout a64687d
conda run -n pysyndna-oracle pip install --no-deps -e ~/dev/pysyndna

# confirm the environment is sound before trusting it (expect 106 passed)
cd ~/dev/pysyndna && conda run -n pysyndna-oracle python -m pytest -q

conda run -n pysyndna-oracle python \
    ~/dev/duckdb-miint-localdocs/gen_syndna_oracle.py \
    ~/dev/duckdb-miint/data/syndna
```

The generator is idempotent: re-running produces byte-identical files, so `git diff`
after a regeneration should be empty. Floats are written with Python `repr`, the
shortest round-trippable form.

Set A inputs are not transcribed by hand — the generator imports pysyndna's own test-data
classes (`FitSyndnaModelsTestData`, `TestCalcCellCountsData`, `TestQuantOrfsData`), so the
committed inputs are provably the ones pysyndna tests against. The generator then asserts
its outputs against pysyndna's own published expectations before writing anything.

## Two conversions applied on the way in

Both are cosmetic — they change no value, and pysyndna produces identical results either
way — but mistaking them for discrepancies would waste an afternoon:

1. **Coverage is a fraction here, a percent in pysyndna.** pysyndna's fixture coverages
   are percents (`92.597…`) compared against `min_coverage = 10`. These files store
   `coverage / 100`, to be compared against `min_coverage = 0.1`. miint is
   fraction-only, which also lets `genome_coverage_per_sample` output feed straight in.
2. **`sample_name` → `sample_id`, `ogu_id` / `syndna_id` / `ogu_orf_id` → `feature_id`.**
   pysyndna uses Qiita metadata names; miint uses `sample_id` / `feature_id` everywhere.
   Non-identifier parameter columns keep pysyndna's names verbatim.

Count and oracle tables are stored **as dense as pysyndna's own output is** — which is not
the same as "every cell". Two different things happen in pysyndna, and the distinction is
what the parity idiom relies on:

- A feature that survives the coverage filter in *at least one* sample gets an explicit
  `0.0` for the samples where it failed, because `calc_ogu_cell_counts_biom` pivots and
  then `fillna(0)`. `cellsb_oracle.csv`'s `(f_only1, s2)` is exactly this.
- A feature that fails in *every* sample (or a sample that fails everywhere) never
  reaches the pivot at all and is **genuinely absent** — not zeroed — on pysyndna's side
  too. `f_neither` and `sallcov` are these.

miint omits zero-valued rows in both situations (long-form sparse invariant). Because the
first case is dense in the golden, the parity tests can assert
`COALESCE(mine.value, 0) = gold.value` and so *prove* equivalence rather than assume it;
because the second case is absent on both sides, a `FULL JOIN` still catches any row
miint invents.

## Set A — pysyndna's own fixtures

Traceable to `docs/absolute_quant_example.xlsx` in the pysyndna repo and to the
Zaramela `SynDNA_saliva_samples_analysis.ipynb` notebook. This is the primary parity set,
because its provenance is external to both projects.

| File | Rows | Contents |
|---|---|---|
| `fit_counts.csv` | 20 | `(sample_id, feature_id, value)` — 2 samples × 10 synDNAs |
| `fit_concentrations.csv` | 10 | `(feature_id, syndna_indiv_ng_ul)` |
| `fit_params.csv` | 2 | `(sample_id, mass_syndna_input_ng)` |
| `fit_models_oracle.csv` | 2 | fitted models; `min_syndna_counts = 50`, contributing fraction `1.0` |
| `cells_counts.csv` | 26 | 2 samples × 13 OGUs, dense |
| `cells_models.csv` | 2 | the models applied |
| `cells_coverage.csv` | 26 | `(sample_id, feature_id, coverage)`, fractions |
| `cells_lengths.csv` | 13 | `(feature_id, ogu_len_in_bp)` |
| `cells_params.csv` | 2 | all four metrics' denominators present |
| `cells_oracle.csv` | 88 | `(metric, sample_id, feature_id, value)`, all 4 metrics, `min_coverage = 0.1`, `min_rsquared = 0.8` |
| `orf_counts.csv` | 20 | 2 samples × 10 OGU+ORFs, dense |
| `orf_coords.csv` | 10 | `(feature_id, ogu_orf_start, ogu_orf_end)`; 4 are reverse-strand (`start > end`) |
| `orf_params.csv` | 2 | |
| `orf_oracle.csv` | 20 | copies of ORF ssRNA per gram of sample |

*Neisseria subflava* (7.9% coverage) and *Haemophilus influenzae* (1.4%) fall below
`min_coverage` in every sample and are absent from `cells_oracle.csv` — pysyndna's own
fixtures engineer them that way so the filter is always exercised.

## Set B — synthetic edge cases

Hand-chosen integers, no RNG, so there is no seed to record. Each fit sample targets one
branch of `scipy.stats.linregress`; each cells/ORF case targets one filter.

| File | Rows | Contents |
|---|---|---|
| `fitb_counts.csv` | 88 | 11 samples × 8 synDNAs |
| `fitb_concentrations.csv` | 8 | `b5`/`b7`/`b8` share a concentration, mirroring the real pools (`p126` and `p266` are both 1 ng/µL) |
| `fitb_params.csv` | 11 | includes one NULL and one negative mass |
| `fitb_models_oracle.csv` | 5 | only the samples that fit, at `min_syndna_counts = 50` |
| `fitb_boundary_models_oracle.csv` | 8 | same inputs at `min_syndna_counts = 20` — see "Strict-< boundaries" |
| `cellsb_*.csv` | 4–33 | 5 samples × 4 OGUs |
| `cellsb_boundary_oracle.csv` | 45 | same inputs at `min_coverage = 0.02`, `min_rsquared = 0.25` |
| `orfb_*.csv` | 3–9 | 3 samples × 3 ORFs |

## The Student-t survival function grid

`studentt_sf_oracle.csv` — `(t, df, sf)`, 521 rows — stands apart from Sets A and B.
pysyndna is not involved: it comes straight from `scipy.stats.t.sf` (scipy is
BSD-3-Clause, so it is a license-compatible reference).

`pvalue` is the one regression field that cannot be built from DuckDB's `regr_*`
aggregates — it needs a Student-t survival function, i.e. a regularized incomplete beta,
which miint implements itself. The identity used is

```
sf(t, df) = 0.5 * I_{df/(df+t^2)}(df/2, 1/2)
```

which reproduces `scipy.stats.t.sf` bit-for-bit, so matching this grid is matching what
scipy actually computes rather than approximating it.

The grid exists because the fit oracles between them hold only a handful of `pvalue`s,
all at `df <= 6` — far too thin to validate new numerics. It spans `df` 1…1000 and `t`
from 0 through 1e10, in both signs, plus one extra row per `df` at
`t = sqrt(df / 2e-20)`: the value scipy's `TINY = 1e-20` guard produces at a perfect fit
(`r = ±1`). Those are the hardest inputs the function will ever see, and the reason
`TINY` exists — without it the t statistic divides by zero. At `df = 1` it gives
`sf ≈ 4.5e-11`, i.e. `pvalue ≈ 9.0e-11`, a genuine value rather than a floor.

Rows are unique on `(t, df)`. That needs saying because at `df = 2` the perfect-fit
extreme is `sqrt(2/2e-20) = 1e10` **exactly** and collides with the `1e10` grid point;
the generator deduplicates, and `absquant_fixtures.test` asserts the uniqueness.

`fitb` sample intent — six of the eleven deliberately produce **no** model:

| Sample | Outcome |
|---|---|
| `sA`, `sB` | ordinary noisy fits (non-zero `stderr` / `pvalue`) |
| `sPERFECT` | `rvalue == 1.0` exactly — pins scipy's `TINY = 1e-20` guard and the `[-1,1]` clamp |
| `sZERO` | one zero count, so that row is excluded and the fit runs on `n-1` points |
| `sTWO` | exactly 2 usable points — scipy's `n == 2` special case (`pvalue` 0.0, both stderrs hard-set to 0.0) |
| `sCONSTY` | constant y → `rvalue` is NaN → sample dropped |
| `sONE` | 1 usable point → every field NaN → sample dropped |
| `sNONE` | 0 usable points → `linregress` never called → sample dropped |
| `sFLAT` | constant x → scipy raises → sample dropped |
| `sNULL` | NULL `mass_syndna_input_ng` → filtered before fitting |
| `sNEG` | negative `mass_syndna_input_ng` → filtered before fitting |

`b6` totals 20 reads across all samples, so `min_syndna_counts = 50` drops that feature
globally — and drops it *before* any sample filtering, which is the order pysyndna's code
implements (its docstring says otherwise; the code is authoritative).

`cellsb` covers the cases the Set A fixture cannot:

- `f_only1` passes coverage in `s1` but not `s2`, so `(f_only1, s2)` is a **dense zero** —
  the case that proves omitting zero rows loses nothing.
- `f_zero` has a zero read count in `s2` yet passes coverage; the result is `0`. The zero
  does **not** reach `log10` on either side: `calc_ogu_cell_counts_biom` converts with
  `to_dataframe(dense=False)`, whose fill value is NaN, so biom's unstored zero arrives as
  NaN, stays NaN down the chain, and is zeroed by the closing `fillna(0)`. miint's reader
  drops the zero-valued cell outright. The golden `0.0` is matched by `COALESCE`, not by
  reproducing an infinity.
- `f_neither` fails coverage everywhere and is absent entirely.
- `slowr2` has `r² = 0.25 < 0.8` and is dropped from **every** metric.
- `snull` has NULL in `calc_mass_sample_aliquot_input_g` only, so it is dropped from
  `cells_per_g_of_sample` and **survives the other three metrics** — the NaN filter
  applies to the requested metric's denominator, not to all of them.
- `sallcov` fails coverage on all four of its OGUs, so the sample never reaches the
  per-sample calculation. That is a different code path from `f_neither` (one OGU failing
  across all samples), and one an implementation could plausibly crash on.

## Strict-`<` boundaries

All three of pysyndna's thresholds compare with a strict `<`, so a value sitting exactly
*on* the threshold must pass:

| Threshold | Source | Boundary fixture |
|---|---|---|
| `min_syndna_counts` | `fit_syndna_models.py:493` — `sum(axis=1) < min_sample_counts` | `b6` totals exactly 20 reads; `fitb_boundary_models_oracle.csv` was generated at `min_syndna_counts = 20`, so `b6` is retained — which in turn lets `sCONSTY`, `sFLAT` and `sONE` fit, hence 8 models rather than 5 |
| `min_coverage` | `calc_cell_counts.py:602` — `coverage < min_coverage` | `f_only1` has coverage exactly `0.02` in `s2`; `cellsb_boundary_oracle.csv` used `min_coverage = 0.02`, so that cell holds a real value instead of the dense `0.0` it has otherwise |
| `min_rsquared` | `calc_cell_counts.py:680` — `r_squared < min_rsquared` | `slowr2` has `rvalue = 0.5`, so `r²` is exactly `0.25`; the boundary oracle used `min_rsquared = 0.25`, so the sample survives the gate that drops it everywhere else |

These exist because the ordinary fixtures only ever land *clearly* on one side of each
threshold. Without them, a reimplementation that wrote `<=` for any of the three would
satisfy every other assertion in the suite and diverge only on data nobody committed.

The two float thresholds are exact, not lucky: `2.0/100.0` and the literal `0.02` are
both correctly rounded to the same double, and `0.5² = 0.25` is exactly representable.

## Tolerances

`test/sql/absquant_fixtures.test` checks only that these files load and are structurally
intact. The numeric parity tests live with the functions they test and use the
codebase's standard idiom — a `FULL JOIN` plus
`abs(mine - gold) > tol * (1.0 + abs(gold))`, asserted to return 0 rows:

| Golden | Tolerance | Why |
|---|---|---|
| `fit_models_oracle`, `fitb_models_oracle`, `fitb_boundary_models_oracle` | `1e-11` | pysyndna's public API truncates to 12 decimal places, so the golden itself carries only ~1e-12 absolute precision |
| `cells_oracle`, `cellsb_oracle`, `cellsb_boundary_oracle`, `orf_oracle`, `orfb_oracle` | `1e-9` | untruncated doubles; values reach ~1e13, so the bound must be relative |

## Two pysyndna behaviors deliberately not reproduced

Both are divisions by a zero that pysyndna's parameter screen structurally cannot catch:
zero is neither NaN nor negative, so it passes the `< 0` test and reaches the arithmetic.

- A sample whose `total_biological_reads_r1r2` is zero makes the read-fraction step divide
  by zero.
- A sample whose `calc_mass_sample_aliquot_input_g` is zero makes the final per-gram step
  divide by zero.

Either way pysyndna emits `inf` for every feature in that sample, with no warning. No
fixture here encodes them, because miint rejects such input with an error instead — an
infinite copy count is not a usable answer and must not reach a feature table silently.
