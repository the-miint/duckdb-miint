#!/usr/bin/env python3
"""Generate the committed SciPy parity goldens for miint's ks_2samp (issue #218).

    conda run -n <env-with-scipy> python test/scripts/generate_ks2samp_oracle.py

Writes, under data/ks2samp/:

  ks_2samp_fixture.csv       case_id, sample, idx, value   -- the two samples per case
  ks_2samp_oracle.csv        case_id, n1, n2, method, statistic, pvalue
  ks_2samp_large_oracle.csv  case_id, n1, off1, stride1, n2, off2, stride2, method, statistic, pvalue

Consumed by test/sql/ks_2samp_parity.test, which needs no SciPy installed.

This script IS committed on purpose. The small grid is drawn from a seeded NumPy
generator, and a seeded draw sequence cannot be reconstructed from a prose
description of the grid -- so without the script a reviewer facing a SciPy bump has
no way to tell a real regression from a different set of random draws. Running it
against SciPy 1.18.0 must reproduce the committed CSVs byte for byte; if it does
not, the difference is the finding.

The large-n grid deliberately carries NO fixture rows. Its inputs are arithmetic
ranges that SQL can rebuild exactly (`range(n) * stride + offset`, cast to
DOUBLE -- every integer below 2^53 is exact in a double), which is what lets the
parity check reach n = 10000 without committing 40000 rows of CSV.

DO NOT adjust an expected value to make a test pass. A disagreement means the
implementation changed.
"""
import csv
import os
import warnings

import numpy as np
import scipy
from scipy.stats import ks_2samp

OUT = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "data", "ks2samp")

# --------------------------------------------------------------------------------
# Small grid. Every case is an input class that could hide a real bug rather than a
# filler draw -- see data/ks2samp/README.md for what each one is for. The RNG is
# consumed in exactly this order, so do not reorder or insert cases above an
# existing one without regenerating and re-reviewing every later case.
# --------------------------------------------------------------------------------
rng = np.random.default_rng(20260817)


def curve(n, lo, hi, jitter):
    """A cumulative-coverage-shaped sample: monotone, in [0,1], rounded so it ties."""
    base = np.linspace(lo, hi, n) + rng.normal(0, jitter, n)
    return np.round(np.clip(np.sort(base), 0.0, 1.0), 3)


cases = {}
cases["identical_small"] = ([1.0, 2.0, 3.0], [1.0, 2.0, 3.0])
cases["identical_tied"] = ([1.0, 1.0, 2.0, 2.0, 3.0], [1.0, 1.0, 2.0, 2.0, 3.0])
cases["all_same_value"] = ([5.0] * 4, [5.0] * 6)
cases["disjoint_2_2"] = ([1.0, 2.0], [10.0, 20.0])
cases["disjoint_5_5"] = ([1.0, 2.0, 3.0, 4.0, 5.0], [6.0, 7.0, 8.0, 9.0, 10.0])
cases["disjoint_50_50"] = (list(np.arange(50.0)), list(np.arange(50.0) + 1000.0))
cases["disjoint_100_100"] = (list(np.arange(100.0)), list(np.arange(100.0) + 1e6))
cases["h1_overshoot_5_5"] = ([1.0, 2.0, 3.0, 4.0, 5.0], [1.5, 2.5, 3.5, 4.5, 5.5])
cases["single_1_1_equal"] = ([3.0], [3.0])
cases["single_1_1_diff"] = ([3.0], [4.0])
cases["single_1_5"] = ([3.0], [1.0, 2.0, 4.0, 5.0, 6.0])
cases["normal_20_20"] = (list(rng.normal(size=20)), list(rng.normal(size=20)))
cases["normal_20_30"] = (list(rng.normal(size=20)), list(rng.normal(size=30)))
cases["shifted_20_20"] = (list(rng.normal(size=20)), list(rng.normal(loc=1.5, size=20)))
cases["coprime_7_13"] = (list(rng.normal(size=7)), list(rng.normal(size=13)))
cases["coprime_3_11"] = (list(rng.normal(size=3)), list(rng.normal(size=11)))
cases["common_factor_6_9"] = (list(rng.normal(size=6)), list(rng.normal(size=9)))
cases["common_factor_12_30"] = (list(rng.normal(size=12)), list(rng.normal(size=30)))
cases["lopsided_2_97"] = (list(rng.normal(size=2)), list(rng.normal(size=97)))
cases["lopsided_97_2"] = (list(rng.normal(size=97)), list(rng.normal(size=2)))
cases["ties_alphabet4_15_15"] = (list(rng.integers(0, 4, 15).astype(float)), list(rng.integers(0, 4, 15).astype(float)))
cases["ties_alphabet4_9_21"] = (list(rng.integers(0, 4, 9).astype(float)), list(rng.integers(0, 4, 21).astype(float)))
cases["ties_two_values_20_20"] = (
    list(rng.integers(0, 2, 20).astype(float)),
    list(rng.integers(0, 2, 20).astype(float)),
)
cases["extreme_ties_20_20"] = ([0.0] * 20, [1.0] * 20)
cases["one_side_constant"] = ([1.0, 1.0, 1.0, 1.0], [0.0, 1.0, 2.0, 3.0])
cases["negatives_and_zero"] = ([-3.5, -1.0, 0.0, 0.0, 2.25], [-2.0, 0.0, 1.5, 4.0])
cases["curve_30_30"] = (list(curve(30, 0.05, 0.95, 0.02)), list(curve(30, 0.10, 0.99, 0.02)))
cases["curve_12_47"] = (list(curve(12, 0.02, 0.60, 0.03)), list(curve(47, 0.05, 0.98, 0.03)))
cases["curve_dominating"] = (list(curve(25, 0.01, 0.40, 0.01)), list(curve(25, 0.55, 0.99, 0.01)))
cases["curve_identical"] = (lambda c: (list(c), list(c)))(curve(18, 0.03, 0.92, 0.02))
cases["moderate_200_200"] = (list(rng.normal(size=200)), list(rng.normal(loc=0.35, size=200)))
cases["moderate_150_211"] = (list(rng.normal(size=150)), list(rng.normal(loc=0.2, size=211)))

# --------------------------------------------------------------------------------
# Large-n grid. Spec is (n1, off1, stride1, n2, off2, stride2), meaning
# sample = range(n) * stride + off, as DOUBLEs. No RNG and no fixture rows -- SQL
# rebuilds these exactly.
#
# Without them the committed grid topped out at n = 211 while ks_2samp answers up to
# 10000, so 98% of the answerable range had no SciPy evidence behind the docs' parity
# claim. Offsets are chosen to keep D small, which keeps the lattice band narrow and
# each case in the low milliseconds; a large D at n = 10000 is the ~0.5s worst case and
# is covered by the C++ closed-form test instead. The stride is what allows an UNEQUAL
# pair to still have a small D -- range(3000)*3 spans the same support as range(9000).
#
# Every case must keep its p-value in NORMAL double range, and the loop below raises if
# one does not. The reason is NOT that a subnormal result is inaccurate -- miint is
# bit-identical to SciPy at 6.6e-318 for equal-sized disjoint samples. It is that a
# large n together with a large D underflows the sweep's INTERMEDIATE probabilities,
# and there the two implementations diverge: 9000 vs 5000 at D = 0.411 gives 1.86e-321
# against SciPy's 8.90e-321, a factor of ~4.8. A fully underflowed 0.0 is excluded for a
# second reason: it would make the relative-tolerance parity check degenerate into
# comparing two zeros.
large = {
    "gen_10000_10000_shift5": (10000, 0, 1, 10000, 5, 1),
    "gen_10000_10000_shift500": (10000, 0, 1, 10000, 500, 1),
    "gen_10000_10000_shift1500": (10000, 0, 1, 10000, 1500, 1),
    "gen_5000_5000_shift50": (5000, 0, 1, 5000, 50, 1),
    "gen_9000_3000_stride3": (9000, 0, 1, 3000, 0, 3),
}

fixture_rows, oracle_rows, large_rows, fallbacks = [], [], [], []

for case_id in sorted(cases):
    a, b = cases[case_id]
    a = [float(x) for x in a]
    b = [float(x) for x in b]
    for idx, v in enumerate(a):
        fixture_rows.append((case_id, "a", idx, repr(v)))
    for idx, v in enumerate(b):
        fixture_rows.append((case_id, "b", idx, repr(v)))
    for method in ("auto", "exact"):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            r = ks_2samp(a, b, method=method)
            if any("unsuccessful" in str(x.message) for x in w):
                fallbacks.append((case_id, method))
        oracle_rows.append((case_id, len(a), len(b), method, repr(float(r.statistic)), repr(float(r.pvalue))))

SMALLEST_NORMAL = 2.2250738585072014e-308
for case_id in sorted(large):
    n1, off1, s1, n2, off2, s2 = large[case_id]
    a = np.arange(n1, dtype=float) * s1 + off1
    b = np.arange(n2, dtype=float) * s2 + off2
    for method in ("auto", "exact"):
        r = ks_2samp(a, b, method=method)
        pv_ = float(r.pvalue)
        if pv_ < SMALLEST_NORMAL:
            raise SystemExit(
                f"{case_id}/{method}: p={pv_!r} is subnormal or zero; the relative parity "
                f"check cannot discriminate there. Pick a smaller D for this case."
            )
        large_rows.append((case_id, n1, off1, s1, n2, off2, s2, method, repr(float(r.statistic)), repr(pv_)))

os.makedirs(OUT, exist_ok=True)
with open(os.path.join(OUT, "ks_2samp_fixture.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["case_id", "sample", "idx", "value"])
    w.writerows(fixture_rows)
with open(os.path.join(OUT, "ks_2samp_oracle.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["case_id", "n1", "n2", "method", "statistic", "pvalue"])
    w.writerows(oracle_rows)
with open(os.path.join(OUT, "ks_2samp_large_oracle.csv"), "w", newline="") as fh:
    w = csv.writer(fh)
    w.writerow(["case_id", "n1", "off1", "stride1", "n2", "off2", "stride2", "method", "statistic", "pvalue"])
    w.writerows(large_rows)

print(f"scipy {scipy.__version__}")
print(f"small cases {len(cases)}  fixture rows {len(fixture_rows)}  oracle rows {len(oracle_rows)}")
print(f"large cases {len(large)}  large oracle rows {len(large_rows)}")
print(f"scipy internal exact->asymp fallbacks: {fallbacks}")
pv = sorted(float(r[5]) for r in oracle_rows)
print(f"small p-value range: {pv[0]:.6e} .. {pv[-1]:.6e}")
for row in large_rows:
    if row[7] == "auto":
        print(f"  {row[0]:28s} n1={row[1]} n2={row[4]} D={float(row[8]):.6g} p={float(row[9]):.6g}")
