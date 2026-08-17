# `ks_2samp` SciPy parity golden

Fixtures and goldens for `test/sql/ks_2samp_parity.test`, which pins miint's
`ks_2samp` against `scipy.stats.ks_2samp` (issue #218).

micov's published KS statistics and p-values came from SciPy, so parity is an
acceptance criterion rather than a nicety: if these numbers move, previously
published results stop reproducing.

## Files

| File | Contents |
|---|---|
| `ks_2samp_fixture.csv` | `case_id, sample, idx, value` — the two samples per small case in long form, so SQL can rebuild them with `list(value)` |
| `ks_2samp_oracle.csv` | `case_id, n1, n2, method, statistic, pvalue` — SciPy's answers for the small grid |
| `ks_2samp_large_oracle.csv` | `case_id, n1, off1, stride1, n2, off2, stride2, method, statistic, pvalue` — the large-n grid, whose inputs are arithmetic ranges rather than fixture rows |

Generated with **SciPy 1.18.0** (conda env `miint-beta-oracle`), NumPy default
RNG seeded at `20260817`, by the committed generator:

```bash
conda run -n miint-beta-oracle python test/scripts/generate_ks2samp_oracle.py
```

The generator **is** committed, and running it must reproduce these CSVs byte for
byte. That matters because the small grid comes from a seeded NumPy generator: a
seeded draw sequence cannot be reconstructed from a prose description of the grid,
so without the script a reviewer facing a SciPy bump could not tell a real
regression from a different set of random draws. (An earlier version of this
directory followed `community_distances_parity.test` in leaving the generator out;
that is the minority pattern here — `test/scripts/` already holds six committed
fixture generators.) The test itself needs no SciPy on the machine.

## Grid

32 cases × 2 methods = 64 oracle rows. Every case is an input class that could
hide a real bug rather than a filler draw:

- **unequal and coprime `(n1, n2)`** — equal sizes hide index errors, and
  `gcd == 1` puts the lattice band's edges between rows
- **a common factor (`gcd > 1`)** — exercises the lcm/gcd reduction
- **lopsided, both orientations** (`2` vs `97` and `97` vs `2`) — a wide band per
  column, and asymmetry in the band index would show up here and nowhere else
- **`n == 1`** — degenerate ECDF
- **heavy ties from a small integer alphabet** — cumulative coverage curves are
  nothing *but* ties, and SciPy emits no warning on them, so a grid of continuous
  draws would pin only the untied code path
- **curve-shaped monotone input** — micov's actual payload: sorted proportions in
  `[0, 1]`, rounded so they tie
- **identical samples** — `D = 0`, so `pvalue` must be exactly `1`
- **disjoint samples** — `D = 1`, where the answer is the closed form
  `2 / C(n1+n2, n1)`. At `100` vs `100` that is `2.2e-59`, which is what pins
  *precision*: computing the p-value as `1 - P(stay inside)` cancels every
  significant digit away at that magnitude
- **`n1 == n2 == 5` with `D == 0.2`** (case `h1_overshoot_5_5`) — see below

The recorded p-values span `2.2e-59` to `1.0`.

## `h1_overshoot_5_5`, and why it is in the grid

This is the one case where SciPy's own exact routine fails. Its equal-sample-size
formula returns `1.0000000000000002`, `_attempt_exact_2kssamp` sees a value
outside `[0, 1]`, reports failure, and `ks_2samp` **silently falls back to its
asymptotic branch** (emitting a `RuntimeWarning`) before clipping the result to
`1.0`. The generator captures that fallback and it fires for both `method='auto'`
and `method='exact'`.

miint lands on `1.0` directly, without needing a fallback, so the two agree. The
case is committed to keep that agreement pinned.

It does **not** pin miint's clamping, and it is worth being precise about that:
miint's accumulated escape mass for `(5, 5, D=0.2)` is *exactly* `1.0`, so the
clamp is inert here. Mutation testing confirmed it — removing the clamp leaves
this case, and the whole parity grid, green. The clamp is pinned instead by
`test/cpp/test_KsTwoSample.cpp`, at `n1=7, n2=13` with `h=9/91`, where the raw sum
really does reach `1.0000000000000002`.

## `method`

Both `'auto'` and `'exact'` are recorded for every case, and they agree
everywhere — as they must, since SciPy's `'auto'` resolves to `'exact'` whenever
`max(n1, n2) <= 10000` and no case here is larger.

There are **no `method='asymp'` rows, and no rows above `max(n1, n2) = 10000`**.
miint raises rather than answering there: SciPy's asymptotic branch is not the
Kolmogorov series but the exact one-sample KS distribution at an effective sample
size, chosen per region among the Durbin-matrix, Pomeranz, and Pelz-Good
expansions, and approximating it would silently disagree with published values.
The ceiling is covered as an *error* test in `test/sql/ks_2samp.test` and
`test/cpp/test_KsTwoSample.cpp` instead of as a parity row. See issue #218.

Note the ceiling means two different things for the two methods, and only one of
them is SciPy's: for `'auto'` it is SciPy's own `MAX_AUTO_N`, while for `'exact'`
SciPy has no cap at all and the limit is miint's own cost bound.

## Large-n grid

Five cases at `n₁` of 5000–10000, with `n₂` down to 3000 (the unequal case is
9000 vs 3000), in `ks_2samp_large_oracle.csv`. They exist because
the small grid tops out at `n = 211` while `ks_2samp` answers to 10000, which left
98% of the answerable range with no SciPy evidence behind the parity claim in
`docs/diversity.md`.

Their inputs carry no fixture rows: each sample is `range(n) * stride + offset` cast
to DOUBLE, which SQL rebuilds exactly because every integer below 2^53 is exact in a
double. The stride is what lets an *unequal* pair still have a small `D` —
`range(3000) * 3` spans the same support as `range(9000)`. Offsets are chosen to keep
`D` small so the lattice band stays narrow and each case runs in milliseconds; the
large-`D` worst case at the ceiling (~0.5 s) is covered by the C++ closed-form test.

**Subnormal and zero p-values are deliberately excluded, and the generator refuses to
emit one.** Not because a subnormal *result* is inaccurate — miint is bit-identical to
SciPy at `6.6e-318` for equal-sized disjoint samples — but because a large `n` combined
with a large `D` underflows the sweep's *intermediate* probabilities, and there the two
implementations diverge: 9000 vs 5000 at `D = 0.411` gives `1.86e-321` against SciPy's
`8.90e-321`, a factor of ~4.8, since both accumulate in probability space in different
orders. A fully underflowed `0.0` is excluded for a second reason: it would make the
relative-tolerance parity check degenerate into comparing two zeros. The smallest
p-value in the grid is `1.7e-98`.

## Regenerating

```bash
conda run -n miint-beta-oracle python test/scripts/generate_ks2samp_oracle.py
git diff --stat data/ks2samp/          # expect NO diff against a pinned SciPy 1.18.0
./build/release/test/unittest "test/sql/ks_2samp_parity.test"
```

Update the SciPy version recorded here if it changes, and review every changed row
rather than accepting the diff wholesale. Do not adjust an expected value to make the
test pass — a disagreement means the implementation changed.
