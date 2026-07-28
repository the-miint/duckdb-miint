#!/usr/bin/env python
"""Generate parity-oracle fixtures + goldens for the `procrustes` table function.

The C++ implementation (src/procrustes_core.cpp) is ported/adapted from SciPy's
`scipy.spatial.procrustes` + `scipy.linalg.orthogonal_procrustes` (BSD-3-Clause;
Copyright (c) 2001-2002 Enthought, Inc.; 2003-, SciPy Developers). This script
uses SciPy directly as the independent numeric oracle for the FULL-overlap case,
and a NumPy port of the "deconstructed"/partial technique (fit the similarity
transform on a paired subset, apply it to all rows) for the PARTIAL case.

The partial technique is the author's own work (Daniel McDonald,
qiime2/q2-diversity#338 `_partial_procrustes` / `_deconstructed_procrustes`),
reimplemented here; its `_deconstructed_procrustes` base derives from SciPy, so
SciPy attribution applies to that math. When `--check-q2` is passed and a q2
environment is active, the NumPy partial is cross-checked against the shipped
`q2_diversity._procrustes._partial_procrustes` to confirm the port is faithful.

Outputs (under data/procrustes/), all deterministic given the fixed seed:
  <case>.reference.csv   sample_id,axis,coordinate   (input ordination A)
  <case>.other.csv       sample_id,axis,coordinate   (input ordination B)
  <case>.pairing.csv     reference_id,other_id       (partial cases only)
  <case>.expected.csv    matrix,sample_id,axis,coordinate  (standardized A + fitted B)
  <case>.m2.csv          m2                          (disparity on the paired points)
  oracle.manifest.txt    scipy/numpy versions + case descriptions
  goldens.sha256         sha256 of every *.expected.csv / *.m2.csv

Axes are 0-indexed to match the `pcoa`/`unifrac_pcoa` output convention.

Usage:
  conda run -n base python test/scripts/gen_procrustes_oracle.py [--check-q2]
Regenerate goldens intentionally (run_tests.sh will re-pin the SHA sidecar):
  MIINT_PROCRUSTES_REGENERATE=1 conda run -n base python test/scripts/gen_procrustes_oracle.py
"""
import argparse
import hashlib
import os
import sys

import numpy as np
from scipy.spatial import procrustes as scipy_procrustes
from scipy.linalg import orthogonal_procrustes

HERE = os.path.dirname(os.path.abspath(__file__))
OUTDIR = os.path.normpath(os.path.join(HERE, "..", "..", "data", "procrustes"))
SEED = 42


def _case_rng(name):
    """Independent, reorder-stable RNG for one case.

    Seeded from (SEED, hash(name)) rather than a single generator threaded
    through every case, so adding, removing, or reordering a case never shifts
    another case's goldens (which would otherwise re-pin unrelated SHAs). The
    name hash is a stable SHA-256 digest, not Python's salted hash(), so the
    stream is identical across runs and machines.
    """
    tag = int.from_bytes(hashlib.sha256(name.encode()).digest()[:8], "little")
    return np.random.default_rng(np.random.SeedSequence([SEED, tag]))


# ── the deconstructed / partial technique (NumPy port; q2#338, SciPy base) ──
def _deconstructed(mtx1, mtx2):
    """Fit a similarity transform (translate, uniform scale, orthogonal R).

    Returns the six parameters instead of applying them, so the fit computed on
    a paired subset can be applied to a larger set. Derived from scipy.procrustes.
    """
    mtx1 = np.asarray(mtx1, dtype=np.float64).copy()
    mtx2 = np.asarray(mtx2, dtype=np.float64).copy()
    t1 = mtx1.mean(0)
    t2 = mtx2.mean(0)
    mtx1 -= t1
    mtx2 -= t2
    n1 = np.linalg.norm(mtx1)
    n2 = np.linalg.norm(mtx2)
    if n1 == 0 or n2 == 0:
        raise ValueError("input matrices must contain >1 unique points")
    mtx1 /= n1
    mtx2 /= n2
    R, s = orthogonal_procrustes(mtx1, mtx2)
    return t1, t2, n1, n2, R, s


def _partial(ref, other, ref_pair_idx, other_pair_idx):
    """Fit on the paired rows, apply the transform to ALL rows of both matrices.

    ref, other: (n, d) arrays. *_pair_idx: integer row indices of the paired
    (anchor) subset. Returns (ref_transformed, other_transformed).
    """
    t1, t2, n1, n2, R, s = _deconstructed(ref[ref_pair_idx], other[other_pair_idx])
    ref_t = (ref - t1) / n1
    other_t = ((other - t2) / n2) @ R.T * s
    return ref_t, other_t


def _m2(ref_t, other_t, ref_pair_idx, other_pair_idx):
    """Disparity M^2 = sum of squared pointwise diffs over the paired points."""
    diff = ref_t[ref_pair_idx] - other_t[other_pair_idx]
    return float(np.sum(diff * diff))


# ── fixture construction ──
def _similarity_image(x, rng, scale, reflect=False, noise=0.0):
    """Return a translated + scaled + rotated (optionally reflected + noised)
    image of x, so procrustes has a real transform to recover."""
    d = x.shape[1]
    # random orthogonal matrix via QR
    Q, _ = np.linalg.qr(rng.standard_normal((d, d)))
    if reflect:
        Q[:, 0] *= -1.0  # force det = -1
    y = scale * (x @ Q) + rng.standard_normal(d) * 3.0
    if noise:
        y = y + rng.standard_normal(x.shape) * noise
    return y


def _write_long(path, ids, mat):
    with open(path, "w") as fh:
        fh.write("sample_id,axis,coordinate\n")
        for i, sid in enumerate(ids):
            for axis in range(mat.shape[1]):
                fh.write("%s,%d,%.17g\n" % (sid, axis, mat[i, axis]))


def _write_expected(path, ref_ids, ref_t, other_ids, other_t):
    with open(path, "w") as fh:
        fh.write("matrix,sample_id,axis,coordinate\n")
        for label, ids, mat in (
            ("reference", ref_ids, ref_t),
            ("other", other_ids, other_t),
        ):
            for i, sid in enumerate(ids):
                for axis in range(mat.shape[1]):
                    fh.write("%s,%s,%d,%.17g\n" % (label, sid, axis, mat[i, axis]))


def _write_m2(path, m2):
    with open(path, "w") as fh:
        fh.write("m2\n%.17g\n" % m2)


def _write_pairing(path, pairs):
    with open(path, "w") as fh:
        fh.write("reference_id,other_id\n")
        for a, b in pairs:
            fh.write("%s,%s\n" % (a, b))


# ── cases ──
def case_full(name, n, d, reflect=False):
    """Full-overlap case: same ids in both; SciPy is the oracle."""
    rng = _case_rng(name)
    ref = rng.standard_normal((n, d))
    other = _similarity_image(ref, rng, scale=2.5, reflect=reflect, noise=0.05)
    ids = ["S%d" % i for i in range(n)]
    # SciPy oracle (both standardized; disparity over all points == all paired).
    mtx1, mtx2, disparity = scipy_procrustes(ref, other)
    _write_long(os.path.join(OUTDIR, "%s.reference.csv" % name), ids, ref)
    _write_long(os.path.join(OUTDIR, "%s.other.csv" % name), ids, other)
    _write_expected(os.path.join(OUTDIR, "%s.expected.csv" % name), ids, mtx1, ids, mtx2)
    _write_m2(os.path.join(OUTDIR, "%s.m2.csv" % name), disparity)
    # Self-check: our deconstructed(all-paired) must equal SciPy's full procrustes.
    idx = np.arange(n)
    ref_t, other_t = _partial(ref, other, idx, idx)
    assert np.allclose(ref_t, mtx1, atol=1e-9), "%s: ref frame != scipy" % name
    assert np.allclose(other_t, mtx2, atol=1e-9), "%s: other frame != scipy" % name
    return "full  n=%d d=%d reflect=%s" % (n, d, reflect)


def case_partial(name, n_anchor, n_extra, d):
    """Partial case: `other` has the anchors (distinct ids via pairing) + extras
    in a different frame; fit on anchors, apply to all."""
    rng = _case_rng(name)
    # reference = anchors only; other = anchors' image + extra points.
    ref = rng.standard_normal((n_anchor, d))
    extra = rng.standard_normal((n_extra, d))
    both = np.vstack([ref, extra])
    other = _similarity_image(both, rng, scale=1.7, noise=0.02)
    ref_ids = ["A%d" % i for i in range(n_anchor)]
    other_ids = ["B%d" % i for i in range(n_anchor)] + ["X%d" % i for i in range(n_extra)]
    pairs = list(zip(ref_ids, other_ids[:n_anchor]))  # A_i <-> B_i (distinct ids)
    ref_idx = np.arange(n_anchor)
    other_idx = np.arange(n_anchor)  # anchors are the first n_anchor rows of `other`
    ref_t, other_t = _partial(ref, other, ref_idx, other_idx)
    m2 = _m2(ref_t, other_t, ref_idx, other_idx)
    _write_long(os.path.join(OUTDIR, "%s.reference.csv" % name), ref_ids, ref)
    _write_long(os.path.join(OUTDIR, "%s.other.csv" % name), other_ids, other)
    _write_pairing(os.path.join(OUTDIR, "%s.pairing.csv" % name), pairs)
    _write_expected(
        os.path.join(OUTDIR, "%s.expected.csv" % name),
        ref_ids,
        ref_t,
        other_ids,
        other_t,
    )
    _write_m2(os.path.join(OUTDIR, "%s.m2.csv" % name), m2)
    _maybe_check_q2(name, ref, other, ref_ids, other_ids, ref_t, other_t)
    return "partial anchors=%d extra=%d d=%d" % (n_anchor, n_extra, d)


_CHECK_Q2 = False


def _maybe_check_q2(name, ref, other, ref_ids, other_ids, ref_t, other_t):
    """Cross-check the NumPy partial against the shipped q2 implementation."""
    if not _CHECK_Q2:
        return
    import pandas as pd
    from q2_diversity._procrustes import _partial_procrustes

    df1 = pd.DataFrame(ref, index=ref_ids)
    df2 = pd.DataFrame(other, index=other_ids)
    q1, q2 = _partial_procrustes(df1, df2, ref_ids, other_ids[: len(ref_ids)])
    assert np.allclose(q1.values, ref_t, atol=1e-9), "%s: q2 ref mismatch" % name
    assert np.allclose(q2.values, other_t, atol=1e-9), "%s: q2 other mismatch" % name
    print("  [q2 cross-check OK] %s" % name)


def main():
    global _CHECK_Q2
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "--check-q2",
        action="store_true",
        help="cross-check NumPy partial vs q2_diversity (run in a q2 env)",
    )
    args = ap.parse_args()
    _CHECK_Q2 = args.check_q2

    os.makedirs(OUTDIR, exist_ok=True)
    # Each case seeds its own RNG from (SEED, case name) — see _case_rng.
    descriptions = []
    descriptions.append(("full_small", case_full("full_small", n=6, d=3)))
    descriptions.append(("full_reflect", case_full("full_reflect", n=8, d=3, reflect=True)))
    descriptions.append(("full_wide", case_full("full_wide", n=12, d=5)))
    descriptions.append(
        (
            "partial_anchored",
            case_partial("partial_anchored", n_anchor=6, n_extra=4, d=3),
        )
    )

    import scipy

    with open(os.path.join(OUTDIR, "oracle.manifest.txt"), "w") as fh:
        fh.write("procrustes parity oracle\n")
        fh.write("seed=%d (per-case: SeedSequence([SEED, sha256(case_name)]))\n" % SEED)
        fh.write("numpy=%s\nscipy=%s\n" % (np.__version__, scipy.__version__))
        fh.write("oracle: scipy.spatial.procrustes (full) + numpy deconstructed partial (q2#338)\n")
        for nm, desc in descriptions:
            fh.write("case %s: %s\n" % (nm, desc))

    # SHA-pin every golden file (expected coords + m2).
    goldens = sorted(f for f in os.listdir(OUTDIR) if f.endswith(".expected.csv") or f.endswith(".m2.csv"))
    with open(os.path.join(OUTDIR, "goldens.sha256"), "w") as fh:
        for f in goldens:
            h = hashlib.sha256(open(os.path.join(OUTDIR, f), "rb").read()).hexdigest()
            fh.write("%s  %s\n" % (h, f))

    print("wrote %d cases to %s" % (len(descriptions), OUTDIR))
    for nm, desc in descriptions:
        print("  %-18s %s" % (nm, desc))


if __name__ == "__main__":
    sys.exit(main())
