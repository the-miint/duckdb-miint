# LBFGS++ — vendoring provenance

Header-only L-BFGS solver (Eigen-backed), used by the `mmvec` MAP fit
(`src/mmvec.cpp`) as the deterministic quasi-Newton optimizer. See
[`THIRD_PARTY_LICENSES.md`](../../THIRD_PARTY_LICENSES.md) for the licence record
and [`docs/internals/embedded-tools.md`](../../docs/internals/embedded-tools.md)
for how it is built.

## Upstream

- Repository: https://github.com/yixuan/LBFGSpp
- Version: **v0.4.0** (released 2025-04-20)
- Tag commit: `c524a407fb85b74807f53de5a3ca2ddbcc164e54`
- Release tarball: `https://github.com/yixuan/LBFGSpp/archive/refs/tags/v0.4.0.tar.gz`
  sha256 `39c4aaebd8b94ccdc98191d51913a31cddd618cc0869d99f07a4b6da83ca6254`
- License: MIT (`LICENSE.md`, vendored verbatim alongside the headers)

Vendored as tracked file copies rather than a submodule, matching
`ext/kseq++` and `ext/concurrentqueue`: header-only, no build step, and a
submodule would add a checkout for nine files.

## What was copied, and what was not

Only the **unbounded** solver is vendored — `mmvec`'s parameters are
unconstrained. Upstream's bounded solver (`LBFGSB.h`, `LBFGSpp/Cauchy.h`,
`LBFGSpp/SubspaceMin.h`) is deliberately omitted, as are the examples, tests,
CMake scripts and docs.

| Vendored file | Upstream path | sha256 |
|---|---|---|
| `LBFGS.h` | `include/LBFGS.h` | `262ec22f5dad37e4906d73602f99c3ff681a41687e6fe81d0d1c4ff7fc085354` |
| `Param.h` | `include/LBFGSpp/Param.h` | `34466883c92c1b994df551d39528ed3d7c2b7bf5501171ceee370020472b3cd0` |
| `BFGSMat.h` | `include/LBFGSpp/BFGSMat.h` | `67babe7b145a832e9cdb80da3d9c8dd17814c604b071d19378a2c8327a5a05dc` |
| `BKLDLT.h` | `include/LBFGSpp/BKLDLT.h` | `7fc07494b8d0b7fdf4940d02f1a37fbdd3b5aca1064516906c8cd6f0b4906cb3` |
| `LineSearchBacktracking.h` | `include/LBFGSpp/LineSearchBacktracking.h` | `e9f4f0f40a9f7ab475a22257cbc33d0b2c40e766b4a0659ebf780e8b7f1cab30` |
| `LineSearchBracketing.h` | `include/LBFGSpp/LineSearchBracketing.h` | `0d0450fce841bde50091cdb6d54da4e1dfd7e97941861a48d52f9bcdc86a7601` |
| `LineSearchNocedalWright.h` | `include/LBFGSpp/LineSearchNocedalWright.h` | `b342533ea242e7ae1cada096c3f294c5ca63f22b5608f8ce6b0961b129c7231e` |
| `LineSearchMoreThuente.h` | `include/LBFGSpp/LineSearchMoreThuente.h` | `014fa5ead60b0e18cb29eb0f0113c6d5b8fd4bf9d6eca2aa9d1d1ef76a5616af` |
| `LICENSE.md` | `LICENSE.md` | `de8b15e992a8c9665eecb89270530a040c53b1ed52aaed9a9ae86fc586e865c4` |

Every file is **unmodified** — the checksums above are upstream's own. Any local
change would break that property, so fixes belong upstream, not here.

## Layout: flattened, and why that needs no include-path change

Upstream splits the headers across two directories (`include/LBFGS.h` plus
`include/LBFGSpp/*.h`); they are vendored flat into `ext/LBFGSpp/`. This is the
same flattening `ext/kseq++/kseq++.hpp` applies to upstream's
`include/kseq++/kseq++.hpp`.

It works without patching any `#include` line and without adding an include
directory, because `ext` is already on the include path
(`include_directories(... ext ...)` in `CMakeLists.txt`):

- our code writes `#include "LBFGSpp/LBFGS.h"` → found under `ext/`;
- `LBFGS.h`'s own `#include "LBFGSpp/Param.h"` is not found next to `LBFGS.h`,
  so it falls back to the include path and resolves to `ext/LBFGSpp/Param.h`;
- `BFGSMat.h`'s same-directory `#include "Param.h"` / `#include "BKLDLT.h"`
  resolve directly.

## Eigen

`BFGSMat.h` includes `<Eigen/Core>` and `<Eigen/LU>`. Eigen is already a
dependency (vcpkg `eigen3`, MPL-2.0), so LBFGS++ adds no new one.

The build sets `EIGEN_MPL2_ONLY` globally (`CMakeLists.txt`) to keep miint's
Modified-BSD licence clean. Verified before anything depended on these headers:
`LBFGSSolver<double, LineSearchMoreThuente>` compiles and minimizes correctly at
`-std=c++20 -Wall -Wextra` against the vendored Eigen with `EIGEN_MPL2_ONLY`
defined. The `<Eigen/LU>` include is not a problem — the transitive include set
reaches no non-MPL2 Eigen file. (In Eigen 3.4.x no header trips the guard at
all: the one LGPL-derived file, `IterativeLinearSolvers/IncompleteLUT.h`, was
relicensed to MPL-2.0 with the original author's permission, and it is not
reached from `<Eigen/LU>` regardless.)

## Notes for callers

Two upstream behaviours are load-bearing for `mmvec`, both documented at the
call site in `src/mmvec.cpp`:

1. **`LBFGSParam::ftol` and `wolfe` are line-search constants**, not convergence
   tolerances — they are the Armijo and Wolfe parameters. SciPy's convergence
   `ftol` corresponds to `past` / `delta` instead.
2. **`minimize()` throws on line-search failure** where SciPy's L-BFGS-B returns
   a non-zero status, and may leave `x` mid-step. Callers that must survive
   non-convergence (real data does not reach a stationary point) need their own
   best-so-far snapshot.
