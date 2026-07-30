#!/usr/bin/env bash
#
# bump-duckdb-version.sh
#
# Retarget the extension at a new DuckDB release: move both submodules and
# rewrite every in-repo spelling of the version, then verify with
# check-duckdb-version.sh.
#
# The edits are purely mechanical but spread across six files, and two of them
# have drifted silently in the past. This script exists so the mechanical part
# stops being a human transcription exercise.
#
# It deliberately does NOT commit. Review `git diff`, build, test, then commit
# yourself — a version bump can break the C++ build (DuckDB's internal C++ API
# is not stable across releases), so a script that committed for you would be
# committing something unverified.
#
# Usage:
#   ./scripts/bump-duckdb-version.sh v1.5.6
#
# What it does NOT do — these live outside this repo, see docs/UPDATING.md:
#   - the publish host's env file (overrides DUCKDB_VERSION; silently rejects
#     every publish if left stale)
#   - the playground console's duckdb-wasm + EXT_VERSION pair
#   - duckdb/community-extensions' description.yml ref

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

WORKFLOW=".github/workflows/MainDistributionPipeline.yml"

die() { printf 'ERROR: %s\n' "$1" >&2; exit 1; }

# ---------------------------------------------------------------------------
# Arguments and preconditions
# ---------------------------------------------------------------------------

VER="${1:-}"
[[ -n "$VER" ]] || die "usage: $0 vX.Y.Z (e.g. v1.5.6)"
[[ "$VER" =~ ^v[0-9]+\.[0-9]+\.[0-9]+$ ]] || die "version must look like v1.5.6, got '$VER'"
BARE="${VER#v}"

# Refuse to run on top of staged work: this script rewrites tracked files, and
# mixing that into someone's staged changes makes the result hard to review.
if ! git diff --cached --quiet; then
    die "you have staged changes; commit or reset them first (this script edits tracked files)"
fi

printf 'Retargeting DuckDB -> %s\n\n' "$VER"

# ---------------------------------------------------------------------------
# Submodules
# ---------------------------------------------------------------------------

printf '== submodules ==\n'

[[ -e duckdb/.git ]] || die "duckdb submodule not initialised (git submodule update --init duckdb)"

git -C duckdb fetch --quiet --tags origin || die "could not fetch tags in duckdb/"
git -C duckdb rev-parse -q --verify "refs/tags/${VER}^{commit}" >/dev/null \
    || die "duckdb has no tag $VER (is the release out yet?)"

if ! git -C duckdb diff --quiet || ! git -C duckdb diff --cached --quiet; then
    die "duckdb/ has local modifications; refusing to move it"
fi
git -C duckdb checkout --quiet "$VER"
printf '  duckdb            -> %s (%s)\n' "$VER" "$(git -C duckdb rev-parse --short=10 HEAD)"

# extension-ci-tools tracks a release-SERIES branch (e.g. v1.5-variegata), which
# upstream keeps identical to the per-patch branch (v1.5.5).
#
# The branch comes from the workflow's `ci_tools_version` input, NOT from
# .gitmodules. `ci_tools_version` is what CI actually checks out, so it is the
# authority; .gitmodules has said `branch = main` while we tracked
# v1.5-variegata, and following that jumps the submodule to DuckDB's development
# tip. Deriving the codename from the version number is not possible either
# (v1.5.5 -> "variegata" is not computable), so read it from the workflow.
if [[ -e extension-ci-tools/.git ]]; then
    CI_BRANCH=$(grep -oE '^[[:space:]]*ci_tools_version:[[:space:]]*\S+' "$WORKFLOW" \
                | head -1 | awk '{print $2}')
    [[ -n "$CI_BRANCH" ]] || die "could not read ci_tools_version from $WORKFLOW"
    git -C extension-ci-tools fetch --quiet origin "$CI_BRANCH" \
        || die "could not fetch $CI_BRANCH in extension-ci-tools/"
    if ! git -C extension-ci-tools diff --quiet || ! git -C extension-ci-tools diff --cached --quiet; then
        die "extension-ci-tools/ has local modifications; refusing to move it"
    fi
    git -C extension-ci-tools checkout --quiet FETCH_HEAD
    printf '  extension-ci-tools -> %s (%s) on %s\n' \
        "$(git -C extension-ci-tools rev-parse --short=10 HEAD)" \
        "$(git -C extension-ci-tools log -1 --format=%s | cut -c1-52)" \
        "$CI_BRANCH"
    printf '                       (verify this branch has been bumped for %s upstream)\n' "$VER"
else
    printf '  extension-ci-tools: not initialised, skipped\n'
fi

# ---------------------------------------------------------------------------
# In-repo version spellings
#
# `must` edits are functional — if the pattern stops matching, the file has been
# restructured and silently skipping it would ship a stale version. Fail loud.
# `prose` edits are comments; a reworded comment is not a bug, so warn only.
# ---------------------------------------------------------------------------

printf '\n== files ==\n'

edit() {
    local mode="$1" file="$2" pattern="$3" replacement="$4" label="$5"
    [[ -f "$file" ]] || { [[ "$mode" == must ]] && die "$file not found"; return; }

    if ! grep -qE "$pattern" "$file"; then
        if [[ "$mode" == must ]]; then
            die "$label: pattern no longer matches in $file — the file was restructured; fix this script (and check-duckdb-version.sh) before relying on it"
        fi
        printf '  warn  %s: pattern not found in %s (reworded?), skipped\n' "$label" "$file"
        return
    fi

    local before after
    before=$(cat "$file")
    perl -0pi -e "s/${pattern}/${replacement}/g" "$file"
    after=$(cat "$file")

    if [[ "$before" == "$after" ]]; then
        printf '  ok    %s already at %s\n' "$label" "$VER"
    else
        printf '  edit  %s\n' "$label"
    fi
}

# duckdb_version inputs (both the distribution build and the code-quality job)
edit must "$WORKFLOW" \
    '(duckdb_version:\s*)v[0-9]+\.[0-9]+\.[0-9]+' \
    "\${1}${VER}" \
    "$WORKFLOW duckdb_version"

# Artifact name in the commented-out verify-wasm job
edit must "$WORKFLOW" \
    'miint-v[0-9]+\.[0-9]+\.[0-9]+-extension-' \
    "miint-${VER}-extension-" \
    "$WORKFLOW wasm artifact name"

edit must "Dockerfile" \
    '(ARG DUCKDB_VERSION=)[0-9]+\.[0-9]+\.[0-9]+' \
    "\${1}${BARE}" \
    "Dockerfile ARG"

# Only the assignment and the example env file. NOT the tag-glob examples
# elsewhere in that script, which illustrate miint's own release tags
# (v[0-9]+*) and have nothing to do with the DuckDB version.
edit must "scripts/cron-publish-extension.sh" \
    '(: "\$\{DUCKDB_VERSION:=)v[0-9]+\.[0-9]+\.[0-9]+(\}")' \
    "\${1}${VER}\${2}" \
    "cron publish default"

edit must "scripts/cron-publish-extension.sh" \
    '(#\s+DUCKDB_VERSION=)v[0-9]+\.[0-9]+\.[0-9]+' \
    "\${1}${VER}" \
    "cron example env file"

edit must "python/pyproject.toml" \
    '(duckdb==)[0-9]+\.[0-9]+\.[0-9]+' \
    "\${1}${BARE}" \
    "python CLI pin"

edit prose "CMakeLists.txt" \
    '(the version DuckDB )v[0-9]+\.[0-9]+\.[0-9]+( pins)' \
    "\${1}${VER}\${2}" \
    "CMakeLists emsdk comment"

# ---------------------------------------------------------------------------
# Verify, then hand back
# ---------------------------------------------------------------------------

printf '\n== verify ==\n'
bash "$SCRIPT_DIR/check-duckdb-version.sh" || die "consistency check failed after bumping — see above"

printf '\n== changed ==\n'
git diff --stat

cat <<EOF

Next:
  1. bash build.sh                       # the C++ API is not stable across releases
  2. bash run_tests.sh                   # SQL suite
  3. ./build/release/extension/miint/tests   # C++ suite
  4. conda run -n duckdb-format make format-check
  5. WASM is only built by CI on tag pushes, so verify it locally:
       bash scripts/build_wasm.sh        # build + unresolved-import check
       bash scripts/test_wasm_load.sh    # true ABI load test
     If build/wasm_shell/ is left over from an older DuckDB, test_wasm_load.sh
     reuses it and fails with a version-lock error that looks like a real
     incompatibility but is not. Clear it first.
  6. Commit, open a PR, and see docs/UPDATING.md for the cross-repo follow-ups
     (publish-host env file, playground version pair, community-extensions ref).
EOF
