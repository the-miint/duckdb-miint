#!/usr/bin/env bash
#
# check-duckdb-version.sh
#
# Assert that every place this repo spells out its target DuckDB version agrees.
#
# The version is duplicated across workflows, the Dockerfile, the publish cron,
# and the Python CLI's dependency pin. Nothing kept them in sync, and they have
# drifted silently before: python/pyproject.toml sat at duckdb==1.4.4 while the
# extension targeted 1.5.4 — a pin that could not possibly work, because the CLI
# does `LOAD '<extension_path>'` and a C++ extension only loads into the exact
# DuckDB version it was built against. Nothing failed loudly, because nothing
# checked. This script is that check.
#
# `duckdb_version` in MainDistributionPipeline.yml is treated as canonical: it is
# what the distribution build actually passes to `make set_duckdb_version`, so it
# is what the shipped binaries are built against.
#
# Usage:
#   ./scripts/check-duckdb-version.sh          # verify; exit 1 on any mismatch
#   ./scripts/check-duckdb-version.sh --quiet  # only print problems
#   ./scripts/check-duckdb-version.sh --print  # print the canonical version; check nothing
#
# Exit codes: 0 = consistent, 1 = mismatch or unreadable input.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

QUIET=false
[[ "${1:-}" == "--quiet" ]] && QUIET=true
PRINT=false
[[ "${1:-}" == "--print" ]] && PRINT=true

WORKFLOW=".github/workflows/MainDistributionPipeline.yml"

problems=()
notes=()

fail()  { problems+=("$1"); }
note()  { notes+=("$1"); }
say()   { [[ "$QUIET" == true ]] || printf '%s\n' "$1"; }

# Compare versions ignoring a leading "v" — some files carry it, some don't
# (Dockerfile and the pip pin are bare, workflows and the cron are v-prefixed).
bare() { printf '%s' "${1#v}"; }

# ---------------------------------------------------------------------------
# Canonical version: the duckdb_version inputs in the distribution workflow.
# Both jobs must already agree with each other before anything else is checked.
# ---------------------------------------------------------------------------

if [[ ! -f "$WORKFLOW" ]]; then
    printf 'ERROR: %s not found (run from the repo, or fix the path)\n' "$WORKFLOW" >&2
    exit 1
fi

mapfile -t wf_versions < <(grep -oE '^[[:space:]]*duckdb_version:[[:space:]]*v?[0-9]+\.[0-9]+\.[0-9]+' "$WORKFLOW" \
                           | grep -oE 'v?[0-9]+\.[0-9]+\.[0-9]+')

# --print: the canonical version alone, for the scripts that build or test against it
# (build_wasm.sh, test_wasm_browser.sh). None of the checks below run.
if $PRINT; then
    if (( ${#wf_versions[@]} == 0 )); then
        printf 'ERROR: no duckdb_version input in %s\n' "$WORKFLOW" >&2
        exit 1
    fi
    printf '%s\n' "${wf_versions[0]}"
    exit 0
fi

# The commit a duckdb release tag names, or empty: the local submodule clone first, then
# ls-remote, which resolves a tag without cloning. ^{} dereferences an annotated tag.
resolve_tag() {
    local tag="$1" sha="" url ls_out
    if [[ -e duckdb/.git ]]; then
        sha=$(git -C duckdb rev-parse -q --verify "refs/tags/${tag}^{commit}" 2>/dev/null || true)
    fi
    if [[ -z "$sha" ]] && [[ "${MIINT_VERSION_CHECK_OFFLINE:-0}" != "1" ]]; then
        url=$(git config -f .gitmodules --get submodule.duckdb.url || printf 'https://github.com/duckdb/duckdb')
        ls_out=$(git ls-remote "$url" "refs/tags/${tag}" "refs/tags/${tag}^{}" 2>/dev/null || true)
        sha=$(awk -v t="refs/tags/${tag}^{}" '$2==t {print $1; exit}' <<<"$ls_out")
        [[ -z "$sha" ]] && sha=$(awk -v t="refs/tags/${tag}" '$2==t {print $1; exit}' <<<"$ls_out")
    fi
    printf '%s' "$sha"
}

# test-with-bowtie2 builds DuckDB as DUCKDB_VERSION (DuckDB v2 takes its version from it, not
# from git tags), and run_tests.sh installs that release's httpfs. It must name the tag the
# duckdb submodule is pinned at, in the pre-release mode below too.
pinned_sha=$(git ls-tree HEAD duckdb 2>/dev/null | awk '{print $3}')
test_job_version=$(grep -oE '^[[:space:]]*DUCKDB_VERSION:[[:space:]]*v[0-9]+\.[0-9]+\.[0-9]+' "$WORKFLOW" \
                   | grep -oE 'v[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)
if [[ -z "$test_job_version" ]]; then
    fail "$WORKFLOW: no DUCKDB_VERSION for test-with-bowtie2's build"
elif [[ "$(resolve_tag "$test_job_version")" != "$pinned_sha" ]]; then
    fail "$WORKFLOW: test-with-bowtie2 builds DuckDB as $test_job_version, but the duckdb submodule is pinned at ${pinned_sha:0:10}, not that tag"
fi

# Pre-release: before the target DuckDB is tagged, both jobs take duckdb_version from
# extension-ci-tools' _submodule_version.yml (the duckdb submodule pin) instead of a
# literal. There is then no release version to hold the other spellings to, so report
# that loudly and stop. Only the all-or-nothing form is accepted: a workflow mixing a
# literal with the pin falls through to the normal checks and fails there.
n_pin=$(grep -cE '^[[:space:]]*duckdb_version:[[:space:]]*\$\{\{[[:space:]]*needs\.duckdb-submodule-version\.outputs\.version[[:space:]]*\}\}' "$WORKFLOW" || true)
if (( ${#wf_versions[@]} == 0 && n_pin == 2 )); then
    pinned_sha=$(git ls-tree HEAD duckdb 2>/dev/null | awk '{print $3}')
    msg="PRE-RELEASE: $WORKFLOW builds the duckdb submodule pin (${pinned_sha:0:10}), not a release tag. The Dockerfile, publish cron, Python CLI pin and wasm artifact name are NOT checked until literal duckdb_version inputs are restored on release day."
    printf '%s\n' "$msg"
    [[ "${GITHUB_ACTIONS:-}" == true ]] && printf '::warning::%s\n' "$msg"
    if (( ${#problems[@]} > 0 )); then
        for p in "${problems[@]}"; do printf 'FAIL: %s\n' "$p" >&2; done
        exit 1
    fi
    exit 0
fi

if (( ${#wf_versions[@]} == 0 )); then
    printf 'ERROR: no duckdb_version found in %s\n' "$WORKFLOW" >&2
    exit 1
fi

CANON="${wf_versions[0]}"
for v in "${wf_versions[@]}"; do
    if [[ "$(bare "$v")" != "$(bare "$CANON")" ]]; then
        fail "$WORKFLOW: duckdb_version inputs disagree with each other ($(printf '%s ' "${wf_versions[@]}"))"
        break
    fi
done

# We expect one per job (duckdb-stable-build, code-quality-check). Fewer means a
# job lost its pin; more means a job was added without being accounted for here.
if (( ${#wf_versions[@]} != 2 )); then
    note "$WORKFLOW: expected 2 duckdb_version inputs, found ${#wf_versions[@]} — if a job was added or removed, update this check"
fi

say "Canonical target (from $WORKFLOW): $CANON"
say ""

# ---------------------------------------------------------------------------
# Every other in-repo spelling.
# ---------------------------------------------------------------------------

# check_one <label> <file> <extract-regex> <v-prefixed: yes|no>
check_one() {
    local label="$1" file="$2" regex="$3" prefixed="$4" found
    if [[ ! -f "$file" ]]; then
        fail "$label: $file not found"
        return
    fi
    found=$(grep -oE "$regex" "$file" | grep -oE 'v?[0-9]+\.[0-9]+\.[0-9]+' | head -1 || true)
    if [[ -z "$found" ]]; then
        fail "$label: no version found in $file (pattern changed?)"
        return
    fi
    if [[ "$(bare "$found")" == "$(bare "$CANON")" ]]; then
        say "  ok    $label = $found"
    else
        fail "$label = $found, expected $( [[ "$prefixed" == yes ]] && printf '%s' "$CANON" || bare "$CANON" ) ($file)"
    fi
}

check_one "Dockerfile ARG" \
    "Dockerfile" '^ARG DUCKDB_VERSION=v?[0-9]+\.[0-9]+\.[0-9]+' no

check_one "cron publish default" \
    "scripts/cron-publish-extension.sh" '^: "\$\{DUCKDB_VERSION:=v?[0-9]+\.[0-9]+\.[0-9]+\}"' yes

check_one "cron example env file" \
    "scripts/cron-publish-extension.sh" '^#[[:space:]]+DUCKDB_VERSION=v?[0-9]+\.[0-9]+\.[0-9]+' yes

check_one "python CLI pin" \
    "python/pyproject.toml" 'duckdb==v?[0-9]+\.[0-9]+\.[0-9]+' no

# The verify-wasm job is currently commented out, but its artifact name still
# has to be right for the day it is re-enabled — a stale name there fails at
# download time, long after the build has burned its minutes.
check_one "wasm artifact name (commented verify-wasm job)" \
    "$WORKFLOW" 'miint-v?[0-9]+\.[0-9]+\.[0-9]+-extension-' yes

# ---------------------------------------------------------------------------
# The duckdb submodule must actually be at the canonical tag.
#
# This is the one check that can legitimately be unavailable: CI checkouts do
# not always fetch submodule tags, and without them the tag cannot be resolved
# to a SHA. When that happens we say so explicitly rather than passing quietly.
# ---------------------------------------------------------------------------

say ""


if [[ -z "$pinned_sha" ]]; then
    fail "could not read the duckdb gitlink from HEAD (not a git checkout?)"
else
    tag_sha=$(resolve_tag "$CANON")
    source_desc="local submodule or ls-remote"
    if [[ -z "$tag_sha" ]]; then
        fail "could not resolve duckdb tag $CANON (no local clone with tags, and ls-remote unavailable or MIINT_VERSION_CHECK_OFFLINE=1); pinned gitlink is ${pinned_sha:0:10}"
    elif [[ "$pinned_sha" == "$tag_sha" ]]; then
        say "  ok    duckdb submodule pinned at $CANON (${pinned_sha:0:10}, via $source_desc)"
    else
        fail "duckdb submodule is pinned at ${pinned_sha:0:10}, but $WORKFLOW targets $CANON (${tag_sha:0:10}, via $source_desc)"
    fi
fi

# ---------------------------------------------------------------------------
# Report.
# ---------------------------------------------------------------------------

if (( ${#notes[@]} > 0 )); then
    printf '\n'
    for n in "${notes[@]}"; do printf 'NOTE: %s\n' "$n"; done
fi

if (( ${#problems[@]} > 0 )); then
    printf '\nFAIL: DuckDB version is inconsistent across the repo\n' >&2
    for p in "${problems[@]}"; do printf '  - %s\n' "$p" >&2; done
    printf '\nFix by hand, or run: ./scripts/bump-duckdb-version.sh %s\n' "$CANON" >&2
    exit 1
fi

say ""
say "PASS: every in-repo DuckDB version reference agrees on $CANON"
