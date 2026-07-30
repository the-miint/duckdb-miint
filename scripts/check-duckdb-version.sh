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
#
# Exit codes: 0 = consistent, 1 = mismatch or unreadable input.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

QUIET=false
[[ "${1:-}" == "--quiet" ]] && QUIET=true

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

# The recorded gitlink is readable straight out of the tree, so this works even
# when the submodule was never initialised (the usual CI case).
pinned_sha=$(git ls-tree HEAD duckdb 2>/dev/null | awk '{print $3}')

if [[ -z "$pinned_sha" ]]; then
    fail "could not read the duckdb gitlink from HEAD (not a git checkout?)"
else
    # Prefer the local submodule clone; fall back to ls-remote, which resolves a
    # tag without cloning. ^{} dereferences annotated tags to their commit.
    tag_sha=""
    source_desc=""
    if [[ -e duckdb/.git ]]; then
        tag_sha=$(git -C duckdb rev-parse -q --verify "refs/tags/${CANON}^{commit}" 2>/dev/null || true)
        [[ -n "$tag_sha" ]] && source_desc="local submodule"
    fi
    if [[ -z "$tag_sha" ]] && [[ "${MIINT_VERSION_CHECK_OFFLINE:-0}" != "1" ]]; then
        url=$(git config -f .gitmodules --get submodule.duckdb.url || printf 'https://github.com/duckdb/duckdb')
        # Ask for both the tag ref and its peeled form: DuckDB's release tags are
        # lightweight (no ^{} entry exists), but an annotated tag would need the
        # peeled SHA to compare against a commit. Prefer peeled when present.
        ls_out=$(git ls-remote "$url" "refs/tags/${CANON}" "refs/tags/${CANON}^{}" 2>/dev/null || true)
        tag_sha=$(awk -v t="refs/tags/${CANON}^{}" '$2==t {print $1; exit}' <<<"$ls_out")
        [[ -z "$tag_sha" ]] && tag_sha=$(awk -v t="refs/tags/${CANON}" '$2==t {print $1; exit}' <<<"$ls_out")
        [[ -n "$tag_sha" ]] && source_desc="ls-remote"
    fi

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
