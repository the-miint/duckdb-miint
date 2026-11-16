#!/usr/bin/env bash
#
# release-community-extension.sh
#
# Point miint's descriptor in duckdb/community-extensions at a release tag, and open,
# update or ready the PR that proposes it. docs/internals/releasing.md is the release
# procedure this script is one step of.
#
# `INSTALL miint FROM community` serves whatever commit the descriptor names, and
# community-extensions rebuilds from it, without asking, on every DuckDB release. So the
# descriptor names only commits that carry a release tag, and this script refuses
# anything else before it proposes a change:
#
#   - The tag must be an annotated vX.Y.Z or vX.Y.Z-rc.N. From DuckDB v2.0, the build
#     embeds the extension version only when the built commit carries such a tag
#     (otherwise the binary reports a commit hash), and git describe prefers an
#     annotated tag over a lightweight one on the same commit.
#   - The tag must be on origin, where community-extensions checks the commit out.
#   - Its commit must be on the release branch for the key: the default branch for
#     `ref`, the v<X.Y>-<key> branch for a frozen line's key (e.g. `variegata`).
#   - The tag push's CI run (MainDistributionPipeline.yml) must have passed. Only tag
#     builds include wasm_eh.
#
# For `ref` it also sets extension.version to the release (X.Y.Z, without any -rc.N).
# A frozen line's key leaves the version alone: the field is shared, and describes the
# latest-stable line.
#
# A release is one PR, on the fork branch miint-vX.Y.Z (miint-<key>-vX.Y.Z for a frozen
# line). A release candidate opens it as a draft; a later candidate pushes to it; the
# final tag on the candidate's commit marks it ready without pushing.
#
# The maintainer's community-extensions checkout is never switched or edited: the
# commit is made in a temporary worktree of it.
#
# Usage:
#   ./scripts/release-community-extension.sh [--dry-run] <tag> [<key>]
#
#   <key>      descriptor key to point at the tag: `ref` (default) or a frozen line's
#              codename, such as `variegata`
#   --dry-run  run every check and print the descriptor change; push nothing
#
# Environment:
#   COMMUNITY_EXTENSIONS_DIR  clone of your community-extensions fork, with the fork as
#                             `origin` and duckdb/community-extensions as `upstream`
#                             (default: ../community-extensions)
#   CE_FORK_OWNER             the fork's GitHub owner (default: read from origin's URL)
#   MIINT_REPO                default: the-miint/duckdb-miint
#
# Needs git and gh (authenticated). Exit codes: 0 = proposed (or nothing to do),
# 1 = refused or failed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

CE="${COMMUNITY_EXTENSIONS_DIR:-$REPO_DIR/../community-extensions}"
MIINT_REPO="${MIINT_REPO:-the-miint/duckdb-miint}"
CE_REPO=duckdb/community-extensions
DESC=extensions/miint/description.yml
WORKFLOW=MainDistributionPipeline.yml

die() {
    printf 'ERROR: %s\n' "$1" >&2
    exit 1
}
say() { printf '%s\n' "$1"; }

usage() {
    sed -n 's/^#   \(\.\/scripts.*\)/\1/p' "$0" >&2
    exit 1
}

DRY_RUN=false
args=()
for a in "$@"; do
    case "$a" in
        --dry-run) DRY_RUN=true ;;
        -*) usage ;;
        *) args+=("$a") ;;
    esac
done
((${#args[@]} >= 1 && ${#args[@]} <= 2)) || usage
TAG="${args[0]}"
KEY="${args[1]:-ref}"

# ---------------------------------------------------------------------------
# The tag
# ---------------------------------------------------------------------------

[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-rc\.[0-9]+)?$ ]] ||
    die "'$TAG' is not a release tag (vX.Y.Z or vX.Y.Z-rc.N). The descriptor names only tagged commits: the build reports its version only when the built commit carries such a tag."

# Fetching the tag without + fails when the local tag differs from origin's, as it must.
git fetch -q origin "refs/tags/$TAG:refs/tags/$TAG" 2> /dev/null ||
    die "tag $TAG is not on origin, or differs from origin's. community-extensions checks the commit out from GitHub: push the tag first."
[[ "$(git cat-file -t "refs/tags/$TAG")" == tag ]] ||
    die "tag $TAG is lightweight; release tags must be annotated (git tag -a). git describe prefers an annotated tag on the same commit, so the binary would report that one."
SHA="$(git rev-parse "refs/tags/$TAG^{commit}")"
SHORT="${SHA:0:7}"
say "tag $TAG: annotated, on origin, at $SHORT"

TAG_V="${TAG#v}"
VERSION="${TAG_V%%-rc.*}"
IS_RC=false
[[ "$TAG" == *-rc.* ]] && IS_RC=true

# ---------------------------------------------------------------------------
# The release branch for the key
# ---------------------------------------------------------------------------

if [[ "$KEY" == ref ]]; then
    BRANCH="$(git ls-remote --symref origin HEAD | awk '$1 == "ref:" && $3 == "HEAD" {sub("^refs/heads/", "", $2); print $2}')"
    [[ -n "$BRANCH" ]] || die "cannot read origin's default branch"
    PR_BRANCH="miint-v$VERSION"
else
    [[ "$KEY" =~ ^[a-z]+$ ]] || die "'$KEY' is not a descriptor key: use ref or a frozen line's codename"
    mapfile -t matches < <(git ls-remote --heads origin | awk '{print $2}' | sed -n "s#^refs/heads/\(v[0-9]*\.[0-9]*-$KEY\)\$#\1#p")
    ((${#matches[@]} == 1)) ||
        die "expected one v<X.Y>-$KEY branch on origin for key '$KEY', found ${#matches[@]}"
    BRANCH="${matches[0]}"
    PR_BRANCH="miint-$KEY-v$VERSION"
fi
git fetch -q origin "+refs/heads/$BRANCH:refs/remotes/origin/$BRANCH"
git merge-base --is-ancestor "$SHA" "refs/remotes/origin/$BRANCH" ||
    die "$TAG ($SHORT) is not on origin/$BRANCH, the release branch for key '$KEY'. Only merged code ships."
say "commit $SHORT: on origin/$BRANCH"

# ---------------------------------------------------------------------------
# The tag push's CI run
# ---------------------------------------------------------------------------

# The conclusion goes last: it is empty until the run completes, and read would otherwise
# shift the fields after it.
run="$(gh run list --repo "$MIINT_REPO" --workflow "$WORKFLOW" --branch "$TAG" --event push --limit 1 \
    --json headSha,status,conclusion,url --jq '.[0] // empty | "\(.headSha) \(.status) \(.url) \(.conclusion)"')"
[[ -n "$run" ]] || die "no CI run of $WORKFLOW for the push of $TAG. Only a tag push builds every platform (wasm_eh included)."
read -r run_sha run_status run_url run_conclusion <<< "$run"
[[ "$run_sha" == "$SHA" ]] || die "the newest CI run for $TAG tested ${run_sha:0:7}, not $SHORT: $run_url"
[[ "$run_status" == completed && "$run_conclusion" == success ]] ||
    die "the CI run for $TAG has not passed (${run_status} ${run_conclusion:-}): $run_url"
say "CI: passed ($run_url)"

# ---------------------------------------------------------------------------
# The descriptor change
# ---------------------------------------------------------------------------

git -C "$CE" rev-parse --git-dir > /dev/null 2>&1 ||
    die "COMMUNITY_EXTENSIONS_DIR ($CE) is not a git clone of community-extensions"
git -C "$CE" remote get-url upstream > /dev/null 2>&1 ||
    die "$CE has no 'upstream' remote; add duckdb/community-extensions as upstream"
FORK_OWNER="${CE_FORK_OWNER:-$(git -C "$CE" remote get-url origin | sed -nE 's#^.*github\.com[:/]([^/]+)/[^/]+$#\1#p')}"
[[ -n "$FORK_OWNER" ]] || die "cannot read the fork's owner from $CE's origin URL; set CE_FORK_OWNER"

git -C "$CE" fetch -q upstream "+refs/heads/main:refs/remotes/upstream/main"
BASE=refs/remotes/upstream/main
PR_BRANCH_EXISTS=false
if git -C "$CE" fetch -q origin "+refs/heads/$PR_BRANCH:refs/remotes/origin/$PR_BRANCH" 2> /dev/null; then
    PR_BRANCH_EXISTS=true
fi

WORK="$(mktemp -d "${TMPDIR:-/tmp}/release-community-extension.XXXXXX")"
WORKTREE=""
cleanup() {
    [[ -n "$WORKTREE" ]] && git -C "$CE" worktree remove --force "$WORKTREE" 2> /dev/null
    rm -rf "$WORK"
}
trap cleanup EXIT

# edit_descriptor IN OUT: point KEY at SHA (and, for ref, the version at VERSION). Each
# key must appear exactly once at the two-space indent of extension: and repo:; deeper
# lines (block literals in docs:) are never touched.
edit_descriptor() {
    local key_count version_count
    key_count="$(grep -c "^  $KEY: " "$1" || true)"
    [[ "$key_count" == 1 ]] || die "$DESC has $key_count '  $KEY:' lines, expected 1"
    if [[ "$KEY" == ref ]]; then
        version_count="$(grep -c '^  version: ' "$1" || true)"
        [[ "$version_count" == 1 ]] || die "$DESC has $version_count '  version:' lines, expected 1"
        sed -e "s/^  ref: .*/  ref: $SHA/" -e "s/^  version: .*/  version: '$VERSION'/" "$1" > "$2"
    else
        sed -e "s/^  $KEY: .*/  $KEY: $SHA/" "$1" > "$2"
    fi
}

git -C "$CE" show "$BASE:$DESC" > "$WORK/upstream.yml"
edit_descriptor "$WORK/upstream.yml" "$WORK/wanted.yml"
if cmp -s "$WORK/upstream.yml" "$WORK/wanted.yml"; then
    say "community-extensions main already points $KEY at $SHORT; nothing to propose"
    exit 0
fi

# A release already in flight continues on its own branch.
[[ "$PR_BRANCH_EXISTS" == true ]] && BASE="refs/remotes/origin/$PR_BRANCH"
git -C "$CE" show "$BASE:$DESC" > "$WORK/current.yml"
edit_descriptor "$WORK/current.yml" "$WORK/wanted.yml"

if [[ "$DRY_RUN" == true ]]; then
    say "dry run: the descriptor change on $FORK_OWNER:$PR_BRANCH would be"
    diff -u --label "a/$DESC" --label "b/$DESC" "$WORK/current.yml" "$WORK/wanted.yml" || true
    exit 0
fi

PUSHED=false
if ! cmp -s "$WORK/current.yml" "$WORK/wanted.yml"; then
    WORKTREE="$WORK/worktree"
    git -C "$CE" worktree add -q --detach "$WORKTREE" "$BASE"
    cp "$WORK/wanted.yml" "$WORKTREE/$DESC"
    git -C "$WORKTREE" commit -q -m "miint: $KEY -> $TAG ($SHORT)" -- "$DESC"
    git -C "$WORKTREE" push -q origin "HEAD:refs/heads/$PR_BRANCH"
    PUSHED=true
    say "pushed $FORK_OWNER:$PR_BRANCH: $KEY -> $SHORT"
fi

# ---------------------------------------------------------------------------
# The PR
# ---------------------------------------------------------------------------

pr="$(gh pr list --repo "$CE_REPO" --head "$PR_BRANCH" --state open --json number,isDraft,url \
    --jq '.[0] // empty | "\(.number) \(.isDraft) \(.url)"')"
if [[ -z "$pr" ]]; then
    body="Points \`repo.$KEY\` at [$TAG](https://github.com/$MIINT_REPO/tree/$TAG) of $MIINT_REPO ($SHA)"
    [[ "$KEY" == ref ]] && body+=" and \`extension.version\` at $VERSION"
    body+="."
    draft=()
    if [[ "$IS_RC" == true ]]; then
        draft=(--draft)
        body+=$'\n\n'"Draft until v$VERSION is tagged on the same commit."
    fi
    url="$(gh pr create --repo "$CE_REPO" --base main --head "$FORK_OWNER:$PR_BRANCH" \
        --title "miint: v$VERSION" --body "$body" "${draft[@]}")"
    say "opened $([[ "$IS_RC" == true ]] && echo 'draft ')PR: $url"
    exit 0
fi

read -r pr_number pr_draft pr_url <<< "$pr"
if [[ "$IS_RC" == true || "$pr_draft" != true ]]; then
    say "PR: $pr_url"
elif [[ "$PUSHED" == true ]]; then
    say "PR: $pr_url stays a draft: $TAG is on a commit no candidate's community-extensions CI has built. Once that CI passes: gh pr ready $pr_number --repo $CE_REPO"
else
    gh pr ready "$pr_number" --repo "$CE_REPO"
    say "PR marked ready: $pr_url"
fi
