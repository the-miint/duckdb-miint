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
#   - The tag must be on origin, which must be the repository the descriptor's
#     repo.github names: community-extensions checks the commit out from there.
#   - Its commit must be on the release branch for the key: the default branch for
#     `ref`, the v<X.Y>-<key> branch for a frozen line's key (e.g. `variegata`).
#   - The tag push's CI run (MainDistributionPipeline.yml) must have passed. Only tag
#     builds include wasm_eh.
#
# For `ref` it also sets extension.version to the release (X.Y.Z, without any -rc.N).
# A frozen line's key leaves the version alone: the field is shared, and describes the
# latest-stable line.
#
# A release is one PR, from the fork branch miint-vX.Y.Z (miint-<key>-vX.Y.Z for a
# frozen line). A release candidate opens it as a draft, or pushes to it and turns it
# back into a draft. The final tag marks it ready once the PR's own CI has passed on the
# commit it names; a final tag on a new commit pushes and leaves it a draft.
#
# The maintainer's community-extensions checkout is never switched or edited: the
# commit is made in a temporary worktree of it. The script only fetches into it.
#
# Usage:
#   ./scripts/release-community-extension.sh [--dry-run] <tag> [<key>]
#
#   <key>      descriptor key to point at the tag: `ref` (default) or a frozen line's
#              codename, such as `variegata`
#   --dry-run  run every check and print the descriptor change; push nothing and
#              change no PR
#
# Environment:
#   COMMUNITY_EXTENSIONS_DIR  clone of your community-extensions fork, with the fork as
#                             `origin` and duckdb/community-extensions as `upstream`
#                             (default: ../community-extensions)
#
# Needs git and gh (authenticated); runs on bash 3.2. Exit codes: 0 = proposed (or
# nothing to do), 1 = refused or failed.

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "$SCRIPT_DIR/.." && pwd)"
cd "$REPO_DIR"

CE="${COMMUNITY_EXTENSIONS_DIR:-$REPO_DIR/../community-extensions}"
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

# owner/name from a GitHub remote URL (https, ssh or scp-like), or nothing. The URL is
# read from the config, not `git remote get-url`, which applies insteadOf rewrites.
github_repo() {
    printf '%s\n' "$1" | sed -E -n -e 's#\.git$##' -e 's#^(https://|ssh://git@|git@)github\.com[:/]([^/]+/[^/]+)$#\2#p'
}

DRY_RUN=false
TAG=""
KEY=""
for a in "$@"; do
    case "$a" in
        --dry-run) DRY_RUN=true ;;
        -*) usage ;;
        *)
            if [[ -z "$TAG" ]]; then
                TAG="$a"
            elif [[ -z "$KEY" ]]; then
                KEY="$a"
            else
                usage
            fi
            ;;
    esac
done
[[ -n "$TAG" ]] || usage
KEY="${KEY:-ref}"

# ---------------------------------------------------------------------------
# The tag
# ---------------------------------------------------------------------------

[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-rc\.[0-9]+)?$ ]] ||
    die "'$TAG' is not a release tag (vX.Y.Z or vX.Y.Z-rc.N). The descriptor names only tagged commits: the build reports its version only when the built commit carries such a tag."
[[ "$KEY" =~ ^[a-z]+$ ]] || die "'$KEY' is not a descriptor key: use ref or a frozen line's codename"

MIINT_REPO="$(github_repo "$(git config --get remote.origin.url || true)")"
[[ -n "$MIINT_REPO" ]] || die "origin is not a GitHub repository; the tag must be where community-extensions fetches it"

remote_tag="$(git ls-remote --tags origin "refs/tags/$TAG")" || die "cannot read origin's tags"
[[ -n "$remote_tag" ]] ||
    die "tag $TAG is not on origin ($MIINT_REPO). community-extensions checks the commit out from GitHub: push the tag first."
# Without +, the fetch fails when the local tag differs from origin's, as it must.
git fetch -q origin "refs/tags/$TAG:refs/tags/$TAG" || die "the local tag $TAG differs from origin's"
[[ "$(git cat-file -t "refs/tags/$TAG")" == tag ]] ||
    die "tag $TAG is lightweight; release tags must be annotated (git tag -a). git describe prefers an annotated tag on the same commit, so the binary would report that one."
SHA="$(git rev-parse "refs/tags/$TAG^{commit}")"
SHORT="${SHA:0:7}"
say "tag $TAG: annotated, on $MIINT_REPO, at $SHORT"

TAG_V="${TAG#v}"
VERSION="${TAG_V%%-rc.*}"
IS_RC=false
[[ "$TAG" == *-rc.* ]] && IS_RC=true

# ---------------------------------------------------------------------------
# The release branch for the key
# ---------------------------------------------------------------------------

if [[ "$KEY" == ref ]]; then
    heads="$(git ls-remote --symref origin HEAD)" || die "cannot read origin's default branch"
    BRANCH="$(printf '%s\n' "$heads" | awk '$1 == "ref:" && $3 == "HEAD" {sub("^refs/heads/", "", $2); print $2}')"
    [[ -n "$BRANCH" ]] || die "cannot read origin's default branch"
    PR_BRANCH="miint-v$VERSION"
else
    heads="$(git ls-remote --heads origin)" || die "cannot read origin's branches"
    BRANCH="$(printf '%s\n' "$heads" | awk -v k="$KEY" '{sub("^refs/heads/", "", $2)} $2 ~ ("^v[0-9]+\\.[0-9]+-" k "$") {print $2}')"
    n="$(printf '%s' "$BRANCH" | grep -c . || true)"
    [[ "$n" == 1 ]] || die "expected one v<X.Y>-$KEY branch on origin for key '$KEY', found $n"
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
# The community-extensions clone and the descriptor change
# ---------------------------------------------------------------------------

git -C "$CE" rev-parse --git-dir > /dev/null 2>&1 ||
    die "COMMUNITY_EXTENSIONS_DIR ($CE) is not a git clone of community-extensions"
[[ "$(github_repo "$(git -C "$CE" config --get remote.upstream.url || true)")" == "$CE_REPO" ]] ||
    die "$CE's 'upstream' remote must be $CE_REPO"
FORK="$(github_repo "$(git -C "$CE" config --get remote.origin.url || true)")"
FORK_OWNER="${FORK%%/*}"
[[ -n "$FORK" && "$FORK" != "$CE_REPO" && "$FORK_OWNER" =~ ^[A-Za-z0-9-]+$ ]] ||
    die "$CE's 'origin' remote must be your GitHub fork of $CE_REPO"

git -C "$CE" fetch -q upstream "+refs/heads/main:refs/remotes/upstream/main"

WORK="$(mktemp -d "${TMPDIR:-/tmp}/release-community-extension.XXXXXX")"
WORKTREE=""
cleanup() {
    if [[ -n "$WORKTREE" ]]; then
        git -C "$CE" worktree remove --force "$WORKTREE" 2> /dev/null || true
    fi
    rm -rf "$WORK"
}
trap cleanup EXIT

# edit_descriptor IN OUT: point KEY at SHA (and, for ref, the version at VERSION). Only a
# two-space-indented key that already holds a commit (ref, variegata, ...) is rewritten,
# so a field such as repo.github or extension.name can never be.
edit_descriptor() {
    local n
    n="$(grep -cE "^  $KEY: [0-9a-f]{40}\$" "$1" || true)"
    [[ "$n" == 1 ]] || die "$DESC has $n '  $KEY: <commit>' lines, expected 1"
    if [[ "$KEY" == ref ]]; then
        n="$(grep -c '^  version: ' "$1" || true)"
        [[ "$n" == 1 ]] || die "$DESC has $n '  version:' lines, expected 1"
        sed -e "s/^  ref: .*/  ref: $SHA/" -e "s/^  version: .*/  version: '$VERSION'/" "$1" > "$2"
    else
        sed -e "s/^  $KEY: .*/  $KEY: $SHA/" "$1" > "$2"
    fi
}

git -C "$CE" show "refs/remotes/upstream/main:$DESC" > "$WORK/upstream.yml"
n="$(grep -cxF "  github: $MIINT_REPO" "$WORK/upstream.yml" || true)"
[[ "$n" == 1 ]] || die "the descriptor's repo.github is not $MIINT_REPO, the repository origin names"
edit_descriptor "$WORK/upstream.yml" "$WORK/wanted.yml"
if cmp -s "$WORK/upstream.yml" "$WORK/wanted.yml"; then
    say "community-extensions main already points $KEY at $SHORT; nothing to propose"
    exit 0
fi

# This release's PR: the open one from this fork's branch. gh cannot filter --head by
# owner, so a same-named branch on another fork is filtered out here.
pr="$(gh pr list --repo "$CE_REPO" --head "$PR_BRANCH" --state open --json number,isDraft,url,headRepositoryOwner \
    --jq "map(select(.headRepositoryOwner.login == \"$FORK_OWNER\")) | .[0] // empty | \"\\(.number) \\(.isDraft) \\(.url)\"")"
fork_branch="$(git -C "$CE" ls-remote --heads origin "refs/heads/$PR_BRANCH")" || die "cannot read $FORK's branches"
if [[ -n "$pr" ]]; then
    git -C "$CE" fetch -q origin "+refs/heads/$PR_BRANCH:refs/remotes/origin/$PR_BRANCH"
    BASE="refs/remotes/origin/$PR_BRANCH"
elif [[ -n "$fork_branch" ]]; then
    die "$FORK has a branch $PR_BRANCH but no open PR from it; it holds an abandoned attempt. Delete it (git -C $CE push origin --delete $PR_BRANCH) and rerun."
else
    BASE=refs/remotes/upstream/main
fi
git -C "$CE" show "$BASE:$DESC" > "$WORK/current.yml"
edit_descriptor "$WORK/current.yml" "$WORK/wanted.yml"

if [[ "$DRY_RUN" == true ]]; then
    say "dry run: the descriptor change on $FORK_OWNER:$PR_BRANCH would be"
    diff -u -L "a/$DESC" -L "b/$DESC" "$WORK/current.yml" "$WORK/wanted.yml" || true
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

if [[ -z "$pr" ]]; then
    body="Points \`repo.$KEY\` at [$TAG](https://github.com/$MIINT_REPO/tree/$TAG) of $MIINT_REPO ($SHA)"
    [[ "$KEY" == ref ]] && body+=" and \`extension.version\` at $VERSION"
    body+="."
    draft=""
    if [[ "$IS_RC" == true ]]; then
        draft=--draft
        body+=$'\n\n'"Draft until v$VERSION is tagged on the same commit."
    fi
    url="$(gh pr create --repo "$CE_REPO" --base main --head "$FORK_OWNER:$PR_BRANCH" \
        --title "miint: v$VERSION" --body "$body" ${draft:+"$draft"})"
    say "opened ${draft:+draft }PR: $url"
    exit 0
fi

read -r pr_number pr_draft pr_url <<< "$pr"
if [[ "$IS_RC" == true ]]; then
    if [[ "$pr_draft" != true ]]; then
        gh pr ready "$pr_number" --repo "$CE_REPO" --undo
        say "PR back to draft, for candidate $TAG: $pr_url"
    else
        say "PR: $pr_url"
    fi
elif [[ "$pr_draft" != true ]]; then
    say "PR: $pr_url (already ready)"
elif [[ "$PUSHED" == true ]]; then
    say "PR: $pr_url stays a draft: $TAG is on a commit its community-extensions CI has not built. Rerun this once that CI passes."
else
    # Ready only the build its CI passed: the PR's head must be the branch that names SHA,
    # and every check on it must have finished without failing.
    # shellcheck disable=SC2016 # $c is a jq variable
    checks="$(gh pr view "$pr_number" --repo "$CE_REPO" --json headRefOid,statusCheckRollup \
        --jq '"\(.headRefOid) \(.statusCheckRollup | length) \([.statusCheckRollup[] | select((.conclusion // .state) as $c | $c == "SUCCESS" or $c == "SKIPPED" or $c == "NEUTRAL")] | length)"')"
    read -r pr_head n_checks n_passed <<< "$checks"
    [[ "$pr_head" == "$(git -C "$CE" rev-parse "refs/remotes/origin/$PR_BRANCH")" ]] ||
        die "PR $pr_url's head is not $FORK_OWNER:$PR_BRANCH as fetched; rerun"
    [[ "$n_checks" -gt 0 && "$n_checks" == "$n_passed" ]] ||
        die "PR $pr_url's community-extensions CI has not passed ($n_passed of $n_checks checks); rerun once it has"
    gh pr ready "$pr_number" --repo "$CE_REPO"
    say "PR marked ready: $pr_url"
fi
