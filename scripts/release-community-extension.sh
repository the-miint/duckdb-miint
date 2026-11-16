#!/usr/bin/env bash
#
# release-community-extension.sh
#
# Point miint's descriptor in duckdb/community-extensions at a release tag, and open,
# update or ready the PR that proposes it, from your fork. One step of the release
# procedure in docs/internals/releasing.md, which also says what this script refuses
# and why.
#
# Usage:
#   ./scripts/release-community-extension.sh [--dry-run] <tag> [<key>]
#
#   <key>      the descriptor key to point at the tag: `ref` (default) or a frozen
#              release series' codename, such as `variegata`
#   --dry-run  run every check up to the descriptor change, print the change (or that
#              there is nothing to push) and the PR's state; push nothing and change no
#              PR. The check of the PR's own CI, which only readying a PR needs, is not
#              run.
#
# Environment:
#   COMMUNITY_EXTENSIONS_DIR  clone of your community-extensions fork, with the fork as
#                             `origin` and duckdb/community-extensions as `upstream`
#                             (default: ../community-extensions). The script adds a
#                             temporary worktree to it and fetches into it; it never
#                             switches or edits your checkout.
#
# Needs git and gh (authenticated). Avoids bash 4 features so that macOS's bash 3.2 can
# run it; tested on bash 5. Exits 0 when it proposed something or found nothing to do,
# non-zero when it refused or failed.

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

# Set a community-extensions PR's description. Through the REST API, not `gh pr edit`:
# gh 2.46's `pr edit` also reads the PR's classic Projects, which GitHub refuses, and fails
# without changing anything.
set_pr_body() {
    gh api -X PATCH "repos/$CE_REPO/pulls/$1" -f body="$2" --silent
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

[[ "$TAG" =~ ^v[0-9]+\.[0-9]+\.[0-9]+(-rc\.[0-9]+)?$ ]] ||
    die "'$TAG' is not a release tag (vX.Y.Z or vX.Y.Z-rc.N). The descriptor names only tagged commits: the build reports its version only when the built commit carries such a tag."
[[ "$KEY" =~ ^[a-z]+$ ]] || die "'$KEY' is not a descriptor key: use ref or a frozen series' codename"
TAG_V="${TAG#v}"
VERSION="${TAG_V%%-rc.*}"
IS_RC=false
[[ "$TAG" == *-rc.* ]] && IS_RC=true

# ---------------------------------------------------------------------------
# The repositories: origin is the repository the descriptor builds from, and the
# community-extensions clone has upstream and your fork of it.
# ---------------------------------------------------------------------------

MIINT_REPO="$(github_repo "$(git config --get remote.origin.url || true)")"
[[ -n "$MIINT_REPO" ]] || die "origin is not a GitHub repository; the tag must be where community-extensions fetches it"

git -C "$CE" rev-parse --git-dir > /dev/null 2>&1 ||
    die "COMMUNITY_EXTENSIONS_DIR ($CE) is not a git clone of community-extensions"
[[ "$(github_repo "$(git -C "$CE" config --get remote.upstream.url || true)")" == "$CE_REPO" ]] ||
    die "$CE's 'upstream' remote must be $CE_REPO"
FORK="$(github_repo "$(git -C "$CE" config --get remote.origin.url || true)")"
[[ -n "$FORK" ]] || die "$CE's 'origin' remote must be your GitHub fork of $CE_REPO"
# GitHub names the fork's parent, and the owner's login as GitHub spells it, which is
# what PRs report: the URL may spell it in another case.
fork_info="$(gh api "repos/$FORK" --jq '"\(.owner.login) \(.parent.full_name // "-")"')" ||
    die "cannot read $FORK from GitHub"
read -r FORK_OWNER fork_parent <<< "$fork_info"
[[ "$fork_parent" == "$CE_REPO" ]] || die "$CE's 'origin' ($FORK) is not a fork of $CE_REPO"

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
git -C "$CE" show "refs/remotes/upstream/main:$DESC" > "$WORK/upstream.yml"
n="$(grep -cxF "  github: $MIINT_REPO" "$WORK/upstream.yml" || true)"
[[ "$n" == 1 ]] ||
    die "the descriptor's repo.github is not $MIINT_REPO, the repository origin names; community-extensions checks the commit out from repo.github"

# ---------------------------------------------------------------------------
# The tag
# ---------------------------------------------------------------------------

remote_tags="$(git ls-remote --tags origin)" || die "cannot read origin's tags"
printf '%s\n' "$remote_tags" | awk -v t="refs/tags/$TAG" '$2 == t {f = 1} END {exit !f}' ||
    die "tag $TAG is not on origin ($MIINT_REPO). community-extensions checks the commit out from GitHub: push the tag first."
git fetch -q origin "refs/tags/$TAG:refs/tags/$TAG" ||
    die "the local tag $TAG differs from origin's (a fetch without + refuses to move it)"
[[ "$(git cat-file -t "refs/tags/$TAG")" == tag ]] ||
    die "tag $TAG is lightweight; release tags are annotated (git tag -a): the message is the release's summary, and git describe ranks any annotated tag on the commit above it."
SHA="$(git rev-parse "refs/tags/$TAG^{commit}")"
SHORT="${SHA:0:7}"

# The version a DuckDB 2.0 build of SHA reports is `git describe --tags --match "v*.*.*"`
# in a checkout of origin: among the release tags on SHA, annotated beats lightweight,
# then the newest wins. Ask git the same question, over origin's tags on SHA only.
match_args=""
for t in $(printf '%s\n' "$remote_tags" | awk -v s="$SHA" '$1 == s {sub("^refs/tags/", "", $2); sub("\\^\\{\\}$", "", $2); print $2}' | sort -u); do
    [[ "$t" == v*.*.* ]] || continue
    git fetch -q origin "refs/tags/$t:refs/tags/$t" || die "the local tag $t differs from origin's"
    match_args="$match_args --match $t"
done
# shellcheck disable=SC2086 # one --match per tag; release tag names hold no spaces
reported="$(git describe --tags --exact-match $match_args "$SHA")"
[[ "$reported" == "$TAG" ]] ||
    die "a build of $SHORT would report $reported, not $TAG: git describe prefers an annotated tag, then the newest. Tag the final after the candidate."
say "tag $TAG: annotated, on $MIINT_REPO, at $SHORT, and what its build reports"

# ---------------------------------------------------------------------------
# The release branch for the key
# ---------------------------------------------------------------------------

if [[ "$KEY" == ref ]]; then
    heads="$(git ls-remote --symref origin HEAD)" || die "cannot read origin's default branch"
    BRANCH="$(printf '%s\n' "$heads" | awk '$1 == "ref:" && $3 == "HEAD" {sub("^refs/heads/", "", $2); print $2}')"
    [[ -n "$BRANCH" ]] || die "cannot read origin's default branch"
    PR_BRANCH="miint-v$VERSION"
    CE_WORKFLOW=build.yml
else
    heads="$(git ls-remote --heads origin)" || die "cannot read origin's branches"
    BRANCH="$(printf '%s\n' "$heads" | awk -v k="$KEY" '{sub("^refs/heads/", "", $2)} $2 ~ ("^v[0-9]+\\.[0-9]+-" k "$") {print $2}')"
    n="$(printf '%s' "$BRANCH" | grep -c . || true)"
    [[ "$n" == 1 ]] || die "expected one v<X.Y>-$KEY branch on origin for key '$KEY', found $n"
    PR_BRANCH="miint-$KEY-v$VERSION"
    CE_WORKFLOW="build_$KEY.yml"
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
[[ -n "$run" ]] || die "no CI run of $WORKFLOW for the push of $TAG. Only a tag push builds wasm_eh."
read -r run_sha run_status run_url run_conclusion <<< "$run"
[[ "$run_sha" == "$SHA" ]] || die "the newest CI run for $TAG tested ${run_sha:0:7}, not $SHORT: $run_url"
[[ "$run_status" == completed && "$run_conclusion" == success ]] ||
    die "the CI run for $TAG has not passed (${run_status} ${run_conclusion:-}): $run_url"
say "CI: passed ($run_url)"

# ---------------------------------------------------------------------------
# The descriptor change
# ---------------------------------------------------------------------------

# edit_descriptor IN OUT: point KEY at SHA (and, for ref, the version at VERSION). Only a
# two-space-indented key that already holds a commit (ref, variegata, ...) is rewritten,
# so a field such as repo.github or extension.name never is.
edit_descriptor() {
    local n
    n="$(grep -cE "^  $KEY: [0-9a-f]{40}\$" "$1" || true)"
    [[ "$n" == 1 ]] || die "$DESC has $n '  $KEY: <commit>' lines, expected 1"
    if [[ "$KEY" == ref ]]; then
        n="$(grep -c '^  version: ' "$1" || true)"
        [[ "$n" == 1 ]] || die "$DESC has $n '  version:' lines, expected 1"
        sed -E -e "s/^  ref: [0-9a-f]{40}\$/  ref: $SHA/" -e "s/^  version: .*/  version: '$VERSION'/" "$1" > "$2"
    else
        sed -E -e "s/^  $KEY: [0-9a-f]{40}\$/  $KEY: $SHA/" "$1" > "$2"
    fi
}

edit_descriptor "$WORK/upstream.yml" "$WORK/wanted.yml"
if cmp -s "$WORK/upstream.yml" "$WORK/wanted.yml"; then
    say "community-extensions main already points $KEY at $SHORT; nothing to propose"
    exit 0
fi

# This release's PR: the open one from this fork's branch. gh matches --head by branch
# name on any fork, so another fork's same-named branch is filtered out here.
pr="$(gh pr list --repo "$CE_REPO" --head "$PR_BRANCH" --state open --json number,isDraft,url,headRepositoryOwner \
    --jq "map(select(.headRepositoryOwner.login == \"$FORK_OWNER\")) | .[0] // empty | \"\\(.number) \\(.isDraft) \\(.url)\"")"
fork_branch="$(git -C "$CE" ls-remote --heads origin "refs/heads/$PR_BRANCH")" || die "cannot read $FORK's branches"
if [[ -n "$pr" ]]; then
    git -C "$CE" fetch -q origin "+refs/heads/$PR_BRANCH:refs/remotes/origin/$PR_BRANCH"
    BASE="refs/remotes/origin/$PR_BRANCH"
elif [[ -n "$fork_branch" ]]; then
    die "$FORK has a branch $PR_BRANCH with no open PR from it, left by an attempt that did not finish. Delete it (git -C $CE push origin --delete $PR_BRANCH) and rerun."
else
    BASE=refs/remotes/upstream/main
fi
git -C "$CE" show "$BASE:$DESC" > "$WORK/current.yml"
edit_descriptor "$WORK/current.yml" "$WORK/wanted.yml"
WILL_PUSH=false
cmp -s "$WORK/current.yml" "$WORK/wanted.yml" || WILL_PUSH=true
pr_number="" pr_draft="" pr_url=""
[[ -n "$pr" ]] && read -r pr_number pr_draft pr_url <<< "$pr"

if [[ "$DRY_RUN" == true ]]; then
    if [[ "$WILL_PUSH" == true ]]; then
        say "dry run: the descriptor change on $FORK_OWNER:$PR_BRANCH would be"
        diff -u -L "a/$DESC" -L "b/$DESC" "$WORK/current.yml" "$WORK/wanted.yml" || true
    else
        say "dry run: $FORK_OWNER:$PR_BRANCH already points $KEY at $SHORT; nothing to push"
    fi
    if [[ -n "$pr" ]]; then
        say "PR: $pr_url ($([[ "$pr_draft" == true ]] && echo draft || echo ready))"
    fi
    exit 0
fi

# ---------------------------------------------------------------------------
# The push, and the PR. Only the final tag, through the check of the PR's own CI below,
# readies a PR; a ready PR goes back to draft before anything is pushed to it.
# ---------------------------------------------------------------------------

if [[ -n "$pr" && "$pr_draft" != true && ("$IS_RC" == true || "$WILL_PUSH" == true) ]]; then
    gh pr ready "$pr_number" --repo "$CE_REPO" --undo
    pr_draft=true
    say "PR back to draft: $pr_url"
fi

body="Points \`repo.$KEY\` at [$TAG](https://github.com/$MIINT_REPO/tree/$TAG) of $MIINT_REPO ($SHA)"
[[ "$KEY" == ref ]] && body+=" and \`extension.version\` at $VERSION"
body+="."
if [[ "$IS_RC" == true ]]; then
    draft_body="$body"$'\n\n'"Draft until v$VERSION is tagged on the same commit."
else
    draft_body="$body"$'\n\n'"Draft until this PR's CI passes."
fi
RERUN_CMD="./scripts/release-community-extension.sh v$VERSION"
[[ "$KEY" != ref ]] && RERUN_CMD+=" $KEY"
RERUN="rerun $RERUN_CMD once the PR's community-extensions build has passed, to mark it ready"

if [[ "$WILL_PUSH" == true ]]; then
    WORKTREE="$WORK/worktree"
    git -C "$CE" worktree add -q --detach "$WORKTREE" "$BASE"
    cp "$WORK/wanted.yml" "$WORKTREE/$DESC"
    git -C "$WORKTREE" commit -q -m "miint: $KEY -> $TAG ($SHORT)" -- "$DESC"
    git -C "$WORKTREE" push -q origin "HEAD:refs/heads/$PR_BRANCH"
    say "pushed $FORK_OWNER:$PR_BRANCH: $KEY -> $SHORT"
    if [[ -n "$pr" ]]; then
        set_pr_body "$pr_number" "$draft_body"
    fi
fi

if [[ -z "$pr" ]]; then
    url="$(gh pr create --repo "$CE_REPO" --base main --head "$FORK_OWNER:$PR_BRANCH" \
        --title "miint: v$VERSION" --body "$draft_body" --draft)"
    if [[ "$IS_RC" == true ]]; then
        say "opened draft PR: $url"
    else
        say "opened draft PR: $url; $RERUN"
    fi
    exit 0
fi

if [[ "$IS_RC" == true ]]; then
    say "PR: $pr_url (draft, for candidate $TAG)"
elif [[ "$WILL_PUSH" == true ]]; then
    say "PR: $pr_url stays a draft: $TAG is on a commit its CI has not built; $RERUN"
elif [[ "$pr_draft" != true ]]; then
    say "PR: $pr_url (already ready)"
else
    # Ready only the build its CI passed: the PR's head must be the fork branch that names
    # SHA, and the key's community-extensions workflow (build.yml for ref, build_<key>.yml
    # for a frozen series) must have built that commit and passed. A run whose key the PR
    # did not change skips its build and still concludes success; choosing the key's own
    # workflow avoids that run, and requiring a passed platform job (".../build /
    # <platform>", not the matrix job) keeps the check true if community-extensions moves
    # a key to another workflow.
    pr_head="$(gh pr view "$pr_number" --repo "$CE_REPO" --json headRefOid --jq .headRefOid)"
    [[ "$pr_head" == "$(git -C "$CE" rev-parse "refs/remotes/origin/$PR_BRANCH")" ]] ||
        die "PR $pr_url's head is not $FORK_OWNER:$PR_BRANCH as fetched; rerun"
    ce_run="$(gh run list --repo "$CE_REPO" --workflow "$CE_WORKFLOW" --commit "$pr_head" --limit 1 \
        --json databaseId,status,conclusion,url --jq '.[0] // empty | "\(.databaseId) \(.status) \(.url) \(.conclusion)"')"
    [[ -n "$ce_run" ]] || die "PR $pr_url has no $CE_WORKFLOW run on ${pr_head:0:7} yet; $RERUN"
    read -r ce_id ce_status ce_url ce_conclusion <<< "$ce_run"
    [[ "$ce_status" == completed && "$ce_conclusion" == success ]] ||
        die "PR $pr_url's CI has not passed (${ce_status} ${ce_conclusion:-}): $ce_url; $RERUN"
    built="$(gh run view "$ce_id" --repo "$CE_REPO" --json jobs \
        --jq '[.jobs[] | select((.name | test("(^| / )build / ")) and (.name | test("/ Generate matrix$") | not) and .conclusion == "success")] | length')"
    [[ "$built" -gt 0 ]] ||
        die "PR $pr_url's $CE_WORKFLOW run has no passed platform build job (a '.../build / <platform>' job): $ce_url"
    set_pr_body "$pr_number" "$body"
    gh pr ready "$pr_number" --repo "$CE_REPO"
    say "PR marked ready: $pr_url"
fi
