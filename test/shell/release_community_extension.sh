#!/bin/bash
# shellcheck disable=SC2016 # $names inside single-quoted jq programs are jq variables
# Scenario tests for scripts/release-community-extension.sh, which points miint's
# descriptor in duckdb/community-extensions at a release tag and opens, updates or
# readies the PR from the maintainer's fork. Each refusal in docs/internals/releasing.md's
# table has a scenario that triggers it and checks that nothing was proposed and what the
# script said; each PR transition (open, push to, back to draft, ready, stay a draft) has
# a scenario that checks the fork branch and the PR afterwards.
#
# Every remote is a local bare repository behind a GitHub URL (url.<path>.insteadOf),
# and gh is a stub (release_community_extension_stubs/), so nothing leaves the machine.
# Needs jq; skips otherwise. The scenarios' temporary directory is removed when all
# pass and kept otherwise.
#
# Usage: bash test/shell/release_community_extension.sh [script-under-test]

set -e

# A git hook or wrapper that exports these would point every git call below at the
# developer's repository instead of the throwaway ones.
unset GIT_DIR GIT_WORK_TREE GIT_INDEX_FILE GIT_OBJECT_DIRECTORY GIT_COMMON_DIR

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCRIPT="${1:-$HERE/../../scripts/release-community-extension.sh}"
FAILED=0
ALL_PASSED=0

if ! command -v jq > /dev/null 2>&1; then
    echo "release-community-extension tests need jq. Skipping."
    exit 0
fi

export PATH="$HERE/release_community_extension_stubs:$PATH"
TMP_ROOT="$(mktemp -d "${TMPDIR:-/tmp}/release_community_extension.XXXXXX")"
trap 'if [ "$ALL_PASSED" = 1 ]; then rm -rf "$TMP_ROOT"; fi' EXIT
SCENARIO=0

# The developer's global hooks, signing and identity must not reach the throwaway
# repositories, for the git calls here and for the ones the script makes.
export GIT_CONFIG_GLOBAL="$TMP_ROOT/gitconfig"
export GIT_CONFIG_NOSYSTEM=1
git config -f "$GIT_CONFIG_GLOBAL" user.name t
git config -f "$GIT_CONFIG_GLOBAL" user.email t@example.org
git config -f "$GIT_CONFIG_GLOBAL" core.hooksPath /dev/null
git config -f "$GIT_CONFIG_GLOBAL" advice.detachedHead false

DESC=extensions/miint/description.yml
OLD_REF=0000000000000000000000000000000000000000
OLD_VARIEGATA=1111111111111111111111111111111111111111
MIINT_URL=https://github.com/the-miint/duckdb-miint.git
FORK_URL=https://github.com/tester/community-extensions.git
UPSTREAM_URL=https://github.com/duckdb/community-extensions.git

commit() {
    printf '%s\n' "$1" > "$SD/miint/f"
    git -C "$SD/miint" add f
    git -C "$SD/miint" commit -q -m "$1"
    git -C "$SD/miint" rev-parse HEAD
}

# A fresh pair of repositories per scenario:
#   miint: v1.5-variegata (L1, L2) and v2.0-cyanoptera (L1, T1, T2), the remote's
#     default branch, plus an unmerged branch (F1 off T2); tags pushed by the scenario.
#   community-extensions: an upstream and a fork, both bare, and the maintainer's clone
#     of the fork, on a branch of their own with an uncommitted edit.
new_scenario() {
    SCENARIO=$((SCENARIO + 1))
    SD="$TMP_ROOT/$SCENARIO"
    export SCENARIO_DIR="$SD"
    mkdir -p "$SD"
    printf '[]' > "$SD/runs.json"
    printf '[]' > "$SD/prs.json"
    printf '{"tester/community-extensions": {"owner": {"login": "tester"}, "parent": {"full_name": "duckdb/community-extensions"}}}' > "$SD/fork.json"
    : > "$SD/gh.calls"
    TAG_SECOND=0

    git init -q --bare -b v2.0-cyanoptera "$SD/miint.git"
    git init -q -b v1.5-variegata "$SD/miint"
    git -C "$SD/miint" remote add origin "$MIINT_URL"
    git -C "$SD/miint" config "url.$SD/miint.git.insteadOf" "$MIINT_URL"
    commit L1 > /dev/null
    L2="$(commit L2)"
    git -C "$SD/miint" checkout -q -b v2.0-cyanoptera HEAD~1
    commit T1 > /dev/null
    T2="$(commit T2)"
    git -C "$SD/miint" checkout -q -b feature
    F1="$(commit F1)"
    git -C "$SD/miint" checkout -q v2.0-cyanoptera
    git -C "$SD/miint" push -q origin v1.5-variegata v2.0-cyanoptera feature
    mkdir -p "$SD/miint/scripts"
    cp "$SCRIPT" "$SD/miint/scripts/release-community-extension.sh"

    git init -q --bare -b main "$SD/ce-upstream.git"
    git clone -q "$SD/ce-upstream.git" "$SD/ce-seed" 2> /dev/null
    mkdir -p "$SD/ce-seed/extensions/miint"
    cat > "$SD/ce-seed/$DESC" << EOF
extension:
  name: miint
  description: test descriptor
  version: '1.0.1'
  language: C++
  build: cmake
  license: BSD-3-Clause
  maintainers:
    - tester

repo:
  github: the-miint/duckdb-miint
  ref: $OLD_REF
  variegata: $OLD_VARIEGATA

docs:
  hello_world: |
    SELECT 1;
  extended_description: |
    Indented lines that only look like keys:
      version: '9.9.9'
      ref: $OLD_REF
EOF
    printf 'readme\n' > "$SD/ce-seed/README.md"
    git -C "$SD/ce-seed" add .
    git -C "$SD/ce-seed" commit -q -m seed
    git -C "$SD/ce-seed" push -q origin main
    git clone -q --bare "$SD/ce-upstream.git" "$SD/ce-fork.git"
    git clone -q "$SD/ce-fork.git" "$SD/ce"
    git -C "$SD/ce" remote set-url origin "$FORK_URL"
    git -C "$SD/ce" config "url.$SD/ce-fork.git.insteadOf" "$FORK_URL"
    git -C "$SD/ce" remote add upstream "$UPSTREAM_URL"
    git -C "$SD/ce" config "url.$SD/ce-upstream.git.insteadOf" "$UPSTREAM_URL"
    git -C "$SD/ce" checkout -q -b maintainers-own-work
    printf 'uncommitted\n' >> "$SD/ce/README.md"
}

# tag NAME SHA [lightweight]: create the tag in the miint clone and push it. Annotated
# tags are dated one second apart, in the order created, as git describe needs.
tag() {
    if [ "${3:-}" = lightweight ]; then
        git -C "$SD/miint" tag "$1" "$2"
    else
        TAG_SECOND=$((TAG_SECOND + 1))
        GIT_COMMITTER_DATE="2026-01-01T00:00:$(printf '%02d' "$TAG_SECOND")" \
            git -C "$SD/miint" tag -a "$1" -m "Release $1" "$2"
    fi
    git -C "$SD/miint" push -q origin "refs/tags/$1"
}

# add_run_in REPO WORKFLOW EVENT BRANCH SHA STATUS CONCLUSION [JOBS]: a CI run as gh
# lists it, with the jobs `gh run view` reports (a JSON array, default none).
add_run_in() {
    local id
    id=$(($(jq length "$SD/runs.json") + 1))
    jq -c --argjson id "$id" --arg r "$1" --arg w "$2" --arg e "$3" --arg b "$4" --arg s "$5" --arg st "$6" --arg c "$7" \
        --argjson j "${8:-[]}" \
        '. + [{databaseId: $id, repo: $r, workflow: $w, event: $e, headBranch: $b, headSha: $s, status: $st, conclusion: $c,
               jobs: $j, url: ("https://example/run/" + ($id | tostring))}]' "$SD/runs.json" > "$SD/runs.new"
    mv "$SD/runs.new" "$SD/runs.json"
}
# add_run TAG SHA STATUS CONCLUSION: the tag push's run in miint's CI.
add_run() { add_run_in the-miint/duckdb-miint MainDistributionPipeline.yml push "$1" "$2" "$3" "$4"; }
# add_ce_run BRANCH STATUS CONCLUSION [WORKFLOW] [skipped]: a community-extensions run of
# WORKFLOW (default build.yml) on the fork branch's head. Its build jobs ran, named as
# GitHub names a called workflow's jobs, unless "skipped": then its build job was skipped,
# as for a key the PR did not change, and the run still concludes as given.
add_ce_run() {
    local wf="${4:-build.yml}" prefix="" jobs
    [ "$wf" = build.yml ] || prefix="${wf%.yml} / "
    if [ "${5:-}" = skipped ]; then
        jobs="[{\"name\": \"${prefix}build\", \"conclusion\": \"skipped\"}]"
    else
        jobs="[{\"name\": \"${prefix}build / linux_amd64\", \"conclusion\": \"$3\"}]"
    fi
    add_run_in duckdb/community-extensions "$wf" pull_request "$1" "$(fork_branch "$1")" "$2" "$3" "$jobs"
}

prs() {
    jq -c "$@" "$SD/prs.json" > "$SD/prs.new"
    mv "$SD/prs.new" "$SD/prs.json"
}
set_pr() { prs --argjson n "$1" "map(if .number == \$n then $2 else . end)"; }
# add_foreign_pr BRANCH OWNER: an open draft PR from someone else's fork.
add_foreign_pr() {
    prs --arg h "$1" --arg o "$2" \
        '. + [{number: (length + 1), headRefName: $h, headRepositoryOwner: {login: $o}, isDraft: true, state: "OPEN", url: "https://example/pull/foreign"}]'
}

invoke() {
    EXIT=0
    (cd "$SD/miint" && COMMUNITY_EXTENSIONS_DIR="$SD/ce" bash scripts/release-community-extension.sh "$@") > "$SD/out" 2>&1 || EXIT=$?
    OUT="$(cat "$SD/out")"
}

fork_branch() { git --git-dir="$SD/ce-fork.git" rev-parse -q --verify "refs/heads/$1" 2> /dev/null || true; }
fork_desc() { git --git-dir="$SD/ce-fork.git" show "$1:$DESC"; }
upstream_desc() { git --git-dir="$SD/ce-upstream.git" show "main:$DESC"; }
# The number of descriptor lines a fork branch changes against upstream main.
changed_lines() { diff <(upstream_desc) <(fork_desc "$1") | grep -c '^>' || true; }
pr_count() { jq length "$SD/prs.json"; }
pr_at() { jq -r ".[$1].$2" "$SD/prs.json"; }
# gh calls that change a PR: a description is set with `api -X PATCH`.
pr_changes() { grep -cE '^(pr (create|ready|edit)|api -X PATCH)' "$SD/gh.calls" || true; }
has_line() { fork_desc "$1" | grep -qxF -- "$2"; }
said() { printf '%s\n' "$OUT" | grep -qF -- "$1"; }
# The maintainer's community-extensions checkout is exactly as they left it.
ce_untouched() {
    [ "$(git -C "$SD/ce" branch --show-current)" = maintainers-own-work ] &&
        [ "$(git -C "$SD/ce" status --porcelain)" = " M README.md" ] &&
        [ "$(git -C "$SD/ce" worktree list | wc -l)" -eq 1 ]
}
nothing_proposed() { [ -z "$(git --git-dir="$SD/ce-fork.git" for-each-ref refs/heads/miint-*)" ] && [ "$(pr_changes)" -eq 0 ]; }

report() {
    local name="$1" why="$2" ok="$3"
    echo "Running: $name"
    if [ "$ok" = 1 ]; then
        echo "  PASS"
    else
        echo "  FAIL"
        echo "  Why it matters: $why"
        echo "  Exit: $EXIT"
        echo "  Output: $(printf '%s' "$OUT" | tr '\n' '|')"
        echo "  gh calls: $(tr '\n' '|' < "$SD/gh.calls")"
        echo "  Scenario: $SD"
        FAILED=1
    fi
}

# A release candidate v1.1.0-rc.1 on T2, its CI passed, proposed.
propose_rc1() {
    tag v1.1.0-rc.1 "$T2"
    add_run v1.1.0-rc.1 "$T2" completed success
    invoke v1.1.0-rc.1
}
# A new commit T3 on the default branch.
new_commit() {
    T3="$(commit T3)"
    git -C "$SD/miint" push -q origin v2.0-cyanoptera
}

# --- PR transitions ----------------------------------------------------------

new_scenario
propose_rc1
ok=0
B=miint-v1.1.0
if [ "$EXIT" = 0 ] && [ -n "$(fork_branch $B)" ] && has_line $B "  ref: $T2" &&
    has_line $B "  version: '1.1.0'" && has_line $B "  variegata: $OLD_VARIEGATA" &&
    [ "$(changed_lines $B)" = 2 ] && [ "$(pr_count)" = 1 ] && [ "$(pr_at 0 isDraft)" = true ] &&
    [ "$(pr_at 0 headRefName)" = $B ] && [ "$(pr_at 0 headRepositoryOwner.login)" = tester ] && ce_untouched; then
    ok=1
fi
report "a release candidate opens a draft PR pointing ref at the tagged commit" \
    "community-extensions' CI checks a candidate before anything can ship: the PR must name the tagged commit, change only ref and version, stay a draft so it cannot be merged early, and leave the maintainer's own checkout alone" "$ok"

new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
add_ce_run miint-v1.1.0 completed success
tag v1.1.0 "$T2"
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" = 0 ] && [ "$(fork_branch miint-v1.1.0)" = "$before" ] && [ "$(pr_count)" = 1 ] &&
    [ "$(pr_at 0 isDraft)" = false ] && pr_at 0 body | grep -qF "tree/v1.1.0)" && ce_untouched; then
    ok=1
fi
report "the final tag on the candidate's commit readies the PR, naming the final tag, without pushing to it" \
    "the final tag names the same build the candidate did; pushing again would rebuild every platform on community-extensions for nothing, and the merged PR should name the release, not the candidate" "$ok"

for ce in "in_progress " "completed failure" "none" "completed success skipped" "completed success matrix-only"; do
    new_scenario
    propose_rc1
    # shellcheck disable=SC2086 # the status and conclusion are two words
    case "$ce" in
        none) ;;
        *skipped) add_ce_run miint-v1.1.0 completed success build.yml skipped ;;
        *matrix-only)
            add_run_in duckdb/community-extensions build.yml pull_request miint-v1.1.0 "$(fork_branch miint-v1.1.0)" completed success \
                '[{"name": "build / Generate matrix", "conclusion": "success"}, {"name": "build / linux_amd64", "conclusion": "skipped"}]'
            ;;
        *) add_ce_run miint-v1.1.0 $ce ;;
    esac
    tag v1.1.0 "$T2"
    add_run v1.1.0 "$T2" completed success
    invoke v1.1.0
    ok=0
    hint="rerun ./scripts/release-community-extension.sh v1.1.0 once"
    case "$ce" in
        none) expected="has no build.yml run" ;;
        *skipped | *matrix-only) expected="has no passed platform build job" hint="" ;;
        *) expected="has not passed" ;;
    esac
    if [ "$EXIT" != 0 ] && [ "$(pr_at 0 isDraft)" = true ] && said "$expected" && { [ -z "$hint" ] || said "$hint"; }; then ok=1; fi
    report "the final tag leaves the PR a draft until its community-extensions build passes (build: $ce)" \
        "readying the PR offers the build to DuckDB's maintainers to merge; a build still running, failed, never run, or skipped must not be offered" "$ok"
done

new_scenario
propose_rc1
add_ce_run miint-v1.1.0 completed success
new_commit
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
tag v1.1.0 "$T3"
add_run v1.1.0 "$T3" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && [ "$(pr_at 0 isDraft)" = true ] && said "has no build.yml run on"; then ok=1; fi
report "a passed build of an earlier candidate does not ready the PR for a later one" \
    "the build that passed must be the build of the commit the PR now names" "$ok"

new_scenario
git -C "$SD/ce" remote set-url origin https://github.com/Tester/community-extensions.git
git -C "$SD/ce" config "url.$SD/ce-fork.git.insteadOf" https://github.com/Tester/community-extensions.git
printf '{"Tester/community-extensions": {"owner": {"login": "tester"}, "parent": {"full_name": "duckdb/community-extensions"}}}' > "$SD/fork.json"
propose_rc1
add_ce_run miint-v1.1.0 completed success
tag v1.1.0 "$T2"
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" = 0 ] && [ "$(pr_count)" = 1 ] && [ "$(pr_at 0 isDraft)" = false ]; then ok=1; fi
report "a fork URL that spells the owner in another case still finds the release's PR" \
    "GitHub reports the owner's login as it spells it; matching the URL's spelling would miss the open PR, and then refuse its branch as abandoned" "$ok"

new_scenario
propose_rc1
tag v1.1.0 "$T2"
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && said "no CI run" && [ "$(pr_changes)" = 1 ] && [ "$(pr_at 0 isDraft)" = true ]; then ok=1; fi
report "a final tag with no CI run of its own is refused, though its candidate's run passed on the same commit" \
    "each tag's own push builds what that tag's binary reports; a candidate's run built a binary that reports the candidate" "$ok"

new_scenario
propose_rc1
add_ce_run miint-v1.1.0 completed success
set_pr 1 '.headRefOid = "ffffffffffffffffffffffffffffffffffffffff"'
tag v1.1.0 "$T2"
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && [ "$(pr_at 0 isDraft)" = true ] && said "head is not"; then ok=1; fi
report "the final tag does not ready a PR whose head is not the fork branch it checked" \
    "the CI that passed is the CI of the commit the script checked; a PR head that moved since names something else" "$ok"

new_scenario
propose_rc1
new_commit
tag v1.1.0 "$T3"
add_run v1.1.0 "$T3" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" = 0 ] && has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_count)" = 1 ] &&
    [ "$(pr_at 0 isDraft)" = true ] && said "stays a draft"; then
    ok=1
fi
report "a final tag on a new commit pushes to the PR but leaves it a draft" \
    "a PR is marked ready only for a build community-extensions' CI has already passed; the new commit's build has only just started" "$ok"

new_scenario
propose_rc1
set_pr 1 '.isDraft = false'
new_commit
tag v1.1.0 "$T3"
add_run v1.1.0 "$T3" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" = 0 ] && has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_at 0 isDraft)" = true ]; then ok=1; fi
report "a final tag on a new commit turns a ready PR back into a draft" \
    "a ready PR can be merged, and the new commit's build has not been checked" "$ok"

new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
new_commit
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
after="$(fork_branch miint-v1.1.0)"
ok=0
if [ "$EXIT" = 0 ] && [ -n "$after" ] && [ "$after" != "$before" ] &&
    git --git-dir="$SD/ce-fork.git" merge-base --is-ancestor "$before" "$after" &&
    has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_count)" = 1 ] && [ "$(pr_at 0 isDraft)" = true ] &&
    pr_at 0 body | grep -qF "tree/v1.1.0-rc.2)"; then
    ok=1
fi
report "a second candidate on a new commit updates the same draft PR, and its body" \
    "one release is one PR: the maintainers review a single thread, and its CI always tests the newest candidate" "$ok"

new_scenario
propose_rc1
set_pr 1 '.isDraft = false'
new_commit
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
ok=0
if [ "$EXIT" = 0 ] && has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_at 0 isDraft)" = true ]; then ok=1; fi
report "a candidate pushed to a ready PR turns it back into a draft" \
    "a ready PR can be merged, and a candidate's commit must not ship as the release that extension.version names" "$ok"

new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
set_pr 1 '.isDraft = false | .undoFails = true'
new_commit
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
ok=0
if [ "$EXIT" != 0 ] && [ "$(fork_branch miint-v1.1.0)" = "$before" ]; then ok=1; fi
report "a PR that cannot be turned back into a draft is not pushed to" \
    "pushing first and failing to undo would leave a mergeable PR naming a commit nobody has checked" "$ok"

new_scenario
tag v1.0.2 "$L2"
add_run v1.0.2 "$L2" completed success
invoke v1.0.2 variegata
first_draft="$(pr_at 0 isDraft)"
first_said_key=0
said "release-community-extension.sh v1.0.2 variegata once" && first_said_key=1
add_ce_run miint-variegata-v1.0.2 completed success build_variegata.yml
invoke v1.0.2 variegata
ok=0
B=miint-variegata-v1.0.2
if [ "$EXIT" = 0 ] && [ "$first_draft" = true ] && [ "$first_said_key" = 1 ] && [ "$(pr_at 0 isDraft)" = false ] && has_line $B "  variegata: $L2" &&
    has_line $B "  ref: $OLD_REF" && has_line $B "  version: '1.0.1'" && [ "$(changed_lines $B)" = 1 ] &&
    [ "$(pr_count)" = 1 ]; then
    ok=1
fi
report "a frozen series' release changes only its own key, and opens as a draft until its build passes" \
    "community-extensions rebuilds only the keys a PR changes, and extension.version describes the latest-stable series; and a release with no candidate gets the same CI check before it is offered to merge" "$ok"

new_scenario
tag v1.0.2 "$L2"
add_run v1.0.2 "$L2" completed success
invoke v1.0.2 variegata
add_ce_run miint-variegata-v1.0.2 completed success build.yml skipped
invoke v1.0.2 variegata
ok=0
if [ "$EXIT" != 0 ] && [ "$(pr_at 0 isDraft)" = true ] && said "has no build_variegata.yml run"; then ok=1; fi
report "a frozen series' PR is readied only by its own series' build" \
    "community-extensions builds a codename key in build_<codename>.yml; build.yml's run on that PR builds nothing and still passes" "$ok"

new_scenario
add_foreign_pr miint-v1.1.0 someone-else
propose_rc1
ok=0
if [ "$EXIT" = 0 ] && [ "$(pr_count)" = 2 ] && [ "$(pr_at 1 headRepositoryOwner.login)" = tester ] &&
    [ "$(pr_at 1 isDraft)" = true ] && [ "$(pr_at 0 isDraft)" = true ] && [ -n "$(fork_branch miint-v1.1.0)" ]; then
    ok=1
fi
report "another fork's open PR from a same-named branch is not taken for this release's" \
    "gh matches --head by branch name on any fork; acting on that PR would ready or report someone else's change, and never propose ours" "$ok"

new_scenario
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke --dry-run v1.1.0-rc.1
ok=0
if [ "$EXIT" = 0 ] && nothing_proposed && ce_untouched &&
    printf '%s\n' "$OUT" | grep -qxF "+  ref: $T2" && printf '%s\n' "$OUT" | grep -qxF "+  version: '1.1.0'"; then
    ok=1
fi
report "--dry-run shows the descriptor change and proposes nothing" \
    "the maintainer can see exactly what will be proposed before anything reaches a shared repository" "$ok"

new_scenario
propose_rc1
tag v1.1.0 "$T2"
add_run v1.1.0 "$T2" completed success
invoke --dry-run v1.1.0
ok=0
if [ "$EXIT" = 0 ] && [ "$(pr_changes)" = 1 ] && said "nothing to push" && said "(draft)"; then ok=1; fi
report "--dry-run with nothing to push says so, and the PR's state" \
    "before readying, the maintainer checks the run will not push; an empty diff alone does not say which PR it would act on" "$ok"

# --- Refusals ----------------------------------------------------------------

new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
set_pr 1 '.state = "CLOSED"'
new_commit
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
ok=0
if [ "$EXIT" != 0 ] && [ "$(fork_branch miint-v1.1.0)" = "$before" ] && [ "$(pr_count)" = 1 ] &&
    said "push origin --delete miint-v1.1.0"; then
    ok=1
fi
report "a fork branch with no open PR from it is refused, with the command that clears it" \
    "that branch holds an unfinished attempt on an old upstream; building on it would propose a PR that conflicts with what upstream has since merged" "$ok"

new_scenario
add_run v1.1.0-rc.1 "$T2" completed success
invoke "$T2"
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'vX.Y.Z'; then ok=1; fi
report "a commit hash is refused" \
    "community-extensions builds whatever ref names, and only a commit carrying a vX.Y.Z tag makes the shipped binary report its version instead of a hash" "$ok"

new_scenario
tag v1.1.0 "$T2" lightweight
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'lightweight'; then ok=1; fi
report "a lightweight tag is refused" \
    "a release's tag message is its summary, and an annotated tag records who tagged it and when; a lightweight tag has none of them" "$ok"

new_scenario
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
git -C "$SD/miint" tag -d v1.1.0-rc.1 > /dev/null
git -C "$SD/miint" tag -a v1.1.0-rc.1 -m "moved" "$F1"
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "differs from origin's"; then ok=1; fi
report "a local tag that differs from origin's is refused" \
    "the checks would read a commit other than the one community-extensions fetches" "$ok"

new_scenario
tag v1.1.0 "$T2"
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0 "$T2" completed success
git -C "$SD/miint" tag -d v1.1.0-rc.1 > /dev/null
git -C "$SD/miint" tag -a v1.1.0-rc.1 -m "moved" "$F1"
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "local tag v1.1.0-rc.1 differs"; then ok=1; fi
report "a local copy of another release tag on the commit that differs from origin's is refused" \
    "git describe would then see another set of tags than community-extensions' checkout; here origin's build reports the candidate, which the moved local copy would hide" "$ok"

new_scenario
tag v1.1.0 "$T2"
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "would report v1.1.0-rc.1"; then ok=1; fi
report "a final tag older than a candidate on the same commit is refused" \
    "git describe takes the newest annotated tag, so that build would report the candidate's version" "$ok"

new_scenario
git -C "$SD/miint" tag -a v1.1.0-rc.1 -m "Release v1.1.0-rc.1" "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'push the tag first'; then ok=1; fi
report "a tag that was never pushed is refused" \
    "community-extensions checks the commit out from GitHub; a tag that exists only locally is not there for the build to report" "$ok"

new_scenario
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
git -C "$SD/miint" config --unset "url.$SD/miint.git.insteadOf"
git -C "$SD/miint" config "url.$SD/nowhere.git.insteadOf" "$MIINT_URL"
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "cannot read origin" && ! said 'push the tag first'; then ok=1; fi
report "an unreachable origin is reported as such" \
    "telling the maintainer to push a tag that is already pushed sends them after the wrong problem" "$ok"

new_scenario
git -C "$SD/miint" remote set-url origin https://github.com/someone/duckdb-miint.git
git -C "$SD/miint" config "url.$SD/miint.git.insteadOf" https://github.com/someone/duckdb-miint.git
tag v1.1.0-rc.1 "$T2"
add_run_in someone/duckdb-miint MainDistributionPipeline.yml push v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'repo.github'; then ok=1; fi
report "an origin other than the descriptor's repository is refused" \
    "community-extensions checks the commit out from repo.github; a tag and branch found on a personal fork may not exist there" "$ok"

new_scenario
git -C "$SD/ce" remote set-url upstream https://github.com/someone/community-extensions.git
git -C "$SD/ce" config "url.$SD/ce-upstream.git.insteadOf" https://github.com/someone/community-extensions.git
propose_rc1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "'upstream' remote must be duckdb/community-extensions"; then ok=1; fi
report "an upstream other than duckdb/community-extensions is refused" \
    "the change is proposed against upstream's descriptor; any other repository's is not the one DuckDB builds from" "$ok"

new_scenario
printf '{"tester/community-extensions": {"owner": {"login": "tester"}, "parent": null}}' > "$SD/fork.json"
propose_rc1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "is not a fork of duckdb/community-extensions"; then ok=1; fi
report "an origin that is not a fork of duckdb/community-extensions is refused" \
    "the script pushes to origin before it opens the PR; pushing release branches to an unrelated repository is a change nobody asked for" "$ok"

new_scenario
tag v1.2.0 "$F1"
add_run v1.2.0 "$F1" completed success
invoke v1.2.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'v2.0-cyanoptera'; then ok=1; fi
report "a tag off the default branch is refused" \
    "only code merged to the release branch ships; a tag on an unmerged branch would publish code no review has passed" "$ok"

for state in "completed failure" "in_progress " "none" "other-commit"; do
    new_scenario
    tag v1.1.0-rc.1 "$T2"
    # shellcheck disable=SC2086 # the status and conclusion are two words
    case "$state" in
        none) ;;
        other-commit) add_run v1.1.0-rc.1 "$L2" completed success ;;
        *) add_run v1.1.0-rc.1 "$T2" $state ;;
    esac
    invoke v1.1.0-rc.1
    ok=0
    case "$state" in
        none) expected="no CI run" ;;
        other-commit) expected="tested ${L2:0:7}" ;;
        *) expected="has not passed" ;;
    esac
    if [ "$EXIT" != 0 ] && nothing_proposed && said "$expected"; then ok=1; fi
    report "a tag whose CI run has not passed on its commit is refused ($state)" \
        "wasm_eh builds only on a tag push; a tag whose run has not passed, or tested another commit, may not build on community-extensions either" "$ok"
done

new_scenario
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1 variegata
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'v1.5-variegata'; then ok=1; fi
report "a frozen series' key needs a tag on that series" \
    "the key decides which DuckDB version builds the commit, and a commit from the DuckDB 2.0 series does not build against DuckDB 1.5" "$ok"

new_scenario
git -C "$SD/miint" push -q origin "$T2:refs/heads/v9.9-github"
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1 github
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said "'  github: <commit>'"; then ok=1; fi
report "a key that does not hold a commit is refused" \
    "rewriting repo.github, or an extension: field, with a commit hash would propose a descriptor community-extensions cannot build" "$ok"

# ---------------------------------------------------------------------------
if [ "$FAILED" = 0 ]; then
    ALL_PASSED=1
    echo "All release-community-extension scenarios passed."
else
    echo "Some release-community-extension scenarios failed; kept $TMP_ROOT"
    exit 1
fi
