#!/bin/bash
# shellcheck disable=SC2016 # $names inside single-quoted jq programs are jq variables
# Scenario tests for scripts/release-community-extension.sh, which points miint's
# descriptor in duckdb/community-extensions at a release tag and opens, updates or
# readies the PR. Covers what it refuses (a commit hash, a lightweight tag, a tag that
# was never pushed, an unreachable origin, an origin that is not the descriptor's
# repository, a tag off the release branch, a tag whose CI has not passed, a frozen-line
# key with a tag from another line, a key that does not hold a commit, a fork branch
# left from an abandoned PR), what it does to the PR for each tag (a candidate, a second
# candidate, the final tag before and after the PR's CI passes, a final tag on a new
# commit), and that it ignores another fork's same-named branch.
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
    : > "$SD/gh.calls"

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

# tag NAME SHA [lightweight]: create the tag in the miint clone and push it.
tag() {
    if [ "${3:-}" = lightweight ]; then
        git -C "$SD/miint" tag "$1" "$2"
    else
        git -C "$SD/miint" tag -a "$1" -m "Release $1" "$2"
    fi
    git -C "$SD/miint" push -q origin "refs/tags/$1"
}

# add_run TAG SHA STATUS CONCLUSION: the tag push's CI run, as gh run list reports it.
add_run() {
    local id
    id=$(($(jq length "$SD/runs.json") + 1))
    jq -c --argjson id "$id" --arg b "$1" --arg s "$2" --arg st "$3" --arg c "$4" \
        '. + [{databaseId: $id, event: "push", headBranch: $b, headSha: $s, status: $st, conclusion: $c, url: ("https://example/run/" + ($id | tostring))}]' \
        "$SD/runs.json" > "$SD/runs.new"
    mv "$SD/runs.new" "$SD/runs.json"
}

prs() {
    jq -c "$@" "$SD/prs.json" > "$SD/prs.new"
    mv "$SD/prs.new" "$SD/prs.json"
}
# set_checks PR CONCLUSION...: the PR's community-extensions checks; "" is still running.
set_checks() {
    local n="$1" checks
    shift
    checks="$(printf '%s\n' "$@" | jq -R '{conclusion: .}' | jq -sc '.')"
    prs --argjson n "$n" --argjson c "$checks" 'map(if .number == $n then .checks = $c else . end)'
}
set_pr() { prs --argjson n "$1" "map(if .number == \$n then $2 else . end)"; }
# add_foreign_pr BRANCH OWNER: an open draft PR from someone else's fork.
add_foreign_pr() {
    prs --arg h "$1" --arg o "$2" \
        '. + [{number: (length + 1), headRefName: $h, headRepositoryOwner: {login: $o}, isDraft: true, state: "OPEN", url: "https://example/pull/foreign", checks: []}]'
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
# gh calls that change a PR.
pr_changes() { grep -cE '^pr (create|ready)' "$SD/gh.calls" || true; }
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

# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
set_checks 1 SUCCESS SKIPPED
tag v1.1.0 "$T2"
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" = 0 ] && [ "$(fork_branch miint-v1.1.0)" = "$before" ] && [ "$(pr_count)" = 1 ] &&
    [ "$(pr_at 0 isDraft)" = false ] && ce_untouched; then
    ok=1
fi
report "the final tag on the candidate's commit readies the PR without pushing to it" \
    "the final tag names the same build the candidate did; pushing to the PR again would rebuild every platform on community-extensions for nothing" "$ok"

# ---------------------------------------------------------------------------
for checks in "SUCCESS ''" "none"; do
    new_scenario
    propose_rc1
    if [ "$checks" != none ]; then set_checks 1 SUCCESS ""; fi
    tag v1.1.0 "$T2"
    add_run v1.1.0 "$T2" completed success
    invoke v1.1.0
    ok=0
    if [ "$EXIT" != 0 ] && [ "$(pr_at 0 isDraft)" = true ] && said "CI has not passed"; then ok=1; fi
    report "the final tag leaves the PR a draft until its community-extensions CI passes (checks: $checks)" \
        "readying the PR offers the build to DuckDB's maintainers to merge; a build whose own CI is still running or never ran must not be offered" "$ok"
done

# ---------------------------------------------------------------------------
new_scenario
propose_rc1
T3="$(commit T3)"
git -C "$SD/miint" push -q origin v2.0-cyanoptera
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

# ---------------------------------------------------------------------------
new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
T3="$(commit T3)"
git -C "$SD/miint" push -q origin v2.0-cyanoptera
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
after="$(fork_branch miint-v1.1.0)"
ok=0
if [ "$EXIT" = 0 ] && [ -n "$after" ] && [ "$after" != "$before" ] &&
    git --git-dir="$SD/ce-fork.git" merge-base --is-ancestor "$before" "$after" &&
    has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_count)" = 1 ] && [ "$(pr_at 0 isDraft)" = true ]; then
    ok=1
fi
report "a second candidate on a new commit updates the same draft PR" \
    "one release is one PR: the maintainers review a single thread, and its CI always tests the newest candidate" "$ok"

# ---------------------------------------------------------------------------
new_scenario
propose_rc1
set_pr 1 '.isDraft = false'
T3="$(commit T3)"
git -C "$SD/miint" push -q origin v2.0-cyanoptera
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
ok=0
if [ "$EXIT" = 0 ] && has_line miint-v1.1.0 "  ref: $T3" && [ "$(pr_at 0 isDraft)" = true ]; then ok=1; fi
report "a candidate pushed to a ready PR turns it back into a draft" \
    "a ready PR can be merged, and a candidate's commit must not ship as the release that extension.version names" "$ok"

# ---------------------------------------------------------------------------
new_scenario
propose_rc1
before="$(fork_branch miint-v1.1.0)"
set_pr 1 '.state = "CLOSED"'
T3="$(commit T3)"
git -C "$SD/miint" push -q origin v2.0-cyanoptera
tag v1.1.0-rc.2 "$T3"
add_run v1.1.0-rc.2 "$T3" completed success
invoke v1.1.0-rc.2
ok=0
if [ "$EXIT" != 0 ] && [ "$(fork_branch miint-v1.1.0)" = "$before" ] && [ "$(pr_count)" = 1 ] &&
    said "push origin --delete miint-v1.1.0"; then
    ok=1
fi
report "a fork branch left by a closed PR is refused, with the command that clears it" \
    "that branch holds an abandoned descriptor on an old upstream; building on it would propose a PR that conflicts with what upstream has since merged" "$ok"

# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
new_scenario
add_run v1.1.0-rc.1 "$T2" completed success
invoke "$T2"
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'vX.Y.Z'; then ok=1; fi
report "a commit hash is refused" \
    "community-extensions builds whatever ref names, and only a commit carrying a vX.Y.Z tag makes the shipped binary report its version instead of a hash" "$ok"

# ---------------------------------------------------------------------------
new_scenario
tag v1.1.0-rc.1 "$T2"
tag v1.1.0 "$T2" lightweight
add_run v1.1.0 "$T2" completed success
invoke v1.1.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'annotated'; then ok=1; fi
report "a lightweight tag is refused" \
    "git describe prefers an annotated tag on the same commit, so a lightweight final tag beside an annotated candidate ships a binary that reports the candidate's version" "$ok"

# ---------------------------------------------------------------------------
new_scenario
git -C "$SD/miint" tag -a v1.1.0-rc.1 -m "Release v1.1.0-rc.1" "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'push the tag first'; then ok=1; fi
report "a tag that was never pushed is refused" \
    "community-extensions checks the commit out from GitHub; a tag that exists only locally is not there for the build to report" "$ok"

# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
new_scenario
git -C "$SD/miint" remote set-url origin https://github.com/someone/duckdb-miint.git
git -C "$SD/miint" config "url.$SD/miint.git.insteadOf" https://github.com/someone/duckdb-miint.git
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'repo.github'; then ok=1; fi
report "an origin other than the descriptor's repository is refused" \
    "community-extensions checks the commit out from repo.github; a tag and branch found on a personal fork may not exist there" "$ok"

# ---------------------------------------------------------------------------
new_scenario
tag v1.2.0 "$F1"
add_run v1.2.0 "$F1" completed success
invoke v1.2.0
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'v2.0-cyanoptera'; then ok=1; fi
report "a tag off the default branch is refused" \
    "only code merged to the release branch ships; a tag on an unmerged branch would publish code no review has passed" "$ok"

# ---------------------------------------------------------------------------
for state in "completed failure" "in_progress " "none"; do
    new_scenario
    tag v1.1.0-rc.1 "$T2"
    if [ "$state" != none ]; then
        # shellcheck disable=SC2086 # the status and conclusion are two words on purpose
        add_run v1.1.0-rc.1 "$T2" $state
    fi
    invoke v1.1.0-rc.1
    ok=0
    if [ "$EXIT" != 0 ] && nothing_proposed && said 'CI'; then ok=1; fi
    report "a tag whose CI run has not passed is refused ($state)" \
        "wasm_eh and the full platform list build only on a tag push; a tag whose run has not passed may not build on community-extensions either" "$ok"
done

# ---------------------------------------------------------------------------
new_scenario
tag v1.0.2 "$L2"
add_run v1.0.2 "$L2" completed success
invoke v1.0.2 variegata
ok=0
B=miint-variegata-v1.0.2
if [ "$EXIT" = 0 ] && [ -n "$(fork_branch $B)" ] && has_line $B "  variegata: $L2" &&
    has_line $B "  ref: $OLD_REF" && has_line $B "  version: '1.0.1'" && [ "$(changed_lines $B)" = 1 ] &&
    [ "$(pr_count)" = 1 ] && [ "$(pr_at 0 isDraft)" = false ]; then
    ok=1
fi
report "a frozen-line release changes only its own key" \
    "community-extensions rebuilds only the keys a PR changes, and extension.version describes the latest-stable line, so touching ref or version would rebuild or mislabel the current release" "$ok"

# ---------------------------------------------------------------------------
new_scenario
tag v1.1.0-rc.1 "$T2"
add_run v1.1.0-rc.1 "$T2" completed success
invoke v1.1.0-rc.1 variegata
ok=0
if [ "$EXIT" != 0 ] && nothing_proposed && said 'v1.5-variegata'; then ok=1; fi
report "a frozen-line key needs a tag on that line" \
    "the key decides which DuckDB version builds the commit, and a commit from the DuckDB 2.0 line does not build against DuckDB 1.5" "$ok"

# ---------------------------------------------------------------------------
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

# ---------------------------------------------------------------------------
if [ "$FAILED" = 0 ]; then
    ALL_PASSED=1
    echo "All release-community-extension scenarios passed."
else
    echo "Some release-community-extension scenarios failed; kept $TMP_ROOT"
    exit 1
fi
