#!/usr/bin/env bash
#
# cron-publish-extension.sh
#
# Poll GitHub Actions for the latest passing build of duckdb-miint on a
# given branch, download the extension artifacts, gzip them, and publish
# them to a custom DuckDB extension repository via rsync. Designed to run
# from cron on a host other than the builder.
#
# Deployment is atomic via a symlink swap: each release lives under
# $DEST_BASE/releases/<duckdb_version>-run<run_id>/ and
# $DEST_BASE/<duckdb_version> is a symlink that flips to the new release
# in a single rename(2) call. Old releases are left in place; prune them
# manually if disk space becomes an issue.
#
# Idempotent: re-running without a new passing build is a no-op deploy,
# but always emails a status. Concurrent runs are prevented via flock.
# A FAIL email is sent if anything in the pipeline errors out.
#
# Prereqs on the cron host (none scripted — all assumed ready):
#   - `gh` authenticated (check with: gh auth status)
#   - SSH public key trusted by $DEST_USER@$DEST_HOST (no-password rsync/ssh)
#   - `mail` configured to deliver to $NOTIFY_EMAIL
#   - gzip, rsync, ssh, flock, jq, find, mktemp (standard)
#
# First-deploy note: if $DEST_BASE/<duckdb_version> already exists as a
# plain directory, the symlink swap will fail (rename(2) refuses to
# replace a non-empty directory with a non-directory). Rename or move
# that directory out of the way before the first cron run.
#
# Configure via env vars or by passing an env file as $1. All vars marked
# REQUIRED must be set; others have sensible defaults. See bottom of file
# for an example env file.
#
# Example crontab (daily at 03:17 local time):
#   17 3 * * * /path/to/repo/scripts/cron-publish-extension.sh /etc/miint-publish.env

set -euo pipefail

# ----- config -----

if [[ -n "${1:-}" ]]; then
    # shellcheck disable=SC1090
    source "$1"
fi

: "${REPO:=the-miint/duckdb-miint}"
: "${BRANCH:=v1.5-variegata}"
: "${WORKFLOW:=MainDistributionPipeline.yml}"

# DUCKDB_VERSION must match both the duckdb submodule pin and the
# duckdb_version input in .github/workflows/MainDistributionPipeline.yml.
# When the extension targets a new DuckDB version, bump all three together.
: "${DUCKDB_VERSION:=v1.5.2}"

: "${DEST_USER:?DEST_USER is required}"
: "${DEST_HOST:?DEST_HOST is required}"
: "${DEST_BASE:?DEST_BASE is required (absolute path on \$DEST_HOST)}"
: "${NOTIFY_EMAIL:?NOTIFY_EMAIL is required}"

: "${STATE_DIR:=${HOME}/.local/state/miint-publish}"
STATE_FILE="${STATE_DIR}/last_run_id"
LOCK_FILE="${STATE_DIR}/lock"
WORK_DIR="${STATE_DIR}/work"
LOG_FILE="${STATE_DIR}/last_run.log"

ARTIFACT_PREFIX="miint-${DUCKDB_VERSION}-extension-"
HOSTNAME_SHORT="$(hostname -s)"

mkdir -p "$STATE_DIR" "$WORK_DIR" "$WORK_DIR/.empty"

# Log everything to $LOG_FILE (truncated each run). Cron swallows stderr
# unless captured; this keeps a record for the FAIL email tail.
: > "$LOG_FILE"
exec >>"$LOG_FILE" 2>&1

log() { printf '[%s] %s\n' "$(date -u +%FT%TZ)" "$*"; }

send_mail() {
    local subject="$1" body="$2"
    printf '%s\n' "$body" | mail -s "$subject" "$NOTIFY_EMAIL"
}

fail() {
    local msg="$1"
    log "FAIL: $msg"
    local subject="[FAIL] miint publish on ${HOSTNAME_SHORT} — ${msg}"
    local body
    body=$(cat <<EOF
Host:          $HOSTNAME_SHORT
Repo:          $REPO
Branch:        $BRANCH
Workflow:      $WORKFLOW
DuckDB ver:    $DUCKDB_VERSION
Destination:   ${DEST_USER}@${DEST_HOST}:${DEST_BASE}

--- last log lines ---
$(tail -n 200 "$LOG_FILE")
EOF
)
    send_mail "$subject" "$body"
    exit 1
}

trap 'fail "unexpected error at line $LINENO"' ERR

# ----- lock -----

exec 9>"$LOCK_FILE"
if ! flock -n 9; then
    log "another instance is holding the lock; exiting quietly"
    exit 0
fi

# ----- reset fixed work dirs (rsync --delete trick; no rm) -----

reset_dir() {
    local d="$1"
    mkdir -p "$d"
    rsync -a --delete "$WORK_DIR/.empty/" "$d/"
}

ART_DIR="$WORK_DIR/artifacts"
STAGE_DIR="$WORK_DIR/stage"
reset_dir "$ART_DIR"
reset_dir "$STAGE_DIR"

# ----- query latest passing run -----

log "Querying latest passing run of $WORKFLOW on $REPO@$BRANCH"
RUN_JSON=$(gh run list \
    --repo "$REPO" \
    --workflow "$WORKFLOW" \
    --branch "$BRANCH" \
    --status success \
    --limit 1 \
    --json databaseId,headSha,displayTitle,createdAt,url)

RUN_ID=$(jq -r '.[0].databaseId // empty' <<<"$RUN_JSON")
RUN_SHA=$(jq -r '.[0].headSha // empty' <<<"$RUN_JSON")
RUN_TITLE=$(jq -r '.[0].displayTitle // empty' <<<"$RUN_JSON")
RUN_TIME=$(jq -r '.[0].createdAt // empty' <<<"$RUN_JSON")
RUN_URL=$(jq -r '.[0].url // empty' <<<"$RUN_JSON")

[[ -n "$RUN_ID" ]] || fail "no passing run found for $WORKFLOW on $BRANCH"

LAST_RUN_ID=""
[[ -f "$STATE_FILE" ]] && LAST_RUN_ID=$(cat "$STATE_FILE")

log "Latest passing run: $RUN_ID ($RUN_SHA)"
log "Last deployed run:  ${LAST_RUN_ID:-<none>}"

# ----- no-op path: heartbeat SUCCEED -----

if [[ "$RUN_ID" == "$LAST_RUN_ID" ]]; then
    subject="[SUCCEED] miint publish on ${HOSTNAME_SHORT} — no new build (run ${RUN_ID})"
    body=$(cat <<EOF
Host:          $HOSTNAME_SHORT
Repo:          $REPO @ $BRANCH
Workflow:      $WORKFLOW
DuckDB ver:    $DUCKDB_VERSION

No new passing build since last deploy.

Latest passing run (already deployed):
  id:    $RUN_ID
  sha:   $RUN_SHA
  title: $RUN_TITLE
  time:  $RUN_TIME
  url:   $RUN_URL
EOF
)
    send_mail "$subject" "$body"
    log "No new build; exit."
    exit 0
fi

# ----- download artifacts -----

log "Downloading artifacts for run $RUN_ID into $ART_DIR"
gh run download "$RUN_ID" \
    --repo "$REPO" \
    --dir "$ART_DIR" \
    --pattern "${ARTIFACT_PREFIX}*"

# ----- stage: gzip + lay out release tree -----

RELEASE_NAME="${DUCKDB_VERSION}-run${RUN_ID}"
RELEASE_DIR="$STAGE_DIR/$RELEASE_NAME"
mkdir -p "$RELEASE_DIR"

deployed=()
while IFS= read -r -d '' art; do
    name=$(basename "$art")
    # Strip the prefix; remainder is the platform directory name that the
    # DuckDB client's Platform() resolves to at INSTALL time.
    platform="${name#"$ARTIFACT_PREFIX"}"
    if [[ "$platform" == "$name" ]]; then
        log "  skip $name (no prefix match)"
        continue
    fi

    # Accept .duckdb_extension (CI default) or already-gzipped .duckdb_extension.gz
    ext_path=$(find "$art" -type f -name '*.duckdb_extension' -print -quit)
    if [[ -z "$ext_path" ]]; then
        ext_path=$(find "$art" -type f -name '*.duckdb_extension.gz' -print -quit)
    fi
    [[ -n "$ext_path" ]] || fail "artifact $name contains no .duckdb_extension(.gz) file"

    out_dir="$RELEASE_DIR/$platform"
    mkdir -p "$out_dir"
    out_file="$out_dir/miint.duckdb_extension.gz"

    if [[ "$ext_path" == *.gz ]]; then
        cp "$ext_path" "$out_file"
    else
        gzip -9 -c "$ext_path" > "$out_file"
    fi

    size=$(stat -c %s "$out_file")
    deployed+=("$platform ($(numfmt --to=iec --suffix=B "$size"))")
    log "  staged $platform <- $(basename "$ext_path") ($size bytes gz)"
done < <(find "$ART_DIR" -mindepth 1 -maxdepth 1 -type d -print0)

[[ ${#deployed[@]} -gt 0 ]] || fail "no artifacts matched prefix $ARTIFACT_PREFIX"

# ----- push to remote + atomic symlink swap -----

REMOTE_RELEASES="$DEST_BASE/releases"
REMOTE_RELEASE="$REMOTE_RELEASES/$RELEASE_NAME"

log "Ensuring $REMOTE_RELEASES exists on $DEST_HOST"
# shellcheck disable=SC2029  # intentional local expansion of $REMOTE_RELEASES
ssh "${DEST_USER}@${DEST_HOST}" "mkdir -p '$REMOTE_RELEASES'"

log "Syncing $RELEASE_DIR/ -> ${DEST_USER}@${DEST_HOST}:${REMOTE_RELEASE}/"
rsync -az --chmod=D755,F644 \
    "$RELEASE_DIR/" \
    "${DEST_USER}@${DEST_HOST}:${REMOTE_RELEASE}/"

log "Flipping $DEST_BASE/$DUCKDB_VERSION -> releases/$RELEASE_NAME"
# ln -sfn + mv -Tf is the atomic symlink-swap idiom: ln stages a new link
# at a temp name, mv -T uses rename(2) to replace the live link in one step.
# shellcheck disable=SC2029  # intentional local expansion of $DEST_BASE etc.
ssh "${DEST_USER}@${DEST_HOST}" "
    set -euo pipefail
    cd '$DEST_BASE'
    ln -sfn 'releases/$RELEASE_NAME' '${DUCKDB_VERSION}.tmp'
    mv -Tf '${DUCKDB_VERSION}.tmp' '$DUCKDB_VERSION'
"

# ----- record state + success email -----

echo "$RUN_ID" > "$STATE_FILE"

subject="[SUCCEED] miint publish on ${HOSTNAME_SHORT} — deployed $DUCKDB_VERSION run ${RUN_ID}"
body=$(cat <<EOF
Host:          $HOSTNAME_SHORT
Repo:          $REPO @ $BRANCH
Workflow:      $WORKFLOW
DuckDB ver:    $DUCKDB_VERSION

Deployed run:
  id:    $RUN_ID
  sha:   $RUN_SHA
  title: $RUN_TITLE
  time:  $RUN_TIME
  url:   $RUN_URL

Remote release dir: ${DEST_USER}@${DEST_HOST}:${REMOTE_RELEASE}
Live symlink:       ${DEST_USER}@${DEST_HOST}:${DEST_BASE}/${DUCKDB_VERSION}

Previous run id:    ${LAST_RUN_ID:-<none>}

Platforms published (${#deployed[@]}):
$(printf '  - %s\n' "${deployed[@]}")
EOF
)
send_mail "$subject" "$body"
log "Done."

# ──────────────────────────────────────────────────────────────────────
# Example env file (save as e.g. /etc/miint-publish.env, chmod 600):
#
#   REPO=the-miint/duckdb-miint
#   BRANCH=v1.5-variegata
#   WORKFLOW=MainDistributionPipeline.yml
#   DUCKDB_VERSION=v1.5.2
#   DEST_USER=deploy
#   DEST_HOST=repo.example.com
#   DEST_BASE=/var/www/miint-ext
#   NOTIFY_EMAIL=releases@example.com
#   # optional:
#   # STATE_DIR=/var/lib/miint-publish
#
# After install:
#   chmod +x scripts/cron-publish-extension.sh
#   scripts/cron-publish-extension.sh /etc/miint-publish.env   # run once to verify
#   crontab -e
#     17 3 * * * /path/to/repo/scripts/cron-publish-extension.sh /etc/miint-publish.env
# ──────────────────────────────────────────────────────────────────────
