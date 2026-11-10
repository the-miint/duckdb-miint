#!/bin/bash
# A GTDB release download killed partway leaves no file at its final name in the cache.
# The readers take a file at its final name as complete (src/read_gtdb.cpp,
# DownloadToCache), so a cut one there would be read as the release until
# refresh := true. A download streams to a file of its own and is renamed into place
# only once complete; this kills the process mid-download and checks that.
#
# Requires MIINT_HTTPS_TEST_URL (run_tests.sh sets it): its /slow/ paths send a
# download's first 600 bytes, then a byte every 100 ms (test/scripts/range_http_server.py).

set -u

DUCKDB="./build/release/duckdb"

if [ -z "${MIINT_HTTPS_TEST_URL:-}" ]; then
    echo "MIINT_HTTPS_TEST_URL not set, skipping the GTDB download test"
    exit 0
fi
if ! $DUCKDB -c "LOAD httpfs; SELECT 1;" > /dev/null 2>&1; then
    echo "httpfs not available, skipping the GTDB download test"
    exit 0
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/miint_gtdb_download.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
CACHE="$WORK/cache"
FINAL="$CACHE/232.0/qc_failed_r232.tsv"

read_release() {
    $DUCKDB -csv -noheader -c "LOAD httpfs;
        SET miint_gtdb_url = '$1/gtdb/releases';
        SET miint_gtdb_cache_dir = '$CACHE';
        SELECT count(*) FROM read_gtdb_qc_failed(release := '232');"
}

read_release "$MIINT_HTTPS_TEST_URL/slow" > /dev/null 2>&1 &
PID=$!
# The download has started once its partial file exists.
started=0
for _ in $(seq 1 100); do
    if compgen -G "$CACHE/232.0/*.part" > /dev/null; then
        started=1
        break
    fi
    sleep 0.1
done
kill -9 "$PID" 2> /dev/null
wait "$PID" 2> /dev/null

FAILED=0
if [ "$started" != 1 ]; then
    echo "FAIL: the download never started (no partial file in $CACHE/232.0)"
    FAILED=1
elif [ -e "$FINAL" ]; then
    echo "FAIL: a download killed partway left $FINAL"
    FAILED=1
else
    echo "PASS: a download killed partway left no file at its final name"
fi

# The next read downloads the release whole.
result=$(read_release "$MIINT_HTTPS_TEST_URL" 2>&1)
if [ "$result" = "3" ] && cmp -s "$FINAL" data/gtdb/releases/release232/232.0/auxillary_files/qc_failed_r232.tsv; then
    echo "PASS: the next read downloaded the release"
else
    echo "FAIL: the next read gave '$result', or a cached file unlike the server's"
    FAILED=1
fi

exit $FAILED
