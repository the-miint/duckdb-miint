#!/bin/bash
# read_ena_sequences over Aspera, with test/shell/read_ena_sequences_aspera_stubs/ascp
# standing in for ascp:
#
# - A transfer that works: the rows match the file sent, single- and paired-end,
#   and filepath names the Aspera path.
# - ascp fails before sending any data, as when the server rejects the key:
#   - download_method => 'aspera' raises an error carrying ascp's message,
#     also when DuckDB runs with stdin closed or under a terminal;
#   - download_method => 'auto' (the default) warns and reads over HTTP, and
#     the scan's later runs go straight to HTTP;
#   - an ascp that exits 0 without sending anything fails the same way;
#   - a scan of several runs that all fail raises instead of returning no rows,
#     and so does one run that fails partway, instead of returning part of it.
# - With ASPERA_AVAILABLE, the real ascp neither hangs nor fails silently.
#
# None may hang: a prompt that reaches stdin or a terminal waits forever, so
# every query runs under a hard timeout, and a hang reads as exit 137. Needs
# live ENA for run resolution, and httpfs for the HTTPS reads.

set -u

DUCKDB="./build/release/duckdb"
HERE="$(cd "$(dirname "$0")" && pwd)"
STUBS="$HERE/read_ena_sequences_aspera_stubs"
ROOT="$(cd "$HERE/../.." && pwd)"

if [ -z "${ENA_AVAILABLE:-}" ]; then
    echo "ENA_AVAILABLE not set, skipping Aspera tests"
    exit 0
fi
if ! "$DUCKDB" -c "LOAD httpfs; SELECT 1;" > /dev/null 2>&1; then
    echo "httpfs not available, skipping Aspera tests"
    exit 0
fi
TIMEOUT_CMD="$(command -v timeout || command -v gtimeout || true)"
if [ -z "$TIMEOUT_CMD" ]; then
    echo "ERROR: no timeout command; a regression here hangs instead of failing" >&2
    exit 1
fi
HAVE_SCRIPT=""
if script --version 2>/dev/null | grep -q util-linux; then HAVE_SCRIPT=1; fi

REAL_PATH="$PATH"
export PATH="$STUBS:$PATH"
# Key lookup tries $CONDA_PREFIX/etc first; the placeholder there keeps it from
# downloading the real key, so the test does not depend on that download.
export CONDA_PREFIX="$STUBS"
FAILED=0
ESC="$(printf '\033')"

# run_query SQL [ENV_ARG ...]: sets OUT (stdout + stderr) and RC. The ENV_ARGs
# go to env(1), e.g. VAR=value or -u VAR. CLOSE_STDIN=1 runs DuckDB with stdin
# closed; PTY=1 runs it under a pseudo-terminal (util-linux script).
run_query() {
    local sql="LOAD httpfs; $1"
    shift
    RC=0
    if [ -n "${PTY:-}" ]; then
        local cmd
        cmd="$(printf '%q ' env "$@" "$DUCKDB" -dark-mode -noheader -list -c "$sql")"
        OUT="$("$TIMEOUT_CMD" -s KILL 120 script -qec "$cmd" /dev/null 2>&1 < /dev/null)" || RC=$?
        OUT="$(printf '%s\n' "$OUT" | tr -d '\r')"
    elif [ -n "${CLOSE_STDIN:-}" ]; then
        OUT="$(env "$@" "$TIMEOUT_CMD" -s KILL 120 "$DUCKDB" -noheader -list -c "$sql" 2>&1 0<&-)" || RC=$?
    else
        OUT="$(env "$@" "$TIMEOUT_CMD" -s KILL 120 "$DUCKDB" -noheader -list -c "$sql" 2>&1)" || RC=$?
    fi
    # The shell colours warnings, which leaves escape codes on the result line.
    OUT="$(printf '%s\n' "$OUT" | sed "s/${ESC}\[[0-9;]*m//g")"
}

# error_has TEXT: the raised error, not a warning that quotes it, contains TEXT.
# Under a terminal the CLI puts the message on the line after "IO Error:".
error_has() {
    echo "$OUT" | awk '/^[A-Za-z ]*Error:/ { print; getline; print }' | grep -qF "$1"
}

# fell_back_runs: how many runs the printed fallback warnings name. Counted on
# stderr because the CLI keeps no warnings for miint_warnings() (its log
# storage prints them), and it prints each one twice.
fell_back_runs() {
    echo "$OUT" | grep "falling back to HTTP" | grep -o "for run '[^']*'" | sort -u | wc -l
}

check() {
    local name="$1" ok="$2"
    if [ "$ok" = 1 ]; then
        echo "  PASS: $name"
    else
        echo "  FAIL: $name (exit $RC)"
        while IFS= read -r line; do echo "    $line"; done <<<"$OUT"
        FAILED=1
    fi
}

R1="$ROOT/data/merge/ltp_50_R1.fq"
R2="$ROOT/data/merge/ltp_50_R2.fq"
DATA=(MIINT_ASCP_STUB_MODE=data "MIINT_ASCP_STUB_DATA_R1=$R1" "MIINT_ASCP_STUB_DATA_R2=$R2")

# Every byte must arrive once, in order, including the first, which is read
# early to tell a transfer that started from one that failed. With the count,
# matching sequence_index on both sides pins each read's position.
echo "Running: an Aspera transfer that works, single-end"
run_query "CREATE TABLE got AS SELECT * FROM read_ena_sequences('ERR1074767', download_method => 'aspera', include_filepath => true);
SELECT (SELECT count(*) FROM got) = 50,
       (SELECT count(*) FROM (SELECT sequence_index, read_id, sequence1, qual1 FROM got
                              EXCEPT SELECT sequence_index, read_id, sequence1, qual1 FROM read_fastx('$R1'))) = 0,
       (SELECT bool_and(filepath LIKE 'fasp.sra.ebi.ac.uk:%') FROM got);" "${DATA[@]}"
ok=0
if [ "$RC" -eq 0 ] && echo "$OUT" | grep -qx "true|true|true"; then ok=1; fi
check "the file's 50 reads, and an Aspera filepath" "$ok"

echo "Running: an Aspera transfer that works, paired-end"
run_query "CREATE TABLE got AS SELECT * FROM read_ena_sequences('ERR10682882', download_method => 'aspera');
SELECT (SELECT count(*) FROM got) = 50,
       (SELECT count(*) FROM (SELECT sequence_index, read_id, sequence1, qual1, sequence2, qual2 FROM got
                              EXCEPT SELECT sequence_index, read_id, sequence1, qual1, sequence2, qual2
                              FROM read_fastx('$R1', sequence2 => '$R2'))) = 0;" "${DATA[@]}"
ok=0
if [ "$RC" -eq 0 ] && echo "$OUT" | grep -qx "true|true"; then ok=1; fi
check "the pair's 50 reads" "$ok"

# A single-run scan raises the run's own error, not the every-run-failed one.
echo "Running: download_method => 'aspera' raises ascp's error"
run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera', max_sequences => 10);"
ok=0
if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "failed to authenticate" && ! error_has "every run failed"; then
    ok=1
fi
check "the run's error names ascp's failure" "$ok"

# filepath names where the rows came from: the HTTPS URL, not the Aspera path.
echo "Running: download_method => 'auto' falls back to HTTP"
run_query "SELECT COUNT(*), bool_and(filepath LIKE 'https://%') FROM read_ena_sequences('ERR1074767', include_filepath => true, max_sequences => 10);"
ok=0
if [ "$RC" -eq 0 ] && echo "$OUT" | grep -qx "10|true" && [ "$(fell_back_runs)" -eq 1 ]; then ok=1; fi
check "10 rows over HTTP, HTTPS filepaths, one warning" "$ok"

# A parent with fd 0 closed hands the lowest free descriptor to whatever the
# child opens; the child's stdin must still be /dev/null.
echo "Running: download_method => 'aspera' with DuckDB's stdin closed"
CLOSE_STDIN=1 run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera', max_sequences => 10);"
ok=0
if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "failed to authenticate"; then ok=1; fi
check "the run's error names ascp's failure" "$ok"

# Under a terminal, ascp asks on /dev/tty instead; it must not have one.
if [ -n "$HAVE_SCRIPT" ]; then
    echo "Running: download_method => 'aspera' under a terminal"
    PTY=1 run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera', max_sequences => 10);"
    ok=0
    if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "failed to authenticate"; then ok=1; fi
    check "the run's error names ascp's failure" "$ok"
else
    echo "  SKIP: under a terminal (needs util-linux script)"
fi

# One thread, so the second run opens after the first has failed over.
echo "Running: download_method => 'auto' over two runs falls back once"
run_query "SET threads = 1; SELECT COUNT(*) FROM read_ena_sequences(['ERR1074767', 'ERR10360000'], max_sequences => 10);"
ok=0
if [ "$RC" -eq 0 ] && echo "$OUT" | grep -qx "20" && [ "$(fell_back_runs)" -eq 1 ]; then ok=1; fi
check "20 rows over HTTP, one warning" "$ok"

echo "Running: an ascp that exits 0 without sending data"
run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera', max_sequences => 10);" \
    MIINT_ASCP_STUB_MODE=empty
ok=0
if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "exited with code 0 before sending any data"; then ok=1; fi
check "the run's error says no data came" "$ok"

# Each run is skipped after its retry, which leaves nothing: an error, not 0 rows.
echo "Running: download_method => 'aspera' over two runs that both fail"
run_query "SELECT COUNT(*) FROM read_ena_sequences(['ERR1074767', 'ERR10360000'], download_method => 'aspera', max_sequences => 10);"
ok=0
if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "every run failed" && error_has "failed to authenticate"; then ok=1; fi
check "the error says every run failed, with ascp's message" "$ok"

# Rows from the first half went downstream before the break; the scan must not
# end as if that were the whole run.
echo "Running: a single-run transfer that breaks partway"
run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera');" \
    MIINT_ASCP_STUB_MODE=truncated "MIINT_ASCP_STUB_DATA_R1=$R1"
ok=0
if [ "$RC" -ne 0 ] && [ "$RC" -ne 137 ] && error_has "every run failed (1): ERR1074767"; then ok=1; fi
check "the error says the run failed" "$ok"

# The real ascp and key, against EBI: rows if the key works, ascp's error if it
# does not, and never a hang, whatever DuckDB's stdin is.
if [ -z "${ASPERA_AVAILABLE:-}" ]; then
    echo "  SKIP: the real ascp (ASPERA_AVAILABLE not set)"
else
    for how in open closed terminal; do
        if [ "$how" = terminal ] && [ -z "$HAVE_SCRIPT" ]; then
            echo "  SKIP: the real ascp under a terminal (needs util-linux script)"
            continue
        fi
        echo "Running: the real ascp, stdin $how"
        CLOSE_STDIN="$([ "$how" = closed ] && echo 1)" PTY="$([ "$how" = terminal ] && echo 1)" \
            run_query "SELECT COUNT(*) FROM read_ena_sequences('ERR1074767', download_method => 'aspera', max_sequences => 10);" \
            -u CONDA_PREFIX "PATH=$REAL_PATH"
        ok=0
        if [ "$RC" -ne 137 ] && { echo "$OUT" | grep -qx "10" || error_has "ascp"; }; then ok=1; fi
        check "10 rows or ascp's error, not a hang" "$ok"
    done
fi

if [ "$FAILED" -ne 0 ]; then
    exit 1
fi
echo "All Aspera tests passed"
