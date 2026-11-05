#!/bin/bash
# predict_genes refuses a gpl-boundary whose prodigal or prodigal-train is not Prodigal 2.6.5 (the version
# COPY FORMAT GFF's DIALECT 'prodigal' names), or does not say, before it starts the daemon; and checks again on
# every call, so a binary swapped in place is seen. (The release gate itself: gpl_boundary_version.sh.) The fake binaries are bash
# scripts answering --version, so this needs predict_genes compiled in, not gpl-boundary installed.

set -u

DUCKDB=./build/release/duckdb
if ! echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'predict_genes';" | $DUCKDB -csv -noheader 2>/dev/null | grep -q 1; then
    echo "predict_genes not compiled in; skipping"
    exit 0
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/miint_predict_genes_prodigal_version.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
FAILED=0

# A gpl-boundary of a release new enough, whose --version prints `json`; anything else (daemon mode) records that it
# ran.
fake_binary() {
    local path="$WORK/$1" json="$2"
    cat > "$path" <<EOF
#!/bin/bash
if [ "\${1:-}" = --version ]; then
    printf '%s\n' '$json'
    exit 0
fi
touch "$WORK/$1.daemon_started"
exit 1
EOF
    chmod +x "$path"
    echo "$path"
}

SQL="CREATE TABLE c AS SELECT 'g' AS genome, 0 AS sequence_index, 'c1' AS read_id, 'ACGTACGT' AS sequence1;
SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);"

run_with() {
    printf '%s\n' "$SQL" | MIINT_GPL_BOUNDARY_PATH="$1" $DUCKDB -no-agent -csv 2>&1
}

expect() {
    local name="$1" out="$2" pattern="$3"
    if echo "$out" | grep -qE -- "$pattern"; then
        echo "  ✓ $name"
    else
        echo "  ✗ $name: expected /$pattern/ in:"
        echo "$out" | sed 's/^/      /'
        FAILED=1
    fi
}

echo "predict_genes Prodigal version:"

other=$(fake_binary other '{"gpl_boundary":"0.6.0","tools":[{"name":"prodigal","version":"2.6.6"},{"name":"prodigal-train","version":"2.6.6"}]}')
out=$(run_with "$other")
expect "another Prodigal is refused, naming both versions" "$out" "whose prodigal is version 2\.6\.5.*reports version 2\.6\.6"
if [ -e "$WORK/other.daemon_started" ]; then
    echo "  ✗ the daemon was started for another Prodigal"
    FAILED=1
else
    echo "  ✓ the daemon is not started for another Prodigal"
fi

unlisted=$(fake_binary unlisted '{"gpl_boundary":"0.5.0","tools":[{"name":"prodigal","version":"2.6.5"}]}')
out=$(run_with "$unlisted")
expect "a prodigal-train with no version is refused" "$out" "whose prodigal-train is version 2\.6\.5.*reports no version for it"

# A Prodigal swapped in after a query that passed is seen by the next query in the same process. The fake's
# --version reads the Prodigal version from a file; the daemon it does not start fails the passing query, after
# the check.
cat > "$WORK/swapped" <<EOF
#!/bin/bash
if [ "\${1:-}" = --version ]; then
    p=\$(cat "$WORK/prodigal")
    printf '{"gpl_boundary":"0.5.0","tools":[{"name":"prodigal","version":"%s"},{"name":"prodigal-train","version":"%s"}]}\n' "\$p" "\$p"
    exit 0
fi
exit 1
EOF
chmod +x "$WORK/swapped"
echo 2.6.5 > "$WORK/prodigal"
# The CLI's .shell runs a command without a shell, so the swap is a script.
printf '#!/bin/bash\necho 2.6.6 > "%s/prodigal"\n' "$WORK" > "$WORK/swap"
chmod +x "$WORK/swap"
out=$(printf '%s\n' "CREATE TABLE c AS SELECT 'g' AS genome, 0 AS sequence_index, 'c1' AS read_id, 'ACGTACGT' AS sequence1;" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" \
    ".shell $WORK/swap" \
    "SELECT 'second query';" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" |
    MIINT_GPL_BOUNDARY_PATH="$WORK/swapped" $DUCKDB -no-agent -csv 2>&1)
before=$(echo "$out" | sed -n '1,/second query/p')
if echo "$before" | grep -q "whose prodigal"; then
    echo "  ✗ Prodigal 2.6.5 is refused:"
    echo "$before" | sed 's/^/      /'
    FAILED=1
else
    echo "  ✓ Prodigal 2.6.5 passes"
fi
expect "a Prodigal swapped in after a passing query is refused" "$(echo "$out" | sed -n '/second query/,$p')" \
    "whose prodigal is version 2\.6\.5.*reports version 2\.6\.6"

exit $FAILED
