#!/bin/bash
# predict_genes refuses a gpl-boundary older than 0.5.0, and one whose release it cannot read, before it starts
# the daemon. The fake binaries are bash scripts answering --version, so this needs predict_genes compiled in, not
# gpl-boundary installed.

set -u

DUCKDB=./build/release/duckdb
if ! echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'predict_genes';" | $DUCKDB -csv -noheader 2>/dev/null | grep -q 1; then
    echo "predict_genes not compiled in; skipping"
    exit 0
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/miint_predict_genes_version.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
FAILED=0

# A fake gpl-boundary whose --version prints `json`; anything else (daemon mode) records that it ran.
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

echo "predict_genes version gate:"

old=$(fake_binary old '{"gpl_boundary":"0.4.2","tools":[]}')
out=$(run_with "$old")
expect "0.4.2 is refused, naming both releases" "$out" "requires gpl-boundary >= 0\.5\.0.*reports version 0\.4\.2"
if [ -e "$WORK/old.daemon_started" ]; then
    echo "  ✗ the daemon was started for a refused release"
    FAILED=1
else
    echo "  ✓ the daemon is not started for a refused release"
fi

garbled=$(fake_binary garbled 'not json')
out=$(run_with "$garbled")
expect "an unreadable release is refused" "$out" "requires gpl-boundary >= 0\.5\.0.*unknown"

# A binary upgraded in place (install_gpl_boundary(true) replaces it at the same path) is seen by the next query in
# the same process: a release too old for the caller is probed again. Its --version reads the release from a file.
cat > "$WORK/upgraded" <<EOF
#!/bin/bash
if [ "\${1:-}" = --version ]; then
    printf '{"gpl_boundary":"%s","tools":[]}\n' "\$(cat "$WORK/release")"
    exit 0
fi
exit 1
EOF
chmod +x "$WORK/upgraded"
echo 0.4.2 > "$WORK/release"
# The CLI's .shell runs a command without a shell, so the upgrade is a script.
printf '#!/bin/bash\necho 0.5.0 > "%s/release"\n' "$WORK" > "$WORK/upgrade"
chmod +x "$WORK/upgrade"
out=$(printf '%s\n' "CREATE TABLE c AS SELECT 'g' AS genome, 0 AS sequence_index, 'c1' AS read_id, 'ACGTACGT' AS sequence1;" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" \
    ".shell $WORK/upgrade" \
    "SELECT 'second query';" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" |
    MIINT_GPL_BOUNDARY_PATH="$WORK/upgraded" $DUCKDB -no-agent -csv 2>&1)
expect "the first query is refused at 0.4.2" "$out" "reports version 0\.4\.2"
after=$(echo "$out" | sed -n '/second query/,$p')
if echo "$after" | grep -q "requires gpl-boundary"; then
    echo "  ✗ the upgraded release is still refused:"
    echo "$after" | sed 's/^/      /'
    FAILED=1
else
    echo "  ✓ the upgraded release passes the gate without a restart"
fi

exit $FAILED
