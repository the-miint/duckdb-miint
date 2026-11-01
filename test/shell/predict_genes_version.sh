#!/bin/bash
# predict_genes refuses a gpl-boundary older than 0.5.0, one whose release it cannot read, and one whose Prodigal is
# not 2.6.5, before it starts the daemon. The fake binaries are bash scripts answering --version, so this needs
# predict_genes compiled in, not gpl-boundary installed.

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

# A 0.5.0 or later whose Prodigal is not the one this miint names in its GFF, or that lists no version for it.
other=$(fake_binary other '{"gpl_boundary":"0.6.0","tools":[{"name":"prodigal","version":"2.6.6"},{"name":"prodigal-train","version":"2.6.6"}]}')
out=$(run_with "$other")
expect "another Prodigal is refused, naming both versions" "$out" "whose prodigal is version 2\.6\.5.*reports version 2\.6\.6"
if [ -e "$WORK/other.daemon_started" ]; then
    echo "  ✗ the daemon was started for another Prodigal"
    FAILED=1
fi
unlisted=$(fake_binary unlisted '{"gpl_boundary":"0.5.0","tools":[{"name":"prodigal","version":"2.6.5"}]}')
out=$(run_with "$unlisted")
expect "a prodigal-train with no version is refused" "$out" "whose prodigal-train is version 2\.6\.5.*reports no version for it"

# A binary upgraded in place (install_gpl_boundary(true) replaces it at the same path) is seen by the next query in
# the same process: a release too old for the caller is probed again, and so is Prodigal's version after a query that
# passed. Its --version reads the release and the Prodigal version from files.
cat > "$WORK/upgraded" <<EOF
#!/bin/bash
if [ "\${1:-}" = --version ]; then
    p=\$(cat "$WORK/prodigal")
    printf '{"gpl_boundary":"%s","tools":[{"name":"prodigal","version":"%s"},{"name":"prodigal-train","version":"%s"}]}\n' "\$(cat "$WORK/release")" "\$p" "\$p"
    exit 0
fi
exit 1
EOF
chmod +x "$WORK/upgraded"
echo 0.4.2 > "$WORK/release"
echo 2.6.5 > "$WORK/prodigal"
# The CLI's .shell runs a command without a shell, so the upgrade is a script.
printf '#!/bin/bash\necho 0.5.0 > "%s/release"\n' "$WORK" > "$WORK/upgrade"
chmod +x "$WORK/upgrade"
printf '#!/bin/bash\necho 2.6.6 > "%s/prodigal"\n' "$WORK" > "$WORK/swap"
chmod +x "$WORK/swap"
out=$(printf '%s\n' "CREATE TABLE c AS SELECT 'g' AS genome, 0 AS sequence_index, 'c1' AS read_id, 'ACGTACGT' AS sequence1;" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" \
    ".shell $WORK/upgrade" \
    "SELECT 'second query';" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" \
    ".shell $WORK/swap" \
    "SELECT 'third query';" \
    "SELECT * FROM predict_genes('c', genome_id := 'genome', meta_mode := true);" |
    MIINT_GPL_BOUNDARY_PATH="$WORK/upgraded" $DUCKDB -no-agent -csv 2>&1)
expect "the first query is refused at 0.4.2" "$out" "reports version 0\.4\.2"
after=$(echo "$out" | sed -n '/second query/,/third query/p')
if echo "$after" | grep -q "requires gpl-boundary"; then
    echo "  ✗ the upgraded release is still refused:"
    echo "$after" | sed 's/^/      /'
    FAILED=1
else
    echo "  ✓ the upgraded release passes the gate without a restart"
fi
expect "a Prodigal swapped in after a passing query is refused" "$(echo "$out" | sed -n '/third query/,$p')" \
    "whose prodigal is version 2\.6\.5.*reports version 2\.6\.6"

exit $FAILED
