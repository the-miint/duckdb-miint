#!/bin/bash
# Failures that are not a genome's: they fail the query instead of becoming a training_failed row.
# - gpl-boundary dying while training a genome (a python3 proxy in front of the real gpl-boundary exits on the first
#   prodigal-train request);
# - a training_path to be written in a directory that cannot be written (skipped as root, who can write it).
# Needs predict_genes compiled in and a gpl-boundary >= 0.5.0 (MIINT_GPL_BOUNDARY_PATH or on PATH).

set -u

DUCKDB=./build/release/duckdb
if ! echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'predict_genes';" | $DUCKDB -csv -noheader 2>/dev/null | grep -q 1; then
    echo "predict_genes not compiled in; skipping"
    exit 0
fi
REAL="${MIINT_GPL_BOUNDARY_PATH:-$(command -v gpl-boundary || true)}"
if [ -z "$REAL" ] || [ ! -x "$REAL" ] || ! command -v python3 > /dev/null; then
    echo "gpl-boundary or python3 not found; skipping"
    exit 0
fi

WORK=$(mktemp -d "${TMPDIR:-/tmp}/miint_predict_genes_daemon.XXXXXX")
trap 'rm -rf "$WORK"' EXIT
FAILED=0

cat > "$WORK/dying.py" <<'EOF'
import os, subprocess, sys, threading
real = os.environ["REAL_GPL_BOUNDARY"]
if len(sys.argv) > 1:  # --version and the like go to the real binary
    os.execv(real, [real] + sys.argv[1:])
child = subprocess.Popen([real], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
def pump():
    for line in child.stdout:
        sys.stdout.buffer.write(line)
        sys.stdout.buffer.flush()
threading.Thread(target=pump, daemon=True).start()
for line in sys.stdin.buffer:
    if b'"prodigal-train"' in line:
        child.kill()
        os._exit(3)
    child.stdin.write(line)
    child.stdin.flush()
child.stdin.close()
child.wait()
EOF
cat > "$WORK/gpl-boundary" <<EOF
#!/bin/bash
exec python3 "$WORK/dying.py" "\$@"
EOF
chmod +x "$WORK/gpl-boundary"

echo "predict_genes when gpl-boundary dies:"
out=$(printf '%s\n' "CREATE TABLE c AS SELECT 'g' AS genome, * FROM read_fastx('data/prodigal/anthus_aco.fas');" \
    "SELECT genome_id, status FROM predict_genes('c', genome_id := 'genome');" |
    REAL_GPL_BOUNDARY="$REAL" MIINT_GPL_BOUNDARY_PATH="$WORK/gpl-boundary" $DUCKDB -no-agent -csv 2>&1)
if echo "$out" | grep -q "training_failed"; then
    echo "  ✗ the failure became a training_failed row:"
    echo "$out" | sed 's/^/      /'
    FAILED=1
elif ! echo "$out" | grep -q "Error"; then
    echo "  ✗ expected an error, got:"
    echo "$out" | sed 's/^/      /'
    FAILED=1
else
    echo "  ✓ the query fails, naming the daemon's failure"
fi

# A training_path to be written in a directory that cannot be written fails the call up front, rather than reading as
# the genome failing to train.
if [ "$(id -u)" = 0 ]; then
    echo "  - running as root, which can write a 0555 directory: skipped"
    exit $FAILED
fi
mkdir "$WORK/read_only"
chmod 555 "$WORK/read_only"
out=$(printf '%s\n' "CREATE TABLE c AS SELECT 'g' AS genome, * FROM read_fastx('data/prodigal/anthus_aco.fas');" \
    "SELECT genome_id, status FROM predict_genes('c', genome_id := 'genome', training_path := '$WORK/read_only/m.trn');" |
    MIINT_GPL_BOUNDARY_PATH="$REAL" $DUCKDB -no-agent -csv 2>&1)
chmod 755 "$WORK/read_only"
expect_error() {
    if echo "$2" | grep -q "$3"; then
        echo "  ✓ $1"
    else
        echo "  ✗ $1: expected /$3/ in:"
        echo "$2" | sed 's/^/      /'
        FAILED=1
    fi
}
expect_error "an unwritable training_path directory is refused" "$out" "is not writable"

exit $FAILED
