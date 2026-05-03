set -e

# Start local HTTP server for HTTPS reader tests (unless already set externally)
HTTP_SERVER_PID=""
if [ -z "$MIINT_HTTPS_TEST_URL" ]; then
    # Find a free port: bind to port 0, read assigned port, close immediately
    HTTPS_TEST_PORT=$(python3 -c "import socket; s=socket.socket(); s.bind(('',0)); print(s.getsockname()[1]); s.close()")
    python3 -m http.server "$HTTPS_TEST_PORT" --directory data/ > /dev/null 2>&1 &
    HTTP_SERVER_PID=$!
    trap "kill $HTTP_SERVER_PID 2>/dev/null || true" EXIT

    # Wait for server to be ready (up to 3 seconds)
    HTTP_SERVER_READY=0
    for i in $(seq 1 15); do
        if curl -s "http://localhost:$HTTPS_TEST_PORT/" > /dev/null 2>&1; then
            HTTP_SERVER_READY=1
            break
        fi
        sleep 0.2
    done

    if [ "$HTTP_SERVER_READY" -ne 1 ]; then
        echo "ERROR: HTTP test server failed to start on port $HTTPS_TEST_PORT"
        exit 1
    fi

    export MIINT_HTTPS_TEST_URL="http://localhost:$HTTPS_TEST_PORT"
fi

# Detect optional external tools
if command -v bowtie2 &> /dev/null; then
    export BOWTIE2_AVAILABLE=1
fi

if command -v ascp &> /dev/null; then
    export ASPERA_AVAILABLE=1
fi

# ENA network reachability. Probes the Portal API with a cheap HEAD-ish call
# so tests that need live ENA access can skip gracefully in offline CI.
# Timeout is tight (3s) to avoid slowing down the common case.
if curl -sSf --max-time 3 -o /dev/null \
        "https://www.ebi.ac.uk/ena/portal/api/search?result=read_run&query=run_accession%3D%22ERR1074767%22&fields=run_accession&limit=1&format=tsv" \
        2>/dev/null; then
    export ENA_AVAILABLE=1
fi

# ENA Webin V2 submission test endpoint. Live-integration tests for
# INSERT INTO ena.{projects,samples,...} require a real Webin account
# (free, registered at https://www.ebi.ac.uk/ena/submit/webin/). When the
# user provides credentials AND the dev/test server responds, we export
# ENA_WEBIN_TEST_AVAILABLE so guarded SQL tests can run.
if [ -n "$ENA_WEBIN_TEST_USER" ] && [ -n "$ENA_WEBIN_TEST_PASSWORD" ]; then
    export ENA_WEBIN_TEST_USER
    export ENA_WEBIN_TEST_PASSWORD
    if curl -sSf --max-time 3 -o /dev/null \
            "https://wwwdev.ebi.ac.uk/ena/submit/webin-v2/swagger-ui/index.html" \
            2>/dev/null; then
        export ENA_WEBIN_TEST_AVAILABLE=1
    fi
fi

# ENA Webin V2 mock server. Provides canned XML receipts for INSERT INTO
# ena.{projects,samples,...} round-trip tests with no live credentials. The
# helper python script is in test/scripts/ena_webin_mock.py.
ENA_MOCK_PID=""
if [ -z "$ENA_WEBIN_MOCK_URL" ] && python3 -c "import http.server" 2>/dev/null; then
    python3 test/scripts/ena_webin_mock.py 0 > /tmp/ena_webin_mock.boot 2>&1 &
    ENA_MOCK_PID=$!
    # Inherit the existing trap (HTTP_SERVER_PID) and add the mock kill
    trap "kill $HTTP_SERVER_PID $ENA_MOCK_PID 2>/dev/null || true" EXIT
    for i in $(seq 1 25); do
        if [ -s /tmp/ena_webin_mock.boot ] && grep -q ENA_WEBIN_MOCK_URL /tmp/ena_webin_mock.boot; then
            export ENA_WEBIN_MOCK_URL="$(grep ENA_WEBIN_MOCK_URL /tmp/ena_webin_mock.boot | head -n1 | cut -d= -f2-)"
            break
        fi
        sleep 0.1
    done
    if [ -z "$ENA_WEBIN_MOCK_URL" ]; then
        echo "Warning: ENA Webin mock failed to start; submission round-trip tests will skip"
        kill $ENA_MOCK_PID 2>/dev/null || true
    fi
fi

if conda run -n massql python3 -c "from massql import msql_engine" 2>/dev/null; then
    export MASSQL_PYTHON_AVAILABLE=1
fi

# Download and cache MassQL test data if not present
CACHE_DIR="data/cache"
mkdir -p "$CACHE_DIR"
MZML_FILE="$CACHE_DIR/bld_plt1_07_120_1.mzML"
EXPECTED_SHA="9be69d251cbf7d9b3c439ccc76dd65f055ddac96b946fc2f565e9a95bfd4ca46"

if [ ! -f "$MZML_FILE" ]; then
    echo "Downloading MassQL test data..."
    if ! wget -q -O "$MZML_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000085944/ccms_peak/raw_data/bld_plt1_07_120_1.mzML"; then
        rm -f "$MZML_FILE"
        echo "Warning: MassQL test data download failed, skipping MassQL tests"
    fi
fi
if echo "$EXPECTED_SHA  $MZML_FILE" | sha256sum -c --quiet 2>/dev/null; then
    export MASSQL_TEST_DATA="$MZML_FILE"
fi

# Download and cache GNPS test data if not present
GNPS_FILE="$CACHE_DIR/GNPS00002_A3_p.mzML"
if [ ! -f "$GNPS_FILE" ]; then
    echo "Downloading GNPS test data..."
    if ! wget -q -O "$GNPS_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?forceDownload=true&file=f.MSV000084494/ccms_peak/raw/GNPS00002_A3_p.mzML"; then
        rm -f "$GNPS_FILE"
        echo "Warning: GNPS test data download failed, skipping GNPS tests"
    fi
fi
if [ -f "$GNPS_FILE" ]; then
    export MASSQL_GNPS_DATA=1
fi

# Download and cache real mzXML test data if not present
MZXML_FILE="$CACHE_DIR/10125_P1_RA12_01_49.mzXML"
MZXML_EXPECTED_SHA="d026938157e5f640568c371ed6a0bca0de34e8c181da1443298671005e7ffe5a"
if [ ! -f "$MZXML_FILE" ]; then
    echo "Downloading mzXML test data..."
    if ! wget -q -O "$MZXML_FILE" "https://massive.ucsd.edu/ProteoSAFe/DownloadResultFile?file=f.MSV000080179%2Fccms_peak%2FmzXMLs%2F10125_P1_RA12_01_49.mzXML&forceDownload=true"; then
        rm -f "$MZXML_FILE"
        echo "Warning: mzXML test data download failed, skipping mzXML real data tests"
    fi
fi
if echo "$MZXML_EXPECTED_SHA  $MZXML_FILE" | sha256sum -c --quiet 2>/dev/null; then
    export MZXML_REAL_DATA="$MZXML_FILE"
fi

# Detect compile-time optional features by querying the built extension
if echo "SELECT * FROM miint_versions() WHERE library = 'HDF5';" | ./build/release/duckdb -csv 2>/dev/null | grep -q HDF5; then
    export HDF5_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'detect_chimera_uchime';" | ./build/release/duckdb -csv 2>/dev/null | grep -q 1; then
    export VSEARCH_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'align_mafft';" | ./build/release/duckdb -csv 2>/dev/null | grep -q 1; then
    export MAFFT_AVAILABLE=1
fi
if echo "SELECT 1 FROM duckdb_functions() WHERE function_name = 'align_sortmerna_rrna';" | ./build/release/duckdb -csv 2>/dev/null | grep -q 1; then
    export SORTMERNA_AVAILABLE=1
fi

# Phase 5: real-data regression oracle.
#
# The BLAST-format output of `sortmerna 4.4.0 --blast '1 cigar qcov'` on the
# in-tree 100k-read environmental amplicons vs. gg_13_8_ref_set.fasta is
# checked in at data/sortmerna/real_oracle.blast.gz (1.9 MB gzipped). The
# submodule SHA that produced it is recorded alongside so we can detect when
# the oracle has gone stale.
#
# To regenerate after bumping ext/sortmerna: remove both files and rerun
# run_tests.sh with MIINT_SORTMERNA_REAL_DATA=1. That rebuilds the oracle
# using the submodule-built native CLI (the exact binary our library links),
# so the oracle and our embedded library always run the same version.
SMR_NATIVE_CLI="./build/release/extension/miint/sortmerna_build/src/sortmerna/sortmerna"
SMR_REAL_QUERY="ext/sortmerna/data/set2_environmental_study_550_amplicon.fasta"
SMR_REAL_REF="ext/sortmerna/data/gg_13_8_ref_set.fasta"
SMR_ORACLE_GZ="data/sortmerna/real_oracle.blast.gz"
SMR_ORACLE_SHA="data/sortmerna/real_oracle.submodule.sha"
SMR_SUBMODULE_HEAD=$(git -C ext/sortmerna rev-parse HEAD 2>/dev/null || echo "unknown")

if [ "${MIINT_SORTMERNA_REAL_DATA:-0}" = "1" ] && [ -x "$SMR_NATIVE_CLI" ] \
        && [ -f "$SMR_REAL_QUERY" ] && [ -f "$SMR_REAL_REF" ]; then
    echo "Regenerating sortmerna real-data oracle via native 4.4.0 CLI..."
    SMR_TMP=$(mktemp -d)
    if "$SMR_NATIVE_CLI" \
            --ref "$SMR_REAL_REF" --reads "$SMR_REAL_QUERY" \
            --workdir "$SMR_TMP" --aligned "$SMR_TMP/aligned" \
            --blast '1 cigar qcov' --num_alignments 1 --threads 4 > "$SMR_TMP/smr.log" 2>&1 \
            && [ -f "$SMR_TMP/aligned.blast" ]; then
        mkdir -p data/sortmerna
        gzip -c "$SMR_TMP/aligned.blast" > "$SMR_ORACLE_GZ"
        echo "$SMR_SUBMODULE_HEAD" > "$SMR_ORACLE_SHA"
        echo "Oracle regenerated at $SMR_ORACLE_GZ (submodule $SMR_SUBMODULE_HEAD)"
        rm -rf "$SMR_TMP"
    else
        # Preserve the CLI log for debugging — the tmpdir is about to go.
        SMR_FAIL_LOG="$(pwd)/smr_oracle_gen_fail.log"
        cp "$SMR_TMP/smr.log" "$SMR_FAIL_LOG" 2>/dev/null || true
        echo "Warning: sortmerna oracle generation failed; log saved to $SMR_FAIL_LOG"
        rm -rf "$SMR_TMP"
    fi
fi

# Export the oracle path only if both the gzipped file and its SHA record
# exist, AND the recorded SHA matches the current submodule HEAD. A mismatch
# means the submodule was bumped without regenerating — skip the test rather
# than silently validate against stale goldens.
if [ -f "$SMR_ORACLE_GZ" ] && [ -f "$SMR_ORACLE_SHA" ]; then
    SMR_RECORDED_SHA=$(cat "$SMR_ORACLE_SHA")
    if [ "$SMR_RECORDED_SHA" = "$SMR_SUBMODULE_HEAD" ]; then
        export SORTMERNA_REAL_ORACLE="$SMR_ORACLE_GZ"
    else
        echo "Warning: sortmerna oracle ($SMR_RECORDED_SHA) does not match submodule HEAD ($SMR_SUBMODULE_HEAD); skipping real-data test"
    fi
fi

make test
./build/release/extension/miint/tests

# Run shell script tests
echo "Running shell script tests..."
for test_script in test/shell/*.sh; do
    if [ -f "$test_script" ]; then
        echo "Running $test_script..."
        if ! bash "$test_script"; then
            echo "Shell test failed: $test_script"
            exit 1
        fi
    fi
done
echo "All shell script tests passed!"
