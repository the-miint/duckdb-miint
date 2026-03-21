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
