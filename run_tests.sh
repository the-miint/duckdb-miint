set -e

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
