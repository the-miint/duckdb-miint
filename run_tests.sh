set -e

# Detect optional external tools
if command -v bowtie2 &> /dev/null; then
    export BOWTIE2_AVAILABLE=1
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
