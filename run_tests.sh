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
