PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Configuration of extension
EXT_NAME=miint
EXT_CONFIG=${PROJ_DIR}extension_config.cmake

# Include the Makefile from extension-ci-tools
include extension-ci-tools/makefiles/duckdb_extension.Makefile

# Override tidy-check to include VCPKG_MANIFEST_FLAGS (needed for HDF5 dependency)
# Run vcpkg install first to ensure dependencies are available
tidy-check:
	@if [ -n "$(VCPKG_TOOLCHAIN_PATH)" ] && [ -f "vcpkg/vcpkg" ]; then \
		echo "Installing vcpkg dependencies..."; \
		./vcpkg/vcpkg install; \
	fi
	mkdir -p ./build/tidy
	cmake $(GENERATOR) $(BUILD_FLAGS) $(EXT_DEBUG_FLAGS) $(VCPKG_MANIFEST_FLAGS) -DDISABLE_UNITY=1 -DCLANG_TIDY=1 -S $(DUCKDB_SRCDIR) -B build/tidy
	cp duckdb/.clang-tidy build/tidy/.clang-tidy
	cd build/tidy && python3 ../../duckdb/scripts/run-clang-tidy.py '$(PROJ_DIR)src/.*/' -header-filter '$(PROJ_DIR)src/.*/' -quiet ${TIDY_THREAD_PARAMETER} ${TIDY_BINARY_PARAMETER} ${TIDY_PERFORM_CHECKS}