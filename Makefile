PROJ_DIR := $(dir $(abspath $(lastword $(MAKEFILE_LIST))))

# Configuration of extension
EXT_NAME=miint
EXT_CONFIG=${PROJ_DIR}extension_config.cmake

# Include the Makefile from extension-ci-tools
include extension-ci-tools/makefiles/duckdb_extension.Makefile

# Extend clean to also remove BUILD_IN_SOURCE artifacts from ext/ directories.
# minimap2 and WFA2-lib build in their source trees, so their .o and .a files
# survive the default `rm -rf build`.
clean: clean-ext

.PHONY: clean-ext
clean-ext:
	-make -C $(PROJ_DIR)ext/minimap2 clean 2>/dev/null
	-make -C $(PROJ_DIR)ext/WFA2-lib clean 2>/dev/null
