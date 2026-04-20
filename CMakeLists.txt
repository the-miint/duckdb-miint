cmake_minimum_required(VERSION 3.15)
include(FetchContent)   # <--- REQUIRED for FetchContent_Declare!

# Set extension name here
set(TARGET_NAME miint)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# from stobealign, seem needed for kseq++
# note: bumped from c++17 to 20
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# DuckDB's extension distribution supports vcpkg. As such, dependencies can be added in ./vcpkg.json and then
# used in cmake with find_package. Feel free to remove or replace with other dependencies.
# Note that it should also be removed from vcpkg.json to prevent needlessly installing it..
find_package(ZLIB REQUIRED)
find_package(Threads REQUIRED)
find_package(expat CONFIG REQUIRED)

# Optional feature flags for WASM-compatible builds
option(MIINT_ENABLE_HDF5 "Build with HDF5/BIOM support" ON)
option(MIINT_ENABLE_BOWTIE2 "Build with Bowtie2 alignment support" ON)
option(MIINT_ENABLE_MAFFT "Build with MAFFT multiple sequence alignment support" ON)
option(MIINT_ENABLE_SORTMERNA "Build with SortMeRNA rRNA alignment support" ON)

# Auto-disable platform-incompatible features.
# Bowtie2 requires POSIX fork/exec (unavailable on WASM and Windows).
# MAFFT uses mkdtemp and other POSIX APIs (unavailable on Windows; segfaults on MinGW).
# HDF5 C++ static class members (e.g. LinkCreatPropList::DEFAULT) become
# unresolvable GOT.mem imports in WASM, causing "bad export type: undefined"
# at extension load time.
# SortMeRNA library builds on all three OSes (its own CI matrix covers windows/
# linux/macos). Only disable on WASM.
if(EMSCRIPTEN)
    set(MIINT_ENABLE_BOWTIE2 OFF)
    set(MIINT_ENABLE_HDF5 OFF)
    set(MIINT_ENABLE_SORTMERNA OFF)
    message(STATUS "Emscripten detected: disabling Bowtie2, HDF5, and SortMeRNA (WASM-incompatible)")
elseif(WIN32)
    set(MIINT_ENABLE_BOWTIE2 OFF)
    set(MIINT_ENABLE_MAFFT OFF)
    message(STATUS "Windows detected: disabling Bowtie2 and MAFFT (require POSIX APIs)")
endif()

# HDF5 is required for BIOM support, but can be disabled entirely or made optional for clang-tidy
if(MIINT_ENABLE_HDF5)
    if(CLANG_TIDY)
        find_package(HDF5 COMPONENTS CXX)
        if(HDF5_FOUND)
            message(STATUS "HDF5 found for tidy-check")
        else()
            message(STATUS "HDF5 not found for tidy-check - BIOM support will be excluded from analysis")
        endif()
    else()
        find_package(HDF5 REQUIRED COMPONENTS CXX)
    endif()
else()
    set(HDF5_FOUND OFF)
    message(STATUS "HDF5/BIOM support disabled (MIINT_ENABLE_HDF5=OFF)")
endif()

# Detect threading flags for Emscripten wasm_threads builds.
# DuckDB passes -pthread in CMAKE_CXX_FLAGS for wasm_threads but not wasm_eh/wasm_mvp.
# ExternalProject builds (htslib, minimap2, WFA2) use explicit CFLAGS and don't
# inherit CMAKE_CXX_FLAGS, so we propagate -pthread to their compile flags here.
# Without this, object files like kthread.o lack atomics/bulk-memory features and
# the linker rejects them when --shared-memory is used.
# Note: CMAKE_CXX_FLAGS may not be visible here (set via -D on parent cmake),
# so we also check DUCKDB_EXPLICIT_PLATFORM which is reliably propagated.
set(EMSCRIPTEN_THREAD_CFLAGS "")
if(EMSCRIPTEN AND (CMAKE_CXX_FLAGS MATCHES "-pthread" OR DUCKDB_EXPLICIT_PLATFORM STREQUAL "wasm_threads"))
    set(EMSCRIPTEN_THREAD_CFLAGS "-pthread")
    message(STATUS "Emscripten threading detected: propagating -pthread to C library builds")
endif()

# Find zstd for compression support (optional - build will succeed without it but zstd won't work)
find_package(zstd QUIET)
if(zstd_FOUND)
    # Determine which target is available
    if(TARGET zstd::libzstd_static)
        set(ZSTD_TARGET zstd::libzstd_static)
        set(HAVE_LIBZSTD ON)
        message(STATUS "Found zstd: enabling zstd compression support (target: zstd::libzstd_static)")
    elseif(TARGET zstd::libzstd_shared)
        set(ZSTD_TARGET zstd::libzstd_shared)
        set(HAVE_LIBZSTD ON)
        message(STATUS "Found zstd: enabling zstd compression support (target: zstd::libzstd_shared)")
    elseif(TARGET zstd::libzstd)
        set(ZSTD_TARGET zstd::libzstd)
        set(HAVE_LIBZSTD ON)
        message(STATUS "Found zstd: enabling zstd compression support (target: zstd::libzstd)")
    elseif(TARGET zstd)
        set(ZSTD_TARGET zstd)
        set(HAVE_LIBZSTD ON)
        message(STATUS "Found zstd: enabling zstd compression support (target: zstd)")
    else()
        message(WARNING "zstd found but no known target available - disabling zstd support")
        set(HAVE_LIBZSTD OFF)
    endif()
else()
    message(WARNING "zstd not found - building without zstd compression support")
    set(HAVE_LIBZSTD OFF)
endif()

set(EXTENSION_NAME ${TARGET_NAME}_extension)
set(LOADABLE_EXTENSION_NAME ${TARGET_NAME}_loadable_extension)

project(${TARGET_NAME})

include(ExternalProject)

# Capture git describe for submodules at configure time
execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2
    OUTPUT_VARIABLE MINIMAP2_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT MINIMAP2_GIT_VERSION)
    set(MINIMAP2_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib
    OUTPUT_VARIABLE WFA2_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT WFA2_GIT_VERSION)
    set(WFA2_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype
    OUTPUT_VARIABLE RYPE_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT RYPE_GIT_VERSION)
    set(RYPE_GIT_VERSION "unknown")
endif()

message(STATUS "Submodule versions: minimap2=${MINIMAP2_GIT_VERSION} WFA2=${WFA2_GIT_VERSION} rype=${RYPE_GIT_VERSION}")

# Build htslib from release tarball
# Configure with zstd support if available (htslib autodetects via pkg-config)
if(HAVE_LIBZSTD)
    # Set PKG_CONFIG_PATH so htslib configure can find zstd
    set(HTSLIB_PKG_CONFIG_PATH "PKG_CONFIG_PATH=${CMAKE_CURRENT_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET}/lib/pkgconfig")
else()
    set(HTSLIB_PKG_CONFIG_PATH "")
endif()

# Set CFLAGS and LDFLAGS to find ZLIB from vcpkg
if(ZLIB_LIBRARY_RELEASE)
    get_filename_component(ZLIB_LIBRARY_DIR "${ZLIB_LIBRARY_RELEASE}" DIRECTORY)
    message(STATUS "Using ZLIB from ZLIB_LIBRARY_RELEASE: ${ZLIB_LIBRARY_DIR}")
elseif(ZLIB_LIBRARY)
    get_filename_component(ZLIB_LIBRARY_DIR "${ZLIB_LIBRARY}" DIRECTORY)
    message(STATUS "Using ZLIB from ZLIB_LIBRARY: ${ZLIB_LIBRARY_DIR}")
else()
    # Extract first real library path (skip keywords like "optimized"/"debug")
    set(ZLIB_LIBRARY_DIR "")
    foreach(lib ${ZLIB_LIBRARIES})
        if(EXISTS "${lib}" AND NOT ZLIB_LIBRARY_DIR)
            get_filename_component(ZLIB_LIBRARY_DIR "${lib}" DIRECTORY)
            message(STATUS "Using ZLIB from ZLIB_LIBRARIES: ${ZLIB_LIBRARY_DIR}")
        endif()
    endforeach()
endif()

if(NOT ZLIB_LIBRARY_DIR)
    message(WARNING "Could not determine ZLIB library directory. ZLIB_LIBRARY_RELEASE=${ZLIB_LIBRARY_RELEASE}, ZLIB_LIBRARY=${ZLIB_LIBRARY}, ZLIB_LIBRARIES=${ZLIB_LIBRARIES}")
endif()

# Handle cross-compilation for OSX (e.g., building x86_64 on ARM64 host)
set(HTSLIB_HOST_FLAG "")
set(HTSLIB_ARCH_CFLAGS "")
if(EMSCRIPTEN)
    # Emscripten: tell autotools we are cross-compiling.
    # -fPIC IS required: the extension is linked as a WASM side module (-sSIDE_MODULE=2)
    # which needs position-independent code for all linked objects.
    set(HTSLIB_HOST_FLAG "--host=wasm32-unknown-emscripten")
    message(STATUS "HTSlib cross-compiling for Emscripten/WASM")
elseif(APPLE AND OSX_BUILD_ARCH)
    if(OSX_BUILD_ARCH STREQUAL "x86_64")
        set(HTSLIB_HOST_FLAG "--host=x86_64-apple-darwin")
        set(HTSLIB_ARCH_CFLAGS "-arch x86_64")
        message(STATUS "HTSlib cross-compiling for x86_64 on Apple Silicon")
    elseif(OSX_BUILD_ARCH STREQUAL "arm64")
        set(HTSLIB_HOST_FLAG "--host=aarch64-apple-darwin")
        set(HTSLIB_ARCH_CFLAGS "-arch arm64")
        message(STATUS "HTSlib cross-compiling for arm64")
    endif()
endif()

if(EMSCRIPTEN)
    set(HTSLIB_CFLAGS "CFLAGS=-fPIC ${HTSLIB_ARCH_CFLAGS} ${EMSCRIPTEN_THREAD_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
else()
    set(HTSLIB_CFLAGS "CFLAGS=-fPIC ${HTSLIB_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
endif()
if(ZLIB_LIBRARY_DIR)
    set(HTSLIB_LDFLAGS "LDFLAGS=${HTSLIB_ARCH_CFLAGS} -L${ZLIB_LIBRARY_DIR}")
else()
    set(HTSLIB_LDFLAGS "LDFLAGS=${HTSLIB_ARCH_CFLAGS}")
endif()

# Get the actual zlib library file for explicit linking
# This is needed because HTSlib's configure uses -lz which may not find vcpkg's zlib
if(ZLIB_LIBRARY_RELEASE)
    set(ZLIB_LIBRARY_FILE "${ZLIB_LIBRARY_RELEASE}")
elseif(ZLIB_LIBRARY)
    set(ZLIB_LIBRARY_FILE "${ZLIB_LIBRARY}")
else()
    # Find first actual library file in ZLIB_LIBRARIES
    set(ZLIB_LIBRARY_FILE "")
    foreach(lib ${ZLIB_LIBRARIES})
        if(EXISTS "${lib}" AND NOT ZLIB_LIBRARY_FILE)
            set(ZLIB_LIBRARY_FILE "${lib}")
        endif()
    endforeach()
endif()

# Set LIBS to pass explicit library path to configure (needed for OSX with vcpkg)
if(ZLIB_LIBRARY_FILE)
    set(HTSLIB_LIBS "LIBS=${ZLIB_LIBRARY_FILE}")
    message(STATUS "HTSlib will use explicit ZLIB library: ${ZLIB_LIBRARY_FILE}")
else()
    set(HTSLIB_LIBS "")
    message(WARNING "Could not determine explicit ZLIB library file for HTSlib")
endif()

if(WIN32)
    # Windows/MinGW: skip autotools configure (requires full MSYS2 shell which is
    # not available in DuckDB's rtools42 CI environment). Instead, generate config.h
    # and config.mk manually to match our --disable-libcurl --disable-bz2 etc. flags.
    # config.mk sets NONCONFIGURE_OBJS= to prevent hfile_libcurl.o from being built.
    # config.h is generated with only the features we need (zlib only).
    # See: https://github.com/samtools/htslib/blob/develop/.github/workflows/windows-build.yml
    set(HTSLIB_SRC ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1)
    set(HTSLIB_CONFIG_H "${HTSLIB_SRC}/config.h")
    set(HTSLIB_CONFIG_MK "${HTSLIB_SRC}/config.mk")
    # Generate config.h: minimal — zlib only, no curl/bz2/lzma/libdeflate
    file(WRITE ${HTSLIB_CONFIG_H}
        "/* config.h generated by CMake for Windows/MinGW (no autotools) */\n"
        "#define HAVE_FSEEKO 1\n"
    )
    # Generate config.mk: blank NONCONFIGURE_OBJS to prevent hfile_libcurl.o,
    # set LIBS to just -lz, and disable features we don't need.
    file(WRITE ${HTSLIB_CONFIG_MK}
        "# config.mk generated by CMake for Windows/MinGW (no autotools)\n"
        "NONCONFIGURE_OBJS =\n"
        "CFLAGS = -g -Wall -O2 -fvisibility=hidden -I${ZLIB_INCLUDE_DIRS}\n"
        "LDFLAGS = -fvisibility=hidden\n"
        "LIBS = -lz\n"
        "prefix = ${CMAKE_CURRENT_BINARY_DIR}/htslib\n"
    )
    ExternalProject_Add(
        htslib_build
        SOURCE_DIR ${HTSLIB_SRC}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -C ${HTSLIB_SRC} clean libhts.a
        INSTALL_COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib
            COMMAND ${CMAKE_COMMAND} -E copy ${HTSLIB_SRC}/libhts.a ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/htslib/include/htslib
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${HTSLIB_SRC}/htslib ${CMAKE_CURRENT_BINARY_DIR}/htslib/include/htslib
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
    )
else()
    if(EMSCRIPTEN)
        # Emscripten/WASM doesn't support shared libraries, so htslib's default
        # 'make' (which builds both static and shared) and 'make install' (which
        # requires the shared lib) both fail. Build only the static library and
        # manually copy the artifacts.
        set(HTSLIB_BUILD_COMMAND make lib-static)
        # libhts.a is in <BINARY_DIR> (out-of-source build), but headers are
        # read-only source files so we copy them from the source tree directly.
        set(HTSLIB_INSTALL_COMMAND
            ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib
            COMMAND ${CMAKE_COMMAND} -E copy <BINARY_DIR>/libhts.a ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
            COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_CURRENT_BINARY_DIR}/htslib/include/htslib
            COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1/htslib ${CMAKE_CURRENT_BINARY_DIR}/htslib/include/htslib)
    else()
        set(HTSLIB_BUILD_COMMAND make)
        set(HTSLIB_INSTALL_COMMAND make install)
    endif()
    ExternalProject_Add(
        htslib_build
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1
        CONFIGURE_COMMAND env ${HTSLIB_CFLAGS} ${HTSLIB_LDFLAGS} ${HTSLIB_LIBS} ${HTSLIB_PKG_CONFIG_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1/configure
            ${HTSLIB_HOST_FLAG}
            --prefix=${CMAKE_CURRENT_BINARY_DIR}/htslib
            --disable-libcurl
            --disable-bz2
            --disable-lzma
            --without-libdeflate
        BUILD_COMMAND ${HTSLIB_BUILD_COMMAND}
        INSTALL_COMMAND ${HTSLIB_INSTALL_COMMAND}
        BUILD_IN_SOURCE 0
        BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
    )
endif()

# Create imported library
add_library(hts STATIC IMPORTED GLOBAL)
set_target_properties(hts PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
)
set_property(TARGET hts PROPERTY
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/htslib/include)
add_dependencies(hts htslib_build)

# Detect architecture for minimap2 SIMD support
# Use OSX_BUILD_ARCH for cross-compilation, otherwise use host architecture
set(MINIMAP2_ARCH_FLAGS "")
set(MINIMAP2_ARCH_CFLAGS "")
set(MINIMAP2_FPIC "-fPIC")
if(EMSCRIPTEN)
    # Emscripten: use WASM SIMD via SSE2 emulation (minimap2 has SIMDe support).
    # sse2only=1 avoids SSE4.1 dispatch which requires runtime CPU detection.
    set(MINIMAP2_ARCH_FLAGS "sse2only=1")
    set(MINIMAP2_ARCH_CFLAGS "-msimd128 -msse -msse2 ${EMSCRIPTEN_THREAD_CFLAGS}")
    # Keep -fPIC: required for WASM side module linking (-sSIDE_MODULE=2)
    message(STATUS "Minimap2 building for Emscripten/WASM with SIMD via SSE2 emulation")
elseif(APPLE AND OSX_BUILD_ARCH)
    # Cross-compilation: use target architecture, not host
    if(OSX_BUILD_ARCH STREQUAL "arm64")
        set(MINIMAP2_ARCH_FLAGS "aarch64=1")
        set(MINIMAP2_ARCH_CFLAGS "-arch arm64")
        message(STATUS "Minimap2 cross-compiling for arm64")
    elseif(OSX_BUILD_ARCH STREQUAL "x86_64")
        # x86_64 target: don't set aarch64=1, use x86_64 arch flag
        set(MINIMAP2_ARCH_CFLAGS "-arch x86_64")
        message(STATUS "Minimap2 cross-compiling for x86_64 on Apple Silicon")
    endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
    # Native ARM64 build
    set(MINIMAP2_ARCH_FLAGS "aarch64=1")
endif()

# Detect DuckDB's jemalloc (not present on musl/macOS/Windows builds).
# DUCKDB_EXTENSION_JEMALLOC_SHOULD_LINK is set during extension registration,
# before add_subdirectory, so it's available here.
if(DUCKDB_EXTENSION_JEMALLOC_SHOULD_LINK)
    set(MIINT_HAS_JEMALLOC TRUE)
    set(MINIMAP2_ALLOC_FLAGS "-DMIINT_USE_JEMALLOC -include ${CMAKE_CURRENT_SOURCE_DIR}/src/include/mm_alloc.h")
    message(STATUS "DuckDB jemalloc detected — minimap2 will use jemalloc allocator")
else()
    set(MIINT_HAS_JEMALLOC FALSE)
    set(MINIMAP2_ALLOC_FLAGS "")
    message(STATUS "DuckDB jemalloc not detected — minimap2 will use system allocator")
endif()

# Build minimap2 from source (v2.30)
ExternalProject_Add(
    minimap2_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 clean libminimap2.a "CFLAGS=-g -Wall -O2 -Wc++-compat ${MINIMAP2_FPIC} ${MINIMAP2_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}" "CPPFLAGS=-DHAVE_KALLOC ${MINIMAP2_ALLOC_FLAGS}" ${MINIMAP2_ARCH_FLAGS}
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2/libminimap2.a
        ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a
)

# Create imported library for minimap2
add_library(minimap2 STATIC IMPORTED GLOBAL)
set_target_properties(minimap2 PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a
)
add_dependencies(minimap2 minimap2_build)

# Handle cross-compilation arch flags for WFA2
set(WFA2_ARCH_CFLAGS "")
set(WFA2_FPIC "-fPIC")
if(EMSCRIPTEN)
    # Keep -fPIC: required for WASM side module linking (-sSIDE_MODULE=2)
    set(WFA2_ARCH_CFLAGS "${EMSCRIPTEN_THREAD_CFLAGS}")
elseif(APPLE AND OSX_BUILD_ARCH)
    if(OSX_BUILD_ARCH STREQUAL "x86_64")
        set(WFA2_ARCH_CFLAGS "-arch x86_64")
    elseif(OSX_BUILD_ARCH STREQUAL "arm64")
        set(WFA2_ARCH_CFLAGS "-arch arm64")
    endif()
endif()

# Build WFA2-lib from source (v2.3.5) using Makefile (matching minimap2 pattern)
ExternalProject_Add(
    wfa2_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib clean all
        "BUILD_TOOLS=0" "BUILD_EXAMPLES=0" "BUILD_WFA_PARALLEL=0"
        "CC_FLAGS=-Wall -O2 ${WFA2_FPIC} ${WFA2_ARCH_CFLAGS}"
    INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfa.a
        ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfa.a
        COMMAND ${CMAKE_COMMAND} -E copy
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfacpp.a
        ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfacpp.a
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS
        ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfa.a
        ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfacpp.a
)

# Create imported libraries for WFA2 (C++ bindings depend on C core — link order matters)
add_library(wfa2 STATIC IMPORTED GLOBAL)
set_target_properties(wfa2 PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfa.a
)
add_dependencies(wfa2 wfa2_build)

add_library(wfa2cpp STATIC IMPORTED GLOBAL)
set_target_properties(wfa2cpp PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfacpp.a
)
add_dependencies(wfa2cpp wfa2_build)

# vsearch: autotools-based build, not supported on Windows or WASM.
if(WIN32 OR EMSCRIPTEN)
    set(MIINT_ENABLE_VSEARCH OFF)
else()
    set(MIINT_ENABLE_VSEARCH ON)
endif()

if(MIINT_ENABLE_VSEARCH)
    # Build vsearch from source (v2.30.5-miint fork) using autotools.
    # Produces libvsearch.a — a PIC static archive for embedding vsearch's
    # search, clustering, chimera detection, DUST masking, and merge-pairs
    # algorithms into the extension.
    set(VSEARCH_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/vsearch)

    # Handle cross-compilation for vsearch (same pattern as htslib)
    set(VSEARCH_HOST_FLAG "")
    set(VSEARCH_ARCH_CFLAGS "")
    if(APPLE AND OSX_BUILD_ARCH)
        if(OSX_BUILD_ARCH STREQUAL "x86_64")
            set(VSEARCH_HOST_FLAG "--host=x86_64-apple-darwin")
            set(VSEARCH_ARCH_CFLAGS "-arch x86_64")
        elseif(OSX_BUILD_ARCH STREQUAL "arm64")
            set(VSEARCH_HOST_FLAG "--host=aarch64-apple-darwin")
            set(VSEARCH_ARCH_CFLAGS "-arch arm64")
        endif()
    endif()

    ExternalProject_Add(
        vsearch_build
        SOURCE_DIR ${VSEARCH_SRC_DIR}
        # Touch pre-generated autotools files so make doesn't try to regenerate
        # them (git checkout gives all files the same timestamp, which can trigger
        # autotools rebuild rules if make thinks sources are newer).
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E touch
            ${VSEARCH_SRC_DIR}/aclocal.m4
            ${VSEARCH_SRC_DIR}/configure
            ${VSEARCH_SRC_DIR}/Makefile.in
            ${VSEARCH_SRC_DIR}/src/Makefile.in
            ${VSEARCH_SRC_DIR}/man/Makefile.in
            ${VSEARCH_SRC_DIR}/config.h.in
            COMMAND ${CMAKE_COMMAND} -E env "CXXFLAGS=-fPIC ${VSEARCH_ARCH_CFLAGS}" "CFLAGS=-fPIC ${VSEARCH_ARCH_CFLAGS}"
            ${VSEARCH_SRC_DIR}/configure ${VSEARCH_HOST_FLAG}
        BUILD_COMMAND make -C src libvsearch.a
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS ${VSEARCH_SRC_DIR}/src/libvsearch.a
    )

    # Create imported library for vsearch
    add_library(vsearch STATIC IMPORTED GLOBAL)
    set_target_properties(vsearch PROPERTIES
        IMPORTED_LOCATION ${VSEARCH_SRC_DIR}/src/libvsearch.a
    )
    add_dependencies(vsearch vsearch_build)
endif()

# SortMeRNA: embedded library for rRNA filtering and alignment.
# Upstream CMake produces several archives; smr_api is the C-API entry point.
# smr_objs, cmph, alp, build_version are CMake OBJECT libraries whose object
# files get archived into libsmr_api.a via the STATIC→OBJECT link rule.
#
# Like the HDF5 block above: clang-tidy builds don't have vcpkg wired up, so
# RocksDB (the only hard external dependency) may be absent. In that case we
# silently drop sortmerna from the analysis set rather than fail the whole
# tidy pipeline. Regular builds still REQUIRE RocksDB.
if(MIINT_ENABLE_SORTMERNA)
    if(CLANG_TIDY)
        find_package(RocksDB CONFIG QUIET)
        if(NOT RocksDB_FOUND)
            message(STATUS "RocksDB not found for tidy-check - SortMeRNA support will be excluded from analysis")
            set(MIINT_ENABLE_SORTMERNA OFF)
        else()
            message(STATUS "RocksDB found for tidy-check - SortMeRNA support enabled")
        endif()
    else()
        find_package(RocksDB CONFIG REQUIRED)
    endif()
endif()

if(MIINT_ENABLE_SORTMERNA)
    set(SORTMERNA_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sortmerna)
    set(SORTMERNA_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/sortmerna_build)
    set(SORTMERNA_API_LIB ${SORTMERNA_BIN_DIR}/src/smr_api/libsmr_api.a)
    set(SORTMERNA_BUNDLE_LIB ${CMAKE_CURRENT_BINARY_DIR}/sortmerna/libsortmerna_bundle.a)
    set(CONCURRENTQUEUE_HOME ${CMAKE_CURRENT_SOURCE_DIR}/ext/concurrentqueue)

    # macOS cross-compile (GitHub runners are arm64; osx_amd64 builds still need
    # to emit x86_64 objects). Without forwarding CMAKE_OSX_ARCHITECTURES into
    # the sub-CMake, sortmerna's CLI target gets linked for the host arch and
    # trips on vcpkg's RocksDB which is built for the target arch.
    set(SORTMERNA_OSX_ARCH_ARG "")
    if(APPLE AND OSX_BUILD_ARCH)
        set(SORTMERNA_OSX_ARCH_ARG "-DCMAKE_OSX_ARCHITECTURES=${OSX_BUILD_ARCH}")
    endif()

    # alp's sources use the deprecated `register` keyword, which is a hard
    # error on Clang 16+. Suppress it so the library still compiles.
    ExternalProject_Add(
        sortmerna_build
        SOURCE_DIR ${SORTMERNA_SRC_DIR}
        BINARY_DIR ${SORTMERNA_BIN_DIR}
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DCMAKE_CXX_FLAGS=-Wno-register
            -DWITH_TESTS=OFF
            -DCONCURRENTQUEUE_HOME=${CONCURRENTQUEUE_HOME}
            -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
            -DVCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}
            -DVCPKG_MANIFEST_MODE=OFF
            ${SORTMERNA_OSX_ARCH_ARG}
            # Force sortmerna's FindRocksDB branch to use vcpkg's CONFIG file.
            # Without this, its else-branch runs MODULE mode which hits system
            # /usr/lib rocksdb first on Debian/Ubuntu.
            -DROCKSDB_DIST=${CMAKE_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET}/share/rocksdb
            # vcpkg builds static-only; switch sortmerna off its default shared target.
            -DROCKSDB_USE_STATIC_LIBS=ON
        BUILD_BYPRODUCTS ${SORTMERNA_API_LIB}
        INSTALL_COMMAND ""
    )

    # Post-build: bundle cmph/alp/build_version OBJECT library outputs into a
    # single archive so consumers only need to link one -lsortmerna equivalent.
    add_custom_command(
        OUTPUT ${SORTMERNA_BUNDLE_LIB}
        COMMAND ${CMAKE_COMMAND}
            -DSORTMERNA_BIN_DIR=${SORTMERNA_BIN_DIR}
            -DSORTMERNA_API_LIB=${SORTMERNA_API_LIB}
            -DSORTMERNA_BUNDLE_LIB=${SORTMERNA_BUNDLE_LIB}
            -DCMAKE_AR=${CMAKE_AR}
            -P ${CMAKE_CURRENT_SOURCE_DIR}/cmake/bundle_sortmerna.cmake
        DEPENDS sortmerna_build
        COMMENT "Bundling sortmerna object libraries into ${SORTMERNA_BUNDLE_LIB}"
    )
    add_custom_target(sortmerna_bundle DEPENDS ${SORTMERNA_BUNDLE_LIB})

    add_library(sortmerna STATIC IMPORTED GLOBAL)
    set_target_properties(sortmerna PROPERTIES
        IMPORTED_LOCATION ${SORTMERNA_BUNDLE_LIB}
        INTERFACE_INCLUDE_DIRECTORIES ${SORTMERNA_SRC_DIR}/include
    )
    add_dependencies(sortmerna sortmerna_bundle)
endif()

# Build MAFFT PartTree static library for multiple sequence alignment
if(MIINT_ENABLE_MAFFT)
    # Handle cross-compilation arch flags for MAFFT (same pattern as WFA2/minimap2)
    set(MAFFT_ARCH_CFLAGS "")
    set(MAFFT_FPIC "-fPIC")
    if(EMSCRIPTEN)
        # Keep -fPIC: required for WASM side module linking (-sSIDE_MODULE=2)
        set(MAFFT_ARCH_CFLAGS "${EMSCRIPTEN_THREAD_CFLAGS}")
    elseif(APPLE AND OSX_BUILD_ARCH)
        if(OSX_BUILD_ARCH STREQUAL "x86_64")
            set(MAFFT_ARCH_CFLAGS "-arch x86_64")
        elseif(OSX_BUILD_ARCH STREQUAL "arm64")
            set(MAFFT_ARCH_CFLAGS "-arch arm64")
        endif()
    endif()

    ExternalProject_Add(
        mafft_build
        SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core
        CONFIGURE_COMMAND ""
        BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core clean libmafft_parttree.a
            "CC=${CMAKE_C_COMPILER}"
            "CFLAGS=-O3 ${MAFFT_FPIC} ${MAFFT_ARCH_CFLAGS}"
            "ENABLE_MULTITHREAD=-Denablemultithread"
        INSTALL_COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core/libmafft_parttree.a
            ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a
    )

    add_library(mafft_parttree STATIC IMPORTED GLOBAL)
    set_target_properties(mafft_parttree PROPERTIES
        IMPORTED_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a
    )
    add_dependencies(mafft_parttree mafft_build)
endif()

# Build RYpe from Rust source with Arrow support
find_program(CARGO cargo REQUIRED)
message(STATUS "Found cargo: ${CARGO}")

# Determine the static library name based on platform.
# On Windows/MinGW, cargo with --target x86_64-pc-windows-gnu produces librype.a (not .lib).
set(RYPE_LIB_NAME "librype.a")

set(RYPE_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/rype_target)

# Handle cross-compilation for Rust via --target flag.
# When cross-compiling, cargo places output in <target-dir>/<triple>/release/
# instead of <target-dir>/release/.
# Also ensure the Rust target is installed via `rustup target add`.
set(RYPE_CARGO_TARGET_FLAG "")
set(RYPE_LIB_SUBDIR "release")
set(RYPE_CONFIGURE_CMD "")
if(EMSCRIPTEN)
    find_program(RUSTUP rustup)
    set(RYPE_CARGO_TARGET_FLAG "--target" "wasm32-unknown-emscripten")
    set(RYPE_LIB_SUBDIR "wasm32-unknown-emscripten/release")
    if(RUSTUP)
        set(RYPE_CONFIGURE_CMD ${RUSTUP} target add wasm32-unknown-emscripten)
    endif()
    message(STATUS "RYpe cross-compiling for wasm32-unknown-emscripten")
elseif(WIN32)
    # Windows/MinGW: Rust defaults to MSVC target but we link with MinGW ld.exe.
    # Must target x86_64-pc-windows-gnu to produce compatible object files.
    find_program(RUSTUP rustup)
    set(RYPE_CARGO_TARGET_FLAG "--target" "x86_64-pc-windows-gnu")
    set(RYPE_LIB_SUBDIR "x86_64-pc-windows-gnu/release")
    if(RUSTUP)
        set(RYPE_CONFIGURE_CMD ${RUSTUP} target add x86_64-pc-windows-gnu)
    endif()
    message(STATUS "RYpe targeting x86_64-pc-windows-gnu for MinGW compatibility")
elseif(APPLE AND OSX_BUILD_ARCH)
    find_program(RUSTUP rustup)
    if(OSX_BUILD_ARCH STREQUAL "x86_64")
        set(RYPE_CARGO_TARGET_FLAG "--target" "x86_64-apple-darwin")
        set(RYPE_LIB_SUBDIR "x86_64-apple-darwin/release")
        if(RUSTUP)
            set(RYPE_CONFIGURE_CMD ${RUSTUP} target add x86_64-apple-darwin)
        endif()
        message(STATUS "RYpe cross-compiling for x86_64 on Apple Silicon")
    elseif(OSX_BUILD_ARCH STREQUAL "arm64")
        set(RYPE_CARGO_TARGET_FLAG "--target" "aarch64-apple-darwin")
        set(RYPE_LIB_SUBDIR "aarch64-apple-darwin/release")
        if(RUSTUP)
            set(RYPE_CONFIGURE_CMD ${RUSTUP} target add aarch64-apple-darwin)
        endif()
        message(STATUS "RYpe cross-compiling for arm64")
    endif()
endif()

set(RYPE_LIB_PATH ${RYPE_TARGET_DIR}/${RYPE_LIB_SUBDIR}/${RYPE_LIB_NAME})

# Build RYpe with Arrow FFI support for zero-copy data exchange.
# On Emscripten and Windows: --no-default-features disables the "fastx" feature
# (needletail), whose cdylib crate-type triggers linker errors (emcc standalone
# linking on WASM, unresolved ___chkstk_ms on MSVC).
if(EMSCRIPTEN OR WIN32)
    set(RYPE_CARGO_FEATURES "--no-default-features" "--features" "arrow-ffi")
else()
    set(RYPE_CARGO_FEATURES "--features" "arrow-ffi")
endif()
# For WASM side modules (-sSIDE_MODULE=2), all objects need PIC relocations.
# Rust's default for wasm32-unknown-emscripten is static; override to PIC.
if(EMSCRIPTEN)
    set(RYPE_BUILD_ENV ${CMAKE_COMMAND} -E env "RUSTFLAGS=-C relocation-model=pic")
else()
    # Unset (not empty string "") so it expands to zero arguments in BUILD_COMMAND
    set(RYPE_BUILD_ENV)
endif()
ExternalProject_Add(
    rype_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype
    CONFIGURE_COMMAND "${RYPE_CONFIGURE_CMD}"
    BUILD_COMMAND ${RYPE_BUILD_ENV} ${CARGO} rustc --release ${RYPE_CARGO_FEATURES} --lib --crate-type=staticlib
        --target-dir ${RYPE_TARGET_DIR}
        --manifest-path ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype/Cargo.toml
        ${RYPE_CARGO_TARGET_FLAG}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 0
    BUILD_BYPRODUCTS ${RYPE_LIB_PATH}
)

# Create imported library for rype
add_library(rype STATIC IMPORTED GLOBAL)
set_target_properties(rype PROPERTIES
    IMPORTED_LOCATION ${RYPE_LIB_PATH}
)
add_dependencies(rype rype_build)

# adapted kseq++ include use from strobealign
include_directories(src/include ext ext/WFA2-lib ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype)
if(MIINT_ENABLE_MAFFT)
    include_directories(ext/mafft/core)
endif()
if(MIINT_ENABLE_VSEARCH)
    include_directories(${VSEARCH_SRC_DIR}/src)
endif()
set(EXTENSION_SOURCES
    src/SequenceReader.cpp
    src/QualScore.cpp
    src/read_fastx.cpp
    src/miint_extension.cpp
    src/read_alignments.cpp
    src/read_sequences_sam.cpp
    src/SAMReader.cpp
    src/alignment_flag_functions.cpp
    src/alignment_functions.cpp
    src/sequence_functions.cpp
    src/IntervalCompressor.cpp
    src/compress_intervals.cpp
    src/CoverageDepthCalculator.cpp
    src/compute_coverage_depth.cpp
    src/AlignmentSlicer.cpp
    src/alignment_slice.cpp
    src/copy_format_common.cpp
    src/table_function_common.cpp
    src/copy_fastq.cpp
    src/copy_fasta.cpp
    src/copy_sam.cpp
    src/catalog_utils.cpp
    src/reference_table_reader.cpp
    src/placement_table_reader.cpp
    src/NewickTree.cpp
    src/read_jplace_newick.cpp
    src/read_newick.cpp
    src/tree_table_reader.cpp
    src/tree_resolve_placement.cpp
    src/copy_newick.cpp
    src/Minimap2Aligner.cpp
    src/sequence_table_reader.cpp
    src/sequence_data_reader.cpp
    src/align_minimap2.cpp
    src/align_minimap2_sharded.cpp
    src/save_minimap2_index.cpp
    src/ncbi_parser.cpp
    src/ncbi_client.cpp
    src/zip_utils.cpp
    src/read_ncbi_fasta.cpp
    src/read_ncbi.cpp
    src/read_ncbi_annotation.cpp
    src/ena_parser.cpp
    src/ena_client.cpp
    src/read_ena.cpp
    src/read_ena_attributes.cpp
    src/read_ena_sequences.cpp
    src/aspera_utils.cpp
    src/aspera_stream.cpp
    src/SFFReader.cpp
    src/read_sequences_sff.cpp
    src/WFA2Aligner.cpp
    src/align_pairwise_functions.cpp
    src/rype_classify.cpp
    src/rype_extract.cpp
    src/MsBinaryDecoder.cpp
    src/MzMLReader.cpp
    src/read_mzml.cpp
    src/read_mzml_chromatograms.cpp
    src/rype_log_ratio.cpp
    src/formula_parser.cpp
    src/formula_function.cpp
    src/massql_parser.cpp
    src/massql_transpiler.cpp
    src/massql_function.cpp
    src/mzml_peak_pair_function.cpp
    src/mass_tables.cpp
    src/MzXMLReader.cpp
    src/read_mzxml.cpp
    src/remote_file_helper.cpp
    src/duckdb_seq_stream.cpp
    src/hfile_duckdb.cpp
    src/deblur.cpp
    src/deblur_table_function.cpp)

if(MIINT_ENABLE_VSEARCH)
    list(APPEND EXTENSION_SOURCES
        src/VsearchChimeraWrapper.cpp src/VsearchSearchWrapper.cpp src/VsearchClusterWrapper.cpp
        src/mask_dust_function.cpp src/merge_pairs_function.cpp src/search_sequences.cpp
        src/cluster_sequences.cpp src/uchime_common.cpp src/uchime_ref.cpp src/uchime_denovo.cpp)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    list(APPEND EXTENSION_SOURCES src/copy_biom.cpp src/BIOMTable.cpp src/BIOMReader.cpp src/read_biom.cpp)
endif()
if(MIINT_ENABLE_BOWTIE2)
    list(APPEND EXTENSION_SOURCES src/Bowtie2Aligner.cpp src/align_bowtie2.cpp src/align_bowtie2_sharded.cpp)
endif()
if(MIINT_ENABLE_MAFFT)
    list(APPEND EXTENSION_SOURCES src/MafftAligner.cpp src/align_mafft.cpp)
endif()
if(MIINT_ENABLE_SORTMERNA)
    list(APPEND EXTENSION_SOURCES src/SortMeRNAAligner.cpp src/align_sortmerna.cpp src/align_sortmerna_rrna.cpp)
endif()

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})

# WASM builds: DuckDB's build_loadable_extension runs a POST_BUILD emcc command
# to create the .wasm side module. It reads additional libraries to link from
# DUCKDB_EXTENSION_<NAME>_LINKED_LIBS — NOT from target_link_libraries (which
# only records metadata for STATIC library targets on Emscripten). Without this,
# symbols like gzclose become unresolved imports and fail at extension load time
# with "bad export type: undefined".
if(EMSCRIPTEN)
    # Space-separated string of library paths. extension_build_tools.cmake runs
    # separate_arguments() then passes the result to add_custom_command COMMAND,
    # where generator expressions ($<TARGET_FILE:...>) are evaluated at build time.
    set(WASM_LIBS "")
    string(APPEND WASM_LIBS "$<TARGET_FILE:ZLIB::ZLIB> ")
    string(APPEND WASM_LIBS "$<TARGET_FILE:expat::expat> ")
    string(APPEND WASM_LIBS "${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a ")
    string(APPEND WASM_LIBS "${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a ")
    string(APPEND WASM_LIBS "${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfacpp.a ")
    string(APPEND WASM_LIBS "${CMAKE_CURRENT_BINARY_DIR}/wfa2/lib/libwfa.a ")
    string(APPEND WASM_LIBS "${RYPE_LIB_PATH}")
    if(HAVE_LIBZSTD)
        string(APPEND WASM_LIBS " $<TARGET_FILE:${ZSTD_TARGET}>")
    endif()
    if(MIINT_ENABLE_MAFFT)
        string(APPEND WASM_LIBS " ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a")
    endif()
    set(DUCKDB_EXTENSION_MIINT_LINKED_LIBS "${WASM_LIBS}")
    message(STATUS "WASM linked libs: ${DUCKDB_EXTENSION_MIINT_LINKED_LIBS}")
endif()

build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})

# Enable compiler warnings for extension code only (not external libraries)
#if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#    target_compile_options(${EXTENSION_NAME} PRIVATE -Wall -Wextra)
#    target_compile_options(${LOADABLE_EXTENSION_NAME} PRIVATE -Wall -Wextra)
#elseif(MSVC)
#    target_compile_options(${EXTENSION_NAME} PRIVATE /W4)
#    target_compile_options(${LOADABLE_EXTENSION_NAME} PRIVATE /W4)
#endif()

target_compile_definitions(${EXTENSION_NAME} PRIVATE RYPE_ARROW
    MIINT_STATIC_BUILD
    MINIMAP2_GIT_VERSION="${MINIMAP2_GIT_VERSION}"
    WFA2_GIT_VERSION="${WFA2_GIT_VERSION}"
    RYPE_GIT_VERSION="${RYPE_GIT_VERSION}"
    $<$<BOOL:${MIINT_HAS_JEMALLOC}>:MIINT_USE_JEMALLOC>)
if(MIINT_ENABLE_VSEARCH)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_VSEARCH)
endif()
if(MIINT_ENABLE_BOWTIE2)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_BOWTIE2)
endif()
if(MIINT_ENABLE_MAFFT)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_MAFFT)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_SORTMERNA)
endif()
target_link_libraries(${EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts
                      minimap2
                      wfa2cpp
                      wfa2
                      rype
                      expat::expat
)
if(MIINT_ENABLE_MAFFT)
    target_link_libraries(${EXTENSION_NAME} mafft_parttree)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(${EXTENSION_NAME} vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(${EXTENSION_NAME} sortmerna RocksDB::rocksdb)
endif()
if(APPLE)
    target_link_libraries(${EXTENSION_NAME} "-framework CoreFoundation")
elseif(WIN32)
    # HTSlib: -lregex for POSIX regex (regcomp/regexec/regfree).
    # Rust x86_64-pc-windows-gnu stdlib system deps (from rustc_target spec and
    # std #[link] attributes): ntdll, userenv, bcrypt, advapi32, ws2_32.
    # MinGW driver auto-adds: kernel32, user32, msvcrt, mingwex, mingw32, gcc.
    target_link_libraries(${EXTENSION_NAME}
        regex ws2_32 ntdll userenv bcrypt advapi32)
elseif(UNIX)
    target_link_libraries(${EXTENSION_NAME} dl m)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_link_libraries(${EXTENSION_NAME}
                          hdf5::hdf5-static
                          hdf5::hdf5_hl-static
                          hdf5::hdf5_cpp-static
                          hdf5::hdf5_hl_cpp-static
    )
    target_link_libraries(${LOADABLE_EXTENSION_NAME}
                          hdf5::hdf5-static
                          hdf5::hdf5_hl-static
                          hdf5::hdf5_cpp-static
                          hdf5::hdf5_hl_cpp-static
    )
endif()
if(HAVE_LIBZSTD)
    target_link_libraries(${EXTENSION_NAME} ${ZSTD_TARGET})
endif()

target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE RYPE_ARROW
    MIINT_STATIC_BUILD
    MINIMAP2_GIT_VERSION="${MINIMAP2_GIT_VERSION}"
    WFA2_GIT_VERSION="${WFA2_GIT_VERSION}"
    RYPE_GIT_VERSION="${RYPE_GIT_VERSION}"
    $<$<BOOL:${MIINT_HAS_JEMALLOC}>:MIINT_USE_JEMALLOC>)
if(MIINT_ENABLE_VSEARCH)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_VSEARCH)
endif()
if(MIINT_ENABLE_BOWTIE2)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_BOWTIE2)
endif()
if(MIINT_ENABLE_MAFFT)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_MAFFT)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_SORTMERNA)
endif()
target_link_libraries(${LOADABLE_EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts
                      minimap2
                      wfa2cpp
                      wfa2
                      rype
                      expat::expat
)
if(MIINT_ENABLE_MAFFT)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} mafft_parttree)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} sortmerna RocksDB::rocksdb)
endif()
if(APPLE)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} "-framework CoreFoundation")
elseif(WIN32)
    target_link_libraries(${LOADABLE_EXTENSION_NAME}
        regex ws2_32 ntdll userenv bcrypt advapi32)
elseif(UNIX)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} dl m)
endif()
if(HAVE_LIBZSTD)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} ${ZSTD_TARGET})
endif()
# The loadable extension is dlopen'd at runtime. When the host DuckDB is loaded
# via Python (RTLD_LOCAL), duckdb_je_* symbols from the host are not visible.
# Link stubs that forward to system malloc so the .so is self-contained.
if(MIINT_HAS_JEMALLOC)
    target_sources(${LOADABLE_EXTENSION_NAME} PRIVATE src/mm_alloc_stubs.c)
endif()

install(
  TARGETS ${EXTENSION_NAME}
  EXPORT "${DUCKDB_EXPORT_SET}"
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}"
  ARCHIVE DESTINATION "${INSTALL_LIB_DIR}")

FetchContent_Declare(
  catch2
  GIT_REPOSITORY https://github.com/catchorg/Catch2.git
  GIT_TAG v3.4.0
)
FetchContent_MakeAvailable(catch2)
set(TEST_SOURCES
    src/SequenceReader.cpp
    test/cpp/test_SequenceReader.cpp
    src/QualScore.cpp
    test/cpp/test_QualScore.cpp
    src/SAMReader.cpp
    test/cpp/test_SAMReader.cpp
    test/cpp/test_AlignmentFunctions.cpp
    src/IntervalCompressor.cpp
    test/cpp/test_IntervalCompressor.cpp
    src/CoverageDepthCalculator.cpp
    test/cpp/test_CoverageDepthCalculator.cpp
    src/AlignmentSlicer.cpp
    test/cpp/test_AlignmentSlicer.cpp
    src/NewickTree.cpp
    test/cpp/test_NewickParser.cpp
    test/cpp/test_InsertFullyResolved.cpp
    src/Minimap2Aligner.cpp
    test/cpp/test_Minimap2Aligner.cpp
    src/ncbi_parser.cpp
    src/zip_utils.cpp
    duckdb/third_party/miniz/miniz.cpp  # For ExtractFromZip tests
    test/cpp/test_ncbi_client.cpp
    src/ena_parser.cpp
    test/cpp/test_ENAClient.cpp
    src/aspera_utils.cpp
    test/cpp/test_AsperaUtils.cpp
    src/SFFReader.cpp
    test/cpp/test_SFFReader.cpp
    src/WFA2Aligner.cpp
    test/cpp/test_WFA2Aligner.cpp
    src/MsBinaryDecoder.cpp
    test/cpp/test_MsBinaryDecoder.cpp
    src/MzMLReader.cpp
    test/cpp/test_MzMLReader.cpp
    src/MzXMLReader.cpp
    test/cpp/test_MzXMLReader.cpp
    src/formula_parser.cpp
    src/massql_parser.cpp
    src/massql_transpiler.cpp
    src/mass_tables.cpp
    test/cpp/test_MassQLParser.cpp
    src/remote_file_helper.cpp
    test/cpp/test_RemoteFileHelper.cpp
    test/cpp/test_hfile_duckdb.cpp
    src/deblur.cpp
    test/cpp/test_deblur.cpp
)
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    list(APPEND TEST_SOURCES src/BIOMTable.cpp src/BIOMReader.cpp
         test/cpp/test_BIOMReader.cpp test/cpp/test_BIOMTable.cpp)
endif()
if(MIINT_ENABLE_VSEARCH)
    list(APPEND TEST_SOURCES
        src/VsearchChimeraWrapper.cpp test/cpp/test_VsearchChimeraWrapper.cpp
        src/VsearchSearchWrapper.cpp test/cpp/test_VsearchSearchWrapper.cpp
        src/VsearchClusterWrapper.cpp test/cpp/test_VsearchClusterWrapper.cpp)
endif()
if(MIINT_ENABLE_BOWTIE2)
    list(APPEND TEST_SOURCES src/Bowtie2Aligner.cpp test/cpp/test_Bowtie2Aligner.cpp)
endif()
if(MIINT_ENABLE_MAFFT)
    list(APPEND TEST_SOURCES src/MafftAligner.cpp test/cpp/test_MafftAligner.cpp)
endif()
if(MIINT_ENABLE_SORTMERNA)
    list(APPEND TEST_SOURCES
        src/SortMeRNAAligner.cpp
        test/cpp/test_sortmerna_build.cpp
        test/cpp/test_SortMeRNAAligner.cpp)
endif()

add_executable(tests ${TEST_SOURCES})

target_link_libraries(tests
  Catch2::Catch2WithMain
  Threads::Threads
  ZLIB::ZLIB
  hts
  minimap2
  wfa2cpp
  wfa2
  expat::expat
)
if(MIINT_ENABLE_MAFFT)
    target_link_libraries(tests mafft_parttree)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(tests vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(tests sortmerna RocksDB::rocksdb)
    target_compile_definitions(tests PRIVATE MIINT_HAS_SORTMERNA)
endif()
if(WIN32)
    # Tests target only needs HTSlib deps (no Rust code linked here)
    target_link_libraries(tests regex ws2_32)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_link_libraries(tests
                          hdf5::hdf5-static
                          hdf5::hdf5_hl-static
                          hdf5::hdf5_cpp-static
                          hdf5::hdf5_hl_cpp-static)
endif()
if(HAVE_LIBZSTD)
    target_link_libraries(tests ${ZSTD_TARGET})
endif()
target_include_directories(tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_BINARY_DIR}/htslib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/miniz  # For ExtractFromZip tests
)
if(MIINT_HAS_JEMALLOC)
    target_sources(tests PRIVATE src/mm_alloc_stubs.c)
    target_compile_definitions(tests PRIVATE MIINT_USE_JEMALLOC)
endif()
