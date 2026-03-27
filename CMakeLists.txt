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

# Auto-disable platform-incompatible features.
# Bowtie2 requires POSIX fork/exec (unavailable on WASM and Windows).
if(EMSCRIPTEN)
    set(MIINT_ENABLE_BOWTIE2 OFF)
    message(STATUS "Emscripten detected: disabling Bowtie2 (WASM-incompatible)")
elseif(WIN32)
    set(MIINT_ENABLE_BOWTIE2 OFF)
    message(STATUS "Windows detected: disabling Bowtie2 (requires POSIX fork/exec)")
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
    # Emscripten: tell autotools we are cross-compiling, drop -fPIC (not meaningful)
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
    set(HTSLIB_CFLAGS "CFLAGS=${HTSLIB_ARCH_CFLAGS} ${EMSCRIPTEN_THREAD_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
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
    set(MINIMAP2_FPIC "")
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
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2/libminimap2.a
)

# Create imported library for minimap2
add_library(minimap2 STATIC IMPORTED GLOBAL)
set_target_properties(minimap2 PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2/libminimap2.a
)
add_dependencies(minimap2 minimap2_build)

# Handle cross-compilation arch flags for WFA2
set(WFA2_ARCH_CFLAGS "")
set(WFA2_FPIC "-fPIC")
if(EMSCRIPTEN)
    set(WFA2_FPIC "")
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
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 1
    BUILD_BYPRODUCTS
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfa.a
        ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfacpp.a
)

# Create imported libraries for WFA2 (C++ bindings depend on C core — link order matters)
add_library(wfa2 STATIC IMPORTED GLOBAL)
set_target_properties(wfa2 PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfa.a
)
add_dependencies(wfa2 wfa2_build)

add_library(wfa2cpp STATIC IMPORTED GLOBAL)
set_target_properties(wfa2cpp PROPERTIES
    IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib/lib/libwfacpp.a
)
add_dependencies(wfa2cpp wfa2_build)

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
ExternalProject_Add(
    rype_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype
    CONFIGURE_COMMAND "${RYPE_CONFIGURE_CMD}"
    BUILD_COMMAND ${CARGO} rustc --release ${RYPE_CARGO_FEATURES} --lib --crate-type=staticlib
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
    src/copy_format_common.cpp
    src/table_function_common.cpp
    src/copy_fastq.cpp
    src/copy_fasta.cpp
    src/copy_sam.cpp
    src/reference_table_reader.cpp
    src/placement_table_reader.cpp
    src/NewickTree.cpp
    src/read_newick.cpp
    src/copy_newick.cpp
    src/Minimap2Aligner.cpp
    src/sequence_table_reader.cpp
    src/align_minimap2.cpp
    src/align_minimap2_sharded.cpp
    src/save_minimap2_index.cpp
    src/ncbi_parser.cpp
    src/ncbi_client.cpp
    src/zip_utils.cpp
    src/read_ncbi_fasta.cpp
    src/read_ncbi.cpp
    src/read_ncbi_annotation.cpp
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
    src/mass_tables.cpp
    src/MzXMLReader.cpp
    src/read_mzxml.cpp
    src/remote_file_helper.cpp
    src/duckdb_seq_stream.cpp
    src/hfile_duckdb.cpp
    src/KmerIndex.cpp
    src/ChimeraDetector.cpp
    src/uchime_ref.cpp
    src/uchime_denovo.cpp)

if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    list(APPEND EXTENSION_SOURCES src/copy_biom.cpp src/BIOMTable.cpp src/BIOMReader.cpp src/read_biom.cpp)
endif()
if(MIINT_ENABLE_BOWTIE2)
    list(APPEND EXTENSION_SOURCES src/Bowtie2Aligner.cpp src/align_bowtie2.cpp src/align_bowtie2_sharded.cpp)
endif()

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
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
if(MIINT_ENABLE_BOWTIE2)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_BOWTIE2)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
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
if(MIINT_ENABLE_BOWTIE2)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_BOWTIE2)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
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
    src/NewickTree.cpp
    test/cpp/test_NewickParser.cpp
    test/cpp/test_InsertFullyResolved.cpp
    src/Minimap2Aligner.cpp
    test/cpp/test_Minimap2Aligner.cpp
    src/ncbi_parser.cpp
    src/zip_utils.cpp
    duckdb/third_party/miniz/miniz.cpp  # For ExtractFromZip tests
    test/cpp/test_ncbi_client.cpp
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
    src/KmerIndex.cpp
    test/cpp/test_KmerIndex.cpp
    src/ChimeraDetector.cpp
    test/cpp/test_ChimeraDetector.cpp
)
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    list(APPEND TEST_SOURCES src/BIOMTable.cpp src/BIOMReader.cpp
         test/cpp/test_BIOMReader.cpp test/cpp/test_BIOMTable.cpp)
endif()
if(MIINT_ENABLE_BOWTIE2)
    list(APPEND TEST_SOURCES src/Bowtie2Aligner.cpp test/cpp/test_Bowtie2Aligner.cpp)
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
