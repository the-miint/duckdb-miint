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

# libcurl powers the streaming HTTPS-PUT / FTP / FTPS transport in
# `ena_upload_reads`. Disabled on Apple because vcpkg's libcurl pulls in
# OpenSSL's libcrypto, which exports `MD5_*` / `SHA1_*` symbols that
# vsearch's vendored copies also export. GNU ld supports
# `--allow-multiple-definition` to resolve this; Apple's ld does not.
# Since vsearch is first-class on macOS and not all users need the
# streaming transport (Aspera + file:// still work), we drop libcurl on
# Apple rather than drop vsearch.
option(MIINT_ENABLE_CURL "Build with libcurl streaming-upload support" ON)
if(APPLE)
    set(MIINT_ENABLE_CURL OFF)
    message(STATUS "Apple detected: disabling libcurl streaming transport (vsearch+OpenSSL symbol clash)")
endif()
# vcpkg's libcurl port for wasm32-emscripten is not built with -pthread,
# so its object files lack atomics/bulk-memory features. wasm_threads links
# with --shared-memory which requires those features, so wasm-ld rejects the
# port. wasm_eh / wasm_mvp don't use shared memory and link fine. Rebuilding
# vcpkg libcurl with pthread is a port-level change; cheaper to drop the
# streaming transport on wasm_threads (file:// + Aspera paths still work).
if(EMSCRIPTEN AND (CMAKE_CXX_FLAGS MATCHES "-pthread" OR DUCKDB_EXPLICIT_PLATFORM STREQUAL "wasm_threads"))
    set(MIINT_ENABLE_CURL OFF)
    message(STATUS "Emscripten wasm_threads detected: disabling libcurl (vcpkg port not built with -pthread)")
endif()
if(MIINT_ENABLE_CURL)
    # Try VCPKG-style CONFIG package first (production builds), fall back
    # to the system FindCURL Module (apt-installed libcurl4-openssl-dev,
    # macOS Homebrew, etc.). Soft-disable if neither is available so the
    # tidy / no-vcpkg CI lanes don't fail at configure time before tidy
    # has a chance to run.
    find_package(CURL CONFIG QUIET)
    if(NOT CURL_FOUND)
        find_package(CURL QUIET)
    endif()
    if(NOT CURL_FOUND)
        message(STATUS "libcurl not found (no VCPKG package, no system FindCURL); disabling MIINT_ENABLE_CURL")
        set(MIINT_ENABLE_CURL OFF)
    endif()
endif()

# Optional feature flags for WASM-compatible builds
option(MIINT_ENABLE_HDF5 "Build with HDF5/BIOM support" ON)
option(MIINT_ENABLE_MAFFT "Build with MAFFT multiple sequence alignment support" ON)
option(MIINT_ENABLE_ABPOA "Build with abPOA partial-order alignment support" ON)
option(MIINT_ENABLE_SORTMERNA "Build with SortMeRNA rRNA alignment support" ON)
option(MIINT_ENABLE_SYLPH "Build with sylph FracMinHash relative-abundance profiling" ON)
option(MIINT_ENABLE_GPL_BOUNDARY "Build with gpl-boundary subsystem (process pipes + Arrow IPC over POSIX shm)" ON)
option(MIINT_ENABLE_UNIFRAC "Build with UniFrac (PCoA / PERMANOVA / Faith PD) support" ON)

# libdeflate: Linux-only accelerator for htslib bgzf gzip (FASTQ/FASTA/BAM COPY).
# Default ON so the shipped Linux binary actually gets the win (~3x over zlib-parallel bgzf,
# ~37x over the original single-threaded gzip writer; measured 2026-06-03). Hard-gated to
# Linux below: the Windows htslib build is a hand-written config.h/config.mk (no autotools
# autodetect) and wasm_eh is single-threaded in-browser, so neither benefits — the gate
# force-disables them. A missing vcpkg port also soft-disables (find_package below), so
# default-ON degrades gracefully everywhere it can't apply.
option(MIINT_ENABLE_LIBDEFLATE "Build htslib bgzf with libdeflate (Linux only)" ON)
if(MIINT_ENABLE_LIBDEFLATE AND NOT CMAKE_SYSTEM_NAME STREQUAL "Linux")
    message(STATUS "MIINT_ENABLE_LIBDEFLATE only supported on Linux; disabling here")
    set(MIINT_ENABLE_LIBDEFLATE OFF)
endif()
if(MIINT_ENABLE_LIBDEFLATE)
    find_package(libdeflate CONFIG QUIET)
    if(TARGET libdeflate::libdeflate_static)
        set(LIBDEFLATE_TARGET libdeflate::libdeflate_static)
    elseif(TARGET libdeflate::libdeflate_shared)
        set(LIBDEFLATE_TARGET libdeflate::libdeflate_shared)
    else()
        message(WARNING "MIINT_ENABLE_LIBDEFLATE=ON but no libdeflate package found; disabling")
        set(MIINT_ENABLE_LIBDEFLATE OFF)
    endif()
endif()

# libssu_inmem.a / libskbb_inmem.a are compiled with -fopenmp; the final
# extension must resolve their OpenMP runtime symbols (__kmpc_*, omp_*).
#
# Linux: find_package(OpenMP) discovers libgomp via system paths.
# macOS: Apple Clang doesn't ship libomp. Probe Homebrew's libomp (set up
#   via the `omp` extra_toolchain hook in our workflow) and link libomp.a
#   *statically* (positional arg, not -lomp) so the released
#   .duckdb_extension doesn't carry an RPATH to /opt/homebrew/opt/libomp.
#   Mirrors GPL-boundary's build.rs locate_macos_libomp().
#
# Exposes a single INTERFACE target `miint_openmp` so link sites below
# don't need to know which platform-specific dance produced it.
if(MIINT_ENABLE_UNIFRAC AND NOT EMSCRIPTEN)
    if(APPLE)
        execute_process(
            COMMAND brew --prefix libomp
            OUTPUT_VARIABLE BREW_LIBOMP_PREFIX
            OUTPUT_STRIP_TRAILING_WHITESPACE
            ERROR_QUIET)
        set(MIINT_LIBOMP_PREFIX "")
        foreach(candidate
                "${BREW_LIBOMP_PREFIX}"
                "/opt/homebrew/opt/libomp"
                "/usr/local/opt/libomp")
            if(candidate AND
               EXISTS "${candidate}/include/omp.h" AND
               EXISTS "${candidate}/lib/libomp.a")
                set(MIINT_LIBOMP_PREFIX "${candidate}")
                break()
            endif()
        endforeach()
        if(NOT MIINT_LIBOMP_PREFIX)
            message(FATAL_ERROR
                "OpenMP static archive (libomp.a) not found on macOS.\n"
                "  Install with: brew install libomp\n"
                "  Required when MIINT_ENABLE_UNIFRAC=ON.")
        endif()
        message(STATUS "Found Homebrew libomp at: ${MIINT_LIBOMP_PREFIX}")
        set(MIINT_LIBOMP_INCLUDE_DIR "${MIINT_LIBOMP_PREFIX}/include")
        # IMPORTED INTERFACE (not a regular INTERFACE library): exempts it
        # from install(EXPORT)'s "all linked targets must be in the export
        # set" check, mirroring how OpenMP::OpenMP_CXX itself behaves.
        add_library(miint_openmp INTERFACE IMPORTED GLOBAL)
        set_property(TARGET miint_openmp PROPERTY
            INTERFACE_INCLUDE_DIRECTORIES "${MIINT_LIBOMP_INCLUDE_DIR}")
        # Positional .a (not -lomp): Homebrew ships both libomp.a and
        # libomp.dylib; -lomp prefers the dylib, baking a Homebrew path
        # into the extension's RPATH.
        set_property(TARGET miint_openmp PROPERTY
            INTERFACE_LINK_LIBRARIES "${MIINT_LIBOMP_PREFIX}/lib/libomp.a")
    else()
        find_package(OpenMP REQUIRED COMPONENTS CXX)
        add_library(miint_openmp INTERFACE IMPORTED GLOBAL)
        set_property(TARGET miint_openmp PROPERTY
            INTERFACE_LINK_LIBRARIES OpenMP::OpenMP_CXX)
    endif()
endif()

# Auto-disable platform-incompatible features.
# MAFFT uses mkdtemp and other POSIX APIs (unavailable on Windows; segfaults on MinGW).
# HDF5 C++ static class members (e.g. LinkCreatPropList::DEFAULT) become
# unresolvable GOT.mem imports in WASM, causing "bad export type: undefined"
# at extension load time.
# SortMeRNA transitively pulls in cmph, whose cmph_time.h assumes POSIX
# <sys/time.h> / gettimeofday — MinGW headers don't provide these, so the
# sub-build fails with "storage size of 'e_time' isn't known". Upstream
# sortmerna CI tests MSVC, not MinGW; we disable it here rather than patching
# cmph. Still enabled on Linux/macOS/MSVC.
if(EMSCRIPTEN)
    set(MIINT_ENABLE_HDF5 OFF)
    set(MIINT_ENABLE_SORTMERNA OFF)
    set(MIINT_ENABLE_SYLPH OFF)
    set(MIINT_ENABLE_GPL_BOUNDARY OFF)
    message(STATUS "Emscripten detected: disabling HDF5, SortMeRNA, sylph, and gpl-boundary (WASM-incompatible)")
elseif(WIN32)
    set(MIINT_ENABLE_MAFFT OFF)
    set(MIINT_ENABLE_ABPOA OFF)
    set(MIINT_ENABLE_SYLPH OFF)
    set(MIINT_ENABLE_GPL_BOUNDARY OFF)
    set(MIINT_ENABLE_UNIFRAC OFF)
    if(MINGW)
        set(MIINT_ENABLE_SORTMERNA OFF)
        message(STATUS "MinGW detected: disabling MAFFT, SortMeRNA, sylph, gpl-boundary, and UniFrac (cmph + gpl-boundary rely on POSIX APIs; libssu inmem assumes POSIX)")
    else()
        message(STATUS "Windows detected: disabling MAFFT, sylph, gpl-boundary, and UniFrac (require POSIX APIs)")
    endif()
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

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/vsearch
    OUTPUT_VARIABLE VSEARCH_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT VSEARCH_GIT_VERSION)
    set(VSEARCH_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft
    OUTPUT_VARIABLE MAFFT_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT MAFFT_GIT_VERSION)
    set(MAFFT_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/abpoa
    OUTPUT_VARIABLE ABPOA_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT ABPOA_GIT_VERSION)
    set(ABPOA_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/unifrac-binaries
    OUTPUT_VARIABLE UNIFRAC_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT UNIFRAC_GIT_VERSION)
    set(UNIFRAC_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/scikit-bio-binaries
    OUTPUT_VARIABLE SKBB_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT SKBB_GIT_VERSION)
    set(SKBB_GIT_VERSION "unknown")
endif()

execute_process(COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/ext/sylph
    OUTPUT_VARIABLE SYLPH_GIT_VERSION OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
if(NOT SYLPH_GIT_VERSION)
    set(SYLPH_GIT_VERSION "unknown")
endif()

message(STATUS "Submodule versions: minimap2=${MINIMAP2_GIT_VERSION} WFA2=${WFA2_GIT_VERSION} rype=${RYPE_GIT_VERSION} vsearch=${VSEARCH_GIT_VERSION} mafft=${MAFFT_GIT_VERSION} abpoa=${ABPOA_GIT_VERSION} sylph=${SYLPH_GIT_VERSION} unifrac=${UNIFRAC_GIT_VERSION} skbb=${SKBB_GIT_VERSION}")

# Build htslib from release tarball
# Configure with zstd support if available (htslib autodetects via pkg-config)
if(HAVE_LIBZSTD OR MIINT_ENABLE_LIBDEFLATE)
    # Set PKG_CONFIG_PATH so htslib configure can find zstd / libdeflate
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

# -O2: passing CFLAGS to htslib's configure overrides its Makefile default
# (-g -Wall -O2), so an explicit -O is required or htslib builds at -O0. -O2
# matches htslib/samtools upstream's shipping optimization level.
if(EMSCRIPTEN)
    set(HTSLIB_CFLAGS "CFLAGS=-O2 -fPIC ${HTSLIB_ARCH_CFLAGS} ${EMSCRIPTEN_THREAD_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
else()
    set(HTSLIB_CFLAGS "CFLAGS=-O2 -fPIC ${HTSLIB_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
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

# libdeflate (Linux opt-in): point htslib's configure at the vcpkg-installed headers
# and lib. htslib probes libdeflate via AC_CHECK_HEADER/AC_CHECK_LIB (not pkg-config),
# so feed -I/-L directly rather than trusting libdeflate.pc. --with-libdeflate makes it
# mandatory: configure errors loudly if not found (Rule 10) instead of silently using zlib.
if(MIINT_ENABLE_LIBDEFLATE)
    set(HTSLIB_VCPKG_ROOT ${CMAKE_CURRENT_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET})
    set(HTSLIB_CFLAGS  "${HTSLIB_CFLAGS} -I${HTSLIB_VCPKG_ROOT}/include")
    set(HTSLIB_LDFLAGS "${HTSLIB_LDFLAGS} -L${HTSLIB_VCPKG_ROOT}/lib")
    set(HTSLIB_LIBDEFLATE_FLAG --with-libdeflate)
    message(STATUS "HTSlib will build WITH libdeflate from ${HTSLIB_VCPKG_ROOT}")
else()
    set(HTSLIB_LIBDEFLATE_FLAG --without-libdeflate)
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
        BUILD_COMMAND make -C ${HTSLIB_SRC} clean
            COMMAND make -C ${HTSLIB_SRC} libhts.a
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
            ${HTSLIB_LIBDEFLATE_FLAG}
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
    BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 clean
        COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 libminimap2.a "CFLAGS=-g -Wall -O2 -Wc++-compat ${MINIMAP2_FPIC} ${MINIMAP2_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}" "CPPFLAGS=-DHAVE_KALLOC ${MINIMAP2_ALLOC_FLAGS}" ${MINIMAP2_ARCH_FLAGS}
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
    BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib clean
        COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/WFA2-lib all
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
        # vsearch upstream gitignores `configure`, `Makefile.in`, etc. — they are
        # generated on demand by autogen.sh (autoreconf --force --install). A
        # fresh CI checkout has none of them, so we must regenerate here before
        # invoking configure.
        CONFIGURE_COMMAND ${CMAKE_COMMAND} -E chdir ${VSEARCH_SRC_DIR} sh ./autogen.sh
            # -O3: autoconf treats a user-supplied CXXFLAGS/CFLAGS as authoritative and
            # drops its default -O2, and vsearch's Makefile.am does not hardcode an
            # optimization level — so without an explicit -O here the whole library
            # (including the SIMD Needleman-Wunsch aligner) builds at -O0, ~4x slower.
            COMMAND ${CMAKE_COMMAND} -E env "CXXFLAGS=-O3 -fPIC ${VSEARCH_ARCH_CFLAGS}" "CFLAGS=-O3 -fPIC ${VSEARCH_ARCH_CFLAGS}"
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

# unifrac-binaries + scikit-bio-binaries: in-memory static archives for
# UniFrac (PCoA / PERMANOVA / Faith PD). Native uses `make inmem_static`
# (libssu_inmem.a / libskbb_inmem.a — Eigen-backed, no HDF5/lz4/BLAS),
# emscripten uses `make wasm` (libssu_wasm.a / libskbb_wasm.a — same
# subset, compiled with em++). Disabled on Windows.
if(MIINT_ENABLE_UNIFRAC)
    set(SKBB_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/scikit-bio-binaries)
    set(UNIFRAC_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/unifrac-binaries)

    if(EMSCRIPTEN)
        set(SKBB_LIB_NAME libskbb_wasm.a)
        set(UNIFRAC_LIB_NAME libssu_wasm.a)
        set(SKBB_MAKE_TARGET wasm)
        set(UNIFRAC_MAKE_TARGET wasm)
    else()
        set(SKBB_LIB_NAME libskbb_inmem.a)
        set(UNIFRAC_LIB_NAME libssu_inmem.a)
        set(SKBB_MAKE_TARGET inmem_static)
        set(UNIFRAC_MAKE_TARGET inmem_static)
    endif()

    set(SKBB_LIB_PATH ${SKBB_SRC_DIR}/src/${SKBB_LIB_NAME})
    set(UNIFRAC_LIB_PATH ${UNIFRAC_SRC_DIR}/src/${UNIFRAC_LIB_NAME})

    # skbb's inmem and wasm builds both depend on Eigen 3.4.0 in
    # .wasm-cache/eigen/include/. fetch_eigen.sh is idempotent — cache hit
    # returns immediately. Requires network on first configure (mirrors
    # the upstream WASM build's expectation).
    #
    # WASM PIC patch: both libskbb_wasm.a and libssu_wasm.a are linked
    # into the duckdb-miint extension as a WASM side module
    # (-sSIDE_MODULE=2), which requires every translation unit to be
    # position-independent. Upstream's wasm/emscripten_build.mk does NOT
    # include -fPIC in WASM_CXXFLAGS (the inmem path does, but the WASM
    # path was authored standalone and never linked into a side module
    # before). Patch -fPIC in idempotently before the WASM build. File
    # an upstream PR to land this in skbb / unifrac-binaries' WASM build
    # rules and drop these patches when both submodules ship the fix.
    #
    # Gated on EMSCRIPTEN to keep native builds from dirtying the
    # submodule worktree (the patch is a no-op on native, but the
    # touched-by-sed file still shows up in `git status` of the
    # submodule).
    # The patch script grep-guards the sed, then re-greps to fail loud if
    # the sed didn't actually match — sed -i exits 0 whether or not it
    # substituted, so a submodule advance that renames the assignment
    # would otherwise produce a cryptic R_WASM_MEMORY_ADDR_SLEB link
    # error downstream instead of a clear "patch didn't take" message.
    #
    # Written to a file rather than passed inline because CMake splits
    # `set(VAR "a;b;c")`-style strings on `;` when expanding them into a
    # command list, and shell flow control here needs semicolons.
    if(EMSCRIPTEN)
        set(WASM_PIC_PATCH_SH ${CMAKE_CURRENT_BINARY_DIR}/wasm_pic_patch.sh)
        file(WRITE ${WASM_PIC_PATCH_SH}
"#!/bin/sh
# Idempotent -fPIC patch for libssu/libskbb's WASM emscripten_build.mk.
# Generated by CMakeLists.txt; do not edit.
set -e
if grep -q '^WASM_CXXFLAGS := -fPIC' \"$1\"; then
    exit 0
fi
sed -i 's|^WASM_CXXFLAGS := |WASM_CXXFLAGS := -fPIC |' \"$1\"
if ! grep -q '^WASM_CXXFLAGS := -fPIC' \"$1\"; then
    echo \"ERROR: WASM -fPIC patch did not take in $1 — submodule advanced past the patch's assumptions; update CMakeLists.txt or land the fix upstream\" >&2
    exit 1
fi
")
        set(SKBB_WASM_PATCH sh ${WASM_PIC_PATCH_SH} ${SKBB_SRC_DIR}/src/wasm/emscripten_build.mk)
        set(UNIFRAC_WASM_PATCH sh ${WASM_PIC_PATCH_SH} ${UNIFRAC_SRC_DIR}/src/wasm/emscripten_build.mk)
    else()
        set(SKBB_WASM_PATCH ${CMAKE_COMMAND} -E true)
        set(UNIFRAC_WASM_PATCH ${CMAKE_COMMAND} -E true)
    endif()

    # ext/{scikit-bio-binaries,unifrac-binaries}/src/inmem_build.mk default
    # INMEM_MPFLAG=-fopenmp. Apple Clang rejects bare -fopenmp outright
    # ("unsupported option '-fopenmp'"); it requires -Xpreprocessor -fopenmp
    # to enable the pragmas at preprocess time, then libomp gets linked
    # separately. Also needs -I for libomp's headers since Apple's clang
    # doesn't auto-search Homebrew paths. Mirrors GPL-boundary build.rs.
    if(APPLE)
        set(MIINT_INMEM_MPFLAG_ARG "INMEM_MPFLAG=-Xpreprocessor -fopenmp -I${MIINT_LIBOMP_INCLUDE_DIR}")
    else()
        set(MIINT_INMEM_MPFLAG_ARG "")
    endif()

    ExternalProject_Add(
        skbb_build
        SOURCE_DIR ${SKBB_SRC_DIR}
        CONFIGURE_COMMAND bash ${SKBB_SRC_DIR}/scripts/fetch_eigen.sh
        BUILD_COMMAND ${SKBB_WASM_PATCH}
        COMMAND make -C ${SKBB_SRC_DIR}/src ${SKBB_MAKE_TARGET} ${MIINT_INMEM_MPFLAG_ARG}
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS ${SKBB_LIB_PATH}
    )

    # unifrac's inmem build expects SKBB_DIR pointing at the skbb source
    # tree (defaults to ../../scikit-bio-binaries which is wrong for our
    # layout). The wasm path uses the same variable.
    ExternalProject_Add(
        unifrac_build
        SOURCE_DIR ${UNIFRAC_SRC_DIR}
        CONFIGURE_COMMAND ""
        BUILD_COMMAND ${UNIFRAC_WASM_PATCH}
        COMMAND make -C ${UNIFRAC_SRC_DIR}/src ${UNIFRAC_MAKE_TARGET} SKBB_DIR=${SKBB_SRC_DIR} ${MIINT_INMEM_MPFLAG_ARG}
        INSTALL_COMMAND ""
        BUILD_IN_SOURCE 1
        BUILD_BYPRODUCTS ${UNIFRAC_LIB_PATH}
        DEPENDS skbb_build
    )

    add_library(libskbb STATIC IMPORTED GLOBAL)
    set_target_properties(libskbb PROPERTIES IMPORTED_LOCATION ${SKBB_LIB_PATH})
    add_dependencies(libskbb skbb_build)

    add_library(libssu STATIC IMPORTED GLOBAL)
    set_target_properties(libssu PROPERTIES IMPORTED_LOCATION ${UNIFRAC_LIB_PATH})
    add_dependencies(libssu unifrac_build)
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
    # The other -Wno-* flags silence cosmetic warnings from upstream
    # sortmerna/cmph/alp code (type-punned pointer casts, intentional
    # switch fall-through without [[fallthrough]], unchecked fread
    # return values). All would require upstream patches to fix properly.
    set(SORTMERNA_WARNING_SUPPRESSIONS "-Wno-register -Wno-strict-aliasing -Wno-implicit-fallthrough -Wno-unused-result")
    ExternalProject_Add(
        sortmerna_build
        SOURCE_DIR ${SORTMERNA_SRC_DIR}
        BINARY_DIR ${SORTMERNA_BIN_DIR}
        CMAKE_ARGS
            -DCMAKE_BUILD_TYPE=Release
            -DCMAKE_POSITION_INDEPENDENT_CODE=ON
            -DCMAKE_CXX_FLAGS=${SORTMERNA_WARNING_SUPPRESSIONS}
            -DCMAKE_C_FLAGS=${SORTMERNA_WARNING_SUPPRESSIONS}
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

# rype + sylph are independent Rust crates; each emits a staticlib that
# statically embeds its own copy of the Rust std runtime (rust_eh_personality,
# EMPTY_PANIC, ...). Linking both as separate .a archives produces
# duplicate-symbol errors on every supported linker (GNU ld, Apple ld-prime,
# Mach-O classic, lld) — see the Rust Reference (Linkage chapter).
#
# Fix: an umbrella `miint_rust_glue` crate (ext/miint-rust-glue/) lists both
# as path dependencies, cargo compiles them through a shared std, and the
# resulting libmiint_rust_glue.a contains a single deduplicated copy of std
# plus the #[no_mangle] C ABI surface from both crates. This is the canonical
# pattern (Mozilla's "rul" super-crate for Firefox/Gecko).
find_program(CARGO cargo REQUIRED)
message(STATUS "Found cargo: ${CARGO}")

set(GLUE_LIB_NAME "libmiint_rust_glue.a")
set(GLUE_TARGET_DIR ${CMAKE_CURRENT_BINARY_DIR}/miint_rust_glue_target)

# Cross-compile flags mirror the original per-crate logic.
set(GLUE_CARGO_TARGET_FLAG "")
set(GLUE_LIB_SUBDIR "release")
set(GLUE_CONFIGURE_CMD "")
if(EMSCRIPTEN)
    find_program(RUSTUP rustup)
    set(GLUE_CARGO_TARGET_FLAG "--target" "wasm32-unknown-emscripten")
    set(GLUE_LIB_SUBDIR "wasm32-unknown-emscripten/release")
    if(RUSTUP)
        set(GLUE_CONFIGURE_CMD ${RUSTUP} target add wasm32-unknown-emscripten)
    endif()
    message(STATUS "Rust glue cross-compiling for wasm32-unknown-emscripten")
elseif(WIN32)
    find_program(RUSTUP rustup)
    set(GLUE_CARGO_TARGET_FLAG "--target" "x86_64-pc-windows-gnu")
    set(GLUE_LIB_SUBDIR "x86_64-pc-windows-gnu/release")
    if(RUSTUP)
        set(GLUE_CONFIGURE_CMD ${RUSTUP} target add x86_64-pc-windows-gnu)
    endif()
    message(STATUS "Rust glue targeting x86_64-pc-windows-gnu for MinGW compatibility")
elseif(APPLE AND OSX_BUILD_ARCH)
    find_program(RUSTUP rustup)
    if(OSX_BUILD_ARCH STREQUAL "x86_64")
        set(GLUE_CARGO_TARGET_FLAG "--target" "x86_64-apple-darwin")
        set(GLUE_LIB_SUBDIR "x86_64-apple-darwin/release")
        if(RUSTUP)
            set(GLUE_CONFIGURE_CMD ${RUSTUP} target add x86_64-apple-darwin)
        endif()
    elseif(OSX_BUILD_ARCH STREQUAL "arm64")
        set(GLUE_CARGO_TARGET_FLAG "--target" "aarch64-apple-darwin")
        set(GLUE_LIB_SUBDIR "aarch64-apple-darwin/release")
        if(RUSTUP)
            set(GLUE_CONFIGURE_CMD ${RUSTUP} target add aarch64-apple-darwin)
        endif()
    endif()
endif()

# Feature selection: rype's fastx (needletail) is off on WASM and Windows
# (matches the original rype_build wiring). sylph is off on WASM/Windows
# (matches MIINT_ENABLE_SYLPH gating).
set(GLUE_FEATURES "")
if(NOT EMSCRIPTEN AND NOT WIN32)
    list(APPEND GLUE_FEATURES "rype-fastx")
endif()
if(MIINT_ENABLE_SYLPH)
    list(APPEND GLUE_FEATURES "with-sylph")
endif()
string(REPLACE ";" "," GLUE_FEATURES_CSV "${GLUE_FEATURES}")
if(GLUE_FEATURES_CSV)
    set(GLUE_FEATURES_ARGS "--features" "${GLUE_FEATURES_CSV}")
else()
    set(GLUE_FEATURES_ARGS "")
endif()

# Per-target RUSTFLAGS for the umbrella's cargo build. Cargo's RUSTFLAGS env
# var overrides .cargo/config.toml's rustflags rather than merging, so all
# platform-specific flags must accumulate here.
set(GLUE_RUSTFLAGS_LIST "")

if(EMSCRIPTEN)
    # PIC for WASM side modules. The umbrella itself is a staticlib (no link
    # step), but rype's [lib] declares crate-type = ["cdylib", "rlib"] and
    # cargo unconditionally builds the cdylib for the wasm32-unknown-emscripten
    # path-dep. -sSIDE_MODULE=2 tells wasm-ld that rype.wasm is a side module
    # rather than a standalone wasm — so no `main` symbol is required. We
    # never actually link rype.wasm into anything; only libmiint_rust_glue.a
    # goes into the extension. This just unblocks the cargo build.
    list(APPEND GLUE_RUSTFLAGS_LIST "-C" "relocation-model=pic")
    list(APPEND GLUE_RUSTFLAGS_LIST "-C" "link-arg=-sSIDE_MODULE=2")
endif()

if(WIN32 AND MINGW)
    # rtools42's MinGW64 doesn't ship libgcc_eh.a (modern GCC folded EH into
    # libgcc), but rustc's windows-gnu linker spec still passes -lgcc_eh. Same
    # situation as the WASM case above — cargo builds rype's cdylib
    # (rype.dll) even though we never link it. Stub an empty archive in a -L
    # search path so the cdylib link step finds something and proceeds.
    set(GLUE_MINGW_STUB_DIR ${CMAKE_CURRENT_BINARY_DIR}/mingw_libgcc_eh_stub)
    file(MAKE_DIRECTORY ${GLUE_MINGW_STUB_DIR})
    if(NOT EXISTS ${GLUE_MINGW_STUB_DIR}/libgcc_eh.a)
        execute_process(COMMAND ${CMAKE_AR} rcs ${GLUE_MINGW_STUB_DIR}/libgcc_eh.a)
    endif()
    list(APPEND GLUE_RUSTFLAGS_LIST "-C" "link-arg=-L${GLUE_MINGW_STUB_DIR}")
endif()

if(GLUE_RUSTFLAGS_LIST)
    string(JOIN " " GLUE_RUSTFLAGS_STR ${GLUE_RUSTFLAGS_LIST})
    set(GLUE_BUILD_ENV ${CMAKE_COMMAND} -E env "RUSTFLAGS=${GLUE_RUSTFLAGS_STR}")
else()
    set(GLUE_BUILD_ENV)
endif()

set(GLUE_LIB_PATH ${GLUE_TARGET_DIR}/${GLUE_LIB_SUBDIR}/${GLUE_LIB_NAME})

ExternalProject_Add(
    miint_rust_glue_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/miint-rust-glue
    CONFIGURE_COMMAND "${GLUE_CONFIGURE_CMD}"
    BUILD_COMMAND ${GLUE_BUILD_ENV} ${CARGO} build --release
        ${GLUE_FEATURES_ARGS}
        --target-dir ${GLUE_TARGET_DIR}
        --manifest-path ${CMAKE_CURRENT_SOURCE_DIR}/ext/miint-rust-glue/Cargo.toml
        ${GLUE_CARGO_TARGET_FLAG}
    INSTALL_COMMAND ""
    BUILD_IN_SOURCE 0
    # Run the cargo build on EVERY build, deferring change detection to cargo's own
    # fingerprinting. ExternalProject does not track SOURCE_DIR (or path-dependency,
    # e.g. ext/rype / ext/sylph) file contents, so without this an incremental build
    # after a submodule bump or a local rype/sylph edit would silently relink the
    # STALE libmiint_rust_glue.a — the dev gets a binary that doesn't match their
    # source and has no way to notice. Cargo's no-op when nothing changed is ~1-3s.
    # (Clean builds — CI's fresh ubuntu-latest runners, or `make clean` locally —
    # were always correct; this only closes the incremental-build footgun.)
    BUILD_ALWAYS 1
    BUILD_BYPRODUCTS ${GLUE_LIB_PATH}
)

# `rype` is the sole IMPORTED STATIC target for the merged archive. The
# archive contains both rype and sylph symbols, so linking `rype` provides
# everything — there is no separate `sylph` link target.
#
# Why no second target: declaring an additional IMPORTED STATIC target
# pointing at the same archive caused libmiint_rust_glue.a to appear twice on
# the link line. Combined with -Wl,--allow-multiple-definition + --gc-sections,
# each pass resolved the Rust runtime symbols (__rust_alloc /
# rust_eh_personality / EMPTY_PANIC / ...) from a different instance and
# pruned the other; allocations made through one resolution freed via the
# other corrupted the heap, producing non-deterministic SIGSEGVs in unrelated
# code paths (mzXML reader, abPOA consensus, ...) on Linux glibc CI runners.
# Local Red Hat builds happened to dodge it.
add_library(rype STATIC IMPORTED GLOBAL)
set_target_properties(rype PROPERTIES IMPORTED_LOCATION ${GLUE_LIB_PATH})
add_dependencies(rype miint_rust_glue_build)

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
        BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core clean
            COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core libmafft_parttree.a
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

# Build abPOA partial-order alignment library via add_subdirectory
# (the fork's CMakeLists.txt is subdirectory-friendly: no global flag
# pollution, gated install targets, multi-ISA SIMD dispatch + WASM)
if(MIINT_ENABLE_ABPOA)
    if(APPLE AND OSX_BUILD_ARCH)
        set(CMAKE_OSX_ARCHITECTURES ${OSX_BUILD_ARCH})
    endif()
    add_subdirectory(ext/abpoa EXCLUDE_FROM_ALL)
    if(APPLE AND OSX_BUILD_ARCH)
        unset(CMAKE_OSX_ARCHITECTURES)
    endif()
endif()

# Backward-compat alias: the WASM linked-libs block below still references
# RYPE_LIB_PATH. rype now lives inside libmiint_rust_glue.a (built above).
set(RYPE_LIB_PATH ${GLUE_LIB_PATH})

# adapted kseq++ include use from strobealign
include_directories(src/include third_party/microtar ext ext/WFA2-lib ${CMAKE_CURRENT_SOURCE_DIR}/ext/rype)
if(MIINT_ENABLE_SYLPH)
    include_directories(${CMAKE_CURRENT_SOURCE_DIR}/ext/sylph)
endif()
if(MIINT_ENABLE_MAFFT)
    include_directories(ext/mafft/core)
endif()
if(MIINT_ENABLE_ABPOA)
    include_directories(ext/abpoa/include)
endif()
if(MIINT_ENABLE_VSEARCH)
    include_directories(${VSEARCH_SRC_DIR}/src)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    add_subdirectory(third_party/nanoarrow_ipc)
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
    src/qc_functions.cpp
    src/qc_algorithms.cpp
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
    src/shear_tree.cpp
    src/copy_newick.cpp
    src/Minimap2Aligner.cpp
    src/sequence_table_reader.cpp
    src/sequence_data_reader.cpp
    src/align_minimap2.cpp
    src/align_minimap2_sharded.cpp
    src/save_minimap2_index.cpp
    src/ncbi_parser.cpp
    src/ncbi_client.cpp
    src/taxdump_parser.cpp
    src/taxdump_archive.cpp
    src/taxonomy_lineage.cpp
    third_party/microtar/microtar.c
    src/blast_parser.cpp
    src/blast_client.cpp
    src/blast_search.cpp
    src/zip_utils.cpp
    src/read_ncbi_fasta.cpp
    src/read_ncbi.cpp
    src/read_ncbi_annotation.cpp
    src/read_ncbi_taxdump.cpp
    src/read_ncbi_lineage.cpp
    src/miint_cache.cpp
    src/ensure_httpfs.cpp
    src/ena_parser.cpp
    src/ena_client.cpp
    src/ena_run_info_extractor.cpp
    src/ena_resolver_cache.cpp
    src/ena_search_field_registry.cpp
    src/ena_attributes_filter.cpp
    src/ena_secret.cpp
    src/ena_envelope_builder.cpp
    src/ena_receipt_parser.cpp
    src/ena_lifecycle_submit.cpp
    src/ena_lifecycle_functions.cpp
    src/ena_lifecycle_delete.cpp
    src/ena_modify_functions.cpp
    src/ena_submission_log.cpp
    src/ena_storage.cpp
    src/ena_insert_common.cpp
    src/ena_alias_check.cpp
    src/ena_reports_client.cpp
    src/ena_checklist.cpp
    src/ena_projects_insert.cpp
    src/ena_projects_insert_op.cpp
    src/ena_samples_insert.cpp
    src/ena_samples_insert_op.cpp
    src/ena_experiments_insert.cpp
    src/ena_experiments_insert_op.cpp
    src/ena_runs_insert.cpp
    src/ena_runs_insert_op.cpp
    src/fastq_encoder.cpp
    src/ena_upload_helpers.cpp
    src/ena_upload_reads.cpp
    src/ena_lifecycle.cpp
    src/miint_log.cpp
    src/per_run_reader.cpp
    src/read_ena.cpp
    src/read_ena_attributes.cpp
    src/read_ena_searchable_fields.cpp
    src/read_ena_sequences.cpp
    src/aspera_utils.cpp
    src/aspera_stream.cpp
    src/aspera_send.cpp
    src/SFFReader.cpp
    src/read_sequences_sff.cpp
    src/WFA2Aligner.cpp
    src/align_pairwise_wfa2_functions.cpp
    src/KSW2Aligner.cpp
    src/align_pairwise_ksw2_functions.cpp
    src/align_pairwise_ksw2_dual_affine_functions.cpp
    src/align_pairwise_ksw2_splice_functions.cpp
    src/rype_classify.cpp
    src/rype_extract.cpp
    src/rype_index_create.cpp
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
    src/woltka_ogu_function.cpp
    src/mzml_peak_pair_function.cpp
    src/mass_tables.cpp
    src/MzXMLReader.cpp
    src/read_mzxml.cpp
    src/remote_file_helper.cpp
    src/duckdb_seq_stream.cpp
    src/hfile_duckdb.cpp
    src/deblur.cpp
    src/deblur_table_function.cpp
    src/per_sample_table_function.cpp
    src/id_column_codec.cpp
    src/extract_linked_amplicon.cpp
    src/match_short_barcodes.cpp
    src/compute_pileup.cpp
    src/PileupWalker.cpp
    src/compute_msa_consensus.cpp
    src/MsaConsensus.cpp)

if(MIINT_ENABLE_VSEARCH)
    list(APPEND EXTENSION_SOURCES
        src/VsearchChimeraWrapper.cpp src/VsearchSearchWrapper.cpp src/VsearchClusterWrapper.cpp
        src/mask_dust_function.cpp src/merge_pairs_function.cpp src/search_sequences.cpp
        src/cluster_sequences.cpp src/uchime_common.cpp src/uchime_ref.cpp src/uchime_denovo.cpp)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    list(APPEND EXTENSION_SOURCES src/copy_biom.cpp src/BIOMTable.cpp src/BIOMReader.cpp src/read_biom.cpp)
endif()
if(MIINT_ENABLE_MAFFT)
    list(APPEND EXTENSION_SOURCES src/MafftAligner.cpp src/align_mafft.cpp)
endif()
if(MIINT_ENABLE_ABPOA)
    list(APPEND EXTENSION_SOURCES src/AbpoaAligner.cpp src/align_abpoa.cpp src/consensus_abpoa.cpp)
endif()
if(MIINT_ENABLE_SORTMERNA)
    list(APPEND EXTENSION_SOURCES src/SortMeRNAAligner.cpp src/align_sortmerna.cpp src/align_sortmerna_rrna.cpp)
endif()
if(MIINT_ENABLE_SYLPH)
    list(APPEND EXTENSION_SOURCES src/SylphDatabase.cpp src/sylph_profile.cpp src/sylph_index_create.cpp)
endif()
if(MIINT_ENABLE_CURL)
    list(APPEND EXTENSION_SOURCES src/curl_send.cpp)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    list(APPEND EXTENSION_SOURCES
        src/gpl_boundary/process.cpp
        src/gpl_boundary/session.cpp
        src/gpl_boundary/shm.cpp
        src/gpl_boundary/arrow_ipc.cpp
        src/install_gpl_boundary.cpp
        src/phylogeny_fasttree.cpp
        src/align_bowtie2_daemon_common.cpp
        src/align_bowtie2.cpp
        src/align_bowtie2_sharded.cpp
        src/save_bowtie2_index.cpp)
endif()
if(MIINT_ENABLE_UNIFRAC)
    list(APPEND EXTENSION_SOURCES
        src/unifrac_support_biom.cpp
        src/unifrac_bptree.cpp
        src/unifrac_metadata.cpp
        src/unifrac_distance.cpp
        src/unifrac_omp_scope.cpp
        src/unifrac_subsample_bridge.cpp
        src/unifrac_function_common.cpp
        src/unifrac_pcoa_function.cpp
        src/unifrac_permanova_function.cpp
        src/unifrac_faith_pd_function.cpp)
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
    if(MIINT_ENABLE_ABPOA)
        string(APPEND WASM_LIBS " $<TARGET_FILE:abpoa>")
    endif()
    if(MIINT_ENABLE_UNIFRAC)
        string(APPEND WASM_LIBS " ${UNIFRAC_LIB_PATH}")
        string(APPEND WASM_LIBS " ${SKBB_LIB_PATH}")
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
    MAFFT_GIT_VERSION="${MAFFT_GIT_VERSION}"
    ABPOA_GIT_VERSION="${ABPOA_GIT_VERSION}"
    $<$<BOOL:${MIINT_HAS_JEMALLOC}>:MIINT_USE_JEMALLOC>
    HAVE_KALLOC)
if(MIINT_ENABLE_VSEARCH)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_VSEARCH
        VSEARCH_GIT_VERSION="${VSEARCH_GIT_VERSION}")
endif()
if(MIINT_ENABLE_MAFFT)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_MAFFT)
endif()
if(MIINT_ENABLE_ABPOA)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_ABPOA)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
endif()
if(MIINT_ENABLE_LIBDEFLATE)
    # Include dirs propagate from the linked imported target's interface.
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_LIBDEFLATE)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_SORTMERNA)
endif()
if(MIINT_ENABLE_SYLPH)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_SYLPH
        SYLPH_GIT_VERSION="${SYLPH_GIT_VERSION}")
endif()
if(MIINT_ENABLE_CURL)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_CURL)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_GPL_BOUNDARY)
endif()
if(MIINT_ENABLE_UNIFRAC)
    target_compile_definitions(${EXTENSION_NAME} PRIVATE MIINT_HAS_UNIFRAC
        UNIFRAC_GIT_VERSION="${UNIFRAC_GIT_VERSION}"
        SKBB_GIT_VERSION="${SKBB_GIT_VERSION}")
    target_include_directories(${EXTENSION_NAME} PRIVATE
        ${UNIFRAC_SRC_DIR}/src
        ${SKBB_SRC_DIR}/src/extern)
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
if(MIINT_ENABLE_ABPOA)
    target_link_libraries(${EXTENSION_NAME} $<TARGET_FILE:abpoa>)
    add_dependencies(${EXTENSION_NAME} abpoa)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(${EXTENSION_NAME} vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(${EXTENSION_NAME} sortmerna RocksDB::rocksdb)
endif()
# MIINT_ENABLE_SYLPH: sylph's symbols ship inside libmiint_rust_glue.a alongside
# rype's, so the existing `rype` link above already provides them. No separate
# target_link_libraries(... sylph) is needed; adding one would double-list the
# archive on the link line — see the comment near the rype target declaration.
if(MIINT_ENABLE_CURL)
    target_link_libraries(${EXTENSION_NAME} CURL::libcurl)
endif()
if(MIINT_ENABLE_UNIFRAC)
    # libssu (unifrac-binaries) references libskbb (skbio) symbols, so libssu
    # must precede libskbb in link order for static-archive resolution.
    # Native libssu_inmem.a is compiled with -fopenmp (see ext/unifrac-binaries/
    # src/inmem_build.mk INMEM_MPFLAG); libssu_wasm.a is OpenMP-free per
    # ext/unifrac-binaries/src/wasm/emscripten_build.mk header. Only link
    # libgomp on non-EMSCRIPTEN.
    target_link_libraries(${EXTENSION_NAME} libssu libskbb)
    if(NOT EMSCRIPTEN)
        target_link_libraries(${EXTENSION_NAME} miint_openmp)
    endif()
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    # Object lib: feed its compiled .o's into the extension target without
    # registering a link dependency (would require nanoarrow_ipc to also be
    # in the DuckDB export set, which it has no business being in).
    target_sources(${EXTENSION_NAME} PRIVATE $<TARGET_OBJECTS:miint_nanoarrow_ipc>)
    target_include_directories(${EXTENSION_NAME} PRIVATE
        $<TARGET_PROPERTY:miint_nanoarrow_ipc,INTERFACE_INCLUDE_DIRECTORIES>)
    # yyjson — bundled in DuckDB as a static lib (target `duckdb_yyjson`) at
    # `duckdb/third_party/yyjson/CMakeLists.txt`. We just consume its include
    # directories; the symbols are linked transitively via the duckdb library.
    target_include_directories(${EXTENSION_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/include)
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
    if(MIINT_ENABLE_GPL_BOUNDARY AND NOT APPLE)
        # gpl_boundary/shm.cpp uses shm_open/shm_unlink. On glibc these live in
        # librt (older RHEL/CentOS toolchains require explicit linkage; newer
        # glibc folded them into libc but we link defensively). macOS exposes
        # them in libSystem so we skip the flag there.
        target_link_libraries(${EXTENSION_NAME} rt)
    endif()
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_link_libraries(${EXTENSION_NAME}
                          hdf5::hdf5-static
                          hdf5::hdf5_cpp-static
    )
    target_link_libraries(${LOADABLE_EXTENSION_NAME}
                          hdf5::hdf5-static
                          hdf5::hdf5_cpp-static
    )
endif()
if(HAVE_LIBZSTD)
    target_link_libraries(${EXTENSION_NAME} ${ZSTD_TARGET})
endif()
# libdeflate appended AFTER the hts-containing block above, so -ldeflate follows -lhts.
if(MIINT_ENABLE_LIBDEFLATE)
    target_link_libraries(${EXTENSION_NAME} ${LIBDEFLATE_TARGET})
endif()

target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE RYPE_ARROW
    MIINT_STATIC_BUILD
    MINIMAP2_GIT_VERSION="${MINIMAP2_GIT_VERSION}"
    WFA2_GIT_VERSION="${WFA2_GIT_VERSION}"
    RYPE_GIT_VERSION="${RYPE_GIT_VERSION}"
    MAFFT_GIT_VERSION="${MAFFT_GIT_VERSION}"
    ABPOA_GIT_VERSION="${ABPOA_GIT_VERSION}"
    $<$<BOOL:${MIINT_HAS_JEMALLOC}>:MIINT_USE_JEMALLOC>
    HAVE_KALLOC)
if(MIINT_ENABLE_VSEARCH)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_VSEARCH
        VSEARCH_GIT_VERSION="${VSEARCH_GIT_VERSION}")
endif()
if(MIINT_ENABLE_MAFFT)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_MAFFT)
endif()
if(MIINT_ENABLE_ABPOA)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_ABPOA)
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_HDF5)
endif()
if(MIINT_ENABLE_LIBDEFLATE)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_LIBDEFLATE)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_SORTMERNA)
endif()
if(MIINT_ENABLE_SYLPH)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_SYLPH
        SYLPH_GIT_VERSION="${SYLPH_GIT_VERSION}")
endif()
if(MIINT_ENABLE_CURL)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_CURL)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_GPL_BOUNDARY)
endif()
if(MIINT_ENABLE_UNIFRAC)
    target_compile_definitions(${LOADABLE_EXTENSION_NAME} PRIVATE MIINT_HAS_UNIFRAC
        UNIFRAC_GIT_VERSION="${UNIFRAC_GIT_VERSION}"
        SKBB_GIT_VERSION="${SKBB_GIT_VERSION}")
    target_include_directories(${LOADABLE_EXTENSION_NAME} PRIVATE
        ${UNIFRAC_SRC_DIR}/src
        ${SKBB_SRC_DIR}/src/extern)
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
if(MIINT_ENABLE_ABPOA)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} $<TARGET_FILE:abpoa>)
    add_dependencies(${LOADABLE_EXTENSION_NAME} abpoa)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} sortmerna RocksDB::rocksdb)
endif()
# sylph's symbols come in via rype's archive — see comment in EXTENSION_NAME
# block above.
if(MIINT_ENABLE_CURL)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} CURL::libcurl)
endif()
if(MIINT_ENABLE_UNIFRAC)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} libssu libskbb)
    if(NOT EMSCRIPTEN)
        target_link_libraries(${LOADABLE_EXTENSION_NAME} miint_openmp)
    endif()
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    target_sources(${LOADABLE_EXTENSION_NAME} PRIVATE $<TARGET_OBJECTS:miint_nanoarrow_ipc>)
    target_include_directories(${LOADABLE_EXTENSION_NAME} PRIVATE
        $<TARGET_PROPERTY:miint_nanoarrow_ipc,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(${LOADABLE_EXTENSION_NAME} PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/include)
endif()
if(APPLE)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} "-framework CoreFoundation")
elseif(WIN32)
    target_link_libraries(${LOADABLE_EXTENSION_NAME}
        regex ws2_32 ntdll userenv bcrypt advapi32)
elseif(UNIX)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} dl m)
    if(MIINT_ENABLE_GPL_BOUNDARY AND NOT APPLE)
        target_link_libraries(${LOADABLE_EXTENSION_NAME} rt)
    endif()
endif()
if(HAVE_LIBZSTD)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} ${ZSTD_TARGET})
endif()
if(MIINT_ENABLE_LIBDEFLATE)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} ${LIBDEFLATE_TARGET})
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
    test/cpp/test_DuckDBSeqStream.cpp
    test/cpp/test_CopyFileHandlePipe.cpp
    src/QualScore.cpp
    test/cpp/test_QualScore.cpp
    src/qc_algorithms.cpp
    test/cpp/test_qc_algorithms.cpp
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
    test/cpp/test_ShearTree.cpp
    src/Minimap2Aligner.cpp
    test/cpp/test_Minimap2Aligner.cpp
    test/cpp/test_AlignBowtie2Sharded.cpp
    test/cpp/test_ShardProgress.cpp
    src/ncbi_parser.cpp
    src/zip_utils.cpp
    duckdb/third_party/miniz/miniz.cpp  # For ExtractFromZip tests
    test/cpp/test_ncbi_client.cpp
    src/blast_parser.cpp
    test/cpp/test_blast_parser.cpp
    src/taxdump_parser.cpp
    test/cpp/test_taxdump_parser.cpp
    src/taxdump_archive.cpp
    third_party/microtar/microtar.c
    test/cpp/test_taxdump_archive.cpp
    src/taxonomy_lineage.cpp
    test/cpp/test_taxonomy_lineage.cpp
    src/ena_parser.cpp
    src/ena_run_info_extractor.cpp
    src/ena_resolver_cache.cpp
    src/ena_search_field_registry.cpp
    src/ena_attributes_filter.cpp
    test/cpp/test_ENAClient.cpp
    test/cpp/test_ENAResolverCache.cpp
    test/cpp/test_ena_search_field_registry.cpp
    test/cpp/test_ena_attributes_filter.cpp
    test/cpp/test_ena_attributes_pushdown.cpp
    test/cpp/test_ena_secret.cpp
    src/ena_envelope_builder.cpp
    test/cpp/test_ena_envelope.cpp
    src/ena_receipt_parser.cpp
    test/cpp/test_ena_receipt.cpp
    src/ena_lifecycle_submit.cpp
    test/cpp/test_ena_lifecycle_submit.cpp
    src/ena_projects_insert.cpp
    test/cpp/test_ena_projects_insert.cpp
    src/ena_samples_insert.cpp
    test/cpp/test_ena_samples_insert.cpp
    src/ena_experiments_insert.cpp
    test/cpp/test_ena_experiments_insert.cpp
    src/ena_runs_insert.cpp
    test/cpp/test_ena_runs_insert.cpp
    src/ena_alias_check.cpp
    test/cpp/test_ena_alias_check.cpp
    src/ena_reports_client.cpp
    test/cpp/test_ena_reports_client.cpp
    src/ena_checklist.cpp
    test/cpp/test_ena_checklist.cpp
    src/fastq_encoder.cpp
    test/cpp/test_fastq_encoder.cpp
    src/ena_upload_helpers.cpp
    test/cpp/test_ena_upload_reads.cpp
    src/aspera_utils.cpp
    test/cpp/test_AsperaUtils.cpp
    src/SFFReader.cpp
    test/cpp/test_SFFReader.cpp
    src/WFA2Aligner.cpp
    test/cpp/test_WFA2Aligner.cpp
    src/KSW2Aligner.cpp
    test/cpp/test_KSW2Aligner.cpp
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
    src/id_column_codec.cpp
    test/cpp/test_id_column_codec.cpp
    src/PileupWalker.cpp
    test/cpp/test_compute_pileup.cpp
    src/MsaConsensus.cpp
    test/cpp/test_consensus_msa.cpp
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
if(MIINT_ENABLE_MAFFT)
    list(APPEND TEST_SOURCES src/MafftAligner.cpp test/cpp/test_MafftAligner.cpp)
endif()
if(MIINT_ENABLE_ABPOA)
    list(APPEND TEST_SOURCES src/AbpoaAligner.cpp test/cpp/test_AbpoaAligner.cpp)
endif()
if(MIINT_ENABLE_SORTMERNA)
    list(APPEND TEST_SOURCES
        src/SortMeRNAAligner.cpp
        test/cpp/test_sortmerna_build.cpp
        test/cpp/test_SortMeRNAAligner.cpp)
endif()
if(MIINT_ENABLE_SYLPH)
    # Tests only exercise SylphDatabaseHandle + the FFI default getters from
    # libsylph.a; sylph_profile.cpp itself isn't linked into the tests binary
    # because it pulls in DuckDB-internal symbols (ArrowArrayWrapper) +
    # miint's PerSampleGlobalState that the tests target doesn't otherwise
    # depend on. The behavioral coverage for sylph_profile lives in the SQL
    # tests, which run against the full extension binary.
    list(APPEND TEST_SOURCES
        src/SylphDatabase.cpp
        test/cpp/test_SylphDatabase.cpp)
endif()
if(MIINT_ENABLE_CURL)
    list(APPEND TEST_SOURCES src/curl_send.cpp test/cpp/test_curl_send.cpp)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    list(APPEND TEST_SOURCES
        src/gpl_boundary/process.cpp
        src/gpl_boundary/session.cpp
        src/gpl_boundary/shm.cpp
        src/gpl_boundary/arrow_ipc.cpp
        test/cpp/test_gpl_boundary_process.cpp
        test/cpp/test_gpl_boundary_session.cpp
        test/cpp/test_gpl_boundary_shm.cpp
        test/cpp/test_gpl_boundary_arrow_ipc.cpp)
endif()
if(MIINT_ENABLE_UNIFRAC)
    # NewickTree.cpp already in TEST_SOURCES (line 1265) for the existing
    # NewickParser tests; don't duplicate.
    list(APPEND TEST_SOURCES
        src/unifrac_support_biom.cpp
        src/unifrac_bptree.cpp
        src/unifrac_metadata.cpp
        src/unifrac_distance.cpp
        src/unifrac_omp_scope.cpp
        src/unifrac_subsample_bridge.cpp
        test/cpp/test_UnifracSupportBiom.cpp
        test/cpp/test_UnifracBptree.cpp
        test/cpp/test_UnifracMetadata.cpp
        test/cpp/test_UnifracDistance.cpp
        test/cpp/test_UnifracSubsampleBridge.cpp)
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
if(MIINT_ENABLE_CURL)
    target_link_libraries(tests CURL::libcurl)
    target_compile_definitions(tests PRIVATE MIINT_HAS_CURL)
endif()
if(MIINT_ENABLE_MAFFT)
    target_link_libraries(tests mafft_parttree)
endif()
if(MIINT_ENABLE_ABPOA)
    target_link_libraries(tests $<TARGET_FILE:abpoa>)
    add_dependencies(tests abpoa)
    target_compile_definitions(tests PRIVATE MIINT_HAS_ABPOA)
endif()
if(MIINT_ENABLE_VSEARCH)
    target_link_libraries(tests vsearch)
endif()
if(MIINT_ENABLE_SORTMERNA)
    target_link_libraries(tests sortmerna RocksDB::rocksdb)
    target_compile_definitions(tests PRIVATE MIINT_HAS_SORTMERNA)
endif()
if(MIINT_ENABLE_SYLPH)
    # Link rype (the umbrella archive), which transitively provides sylph
    # symbols — see comment near the rype IMPORTED target declaration.
    target_link_libraries(tests rype)
    target_compile_definitions(tests PRIVATE MIINT_HAS_SYLPH)
endif()
if(MIINT_ENABLE_GPL_BOUNDARY)
    target_sources(tests PRIVATE $<TARGET_OBJECTS:miint_nanoarrow_ipc>)
    target_include_directories(tests PRIVATE
        $<TARGET_PROPERTY:miint_nanoarrow_ipc,INTERFACE_INCLUDE_DIRECTORIES>)
    target_include_directories(tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/include)
    # The catch2 tests binary doesn't link against libduckdb (it's a
    # standalone executable), so it needs its own copy of yyjson.cpp.
    target_sources(tests PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/yyjson.cpp)
    target_compile_definitions(tests PRIVATE MIINT_HAS_GPL_BOUNDARY)
    if(UNIX AND NOT APPLE)
        # shm_open/shm_unlink in glibc — see comment on the EXTENSION_NAME side.
        target_link_libraries(tests rt)
    endif()
endif()
if(WIN32)
    # Tests target only needs HTSlib deps (no Rust code linked here)
    target_link_libraries(tests regex ws2_32)
endif()
if(APPLE)
    # libmiint_rust_glue.a transitively pulls rype's chrono → iana_time_zone,
    # which calls CoreFoundation symbols (CFRelease, CFTimeZoneCopySystem, ...).
    # The extension targets already link CoreFoundation; the tests target
    # historically didn't, because pre-glue it only linked libsylph.a.
    target_link_libraries(tests "-framework CoreFoundation")
endif()
if(MIINT_ENABLE_HDF5 AND HDF5_FOUND)
    target_link_libraries(tests
                          hdf5::hdf5-static
                          hdf5::hdf5_cpp-static)
endif()
# vsearch ships its own MD5/SHA1 symbols (libvsearch_core_a-md5.o,
# libvsearch_core_a-sha1.o), and OpenSSL's libcrypto — pulled in transitively
# when we link CURL::libcurl — exports the same names. GNU ld's
# --allow-multiple-definition picks the first-defined symbol; since each
# archive uses its own copy at compile time (intra-.o calls stay within the
# archive), runtime behaviour is preserved. libcurl itself uses `EVP_Digest*`
# rather than the colliding `MD5_*` C API, so the choice is benign.
#
# This workaround is only needed when both vsearch AND libcurl are linked.
# Apple disables MIINT_ENABLE_CURL (see top of file) precisely because
# Apple's ld lacks an equivalent flag and vsearch must remain first-class on
# macOS; Windows auto-disables vsearch. So we only need this on Linux.
if(MIINT_ENABLE_VSEARCH AND MIINT_ENABLE_CURL AND UNIX AND NOT APPLE)
    target_link_options(tests PRIVATE "-Wl,--allow-multiple-definition")
    # The static `${EXTENSION_NAME}` archive is consumed by DuckDB-built
    # binaries (duckdb shell, unittest, plan_serializer) that we don't
    # control. INTERFACE makes the link option propagate to whichever
    # executable links against this archive, so the duplicate-symbol
    # tolerance follows the static lib instead of being a per-target chore.
    target_link_options(${EXTENSION_NAME} INTERFACE "-Wl,--allow-multiple-definition")
    target_link_options(${LOADABLE_EXTENSION_NAME} PRIVATE "-Wl,--allow-multiple-definition")
endif()

# Windows MinGW: duckdb's LogicalType static members (VARCHAR, DOUBLE, BIGINT, ...)
# are emitted into COMDAT sections in every translation unit that references them.
# GNU ld on Linux folds these by default; rtools42's MinGW ld on Windows does not,
# so the duckdb.dll link sees duplicate `duckdb::LogicalType::VARCHAR` symbols
# between libmiint_extension.a's TUs and duckdb_common's TUs. Pre-sylph-merge this
# was masked by --allow-multiple-definition being applied for the unrelated
# vsearch+curl MD5/SHA1 collision; that flag is now gated on
# `MIINT_ENABLE_VSEARCH AND MIINT_ENABLE_CURL`, both of which are OFF on Windows.
# Re-add the flag for MinGW with the actual C++ COMDAT reason.
if(WIN32 AND MINGW)
    target_link_options(tests PRIVATE "-Wl,--allow-multiple-definition")
    target_link_options(${EXTENSION_NAME} INTERFACE "-Wl,--allow-multiple-definition")
    target_link_options(${LOADABLE_EXTENSION_NAME} PRIVATE "-Wl,--allow-multiple-definition")
endif()

if(HAVE_LIBZSTD)
    target_link_libraries(tests ${ZSTD_TARGET})
endif()
if(MIINT_ENABLE_LIBDEFLATE)
    target_link_libraries(tests ${LIBDEFLATE_TARGET})
endif()
if(MIINT_ENABLE_UNIFRAC)
    target_include_directories(tests PRIVATE
        ${UNIFRAC_SRC_DIR}/src
        ${SKBB_SRC_DIR}/src/extern)
    target_compile_definitions(tests PRIVATE MIINT_HAS_UNIFRAC)
    # Tests call into libssu directly for the BP round-trip check
    # (load_bptree_opaque, convert_bptree_opaque, get_bptree_opaque_els).
    # libssu depends on libskbb symbols, and the native archive uses OpenMP.
    target_link_libraries(tests libssu libskbb)
    if(NOT EMSCRIPTEN)
        target_link_libraries(tests miint_openmp)
    endif()
endif()
target_include_directories(tests PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}/test/cpp/include  # Catch2 helper headers
    ${CMAKE_CURRENT_BINARY_DIR}/htslib/include
    ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/miniz  # For ExtractFromZip tests
)
if(MIINT_HAS_JEMALLOC)
    target_sources(tests PRIVATE src/mm_alloc_stubs.c)
    target_compile_definitions(tests PRIVATE MIINT_USE_JEMALLOC)
endif()
target_compile_definitions(tests PRIVATE HAVE_KALLOC)
