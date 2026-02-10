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

# HDF5 is required for BIOM support, but make it optional for clang-tidy runs
# In CI, the tidy-check doesn't have HDF5 available, but it doesn't need to link - just parse
if(CLANG_TIDY)
    find_package(HDF5 COMPONENTS CXX)
    if(HDF5_FOUND)
        message(STATUS "HDF5 found for tidy-check")
    else()
        message(STATUS "HDF5 not found for tidy-check - BIOM support will be excluded from analysis")
        # Create dummy targets so the build can configure
        add_library(HDF5::HDF5 INTERFACE IMPORTED)
    endif()
else()
    find_package(HDF5 REQUIRED COMPONENTS CXX)
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
if(APPLE AND OSX_BUILD_ARCH)
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

set(HTSLIB_CFLAGS "CFLAGS=-fPIC ${HTSLIB_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}")
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

ExternalProject_Add(
    htslib_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1
    CONFIGURE_COMMAND env ${HTSLIB_CFLAGS} ${HTSLIB_LDFLAGS} ${HTSLIB_LIBS} ${HTSLIB_PKG_CONFIG_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1/configure
        ${HTSLIB_HOST_FLAG}
        --prefix=${CMAKE_CURRENT_BINARY_DIR}/htslib
        --disable-libcurl
        --disable-bz2
        --disable-lzma
    BUILD_COMMAND make
    INSTALL_COMMAND make install
    BUILD_IN_SOURCE 0
    BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/htslib/lib/libhts.a
)

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
if(APPLE AND OSX_BUILD_ARCH)
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

# Build minimap2 from source (v2.30)
ExternalProject_Add(
    minimap2_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2
    CONFIGURE_COMMAND ""
    BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 libminimap2.a "CFLAGS=-g -Wall -O2 -Wc++-compat -fPIC ${MINIMAP2_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}" "CPPFLAGS=-DHAVE_KALLOC" ${MINIMAP2_ARCH_FLAGS}
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
if(APPLE AND OSX_BUILD_ARCH)
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
        "CC_FLAGS=-Wall -O2 -fPIC ${WFA2_ARCH_CFLAGS}"
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

# adapted kseq++ include use from strobealign
include_directories(src/include ext ext/WFA2-lib)
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
    src/copy_biom.cpp
    src/copy_fastq.cpp
    src/copy_fasta.cpp
    src/copy_sam.cpp
    src/BIOMTable.cpp
    src/BIOMReader.cpp
    src/read_biom.cpp
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
    src/Bowtie2Aligner.cpp
    src/align_bowtie2.cpp
    src/align_bowtie2_sharded.cpp
    src/ncbi_parser.cpp
    src/ncbi_client.cpp
    src/zip_utils.cpp
    src/read_ncbi_fasta.cpp
    src/read_ncbi.cpp
    src/read_ncbi_annotation.cpp
    src/SFFReader.cpp
    src/read_sequences_sff.cpp
    src/WFA2Aligner.cpp
    src/align_pairwise_functions.cpp)

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

target_link_libraries(${EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts
                      minimap2
                      wfa2cpp
                      wfa2
)
if(HDF5_FOUND)
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

target_link_libraries(${LOADABLE_EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts
                      minimap2
                      wfa2cpp
                      wfa2
)

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
    src/BIOMTable.cpp
    src/BIOMReader.cpp
    test/cpp/test_BIOMReader.cpp
    test/cpp/test_BIOMTable.cpp
    src/NewickTree.cpp
    test/cpp/test_NewickParser.cpp
    test/cpp/test_InsertFullyResolved.cpp
    src/Minimap2Aligner.cpp
    test/cpp/test_Minimap2Aligner.cpp
    src/Bowtie2Aligner.cpp
    test/cpp/test_Bowtie2Aligner.cpp
    src/ncbi_parser.cpp
    src/zip_utils.cpp
    duckdb/third_party/miniz/miniz.cpp  # For ExtractFromZip tests
    test/cpp/test_ncbi_client.cpp
    src/SFFReader.cpp
    test/cpp/test_SFFReader.cpp
    src/WFA2Aligner.cpp
    test/cpp/test_WFA2Aligner.cpp
)

add_executable(tests ${TEST_SOURCES})

target_link_libraries(tests
  Catch2::Catch2WithMain
  Threads::Threads
  ZLIB::ZLIB
  hts
  minimap2
  wfa2cpp
  wfa2
)
if(HDF5_FOUND)
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
