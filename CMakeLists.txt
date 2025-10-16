cmake_minimum_required(VERSION 3.10)
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

ExternalProject_Add(
    htslib_build
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1
    CONFIGURE_COMMAND env CFLAGS=-fPIC ${HTSLIB_PKG_CONFIG_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/ext/htslib-1.22.1/configure
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

# adapted kseq++ include use from strobealign
include_directories(src/include ext)
set(EXTENSION_SOURCES
    src/SequenceReader.cpp
    src/QualScore.cpp
    src/read_fastx.cpp
    src/miint_extension.cpp
    src/QualScore.cpp
    src/read_sam.cpp
    src/SAMReader.cpp
    src/alignment_flag_functions.cpp
    src/alignment_functions.cpp
    src/IntervalCompressor.cpp
    src/compress_intervals.cpp
    src/copy_format_common.cpp
    src/copy_fastq.cpp
    src/copy_fasta.cpp
    src/copy_sam.cpp
    src/reference_table_reader.cpp)

build_static_extension(${TARGET_NAME} ${EXTENSION_SOURCES})
build_loadable_extension(${TARGET_NAME} " " ${EXTENSION_SOURCES})


target_link_libraries(${EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts)
if(HAVE_LIBZSTD)
    target_link_libraries(${EXTENSION_NAME} ${ZSTD_TARGET})
endif()

target_link_libraries(${LOADABLE_EXTENSION_NAME}
                      Threads::Threads
                      ZLIB::ZLIB
                      hts)
if(HAVE_LIBZSTD)
    target_link_libraries(${LOADABLE_EXTENSION_NAME} ${ZSTD_TARGET})
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
)

add_executable(tests ${TEST_SOURCES})

target_link_libraries(tests
  Catch2::Catch2WithMain
  Threads::Threads
  ZLIB::ZLIB
  hts
)
if(HAVE_LIBZSTD)
    target_link_libraries(tests ${ZSTD_TARGET})
endif()
target_include_directories(tests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include ${CMAKE_CURRENT_BINARY_DIR}/htslib/include)
