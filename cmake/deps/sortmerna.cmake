# SortMeRNA: embedded library for rRNA filtering and alignment. Upstream CMake
# produces several static archives; smr_api is the C-API entry point and holds
# the smr_objs objects (rapidgzip, BBHash and cxxopts are header-only and
# compiled into them). build_version and alp are separate archives, and the
# parasail SIMD aligner is built here from its submodule under
# ext/sortmerna/3rdparty.
#
# Upstream also builds a self-contained libsmr_api.a under src/smr_api/bundle/
# that merges in RocksDB and zlib. We do not link it: miint links its own vcpkg
# RocksDB and zlib, and the bundled copies would define those symbols twice.
# rapidgzip's zlibstatic is left out for the same reason, so its zlib calls
# resolve to the one vcpkg zlib.
#
# Off on MinGW (untested there: RocksDB and parasail have never been built for
# it) and under Emscripten (the vcpkg RocksDB port is not built for wasm32, and
# the alignment pipeline always starts worker threads).
#
# Like the HDF5 deps file: clang-tidy builds don't have vcpkg wired up, so
# RocksDB (the only hard external dependency) may be absent. In that case we
# silently drop sortmerna from the analysis set rather than fail the whole tidy
# pipeline. Regular builds still REQUIRE RocksDB.
if(CLANG_TIDY)
  find_package(RocksDB CONFIG QUIET)
  if(NOT RocksDB_FOUND)
    message(
      STATUS
        "RocksDB not found for tidy-check - SortMeRNA support will be excluded from analysis"
    )
    set(MIINT_ENABLE_SORTMERNA OFF)
    return()
  endif()
  message(STATUS "RocksDB found for tidy-check - SortMeRNA support enabled")
else()
  find_package(RocksDB CONFIG REQUIRED)
endif()

set(SORTMERNA_SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/sortmerna)
set(SORTMERNA_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/sortmerna_build)
set(SORTMERNA_API_LIB ${SORTMERNA_BIN_DIR}/src/smr_api/libsmr_api.a)
set(SORTMERNA_VERSION_LIB ${SORTMERNA_BIN_DIR}/libbuild_version.a)
set(SORTMERNA_ALP_LIB ${SORTMERNA_BIN_DIR}/3rdparty/alp/libalp.a)
set(PARASAIL_SRC_DIR ${SORTMERNA_SRC_DIR}/3rdparty/parasail)
set(PARASAIL_BIN_DIR ${CMAKE_CURRENT_BINARY_DIR}/parasail_build)
set(PARASAIL_LIB_FILE ${PARASAIL_BIN_DIR}/libparasail.a)

# macOS: forward OSX_BUILD_ARCH as CMAKE_OSX_ARCHITECTURES, or the sub-CMake
# links sortmerna's CLI target for its default arch and trips on vcpkg's
# RocksDB, which is built for the target arch.
set(SORTMERNA_OSX_ARCH_ARG "")
if(APPLE AND OSX_BUILD_ARCH)
  set(SORTMERNA_OSX_ARCH_ARG "-DCMAKE_OSX_ARCHITECTURES=${OSX_BUILD_ARCH}")
endif()

# macOS: the smr_api bundle step (an ALL target, though we don't link its
# output) merges archives with Apple's `libtool -static`, found by name on PATH.
# A GNU libtool earlier on PATH (conda, which also supplies autoreconf for
# vsearch) has no -static and fails the build, so name Apple's explicitly.
set(SORTMERNA_LIBTOOL_ARG "")
if(APPLE)
  set(SORTMERNA_LIBTOOL_ARG "-DSMR_API_LIBTOOL=/usr/bin/libtool")
endif()

# Point both sub-builds' FindZLIB at the vcpkg zlib miint links (as htslib.cmake
# does). Left to search, they find the SDK/system zlib -- a forwarded vcpkg
# toolchain has no vcpkg_installed tree in the sub-build -- and v8's smr_api
# bundle step then fails configure on a zlib that is not a static archive.
# ZLIB_LIBRARY_RELEASE is what ZLIB::ZLIB's location comes from, and a build
# tree configured before this change has the system zlib cached there.
# system_deps.cmake warns and leaves ZLIB_LIBRARY_FILE empty when it cannot name
# the archive; the sub-builds then search on their own rather than cache "".
set(SORTMERNA_ZLIB_ARGS "")
if(ZLIB_LIBRARY_FILE)
  set(SORTMERNA_ZLIB_ARGS
      -DZLIB_LIBRARY=${ZLIB_LIBRARY_FILE}
      -DZLIB_LIBRARY_RELEASE=${ZLIB_LIBRARY_FILE}
      -DZLIB_INCLUDE_DIR=${ZLIB_INCLUDE_DIR})
endif()

# Only the library target: parasail's default targets also build its aligner app
# (OpenMP) and test programs. It requires CMake 3.2, which CMake 4 refuses
# without a policy floor. Its SIMD variants are chosen at runtime, so one build
# serves every CPU of the target architecture.
ExternalProject_Add(
  parasail_build
  SOURCE_DIR ${PARASAIL_SRC_DIR}
  BINARY_DIR ${PARASAIL_BIN_DIR}
  CMAKE_ARGS -DCMAKE_BUILD_TYPE=Release
             -DBUILD_SHARED_LIBS=OFF
             -DCMAKE_POSITION_INDEPENDENT_CODE=ON
             -DWITH_OPENMP=OFF
             -DCMAKE_POLICY_VERSION_MINIMUM=3.5
             ${SORTMERNA_ZLIB_ARGS}
             ${SORTMERNA_OSX_ARCH_ARG}
  BUILD_COMMAND ${CMAKE_COMMAND} --build <BINARY_DIR> --target parasail
  BUILD_BYPRODUCTS ${PARASAIL_LIB_FILE}
  INSTALL_COMMAND "")

# alp's sources use the deprecated `register` keyword, which is a hard error on
# Clang 16+. Suppress it so the library still compiles. The other -Wno-* flags
# silence cosmetic warnings from upstream sortmerna/alp code (type-punned
# pointer casts, intentional switch fall-through without [[fallthrough]],
# unchecked fread return values). All would require upstream patches to fix
# properly.
set(SORTMERNA_WARNING_SUPPRESSIONS
    "-Wno-register -Wno-strict-aliasing -Wno-implicit-fallthrough -Wno-unused-result"
)
# All targets are built, not just smr_api: run_tests.sh regenerates the
# real-data oracle with the native sortmerna CLI from this build tree.
ExternalProject_Add(
  sortmerna_build
  SOURCE_DIR ${SORTMERNA_SRC_DIR}
  BINARY_DIR ${SORTMERNA_BIN_DIR}
  DEPENDS parasail_build
  CMAKE_ARGS
    -DCMAKE_BUILD_TYPE=Release
    -DCMAKE_POSITION_INDEPENDENT_CODE=ON
    -DCMAKE_CXX_FLAGS=${SORTMERNA_WARNING_SUPPRESSIONS}
    -DCMAKE_C_FLAGS=${SORTMERNA_WARNING_SUPPRESSIONS}
    -DWITH_TESTS=OFF
    -DCMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}
    -DVCPKG_TARGET_TRIPLET=${VCPKG_TARGET_TRIPLET}
    -DVCPKG_MANIFEST_MODE=OFF
    ${SORTMERNA_ZLIB_ARGS}
    ${SORTMERNA_OSX_ARCH_ARG}
    ${SORTMERNA_LIBTOOL_ARG}
    # Force sortmerna's FindRocksDB branch to use vcpkg's CONFIG file. Without
    # this, its else-branch runs MODULE mode which hits system /usr/lib rocksdb
    # first on Debian/Ubuntu.
    -DROCKSDB_DIST=${CMAKE_BINARY_DIR}/vcpkg_installed/${VCPKG_TARGET_TRIPLET}/share/rocksdb
    # vcpkg builds static-only; switch sortmerna off its default shared target.
    -DROCKSDB_USE_STATIC_LIBS=ON
    # Pre-seeding the cache entries skips sortmerna's find_library/find_path,
    # which would otherwise look in a setup.py install prefix.
    -DPARASAIL_LIB=${PARASAIL_LIB_FILE}
    -DPARASAIL_INCLUDE_DIR=${PARASAIL_SRC_DIR}
  BUILD_BYPRODUCTS ${SORTMERNA_API_LIB} ${SORTMERNA_VERSION_LIB}
                   ${SORTMERNA_ALP_LIB}
  INSTALL_COMMAND "")

# Imported archives in link order: smr_api calls into build_version, alp and
# parasail.
foreach(lib api version alp)
  string(TOUPPER ${lib} lib_upper)
  add_library(sortmerna_${lib} STATIC IMPORTED GLOBAL)
  set_target_properties(
    sortmerna_${lib} PROPERTIES IMPORTED_LOCATION ${SORTMERNA_${lib_upper}_LIB})
  add_dependencies(sortmerna_${lib} sortmerna_build)
endforeach()
set_target_properties(sortmerna_api PROPERTIES INTERFACE_INCLUDE_DIRECTORIES
                                               ${SORTMERNA_SRC_DIR}/include)

add_library(sortmerna_parasail STATIC IMPORTED GLOBAL)
set_target_properties(sortmerna_parasail PROPERTIES IMPORTED_LOCATION
                                                    ${PARASAIL_LIB_FILE})
add_dependencies(sortmerna_parasail parasail_build)

list(APPEND MIINT_DEFINITIONS MIINT_HAS_SORTMERNA)
list(APPEND MIINT_EXTENSION_SOURCES src/SortMeRNAAligner.cpp
     src/align_sortmerna.cpp src/align_sortmerna_rrna.cpp)
list(APPEND MIINT_TEST_SOURCES src/SortMeRNAAligner.cpp
     test/cpp/test_sortmerna_build.cpp test/cpp/test_SortMeRNAAligner.cpp)
list(
  APPEND
  MIINT_LINK_LIBRARIES
  sortmerna_api
  sortmerna_version
  sortmerna_alp
  sortmerna_parasail
  RocksDB::rocksdb)
