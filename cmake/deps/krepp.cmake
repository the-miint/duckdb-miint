# krepp phylogenetic placement (commit 31205033; VERSION string v0.10.0,
# untagged). Header-only dependencies come from vcpkg rather than krepp's own
# submodules. ext/krepp/src/krepp.cpp is omitted: it holds main() and the CLI11
# argument layer. KreppPlacer therefore reimplements krepp's index discovery,
# which lives in that same file and will drift if krepp changes its index
# layout. Off on Windows and under Emscripten.
#
# Like the HDF5 and SortMeRNA deps files: clang-tidy builds do not have vcpkg
# wired up, so these headers can be absent there. Dropping krepp from the
# analysis set is the house response - a FATAL_ERROR here fails the whole tidy
# pipeline for every other extension too, which is what it did on the first CI
# run of this feature. Regular builds still require both headers.
find_path(PHMAP_INCLUDE_DIR parallel_hashmap/phmap.h)
find_path(BOOST_MATH_INCLUDE_DIR boost/math/tools/minima.hpp)
if(CLANG_TIDY)
  if(NOT PHMAP_INCLUDE_DIR OR NOT BOOST_MATH_INCLUDE_DIR)
    message(
      STATUS
        "parallel-hashmap/boost-math not found for tidy-check - krepp placement support will be "
        "excluded from analysis")
    set(MIINT_ENABLE_KREPP OFF)
    return()
  endif()
  message(
    STATUS
      "parallel-hashmap and boost-math found for tidy-check - krepp placement support enabled"
  )
elseif(NOT PHMAP_INCLUDE_DIR OR NOT BOOST_MATH_INCLUDE_DIR)
  message(
    FATAL_ERROR
      "krepp requires parallel-hashmap and boost-math headers (vcpkg). "
      "Configure with -DMIINT_ENABLE_KREPP=OFF to build without placement support."
  )
endif()

set(KREPP_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/krepp)
add_library(
  krepp STATIC
  ${KREPP_DIR}/src/common.cpp
  ${KREPP_DIR}/src/MurmurHash3.cpp
  ${KREPP_DIR}/src/lshf.cpp
  ${KREPP_DIR}/src/phytree.cpp
  ${KREPP_DIR}/src/rqseq.cpp
  ${KREPP_DIR}/src/index.cpp
  ${KREPP_DIR}/src/sketch.cpp
  ${KREPP_DIR}/src/query.cpp
  ${KREPP_DIR}/src/seek.cpp
  ${KREPP_DIR}/src/record.cpp
  ${KREPP_DIR}/src/table.cpp)
set_target_properties(krepp PROPERTIES POSITION_INDEPENDENT_CODE ON)
# krepp and minimap2 both define a global `seq_nt4_table`, and they are not the
# same table: minimap2's is `unsigned char [256]`, krepp's is `const unsigned
# char [128]`. Linking both is a duplicate-symbol error, and letting one win
# would be worse than the error - minimap2 indexes bytes above 127, which would
# read off the end of krepp's array. Renaming krepp's copy keeps the two apart.
# This line covers krepp's own sources; the single miint TU that includes
# krepp's headers gets the same define via a source-file property below -
# deliberately not target-wide, for the reason recorded there. PRIVATE, not
# PUBLIC: krepp's own TUs still get both, but neither becomes a usage
# requirement. The three link sites use $<TARGET_FILE:krepp>, which carries no
# usage requirements, so PUBLIC is inert today - but one target_link_libraries(x
# krepp) anywhere would re-arm exactly the trap the comment above describes, and
# put krepp's kseq.h/kvec.h ahead of minimap2's for that target. krepp gates
# every OpenMP region on `defined(_OPENMP) && _WOPENMP == 1`, so the define is
# only half of it: the compile flag has to actually turn OpenMP on. Verified
# rather than assumed, because if it silently did not, krepp_index_create would
# advertise the feature, accept a thread count, and still build on one core.
set(MIINT_KREPP_OPENMP_OK OFF)
if(MIINT_OPENMP_AVAILABLE)
  if(APPLE)
    set(MIINT_KREPP_OPENMP_FLAG -Xpreprocessor -fopenmp)
  else()
    set(MIINT_KREPP_OPENMP_FLAG ${OpenMP_CXX_FLAGS})
  endif()
  include(CheckCXXSourceCompiles)
  set(MIINT_SAVED_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS}")
  string(REPLACE ";" " " CMAKE_REQUIRED_FLAGS "${MIINT_KREPP_OPENMP_FLAG}")
  check_cxx_source_compiles(
    "#ifndef _OPENMP\n#error OpenMP is not enabled\n#endif\nint main() { return 0; }"
    MIINT_KREPP_OPENMP_FLAG_DEFINES_OPENMP)
  set(CMAKE_REQUIRED_FLAGS "${MIINT_SAVED_REQUIRED_FLAGS}")
  if(MIINT_KREPP_OPENMP_FLAG_DEFINES_OPENMP)
    set(MIINT_KREPP_OPENMP_OK ON)
  else()
    message(
      STATUS
        "OpenMP flag does not define _OPENMP - krepp index builds stay single-threaded."
    )
  endif()
endif()
# _WOPENMP feeds both the archive and the two miint TUs that include krepp's
# headers. In those headers it gates only `#include \"omp.h\"` - none of them
# declares an omp type or holds a pragma - so a mismatch would be a preprocessor
# inconsistency, not an ABI one. One variable regardless, since the same value
# decides MIINT_KREPP_OPENMP below.
if(MIINT_KREPP_OPENMP_OK)
  set(MIINT_KREPP_WOPENMP 1)
else()
  set(MIINT_KREPP_WOPENMP 0)
endif()
# krepp's src/sdust.h, compiled into rqseq.cpp, is a copy of minimap2's sdust
# with the same global names, so a link that pulls in both archives' copies
# (tests_duckdb's does) fails with duplicate symbols. No miint TU includes
# sdust.h, so unlike seq_nt4_table the rename needs no source-file property. See
# docs/internals/embedded-tools.md, krepp.
target_compile_definitions(
  krepp
  PRIVATE _WOPENMP=${MIINT_KREPP_WOPENMP}
          _WLCURL=0
          seq_nt4_table=krepp_seq_nt4_table
          sdust=krepp_sdust
          sdust_core=krepp_sdust_core
          sdust_buf_init=krepp_sdust_buf_init
          sdust_buf_destroy=krepp_sdust_buf_destroy
          seq2nt4=krepp_seq2nt4)
if(MIINT_KREPP_OPENMP_OK)
  # Apple Clang rejects a bare -fopenmp and wants the preprocessor form, with
  # libomp linked separately - the same choice the skbb/unifrac sub-makes make
  # for INMEM_MPFLAG in unifrac.cmake. Off Darwin the flag arrives with
  # OpenMP::OpenMP_CXX's interface compile options through the link below, so it
  # is not repeated here: FindOpenMP spells some candidates with a space
  # ("-Xclang -fopenmp"), and target_compile_options would pass that as a single
  # shell-quoted argument.
  if(APPLE)
    target_compile_options(krepp PRIVATE ${MIINT_KREPP_OPENMP_FLAG})
  endif()
  # For omp.h on the include path. The runtime itself is linked at each
  # extension target, because those link krepp by archive path rather than by
  # target name and so inherit nothing from here.
  target_link_libraries(krepp PRIVATE miint_openmp)
endif()
target_include_directories(krepp PRIVATE ${KREPP_DIR}/src ${PHMAP_INCLUDE_DIR}
                                         ${BOOST_MATH_INCLUDE_DIR})
# PUBLIC here, unlike the two above, is deliberate: zlib is a real link
# dependency any consumer of this archive needs, whereas the rename and the klib
# include dir are things a consumer must NOT inherit.
target_link_libraries(krepp PUBLIC ZLIB::ZLIB)

# Scoped to the same single TU as the defines below, and for the same reason.
# ${KREPP_DIR}/src ships a klib set - kseq.h, kvec.h, kalloc.h, kdq.h, ketopt.h,
# sdust.h - whose basenames collide with minimap2's. On the directory-level
# include path those would shadow for every target in this file. Nothing
# includes them bare today, so this is latent rather than a live bug - but the
# seq_nt4_table breakage below was latent in exactly the same way until a header
# happened to line up, so scope it now.
set_source_files_properties(
  src/KreppPlacer.cpp src/KreppIndexBuilder.cpp
  PROPERTIES INCLUDE_DIRECTORIES
             "${KREPP_DIR}/src;${PHMAP_INCLUDE_DIR};${BOOST_MATH_INCLUDE_DIR}")
# These three belong to the ONE translation unit that includes krepp's headers,
# not to the targets. Applied target-wide they also rewrote the declaration in
# src/KSW2Aligner.cpp, which declares minimap2's `unsigned char
# seq_nt4_table[256]` and indexes it over the full byte range - so it bound
# krepp's `const unsigned char [128]` instead and read up to 128 bytes past the
# end for any input byte >= 128. That is precisely the out-of-bounds the rename
# exists to prevent, pointed the wrong way. Source-file properties are
# directory-scoped, so this reaches miint_extension, miint_loadable_extension
# and tests alike.
set_source_files_properties(
  src/KreppPlacer.cpp src/KreppIndexBuilder.cpp
  PROPERTIES
    COMPILE_DEFINITIONS
    "_WOPENMP=${MIINT_KREPP_WOPENMP};_WLCURL=0;seq_nt4_table=krepp_seq_nt4_table;MIINT_KREPP_OPENMP=${MIINT_KREPP_WOPENMP}"
)

miint_git_version(KREPP ext/krepp)
list(APPEND MIINT_DEFINITIONS MIINT_HAS_KREPP)
list(APPEND MIINT_EXTENSION_SOURCES src/KreppPlacer.cpp src/place_krepp.cpp
     src/KreppIndexBuilder.cpp src/krepp_index_create.cpp)
list(APPEND MIINT_TEST_SOURCES src/KreppPlacer.cpp
     test/cpp/test_KreppPlacer.cpp src/KreppIndexBuilder.cpp
     test/cpp/test_KreppIndexBuilder.cpp)
# Link the archive by path, not by target name: miint_extension goes into
# DuckDB's install(EXPORT) set, which rejects link dependencies on targets
# outside it. Same treatment as abpoa.
list(APPEND MIINT_LINK_LIBRARIES $<TARGET_FILE:krepp>)
list(APPEND MIINT_TARGET_DEPENDENCIES krepp)
# The archive's OpenMP runtime: linked here only when UniFrac does not already
# link it.
if(MIINT_KREPP_OPENMP_OK AND NOT MIINT_ENABLE_UNIFRAC)
  list(APPEND MIINT_LINK_LIBRARIES miint_openmp)
endif()
