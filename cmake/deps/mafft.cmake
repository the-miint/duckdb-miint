# MAFFT PartTree static library for multiple sequence alignment.
#
# Off on Windows: MAFFT uses mkdtemp and other POSIX APIs, and segfaults on
# MinGW. Off under Emscripten: its local arrays are sized by M (mltaln.h:
# 500000, the maximum sequence count), giving 2-4 MB stack frames, and
# duckdb-wasm runs with a 1 MB stack (-sSTACK_SIZE=1048576) that has no guard
# page. Entering such a frame writes past the stack into the heap, and the
# damage surfaces later as "memory access out of bounds" in free() --
# align_mafft on the 36-sequence data/mafft/sample.fa. The frames are MAFFT's
# own (M is not overridable) and the stack is the runtime's, so the fix would be
# MAFFT-side.

# Keep -fPIC on every platform: required for WASM side module linking
# (-sSIDE_MODULE=2) too.
set(MAFFT_FPIC "-fPIC")
set(MAFFT_ARCH_CFLAGS "${MIINT_OSX_ARCH_CFLAGS} ${EMSCRIPTEN_THREAD_CFLAGS}")

ExternalProject_Add(
  mafft_build
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core clean
  COMMAND
    make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core -j${MIINT_SUBBUILD_JOBS}
    libmafft_parttree.a ${MIINT_SUBBUILD_COMPILER_VARS}
    "CFLAGS=-O3 ${MAFFT_FPIC} ${MAFFT_ARCH_CFLAGS}"
    "ENABLE_MULTITHREAD=-Denablemultithread"
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/mafft/core/libmafft_parttree.a
    ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a)

add_library(mafft_parttree STATIC IMPORTED GLOBAL)
set_target_properties(
  mafft_parttree
  PROPERTIES IMPORTED_LOCATION
             ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a)
add_dependencies(mafft_parttree mafft_build)

miint_git_version(MAFFT ext/mafft)
list(APPEND MIINT_DEFINITIONS MIINT_HAS_MAFFT)
list(APPEND MIINT_INCLUDE_DIRS ext/mafft/core)
list(APPEND MIINT_EXTENSION_SOURCES src/MafftAligner.cpp src/align_mafft.cpp)
list(APPEND MIINT_TEST_SOURCES src/MafftAligner.cpp
     test/cpp/test_MafftAligner.cpp)
list(APPEND MIINT_LINK_LIBRARIES mafft_parttree)
list(APPEND MIINT_WASM_LIBS
     ${CMAKE_CURRENT_BINARY_DIR}/mafft/lib/libmafft_parttree.a)
