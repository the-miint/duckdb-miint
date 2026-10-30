# minimap2 (v2.31): read alignment and indexing. Always built. Its ksw2 kernels
# also back the align_pairwise_ksw2* functions.

# Detect architecture for minimap2 SIMD support: OSX_BUILD_ARCH on macOS,
# otherwise the host architecture
set(MINIMAP2_ARCH_FLAGS "")
set(MINIMAP2_ARCH_CFLAGS "")
set(MINIMAP2_FPIC "-fPIC")
if(EMSCRIPTEN)
  # Emscripten: use WASM SIMD via Emscripten's SSE2 emulation. sse2only=1 avoids
  # SSE4.1 dispatch which requires runtime CPU detection.
  set(MINIMAP2_ARCH_FLAGS "sse2only=1")
  set(MINIMAP2_ARCH_CFLAGS "-msimd128 -msse -msse2 ${EMSCRIPTEN_THREAD_CFLAGS}")
  # Keep -fPIC: required for WASM side module linking (-sSIDE_MODULE=2)
  message(
    STATUS "Minimap2 building for Emscripten/WASM with SIMD via SSE2 emulation")
elseif(APPLE AND OSX_BUILD_ARCH)
  # Cross-compilation: use target architecture, not host
  set(MINIMAP2_ARCH_CFLAGS "${MIINT_OSX_ARCH_CFLAGS}")
  if(OSX_BUILD_ARCH STREQUAL "arm64")
    set(MINIMAP2_ARCH_FLAGS "aarch64=1")
  endif()
  message(STATUS "Minimap2 cross-compiling (${MIINT_OSX_ARCH_CFLAGS})")
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
  # Native ARM64 build
  set(MINIMAP2_ARCH_FLAGS "aarch64=1")
endif()

# Route minimap2's allocations through DuckDB's jemalloc when DuckDB builds it
# (64-bit Linux; not musl, macOS, Windows or WASM). DuckDB makes that choice
# after it adds the extensions, so no variable can tell us here; TARGET_EXISTS
# is evaluated at generate time, when duckdb_jemalloc exists or does not.
# CMakeLists.txt adds src/mm_alloc_stubs.c under the same condition.
set(MIINT_JEMALLOC "$<TARGET_EXISTS:duckdb_jemalloc>")
set(MINIMAP2_ALLOC_FLAGS
    "$<${MIINT_JEMALLOC}:-DMIINT_USE_JEMALLOC -include ${CMAKE_CURRENT_SOURCE_DIR}/src/include/mm_alloc.h>"
)
list(APPEND MIINT_DEFINITIONS "$<${MIINT_JEMALLOC}:MIINT_USE_JEMALLOC>")

ExternalProject_Add(
  minimap2_build
  SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2
  CONFIGURE_COMMAND ""
  BUILD_COMMAND make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 clean
  COMMAND
    make -C ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2 libminimap2.a
    "CFLAGS=-g -Wall -O2 -Wc++-compat ${MINIMAP2_FPIC} ${MINIMAP2_ARCH_CFLAGS} -I${ZLIB_INCLUDE_DIRS}"
    "CPPFLAGS=-DHAVE_KALLOC ${MINIMAP2_ALLOC_FLAGS}" ${MINIMAP2_ARCH_FLAGS}
  INSTALL_COMMAND
    ${CMAKE_COMMAND} -E copy
    ${CMAKE_CURRENT_SOURCE_DIR}/ext/minimap2/libminimap2.a
    ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a
  BUILD_IN_SOURCE 1
  BUILD_BYPRODUCTS ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a)

add_library(minimap2 STATIC IMPORTED GLOBAL)
set_target_properties(
  minimap2 PROPERTIES IMPORTED_LOCATION
                      ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a)
add_dependencies(minimap2 minimap2_build)

miint_git_version(MINIMAP2 ext/minimap2)
# ksw2.h routes kmalloc/kfree through kalloc only under HAVE_KALLOC (otherwise
# to plain malloc), so miint's TUs must agree with the archive's CPPFLAGS above.
list(APPEND MIINT_DEFINITIONS HAVE_KALLOC)
list(APPEND MIINT_LINK_LIBRARIES minimap2)
list(APPEND MIINT_WASM_LIBS
     ${CMAKE_CURRENT_BINARY_DIR}/minimap2/lib/libminimap2.a)
