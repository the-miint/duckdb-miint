# The gpl-boundary process host: miint talks to GPL tools (bowtie2, FastTree)
# over pipes and Arrow IPC in POSIX shared memory, so no GPL code is linked in.
# Off on Windows and under Emscripten. skb_arrow.cmake, included at the end,
# reuses this transport.
add_subdirectory(third_party/nanoarrow_ipc)

list(APPEND MIINT_DEFINITIONS MIINT_HAS_GPL_BOUNDARY)
# nanoarrow_ipc is an OBJECT library: its .o's go into each target without a
# link dependency, which would require nanoarrow_ipc to also be in the DuckDB
# export set, which it has no business being in.
list(APPEND MIINT_OBJECT_SOURCES $<TARGET_OBJECTS:miint_nanoarrow_ipc>)
# yyjson is bundled in DuckDB as a static lib (target `duckdb_yyjson`) at
# `duckdb/third_party/yyjson/CMakeLists.txt`. We just consume its include
# directories; the extension gets the symbols transitively via the duckdb
# library.
list(APPEND MIINT_INCLUDE_DIRS
     $<TARGET_PROPERTY:miint_nanoarrow_ipc,INTERFACE_INCLUDE_DIRECTORIES>
     ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/include)
list(
  APPEND
  MIINT_EXTENSION_SOURCES
  src/gpl_boundary/process.cpp
  src/gpl_boundary/session.cpp
  src/gpl_boundary/shm.cpp
  src/gpl_boundary/arrow_ipc.cpp
  src/install_gpl_boundary.cpp
  src/gpl_boundary_version.cpp
  src/phylogeny_fasttree.cpp
  src/predict_genes.cpp
  src/align_bowtie2_daemon_common.cpp
  src/align_bowtie2.cpp
  src/align_bowtie2_sharded.cpp
  src/save_bowtie2_index.cpp)
# The catch2 tests binary doesn't link against libduckdb (it's a standalone
# executable), so it needs its own copy of yyjson.cpp.
list(
  APPEND
  MIINT_TEST_SOURCES
  src/gpl_boundary/process.cpp
  src/gpl_boundary/session.cpp
  src/gpl_boundary/shm.cpp
  src/gpl_boundary/arrow_ipc.cpp
  test/cpp/test_gpl_boundary_process.cpp
  test/cpp/test_gpl_boundary_session.cpp
  test/cpp/test_gpl_boundary_shm.cpp
  test/cpp/test_gpl_boundary_arrow_ipc.cpp
  ${CMAKE_CURRENT_SOURCE_DIR}/duckdb/third_party/yyjson/yyjson.cpp)
if(UNIX AND NOT APPLE)
  # gpl_boundary/shm.cpp uses shm_open/shm_unlink. On glibc these live in librt
  # (older RHEL/CentOS toolchains require explicit linkage; newer glibc folded
  # them into libc but we link defensively). macOS exposes them in libSystem so
  # we skip the flag there.
  list(APPEND MIINT_LINK_LIBRARIES rt)
endif()

include(${CMAKE_CURRENT_LIST_DIR}/skb_arrow.cmake)
