# Helpers shared by CMakeLists.txt and cmake/deps/*.cmake. The registration
# lists they fill are described in docs/internals/build-system.md.

# Record `git describe` of the submodule at <dir> as <name>_GIT_VERSION="..." in
# MIINT_DEFINITIONS, for miint_versions() to report. "unknown" outside a git
# checkout (e.g. a source tarball).
function(miint_git_version name dir)
  execute_process(
    COMMAND git describe --tags --dirty --always
    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/${dir}
    OUTPUT_VARIABLE version
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT version)
    set(version "unknown")
  endif()
  list(APPEND MIINT_DEFINITIONS "${name}_GIT_VERSION=\"${version}\"")
  list(APPEND MIINT_SUBMODULE_VERSIONS "${name}=${version}")
  set(MIINT_DEFINITIONS
      "${MIINT_DEFINITIONS}"
      PARENT_SCOPE)
  set(MIINT_SUBMODULE_VERSIONS
      "${MIINT_SUBMODULE_VERSIONS}"
      PARENT_SCOPE)
endfunction()

# Apply the registration lists to one of the three targets: the static
# extension, the loadable extension, or `tests`. Called once per target, so the
# three cannot drift apart; what is specific to one target is added next to its
# call in CMakeLists.txt.
#
# Uses the plain target_link_libraries signature (no PRIVATE/PUBLIC): DuckDB's
# build_*_extension macros use it on the extension targets, and CMake rejects
# mixing the two signatures on one target.
function(miint_configure_target target)
  target_sources(${target} PRIVATE ${MIINT_OBJECT_SOURCES})
  target_compile_definitions(${target} PRIVATE ${MIINT_DEFINITIONS})
  target_link_libraries(${target} ${MIINT_LINK_LIBRARIES})
  if(MIINT_TARGET_DEPENDENCIES)
    add_dependencies(${target} ${MIINT_TARGET_DEPENDENCIES})
  endif()
endfunction()

# Cross-compilation flags for the C libraries built with their own Makefiles or
# autotools, which do not see CMAKE_C_FLAGS. At most one platform branch
# applies.
set(MIINT_OSX_ARCH_CFLAGS "") # -arch for an osx cross build
set(MIINT_AUTOTOOLS_HOST "") # --host= for an autotools configure
if(EMSCRIPTEN)
  set(MIINT_AUTOTOOLS_HOST "--host=wasm32-unknown-emscripten")
elseif(APPLE AND OSX_BUILD_ARCH STREQUAL "arm64")
  set(MIINT_OSX_ARCH_CFLAGS "-arch arm64")
  set(MIINT_AUTOTOOLS_HOST "--host=aarch64-apple-darwin")
endif()

# Detect threading flags for Emscripten wasm_threads builds. DuckDB passes
# -pthread in CMAKE_CXX_FLAGS for wasm_threads but not wasm_eh/wasm_mvp.
# ExternalProject builds (htslib, minimap2, WFA2) use explicit CFLAGS and don't
# inherit CMAKE_CXX_FLAGS, so we propagate -pthread to their compile flags here.
# Without this, object files like kthread.o lack atomics/bulk-memory features
# and the linker rejects them when --shared-memory is used. Note:
# CMAKE_CXX_FLAGS may not be visible here (set via -D on parent cmake), so we
# also check DUCKDB_EXPLICIT_PLATFORM which is reliably propagated.
set(EMSCRIPTEN_THREAD_CFLAGS "")
if(EMSCRIPTEN AND (CMAKE_CXX_FLAGS MATCHES "-pthread"
                   OR DUCKDB_EXPLICIT_PLATFORM STREQUAL "wasm_threads"))
  set(EMSCRIPTEN_THREAD_CFLAGS "-pthread")
  message(
    STATUS
      "Emscripten threading detected: propagating -pthread to C library builds")
endif()

# Job count and compiler launcher for the ExternalProject sub-builds in
# cmake/deps/. They run their own make or CMake outside this build's graph, so
# they get neither its parallelism nor the launcher (ccache) DuckDB's CMake
# detects, and without them each compiles one file at a time, uncached.
#
# Jobs: CMAKE_BUILD_PARALLEL_LEVEL if set, else the CPUs this process may run on
# (`nproc` honours taskset and cpuset limits), else CMake's logical core count
# (macOS has no `nproc`). Each sub-build gets the whole count: several start
# together, but only for seconds, and a share of it would leave the longest
# (htslib) short once the others finish. No `make -l`: on a small runner the
# main build's own jobs keep the load above the core count, so it would put
# every sub-build back on one job.
if("$ENV{CMAKE_BUILD_PARALLEL_LEVEL}" MATCHES "^[1-9][0-9]*$")
  set(MIINT_SUBBUILD_JOBS "$ENV{CMAKE_BUILD_PARALLEL_LEVEL}")
else()
  execute_process(
    COMMAND nproc
    OUTPUT_VARIABLE MIINT_SUBBUILD_JOBS
    OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
  if(NOT MIINT_SUBBUILD_JOBS MATCHES "^[1-9][0-9]*$")
    cmake_host_system_information(RESULT MIINT_SUBBUILD_JOBS
                                  QUERY NUMBER_OF_LOGICAL_CORES)
  endif()
endif()

# Launcher: MIINT_SUBBUILD_COMPILER_VARS ("CC=<launcher> <compiler>", "CXX=...")
# for make and autotools sub-builds, MIINT_SUBBUILD_LAUNCHER_ARGS for CMake
# ones, MIINT_SUBBUILD_CC for a sub-build that already names its C compiler.
# Without a launcher the first two are empty and MIINT_SUBBUILD_CC is the bare
# compiler, so those sub-builds keep the compiler they chose before.
set(MIINT_SUBBUILD_CC "${CMAKE_C_COMPILER}")
set(MIINT_SUBBUILD_COMPILER_VARS "")
set(MIINT_SUBBUILD_LAUNCHER_ARGS "")
if(CMAKE_C_COMPILER_LAUNCHER)
  set(MIINT_SUBBUILD_CC "${CMAKE_C_COMPILER_LAUNCHER} ${CMAKE_C_COMPILER}")
  list(APPEND MIINT_SUBBUILD_COMPILER_VARS "CC=${MIINT_SUBBUILD_CC}")
  list(APPEND MIINT_SUBBUILD_LAUNCHER_ARGS
       -DCMAKE_C_COMPILER_LAUNCHER=${CMAKE_C_COMPILER_LAUNCHER})
endif()
if(CMAKE_CXX_COMPILER_LAUNCHER)
  list(APPEND MIINT_SUBBUILD_COMPILER_VARS
       "CXX=${CMAKE_CXX_COMPILER_LAUNCHER} ${CMAKE_CXX_COMPILER}")
  list(APPEND MIINT_SUBBUILD_LAUNCHER_ARGS
       -DCMAKE_CXX_COMPILER_LAUNCHER=${CMAKE_CXX_COMPILER_LAUNCHER})
endif()

# Fail on a source listed twice in <list_name>. The likely cause is a feature's
# source also added to the core list in CMakeLists.txt, which would compile it
# even with that feature OFF.
function(miint_check_unique list_name)
  set(seen "")
  foreach(item ${${list_name}})
    if(item IN_LIST seen)
      message(FATAL_ERROR "${item} is listed twice in ${list_name}")
    endif()
    list(APPEND seen ${item})
  endforeach()
endfunction()
