# HMMER 3.4 and Easel 0.49 (the-miint forks; Easel is a submodule inside
# ext/hmmer): profile HMMs. Compiled here from a curated source list, with the
# two config headers generated from the forks' own templates, instead of running
# HMMER's autotools. Off on Windows: never built there, and the generated config
# headers declare a POSIX system (unistd.h, fseeko, times).
#
# The forks' patch series, and why each patch exists, is in
# ext/hmmer/miint/README.md (Easel's in ext/hmmer/easel/miint/README.md).
set(HMMER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ext/hmmer)
set(EASEL_DIR ${HMMER_DIR}/easel)
if(NOT EXISTS ${EASEL_DIR}/easel.h)
  message(
    FATAL_ERROR
      "ext/hmmer/easel is empty: HMMER's Easel is a nested submodule. Run "
      "git submodule update --init --recursive")
endif()

# One SIMD implementation, chosen by target: SSE on x86-64, NEON on arm64, and
# under Emscripten the same SSE code on its SSE emulation (lowered to wasm
# SIMD128), as minimap2 does. HMMER 3.4 has no scalar implementation, and miint
# does not build its VMX (POWER) one.
if(EMSCRIPTEN)
  set(HMMER_SIMD sse)
elseif(APPLE AND OSX_BUILD_ARCH)
  if(OSX_BUILD_ARCH STREQUAL "arm64")
    set(HMMER_SIMD neon)
  else()
    set(HMMER_SIMD sse)
  endif()
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "aarch64|arm64|ARM64")
  set(HMMER_SIMD neon)
elseif(CMAKE_SYSTEM_PROCESSOR MATCHES "x86_64|AMD64|amd64")
  set(HMMER_SIMD sse)
else()
  message(
    FATAL_ERROR
      "miint builds HMMER's SSE or NEON implementation; neither runs on "
      "${CMAKE_SYSTEM_PROCESSOR}. Configure with -DMIINT_ENABLE_HMMER=OFF.")
endif()

# Generate <out> from the autoconf template <in>: each NAME=VALUE in ARGN turns
# the template's "#undef NAME" line into "#define NAME VALUE". Every other line,
# including every other #undef, is kept as the fork has it, so fixed settings
# such as p7_MAXM come from the fork and cannot drift.
function(miint_hmmer_config in out)
  file(READ ${in} text)
  foreach(def ${ARGN})
    string(REGEX MATCH "^[^=]+" name "${def}")
    string(REGEX REPLACE "^[^=]+=" "" value "${def}")
    string(REGEX MATCH "\n#undef ${name}([ \t][^\n]*)?\n" found "${text}")
    if(NOT found)
      message(FATAL_ERROR "${in} has no '#undef ${name}' line")
    endif()
    string(REGEX REPLACE "\n#undef ${name}([ \t][^\n]*)?\n"
                         "\n#define ${name} ${value}\n" text "${text}")
  endforeach()
  # Written only when its text changes (docs/internals/build-system.md, Rules).
  set(old "")
  if(EXISTS ${out})
    file(READ ${out} old)
  endif()
  if(NOT "${old}" STREQUAL "${text}")
    file(WRITE ${out} "${text}")
  endif()
  set_property(
    DIRECTORY
    APPEND
    PROPERTY CMAKE_CONFIGURE_DEPENDS ${in})
endfunction()

# Only macros the sources below read are defined. Which are never defined, and
# why, is in docs/internals/embedded-tools.md (HMMER, generated config headers);
# test/cpp/test_HmmerConfig.cpp fails to compile if one of them is.
if(HMMER_SIMD STREQUAL "neon")
  set(HMMER_SIMD_DEFINES eslENABLE_NEON=1)
  set(EASEL_SIMD_DEFINES eslENABLE_NEON=1 eslHAVE_NEON_AARCH64=1)
else()
  set(HMMER_SIMD_DEFINES eslENABLE_SSE=1)
  set(EASEL_SIMD_DEFINES eslENABLE_SSE=1)
endif()
set(HMMER_CONFIG_DIR ${CMAKE_CURRENT_BINARY_DIR}/hmmer/include)
file(MAKE_DIRECTORY ${HMMER_CONFIG_DIR})
miint_hmmer_config(
  ${HMMER_DIR}/src/p7_config.h.in
  ${HMMER_CONFIG_DIR}/p7_config.h
  "HMMER_VERSION=\"3.4\""
  "HMMER_DATE=\"Aug 2023\""
  "HMMER_COPYRIGHT=\"Copyright (C) 2023 Howard Hughes Medical Institute.\""
  "HMMER_LICENSE=\"Freely distributed under the BSD open source license.\""
  "HMMER_URL=\"http://hmmer.org/\""
  ${HMMER_SIMD_DEFINES})
miint_hmmer_config(
  ${EASEL_DIR}/esl_config.h.in
  ${HMMER_CONFIG_DIR}/esl_config.h
  "EASEL_VERSION=\"0.49\""
  "EASEL_DATE=\"Aug 2023\""
  "EASEL_COPYRIGHT=\"Copyright (C) 2023 Howard Hughes Medical Institute.\""
  "EASEL_LICENSE=\"Freely distributed under the BSD open source license.\""
  eslDEBUGLEVEL=0
  ${EASEL_SIMD_DEFINES}
  HAVE_INTTYPES_H=1
  HAVE_STDINT_H=1
  HAVE_UNISTD_H=1
  HAVE_SYS_TYPES_H=1
  HAVE_STRINGS_H=1
  HAVE_FUNC_ATTRIBUTE_NORETURN=1
  HAVE_FUNC_ATTRIBUTE_FORMAT=1
  HAVE_GETPID=1
  HAVE_STRCASECMP=1
  HAVE_STRSEP=1
  HAVE_TIMES=1
  HAVE_FSEEKO=1)

# The link closure of what miint calls, and nothing else; how it was found is in
# docs/internals/embedded-tools.md (HMMER, Build). A missing file fails the
# native link, and the wasm import check (scripts/build_wasm.sh).
set(HMMER_IMPL_DIR ${HMMER_DIR}/src/impl_${HMMER_SIMD})
add_library(
  hmmer STATIC
  ${HMMER_DIR}/src/emit.c
  ${HMMER_DIR}/src/errors.c
  ${HMMER_DIR}/src/evalues.c
  ${HMMER_DIR}/src/generic_fwdback.c
  ${HMMER_DIR}/src/hmmer.c
  ${HMMER_DIR}/src/logsum.c
  ${HMMER_DIR}/src/modelconfig.c
  ${HMMER_DIR}/src/modelstats.c
  ${HMMER_DIR}/src/p7_bg.c
  ${HMMER_DIR}/src/p7_gmx.c
  ${HMMER_DIR}/src/p7_hmm.c
  ${HMMER_DIR}/src/p7_hmmfile.c
  ${HMMER_DIR}/src/p7_hmmwindow.c
  ${HMMER_DIR}/src/p7_profile.c
  ${HMMER_DIR}/src/p7_trace.c
  ${HMMER_IMPL_DIR}/fwdback.c
  ${HMMER_IMPL_DIR}/msvfilter.c
  ${HMMER_IMPL_DIR}/p7_omx.c
  ${HMMER_IMPL_DIR}/p7_oprofile.c
  ${HMMER_IMPL_DIR}/ssvfilter.c
  ${HMMER_IMPL_DIR}/vitfilter.c
  ${EASEL_DIR}/easel.c
  ${EASEL_DIR}/esl_alphabet.c
  ${EASEL_DIR}/esl_arr2.c
  ${EASEL_DIR}/esl_dirichlet.c
  ${EASEL_DIR}/esl_fileparser.c
  ${EASEL_DIR}/esl_getopts.c
  ${EASEL_DIR}/esl_gumbel.c
  ${EASEL_DIR}/esl_hmm.c
  ${EASEL_DIR}/esl_mem.c
  ${EASEL_DIR}/esl_minimizer.c
  ${EASEL_DIR}/esl_rand64.c
  ${EASEL_DIR}/esl_random.c
  ${EASEL_DIR}/esl_randomseq.c
  ${EASEL_DIR}/esl_sq.c
  ${EASEL_DIR}/esl_${HMMER_SIMD}.c
  ${EASEL_DIR}/esl_ssi.c
  ${EASEL_DIR}/esl_stats.c
  ${EASEL_DIR}/esl_vectorops.c)
set_target_properties(hmmer PROPERTIES POSITION_INDEPENDENT_CODE ON)
# The generated headers first. HMMER's headers include <p7_config.h> with angle
# brackets, so the include path decides which copy they get.
set(HMMER_INCLUDE_DIRS ${HMMER_CONFIG_DIR} ${HMMER_DIR}/src ${EASEL_DIR})
target_include_directories(hmmer PRIVATE ${HMMER_INCLUDE_DIRS})
# Vendored C, compiled as upstream ships it: its warnings are upstream's.
target_compile_options(hmmer PRIVATE -w)
if(EMSCRIPTEN)
  set(HMMER_WASM_SIMD_FLAGS -msimd128 -msse -msse2)
  target_compile_options(hmmer PRIVATE ${HMMER_WASM_SIMD_FLAGS}
                                       ${EMSCRIPTEN_THREAD_CFLAGS})
  message(STATUS "HMMER building for Emscripten/WASM with SSE emulation")
endif()
message(STATUS "HMMER SIMD implementation: ${HMMER_SIMD}")

# miint's translation units that include HMMER's headers. hmmer.h includes the
# SIMD implementation's header, so under Emscripten these need the library's
# SIMD flags too. The include path is scoped to these files, not added for every
# target: HMMER's src/ directory holds generically named headers.
set(HMMER_MIINT_SOURCES src/HmmerProfile.cpp test/cpp/test_HmmerConfig.cpp
                        test/cpp/test_HmmerProfile.cpp)
set_source_files_properties(
  ${HMMER_MIINT_SOURCES} PROPERTIES INCLUDE_DIRECTORIES "${HMMER_INCLUDE_DIRS}")
if(EMSCRIPTEN)
  set_source_files_properties(
    ${HMMER_MIINT_SOURCES} PROPERTIES COMPILE_OPTIONS
                                      "${HMMER_WASM_SIMD_FLAGS}")
endif()

miint_git_version(HMMER ext/hmmer)
miint_git_version(EASEL ext/hmmer/easel)
list(APPEND MIINT_DEFINITIONS MIINT_HAS_HMMER)
list(APPEND MIINT_EXTENSION_SOURCES src/HmmerProfile.cpp src/read_hmm.cpp)
list(APPEND MIINT_TEST_SOURCES src/HmmerProfile.cpp
     test/cpp/test_HmmerConfig.cpp test/cpp/test_HmmerProfile.cpp)
# Linked by archive path, not target name: miint_extension goes into DuckDB's
# install(EXPORT) set, which rejects link dependencies on targets outside it.
list(APPEND MIINT_LINK_LIBRARIES $<TARGET_FILE:hmmer>)
list(APPEND MIINT_TARGET_DEPENDENCIES hmmer)
list(APPEND MIINT_WASM_LIBS $<TARGET_FILE:hmmer>)
