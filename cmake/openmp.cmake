# The miint_openmp INTERFACE target: the OpenMP runtime for UniFrac's archives
# (required) and krepp's index builds (optional). Sets MIINT_OPENMP_AVAILABLE.

# libssu_inmem.a / libskbb_inmem.a are compiled with -fopenmp; the final
# extension must resolve their OpenMP runtime symbols (__kmpc_*, omp_*).
#
# Linux: find_package(OpenMP) discovers libgomp via system paths. macOS: Apple
# Clang doesn't ship libomp. Probe Homebrew's libomp (set up via the `omp`
# extra_toolchain hook in our workflow) and link libomp.a *statically*
# (positional arg, not -lomp) so the released .duckdb_extension doesn't carry an
# RPATH to /opt/homebrew/opt/libomp. Mirrors GPL-boundary's build.rs
# locate_macos_libomp().
#
# Exposes a single INTERFACE target `miint_openmp` so link sites don't need to
# know which platform-specific dance produced it.
#
# MIINT_OPENMP_STATIC_LIBOMP is set below whenever libomp.a is the runtime
# linked, and passed to the sources as a definition: a static libomp has to be
# kept from shutting down (docs/internals/embedded-tools.md, "The OpenMP
# runtime").
#
# krepp wants it too, to build indexes on more than one core. Unlike UniFrac it
# can do without: compiled with _WOPENMP=0 its pragmas are ignored and the build
# stays serial, which is what every release before this did. So a missing
# runtime is fatal only when UniFrac is the one asking, and
# MIINT_OPENMP_AVAILABLE is what cmake/deps/krepp.cmake consults.
set(MIINT_OPENMP_STATIC_LIBOMP OFF)
if((MIINT_ENABLE_UNIFRAC OR MIINT_ENABLE_KREPP) AND NOT EMSCRIPTEN)
  set(MIINT_OPENMP_AVAILABLE ON)
  if(APPLE)
    # Opt-in override for hosts without Homebrew: point at a prefix holding
    # include/omp.h + lib/libomp.{a,dylib} (e.g. a conda env's $CONDA_PREFIX).
    # Set once via -DMIINT_LIBOMP_PREFIX=... and it persists in the cache. Empty
    # (default) keeps the Homebrew auto-probe below untouched, so CI and
    # distribution builds are entirely unaffected.
    set(MIINT_LIBOMP_PREFIX
        ""
        CACHE
          PATH
          "macOS: explicit libomp prefix (include/omp.h + lib/libomp.{a,dylib}); empty = auto-probe Homebrew"
    )
    # IMPORTED INTERFACE (not a regular INTERFACE library): exempts it from
    # install(EXPORT)'s "all linked targets must be in the export set" check,
    # mirroring how OpenMP::OpenMP_CXX itself behaves.
    add_library(miint_openmp INTERFACE IMPORTED GLOBAL)
    if(MIINT_LIBOMP_PREFIX)
      # MIINT_LIBOMP_PREFIX is a cache entry, so it outlives the configure that
      # set it. Now that krepp enters this block too, a prefix that has since
      # gone stale must not abort a build that never asked for one - only
      # UniFrac, which cannot proceed without it, still fails.
      if(NOT EXISTS "${MIINT_LIBOMP_PREFIX}/include/omp.h")
        if(MIINT_ENABLE_UNIFRAC)
          message(
            FATAL_ERROR
              "MIINT_LIBOMP_PREFIX='${MIINT_LIBOMP_PREFIX}' has no include/omp.h"
          )
        endif()
        message(
          STATUS
            "MIINT_LIBOMP_PREFIX='${MIINT_LIBOMP_PREFIX}' has no include/omp.h - "
            "krepp index builds stay single-threaded.")
        set(MIINT_OPENMP_AVAILABLE OFF)
      endif()
      # Also consumed by the skbb/unifrac inmem sub-makes' INMEM_MPFLAG
      # (cmake/deps/unifrac.cmake).
      set(MIINT_LIBOMP_INCLUDE_DIR "${MIINT_LIBOMP_PREFIX}/include")
      set_property(TARGET miint_openmp PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                "${MIINT_LIBOMP_INCLUDE_DIR}")
      if(EXISTS "${MIINT_LIBOMP_PREFIX}/lib/libomp.a")
        # Static: link positionally, no RPATH (distribution-shaped).
        message(
          STATUS
            "Using static libomp from MIINT_LIBOMP_PREFIX: ${MIINT_LIBOMP_PREFIX}"
        )
        set_property(
          TARGET miint_openmp PROPERTY INTERFACE_LINK_LIBRARIES
                                       "${MIINT_LIBOMP_PREFIX}/lib/libomp.a")
        set(MIINT_OPENMP_STATIC_LIBOMP ON)
      elseif(EXISTS "${MIINT_LIBOMP_PREFIX}/lib/libomp.dylib")
        # Dynamic (conda-forge llvm-openmp ships no .a): link the dylib and
        # record an RPATH to it. Fine for a local build; not the shape we
        # distribute — hence this path is only ever reached by opt-in.
        message(
          STATUS
            "Using dynamic libomp from MIINT_LIBOMP_PREFIX: ${MIINT_LIBOMP_PREFIX} (adds an RPATH; local builds only)"
        )
        set_property(
          TARGET miint_openmp
          PROPERTY INTERFACE_LINK_LIBRARIES
                   "${MIINT_LIBOMP_PREFIX}/lib/libomp.dylib")
        set_property(
          TARGET miint_openmp PROPERTY INTERFACE_LINK_OPTIONS
                                       "-Wl,-rpath,${MIINT_LIBOMP_PREFIX}/lib")
      else()
        if(MIINT_ENABLE_UNIFRAC)
          message(
            FATAL_ERROR
              "MIINT_LIBOMP_PREFIX='${MIINT_LIBOMP_PREFIX}' has neither lib/libomp.a nor lib/libomp.dylib"
          )
        endif()
        message(
          STATUS
            "MIINT_LIBOMP_PREFIX='${MIINT_LIBOMP_PREFIX}' has no libomp library - "
            "krepp index builds stay single-threaded.")
        set(MIINT_OPENMP_AVAILABLE OFF)
      endif()
    else()
      # Default: probe Homebrew for a STATIC libomp.a (the CI path, unchanged).
      execute_process(
        COMMAND brew --prefix libomp
        OUTPUT_VARIABLE BREW_LIBOMP_PREFIX
        OUTPUT_STRIP_TRAILING_WHITESPACE ERROR_QUIET)
      set(MIINT_LIBOMP_FOUND_PREFIX "")
      foreach(candidate "${BREW_LIBOMP_PREFIX}" "/opt/homebrew/opt/libomp"
                        "/usr/local/opt/libomp")
        if(candidate
           AND EXISTS "${candidate}/include/omp.h"
           AND EXISTS "${candidate}/lib/libomp.a")
          set(MIINT_LIBOMP_FOUND_PREFIX "${candidate}")
          break()
        endif()
      endforeach()
      if(NOT MIINT_LIBOMP_FOUND_PREFIX AND MIINT_ENABLE_UNIFRAC)
        message(
          FATAL_ERROR
            "OpenMP static archive (libomp.a) not found on macOS.\n"
            "  Install with: brew install libomp\n"
            "  Or set -DMIINT_LIBOMP_PREFIX=/prefix (include/omp.h + lib/libomp.{a,dylib}).\n"
            "  Required when MIINT_ENABLE_UNIFRAC=ON.")
      elseif(NOT MIINT_LIBOMP_FOUND_PREFIX)
        message(
          STATUS
            "No libomp on macOS - krepp index builds stay single-threaded. "
            "brew install libomp, or -DMIINT_LIBOMP_PREFIX=/prefix, to build them on more than one core."
        )
        set(MIINT_OPENMP_AVAILABLE OFF)
      endif()
      if(MIINT_OPENMP_AVAILABLE)
        message(STATUS "Found Homebrew libomp at: ${MIINT_LIBOMP_FOUND_PREFIX}")
        # Also consumed by the skbb/unifrac inmem sub-makes' INMEM_MPFLAG
        # (cmake/deps/unifrac.cmake).
        set(MIINT_LIBOMP_INCLUDE_DIR "${MIINT_LIBOMP_FOUND_PREFIX}/include")
        set_property(TARGET miint_openmp PROPERTY INTERFACE_INCLUDE_DIRECTORIES
                                                  "${MIINT_LIBOMP_INCLUDE_DIR}")
        # Positional .a (not -lomp): Homebrew ships both libomp.a and
        # libomp.dylib; -lomp prefers the dylib, baking a Homebrew path into the
        # extension's RPATH.
        set_property(
          TARGET miint_openmp
          PROPERTY INTERFACE_LINK_LIBRARIES
                   "${MIINT_LIBOMP_FOUND_PREFIX}/lib/libomp.a")
        set(MIINT_OPENMP_STATIC_LIBOMP ON)
      endif()
    endif()
  else()
    if(MIINT_ENABLE_UNIFRAC)
      find_package(OpenMP REQUIRED COMPONENTS CXX)
    else()
      find_package(OpenMP QUIET COMPONENTS CXX)
    endif()
    if(OpenMP_CXX_FOUND)
      add_library(miint_openmp INTERFACE IMPORTED GLOBAL)
      set_property(TARGET miint_openmp PROPERTY INTERFACE_LINK_LIBRARIES
                                                OpenMP::OpenMP_CXX)
    else()
      message(
        STATUS "OpenMP not found - krepp index builds stay single-threaded.")
      set(MIINT_OPENMP_AVAILABLE OFF)
    endif()
  endif()
endif()
# A prefix can hold libomp.a without omp.h; then nothing is linked after all.
if(NOT MIINT_OPENMP_AVAILABLE)
  set(MIINT_OPENMP_STATIC_LIBOMP OFF)
endif()
if(MIINT_OPENMP_STATIC_LIBOMP)
  list(APPEND MIINT_DEFINITIONS MIINT_OPENMP_STATIC_LIBOMP)
endif()
