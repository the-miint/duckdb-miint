set(VCPKG_ENV_PASSTHROUGH_UNTRACKED EMSCRIPTEN_ROOT EMSDK PATH)

if(NOT DEFINED ENV{EMSCRIPTEN_ROOT})
   find_path(EMSCRIPTEN_ROOT "emcc")
else()
   set(EMSCRIPTEN_ROOT "$ENV{EMSCRIPTEN_ROOT}")
endif()

if(NOT EMSCRIPTEN_ROOT)
   if(NOT DEFINED ENV{EMSDK})
      message(FATAL_ERROR "The emcc compiler not found in PATH")
   endif()
   set(EMSCRIPTEN_ROOT "$ENV{EMSDK}/upstream/emscripten")
endif()

if(NOT EXISTS "${EMSCRIPTEN_ROOT}/cmake/Modules/Platform/Emscripten.cmake")
   message(FATAL_ERROR "Emscripten.cmake toolchain file not found")
endif()

set(VCPKG_TARGET_ARCHITECTURE wasm32)
set(VCPKG_CRT_LINKAGE dynamic)
set(VCPKG_LIBRARY_LINKAGE static)
set(VCPKG_CMAKE_SYSTEM_NAME Emscripten)
set(VCPKG_CHAINLOAD_TOOLCHAIN_FILE "${EMSCRIPTEN_ROOT}/cmake/Modules/Platform/Emscripten.cmake")

# Emscripten's fenv.h stubs out floating-point exception support (FE_ALL_EXCEPT=0)
# but omits individual exception flag macros like FE_INVALID. HDF5 1.14.x uses
# FE_INVALID unconditionally (fixed upstream in HDF5 2.0.0 via HDFGroup/hdf5#4952).
# Define it as 0 (no-op) to match Emscripten's no-FP-exception semantics.
# Injected via VCPKG_CMAKE_CONFIGURE_OPTIONS because VCPKG_C_FLAGS doesn't propagate
# for Emscripten (vcpkg has no emscripten toolchain in scripts/toolchains/).
set(VCPKG_CMAKE_CONFIGURE_OPTIONS "-DCMAKE_C_FLAGS=-DFE_INVALID=0" "-DCMAKE_CXX_FLAGS=-DFE_INVALID=0")
