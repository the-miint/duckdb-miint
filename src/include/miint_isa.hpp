#pragma once
//
// Which x86 instruction set may miint's dispatched kernels use on THIS machine?
//
// Deliberately subsystem-agnostic. mmvec is the first caller, not the intended
// only one: Eigen fixes its packet width at compile time, so any future dense
// kernel that wants a wider one needs the same build-time-variants plus
// runtime-selection arrangement, and it should not have to re-derive it. The
// companion CMake helper is `miint_add_isa_variants()`.
//
// The pattern -- build every variant, choose at load time from CPUID -- is the
// one htslib's htscodecs already uses, and the shipped CI artifact proves it
// survives the distribution pipeline: it carries AVX2 and AVX-512 rANS kernels,
// abPOA's AVX-512 aligner, and OpenSSL's AVX-512 AES. DuckDB's platform string
// (`os_arch`, plus at most a `_musl`/`_android`/`_mingw` postfix) has no
// CPU-feature component, so a per-ISA download is not expressible; runtime
// selection inside one binary is the only mechanism that reaches users.
//
// Non-x86 targets (arm64, wasm) always report Baseline: NEON and wasm SIMD are
// unconditional there, so Eigen already vectorizes and there is nothing to pick.

namespace miint {

//! Instruction-set tiers miint compiles kernels for. Ordered, so `<=` means
//! "no wider than", which is what the clamping in DetectIsa() relies on.
enum class IsaLevel : int {
	Baseline = 0, //!< whatever the target's default is (SSE2 on x86-64)
	Avx2 = 1,     //!< AVX2 + FMA
	Avx512 = 2,   //!< AVX-512 F/DQ/VL/BW + FMA
};

//! The level this process should use, resolved ONCE into a function-local static.
//!
//! Honours the `MIINT_SIMD` environment variable -- `baseline`, `avx2`, `avx512`
//! or `auto` (the default). A request the CPU or the build cannot satisfy is
//! CLAMPED down to the best available rather than honoured, so forcing `avx512`
//! on a machine without it yields a slower answer, never SIGILL. An unrecognised
//! value is treated as `auto`.
//!
//! Deliberately an environment variable rather than a DuckDB setting: it has to
//! be readable from the C++ unit tests, which never construct a database.
IsaLevel DetectIsa();

//! The widest level this BUILD contains variants for, ignoring the CPU and the
//! override. Baseline unless the toolchain accepted the wider flags at configure
//! time. Exposed for tests and for reporting.
IsaLevel BuiltIsaCeiling();

//! The widest level THIS CPU supports, ignoring the build and the override.
IsaLevel CpuIsaCeiling();

//! Lower-case name (`baseline` / `avx2` / `avx512`), suitable for the `message`
//! column of a fit and for test assertions. Never null.
const char *IsaLevelName(IsaLevel level);

} // namespace miint
