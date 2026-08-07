#include "miint_isa.hpp"

#include <cstdlib>
#include <cstring>

namespace miint {

namespace {

//! Does the CPU support everything the AVX2 variant is compiled with?
//!
//! `__builtin_cpu_supports` rather than raw `__cpuid_count` because it also
//! consults XCR0 -- the OS has to have enabled the wider register state, and a
//! CPUID feature bit alone does not tell you that. Getting that wrong is exactly
//! how CPU dispatch turns into SIGILL on an old kernel or inside a hypervisor
//! that masks AVX state. GCC and Clang both provide it on x86.
bool CpuHasAvx2() {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
	return __builtin_cpu_supports("avx2") && __builtin_cpu_supports("fma");
#else
	return false;
#endif
}

//! The AVX-512 variant is compiled with F/DQ/VL/BW, so every one of them must be
//! present -- not just the `avx512f` that a coarser check would settle for. Skylake-X
//! era parts have all four; Knights Landing has F but neither DQ nor VL, and would
//! fault on the first instruction from either.
bool CpuHasAvx512() {
#if defined(__x86_64__) || defined(_M_X64) || defined(__i386__)
	return __builtin_cpu_supports("avx512f") && __builtin_cpu_supports("avx512dq") &&
	       __builtin_cpu_supports("avx512vl") && __builtin_cpu_supports("avx512bw") && __builtin_cpu_supports("fma");
#else
	return false;
#endif
}

//! Parse MIINT_SIMD. Unset or unrecognised means auto, which is represented by
//! the build ceiling -- i.e. "no cap beyond what exists".
IsaLevel RequestedCeiling() {
	const char *env = std::getenv("MIINT_SIMD");
	if (env == nullptr) {
		return IsaLevel::Avx512;
	}
	if (std::strcmp(env, "baseline") == 0) {
		return IsaLevel::Baseline;
	}
	if (std::strcmp(env, "avx2") == 0) {
		return IsaLevel::Avx2;
	}
	if (std::strcmp(env, "avx512") == 0) {
		return IsaLevel::Avx512;
	}
	// "auto" and anything unrecognised. Deliberately not an error: this is read
	// on a hot-ish path with no way to report, and a typo silently disabling
	// dispatch would be a worse failure than a typo being ignored.
	return IsaLevel::Avx512;
}

IsaLevel Min(IsaLevel a, IsaLevel b) {
	return static_cast<int>(a) < static_cast<int>(b) ? a : b;
}

} // namespace

IsaLevel BuiltIsaCeiling() {
#if defined(MIINT_HAS_AVX512)
	return IsaLevel::Avx512;
#elif defined(MIINT_HAS_AVX2)
	return IsaLevel::Avx2;
#else
	return IsaLevel::Baseline;
#endif
}

IsaLevel CpuIsaCeiling() {
	if (CpuHasAvx512()) {
		return IsaLevel::Avx512;
	}
	if (CpuHasAvx2()) {
		return IsaLevel::Avx2;
	}
	return IsaLevel::Baseline;
}

IsaLevel DetectIsa() {
	// Resolved once per process. CPUID is not expensive, but the objective calls
	// through this selection once per evaluation -- 1072 times in a cystic-fibrosis
	// L-BFGS fit, 196000 in an Adam one -- and a static keeps that a load.
	static const IsaLevel level = Min(RequestedCeiling(), Min(BuiltIsaCeiling(), CpuIsaCeiling()));
	return level;
}

const char *IsaLevelName(IsaLevel level) {
	switch (level) {
	case IsaLevel::Avx512:
		return "avx512";
	case IsaLevel::Avx2:
		return "avx2";
	case IsaLevel::Baseline:
		break;
	}
	return "baseline";
}

} // namespace miint
