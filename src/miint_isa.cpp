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

IsaLevel Min(IsaLevel a, IsaLevel b) {
	return static_cast<int>(a) < static_cast<int>(b) ? a : b;
}

} // namespace

IsaLevel ResolveIsa(const char *request, IsaLevel built, IsaLevel cpu) {
	// Unset, "auto", and anything unrecognised all mean "no cap beyond what
	// exists". Unrecognised is deliberately not an error: there is nowhere to
	// report it from, and a typo silently disabling dispatch would be a worse
	// failure than a typo being ignored.
	IsaLevel requested = IsaLevel::Avx512;
	if (request != nullptr) {
		if (std::strcmp(request, "baseline") == 0) {
			requested = IsaLevel::Baseline;
		} else if (std::strcmp(request, "avx2") == 0) {
			requested = IsaLevel::Avx2;
		} else if (std::strcmp(request, "avx512") == 0) {
			requested = IsaLevel::Avx512;
		}
	}
	// CLAMPING, not honouring, is what keeps a forced request from becoming a
	// SIGILL: asking for avx512 on a machine or a build without it yields the
	// best available variant instead of an illegal instruction.
	return Min(requested, Min(built, cpu));
}

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
	static const IsaLevel level = ResolveIsa(std::getenv("MIINT_SIMD"), BuiltIsaCeiling(), CpuIsaCeiling());
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
