#include <catch2/catch_test_macros.hpp>

#include <cstdint>
#include <string>
#include <vector>

#include "miint_aligned_vector.hpp"
#include "miint_isa.hpp"

using miint::AlignedVector;
using miint::IsaLevel;
using miint::kSimdAlignment;

namespace {

bool IsAligned(const void *p) {
	return reinterpret_cast<uintptr_t>(p) % kSimdAlignment == 0;
}

} // namespace

// The clamping is the whole safety argument for shipping runtime dispatch: a
// request the CPU or the build cannot satisfy has to come back as the best
// AVAILABLE variant, because honouring it would mean executing an instruction the
// machine does not have. That path cannot be reached through DetectIsa() -- it
// memoizes, so a process tests one combination, and the combination that matters
// is by definition the one this machine cannot run. Hence ResolveIsa is pure.
TEST_CASE("ResolveIsa never returns more than the build and the CPU can do", "[isa]") {
	SECTION("a request beyond the BUILD is clamped, not honoured") {
		REQUIRE(miint::ResolveIsa("avx512", IsaLevel::Baseline, IsaLevel::Avx512) == IsaLevel::Baseline);
		REQUIRE(miint::ResolveIsa("avx512", IsaLevel::Avx2, IsaLevel::Avx512) == IsaLevel::Avx2);
		REQUIRE(miint::ResolveIsa("avx2", IsaLevel::Baseline, IsaLevel::Avx512) == IsaLevel::Baseline);
	}

	SECTION("a request beyond the CPU is clamped, not honoured") {
		REQUIRE(miint::ResolveIsa("avx512", IsaLevel::Avx512, IsaLevel::Baseline) == IsaLevel::Baseline);
		REQUIRE(miint::ResolveIsa("avx512", IsaLevel::Avx512, IsaLevel::Avx2) == IsaLevel::Avx2);
		REQUIRE(miint::ResolveIsa("avx2", IsaLevel::Avx512, IsaLevel::Baseline) == IsaLevel::Baseline);
	}

	SECTION("a request BELOW what is available is honoured -- that is the escape hatch") {
		// The reproducibility override. If this stopped working, someone who needs
		// bit-identical results across machines would silently stop getting them.
		REQUIRE(miint::ResolveIsa("baseline", IsaLevel::Avx512, IsaLevel::Avx512) == IsaLevel::Baseline);
		REQUIRE(miint::ResolveIsa("avx2", IsaLevel::Avx512, IsaLevel::Avx512) == IsaLevel::Avx2);
	}

	SECTION("unset, auto and unrecognised all mean 'the best available'") {
		for (const char *request : {static_cast<const char *>(nullptr), "auto", "", "AVX512", "nonsense"}) {
			REQUIRE(miint::ResolveIsa(request, IsaLevel::Avx2, IsaLevel::Avx512) == IsaLevel::Avx2);
			REQUIRE(miint::ResolveIsa(request, IsaLevel::Avx512, IsaLevel::Baseline) == IsaLevel::Baseline);
		}
		// Note "AVX512" above: the parse is case-SENSITIVE, so a mis-cased value
		// falls through to auto rather than being honoured. Asserted so the
		// behaviour is a decision rather than an accident.
	}
}

TEST_CASE("DetectIsa respects both ceilings on this machine", "[isa]") {
	const IsaLevel got = miint::DetectIsa();
	REQUIRE(static_cast<int>(got) <= static_cast<int>(miint::BuiltIsaCeiling()));
	REQUIRE(static_cast<int>(got) <= static_cast<int>(miint::CpuIsaCeiling()));
}

TEST_CASE("IsaLevelName round-trips every level", "[isa]") {
	REQUIRE(std::string(miint::IsaLevelName(IsaLevel::Baseline)) == "baseline");
	REQUIRE(std::string(miint::IsaLevelName(IsaLevel::Avx2)) == "avx2");
	REQUIRE(std::string(miint::IsaLevelName(IsaLevel::Avx512)) == "avx512");
	// Every name must be one MIINT_SIMD accepts, or the value reported back to a
	// user would not be one they could set.
	for (const IsaLevel level : {IsaLevel::Baseline, IsaLevel::Avx2, IsaLevel::Avx512}) {
		REQUIRE(miint::ResolveIsa(miint::IsaLevelName(level), IsaLevel::Avx512, IsaLevel::Avx512) == level);
	}
}

// Not a tuning assertion. Eigen peels leading scalars off an under-aligned buffer
// onto a different arithmetic path, and the count comes from the runtime address,
// so an unaligned buffer makes a fit's result depend on where the allocator put
// its scratch. See miint_aligned_vector.hpp for the measurement.
TEST_CASE("AlignedVector storage starts on a SIMD boundary", "[isa]") {
	// Sizes chosen to span the allocator's small/large paths, and deliberately not
	// multiples of the alignment.
	for (const size_t n : {1u, 3u, 7u, 100u, 1001u, 100003u}) {
		AlignedVector<double> v(n);
		REQUIRE(IsAligned(v.data()));
		v.resize(n * 2 + 1);
		REQUIRE(IsAligned(v.data()));
	}
}

TEST_CASE("AlignedVector stays aligned as it grows, and keeps its contents", "[isa]") {
	AlignedVector<double> v;
	for (size_t n = 1; n < 5000; n = n * 3 + 1) {
		v.push_back(static_cast<double>(n));
		v.resize(n);
		REQUIRE(IsAligned(v.data()));
		REQUIRE(v.size() == n);
	}
	AlignedVector<double> a(10, 2.5);
	const AlignedVector<double> b = a; // copy must also land aligned
	REQUIRE(IsAligned(b.data()));
	REQUIRE(b[9] == 2.5);
}
