#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <stdexcept>
#include <vector>

#include "simulate_resemblance.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;

// These pin the DETERMINISTIC parts of the gradient model with exact numeric
// expectations, so a broken formula (exponent sign, normalization constant,
// linspace spacing) fails here regardless of the (non-portable) RNG stream.
// The random/statistical fidelity is checked separately by the SQL parity test.

TEST_CASE("Linspace matches numpy semantics", "[simulate]") {
	// n == 1 collapses to {lo}.
	auto one = miint::Linspace(0.3, 0.9, 1);
	REQUIRE(one.size() == 1);
	CHECK(one[0] == Approx(0.3));

	// n >= 2 includes both endpoints, evenly spaced.
	auto two = miint::Linspace(0.1, 0.9, 2);
	REQUIRE(two.size() == 2);
	CHECK(two[0] == Approx(0.1));
	CHECK(two[1] == Approx(0.9));

	auto three = miint::Linspace(0.1, 0.9, 3);
	REQUIRE(three.size() == 3);
	CHECK(three[0] == Approx(0.1));
	CHECK(three[1] == Approx(0.5));
	CHECK(three[2] == Approx(0.9));

	auto five = miint::Linspace(0.1, 0.9, 5);
	REQUIRE(five.size() == 5);
	CHECK(five[0] == Approx(0.1));
	CHECK(five[1] == Approx(0.3));
	CHECK(five[2] == Approx(0.5));
	CHECK(five[3] == Approx(0.7));
	CHECK(five[4] == Approx(0.9));
}

TEST_CASE("BuildTrueDistribution pins the Gaussian response and normalization", "[simulate]") {
	// One sample at position 0.5; two species with optima 0.5 (on-peak) and 0.3
	// (|delta|=0.2 => exponent -2 at sp_width 0.1). Expected values computed
	// independently (norm = 1/(0.1*sqrt(2*pi)) = 3.989422804014327).
	const std::vector<double> abundances = {1.0, 2.0};
	const std::vector<double> positions = {0.5};
	const std::vector<double> optima = {0.5, 0.3};
	auto tru = miint::BuildTrueDistribution(abundances, positions, optima, 0.1);

	REQUIRE(tru.size() == 1);
	REQUIRE(tru[0].size() == 2);
	// On-peak: abund * 1 * norm.
	CHECK(tru[0][0] == Approx(3.989422804014327).epsilon(1e-12));
	// Off-peak: 2 * exp(-2) * norm. A flipped exponent sign would give 2*exp(+2)*norm ~ 58.9.
	CHECK(tru[0][1] == Approx(1.079819330263761).epsilon(1e-12));
	// Off-peak strictly below on-peak (basic monotonicity of the response curve).
	CHECK(tru[0][1] < tru[0][0]);
}

TEST_CASE("SimulateGradient is deterministic for a fixed seed", "[simulate]") {
	const std::vector<double> abund = {0.4, 0.25, 0.2, 0.1, 0.05};
	auto a = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "*sample", 0.1, 0.9, 123);
	auto b = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "*sample", 0.1, 0.9, 123);
	REQUIRE(a.size() == b.size());
	CHECK(a.sample_id == b.sample_id);
	CHECK(a.otu_id == b.otu_id);
	CHECK(a.count == b.count);
	CHECK(a.ground_truth == b.ground_truth);
}

TEST_CASE("SimulateGradient enforces per-sample depth and nonzero-only COO", "[simulate]") {
	const std::vector<double> abund = {0.4, 0.25, 0.2, 0.1, 0.05};
	auto coo = miint::SimulateGradient(abund, 6, 500, 0.1, 0.0, "*sample", 0.1, 0.9, 7);

	// All emitted counts are >= 1 (COO stores no zeros).
	for (auto c : coo.count) {
		CHECK(c >= 1);
	}
	// Each of the 6 samples sums to exactly 500 reads (multinomial invariant).
	std::vector<int64_t> per_sample(6, 0);
	for (size_t i = 0; i < coo.size(); i++) {
		REQUIRE(coo.sample_id[i] >= 0);
		REQUIRE(coo.sample_id[i] < 6);
		per_sample[static_cast<size_t>(coo.sample_id[i])] += coo.count[i];
	}
	for (int s = 0; s < 6; s++) {
		CHECK(per_sample[static_cast<size_t>(s)] == 500);
	}
}

TEST_CASE("SimulateGradient '+species' noise scales the width per species", "[simulate]") {
	const std::vector<double> abund = {0.4, 0.25, 0.2, 0.1, 0.05};

	// '+species' draws each species' noise width from that species' own abundance
	// (x += N(0, noise*x) == x*(1 + noise*Z)), so the perturbation is a per-element
	// RELATIVE jitter. That makes the whole pipeline invariant to a uniform rescale
	// of the abundance vector: every step downstream (floor-shift, mean-rescale,
	// multinomial normalization) is homogeneous in the row, so scaling the input by
	// k must reproduce the identical count matrix at a fixed seed.
	std::vector<double> scaled(abund.size());
	for (size_t i = 0; i < abund.size(); i++) {
		scaled[i] = abund[i] * 1000.0;
	}

	auto base = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "+species", 0.1, 0.9, 123);
	auto big = miint::SimulateGradient(scaled, 6, 500, 0.1, 0.5, "+species", 0.1, 0.9, 123);
	REQUIRE(base.size() == big.size());
	CHECK(base.sample_id == big.sample_id);
	CHECK(base.otu_id == big.otu_id);
	CHECK(base.count == big.count);

	// '*sample' takes its width from the SAMPLE total instead, so the same rescale
	// changes the noise magnitude and the realization must differ. This is the
	// behavioural difference between the two widths, not just a different RNG path.
	auto s_base = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "*sample", 0.1, 0.9, 123);
	auto s_big = miint::SimulateGradient(scaled, 6, 500, 0.1, 0.5, "*sample", 0.1, 0.9, 123);
	CHECK(s_base.count != s_big.count);

	// '+species' is a genuinely different model from BOTH sample-total widths at the
	// same seed. The '+sample' comparison matters most: it is additive too and is
	// likewise scale-equivariant, so the rescale check above cannot tell them apart
	// on its own.
	CHECK(base.count != s_base.count);
	auto add_sample = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "+sample", 0.1, 0.9, 123);
	CHECK(base.count != add_sample.count);

	// ... and deterministic in its own right.
	auto again = miint::SimulateGradient(abund, 6, 500, 0.1, 0.5, "+species", 0.1, 0.9, 123);
	CHECK(base.count == again.count);
}

TEST_CASE("SimulateGradient '+species' gives a zero-abundance species zero width", "[simulate]") {
	// A species with zero abundance has noise width noise*0 == 0 under '+species',
	// so it is never perturbed and stays absent. Under '+sample' the width comes
	// from the sample total, so the same species IS jittered off zero and (after the
	// floor-shift lifts the row) becomes observable. This pins the per-species
	// width semantics on behaviour a caller can see, not on an RNG stream.
	const std::vector<double> abund = {0.5, 0.5, 0.0};
	const int32_t nsamp = 4;
	const std::vector<int64_t> seeds = {11, 42, 7, 123, 2024};

	// '+species' side is an analytic guarantee, not a seed accident: species 2 takes
	// width noise*0 == 0, and at noise 0.05 a positive species only goes negative on
	// a |z| > 20 draw, so the row minimum IS species 2 and the floor-shift is a
	// no-op. It must therefore be absent for EVERY seed.
	for (int64_t seed : seeds) {
		auto sp = miint::SimulateGradient(abund, nsamp, 400, 0.1, 0.05, "+species", 0.1, 0.9, seed);
		for (size_t i = 0; i < sp.size(); i++) {
			CHECK(sp.otu_id[i] != 2);
		}
	}

	// '+sample' side is probabilistic per seed (species 2 is jittered off zero, then
	// has to win at least one of 4 x 400 multinomial draws), so requiring it at one
	// fixed seed would be brittle against any change in how many draws are consumed
	// upstream. Require it on a majority of seeds instead: that is what distinguishes
	// the sample-total width from the per-species width, without pinning an RNG path.
	int emitted_in = 0;
	for (int64_t seed : seeds) {
		auto sa = miint::SimulateGradient(abund, nsamp, 400, 0.1, 0.05, "+sample", 0.1, 0.9, seed);
		for (size_t i = 0; i < sa.size(); i++) {
			if (sa.otu_id[i] == 2) {
				emitted_in++;
				break;
			}
		}
	}
	CHECK(emitted_in * 2 > static_cast<int>(seeds.size()));
}

TEST_CASE("SimulateGradient throws on a degenerate empty sample", "[simulate]") {
	// A single species under multiplicative noise: the lone element is its own
	// minimum, so the floor-subtract zeroes the whole row -> empty sample.
	const std::vector<double> one_sp = {1.0};
	CHECK_THROWS_WITH(miint::SimulateGradient(one_sp, 3, 100, 0.1, 0.5, "*sample", 0.1, 0.9, 7),
	                  ContainsSubstring("empty sample"));
}
