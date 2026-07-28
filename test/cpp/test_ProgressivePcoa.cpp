#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include <cmath>
#include <cstdint>
#include <random>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <vector>

#include "ordination.h" // skbb_pcoa_fsvd_fp32

#include "procrustes_core.hpp" // FullDisparity — frame-invariant M^2 for the oracle comparison
#include "progressive_pcoa_core.hpp"

using Catch::Approx;
using miint::procrustes::FullDisparity;
using miint::progressive::DistanceBlock;
using miint::progressive::ProgressivePcoaResult;
using miint::progressive::RunProgressivePcoa;

namespace {

// A synthetic Euclidean oracle: known points in R^d_true, from which an exact
// Euclidean distance matrix is built. Classical PCoA of a Euclidean distance
// matrix recovers the original configuration (up to a similarity transform), so
// the full PCoA is a ground truth the progressive result must match to within a
// procrustes disparity M^2 — the same correctness criterion as the origin PoC
// (2025.07.30-progressive-pcoa/poc.py), which scores progressive vs full via
// q2 procrustes_analysis.
struct Oracle {
	uint32_t n;
	uint32_t d_true;
	std::vector<std::string> ids; // "s0000".. — lexicographically stable
	std::vector<double> points;   // n * d_true row-major
	std::vector<float> dm;        // n * n row-major fp32 Euclidean distances
	std::unordered_map<std::string, uint32_t> index;
};

Oracle MakeEuclideanOracle(uint32_t n, uint32_t d_true, uint64_t seed) {
	Oracle o;
	o.n = n;
	o.d_true = d_true;
	std::mt19937_64 rng(seed);
	std::uniform_real_distribution<double> u(-5.0, 5.0);
	o.points.resize(static_cast<size_t>(n) * d_true);
	o.ids.reserve(n);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t k = 0; k < d_true; ++k) {
			o.points[i * d_true + k] = u(rng);
		}
		char buf[16];
		std::snprintf(buf, sizeof(buf), "s%04u", i);
		o.ids.emplace_back(buf);
		o.index.emplace(o.ids.back(), i);
	}
	o.dm.resize(static_cast<size_t>(n) * n, 0.0f);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j) {
			double acc = 0.0;
			for (uint32_t k = 0; k < d_true; ++k) {
				const double diff = o.points[i * d_true + k] - o.points[j * d_true + k];
				acc += diff * diff;
			}
			const float dist = static_cast<float>(std::sqrt(acc));
			o.dm[i * n + j] = dist;
			o.dm[j * n + i] = dist;
		}
	}
	return o;
}

// Full-matrix PCoA (the ground truth), returned as an n × n_dims row-major double
// matrix in the oracle's id order. Mirrors RunPcoaOnMatrix's skbb call/layout.
std::vector<double> FullPcoa(const Oracle &o, uint32_t n_dims, int seed) {
	std::vector<float> eig(n_dims), samples(static_cast<size_t>(o.n) * n_dims), prop(n_dims);
	skbb_pcoa_fsvd_fp32(o.n, o.dm.data(), n_dims, seed, eig.data(), samples.data(), prop.data());
	std::vector<double> out(static_cast<size_t>(o.n) * n_dims);
	for (size_t i = 0; i < out.size(); ++i) {
		out[i] = static_cast<double>(samples[i]);
	}
	return out;
}

// Slice the oracle's precomputed full DM into a dense block over exactly the
// requested ids (in the requested order). Stands in for the real block sources:
// a COO slice (progressive_pcoa_from_distances) or an on-the-fly UniFrac compute
// (progressive_pcoa_from_unifrac).
DistanceBlock SliceBlock(const Oracle &o, const std::vector<std::string> &requested) {
	DistanceBlock block;
	block.ids = requested;
	const uint32_t m = static_cast<uint32_t>(requested.size());
	block.matrix.resize(static_cast<size_t>(m) * m, 0.0f);
	for (uint32_t r = 0; r < m; ++r) {
		const uint32_t ri = o.index.at(requested[r]);
		for (uint32_t c = 0; c < m; ++c) {
			const uint32_t ci = o.index.at(requested[c]);
			block.matrix[static_cast<size_t>(r) * m + c] = o.dm[static_cast<size_t>(ri) * o.n + ci];
		}
	}
	return block;
}

// Reassemble a progressive result into an n × d row-major matrix in the oracle's
// id order (so it lines up row-for-row with the full-PCoA ground truth), while
// asserting every (sample, axis) appears exactly once.
std::vector<double> AssembleInOracleOrder(const Oracle &o, const ProgressivePcoaResult &result, uint32_t d) {
	std::vector<double> out(static_cast<size_t>(o.n) * d);
	std::vector<int> filled(static_cast<size_t>(o.n) * d, 0);
	for (const auto &c : result.coords) {
		auto it = o.index.find(c.sample_id);
		REQUIRE(it != o.index.end());
		REQUIRE(c.axis >= 0);
		REQUIRE(static_cast<uint32_t>(c.axis) < d);
		const size_t pos = static_cast<size_t>(it->second) * d + static_cast<uint32_t>(c.axis);
		REQUIRE(filled[pos] == 0); // no duplicate (sample, axis)
		filled[pos] = 1;
		out[pos] = c.coordinate;
	}
	for (int f : filled) {
		REQUIRE(f == 1);
	}
	return out;
}

} // namespace

TEST_CASE("progressive PCoA matches full PCoA on a Euclidean oracle", "[progressive]") {
	const uint32_t n = 60, d_true = 3, n_dims = 3;
	const int seed = 42;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/12345);

	// Ground truth: PCoA over the whole distance matrix.
	const std::vector<double> full = FullPcoa(oracle, n_dims, seed);

	// Anchors = first 15 ids; the rest stream in batches of 20 (→ 20, 20, 5).
	const uint32_t a = 15, batch_size = 20;
	std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.begin() + a);
	std::vector<std::string> remaining(oracle.ids.begin() + a, oracle.ids.end());

	const auto provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};
	const ProgressivePcoaResult result =
	    RunProgressivePcoa(anchors, remaining, n_dims, batch_size, seed, /*n_threads=*/1, provider);

	// Every sample gets every axis exactly once.
	REQUIRE(result.d == n_dims);
	REQUIRE(result.coords.size() == static_cast<size_t>(n) * n_dims);
	REQUIRE(result.eigvals.size() == n_dims);
	REQUIRE(result.proportion_explained.size() == n_dims);

	// Reassemble into oracle id order so it lines up row-for-row with `full`.
	const std::vector<double> prog = AssembleInOracleOrder(oracle, result, n_dims);

	// The two ordinations describe the same points, so they must coincide up to a
	// similarity transform: the procrustes disparity is ~0 (only FSVD numerical
	// error) for an exactly-Euclidean oracle. Observed here ~1.6e-13; the 1e-6
	// bound leaves generous headroom for cross-platform fp32/FSVD rounding while
	// still asserting "progressive reproduces full PCoA to numerical precision".
	const double m2 = FullDisparity(full.data(), prog.data(), n, n_dims);
	INFO("progressive-vs-full procrustes M^2 = " << m2);
	REQUIRE(m2 < 1e-6);
}

TEST_CASE("progressive PCoA anchors are batch-size invariant; quality does not compound", "[progressive]") {
	// The reference frame is fixed by the anchors alone, and each batch is aligned
	// to it independently — so the anchor coordinates must be identical regardless
	// of how the remaining samples are batched, and overall accuracy must not
	// degrade with more (smaller) batches. This pins the "errors don't compound"
	// design claim.
	const uint32_t n = 80, d_true = 3, n_dims = 3, a = 20;
	const int seed = 7;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/999);
	const std::vector<double> full = FullPcoa(oracle, n_dims, seed);
	const std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.begin() + a);
	const std::vector<std::string> remaining(oracle.ids.begin() + a, oracle.ids.end());
	const auto provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};

	// Many small batches vs a single batch covering all 60 remaining samples.
	const ProgressivePcoaResult many =
	    RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/12, seed, /*n_threads=*/1, provider);
	const ProgressivePcoaResult one =
	    RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/1000, seed, /*n_threads=*/1, provider);

	const std::vector<double> prog_many = AssembleInOracleOrder(oracle, many, n_dims);
	const std::vector<double> prog_one = AssembleInOracleOrder(oracle, one, n_dims);

	// Both track the ground truth tightly, independent of batch size.
	REQUIRE(FullDisparity(full.data(), prog_many.data(), n, n_dims) < 1e-6);
	REQUIRE(FullDisparity(full.data(), prog_one.data(), n, n_dims) < 1e-6);

	// Anchor coordinates are byte-for-byte independent of batching: same anchors +
	// same seed → same reference PCoA → same standardized anchor frame.
	for (uint32_t i = 0; i < a; ++i) {
		const size_t base = static_cast<size_t>(i) * n_dims;
		for (uint32_t k = 0; k < n_dims; ++k) {
			REQUIRE(prog_many[base + k] == Approx(prog_one[base + k]).margin(1e-12));
		}
	}

	// eigenvalues/proportions come from the anchor reference PCoA — also batch-invariant.
	REQUIRE(many.eigvals.size() == n_dims);
	for (uint32_t k = 0; k < n_dims; ++k) {
		REQUIRE(many.eigvals[k] == Approx(one.eigvals[k]).margin(1e-9));
		REQUIRE(many.proportion_explained[k] == Approx(one.proportion_explained[k]).margin(1e-9));
	}
}

TEST_CASE("progressive PCoA with no remaining samples emits only the reference anchors", "[progressive]") {
	// Degenerate case: everything is an anchor. The result is just the standardized
	// reference ordination (the self-fit path), with no batch phase.
	const uint32_t n = 20, d_true = 3, n_dims = 3;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/3);
	const std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.end());
	const auto provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};

	const ProgressivePcoaResult result =
	    RunProgressivePcoa(anchors, {}, n_dims, /*batch_size=*/10, /*seed=*/1, /*n_threads=*/1, provider);
	REQUIRE(result.d == n_dims);
	REQUIRE(result.coords.size() == static_cast<size_t>(n) * n_dims);
	(void)AssembleInOracleOrder(oracle, result, n_dims); // every (sample, axis) present once
}

TEST_CASE("progressive PCoA fails loud on misuse", "[progressive]") {
	const Oracle oracle = MakeEuclideanOracle(30, 3, /*seed=*/5);
	const auto provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};
	const std::vector<std::string> ids = oracle.ids;

	SECTION("too few anchors for the requested dimensionality") {
		const std::vector<std::string> anchors(ids.begin(), ids.begin() + 3); // n_dims=3 needs >= 4
		const std::vector<std::string> remaining(ids.begin() + 3, ids.end());
		REQUIRE_THROWS_AS(RunProgressivePcoa(anchors, remaining, 3, 10, 1, 1, provider), std::invalid_argument);
	}
	SECTION("a sample is both an anchor and a batch sample") {
		const std::vector<std::string> anchors(ids.begin(), ids.begin() + 6);
		std::vector<std::string> remaining(ids.begin() + 4, ids.end()); // overlaps ids[4], ids[5]
		REQUIRE_THROWS_AS(RunProgressivePcoa(anchors, remaining, 3, 10, 1, 1, provider), std::invalid_argument);
	}
	SECTION("n_dims must be >= 1") {
		const std::vector<std::string> anchors(ids.begin(), ids.begin() + 6);
		const std::vector<std::string> remaining(ids.begin() + 6, ids.end());
		REQUIRE_THROWS_AS(RunProgressivePcoa(anchors, remaining, 0, 10, 1, 1, provider), std::invalid_argument);
	}
	SECTION("batch_size must be >= 1") {
		const std::vector<std::string> anchors(ids.begin(), ids.begin() + 6);
		const std::vector<std::string> remaining(ids.begin() + 6, ids.end());
		REQUIRE_THROWS_AS(RunProgressivePcoa(anchors, remaining, 3, 0, 1, 1, provider), std::invalid_argument);
	}
	SECTION("n_threads must be >= 1") {
		const std::vector<std::string> anchors(ids.begin(), ids.begin() + 6);
		const std::vector<std::string> remaining(ids.begin() + 6, ids.end());
		REQUIRE_THROWS_AS(RunProgressivePcoa(anchors, remaining, 3, 10, 1, 0, provider), std::invalid_argument);
	}
}
