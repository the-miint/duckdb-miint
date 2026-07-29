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

TEST_CASE("progressive PCoA reports each batch's anchor-overlap disparity", "[progressive]") {
	// WHY this exists: at the scale progressive PCoA is for, a user cannot check the
	// result against a full PCoA — so the run must carry its own quality evidence.
	// Each batch is placed into the reference frame by a procrustes fit on the
	// anchor overlap, and that fit's disparity is already computed internally. It is
	// the per-batch measure of "did this batch's local geometry agree with the
	// reference frame", and reporting it turns an unverifiable run into an
	// auditable one. On a Euclidean oracle every batch's geometry is exactly
	// consistent, so every disparity must be ~0.
	const uint32_t n = 80, d_true = 3, n_dims = 3, a = 20;
	const int seed = 7;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/4242);
	const std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.begin() + a);
	const std::vector<std::string> remaining(oracle.ids.begin() + a, oracle.ids.end());
	const auto provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};

	const ProgressivePcoaResult r =
	    RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/15, seed, /*n_threads=*/1, provider);

	// 60 remaining samples in batches of 15 → exactly 4 batches, reported in order.
	REQUIRE(r.batches.size() == 4);
	uint32_t counted = 0;
	for (size_t b = 0; b < r.batches.size(); ++b) {
		REQUIRE(r.batches[b].batch == static_cast<int32_t>(b));
		REQUIRE(r.batches[b].n_samples == 15);
		// Exact Euclidean geometry → the anchor overlap fits the reference frame
		// essentially perfectly. A nonzero value here is the signal a user acts on.
		REQUIRE(r.batches[b].anchor_m2 < 1e-9);
		counted += r.batches[b].n_samples;
	}
	REQUIRE(counted == remaining.size());
}

TEST_CASE("progressive PCoA batch disparity rises when a batch disagrees with the frame", "[progressive]") {
	// The complement of the test above: the diagnostic must be able to FAIL. A
	// disparity that is always ~0 would be decoration, not evidence. Here one
	// batch's block is served with its non-anchor distances corrupted (scrambled
	// against the anchors), so its local geometry genuinely contradicts the
	// reference frame — that batch's disparity must stand out from the rest.
	const uint32_t n = 80, d_true = 3, n_dims = 3, a = 20;
	const int seed = 7;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/4242);
	const std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.begin() + a);
	const std::vector<std::string> remaining(oracle.ids.begin() + a, oracle.ids.end());

	// Corrupt only the LAST batch's request (its first non-anchor id identifies it):
	// rotate the anchor block's rows so the anchors' mutual geometry is wrong for
	// that block alone.
	const std::string poisoned_first_id = remaining[45]; // batch 3 of 4 with batch_size=15
	const auto provider = [&](const std::vector<std::string> &req) {
		DistanceBlock blk = SliceBlock(oracle, req);
		if (!req.empty() && req[0] == poisoned_first_id) {
			const uint32_t m = static_cast<uint32_t>(req.size());
			// Swap two anchor rows/cols — a valid symmetric matrix, but one whose
			// anchor configuration no longer matches the reference block's.
			const uint32_t i = m - 1, j = m - 2;
			for (uint32_t c = 0; c < m; ++c) {
				std::swap(blk.matrix[static_cast<size_t>(i) * m + c], blk.matrix[static_cast<size_t>(j) * m + c]);
			}
			for (uint32_t rr = 0; rr < m; ++rr) {
				std::swap(blk.matrix[static_cast<size_t>(rr) * m + i], blk.matrix[static_cast<size_t>(rr) * m + j]);
			}
		}
		return blk;
	};

	const ProgressivePcoaResult r =
	    RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/15, seed, /*n_threads=*/1, provider);
	REQUIRE(r.batches.size() == 4);
	// The clean batches stay ~0; the poisoned one is orders of magnitude worse.
	REQUIRE(r.batches[0].anchor_m2 < 1e-9);
	REQUIRE(r.batches[1].anchor_m2 < 1e-9);
	REQUIRE(r.batches[2].anchor_m2 < 1e-9);
	REQUIRE(r.batches[3].anchor_m2 > 1e-6);
}

TEST_CASE("progressive PCoA announces each wave's requests and is unchanged by wave width", "[progressive]") {
	// WHY: the per-block provider for a stored relation costs one full scan of the
	// relation per batch — O(B · pairs). Waves exist so a provider can satisfy many
	// blocks from ONE pass. Two properties make that safe, and both are pinned here:
	//   1. every batch request is announced before it is fetched, in fetch order, so
	//      a provider may serve batch blocks purely from its wave cache;
	//   2. wave width is an I/O choice ONLY — coordinates are identical for any W,
	//      so tuning it can never change a scientific result.
	const uint32_t n = 80, d_true = 3, n_dims = 3, a = 20;
	const int seed = 7;
	const Oracle oracle = MakeEuclideanOracle(n, d_true, /*seed=*/31337);
	const std::vector<std::string> anchors(oracle.ids.begin(), oracle.ids.begin() + a);
	const std::vector<std::string> remaining(oracle.ids.begin() + a, oracle.ids.end());
	const std::vector<std::string> anchors_only = anchors;

	// A deliberately strict provider: batch blocks may ONLY come from the wave
	// cache. Anything else means the core asked for a block it never announced.
	std::vector<std::vector<std::string>> announced;
	size_t prefetch_calls = 0;
	std::vector<std::pair<std::vector<std::string>, DistanceBlock>> cache;
	const auto prefetch = [&](const std::vector<std::vector<std::string>> &requests) {
		++prefetch_calls;
		cache.clear();
		for (const auto &req : requests) {
			announced.push_back(req);
			cache.emplace_back(req, SliceBlock(oracle, req));
		}
	};
	const auto cache_only_provider = [&](const std::vector<std::string> &req) {
		for (const auto &entry : cache) {
			if (entry.first == req) {
				return entry.second;
			}
		}
		// The anchors-alone reference block is fetched once, before any wave.
		if (req == anchors_only) {
			return SliceBlock(oracle, req);
		}
		throw std::invalid_argument("provider asked for an unannounced block");
	};

	const auto plain_provider = [&oracle](const std::vector<std::string> &req) {
		return SliceBlock(oracle, req);
	};

	// 60 remaining in batches of 12 → 5 batches. W=2 → waves of 2,2,1.
	const ProgressivePcoaResult waved = RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/12, seed,
	                                                       /*n_threads=*/1, cache_only_provider, prefetch,
	                                                       /*wave_batches=*/2);
	REQUIRE(prefetch_calls == 3);
	REQUIRE(announced.size() == 5); // every batch announced exactly once...
	// ...and in fetch order: request k starts with the (k*12)th remaining sample.
	for (size_t k = 0; k < announced.size(); ++k) {
		REQUIRE(announced[k].front() == remaining[k * 12]);
	}

	// Same run with no waves at all (one block at a time, no prefetch) must produce
	// identical coordinates — W is not a scientific parameter.
	const ProgressivePcoaResult plain =
	    RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/12, seed, /*n_threads=*/1, plain_provider);
	const std::vector<double> waved_m = AssembleInOracleOrder(oracle, waved, n_dims);
	const std::vector<double> plain_m = AssembleInOracleOrder(oracle, plain, n_dims);
	REQUIRE(waved_m.size() == plain_m.size());
	for (size_t i = 0; i < plain_m.size(); ++i) {
		REQUIRE(waved_m[i] == Approx(plain_m[i]).margin(1e-12));
	}
	// Diagnostics too: same batches, same disparities.
	REQUIRE(waved.batches.size() == plain.batches.size());
	for (size_t b = 0; b < plain.batches.size(); ++b) {
		REQUIRE(waved.batches[b].batch == plain.batches[b].batch);
		REQUIRE(waved.batches[b].n_samples == plain.batches[b].n_samples);
		REQUIRE(waved.batches[b].anchor_m2 == Approx(plain.batches[b].anchor_m2).margin(1e-12));
	}

	// A wave wider than the batch count is legal: one announcement, everything in it.
	prefetch_calls = 0;
	announced.clear();
	const ProgressivePcoaResult one_wave = RunProgressivePcoa(anchors, remaining, n_dims, /*batch_size=*/12, seed,
	                                                          /*n_threads=*/1, cache_only_provider, prefetch,
	                                                          /*wave_batches=*/999);
	REQUIRE(prefetch_calls == 1);
	REQUIRE(announced.size() == 5);
	const std::vector<double> one_wave_m = AssembleInOracleOrder(oracle, one_wave, n_dims);
	for (size_t i = 0; i < plain_m.size(); ++i) {
		REQUIRE(one_wave_m[i] == Approx(plain_m[i]).margin(1e-12));
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
