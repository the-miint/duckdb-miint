#include <catch2/catch_test_macros.hpp>

#include <chrono>
#include <cmath>
#include <condition_variable>
#include <cstdint>
#include <mutex>
#include <random>
#include <set>
#include <thread>
#include <vector>

#ifndef __EMSCRIPTEN__
#include <omp.h>
#endif

#include "distance.h"   // skbb_permanova_fp32
#include "ordination.h" // skbb_pcoa_fsvd_fp32

#include "unifrac_omp_scope.hpp"

using miint::unifrac::OmpThreadPin;
using miint::unifrac::SkbbCallScope;

namespace {

// An exact Euclidean distance matrix over n random points — a well-conditioned
// input for both PCoA and PERMANOVA, built the same way test_ProgressivePcoa's
// oracle builds its own.
std::vector<float> EuclideanDm(uint32_t n, uint32_t d_true, uint64_t seed) {
	std::mt19937_64 rng(seed);
	std::uniform_real_distribution<double> u(-5.0, 5.0);
	std::vector<double> pts(static_cast<size_t>(n) * d_true);
	for (auto &v : pts) {
		v = u(rng);
	}
	std::vector<float> dm(static_cast<size_t>(n) * n, 0.0f);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j) {
			double acc = 0.0;
			for (uint32_t k = 0; k < d_true; ++k) {
				const double diff = pts[i * d_true + k] - pts[j * d_true + k];
				acc += diff * diff;
			}
			const float dist = static_cast<float>(std::sqrt(acc));
			dm[static_cast<size_t>(i) * n + j] = dist;
			dm[static_cast<size_t>(j) * n + i] = dist;
		}
	}
	return dm;
}

// Run `body` on `n_threads` threads that are all inside it at once, so a race
// between them has a chance to happen rather than being serialized by luck.
// Threads rendezvous with a bounded wait: a timeout degrades the test to "less
// overlap than intended", never to a hang.
void RunConcurrently(size_t n_threads, const std::function<void(size_t)> &body) {
	std::mutex mu;
	std::condition_variable cv;
	size_t arrived = 0;
	std::vector<std::thread> pool;
	pool.reserve(n_threads);
	for (size_t t = 0; t < n_threads; ++t) {
		pool.emplace_back([&, t]() {
			{
				std::unique_lock<std::mutex> lk(mu);
				++arrived;
				cv.notify_all();
				cv.wait_for(lk, std::chrono::seconds(5), [&] { return arrived == n_threads; });
			}
			body(t);
		});
	}
	for (auto &th : pool) {
		th.join();
	}
}

} // namespace

TEST_CASE("OmpThreadPin pins the calling thread, not the process", "[omp]") {
	// WHY THIS IS THE FOUNDATION: every concurrent-compute decision here rests on
	// omp_set_num_threads affecting only the calling thread (OpenMP specifies it as
	// the current task's nthreads-var ICV). If the runtime we actually link treated
	// it as process-global, two concurrent pinned computes would silently overwrite
	// each other's fan-out — no crash, no wrong answer, just unexplained thread
	// counts. So pin two live threads to different values and require each to see
	// its own.
	int observed_one = -1, observed_four = -1;
	std::thread a([&]() {
		OmpThreadPin pin(1);
		std::this_thread::sleep_for(std::chrono::milliseconds(50)); // overlap with b
		observed_one = omp_get_max_threads();
	});
	std::thread b([&]() {
		OmpThreadPin pin(4);
		std::this_thread::sleep_for(std::chrono::milliseconds(50));
		observed_four = omp_get_max_threads();
	});
	a.join();
	b.join();
	REQUIRE(observed_one == 1);
	REQUIRE(observed_four == 4);
}

TEST_CASE("SkbbCallScope hands out a seed skbb can use without its global RNG", "[omp]") {
	// The whole reason concurrent skbb calls are safe is that a NON-NEGATIVE seed
	// keeps skbb on a per-call generator; seed < 0 sends it to a process-global
	// mt19937 (scikit-bio-binaries util/rand.cpp) that concurrent callers would
	// race on. This scope is what makes that unmissable: you cannot obtain the
	// thread pin without also obtaining a usable seed.
	SECTION("an explicit seed is passed through untouched, so seeded runs stay reproducible") {
		SkbbCallScope scope(1, 42);
		REQUIRE(scope.seed() == 42);
		REQUIRE(scope.seed() >= 0);
	}
	SECTION("an unseeded call gets a non-negative seed of its own, and they vary") {
		// Varying is the point: "no seed" must keep meaning "not reproducible", not
		// silently become a fixed seed.
		std::set<int> seeds;
		for (int i = 0; i < 8; ++i) {
			SkbbCallScope scope(1, -1);
			REQUIRE(scope.seed() >= 0);
			seeds.insert(scope.seed());
		}
		REQUIRE(seeds.size() > 1);
	}
	SECTION("it pins the calling thread's fan-out like a bare pin does") {
		SkbbCallScope scope(3, 1);
		REQUIRE(omp_get_max_threads() == 3);
	}
}

TEST_CASE("concurrent seeded skbb PCoA reproduces the serial result exactly", "[omp]") {
	// The capability claim: several ordinations may run at once in one process
	// without a process-wide lock. Bit-identical, not merely close — anything less
	// would mean the concurrent calls perturbed each other's numerics, which is
	// exactly the failure a shared RNG or a fought-over thread count would cause.
	const uint32_t n = 120, n_dims = 5;
	const int seed = 7;
	const std::vector<float> dm = EuclideanDm(n, 3, /*seed=*/8675309);

	const auto pcoa_once = [&](std::vector<float> &eig, std::vector<float> &samples, std::vector<float> &prop) {
		eig.assign(n_dims, 0.0f);
		samples.assign(static_cast<size_t>(n) * n_dims, 0.0f);
		prop.assign(n_dims, 0.0f);
		SkbbCallScope scope(1, seed);
		skbb_pcoa_fsvd_fp32(n, dm.data(), n_dims, scope.seed(), eig.data(), samples.data(), prop.data());
	};

	std::vector<float> base_eig, base_samples, base_prop;
	pcoa_once(base_eig, base_samples, base_prop);

	const size_t workers = 8;
	std::vector<std::vector<float>> eig(workers), samples(workers), prop(workers);
	RunConcurrently(workers, [&](size_t t) { pcoa_once(eig[t], samples[t], prop[t]); });

	for (size_t t = 0; t < workers; ++t) {
		REQUIRE(eig[t] == base_eig);
		REQUIRE(prop[t] == base_prop);
		REQUIRE(samples[t] == base_samples);
	}
}

TEST_CASE("concurrent seeded skbb PERMANOVA reproduces the serial result exactly", "[omp]") {
	// Same claim for the other skbb entry point we call. PERMANOVA permutes its
	// groupings from a per-call-seeded generator array, so identical seeds must give
	// identical F and p no matter how many run at once.
	const uint32_t n = 90;
	const unsigned int n_perm = 99;
	const int seed = 11;
	const std::vector<float> dm = EuclideanDm(n, 3, /*seed=*/24601);
	std::vector<unsigned int> grouping(n);
	for (uint32_t i = 0; i < n; ++i) {
		grouping[i] = i % 3; // three balanced groups
	}

	const auto permanova_once = [&](float &f_stat, float &p_value) {
		f_stat = 0.0f;
		p_value = 0.0f;
		SkbbCallScope scope(1, seed);
		skbb_permanova_fp32(n, dm.data(), grouping.data(), n_perm, scope.seed(), &f_stat, &p_value);
	};

	float base_f = 0.0f, base_p = 0.0f;
	permanova_once(base_f, base_p);

	const size_t workers = 8;
	std::vector<float> f(workers, 0.0f), p(workers, 0.0f);
	RunConcurrently(workers, [&](size_t t) { permanova_once(f[t], p[t]); });

	for (size_t t = 0; t < workers; ++t) {
		REQUIRE(f[t] == base_f);
		REQUIRE(p[t] == base_p);
	}
}
