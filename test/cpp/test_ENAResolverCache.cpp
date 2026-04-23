#include <catch2/catch_test_macros.hpp>
#include "ena_resolver_cache.hpp"
#include <atomic>
#include <thread>
#include <vector>

static miint::ENARunInfo MakeRun(const std::string &acc) {
	miint::ENARunInfo r;
	r.run_accession = acc;
	return r;
}

TEST_CASE("ENAResolverCache: miss then hit", "[ena][cache]") {
	miint::ENAResolverCache cache(4);
	miint::ENAResolverCache::Key key {"ERR1", "auto"};

	std::vector<miint::ENARunInfo> out;
	CHECK_FALSE(cache.Get(key, out));

	std::vector<miint::ENARunInfo> value;
	value.push_back(MakeRun("ERR1"));
	cache.Put(key, value);

	CHECK(cache.Get(key, out));
	REQUIRE(out.size() == 1);
	CHECK(out[0].run_accession == "ERR1");
}

TEST_CASE("ENAResolverCache: LRU eviction at capacity", "[ena][cache]") {
	miint::ENAResolverCache cache(2);

	cache.Put({"ERR1", "auto"}, {MakeRun("ERR1")});
	cache.Put({"ERR2", "auto"}, {MakeRun("ERR2")});
	CHECK(cache.Size() == 2);

	// Touch ERR1 to make it MRU; ERR2 becomes LRU
	std::vector<miint::ENARunInfo> out;
	CHECK(cache.Get({"ERR1", "auto"}, out));

	// Insert ERR3 → ERR2 should be evicted
	cache.Put({"ERR3", "auto"}, {MakeRun("ERR3")});
	CHECK(cache.Size() == 2);
	CHECK(cache.Get({"ERR1", "auto"}, out));
	CHECK(cache.Get({"ERR3", "auto"}, out));
	CHECK_FALSE(cache.Get({"ERR2", "auto"}, out));
}

TEST_CASE("ENAResolverCache: key includes prefer_format", "[ena][cache]") {
	miint::ENAResolverCache cache(4);
	cache.Put({"ERR1", "auto"}, {MakeRun("ERR1-auto")});
	cache.Put({"ERR1", "sff"}, {MakeRun("ERR1-sff")});

	std::vector<miint::ENARunInfo> out;
	REQUIRE(cache.Get({"ERR1", "auto"}, out));
	CHECK(out[0].run_accession == "ERR1-auto");

	REQUIRE(cache.Get({"ERR1", "sff"}, out));
	CHECK(out[0].run_accession == "ERR1-sff");

	CHECK_FALSE(cache.Get({"ERR1", "fastq"}, out));
}

TEST_CASE("ENAResolverCache: negative caching (empty vector) returns hit", "[ena][cache]") {
	miint::ENAResolverCache cache(4);
	miint::ENAResolverCache::Key key {"ERR_MISSING", "auto"};

	cache.Put(key, std::vector<miint::ENARunInfo> {}); // negative cache entry

	std::vector<miint::ENARunInfo> out;
	out.push_back(MakeRun("stale")); // ensure Get overwrites
	CHECK(cache.Get(key, out));
	CHECK(out.empty());
}

TEST_CASE("ENAResolverCache: replace existing key updates value and bumps MRU", "[ena][cache]") {
	miint::ENAResolverCache cache(2);
	cache.Put({"ERR1", "auto"}, {MakeRun("v1")});
	cache.Put({"ERR2", "auto"}, {MakeRun("ERR2")});

	// Overwrite ERR1 with v2; this bumps ERR1 to MRU (so ERR2 is LRU)
	cache.Put({"ERR1", "auto"}, {MakeRun("v2")});

	// Evict ERR2 by inserting ERR3
	cache.Put({"ERR3", "auto"}, {MakeRun("ERR3")});

	std::vector<miint::ENARunInfo> out;
	REQUIRE(cache.Get({"ERR1", "auto"}, out));
	CHECK(out[0].run_accession == "v2");
	CHECK(cache.Get({"ERR3", "auto"}, out));
	CHECK_FALSE(cache.Get({"ERR2", "auto"}, out));
}

TEST_CASE("ENAResolverCache: capacity=0 disables caching", "[ena][cache]") {
	miint::ENAResolverCache cache(0);
	cache.Put({"ERR1", "auto"}, {MakeRun("ERR1")});
	CHECK(cache.Size() == 0);
	std::vector<miint::ENARunInfo> out;
	CHECK_FALSE(cache.Get({"ERR1", "auto"}, out));
}

TEST_CASE("ENAResolverCache: concurrent Put/Get (thread safety sanity)", "[ena][cache]") {
	miint::ENAResolverCache cache(128);
	constexpr int kThreads = 8;
	constexpr int kIters = 2000;

	std::atomic<int> failures {0};
	std::vector<std::thread> ts;
	for (int t = 0; t < kThreads; t++) {
		ts.emplace_back([&, t]() {
			for (int i = 0; i < kIters; i++) {
				std::string acc = "ERR" + std::to_string((t * kIters + i) % 50); // 50 distinct keys
				miint::ENAResolverCache::Key key {acc, "auto"};
				std::vector<miint::ENARunInfo> out;
				if (cache.Get(key, out)) {
					if (!out.empty() && out[0].run_accession != acc) {
						failures.fetch_add(1);
					}
				} else {
					cache.Put(key, {MakeRun(acc)});
				}
			}
		});
	}
	for (auto &th : ts) {
		th.join();
	}
	CHECK(failures.load() == 0);
}

TEST_CASE("ENAResolverCache: adversarial eviction under contention", "[ena][cache]") {
	// Capacity < live-key count forces constant eviction. 8 threads hammering 8 keys
	// against a 4-slot cache stresses the Put eviction path (order.pop_back +
	// map.erase + order.push_front + map.emplace) under concurrent access.
	miint::ENAResolverCache cache(4);
	constexpr int kThreads = 8;
	constexpr int kKeys = 8;
	constexpr int kIters = 5000;

	std::atomic<int> corruption {0};
	std::vector<std::thread> ts;
	for (int t = 0; t < kThreads; t++) {
		ts.emplace_back([&, t]() {
			for (int i = 0; i < kIters; i++) {
				std::string acc = "ACC" + std::to_string((t + i) % kKeys);
				miint::ENAResolverCache::Key key {acc, "auto"};
				std::vector<miint::ENARunInfo> out;
				if (cache.Get(key, out)) {
					// A hit must return a value matching its key (no cross-key contamination).
					if (!out.empty() && out[0].run_accession != acc) {
						corruption.fetch_add(1);
					}
				} else {
					cache.Put(key, {MakeRun(acc)});
				}
			}
		});
	}
	for (auto &th : ts) {
		th.join();
	}
	CHECK(corruption.load() == 0);
	// Cache must remain within capacity regardless of the race pattern.
	CHECK(cache.Size() <= 4);
}
