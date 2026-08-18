#include <catch2/catch_test_macros.hpp>
#include <IntervalCompressor.hpp>

using namespace miint;

TEST_CASE("IntervalCompressor - empty state", "[interval_compressor]") {
	IntervalCompressor compressor;
	REQUIRE(compressor.Empty());
	REQUIRE(compressor.Size() == 0);
}

TEST_CASE("IntervalCompressor - single interval", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Compress();

	REQUIRE_FALSE(compressor.Empty());
	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 20);
}

TEST_CASE("IntervalCompressor - non-overlapping intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Add(100, 120);
	compressor.Add(200, 220);
	compressor.Compress();

	REQUIRE(compressor.Size() == 3);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 20);
	REQUIRE(compressor.starts[1] == 100);
	REQUIRE(compressor.stops[1] == 120);
	REQUIRE(compressor.starts[2] == 200);
	REQUIRE(compressor.stops[2] == 220);
}

TEST_CASE("IntervalCompressor - fully overlapping intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 110);
	compressor.Add(100, 220);
	compressor.Add(200, 300);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 300);
}

TEST_CASE("IntervalCompressor - partially overlapping intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 50);
	compressor.Add(40, 80);
	compressor.Add(100, 150);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 80);
	REQUIRE(compressor.starts[1] == 100);
	REQUIRE(compressor.stops[1] == 150);
}

TEST_CASE("IntervalCompressor - touching intervals (should merge)", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Add(20, 30);
	compressor.Add(30, 40);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 40);
}

TEST_CASE("IntervalCompressor - unsorted input", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(200, 220);
	compressor.Add(10, 20);
	compressor.Add(100, 120);
	compressor.Compress();

	REQUIRE(compressor.Size() == 3);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 20);
	REQUIRE(compressor.starts[1] == 100);
	REQUIRE(compressor.stops[1] == 120);
	REQUIRE(compressor.starts[2] == 200);
	REQUIRE(compressor.stops[2] == 220);
}

TEST_CASE("IntervalCompressor - nested intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 100);
	compressor.Add(20, 30);
	compressor.Add(40, 50);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 100);
}

TEST_CASE("IntervalCompressor - complex overlapping", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 30);
	compressor.Add(20, 50);
	compressor.Add(40, 60);
	compressor.Add(70, 80);
	compressor.Add(75, 90);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 60);
	REQUIRE(compressor.starts[1] == 70);
	REQUIRE(compressor.stops[1] == 90);
}

TEST_CASE("IntervalCompressor - multiple compressions", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Add(15, 25);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 25);

	compressor.Add(30, 40);
	compressor.Add(35, 45);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 25);
	REQUIRE(compressor.starts[1] == 30);
	REQUIRE(compressor.stops[1] == 45);
}

TEST_CASE("IntervalCompressor - identical intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Add(10, 20);
	compressor.Add(10, 20);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 20);
}

TEST_CASE("IntervalCompressor - zero-length intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 10);
	compressor.Add(20, 20);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 10);
	REQUIRE(compressor.starts[1] == 20);
	REQUIRE(compressor.stops[1] == 20);
}

TEST_CASE("IntervalCompressor - merge overlapping with touching", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(10, 20);
	compressor.Add(15, 30);
	compressor.Add(30, 40);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 10);
	REQUIRE(compressor.stops[0] == 40);
}

TEST_CASE("IntervalCompressor - negative coordinates", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(-100, -50);
	compressor.Add(-60, -40);
	compressor.Add(10, 20);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == -100);
	REQUIRE(compressor.stops[0] == -40);
	REQUIRE(compressor.starts[1] == 10);
	REQUIRE(compressor.stops[1] == 20);
}

TEST_CASE("IntervalCompressor - very large coordinates", "[interval_compressor]") {
	IntervalCompressor compressor;
	constexpr int64_t large = 9223372036854775000LL;
	compressor.Add(large, large + 100);
	compressor.Add(large + 50, large + 200);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == large);
	REQUIRE(compressor.stops[0] == large + 200);
}

TEST_CASE("IntervalCompressor - intervals starting at zero", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(0, 10);
	compressor.Add(5, 15);
	compressor.Add(20, 30);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
	REQUIRE(compressor.starts[0] == 0);
	REQUIRE(compressor.stops[0] == 15);
	REQUIRE(compressor.starts[1] == 20);
	REQUIRE(compressor.stops[1] == 30);
}

TEST_CASE("IntervalCompressor - inverted intervals (stop < start)", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(20, 10);
	compressor.Add(100, 50);
	compressor.Compress();

	REQUIRE(compressor.Size() == 2);
}

TEST_CASE("IntervalCompressor - many intervals merging into one", "[interval_compressor]") {
	IntervalCompressor compressor;
	for (int i = 0; i < 100; i++) {
		compressor.Add(i * 10, i * 10 + 15);
	}
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 0);
	REQUIRE(compressor.stops[0] == 1005);
}

TEST_CASE("IntervalCompressor - many non-overlapping intervals", "[interval_compressor]") {
	IntervalCompressor compressor;
	for (int i = 0; i < 100; i++) {
		compressor.Add(i * 100, i * 100 + 10);
	}
	compressor.Compress();

	REQUIRE(compressor.Size() == 100);
	REQUIRE(compressor.starts[0] == 0);
	REQUIRE(compressor.stops[0] == 10);
	REQUIRE(compressor.starts[99] == 9900);
	REQUIRE(compressor.stops[99] == 9910);
}

TEST_CASE("IntervalCompressor - automatic compression at threshold", "[interval_compressor]") {
	IntervalCompressor compressor;

	for (int i = 0; i < 999999; i++) {
		compressor.Add(i, i + 10);
	}

	REQUIRE(compressor.Size() == 999999);

	compressor.Add(999999, 1000009);

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == 0);
	REQUIRE(compressor.stops[0] == 1000009);
}

TEST_CASE("IntervalCompressor - compression is amortized, not per-row", "[interval_compressor]") {
	// The input shape Compress() cannot reclaim anything from: every interval disjoint.
	// The threshold is a FLOOR THAT MUST GROW, because a fixed one re-triggers on the very
	// next Add() once compression fails to get the state below it -- every subsequent Add()
	// then re-sorts the whole state, which is O(n^2 log n) overall (measured ~10 ms per row
	// past 1e6). The test above covers the opposite shape, where everything merges to one.
	//
	// Asserted as a bound on compression COUNT so it is deterministic rather than a flaky
	// wall-clock limit. A fixed floor gives 101 here; a doubling floor gives 1.
	const size_t n = IntervalCompressor::COMPRESS_THRESHOLD + 100;

	IntervalCompressor compressor;
	for (size_t i = 0; i < n; i++) {
		const int64_t start = static_cast<int64_t>(i) * 10;
		compressor.Add(start, start + 5);
	}

	REQUIRE(compressor.CompressionCount() <= 3);

	// The policy must not change the answer: n disjoint intervals stay n after a final
	// Compress(), and the total width is exactly 5 per interval.
	compressor.Compress();
	REQUIRE(compressor.Size() == n);
	int64_t total = 0;
	for (size_t i = 0; i < compressor.Size(); i++) {
		total += compressor.stops[i] - compressor.starts[i];
	}
	REQUIRE(total == static_cast<int64_t>(n) * 5);
}

TEST_CASE("IntervalCompressor - a loop of Add() bounds the state on its own", "[interval_compressor]") {
	// This is exactly what compress_intervals' Combine does: replay the source state's
	// intervals into the target through Add(). Add() checks the growing compression floor on
	// EVERY push, so when the loop ends the state is already bounded -- which is why Combine
	// must NOT follow it with an unconditional Compress(). Doing so reclaimed nothing on
	// disjoint input and cost O(m log m) per combine -- about 25% of the aggregate on 500k
	// rows / 50 groups, at either thread count -- and it could compress twice in a row when
	// the loop had just crossed the floor itself.
	//
	// The DuckDB glue is not reachable from this binary (it links no libduckdb), so this
	// pins the class-level invariant the removal rests on rather than the Combine call.
	IntervalCompressor target;
	IntervalCompressor source;

	// Disjoint, so compression can reclaim nothing -- the shape that made the old policy
	// pathological, and the shape compress_intervals is normally fed.
	for (int i = 0; i < 5000; i++) {
		source.Add(i * 10, i * 10 + 5);
	}
	REQUIRE(source.CompressionCount() == 0);

	for (size_t i = 0; i < source.Size(); i++) {
		target.Add(source.starts[i], source.stops[i]);
	}

	// Well under the floor, so no compression was needed and none happened.
	REQUIRE(target.CompressionCount() == 0);
	REQUIRE(target.Size() == 5000);

	// And the answer is intact once the caller does compress, which Finalize always does.
	target.Compress();
	REQUIRE(target.Size() == 5000);
	int64_t total = 0;
	for (size_t i = 0; i < target.Size(); i++) {
		total += target.stops[i] - target.starts[i];
	}
	REQUIRE(total == 5000 * 5);
}

TEST_CASE("IntervalCompressor - repeated merges stay bounded without an external Compress", "[interval_compressor]") {
	// Many combines, as a parallel GROUP BY produces. The state must stay bounded from
	// Add()'s floor alone, with the caller never compressing: this is what makes dropping
	// Combine's Compress() safe rather than merely faster.
	IntervalCompressor target;
	for (int round = 0; round < 40; round++) {
		for (int i = 0; i < 30000; i++) {
			const int64_t base = (static_cast<int64_t>(round) * 30000 + i) * 10;
			target.Add(base, base + 5);
		}
	}

	// 1.2M disjoint intervals: the floor is crossed once at 1e6 and then doubles past the
	// total, so the whole run costs a small constant number of compressions.
	REQUIRE(target.CompressionCount() <= 3);

	target.Compress();
	REQUIRE(target.Size() == 40 * 30000);
	int64_t total = 0;
	for (size_t i = 0; i < target.Size(); i++) {
		total += target.stops[i] - target.starts[i];
	}
	REQUIRE(total == 40LL * 30000LL * 5LL);
}

TEST_CASE("IntervalCompressor - mix of positive and negative", "[interval_compressor]") {
	IntervalCompressor compressor;
	compressor.Add(-50, -10);
	compressor.Add(-20, 20);
	compressor.Add(10, 50);
	compressor.Compress();

	REQUIRE(compressor.Size() == 1);
	REQUIRE(compressor.starts[0] == -50);
	REQUIRE(compressor.stops[0] == 50);
}
