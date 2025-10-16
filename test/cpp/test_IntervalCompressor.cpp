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
