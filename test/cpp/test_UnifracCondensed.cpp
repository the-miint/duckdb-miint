#include <catch2/catch_test_macros.hpp>

#include "unifrac_condensed.hpp"

#include <cstdint>
#include <utility>
#include <vector>

using miint::unifrac::CondensedCursor;

namespace {

// Drive a cursor to exhaustion, collecting every (i, j) it yields.
std::vector<std::pair<uint32_t, uint32_t>> Drain(uint32_t n) {
	std::vector<std::pair<uint32_t, uint32_t>> pairs;
	auto c = CondensedCursor::Begin(n);
	while (!c.done()) {
		pairs.emplace_back(c.i, c.j);
		c.advance();
	}
	return pairs;
}

// Independent reference: the strict upper triangle in row-major order.
std::vector<std::pair<uint32_t, uint32_t>> Reference(uint32_t n) {
	std::vector<std::pair<uint32_t, uint32_t>> pairs;
	for (uint32_t i = 0; i + 1 < n; ++i) {
		for (uint32_t j = i + 1; j < n; ++j) {
			pairs.emplace_back(i, j);
		}
	}
	return pairs;
}

} // namespace

TEST_CASE("CondensedCursor yields nothing for n < 2", "[unifrac][condensed]") {
	// A distance "matrix" with 0 or 1 samples has no unordered pairs. done()
	// must be true from the start so the streaming loop emits an empty result
	// rather than underflowing on n - 1.
	REQUIRE(CondensedCursor::Begin(0).done());
	REQUIRE(CondensedCursor::Begin(1).done());
	REQUIRE(Drain(0).empty());
	REQUIRE(Drain(1).empty());
}

TEST_CASE("CondensedCursor enumerates the strict upper triangle in order", "[unifrac][condensed]") {
	// The emitted sequence is exactly (0,1),(0,2),...,(0,n-1),(1,2),... — one
	// unordered pair each, diagonal excluded — which is what makes sample_a <
	// sample_b hold for lexicographically sorted ids.
	for (uint32_t n : {2u, 3u, 5u, 10u}) {
		const auto got = Drain(n);
		const auto want = Reference(n);
		REQUIRE(got == want);
		// n*(n-1)/2 pairs, all with i < j.
		REQUIRE(got.size() == static_cast<size_t>(n) * (n - 1) / 2);
		for (const auto &p : got) {
			REQUIRE(p.first < p.second);
		}
	}
}

TEST_CASE("CondensedCursor is resumable across a boundary by copying", "[unifrac][condensed]") {
	// The execute path snapshots the cursor between output chunks. Copying the
	// value type mid-stream and continuing must produce the identical tail —
	// this is the invariant that keeps multi-chunk streaming correct.
	const uint32_t n = 12; // 66 pairs
	auto full = Drain(n);

	auto c = CondensedCursor::Begin(n);
	std::vector<std::pair<uint32_t, uint32_t>> rejoined;
	// Emit 20 pairs, snapshot, then resume from the copy.
	for (int k = 0; k < 20 && !c.done(); ++k) {
		rejoined.emplace_back(c.i, c.j);
		c.advance();
	}
	CondensedCursor snapshot = c; // copy
	while (!snapshot.done()) {
		rejoined.emplace_back(snapshot.i, snapshot.j);
		snapshot.advance();
	}
	REQUIRE(rejoined == full);
}
