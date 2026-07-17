#pragma once

#include <cstdint>

namespace miint::unifrac {

// Row-major strict-upper-triangle iterator over an n×n symmetric matrix.
//
// Yields index pairs (i, j) with 0 <= i < j < n, in order:
//   (0,1), (0,2), ..., (0,n-1), (1,2), ..., (n-2,n-1)
// i.e. the condensed form of a symmetric distance matrix (each unordered pair
// once, diagonal excluded). Total n*(n-1)/2 pairs; for n < 2 there are none.
//
// A plain value type with no heap state, so it can be snapshotted and resumed
// across DuckDB execute-chunk boundaries just by copying it — which is exactly
// how unifrac_distances streams the upper triangle without materializing all
// O(n²) rows at once.
struct CondensedCursor {
	uint32_t n = 0;
	uint32_t i = 0;
	uint32_t j = 1;

	static CondensedCursor Begin(uint32_t n_samples) {
		CondensedCursor c;
		c.n = n_samples;
		c.i = 0;
		c.j = 1;
		return c;
	}

	// True when no pairs remain. `i + 1 >= n` (rather than `i >= n - 1`) avoids
	// unsigned underflow for n == 0.
	bool done() const {
		return n < 2 || i + 1 >= n;
	}

	// Advance to the next (i, j). Precondition: !done().
	void advance() {
		++j;
		if (j >= n) {
			++i;
			j = i + 1;
		}
	}
};

} // namespace miint::unifrac
