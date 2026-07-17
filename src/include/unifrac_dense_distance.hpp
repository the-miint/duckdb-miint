#pragma once
//
// Pure (DuckDB-free) construction of a dense UniFrac-style distance matrix from
// a condensed COO distance table (the `sample_a, sample_b, distance` shape that
// `unifrac_distances` emits, or any precomputed Bray-Curtis/Jaccard relation).
//
// Kept header-only and dependency-free so the Catch2 unit-test target links it
// without the duckdb library — mirroring unifrac_condensed.hpp. The DuckDB-aware
// reader (ReadDistanceTable, unifrac_function_common.cpp) resolves ids to
// indices, then delegates the fill + validation here.

#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <string>
#include <vector>

namespace miint::unifrac {

// A single off-diagonal distance, addressed by indices into a caller-owned,
// lexicographically-sorted `ids` vector. `distance` is carried as double (the
// SQL wire type) and narrowed to fp32 on storage to match libssu's matrices.
struct DistanceEntry {
	uint32_t i;
	uint32_t j;
	double distance;
};

// Build a dense n×n symmetric row-major fp32 distance matrix (zero diagonal)
// from resolved off-diagonal entries. Every unordered pair must be provided
// exactly once; identical duplicates (including the reversed (b,a) orientation)
// are tolerated as no-ops. `ids` (size == n) is used only for error messages.
//
// `ids` must have exactly `n` entries and every `DistanceEntry` index must be
// < n — both are enforced (this header is DuckDB-free and callable directly, so
// a bad index must become a clean throw, never an out-of-bounds access).
//
// Throws std::invalid_argument on:
//   - n < 2                    — no unordered pairs; ordination/PERMANOVA undefined
//   - ids.size() != n          — caller contract violation
//   - n too large for size_t   — n*n would overflow the allocation (32-bit size_t)
//   - an index >= n            — caller contract violation
//   - a non-finite distance    — NaN/Inf poisons the downstream fp32 routines
//   - distance < 0             — nonsensical distance
//   - a nonzero self-distance  — the diagonal is definitionally zero
//   - a duplicate unordered pair set to a *different* fp32 value — ambiguous input
//   - an incomplete matrix     — names the first missing (ids[i], ids[j]) in
//                                row-major upper-triangle order
//
// Self-pairs (i == j) with distance 0 are ignored: the diagonal is
// definitionally zero and does not count toward completeness.
inline std::vector<float> BuildDenseDistanceMatrix(const std::vector<DistanceEntry> &entries, uint32_t n,
                                                   const std::vector<std::string> &ids) {
	if (n < 2) {
		throw std::invalid_argument("a distance matrix requires at least 2 distinct samples (got " + std::to_string(n) +
		                            ")");
	}
	if (ids.size() != static_cast<size_t>(n)) {
		throw std::invalid_argument("id vector size (" + std::to_string(ids.size()) + ") does not match n (" +
		                            std::to_string(n) + ")");
	}
	// Guard n*n against wraparound on a 32-bit size_t (wasm32) before it drives
	// the allocation size and every cell index. On 64-bit this never fires
	// (uint32 squared fits in size_t); on 32-bit it turns a would-be silent
	// empty allocation + out-of-bounds write into a clean throw.
	const uint64_t nn64 = static_cast<uint64_t>(n) * static_cast<uint64_t>(n);
	if (nn64 > static_cast<uint64_t>(SIZE_MAX)) {
		throw std::invalid_argument("distance matrix too large: " + std::to_string(n) + " samples");
	}
	const size_t nn = static_cast<size_t>(nn64);
	std::vector<float> matrix(nn, 0.0f);
	// Separate "filled" bitmap: a genuine off-diagonal distance of 0.0 is
	// indistinguishable from an unwritten cell in `matrix` alone, so we cannot
	// use a value sentinel for the completeness / conflict checks.
	std::vector<char> filled(nn, 0);
	uint64_t pairs_filled = 0;

	for (const auto &e : entries) {
		if (e.i >= n || e.j >= n) {
			throw std::invalid_argument("distance entry index out of range for n=" + std::to_string(n));
		}
		if (!std::isfinite(e.distance)) {
			throw std::invalid_argument("non-finite distance between '" + ids[e.i] + "' and '" + ids[e.j] + "'");
		}
		if (e.distance < 0.0) {
			throw std::invalid_argument("negative distance " + std::to_string(e.distance) + " between '" + ids[e.i] +
			                            "' and '" + ids[e.j] + "'");
		}
		if (e.i == e.j) {
			// The diagonal is definitionally zero; a nonzero self-distance is a
			// contradiction (e.g. a botched join counting a sample against
			// itself) and must be surfaced, not silently dropped.
			if (e.distance != 0.0) {
				throw std::invalid_argument("nonzero self-distance " + std::to_string(e.distance) + " for sample '" +
				                            ids[e.i] + "'");
			}
			continue;
		}
		const size_t ij = static_cast<size_t>(e.i) * n + e.j;
		const size_t ji = static_cast<size_t>(e.j) * n + e.i;
		const float v = static_cast<float>(e.distance);
		if (filled[ij]) {
			if (matrix[ij] != v) {
				throw std::invalid_argument("conflicting distances for pair ('" + ids[e.i] + "', '" + ids[e.j] + "')");
			}
			continue; // identical duplicate — no-op
		}
		matrix[ij] = v;
		matrix[ji] = v;
		filled[ij] = 1;
		filled[ji] = 1;
		++pairs_filled;
	}

	const uint64_t expected = static_cast<uint64_t>(n) * static_cast<uint64_t>(n - 1) / 2;
	if (pairs_filled < expected) {
		for (uint32_t i = 0; i + 1 < n; ++i) {
			for (uint32_t j = i + 1; j < n; ++j) {
				if (!filled[static_cast<size_t>(i) * n + j]) {
					throw std::invalid_argument("incomplete distance matrix: no distance provided for pair ('" +
					                            ids[i] + "', '" + ids[j] + "')");
				}
			}
		}
	}
	return matrix;
}

} // namespace miint::unifrac
