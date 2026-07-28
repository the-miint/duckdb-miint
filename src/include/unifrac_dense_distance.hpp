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

// Incremental builder for a dense n×n symmetric row-major fp32 distance matrix
// (zero diagonal). Cells are written as entries arrive, so a streaming reader
// holds nothing per-row: the matrix + fill bitmap are the only O(n²) state.
// `BuildDenseDistanceMatrix` below is a thin wrapper over this for callers that
// already hold every entry in a vector.
//
// `ids` (size == n) is used only for error messages, and is referenced (not
// copied) — it must outlive the builder.
//
// Throws std::invalid_argument from the constructor on:
//   - n < 2                    — no unordered pairs; ordination/PERMANOVA undefined
//   - ids.size() != n          — caller contract violation
//   - n too large for size_t   — n*n would overflow the allocation (32-bit size_t)
// ...from Add() on:
//   - an index >= n            — caller contract violation
//   - a non-finite distance    — NaN/Inf poisons the downstream fp32 routines
//   - distance < 0             — nonsensical distance
//   - a nonzero self-distance  — the diagonal is definitionally zero
//   - a duplicate unordered pair set to a *different* fp32 value — ambiguous input
// ...and from Finish() on:
//   - an incomplete matrix     — names the first missing (ids[i], ids[j]) in
//                                row-major upper-triangle order
//
// The split is deliberate: everything knowable from one entry fails on that
// entry (so a streaming caller can attribute it to the row in hand), while
// completeness is only knowable once the stream ends.
//
// Self-pairs (i == j) with distance 0 are ignored: the diagonal is
// definitionally zero and does not count toward completeness.
class DenseDistanceMatrixBuilder {
public:
	DenseDistanceMatrixBuilder(uint32_t n, const std::vector<std::string> &ids) : n_(n), ids_(ids) {
		if (n < 2) {
			throw std::invalid_argument("a distance matrix requires at least 2 distinct samples (got " +
			                            std::to_string(n) + ")");
		}
		if (ids.size() != static_cast<size_t>(n)) {
			throw std::invalid_argument("id vector size (" + std::to_string(ids.size()) + ") does not match n (" +
			                            std::to_string(n) + ")");
		}
		// Guard n*n against wraparound on a 32-bit size_t (wasm32) before it
		// drives the allocation size and every cell index. On 64-bit this never
		// fires (uint32 squared fits in size_t); on 32-bit it turns a would-be
		// silent empty allocation + out-of-bounds write into a clean throw.
		const uint64_t nn64 = static_cast<uint64_t>(n) * static_cast<uint64_t>(n);
		if (nn64 > static_cast<uint64_t>(SIZE_MAX)) {
			throw std::invalid_argument("distance matrix too large: " + std::to_string(n) + " samples");
		}
		const size_t nn = static_cast<size_t>(nn64);
		matrix_.assign(nn, 0.0f);
		// Separate "filled" bitmap: a genuine off-diagonal distance of 0.0 is
		// indistinguishable from an unwritten cell in `matrix_` alone, so we
		// cannot use a value sentinel for the completeness / conflict checks.
		filled_.assign(nn, 0);
	}

	void Add(uint32_t i, uint32_t j, double distance) {
		if (i >= n_ || j >= n_) {
			throw std::invalid_argument("distance entry index out of range for n=" + std::to_string(n_));
		}
		if (!std::isfinite(distance)) {
			throw std::invalid_argument("non-finite distance between '" + ids_[i] + "' and '" + ids_[j] + "'");
		}
		if (distance < 0.0) {
			throw std::invalid_argument("negative distance " + std::to_string(distance) + " between '" + ids_[i] +
			                            "' and '" + ids_[j] + "'");
		}
		if (i == j) {
			// The diagonal is definitionally zero; a nonzero self-distance is a
			// contradiction (e.g. a botched join counting a sample against
			// itself) and must be surfaced, not silently dropped.
			if (distance != 0.0) {
				throw std::invalid_argument("nonzero self-distance " + std::to_string(distance) + " for sample '" +
				                            ids_[i] + "'");
			}
			return;
		}
		const size_t ij = static_cast<size_t>(i) * n_ + j;
		const size_t ji = static_cast<size_t>(j) * n_ + i;
		const float v = static_cast<float>(distance);
		if (filled_[ij]) {
			if (matrix_[ij] != v) {
				throw std::invalid_argument("conflicting distances for pair ('" + ids_[i] + "', '" + ids_[j] + "')");
			}
			return; // identical duplicate — no-op
		}
		matrix_[ij] = v;
		matrix_[ji] = v;
		filled_[ij] = 1;
		filled_[ji] = 1;
		++pairs_filled_;
	}

	// Validates completeness and hands over the matrix. The fill bitmap is
	// released here, so the caller is left holding only the matrix itself.
	std::vector<float> Finish() {
		const uint64_t expected = static_cast<uint64_t>(n_) * static_cast<uint64_t>(n_ - 1) / 2;
		if (pairs_filled_ < expected) {
			for (uint32_t i = 0; i + 1 < n_; ++i) {
				for (uint32_t j = i + 1; j < n_; ++j) {
					if (!filled_[static_cast<size_t>(i) * n_ + j]) {
						throw std::invalid_argument("incomplete distance matrix: no distance provided for pair ('" +
						                            ids_[i] + "', '" + ids_[j] + "')");
					}
				}
			}
		}
		filled_.clear();
		filled_.shrink_to_fit();
		return std::move(matrix_);
	}

private:
	uint32_t n_;
	const std::vector<std::string> &ids_;
	std::vector<float> matrix_;
	std::vector<char> filled_;
	uint64_t pairs_filled_ = 0;
};

// Build a dense n×n symmetric row-major fp32 distance matrix (zero diagonal)
// from resolved off-diagonal entries. Every unordered pair must be provided
// exactly once; identical duplicates (including the reversed (b,a) orientation)
// are tolerated as no-ops. `ids` (size == n) is used only for error messages.
//
// A wrapper over DenseDistanceMatrixBuilder — see it for the full error
// contract. Prefer the builder when entries can be streamed: this form requires
// the caller to hold all of them, which is O(n²) heap on top of the matrix.
inline std::vector<float> BuildDenseDistanceMatrix(const std::vector<DistanceEntry> &entries, uint32_t n,
                                                   const std::vector<std::string> &ids) {
	DenseDistanceMatrixBuilder builder(n, ids);
	for (const auto &e : entries) {
		builder.Add(e.i, e.j, e.distance);
	}
	return builder.Finish();
}

} // namespace miint::unifrac
