#include "pick_anchors.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <unordered_map>
#include <utility>

namespace miint {

namespace {

// Shared argument checks. Both selectors are SQL-facing through pick_anchors, so
// the messages are prefixed with that name (the cluster_kmeans convention).
void ValidateCloud(const std::vector<double> &points, uint32_t n_points, uint32_t n_dims, uint32_t k) {
	if (n_points == 0) {
		throw std::invalid_argument("pick_anchors: no points");
	}
	if (n_dims == 0) {
		throw std::invalid_argument("pick_anchors: n_dims must be >= 1");
	}
	if (k < 1) {
		throw std::invalid_argument("pick_anchors: k must be >= 1 (got " + std::to_string(k) + ")");
	}
	if (k > n_points) {
		throw std::invalid_argument("pick_anchors: k (" + std::to_string(k) + ") must be <= number of points (" +
		                            std::to_string(n_points) + ")");
	}
	if (points.size() != static_cast<size_t>(n_points) * n_dims) {
		throw std::invalid_argument("pick_anchors: points size (" + std::to_string(points.size()) +
		                            ") does not match n_points * n_dims (" + std::to_string(n_points) + " * " +
		                            std::to_string(n_dims) + ")");
	}
	// Both rules order coordinates, and a NaN breaks the strict weak ordering that
	// std::sort requires -- undefined behavior, not merely a wrong answer. Infinite
	// coordinates would poison the centroid and every squared distance. Neither can
	// come out of a real ordination, so this is a fail-loud check on a broken input
	// rather than a supported mode.
	for (const double v : points) {
		if (!std::isfinite(v)) {
			throw std::invalid_argument("pick_anchors: coordinates must all be finite");
		}
	}
}

// FNV-1a 64 over the id bytes, finalized with splitmix64 mixed with the seed.
// Hand-rolled rather than std::hash because std::hash is implementation-defined:
// an anchor set has to reproduce across platforms and standard libraries. The
// splitmix64 finalizer avalanches, so adjacent seeds give unrelated orders.
uint64_t SaltedIdHash(const std::string &id, int64_t seed) {
	uint64_t h = 14695981039346656037ull; // FNV-1a offset basis
	for (const char c : id) {
		h ^= static_cast<uint64_t>(static_cast<unsigned char>(c));
		h *= 1099511628211ull; // FNV-1a prime
	}
	h ^= static_cast<uint64_t>(seed);
	h += 0x9e3779b97f4a7c15ull;
	h = (h ^ (h >> 30)) * 0xbf58476d1ce4e5b9ull;
	h = (h ^ (h >> 27)) * 0x94d049bb133111ebull;
	return h ^ (h >> 31);
}

// DuckDB's ntile(b) semantics over n ordered rows: the first (n % b) buckets take
// ceil(n/b) rows and the rest take floor(n/b); when b > n every row is its own
// bucket. Reproduced exactly so the compiled selector bins the same way as the
// SQL prototype the 1.74x measurement came from.
uint32_t NtileBucket(uint32_t position, uint32_t n_points, uint32_t n_bins) {
	const uint32_t q = n_points / n_bins;
	const uint32_t r = n_points % n_bins;
	if (q == 0) {
		// n_points < n_bins, so position < n_points < n_bins: still a valid bin index.
		return position;
	}
	const uint32_t big_rows = r * (q + 1);
	return position < big_rows ? position / (q + 1) : r + (position - big_rows) / q;
}

} // namespace

std::vector<uint32_t> SelectStratified(const std::vector<double> &points, const std::vector<std::string> &ids,
                                       uint32_t n_points, uint32_t n_dims, uint32_t k, uint32_t n_bins, int64_t seed) {
	ValidateCloud(points, n_points, n_dims, k);
	if (ids.size() != static_cast<size_t>(n_points)) {
		throw std::invalid_argument("pick_anchors: ids size (" + std::to_string(ids.size()) +
		                            ") does not match n_points (" + std::to_string(n_points) + ")");
	}
	if (n_bins < 1) {
		throw std::invalid_argument("pick_anchors: n_bins must be >= 1 (got " + std::to_string(n_bins) + ")");
	}
	// The stratum id is the bin tuple in base n_bins, so the grid has to be
	// indexable. Refuse rather than let it wrap and alias unrelated cells together.
	uint64_t cells = 1;
	for (uint32_t j = 0; j < n_dims; ++j) {
		if (cells > std::numeric_limits<uint64_t>::max() / n_bins) {
			throw std::invalid_argument("pick_anchors: n_bins^n_dims strata (" + std::to_string(n_bins) + "^" +
			                            std::to_string(n_dims) +
			                            ") exceeds what a 64-bit cell index can hold; lower n_bins or n_dims");
		}
		cells *= n_bins;
	}

	// Stratum = the equal-frequency grid cell, accumulated one axis at a time in
	// base n_bins. stable_sort so coordinate ties fall in ascending row order,
	// which makes the binning a function of the data alone.
	std::vector<uint64_t> stratum(n_points, 0);
	std::vector<uint32_t> order(n_points);
	for (uint32_t j = 0; j < n_dims; ++j) {
		std::iota(order.begin(), order.end(), 0u);
		std::stable_sort(order.begin(), order.end(), [&](uint32_t a, uint32_t b) {
			return points[static_cast<size_t>(a) * n_dims + j] < points[static_cast<size_t>(b) * n_dims + j];
		});
		for (uint32_t p = 0; p < n_points; ++p) {
			stratum[order[p]] = stratum[order[p]] * n_bins + NtileBucket(p, n_points, n_bins);
		}
	}

	std::unordered_map<uint64_t, std::vector<uint32_t>> members;
	for (uint32_t i = 0; i < n_points; ++i) {
		members[stratum[i]].push_back(i);
	}

	// Within a stratum: rank by salted hash, carrying NO geometric preference --
	// that unbiasedness is the whole mechanism (see the header). Across strata:
	// order by rank/size, so taking a prefix draws proportionally to stratum size.
	struct Keyed {
		double key;
		uint64_t hash;
		uint32_t index;
	};
	std::vector<Keyed> keyed;
	keyed.reserve(n_points);
	std::vector<std::pair<uint64_t, uint32_t>> ranked;
	for (auto &entry : members) {
		auto &mem = entry.second;
		ranked.clear();
		ranked.reserve(mem.size());
		for (const auto i : mem) {
			ranked.emplace_back(SaltedIdHash(ids[i], seed), i);
		}
		// Index breaks a hash collision, so the order is total and reproducible.
		std::sort(ranked.begin(), ranked.end());
		const double size = static_cast<double>(mem.size());
		for (size_t rank = 0; rank < ranked.size(); ++rank) {
			keyed.push_back({static_cast<double>(rank + 1) / size, ranked[rank].first, ranked[rank].second});
		}
	}

	// Two strata of equal size hold identical keys; the hash breaks that tie, so
	// which stratum yields first is unbiased too, and the index makes it total.
	std::sort(keyed.begin(), keyed.end(), [](const Keyed &a, const Keyed &b) {
		if (a.key != b.key) {
			return a.key < b.key;
		}
		if (a.hash != b.hash) {
			return a.hash < b.hash;
		}
		return a.index < b.index;
	});

	std::vector<uint32_t> chosen;
	chosen.reserve(k);
	for (uint32_t i = 0; i < k; ++i) {
		chosen.push_back(keyed[i].index);
	}
	return chosen;
}

std::vector<uint32_t> SelectFarthestPoint(const std::vector<double> &points, uint32_t n_points, uint32_t n_dims,
                                          uint32_t k) {
	ValidateCloud(points, n_points, n_dims, k);

	const auto squared_distance = [&](uint32_t a, uint32_t b) {
		double total = 0.0;
		for (uint32_t j = 0; j < n_dims; ++j) {
			const double delta =
			    points[static_cast<size_t>(a) * n_dims + j] - points[static_cast<size_t>(b) * n_dims + j];
			total += delta * delta;
		}
		return total;
	};

	// Rank 0 = farthest from the centroid, which in Euclidean space is exactly the
	// largest total squared distance to all other points (see the header) -- the
	// same "most peripheral" start a distance-matrix implementation reaches by row
	// sums, at O(n*d) instead of O(n^2*d).
	std::vector<double> centroid(n_dims, 0.0);
	for (uint32_t i = 0; i < n_points; ++i) {
		for (uint32_t j = 0; j < n_dims; ++j) {
			centroid[j] += points[static_cast<size_t>(i) * n_dims + j];
		}
	}
	for (auto &value : centroid) {
		value /= static_cast<double>(n_points);
	}
	uint32_t first = 0;
	double best_from_centroid = -1.0;
	for (uint32_t i = 0; i < n_points; ++i) {
		double total = 0.0;
		for (uint32_t j = 0; j < n_dims; ++j) {
			const double delta = points[static_cast<size_t>(i) * n_dims + j] - centroid[j];
			total += delta * delta;
		}
		if (total > best_from_centroid) { // strict >: ties keep the lower index
			best_from_centroid = total;
			first = i;
		}
	}

	// `min_d2[i]` = squared distance from i to its nearest already-chosen point.
	// Chosen points are set to -1 so they can never be picked twice.
	std::vector<double> min_d2(n_points, std::numeric_limits<double>::infinity());
	std::vector<uint32_t> chosen;
	chosen.reserve(k);
	const auto take = [&](uint32_t pick) {
		chosen.push_back(pick);
		min_d2[pick] = -1.0;
		for (uint32_t i = 0; i < n_points; ++i) {
			if (min_d2[i] < 0.0) {
				continue;
			}
			const double d2 = squared_distance(i, pick);
			if (d2 < min_d2[i]) {
				min_d2[i] = d2;
			}
		}
	};
	take(first);

	while (chosen.size() < k) {
		uint32_t best_i = 0;
		double best_d2 = -1.0;
		for (uint32_t i = 0; i < n_points; ++i) {
			if (min_d2[i] > best_d2) { // strict >: ties keep the lower index
				best_d2 = min_d2[i];
				best_i = i;
			}
		}
		// best_d2 < 0 would mean every point is already chosen, which k <= n_points
		// makes unreachable.
		take(best_i);
	}
	return chosen;
}

} // namespace miint
