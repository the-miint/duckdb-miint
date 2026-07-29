#include "cluster_upgma.hpp"

#include <limits>
#include <stdexcept>
#include <string>

namespace miint {

std::vector<UpgmaNode> UpgmaAverageLinkage(const std::vector<double> &dist, uint32_t n) {
	if (n < 2) {
		throw std::invalid_argument("cluster_upgma requires at least 2 samples (got " + std::to_string(n) + ")");
	}
	if (dist.size() != static_cast<size_t>(n) * static_cast<size_t>(n)) {
		throw std::invalid_argument("cluster_upgma: distance matrix size (" + std::to_string(dist.size()) +
		                            ") does not match n*n (n=" + std::to_string(n) + ")");
	}

	const uint32_t total = 2 * n - 1;
	std::vector<UpgmaNode> nodes(total);
	for (uint32_t i = 0; i < total; ++i) {
		nodes[i] = {i < n ? static_cast<int32_t>(i) : -1, -1, 0.0, 0.0};
	}

	// Cluster sizes (# original leaves) for the size-weighted (arithmetic-mean)
	// distance update that makes this UPGMA rather than WPGMA.
	std::vector<double> csize(total, 0.0);
	for (uint32_t i = 0; i < n; ++i) {
		csize[i] = 1.0;
	}

	// Working inter-cluster distances, indexed by node id (grows as merges add
	// internal nodes). total x total is small for the paper's scale (n <= ~90).
	std::vector<double> D(static_cast<size_t>(total) * total, 0.0);
	for (uint32_t i = 0; i < n; ++i) {
		for (uint32_t j = 0; j < n; ++j) {
			D[static_cast<size_t>(i) * total + j] = dist[static_cast<size_t>(i) * n + j];
		}
	}

	std::vector<char> active(total, 0);
	for (uint32_t i = 0; i < n; ++i) {
		active[i] = 1;
	}

	uint32_t next = n;
	for (uint32_t m = 0; m + 1 < n; ++m) {
		// Closest active pair (ties -> lowest (i, j), deterministic).
		double best = std::numeric_limits<double>::infinity();
		uint32_t ba = 0, bb = 0;
		for (uint32_t i = 0; i < next; ++i) {
			if (!active[i]) {
				continue;
			}
			for (uint32_t j = i + 1; j < next; ++j) {
				if (!active[j]) {
					continue;
				}
				const double dd = D[static_cast<size_t>(i) * total + j];
				if (dd < best) {
					best = dd;
					ba = i;
					bb = j;
				}
			}
		}

		const uint32_t node = next++;
		const double h = best / 2.0;
		csize[node] = csize[ba] + csize[bb];
		nodes[node].height = h;
		nodes[ba].parent = static_cast<int32_t>(node);
		nodes[ba].branch_length = h - nodes[ba].height;
		nodes[bb].parent = static_cast<int32_t>(node);
		nodes[bb].branch_length = h - nodes[bb].height;

		// Size-weighted arithmetic mean of the two merged clusters' distances to
		// every other active cluster.
		for (uint32_t c = 0; c < node; ++c) {
			if (!active[c] || c == ba || c == bb) {
				continue;
			}
			const double dnew = (csize[ba] * D[static_cast<size_t>(ba) * total + c] +
			                     csize[bb] * D[static_cast<size_t>(bb) * total + c]) /
			                    (csize[ba] + csize[bb]);
			D[static_cast<size_t>(node) * total + c] = dnew;
			D[static_cast<size_t>(c) * total + node] = dnew;
		}

		active[ba] = 0;
		active[bb] = 0;
		active[node] = 1;
	}

	// The last-created node (total-1) is the root; its parent stays -1.
	return nodes;
}

} // namespace miint
