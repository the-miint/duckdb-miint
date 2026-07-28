#include "simulate_resemblance.hpp"

#include <algorithm>
#include <cmath>
#include <random>
#include <stdexcept>
#include <string>

namespace miint {

namespace {
constexpr double TWO_PI = 6.28318530717958647692;

double VecMean(const std::vector<double> &v) {
	double s = 0.0;
	for (double x : v) {
		s += x;
	}
	return s / static_cast<double>(v.size());
}

double VecSum(const std::vector<double> &v) {
	double s = 0.0;
	for (double x : v) {
		s += x;
	}
	return s;
}

//! Draw `n` reads from a categorical distribution over `weights` (a multinomial
//! sample), writing per-category counts into `counts` (resized/zeroed here).
//! Shared by the gradient and cluster simulators. `weights` must be non-negative
//! with a positive sum (callers guard the empty-sample case with a clear error).
void MultinomialDraw(std::mt19937_64 &rng, const std::vector<double> &weights, int64_t n,
                     std::vector<int64_t> &counts) {
	counts.assign(weights.size(), 0);
	std::discrete_distribution<size_t> draw(weights.begin(), weights.end());
	for (int64_t k = 0; k < n; k++) {
		counts[draw(rng)]++;
	}
}

//! One perturbation step of ord_survey `_perturb_node` (assumes node.sum()==1,
//! preserves it). Multiplicative: x *= N(1, mag); additive: x += N(0, mag). Then
//! 'clip' floors negatives at 0, 'rescale' shifts the floor to 0; finally renorm
//! so the result sums to 1. Throws if the result would be all-zero (degenerate).
std::vector<double> PerturbNode(std::mt19937_64 &rng, const std::vector<double> &node, double mag,
                                const std::string &noise_type, const std::string &normalization) {
	std::vector<double> out(node.size());
	const bool multiplicative = (noise_type == "*sample");
	std::normal_distribution<double> perturb(multiplicative ? 1.0 : 0.0, mag);
	for (size_t s = 0; s < node.size(); s++) {
		const double z = perturb(rng);
		out[s] = multiplicative ? node[s] * z : node[s] + z;
	}
	if (normalization == "clip") {
		for (double &x : out) {
			if (x < 0.0) {
				x = 0.0;
			}
		}
	} else { // "rescale": shift the minimum to 0
		const double floor_val = *std::min_element(out.begin(), out.end());
		for (double &x : out) {
			x -= floor_val;
		}
	}
	const double sum = VecSum(out);
	// !(sum > 0) rather than (sum <= 0) so a NaN sum (e.g. from a non-finite node)
	// is also rejected -- NaN comparisons are all false, so `sum <= 0` would pass it.
	if (!(sum > 0.0)) {
		throw std::invalid_argument("simulate_cluster_otus: degenerate perturbation produced an all-zero node; "
		                            "reduce cluster_spacing/sample_spacing");
	}
	for (double &x : out) {
		x /= sum;
	}
	return out;
}
} // namespace

std::vector<double> Linspace(double lo, double hi, int32_t n) {
	std::vector<double> pos(static_cast<size_t>(n));
	if (n == 1) {
		pos[0] = lo;
	} else {
		for (int32_t i = 0; i < n; i++) {
			pos[static_cast<size_t>(i)] = lo + (hi - lo) * i / (n - 1);
		}
	}
	return pos;
}

std::vector<std::vector<double>> BuildTrueDistribution(const std::vector<double> &abundances,
                                                       const std::vector<double> &positions,
                                                       const std::vector<double> &optima, double sp_width) {
	const size_t num_sp = abundances.size();
	const double norm = 1.0 / (sp_width * std::sqrt(TWO_PI));
	const double two_var = 2.0 * sp_width * sp_width;
	std::vector<std::vector<double>> tru(positions.size(), std::vector<double>(num_sp));
	for (size_t i = 0; i < positions.size(); i++) {
		auto &row = tru[i];
		for (size_t s = 0; s < num_sp; s++) {
			const double d = positions[i] - optima[s];
			row[s] = abundances[s] * std::exp(-(d * d) / two_var) * norm;
		}
	}
	return tru;
}

SimulationCOO SimulateGradient(const std::vector<double> &abundances, int32_t num_samples, int64_t seqs_per_sample,
                               double sp_width, double noise, const std::string &noise_type, double range_lo,
                               double range_hi, int64_t seed) {
	const size_t num_sp = abundances.size();
	const int32_t n = num_samples;

	std::mt19937_64 rng(seed >= 0 ? static_cast<uint64_t>(seed) : std::random_device {}());

	// Evenly-spaced sample positions along the gradient.
	const std::vector<double> pos = Linspace(range_lo, range_hi, n);

	// Per-species Gaussian response curves with random optima in [0,1).
	std::uniform_real_distribution<double> unit(0.0, 1.0);
	std::vector<double> optima(num_sp);
	for (size_t s = 0; s < num_sp; s++) {
		optima[s] = unit(rng);
	}
	std::vector<std::vector<double>> tru = BuildTrueDistribution(abundances, pos, optima, sp_width);

	// Optional perturbation. The noise width is drawn from either the species' own
	// abundance or the sample total:
	//   '+species' : x += N(0, noise*x)       -- per-species width (default)
	//   '*sample'  : x *= N(1, noise*rowsum)  -- shared sample-total width
	//   '+sample'  : x += N(0, noise*rowsum)  -- shared sample-total width
	// Then shift floor to 0 and rescale to preserve the sample's original mean
	// abundance (matches ord_survey).
	if (noise != 0.0) {
		const bool per_species = (noise_type == "+species");
		const bool multiplicative = (noise_type == "*sample");
		for (int32_t i = 0; i < n; i++) {
			auto &row = tru[static_cast<size_t>(i)];
			const double original_mean = VecMean(row);
			if (per_species) {
				// Width scales with each element, so draw a standard normal and
				// scale per species. One distribution object keeps the engine
				// consumption uniform across species.
				std::normal_distribution<double> unit_normal(0.0, 1.0);
				for (size_t s = 0; s < num_sp; s++) {
					row[s] += unit_normal(rng) * noise * row[s];
				}
			} else {
				const double scale = noise * VecSum(row); // rowsum BEFORE perturbation
				std::normal_distribution<double> perturb(multiplicative ? 1.0 : 0.0, scale);
				for (size_t s = 0; s < num_sp; s++) {
					const double z = perturb(rng);
					if (multiplicative) {
						row[s] *= z;
					} else {
						row[s] += z;
					}
				}
			}
			const double floor_val = *std::min_element(row.begin(), row.end());
			for (size_t s = 0; s < num_sp; s++) {
				row[s] -= floor_val;
			}
			const double new_mean = VecMean(row);
			if (new_mean != 0.0) {
				const double factor = original_mean / new_mean;
				for (size_t s = 0; s < num_sp; s++) {
					row[s] *= factor;
				}
			}
		}
	}

	// Multinomial draw of seqs_per_sample reads from each sample's abundance
	// profile; emit only nonzero cells (globally-absent OTUs never appear).
	SimulationCOO coo;
	std::vector<int64_t> counts;
	for (int32_t i = 0; i < n; i++) {
		auto &row = tru[static_cast<size_t>(i)];
		if (VecSum(row) <= 0.0) {
			throw std::invalid_argument("simulate_gradient_otus: empty sample at index " + std::to_string(i) +
			                            " (all-zero abundance); check abundances / noise settings");
		}
		MultinomialDraw(rng, row, seqs_per_sample, counts);
		for (size_t s = 0; s < num_sp; s++) {
			if (counts[s] > 0) {
				coo.sample_id.push_back(i);
				coo.otu_id.push_back(static_cast<int32_t>(s));
				coo.count.push_back(counts[s]);
				coo.ground_truth.push_back(pos[static_cast<size_t>(i)]);
			}
		}
	}
	return coo;
}

SimulationCOO SimulateCluster(const std::vector<double> &abundances, const std::vector<int32_t> &cluster_sizes,
                              int64_t seqs_per_sample, double cluster_spacing, double sample_spacing,
                              const std::string &noise_type, const std::string &normalization, int64_t seed) {
	const size_t num_sp = abundances.size();

	std::mt19937_64 rng(seed >= 0 ? static_cast<uint64_t>(seed) : std::random_device {}());

	// Base ("root") node: relative abundances summing to 1 (perturbation is
	// scale-invariant, but this matches _perturb_node's stated precondition).
	// Guard before dividing: an all-zero (or non-finite) sum would make root NaN,
	// which the SQL wrapper's ParseAbundancesArg blocks -- but this pure core is
	// public/independently-tested, so it must fail loud on its own.
	std::vector<double> root(abundances);
	const double root_sum = VecSum(root);
	if (!(root_sum > 0.0)) {
		throw std::invalid_argument("simulate_cluster_otus: abundances must sum to a positive, finite value");
	}
	for (double &x : root) {
		x /= root_sum;
	}

	SimulationCOO coo;
	std::vector<int64_t> counts;
	int32_t sample_id = 0;
	for (size_t k = 0; k < cluster_sizes.size(); k++) {
		// One centroid per cluster (perturb the root by cluster_spacing), then each
		// within-cluster sample perturbs that centroid by sample_spacing.
		const std::vector<double> centroid = PerturbNode(rng, root, cluster_spacing, noise_type, normalization);
		for (int32_t j = 0; j < cluster_sizes[k]; j++) {
			const std::vector<double> sample_node =
			    PerturbNode(rng, centroid, sample_spacing, noise_type, normalization);
			MultinomialDraw(rng, sample_node, seqs_per_sample, counts);
			for (size_t s = 0; s < num_sp; s++) {
				if (counts[s] > 0) {
					coo.sample_id.push_back(sample_id);
					coo.otu_id.push_back(static_cast<int32_t>(s));
					coo.count.push_back(counts[s]);
					coo.ground_truth.push_back(static_cast<double>(k));
				}
			}
			sample_id++;
		}
	}
	return coo;
}

} // namespace miint
