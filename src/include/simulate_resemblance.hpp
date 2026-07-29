#pragma once

#include "duckdb.hpp"

#include <cstdint>
#include <string>
#include <vector>

namespace duckdb {
class ExtensionLoader;
}

namespace miint {

//! Long/COO result of a resemblance simulation: parallel columns, one entry per
//! nonzero (sample, OTU) cell. `ground_truth` carries the per-sample gradient
//! position (gradient model) or cluster id (cluster model).
struct SimulationCOO {
	std::vector<int32_t> sample_id;
	std::vector<int32_t> otu_id;
	std::vector<int64_t> count;
	std::vector<double> ground_truth;

	size_t size() const {
		return sample_id.size();
	}
};

//! Evenly-spaced points on [lo, hi] with numpy.linspace semantics: n==1 => {lo};
//! n>=2 => lo + (hi-lo)*i/(n-1) for i in [0,n). Deterministic (RNG-free) so the
//! gradient geometry can be pinned by unit tests.
std::vector<double> Linspace(double lo, double hi, int32_t n);

//! Deterministic species-response matrix (rows = samples, cols = species): each
//! entry = abundances[s] * exp(-(positions[i]-optima[s])^2 / (2*sp_width^2))
//! / (sp_width*sqrt(2*pi)). Split out from SimulateGradient (which draws `optima`
//! at random) so the exponent/normalization math is unit-testable without RNG.
//! Requires optima.size() == abundances.size().
std::vector<std::vector<double>> BuildTrueDistribution(const std::vector<double> &abundances,
                                                       const std::vector<double> &positions,
                                                       const std::vector<double> &optima, double sp_width);

//! Non-phylogenetic gradient simulator (Kuczynski et al. 2010). Mirrors the model
//! in ord_survey `generate_1d_gradient_data`: per-species Gaussian response curves
//! with random optima along [0,1], sampled at `num_samples` evenly-spaced positions
//! in [range_lo, range_hi], optional noise, then a multinomial draw of
//! `seqs_per_sample` reads per sample. Globally-absent OTUs never appear (COO is
//! nonzero-only). `ground_truth` = the sample's gradient position.
//!
//! noise_type selects where each perturbation's width comes from:
//!   "+species" : x += N(0, noise*x)       -- that species' own abundance
//!   "*sample"  : x *= N(1, noise*rowsum)  -- the sample total, shared by all species
//!   "+sample"  : x += N(0, noise*rowsum)  -- the sample total, shared by all species
//! Negatives are then floored to 0 by subtracting the row minimum, and the row is
//! rescaled to its original mean.
//!
//! seed < 0 => nondeterministic (random_device); seed >= 0 => reproducible.
SimulationCOO SimulateGradient(const std::vector<double> &abundances, int32_t num_samples, int64_t seqs_per_sample,
                               double sp_width, double noise, const std::string &noise_type, double range_lo,
                               double range_hi, int64_t seed);

//! Non-phylogenetic cluster simulator (Kuczynski et al. 2010). Mirrors the model
//! in ord_survey `ClusterSimulation`: from a base abundance vector, build one
//! centroid per cluster by perturbing the base (magnitude `cluster_spacing`), then
//! draw each within-cluster sample by perturbing its centroid (magnitude
//! `sample_spacing`); each perturbation is `_perturb_node` (see the .cpp), then a
//! multinomial draw of `seqs_per_sample` reads. `cluster_sizes[k]` samples are
//! emitted for cluster k. Globally-absent OTUs never appear (COO is nonzero-only).
//! `ground_truth` = the sample's integer cluster id (as a double).
//!
//! noise_type in {"*sample","+sample"}, normalization in {"clip","rescale"}.
//! seed < 0 => nondeterministic; seed >= 0 => reproducible.
SimulationCOO SimulateCluster(const std::vector<double> &abundances, const std::vector<int32_t> &cluster_sizes,
                              int64_t seqs_per_sample, double cluster_spacing, double sample_spacing,
                              const std::string &noise_type, const std::string &normalization, int64_t seed);

} // namespace miint

namespace duckdb {
//! Registers the resemblance-simulator table functions (simulate_gradient_otus,
//! and later simulate_cluster_otus) into the extension catalog.
void RegisterSimulateResemblance(ExtensionLoader &loader);
} // namespace duckdb
