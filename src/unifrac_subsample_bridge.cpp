#include "unifrac_subsample_bridge.hpp"

#include <stdexcept>
#include <string>
#include <vector>

#include "unifrac_libssu.hpp"

namespace miint::unifrac {

namespace {

const char *ComputeStatusName(ComputeStatus s) {
	switch (s) {
	case okay:
		return "okay";
	case tree_missing:
		return "tree_missing";
	case table_missing:
		return "table_missing";
	case table_empty:
		return "table_empty";
	case unknown_method:
		return "unknown_method";
	case table_and_tree_do_not_overlap:
		return "table_and_tree_do_not_overlap";
	case output_error:
		return "output_error";
	case invalid_method:
		return "invalid_method";
	case grouping_missing:
		return "grouping_missing";
	}
	return "unknown";
}

// RAII over opaque_biom_inmem_t* so the destroy_subsampled_inmem call
// happens even if downstream materialization throws.
class OpaqueBiomHandle {
public:
	~OpaqueBiomHandle() {
		if (ptr_ != nullptr) {
			destroy_subsampled_inmem(&ptr_);
		}
	}
	opaque_biom_inmem_t **out() {
		return &ptr_;
	}
	opaque_biom_inmem_t *get() const {
		return ptr_;
	}

	OpaqueBiomHandle() = default;
	OpaqueBiomHandle(const OpaqueBiomHandle &) = delete;
	OpaqueBiomHandle &operator=(const OpaqueBiomHandle &) = delete;

private:
	opaque_biom_inmem_t *ptr_ = nullptr;
};

} // namespace

UnifracSupportBiomView BridgeSubsample(const UnifracSupportBiomView &biom_view, uint32_t subsample_depth,
                                       bool with_replacement, int seed) {
	OpaqueBiomHandle handle;
	ComputeStatus status =
	    subsample_table_inmem_seeded(biom_view.support_biom(), subsample_depth, with_replacement, seed, handle.out());
	if (status != okay) {
		throw std::runtime_error(std::string("subsample_table_inmem_seeded returned status ") +
		                         ComputeStatusName(status));
	}

	const auto *opaque = handle.get();
	const unsigned int n_samples = subsampled_n_samples(opaque);
	const unsigned int n_obs = subsampled_n_obs(opaque);
	if (n_samples == 0 || n_obs == 0) {
		throw std::invalid_argument(
		    "BridgeSubsample: subsampled table is empty at subsample_depth=" + std::to_string(subsample_depth) +
		    " (every sample's total count fell below the depth, or every OTU is zero after "
		    "rarefaction); pick a smaller depth");
	}

	// subsampled_get_obs_data fills a dense per-OTU vector indexed by
	// subsampled sample order. Walk obs × sample, emit non-zero cells.
	std::vector<double> obs_buf(n_samples);
	std::vector<CooRow> rows;
	rows.reserve(static_cast<size_t>(n_obs) * static_cast<size_t>(n_samples) / 4); // sparse heuristic
	for (unsigned int o = 0; o < n_obs; ++o) {
		const char *obs_id_cstr = subsampled_get_obs_id(opaque, o);
		if (obs_id_cstr == nullptr) {
			throw std::runtime_error("BridgeSubsample: subsampled_get_obs_id returned null at index " +
			                         std::to_string(o));
		}
		const std::string obs_id = obs_id_cstr;
		// Invariant assertion: `obs_id` came from the same hash-indexed obs
		// map that get_obs_data looks up against. A `false` return here would
		// mean the opaque table is internally inconsistent (libssu bug) — we
		// throw rather than silently mis-emit zero counts.
		if (!subsampled_get_obs_data(opaque, obs_id_cstr, obs_buf.data())) {
			throw std::runtime_error("BridgeSubsample: subsampled_get_obs_data lookup failed for obs '" + obs_id +
			                         "' (opaque table internal inconsistency)");
		}
		for (unsigned int s = 0; s < n_samples; ++s) {
			const double v = obs_buf[s];
			if (v == 0.0) {
				continue;
			}
			const char *sample_id_cstr = subsampled_get_sample_id(opaque, s);
			if (sample_id_cstr == nullptr) {
				throw std::runtime_error("BridgeSubsample: subsampled_get_sample_id returned null at index " +
				                         std::to_string(s));
			}
			rows.push_back({sample_id_cstr, obs_id, v});
		}
	}

	if (rows.empty()) {
		throw std::invalid_argument("BridgeSubsample: subsampled table produced no non-zero cells");
	}

	// UnifracSupportBiomView::FromCoo re-sorts to canonical lex order;
	// after bridging, callers cannot assume the subsampled sample order
	// matches the opaque handle's enumeration. That's intentional —
	// faith_pd_inmem's r_vec output indexes by the new view's order.
	return UnifracSupportBiomView::FromCoo(std::move(rows));
}

} // namespace miint::unifrac
