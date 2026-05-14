#pragma once

#include <cstdint>
#include <string>
#include <vector>

#include "unifrac_libssu.hpp" // guarded wrapper over unifrac api.hpp (support_biom_t)

namespace miint::unifrac {

struct CooRow {
	std::string sample_id;
	std::string feature_id;
	double count;
};

// Builds a libssu-compatible CSR (obs-major) view over heap-owned storage.
//
// libssu's support_biom_t is documented at api.hpp:223 ("CSR-encoded input
// table") and consumed by biom_inmem.cpp:144-156 which iterates over
// [0, n_obs) and reads sample column indices from `indices`. Layout:
//   indptr   length n_obs+1; indptr[i]..indptr[i+1] indexes obs i's slice
//   indices  length nnz;     sample column index for each (obs, sample) cell
//   data     length nnz;     count value for each (obs, sample) cell
//
// Sums duplicate (sample, feature) rows; drops cells that aggregate to zero;
// canonicalizes sample_ids and feature_ids by lexicographic order so that
// downstream seeded computations produce identical results across input
// shuffles.
//
// Lifetime: support_biom() returns a pointer whose char**/uint32_t*/double*
// fields are valid for the lifetime of the view. Move-only.
class UnifracSupportBiomView {
public:
	// Throws std::invalid_argument on empty input.
	static UnifracSupportBiomView FromCoo(std::vector<CooRow> rows);

	const support_biom_t *support_biom() const {
		return &biom_;
	}

	UnifracSupportBiomView(UnifracSupportBiomView &&) noexcept;
	UnifracSupportBiomView &operator=(UnifracSupportBiomView &&) noexcept;
	UnifracSupportBiomView(const UnifracSupportBiomView &) = delete;
	UnifracSupportBiomView &operator=(const UnifracSupportBiomView &) = delete;
	~UnifracSupportBiomView() = default;

private:
	UnifracSupportBiomView() = default;
	void wire_biom_ptrs();

	std::vector<std::string> sample_ids_;
	std::vector<std::string> feature_ids_;
	std::vector<uint32_t> indices_;
	std::vector<uint32_t> indptr_;
	std::vector<double> data_;

	std::vector<char *> sample_id_ptrs_;
	std::vector<char *> feature_id_ptrs_;
	support_biom_t biom_ {};
};

} // namespace miint::unifrac
