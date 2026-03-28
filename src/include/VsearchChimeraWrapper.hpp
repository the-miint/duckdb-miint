#pragma once

#include "uchime_common.hpp"

#include <cstdint>
#include <memory>
#include <mutex>
#include <string>
#include <vector>

namespace miint {

// Wrapper around vsearch's chimera detection library API.
// Provides the same interface as the previous ChimeraDetector class
// but delegates all computation to vsearch's battle-tested implementation.
//
// Thread safety:
// - set_reference() must be called from a single thread (acquires global mutex)
// - detect() is safe for concurrent calls with different wrapper instances
//   (each instance has its own chimera_info_s; global DB/index is read-only)
// - detect_denovo() is single-threaded (modifies global DB incrementally)
class VsearchChimeraWrapper {
public:
	explicit VsearchChimeraWrapper(const UchimeParams &params = UchimeParams {});
	~VsearchChimeraWrapper();

	// Non-copyable (wraps global state)
	VsearchChimeraWrapper(const VsearchChimeraWrapper &) = delete;
	VsearchChimeraWrapper &operator=(const VsearchChimeraWrapper &) = delete;
	VsearchChimeraWrapper(VsearchChimeraWrapper &&) noexcept;
	VsearchChimeraWrapper &operator=(VsearchChimeraWrapper &&) noexcept;

	// Load reference sequences into vsearch's global DB.
	// Applies DUST masking and builds k-mer index.
	void set_reference(const std::vector<std::string> &labels, const std::vector<std::string> &sequences);

	// Add a single sequence to the reference (for de novo incremental mode).
	void add_to_reference(const std::string &label, const std::string &sequence, int64_t abundance = 0);

	// Index a sequence that was previously added via add_to_reference.
	// Call after detect_denovo classifies a sequence as non-chimeric.
	void index_sequence(uint64_t seqno);

	// Detect chimera for a single query.
	UchimeResult detect(const std::string &query_label, const std::string &query_sequence);

	// Detect chimera with abundance skew filtering (de novo mode).
	UchimeResult detect_denovo(const std::string &query_label, const std::string &query_sequence,
	                           int64_t query_abundance);

	const std::vector<std::string> &ref_labels() const {
		return ref_labels_;
	}
	const std::vector<std::string> &ref_sequences() const {
		return ref_sequences_;
	}

private:
	struct Impl; // Pimpl to hide vsearch types from the header
	std::unique_ptr<Impl> impl_;

	UchimeParams params_;
	std::vector<std::string> ref_labels_;
	std::vector<std::string> ref_sequences_;
	std::vector<int64_t> ref_abundances_;
	bool initialized_ = false;

	void init_vsearch_globals();
	static UchimeResult convert_result(const void *vsearch_result); // chimera_result_s*
};

} // namespace miint
