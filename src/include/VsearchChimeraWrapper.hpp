#pragma once

#include "uchime_common.hpp"

#include <cassert>
#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Wrapper around vsearch's chimera detection library API.
// Delegates all computation to vsearch's battle-tested implementation.
//
// Thread safety:
// - set_reference() / prepare_denovo(): single-threaded only
// - DetectHandle::detect(): thread-safe (each handle has its own chimera_info_s;
//   global DB/index is read-only after set_reference/prepare_denovo)
// - detect_denovo(): single-threaded (uses internal handle)
// - index_sequence(): single-threaded (modifies global k-mer index)
class VsearchChimeraWrapper {
public:
	// Per-thread detection handle. Wraps a chimera_info_s for thread-safe
	// concurrent detection against the shared read-only reference DB.
	//
	// LIFETIME: Must not outlive the VsearchChimeraWrapper that created it.
	// DuckDB enforces this: LocalState (holds handles) is destroyed before
	// GlobalState (holds wrapper). Using a handle after the wrapper is
	// destroyed is undefined behavior (reads freed vsearch global DB).
	class DetectHandle {
	public:
		~DetectHandle();
		DetectHandle(DetectHandle &&) noexcept;
		DetectHandle &operator=(DetectHandle &&) noexcept;

		DetectHandle(const DetectHandle &) = delete;
		DetectHandle &operator=(const DetectHandle &) = delete;

		UchimeResult detect(const std::string &query_label, const std::string &query_sequence);

	private:
		friend class VsearchChimeraWrapper;
		DetectHandle(); // Only created by VsearchChimeraWrapper::create_detect_handle()
		struct Impl;
		std::unique_ptr<Impl> impl_;
	};

	explicit VsearchChimeraWrapper(const UchimeParams &params = UchimeParams {});
	~VsearchChimeraWrapper();

	VsearchChimeraWrapper(const VsearchChimeraWrapper &) = delete;
	VsearchChimeraWrapper &operator=(const VsearchChimeraWrapper &) = delete;
	VsearchChimeraWrapper(VsearchChimeraWrapper &&) noexcept;
	VsearchChimeraWrapper &operator=(VsearchChimeraWrapper &&) noexcept;

	// Load reference sequences into vsearch's global DB (uchime_ref mode).
	// Acquires vsearch session mutex. Applies DUST masking, builds full k-mer
	// index, and initializes chimera detection infrastructure. After this,
	// create per-thread DetectHandles for concurrent detection.
	void set_reference(const std::vector<std::string> &labels, const std::vector<std::string> &sequences);

	// Load all sequences into vsearch's global DB (uchime_denovo mode).
	// Acquires vsearch session mutex. Applies DUST masking and prepares the
	// k-mer index structure, but does NOT index any sequences. Caller must
	// incrementally index non-chimeras via index_sequence().
	// Sequences must be sorted by abundance descending. When two sequences
	// have equal abundance, the tie-break is caller-determined (SQL ORDER BY).
	void prepare_denovo(const std::vector<std::string> &labels, const std::vector<std::string> &sequences,
	                    const std::vector<int64_t> &abundances);

	// Create a per-thread detection handle for parallel detection.
	// Must be called AFTER set_reference() or prepare_denovo().
	DetectHandle create_detect_handle();

	// Index a single sequence in the k-mer index (de novo mode).
	void index_sequence(uint64_t seqno);

	// Detect chimeras for a batch of queries using vsearch's internal thread pool.
	// Internally parallelizes across opt_threads. Results are appended to output.
	// Must be called AFTER set_reference(). Not thread-safe — call from one thread.
	// Reference mode only (not for de novo).
	void detect_batch(const std::vector<std::string> &query_labels, const std::vector<std::string> &query_sequences,
	                  std::vector<UchimeResult> &output);

	// Detect chimera with abundance skew filtering (de novo mode, single-threaded).
	UchimeResult detect_denovo(const std::string &query_label, const std::string &query_sequence,
	                           int64_t query_abundance);

	const std::vector<std::string> &ref_labels() const {
		return ref_labels_;
	}
	const std::vector<std::string> &ref_sequences() const {
		return ref_sequences_;
	}

private:
	UchimeParams params_;
	std::vector<std::string> ref_labels_;
	std::vector<std::string> ref_sequences_;
	std::vector<int64_t> ref_abundances_;

	// Vsearch state — all managed together for correct teardown ordering.
	// These are only non-null/true between a successful set_reference/prepare_denovo
	// and destruction.
	struct VsearchState;
	std::unique_ptr<VsearchState> state_;

	void teardown();
	void init_common(bool denovo);
	static void load_sequences_into_db(const std::vector<std::string> &labels,
	                                   const std::vector<std::string> &sequences,
	                                   const std::vector<int64_t> &abundances);
	static UchimeResult convert_result(const void *vsearch_result);
};

} // namespace miint
