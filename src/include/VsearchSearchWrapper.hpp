#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <vector>

namespace miint {

// Parameters for global search.
struct SearchParams {
	double id = 0.0;     // minimum identity threshold (0.0-1.0)
	int maxaccepts = 1;  // max accepted hits per query
	int maxrejects = 32; // max rejected targets before stopping
};

// Single search hit result.
struct SearchResult {
	std::string query_id;
	std::string target_id;
	double identity = 0.0; // percent identity (0-100)
	int matches = 0;
	int mismatches = 0;
	int gaps = 0;
	int alignment_length = 0;
	int query_length = 0;
	int target_length = 0;
	bool accepted = false;
};

// Wrapper around vsearch's global search library API.
// Delegates all computation to vsearch's battle-tested implementation.
//
// Thread safety:
// - set_database(): single-threaded only
// - SearchHandle::search(): not thread-safe across handles (vsearch permits
//   only one active search session per process; do not share a handle across
//   threads). Use search_batch() for parallel queries.
class VsearchSearchWrapper {
public:
	// Search handle. Wraps a search_session_s for per-query search against
	// the shared reference DB. Only one handle may be active per process
	// (vsearch enforces one active session); for parallel queries, use
	// search_batch() which internally parallelizes across opt_threads.
	//
	// LIFETIME: Must not outlive the VsearchSearchWrapper that created it.
	class SearchHandle {
	public:
		~SearchHandle();
		SearchHandle(SearchHandle &&) noexcept;
		SearchHandle &operator=(SearchHandle &&) noexcept;

		SearchHandle(const SearchHandle &) = delete;
		SearchHandle &operator=(const SearchHandle &) = delete;

		// Search a single query against the database.
		// Returns 0..maxaccepts results ordered by identity (descending).
		std::vector<SearchResult> search(const std::string &query_label, const std::string &query_sequence,
		                                 int query_size = 1);

	private:
		friend class VsearchSearchWrapper;
		SearchHandle(int maxaccepts);
		struct Impl;
		std::unique_ptr<Impl> impl_;
	};

	explicit VsearchSearchWrapper(const SearchParams &params = SearchParams {});
	~VsearchSearchWrapper();

	VsearchSearchWrapper(const VsearchSearchWrapper &) = delete;
	VsearchSearchWrapper &operator=(const VsearchSearchWrapper &) = delete;
	VsearchSearchWrapper(VsearchSearchWrapper &&) noexcept;
	VsearchSearchWrapper &operator=(VsearchSearchWrapper &&) noexcept;

	// Load reference sequences into vsearch's global DB and build index.
	// Acquires vsearch session mutex. Applies DUST masking and builds full
	// k-mer index. After this, create a SearchHandle for per-query search.
	void set_database(const std::vector<std::string> &labels, const std::vector<std::string> &sequences);

	// Create a search handle. Must be called AFTER set_database().
	// Only one handle may be active per process at a time.
	SearchHandle create_search_handle();

	// Search a batch of queries using vsearch's internal thread pool.
	// Internally parallelizes across opt_threads. Results are appended to output.
	// Must be called AFTER set_database(). Not thread-safe — call from one thread.
	void search_batch(const std::vector<std::string> &query_labels, const std::vector<std::string> &query_sequences,
	                  std::vector<SearchResult> &output);

private:
	SearchParams params_;

	struct VsearchState;
	std::unique_ptr<VsearchState> state_;
};

} // namespace miint
