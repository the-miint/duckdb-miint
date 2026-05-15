#pragma once

#include <mutex>

namespace miint::unifrac {

// RAII scope that pins libssu / scikit-bio-binaries OpenMP fan-out to a
// caller-chosen thread count for the duration of one upstream compute call.
//
// libssu and skbb both use raw `#pragma omp parallel for` (see api.cpp:34-36
// upstream: "Threading is now full controlled by OpenMP. Any threads variable
// is really referring to n_substeps."). Without explicit control, OpenMP
// defaults to all detected cores and ignores DuckDB's `SET threads = N`.
// The duckdb-miint contract is the opposite: every compute call honors
// DuckDB's TaskScheduler thread count (or a per-call override).
//
// Thread safety: `omp_set_num_threads` mutates a process-global setting, so
// concurrent libssu/skbb calls from different DuckDB queries would race.
// This class holds a process-wide static mutex for its entire lifetime —
// callers therefore also get the serialization that libssu's RNG state
// requires (this scope replaces the previous standalone GlobalRngMutex in
// unifrac_distance.cpp).
//
// Construction sets the OpenMP thread count to `n_threads` after saving the
// previous value; destruction restores it. `n_threads` must be >= 1.
class OmpThreadScope {
public:
	explicit OmpThreadScope(int n_threads);
	~OmpThreadScope();

	OmpThreadScope(const OmpThreadScope &) = delete;
	OmpThreadScope &operator=(const OmpThreadScope &) = delete;
	OmpThreadScope(OmpThreadScope &&) = delete;
	OmpThreadScope &operator=(OmpThreadScope &&) = delete;

private:
	std::unique_lock<std::mutex> lock_;
	int prev_threads_;
};

} // namespace miint::unifrac
