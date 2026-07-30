#pragma once

#include <mutex>

namespace miint::unifrac {

class OmpThreadScope;

// Proof that some thread of control already holds an OmpThreadScope, and with it
// the process-wide libssu/skbb lock. Only OmpThreadScope can mint one (private
// constructor + friendship), so OmpThreadPin below cannot be constructed by a
// caller that merely believes the lock is held — the unsafe case is unwritable
// rather than merely discouraged.
class OmpScopeHeld {
private:
	friend class OmpThreadScope;
	OmpScopeHeld() = default;
};

// RAII pin of ONE thread's OpenMP fan-out, WITHOUT taking the process-wide lock.
//
// `omp_set_num_threads` writes the calling thread's ICV, so a worker thread
// spawned by an OmpThreadScope holder does not inherit that scope's pin — left
// alone it would fan out to every core, once per worker. Each worker therefore
// needs its own pin, and taking the scope instead would just serialize them.
//
// WHAT MAY RUN UNDER A BARE PIN — skbb's ordination (`skbb_pcoa_fsvd_*`) called
// with a seed >= 0. Its randomization then runs off a generator seeded per call,
// touching no shared state. A seed < 0 draws from skbb's process-global RNG
// (scikit-bio-binaries src/util/rand.cpp) and is NOT thread-safe: callers that
// fan out must substitute non-negative per-call seeds.
//
// WHAT MAY NOT — libssu's UniFrac compute (`UnifracDistanceMatrix::Compute`).
// Its global RNG/state is precisely what the OmpThreadScope mutex serializes, so
// those calls stay one at a time no matter how the surrounding work is arranged.
//
// `n_threads` must be >= 1. Construction saves the calling thread's previous
// value; destruction restores it.
class OmpThreadPin {
public:
	OmpThreadPin(int n_threads, OmpScopeHeld held);
	~OmpThreadPin();

	OmpThreadPin(const OmpThreadPin &) = delete;
	OmpThreadPin &operator=(const OmpThreadPin &) = delete;
	OmpThreadPin(OmpThreadPin &&) = delete;
	OmpThreadPin &operator=(OmpThreadPin &&) = delete;

private:
	int prev_threads_;
};

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
// Thread safety: this class holds a process-wide static mutex for its entire
// lifetime, so callers get the serialization libssu's RNG state requires (this
// scope replaces the previous standalone GlobalRngMutex in unifrac_distance.cpp)
// on top of the thread-count pin.
//
// A holder that fans work out to its own worker threads gives each worker an
// OmpThreadPin(n, held()) — one scope for the group, one pin per thread. See
// OmpThreadPin for what is safe to run that way.
class OmpThreadScope {
public:
	explicit OmpThreadScope(int n_threads);
	~OmpThreadScope();

	OmpThreadScope(const OmpThreadScope &) = delete;
	OmpThreadScope &operator=(const OmpThreadScope &) = delete;
	OmpThreadScope(OmpThreadScope &&) = delete;
	OmpThreadScope &operator=(OmpThreadScope &&) = delete;

	// Hand worker threads the proof they need to pin themselves. Cheap and
	// copyable; it carries no state, only the guarantee that this scope exists.
	// The bearer must not outlive the scope.
	OmpScopeHeld held() const {
		return OmpScopeHeld();
	}

private:
	std::unique_lock<std::mutex> lock_;
	OmpThreadPin pin_; // destroyed before lock_ is released
};

} // namespace miint::unifrac
