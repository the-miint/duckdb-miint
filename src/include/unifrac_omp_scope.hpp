#pragma once

#include <mutex>

namespace miint::unifrac {

// RAII pin of ONE thread's OpenMP fan-out. Takes no lock.
//
// `omp_set_num_threads` writes the calling thread's ICV (OpenMP specifies it as
// the current task's nthreads-var), so each thread that will enter an
// OpenMP-parallel upstream call needs its own pin — a thread spawned by a pin
// holder does not inherit one, and left alone would fan out to every core.
// test_SkbbConcurrency pins this "per-thread, not per-process" property directly,
// because everything below depends on it.
//
// Use this ALONE only for upstream calls that are re-entrant. In this codebase
// that means scikit-bio-binaries' ordination and PERMANOVA — and only with a
// non-negative seed, which is what SkbbCallScope exists to guarantee. libssu's
// UniFrac/Faith PD are NOT re-entrant (see OmpThreadScope); they are reached only
// through wrappers that take the scope themselves.
//
// `n_threads` must be >= 1. Construction saves the calling thread's previous
// value; destruction restores it.
class OmpThreadPin {
public:
	explicit OmpThreadPin(int n_threads);
	~OmpThreadPin();

	OmpThreadPin(const OmpThreadPin &) = delete;
	OmpThreadPin &operator=(const OmpThreadPin &) = delete;
	OmpThreadPin(OmpThreadPin &&) = delete;
	OmpThreadPin &operator=(OmpThreadPin &&) = delete;

private:
	int prev_threads_;
};

// One scikit-bio-binaries compute call (`skbb_pcoa_fsvd_*`, `skbb_permanova_*`),
// runnable CONCURRENTLY with other such calls in the same process — no
// process-wide lock.
//
// Two things make that true, and this class supplies both so a caller cannot get
// one without the other:
//
//   1. the calling thread's OpenMP fan-out is pinned, so concurrent callers do
//      not fight over the thread count;
//   2. `seed()` is ALWAYS >= 0. That is the whole ballgame: with a non-negative
//      seed, skbb seeds a generator per call (ordination/principal_coordinate_
//      analysis.cpp centered_randomize_T, distance/permanova.cpp's
//      RandomGeneratorArray) and touches no shared state. A NEGATIVE seed instead
//      draws from skbb's process-global mt19937 (util/rand.cpp) — an
//      unsynchronized read-modify-write that concurrent callers would race on.
//
// A caller that has no seed (the `seed := -1` default of every miint function
// that reaches skbb) passes it through and gets one derived here, per call, from
// a thread-local generator seeded once from random_device. So "no seed" keeps
// meaning "not reproducible" exactly as before, while never touching skbb's
// generator. An explicit seed >= 0 is passed through untouched, which is what
// keeps seeded results reproducible and comparable across thread counts.
//
// NOT for libssu (UniFrac, Faith PD) — see OmpThreadScope.
class SkbbCallScope {
public:
	SkbbCallScope(int n_threads, int seed);

	SkbbCallScope(const SkbbCallScope &) = delete;
	SkbbCallScope &operator=(const SkbbCallScope &) = delete;
	SkbbCallScope(SkbbCallScope &&) = delete;
	SkbbCallScope &operator=(SkbbCallScope &&) = delete;

	// The seed to hand skbb — never negative. Pass THIS, not the caller's seed.
	int seed() const {
		return seed_;
	}

private:
	OmpThreadPin pin_;
	int seed_;
};

// RAII scope for one libssu compute call: pins OpenMP fan-out AND holds a
// process-wide mutex for its entire lifetime.
//
// libssu and skbb both use raw `#pragma omp parallel for` (see api.cpp:34-36
// upstream: "Threading is now full controlled by OpenMP. Any threads variable
// is really referring to n_substeps."). Without explicit control, OpenMP
// defaults to all detected cores and ignores DuckDB's `SET threads = N`.
// The duckdb-miint contract is the opposite: every compute call honors
// DuckDB's TaskScheduler thread count (or a per-call override).
//
// The mutex is required because libssu is NOT re-entrant, and not merely for RNG
// reasons: `su::process_stripes` brackets every compute with
// register_report_status()/remove_report_status(), which calloc and then
// free+NULL a process-global `report_status` that the compute's inner loop
// dereferences without a NULL check (ext/unifrac-binaries/src/unifrac_internal.cpp),
// and re-init/destroys a static pthread mutex around it. Two concurrent computes
// are therefore a use-after-free, not a numerical inconvenience. It also covers
// `ssu_set_random_seed`, since libssu's entry points take no per-call seed.
// Upstream has been asked to fix the report_status lifecycle
// (localdocs/unifrac-binaries-concurrency-issue.md); until that lands, every
// libssu call goes through this scope.
//
// Prefer SkbbCallScope for skbb ordination/PERMANOVA — those ARE re-entrant with
// a non-negative seed and do not need this lock.
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
	OmpThreadPin pin_; // destroyed before lock_ is released
};

} // namespace miint::unifrac
