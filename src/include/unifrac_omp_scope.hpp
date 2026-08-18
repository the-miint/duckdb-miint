#pragma once

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
// Use this ALONE for upstream calls that draw no randomness at all — Faith PD's
// `faith_pd_inmem` is the one such caller here. Anything that can draw needs a
// non-negative seed as well, which is what ComputeCallScope exists to guarantee.
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

// One compute call into an OpenMP-parallel native library — libssu's
// `one_off_matrix_inmem_*_v4` or scikit-bio-binaries' `skbb_pcoa_fsvd_*` /
// `skbb_permanova_*` — runnable CONCURRENTLY with other such calls in the same
// process, with no process-wide lock.
//
// Two things make that true, and this class supplies both so a caller cannot get
// one without the other:
//
//   1. the calling thread's OpenMP fan-out is pinned, so concurrent callers do
//      not fight over the thread count;
//   2. `seed()` is ALWAYS >= 0. That is the whole ballgame. With a non-negative
//      seed each library seeds a generator per call and touches no shared state
//      (skbb: ordination/principal_coordinate_analysis.cpp centered_randomize_T
//      and distance/permanova.cpp's RandomGeneratorArray; libssu: a per-call
//      `su::biom_subsampled`). A NEGATIVE seed instead routes both to a
//      process-global mt19937 — an unsynchronized read-modify-write that
//      concurrent callers race on. Upstream measured what that costs: 4 threads
//      seeding 42 made 52 of 161 results disagree with the serial answer.
//
// A caller that has no seed (the `seed := -1` default of every miint function
// reaching these libraries) passes it through and gets one derived here, per
// call, from a thread-local generator seeded once from random_device. So "no
// seed" keeps meaning "not reproducible" exactly as before, while never touching
// either global generator. An explicit seed >= 0 is passed through untouched,
// which is what keeps seeded results reproducible.
//
// Never reseed a global instead: libssu's `ssu_set_random_seed` chains TWO
// generators (it reseeds libssu's mt19937, then draws once from it to seed
// skbb's), so one such call perturbs every later unseeded draw in the process.
// `test_UnifracDistance.cpp` pins that Compute does not do this.
class ComputeCallScope {
public:
	ComputeCallScope(int n_threads, int seed);

	ComputeCallScope(const ComputeCallScope &) = delete;
	ComputeCallScope &operator=(const ComputeCallScope &) = delete;
	ComputeCallScope(ComputeCallScope &&) = delete;
	ComputeCallScope &operator=(ComputeCallScope &&) = delete;

	// The seed to hand the library — never negative. Pass THIS, not the caller's.
	int seed() const {
		return seed_;
	}

private:
	OmpThreadPin pin_;
	int seed_;
};

} // namespace miint::unifrac
