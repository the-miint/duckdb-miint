#include "KreppIndexBuilder.hpp"

#include <mutex>

// krepp's own headers. krepp.hpp is deliberately not included: it pulls in
// CLI11 and the command-line layer, neither of which applies here. Everything
// the builder needs moved into index.hpp/cpp upstream (bo1929/krepp#13) for
// exactly this reason.
#include "common.hpp"
#include "index.hpp"

namespace {

#if defined(MIINT_KREPP_OPENMP) && MIINT_KREPP_OPENMP == 1
// Declared rather than pulled in from <omp.h>: this TU is compiled with
// _WOPENMP=1 but WITHOUT the OpenMP flag (see cmake/deps/krepp.cmake - the flag is
// scoped to the krepp archive), so _OPENMP is undefined here and omp.h is not
// on the include path. These two symbols come from the libomp the krepp target
// links, and their signatures are fixed by the OpenMP ABI.
extern "C" int omp_get_max_threads(void);
extern "C" void omp_set_num_threads(int);

// krepp's build_index() calls omp_set_num_threads(num_threads). That writes the
// *encountering thread's* nthreads-var, which is thread state, not krepp state:
// it outlives the build and the lock. Measured with a probe against the real
// libomp - a thread's max_threads went 10 -> 8 and stayed 8 after the call
// returned. UniFrac ships in the same binary (both MIINT_ENABLE_UNIFRAC and
// MIINT_ENABLE_KREPP default ON) and 29 of its 30 parallel regions carry no
// num_threads clause, so they take that per-thread default. A later UniFrac
// operator landing on the DuckDB worker that just built an index would inherit
// whatever the index build asked for. Put it back.
struct RestoreOmpThreadCount {
	int saved;
	RestoreOmpThreadCount() : saved(omp_get_max_threads()) {
	}
	~RestoreOmpThreadCount() {
		omp_set_num_threads(saved);
	}
	RestoreOmpThreadCount(const RestoreOmpThreadCount &) = delete;
	RestoreOmpThreadCount &operator=(const RestoreOmpThreadCount &) = delete;
};
#endif

} // namespace

namespace miint {

void InstallKreppErrorHandler() {
	static std::once_flag installed;
	std::call_once(installed, [] {
		set_error_handler([](const std::string &message, int code) { throw KreppFatalError(message, code); });
	});
}

bool KreppIndexThreadsSupported() {
#if defined(MIINT_KREPP_OPENMP) && MIINT_KREPP_OPENMP == 1
	return true;
#else
	return false;
#endif
}

const char *KreppVersionString() {
	return VERSION;
}

void BuildKreppIndex(const KreppIndexOptions &options) {
	// krepp keeps its thread count in a plain global - `uint32_t num_threads` in
	// common.cpp - written just below and read much later, inside build_index().
	// Two krepp_index_create queries running at once would race on that one
	// word, and a build would then run at the other query's thread count. That
	// is not purely a performance matter: krepp merges each node's children in
	// task-completion order, so the colour arrays (cmer's subset column, and
	// crecord) do depend on the thread count - measured differing between a
	// serial and a threaded build of the same corpus. The k-mers indexed and
	// their buckets do not (krepp_index_create_threads.test compares inc-m4r1-frac
	// byte for byte across 1 and 4 threads). It is a data race either way.
	//
	// The whole build is serialised rather than just the write, because the read
	// is deep inside build_index(): a lock around set_num_threads alone would
	// not make the write-then-read pair atomic. Queueing costs little, since
	// each build is already parallel across its own OpenMP team.
	static std::mutex build_mutex;
	const std::lock_guard<std::mutex> build_guard(build_mutex);

#if defined(MIINT_KREPP_OPENMP) && MIINT_KREPP_OPENMP == 1
	// Destructor runs on the throwing path too, which matters: BuildKreppIndex
	// exits by exception for anything krepp rejects.
	const RestoreOmpThreadCount omp_guard;
#endif

	// Before the build and on this thread: krepp reads the num_threads global
	// inside index_files()/index_sequences() to size its team. set_num_threads
	// clamps 0, which omp_set_num_threads does not accept.
	set_num_threads(options.threads);

	// set_lshf() draws the LSH positions from `gen`, a thread_local mt19937 that
	// krepp's CLI seeds only under --seed. A fresh krepp process therefore starts
	// from mt19937's default state, so the same input gives the same index. A
	// long-lived one does not: the second build resumes the stream where the
	// first left it and picks a different hash function for identical
	// references. Seeding here restores the CLI's property. The first build in a
	// process is unaffected - it already started from this state.
	gen.seed(std::mt19937::default_seed);

	IndexConfig config;
	config.input = options.input_map_path;
	config.index_dir = options.index_dir;
	config.nwk_path = options.newick_path;
	config.k = options.k;
	config.w = options.w;
	config.h = options.h;
	config.m = options.m;
	config.r = options.r;
	config.frac = options.frac;
	config.sdust_t = options.sdust_t;
	config.sdust_w = options.sdust_w;

	// The order krepp's main uses for its `index` subcommand, unchanged. It is
	// not interchangeable: set_nrows/set_lshf seed the hash function the whole
	// build hangs off, read_input_file decides per-sequence vs per-file mode,
	// and obtain_build_tree consults that decision.
	IndexMultiple index(config);
	index.set_nrows();
	index.set_lshf();
	index.read_input_file();
	index.obtain_build_tree();
	index.build_index();
	index.save_index();
}

} // namespace miint
