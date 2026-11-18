#pragma once

/*
 * KreppIndexBuilder - C++ wrapper for krepp index construction.
 *
 * The mirror of KreppPlacer: that one loads an index krepp's CLI wrote, this
 * one writes an index KreppPlacer can load, from the same linked krepp. Same
 * PIMPL discipline for the same reason - krepp's headers reach exactly one
 * translation unit, because their klib basenames collide with minimap2's and
 * their seq_nt4_table is a different table of a different length.
 *
 * INPUT IS FILES, NOT MEMORY, AND THAT IS KREPP'S SHAPE
 * ====================================================
 * IndexMultiple::read_input_file takes a path and nothing else. Handed a FASTA
 * it switches to "per sequence" mode, where every record is its own reference -
 * and obtain_build_tree then REFUSES a guide tree in that mode outright
 * (index.cpp: "A guide tree (-t) is incompatible with a per sequence
 * indexing."). So an index with a backbone can only be built the other way:
 * a TSV mapping each reference name to a file of its own.
 *
 * That is why the caller materialises one FASTA per reference. It is not a
 * convenience; it is the only shape that keeps the backbone.
 *
 * ERRORS
 * ======
 * krepp reports fatal errors with error_exit, which is std::exit - a malformed
 * tree would take the DuckDB process down with no SQL error at all.
 * InstallKreppErrorHandler redirects those into KreppFatalError.
 *
 * A throwing handler is only sound because krepp catches at the boundary of
 * every OpenMP region in its index build path and rethrows on the calling
 * thread (bo1929/krepp#14). Before that an exception could not leave a
 * structured block at all, so this build had to compile krepp with _WOPENMP=0
 * to stay correct - which is why there was no thread count to offer.
 *
 * That guarantee is scoped to the BUILD. krepp's index load path still has
 * unguarded critical regions, and KreppPlacer calls into it. The throws
 * reachable there are allocation failures rather than bad input, but an
 * error_exit raised inside one would still terminate the process.
 */

#include <cstdint>
#include <optional>
#include <stdexcept>
#include <string>

namespace miint {

// What krepp's error_exit becomes once InstallKreppErrorHandler has run.
// `code` is the exit status krepp would have used; it is EXIT_FAILURE for
// every call site in the library sources, and carried anyway so a caller can
// tell them apart if that ever stops being true.
class KreppFatalError : public std::runtime_error {
public:
	KreppFatalError(const std::string &message, int code) : std::runtime_error(message), code_(code) {
	}
	int code() const noexcept {
		return code_;
	}

private:
	int code_;
};

// Install the handler that turns krepp's fatal errors into KreppFatalError.
// Idempotent, and safe to call from several threads.
//
// Called once at extension load rather than lazily at first use, which is what
// krepp's contract asks for: error_exit reads the handler without a lock, so
// installing it while another thread might already be inside krepp is a data
// race. Loading the extension is the one moment when no krepp call can be in
// flight.
void InstallKreppErrorHandler();

// Everything krepp's IndexConfig needs, minus the CLI. Defaults match krepp's
// own; `w` and `h` are optional so "not given" (derive from k) stays distinct
// from "given as 0" (invalid, and krepp says so).
struct KreppIndexOptions {
	std::string index_dir;      // directory to create and write the index into
	std::string input_map_path; // TSV: reference name <TAB> path to its FASTA
	std::string newick_path;    // backbone tree; empty means krepp generates one
	uint8_t k = 29;
	std::optional<uint8_t> w; // unset: k + 6
	std::optional<uint8_t> h; // unset: k - 16
	uint32_t m = 4;
	uint32_t r = 1;
	bool frac = true;
	uint32_t sdust_t = 0;
	uint32_t sdust_w = 0;
	// 1 is what every release before this did. Above 1 only does anything when
	// krepp was compiled with its OpenMP regions in - ask
	// KreppIndexThreadsSupported() rather than assuming.
	uint32_t threads = 1;
};

// Whether a `threads` above 1 can do anything: true when krepp was built with
// an OpenMP runtime available. A build without one takes threads = 1 and
// refuses more, rather than accepting a thread count and quietly ignoring it.
bool KreppIndexThreadsSupported();

// krepp's own VERSION string (ext/krepp/src/common.hpp), which is what krepp
// writes as `krepp version:` into each partial's metadata<suffix>.txt. Not the
// same as KREPP_GIT_VERSION, which is `git describe` of the submodule and
// agrees with it only on a tagged checkout. Exposed through this TU because
// krepp's headers are on the include path of the krepp sources alone.
const char *KreppVersionString();

// Build the index. Runs krepp's own index pipeline in the order krepp's main
// runs it, and writes the same files to `index_dir`.
//
// Throws KreppFatalError for anything krepp itself rejects - an unreadable
// input, a malformed tree, an out-of-range parameter. Runs on
// options.threads threads.
void BuildKreppIndex(const KreppIndexOptions &options);

} // namespace miint
