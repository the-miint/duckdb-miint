#pragma once

/*
 * KreppPlacer - C++ wrapper for krepp phylogenetic placement.
 *
 * krepp answers "where on this backbone tree does each read belong?" by
 * matching k-mers against a prebuilt index and maximising a pseudo-likelihood
 * over candidate edges.
 *
 * OWNERSHIP MODEL:
 * ===============
 * SharedKreppIndex holds the loaded index and the backbone tree. Both are
 * immutable once constructed and the query path only reads them, so a single
 * instance is shared by every thread rather than loaded per thread - indexes
 * routinely run to tens of gigabytes.
 *
 * KreppPlacer holds the per-thread scratch state and refers to the shared
 * index. Construct one per thread.
 *
 * Edge numbering follows the backbone tree. When the supplied Newick carries
 * jplace-style {N} decorations those numbers are reported verbatim, which is
 * what lets placements accumulate against a phylogeny whose edge IDs are
 * assigned elsewhere and must stay stable.
 */

#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <string>
#include <vector>

namespace miint {

// Tuning parameters for a placement run; defaults match krepp's own.
struct KreppConfig {
	uint32_t hdist_th = 4; // max Hamming distance for a k-mer to match
	uint32_t tau = 2;      // Hamming distance threshold for the placement filter
	double chisq = 2.706;  // chi-square cutoff for distinguishability (alpha=90%)
	bool multi = true;     // report every candidate edge, not just the best
	bool filter = true;    // drop placements without enough k-mer support
};

// One query sequence to place.
struct KreppQuery {
	std::string id;
	std::string sequence;
};

// One placement of one query onto one edge of the backbone tree.
//
// Column names deliberately follow read_jplace, but the two are not
// interchangeable: read_jplace reports only the best placement per fragment
// and carries a filepath, while this carries krepp's `distance` and, with
// multi enabled, one row per candidate edge.
struct KreppPlacement {
	std::string fragment;
	int64_t edge_num = 0;
	double likelihood = 0.0;
	double like_weight_ratio = 0.0;
	double distal_length = 0.0;
	double pendant_length = 0.0;
	double distance = 0.0;
};

// A loaded krepp index together with the backbone tree placements are
// reported against. Immutable after construction.
class SharedKreppIndex {
public:
	// Loads every partial index under index_dir. If newick_path is non-empty
	// that tree becomes the backbone, overriding any tree stored in the index;
	// otherwise the index must carry its own.
	SharedKreppIndex(const std::string &index_dir, const std::string &newick_path);
	~SharedKreppIndex();

	SharedKreppIndex(const SharedKreppIndex &) = delete;
	SharedKreppIndex &operator=(const SharedKreppIndex &) = delete;

	// k-mer length the index was built with. Queries shorter than this
	// cannot produce a placement.
	uint32_t kmer_length() const;

	// Labeled tips in the backbone tree, and how many of those the index also
	// knows. Both 0 when the index supplied its own backbone. A tip the index
	// does not know contributes no placements and krepp says nothing about it,
	// so a caller that wants to distinguish "nothing placed" from "wrong tree"
	// has to compare these. All-or-nothing mismatch is rejected in the
	// constructor; a partial one is the caller's to judge.
	size_t backbone_tips_total() const;
	size_t backbone_tips_matched() const;

private:
	friend class KreppPlacer;
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

// Places query sequences against a SharedKreppIndex. Not thread-safe;
// construct one per thread.
class KreppPlacer {
public:
	KreppPlacer(std::shared_ptr<SharedKreppIndex> index, const KreppConfig &config);
	~KreppPlacer();

	KreppPlacer(const KreppPlacer &) = delete;
	KreppPlacer &operator=(const KreppPlacer &) = delete;

	// Appends placements for `queries` to `out` and returns how many queries
	// were skipped for being shorter than the index k-mer length.
	//
	// All or nothing: on failure this throws and `out` is left exactly as it
	// was, so a caller accumulating across batches never sees half of one.
	//
	// Rows are not one-to-one with queries in either direction: a query whose
	// placements are all rejected by the filter contributes none, and with
	// `multi` enabled a query contributes one row per candidate edge, with
	// like_weight_ratio summing to 1 across them.
	//
	// A sequence too short to place is skipped and counted; one containing a
	// character outside the nucleotide alphabet throws. The asymmetry is
	// deliberate - a short read is ordinary input, a non-nucleotide byte means
	// the column does not hold what the caller thinks it holds.
	size_t place(const std::vector<KreppQuery> &queries, std::vector<KreppPlacement> &out);

private:
	struct Impl;
	std::unique_ptr<Impl> impl_;
};

// Exposed for testing. Not part of the supported surface.
namespace krepp_detail {

// IUPAC nucleotide codes, both cases. krepp reads every other byte as N
// (seq_nt4_table, ext/krepp/src/common.cpp:11), so a stray character would be
// placed as an ambiguous base rather than reported. Shared with the index builder, which needs the same rule for
// reference sequences: a stray newline followed by '>' in a sequence would
// otherwise open a second FASTA record inside a reference's file, and krepp
// folds every record in that file into the same leaf.
extern const char *const kNucleotideAlphabet;

// Rejects the Newick shapes krepp's split_nwk refuses, before krepp ever sees
// the file. krepp reports every one of them by calling error_exit, which is
// std::exit, so without this a perfectly ordinary tree - one indented, or with
// CRLF line endings, or carrying a [&R] marker - would terminate the host
// process instead of raising.
//
// Deliberately stricter than krepp on whitespace: krepp errors on it only when
// a token is already accumulating, and otherwise folds the character into the
// following taxon name, which then fails to match the index and is dropped
// silently. Rejecting covers the quiet case as well as the fatal one.
//
// Exposed for testing; the shapes it rejects are the point, and they need no
// index to exercise. Throws miint::InvalidInputException naming `path`.
void ValidateNewickLexically(const std::string &newick_text, const std::string &path);

// Checks that `index_dir` holds exactly one krepp index, and returns its
// partials keyed by filename suffix. krepp keeps this discovery in krepp.cpp
// beside main(), which this build excludes; reimplementing it is also the only
// chance to convert these failures into exceptions, since krepp reports a
// missing file or an incompatible pair with error_exit.
//
// Decides on filenames alone and opens nothing - which is what makes it
// testable, and what lets it run before krepp can exit the process.
// Throws std::runtime_error naming `index_dir`.
std::map<std::string, std::set<std::string>> ValidateIndexLayout(const std::string &index_dir);

// The part of a partial's filename suffix that must AGREE across the partials of
// one index, with the part that is expected to differ removed.
//
// krepp builds the suffix as "-m<M>r<R>" + ("-frac"|"-no_frac")
// (ext/krepp/src/index.cpp:406-408). Index::load_partial_index reads k, w, h, m,
// r and frac out of each partial's metadata, and LSHF::check_compatible requires
// m, h, k, frac and the LSH positions to match, plus r when frac is true
// (ext/krepp/src/lshf.cpp:161-168). So what has to agree depends on frac:
//   - frac := false: r is the residue that says WHICH partial this is, and is
//     removed. Splitting the residues of one index across jobs and pointing
//     krepp at the collected directory is the supported way to build a large
//     index, so treating a differing r as a different index refuses layouts that
//     work.
//   - frac := true: r is a cumulative threshold (ext/krepp/src/rqseq.cpp:126)
//     and is kept, so a frac := true index is a single partial.
//
// Returns the suffix unchanged if it does not have that shape, so an
// unrecognised name groups only with itself rather than being merged.
std::string PartialHashConfig(const std::string &suffix);

// The filename suffix krepp gives every file of one partial. The single place
// this is built; ext/krepp/src/index.cpp:406-408 is what it mirrors.
std::string PartialSuffix(uint32_t m, uint32_t r, bool frac);

// Reads k, w and h out of a partial's BINARY `metadata<suffix>`.
//
// Not the `metadata<suffix>.txt` sidecar. The sidecar is advisory: krepp writes
// it after the binary (ext/krepp/src/index.cpp:567-577), synthesises it when it
// is missing (:186-199), and never consults it for compatibility. The binary is
// mandatory - load_partial_index error_exits without it (:122-125) - it is what
// LSHF::check_compatible is fed from (:127-140), and the completeness check above
// already requires it, so a complete partial always has one. Reading the
// sidecar instead meant a partial that had lost its .txt was silently skipped,
// which mattered most for `w`: krepp compares m, h, k, frac, r under
// frac := true, and the positions, but never w (ext/krepp/src/lshf.cpp:161-168),
// so this is the only w check there
// is, and skipping it produced a wrong placement result with no error at all.
//
// The layout is fixed and carries no version marker, so m, r and frac are read
// back too and checked against `expect_suffix`, which the filename already
// gave us. Those three are the format's own witness: if they disagree, the
// layout has moved and k/w/h cannot be trusted either.
//
// Throws std::runtime_error if the file cannot be read or the witness fails.
void ReadPartialConfig(const std::string &path, const std::string &expect_suffix, int32_t &k, int32_t &w, int32_t &h);

} // namespace krepp_detail

} // namespace miint
