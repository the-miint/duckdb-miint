#include "KreppPlacer.hpp"

#include "NewickTree.hpp"
#include "alignment_functions_internal.hpp"

#include <algorithm>
#include <cctype>
#include <cstdint>
#include <deque>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <system_error>

// krepp's own headers. krepp.hpp is deliberately not included: it pulls in
// CLI11 and the command-line layer, neither of which applies here.
#include "index.hpp"
#include "phytree.hpp"
#include "query.hpp"
#include "rqseq.hpp"

// normalize_rna lives here for the vsearch wrappers but has no vsearch
// dependency; see the U/u note in place() for why krepp needs it too.
#include "vsearch_utils.hpp"

namespace miint {

namespace {

// The file kinds that make up one partial index, and the two complete sets.
const std::set<std::string> kIndexFileTypes = {"cmer", "crecord", "inc", "metadata", "tree", "reflist"};
const std::set<std::string> kWithBackbone = {"cmer", "crecord", "inc", "metadata", "tree"};
const std::set<std::string> kWithoutBackbone = {"cmer", "crecord", "inc", "metadata", "reflist"};

// IUPAC nucleotide codes, both cases. krepp reads every other byte as N
// (seq_nt4_table, ext/krepp/src/common.cpp:11), so a stray character would be
// placed as an ambiguous base rather than reported.
//
// Groups the files in index_dir by their "-m4r1-frac" style suffix. Mirrors how
// krepp names the pieces of a partial index.
std::map<std::string, std::set<std::string>> DiscoverPartials(const std::string &index_dir) {
	std::map<std::string, std::set<std::string>> suffix_to_types;
	for (const auto &entry : std::filesystem::directory_iterator(index_dir)) {
		// Only the extensionless binaries participate; metadata-*.txt is prose.
		if (!entry.path().extension().empty()) {
			continue;
		}
		const std::string filename = entry.path().filename().string();
		const size_t first = filename.find('-', 0);
		if (first == std::string::npos) {
			continue;
		}
		const std::string file_type = filename.substr(0, first);
		if (kIndexFileTypes.find(file_type) == kIndexFileTypes.end()) {
			continue;
		}
		suffix_to_types[filename.substr(first)].insert(file_type);
	}
	return suffix_to_types;
}

} // namespace

namespace krepp_detail {

std::string PartialHashConfig(const std::string &suffix) {
	// "-m<M>r<R><rest>" -> "-m<M><rest>", except that a "-frac" suffix keeps its r
	// (see the header). Anything else is returned unchanged.
	if (suffix.size() < 2 || suffix[0] != '-' || suffix[1] != 'm') {
		return suffix;
	}
	size_t i = 2;
	while (i < suffix.size() && std::isdigit(static_cast<unsigned char>(suffix[i]))) {
		i++;
	}
	if (i == 2 || i >= suffix.size() || suffix[i] != 'r') {
		return suffix;
	}
	const size_t r_begin = i++;
	while (i < suffix.size() && std::isdigit(static_cast<unsigned char>(suffix[i]))) {
		i++;
	}
	if (i == r_begin + 1) {
		return suffix; // an 'r' with no digits after it is not a residue
	}
	if (suffix.compare(i, std::string::npos, "-frac") == 0) {
		return suffix;
	}
	return suffix.substr(0, r_begin) + suffix.substr(i);
}

// U and u are accepted only because place() rewrites them to T and t before
// krepp sees them; krepp's own table maps them to the ambiguous code, so an RNA
// read would otherwise validate here and then place nothing.
const char *const kNucleotideAlphabet = "ACGTURYKMSWBDHVNacgturykmswbdhvn";

void ValidateNewickLexically(const std::string &newick_text, const std::string &path) {
	// krepp's split_nwk runs before any parsing and rejects several shapes
	// that miint's parser accepts. Mirrored here because krepp reports them
	// through error_exit, which is std::exit - measured, each of these killed
	// the DuckDB process outright with no SQL error at all: a trailing blank
	// line, CRLF line endings, a `[&R]` rooted marker, an inline `[comment]`,
	// an indented tree, and a two-tree file.
	//
	// The whitespace rule below is the exception and is broader than krepp's
	// on purpose; see the comment on it.
	//
	// The terminating-';' rule is krepp's own, quirk included: it pops
	// exactly ONE trailing newline and then insists on ';'. So ";\n" passes
	// and ";\n\n" or ";\r\n" does not.
	std::string tail = newick_text;
	if (!tail.empty() && tail.back() == '\n') {
		tail.pop_back();
	}
	if (tail.empty() || tail.back() != ';') {
		throw miint::InvalidInputException(
		    "Backbone tree " + path +
		    " must end with ';' followed by at most one newline; krepp rejects anything else. "
		    "A trailing blank line or CRLF line endings are the usual causes.");
	}
	// One pass over the body, quotes tracked so a quoted label may contain
	// anything. krepp treats both ' and " as quote characters.
	//
	// Deliberately simpler than krepp's quote handling rather than strictly
	// stricter: krepp treats a doubled quote as an escape and clears its
	// state, where this toggles, so the two diverge in both directions on
	// pathological input. NewickTree::parse has already rejected anything
	// with an unterminated quote by this point, which is what keeps the
	// difference theoretical.
	bool quoted = false;
	for (size_t i = 0; i < tail.size(); ++i) {
		const char c = tail[i];
		if (c == '\'' || c == '"') {
			quoted = !quoted;
			continue;
		}
		if (quoted) {
			continue;
		}
		if (c == '[' || c == ']') {
			throw miint::InvalidInputException(
			    "Backbone tree " + path +
			    " contains an unquoted '[' or ']'; krepp's Newick reader has no comment support, so "
			    "rooted markers like [&R] and inline comments must be removed.");
		}
		if (c == ';' && i + 1 != tail.size()) {
			throw miint::InvalidInputException(
			    "Backbone tree " + path +
			    " contains a ';' before the end of the file; krepp accepts exactly one tree per file.");
		}
		// Any unquoted whitespace at all. krepp has no case where interior
		// whitespace is harmless: where it lands next to a token it exits,
		// and where it lands at the start of one it is absorbed INTO the
		// label - "(A:1, B:1)" names the second tip " B", which then fails
		// to match the index and is dropped by map_to_qtree without a word.
		// Rejecting covers the silent case as well as the fatal one.
		if (c == ' ' || c == '\t' || c == '\n' || c == '\r') {
			throw miint::InvalidInputException(
			    "Backbone tree " + path +
			    " contains whitespace outside a quoted label; krepp cannot read an indented or "
			    "line-wrapped Newick, and would silently fold a space into the adjacent taxon name. "
			    "Write the tree on one line.");
		}
	}
}

std::string PartialSuffix(uint32_t m, uint32_t r, bool frac) {
	return "-m" + std::to_string(m) + "r" + std::to_string(r) + (frac ? "-frac" : "-no_frac");
}

void ReadPartialConfig(const std::string &path, const std::string &expect_suffix, int32_t &k, int32_t &w, int32_t &h) {
	std::ifstream in(path, std::ifstream::binary);
	if (!in) {
		throw std::runtime_error("krepp metadata " + path + " could not be opened");
	}
	// Field order and widths are krepp's, read the same way it reads them
	// (ext/krepp/src/index.cpp:127-135). Only the first six are needed; nrows and
	// the position vectors follow and are left alone.
	uint8_t k_raw = 0, w_raw = 0, h_raw = 0;
	uint32_t m_raw = 0, r_raw = 0;
	bool frac_raw = false;
	in.read(reinterpret_cast<char *>(&k_raw), sizeof(uint8_t));
	in.read(reinterpret_cast<char *>(&w_raw), sizeof(uint8_t));
	in.read(reinterpret_cast<char *>(&h_raw), sizeof(uint8_t));
	in.read(reinterpret_cast<char *>(&m_raw), sizeof(uint32_t));
	in.read(reinterpret_cast<char *>(&r_raw), sizeof(uint32_t));
	in.read(reinterpret_cast<char *>(&frac_raw), sizeof(bool));
	if (!in) {
		throw std::runtime_error("krepp metadata " + path +
		                         " is shorter than krepp's own header; the file is "
		                         "truncated or was not written by krepp");
	}
	// The format's own witness. m, r and frac are the three fields whose value
	// is already known from the filename, so a disagreement means the layout has
	// moved underneath us and k, w and h are being read from the wrong offsets.
	const std::string found_suffix = PartialSuffix(m_raw, r_raw, frac_raw);
	if (found_suffix != expect_suffix) {
		throw std::runtime_error("krepp metadata " + path + " describes '" + found_suffix +
		                         "' but its filename says '" + expect_suffix +
		                         "'; krepp's metadata layout has changed");
	}
	k = static_cast<int32_t>(k_raw);
	w = static_cast<int32_t>(w_raw);
	h = static_cast<int32_t>(h_raw);
}

std::map<std::string, std::set<std::string>> ValidateIndexLayout(const std::string &index_dir) {
	// `ec` is reported rather than discarded, but only when it says something
	// new: is_directory returns false both for "not there" and for "there but
	// unreadable", and diagnosing a permission problem as a missing path sends
	// the caller looking in the wrong place. libc++ sets ec to ENOENT for the
	// ordinary missing case, so that one is filtered back out.
	std::error_code ec;
	if (!std::filesystem::is_directory(index_dir, ec)) {
		if (ec && ec != std::errc::no_such_file_or_directory) {
			throw std::runtime_error("krepp index directory is not readable: " + index_dir + " (" + ec.message() + ")");
		}
		// is_directory is false-with-ec-clear for a path that exists but is a
		// regular file, so "does not exist" would be the wrong diagnosis - a
		// different mistake with a different fix.
		if (std::filesystem::exists(index_dir)) {
			throw std::runtime_error("krepp index path is not a directory: " + index_dir);
		}
		throw std::runtime_error("krepp index directory does not exist: " + index_dir);
	}

	std::map<std::string, std::set<std::string>> partials = DiscoverPartials(index_dir);
	if (partials.empty()) {
		throw std::runtime_error("No krepp index found in directory: " + index_dir);
	}
	for (const auto &entry : partials) {
		const bool complete =
		    std::includes(entry.second.begin(), entry.second.end(), kWithBackbone.begin(), kWithBackbone.end()) ||
		    std::includes(entry.second.begin(), entry.second.end(), kWithoutBackbone.begin(), kWithoutBackbone.end());
		if (!complete) {
			throw std::runtime_error("krepp index '" + entry.first + "' in " + index_dir +
			                         " is missing one or more files");
		}
	}
	// Several partials of one index are normal, and krepp loads them all
	// (ext/krepp/src/api.cpp:44-93 collects every distinct suffix and calls
	// load_partial_index on each). Two *different* indexes in one directory are
	// not: `krepp index` never clears what is already there, so re-indexing with
	// new parameters leaves both behind. krepp discovers that only inside its
	// partial loaders (Index::generate_partial_tree, load_partial_tree and
	// load_partial_index), which report it with error_exit (index.cpp:78, :116,
	// :157). InstallKreppErrorHandler turns that into an exception rather
	// than a std::exit, so it surfaces as a DuckDB error - but only after the
	// files have been opened and read. It is decidable here, from the filenames,
	// before anything opens.
	//
	// What has to agree is the hash configuration - see PartialHashConfig, which
	// keeps the residue only under frac := true. An earlier version keyed on
	// `-m4r1` for every partial, which folded the residue into the identity and so
	// rejected the ordinary frac := false multi-partial layout: three residues of
	// one index read as three different indexes.
	std::set<std::string> hash_configs;
	for (const auto &entry : partials) {
		hash_configs.insert(PartialHashConfig(entry.first));
	}
	if (hash_configs.size() > 1) {
		std::string listed;
		for (const auto &config : hash_configs) {
			listed += (listed.empty() ? "" : ", ") + config.substr(1);
		}
		throw std::runtime_error("Directory " + index_dir +
		                         " holds krepp indexes built with different hash "
		                         "configurations (" +
		                         listed +
		                         "); krepp would try to load them as one index. Partials of one index share m "
		                         "and frac and differ only in r, which frac := true does not allow: a frac := true "
		                         "index is a single partial. Keep one index per directory.");
	}

	// k, w and h are in no filename, so two builds that differ in them are
	// indistinguishable above. That matters in both directions. krepp catches k
	// and h itself inside LSHF::check_compatible, but only once the files are
	// open - and it never compares w across partials (ext/krepp/src/lshf.cpp:161-168),
	// checking each partial's w only against its own k (LSHF::check_configuration,
	// index.cpp:137, lshf.cpp:18), so a differing w is
	// caught nowhere and simply leaves one index holding two different sets of
	// minimizers. Compare all three here, from the binary metadata krepp itself
	// reads, before anything is loaded.
	const std::string *named = nullptr;
	int32_t first_k = 0, first_w = 0, first_h = 0;
	for (const auto &entry : partials) {
		int32_t k = 0, w = 0, h = 0;
		// Every partial, unconditionally. The binary metadata is mandatory and
		// the completeness check above already required it, so there is no
		// absent case to tolerate - and tolerating one would delete the only
		// check on w that exists anywhere.
		ReadPartialConfig(index_dir + "/metadata" + entry.first, entry.first, k, w, h);
		if (named == nullptr) {
			named = &entry.first;
			first_k = k;
			first_w = w;
			first_h = h;
			continue;
		}
		if (k == first_k && w == first_w && h == first_h) {
			continue;
		}
		throw std::runtime_error("Directory " + index_dir + " holds partials that disagree on k, w or h: '" + *named +
		                         "' has k=" + std::to_string(first_k) + " w=" + std::to_string(first_w) + " h=" +
		                         std::to_string(first_h) + ", '" + entry.first + "' has k=" + std::to_string(k) +
		                         " w=" + std::to_string(w) + " h=" + std::to_string(h) +
		                         ". Partials of one index must agree on all three; these are two different indexes.");
	}
	return partials;
}

} // namespace krepp_detail

struct SharedKreppIndex::Impl {
	index_sptr_t index;
	tree_sptr_t backbone;
	uint32_t k = 0;
	// Backbone tip overlap with the index, both 0 when no newick_path was given.
	size_t backbone_tips_total = 0;
	size_t backbone_tips_matched = 0;
};

SharedKreppIndex::SharedKreppIndex(const std::string &index_dir, const std::string &newick_path) : impl_(new Impl()) {
	// krepp reports failures through error_exit, which InstallKreppErrorHandler
	// turns into a thrown KreppFatalError rather than std::exit. Paths and the
	// index file set are still checked here, so the common mistakes are refused
	// before anything opens and with a message naming the file. Errors in the
	// *content* of an index - partials built against different trees or with
	// incompatible hash configurations - reach krepp's error_exit mid-load
	// (ext/krepp/src/index.cpp:78, :116, :157) and surface as that error.
	const auto partials = krepp_detail::ValidateIndexLayout(index_dir);
	std::error_code ec;

	// Parse the backbone with miint's own reader first. krepp's parser reports
	// malformed Newick through error_exit, so anything it would reject at that
	// level is caught here as an exception instead.
	//
	// Parsing is not enough on its own. krepp's Tree::load rejects two further
	// shapes, both by error_exit, and both of which miint's parser accepts:
	//   - a node with exactly one child (phytree.cpp:171-173), which
	//     NewickTree deliberately preserves; and
	//   - two nodes sharing a name (check_unique_labels, phytree.cpp:468-480).
	//     This one is not exotic. krepp's parser puts the token after ')' into
	//     Node::name, so an internal node's bootstrap support becomes its name,
	//     and any tree with the same support value on two clades - which is to
	//     say most trees carrying support - is a duplicate.
	// Both were confirmed to kill the process: '((A:1,B:1)100:1,(C:1,D:1)100:1);'
	// and '((A:1):1,B:1);' each exited the shell with no SQL error at all.
	// So they are checked here, where a throw is still possible.
	// Held open across both reads below: krepp's Tree::load takes an ifstream,
	// so the validated bytes and the parsed bytes have to come from one handle
	// rather than two opens of the same path.
	std::ifstream tree_stream;
	if (!newick_path.empty()) {
		if (!std::filesystem::is_regular_file(newick_path, ec)) {
			if (ec && ec != std::errc::no_such_file_or_directory) {
				throw std::runtime_error("Backbone tree file is not readable: " + newick_path + " (" + ec.message() +
				                         ")");
			}
			if (std::filesystem::exists(newick_path)) {
				throw std::runtime_error("Backbone tree path is not a regular file: " + newick_path);
			}
			throw std::runtime_error("Backbone tree file does not exist: " + newick_path);
		}
		tree_stream.open(newick_path);
		if (!tree_stream.good()) {
			throw std::runtime_error("Failed to open backbone tree: " + newick_path);
		}
		std::string newick_text((std::istreambuf_iterator<char>(tree_stream)), std::istreambuf_iterator<char>());
		// NewickTree::parse throws std::runtime_error, which the table function
		// maps to IOException - misleading for what is a malformed input, not an
		// I/O failure. Rethrow as an input error naming the file.
		//
		// Note this parser is narrower than krepp's on one point: it reads a
		// jplace {N} decoration only after the branch length (":1{1}", the form
		// the jplace spec uses), while krepp also accepts it straight after the
		// closing paren ("){1}:1"). A backbone in that second form is rejected
		// here even though krepp would read it. Widening NewickTree is the fix,
		// but it is shared with read_newick and every other phylogeny function,
		// so it is not this wrapper's call to make.
		std::unique_ptr<NewickTree> parsed;
		try {
			parsed = std::make_unique<NewickTree>(NewickTree::parse(newick_text));
		} catch (const std::exception &e) {
			throw miint::InvalidInputException("Backbone tree " + newick_path + " could not be parsed: " + e.what());
		}
		const NewickTree &tree = *parsed;

		krepp_detail::ValidateNewickLexically(newick_text, newick_path);
		std::set<std::string> seen_names;
		std::set<int64_t> seen_edge_ids;
		size_t decorated = 0;
		for (uint32_t node = 0; node < tree.num_nodes(); ++node) {
			// jplace {N} edge decorations. krepp requires every node decorated or
			// none (Tree::load, ext/krepp/src/phytree.cpp:461-464) and, since
			// v0.9.1, that the numbers are unique (:455-458). Both are error_exit.
			//
			// Safe to check with miint's parser because the one form it reads -
			// {N} after the branch length, the jplace spec's own - is a form krepp
			// reads identically. Verified on a fully decorated tree: miint reports
			// every node decorated and krepp accepts it; on a partially decorated
			// one miint reports 3 of 7 and krepp exits.
			const auto edge_id = tree.edge_id(node);
			if (edge_id.has_value()) {
				++decorated;
				// krepp narrows the decoration to uint32 (`static_cast<se_t>(std::atol(...))`,
				// ext/krepp/src/phytree.cpp:236, with se_t = uint32_t at common.hpp:69) while
				// miint's parser reads it as int64. Without this bound {-5} would be reported
				// back as 4294967291 - breaking the contract that these numbers are echoed
				// verbatim - and {-1} would collide with {4294967295} inside krepp while
				// looking distinct to the uniqueness check below, reaching the error_exit at
				// phytree.cpp:456.
				if (*edge_id < 0 || *edge_id > static_cast<int64_t>(std::numeric_limits<uint32_t>::max())) {
					throw miint::InvalidInputException("Backbone tree " + newick_path + " uses edge number " +
					                                   std::to_string(*edge_id) +
					                                   ", which is outside the range krepp represents (0 to " +
					                                   std::to_string(std::numeric_limits<uint32_t>::max()) +
					                                   "); krepp would silently wrap it to a different number.");
				}
				if (!seen_edge_ids.insert(*edge_id).second) {
					throw miint::InvalidInputException(
					    "Backbone tree " + newick_path + " uses edge number " + std::to_string(*edge_id) +
					    " on more than one node; krepp requires {N} decorations to be unique.");
				}
			}
			if (tree.children(node).size() == 1) {
				throw miint::InvalidInputException(
				    "Backbone tree " + newick_path +
				    " contains a unifurcation (a node with exactly one child); krepp rejects these. "
				    "Suppress unifurcations before passing the tree.");
			}
			// Empty names are skipped, matching krepp's own check.
			const std::string &node_name = tree.name(node);
			if (node_name.empty()) {
				continue;
			}
			if (!seen_names.insert(node_name).second) {
				throw miint::InvalidInputException(
				    "Backbone tree " + newick_path + " has more than one node named '" + node_name +
				    "'; krepp requires every node name to be unique. Note that krepp reads an internal node's "
				    "bootstrap support value as its name, so a tree carrying support values will collide "
				    "wherever two clades share one - strip them before passing the tree.");
			}
		}
		if (decorated > 0 && decorated != tree.num_nodes()) {
			throw miint::InvalidInputException(
			    "Backbone tree " + newick_path + " is only partially decorated with {N} edge numbers (" +
			    std::to_string(decorated) + " of " + std::to_string(tree.num_nodes()) +
			    " nodes); krepp requires all of them or none, so that one column cannot mix the caller's edge "
			    "numbers with krepp's internal ones.");
		}
		tree_stream.clear();
		tree_stream.seekg(0);
	}

	impl_->index = std::make_shared<Index>(index_dir);
	for (const auto &entry : partials) {
		const bool with_backbone =
		    std::includes(entry.second.begin(), entry.second.end(), kWithBackbone.begin(), kWithBackbone.end());
		if (with_backbone) {
			impl_->index->load_partial_tree(entry.first);
		} else {
			impl_->index->generate_partial_tree(entry.first);
		}
		impl_->index->load_partial_index(entry.first);
	}
	impl_->index->make_rho_partial();
	impl_->k = impl_->index->get_lshf()->get_k();

	if (!newick_path.empty()) {
		// Collect the index's leaf names BEFORE map_to_qtree, which reassigns the
		// index tree's root to the backbone's (phytree.cpp:511) and so cannot be
		// traversed for them afterwards.
		std::set<std::string> index_leaf_names;
		tree_sptr_t index_tree = impl_->index->get_tree();
		index_tree->reset_traversal();
		while (node_sptr_t node = index_tree->next_post_order()) {
			if (node->check_leaf()) {
				index_leaf_names.insert(node->get_name());
			}
		}

		impl_->backbone = std::make_shared<Tree>();
		impl_->backbone->load(tree_stream);

		// map_to_qtree silently skips backbone leaves absent from the index - its
		// only diagnostic is commented out (phytree.cpp:519-525). Every skipped
		// leaf leaves eff_nchildren at 0, and collect_placements then rejects
		// every candidate on the `get_nchildren() != get_eff_nchildren()` test
		// (query.cpp:339). A backbone whose labels do not match the index is
		// therefore not an error anywhere in krepp: it is zero rows, no warning,
		// indistinguishable from "nothing placed". Count the overlap ourselves.
		impl_->backbone->reset_traversal();
		while (node_sptr_t node = impl_->backbone->next_post_order()) {
			if (node->check_leaf() && node->is_labeled()) {
				++impl_->backbone_tips_total;
				if (index_leaf_names.count(node->get_name())) {
					++impl_->backbone_tips_matched;
				}
			}
		}
		// Nothing matched because nothing was named. Caught separately because
		// `backbone_tips_total == 0` is also how "no newick_path was given" is
		// represented to callers, and because the fix a caller needs differs:
		// label the tree, rather than relabel it.
		if (impl_->backbone_tips_total == 0) {
			throw miint::InvalidInputException(
			    "Backbone tree " + newick_path +
			    " has no labeled tips. krepp matches backbone leaves to the index by name, so an unlabeled tree "
			    "can place nothing and would return zero rows rather than fail.");
		}
		if (impl_->backbone_tips_matched == 0) {
			throw miint::InvalidInputException(
			    "Backbone tree " + newick_path + " shares no tip labels with krepp index '" + index_dir + "' (" +
			    std::to_string(impl_->backbone_tips_total) +
			    " tips, none matched). Every placement would be discarded, so this would return zero rows "
			    "rather than fail. Check that the tree's tip labels are the index's reference names.");
		}

		impl_->backbone->reset_traversal();
		index_tree->map_to_qtree(impl_->backbone);
	} else if (impl_->index->check_wbackbone()) {
		impl_->backbone = impl_->index->get_tree();
	} else {
		throw miint::InvalidInputException("krepp index '" + index_dir +
		                                   "' has no backbone tree; supply one with the newick_path parameter");
	}
}

SharedKreppIndex::~SharedKreppIndex() = default;

uint32_t SharedKreppIndex::kmer_length() const {
	return impl_->k;
}

size_t SharedKreppIndex::backbone_tips_total() const {
	return impl_->backbone_tips_total;
}

size_t SharedKreppIndex::backbone_tips_matched() const {
	return impl_->backbone_tips_matched;
}

struct KreppPlacer::Impl {
	std::shared_ptr<SharedKreppIndex> shared;
	KreppConfig config;
	// Built once per placer, and a placer is per-thread (LocalState owns it), so
	// these are never shared across workers. See the note in place() for why
	// reuse across batches is sound and why the QSeq exists at all.
	qseq_sptr_t qs;
	std::shared_ptr<IBatch> batch;
};

KreppPlacer::KreppPlacer(std::shared_ptr<SharedKreppIndex> index, const KreppConfig &config) : impl_(new Impl()) {
	impl_->shared = std::move(index);
	impl_->config = config;
	// gzopen of /dev/null is the one error_exit left on a non-index path
	// (rqseq.cpp:168-170, reached only under EMFILE). Constructing here rather
	// than per batch narrows that window from once per batch per worker to once
	// per worker, and takes the per-batch HDistHistLLH construction with it.
	impl_->qs = std::make_shared<QSeq>("/dev/null");
	impl_->batch = std::make_shared<IBatch>(impl_->shared->impl_->index, impl_->qs, impl_->config.hdist_th,
	                                        impl_->config.chisq, std::numeric_limits<double>::quiet_NaN(),
	                                        impl_->config.tau, !impl_->config.filter, impl_->config.multi,
	                                        /*summarize=*/false, /*p_value=*/false);
}

KreppPlacer::~KreppPlacer() = default;

size_t KreppPlacer::place(const std::vector<KreppQuery> &queries, std::vector<KreppPlacement> &out) {
	if (queries.empty()) {
		return 0;
	}

	// krepp reads any byte outside ACGT as N (nt4_of, ext/krepp/src/common.hpp:241),
	// so a stray character would be searched as an ambiguous base rather than
	// reported. Sequences here come from a SQL column rather than a curated
	// file, so reject instead.
	//
	// Anything shorter than k cannot place; search_mers computes `len - k + 1`
	// into an unsigned type, so passing those through would underflow.
	const uint32_t k = impl_->shared->kmer_length();
	size_t skipped_short = 0;
	std::vector<const std::string *> sequences;
	std::vector<const std::string *> ids;
	// Backing store for the RNA rewrite below. A deque, not a vector: `sequences`
	// holds pointers into it, and a vector would invalidate them on regrowth.
	std::deque<std::string> rewritten;
	sequences.reserve(queries.size());
	ids.reserve(queries.size());
	for (const KreppQuery &query : queries) {
		const size_t bad = query.sequence.find_first_not_of(krepp_detail::kNucleotideAlphabet);
		if (bad != std::string::npos) {
			throw miint::InvalidInputException(
			    "Query sequence for '" + query.id + "' contains a character that is not a nucleotide code: byte " +
			    std::to_string(static_cast<int>(static_cast<unsigned char>(query.sequence[bad]))) + " at offset " +
			    std::to_string(bad));
		}
		if (query.sequence.size() < k) {
			++skipped_short;
			continue;
		}
		// krepp maps U and u to the ambiguous code (ext/krepp/src/common.cpp:13-15),
		// which breaks the k-mer run at every uracil - an RNA read would pass the
		// check above and then place nothing, with no way to tell that apart from
		// "this read does not belong on this tree". The vsearch wrappers normalize
		// for the same reason. Only reads that actually carry a U are copied;
		// DNA is used in place.
		if (query.sequence.find_first_of("Uu") != std::string::npos) {
			rewritten.push_back(normalize_rna(query.sequence));
			sequences.push_back(&rewritten.back());
		} else {
			sequences.push_back(&query.sequence);
		}
		ids.push_back(&query.id);
	}
	if (sequences.empty()) {
		return skipped_short;
	}

	// Placements are staged here rather than appended straight to `out`: the
	// loop below can throw part way through a batch, and a caller accumulating
	// across batches must not be handed half of one.
	std::vector<KreppPlacement> staged;

	// krepp's own query loop (IBatch::place_sequences) reads sequences out of a
	// QSeq, formats each result as jplace text and writes it to a stream. We use
	// the three public steps underneath it instead - search_mers,
	// summarize_matches, collect_placements - which take the sequence as a
	// pointer and hand back `placement_t` structs. That skips the text entirely:
	// no FASTA to generate, no pipe to feed it through, and no jplace to parse
	// back. It also removes the whole class of field-order and precision bugs,
	// since the doubles are read straight out of the struct.
	//
	// The one wart is the QSeq that IBatch's constructor takes; it is built once
	// in this class's constructor, which is where the reasoning for it lives.
	// Both live on Impl, built once in the constructor. IBatch's constructor
	// reads `cbatch_size` off the QSeq and swaps away two empty vectors; it does
	// not store it, so the QSeq is inert from that point and one IBatch is good
	// for every batch this placer sees.
	IBatch &batch = *impl_->batch;

	vec<placement_t> placements;
	for (size_t i = 0; i < sequences.size(); ++i) {
		const std::string &sequence = *sequences[i];
		const uint64_t length = sequence.size();

		// Fresh per query, exactly as place_sequences does. summarize_matches
		// resets the rest of IBatch's per-query state (nd_closest, mi_closest,
		// node_to_minfo), which is what makes reusing one IBatch correct.
		imers_sptr_t imers_or = std::make_shared<IMers>(impl_->shared->impl_->index, length, impl_->config.hdist_th);
		imers_sptr_t imers_rc = std::make_shared<IMers>(impl_->shared->impl_->index, length, impl_->config.hdist_th);

		batch.search_mers(sequence.data(), length, imers_or, imers_rc);
		batch.summarize_matches(imers_or, imers_rc);
		if (!batch.collect_placements(placements)) {
			continue;
		}
		for (const placement_t &placement : placements) {
			KreppPlacement result;
			result.fragment = *ids[i];
			result.edge_num = static_cast<int64_t>(placement.edge_num);
			result.pendant_length = placement.pendant_length;
			result.distal_length = placement.distal_length;
			result.likelihood = placement.likelihood;
			result.like_weight_ratio = placement.like_weight_ratio;
			result.distance = placement.distance;
			staged.push_back(std::move(result));
		}
	}

	out.insert(out.end(), std::make_move_iterator(staged.begin()), std::make_move_iterator(staged.end()));
	return skipped_short;
}

} // namespace miint
