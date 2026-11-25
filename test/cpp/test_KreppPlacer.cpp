#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cstdint>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "KreppPlacer.hpp"

using miint::krepp_detail::ValidateIndexLayout;
using miint::krepp_detail::ValidateNewickLexically;

// krepp's Newick reader rejects several shapes that miint's own parser accepts,
// and it rejects them by calling error_exit - std::exit, which would take the
// DuckDB process down with no SQL error at all. Every shape below was confirmed
// to do exactly that before this check existed.
//
// These tests exist because that check is the only thing standing between a
// user's tree and a dead server, and because every shape here is an ORDINARY
// file: trees that other phylogenetics tools write by default. They need no
// index, which is what makes them unit-testable at all - the guard runs inside
// SharedKreppIndex's constructor, and everything else there needs 69 MB of
// index on disk.

namespace {

// The shape of a tree krepp accepts: one line, no comments, no support values,
// exactly one trailing newline.
const char *const kValid = "((A:1,B:1):1,(C:1,D:1):1);\n";

} // namespace

TEST_CASE("ValidateNewickLexically accepts a well-formed tree", "[krepp]") {
	// Anti-vacuity for every rejection below: if this threw, the tests would all
	// pass without discriminating anything.
	CHECK_NOTHROW(ValidateNewickLexically(kValid, "t.nwk"));
	// krepp pops exactly one trailing newline, so a file with none is fine too.
	CHECK_NOTHROW(ValidateNewickLexically("((A:1,B:1):1,(C:1,D:1):1);", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically rejects trailing content after the ';'", "[krepp]") {
	// krepp pops ONE trailing newline and then insists the last character is
	// ';'. A blank line at end of file, or CRLF endings, leaves something else
	// there. Both are what an editor or a Windows-authored file produces.
	CHECK_THROWS(ValidateNewickLexically("((A:1,B:1):1,(C:1,D:1):1);\n\n", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("((A:1,B:1):1,(C:1,D:1):1);\r\n", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("((A:1,B:1):1,(C:1,D:1):1); ", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically rejects unquoted brackets", "[krepp]") {
	// krepp's reader has no comment support, so a [&R] rooted marker or an
	// inline comment - both standard Newick - are fatal.
	CHECK_THROWS(ValidateNewickLexically("[&R] ((A:1,B:1):1,(C:1,D:1):1);\n", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("((A:1,B:1)[c]:1,(C:1,D:1):1);\n", "t.nwk"));
	// Inside a quoted label they are ordinary characters and must pass.
	CHECK_NOTHROW(ValidateNewickLexically("(('A[1]':1,B:1):1,(C:1,D:1):1);\n", "t.nwk"));
	// krepp treats " as a quote character too, not just '.
	CHECK_NOTHROW(ValidateNewickLexically("((\"A[1]\":1,B:1):1,(C:1,D:1):1);\n", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically rejects a second tree in the file", "[krepp]") {
	// A ';' anywhere but the very end means more than one tree - a file of
	// bootstrap replicates, most often. krepp reads the first and exits.
	CHECK_THROWS(ValidateNewickLexically("(A:1,B:1);\n(C:1,D:1);\n", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically rejects unquoted whitespace", "[krepp]") {
	// This one is deliberately stricter than krepp. krepp errors on whitespace
	// only when a token is already accumulating; otherwise it folds the
	// character into the FOLLOWING label (ext/krepp/src/phytree.cpp:141-144), so
	// "(A:1, B:1)" silently yields a tip named " B" which then fails to match
	// the index and is dropped without a word (phytree.cpp:507-514). Rejecting
	// covers the quiet corruption as well as the crash.
	CHECK_THROWS(ValidateNewickLexically("(\n  (A:1,B:1):1,\n  (C:1,D:1):1\n);\n", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("((A:1, B:1):1,(C:1,D:1):1);\n", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("((A:1,B:1):1,\t(C:1,D:1):1);\n", "t.nwk"));
	// Quoted labels may contain spaces; krepp handles those.
	CHECK_NOTHROW(ValidateNewickLexically("(('Homo sapiens':1,B:1):1,(C:1,D:1):1);\n", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically rejects an empty tree", "[krepp]") {
	// krepp error_exits on an empty nwk_str before anything else.
	CHECK_THROWS(ValidateNewickLexically("", "t.nwk"));
	CHECK_THROWS(ValidateNewickLexically("\n", "t.nwk"));
}

TEST_CASE("ValidateNewickLexically names the offending file", "[krepp]") {
	// The message has to say which tree, because place_krepp takes the path from
	// a named parameter and the user may be passing several.
	try {
		ValidateNewickLexically("(A:1,B:1)", "/data/backbone.nwk");
		FAIL("expected a throw");
	} catch (const std::exception &e) {
		CHECK(std::string(e.what()).find("/data/backbone.nwk") != std::string::npos);
	}
}

// Index discovery. krepp keeps this in krepp.cpp next to main(), which this
// build excludes, so KreppPlacer reimplements it - and the reimplementation is
// the only place that can convert a discovery failure into an exception. Every
// case below reaches krepp as error_exit(), i.e. std::exit, and would kill the
// DuckDB process instead of failing the query.
//
// ValidateIndexLayout decides on filenames alone and opens nothing, which is
// what makes these testable: empty files are enough, and the legitimate case
// can be asserted to PASS without loading an index that is not there. Going
// through SharedKreppIndex instead would reach krepp's loader, and a
// well-formed-but-empty tree file exits the process - measured.
namespace {

// The five extensionless files that make up one complete krepp index with a
// backbone tree. `suffix` is krepp's own: "-m<M>r<R>" (the hash configuration,
// which every partial of one index shares) followed by a partial id.
// Writes one complete partial. `metadata<suffix>` is a REAL binary header in
// krepp's own layout (ext/krepp/src/index.cpp:127-140, written by
// save_configuration at :567-570), not a stub: ValidateIndexLayout reads k, w
// and h back out of that file exactly the way krepp does, so a placeholder
// would fail for the wrong reason.
void WriteIndexFiles(const std::filesystem::path &dir, uint32_t m, uint32_t r, bool frac, uint8_t k = 25,
                     uint8_t w = 31, uint8_t h = 9) {
	std::filesystem::create_directories(dir);
	const std::string suffix = miint::krepp_detail::PartialSuffix(m, r, frac);
	for (const char *type : {"cmer", "crecord", "inc", "tree"}) {
		std::ofstream(dir / (std::string(type) + suffix)).put('\0');
	}
	std::ofstream meta(dir / ("metadata" + suffix), std::ofstream::binary);
	const uint32_t nrows = 1;
	meta.write(reinterpret_cast<const char *>(&k), sizeof(uint8_t));
	meta.write(reinterpret_cast<const char *>(&w), sizeof(uint8_t));
	meta.write(reinterpret_cast<const char *>(&h), sizeof(uint8_t));
	meta.write(reinterpret_cast<const char *>(&m), sizeof(uint32_t));
	meta.write(reinterpret_cast<const char *>(&r), sizeof(uint32_t));
	meta.write(reinterpret_cast<const char *>(&frac), sizeof(bool));
	meta.write(reinterpret_cast<const char *>(&nrows), sizeof(uint32_t));
	const std::vector<uint8_t> ppos(h, 0), npos(static_cast<size_t>(k - h), 0);
	meta.write(reinterpret_cast<const char *>(ppos.data()), static_cast<std::streamsize>(ppos.size()));
	meta.write(reinterpret_cast<const char *>(npos.data()), static_cast<std::streamsize>(npos.size()));
}

// A directory nothing else in this file shares, cleared on the way in so a
// previous run cannot leak into this one.
std::filesystem::path FreshDir(const std::string &name) {
	const std::filesystem::path dir = std::filesystem::temp_directory_path() / ("miint_krepp_" + name);
	std::filesystem::remove_all(dir);
	return dir;
}

} // namespace

TEST_CASE("ValidateIndexLayout rejects two hash configurations in one directory", "[krepp]") {
	// `krepp index` writes into an existing directory without clearing it, so
	// re-indexing with different -h/-w leaves both file sets behind. krepp loads
	// every complete group it finds and only notices the mismatch inside its
	// partial loaders (ext/krepp/src/index.cpp:78, :116, :157), which report it with
	// error_exit. The suffix says it first.
	// frac := false, so m is the only thing separating the two: under frac := true
	// their different r would reject them on its own.
	const std::filesystem::path dir = FreshDir("twocfg");
	WriteIndexFiles(dir, 4, 1, false);
	WriteIndexFiles(dir, 8, 2, false);

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()),
	                    Catch::Matchers::ContainsSubstring("different hash configurations"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout accepts several partials of one index", "[krepp]") {
	// The counterpart, and the reason the check keys on the hash configuration
	// rather than the whole suffix: multiple partials are the normal layout for
	// a large index, and they differ precisely in the residue. Keying on the
	// whole suffix would reject this, which is a worse failure than the one
	// being prevented - it would refuse indexes that work.
	//
	// These are the real shapes. krepp builds the suffix as "-m<M>r<R>" plus
	// "-frac" or "-no_frac" (ext/krepp/src/index.cpp:406-408), and one residue
	// per job into a shared directory is how a large index is meant to be
	// built - with frac := false, the only setting under which krepp loads
	// several partials together (frac := true is the next test). An earlier
	// version of this test used "-m4r1-frac" against "-m4r1-frac2" - a shape
	// krepp never emits - and so passed while the code rejected every genuine
	// multi-partial index.
	const std::filesystem::path dir = FreshDir("onecfg");
	WriteIndexFiles(dir, 4, 1, false);
	WriteIndexFiles(dir, 4, 2, false);
	WriteIndexFiles(dir, 4, 3, false);

	REQUIRE_NOTHROW(ValidateIndexLayout(dir.string()));
	REQUIRE(ValidateIndexLayout(dir.string()).size() == 3);
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout rejects several frac := true partials", "[krepp]") {
	// Under frac := true, r is a cumulative threshold rather than a residue
	// (ext/krepp/src/rqseq.cpp:126), and krepp refuses to load two partials that
	// differ in it (LSHF::check_compatible, ext/krepp/src/lshf.cpp:161-168), which
	// it reports only once the files are open. Two frac := true partials of one m
	// are two indexes, and the layout check says so before anything is opened.
	const std::filesystem::path dir = FreshDir("twofrac");
	WriteIndexFiles(dir, 4, 1, true);
	WriteIndexFiles(dir, 4, 2, true);

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()),
	                    Catch::Matchers::ContainsSubstring("different hash configurations (m4r1-frac, m4r2-frac)"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout rejects partials that disagree on w", "[krepp]") {
	// w gets its own case because this check is the ONLY one anywhere. krepp
	// compares m, h, k, frac, r under frac := true, and the two position vectors,
	// and never w (ext/krepp/src/lshf.cpp:161-168), and load_partial_index checks each
	// partial's w only against its own k (LSHF::check_configuration, index.cpp:137).
	// So unlike k and h, a differing w
	// produces no error at any later point - just one index quietly holding two
	// different sets of minimizers, and a placement result that is wrong rather
	// than absent. Measured on a real pair: 148 rows against 140, with 13
	// placements present in the consistent index and missing from the mixed one.
	const std::filesystem::path dir = FreshDir("mixedw");
	WriteIndexFiles(dir, 4, 1, false, 25, 40, 9);
	WriteIndexFiles(dir, 4, 2, false, 25, 31, 9);

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()), Catch::Matchers::ContainsSubstring("disagree on k, w or h"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout rejects partials that disagree on k or h", "[krepp]") {
	// The loud half of the same defect. krepp does catch these, from inside
	// LSHF::check_compatible - but only once the files are open, which for a
	// real index means after the build that produced them. Deciding it here
	// keeps the failure at the point where the mistake is still fixable.
	const std::filesystem::path dir = FreshDir("mixedh");
	WriteIndexFiles(dir, 4, 1, false, 25, 31, 10);
	WriteIndexFiles(dir, 4, 2, false, 25, 31, 9);

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()), Catch::Matchers::ContainsSubstring("disagree on k, w or h"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout rejects metadata that contradicts its own filename", "[krepp]") {
	// The binary header has no version marker, so reading k, w and h from fixed
	// offsets is only safe while the layout holds. m, r and frac are the three
	// fields whose values the filename already states, which makes them the
	// format's own witness: writing a header that says m=8 into a file named
	// -m4r1-no_frac is indistinguishable from krepp having moved the fields,
	// and either way k/w/h would be read from the wrong place.
	const std::filesystem::path dir = FreshDir("badwitness");
	WriteIndexFiles(dir, 4, 1, false);
	// Same file set, but the header inside describes a different partial.
	std::filesystem::remove(dir / "metadata-m4r1-no_frac");
	const std::filesystem::path other = FreshDir("badwitness_src");
	WriteIndexFiles(other, 8, 1, false);
	std::filesystem::copy_file(other / "metadata-m8r1-no_frac", dir / "metadata-m4r1-no_frac");

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()),
	                    Catch::Matchers::ContainsSubstring("metadata layout has changed"));
	std::filesystem::remove_all(dir);
	std::filesystem::remove_all(other);
}

TEST_CASE("PartialSuffix builds the shape krepp writes", "[krepp]") {
	// The one place the suffix is constructed, mirroring
	// ext/krepp/src/index.cpp:406-408. PartialHashConfig parses what this emits
	// and ReadPartialConfig checks a header against it, so all three move
	// together or none of them do.
	CHECK(miint::krepp_detail::PartialSuffix(4, 1, true) == "-m4r1-frac");
	CHECK(miint::krepp_detail::PartialSuffix(4, 1, false) == "-m4r1-no_frac");
	CHECK(miint::krepp_detail::PartialSuffix(64, 10, true) == "-m64r10-frac");
}

TEST_CASE("PartialHashConfig keeps m and frac and drops a frac := false residue", "[krepp]") {
	// What must agree across partials, and what may differ. Under frac := false
	// krepp compares the hash configuration without r (LSHF::check_compatible,
	// ext/krepp/src/lshf.cpp:161-168); frac := true is the next test.
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-no_frac") ==
	      miint::krepp_detail::PartialHashConfig("-m4r2-no_frac"));
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-no_frac") ==
	      miint::krepp_detail::PartialHashConfig("-m4r16-no_frac"));
	// m, and the frac flag, are part of the identity. The m checks use -no_frac,
	// the one shape whose residue is stripped, so the stripping is what they test.
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-no_frac") !=
	      miint::krepp_detail::PartialHashConfig("-m8r1-no_frac"));
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-frac") !=
	      miint::krepp_detail::PartialHashConfig("-m4r1-no_frac"));
	// Multi-digit m must not be confused with the residue.
	CHECK(miint::krepp_detail::PartialHashConfig("-m64r1-no_frac") !=
	      miint::krepp_detail::PartialHashConfig("-m6r1-no_frac"));
	// ...and must still group with its OWN residues. The line above only rules
	// out a false merge; without this one a false split goes unnoticed. Parsing
	// a single digit of m instead of the whole run passes every other assertion
	// in this file and in the SQL tests - all of which use a one-digit m - and
	// then rejects a real two-residue m := 64 build as two different indexes.
	CHECK(miint::krepp_detail::PartialHashConfig("-m64r1-no_frac") ==
	      miint::krepp_detail::PartialHashConfig("-m64r2-no_frac"));
	// Anything that is not that shape groups only with itself, rather than
	// being merged with a partial it has nothing to do with.
	CHECK(miint::krepp_detail::PartialHashConfig("-nonsense") == "-nonsense");
	CHECK(miint::krepp_detail::PartialHashConfig("-mr1-frac") == "-mr1-frac");
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r-frac") == "-m4r-frac");
}

TEST_CASE("PartialHashConfig keeps the residue of a frac := true partial", "[krepp]") {
	// What must agree depends on frac. Under frac := false, r says which partial
	// this is and is dropped. Under frac := true it is a threshold that krepp
	// compares across partials (ext/krepp/src/lshf.cpp:161-168), so it stays part
	// of the identity.
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-frac") == "-m4r1-frac");
	CHECK(miint::krepp_detail::PartialHashConfig("-m4r1-frac") != miint::krepp_detail::PartialHashConfig("-m4r2-frac"));
	CHECK(miint::krepp_detail::PartialHashConfig("-m64r1-frac") !=
	      miint::krepp_detail::PartialHashConfig("-m64r2-frac"));
}

TEST_CASE("ValidateIndexLayout rejects an incomplete index", "[krepp]") {
	// One file short of a complete group. krepp's own response is
	// error_exit("There is a partial index with a missing file!").
	const std::filesystem::path dir = FreshDir("incomplete");
	std::filesystem::create_directories(dir);
	for (const char *type : {"cmer", "crecord", "inc"}) {
		std::ofstream(dir / (std::string(type) + "-m4r1-frac")).put('\0');
	}

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()),
	                    Catch::Matchers::ContainsSubstring("missing one or more files"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout rejects an index path that is not a directory", "[krepp]") {
	// The counterpart to the missing-directory case: is_directory returns false
	// with the error code CLEAR for a regular file, so this used to be reported
	// as "does not exist". Pointing index_path at a file is a different mistake
	// from mistyping it and deserves a different message.
	const std::filesystem::path file = FreshDir("notadir");
	std::filesystem::create_directories(file.parent_path());
	std::ofstream(file).put('\0');

	REQUIRE_THROWS_WITH(ValidateIndexLayout(file.string()), Catch::Matchers::ContainsSubstring("is not a directory"));
	std::filesystem::remove_all(file);
}

TEST_CASE("ValidateIndexLayout distinguishes a missing directory from an empty one", "[krepp]") {
	// Two different mistakes with two different fixes: a typo in the path, and a
	// path that is right but was never indexed. krepp conflates neither because
	// krepp never gets this far.
	const std::filesystem::path missing = FreshDir("missing");
	REQUIRE_THROWS_WITH(ValidateIndexLayout(missing.string()), Catch::Matchers::ContainsSubstring("does not exist"));

	const std::filesystem::path empty = FreshDir("empty");
	std::filesystem::create_directories(empty);
	REQUIRE_THROWS_WITH(ValidateIndexLayout(empty.string()),
	                    Catch::Matchers::ContainsSubstring("No krepp index found"));
	std::filesystem::remove_all(empty);
}
