#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>

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
	// the index and is dropped without a word (phytree.cpp:488-495). Rejecting
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
void WriteIndexFiles(const std::filesystem::path &dir, const std::string &suffix) {
	std::filesystem::create_directories(dir);
	for (const char *type : {"cmer", "crecord", "inc", "metadata", "tree"}) {
		std::ofstream(dir / (std::string(type) + suffix)).put('\0');
	}
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
	// every complete group it finds and only notices the mismatch inside
	// Index::check_compatible, which reports it with error_exit. The suffix says
	// it first.
	const std::filesystem::path dir = FreshDir("twocfg");
	WriteIndexFiles(dir, "-m4r1-frac");
	WriteIndexFiles(dir, "-m8r2-frac");

	REQUIRE_THROWS_WITH(ValidateIndexLayout(dir.string()),
	                    Catch::Matchers::ContainsSubstring("different hash configurations"));
	std::filesystem::remove_all(dir);
}

TEST_CASE("ValidateIndexLayout accepts several partials of one index", "[krepp]") {
	// The counterpart, and the reason the check keys on the hash configuration
	// rather than the whole suffix: multiple partials are the normal layout for
	// a large index, and they differ precisely in the part after it. Keying on
	// the whole suffix would reject this, which is a worse failure than the one
	// being prevented - it would refuse indexes that work.
	const std::filesystem::path dir = FreshDir("onecfg");
	WriteIndexFiles(dir, "-m4r1-frac");
	WriteIndexFiles(dir, "-m4r1-frac2");

	REQUIRE_NOTHROW(ValidateIndexLayout(dir.string()));
	REQUIRE(ValidateIndexLayout(dir.string()).size() == 2);
	std::filesystem::remove_all(dir);
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
