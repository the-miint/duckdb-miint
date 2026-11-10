#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "GtdbFiles.hpp"

#include <optional>
#include <stdexcept>
#include <string>
#include <vector>

using namespace miint::gtdb;
using Catch::Matchers::ContainsSubstring;

namespace {

Release R(const char *major, const char *minor) {
	return Release {major, minor};
}

} // namespace

// ---------------------------------------------------------------------------
// Releases
// ---------------------------------------------------------------------------

// A bare major number is the release's first sub-release, so `release := '232'` and
// '232.0' name the same files and the same cache directory.
TEST_CASE("ParseRelease gives a bare major number its .0", "[gtdb]") {
	CHECK(ParseRelease("232").ToString() == "232.0");
	CHECK(ParseRelease("232.0").ToString() == "232.0");
	CHECK(ParseRelease("214.1").ToString() == "214.1");
	CHECK(ParseRelease("95").ToString() == "95.0");
	const auto r = ParseRelease("214.1");
	CHECK(r.major == "214");
	CHECK(r.minor == "1");
}

// The release becomes part of a URL and a cache path, so anything but GTDB's numbers (no
// leading zero) and one dot is refused rather than passed on: '0232' would be looked for
// where GTDB serves nothing.
TEST_CASE("ParseRelease refuses anything but GTDB release numbers", "[gtdb]") {
	for (const char *bad : {"", "r232", "release232", "232.0.1", "-1", "+232", "232.", ".1", " 232", "232 ", "1e3",
	                        "23 2", "232/../x", "232.a", "0232", "232.01", "00"}) {
		INFO("release: '" << bad << "'");
		CHECK_THROWS_AS(ParseRelease(bad), std::invalid_argument);
		CHECK_THROWS_WITH(ParseRelease(bad), ContainsSubstring("'" + std::string(bad) + "'") &&
		                                         ContainsSubstring("232") && ContainsSubstring("214.1"));
	}
}

// ---------------------------------------------------------------------------
// Where a release's files are served and cached
// ---------------------------------------------------------------------------

// The paths GTDB serves, as probed on data.gtdb.ecogenomic.org: the major number names
// the release directory and the files, and qc_failed sits under auxillary_files/.
TEST_CASE("ServedPath follows GTDB's layout", "[gtdb]") {
	CHECK(ServedPath(ReleaseFile::BAC120_METADATA, R("232", "0")) == "release232/232.0/bac120_metadata_r232.tsv.gz");
	CHECK(ServedPath(ReleaseFile::AR53_METADATA, R("214", "1")) == "release214/214.1/ar53_metadata_r214.tsv.gz");
	CHECK(ServedPath(ReleaseFile::QC_FAILED, R("232", "0")) == "release232/232.0/auxillary_files/qc_failed_r232.tsv");
	CHECK(ReleaseFileName(ReleaseFile::QC_FAILED, R("220", "0")) == "qc_failed_r220.tsv");
}

// Every file of a release sits in one directory named for the full release, so the
// directory can be passed as a source, and 214.0 and 214.1 never share one.
TEST_CASE("CachedPath keeps a release's files side by side", "[gtdb]") {
	CHECK(CachedPath(ReleaseFile::BAC120_METADATA, R("232", "0")) == "232.0/bac120_metadata_r232.tsv.gz");
	CHECK(CachedPath(ReleaseFile::QC_FAILED, R("232", "0")) == "232.0/qc_failed_r232.tsv");
	CHECK(CachedPath(ReleaseFile::AR53_METADATA, R("214", "1")) == "214.1/ar53_metadata_r214.tsv.gz");
}

// ---------------------------------------------------------------------------
// File names
// ---------------------------------------------------------------------------

// The domain comes from the file's name, so the name must be found in a path, a Windows
// path and a URL alike, of any scheme; a URL's query is not part of it.
TEST_CASE("FileName finds the name in paths and URLs", "[gtdb]") {
	CHECK(FileName("bac120_metadata_r232.tsv.gz") == "bac120_metadata_r232.tsv.gz");
	CHECK(FileName("/data/gtdb/ar53_metadata_r232.tsv") == "ar53_metadata_r232.tsv");
	CHECK(FileName("C:\\gtdb\\ar53_metadata_r232.tsv") == "ar53_metadata_r232.tsv");
	CHECK(FileName("https://h/releases/ar53_metadata_r232.tsv.gz?token=x/y") == "ar53_metadata_r232.tsv.gz");
	CHECK(FileName("http://h/ar53_metadata_r232.tsv#part") == "ar53_metadata_r232.tsv");
	CHECK(FileName("s3://bucket/gtdb/ar53_metadata_r232.tsv.gz?versionId=3") == "ar53_metadata_r232.tsv.gz");
	CHECK(FileName("HTTPS://h/bac120_metadata_r232.tsv.gz?x") == "bac120_metadata_r232.tsv.gz");
	// Only a URL has a query: a local name may hold '?' or '#'.
	CHECK(FileName("/data/odd?name#1.tsv") == "odd?name#1.tsv");
	CHECK(FileName("/data/dir/") == "");
}

// GTDB names its metadata files for the marker set: bac120 for bacteria, ar53 (ar122
// before release 207) for archaea. The match is exact, as GTDB writes the names.
TEST_CASE("MetadataDomain reads the domain from the file name", "[gtdb]") {
	CHECK(MetadataDomain("bac120_metadata_r232.tsv.gz") == std::optional<std::string>("Bacteria"));
	CHECK(MetadataDomain("ar53_metadata_r232.tsv") == std::optional<std::string>("Archaea"));
	CHECK(MetadataDomain("ar122_metadata_r202.tsv") == std::optional<std::string>("Archaea"));
	CHECK(MetadataDomain("bac120_renamed.tsv") == std::optional<std::string>("Bacteria"));
	for (const char *other : {"qc_failed_r232.tsv", "Bac120_metadata_r232.tsv", "xbac120_metadata.tsv",
	                          "bac120metadata.tsv", "ar5_metadata.tsv", "ar53", "", "metadata.tsv"}) {
		INFO("name: '" << other << "'");
		CHECK_FALSE(MetadataDomain(other).has_value());
	}
}

// ---------------------------------------------------------------------------
// Directories
// ---------------------------------------------------------------------------

// A directory is read as a release: the cache's layout, or GTDB's files saved side by
// side. Other files are ignored, among them a download still in progress (.part).
TEST_CASE("MatchMetadataFiles picks a release's metadata files", "[gtdb]") {
	const std::vector<std::string> cache = {"qc_failed_r232.tsv", "ar53_metadata_r232.tsv.gz",
	                                        "bac120_metadata_r232.tsv.gz.4f2a.part", "README",
	                                        "bac120_metadata_r232.tsv.gz"};
	CHECK(MatchMetadataFiles(cache) ==
	      std::vector<std::string> {"bac120_metadata_r232.tsv.gz", "ar53_metadata_r232.tsv.gz"});
	CHECK(MatchMetadataFiles({"ar53_metadata_r232.tsv"}) == std::vector<std::string> {"ar53_metadata_r232.tsv"});
	CHECK(MatchMetadataFiles({"bac120_metadata_r95.tsv", "x.tsv"}) ==
	      std::vector<std::string> {"bac120_metadata_r95.tsv"});
	CHECK(MatchMetadataFiles({"ar122_metadata_r202.tsv", "bac120_metadata_r202.tsv"}) ==
	      std::vector<std::string> {"bac120_metadata_r202.tsv", "ar122_metadata_r202.tsv"});
}

// A directory that does not name one release's files is refused, not guessed at: two
// bacterial files would read the same genomes twice, and two releases would mix
// taxonomies.
TEST_CASE("MatchMetadataFiles refuses a directory that is not one release", "[gtdb]") {
	using V = std::vector<std::string>;
	CHECK_THROWS_WITH(MatchMetadataFiles(V {}), ContainsSubstring("bac120_metadata_r"));
	CHECK_THROWS_WITH(MatchMetadataFiles(V {"qc_failed_r232.tsv", "bac120_metadata.tsv", "bac120_metadata_r.tsv",
	                                        "bac120_metadata_rx.tsv", "bac120_metadata_r232.csv",
	                                        "bac120_taxonomy_r232.tsv", "bac120_metadata_r232.tsv.gz.part"}),
	                  ContainsSubstring("bac120_metadata_r") && ContainsSubstring("ar53_metadata_r"));
	CHECK_THROWS_WITH(MatchMetadataFiles(V {"bac120_metadata_r226.tsv.gz", "bac120_metadata_r232.tsv.gz"}),
	                  ContainsSubstring("bac120_metadata_r226.tsv.gz") &&
	                      ContainsSubstring("bac120_metadata_r232.tsv.gz"));
	CHECK_THROWS_WITH(MatchMetadataFiles(V {"bac120_metadata_r232.tsv", "bac120_metadata_r232.tsv.gz"}),
	                  ContainsSubstring("bac120_metadata_r232.tsv") &&
	                      ContainsSubstring("bac120_metadata_r232.tsv.gz"));
	CHECK_THROWS_WITH(MatchMetadataFiles(V {"ar53_metadata_r232.tsv.gz", "ar122_metadata_r202.tsv"}),
	                  ContainsSubstring("ar53_metadata_r232.tsv.gz") && ContainsSubstring("ar122_metadata_r202.tsv"));
	CHECK_THROWS_WITH(MatchMetadataFiles(V {"bac120_metadata_r226.tsv.gz", "ar53_metadata_r232.tsv.gz"}),
	                  ContainsSubstring("release") && ContainsSubstring("bac120_metadata_r226.tsv.gz") &&
	                      ContainsSubstring("ar53_metadata_r232.tsv.gz"));
}

TEST_CASE("MatchQcFailedFile wants exactly one qc_failed file", "[gtdb]") {
	using V = std::vector<std::string>;
	CHECK(MatchQcFailedFile(V {"bac120_metadata_r232.tsv.gz", "qc_failed_r232.tsv"}) == "qc_failed_r232.tsv");
	CHECK(MatchQcFailedFile(V {"qc_failed_r232.tsv.gz", "qc_failed_r232.tsv.gz.1.part"}) == "qc_failed_r232.tsv.gz");
	CHECK_THROWS_WITH(MatchQcFailedFile(V {"bac120_metadata_r232.tsv.gz", "qc_failed.tsv"}),
	                  ContainsSubstring("qc_failed_r"));
	CHECK_THROWS_WITH(MatchQcFailedFile(V {"qc_failed_r226.tsv", "qc_failed_r232.tsv"}),
	                  ContainsSubstring("qc_failed_r226.tsv") && ContainsSubstring("qc_failed_r232.tsv"));
}

// ---------------------------------------------------------------------------
// Header lines
// ---------------------------------------------------------------------------

// The header names the columns as GTDB writes them: qc_failed's hold spaces and "(%)".
TEST_CASE("SplitHeader keeps the names as written", "[gtdb]") {
	CHECK(SplitHeader("accession\tambiguous_bases\tgtdb_taxonomy") ==
	      std::vector<std::string> {"accession", "ambiguous_bases", "gtdb_taxonomy"});
	CHECK(SplitHeader("Accession\tNCBI species\tCM1 completeness (%)") ==
	      std::vector<std::string> {"Accession", "NCBI species", "CM1 completeness (%)"});
	CHECK(SplitHeader("one") == std::vector<std::string> {"one"});
}

// A nameless column cannot be selected, and DuckDB treats names that differ only in case
// as the same column, so either is refused at the header instead of failing later.
TEST_CASE("SplitHeader refuses empty and repeated names", "[gtdb]") {
	CHECK_THROWS_WITH(SplitHeader(""), ContainsSubstring("column 1"));
	CHECK_THROWS_WITH(SplitHeader("a\t\tb"), ContainsSubstring("column 2"));
	CHECK_THROWS_WITH(SplitHeader("a\tb\t"), ContainsSubstring("column 3"));
	CHECK_THROWS_WITH(SplitHeader("a\tb\ta"), ContainsSubstring("'a'") && ContainsSubstring("column 3"));
	CHECK_THROWS_WITH(SplitHeader("Accession\taccession"),
	                  ContainsSubstring("'accession'") && ContainsSubstring("column 2"));
}

// The header is split on tabs, not read as CSV, while the rows are read as CSV, where a
// '"' quotes. GTDB never quotes a name, so a name holding one is refused rather than read
// one way here and another by the row reader.
TEST_CASE("SplitHeader refuses names holding a double quote", "[gtdb]") {
	CHECK_THROWS_WITH(SplitHeader("accession\t\"x\""),
	                  ContainsSubstring("column 2") && ContainsSubstring("double quote"));
	CHECK_THROWS_WITH(SplitHeader("say \"hi\""), ContainsSubstring("column 1") && ContainsSubstring("double quote"));
}
