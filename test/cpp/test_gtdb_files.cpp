#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "GtdbFiles.hpp"
#include "taxdump_archive.hpp"

#include <zlib.h>

#include <cstring>
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

// gzip framing of `text`, as GTDB's .tsv.gz files are.
std::string Gzip(const std::string &text) {
	z_stream zs;
	std::memset(&zs, 0, sizeof(zs));
	REQUIRE(deflateInit2(&zs, Z_DEFAULT_COMPRESSION, Z_DEFLATED, 16 + MAX_WBITS, 8, Z_DEFAULT_STRATEGY) == Z_OK);
	std::string out(deflateBound(&zs, static_cast<uLong>(text.size())), '\0');
	zs.next_in = reinterpret_cast<Bytef *>(const_cast<char *>(text.data()));
	zs.avail_in = static_cast<uInt>(text.size());
	zs.next_out = reinterpret_cast<Bytef *>(&out[0]);
	zs.avail_out = static_cast<uInt>(out.size());
	REQUIRE(deflate(&zs, Z_FINISH) == Z_STREAM_END);
	out.resize(zs.total_out);
	deflateEnd(&zs);
	return out;
}

// The first line the slow way, as the oracle for HeaderLineReader: the whole input,
// gunzipped whole by the taxdump reader's separate code, then cut at the first newline.
std::optional<std::string> NaiveFirstLine(const std::string &bytes) {
	std::string text = bytes.compare(0, 2, "\x1f\x8b") == 0 ? miint::TaxdumpArchive::Gunzip(bytes) : bytes;
	const auto newline = text.find('\n');
	if (newline == std::string::npos) {
		return std::nullopt;
	}
	std::string line = text.substr(0, newline);
	if (!line.empty() && line.back() == '\r') {
		line.pop_back();
	}
	return line;
}

// Feed `bytes` to a reader `chunk` bytes at a time, as a download hands them over, until
// it has the line; then end the input. Returns the line, and how many bytes were fed.
struct Read {
	std::string line;
	size_t fed;
	bool gzip;
};
Read ReadFirstLine(const std::string &bytes, size_t chunk, size_t max_line = HeaderLineReader::DEFAULT_MAX_LINE) {
	HeaderLineReader reader(max_line);
	size_t fed = 0;
	while (fed < bytes.size()) {
		const size_t n = std::min(chunk, bytes.size() - fed);
		const bool done = reader.Feed(bytes.data() + fed, n);
		fed += n;
		if (done) {
			break;
		}
	}
	auto line = reader.Finish();
	return Read {line, fed, reader.IsGzip()};
}

std::string Repeat(const std::string &s, size_t n) {
	std::string out;
	for (size_t i = 0; i < n; ++i) {
		out += s;
	}
	return out;
}

// A header shaped like GTDB's metadata header, and rows after it.
std::string MetadataLike(size_t rows) {
	std::string header = "accession";
	for (int i = 2; i <= 113; ++i) {
		header += "\tcolumn_" + std::to_string(i);
	}
	std::string text = header + "\n";
	for (size_t r = 0; r < rows; ++r) {
		text += "RS_GCF_" + std::to_string(r);
		for (int i = 2; i <= 113; ++i) {
			text += "\t" + std::to_string((r * 7919 + i * 104729) % 1000003);
		}
		text += "\n";
	}
	return text;
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
TEST_CASE("ParseRelease refuses anything but digits and one dot", "[gtdb]") {
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
TEST_CASE("ReleaseFileUrl follows GTDB's layout", "[gtdb]") {
	const std::string base = "https://data.gtdb.ecogenomic.org/releases";
	CHECK(ReleaseFileUrl(base, ReleaseFile::BAC120_METADATA, R("232", "0")) ==
	      base + "/release232/232.0/bac120_metadata_r232.tsv.gz");
	CHECK(ReleaseFileUrl(base, ReleaseFile::AR53_METADATA, R("214", "1")) ==
	      base + "/release214/214.1/ar53_metadata_r214.tsv.gz");
	CHECK(ReleaseFileUrl(base, ReleaseFile::QC_FAILED, R("232", "0")) ==
	      base + "/release232/232.0/auxillary_files/qc_failed_r232.tsv");
	CHECK(ReleaseFileName(ReleaseFile::QC_FAILED, R("220", "0")) == "qc_failed_r220.tsv");
}

// Every file of a release sits in one directory named for the full release, so the
// directory can be passed as a source, and 214.0 and 214.1 never share one.
TEST_CASE("ReleaseCachePath keeps a release's files side by side", "[gtdb]") {
	CHECK(ReleaseCachePath("/c/gtdb", ReleaseFile::BAC120_METADATA, R("232", "0")) ==
	      "/c/gtdb/232.0/bac120_metadata_r232.tsv.gz");
	CHECK(ReleaseCachePath("/c/gtdb", ReleaseFile::QC_FAILED, R("232", "0")) == "/c/gtdb/232.0/qc_failed_r232.tsv");
	CHECK(ReleaseCachePath("/c/gtdb/", ReleaseFile::AR53_METADATA, R("214", "1")) ==
	      "/c/gtdb/214.1/ar53_metadata_r214.tsv.gz");
}

// ---------------------------------------------------------------------------
// File names
// ---------------------------------------------------------------------------

// The domain comes from the file's name, so the name must be found in a path, a Windows
// path and a URL alike; a URL's query is not part of it.
TEST_CASE("FileName finds the name in paths and URLs", "[gtdb]") {
	CHECK(FileName("bac120_metadata_r232.tsv.gz") == "bac120_metadata_r232.tsv.gz");
	CHECK(FileName("/data/gtdb/ar53_metadata_r232.tsv") == "ar53_metadata_r232.tsv");
	CHECK(FileName("C:\\gtdb\\ar53_metadata_r232.tsv") == "ar53_metadata_r232.tsv");
	CHECK(FileName("https://h/releases/ar53_metadata_r232.tsv.gz?token=x/y") == "ar53_metadata_r232.tsv.gz");
	CHECK(FileName("http://h/ar53_metadata_r232.tsv#part") == "ar53_metadata_r232.tsv");
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

// The reader must agree with reading the whole input and cutting at the first newline,
// however the bytes arrive: all at once, a byte at a time (gzip's magic split across
// calls), or in odd-sized pieces.
TEST_CASE("HeaderLineReader agrees with reading the whole input", "[gtdb]") {
	const std::vector<std::string> texts = {
	    "a\tb\nrow1\trow2\n",
	    "a\tb\r\nrow1\trow2\r\n",
	    "\nnext\n",
	    "\r\n",
	    "only\n",
	    "x\ry\tz\n",
	    // gzip's first magic byte alone does not make a file gzip.
	    "\x1f\tb\nrow\n",
	    MetadataLike(3),
	    Repeat("w", 100000) + "\n" + "row\n",
	};
	for (const auto &text : texts) {
		for (const bool gz : {false, true}) {
			const auto bytes = gz ? Gzip(text) : text;
			const auto expected = NaiveFirstLine(bytes);
			REQUIRE(expected.has_value());
			for (const size_t chunk : {size_t(1), size_t(2), size_t(7), size_t(4096), bytes.size()}) {
				INFO("gzip " << gz << ", chunk " << chunk << ", text starts '" << text.substr(0, 20) << "'");
				const auto read = ReadFirstLine(bytes, chunk);
				CHECK(read.line == *expected);
				CHECK(read.gzip == gz);
			}
		}
	}
}

// The point of reading the header this way: it stops at the line, so binding a query
// does not read or download the file. A plain file stops at the newline's byte; a gzip
// stops within the compressed bytes that hold it.
TEST_CASE("HeaderLineReader stops at the first line", "[gtdb]") {
	const auto text = MetadataLike(20000);
	const auto plain = ReadFirstLine(text, 1);
	CHECK(plain.fed == text.find('\n') + 1);

	const auto gz = Gzip(text);
	REQUIRE(gz.size() > 1000000);
	const auto read = ReadFirstLine(gz, 1024);
	CHECK(read.fed < 65536);
	CHECK(read.line == text.substr(0, text.find('\n')));
}

// A download stopped once the header arrived is a gzip stream cut short, and that is
// fine; one cut before the header ends is refused, as is a file with no newline at all.
TEST_CASE("HeaderLineReader refuses input that ends before the line", "[gtdb]") {
	const auto text = MetadataLike(100);
	const auto gz = Gzip(text);
	CHECK(ReadFirstLine(gz.substr(0, gz.size() / 2), 64).line == text.substr(0, text.find('\n')));

	CHECK_THROWS_WITH(ReadFirstLine(Gzip(text).substr(0, 40), 7),
	                  ContainsSubstring("first line") && ContainsSubstring("cut short"));
	CHECK_THROWS_WITH(ReadFirstLine("a\tb", 1), ContainsSubstring("first line"));
	CHECK_THROWS_WITH(ReadFirstLine("", 1), ContainsSubstring("first line"));
	CHECK_THROWS_WITH(ReadFirstLine("\x1f", 1), ContainsSubstring("first line"));
}

// A gzip file may hold several members one after another, as bgzip writes it, and
// DuckDB's CSV reader reads the rows across them: so a header line that spans members
// is read whole, however the bytes arrive.
TEST_CASE("HeaderLineReader reads a line across gzip members", "[gtdb]") {
	const auto bytes = Gzip("a\tb") + Gzip("\tc") + Gzip("\nrow\n");
	for (const size_t chunk : {size_t(1), size_t(7), bytes.size()}) {
		INFO("chunk " << chunk);
		const auto read = ReadFirstLine(bytes, chunk);
		CHECK(read.line == "a\tb\tc");
		CHECK(read.gzip);
	}
}

// A whole gzip stream with no newline in it is refused, saying it ended rather than was
// cut short; and what follows a member must be another member.
TEST_CASE("HeaderLineReader refuses a gzip stream that ends before the line", "[gtdb]") {
	CHECK_THROWS_WITH(ReadFirstLine(Gzip("a\tb"), 1),
	                  ContainsSubstring("first line") && ContainsSubstring("ends before"));
	CHECK_THROWS_WITH(ReadFirstLine(Gzip("a\tb") + Gzip("\tc"), 5),
	                  ContainsSubstring("first line") && ContainsSubstring("ends before"));
	CHECK_THROWS_WITH(ReadFirstLine(Gzip("a\tb") + Repeat("\n", 1000), 64), ContainsSubstring("corrupt gzip"));
}

TEST_CASE("HeaderLineReader refuses a corrupt gzip stream", "[gtdb]") {
	auto gz = Gzip(MetadataLike(10));
	gz[12] = static_cast<char>(gz[12] ^ 0x5a);
	gz[13] = static_cast<char>(gz[13] ^ 0xa5);
	CHECK_THROWS_WITH(ReadFirstLine(gz, 4096), ContainsSubstring("gzip"));
	CHECK_THROWS_WITH(ReadFirstLine(std::string("\x1f\x8b", 2) + Repeat("z", 100), 4096), ContainsSubstring("gzip"));
}

// A file whose first line never ends (a binary file, or the wrong file) must not be read
// whole into memory: past the cap the reader refuses, and says so while bytes still
// arrive.
TEST_CASE("HeaderLineReader caps the line's length", "[gtdb]") {
	for (const bool gz : {false, true}) {
		INFO("gzip " << gz);
		const auto at_cap = Repeat("0123456789", 10) + "\nrow\n";
		const auto past_cap = Repeat("0123456789", 10) + "x\nrow\n";
		const auto wrap = [gz](const std::string &s) {
			return gz ? Gzip(s) : s;
		};
		// A byte at a time, the '\r' arrives before the '\n' that makes it the line's end.
		for (const size_t chunk : {size_t(1), size_t(3)}) {
			INFO("chunk " << chunk);
			CHECK(ReadFirstLine(wrap(at_cap), chunk, 100).line == Repeat("0123456789", 10));
			CHECK(ReadFirstLine(wrap(Repeat("0123456789", 10) + "\r\n"), chunk, 100).line == Repeat("0123456789", 10));
			CHECK_THROWS_WITH(ReadFirstLine(wrap(past_cap), chunk, 100), ContainsSubstring("100 bytes"));
		}

		const auto endless = wrap(Repeat("0123456789", 100000));
		HeaderLineReader reader(100);
		size_t fed = 0;
		bool threw = false;
		try {
			while (fed < endless.size()) {
				const size_t n = std::min<size_t>(64, endless.size() - fed);
				reader.Feed(endless.data() + fed, n);
				fed += n;
			}
		} catch (const std::runtime_error &e) {
			threw = true;
			CHECK_THAT(e.what(), ContainsSubstring("100 bytes"));
		}
		CHECK(threw);
		CHECK(fed < endless.size());
	}
}

// A retried download starts the body over, so the caller replaces its reader with a
// fresh one; the moved-to reader starts from nothing.
TEST_CASE("HeaderLineReader starts over when replaced", "[gtdb]") {
	const auto gz = Gzip("first\tline\nrow\n");
	HeaderLineReader reader;
	reader.Feed(gz.data(), 5);
	reader = HeaderLineReader();
	REQUIRE(reader.Feed("plain\n", 6));
	CHECK(reader.Finish() == "plain");
	CHECK_FALSE(reader.IsGzip());
}
