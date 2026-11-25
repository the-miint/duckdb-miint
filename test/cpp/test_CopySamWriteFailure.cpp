// COPY FORMAT SAM/BAM/UBAM when a write fails, as on a full disk. A failed close must fail the COPY and leave no file
// (why htslib reports a failure there: docs/internals/embedded-tools.md, htslib): a writer that dropped sam_close's
// result would report success over a truncated file.
//
// The failure is made with a file-size limit (RLIMIT_FSIZE, with SIGXFSZ ignored, so a write past it fails with
// EFBIG), as test_CopyBiomWriteFailure.cpp does. Set here, in the process, the limit can sit one byte under the size
// the file needs, which no SQL test can arrange.
//
// POSIX only (CMakeLists.txt).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <chrono>
#include <csignal>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <string>
#include <sys/resource.h>

namespace {

using namespace duckdb;
using Catch::Matchers::ContainsSubstring;

// Lowers the soft file-size limit to `bytes` and ignores SIGXFSZ until destroyed, then restores both.
class FileSizeLimit {
public:
	explicit FileSizeLimit(rlim_t bytes) {
		REQUIRE(getrlimit(RLIMIT_FSIZE, &saved_limit) == 0);
		auto lowered = saved_limit;
		lowered.rlim_cur = bytes;
		saved_handler = std::signal(SIGXFSZ, SIG_IGN);
		REQUIRE(setrlimit(RLIMIT_FSIZE, &lowered) == 0);
	}
	~FileSizeLimit() {
		setrlimit(RLIMIT_FSIZE, &saved_limit);
		std::signal(SIGXFSZ, saved_handler);
	}

private:
	struct rlimit saved_limit {};
	void (*saved_handler)(int) = SIG_DFL;
};

// One output format: the table it writes, the COPY options and the name its errors give.
struct Output {
	std::string name;      // as error messages name the format
	std::string extension; // of the target file
	std::string table;     // holding the columns the format writes
	std::string options;
};

const Output OUTPUTS[] = {
    {"BAM", ".bam", "aln", "FORMAT BAM, REFERENCE_LENGTHS 'ref_lens'"},
    {"SAM", ".sam.gz", "aln", "FORMAT SAM, COMPRESSION gzip, REFERENCE_LENGTHS 'ref_lens'"},
    {"SAM", ".sam", "aln", "FORMAT SAM, REFERENCE_LENGTHS 'ref_lens'"},
    {"UBAM", ".bam", "reads", "FORMAT UBAM"},
};

struct SamWriteFixture {
	DuckDB db;
	Connection con;
	std::filesystem::path dir;

	SamWriteFixture() : db(nullptr), con(db) {
		dir = std::filesystem::temp_directory_path() /
		      ("miint_copy_sam_write_failure_" +
		       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		for (auto sql :
		     {"CREATE TABLE ref_lens AS SELECT 'ref' || k AS name, 1000000::BIGINT AS length FROM range(4) t(k)",
		      "CREATE TABLE aln AS SELECT 'r' || i AS read_id, (CASE WHEN i % 2 = 0 THEN 16 ELSE 0 END)::USMALLINT "
		      "AS flags, 'ref' || (i % 4) AS reference, (1 + i % 900000)::BIGINT AS position, 60::UTINYINT AS "
		      "mapq, '50M' AS cigar, '*' AS mate_reference, 0::BIGINT AS mate_position, 0::BIGINT AS "
		      "template_length, (i % 100)::BIGINT AS tag_as FROM range(20000) t(i)",
		      // Reads that compress little, so that the uBAM, like the others, is many times 64 KiB.
		      "CREATE TABLE reads AS SELECT 'r' || i AS read_id, left(translate(md5(i::VARCHAR), "
		      "'0123456789abcdef', 'ACGTACGTACGTACGT') || translate(md5((i + 1)::VARCHAR), '0123456789abcdef', "
		      "'ACGTACGTACGTACGT'), 50) AS sequence1, list_transform(range(50), lambda j: ((i * 7 + j * 13) % "
		      "41)::UTINYINT) AS qual1 FROM range(20000) t(i)"}) {
			auto created = con.Query(sql);
			INFO(created->GetError());
			REQUIRE_FALSE(created->HasError());
		}
	}

	~SamWriteFixture() {
		std::error_code ignored;
		std::filesystem::remove_all(dir, ignored);
	}

	void SetThreads(int threads) {
		auto set = con.Query("SET threads = " + std::to_string(threads));
		REQUIRE_FALSE(set->HasError());
	}

	// The limit, if any, is lifted before this returns, so Catch's own output is never limited.
	unique_ptr<QueryResult> Copy(const std::filesystem::path &target, const Output &output,
	                             rlim_t limit = RLIM_INFINITY) {
		auto sql = "COPY " + output.table + " TO '" + target.string() + "' (" + output.options + ")";
		if (limit == RLIM_INFINITY) {
			return con.Query(sql);
		}
		FileSizeLimit lowered(limit);
		return con.Query(sql);
	}

	// The message names the format and the file; DuckDB removes a file its failed COPY created.
	void RequireWriteRefused(const unique_ptr<QueryResult> &result, const std::filesystem::path &target,
	                         const Output &output) {
		REQUIRE(result->HasError());
		REQUIRE_THAT(result->GetError(),
		             ContainsSubstring("COPY FORMAT " + output.name + ": cannot write '" + target.string() + "'"));
		REQUIRE_FALSE(std::filesystem::exists(target));
	}
};

std::string ReadFile(const std::filesystem::path &path) {
	std::ifstream in(path, std::ios::binary);
	return std::string(std::istreambuf_iterator<char>(in), std::istreambuf_iterator<char>());
}

} // namespace

TEST_CASE("COPY FORMAT SAM/BAM/UBAM reports a write that fails when the file closes", "[copy_sam]") {
	SamWriteFixture fixture;
	// One thread, and several: with several, BGZF blocks are compressed and written on htslib's own threads, which
	// hold a failed write until the file is closed.
	for (int threads : {1, 4}) {
		fixture.SetThreads(threads);
		for (const auto &output : OUTPUTS) {
			auto stem = output.name + "_t" + std::to_string(threads);
			auto full = fixture.dir / ("full_" + stem + output.extension);
			auto written = fixture.Copy(full, output);
			INFO(written->GetError());
			REQUIRE_FALSE(written->HasError());
			auto size = std::filesystem::file_size(full);

			// One byte short: everything fits but the last write, which only the close makes.
			INFO(output.options << " at threads=" << threads << ", file size " << size);
			auto target = fixture.dir / ("short_" + stem + output.extension);
			fixture.RequireWriteRefused(fixture.Copy(target, output, size - 1), target, output);
		}
	}
}

TEST_CASE("COPY FORMAT SAM/BAM/UBAM reports a write that fails partway", "[copy_sam]") {
	SamWriteFixture fixture;
	for (int threads : {1, 4}) {
		fixture.SetThreads(threads);
		for (const auto &output : OUTPUTS) {
			// 64 KiB holds the header and a little of the records.
			auto target = fixture.dir / ("partway_" + output.name + "_t" + std::to_string(threads) + output.extension);
			INFO(output.options << " at threads=" << threads);
			fixture.RequireWriteRefused(fixture.Copy(target, output, 64 * 1024), target, output);
		}
	}

	// htslib is left usable: the same COPY with no limit writes a file that reads back whole.
	auto again = fixture.dir / "again.bam";
	REQUIRE_FALSE(fixture.Copy(again, OUTPUTS[0])->HasError());
	auto read = fixture.con.Query("SELECT count(*) FROM read_alignments('" + again.string() + "')");
	INFO(read->GetError());
	REQUIRE_FALSE(read->HasError());
	auto chunk = read->Fetch();
	REQUIRE(chunk);
	REQUIRE(chunk->GetValue(0, 0).GetValue<int64_t>() == 20000);
}

TEST_CASE("COPY FORMAT SAM/BAM/UBAM over an existing file: a failed write leaves the file as it was", "[copy_sam]") {
	SamWriteFixture fixture;
	for (const auto &output : OUTPUTS) {
		INFO(output.options);
		auto target = fixture.dir / ("existing_" + output.name + output.extension);
		auto written = fixture.Copy(target, output);
		INFO(written->GetError());
		REQUIRE_FALSE(written->HasError());
		const auto before = ReadFile(target);

		// DuckDB writes over an existing file through tmp_<name>, so the error names that file
		// (docs/internals/duckdb-engine-notes.md, "A failed COPY leaves no output").
		auto tmp = fixture.dir / ("tmp_existing_" + output.name + output.extension);
		auto result = fixture.Copy(target, output, 64 * 1024);
		REQUIRE(result->HasError());
		REQUIRE_THAT(result->GetError(),
		             ContainsSubstring("COPY FORMAT " + output.name + ": cannot write '" + tmp.string() + "'"));
		REQUIRE(ReadFile(target) == before);
		REQUIRE_FALSE(std::filesystem::exists(tmp));
	}
}
