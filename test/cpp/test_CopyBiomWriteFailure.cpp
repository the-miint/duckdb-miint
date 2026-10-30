// COPY FORMAT BIOM when its write fails partway, as on a full disk. HDF5 writes a compressed dataset's chunks
// when the dataset closes, and the metadata it holds when the file closes, so a failure in either close must
// reach the COPY. Left to the HDF5 objects' destructors, it was dropped: the COPY reported success, left a
// file read_biom could not open, and the process crashed at exit.
//
// The failure is made with a file-size limit (RLIMIT_FSIZE, with SIGXFSZ ignored, so a write past it fails
// with EFBIG). Set here, in the process, the limit can sit one byte under the size the file needs, which no
// SQL test can arrange.
//
// POSIX only; built only with HDF5 (CMakeLists.txt).

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

#include <chrono>
#include <csignal>
#include <filesystem>
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

struct BiomWriteFixture {
	DuckDB db;
	Connection con;
	std::filesystem::path dir;

	BiomWriteFixture() : db(nullptr), con(db) {
		dir = std::filesystem::temp_directory_path() /
		      ("miint_copy_biom_write_failure_" +
		       std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()));
		std::filesystem::create_directories(dir);
		auto created = con.Query("CREATE TABLE biom_rows AS SELECT 'F' || (r % 50) AS feature_id, 'S' || (r // 50) "
		                         "AS sample_id, (r + 1)::DOUBLE AS value FROM range(2000) t(r)");
		REQUIRE_FALSE(created->HasError());
	}

	~BiomWriteFixture() {
		std::error_code ignored;
		std::filesystem::remove_all(dir, ignored);
	}

	// The limit, if any, is lifted before this returns, so Catch's own output is never limited.
	unique_ptr<QueryResult> Copy(const std::filesystem::path &target, const std::string &compression,
	                             rlim_t limit = RLIM_INFINITY) {
		auto sql = "COPY biom_rows TO '" + target.string() + "' (FORMAT BIOM, COMPRESSION '" + compression + "')";
		if (limit == RLIM_INFINITY) {
			return con.Query(sql);
		}
		FileSizeLimit lowered(limit);
		return con.Query(sql);
	}

	// The message names the file, the HDF5 call that failed, and the OS error
	void RequireWriteRefused(const unique_ptr<QueryResult> &result, const std::filesystem::path &target,
	                         const std::string &failed_call) {
		REQUIRE(result->HasError());
		REQUIRE_THAT(result->GetError(), ContainsSubstring("COPY FORMAT BIOM: cannot write '" + target.string() +
		                                                   "': " + failed_call + " failed: File too large"));
		// DuckDB removes a file its failed COPY created
		REQUIRE_FALSE(std::filesystem::exists(target));
	}
};

} // namespace

TEST_CASE("COPY FORMAT BIOM reports a write that fails when a dataset closes", "[biom][copy]") {
	BiomWriteFixture fixture;
	for (std::string compression : {"gzip", "none"}) {
		// HDF5's overhead alone makes a BIOM file larger than 16 KiB, and this table's data is held until each
		// dataset closes
		auto target = fixture.dir / ("dataset_close_" + compression + ".biom");
		fixture.RequireWriteRefused(fixture.Copy(target, compression, 16 * 1024), target, "H5Dclose");
	}
}

TEST_CASE("COPY FORMAT BIOM reports a write that fails when the file closes", "[biom][copy]") {
	BiomWriteFixture fixture;
	auto full = fixture.dir / "full.biom";
	REQUIRE_FALSE(fixture.Copy(full, "gzip")->HasError());
	auto size = std::filesystem::file_size(full);

	// One byte short. With gzip, every dataset's data fits, and the write that fails is the last one, of the
	// metadata the file's close writes. ('none' puts its data at the end of the file, so there the write that
	// fails is a dataset's close, the case above.)
	auto target = fixture.dir / "file_close.biom";
	fixture.RequireWriteRefused(fixture.Copy(target, "gzip", size - 1), target, "H5Fclose");

	// HDF5 is left usable: the same COPY with no limit writes a file that reads back whole
	auto again = fixture.dir / "again.biom";
	REQUIRE_FALSE(fixture.Copy(again, "gzip")->HasError());
	auto read = fixture.con.Query("SELECT count(*) FROM read_biom('" + again.string() + "')");
	REQUIRE_FALSE(read->HasError());
	auto chunk = read->Fetch();
	REQUIRE(chunk);
	REQUIRE(chunk->GetValue(0, 0).GetValue<int64_t>() == 2000);
}
