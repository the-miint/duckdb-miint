#include <catch2/catch_test_macros.hpp>
#include "aspera_utils.hpp"

TEST_CASE("AsperaUtils path parsing", "[aspera]") {
	SECTION("Single path") {
		auto paths = miint::AsperaUtils::ParseAsperaPaths(
		    "fasp.sra.ebi.ac.uk:/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz");
		REQUIRE(paths.size() == 1);
		CHECK(paths[0].host == "fasp.sra.ebi.ac.uk");
		CHECK(paths[0].remote_path == "/vol1/fastq/ERR107/007/ERR1074767/ERR1074767.fastq.gz");
	}

	SECTION("Paired-end paths (semicolon-separated)") {
		auto paths = miint::AsperaUtils::ParseAsperaPaths(
		    "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR134/085/SRR13425885/SRR13425885_1.fastq.gz;"
		    "fasp.sra.ebi.ac.uk:/vol1/fastq/SRR134/085/SRR13425885/SRR13425885_2.fastq.gz");
		REQUIRE(paths.size() == 2);
		CHECK(paths[0].host == "fasp.sra.ebi.ac.uk");
		CHECK(paths[0].remote_path == "/vol1/fastq/SRR134/085/SRR13425885/SRR13425885_1.fastq.gz");
		CHECK(paths[1].host == "fasp.sra.ebi.ac.uk");
		CHECK(paths[1].remote_path == "/vol1/fastq/SRR134/085/SRR13425885/SRR13425885_2.fastq.gz");
	}

	SECTION("Empty string") {
		auto paths = miint::AsperaUtils::ParseAsperaPaths("");
		CHECK(paths.empty());
	}

	SECTION("Path without host prefix") {
		auto paths = miint::AsperaUtils::ParseAsperaPaths("vol1/fastq/ERR107/ERR1074767.fastq.gz");
		REQUIRE(paths.size() == 1);
		CHECK(paths[0].host == "fasp.sra.ebi.ac.uk");
		CHECK(paths[0].remote_path == "/vol1/fastq/ERR107/ERR1074767.fastq.gz");
	}

	SECTION("Multiple semicolons with empty segments") {
		auto paths = miint::AsperaUtils::ParseAsperaPaths("fasp.sra.ebi.ac.uk:/path1.gz;;fasp.sra.ebi.ac.uk:/path2.gz");
		REQUIRE(paths.size() == 2);
		CHECK(paths[0].remote_path == "/path1.gz");
		CHECK(paths[1].remote_path == "/path2.gz");
	}
}

TEST_CASE("AsperaUtils config building", "[aspera]") {
	auto config = miint::AsperaUtils::BuildConfig("/usr/bin/ascp", "/home/user/.aspera/key.openssh");
	CHECK(config.ascp_path == "/usr/bin/ascp");
	CHECK(config.key_path == "/home/user/.aspera/key.openssh");
	CHECK(config.host == "fasp.sra.ebi.ac.uk");
	CHECK(config.user == "era-fasp");
	CHECK(config.port == 33001);
	CHECK(config.max_rate == "300m");
}
