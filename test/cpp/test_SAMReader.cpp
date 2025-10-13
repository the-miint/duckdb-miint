#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <string>
#include <fstream>
#include <filesystem>
#include <iostream>
#include "SAMReader.hpp"

std::string header_reference_line(const std::string &id, const std::string &length) {
	return "@SQ\tSN:" + id + "\tLN:" + length + "\n";
}
// Utilities for test data
std::string unpaired_record(const std::string &id, const std::string &ref, const std::string &flags,
                            const std::string &cigar, const std::string &fields = "") {
	return id + "\t" + flags + "\t" + ref + "\t" + "1\t0\t" + cigar + "\t*\t0\t0\tATGC\tIIII\n";
}

// TODO: move to utilities
void write_temp_file2(const std::string &path, const std::vector<std::string> &records) {
	std::ofstream out(path);
	for (const auto &r : records) {
		out << r;
	}
}

TEST_CASE("SAMReader with header unpaired no flags", "[SAMReader]") {
	auto path = "test.sam";
	write_temp_file2(path, {header_reference_line("g1", "1000"), unpaired_record("r1", "g1", "0", "4M"),
	                        unpaired_record("r2", "g1", "256", "4M")});

	miint::SAMReader reader(path);
	auto batch = reader.read(5);

	REQUIRE((batch.size() == 2));
	REQUIRE((batch[0].read_id == "r1"));
	REQUIRE((batch[0].flags == 0));
	REQUIRE((batch[0].reference == "g1"));
	REQUIRE((batch[0].position == 1));
	REQUIRE((batch[0].mapq == 0));
	REQUIRE((batch[0].cigar == "4M"));
	REQUIRE((batch[0].mate_reference == "*"));
	REQUIRE((batch[0].mate_position == 0));
	REQUIRE((batch[0].template_length == 0));
	REQUIRE((batch[0].tag_as == -1));
	REQUIRE((batch[0].tag_xs == -1));
	REQUIRE((batch[0].tag_ys == -1));
	REQUIRE((batch[0].tag_xn == -1));
	REQUIRE((batch[0].tag_xm == -1));
	REQUIRE((batch[0].tag_xo == -1));
	REQUIRE((batch[0].tag_xg == -1));
	REQUIRE((batch[0].tag_nm == -1));
	REQUIRE((batch[0].tag_yt == ""));
	REQUIRE((batch[0].tag_md == ""));
	REQUIRE((batch[0].tag_sa == ""));

	REQUIRE((batch[1].read_id == "r2"));
	REQUIRE((batch[1].flags == 256));
	REQUIRE((batch[1].reference == "g1"));
	REQUIRE((batch[1].position == 1));
	REQUIRE((batch[1].mapq == 0));
	REQUIRE((batch[1].cigar == "4M"));
	REQUIRE((batch[1].mate_reference == "*"));
	REQUIRE((batch[1].mate_position == 0));
	REQUIRE((batch[1].template_length == 0));
	REQUIRE((batch[1].tag_as == -1));
	REQUIRE((batch[1].tag_xs == -1));
	REQUIRE((batch[1].tag_ys == -1));
	REQUIRE((batch[1].tag_xn == -1));
	REQUIRE((batch[1].tag_xm == -1));
	REQUIRE((batch[1].tag_xo == -1));
	REQUIRE((batch[1].tag_xg == -1));
	REQUIRE((batch[1].tag_nm == -1));
	REQUIRE((batch[1].tag_yt == ""));
	REQUIRE((batch[1].tag_md == ""));
	REQUIRE((batch[1].tag_sa == ""));
	std::filesystem::remove(path);
}
