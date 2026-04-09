#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "VsearchChimeraWrapper.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

// Simple FASTA parser for test data.
static void read_fasta(const std::string &path, std::vector<std::string> &labels, std::vector<std::string> &sequences) {
	std::ifstream in(path);
	REQUIRE(in.good());
	std::string line;
	std::string current_label;
	std::string current_seq;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		if (line[0] == '>') {
			if (!current_label.empty()) {
				labels.push_back(current_label);
				sequences.push_back(current_seq);
			}
			current_label = line.substr(1);
			current_seq.clear();
		} else {
			current_seq += line;
		}
	}
	if (!current_label.empty()) {
		labels.push_back(current_label);
		sequences.push_back(current_seq);
	}
}

// Parse expected_ref.tsv: tab-separated, columns are:
// score, query_label, parent_a, parent_b, closest_parent,
// id_QM, id_QA, id_QB, id_AB, id_QT,
// left_yes, left_no, left_abstain, right_yes, right_no, right_abstain,
// divergence, flag
struct ExpectedRow {
	std::string query_label;
	std::string flag;
	double score;
	std::string parent_a;
	std::string parent_b;
};

static std::unordered_map<std::string, ExpectedRow> read_expected(const std::string &path) {
	std::ifstream in(path);
	REQUIRE(in.good());
	std::unordered_map<std::string, ExpectedRow> rows;
	std::string line;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		std::istringstream ss(line);
		std::string tok;
		ExpectedRow row;
		std::getline(ss, tok, '\t');
		row.score = std::stod(tok);
		std::getline(ss, row.query_label, '\t');
		std::getline(ss, row.parent_a, '\t');
		std::getline(ss, row.parent_b, '\t');
		// skip columns 5-17 (closest_parent through divergence)
		for (int i = 0; i < 13; i++) {
			std::getline(ss, tok, '\t');
		}
		// column 18: flag (last column, no trailing tab)
		std::getline(ss, row.flag, '\t');
		// Strip trailing whitespace (newline)
		while (!row.flag.empty() && (row.flag.back() == '\n' || row.flag.back() == '\r')) {
			row.flag.pop_back();
		}
		rows[row.query_label] = row;
	}
	return rows;
}

TEST_CASE("VsearchChimeraWrapper synthetic uchime_ref", "[VsearchChimera]") {
	std::vector<std::string> ref_labels, ref_seqs;
	read_fasta("data/uchime/chimera_ref.fasta", ref_labels, ref_seqs);
	REQUIRE(ref_labels.size() == 6);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/uchime/chimera_queries.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 8);

	auto expected = read_expected("data/uchime/expected_ref.tsv");
	REQUIRE(expected.size() == 8);

	miint::VsearchChimeraWrapper wrapper;
	wrapper.set_reference(ref_labels, ref_seqs);
	auto handle = wrapper.create_detect_handle();

	for (size_t i = 0; i < query_labels.size(); i++) {
		auto result = handle.detect(query_labels[i], query_seqs[i]);

		auto it = expected.find(query_labels[i]);
		REQUIRE(it != expected.end());
		auto &exp = it->second;

		INFO("Query: " << query_labels[i]);
		REQUIRE(result.flag == exp.flag);

		if (exp.flag == "Y") {
			// For chimeras, verify parents match
			CHECK(result.parent_a_label == exp.parent_a);
			CHECK(result.parent_b_label == exp.parent_b);
			// Score should be close (exact match expected since we ARE vsearch)
			CHECK_THAT(result.score, Catch::Matchers::WithinRel(exp.score, 0.01));
		}
	}
}

TEST_CASE("VsearchChimeraWrapper detect_batch lifecycle", "[VsearchChimera]") {
	// Regression test: detect_batch manages its own chimera session lifecycle
	// internally (chimera_detect_batch calls chimera_session_init/cleanup).
	// Previously, set_reference() eagerly initialized the chimera session,
	// causing double mutex destroy on wrapper destruction.
	std::vector<std::string> ref_labels, ref_seqs;
	read_fasta("data/uchime/chimera_ref.fasta", ref_labels, ref_seqs);
	REQUIRE(ref_labels.size() == 6);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/uchime/chimera_queries.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 8);

	miint::VsearchChimeraWrapper wrapper;
	wrapper.set_reference(ref_labels, ref_seqs);

	// First batch
	std::vector<miint::UchimeResult> results1;
	wrapper.detect_batch(query_labels, query_seqs, results1);
	REQUIRE(results1.size() == 8);

	// Second batch — exercises repeated session init/cleanup
	std::vector<miint::UchimeResult> results2;
	wrapper.detect_batch(query_labels, query_seqs, results2);
	REQUIRE(results2.size() == 8);

	// Results should be identical
	for (size_t i = 0; i < results1.size(); i++) {
		CHECK(results1[i].flag == results2[i].flag);
		CHECK(results1[i].query_label == results2[i].query_label);
	}

	// Wrapper destruction should not crash (no double mutex destroy)
}

TEST_CASE("VsearchChimeraWrapper real-world LTP chimera", "[VsearchChimera]") {
	std::vector<std::string> ref_labels, ref_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", ref_labels, ref_seqs);
	REQUIRE(ref_labels.size() == 500);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/uchime/EU778281_chimera.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 1);

	miint::VsearchChimeraWrapper wrapper;
	wrapper.set_reference(ref_labels, ref_seqs);
	auto handle = wrapper.create_detect_handle();

	auto result = handle.detect(query_labels[0], query_seqs[0]);
	// EU778281 is a published chimera, but its actual parents (MK567960/JAOQKB000000000)
	// are NOT in the 500-ref subset. vsearch correctly classifies it as N here.
	// This test verifies the wrapper produces the same result as native vsearch.
	REQUIRE(result.flag == "N");
}
