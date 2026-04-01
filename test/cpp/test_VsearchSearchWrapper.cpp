#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "VsearchSearchWrapper.hpp"

#include <fstream>
#include <sstream>
#include <string>
#include <thread>
#include <unordered_map>
#include <unordered_set>
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

// Parse expected search TSV: query target id alnlen mism gaps ql tl
struct ExpectedHit {
	std::string query;
	std::string target;
	double identity;
	int alnlen;
	int mismatches;
	int gaps;
	int query_length;
	int target_length;
};

static std::vector<ExpectedHit> read_expected_search(const std::string &path) {
	std::ifstream in(path);
	REQUIRE(in.good());
	std::vector<ExpectedHit> rows;
	std::string line;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		std::istringstream ss(line);
		ExpectedHit h;
		std::string tok;
		std::getline(ss, h.query, '\t');
		std::getline(ss, h.target, '\t');
		std::getline(ss, tok, '\t');
		h.identity = std::stod(tok);
		std::getline(ss, tok, '\t');
		h.alnlen = std::stoi(tok);
		std::getline(ss, tok, '\t');
		h.mismatches = std::stoi(tok);
		std::getline(ss, tok, '\t');
		h.gaps = std::stoi(tok);
		std::getline(ss, tok, '\t');
		h.query_length = std::stoi(tok);
		std::getline(ss, tok, '\t');
		h.target_length = std::stoi(tok);
		rows.push_back(h);
	}
	return rows;
}

TEST_CASE("VsearchSearchWrapper synthetic exact match", "[VsearchSearch]") {
	// DB with 3 distinct sequences; query is identical to one of them.
	std::vector<std::string> db_labels = {"seqA", "seqB", "seqC"};
	std::vector<std::string> db_seqs = {
	    "ATCGATCGATCGATCGATCGATCG",
	    "GGCCTTAAGGCCTTAAGGCCTTAA",
	    "AAAAAAAACCCCCCCCGGGGGGGG",
	};

	miint::SearchParams params;
	params.id = 0.97;
	params.maxaccepts = 1;

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	auto results = handle.search("queryA", "ATCGATCGATCGATCGATCGATCG");
	REQUIRE(results.size() == 1);
	CHECK(results[0].target_id == "seqA");
	CHECK(results[0].identity == 100.0);
	CHECK(results[0].accepted == true);
	CHECK(results[0].matches == 24);
	CHECK(results[0].mismatches == 0);
	CHECK(results[0].gaps == 0);
}

TEST_CASE("VsearchSearchWrapper no match below threshold", "[VsearchSearch]") {
	std::vector<std::string> db_labels = {"seqA"};
	std::vector<std::string> db_seqs = {"ATCGATCGATCGATCGATCGATCG"};

	miint::SearchParams params;
	params.id = 0.99; // Very high threshold

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	// Query with many mismatches — should not pass 99% threshold
	auto results = handle.search("queryX", "TTTTTTTTTTTTTTTTTTTTTTTT");
	CHECK(results.empty());
}

TEST_CASE("VsearchSearchWrapper RNA normalization", "[VsearchSearch]") {
	// DB has DNA, query has RNA (U instead of T)
	std::vector<std::string> db_labels = {"dna_seq"};
	std::vector<std::string> db_seqs = {"ATCGATCGATCGATCGATCGATCG"};

	miint::SearchParams params;
	params.id = 0.97;

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	auto results = handle.search("rna_query", "AUCGAUCGAUCGAUCGAUCGAUCG");
	REQUIRE(results.size() == 1);
	CHECK(results[0].identity == 100.0);
}

TEST_CASE("VsearchSearchWrapper label too long", "[VsearchSearch]") {
	std::vector<std::string> db_labels = {"seqA"};
	std::vector<std::string> db_seqs = {"ATCGATCGATCGATCGATCGATCG"};

	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	std::string long_label(1024, 'X');
	CHECK_THROWS_AS(handle.search(long_label, "ATCGATCG"), std::runtime_error);
}

TEST_CASE("VsearchSearchWrapper empty DB throws", "[VsearchSearch]") {
	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);

	CHECK_THROWS_AS(wrapper.set_database({}, {}), std::runtime_error);
}

TEST_CASE("VsearchSearchWrapper db label too long throws", "[VsearchSearch]") {
	std::string long_label(1024, 'X');
	std::vector<std::string> db_labels = {long_label};
	std::vector<std::string> db_seqs = {"ATCGATCGATCGATCGATCGATCG"};

	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);

	CHECK_THROWS_AS(wrapper.set_database(db_labels, db_seqs), std::runtime_error);
}

TEST_CASE("VsearchSearchWrapper empty query sequence returns empty", "[VsearchSearch]") {
	std::vector<std::string> db_labels = {"seqA"};
	std::vector<std::string> db_seqs = {"ATCGATCGATCGATCGATCGATCG"};

	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	auto results = handle.search("emptyQ", "");
	CHECK(results.empty());
}

TEST_CASE("VsearchSearchWrapper uninitialized handle throws", "[VsearchSearch]") {
	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);

	// create_search_handle before set_database should throw
	CHECK_THROWS_AS(wrapper.create_search_handle(), std::runtime_error);
}

TEST_CASE("VsearchSearchWrapper reinitialize with set_database", "[VsearchSearch]") {
	// Calling set_database twice should work — second call replaces the first session.
	std::vector<std::string> db_labels1 = {"seqA"};
	std::vector<std::string> db_seqs1 = {"ATCGATCGATCGATCGATCGATCG"};

	// Use a completely different sequence family for DB2 so there's no cross-match at 97%.
	std::vector<std::string> db_labels2 = {"seqX"};
	std::vector<std::string> db_seqs2 = {"CCCCCCCCCCCCCCCCCCCCCCCC"};

	miint::SearchParams params;
	params.id = 0.97;
	miint::VsearchSearchWrapper wrapper(params);

	// First database — handle scoped so it's destroyed before re-init.
	{
		wrapper.set_database(db_labels1, db_seqs1);
		auto handle1 = wrapper.create_search_handle();
		auto r1 = handle1.search("q", "ATCGATCGATCGATCGATCGATCG");
		REQUIRE(r1.size() == 1);
		CHECK(r1[0].target_id == "seqA");
	}

	// Replace database — releases old session, acquires new one.
	wrapper.set_database(db_labels2, db_seqs2);
	auto handle2 = wrapper.create_search_handle();
	auto r2 = handle2.search("q", "CCCCCCCCCCCCCCCCCCCCCCCC");
	REQUIRE(r2.size() == 1);
	CHECK(r2[0].target_id == "seqX");

	// Old DB sequence should not match in new DB (ATCG vs all-C at 97%)
	auto r3 = handle2.search("q", "ATCGATCGATCGATCGATCGATCG");
	CHECK(r3.empty());
}

TEST_CASE("VsearchSearchWrapper parallel search", "[VsearchSearch]") {
	// Verify thread-safe concurrent searching with multiple handles.
	std::vector<std::string> db_labels, db_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", db_labels, db_seqs);
	REQUIRE(db_labels.size() == 500);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/search/ltp_query_50.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 50);

	miint::SearchParams params;
	params.id = 0.97;
	params.maxaccepts = 1;

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);

	constexpr int NUM_THREADS = 4;

	// Create handles on main thread — search_init sets file-statics
	// (seqcount, tophits) that are not thread-safe to set concurrently.
	std::vector<miint::VsearchSearchWrapper::SearchHandle> handles;
	for (int t = 0; t < NUM_THREADS; t++) {
		handles.push_back(wrapper.create_search_handle());
	}

	// Each thread searches a slice of queries and stores results.
	std::vector<std::vector<miint::SearchResult>> thread_results(NUM_THREADS);
	std::vector<std::thread> threads;

	for (int t = 0; t < NUM_THREADS; t++) {
		threads.emplace_back([&, t]() {
			size_t start = t * (query_labels.size() / NUM_THREADS);
			size_t end = (t == NUM_THREADS - 1) ? query_labels.size() : (t + 1) * (query_labels.size() / NUM_THREADS);

			for (size_t i = start; i < end; i++) {
				auto results = handles[t].search(query_labels[i], query_seqs[i]);
				for (auto &r : results) {
					thread_results[t].push_back(std::move(r));
				}
			}
		});
	}

	for (auto &th : threads) {
		th.join();
	}

	// Collect all results
	int total_hits = 0;
	for (auto &tr : thread_results) {
		total_hits += tr.size();
	}

	// All 50 queries should find a self-match at 97%
	REQUIRE(total_hits == 50);

	// Verify all results are self-matches at 100%
	for (auto &tr : thread_results) {
		for (auto &r : tr) {
			CHECK(r.query_id == r.target_id);
			CHECK(r.identity == 100.0);
		}
	}
}

TEST_CASE("VsearchSearchWrapper real-world LTP 97%", "[VsearchSearch]") {
	// Ground truth: vsearch --usearch_global data/search/ltp_query_50.fasta
	//   --db data/uchime/ltp_subset_500.fasta --id 0.97
	//   --userout data/search/expected_97.tsv
	//   --userfields query+target+id+alnlen+mism+gaps+ql+tl
	// Run date: 2026-03-30, vsearch v2.30.5

	std::vector<std::string> db_labels, db_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", db_labels, db_seqs);
	REQUIRE(db_labels.size() == 500);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/search/ltp_query_50.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 50);

	auto expected = read_expected_search("data/search/expected_97.tsv");
	REQUIRE(expected.size() == 50);

	// Build expected lookup: query -> expected hit
	std::unordered_map<std::string, ExpectedHit> expected_map;
	for (auto &e : expected) {
		expected_map[e.query] = e;
	}

	miint::SearchParams params;
	params.id = 0.97;
	params.maxaccepts = 1;
	params.maxrejects = 32;

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	int matched = 0;
	for (size_t i = 0; i < query_labels.size(); i++) {
		auto results = handle.search(query_labels[i], query_seqs[i]);
		auto it = expected_map.find(query_labels[i]);
		REQUIRE(it != expected_map.end());
		auto &exp = it->second;

		INFO("Query: " << query_labels[i]);
		REQUIRE(results.size() == 1);
		CHECK(results[0].target_id == exp.target);
		CHECK_THAT(results[0].identity, Catch::Matchers::WithinAbs(exp.identity, 0.01));
		CHECK(results[0].alignment_length == exp.alnlen);
		CHECK(results[0].mismatches == exp.mismatches);
		CHECK(results[0].gaps == exp.gaps);
		matched++;
	}
	REQUIRE(matched == 50);
}

TEST_CASE("VsearchSearchWrapper real-world LTP 90% maxaccepts=3", "[VsearchSearch]") {
	// Ground truth: vsearch --usearch_global data/search/ltp_query_50.fasta
	//   --db data/uchime/ltp_subset_500.fasta --id 0.90 --maxaccepts 3
	//   --userout data/search/expected_90_max3.tsv
	//   --userfields query+target+id+alnlen+mism+gaps+ql+tl
	// Run date: 2026-03-30, vsearch v2.30.5

	std::vector<std::string> db_labels, db_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", db_labels, db_seqs);
	REQUIRE(db_labels.size() == 500);

	std::vector<std::string> query_labels, query_seqs;
	read_fasta("data/search/ltp_query_50.fasta", query_labels, query_seqs);
	REQUIRE(query_labels.size() == 50);

	auto expected = read_expected_search("data/search/expected_90_max3.tsv");
	REQUIRE(expected.size() == 150);

	// Build expected lookup: query -> list of hits (ordered)
	std::unordered_map<std::string, std::vector<ExpectedHit>> expected_map;
	for (auto &e : expected) {
		expected_map[e.query].push_back(e);
	}

	miint::SearchParams params;
	params.id = 0.90;
	params.maxaccepts = 3;
	params.maxrejects = 32;

	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(db_labels, db_seqs);
	auto handle = wrapper.create_search_handle();

	int total_hits = 0;
	for (size_t i = 0; i < query_labels.size(); i++) {
		auto results = handle.search(query_labels[i], query_seqs[i]);
		auto it = expected_map.find(query_labels[i]);
		REQUIRE(it != expected_map.end());
		auto &exp_hits = it->second;

		INFO("Query: " << query_labels[i]);
		REQUIRE(results.size() == exp_hits.size());

		for (size_t j = 0; j < results.size(); j++) {
			INFO("  Hit " << j << ": " << results[j].target_id);
			CHECK(results[j].target_id == exp_hits[j].target);
			// vsearch CLI --userfields id rounds to 1 decimal place, so
			// expected values are only 0.1% precision. Tolerance ±0.1.
			CHECK_THAT(results[j].identity, Catch::Matchers::WithinAbs(exp_hits[j].identity, 0.1));
		}
		total_hits += results.size();
	}
	REQUIRE(total_hits == 150);
}
