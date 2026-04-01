#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "VsearchClusterWrapper.hpp"

#include <fstream>
#include <sstream>
#include <string>
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

// Parse vsearch UC file.  Returns map: query_label -> {cluster_id, is_centroid, centroid_label, identity}.
struct UCRecord {
	int cluster_id;
	bool is_centroid;
	std::string centroid_label;
	double identity;
};

static std::unordered_map<std::string, UCRecord> read_uc(const std::string &path) {
	std::ifstream in(path);
	REQUIRE(in.good());
	std::unordered_map<std::string, UCRecord> records;
	std::string line;
	while (std::getline(in, line)) {
		if (line.empty()) {
			continue;
		}
		// Only parse S (centroid) and H (hit) records
		if (line[0] != 'S' && line[0] != 'H') {
			continue;
		}
		std::istringstream ss(line);
		std::string type, cluster_str, len_str, id_str, strand, col5, col6, cigar, query_label, target_label;
		std::getline(ss, type, '\t');
		std::getline(ss, cluster_str, '\t');
		std::getline(ss, len_str, '\t');
		std::getline(ss, id_str, '\t');
		std::getline(ss, strand, '\t');
		std::getline(ss, col5, '\t');
		std::getline(ss, col6, '\t');
		std::getline(ss, cigar, '\t');
		std::getline(ss, query_label, '\t');
		std::getline(ss, target_label, '\t');

		UCRecord rec;
		rec.cluster_id = std::stoi(cluster_str);
		rec.is_centroid = (type == "S");
		rec.centroid_label = rec.is_centroid ? query_label : target_label;
		rec.identity = rec.is_centroid ? 100.0 : std::stod(id_str);
		records[query_label] = rec;
	}
	return records;
}

TEST_CASE("VsearchClusterWrapper identical sequences", "[VsearchCluster]") {
	// Three identical sequences should all end up in the same cluster.
	std::vector<std::string> labels = {"seq1", "seq2", "seq3"};
	std::vector<std::string> seqs = {
	    "ATCGATCGATCGATCGATCGATCG",
	    "ATCGATCGATCGATCGATCGATCG",
	    "ATCGATCGATCGATCGATCGATCG",
	};

	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	auto results = wrapper.cluster_all();

	REQUIRE(results.size() == 3);

	// First sequence is the centroid
	CHECK(results[0].is_centroid == true);
	CHECK(results[0].cluster_id == 0);
	CHECK(results[0].identity == 100.0);

	// Others join the same cluster
	CHECK(results[1].is_centroid == false);
	CHECK(results[1].cluster_id == 0);
	CHECK(results[1].identity == 100.0);
	CHECK(results[1].centroid_id == "seq1");

	CHECK(results[2].is_centroid == false);
	CHECK(results[2].cluster_id == 0);
	CHECK(results[2].centroid_id == "seq1");
}

TEST_CASE("VsearchClusterWrapper distinct sequences", "[VsearchCluster]") {
	// Three very different sequences should each be their own centroid at 97%.
	std::vector<std::string> labels = {"seqA", "seqB", "seqC"};
	std::vector<std::string> seqs = {
	    "ATCGATCGATCGATCGATCGATCG",
	    "GGCCTTAAGGCCTTAAGGCCTTAA",
	    "AAAAAAAACCCCCCCCGGGGGGGG",
	};

	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	auto results = wrapper.cluster_all();

	REQUIRE(results.size() == 3);
	CHECK(results[0].is_centroid == true);
	CHECK(results[0].cluster_id == 0);
	CHECK(results[1].is_centroid == true);
	CHECK(results[1].cluster_id == 1);
	CHECK(results[2].is_centroid == true);
	CHECK(results[2].cluster_id == 2);
}

TEST_CASE("VsearchClusterWrapper empty input throws", "[VsearchCluster]") {
	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	CHECK_THROWS_AS(wrapper.set_sequences({}, {}), std::runtime_error);
}

TEST_CASE("VsearchClusterWrapper label too long throws", "[VsearchCluster]") {
	std::string long_label(1024, 'X');
	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	CHECK_THROWS_AS(wrapper.set_sequences({long_label}, {"ATCGATCG"}), std::runtime_error);
}

TEST_CASE("VsearchClusterWrapper cluster_all before set_sequences throws", "[VsearchCluster]") {
	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	CHECK_THROWS_AS(wrapper.cluster_all(), std::runtime_error);
}

TEST_CASE("VsearchClusterWrapper RNA normalization", "[VsearchCluster]") {
	// RNA (U) should be treated as DNA (T)
	std::vector<std::string> labels = {"dna", "rna"};
	std::vector<std::string> seqs = {
	    "ATCGATCGATCGATCGATCGATCG",
	    "AUCGAUCGAUCGAUCGAUCGAUCG",
	};

	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	auto results = wrapper.cluster_all();

	REQUIRE(results.size() == 2);
	CHECK(results[0].is_centroid == true);
	CHECK(results[1].is_centroid == false);
	CHECK(results[1].cluster_id == 0);
	CHECK(results[1].identity == 100.0);
}

TEST_CASE("VsearchClusterWrapper reinitialize", "[VsearchCluster]") {
	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);

	// First run
	wrapper.set_sequences({"s1", "s2"}, {"ATCGATCGATCGATCGATCGATCG", "ATCGATCGATCGATCGATCGATCG"});
	auto r1 = wrapper.cluster_all();
	REQUIRE(r1.size() == 2);

	// Second run with different data — should work cleanly
	wrapper.set_sequences({"a", "b", "c"},
	                      {"GGCCTTAAGGCCTTAAGGCCTTAA", "AAAAAAAACCCCCCCCGGGGGGGG", "TTTTTTTTTTTTTTTTTTTTTTTT"});
	auto r2 = wrapper.cluster_all();
	REQUIRE(r2.size() == 3);
	CHECK(r2[0].is_centroid == true);
	CHECK(r2[1].is_centroid == true);
	CHECK(r2[2].is_centroid == true);
}

TEST_CASE("VsearchClusterWrapper strand both vs plus", "[VsearchCluster]") {
	// A sequence and its reverse complement should be 2 clusters with plus-only
	// but 1 cluster with both-strand search.
	std::vector<std::string> labels = {"fwd", "rc_of_fwd"};
	std::vector<std::string> seqs = {
	    "AAAAAAAAAAAAAAAAAAAAAAACCCCCCCCCCCCCCCCCCCCCCC",
	    "GGGGGGGGGGGGGGGGGGGGGGGTTTTTTTTTTTTTTTTTTTTTTTT",
	};

	// Plus-only: no RC match, 2 clusters
	{
		miint::ClusterParams params;
		params.id = 0.90;
		params.strand = 1;
		miint::VsearchClusterWrapper wrapper(params);
		wrapper.set_sequences(labels, seqs);
		auto results = wrapper.cluster_all();
		REQUIRE(results.size() == 2);
		CHECK(results[0].is_centroid == true);
		CHECK(results[1].is_centroid == true);
	}

	// Both strands: RC match found, 1 cluster
	{
		miint::ClusterParams params;
		params.id = 0.90;
		params.strand = 2;
		miint::VsearchClusterWrapper wrapper(params);
		wrapper.set_sequences(labels, seqs);
		auto results = wrapper.cluster_all();
		REQUIRE(results.size() == 2);
		CHECK(results[0].is_centroid == true);
		CHECK(results[1].is_centroid == false);
		CHECK(results[1].cluster_id == 0);
	}
}

TEST_CASE("VsearchClusterWrapper real-world LTP 97%", "[VsearchCluster]") {
	// Ground truth: vsearch --cluster_fast data/cluster/ltp_200.fasta --id 0.97
	//   --uc data/cluster/expected_97.uc
	// Run date: 2026-03-31, vsearch v2.30.5

	std::vector<std::string> labels, seqs;
	read_fasta("data/cluster/ltp_200_sorted.fasta", labels, seqs);
	REQUIRE(labels.size() == 200);

	auto expected = read_uc("data/cluster/expected_97.uc");
	REQUIRE(expected.size() == 200);

	miint::ClusterParams params;
	params.id = 0.97;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	auto results = wrapper.cluster_all();

	REQUIRE(results.size() == 200);

	// Count centroids
	int centroid_count = 0;
	for (auto &r : results) {
		if (r.is_centroid) {
			centroid_count++;
		}
	}
	CHECK(centroid_count == 24);

	// Verify each sequence's cluster assignment matches vsearch ground truth.
	for (size_t i = 0; i < results.size(); i++) {
		auto it = expected.find(results[i].read_id);
		INFO("Sequence: " << results[i].read_id);
		REQUIRE(it != expected.end());
		auto &exp = it->second;

		CHECK(results[i].is_centroid == exp.is_centroid);
		CHECK(results[i].cluster_id == exp.cluster_id);
		// vsearch CLI rounds identity to 1 decimal place
		CHECK_THAT(results[i].identity, Catch::Matchers::WithinAbs(exp.identity, 0.1));
	}
}

TEST_CASE("VsearchClusterWrapper real-world LTP 99%", "[VsearchCluster]") {
	// Ground truth: vsearch --cluster_fast data/cluster/ltp_200.fasta --id 0.99
	//   --uc data/cluster/expected_99.uc
	// Run date: 2026-03-31, vsearch v2.30.5

	std::vector<std::string> labels, seqs;
	read_fasta("data/cluster/ltp_200_sorted.fasta", labels, seqs);
	REQUIRE(labels.size() == 200);

	auto expected = read_uc("data/cluster/expected_99.uc");
	REQUIRE(expected.size() == 200);

	miint::ClusterParams params;
	params.id = 0.99;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	auto results = wrapper.cluster_all();

	REQUIRE(results.size() == 200);

	int centroid_count = 0;
	for (auto &r : results) {
		if (r.is_centroid) {
			centroid_count++;
		}
	}
	CHECK(centroid_count == 114);

	for (size_t i = 0; i < results.size(); i++) {
		auto it = expected.find(results[i].read_id);
		INFO("Sequence: " << results[i].read_id);
		REQUIRE(it != expected.end());
		auto &exp = it->second;

		CHECK(results[i].is_centroid == exp.is_centroid);
		CHECK(results[i].cluster_id == exp.cluster_id);
		CHECK_THAT(results[i].identity, Catch::Matchers::WithinAbs(exp.identity, 0.1));
	}
}
