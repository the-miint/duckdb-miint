#include <catch2/catch_test_macros.hpp>

#include "VsearchChimeraWrapper.hpp"
#include "VsearchClusterWrapper.hpp"
#include "VsearchSearchWrapper.hpp"
#include "vsearch_serial.hpp"

#include <algorithm>
#include <cctype>
#include <fstream>
#include <string>
#include <tuple>
#include <vector>

// The serial stand-ins in vsearch_serial.hpp are the only vsearch code paths a
// thread-less wasm build can run. What makes them a substitute rather than an
// approximation is that they return exactly what the threaded entry points they
// replace return, so each test runs the same real-world LTP input through both
// and requires identical results, field for field. Each threaded run is also
// anchored to the value its SQL test commits to, so the two paths cannot agree
// by both being empty.

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

template <class T>
static std::vector<T> slice(const std::vector<T> &v, size_t begin, size_t end) {
	return std::vector<T>(v.begin() + begin, v.begin() + end);
}

// ---------------------------------------------------------------------------
// DustAllSerial
// ---------------------------------------------------------------------------

static std::vector<std::string> mask_database(const std::vector<std::string> &labels,
                                              const std::vector<std::string> &seqs, void (*mask)()) {
	vsearch_init_defaults();
	vsearch_apply_defaults_fixups();
	db_init();
	for (size_t i = 0; i < labels.size(); i++) {
		db_add(false, labels[i].c_str(), seqs[i].c_str(), nullptr, labels[i].size(), seqs[i].size(), 1);
	}
	mask();
	std::vector<std::string> masked;
	for (uint64_t i = 0; i < db_getsequencecount(); i++) {
		masked.emplace_back(db_getsequence(i), db_getsequencelen(i));
	}
	db_free();
	vsearch_session_end();
	return masked;
}

TEST_CASE("DustAllSerial masks the database exactly as dust_all does", "[VsearchSerial]") {
	std::vector<std::string> labels, seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", labels, seqs);
	REQUIRE(labels.size() == 500);
	// A 16S sequence with a poly-A run spliced in: DUST must soft-mask that run,
	// which is what lets this test fail if a stand-in skipped or misapplied masking
	// (on clean 16S alone, "masked nothing" and "masked correctly" look the same).
	labels.emplace_back("low_complexity");
	seqs.push_back(seqs[0].substr(0, 200) + std::string(80, 'A') + seqs[1].substr(0, 200));

	const auto threaded = mask_database(labels, seqs, dust_all);
	const auto serial = mask_database(labels, seqs, miint::DustAllSerial);

	REQUIRE(serial.size() == labels.size());
	REQUIRE(serial == threaded);
	const auto &spliced = serial.back();
	REQUIRE(std::any_of(spliced.begin(), spliced.end(), [](char c) { return std::islower(c) != 0; }));
}

// ---------------------------------------------------------------------------
// search_batch
// ---------------------------------------------------------------------------

static std::vector<miint::SearchResult> search_ltp(bool serial) {
	std::vector<std::string> ref_labels, ref_seqs, q_labels, q_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", ref_labels, ref_seqs);
	read_fasta("data/search/ltp_query_50.fasta", q_labels, q_seqs);
	REQUIRE(q_labels.size() == 50);

	miint::SearchParams params;
	params.id = 0.90;
	params.maxaccepts = 3;
	params.threads = 4;
	params.serial = serial;
	miint::VsearchSearchWrapper wrapper(params);
	wrapper.set_database(ref_labels, ref_seqs);

	// Two calls, the way search_sequences_vsearch feeds the wrapper chunk by chunk.
	std::vector<miint::SearchResult> out;
	wrapper.search_batch(slice(q_labels, 0, 25), slice(q_seqs, 0, 25), out);
	wrapper.search_batch(slice(q_labels, 25, 50), slice(q_seqs, 25, 50), out);
	return out;
}

TEST_CASE("serial search_batch returns what the threaded search_batch returns", "[VsearchSerial]") {
	const auto threaded = search_ltp(false);
	const auto serial = search_ltp(true);

	// search_sequences_realworld.test: id 0.90 with maxaccepts 3 yields 150 hits.
	REQUIRE(threaded.size() == 150);
	REQUIRE(serial.size() == threaded.size());
	for (size_t i = 0; i < threaded.size(); i++) {
		INFO("hit " << i << ": " << threaded[i].query_id << " -> " << threaded[i].target_id);
		CHECK(serial[i].query_id == threaded[i].query_id);
		CHECK(serial[i].target_id == threaded[i].target_id);
		CHECK(serial[i].identity == threaded[i].identity);
		CHECK(serial[i].matches == threaded[i].matches);
		CHECK(serial[i].mismatches == threaded[i].mismatches);
		CHECK(serial[i].gaps == threaded[i].gaps);
		CHECK(serial[i].alignment_length == threaded[i].alignment_length);
		CHECK(serial[i].query_length == threaded[i].query_length);
		CHECK(serial[i].target_length == threaded[i].target_length);
		CHECK(serial[i].accepted == threaded[i].accepted);
	}
}

// ---------------------------------------------------------------------------
// chimera_detect_batch (uchime_ref)
// ---------------------------------------------------------------------------

// The option globals chimera_session_init() overwrites. chimera_detect_batch()
// puts them back afterwards, so the stand-in must too. Nothing in the wrapper
// reads them in between (each batch re-initializes them to fixed values), so
// reading them back directly is the only way a skipped restore can show.
static auto chimera_globals() {
	return std::make_tuple(opt_maxaccepts, opt_maxrejects, opt_id, opt_weak_id, opt_self, opt_selfid, opt_maxsizeratio);
}
using ChimeraGlobals = decltype(chimera_globals());

static std::vector<miint::UchimeResult> chimera_ltp(bool serial, ChimeraGlobals &globals_after) {
	std::vector<std::string> ref_labels, ref_seqs, q_labels, q_seqs;
	read_fasta("data/uchime/ltp_subset_500.fasta", ref_labels, ref_seqs);
	read_fasta("data/uchime/ltp_chimera_queries.fasta", q_labels, q_seqs);
	REQUIRE(q_labels.size() == 9);

	miint::UchimeParams params;
	params.threads = 4;
	params.serial = serial;
	miint::VsearchChimeraWrapper wrapper(params);
	wrapper.set_reference(ref_labels, ref_seqs);

	// Two calls, the way detect_chimera_uchime feeds the wrapper chunk by chunk, so
	// each call's session setup and teardown run more than once.
	std::vector<miint::UchimeResult> out;
	wrapper.detect_batch(slice(q_labels, 0, 5), slice(q_seqs, 0, 5), out);
	wrapper.detect_batch(slice(q_labels, 5, 9), slice(q_seqs, 5, 9), out);
	globals_after = chimera_globals(); // while this wrapper's vsearch session is still live
	return out;
}

TEST_CASE("serial detect_batch returns what the threaded detect_batch returns", "[VsearchSerial]") {
	ChimeraGlobals threaded_globals, serial_globals;
	const auto threaded = chimera_ltp(false, threaded_globals);
	const auto serial = chimera_ltp(true, serial_globals);
	REQUIRE(serial_globals == threaded_globals);

	// uchime_ref_realworld.test: the five chimera_* queries are flagged Y.
	REQUIRE(threaded.size() == 9);
	REQUIRE(std::count_if(threaded.begin(), threaded.end(),
	                      [](const miint::UchimeResult &r) { return r.flag == "Y"; }) == 5);
	REQUIRE(serial.size() == threaded.size());
	for (size_t i = 0; i < threaded.size(); i++) {
		INFO("query " << threaded[i].query_label);
		CHECK(serial[i].query_label == threaded[i].query_label);
		CHECK(serial[i].flag == threaded[i].flag);
		CHECK(serial[i].score == threaded[i].score);
		CHECK(serial[i].parent_a_label == threaded[i].parent_a_label);
		CHECK(serial[i].parent_b_label == threaded[i].parent_b_label);
		CHECK(serial[i].closest_parent_label == threaded[i].closest_parent_label);
		CHECK(serial[i].id_query_model == threaded[i].id_query_model);
		CHECK(serial[i].id_query_a == threaded[i].id_query_a);
		CHECK(serial[i].id_query_b == threaded[i].id_query_b);
		CHECK(serial[i].id_a_b == threaded[i].id_a_b);
		CHECK(serial[i].id_query_top == threaded[i].id_query_top);
		CHECK(serial[i].left_yes == threaded[i].left_yes);
		CHECK(serial[i].left_no == threaded[i].left_no);
		CHECK(serial[i].left_abstain == threaded[i].left_abstain);
		CHECK(serial[i].right_yes == threaded[i].right_yes);
		CHECK(serial[i].right_no == threaded[i].right_no);
		CHECK(serial[i].right_abstain == threaded[i].right_abstain);
		CHECK(serial[i].divergence == threaded[i].divergence);
	}
}

// ---------------------------------------------------------------------------
// cluster_assign_batch
// ---------------------------------------------------------------------------

static std::vector<miint::ClusterResult> cluster_ltp(bool serial) {
	std::vector<std::string> labels, seqs;
	read_fasta("data/cluster/ltp_200_sorted.fasta", labels, seqs);
	REQUIRE(labels.size() == 200);

	miint::ClusterParams params;
	params.id = 0.97;
	params.threads = 4;
	params.serial = serial;
	miint::VsearchClusterWrapper wrapper(params);
	wrapper.set_sequences(labels, seqs);
	return wrapper.cluster_all();
}

TEST_CASE("serial cluster_all returns what the threaded cluster_all returns", "[VsearchSerial]") {
	const auto threaded = cluster_ltp(false);
	const auto serial = cluster_ltp(true);

	// cluster_sequences_realworld.test / vsearch --cluster_fast: 24 centroids at 97%.
	REQUIRE(threaded.size() == 200);
	REQUIRE(std::count_if(threaded.begin(), threaded.end(),
	                      [](const miint::ClusterResult &r) { return r.is_centroid; }) == 24);
	REQUIRE(serial.size() == threaded.size());
	for (size_t i = 0; i < threaded.size(); i++) {
		INFO("sequence " << threaded[i].read_id);
		CHECK(serial[i].read_id == threaded[i].read_id);
		CHECK(serial[i].is_centroid == threaded[i].is_centroid);
		CHECK(serial[i].cluster_id == threaded[i].cluster_id);
		CHECK(serial[i].centroid_id == threaded[i].centroid_id);
		CHECK(serial[i].identity == threaded[i].identity);
		CHECK(serial[i].cigar == threaded[i].cigar);
		CHECK(serial[i].cigar_truncated == threaded[i].cigar_truncated);
	}
}
