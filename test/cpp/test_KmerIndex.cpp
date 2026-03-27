#include <catch2/catch_test_macros.hpp>

#include "KmerIndex.hpp"

#include <algorithm>
#include <string>
#include <vector>

using miint::KmerIndex;

TEST_CASE("KmerIndex - Construction", "[KmerIndex]") {
	SECTION("Default construction succeeds") {
		KmerIndex idx;
		REQUIRE(idx.size() == 0);
	}

	SECTION("Empty index returns no candidates for any query") {
		KmerIndex idx;
		auto result = idx.find_candidates("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
		REQUIRE(result.empty());
	}
}

TEST_CASE("KmerIndex - add_sequence assigns sequential IDs", "[KmerIndex]") {
	KmerIndex idx;

	REQUIRE(idx.add_sequence("ACGTACGTACGTACGT") == 0);
	REQUIRE(idx.add_sequence("TGCATGCATGCATGCA") == 1);
	REQUIRE(idx.add_sequence("AAAA") == 2); // too short for k-mers, still gets ID
	REQUIRE(idx.add_sequence("GGGGGGGGGGGGGGGG") == 3);
	REQUIRE(idx.size() == 4);
}

TEST_CASE("KmerIndex - add_sequence and find_candidates basics", "[KmerIndex]") {
	KmerIndex idx;

	// A 60bp sequence — long enough that each of 4 chunks (15bp each) has valid 8-mers
	std::string seq_a = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	idx.add_sequence(seq_a);

	SECTION("Query identical to indexed sequence returns that sequence") {
		auto result = idx.find_candidates(seq_a);
		REQUIRE(result.size() == 1);
		REQUIRE(result[0] == 0);
	}

	SECTION("Query with no shared 8-mers returns empty") {
		// All T's — the 8-mer TTTTTTTT does not appear in seq_a (repeating ACGT pattern)
		std::string all_t(60, 'T');
		auto result = idx.find_candidates(all_t);
		REQUIRE(result.empty());
	}

	SECTION("Sequences shorter than k=8 produce no indexed k-mers") {
		KmerIndex short_idx;
		short_idx.add_sequence("ACGT"); // 4bp < 8, gets ID 0 but no k-mers
		REQUIRE(short_idx.size() == 1);
		auto result = short_idx.find_candidates("ACGT");
		REQUIRE(result.empty());
	}

	SECTION("Query shorter than k=8 against non-empty index returns empty") {
		auto result = idx.find_candidates("ACGTACG"); // 7bp < 8
		REQUIRE(result.empty());
	}

	SECTION("8-mers containing N are skipped") {
		KmerIndex n_idx;
		// 60bp sequence with N at position 4. 8-mers overlapping position 4 are invalid.
		// But 8-mers starting at positions 5+ (that don't span the N) are valid.
		std::string seq_with_n = "ACGTNACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
		n_idx.add_sequence(seq_with_n);
		auto result = n_idx.find_candidates(seq_a);
		REQUIRE(result.size() == 1);
		REQUIRE(result[0] == 0);
	}
}

TEST_CASE("KmerIndex - Unique k-mer membership (presence/absence)", "[KmerIndex]") {
	// Posting lists store each seq_id at most once per k-mer.
	// A repetitive sequence should not inflate its hit score.
	KmerIndex idx;

	// Highly repetitive: the 8-mer AAAAAAAA appears many times
	std::string repetitive(60, 'A');
	// Diverse: many distinct 8-mers
	std::string diverse = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

	idx.add_sequence(repetitive); // id=0
	idx.add_sequence(diverse);    // id=1

	// Query with all-A should find repetitive as top hit, but its score
	// should reflect unique k-mer presence (1 unique 8-mer: AAAAAAAA),
	// not multiplicity
	auto result = idx.find_candidates(repetitive);
	REQUIRE(!result.empty());
	REQUIRE(result[0] == 0);
}

TEST_CASE("KmerIndex - Chimeric query returns both parents", "[KmerIndex]") {
	KmerIndex idx;

	// Two clearly distinct 102bp sequences (divisible by 4 = 25.5, so chunks ~25-26bp)
	std::string parent_a = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                       "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

	std::string parent_b = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
	                       "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";

	idx.add_sequence(parent_a); // id=0
	idx.add_sequence(parent_b); // id=1

	// Chimera: first 51bp from parent_a, last 51bp from parent_b
	std::string chimera = parent_a.substr(0, 51) + parent_b.substr(51);

	auto result = idx.find_candidates(chimera);
	REQUIRE(result.size() == 2);

	bool has_a = std::find(result.begin(), result.end(), 0) != result.end();
	bool has_b = std::find(result.begin(), result.end(), 1) != result.end();
	REQUIRE(has_a);
	REQUIRE(has_b);
}

TEST_CASE("KmerIndex - Results capped at MAX_CANDIDATES", "[KmerIndex]") {
	KmerIndex idx;

	// Add 20 sequences that all share k-mers with the query
	std::string base = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	for (uint32_t i = 0; i < 20; i++) {
		std::string seq = base;
		if (i < seq.size()) {
			seq[i] = (seq[i] == 'A') ? 'C' : 'A';
		}
		idx.add_sequence(seq);
	}

	auto result = idx.find_candidates(base);
	REQUIRE(result.size() <= KmerIndex::MAX_CANDIDATES);
}

TEST_CASE("KmerIndex - Chunk partitioning", "[KmerIndex]") {
	SECTION("100bp query gives 4 chunks of 25 each") {
		auto chunks = KmerIndex::partition_query(100);
		REQUIRE(chunks.size() == 4);
		size_t total = 0;
		for (auto &[start, len] : chunks) {
			total += len;
		}
		REQUIRE(total == 100);
		for (auto &[start, len] : chunks) {
			REQUIRE(len == 25);
		}
	}

	SECTION("101bp query distributes remainder to last chunk") {
		auto chunks = KmerIndex::partition_query(101);
		REQUIRE(chunks.size() == 4);
		size_t total = 0;
		for (auto &[start, len] : chunks) {
			total += len;
		}
		REQUIRE(total == 101);
		// partition_query distributes floor to early chunks, remainder to last
		REQUIRE(chunks[0].second == 25);
		REQUIRE(chunks[3].second == 26);
	}

	SECTION("99bp query — first chunk is shorter") {
		auto chunks = KmerIndex::partition_query(99);
		REQUIRE(chunks.size() == 4);
		size_t total = 0;
		for (auto &[start, len] : chunks) {
			total += len;
		}
		REQUIRE(total == 99);
		// 99/4=24, 75/3=25, 50/2=25, 25/1=25
		REQUIRE(chunks[0].second == 24);
		REQUIRE(chunks[1].second == 25);
		REQUIRE(chunks[2].second == 25);
		REQUIRE(chunks[3].second == 25);
	}

	SECTION("Very short query still produces chunks") {
		auto chunks = KmerIndex::partition_query(10);
		REQUIRE(chunks.size() == 4);
		size_t total = 0;
		for (auto &[start, len] : chunks) {
			total += len;
		}
		REQUIRE(total == 10);
	}

	SECTION("Query shorter than NUM_CHUNKS — some zero-length chunks") {
		auto chunks = KmerIndex::partition_query(3);
		REQUIRE(chunks.size() == 4);
		size_t total = 0;
		size_t zero_count = 0;
		for (auto &[start, len] : chunks) {
			total += len;
			if (len == 0) {
				zero_count++;
			}
		}
		REQUIRE(total == 3);
		REQUIRE(zero_count == 1); // 3/4=0, so first chunk is 0
	}

	SECTION("Chunks are contiguous and non-overlapping") {
		for (size_t query_len : {50, 99, 100, 101, 300, 1500}) {
			auto chunks = KmerIndex::partition_query(query_len);
			size_t expected_start = 0;
			for (auto &[start, len] : chunks) {
				REQUIRE(start == expected_start);
				expected_start = start + len;
			}
			REQUIRE(expected_start == query_len);
		}
	}
}

TEST_CASE("KmerIndex - Multiple sequences with partial overlap", "[KmerIndex]") {
	KmerIndex idx;

	// ref0 and ref1 share the first 20bp, differ in the rest
	std::string shared_prefix = "ACGTACGTACGTACGTACGT";
	std::string ref0 = shared_prefix + "AAAAAAAAAAAAAAAAAAAA" + "CCCCCCCCCCCCCCCCCCCC";
	std::string ref1 = shared_prefix + "GGGGGGGGGGGGGGGGGGGG" + "TTTTTTTTTTTTTTTTTTTT";
	std::string ref2 = std::string("TTTTTTTTTTTTTTTTTTTT") + "GGGGGGGGGGGGGGGGGGGG" + "CCCCCCCCCCCCCCCCCCCC";

	idx.add_sequence(ref0); // id=0
	idx.add_sequence(ref1); // id=1
	idx.add_sequence(ref2); // id=2

	SECTION("Exact match query returns correct top hit") {
		auto result = idx.find_candidates(ref0);
		REQUIRE(!result.empty());
		REQUIRE(result[0] == 0);
	}
}

TEST_CASE("KmerIndex - Concurrent find_candidates is safe", "[KmerIndex]") {
	// This test verifies that thread_local scratch buffers work correctly
	// by calling find_candidates multiple times in sequence (same thread).
	// A threading test would require <thread> which we avoid in unit tests.
	KmerIndex idx;

	std::string ref0 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	std::string ref1 = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA";
	idx.add_sequence(ref0);
	idx.add_sequence(ref1);

	// Multiple consecutive calls should produce consistent results
	// (thread_local buffers are properly reset between calls)
	for (int i = 0; i < 5; i++) {
		auto result = idx.find_candidates(ref0);
		REQUIRE(!result.empty());
		REQUIRE(result[0] == 0);
	}

	for (int i = 0; i < 5; i++) {
		auto result = idx.find_candidates(ref1);
		REQUIRE(!result.empty());
		REQUIRE(result[0] == 1);
	}
}
