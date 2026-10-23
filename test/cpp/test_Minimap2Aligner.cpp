#include <catch2/catch_test_macros.hpp>
#include "Minimap2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include "minimap2_part_cursor.hpp"
#include "sequence_utils.hpp"
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <iterator>
#include <memory>
#include <mutex>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

// Helper to sum query-consuming CIGAR operations (M, I, S, =, X)
// D, N, H, P do not consume query bases
static int cigar_query_consumed(const std::string &cigar) {
	int total = 0;
	int num = 0;
	for (char c : cigar) {
		if (c >= '0' && c <= '9') {
			num = num * 10 + (c - '0');
		} else {
			if (c == 'M' || c == 'I' || c == 'S' || c == '=' || c == 'X') {
				total += num;
			}
			num = 0;
		}
	}
	return total;
}

using namespace miint;

// Helper to create a single unpaired query batch
static SequenceRecordBatch make_query_batch(const std::string &read_id, const std::string &sequence) {
	SequenceRecordBatch batch(false); // unpaired
	batch.read_ids.push_back(read_id);
	batch.comments.push_back("");
	batch.sequences1.push_back(sequence);
	batch.quals1.push_back(QualScore(""));
	return batch;
}

// Helper to create a paired query batch
static SequenceRecordBatch make_paired_query_batch(const std::string &read_id, const std::string &seq1,
                                                   const std::string &seq2) {
	SequenceRecordBatch batch(true); // paired
	batch.read_ids.push_back(read_id);
	batch.comments.push_back("");
	batch.sequences1.push_back(seq1);
	batch.sequences2.push_back(seq2);
	batch.quals1.push_back(QualScore(""));
	batch.quals2.push_back(QualScore(""));
	return batch;
}

TEST_CASE("Minimap2Aligner construction with default config", "[Minimap2Aligner]") {
	Minimap2Config config;
	REQUIRE(config.preset == "sr");
	REQUIRE(config.max_secondary == 5);
	REQUIRE(config.eqx == true);
	REQUIRE(config.k == 0); // 0 means use preset default
	REQUIRE(config.w == 0);

	Minimap2Aligner aligner(config);
	// Should not throw
}

TEST_CASE("Minimap2Aligner build index from subjects", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"ref1", "ACGTACGTACGTACGTACGTACGTACGTACGT"});
	subjects.push_back({"ref2", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"});

	REQUIRE_NOTHROW(aligner.build_index(subjects));
}

TEST_CASE("Minimap2Aligner single-end alignment - exact match", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// Create a subject sequence (must be longer than default k-mer size of 21)
	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	// Create a query that exactly matches part of the subject (50bp)
	auto queries = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	REQUIRE(batch.read_ids[0] == "query1");
	REQUIRE(batch.references[0] == "reference");
	REQUIRE(batch.positions[0] == 1); // 1-based position
	// With EQX mode, CIGAR should contain = for matches
	INFO("CIGAR: " << batch.cigars[0]);
	REQUIRE(batch.cigars[0].find('=') != std::string::npos);
}

TEST_CASE("Minimap2Aligner single-end alignment - with mismatch", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// 100bp reference - unique sequence with recognizable pattern
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	// Query: same as first 52bp of reference but with single mismatch near position 40
	// (far enough from the ends to have good minimizer anchors on both sides)
	// Original: ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT**AC**GTACGTACGT
	// Changed:  ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT**TT**GTACGTACGT
	std::string query = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTTTGTACGTACGT";
	auto queries = make_query_batch("query_mm", query);

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	REQUIRE(batch.read_ids[0] == "query_mm");
	// With EQX mode, CIGAR should contain X for mismatches
	INFO("CIGAR: " << batch.cigars[0]);
	REQUIRE(batch.cigars[0].find('X') != std::string::npos);
}

TEST_CASE("Minimap2Aligner unmapped query", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	// Completely different sequence (50bp of Ts)
	auto queries = make_query_batch("unmapped_query", "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT");

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// May have zero alignments or unmapped entry
	if (batch.size() > 0) {
		// If there's an entry, check the unmapped flag
		INFO("Flags: " << batch.flags[0]);
		bool is_unmapped = (batch.flags[0] & 0x4) != 0;
		if (!is_unmapped) {
			// If it mapped somehow, check mapq is low
			REQUIRE(batch.mapqs[0] < 10);
		}
	}
}

TEST_CASE("Minimap2Aligner max_secondary limits alignments", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.max_secondary = 2; // Only allow 2 secondary alignments
	Minimap2Aligner aligner(config);

	// Create a single reference with multiple similar regions (the query will map to
	// multiple positions within the same reference, creating secondary alignments)
	// 400bp reference with the query sequence repeated multiple times with gaps
	std::string query_part = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
	std::string spacer = "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
	std::string ref_seq =
	    query_part + spacer + query_part + spacer + query_part + spacer + query_part + spacer + query_part;
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	auto queries = make_query_batch("multi_hit", query_part);

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// Should have at most max_secondary + 1 (primary) = 3 alignments
	// minimap2's behavior with repetitive sequences may vary, so we use <=
	INFO("Number of alignments: " << batch.size());
	REQUIRE(batch.size() >= 1);                        // At least one alignment
	REQUIRE(batch.size() <= config.max_secondary + 1); // At most primary + secondaries
}

TEST_CASE("Minimap2Aligner SAM fields are populated correctly", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	auto queries = make_query_batch("test_query", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);

	// Check all required fields are populated
	REQUIRE(!batch.read_ids[0].empty());
	REQUIRE(!batch.references[0].empty());
	REQUIRE(batch.positions[0] > 0); // 1-based
	REQUIRE(batch.stop_positions[0] >= batch.positions[0]);
	REQUIRE(!batch.cigars[0].empty());
	// mapq should be reasonable (0-60)
	REQUIRE(batch.mapqs[0] <= 60);

	// For single-end, mate fields should be default/unmapped
	REQUIRE(batch.mate_references[0] == "*");
	REQUIRE(batch.mate_positions[0] == 0);
	REQUIRE(batch.template_lengths[0] == 0);

	// tag_as should have alignment score
	REQUIRE(batch.tag_as_values[0] > 0);
}

TEST_CASE("Minimap2Aligner different presets", "[Minimap2Aligner]") {
	std::vector<std::string> presets = {"sr", "map-ont", "map-pb"};

	// Use a longer reference (200bp) to support all presets
	// map-pb uses k=19, w=10; map-ont uses k=15, w=10; sr uses k=21, w=11
	std::string long_ref = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                       "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
	                       "ATATATATATATATATATATATATATATATATATATATATATATATATATAT"
	                       "GCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGC";

	for (const auto &preset : presets) {
		INFO("Testing preset: " << preset);
		Minimap2Config config;
		config.preset = preset;
		Minimap2Aligner aligner(config);

		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref", long_ref});
		REQUIRE_NOTHROW(aligner.build_index(subjects));
	}
}

TEST_CASE("Minimap2Aligner paired-end alignment", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// Create a longer reference for paired-end with unique 50bp regions separated by spacers
	// sr preset uses k=21, so we need sequences longer than 21bp
	std::string read1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"; // 52bp
	std::string read2_seq = "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"; // 52bp
	std::string spacer = std::string(200, 'N');                                     // N's won't match
	std::string long_ref = read1_seq + spacer + read2_seq;
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", long_ref});
	aligner.build_index(subjects);

	auto queries = make_paired_query_batch("paired_query",
	                                       read1_seq, // Matches at position 1
	                                       read2_seq  // Matches at position ~254
	);

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// Should have at least 2 records (one for each mate)
	INFO("Number of alignments: " << batch.size());
	REQUIRE(batch.size() >= 2);

	// Check paired flags
	bool found_read1 = false;
	bool found_read2 = false;
	for (size_t i = 0; i < batch.size(); i++) {
		uint16_t flags = batch.flags[i];
		bool is_paired = (flags & 0x1) != 0;
		bool is_read1 = (flags & 0x40) != 0;
		bool is_read2 = (flags & 0x80) != 0;

		if (is_paired) {
			if (is_read1)
				found_read1 = true;
			if (is_read2)
				found_read2 = true;
		}
	}
	REQUIRE(found_read1);
	REQUIRE(found_read2);
}

TEST_CASE("Minimap2Aligner build_single_index for per-subject mode", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	AlignmentSubject subject {"single_ref", ref_seq};

	REQUIRE_NOTHROW(aligner.build_single_index(subject));

	// Align against single-subject index with 52bp query
	auto queries = make_query_batch("query", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	REQUIRE(batch.references[0] == "single_ref");
}

TEST_CASE("AlignmentSubject length computed from sequence", "[Minimap2Aligner]") {
	AlignmentSubject subject {"test", "ACGTACGT"};
	REQUIRE(subject.length() == 8);

	AlignmentSubject empty_subject {"empty", ""};
	REQUIRE(empty_subject.length() == 0);
}

TEST_CASE("Minimap2Aligner CIGAR string generation", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.eqx = true; // Use =/X instead of M
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"ref", ref_seq});
	aligner.build_index(subjects);

	SECTION("exact match produces = operations") {
		// 52bp query that matches the reference
		auto queries = make_query_batch("exact", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

		SAMRecordBatch batch;
		aligner.align(queries, batch);

		REQUIRE(batch.size() >= 1);
		INFO("CIGAR: " << batch.cigars[0]);
		// Should have = for sequence matches
		REQUIRE(batch.cigars[0].find('=') != std::string::npos);
		// Should not have M (since EQX is enabled)
		REQUIRE(batch.cigars[0].find('M') == std::string::npos);
	}
}

TEST_CASE("Minimap2Aligner empty input handling", "[Minimap2Aligner]") {
	Minimap2Config config;
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"ref", ref_seq});
	aligner.build_index(subjects);

	SECTION("empty query batch produces no results") {
		SequenceRecordBatch queries(false);
		SAMRecordBatch batch;
		aligner.align(queries, batch);
		REQUIRE(batch.size() == 0);
	}
}

TEST_CASE("Minimap2Aligner multiple queries in batch", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);

	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	aligner.build_index(subjects);

	// 52bp query sequence
	std::string query_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

	// Create batch with multiple queries
	SequenceRecordBatch queries(false);
	queries.read_ids.push_back("query1");
	queries.read_ids.push_back("query2");
	queries.read_ids.push_back("query3");
	queries.comments.push_back("");
	queries.comments.push_back("");
	queries.comments.push_back("");
	queries.sequences1.push_back(query_seq);
	queries.sequences1.push_back(query_seq);
	queries.sequences1.push_back(query_seq);
	queries.quals1.push_back(QualScore(""));
	queries.quals1.push_back(QualScore(""));
	queries.quals1.push_back(QualScore(""));

	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// Should have at least 3 alignments (one per query)
	REQUIRE(batch.size() >= 3);

	// Check that all query IDs appear in results
	std::set<std::string> found_ids;
	for (const auto &id : batch.read_ids) {
		found_ids.insert(id);
	}
	REQUIRE(found_ids.count("query1") > 0);
	REQUIRE(found_ids.count("query2") > 0);
	REQUIRE(found_ids.count("query3") > 0);
}

// =============================================================================
// SharedMinimap2Index and shared index support tests
// =============================================================================

TEST_CASE("InitOptions produces valid options", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.eqx = true;
	config.max_secondary = 3;

	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt);

	// Check preset-specific values for sr
	REQUIRE(iopt.k == 21);
	REQUIRE(iopt.w == 11);
	// Check CIGAR output enabled
	REQUIRE((mopt.flag & MM_F_CIGAR) != 0);
	// Check EQX enabled
	REQUIRE((mopt.flag & MM_F_EQX) != 0);
	// Check MD tag output
	REQUIRE((mopt.flag & MM_F_OUT_MD) != 0);
	// Check best_n = max_secondary (matches minimap2 CLI default behavior)
	REQUIRE(mopt.best_n == 3);
}

TEST_CASE("InitOptions with custom k and w", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 15;
	config.w = 10;

	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt);

	REQUIRE(iopt.k == 15);
	REQUIRE(iopt.w == 10);
}

TEST_CASE("LoadIndexFromFile loads test .mmi correctly", "[Minimap2Aligner]") {
	// First build and save a test index
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	Minimap2Aligner builder(config);

	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"test_ref", ref_seq});
	builder.build_index(subjects);
	builder.save_index("data/shards/test_load_helper.mmi");

	// Now load it with the static helper
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt);

	mm_idx_t *idx = nullptr;
	std::vector<std::string> names;
	REQUIRE_NOTHROW(Minimap2Aligner::LoadIndexFromFile("data/shards/test_load_helper.mmi", iopt, idx, names));

	REQUIRE(idx != nullptr);
	REQUIRE(names.size() == 1);
	REQUIRE(names[0] == "test_ref");

	mm_idx_destroy(idx);
}

TEST_CASE("SharedMinimap2Index construction and accessors", "[Minimap2Aligner]") {
	// Build and save a test index first
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	Minimap2Aligner builder(config);

	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"shared_ref", ref_seq});
	builder.build_index(subjects);
	builder.save_index("data/shards/test_shared_idx.mmi");

	// Construct SharedMinimap2Index
	SharedMinimap2Index shared_idx("data/shards/test_shared_idx.mmi", config);

	REQUIRE(shared_idx.index() != nullptr);
	REQUIRE(shared_idx.subject_names().size() == 1);
	REQUIRE(shared_idx.subject_names()[0] == "shared_ref");
	// mapopt should have CIGAR flag set
	REQUIRE((shared_idx.mapopt().flag & MM_F_CIGAR) != 0);
}

// Builds a genuine multi-part .mmi fixture by dumping two independently-built
// single-part indexes and concatenating their bytes. mm_idx_dump (called once
// by save_index per file) writes a self-contained MM_IDX_MAGIC-prefixed block;
// minimap2's own `-I <batch>` CLI flag produces a multi-part file by calling
// mm_idx_dump repeatedly into the SAME fp, so byte-concatenating two
// independently-dumped single-part files reproduces that exact on-disk layout.
// Returns the path to the multi-part file; caller owns cleanup of all three
// paths (part1_path, part2_path, and the returned multi-part path).
static std::string build_multipart_mmi_fixture(const std::string &part1_path, const std::string &part2_path,
                                               const std::string &multipart_path) {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	std::string ref_seq1 = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                       "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::string ref_seq2 = "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"
	                       "AAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAACCCCGGGGTTTTAAAA";

	Minimap2Aligner builder1(config);
	builder1.build_index({{"part1_ref", ref_seq1}});
	builder1.save_index(part1_path);

	Minimap2Aligner builder2(config);
	builder2.build_index({{"part2_ref", ref_seq2}});
	builder2.save_index(part2_path);

	std::ifstream in1(part1_path, std::ios::binary);
	std::ifstream in2(part2_path, std::ios::binary);
	std::ofstream out(multipart_path, std::ios::binary | std::ios::trunc);
	REQUIRE(in1.good());
	REQUIRE(in2.good());
	REQUIRE(out.good());
	out << in1.rdbuf();
	out << in2.rdbuf();
	out.close();

	return multipart_path;
}

TEST_CASE("Minimap2IndexReader reads a multi-part index part by part", "[Minimap2Aligner]") {
	const std::string part1 = "data/shards/test_multipart_reader_p1.mmi";
	const std::string part2 = "data/shards/test_multipart_reader_p2.mmi";
	const std::string multipart = "data/shards/test_multipart_reader.mmi";
	build_multipart_mmi_fixture(part1, part2, multipart);

	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	Minimap2IndexReader reader(multipart, config);

	auto p1 = reader.ReadNextPart();
	REQUIRE(p1 != nullptr);
	REQUIRE(p1->subject_names().size() == 1);
	REQUIRE(p1->subject_names()[0] == "part1_ref");
	// A genuinely multi-part fixture: the reader must NOT be at eof after part 1,
	// or this test (and the design it verifies) would be vacuous.
	REQUIRE_FALSE(reader.AtEof());

	auto p2 = reader.ReadNextPart();
	REQUIRE(p2 != nullptr);
	REQUIRE(p2->subject_names().size() == 1);
	REQUIRE(p2->subject_names()[0] == "part2_ref");
	REQUIRE(reader.AtEof());

	auto p3 = reader.ReadNextPart();
	REQUIRE(p3 == nullptr);

	std::remove(part1.c_str());
	std::remove(part2.c_str());
	std::remove(multipart.c_str());
}

TEST_CASE("LoadIndexFromFile throws loud on a multi-part index instead of silently truncating", "[Minimap2Aligner]") {
	// Regression for the bug this streaming feature fixes: before
	// Minimap2IndexReader existed, LoadIndexFromFile (and therefore
	// SharedMinimap2Index(path, config), which align_minimap2_sharded used to
	// load every shard) silently returned only the first part of a multi-part
	// .mmi with no error, dropping every reference in later parts. Both table
	// functions now go through Minimap2PartCursor instead; this guards any
	// other caller of the single-index loader.
	const std::string part1 = "data/shards/test_multipart_loadfromfile_p1.mmi";
	const std::string part2 = "data/shards/test_multipart_loadfromfile_p2.mmi";
	const std::string multipart = "data/shards/test_multipart_loadfromfile.mmi";
	build_multipart_mmi_fixture(part1, part2, multipart);

	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt);

	mm_idx_t *idx = nullptr;
	std::vector<std::string> names;
	REQUIRE_THROWS_AS(Minimap2Aligner::LoadIndexFromFile(multipart, iopt, idx, names), std::runtime_error);

	// SharedMinimap2Index(path, config) goes through LoadIndexFromFile and must
	// throw too.
	REQUIRE_THROWS_AS(SharedMinimap2Index(multipart, config), std::runtime_error);

	std::remove(part1.c_str());
	std::remove(part2.c_str());
	std::remove(multipart.c_str());
}

TEST_CASE("Minimap2PartCursor hands a single-part index straight through", "[Minimap2Aligner]") {
	const std::string part1 = "data/shards/test_part_cursor_single_p1.mmi";
	const std::string part2 = "data/shards/test_part_cursor_single_p2.mmi";
	const std::string multipart = "data/shards/test_part_cursor_single.mmi";
	build_multipart_mmi_fixture(part1, part2, multipart);

	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	Minimap2PartCursor cursor(part1, config, nullptr);
	REQUIRE_FALSE(cursor.IsMultiPart());
	auto idx = cursor.ReleaseSinglePart();
	REQUIRE(idx != nullptr);
	REQUIRE(idx->subject_names() == std::vector<std::string> {"part1_ref"});

	std::remove(part1.c_str());
	std::remove(part2.c_str());
	std::remove(multipart.c_str());
}

// The coordination protocol align_minimap2 and align_minimap2_sharded both rely
// on: every worker must see every part exactly once, exactly one thread loads
// each next part, the freed-memory hook fires on every detach, and the cursor
// reports exhaustion only once the reader has no part left.
TEST_CASE("Minimap2PartCursor walks concurrent workers through every part exactly once", "[Minimap2Aligner]") {
	const std::string part1 = "data/shards/test_part_cursor_multi_p1.mmi";
	const std::string part2 = "data/shards/test_part_cursor_multi_p2.mmi";
	const std::string multipart = "data/shards/test_part_cursor_multi.mmi";
	build_multipart_mmi_fixture(part1, part2, multipart);

	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	std::atomic<int> flushes {0};
	Minimap2PartCursor cursor(multipart, config, [&]() { flushes++; });
	REQUIRE(cursor.IsMultiPart());

	// A query that maps in part 1 only and one that maps in part 2 only: what a
	// worker aligns against each part tells us which part it was attached to.
	SequenceRecordBatch queries;
	queries.is_paired = false;
	queries.read_ids = {"q_part1", "q_part2"};
	queries.comments = {"", ""};
	queries.sequences1 = {"ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
	                      "TTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"};
	queries.quals1 = {{}, {}};

	std::atomic<int> publishes {0};
	std::atomic<int> attached_advances {0};
	constexpr int kWorkers = 4;
	std::vector<std::vector<std::string>> seen(kWorkers); // per worker: mapped reference per part visited
	std::vector<std::thread> workers;
	for (int w = 0; w < kWorkers; w++) {
		workers.emplace_back([&, w]() {
			Minimap2Aligner aligner(config);
			Minimap2PartCursor::Attachment att;
			while (true) {
				// Claim through the cursor exactly as the table functions do; the
				// lambda stands in for their per-part work claim.
				cursor.WithCurrentPart(att, aligner, []() {});
				// att.attached is false only while a leader is mid-transition; the
				// table functions then find their work source exhausted and go
				// straight to Advance, which is what this does too.
				if (att.attached) {
					SAMRecordBatch out;
					aligner.align(queries, out);
					// Distinct references only: secondary/supplementary records repeat
					// the same reference, and each part holds exactly one.
					std::set<std::string> refs;
					for (size_t i = 0; i < out.size(); i++) {
						if ((out.flags[i] & 0x4) == 0) {
							refs.insert(out.references[i]);
						}
					}
					seen[w].insert(seen[w].end(), refs.begin(), refs.end());
				}
				if (att.attached) {
					attached_advances++;
				}
				if (!cursor.Advance(att, aligner, nullptr, [&]() { publishes++; })) {
					break;
				}
			}
		});
	}
	for (auto &t : workers) {
		t.join();
	}

	// One transition (part 1 -> part 2), performed by exactly one leader.
	REQUIRE(publishes.load() == 1);
	// Parts are visited in file order and never twice by the same worker. A
	// worker may legitimately miss a part (one that started after the first
	// transition, or that arrived while a leader was mid-transition) -- the
	// table functions share their work source across workers, so coverage is a
	// property of the group, not of each thread. Between them the workers must
	// have aligned against both parts.
	std::set<std::string> union_seen;
	for (int w = 0; w < kWorkers; w++) {
		REQUIRE(seen[w].size() <= 2);
		if (seen[w].size() == 2) {
			REQUIRE(seen[w] == std::vector<std::string> {"part1_ref", "part2_ref"});
		}
		union_seen.insert(seen[w].begin(), seen[w].end());
	}
	REQUIRE(union_seen == std::set<std::string> {"part1_ref", "part2_ref"});
	// Every thread flushes its OWN detach, not just the leader's: whichever
	// thread drops the last reference to a part is the one that actually frees
	// it, and that is rarely the leader. So the hook must fire at least once per
	// Advance that had a part attached to give up, on top of the leader's own
	// post-reset flush per transition. Asserted as a floor, not an exact count,
	// so tightening when the flush fires stays a free change.
	REQUIRE(flushes.load() >= attached_advances.load());
	REQUIRE(flushes.load() >= publishes.load() + 1);

	// Exhausted stays exhausted, without touching the reader again. A thread
	// arriving after the last part is gone attaches to nothing and is sent
	// straight back out by Advance.
	Minimap2Aligner late(config);
	Minimap2PartCursor::Attachment late_att;
	cursor.WithCurrentPart(late_att, late, []() {});
	REQUIRE_FALSE(late_att.attached);
	REQUIRE_FALSE(cursor.Advance(late_att, late, nullptr, nullptr));

	std::remove(part1.c_str());
	std::remove(part2.c_str());
	std::remove(multipart.c_str());
}

TEST_CASE("Two aligners sharing same SharedMinimap2Index produce identical results", "[Minimap2Aligner]") {
	// Build and save a test index
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	config.max_secondary = 0;
	Minimap2Aligner builder(config);

	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	builder.build_index(subjects);
	builder.save_index("data/shards/test_shared_align.mmi");

	// Create shared index
	auto shared_idx = std::make_shared<SharedMinimap2Index>("data/shards/test_shared_align.mmi", config);

	// Create two aligners and attach the same shared index
	Minimap2Aligner aligner1(config);
	Minimap2Aligner aligner2(config);
	aligner1.attach_shared_index(shared_idx);
	aligner2.attach_shared_index(shared_idx);

	auto queries = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch1, batch2;
	aligner1.align(queries, batch1);
	aligner2.align(queries, batch2);

	REQUIRE(batch1.size() == batch2.size());
	REQUIRE(batch1.size() >= 1);
	REQUIRE(batch1.references[0] == batch2.references[0]);
	REQUIRE(batch1.positions[0] == batch2.positions[0]);
	REQUIRE(batch1.cigars[0] == batch2.cigars[0]);
}

TEST_CASE("attach_shared_index clears owned index", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	Minimap2Aligner aligner(config);

	// Build an owned index
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"owned_ref", ref_seq});
	aligner.build_index(subjects);

	// Align with owned index
	auto queries = make_query_batch("q1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
	SAMRecordBatch batch1;
	aligner.align(queries, batch1);
	REQUIRE(batch1.size() >= 1);
	REQUIRE(batch1.references[0] == "owned_ref");

	// Save and create shared index with different ref name
	aligner.save_index("data/shards/test_owned_clear.mmi");
	auto shared_idx = std::make_shared<SharedMinimap2Index>("data/shards/test_owned_clear.mmi", config);

	// Attach shared index - should clear owned
	aligner.attach_shared_index(shared_idx);

	SAMRecordBatch batch2;
	aligner.align(queries, batch2);
	REQUIRE(batch2.size() >= 1);
	// After attach, ref name comes from shared index (same .mmi file, so same name)
	REQUIRE(batch2.references[0] == "owned_ref");
}

TEST_CASE("load_index clears shared index", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	// Build and save two indexes with different ref names
	Minimap2Aligner builder(config);
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";

	std::vector<AlignmentSubject> subjects1;
	subjects1.push_back({"ref_shared", ref_seq});
	builder.build_index(subjects1);
	builder.save_index("data/shards/test_load_clear_shared.mmi");

	std::vector<AlignmentSubject> subjects2;
	subjects2.push_back({"ref_owned", ref_seq});
	builder.build_index(subjects2);
	builder.save_index("data/shards/test_load_clear_owned.mmi");

	// Aligner starts with shared index
	Minimap2Aligner aligner(config);
	auto shared_idx = std::make_shared<SharedMinimap2Index>("data/shards/test_load_clear_shared.mmi", config);
	aligner.attach_shared_index(shared_idx);

	auto queries = make_query_batch("q1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
	SAMRecordBatch batch1;
	aligner.align(queries, batch1);
	REQUIRE(batch1.size() >= 1);
	REQUIRE(batch1.references[0] == "ref_shared");

	// load_index should clear shared and use owned
	aligner.load_index("data/shards/test_load_clear_owned.mmi");

	SAMRecordBatch batch2;
	aligner.align(queries, batch2);
	REQUIRE(batch2.size() >= 1);
	REQUIRE(batch2.references[0] == "ref_owned");
}

TEST_CASE("detach then re-attach different shared index", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;

	// Build two indexes
	Minimap2Aligner builder(config);
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";

	std::vector<AlignmentSubject> subjects1;
	subjects1.push_back({"ref_A", ref_seq});
	builder.build_index(subjects1);
	builder.save_index("data/shards/test_reattach_A.mmi");

	std::vector<AlignmentSubject> subjects2;
	subjects2.push_back({"ref_B", ref_seq});
	builder.build_index(subjects2);
	builder.save_index("data/shards/test_reattach_B.mmi");

	auto shared_A = std::make_shared<SharedMinimap2Index>("data/shards/test_reattach_A.mmi", config);
	auto shared_B = std::make_shared<SharedMinimap2Index>("data/shards/test_reattach_B.mmi", config);

	Minimap2Aligner aligner(config);
	auto queries = make_query_batch("q1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	// Attach A
	aligner.attach_shared_index(shared_A);
	SAMRecordBatch batch1;
	aligner.align(queries, batch1);
	REQUIRE(batch1.size() >= 1);
	REQUIRE(batch1.references[0] == "ref_A");

	// Detach
	aligner.detach_shared_index();

	// Attach B
	aligner.attach_shared_index(shared_B);
	SAMRecordBatch batch2;
	aligner.align(queries, batch2);
	REQUIRE(batch2.size() >= 1);
	REQUIRE(batch2.references[0] == "ref_B");
}

TEST_CASE("Concurrent alignment on shared index from two threads", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	config.max_secondary = 0;
	Minimap2Aligner builder(config);

	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});
	builder.build_index(subjects);
	builder.save_index("data/shards/test_concurrent.mmi");

	auto shared_idx = std::make_shared<SharedMinimap2Index>("data/shards/test_concurrent.mmi", config);

	auto queries = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

	SAMRecordBatch batch1, batch2;
	bool ok1 = false, ok2 = false;

	// Run two alignments concurrently on separate aligners sharing one index
	std::thread t1([&]() {
		Minimap2Aligner aligner(config);
		aligner.attach_shared_index(shared_idx);
		aligner.align(queries, batch1);
		ok1 = (batch1.size() >= 1 && batch1.references[0] == "reference");
	});
	std::thread t2([&]() {
		Minimap2Aligner aligner(config);
		aligner.attach_shared_index(shared_idx);
		aligner.align(queries, batch2);
		ok2 = (batch2.size() >= 1 && batch2.references[0] == "reference");
	});

	t1.join();
	t2.join();

	REQUIRE(ok1);
	REQUIRE(ok2);
	REQUIRE(batch1.size() == batch2.size());
	REQUIRE(batch1.positions[0] == batch2.positions[0]);
	REQUIRE(batch1.cigars[0] == batch2.cigars[0]);
}

// Clean up temporary .mmi files created by tests
TEST_CASE("Cleanup temp .mmi files", "[Minimap2Aligner]") {
	std::vector<std::string> temp_files = {
	    "data/shards/test_load_helper.mmi", "data/shards/test_shared_idx.mmi", "data/shards/test_shared_align.mmi",
	    "data/shards/test_owned_clear.mmi", "data/shards/test_load_clear_shared.mmi",
	    "data/shards/test_load_clear_owned.mmi", "data/shards/test_reattach_A.mmi", "data/shards/test_reattach_B.mmi",
	    "data/shards/test_concurrent.mmi",
	    // Multi-part fixtures. Each test removes its own three paths on the way
	    // out; these repeats catch the set a failed REQUIRE aborted past.
	    "data/shards/test_multipart_reader_p1.mmi", "data/shards/test_multipart_reader_p2.mmi",
	    "data/shards/test_multipart_reader.mmi", "data/shards/test_multipart_loadfromfile_p1.mmi",
	    "data/shards/test_multipart_loadfromfile_p2.mmi", "data/shards/test_multipart_loadfromfile.mmi",
	    "data/shards/test_part_cursor_single_p1.mmi", "data/shards/test_part_cursor_single_p2.mmi",
	    "data/shards/test_part_cursor_single.mmi", "data/shards/test_part_cursor_multi_p1.mmi",
	    "data/shards/test_part_cursor_multi_p2.mmi", "data/shards/test_part_cursor_multi.mmi"};
	for (const auto &path : temp_files) {
		std::remove(path.c_str());
	}
}

TEST_CASE("Minimap2Aligner CIGAR includes soft clips for partial query alignment", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.eqx = true;
	Minimap2Aligner aligner(config);

	// Reference sequence
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	aligner.build_index({{"reference", ref_seq}});

	// Query: 25bp of non-matching poly-T + full reference sequence + 25bp of non-matching poly-T
	// The reference portion should align perfectly; the poly-T flanks should be soft-clipped
	std::string junk(25, 'T');
	std::string query = junk + ref_seq + junk;

	auto queries = make_query_batch("clipped_query", query);
	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	INFO("CIGAR: " << batch.cigars[0]);

	// The CIGAR must contain soft clips (S) for the non-matching flanking regions
	REQUIRE(batch.cigars[0].find('S') != std::string::npos);

	// Sum of query-consuming CIGAR operations (M/I/S/=/X) must equal query length
	REQUIRE(cigar_query_consumed(batch.cigars[0]) == static_cast<int>(query.length()));
}

TEST_CASE("Minimap2Aligner CIGAR soft clips with map-hifi preset", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "map-hifi";
	config.eqx = true;
	Minimap2Aligner aligner(config);

	// 300+ bp reference required because map-hifi sets a=1 and min_dp_max=200,
	// so alignment score must reach 200 which needs at least 200bp of exact match
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
	                      "AATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATTCCGGAATT"
	                      "TGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCATGCA"
	                      "GCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTAGCTA"
	                      "CGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGATCGAT";
	aligner.build_index({{"reference", ref_seq}});

	// Query: 30bp non-matching + full reference + 30bp non-matching
	std::string junk(30, 'T');
	std::string query = junk + ref_seq + junk;

	auto queries = make_query_batch("hifi_clipped", query);
	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	INFO("CIGAR: " << batch.cigars[0]);

	REQUIRE(batch.cigars[0].find('S') != std::string::npos);
	REQUIRE(cigar_query_consumed(batch.cigars[0]) == static_cast<int>(query.length()));
}

TEST_CASE("Minimap2Aligner CIGAR soft clips on reverse strand", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.eqx = true;
	Minimap2Aligner aligner(config);

	// Reference sequence
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	aligner.build_index({{"reference", ref_seq}});

	// Query: poly-T flank + reverse complement of reference + poly-T flank
	// This forces alignment on the reverse strand (reg->rev = 1),
	// exercising the clip_front/clip_back strand swap logic
	std::string junk(25, 'T');
	std::string query = junk + miint::dna_reverse_complement(ref_seq) + junk;

	auto queries = make_query_batch("rev_clipped", query);
	SAMRecordBatch batch;
	aligner.align(queries, batch);

	REQUIRE(batch.size() >= 1);
	INFO("CIGAR: " << batch.cigars[0]);
	INFO("Flags: " << batch.flags[0]);

	// Must be on reverse strand
	REQUIRE((batch.flags[0] & 0x10) != 0);

	// CIGAR must contain soft clips
	REQUIRE(batch.cigars[0].find('S') != std::string::npos);

	// Query-consuming CIGAR operations must equal query length
	REQUIRE(cigar_query_consumed(batch.cigars[0]) == static_cast<int>(query.length()));
}

// =============================================================================
// Coverage pre-filter tests
// =============================================================================

TEST_CASE("Minimap2Config min_chain_coverage defaults to 0.0", "[Minimap2Aligner]") {
	Minimap2Config config;
	REQUIRE(config.min_chain_coverage == 0.0f);
}

TEST_CASE("InitOptions passes min_chain_coverage to mopt", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	config.min_chain_coverage = 0.75f;

	mm_idxopt_t iopt;
	mm_mapopt_t mopt;
	Minimap2Aligner::InitOptions(config, iopt, mopt);

	REQUIRE(mopt.min_chain_coverage == 0.75f);
}

TEST_CASE("min_chain_coverage=0.0 produces identical results to default", "[Minimap2Aligner]") {
	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});

	// 50bp exact-match query
	std::string query_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

	// Align with default config (min_chain_coverage=0.0)
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);
	aligner.build_index(subjects);

	auto queries = make_query_batch("q1", query_seq);
	SAMRecordBatch batch_default;
	aligner.align(queries, batch_default);

	// Align with explicit min_chain_coverage=0.0
	Minimap2Config config2;
	config2.preset = "sr";
	config2.min_chain_coverage = 0.0f;
	Minimap2Aligner aligner2(config2);
	aligner2.build_index(subjects);

	auto queries2 = make_query_batch("q1", query_seq);
	SAMRecordBatch batch_explicit;
	aligner2.align(queries2, batch_explicit);

	REQUIRE(batch_default.size() == batch_explicit.size());
	for (size_t i = 0; i < batch_default.size(); i++) {
		REQUIRE(batch_default.cigars[i] == batch_explicit.cigars[i]);
		REQUIRE(batch_default.positions[i] == batch_explicit.positions[i]);
	}
}

TEST_CASE("min_chain_coverage=0.99 filters unmappable queries", "[Minimap2Aligner]") {
	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});

	// Random query that does NOT match the reference
	std::string query_seq = "TTTTTTTTTTTTTTTTTTTTTTAAAAAAAAAAAAAAAAAAAACCCCCCCCCC";

	// Without filter: may produce low-quality mapping or unmapped
	Minimap2Config config_base;
	config_base.preset = "sr";
	Minimap2Aligner aligner_base(config_base);
	aligner_base.build_index(subjects);

	auto queries_base = make_query_batch("q_random", query_seq);
	SAMRecordBatch batch_base;
	aligner_base.align(queries_base, batch_base);

	// With filter at 0.99: random sequence has no seed hits → no chains → 0 coverage → filtered
	Minimap2Config config;
	config.preset = "sr";
	config.min_chain_coverage = 0.99f;
	Minimap2Aligner aligner(config);
	aligner.build_index(subjects);

	auto queries = make_query_batch("q_random", query_seq);
	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// Random sequence cannot produce chains against this reference, so zero mapped reads
	int mapped_filtered = 0;
	for (size_t i = 0; i < batch.size(); i++) {
		if ((batch.flags[i] & 0x4) == 0)
			mapped_filtered++;
	}
	REQUIRE(mapped_filtered == 0);
}

TEST_CASE("Coverage filter preserves full-length alignments", "[Minimap2Aligner]") {
	// 100bp reference
	std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
	                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
	std::vector<AlignmentSubject> subjects;
	subjects.push_back({"reference", ref_seq});

	// 50bp exact match (100% coverage chain)
	std::string query_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";

	Minimap2Config config;
	config.preset = "sr";
	config.min_chain_coverage = 0.5f;
	Minimap2Aligner aligner(config);
	aligner.build_index(subjects);

	auto queries = make_query_batch("q_full", query_seq);
	SAMRecordBatch batch;
	aligner.align(queries, batch);

	// Full-length alignment should survive even at 0.5 threshold
	REQUIRE(batch.size() >= 1);
	bool has_mapped = false;
	for (size_t i = 0; i < batch.size(); i++) {
		if ((batch.flags[i] & 0x4) == 0) {
			has_mapped = true;
			REQUIRE(batch.references[i] == "reference");
		}
	}
	REQUIRE(has_mapped);
}

// =============================================================================
// any_mapped: the keep/drop decision of a host filter, without SAM output
// =============================================================================

// A filter that drops every read minimap2 maps must decide exactly as
// align_minimap2(...) + ANTI JOIN on read_id does, or switching a pipeline to the
// filter would change which reads survive. align() emits at least one row for a
// query exactly when some region has a valid reference id, so any_mapped is
// checked against align()'s rows rather than against hand-picked expectations.

static std::string random_dna(std::mt19937 &rng, size_t length) {
	static const char bases[] = "ACGT";
	std::uniform_int_distribution<int> pick(0, 3);
	std::string seq(length, 'A');
	for (auto &c : seq) {
		c = bases[pick(rng)];
	}
	return seq;
}

// Substitutes each base with probability `rate`, so reads range from clean hits to
// ones minimap2 can no longer place.
static std::string mutate(std::mt19937 &rng, std::string seq, double rate) {
	static const char bases[] = "ACGT";
	std::uniform_real_distribution<double> coin(0.0, 1.0);
	std::uniform_int_distribution<int> pick(0, 3);
	for (auto &c : seq) {
		if (coin(rng) < rate) {
			c = bases[pick(rng)];
		}
	}
	return seq;
}

struct AnyMappedFixture {
	std::mt19937 rng {20261021};
	std::vector<AlignmentSubject> subjects;
	Minimap2Config config;

	AnyMappedFixture() {
		config.preset = "sr";
		subjects.push_back({"ref_a", random_dna(rng, 6000)});
		subjects.push_back({"ref_b", random_dna(rng, 6000)});
	}

	// A 150-base read: drawn from a reference (either strand, mutated at `rate`), or
	// random when `rate` is negative.
	std::string read(double rate) {
		if (rate < 0) {
			return random_dna(rng, 150);
		}
		const auto &ref = subjects[rng() % subjects.size()].sequence;
		std::string seq = ref.substr(rng() % (ref.size() - 150), 150);
		if (rng() % 2) {
			seq = dna_reverse_complement(seq);
		}
		return mutate(rng, seq, rate);
	}

	// A 150-base read of random flanks around `length` bases from a reference. Short matches
	// still chain, but the aligned score can fall below the DP cutoff, which is where a decision
	// made from chains alone would differ from align()'s.
	std::string partial(size_t length) {
		const auto &ref = subjects[rng() % subjects.size()].sequence;
		const std::string core = ref.substr(rng() % (ref.size() - length), length);
		const size_t left = (150 - length) / 2;
		return random_dna(rng, left) + core + random_dna(rng, 150 - length - left);
	}
};

TEST_CASE("any_mapped is true for a read from the reference and false for a random one", "[Minimap2Aligner]") {
	AnyMappedFixture fx;
	Minimap2Aligner aligner(fx.config);
	aligner.build_index(fx.subjects);

	const std::string hit = fx.subjects[0].sequence.substr(1000, 150);
	const std::string miss = random_dna(fx.rng, 150);
	REQUIRE(aligner.any_mapped(hit, "", "hit"));
	REQUIRE_FALSE(aligner.any_mapped(miss, "", "miss"));
}

TEST_CASE("any_mapped treats a pair as mapped when either mate maps", "[Minimap2Aligner]") {
	AnyMappedFixture fx;
	Minimap2Aligner aligner(fx.config);
	aligner.build_index(fx.subjects);

	const std::string &ref = fx.subjects[1].sequence;
	const std::string r1_hit = ref.substr(2000, 150);
	const std::string r2_hit = dna_reverse_complement(ref.substr(2250, 150));
	const std::string r1_miss = random_dna(fx.rng, 150);
	const std::string r2_miss = random_dna(fx.rng, 150);

	// A host filter drops the whole pair if any mate is host: the mate that does
	// not map must not rescue it.
	REQUIRE(aligner.any_mapped(r1_hit, r2_hit, "both"));
	REQUIRE(aligner.any_mapped(r1_miss, r2_hit, "r2_only"));
	REQUIRE(aligner.any_mapped(r1_hit, r2_miss, "r1_only"));
	REQUIRE_FALSE(aligner.any_mapped(r1_miss, r2_miss, "neither"));
}

TEST_CASE("any_mapped routes an empty mate the way align() does", "[Minimap2Aligner]") {
	AnyMappedFixture fx;
	Minimap2Aligner aligner(fx.config);
	aligner.build_index(fx.subjects);

	// align() takes the single-end path when sequence2 is empty, so an empty mate
	// is a single-end read, not a pair with a missing mate.
	const std::string hit = fx.subjects[0].sequence.substr(300, 150);
	REQUIRE(aligner.any_mapped(hit, "", "se") == aligner.any_mapped(hit, std::string_view(), "se"));
	REQUIRE(aligner.any_mapped(hit, "", "se"));
	// minimap2 is never called on an empty query; align() emits nothing for it.
	REQUIRE_FALSE(aligner.any_mapped("", "", "empty"));
}

TEST_CASE("any_mapped requires an index", "[Minimap2Aligner]") {
	Minimap2Config config;
	config.preset = "sr";
	Minimap2Aligner aligner(config);
	REQUIRE_THROWS_AS(aligner.any_mapped("ACGTACGTACGT", "", "q"), std::runtime_error);
}

TEST_CASE("any_mapped agrees with align() on a mixed single-end batch", "[Minimap2Aligner]") {
	AnyMappedFixture fx;
	Minimap2Aligner aligner(fx.config);
	aligner.build_index(fx.subjects);

	// Mutation rates straddle the point where minimap2 stops placing a read, so the
	// batch holds clear hits, clear misses and marginal reads.
	const std::vector<double> rates = {0.0, 0.02, 0.05, 0.10, 0.15, 0.20, 0.30, -1.0};
	SequenceRecordBatch queries(false);
	for (int i = 0; i < 240; i++) {
		queries.read_ids.push_back("se" + std::to_string(i));
		queries.comments.push_back("");
		queries.sequences1.push_back(fx.read(rates[i % rates.size()]));
		queries.quals1.push_back(QualScore(""));
	}
	for (int i = 0; i < 160; i++) {
		queries.read_ids.push_back("se_partial" + std::to_string(i));
		queries.comments.push_back("");
		queries.sequences1.push_back(fx.partial(16 + (i % 50)));
		queries.quals1.push_back(QualScore(""));
	}

	SAMRecordBatch rows;
	aligner.align(queries, rows);
	const std::set<std::string> aligned(rows.read_ids.begin(), rows.read_ids.end());

	size_t mapped = 0;
	for (size_t i = 0; i < queries.size(); i++) {
		const bool expected = aligned.count(queries.read_ids[i]) > 0;
		INFO("read " << queries.read_ids[i]);
		REQUIRE(aligner.any_mapped(queries.sequences1[i], "", queries.read_ids[i].c_str()) == expected);
		mapped += expected;
	}
	// Both outcomes must be well represented, or agreement proves nothing.
	REQUIRE(mapped >= 40);
	REQUIRE(queries.size() - mapped >= 40);
}

TEST_CASE("any_mapped agrees with align() on a mixed paired batch", "[Minimap2Aligner]") {
	AnyMappedFixture fx;
	Minimap2Aligner aligner(fx.config);
	aligner.build_index(fx.subjects);

	const std::vector<double> rates = {0.0, 0.05, 0.10, 0.20, -1.0};
	SequenceRecordBatch queries(true);
	for (int i = 0; i < 200; i++) {
		queries.read_ids.push_back("pe" + std::to_string(i));
		queries.comments.push_back("");
		queries.sequences1.push_back(fx.read(rates[i % rates.size()]));
		queries.sequences2.push_back(fx.read(rates[(i / rates.size()) % rates.size()]));
		queries.quals1.push_back(QualScore(""));
		queries.quals2.push_back(QualScore(""));
	}
	for (int i = 0; i < 120; i++) {
		queries.read_ids.push_back("pe_partial" + std::to_string(i));
		queries.comments.push_back("");
		queries.sequences1.push_back(fx.partial(16 + (i % 50)));
		queries.sequences2.push_back(i % 2 ? fx.read(-1.0) : fx.partial(16 + ((i * 7) % 50)));
		queries.quals1.push_back(QualScore(""));
		queries.quals2.push_back(QualScore(""));
	}

	SAMRecordBatch rows;
	aligner.align(queries, rows);
	const std::set<std::string> aligned(rows.read_ids.begin(), rows.read_ids.end());

	size_t mapped = 0;
	for (size_t i = 0; i < queries.size(); i++) {
		const bool expected = aligned.count(queries.read_ids[i]) > 0;
		INFO("pair " << queries.read_ids[i]);
		REQUIRE(aligner.any_mapped(queries.sequences1[i], queries.sequences2[i], queries.read_ids[i].c_str()) ==
		        expected);
		mapped += expected;
	}
	REQUIRE(mapped >= 40);
	REQUIRE(queries.size() - mapped >= 20);
}

// A prebuilt index can be read on several threads, which build its buckets' hash tables, and must
// come out the index read on one. mm_idx_dump writes every hash entry in table order, so two loads
// that dump alike hold the same sequences, the same buckets and every entry in the same place: a
// bucket built from the wrong entries, from part of them, or not at all shows as a difference.
// The reference repeats pieces of itself, so that many minimizers occur more than once and the
// buckets' position arrays are not empty.

static std::string dump_index(const mm_idx_t *idx, const std::string &path) {
	FILE *fp = std::fopen(path.c_str(), "wb");
	REQUIRE(fp != nullptr);
	mm_idx_dump(fp, idx);
	std::fclose(fp);
	std::ifstream in(path, std::ios::binary);
	std::string bytes((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
	std::remove(path.c_str());
	return bytes;
}

//! An .mmi's bytes with each bucket's hash entries put in one fixed order.
//! mm_idx_dump writes them in the order of the hash table's slots, which depends
//! on how the table was filled, so two dumps of one index can differ there and
//! nowhere else. The layout is mm_idx_dump's (ext/minimap2/index.c).
static std::string canonical_index(const std::string &bytes) {
	// Thrown rather than REQUIREd: these run once per hash entry.
	const auto check_within = [&](size_t end) {
		if (end > bytes.size()) {
			throw std::runtime_error("the index ends early");
		}
	};
	size_t pos = 0;
	const auto u32_at = [&](size_t at) {
		check_within(at + 4);
		uint32_t value;
		std::memcpy(&value, bytes.data() + at, 4);
		return value;
	};
	const auto take = [&](size_t n) {
		check_within(pos + n);
		const std::string out = bytes.substr(pos, n);
		pos += n;
		return out;
	};
	// Magic, then w, k, b, n_seq and flag.
	const uint32_t b = u32_at(4 + 2 * 4);
	const uint32_t n_seq = u32_at(4 + 3 * 4);
	std::string out = take(4 + 5 * 4);
	for (uint32_t i = 0; i < n_seq; i++) {
		check_within(pos + 1);
		const auto name_length = static_cast<uint8_t>(bytes[pos]);
		out += take(1 + name_length + 4);
	}
	for (uint32_t bucket = 0; bucket < (1u << b); bucket++) {
		const uint32_t n = u32_at(pos);
		out += take(4 + 8 * static_cast<size_t>(n));
		const uint32_t size = u32_at(pos);
		out += take(4);
		std::vector<std::string> entries;
		for (uint32_t e = 0; e < size; e++) {
			entries.push_back(take(16));
		}
		std::sort(entries.begin(), entries.end());
		for (const auto &entry : entries) {
			out += entry;
		}
	}
	return out + bytes.substr(pos);
}

TEST_CASE("A minimap2 index read on several threads is the index read on one", "[Minimap2Aligner]") {
	std::mt19937 rng {20261022};
	Minimap2Config config;
	config.preset = "sr";
	const std::string repeat = random_dna(rng, 20000);
	const std::string ref_a = repeat + random_dna(rng, 150000) + repeat;
	const std::string ref_b = random_dna(rng, 100000) + mutate(rng, repeat, 0.01);

	const std::string path = "data/shards/test_threaded_load.mmi";
	Minimap2Aligner builder(config);
	builder.build_index({{"ref_a", ref_a}, {"ref_b", ref_b}});
	builder.save_index(path);

	Minimap2IndexReader one(path, config, 1);
	Minimap2IndexReader many(path, config, 4);
	auto from_one = one.ReadNextPart();
	auto from_many = many.ReadNextPart();
	REQUIRE(from_one != nullptr);
	REQUIRE(from_many != nullptr);
	// Read to the end of the part, and no further.
	REQUIRE(many.AtEof());
	REQUIRE(from_many->subject_names() == from_one->subject_names());
	const std::string dumped_one = dump_index(from_one->index(), path + ".one");
	const std::string dumped_many = dump_index(from_many->index(), path + ".many");
	REQUIRE(dumped_one.size() > 1000000);
	REQUIRE(dumped_many == dumped_one);
	// Both reads go through the same loader, so each is also held to the file it
	// read, which the builder wrote from its own index without loading anything:
	// equal but for the order of the hash entries, which loading refills.
	std::ifstream saved_in(path, std::ios::binary);
	const std::string saved((std::istreambuf_iterator<char>(saved_in)), std::istreambuf_iterator<char>());
	REQUIRE(canonical_index(dumped_one) == canonical_index(saved));

	std::remove(path.c_str());
}

TEST_CASE("A multi-part minimap2 index read on several threads stops at each part's end", "[Minimap2Aligner]") {
	const std::string part1 = "data/shards/test_threaded_multipart_p1.mmi";
	const std::string part2 = "data/shards/test_threaded_multipart_p2.mmi";
	const std::string multipart = "data/shards/test_threaded_multipart.mmi";
	build_multipart_mmi_fixture(part1, part2, multipart);

	Minimap2Config config;
	config.preset = "sr";
	config.k = 5;
	Minimap2IndexReader reader(multipart, config, 4);
	auto p1 = reader.ReadNextPart();
	REQUIRE(p1 != nullptr);
	REQUIRE(p1->subject_names() == std::vector<std::string> {"part1_ref"});
	REQUIRE_FALSE(reader.AtEof());
	auto p2 = reader.ReadNextPart();
	REQUIRE(p2 != nullptr);
	REQUIRE(p2->subject_names() == std::vector<std::string> {"part2_ref"});
	REQUIRE(reader.AtEof());

	std::remove(part1.c_str());
	std::remove(part2.c_str());
	std::remove(multipart.c_str());
}
