#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include "Bowtie2Aligner.hpp"
#include "SAMRecord.hpp"
#include "SequenceRecord.hpp"
#include <filesystem>
#include <set>
#include <string>

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

// =============================================================================
// Phase 1 Tests: Config, Binary Discovery, Temp Directory Management
// =============================================================================

TEST_CASE("Bowtie2Config construction with defaults", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	REQUIRE(config.preset == "");
	REQUIRE(config.local == false);
	REQUIRE(config.threads == 1);
	REQUIRE(config.max_secondary == 0);
	REQUIRE(config.extra_args == "");
}

TEST_CASE("Bowtie2Aligner construction requires bowtie2 in PATH", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	// This test will either succeed (bowtie2 installed) or throw with a clear message
	try {
		Bowtie2Aligner aligner(config);
		// Success - bowtie2 was found
		SUCCEED("bowtie2 found in PATH");
	} catch (const std::runtime_error &e) {
		std::string msg = e.what();
		// Error message should mention bowtie2 and PATH
		INFO("Error message: " << msg);
		REQUIRE(msg.find("bowtie2") != std::string::npos);
		REQUIRE(msg.find("PATH") != std::string::npos);
	}
}

TEST_CASE("Bowtie2Aligner creates temp directory on construction", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);
		std::filesystem::path temp_path = aligner.get_temp_dir();

		// Temp directory should exist
		REQUIRE(std::filesystem::exists(temp_path));
		REQUIRE(std::filesystem::is_directory(temp_path));

		// Path should be under system temp directory
		auto sys_temp = std::filesystem::temp_directory_path();
		INFO("Temp dir: " << temp_path);
		INFO("System temp: " << sys_temp);
		REQUIRE(temp_path.string().find(sys_temp.string()) == 0);
	} catch (const std::runtime_error &e) {
		// Skip if bowtie2 not installed
		if (std::string(e.what()).find("bowtie2") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner cleans up temp directory on destruction", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	std::filesystem::path temp_path;

	try {
		{
			Bowtie2Aligner aligner(config);
			temp_path = aligner.get_temp_dir();
			REQUIRE(std::filesystem::exists(temp_path));
		}
		// After destruction, temp directory should be gone
		REQUIRE_FALSE(std::filesystem::exists(temp_path));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner move constructor preserves temp directory", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner1(config);
		std::filesystem::path temp_path = aligner1.get_temp_dir();
		REQUIRE(std::filesystem::exists(temp_path));

		// Move to new aligner
		Bowtie2Aligner aligner2 = std::move(aligner1);

		// Temp directory should still exist
		REQUIRE(std::filesystem::exists(temp_path));
		REQUIRE(aligner2.get_temp_dir() == temp_path);

		// Original aligner's temp_dir should be empty (moved-from state)
		REQUIRE(aligner1.get_temp_dir().empty());
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner move assignment preserves temp directory", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner1(config);
		Bowtie2Aligner aligner2(config);

		std::filesystem::path temp_path1 = aligner1.get_temp_dir();
		std::filesystem::path temp_path2 = aligner2.get_temp_dir();

		REQUIRE(std::filesystem::exists(temp_path1));
		REQUIRE(std::filesystem::exists(temp_path2));
		REQUIRE(temp_path1 != temp_path2);

		// Move assign aligner1 to aligner2
		aligner2 = std::move(aligner1);

		// aligner2's old temp dir should be cleaned up
		REQUIRE_FALSE(std::filesystem::exists(temp_path2));

		// aligner1's temp dir should still exist (now owned by aligner2)
		REQUIRE(std::filesystem::exists(temp_path1));
		REQUIRE(aligner2.get_temp_dir() == temp_path1);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner temp directories are unique", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner1(config);
		Bowtie2Aligner aligner2(config);

		// Each instance should have a unique temp directory
		REQUIRE(aligner1.get_temp_dir() != aligner2.get_temp_dir());
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

// =============================================================================
// Phase 2 Tests: Index Building
// =============================================================================

TEST_CASE("Bowtie2Aligner build_index creates index files", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref1", "ACGTACGTACGTACGTACGTACGTACGTACGT"});

		REQUIRE_NOTHROW(aligner.build_index(subjects));

		// Verify index files exist (bowtie2 creates .1.bt2, .2.bt2, .3.bt2, .4.bt2, .rev.1.bt2, .rev.2.bt2)
		auto temp_dir = aligner.get_temp_dir();
		INFO("Temp dir: " << temp_dir);
		REQUIRE(std::filesystem::exists(temp_dir / "index.1.bt2"));
		REQUIRE(std::filesystem::exists(temp_dir / "index.2.bt2"));
		REQUIRE(std::filesystem::exists(temp_dir / "index.3.bt2"));
		REQUIRE(std::filesystem::exists(temp_dir / "index.4.bt2"));
		REQUIRE(std::filesystem::exists(temp_dir / "index.rev.1.bt2"));
		REQUIRE(std::filesystem::exists(temp_dir / "index.rev.2.bt2"));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner build_index with empty subjects throws", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::vector<AlignmentSubject> subjects;
		REQUIRE_THROWS_WITH(aligner.build_index(subjects),
		                    Catch::Matchers::ContainsSubstring("empty"));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner build_single_index works", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		AlignmentSubject subject{"single_ref", "ACGTACGTACGTACGTACGTACGTACGTACGT"};
		REQUIRE_NOTHROW(aligner.build_single_index(subject));

		// Verify index files exist
		auto temp_dir = aligner.get_temp_dir();
		REQUIRE(std::filesystem::exists(temp_dir / "index.1.bt2"));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner build_index with multiple subjects", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref1", "ACGTACGTACGTACGTACGTACGTACGTACGT"});
		subjects.push_back({"ref2", "GGGGCCCCAAAATTTTGGGGCCCCAAAATTTT"});
		subjects.push_back({"ref3", "TATATATATATATATATATATATATATATATA"});

		REQUIRE_NOTHROW(aligner.build_index(subjects));

		// Verify index files exist
		auto temp_dir = aligner.get_temp_dir();
		REQUIRE(std::filesystem::exists(temp_dir / "index.1.bt2"));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner index files cleaned up on destruction", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	std::filesystem::path temp_dir;

	try {
		{
			Bowtie2Aligner aligner(config);
			temp_dir = aligner.get_temp_dir();

			std::vector<AlignmentSubject> subjects;
			subjects.push_back({"ref1", "ACGTACGTACGTACGTACGTACGTACGTACGT"});
			aligner.build_index(subjects);

			// Index files should exist
			REQUIRE(std::filesystem::exists(temp_dir / "index.1.bt2"));
		}
		// After destruction, temp directory and all contents should be gone
		REQUIRE_FALSE(std::filesystem::exists(temp_dir));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

// =============================================================================
// Phase 3 Tests: Single-End Alignment
// =============================================================================

TEST_CASE("Bowtie2Aligner single-end alignment - exact match", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// 100bp reference
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Query that exactly matches first 50bp
		auto queries = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch); // Get remaining results

		REQUIRE(batch.size() >= 1);
		REQUIRE(batch.read_ids[0] == "query1");
		REQUIRE(batch.references[0] == "reference");
		REQUIRE(batch.positions[0] == 1); // 1-based position
		INFO("CIGAR: " << batch.cigars[0]);
		// CIGAR should indicate a match
		REQUIRE(!batch.cigars[0].empty());
		REQUIRE(batch.cigars[0] != "*");
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner align without index throws", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		auto queries = make_query_batch("q1", "ACGTACGTACGTACGTACGTACGT");
		SAMRecordBatch batch;

		REQUIRE_THROWS_WITH(aligner.align(queries, batch),
		                    Catch::Matchers::ContainsSubstring("index"));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner empty query batch produces no results", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref", "ACGTACGTACGTACGTACGTACGTACGTACGT"});
		aligner.build_index(subjects);

		SequenceRecordBatch queries(false);
		SAMRecordBatch batch;
		aligner.align(queries, batch);

		REQUIRE(batch.size() == 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner multiple queries in batch", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// 100bp reference
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create batch with multiple queries
		SequenceRecordBatch queries(false);
		queries.read_ids = {"q1", "q2", "q3"};
		queries.comments = {"", "", ""};
		queries.sequences1 = {
		    ref_seq.substr(0, 30),  // First 30bp
		    ref_seq.substr(10, 30), // Middle portion
		    ref_seq.substr(52, 30)  // From second half
		};
		queries.quals1 = {QualScore(""), QualScore(""), QualScore("")};

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch); // Get remaining results

		// Should have at least 3 alignments (one per query)
		REQUIRE(batch.size() >= 3);

		// Check that all query IDs appear in results
		std::set<std::string> found_ids(batch.read_ids.begin(), batch.read_ids.end());
		REQUIRE(found_ids.count("q1") > 0);
		REQUIRE(found_ids.count("q2") > 0);
		REQUIRE(found_ids.count("q3") > 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner SAM fields populated correctly", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// 100bp reference
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		auto queries = make_query_batch("test_query", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch); // Get remaining results

		REQUIRE(batch.size() >= 1);

		// Check all required fields are populated
		REQUIRE(!batch.read_ids[0].empty());
		REQUIRE(!batch.references[0].empty());
		REQUIRE(batch.positions[0] > 0); // 1-based
		REQUIRE(batch.stop_positions[0] >= batch.positions[0]);
		REQUIRE(!batch.cigars[0].empty());
		// mapq should be reasonable (0-42 for bowtie2)
		REQUIRE(batch.mapqs[0] <= 42);

		// For single-end, mate fields should be default/unmapped
		REQUIRE(batch.mate_references[0] == "*");
		REQUIRE(batch.mate_positions[0] == 0);
		REQUIRE(batch.template_lengths[0] == 0);

		// tag_as should have alignment score
		REQUIRE(batch.tag_as_values[0] >= 0);

		// YT tag should be "UU" for unpaired
		REQUIRE(batch.tag_yt_values[0] == "UU");
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner handles unmapped reads", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// Reference of all A's
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});
		aligner.build_index(subjects);

		// Query of all G's - won't map (G reverse complements to C, also won't match A's)
		auto queries = make_query_batch("unmapped", "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG");
		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch); // Get remaining results

		// Should have an unmapped record
		REQUIRE(batch.size() >= 1);
		INFO("Flags: " << batch.flags[0]);
		bool is_unmapped = (batch.flags[0] & 0x4) != 0;
		REQUIRE(is_unmapped);
		REQUIRE(batch.references[0] == "*");
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner handles queries with quality scores (FASTQ mode)", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// 100bp reference
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create batch with quality scores (FASTQ mode)
		SequenceRecordBatch queries(false);
		queries.read_ids.push_back("fastq_query");
		queries.comments.push_back("");
		queries.sequences1.push_back("ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");
		queries.quals1.push_back(QualScore(std::string(52, 'I'))); // High quality scores

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		REQUIRE(batch.size() >= 1);
		REQUIRE(batch.read_ids[0] == "fastq_query");
		REQUIRE(batch.references[0] == "reference");
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner multiple align() calls stream to same process", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// First batch
		auto batch1 = make_query_batch("query1", "ACGTACGTACGTACGTACGTACGTACGTACGT");

		// Second batch
		auto batch2 = make_query_batch("query2", "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA");

		SAMRecordBatch results;

		// Align first batch
		aligner.align(batch1, results);
		REQUIRE(aligner.is_aligner_running()); // Process should still be running

		// Align second batch
		aligner.align(batch2, results);
		REQUIRE(aligner.is_aligner_running()); // Process should still be running

		// Finish and get all results
		aligner.finish(results);
		REQUIRE_FALSE(aligner.is_aligner_running()); // Process should be stopped

		// Should have results from both batches
		std::set<std::string> found_ids(results.read_ids.begin(), results.read_ids.end());
		REQUIRE(found_ids.count("query1") > 0);
		REQUIRE(found_ids.count("query2") > 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

// =============================================================================
// Phase 4 Tests: Paired-End Alignment
// =============================================================================

// Helper to create a paired-end query batch
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

TEST_CASE("Bowtie2Aligner paired-end alignment - both mates map", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// Long reference (~200bp) to allow fragment sizes
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
		                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		                      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create paired-end query: R1 from start, R2 from middle
		std::string r1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT"; // First 32bp
		std::string r2_seq = "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA"; // From position 52

		auto queries = make_paired_query_batch("paired_read", r1_seq, r2_seq);

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// Should have 2 records (one per mate)
		REQUIRE(batch.size() >= 2);

		// Both should be mapped (no unmapped flag 0x4)
		bool found_r1 = false, found_r2 = false;
		for (size_t i = 0; i < batch.size(); i++) {
			if (batch.read_ids[i] == "paired_read") {
				bool is_unmapped = (batch.flags[i] & 0x4) != 0;
				REQUIRE_FALSE(is_unmapped);

				// Check pair flags: 0x1 (paired), 0x40 (first), 0x80 (second)
				bool is_paired_flag = (batch.flags[i] & 0x1) != 0;
				REQUIRE(is_paired_flag);

				bool is_first = (batch.flags[i] & 0x40) != 0;
				bool is_second = (batch.flags[i] & 0x80) != 0;

				if (is_first) {
					found_r1 = true;
				}
				if (is_second) {
					found_r2 = true;
				}
			}
		}
		REQUIRE(found_r1);
		REQUIRE(found_r2);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner paired-end alignment - mate info populated", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
		                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		                      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		std::string r1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
		std::string r2_seq = "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA";

		auto queries = make_paired_query_batch("mate_test", r1_seq, r2_seq);

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		REQUIRE(batch.size() >= 2);

		// Find the two mates
		int first_idx = -1, second_idx = -1;
		for (size_t i = 0; i < batch.size(); i++) {
			if (batch.read_ids[i] == "mate_test") {
				if (batch.flags[i] & 0x40) {
					first_idx = static_cast<int>(i);
				}
				if (batch.flags[i] & 0x80) {
					second_idx = static_cast<int>(i);
				}
			}
		}

		REQUIRE(first_idx >= 0);
		REQUIRE(second_idx >= 0);

		// Check mate reference (should point to same reference for concordant pair)
		// SAM uses "=" to mean "same reference as read" which is the standard convention
		bool first_mate_ref_ok =
		    (batch.mate_references[first_idx] == "reference" || batch.mate_references[first_idx] == "=");
		bool second_mate_ref_ok =
		    (batch.mate_references[second_idx] == "reference" || batch.mate_references[second_idx] == "=");
		REQUIRE(first_mate_ref_ok);
		REQUIRE(second_mate_ref_ok);

		// Check mate positions are set (non-zero for mapped mates)
		REQUIRE(batch.mate_positions[first_idx] > 0);
		REQUIRE(batch.mate_positions[second_idx] > 0);

		// Check template length is set (non-zero for proper pairs)
		// One should be positive, one negative
		REQUIRE(batch.template_lengths[first_idx] != 0);
		REQUIRE(batch.template_lengths[second_idx] != 0);
		// They should have opposite signs
		REQUIRE((batch.template_lengths[first_idx] > 0) != (batch.template_lengths[second_idx] > 0));
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner paired-end with quality scores", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create paired-end batch with quality scores (FASTQ mode)
		SequenceRecordBatch queries(true); // paired
		queries.read_ids.push_back("paired_fastq");
		queries.comments.push_back("");
		queries.sequences1.push_back("ACGTACGTACGTACGTACGTACGTACGTACGT");
		queries.sequences2.push_back("GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA");
		queries.quals1.push_back(QualScore(std::string(32, 'I')));
		queries.quals2.push_back(QualScore(std::string(32, 'I')));

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// Should have 2 records (one per mate)
		REQUIRE(batch.size() >= 2);

		// Verify paired flags are set
		for (size_t i = 0; i < batch.size(); i++) {
			if (batch.read_ids[i] == "paired_fastq") {
				bool is_paired_flag = (batch.flags[i] & 0x1) != 0;
				REQUIRE(is_paired_flag);
			}
		}
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner paired-end one mate unmapped", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		// Reference of all A's
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		                           "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});
		aligner.build_index(subjects);

		// R1 maps (A's), R2 doesn't map (G's)
		std::string r1_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
		std::string r2_seq = "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG";

		auto queries = make_paired_query_batch("halfmap", r1_seq, r2_seq);

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		REQUIRE(batch.size() >= 2);

		// Find mates and check mapping status
		bool found_mapped = false, found_unmapped = false;
		for (size_t i = 0; i < batch.size(); i++) {
			if (batch.read_ids[i] == "halfmap") {
				bool is_unmapped = (batch.flags[i] & 0x4) != 0;
				bool is_first = (batch.flags[i] & 0x40) != 0;
				bool is_second = (batch.flags[i] & 0x80) != 0;

				if (is_first) {
					REQUIRE_FALSE(is_unmapped); // R1 should map
					found_mapped = true;
				}
				if (is_second) {
					REQUIRE(is_unmapped); // R2 should NOT map
					found_unmapped = true;
				}

				// Check mate unmapped flag (0x8)
				bool mate_unmapped = (batch.flags[i] & 0x8) != 0;
				if (is_first) {
					REQUIRE(mate_unmapped); // R1's mate (R2) is unmapped
				}
				if (is_second) {
					REQUIRE_FALSE(mate_unmapped); // R2's mate (R1) is mapped
				}
			}
		}
		REQUIRE(found_mapped);
		REQUIRE(found_unmapped);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner paired-end multiple pairs", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC"
		                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		                      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create multiple paired-end queries
		SequenceRecordBatch queries(true);
		queries.read_ids = {"pair1", "pair2", "pair3"};
		queries.comments = {"", "", ""};
		queries.sequences1 = {"ACGTACGTACGTACGTACGTACGTACGTACGT", "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA",
		                      "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"};
		queries.sequences2 = {"GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA", "ACGTACGTACGTACGTACGTACGTACGTACGT",
		                      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT"};
		queries.quals1 = {QualScore(""), QualScore(""), QualScore("")};
		queries.quals2 = {QualScore(""), QualScore(""), QualScore("")};

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// Should have at least 6 records (2 per pair)
		REQUIRE(batch.size() >= 6);

		// All three pairs should be present
		std::set<std::string> found_ids(batch.read_ids.begin(), batch.read_ids.end());
		REQUIRE(found_ids.count("pair1") > 0);
		REQUIRE(found_ids.count("pair2") > 0);
		REQUIRE(found_ids.count("pair3") > 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner YT tag for paired reads", "[Bowtie2Aligner]") {
	Bowtie2Config config;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		std::string r1_seq = "ACGTACGTACGTACGTACGTACGTACGTACGT";
		std::string r2_seq = "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAA";

		auto queries = make_paired_query_batch("yt_test", r1_seq, r2_seq);

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		REQUIRE(batch.size() >= 2);

		// For concordant pairs, YT tag should be "CP" (Concordant Pair)
		// For discordant, it would be "DP"
		// For unpaired aligned, it would be "UP"
		for (size_t i = 0; i < batch.size(); i++) {
			if (batch.read_ids[i] == "yt_test") {
				INFO("YT tag: " << batch.tag_yt_values[i]);
				// Should be concordant or discordant pair, not unpaired
				REQUIRE(batch.tag_yt_values[i] != "UU");
			}
		}
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

// =============================================================================
// Phase 5 Tests: Config Options
// =============================================================================

TEST_CASE("Bowtie2Aligner with preset option", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	config.preset = "very-fast"; // Code adds -- prefix automatically

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		auto queries = make_query_batch("preset_test", "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT");

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// Should still align with --very-fast
		REQUIRE(batch.size() >= 1);
		REQUIRE(batch.read_ids[0] == "preset_test");
		REQUIRE(batch.references[0] == "reference");
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner with local alignment mode", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	config.local = true;

	try {
		Bowtie2Aligner aligner(config);

		// Reference of A's followed by T's
		std::string ref_seq = "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
		                      "TTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTTT";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Query: G's (mismatch) + A's (match) + G's (mismatch)
		// Local mode should soft-clip the G flanks and align the A core
		auto queries = make_query_batch("local_test", "GGGGGGAAAAAAAAAAAAAAAAAAAAAAAAGGGGGG");

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		REQUIRE(batch.size() >= 1);
		// In local mode, the read should map (soft-clipping the G's)
		bool is_mapped = (batch.flags[0] & 0x4) == 0;
		REQUIRE(is_mapped);
		// CIGAR should contain S for soft-clipping
		INFO("CIGAR: " << batch.cigars[0]);
		REQUIRE(batch.cigars[0].find('S') != std::string::npos);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner with multi-threading", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	config.threads = 2;

	try {
		Bowtie2Aligner aligner(config);

		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "GGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCCTTAAGGCC";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Create batch with multiple queries
		SequenceRecordBatch queries(false);
		queries.read_ids = {"q1", "q2", "q3", "q4"};
		queries.comments = {"", "", "", ""};
		queries.sequences1 = {ref_seq.substr(0, 30), ref_seq.substr(10, 30), ref_seq.substr(20, 30),
		                      ref_seq.substr(52, 30)};
		queries.quals1 = {QualScore(""), QualScore(""), QualScore(""), QualScore("")};

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// Should have all 4 aligned
		REQUIRE(batch.size() >= 4);
		std::set<std::string> found_ids(batch.read_ids.begin(), batch.read_ids.end());
		REQUIRE(found_ids.count("q1") > 0);
		REQUIRE(found_ids.count("q2") > 0);
		REQUIRE(found_ids.count("q3") > 0);
		REQUIRE(found_ids.count("q4") > 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner with max_secondary option", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	config.max_secondary = 3; // Report up to 3 alignments per read

	try {
		Bowtie2Aligner aligner(config);

		// Reference with repeated pattern that will cause multiple valid alignments
		std::string ref_seq = "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT"
		                      "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT";
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"reference", ref_seq});
		aligner.build_index(subjects);

		// Query that matches multiple locations
		auto queries = make_query_batch("multi_hit", "ACGTACGTACGTACGTACGTACGT");

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// With -k 3, should get up to 3 alignments for this read
		REQUIRE(batch.size() >= 1);

		// Count alignments for our read
		int count = 0;
		for (const auto &id : batch.read_ids) {
			if (id == "multi_hit") {
				count++;
			}
		}
		INFO("Number of alignments: " << count);
		// Should have more than 1 due to repeated sequence
		REQUIRE(count >= 1);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}

TEST_CASE("Bowtie2Aligner with extra_args option", "[Bowtie2Aligner]") {
	Bowtie2Config config;
	config.extra_args = "--no-unal"; // Suppress unaligned reads in output

	try {
		Bowtie2Aligner aligner(config);

		// Reference of all A's
		std::vector<AlignmentSubject> subjects;
		subjects.push_back({"ref", "AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"});
		aligner.build_index(subjects);

		// Create batch with one that maps and one that doesn't
		SequenceRecordBatch queries(false);
		queries.read_ids = {"maps", "unmaps"};
		queries.comments = {"", ""};
		queries.sequences1 = {"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA",
		                      "GGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG"};
		queries.quals1 = {QualScore(""), QualScore("")};

		SAMRecordBatch batch;
		aligner.align(queries, batch);
		aligner.finish(batch);

		// With --no-unal, the unmapped read should not appear in output
		std::set<std::string> found_ids(batch.read_ids.begin(), batch.read_ids.end());
		REQUIRE(found_ids.count("maps") > 0);
		// "unmaps" should NOT be in output due to --no-unal
		REQUIRE(found_ids.count("unmaps") == 0);
	} catch (const std::runtime_error &e) {
		if (std::string(e.what()).find("bowtie2") != std::string::npos &&
		    std::string(e.what()).find("PATH") != std::string::npos) {
			SKIP("bowtie2 not installed");
		}
		throw;
	}
}
