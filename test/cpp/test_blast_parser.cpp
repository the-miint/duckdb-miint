#include <catch2/catch_test_macros.hpp>
#include <catch2/catch_approx.hpp>
#include "blast_parser.hpp"

using namespace miint;
using Catch::Approx;

// --- ParseSubmitResponse ---------------------------------------------------

TEST_CASE("BlastParser submit response parsing", "[blast]") {
	SECTION("Extract RID and RTOE from typical NCBI HTML") {
		// NCBI BLAST returns HTML with QBlastInfo embedded in comments.
		// Silent parse failure means the job ID is lost and we can never poll.
		std::string html = R"(
<!DOCTYPE html>
<html><head><title>BLAST Search</title></head>
<body>
<!--QBlastInfoBegin
    RID = ABC123XY
    RTOE = 42
QBlastInfoEnd-->
<p>Request submitted.</p>
</body></html>
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.rid == "ABC123XY");
		CHECK(result.rtoe == 42);
		CHECK(result.error.empty());
	}

	SECTION("Missing RID returns empty rid") {
		// If NCBI changes its format or returns garbage, we need to detect
		// the missing RID rather than silently proceeding with an empty one.
		std::string html = R"(
<!--QBlastInfoBegin
    RTOE = 10
QBlastInfoEnd-->
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.rid.empty());
		CHECK(result.rtoe == 10);
	}

	SECTION("Missing RTOE defaults to 0") {
		std::string html = R"(
<!--QBlastInfoBegin
    RID = XYZ789
QBlastInfoEnd-->
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.rid == "XYZ789");
		CHECK(result.rtoe == 0);
	}

	SECTION("Error in response is captured with exact boundary") {
		// NCBI returns errors inside a specific HTML pattern. If we miss
		// the error, the user gets a cryptic "missing RID" instead of the
		// actual problem (e.g., invalid database, bad query format).
		std::string html = R"(
<p class="error">Error: Failed to read the Blast query: Nucleotide FASTA provided for protein sequence</p>
<!--QBlastInfoBegin
    RID =
    RTOE = 0
QBlastInfoEnd-->
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.error == "Error: Failed to read the Blast query: Nucleotide FASTA provided for protein sequence");
	}

	SECTION("Error tag with extra attributes before class") {
		// The error tag may have other attributes; extraction must find
		// the closing > of the tag, not an intermediate one.
		std::string html = R"(
<div id="content"><p id="err1" class="error">Invalid database name</p></div>
<!--QBlastInfoBegin
    RID =
    RTOE = 0
QBlastInfoEnd-->
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.error == "Invalid database name");
	}
}

// --- ParseStatusResponse ---------------------------------------------------

TEST_CASE("BlastParser status response parsing", "[blast]") {
	SECTION("WAITING status") {
		std::string text = R"(
QBlastInfoBegin
	Status=WAITING
QBlastInfoEnd
)";
		CHECK(BlastParser::ParseStatusResponse(text) == BlastStatus::WAITING);
	}

	SECTION("READY status") {
		std::string text = R"(
QBlastInfoBegin
	Status=READY
QBlastInfoEnd
)";
		CHECK(BlastParser::ParseStatusResponse(text) == BlastStatus::READY);
	}

	SECTION("UNKNOWN status") {
		std::string text = R"(
QBlastInfoBegin
	Status=UNKNOWN
QBlastInfoEnd
)";
		CHECK(BlastParser::ParseStatusResponse(text) == BlastStatus::UNKNOWN);
	}

	SECTION("Garbage input returns UNKNOWN") {
		// Defensive: if NCBI changes format or returns an error page,
		// we must not interpret it as READY (premature retrieval) or
		// hang in WAITING forever.
		CHECK(BlastParser::ParseStatusResponse("") == BlastStatus::UNKNOWN);
		CHECK(BlastParser::ParseStatusResponse("some random text") == BlastStatus::UNKNOWN);
	}

	SECTION("Case-sensitive matching") {
		// NCBI API always emits uppercase status values. Lowercase must
		// not match — treating "waiting" as WAITING could mask a format change.
		CHECK(BlastParser::ParseStatusResponse("Status=waiting") == BlastStatus::UNKNOWN);
		CHECK(BlastParser::ParseStatusResponse("Status=Ready") == BlastStatus::UNKNOWN);
	}

	SECTION("Status key present with empty value") {
		CHECK(BlastParser::ParseStatusResponse("Status=") == BlastStatus::UNKNOWN);
	}
}

// --- ParseTabularResults ---------------------------------------------------

TEST_CASE("BlastParser tabular results parsing", "[blast]") {
	SECTION("Standard multi-hit output") {
		// outfmt 6: qseqid sseqid pident length mismatch gapopen
		//           qstart qend sstart send evalue bitscore
		// Each field must map to the correct struct member — column
		// misalignment silently corrupts downstream analysis.
		std::string tabular =
		    "query1\tNC_001416.1\t99.5\t200\t1\t0\t1\t200\t1000\t1199\t1.5e-10\t350.2\n"
		    "query1\tNC_001422.1\t85.3\t150\t22\t1\t10\t159\t500\t649\t0.001\t120.0\n"
		    "query2\tNC_001416.1\t97.0\t100\t3\t0\t1\t100\t2000\t2099\t5e-50\t500.0\n";

		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 3);

		CHECK(hits[0].query_id == "query1");
		CHECK(hits[0].subject_id == "NC_001416.1");
		CHECK(hits[0].pct_identity == Approx(99.5));
		CHECK(hits[0].alignment_length == 200);
		CHECK(hits[0].mismatches == 1);
		CHECK(hits[0].gap_opens == 0);
		CHECK(hits[0].query_start == 1);
		CHECK(hits[0].query_end == 200);
		CHECK(hits[0].subject_start == 1000);
		CHECK(hits[0].subject_end == 1199);
		CHECK(hits[0].evalue == Approx(1.5e-10));
		CHECK(hits[0].bit_score == Approx(350.2));

		CHECK(hits[1].query_id == "query1");
		CHECK(hits[1].subject_id == "NC_001422.1");
		CHECK(hits[1].pct_identity == Approx(85.3));

		CHECK(hits[2].query_id == "query2");
		CHECK(hits[2].evalue == Approx(5e-50));
	}

	SECTION("Empty input returns empty vector") {
		auto hits = BlastParser::ParseTabularResults("");
		CHECK(hits.empty());
	}

	SECTION("Comment lines are skipped") {
		// BLAST tabular output may include comment headers starting with #.
		// These must not be parsed as data rows.
		std::string tabular =
		    "# BLASTN 2.14.0+\n"
		    "# Query: query1\n"
		    "# Database: nt\n"
		    "# 1 hits found\n"
		    "query1\tNC_001416.1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-40\t200.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].query_id == "query1");
	}

	SECTION("Trailing newline does not produce extra record") {
		std::string tabular =
		    "query1\tNC_001416.1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-40\t200.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		CHECK(hits.size() == 1);
	}

	SECTION("Multiple hits for same query") {
		std::string tabular =
		    "seq1\thitA\t95.0\t50\t2\t1\t1\t50\t1\t50\t1e-5\t100.0\n"
		    "seq1\thitB\t90.0\t50\t5\t0\t1\t50\t100\t149\t1e-3\t80.0\n"
		    "seq1\thitC\t80.0\t40\t8\t0\t5\t44\t200\t239\t0.5\t50.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 3);
		CHECK(hits[0].subject_id == "hitA");
		CHECK(hits[1].subject_id == "hitB");
		CHECK(hits[2].subject_id == "hitC");
	}

	SECTION("Scientific notation evalue with positive exponent") {
		std::string tabular =
		    "q1\ts1\t50.0\t30\t15\t0\t1\t30\t1\t30\t2.5e+02\t20.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].evalue == Approx(250.0));
	}

	SECTION("Evalue of zero") {
		// Extremely significant hits emit evalue=0 (not 1e-300).
		// Must parse as 0.0, not error.
		std::string tabular =
		    "q1\ts1\t100.0\t500\t0\t0\t1\t500\t1\t500\t0\t1000.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].evalue == 0.0);
	}

	SECTION("Windows line endings (CRLF)") {
		// NCBI web API can return \r\n; trailing \r must not corrupt
		// the last field (bitscore) and cause a parse failure.
		std::string tabular =
		    "q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-40\t200.0\r\n"
		    "q2\ts2\t95.0\t80\t4\t0\t1\t80\t1\t80\t1e-20\t150.0\r\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 2);
		CHECK(hits[0].bit_score == Approx(200.0));
		CHECK(hits[1].bit_score == Approx(150.0));
	}

	SECTION("Extra fields beyond 12 are accepted") {
		// NCBI BLAST API returns extra columns for some programs (e.g.,
		// blastp adds "% positives" as field 13). We use the first 12
		// and ignore the rest — dropping these rows would silently lose
		// all blastp results.
		std::string tabular =
		    "q1\ts1\t99.0\t100\t1\t0\t1\t100\t1\t100\t1e-40\t200.0\t95.0\n"
		    "q2\ts2\t90.0\t60\t6\t0\t1\t60\t1\t60\t1e-10\t100.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 2);
		CHECK(hits[0].query_id == "q1");
		CHECK(hits[0].bit_score == Approx(200.0));
		CHECK(hits[1].query_id == "q2");
	}

	SECTION("Lines with fewer than 12 fields are skipped") {
		std::string tabular =
		    "q1\ts1\t95.0\t80\t4\t0\t1\t80\t1\t80\t1e-20\n"
		    "q2\ts2\t90.0\t60\t6\t0\t1\t60\t1\t60\t1e-10\t100.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].query_id == "q2");
	}

	SECTION("Malformed numeric field does not crash") {
		// A corrupted line must be skipped, not crash the caller
		// with an unhandled std::invalid_argument from stod/stoi.
		std::string tabular =
		    "q1\ts1\tNOTANUMBER\t100\t1\t0\t1\t100\t1\t100\t1e-40\t200.0\n"
		    "q2\ts2\t95.0\t80\t4\t0\t1\t80\t1\t80\t1e-20\t150.0\n";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].query_id == "q2");
	}
}

// --- ValidateProgram -------------------------------------------------------

TEST_CASE("BlastParser program validation", "[blast]") {
	SECTION("Valid programs are accepted") {
		CHECK(BlastParser::ValidateProgram("blastn") == true);
		CHECK(BlastParser::ValidateProgram("blastp") == true);
		CHECK(BlastParser::ValidateProgram("blastx") == true);
		CHECK(BlastParser::ValidateProgram("tblastn") == true);
		CHECK(BlastParser::ValidateProgram("tblastx") == true);
	}

	SECTION("Invalid programs are rejected") {
		// "megablast" is a mode flag, not a program name.
		// Case-sensitive: NCBI API is case-sensitive for program names.
		CHECK(BlastParser::ValidateProgram("megablast") == false);
		CHECK(BlastParser::ValidateProgram("BLASTN") == false);
		CHECK(BlastParser::ValidateProgram("") == false);
		CHECK(BlastParser::ValidateProgram("blast") == false);
		CHECK(BlastParser::ValidateProgram("psiblast") == false);
	}
}

// --- BuildFastaPayload -----------------------------------------------------

TEST_CASE("BlastParser FASTA payload construction", "[blast]") {
	SECTION("Single sequence") {
		auto fasta = BlastParser::BuildFastaPayload({"seq1"}, {"ACGTACGT"});
		CHECK(fasta == ">seq1\nACGTACGT\n");
	}

	SECTION("Multiple sequences") {
		auto fasta = BlastParser::BuildFastaPayload(
		    {"seq1", "seq2", "seq3"},
		    {"ACGT", "TGCA", "AAAA"});
		CHECK(fasta == ">seq1\nACGT\n>seq2\nTGCA\n>seq3\nAAAA\n");
	}

	SECTION("Empty input") {
		auto fasta = BlastParser::BuildFastaPayload({}, {});
		CHECK(fasta.empty());
	}

	SECTION("Mismatched ids and sequences throws") {
		// Caller bug — ids/sequences must be parallel arrays. Out-of-bounds
		// access is UB; we throw instead.
		CHECK_THROWS(BlastParser::BuildFastaPayload({"a", "b"}, {"ACGT"}));
		CHECK_THROWS(BlastParser::BuildFastaPayload({"a"}, {"ACGT", "TGCA"}));
	}
}

// --- BuildSubmitBody -------------------------------------------------------

TEST_CASE("BlastParser submit body construction", "[blast]") {
	SECTION("Standard blastn with megablast") {
		// The POST body is application/x-www-form-urlencoded. Each parameter
		// must be present and correctly formatted — a missing CMD=Put means
		// NCBI interprets the request as a status check, not a submission.
		auto body = BlastParser::BuildSubmitBody(
		    "blastn", "nt", ">q1\nACGT\n", 10.0, 500, true);
		CHECK(body.find("CMD=Put") != std::string::npos);
		CHECK(body.find("PROGRAM=blastn") != std::string::npos);
		CHECK(body.find("DATABASE=nt") != std::string::npos);
		CHECK(body.find("EXPECT=10") != std::string::npos);
		CHECK(body.find("HITLIST_SIZE=500") != std::string::npos);
		CHECK(body.find("MEGABLAST=on") != std::string::npos);
		CHECK(body.find("QUERY=") != std::string::npos);
		CHECK(body.find("%3Eq1") != std::string::npos);
	}

	SECTION("blastp without megablast") {
		// MEGABLAST param must be absent for non-blastn programs —
		// NCBI ignores it but including it is misleading and may
		// cause validation errors in the future.
		auto body = BlastParser::BuildSubmitBody(
		    "blastp", "nr", ">q1\nMKWV\n", 0.001, 10, false);
		CHECK(body.find("PROGRAM=blastp") != std::string::npos);
		CHECK(body.find("DATABASE=nr") != std::string::npos);
		CHECK(body.find("EXPECT=0.001") != std::string::npos);
		CHECK(body.find("HITLIST_SIZE=10") != std::string::npos);
		CHECK(body.find("MEGABLAST") == std::string::npos);
	}

	SECTION("Small evalue formatting") {
		// Scientific notation in the EXPECT param must be parseable by NCBI.
		auto body = BlastParser::BuildSubmitBody(
		    "blastn", "nt", ">q1\nACGT\n", 1e-50, 100, false);
		CHECK(body.find("EXPECT=") != std::string::npos);
		// Must not be "EXPECT=0" (which std::to_string would produce for tiny doubles)
		CHECK(body.find("EXPECT=0&") == std::string::npos);
	}

	SECTION("Newlines in FASTA query are preserved") {
		// The FASTA payload contains newlines. In application/x-www-form-urlencoded,
		// newlines should be URL-encoded as %0A (or left raw and handled by the
		// server). We encode them so the body is valid.
		auto body = BlastParser::BuildSubmitBody(
		    "blastn", "nt", ">q1\nACGT\n", 10.0, 500, false);
		// The QUERY value must contain the FASTA content
		auto query_pos = body.find("QUERY=");
		REQUIRE(query_pos != std::string::npos);
		auto query_value = body.substr(query_pos + 6);
		auto amp_pos = query_value.find('&');
		if (amp_pos != std::string::npos) {
			query_value = query_value.substr(0, amp_pos);
		}
		// Newlines should be URL-encoded
		CHECK(query_value.find("%0A") != std::string::npos);
		CHECK(query_value.find("%3Eq1") != std::string::npos);
	}
}

// --- SplitIntoBatches ------------------------------------------------------

TEST_CASE("BlastParser batch splitting", "[blast]") {
	SECTION("All sequences fit in one batch") {
		// Small sequences well under the limit should produce exactly one batch.
		auto batches = BlastParser::SplitIntoBatches(
		    {"s1", "s2", "s3"},
		    {"ACGT", "TGCA", "AAAA"},
		    1000000);
		REQUIRE(batches.size() == 1);
		CHECK(batches[0].ids.size() == 3);
		CHECK(batches[0].sequences.size() == 3);
		CHECK(batches[0].ids[0] == "s1");
		CHECK(batches[0].sequences[2] == "AAAA");
	}

	SECTION("Sequences split across two batches") {
		// Each FASTA record is ">id\nSEQ\n" = id.size() + seq.size() + 3 bytes.
		// With max_bytes=30, two 10-byte sequences fit but three don't.
		// ">s1\nACGTACGTAC\n" = 16 bytes, ">s2\nTGCATGCATG\n" = 16 bytes
		// 16+16 = 32 > 30, so they split into separate batches.
		auto batches = BlastParser::SplitIntoBatches(
		    {"s1", "s2"},
		    {"ACGTACGTAC", "TGCATGCATG"},
		    20);
		REQUIRE(batches.size() == 2);
		CHECK(batches[0].ids.size() == 1);
		CHECK(batches[0].ids[0] == "s1");
		CHECK(batches[1].ids.size() == 1);
		CHECK(batches[1].ids[0] == "s2");
	}

	SECTION("Single oversized sequence goes in its own batch") {
		// A sequence larger than max_bytes must not be dropped — it gets
		// its own batch even though it exceeds the limit. Dropping data
		// silently is unacceptable.
		auto batches = BlastParser::SplitIntoBatches(
		    {"big", "small"},
		    {std::string(5000, 'A'), "ACGT"},
		    1000);
		REQUIRE(batches.size() == 2);
		CHECK(batches[0].ids[0] == "big");
		CHECK(batches[0].ids.size() == 1);
		CHECK(batches[1].ids[0] == "small");
	}

	SECTION("Empty input returns empty batch list") {
		auto batches = BlastParser::SplitIntoBatches({}, {}, 1000);
		CHECK(batches.empty());
	}

	SECTION("Multiple small sequences pack into one batch") {
		// Verify greedy packing: 100 tiny sequences should fit in one batch
		// when the limit is generous.
		std::vector<std::string> ids, seqs;
		for (int i = 0; i < 100; i++) {
			ids.push_back("s" + std::to_string(i));
			seqs.push_back("ACGT");
		}
		auto batches = BlastParser::SplitIntoBatches(ids, seqs, 100000);
		CHECK(batches.size() == 1);
		CHECK(batches[0].ids.size() == 100);
	}
}

// --- Edge cases: realistic NCBI response parsing ----------------------------

TEST_CASE("BlastParser edge cases", "[blast]") {
	SECTION("Submit response with multi-attribute error tag") {
		// NCBI error pages can have complex HTML. The parser must extract
		// content from the class="error" element regardless of surrounding
		// attributes or tag nesting.
		std::string html = R"(
<html><body>
<div class="content">
<ul><li><p class="error">Error: Query is Empty</p></li></ul>
</div>
<!--QBlastInfoBegin
    RID =
    RTOE = 0
QBlastInfoEnd-->
</body></html>
)";
		auto result = BlastParser::ParseSubmitResponse(html);
		CHECK(result.error == "Error: Query is Empty");
		CHECK(result.rid.empty());
	}

	SECTION("Tabular output with only comment headers (no hits)") {
		// A valid BLAST run that finds zero hits emits only comment lines.
		// Must return empty vector, not error.
		std::string tabular = R"(<p><!--
QBlastInfoBegin
	Status=READY
QBlastInfoEnd
--><p>
<PRE>
# blastn
# Iteration: 0
# Query: random_junk
# RID: ABCDEF123
# Database: core_nt
# 0 hits found
</PRE>
)";
		auto hits = BlastParser::ParseTabularResults(tabular);
		CHECK(hits.empty());
	}

	SECTION("Tabular output wrapped in HTML PRE tags") {
		// Real NCBI BLAST responses wrap tabular data in <p>, HTML comments,
		// and <PRE> tags. The parser must handle this without choking.
		std::string tabular = R"(<p><!--
QBlastInfoBegin
	Status=READY
QBlastInfoEnd
--><p>
<PRE>
# blastn
# Query: q1
# Fields: query acc.ver, subject acc.ver, % identity, alignment length, mismatches, gap opens, q. start, q. end, s. start, s. end, evalue, bit score
# 1 hits found
q1	NC_001416.1	100.000	100	0	0	1	100	1000	1099	1e-40	185
</PRE>
)";
		auto hits = BlastParser::ParseTabularResults(tabular);
		REQUIRE(hits.size() == 1);
		CHECK(hits[0].query_id == "q1");
		CHECK(hits[0].pct_identity == Approx(100.0));
	}

	SECTION("UrlEncode RFC 3986 compliance") {
		// Only unreserved characters pass through. Everything else is
		// percent-encoded. This matters for FASTA headers with pipe
		// delimiters and query parameters with special characters.
		CHECK(BlastParser::UrlEncode("abc123") == "abc123");
		CHECK(BlastParser::UrlEncode("a b") == "a%20b");
		CHECK(BlastParser::UrlEncode(">seq|1") == "%3Eseq%7C1");
		CHECK(BlastParser::UrlEncode("\n") == "%0A");
		CHECK(BlastParser::UrlEncode("a=b&c") == "a%3Db%26c");
		CHECK(BlastParser::UrlEncode("") == "");
	}
}
