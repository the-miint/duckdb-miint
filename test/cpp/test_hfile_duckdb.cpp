#include <catch2/catch_test_macros.hpp>
#include <cstdio>
#include <string>

// HTSlib headers
#include <htslib-1.22.1/htslib/hfile.h>
#include <htslib-1.22.1/htslib/hts.h>
#include <htslib-1.22.1/htslib/sam.h>
#include <htslib-1.22.1/hfile_internal.h>

#include <SAMReader.hpp>

// ===== FILE*-based hFILE backend for testing =====
// Mirrors the production hFILE_duckdb backend but wraps FILE* instead of FileHandle.
// This validates that SAMReader's hFILE constructor works correctly with any
// custom hFILE backend, without requiring DuckDB linkage in the test binary.

struct hFILE_test {
	hFILE base; // must be first member
	FILE *fp;
};

static ssize_t test_hfile_read(hFILE *fpv, void *buffer, size_t nbytes) {
	auto *fp = reinterpret_cast<hFILE_test *>(fpv);
	return static_cast<ssize_t>(fread(buffer, 1, nbytes, fp->fp));
}

static off_t test_hfile_seek(hFILE *fpv, off_t offset, int whence) {
	auto *fp = reinterpret_cast<hFILE_test *>(fpv);
	if (fseek(fp->fp, offset, whence) != 0) {
		return -1;
	}
	return ftell(fp->fp);
}

static off_t test_hfile_seek_noseek(hFILE * /*fpv*/, off_t /*offset*/, int /*whence*/) {
	errno = ESPIPE;
	return -1;
}

static int test_hfile_close(hFILE *fpv) {
	auto *fp = reinterpret_cast<hFILE_test *>(fpv);
	int ret = fclose(fp->fp);
	fp->fp = nullptr;
	return ret;
}

static const struct hFILE_backend test_backend_seekable = {test_hfile_read, nullptr, test_hfile_seek, nullptr,
                                                           test_hfile_close};

static const struct hFILE_backend test_backend_noseek = {test_hfile_read, nullptr, test_hfile_seek_noseek, nullptr,
                                                         test_hfile_close};

static hFILE *test_hfile_open(const char *path, bool seekable = true) {
	FILE *f = fopen(path, "rb");
	if (!f) {
		return nullptr;
	}

	auto *fp = reinterpret_cast<hFILE_test *>(hfile_init(sizeof(hFILE_test), "r", 0));
	if (!fp) {
		fclose(f);
		return nullptr;
	}

	fp->fp = f;
	fp->base.backend = seekable ? &test_backend_seekable : &test_backend_noseek;
	return &fp->base;
}

// ===== Ground truth helper =====

static int count_records_sam_open(const char *path) {
	auto *fp = sam_open(path, "r");
	if (!fp) {
		return -1;
	}
	auto *hdr = sam_hdr_read(fp);
	if (!hdr) {
		sam_close(fp);
		return -1;
	}
	auto *aln = bam_init1();
	int count = 0;
	while (sam_read1(fp, hdr, aln) >= 0) {
		count++;
	}
	bam_destroy1(aln);
	sam_hdr_destroy(hdr);
	sam_close(fp);
	return count;
}

// ===== Test data =====

static const std::string SAM_WITH_HEADER = "data/sam/foo_has_header.sam";
static const std::string BAM_WITH_TAGS = "data/sam/foo_with_tags.bam";

// ===== Tests: raw hFILE plumbing =====

TEST_CASE("hFILE_duckdb: seekable hFILE reads SAM correctly via hts_hopen", "[hfile_duckdb]") {
	int expected = count_records_sam_open(SAM_WITH_HEADER.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(SAM_WITH_HEADER.c_str(), true);
	REQUIRE(hf != nullptr);

	htsFile *htsfp = hts_hopen(hf, SAM_WITH_HEADER.c_str(), "r");
	REQUIRE(htsfp != nullptr);

	sam_hdr_t *hdr = sam_hdr_read(htsfp);
	REQUIRE(hdr != nullptr);

	bam1_t *aln = bam_init1();
	int count = 0;
	while (sam_read1(htsfp, hdr, aln) >= 0) {
		count++;
	}

	CHECK(count == expected);

	bam_destroy1(aln);
	sam_hdr_destroy(hdr);
	hts_close(htsfp);
}

TEST_CASE("hFILE_duckdb: seekable hFILE reads BAM correctly via hts_hopen", "[hfile_duckdb]") {
	int expected = count_records_sam_open(BAM_WITH_TAGS.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(BAM_WITH_TAGS.c_str(), true);
	REQUIRE(hf != nullptr);

	htsFile *htsfp = hts_hopen(hf, BAM_WITH_TAGS.c_str(), "r");
	REQUIRE(htsfp != nullptr);

	sam_hdr_t *hdr = sam_hdr_read(htsfp);
	REQUIRE(hdr != nullptr);

	bam1_t *aln = bam_init1();
	int count = 0;
	while (sam_read1(htsfp, hdr, aln) >= 0) {
		count++;
	}

	CHECK(count == expected);

	bam_destroy1(aln);
	sam_hdr_destroy(hdr);
	hts_close(htsfp);
}

TEST_CASE("hFILE_duckdb: non-seekable hFILE reads SAM correctly", "[hfile_duckdb]") {
	int expected = count_records_sam_open(SAM_WITH_HEADER.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(SAM_WITH_HEADER.c_str(), false);
	REQUIRE(hf != nullptr);

	htsFile *htsfp = hts_hopen(hf, SAM_WITH_HEADER.c_str(), "r");
	REQUIRE(htsfp != nullptr);

	sam_hdr_t *hdr = sam_hdr_read(htsfp);
	REQUIRE(hdr != nullptr);

	bam1_t *aln = bam_init1();
	int count = 0;
	while (sam_read1(htsfp, hdr, aln) >= 0) {
		count++;
	}

	CHECK(count == expected);

	bam_destroy1(aln);
	sam_hdr_destroy(hdr);
	hts_close(htsfp);
}

TEST_CASE("hFILE_duckdb: non-seekable hFILE reads BAM correctly", "[hfile_duckdb]") {
	int expected = count_records_sam_open(BAM_WITH_TAGS.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(BAM_WITH_TAGS.c_str(), false);
	REQUIRE(hf != nullptr);

	htsFile *htsfp = hts_hopen(hf, BAM_WITH_TAGS.c_str(), "r");
	REQUIRE(htsfp != nullptr);

	sam_hdr_t *hdr = sam_hdr_read(htsfp);
	REQUIRE(hdr != nullptr);

	bam1_t *aln = bam_init1();
	int count = 0;
	while (sam_read1(htsfp, hdr, aln) >= 0) {
		count++;
	}

	CHECK(count == expected);

	bam_destroy1(aln);
	sam_hdr_destroy(hdr);
	hts_close(htsfp);
}

// ===== Tests: SAMReader hFILE constructor =====

TEST_CASE("hFILE_duckdb: SAMReader hFILE constructor reads SAM with header", "[hfile_duckdb]") {
	int expected = count_records_sam_open(SAM_WITH_HEADER.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(SAM_WITH_HEADER.c_str(), true);
	REQUIRE(hf != nullptr);

	miint::SAMReader reader(hf, SAM_WITH_HEADER, false);
	auto batch = reader.read(10000);

	CHECK(static_cast<int>(batch.size()) == expected);
}

TEST_CASE("hFILE_duckdb: SAMReader hFILE constructor reads BAM with header", "[hfile_duckdb]") {
	int expected = count_records_sam_open(BAM_WITH_TAGS.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(BAM_WITH_TAGS.c_str(), true);
	REQUIRE(hf != nullptr);

	miint::SAMReader reader(hf, BAM_WITH_TAGS, false);
	auto batch = reader.read(10000);

	CHECK(static_cast<int>(batch.size()) == expected);
}

TEST_CASE("hFILE_duckdb: SAMReader hFILE constructor reads non-seekable BAM", "[hfile_duckdb]") {
	int expected = count_records_sam_open(BAM_WITH_TAGS.c_str());
	REQUIRE(expected > 0);

	hFILE *hf = test_hfile_open(BAM_WITH_TAGS.c_str(), false);
	REQUIRE(hf != nullptr);

	miint::SAMReader reader(hf, BAM_WITH_TAGS, false);
	auto batch = reader.read(10000);

	CHECK(static_cast<int>(batch.size()) == expected);
}

TEST_CASE("hFILE_duckdb: SAMReader hFILE constructor with include_seq_qual", "[hfile_duckdb]") {
	// ubam_no_sq.sam has actual sequence/quality data on all records
	static const std::string SAM_WITH_SEQ = "data/sam/ubam_no_sq.sam";
	hFILE *hf = test_hfile_open(SAM_WITH_SEQ.c_str(), true);
	REQUIRE(hf != nullptr);

	miint::SAMReader reader(hf, SAM_WITH_SEQ, true, /*require_references=*/false);
	auto batch = reader.read(10000);

	REQUIRE(!batch.empty());
	// With include_seq_qual, sequences and quals should be populated
	CHECK(!batch.sequences.empty());
	CHECK(!batch.quals.empty());
}

TEST_CASE("hFILE_duckdb: SAMReader hFILE constructor with headerless SAM", "[hfile_duckdb]") {
	// Use a headerless SAM file with reference_lengths
	const std::string headerless = "data/sam/foo_no_header.sam";
	int expected = count_records_sam_open(SAM_WITH_HEADER.c_str()); // same data, has header for ground truth
	REQUIRE(expected > 0);

	std::unordered_map<std::string, uint64_t> references = {{"chr1", 10000}, {"chr2", 20000}};

	hFILE *hf = test_hfile_open(headerless.c_str(), true);
	REQUIRE(hf != nullptr);

	miint::SAMReader reader(hf, headerless, references, false);
	auto batch = reader.read(10000);

	CHECK(static_cast<int>(batch.size()) == expected);
}
