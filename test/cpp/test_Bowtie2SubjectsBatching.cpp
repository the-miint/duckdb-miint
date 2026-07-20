#include <catch2/catch_test_macros.hpp>

#include <numeric>
#include <string>
#include <vector>

#include "align_bowtie2_daemon_common.hpp"

using duckdb::bt2_daemon::ComputeSubjectBatchRowCounts;

// ComputeSubjectBatchRowCounts decides how to split subjects into Arrow record
// batches so NEITHER the `name` nor the `sequence` STRING column's cumulative
// byte offset buffer crosses INT32_MAX. That overflow is the actual bug: with
// one giant batch, `save_bowtie2_index` on a >2 GB reference threw
// "subjects append failed at row N" the moment the running offset total wrapped
// past 2^31-1. The intent locked down here: any total reference size (each
// contig < cap) splits into byte-bounded batches whose row counts sum back to
// the full subject count — no subject dropped, no batch over the cap unless a
// single oversized row forces it (still safe because contigs are < 2 GB).
// Tiny strings + tiny caps keep the tests allocation-free; the production cap
// is kMaxSubjectBatchBytes.

namespace {

// Build a vector of `count` strings each of the given byte length, filled with
// 'x'. Content is irrelevant to the size-only batching decision.
std::vector<std::string> strings_of_sizes(const std::vector<size_t> &sizes) {
	std::vector<std::string> out;
	out.reserve(sizes.size());
	for (size_t sz : sizes) {
		out.emplace_back(sz, 'x');
	}
	return out;
}

// sum of a row-count vector — the invariant every non-degenerate case checks.
size_t total(const std::vector<size_t> &counts) {
	return std::accumulate(counts.begin(), counts.end(), size_t {0});
}

} // namespace

TEST_CASE("ComputeSubjectBatchRowCounts: empty input yields no batches", "[bowtie2_subjects_batching]") {
	std::vector<std::string> names;
	std::vector<std::string> sequences;
	REQUIRE(ComputeSubjectBatchRowCounts(names, sequences, 100).empty());
}

TEST_CASE("ComputeSubjectBatchRowCounts: everything under the cap is one batch", "[bowtie2_subjects_batching]") {
	auto names = strings_of_sizes({1, 1, 1});
	auto sequences = strings_of_sizes({1, 1, 1});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 100);
	REQUIRE(counts == std::vector<size_t> {3});
	REQUIRE(total(counts) == names.size());
}

TEST_CASE("ComputeSubjectBatchRowCounts: cumulative bytes exactly at the cap stay in one batch",
          "[bowtie2_subjects_batching]") {
	// seq sizes 3 + 4 == 7 == cap. The cap is inclusive, so no split.
	auto names = strings_of_sizes({0, 0});
	auto sequences = strings_of_sizes({3, 4});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 7);
	REQUIRE(counts == std::vector<size_t> {2});
}

TEST_CASE("ComputeSubjectBatchRowCounts: sequence bytes over the cap force a split", "[bowtie2_subjects_batching]") {
	// 4 + 4 == 8 > cap 7, so the second row starts a new batch.
	auto names = strings_of_sizes({0, 0});
	auto sequences = strings_of_sizes({4, 4});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 7);
	REQUIRE(counts == std::vector<size_t> {1, 1});
	REQUIRE(total(counts) == names.size());
}

TEST_CASE("ComputeSubjectBatchRowCounts: multiple rows pack into each batch", "[bowtie2_subjects_batching]") {
	// seq sizes [3,3,3,3], cap 7: 3+3=6 fits, +3=9 overflows → {2,2}.
	auto names = strings_of_sizes({0, 0, 0, 0});
	auto sequences = strings_of_sizes({3, 3, 3, 3});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 7);
	REQUIRE(counts == std::vector<size_t> {2, 2});
	REQUIRE(total(counts) == names.size());
}

TEST_CASE("ComputeSubjectBatchRowCounts: a single row larger than the cap lands alone", "[bowtie2_subjects_batching]") {
	// Row 0's sequence (10) already exceeds cap 7. It must still be emitted, in
	// a batch by itself; row 1 then starts a fresh batch. (Safe in production:
	// contigs are < 2 GB < INT32_MAX, so a lone oversized row still fits int32.)
	auto names = strings_of_sizes({0, 0});
	auto sequences = strings_of_sizes({10, 2});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 7);
	REQUIRE(counts == std::vector<size_t> {1, 1});
	REQUIRE(total(counts) == names.size());
}

TEST_CASE("ComputeSubjectBatchRowCounts: name bytes alone can force a split", "[bowtie2_subjects_batching]") {
	// Sequences are tiny (1 byte) but names are 4 bytes each: 4+4=8 > cap 7, so
	// the name column — not the sequence column — is what triggers the split.
	auto names = strings_of_sizes({4, 4});
	auto sequences = strings_of_sizes({1, 1});
	auto counts = ComputeSubjectBatchRowCounts(names, sequences, 7);
	REQUIRE(counts == std::vector<size_t> {1, 1});
	REQUIRE(total(counts) == names.size());
}
