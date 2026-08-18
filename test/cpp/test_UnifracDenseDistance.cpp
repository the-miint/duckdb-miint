#include <catch2/catch_test_macros.hpp>

#include "unifrac_dense_distance.hpp"

#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

using miint::unifrac::BuildDenseDistanceMatrix;
using miint::unifrac::DenseDistanceMatrixBuilder;
using miint::unifrac::DistanceEntry;

namespace {

// Canonical 3-sample id vector for the fixtures below.
const std::vector<std::string> kIds3 = {"A", "B", "C"};

} // namespace

TEST_CASE("BuildDenseDistanceMatrix fills a symmetric zero-diagonal matrix", "[unifrac][dense_distance]") {
	// Full upper triangle of a 3-sample matrix: (A,B), (A,C), (B,C).
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}};
	auto m = BuildDenseDistanceMatrix(entries, 3, kIds3);
	REQUIRE(m.size() == 9);
	// Diagonal is definitionally zero.
	REQUIRE(m[0 * 3 + 0] == 0.0f);
	REQUIRE(m[1 * 3 + 1] == 0.0f);
	REQUIRE(m[2 * 3 + 2] == 0.0f);
	// Symmetric off-diagonal.
	REQUIRE(m[0 * 3 + 1] == 0.5f);
	REQUIRE(m[1 * 3 + 0] == 0.5f);
	REQUIRE(m[0 * 3 + 2] == 0.25f);
	REQUIRE(m[2 * 3 + 0] == 0.25f);
	REQUIRE(m[1 * 3 + 2] == 0.75f);
	REQUIRE(m[2 * 3 + 1] == 0.75f);
}

TEST_CASE("BuildDenseDistanceMatrix accepts a distance of exactly zero", "[unifrac][dense_distance]") {
	// A genuine off-diagonal zero (identical composition) must be distinguished
	// from an unfilled cell — the completeness check must still pass.
	std::vector<DistanceEntry> entries = {{0, 1, 0.0}, {0, 2, 0.25}, {1, 2, 0.75}};
	auto m = BuildDenseDistanceMatrix(entries, 3, kIds3);
	REQUIRE(m[0 * 3 + 1] == 0.0f);
	REQUIRE(m[1 * 3 + 0] == 0.0f);
}

TEST_CASE("BuildDenseDistanceMatrix ignores self-pairs", "[unifrac][dense_distance]") {
	// (A,A) is the diagonal — definitionally zero — and must not count toward
	// completeness nor overwrite anything.
	std::vector<DistanceEntry> entries = {{0, 0, 0.0}, {0, 1, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}, {2, 2, 0.0}};
	auto m = BuildDenseDistanceMatrix(entries, 3, kIds3);
	REQUIRE(m[0 * 3 + 1] == 0.5f);
	REQUIRE(m[0 * 3 + 0] == 0.0f);
}

TEST_CASE("BuildDenseDistanceMatrix tolerates identical duplicates and (b,a) mirrors", "[unifrac][dense_distance]") {
	// The same unordered pair given twice with the same value — including the
	// reversed (b,a) orientation — is a no-op, not a conflict.
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {1, 0, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}, {2, 1, 0.75}};
	auto m = BuildDenseDistanceMatrix(entries, 3, kIds3);
	REQUIRE(m[0 * 3 + 1] == 0.5f);
	REQUIRE(m[1 * 3 + 2] == 0.75f);
}

TEST_CASE("BuildDenseDistanceMatrix rejects conflicting duplicates", "[unifrac][dense_distance]") {
	// Same unordered pair, two different values → the input is ambiguous and we
	// must fail loud rather than silently keep one.
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {0, 1, 0.6}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix rejects a conflicting (b,a) mirror", "[unifrac][dense_distance]") {
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {1, 0, 0.6}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix rejects negative distances", "[unifrac][dense_distance]") {
	std::vector<DistanceEntry> entries = {{0, 1, -0.1}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix names the FIRST gap in row-major order", "[unifrac][dense_distance]") {
	// n=4 with three gaps — (A,B), (A,C), (A,D) all missing — and three present
	// pairs (B,C), (B,D), (C,D). The contract is to name the first missing pair
	// in row-major upper-triangle order, which is (A,B). A weaker test with a
	// single gap could not tell "finds the first gap" from "finds any gap", so a
	// scan-order regression (column-major, unordered map) would slip through.
	const std::vector<std::string> ids4 = {"A", "B", "C", "D"};
	std::vector<DistanceEntry> entries = {{1, 2, 0.5}, {1, 3, 0.25}, {2, 3, 0.75}};
	try {
		BuildDenseDistanceMatrix(entries, 4, ids4);
		FAIL("expected std::invalid_argument");
	} catch (const std::invalid_argument &e) {
		const std::string msg = e.what();
		REQUIRE(msg.find("'A'") != std::string::npos);
		REQUIRE(msg.find("'B'") != std::string::npos);
		// The first gap is (A,B); C and D must not appear in the message.
		REQUIRE(msg.find("'C'") == std::string::npos);
		REQUIRE(msg.find("'D'") == std::string::npos);
	}
}

TEST_CASE("BuildDenseDistanceMatrix rejects non-finite distances", "[unifrac][dense_distance]") {
	// NaN/Inf would silently poison skbb_pcoa_fsvd_fp32 / skbb_permanova_fp32,
	// and NaN specifically breaks the identical-duplicate check (NaN != NaN), so
	// it must be rejected up front rather than stored.
	const double nan_v = std::numeric_limits<double>::quiet_NaN();
	const double inf_v = std::numeric_limits<double>::infinity();
	std::vector<DistanceEntry> with_nan = {{0, 1, nan_v}, {0, 2, 0.25}, {1, 2, 0.75}};
	std::vector<DistanceEntry> with_inf = {{0, 1, inf_v}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(with_nan, 3, kIds3), std::invalid_argument);
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(with_inf, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix rejects a nonzero self-distance", "[unifrac][dense_distance]") {
	// A zero self-pair is the (ignored) diagonal; a nonzero one is contradictory
	// input (the diagonal is definitionally zero) and must fail loud, not be
	// silently discarded.
	std::vector<DistanceEntry> entries = {{0, 0, 3.7}, {0, 1, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix rejects out-of-range indices", "[unifrac][dense_distance]") {
	// An index >= n must throw cleanly rather than read/write out of bounds —
	// this header is DuckDB-free and callable directly.
	std::vector<DistanceEntry> entries = {{0, 3, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, kIds3), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix rejects an ids/n size mismatch", "[unifrac][dense_distance]") {
	// ids must describe exactly n samples; a mismatch is a caller-contract bug
	// that would otherwise index ids out of bounds in an error message.
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}};
	const std::vector<std::string> ids2 = {"A", "B"};
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(entries, 3, ids2), std::invalid_argument);
}

TEST_CASE("BuildDenseDistanceMatrix requires at least two samples", "[unifrac][dense_distance]") {
	// A 0- or 1-sample matrix has no unordered pairs; ordination and PERMANOVA
	// are both undefined, so this is an error rather than an empty result.
	std::vector<DistanceEntry> none;
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(none, 0, {}), std::invalid_argument);
	REQUIRE_THROWS_AS(BuildDenseDistanceMatrix(none, 1, {"A"}), std::invalid_argument);
}

// ── DenseDistanceMatrixBuilder — the incremental form ────────────────────────
//
// WHY this exists: the batch entry point needs the caller to hold every entry in
// memory before the fill starts, which is O(N^2) heap on top of the matrix — the
// dominant cost for a condensed COO of a large distance table. The builder lets a
// streaming reader write cells as rows arrive and keep nothing per-row. The two
// forms MUST agree exactly (the batch function is a wrapper over the builder), so
// these tests pin equivalence and that every validation fires at the same point
// in the lifecycle.

TEST_CASE("DenseDistanceMatrixBuilder matches the batch result cell for cell", "[unifrac][dense_distance]") {
	// The equivalence that lets ReadDistanceTable stream without changing output.
	std::vector<DistanceEntry> entries = {{0, 1, 0.5}, {0, 2, 0.25}, {1, 2, 0.75}};
	const auto batch = BuildDenseDistanceMatrix(entries, 3, kIds3);

	DenseDistanceMatrixBuilder builder(3, kIds3);
	for (const auto &e : entries) {
		builder.Add(e.i, e.j, e.distance);
	}
	const auto streamed = builder.Finish();

	REQUIRE(streamed.size() == batch.size());
	for (size_t k = 0; k < batch.size(); ++k) {
		REQUIRE(streamed[k] == batch[k]);
	}
}

TEST_CASE("DenseDistanceMatrixBuilder tolerates duplicates and mirrors, ignores self-pairs",
          "[unifrac][dense_distance]") {
	// Same tolerance rules as the batch form: a repeated unordered pair with the
	// same value is a no-op, and a zero self-pair is the (ignored) diagonal.
	DenseDistanceMatrixBuilder builder(3, kIds3);
	builder.Add(0, 1, 0.5);
	builder.Add(1, 0, 0.5); // mirror, identical value
	builder.Add(0, 0, 0.0); // diagonal, ignored
	builder.Add(0, 2, 0.25);
	builder.Add(1, 2, 0.75);
	const auto m = builder.Finish();
	REQUIRE(m[0 * 3 + 1] == 0.5f);
	REQUIRE(m[1 * 3 + 0] == 0.5f);
	REQUIRE(m[0 * 3 + 0] == 0.0f);
}

TEST_CASE("DenseDistanceMatrixBuilder rejects bad values at Add, not at Finish", "[unifrac][dense_distance]") {
	// A streaming reader must fail on the offending ROW, while the row is still
	// the thing being processed — deferring to Finish would lose which input was
	// at fault. Each of these is the per-entry validation the batch form applies.
	const double nan_v = std::numeric_limits<double>::quiet_NaN();
	const double inf_v = std::numeric_limits<double>::infinity();
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		REQUIRE_THROWS_AS(b.Add(0, 3, 0.5), std::invalid_argument); // index >= n
	}
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		REQUIRE_THROWS_AS(b.Add(0, 1, nan_v), std::invalid_argument);
	}
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		REQUIRE_THROWS_AS(b.Add(0, 1, inf_v), std::invalid_argument);
	}
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		REQUIRE_THROWS_AS(b.Add(0, 1, -0.5), std::invalid_argument);
	}
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		REQUIRE_THROWS_AS(b.Add(0, 0, 3.7), std::invalid_argument); // nonzero self-distance
	}
	{
		DenseDistanceMatrixBuilder b(3, kIds3);
		b.Add(0, 1, 0.5);
		REQUIRE_THROWS_AS(b.Add(1, 0, 0.9), std::invalid_argument); // conflicting mirror
	}
}

TEST_CASE("DenseDistanceMatrixBuilder reports an incomplete matrix at Finish", "[unifrac][dense_distance]") {
	// Completeness is only knowable once the stream ends, so unlike the value
	// checks it must surface from Finish — and must name the first missing pair in
	// row-major order, exactly as the batch form does.
	DenseDistanceMatrixBuilder builder(3, kIds3);
	builder.Add(0, 1, 0.5);
	builder.Add(1, 2, 0.75); // (A,C) never provided
	try {
		builder.Finish();
		FAIL("expected an incomplete-matrix throw");
	} catch (const std::invalid_argument &e) {
		const std::string msg = e.what();
		REQUIRE(msg.find("incomplete distance matrix") != std::string::npos);
		REQUIRE(msg.find("'A'") != std::string::npos);
		REQUIRE(msg.find("'C'") != std::string::npos);
	}
}

TEST_CASE("DenseDistanceMatrixBuilder validates its contract at construction", "[unifrac][dense_distance]") {
	// n < 2 and an ids/n mismatch are caller-contract bugs; catching them at
	// construction means the N^2 allocation never happens on a bad contract.
	REQUIRE_THROWS_AS(DenseDistanceMatrixBuilder(1, {"A"}), std::invalid_argument);
	REQUIRE_THROWS_AS(DenseDistanceMatrixBuilder(3, {"A", "B"}), std::invalid_argument);
}
