#ifdef MIINT_HAS_SYLPH

#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>

#include "SylphDatabase.hpp"
#include "sylph.h"

#include <stdexcept>
#include <string>

using Catch::Approx;

namespace {

// Path is relative to the repo root (where the tests binary is executed
// from, matching Minimap2/SortMeRNA tests).
const std::string kTinySyldb = "data/sylph/tiny.syldb";

} // namespace

TEST_CASE("SylphDatabaseHandle loads a real .syldb", "[sylph][database]") {
	miint::SylphDatabaseHandle db(kTinySyldb);
	REQUIRE(db.num_genomes() == 3);
	REQUIRE(db.raw() != nullptr);
}

TEST_CASE("SylphDatabaseHandle throws on missing path", "[sylph][database]") {
	REQUIRE_THROWS_AS(miint::SylphDatabaseHandle("/tmp/__sylph_no_such_file__.syldb"), std::runtime_error);
	// The thrown message should mention the underlying FFI error.
	try {
		miint::SylphDatabaseHandle bad("/tmp/__sylph_no_such_file__.syldb");
	} catch (const std::runtime_error &e) {
		std::string what = e.what();
		REQUIRE(what.find("failed to open") != std::string::npos);
	}
}

TEST_CASE("SylphDatabaseHandle survives multiple instances", "[sylph][database]") {
	// Construct + destruct two independent handles back-to-back to make
	// sure no thread-local error or global state is clobbered.
	{
		miint::SylphDatabaseHandle a(kTinySyldb);
		REQUIRE(a.num_genomes() == 3);
	}
	{
		miint::SylphDatabaseHandle b(kTinySyldb);
		REQUIRE(b.num_genomes() == 3);
	}
}

TEST_CASE("sylph version strings are non-empty", "[sylph][database]") {
	const char *v = sylph_version();
	REQUIRE(v != nullptr);
	REQUIRE(std::string(v).length() > 0);

	const char *fv = sylph_miint_fork_version();
	REQUIRE(fv != nullptr);
	REQUIRE(std::string(fv).find("miint") != std::string::npos);
}

// Pin the source-of-truth refactor: BuildProfileParams seeds Data.profile_params
// from this FFI getter, and Bind() only overrides fields the user explicitly
// passed via named parameters. If anyone reverts DEREP_PROFILE_ANI to 0.975 or
// silently drops the FFI export, this fails immediately at the unit-test layer
// — the SQL tests would only catch it through downstream behavioral drift on
// borderline genomes.
TEST_CASE("sylph_profile_params_default seeds sylph defaults", "[sylph][ffi]") {
	SylphProfileParams pp {};
	REQUIRE(sylph_profile_params_default(&pp) == 0);
	REQUIRE(pp.pseudotax == 1);
	REQUIRE(pp.estimator == 0);
	REQUIRE(pp.estimate_unknown == 0);
	REQUIRE(pp.minimum_ani == Approx(-1.0));   // sentinel: use sylph's MIN_ANI_P_DEF
	REQUIRE(pp.seq_id == Approx(-1.0));        // sentinel: auto
	REQUIRE(pp.redundant_ani == Approx(99.0)); // DEREP_PROFILE_ANI fix
	REQUIRE(pp.min_count_correct == Approx(3.0));
	REQUIRE(pp.min_number_kmers == Approx(50.0));
	REQUIRE(pp.num_threads == 0);
	REQUIRE(sylph_profile_params_default(nullptr) != 0);
}

// Same regression coverage for the sketch-side defaults. The dedup_fpr field
// is the second source-of-truth fix: SylphSketchParams::Default now uses
// DEFAULT_FPR (0.0001) so the FFI getter agrees with sylph CLI's --fpr
// default. Earlier value was 0.0 (exact dedup), which silently changed the
// algorithm class for any caller that hadn't explicitly overridden it.
TEST_CASE("sylph_sketch_params_default seeds sylph defaults", "[sylph][ffi]") {
	SylphSketchParams sp {};
	REQUIRE(sylph_sketch_params_default(&sp) == 0);
	REQUIRE(sp.dedup == 1);
	REQUIRE(sp.dedup_fpr == Approx(0.0001)); // DEFAULT_FPR alignment fix
	REQUIRE(sylph_sketch_params_default(nullptr) != 0);
}

#endif // MIINT_HAS_SYLPH
