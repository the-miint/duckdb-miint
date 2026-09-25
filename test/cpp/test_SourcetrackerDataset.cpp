#include <catch2/catch_approx.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include <cmath>
#include <limits>
#include <set>
#include <string>
#include <utility>
#include <vector>

#include "sourcetracker_dataset.hpp"

using Catch::Approx;
using Catch::Matchers::ContainsSubstring;
using miint::sourcetracker::CheckDepths;
using miint::sourcetracker::Dataset;
using miint::sourcetracker::IngestDataset;
using miint::unifrac::CooRow;
using miint::unifrac::MetadataRow;

// ---------------------------------------------------------------------------
// The dataset layer is the one place where a wrong answer is silent. st3 takes
// feature and sample INDICES; if a cell lands on the wrong sample, or a sample's
// role or environment is attached to its neighbour, the sampler still runs and
// still returns well-formed proportions -- for the wrong sink. So every
// assertion here is against independently known values, and every precondition
// is asserted by the MESSAGE that names the offending id, because that message
// is what the user sees when their metadata and table disagree.
// ---------------------------------------------------------------------------
namespace {

// Wide metadata (sample_id, source_sink, env) as the unpivoted long-form rows
// ReadWideMetadata produces: one row per (sample, variable). A NULL env arrives
// as "" -- that is what the reader does, and what the empty-env check must see.
void AddMeta(std::vector<MetadataRow> &rows, const std::string &sample, const std::string &role,
             const std::string &env) {
	rows.push_back({sample, "source_sink", role});
	rows.push_back({sample, "env", env});
}

// Sinks s_b, s_a; sources src_z (sewage), src_m and src_m2 (seawater); features
// f3, f1, f2. Every id list is deliberately out of lexicographic order, and the
// cells are emitted scrambled, so a first-seen-order dictionary is caught.
std::vector<MetadataRow> ToyMetadata() {
	std::vector<MetadataRow> m;
	AddMeta(m, "src_z", "source", "sewage");
	AddMeta(m, "s_b", "sink", "");
	AddMeta(m, "src_m2", "source", "seawater");
	AddMeta(m, "s_a", "sink", "ignored-for-sinks");
	AddMeta(m, "src_m", "source", "seawater");
	return m;
}

std::vector<CooRow> ToyTable() {
	return {
	    {"src_m2", "f2", 4.0}, {"s_b", "f3", 1.0},    {"src_z", "f1", 10.0}, {"s_a", "f1", 2.0}, {"src_m", "f3", 6.0},
	    {"s_a", "f3", 5.0},    {"src_z", "f2", 20.0}, {"src_m", "f1", 3.0},  {"s_b", "f2", 7.0}, {"src_m2", "f1", 2.0},
	};
}

size_t IndexOf(const std::vector<std::string> &ids, const std::string &id) {
	for (size_t i = 0; i < ids.size(); ++i) {
		if (ids[i] == id) {
			return i;
		}
	}
	FAIL("id not found in dictionary: " << id);
	return 0;
}

} // namespace

TEST_CASE("IngestDataset: dictionaries are lexicographic and independent of input order", "[sourcetracker]") {
	const Dataset ds = IngestDataset(ToyTable(), ToyMetadata());
	REQUIRE(ds.sample_ids == std::vector<std::string> {"s_a", "s_b", "src_m", "src_m2", "src_z"});
	REQUIRE(ds.feature_ids == std::vector<std::string> {"f1", "f2", "f3"});
}

TEST_CASE("IngestDataset: every cell lands on its own (feature, sample) and nothing else", "[sourcetracker]") {
	const auto table = ToyTable();
	const Dataset ds = IngestDataset(table, ToyMetadata());
	REQUIRE(ds.rows.size() == table.size());
	REQUIRE(ds.cols.size() == table.size());
	REQUIRE(ds.vals.size() == table.size());

	for (const auto &cell : table) {
		const auto fi = static_cast<int32_t>(IndexOf(ds.feature_ids, cell.feature_id));
		const auto si = static_cast<int32_t>(IndexOf(ds.sample_ids, cell.sample_id));
		size_t hits = 0;
		for (size_t k = 0; k < ds.rows.size(); ++k) {
			if (ds.rows[k] == fi && ds.cols[k] == si) {
				REQUIRE(ds.vals[k] == Approx(cell.count));
				++hits;
			}
		}
		REQUIRE(hits == 1);
	}

	// st3's tally has one cell per (source, feature), and the MAP column keyed by
	// feature id needs unique keys: both rest on (feature, sample) being unique here.
	std::set<std::pair<int32_t, int32_t>> keys;
	for (size_t k = 0; k < ds.rows.size(); ++k) {
		keys.emplace(ds.rows[k], ds.cols[k]);
	}
	REQUIRE(keys.size() == ds.rows.size());
}

TEST_CASE("IngestDataset: roles and environments are aligned to sample_ids", "[sourcetracker]") {
	const Dataset ds = IngestDataset(ToyTable(), ToyMetadata());
	REQUIRE(ds.is_source.size() == ds.sample_ids.size());
	REQUIRE(ds.envs.size() == ds.sample_ids.size());

	REQUIRE_FALSE(ds.is_source[IndexOf(ds.sample_ids, "s_a")]);
	REQUIRE_FALSE(ds.is_source[IndexOf(ds.sample_ids, "s_b")]);
	REQUIRE(ds.is_source[IndexOf(ds.sample_ids, "src_m")]);
	REQUIRE(ds.is_source[IndexOf(ds.sample_ids, "src_m2")]);
	REQUIRE(ds.is_source[IndexOf(ds.sample_ids, "src_z")]);

	REQUIRE(ds.envs[IndexOf(ds.sample_ids, "src_z")] == std::optional<std::string> {"sewage"});
	REQUIRE(ds.envs[IndexOf(ds.sample_ids, "src_m")] == std::optional<std::string> {"seawater"});
	REQUIRE(ds.envs[IndexOf(ds.sample_ids, "src_m2")] == std::optional<std::string> {"seawater"});
	// A sink's env is not a source environment; whatever the metadata says there is dropped.
	REQUIRE_FALSE(ds.envs[IndexOf(ds.sample_ids, "s_a")].has_value());
	REQUIRE_FALSE(ds.envs[IndexOf(ds.sample_ids, "s_b")].has_value());
}

TEST_CASE("IngestDataset: sample totals are per-sample sums of the cells", "[sourcetracker]") {
	const Dataset ds = IngestDataset(ToyTable(), ToyMetadata());
	REQUIRE(ds.sample_totals.size() == ds.sample_ids.size());
	REQUIRE(ds.sample_totals[IndexOf(ds.sample_ids, "s_a")] == Approx(7.0));
	REQUIRE(ds.sample_totals[IndexOf(ds.sample_ids, "s_b")] == Approx(8.0));
	REQUIRE(ds.sample_totals[IndexOf(ds.sample_ids, "src_m")] == Approx(9.0));
	REQUIRE(ds.sample_totals[IndexOf(ds.sample_ids, "src_m2")] == Approx(6.0));
	REQUIRE(ds.sample_totals[IndexOf(ds.sample_ids, "src_z")] == Approx(30.0));
}

TEST_CASE("IngestDataset: variable names and role values are matched case-insensitively", "[sourcetracker]") {
	// ReadWideMetadata keeps the relation's own column spelling, so the variable
	// may arrive as `Source_Sink`; users write roles in any case.
	std::vector<MetadataRow> m = {
	    {"s_a", "Source_Sink", "SINK"},
	    {"src_z", "Source_Sink", "Source"},
	    {"src_z", "ENV", "sewage"},
	};
	const Dataset ds = IngestDataset({{"s_a", "f1", 1.0}, {"src_z", "f1", 2.0}}, m);
	REQUIRE_FALSE(ds.is_source[IndexOf(ds.sample_ids, "s_a")]);
	REQUIRE(ds.is_source[IndexOf(ds.sample_ids, "src_z")]);
	REQUIRE(ds.envs[IndexOf(ds.sample_ids, "src_z")] == std::optional<std::string> {"sewage"});
}

TEST_CASE("IngestDataset: fails loud, naming the id, on every malformed input", "[sourcetracker]") {
	SECTION("a feature-table sample with no metadata row") {
		auto table = ToyTable();
		table.push_back({"s_ghost", "f1", 1.0});
		REQUIRE_THROWS_WITH(IngestDataset(table, ToyMetadata()), ContainsSubstring("s_ghost"));
		REQUIRE_THROWS_WITH(IngestDataset(table, ToyMetadata()), ContainsSubstring("metadata"));
	}
	SECTION("a metadata sample with no cells in the feature table") {
		// After the reader drops zero cells this is indistinguishable from an
		// all-zero sample, and both must fail: an empty sink has nothing to
		// attribute and an empty source is a silently missing environment member.
		auto meta = ToyMetadata();
		AddMeta(meta, "src_empty", "source", "sewage");
		REQUIRE_THROWS_WITH(IngestDataset(ToyTable(), meta), ContainsSubstring("src_empty"));
		auto meta2 = ToyMetadata();
		AddMeta(meta2, "s_empty", "sink", "");
		REQUIRE_THROWS_WITH(IngestDataset(ToyTable(), meta2), ContainsSubstring("s_empty"));
	}
	SECTION("a duplicated (sample, feature) cell") {
		auto table = ToyTable();
		table.push_back({"s_a", "f1", 3.0});
		REQUIRE_THROWS_WITH(IngestDataset(table, ToyMetadata()), ContainsSubstring("duplicate"));
		REQUIRE_THROWS_WITH(IngestDataset(table, ToyMetadata()), ContainsSubstring("s_a"));
		REQUIRE_THROWS_WITH(IngestDataset(table, ToyMetadata()), ContainsSubstring("f1"));
	}
	SECTION("a negative or non-finite count") {
		auto neg = ToyTable();
		neg.push_back({"s_b", "f1", -1.0});
		REQUIRE_THROWS_WITH(IngestDataset(neg, ToyMetadata()), ContainsSubstring("s_b"));
		REQUIRE_THROWS_WITH(IngestDataset(neg, ToyMetadata()), ContainsSubstring("f1"));
		auto nan = ToyTable();
		nan.push_back({"s_b", "f1", std::numeric_limits<double>::quiet_NaN()});
		REQUIRE_THROWS_WITH(IngestDataset(nan, ToyMetadata()), ContainsSubstring("s_b"));
		auto inf = ToyTable();
		inf.push_back({"s_b", "f1", std::numeric_limits<double>::infinity()});
		REQUIRE_THROWS_WITH(IngestDataset(inf, ToyMetadata()), ContainsSubstring("s_b"));
	}
	SECTION("a fractional count") {
		// st3 stores integer counts and floors what it is given, so a table of
		// relative abundances would be silently truncated toward zero. Refuse it,
		// as SourceTracker2 does.
		auto frac = ToyTable();
		frac.push_back({"s_b", "f1", 2.5});
		REQUIRE_THROWS_WITH(IngestDataset(frac, ToyMetadata()), ContainsSubstring("s_b"));
		REQUIRE_THROWS_WITH(IngestDataset(frac, ToyMetadata()), ContainsSubstring("whole"));
	}
	SECTION("an unknown source_sink value") {
		std::vector<MetadataRow> meta;
		AddMeta(meta, "src_z", "source", "sewage");
		AddMeta(meta, "s_a", "sauce", "");
		const std::vector<CooRow> table = {{"src_z", "f1", 1.0}, {"s_a", "f1", 1.0}};
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("s_a"));
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("sauce"));
	}
	SECTION("a source with an empty (NULL) env") {
		std::vector<MetadataRow> meta;
		AddMeta(meta, "src_z", "source", "");
		AddMeta(meta, "s_a", "sink", "");
		const std::vector<CooRow> table = {{"src_z", "f1", 1.0}, {"s_a", "f1", 1.0}};
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("src_z"));
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("env"));
	}
	SECTION("a sample listed twice in the metadata") {
		auto meta = ToyMetadata();
		AddMeta(meta, "s_a", "sink", "");
		REQUIRE_THROWS_WITH(IngestDataset(ToyTable(), meta), ContainsSubstring("s_a"));
		REQUIRE_THROWS_WITH(IngestDataset(ToyTable(), meta), ContainsSubstring("duplicate"));
	}
	SECTION("a sample with an env row but no source_sink row") {
		auto meta = ToyMetadata();
		meta.push_back({"s_c", "env", "sewage"});
		auto table = ToyTable();
		table.push_back({"s_c", "f1", 1.0});
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("s_c"));
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("source_sink"));
	}
	SECTION("no source samples at all") {
		std::vector<MetadataRow> meta;
		AddMeta(meta, "s_a", "sink", "");
		AddMeta(meta, "s_b", "sink", "");
		const std::vector<CooRow> table = {{"s_a", "f1", 1.0}, {"s_b", "f1", 1.0}};
		REQUIRE_THROWS_WITH(IngestDataset(table, meta), ContainsSubstring("no source"));
	}
}

TEST_CASE("CheckDepths: sink mode checks sinks against the sink depth, and only those", "[sourcetracker]") {
	// Totals: s_a 7, s_b 8 (sinks); src_m 9, src_m2 6, src_z 30 (sources).
	const Dataset ds = IngestDataset(ToyTable(), ToyMetadata());

	SECTION("depth 0 is off") {
		REQUIRE_NOTHROW(CheckDepths(ds, 0, 0, false));
	}
	SECTION("a sink at exactly the depth passes") {
		REQUIRE_NOTHROW(CheckDepths(ds, 0, 7, false));
	}
	SECTION("one sequence deeper than the shallowest sink fails, in SourceTracker2's words") {
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 8, false), ContainsSubstring("rarefaction of sink samples at 8"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 8, false), ContainsSubstring("1 sink sample"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 8, false), ContainsSubstring("shallowest"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 8, false), ContainsSubstring("'s_a'"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 8, false), ContainsSubstring("7 sequences"));
	}
	SECTION("the count of shallow sinks is reported") {
		REQUIRE_THROWS_WITH(CheckDepths(ds, 0, 1000, false), ContainsSubstring("2 sink samples"));
	}
	SECTION("source totals are st3's business in sink mode (environments are collapsed first)") {
		REQUIRE_NOTHROW(CheckDepths(ds, 1000000, 0, false));
	}
}

TEST_CASE("CheckDepths: leave-one-out checks source samples against the source depth, and only those",
          "[sourcetracker]") {
	const Dataset ds = IngestDataset(ToyTable(), ToyMetadata());

	SECTION("a source sample at exactly the depth passes") {
		REQUIRE_NOTHROW(CheckDepths(ds, 6, 0, true));
	}
	SECTION("the shallowest source sample is named") {
		REQUIRE_THROWS_WITH(CheckDepths(ds, 7, 0, true), ContainsSubstring("rarefaction of source samples at 7"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 7, 0, true), ContainsSubstring("'src_m2'"));
		REQUIRE_THROWS_WITH(CheckDepths(ds, 7, 0, true), ContainsSubstring("6 sequences"));
	}
	SECTION("the sink depth is ignored in leave-one-out") {
		REQUIRE_NOTHROW(CheckDepths(ds, 0, 1000000, true));
	}
}
