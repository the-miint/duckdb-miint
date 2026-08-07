// Pure-core tests for the absquant_* absolute-quantification functions.
//
// The numeric parity target is pysyndna (BSD-3-Clause), via the committed
// goldens in data/syndna/. Nothing here invokes pysyndna or scipy -- those ran
// offline to produce the CSVs; see data/syndna/README.md.
//
// Run from the repo root (the CSV paths are relative), which is what
// run_tests.sh and the CI invocation both do.

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

#include "absquant.hpp"

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

using miint::absquant::CellCountsMetric;
using miint::absquant::CellCountsOptions;
using miint::absquant::CellCountsResult;
using miint::absquant::ComputeCellCounts;
using miint::absquant::ComputeOrfCopies;
using miint::absquant::CountObservation;
using miint::absquant::CoverageObservation;
using miint::absquant::DenominatorColumnName;
using miint::absquant::DuplicatedIds;
using miint::absquant::FeatureLength;
using miint::absquant::FeatureTableValue;
using miint::absquant::FitOptions;
using miint::absquant::FitResult;
using miint::absquant::FitSyndnaModels;
using miint::absquant::FormatIdList;
using miint::absquant::IdsMissingFrom;
using miint::absquant::IsUsableSampleParameter;
using miint::absquant::Linregress;
using miint::absquant::LinregressResult;
using miint::absquant::MetricName;
using miint::absquant::OrfCoords;
using miint::absquant::OrfCopiesResult;
using miint::absquant::ParseCellCountsMetric;
using miint::absquant::RegularizedIncompleteBeta;
using miint::absquant::RejectDuplicateCells;
using miint::absquant::SampleCellParams;
using miint::absquant::SampleMass;
using miint::absquant::SampleOrfParams;
using miint::absquant::SampleRegression;
using miint::absquant::StudentTSurvival;
using miint::absquant::SyndnaConcentration;

namespace {

// Minimal CSV reader for the goldens: header row, no quoted fields, empty cell
// means NaN (that is what the generator's _fmt emits for a NaN). Local rather
// than shared because test_VsearchClusterWrapper.cpp already sets the precedent
// of a small fixture-specific parser per test file.
struct Csv {
	std::vector<std::string> header;
	std::vector<std::vector<std::string>> rows;

	size_t Col(const std::string &name) const {
		for (size_t i = 0; i < header.size(); ++i) {
			if (header[i] == name) {
				return i;
			}
		}
		FAIL("column '" << name << "' not found in golden");
		return 0;
	}
};

Csv ReadCsv(const std::string &path) {
	std::ifstream in(path);
	REQUIRE(in.is_open());
	Csv csv;
	std::string line;
	bool first = true;
	while (std::getline(in, line)) {
		if (!line.empty() && line.back() == '\r') {
			line.pop_back();
		}
		if (line.empty()) {
			continue;
		}
		std::vector<std::string> fields;
		std::stringstream ss(line);
		std::string field;
		while (std::getline(ss, field, ',')) {
			fields.push_back(field);
		}
		// A trailing empty field is lost by getline; restore it.
		if (!line.empty() && line.back() == ',') {
			fields.push_back("");
		}
		if (first) {
			csv.header = fields;
			first = false;
		} else {
			csv.rows.push_back(fields);
		}
	}
	return csv;
}

double AsDouble(const std::string &s) {
	if (s.empty()) {
		return std::nan("");
	}
	return std::strtod(s.c_str(), nullptr);
}

// Relative agreement with an absolute floor, so that values which legitimately
// underflow to 0.0 in the far tail compare sanely.
bool CloseEnough(double got, double want, double rel, double abs_floor) {
	if (std::isnan(got) && std::isnan(want)) {
		return true;
	}
	if (std::isnan(got) || std::isnan(want)) {
		return false;
	}
	const double diff = std::fabs(got - want);
	return diff <= abs_floor || diff <= rel * std::fabs(want);
}

// The parity bound data/syndna/README.md specifies and every parity test in the
// codebase uses: abs(mine - gold) <= tol * (1 + abs(gold)). The +1 keeps a
// golden of exactly 0.0 -- which the fit oracles do contain -- from demanding
// bit equality.
bool ParityOk(double got, double want, double tol) {
	return std::fabs(got - want) <= tol * (1.0 + std::fabs(want));
}

// pysyndna's public API truncates to 12 decimal places, so the goldens carry
// only ~1e-12 of absolute precision; 1e-11 is the bound README.md records for
// the fit oracles.
constexpr double kFitTol = 1e-11;

struct FitFixture {
	std::vector<CountObservation> counts;
	std::vector<SyndnaConcentration> concentrations;
	std::vector<SampleMass> masses;
};

// Loads the committed inputs for one fixture family ("fit" = Set A, pysyndna's
// own; "fitb" = Set B, the synthetic edge cases).
FitFixture LoadFitFixture(const std::string &prefix) {
	FitFixture fixture;
	const std::string dir = "data/syndna/";

	const auto counts = ReadCsv(dir + prefix + "_counts.csv");
	const size_t c_sample = counts.Col("sample_id");
	const size_t c_feature = counts.Col("feature_id");
	const size_t c_value = counts.Col("value");
	for (const auto &row : counts.rows) {
		fixture.counts.push_back({row[c_sample], row[c_feature], AsDouble(row[c_value])});
	}

	const auto concentrations = ReadCsv(dir + prefix + "_concentrations.csv");
	const size_t k_feature = concentrations.Col("feature_id");
	const size_t k_ng_ul = concentrations.Col("syndna_indiv_ng_ul");
	for (const auto &row : concentrations.rows) {
		fixture.concentrations.push_back({row[k_feature], AsDouble(row[k_ng_ul])});
	}

	// An empty mass cell is pysyndna's NaN, which is what SQL NULL will arrive
	// as; AsDouble turns it into one.
	const auto params = ReadCsv(dir + prefix + "_params.csv");
	const size_t p_sample = params.Col("sample_id");
	const size_t p_mass = params.Col("mass_syndna_input_ng");
	for (const auto &row : params.rows) {
		fixture.masses.push_back({row[p_sample], AsDouble(row[p_mass])});
	}
	return fixture;
}

// Compares a FitResult against a committed oracle in BOTH directions: every
// golden sample must be reproduced, and no sample we invent may survive. A
// one-directional check would pass a port that silently emitted extra models.
void CheckAgainstOracle(const FitResult &result, const std::string &oracle_path) {
	const auto gold = ReadCsv(oracle_path);
	const size_t g_sample = gold.Col("sample_id");
	const size_t g_slope = gold.Col("slope");
	const size_t g_intercept = gold.Col("intercept");
	const size_t g_rvalue = gold.Col("rvalue");
	const size_t g_pvalue = gold.Col("pvalue");
	const size_t g_stderr = gold.Col("stderr");
	const size_t g_intercept_stderr = gold.Col("intercept_stderr");

	std::map<std::string, LinregressResult> mine;
	for (const auto &model : result.models) {
		mine[model.sample_id] = model.fit;
	}
	// A duplicated sample id would otherwise hide behind the map.
	REQUIRE(mine.size() == result.models.size());
	CHECK(mine.size() == gold.rows.size());

	for (const auto &row : gold.rows) {
		const std::string sample_id = row[g_sample];
		INFO("sample " << sample_id << " of " << oracle_path);
		REQUIRE(mine.count(sample_id) == 1);
		const auto &fit = mine[sample_id];
		CHECK(ParityOk(fit.slope, AsDouble(row[g_slope]), kFitTol));
		CHECK(ParityOk(fit.intercept, AsDouble(row[g_intercept]), kFitTol));
		CHECK(ParityOk(fit.rvalue, AsDouble(row[g_rvalue]), kFitTol));
		CHECK(ParityOk(fit.pvalue, AsDouble(row[g_pvalue]), kFitTol));
		CHECK(ParityOk(fit.stderr_, AsDouble(row[g_stderr]), kFitTol));
		CHECK(ParityOk(fit.intercept_stderr, AsDouble(row[g_intercept_stderr]), kFitTol));
	}

	for (const auto &model : result.models) {
		bool in_gold = false;
		for (const auto &row : gold.rows) {
			if (row[g_sample] == model.sample_id) {
				in_gold = true;
				break;
			}
		}
		INFO("sample " << model.sample_id << " is not in " << oracle_path);
		CHECK(in_gold);
	}
}

bool Contains(const std::vector<std::string> &haystack, const std::string &needle) {
	return std::find(haystack.begin(), haystack.end(), needle) != haystack.end();
}

// The smallest input FitSyndnaModels accepts, for the validation tests to break
// one field at a time. Counts are 10^3/10^2/10^1 against concentrations
// 1/0.1/0.01 summing to the sample's mass, so the untouched fixture fits a
// perfect line -- any test below that still produces a model is unambiguous.
struct TinyFit {
	std::vector<CountObservation> counts = {{"s1", "f1", 1000.0}, {"s1", "f2", 100.0}, {"s1", "f3", 10.0}};
	std::vector<SyndnaConcentration> concentrations = {{"f1", 1.0}, {"f2", 0.1}, {"f3", 0.01}};
	std::vector<SampleMass> masses = {{"s1", 1.11}};
	FitOptions options;

	TinyFit() {
		options.min_syndna_counts = 1;
		options.syndna_contributing_fraction = 1.0;
	}
	FitResult Run() const {
		return FitSyndnaModels(counts, concentrations, masses, options);
	}
};

// data/syndna/README.md's bound for the cells AND the ORF oracles -- one row of
// its table, so one constant here. Looser than the fit oracles' 1e-11 because
// these are untruncated doubles reaching ~1e13, so the bound has to be
// relative.
//
// Both ORF oracles in fact reproduce EXACTLY, all 26 golden cells at a
// difference of zero, because that chain is multiply/divide/fabs with no pow or
// log10 anywhere. The bound stays where it is rather than dropping to 0 for
// them, because exactness is a property of one build and not of the arithmetic:
// GCC contracts a*b/c into an FMA by default, so the last bit is the compiler's
// to choose. It is still worth knowing which side of that line a future ORF
// failure lands on -- anything above about 1e-16 there is a real regression,
// not tolerance.
constexpr double kFeatureTableTol = 1e-9;

// One cell of a DENSE golden, already selected out of whichever oracle file it
// came from.
struct GoldenCell {
	std::string sample_id;
	std::string feature_id;
	double value = 0.0;
};

// Compares a long-form output feature table against a dense golden in BOTH
// directions.
//
// The golden is dense where pysyndna's pivot+fillna made it so, and miint omits
// zero-valued cells, so a golden 0.0 is matched against an ABSENT cell -- that
// is the D10 claim, and comparing this way proves it rather than assuming it.
// The reverse direction still catches any cell miint invents.
//
// Shared by the cells and the ORF oracles because both functions return the
// same sparse shape under the same invariant; only the selection of golden rows
// differs, and that is the caller's business.
void CheckValuesAgainstGolden(const std::vector<FeatureTableValue> &values, const std::vector<GoldenCell> &gold,
                              const std::string &label) {
	std::map<std::pair<std::string, std::string>, double> mine;
	for (const auto &cell : values) {
		mine[{cell.sample_id, cell.feature_id}] = cell.value;
	}
	// A duplicated cell would otherwise hide behind the map.
	REQUIRE(mine.size() == values.size());
	// An empty golden would make every check below vacuous, so a mis-selected
	// oracle fails here rather than passing silently.
	REQUIRE(!gold.empty());

	std::set<std::pair<std::string, std::string>> golden_keys;
	for (const auto &row : gold) {
		const std::pair<std::string, std::string> key {row.sample_id, row.feature_id};
		golden_keys.insert(key);
		INFO("cell (" << key.first << ", " << key.second << ") of " << label);
		const auto found = mine.find(key);
		const double got = (found == mine.end()) ? 0.0 : found->second;
		CHECK(ParityOk(got, row.value, kFeatureTableTol));
	}

	for (const auto &cell : values) {
		INFO("cell (" << cell.sample_id << ", " << cell.feature_id << ") is not in " << label);
		CHECK(golden_keys.count({cell.sample_id, cell.feature_id}) == 1);
	}
}

struct CellsFixture {
	std::vector<CountObservation> counts;
	std::vector<SampleRegression> models;
	std::vector<CoverageObservation> coverage;
	std::vector<FeatureLength> lengths;
	std::vector<SampleCellParams> params;
	CellCountsOptions options;

	CellCountsResult Run() const {
		return ComputeCellCounts(counts, models, coverage, lengths, params, options);
	}
};

// Loads one cell-count fixture family ("cells" = Set A, pysyndna's own;
// "cellsb" = Set B, the synthetic filter cases) for one metric.
//
// Reads the metric's denominator column and ONLY that one, which is what the
// DuckDB wrapper does too -- asking for cells_per_g_of_gdna must not require a
// sample_volume_ul column to exist.
CellsFixture LoadCellsFixture(const std::string &prefix, CellCountsMetric metric = CellCountsMetric::CellsPerGOfGdna) {
	CellsFixture fixture;
	fixture.options.metric = metric;
	const std::string dir = "data/syndna/";

	const auto counts = ReadCsv(dir + prefix + "_counts.csv");
	const size_t c_sample = counts.Col("sample_id");
	const size_t c_feature = counts.Col("feature_id");
	const size_t c_value = counts.Col("value");
	for (const auto &row : counts.rows) {
		// The DuckDB reader drops zero-valued cells before the core ever sees
		// them (the sparse invariant), so the fixture loader must too -- these
		// CSVs are stored as dense as pysyndna's own input.
		const double value = AsDouble(row[c_value]);
		if (value == 0.0) {
			continue;
		}
		fixture.counts.push_back({row[c_sample], row[c_feature], value});
	}

	const auto models = ReadCsv(dir + prefix + "_models.csv");
	const size_t m_sample = models.Col("sample_id");
	const size_t m_slope = models.Col("slope");
	const size_t m_intercept = models.Col("intercept");
	const size_t m_rvalue = models.Col("rvalue");
	for (const auto &row : models.rows) {
		fixture.models.push_back(
		    {row[m_sample], AsDouble(row[m_slope]), AsDouble(row[m_intercept]), AsDouble(row[m_rvalue])});
	}

	const auto coverage = ReadCsv(dir + prefix + "_coverage.csv");
	const size_t v_sample = coverage.Col("sample_id");
	const size_t v_feature = coverage.Col("feature_id");
	const size_t v_coverage = coverage.Col("coverage");
	for (const auto &row : coverage.rows) {
		fixture.coverage.push_back({row[v_sample], row[v_feature], AsDouble(row[v_coverage])});
	}

	const auto lengths = ReadCsv(dir + prefix + "_lengths.csv");
	const size_t l_feature = lengths.Col("feature_id");
	const size_t l_len = lengths.Col("ogu_len_in_bp");
	for (const auto &row : lengths.rows) {
		fixture.lengths.push_back({row[l_feature], AsDouble(row[l_len])});
	}

	const auto params = ReadCsv(dir + prefix + "_params.csv");
	const size_t p_sample = params.Col("sample_id");
	const size_t p_mass = params.Col("sequenced_sample_gdna_mass_ng");
	const size_t p_conc = params.Col("extracted_gdna_concentration_ng_ul");
	const size_t p_vol = params.Col("vol_extracted_elution_ul");
	const char *denominator_column = DenominatorColumnName(metric);
	const size_t p_denominator = denominator_column == nullptr ? 0 : params.Col(denominator_column);
	for (const auto &row : params.rows) {
		fixture.params.push_back({row[p_sample], AsDouble(row[p_mass]), AsDouble(row[p_conc]), AsDouble(row[p_vol]),
		                          denominator_column == nullptr ? 0.0 : AsDouble(row[p_denominator])});
	}
	return fixture;
}

// Compares cell-count output against one metric's rows of a committed oracle.
//
// The cells oracles carry all four metrics in one file, so the selection is
// what distinguishes this from the ORF case; everything else is the shared
// two-directional comparison above.
void CheckCellsAgainstOracle(const CellCountsResult &result, const std::string &oracle_path,
                             const std::string &metric) {
	const auto gold = ReadCsv(oracle_path);
	const size_t g_metric = gold.Col("metric");
	const size_t g_sample = gold.Col("sample_id");
	const size_t g_feature = gold.Col("feature_id");
	const size_t g_value = gold.Col("value");

	std::vector<GoldenCell> cells;
	for (const auto &row : gold.rows) {
		if (row[g_metric] != metric) {
			continue;
		}
		cells.push_back({row[g_sample], row[g_feature], AsDouble(row[g_value])});
	}
	CheckValuesAgainstGolden(result.values, cells, oracle_path + " (" + metric + ")");
}

} // namespace

TEST_CASE("RegularizedIncompleteBeta known values", "[absquant]") {
	// Boundaries are exact by definition, not by approximation.
	CHECK(RegularizedIncompleteBeta(2.0, 3.0, 0.0) == 0.0);
	CHECK(RegularizedIncompleteBeta(2.0, 3.0, 1.0) == 1.0);

	// I_x(1,1) is the identity -- the simplest closed form, and one a wrong
	// continued fraction fails immediately.
	for (double x = 0.0; x <= 1.0; x += 0.125) {
		CHECK(CloseEnough(RegularizedIncompleteBeta(1.0, 1.0, x), x, 1e-14, 1e-15));
	}

	// I_x(a,b) = 1 - I_{1-x}(b,a): the symmetry the reflection branch relies on,
	// so this fails loudly if the branch is wrong on either side.
	const double pairs[][3] = {{0.5, 2.5, 0.3}, {2.5, 0.5, 0.3}, {3.0, 7.0, 0.75}, {0.1, 0.2, 0.05}, {30.0, 0.5, 0.9}};
	for (auto &p : pairs) {
		const double lhs = RegularizedIncompleteBeta(p[0], p[1], p[2]);
		const double rhs = 1.0 - RegularizedIncompleteBeta(p[1], p[0], 1.0 - p[2]);
		CHECK(CloseEnough(lhs, rhs, 1e-12, 1e-14));
	}

	// I_x(a,b) is a CDF: monotone non-decreasing in x, and within [0,1].
	double prev = -1.0;
	for (double x = 0.0; x <= 1.0; x += 0.05) {
		const double v = RegularizedIncompleteBeta(2.5, 4.5, x);
		CHECK(v >= 0.0);
		CHECK(v <= 1.0);
		CHECK(v >= prev);
		prev = v;
	}

	// Out of domain -> NaN, rather than a plausible-looking number.
	CHECK(std::isnan(RegularizedIncompleteBeta(0.0, 1.0, 0.5)));
	CHECK(std::isnan(RegularizedIncompleteBeta(1.0, -1.0, 0.5)));
	CHECK(std::isnan(RegularizedIncompleteBeta(1.0, 1.0, 1.5)));
	CHECK(std::isnan(RegularizedIncompleteBeta(1.0, 1.0, -0.5)));
}

TEST_CASE("StudentTSurvival identities", "[absquant]") {
	// sf(0, df) == 0.5 EXACTLY for every df. This is the one value that must be
	// exact rather than approximate, and the first thing a broken implementation
	// gets wrong.
	for (double df : {1.0, 2.0, 3.0, 7.0, 30.0, 1000.0}) {
		CHECK(StudentTSurvival(0.0, df) == 0.5);
	}

	// sf(-t) == 1 - sf(t). The fit only ever passes |t|, but an asymmetric
	// implementation would be a latent trap for every later caller.
	for (double df : {1.0, 4.0, 25.0}) {
		for (double t : {0.25, 1.0, 2.5, 10.0}) {
			CHECK(CloseEnough(StudentTSurvival(-t, df), 1.0 - StudentTSurvival(t, df), 1e-12, 1e-14));
		}
	}

	// df = 1 is the Cauchy distribution, where sf has the closed form
	// 0.5 - atan(t)/pi. An independent formula, so this does not merely restate
	// the implementation. (Spelled out rather than using M_PI, which is not
	// standard C++ and is absent under some of this project's targets.)
	constexpr double kPi = 3.14159265358979323846;
	for (double t : {0.0, 0.5, 1.0, 3.0, 25.0}) {
		const double cauchy = 0.5 - std::atan(t) / kPi;
		CHECK(CloseEnough(StudentTSurvival(t, 1.0), cauchy, 1e-12, 1e-15));
	}

	// df = 2 also has a closed form: 0.5 * (1 - t/sqrt(2+t^2)).
	for (double t : {0.0, 0.5, 1.0, 3.0, 25.0}) {
		const double closed = 0.5 * (1.0 - t / std::sqrt(2.0 + t * t));
		CHECK(CloseEnough(StudentTSurvival(t, 2.0), closed, 1e-12, 1e-15));
	}

	// Monotone decreasing in t.
	double prev = 1.0;
	for (double t = -5.0; t <= 5.0; t += 0.25) {
		const double v = StudentTSurvival(t, 6.0);
		CHECK(v <= prev + 1e-15);
		prev = v;
	}

	// Invalid df -> NaN rather than a number that would silently become a pvalue.
	CHECK(std::isnan(StudentTSurvival(1.0, 0.0)));
	CHECK(std::isnan(StudentTSurvival(1.0, -3.0)));
	CHECK(std::isnan(StudentTSurvival(std::nan(""), 5.0)));
}

TEST_CASE("StudentTSurvival matches the scipy grid", "[absquant]") {
	// The real test bed. data/syndna/studentt_sf_oracle.csv is 521 points from
	// scipy.stats.t.sf spanning df 1..1000 and |t| from 0 to 2.2e11, including
	// the extreme t values scipy's TINY = 1e-20 guard produces at a perfect fit.
	// The fit goldens alone contain only six pvalues, all at df <= 6, so without
	// this grid the incomplete beta would be almost untested.
	const auto csv = ReadCsv("data/syndna/studentt_sf_oracle.csv");
	REQUIRE(csv.rows.size() == 521);

	const size_t ct = csv.Col("t");
	const size_t cdf = csv.Col("df");
	const size_t csf = csv.Col("sf");

	// Two of the 521 points are EXCLUDED because scipy itself is the inaccurate
	// side there, and pinning them would enshrine scipy's error as our target.
	// Both are df = 1 (Cauchy) at tiny t, where scipy's stdtr loses precision;
	// they are asserted against the exact closed form in the test below instead.
	// Established by arbitrating every grid point against a 40-digit mpmath
	// evaluation: scipy's relative error exceeds 1e-13 at exactly these two and
	// nowhere else.
	//
	//   t=1e-12 df=1  scipy 0.5                 rel err 6.4e-13
	//   t=1e-08 df=1  scipy 0.4999999952568131  rel err 3.1e-09
	//
	// This matters for parity, not just pedantry: pysyndna's pvalue IS scipy's,
	// so scipy is the right oracle everywhere it is accurate -- and a real fit
	// never produces a t that small, so the excluded region is unreachable from
	// absquant_fit_models regardless.
	auto scipy_is_unreliable = [](double t, double df) {
		return df == 1.0 && std::fabs(t) > 0.0 && std::fabs(t) <= 1e-8;
	};

	// Tolerance is MEASURED, not guessed: the worst relative error across the
	// grid is ~2.6e-12, at df = 1000 where the front factor's lgamma difference
	// (lgamma(500.5) - lgamma(500) - lgamma(0.5), each term ~2600) cancels away
	// about 13 of 16 digits. That is six orders of magnitude tighter than the
	// ~1e-6 a pvalue is asserted to, and df = 1000 would require 1002 synDNA
	// features -- the real fixtures have 8 to 10, i.e. df <= 8, where the
	// measured error is below 1e-15.
	constexpr double kGridTol = 1e-11;

	size_t checked = 0;
	size_t tail_underflow = 0;
	size_t excluded = 0;
	double worst_rel = 0.0;
	std::string worst_at;
	for (const auto &row : csv.rows) {
		const double t = AsDouble(row[ct]);
		const double df = AsDouble(row[cdf]);
		const double want = AsDouble(row[csf]);
		const double got = StudentTSurvival(t, df);

		if (want == 0.0) {
			// scipy underflowed the far tail to exactly zero; we must too, or at
			// least be indistinguishable from it.
			CHECK(got >= 0.0);
			CHECK(got < 1e-300);
			++tail_underflow;
			continue;
		}
		if (scipy_is_unreliable(t, df)) {
			++excluded;
			continue;
		}
		const double rel = std::fabs(got - want) / std::fabs(want);
		if (rel > worst_rel) {
			worst_rel = rel;
			worst_at = "t=" + row[ct] + " df=" + row[cdf];
		}
		INFO("t=" << row[ct] << " df=" << row[cdf] << " want=" << row[csf] << " got=" << std::setprecision(17) << got
		          << " rel=" << std::setprecision(3) << rel);
		CHECK(CloseEnough(got, want, kGridTol, 0.0));
		++checked;
	}
	// Guard against the grid silently becoming all-underflow or all-excluded,
	// which would make the loop above assert almost nothing.
	CHECK(checked > 450);
	CHECK(tail_underflow > 0);
	CHECK(excluded == 2);

	// Report the worst case so a regression in accuracy is visible as a number,
	// not just a pass/fail. If this ever climbs near kGridTol, the front factor
	// needs a more careful log-beta rather than a looser bound.
	INFO("worst relative error " << std::setprecision(3) << worst_rel << " at " << worst_at);
	CHECK(worst_rel < kGridTol);
}

TEST_CASE("StudentTSurvival is exact where scipy is not", "[absquant]") {
	// The two grid points excluded above. df = 1 is Cauchy, so the survival
	// function has the closed form 0.5 - atan(t)/pi and there is a right answer
	// to check against rather than an oracle to defer to.
	//
	// scipy's stdtr loses precision here: at t = 1e-8 it returns
	// 0.4999999952568131, which is wrong by 3.1e-9 relative -- confirmed against
	// a 40-digit mpmath evaluation, which agrees with the closed form and with
	// this implementation to ~4e-17.
	//
	// This test exists so nobody later "corrects" the implementation to match
	// scipy at these points and makes it strictly worse. Being more accurate
	// than the oracle is not a defect, but it is surprising, so it is pinned.
	constexpr double kPi = 3.14159265358979323846;
	for (double t : {1e-12, 1e-8}) {
		const double exact = 0.5 - std::atan(t) / kPi;
		INFO("t=" << t << " exact=" << std::setprecision(17) << exact << " got=" << StudentTSurvival(t, 1.0));
		CHECK(CloseEnough(StudentTSurvival(t, 1.0), exact, 1e-15, 0.0));
	}

	// And the value scipy would have given at t = 1e-8, so the divergence is
	// recorded as a concrete number rather than a claim in a comment.
	CHECK_FALSE(CloseEnough(StudentTSurvival(1e-8, 1.0), 0.4999999952568131, 1e-12, 0.0));
}

TEST_CASE("StudentTSurvival at the linregress TINY extremes", "[absquant]") {
	// scipy's linregress adds TINY = 1e-20 inside the t statistic so a perfect
	// fit (r == +/-1) produces a huge but FINITE t instead of dividing by zero.
	// This reproduces that computation and checks the survival function stays
	// finite and positive there -- the property that makes a perfect fit report
	// pvalue ~9.0e-11 at df = 1 rather than NaN.
	const double kTiny = 1e-20;
	for (double df : {1.0, 2.0, 5.0, 30.0}) {
		const double r = 1.0;
		const double t = r * std::sqrt(df / ((1.0 - r + kTiny) * (1.0 + r + kTiny)));
		REQUIRE(std::isfinite(t));
		REQUIRE(t > 1e9);
		const double sf = StudentTSurvival(t, df);
		CHECK(std::isfinite(sf));
		CHECK(sf >= 0.0);
	}

	// df = 1 is the case with a real (non-underflowing) answer; pvalue = 2*sf.
	const double r = 1.0;
	const double t1 = r * std::sqrt(1.0 / ((1.0 - r + kTiny) * (1.0 + r + kTiny)));
	CHECK(CloseEnough(StudentTSurvival(t1, 1.0), 4.501582e-11, 1e-6, 0.0));
}

TEST_CASE("Linregress reproduces scipy's arithmetic on a hand-checked case", "[absquant]") {
	// x = 1,2,3,4  y = 1,3,2,5, chosen so every quantity is a short rational and
	// the expected values can be derived on paper rather than copied out of an
	// implementation:
	//
	//   xmean 2.5   ymean 2.75
	//   sum dx^2 = 5      -> ssxm  = 5/4    = 1.25      (BIASED: /n, not /(n-1))
	//   sum dx dy = 5.5   -> ssxym = 5.5/4  = 1.375
	//   sum dy^2 = 8.75   -> ssym  = 8.75/4 = 2.1875
	//   slope = 1.375/1.25 = 1.1 ; intercept = 2.75 - 1.1*2.5 = 0
	//   stderr^2 = (ssym - ssxym^2/ssxm)/ssxm/df = 0.675/2.5 = 0.27
	//   intercept_stderr = stderr*sqrt(ssxm + xmean^2) = sqrt(0.27*7.5)
	//
	// (Cross-checked against scipy.stats.linregress during development, which
	// returns exactly these values.)
	const std::vector<double> x = {1.0, 2.0, 3.0, 4.0};
	const std::vector<double> y = {1.0, 3.0, 2.0, 5.0};
	const auto fit = Linregress(x, y);
	REQUIRE(fit.ok);
	CHECK(CloseEnough(fit.slope, 1.1, 1e-14, 1e-15));
	CHECK(CloseEnough(fit.intercept, 0.0, 1e-14, 1e-14));
	CHECK(CloseEnough(fit.rvalue, 1.375 / std::sqrt(1.25 * 2.1875), 1e-14, 1e-15));
	CHECK(CloseEnough(fit.stderr_, std::sqrt(0.27), 1e-14, 1e-15));
	CHECK(CloseEnough(fit.intercept_stderr, std::sqrt(2.025), 1e-14, 1e-15));

	// The ratio intercept_stderr/stderr is the sharpest single discriminator in
	// the whole port: it isolates sqrt(ssxm + xmean^2), where ssxm appears
	// un-normalized, so it separates all three plausible formulas at once.
	// Everywhere else the biased/unbiased choice cancels inside a ratio.
	const double ratio = fit.intercept_stderr / fit.stderr_;
	CHECK(CloseEnough(ratio, std::sqrt(1.25 + 6.25), 1e-13, 1e-15));          // scipy: biased ssxm
	CHECK_FALSE(CloseEnough(ratio, std::sqrt(5.0 / 3.0 + 6.25), 1e-6, 0.0));  // unbiased ssxm
	CHECK_FALSE(CloseEnough(ratio, std::sqrt(0.25 + 6.25 / 5.0), 1e-6, 0.0)); // textbook form

	// A real p-value, not a placeholder.
	CHECK(fit.pvalue > 0.0);
	CHECK(fit.pvalue < 1.0);
}

TEST_CASE("Linregress handles the degenerate inputs scipy special-cases", "[absquant]") {
	SECTION("perfect fit stays finite -- this is what TINY buys") {
		// Without scipy's TINY the t statistic divides by zero here and every
		// downstream field becomes NaN, which would discard the sample.
		const std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
		const std::vector<double> y = {3.0, 5.0, 7.0, 9.0, 11.0};
		const auto fit = Linregress(x, y);
		REQUIRE(fit.ok);
		CHECK(CloseEnough(fit.slope, 2.0, 1e-14, 1e-15));
		CHECK(CloseEnough(fit.intercept, 1.0, 1e-14, 1e-14));
		CHECK(CloseEnough(fit.rvalue, 1.0, 1e-15, 1e-15));
		CHECK(std::isfinite(fit.pvalue));
		CHECK(fit.pvalue >= 0.0);
		// The residual form of stderr, rather than sqrt(1 - r^2)*..., is why this
		// is ~1e-17 instead of ~1e-8 and why it does not depend on r landing on
		// exactly 1.0. See the note on Linregress in absquant.hpp.
		CHECK(fit.stderr_ >= 0.0);
		CHECK(fit.stderr_ < 1e-11);
		CHECK(fit.intercept_stderr < 1e-11);
	}

	SECTION("exactly two points: scipy zeroes both standard errors") {
		const auto fit = Linregress({1.0, 2.0}, {5.0, 3.0});
		REQUIRE(fit.ok);
		CHECK(fit.slope == -2.0);
		CHECK(fit.intercept == 7.0);
		CHECK(fit.rvalue == -1.0);
		// Hard-set, not merely small: with n == 2 there is no residual to
		// estimate from, and scipy declines to pretend otherwise.
		CHECK(fit.pvalue == 0.0);
		CHECK(fit.stderr_ == 0.0);
		CHECK(fit.intercept_stderr == 0.0);
	}

	SECTION("two identical y values: pvalue 1.0, but the sample is still lost") {
		// scipy's other n == 2 arm. It can never reach a caller, because equal
		// y means ssym == 0 which makes rvalue NaN and voids the sample -- so
		// this pins the arm AND the fact that it is unobservable.
		const auto fit = Linregress({1.0, 2.0}, {4.0, 4.0});
		CHECK(fit.pvalue == 1.0);
		CHECK(std::isnan(fit.rvalue));
		CHECK_FALSE(fit.ok);
	}

	SECTION("constant x with n > 1: scipy raises, we report no model") {
		// x = 2 makes the arithmetic decide this on its own -- 2+2+2 over 3 is
		// exactly 2, so every deviation is exactly zero and ssxm is exactly zero.
		CHECK_FALSE(Linregress({2.0, 2.0, 2.0}, {1.0, 5.0, 9.0}).ok);

		// x = 0.1 does NOT. 0.1+0.1+0.1 is 0.30000000000000004, whose third is
		// 0.10000000000000002, so each deviation is -1.4e-17 and ssxm lands on
		// 1.9e-34 instead of 0. Drop the explicit rejection and the fit runs to
		// completion on that noise: slope 0.0, rvalue 0.0, stderr 2.4e17 -- a
		// finished model, no NaN anywhere, and complete nonsense.
		CHECK_FALSE(Linregress({0.1, 0.1, 0.1}, {1.0, 5.0, 9.0}).ok);
		// Same input shape, and here the noise invents a slope of 1.33.
		CHECK_FALSE(Linregress({0.1, 0.1, 0.1}, {0.1, 0.2, 0.4}).ok);
	}

	SECTION("constant y: rvalue is NaN even though the line is well defined") {
		// The trap in pysyndna's all-fields-must-be-non-NaN rule: slope and
		// intercept are perfectly good, and the sample is discarded anyway.
		const auto fit = Linregress({1.0, 2.0, 3.0}, {7.0, 7.0, 7.0});
		CHECK(fit.slope == 0.0);
		CHECK(fit.intercept == 7.0);
		CHECK(std::isnan(fit.rvalue));
		CHECK_FALSE(fit.ok);
		// scipy's stderr carries the NaN r through; the residual form does not
		// on its own, because constant y has exactly zero residuals and would
		// report a confident 0.0. LinregressResult is documented as holding what
		// scipy returns, so the propagation is restored explicitly and pinned
		// here -- ok is already false either way, so nothing else would notice.
		CHECK(std::isnan(fit.stderr_));
		CHECK(std::isnan(fit.intercept_stderr));
	}

	SECTION("TINY turns a perfect fit's p-value into a number, not a floor") {
		// x = 0,1,2,3 against y = 2x is the smallest case where rvalue comes out
		// as EXACTLY 1.0 rather than one ulp short: ssxm = 1.25, ssym = 5.0 and
		// ssxym = 2.5 are all exact, and sqrt(6.25) is exactly 2.5.
		//
		// That is the only regime where scipy's TINY = 1e-20 does anything. With
		// it, t = sqrt(df/(2*TINY)) = 1e10 and the p-value is 1e-20; without it
		// the denominator is exactly zero, t is infinite, and the p-value
		// collapses to exactly 0.0. Both are "essentially zero" to a reader,
		// which is why this needs an explicit test -- every golden's perfect-fit
		// pvalue truncates to 0.0 and cannot tell the two apart.
		const auto fit = Linregress({0.0, 1.0, 2.0, 3.0}, {0.0, 2.0, 4.0, 6.0});
		REQUIRE(fit.ok);
		CHECK(fit.rvalue == 1.0);
		CHECK(fit.slope == 2.0);
		CHECK(fit.pvalue > 0.0);
		CHECK(CloseEnough(fit.pvalue, 1e-20, 1e-6, 0.0));
	}

	SECTION("one point, no points, and mismatched lengths yield no model") {
		const auto one = Linregress({3.0}, {4.0});
		CHECK_FALSE(one.ok);
		CHECK(std::isnan(one.slope));
		CHECK(std::isnan(one.rvalue));
		CHECK_FALSE(Linregress({}, {}).ok);
		CHECK_FALSE(Linregress({1.0, 2.0, 3.0}, {1.0, 2.0}).ok);
	}

	SECTION("a NaN anywhere voids the model") {
		CHECK_FALSE(Linregress({1.0, 2.0, 3.0}, {1.0, std::nan(""), 3.0}).ok);
	}

	SECTION("rvalue is clamped into [-1, 1]") {
		// Rounding can push |r| a hair past 1 on near-perfect data; everything
		// downstream takes a square root of 1 - r^2, so the clamp is load-bearing.
		const std::vector<double> x = {1e6, 2e6, 3e6, 4e6, 5e6, 6e6};
		std::vector<double> y;
		for (double v : x) {
			y.push_back(-7.0 * v + 1e9);
		}
		const auto fit = Linregress(x, y);
		REQUIRE(fit.ok);
		CHECK(fit.rvalue >= -1.0);
		CHECK(fit.rvalue <= 1.0);
	}
}

TEST_CASE("Linregress's stderr stays accurate as the fit approaches perfect", "[absquant]") {
	// The property that justifies computing stderr from the residuals rather
	// than from scipy's sqrt((1 - r^2)*ssym/ssxm/df), stated without needing an
	// oracle: scaling every residual by k must scale stderr by exactly k, since
	// the fit is linear and the residual pattern is unchanged.
	//
	// scipy's form cannot satisfy this. 1 - r^2 is cancellation, so it sheds a
	// digit per decade -- at wobble 1e-4 (rvalue 0.999999999923) scipy returns
	// 1.0103704639979164e-05 where the true value is 1.0103717379024465e-05,
	// already wrong in the sixth digit, and by wobble 1e-7 it returns exactly
	// 0.0 for data that visibly has residuals. The 1e-9 bound below is chosen
	// to sit between the two: the residual form holds to ~4e-13 across this
	// range, scipy's drifts to 1.3e-6.
	//
	// This is not pedantry about a fixture. A clean synDNA standard curve is
	// exactly the regime that gets here, so an M6 parity run against pysyndna
	// on good data will show a disagreement that is this port being right.
	const std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
	const double wobble[] = {1.0, -1.0, 0.5, -0.5, 0.25, -0.25, 0.125, -0.125};

	double reference_ratio = 0.0;
	for (double scale : {1e-1, 1e-2, 1e-3, 1e-4}) {
		std::vector<double> y;
		for (size_t i = 0; i < x.size(); ++i) {
			y.push_back(2.0 * x[i] + 1.0 + scale * wobble[i]);
		}
		const auto fit = Linregress(x, y);
		REQUIRE(fit.ok);
		// The fit really is near-perfect and really is imperfect -- otherwise
		// the scaling claim below would be vacuous. rvalue runs from 0.9999231
		// at the loosest wobble to 0.999999999923 at the tightest, i.e. across
		// and past the point where scipy's form starts shedding digits.
		CHECK(fit.rvalue < 1.0);
		CHECK(fit.rvalue > 0.9999);
		CHECK(fit.stderr_ > 0.0);

		const double ratio = fit.stderr_ / scale;
		if (reference_ratio == 0.0) {
			reference_ratio = ratio;
			// Anchors the scale so a uniformly-wrong-by-a-constant-factor
			// implementation cannot pass on the ratios alone.
			CHECK(CloseEnough(ratio, 0.1010371737902853, 1e-12, 0.0));
		} else {
			INFO("wobble " << scale << " gave stderr/scale " << std::setprecision(17) << ratio << " against "
			               << reference_ratio);
			CHECK(CloseEnough(ratio, reference_ratio, 1e-9, 0.0));
		}
	}
}

TEST_CASE("FitSyndnaModels matches pysyndna's fit goldens", "[absquant]") {
	SECTION("Set A -- pysyndna's own fixtures, min_syndna_counts = 50") {
		const auto fixture = LoadFitFixture("fit");
		FitOptions options;
		options.min_syndna_counts = 50;
		options.syndna_contributing_fraction = 1.0;
		const auto result = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, options);
		CheckAgainstOracle(result, "data/syndna/fit_models_oracle.csv");
	}

	SECTION("Set B -- synthetic edge cases, min_syndna_counts = 50") {
		const auto fixture = LoadFitFixture("fitb");
		FitOptions options;
		options.min_syndna_counts = 50;
		options.syndna_contributing_fraction = 1.0;
		const auto result = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, options);
		CheckAgainstOracle(result, "data/syndna/fitb_models_oracle.csv");
	}

	SECTION("Set B at the strict-< boundary, min_syndna_counts = 20") {
		// b6 totals EXACTLY 20 reads. `<` keeps it, `<=` drops it -- and keeping
		// it is what lets sCONSTY, sFLAT and sONE fit at all, so the golden has
		// 8 models here against 5 above. An off-by-one in the comparison changes
		// the sample set, not just a digit.
		const auto fixture = LoadFitFixture("fitb");
		FitOptions options;
		options.min_syndna_counts = 20;
		options.syndna_contributing_fraction = 1.0;
		const auto result = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, options);
		CHECK(result.dropped_syndna_ids.empty());
		CheckAgainstOracle(result, "data/syndna/fitb_boundary_models_oracle.csv");
	}
}

TEST_CASE("FitSyndnaModels accounts for every sample it was given", "[absquant]") {
	// Set B was built so that six of its eleven samples deliberately produce no
	// model, by six different routes. Silence is not an acceptable answer for
	// any of them: a caller must be able to say why a sample vanished, so the
	// three outcome lists have to partition the input.
	const auto fixture = LoadFitFixture("fitb");
	FitOptions options;
	options.min_syndna_counts = 50;
	const auto result = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, options);

	CHECK(result.models.size() == 5);
	for (const char *id : {"sA", "sB", "sPERFECT", "sZERO", "sTWO"}) {
		INFO("expected a model for " << id);
		bool found = false;
		for (const auto &model : result.models) {
			found = found || model.sample_id == id;
		}
		CHECK(found);
	}

	// NULL and negative mass_syndna_input_ng: rejected on the parameters, before
	// any arithmetic.
	CHECK(result.filtered_sample_ids.size() == 2);
	CHECK(Contains(result.filtered_sample_ids, "sNULL"));
	CHECK(Contains(result.filtered_sample_ids, "sNEG"));

	// Reached the fit and failed there: constant y (NaN rvalue), a single usable
	// point (all NaN), no usable points at all, and constant x (scipy raises).
	CHECK(result.unfittable_sample_ids.size() == 4);
	for (const char *id : {"sCONSTY", "sONE", "sNONE", "sFLAT"}) {
		INFO("expected " << id << " to be reported unfittable");
		CHECK(Contains(result.unfittable_sample_ids, id));
	}

	CHECK(result.models.size() + result.filtered_sample_ids.size() + result.unfittable_sample_ids.size() == 11);

	// b6 totals 20 reads across all samples, well under 50.
	REQUIRE(result.dropped_syndna_ids.size() == 1);
	CHECK(result.dropped_syndna_ids[0] == "b6");
}

TEST_CASE("FitSyndnaModels does not rescale the mass-fraction denominator", "[absquant]") {
	// pysyndna sums syndna_indiv_ng_ul over the FULL concentrations relation
	// (_calc_indiv_syndna_weights:199) and min_syndna_counts never touches that
	// frame -- so dropping a low-count synDNA removes its data points and
	// leaves every other synDNA's mass alone.
	//
	// Summing over the survivors instead is the obvious implementation, and the
	// reason it needs its own test is that it leaves the SLOPE untouched: every
	// y moves by the same constant, so only the intercept shifts. A test that
	// checked the fit "looked right" would not notice.
	//
	// f4 carries 2.0 of the 3.11 ng/uL total but only one read, so at
	// min_syndna_counts = 5 it is dropped while dominating the denominator.
	std::vector<SyndnaConcentration> concentrations = {{"f1", 1.0}, {"f2", 0.1}, {"f3", 0.01}, {"f4", 2.0}};
	std::vector<CountObservation> counts = {
	    {"s", "f1", 1000.0}, {"s", "f2", 100.0}, {"s", "f3", 10.0}, {"s", "f4", 1.0}};
	std::vector<SampleMass> masses = {{"s", 3.11}};
	FitOptions options;
	options.min_syndna_counts = 5;

	const auto result = FitSyndnaModels(counts, concentrations, masses, options);
	REQUIRE(result.models.size() == 1);
	REQUIRE(result.dropped_syndna_ids.size() == 1);
	CHECK(result.dropped_syndna_ids[0] == "f4");

	// With the full denominator the masses come back out as exactly the
	// concentrations, so y = log10(conc) = 0, -1, -2 against x = 3, 2, 1.
	const auto &fit = result.models[0].fit;
	CHECK(CloseEnough(fit.slope, 1.0, 1e-12, 1e-13));
	CHECK(CloseEnough(fit.intercept, -3.0, 1e-12, 1e-12));

	// Had the denominator been rescaled to the survivors (1.11 rather than
	// 3.11), scipy returns slope 1.0 -- identical -- and intercept
	// -2.55256258975982. Asserting the miss explicitly keeps the discrimination
	// from being accidental.
	CHECK(std::fabs(fit.intercept - (-2.55256258975982)) > 0.4);
}

TEST_CASE("FitSyndnaModels drops low-count synDNAs before filtering samples", "[absquant]") {
	// pysyndna computes syndnas_to_drop at fit_syndna_models.py:492-495, which
	// runs BEFORE the bad-parameter filter at :516-525. So reads contributed by
	// a sample that is about to be thrown away still count toward whether a
	// synDNA survives. Reversing the two is an easy and invisible mistake: it
	// changes no field of any model, it changes how many points a model is fit
	// from.
	//
	// f3 has 5 reads in the usable sample and 45 in a sample with a negative
	// mass. min_syndna_counts = 50 keeps it only if the doomed sample's reads
	// are counted.
	const std::vector<SyndnaConcentration> concentrations = {{"f1", 1.0}, {"f2", 0.1}, {"f3", 0.01}};
	const std::vector<CountObservation> good = {{"g", "f1", 1000.0}, {"g", "f2", 100.0}, {"g", "f3", 5.0}};
	FitOptions options;
	options.min_syndna_counts = 50;

	std::vector<CountObservation> with_doomed = good;
	with_doomed.push_back({"bad", "f3", 45.0});
	const std::vector<SampleMass> two_masses = {{"g", 1.11}, {"bad", -1.0}};
	const auto kept = FitSyndnaModels(with_doomed, concentrations, two_masses, options);

	REQUIRE(kept.models.size() == 1);
	CHECK(kept.models[0].sample_id == "g");
	CHECK(kept.dropped_syndna_ids.empty());
	REQUIRE(kept.filtered_sample_ids.size() == 1);
	CHECK(kept.filtered_sample_ids[0] == "bad");

	// Three points, and they are not collinear (x = 3, 2, log10(5)), so the fit
	// is imperfect in a way two points never can be.
	CHECK(CloseEnough(kept.models[0].fit.slope, 0.8642454806701557, 1e-12, 0.0));
	CHECK(CloseEnough(kept.models[0].fit.rvalue, 0.9971596598938065, 1e-12, 0.0));
	CHECK(kept.models[0].fit.stderr_ > 1e-6);

	// Remove the doomed sample and f3 falls to 5 reads, below the threshold. The
	// surviving sample is then fit from two points, which is scipy's special
	// case: a perfect line with both standard errors hard-zeroed. That is the
	// answer the wrong ordering produces, and it looks entirely reasonable.
	const std::vector<SampleMass> one_mass = {{"g", 1.11}};
	const auto dropped = FitSyndnaModels(good, concentrations, one_mass, options);
	REQUIRE(dropped.models.size() == 1);
	REQUIRE(dropped.dropped_syndna_ids.size() == 1);
	CHECK(dropped.dropped_syndna_ids[0] == "f3");
	CHECK(CloseEnough(dropped.models[0].fit.slope, 1.0, 1e-12, 1e-13));
	CHECK(dropped.models[0].fit.stderr_ == 0.0);
}

TEST_CASE("FitSyndnaModels is invariant to how the contributing fraction is split", "[absquant]") {
	// syndna_contributing_fraction only ever scales the pool mass, so 2.5x the
	// mass at 40% contributing must give the same models as 1x at 100%. This is
	// the property that says the parameter means what the documentation claims:
	// if it were applied per-synDNA, or after the fraction-of-pool division, or
	// forgotten in one branch, this would break while every golden still passed
	// (all of which were generated at fraction = 1.0).
	const auto fixture = LoadFitFixture("fitb");

	FitOptions plain;
	plain.min_syndna_counts = 50;
	plain.syndna_contributing_fraction = 1.0;
	const auto baseline = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, plain);

	std::vector<SampleMass> heavier = fixture.masses;
	for (auto &mass : heavier) {
		mass.mass_syndna_input_ng *= 2.5;
	}
	FitOptions scaled;
	scaled.min_syndna_counts = 50;
	scaled.syndna_contributing_fraction = 0.4;
	const auto rescaled = FitSyndnaModels(fixture.counts, fixture.concentrations, heavier, scaled);

	REQUIRE(rescaled.models.size() == baseline.models.size());
	REQUIRE(baseline.models.size() == 5);
	for (size_t i = 0; i < baseline.models.size(); ++i) {
		INFO("sample " << baseline.models[i].sample_id);
		REQUIRE(rescaled.models[i].sample_id == baseline.models[i].sample_id);
		const auto &a = baseline.models[i].fit;
		const auto &b = rescaled.models[i].fit;
		CHECK(CloseEnough(b.slope, a.slope, 1e-12, 1e-13));
		CHECK(CloseEnough(b.intercept, a.intercept, 1e-12, 1e-13));
		CHECK(CloseEnough(b.rvalue, a.rvalue, 1e-12, 1e-13));
		CHECK(CloseEnough(b.stderr_, a.stderr_, 1e-9, 1e-13));
		CHECK(CloseEnough(b.intercept_stderr, a.intercept_stderr, 1e-9, 1e-13));
	}

	// NaN stays NaN and negative stays negative under multiplication, so the two
	// rejected samples are rejected in both runs rather than sneaking back in.
	CHECK(rescaled.filtered_sample_ids == baseline.filtered_sample_ids);
}

TEST_CASE("FitSyndnaModels excludes zero counts before taking the log", "[absquant]") {
	// log10(0) is -inf, which would poison the whole sample rather than removing
	// one point. pysyndna masks with `count > 0` before the transform; Set B's
	// sZERO covers this against the golden, and this states the mechanism.
	const std::vector<SyndnaConcentration> concentrations = {{"f1", 1.0}, {"f2", 0.1}, {"f3", 0.01}, {"f4", 0.001}};
	const std::vector<SampleMass> masses = {{"s", 1.111}};
	FitOptions options;
	options.min_syndna_counts = 1;

	const std::vector<CountObservation> with_zero = {
	    {"s", "f1", 1000.0}, {"s", "f2", 100.0}, {"s", "f3", 10.0}, {"s", "f4", 0.0}};
	const auto zeroed = FitSyndnaModels(with_zero, concentrations, masses, options);
	REQUIRE(zeroed.models.size() == 1);
	CHECK(std::isfinite(zeroed.models[0].fit.slope));
	CHECK(std::isfinite(zeroed.models[0].fit.intercept));

	// f4 has zero reads everywhere, so its total is 0 and min_syndna_counts >= 1
	// drops it as well -- the same outcome pysyndna reaches from an all-zero row
	// of its dense table, which is why miint can accept sparse input without
	// needing to tell "zero everywhere" apart from "absent".
	REQUIRE(zeroed.dropped_syndna_ids.size() == 1);
	CHECK(zeroed.dropped_syndna_ids[0] == "f4");

	// Omitting the zero row entirely -- what a sparse long-form relation
	// actually looks like -- must give the identical model.
	const std::vector<CountObservation> without_zero = {{"s", "f1", 1000.0}, {"s", "f2", 100.0}, {"s", "f3", 10.0}};
	const auto omitted = FitSyndnaModels(without_zero, concentrations, masses, options);
	REQUIRE(omitted.models.size() == 1);
	CHECK(omitted.models[0].fit.slope == zeroed.models[0].fit.slope);
	CHECK(omitted.models[0].fit.intercept == zeroed.models[0].fit.intercept);
}

TEST_CASE("FitSyndnaModels does not depend on the order rows arrive in", "[absquant]") {
	// SQL hands rows over in whatever order the scan produced, and none of this
	// is order-invariant for free -- summation order sets the last bits of every
	// mean and covariance. FitSyndnaModels buys the invariance by visiting
	// synDNAs in sorted order within each sample and summing the concentrations
	// in that same order; reversing the input must then not move a single bit,
	// or the same query could disagree with itself across plans.
	//
	// One residue is left uncovered and is deliberate: the per-synDNA count
	// totals are accumulated in row order, so in principle a differently-ordered
	// scan could shift one in its last bit and flip a min_syndna_counts decision
	// sitting exactly on the threshold. Read counts are non-negative integers,
	// and sums of those are exact in a double out to 2^53, so the case cannot
	// arise from real input.
	auto fixture = LoadFitFixture("fitb");
	FitOptions options;
	options.min_syndna_counts = 50;
	const auto forward = FitSyndnaModels(fixture.counts, fixture.concentrations, fixture.masses, options);

	std::vector<CountObservation> reversed_counts(fixture.counts.rbegin(), fixture.counts.rend());
	std::vector<SyndnaConcentration> reversed_concentrations(fixture.concentrations.rbegin(),
	                                                         fixture.concentrations.rend());
	const auto backward = FitSyndnaModels(reversed_counts, reversed_concentrations, fixture.masses, options);

	REQUIRE(backward.models.size() == forward.models.size());
	std::map<std::string, LinregressResult> by_id;
	for (const auto &model : backward.models) {
		by_id[model.sample_id] = model.fit;
	}
	for (const auto &model : forward.models) {
		INFO("sample " << model.sample_id);
		REQUIRE(by_id.count(model.sample_id) == 1);
		const auto &other = by_id[model.sample_id];
		// Bit equality, not a tolerance: this is a determinism claim.
		CHECK(other.slope == model.fit.slope);
		CHECK(other.intercept == model.fit.intercept);
		CHECK(other.rvalue == model.fit.rvalue);
		CHECK(other.pvalue == model.fit.pvalue);
		CHECK(other.stderr_ == model.fit.stderr_);
		CHECK(other.intercept_stderr == model.fit.intercept_stderr);
	}
}

TEST_CASE("IsUsableSampleParameter keeps zero and rejects NULL or negative", "[absquant]") {
	// pysyndna's filter tests isna() | lt(0) (util.py:258-266), so ZERO PASSES.
	// That is not a rounding detail: a zero-mass sample survives the parameter
	// filter and is then killed downstream by log10(0), which puts it in a
	// different diagnostic bucket than a NULL-mass sample. Collapsing the two
	// would tell a user "bad parameter" about a parameter that was fine.
	CHECK(IsUsableSampleParameter(0.0));
	CHECK(IsUsableSampleParameter(0.25));
	CHECK(IsUsableSampleParameter(1e300));
	CHECK_FALSE(IsUsableSampleParameter(std::nan("")));
	CHECK_FALSE(IsUsableSampleParameter(-1e-300));
	CHECK_FALSE(IsUsableSampleParameter(-0.3));

	// Infinity is the input that looks usable to every naive test: it is not
	// NaN and it is >= 0, so pysyndna accepts it. It then makes every log10 in
	// the sample infinite, whose differences are NaN, and the sample is
	// discarded as "unfittable" -- reported, but as though a fit had been tried
	// and failed rather than as a parameter that was never usable. This helper
	// is shared with M3-M5, so the gap would otherwise reach four functions.
	CHECK_FALSE(IsUsableSampleParameter(std::numeric_limits<double>::infinity()));
	CHECK_FALSE(IsUsableSampleParameter(-std::numeric_limits<double>::infinity()));
}

TEST_CASE("The shared id predicates", "[absquant]") {
	// M3-M5 each read a long-form counts relation plus keyed reference
	// relations, so these three are written once here rather than four times.
	SECTION("IdsMissingFrom deduplicates and sorts") {
		const std::vector<std::string> subject = {"c", "a", "b", "a", "c"};
		const std::vector<std::string> reference = {"b", "b"};
		const std::vector<std::string> expected = {"a", "c"};
		CHECK(IdsMissingFrom(subject, reference) == expected);
		CHECK(IdsMissingFrom(subject, subject).empty());
		CHECK(IdsMissingFrom({}, reference).empty());
		// Everything missing, not "nothing to compare against".
		CHECK(IdsMissingFrom({"z"}, {}) == std::vector<std::string> {"z"});
	}

	SECTION("DuplicatedIds reports each repeat once") {
		CHECK(DuplicatedIds({"a", "b", "a", "c", "a"}) == std::vector<std::string> {"a"});
		CHECK(DuplicatedIds({"b", "a", "b", "a"}) == std::vector<std::string> {"a", "b"});
		CHECK(DuplicatedIds({"a", "b", "c"}).empty());
		CHECK(DuplicatedIds({}).empty());
	}

	SECTION("RejectDuplicateCells keys on the PAIR and names the relation") {
		// Repeating either id alone is legal -- that is what a long-form table
		// looks like. Only the pair repeating is an error, so a check keyed on
		// one column would reject every well-formed relation.
		const std::vector<CountObservation> distinct = {{"s1", "f1", 1.0}, {"s1", "f2", 2.0}, {"s2", "f1", 3.0}};
		CHECK_NOTHROW(RejectDuplicateCells(distinct, "the counts relation"));
		CHECK_NOTHROW(RejectDuplicateCells({}, "the counts relation"));

		const std::vector<CountObservation> repeated = {{"s1", "f1", 1.0}, {"s2", "f2", 2.0}, {"s1", "f1", 3.0}};
		CHECK_THROWS_AS(RejectDuplicateCells(repeated, "the counts relation"), std::invalid_argument);

		// The relation name is the only thing the three call sites vary, so a
		// splice that dropped it would still throw and still pass a bare
		// CHECK_THROWS. Both spellings are pinned.
		try {
			RejectDuplicateCells(repeated, "the synDNA counts relation");
			FAIL("expected a throw");
		} catch (const std::invalid_argument &e) {
			CHECK(std::string(e.what()) ==
			      "the synDNA counts relation has more than one row for the same (sample_id, feature_id): 's1 / f1'");
		}
		try {
			RejectDuplicateCells(repeated, "the counts relation");
			FAIL("expected a throw");
		} catch (const std::invalid_argument &e) {
			CHECK(std::string(e.what()) ==
			      "the counts relation has more than one row for the same (sample_id, feature_id): 's1 / f1'");
		}

		// Three occurrences of one cell arrive as two entries; FormatIdList
		// renders it once.
		const std::vector<CountObservation> thrice = {{"s1", "f1", 1.0}, {"s1", "f1", 2.0}, {"s1", "f1", 3.0}};
		try {
			RejectDuplicateCells(thrice, "the counts relation");
			FAIL("expected a throw");
		} catch (const std::invalid_argument &e) {
			CHECK(std::string(e.what()) ==
			      "the counts relation has more than one row for the same (sample_id, feature_id): 's1 / f1'");
		}
	}

	SECTION("FormatIdList sorts, deduplicates, quotes, and caps at ten") {
		CHECK(FormatIdList({"b", "a"}) == "'a', 'b'");
		CHECK(FormatIdList({}).empty());

		// Callers that collect offenders row by row hand over one entry per
		// offending ROW, so a cell repeated three times arrives twice. Without
		// dedup that renders as "'x', 'x'" and eats two of the ten slots.
		CHECK(FormatIdList({"x", "x", "x"}) == "'x'");
		CHECK(FormatIdList({"b", "a", "b"}) == "'a', 'b'");

		// A relation with thousands of bad ids must not produce a thousand-id
		// error string, so the tail is a count rather than the rest of the list.
		std::vector<std::string> many;
		for (int i = 1; i <= 13; ++i) {
			many.push_back(i < 10 ? "id0" + std::to_string(i) : "id" + std::to_string(i));
		}
		const std::string rendered = FormatIdList(many);
		CHECK(rendered.find("'id01'") == 0);
		CHECK(rendered.find("'id10'") != std::string::npos);
		CHECK(rendered.find("'id11'") == std::string::npos);
		CHECK(rendered.find(", ... (3 more)") != std::string::npos);
	}
}

TEST_CASE("FitSyndnaModels accepts its own minimal fixture", "[absquant]") {
	// The control for every rejection test below: if this ever stops producing
	// a model, those tests would be passing for the wrong reason.
	const TinyFit fixture;
	const auto result = fixture.Run();
	REQUIRE(result.models.size() == 1);
	CHECK(result.models[0].sample_id == "s1");
	CHECK(CloseEnough(result.models[0].fit.slope, 1.0, 1e-12, 1e-13));
	CHECK(result.dropped_syndna_ids.empty());
	CHECK(result.filtered_sample_ids.empty());
	CHECK(result.unfittable_sample_ids.empty());
	CHECK(result.samples_without_counts.empty());
}

TEST_CASE("FitSyndnaModels rejects out-of-range options", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;

	SECTION("syndna_contributing_fraction must be in (0, 1]") {
		TinyFit fixture;
		fixture.options.syndna_contributing_fraction = 1.0;
		CHECK(fixture.Run().models.size() == 1); // the closed upper bound is legal
		fixture.options.syndna_contributing_fraction = 0.4;
		CHECK(fixture.Run().models.size() == 1);

		for (double bad : {0.0, -0.5, 1.0000001, 2.0}) {
			fixture.options.syndna_contributing_fraction = bad;
			INFO("fraction " << bad);
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		}
		// NaN must be rejected too, which a naive `<= 0 || > 1` test would let
		// through -- every comparison against NaN is false.
		fixture.options.syndna_contributing_fraction = std::nan("");
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("syndna_contributing_fraction"));
	}

	SECTION("min_syndna_counts must be at least 1") {
		// Not cosmetic: >= 1 is what makes a synDNA with zero total reads always
		// dropped, which is why miint can accept sparse input that cannot tell
		// "zero everywhere" from "absent". At 0 that guarantee disappears.
		TinyFit fixture;
		fixture.options.min_syndna_counts = 1;
		CHECK(fixture.Run().models.size() == 1);
		for (int64_t bad : {int64_t(0), int64_t(-1)}) {
			fixture.options.min_syndna_counts = bad;
			INFO("min_syndna_counts " << bad);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("min_syndna_counts"));
		}
	}
}

TEST_CASE("FitSyndnaModels rejects malformed relations", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;

	SECTION("a repeated (sample_id, feature_id) cell") {
		// Left unchecked this becomes a second data point and silently
		// double-weights that synDNA. pysyndna cannot express it at all -- its
		// input is a dense table whose cells are unique by construction.
		TinyFit fixture;
		fixture.counts.push_back({"s1", "f2", 100.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f2'"));
		// The same feature in a DIFFERENT sample is ordinary input, not a repeat.
		TinyFit two_samples;
		two_samples.counts.push_back({"s2", "f2", 100.0});
		two_samples.masses.push_back({"s2", 1.11});
		CHECK_NOTHROW(two_samples.Run());
	}

	SECTION("a repeated concentrations feature_id") {
		// Not merely an ambiguous value to tie-break: pysyndna's merge fans the
		// count rows out once per duplicate key AND double-counts the repeat
		// inside the mass-fraction denominator, so there is no "compatible"
		// answer to pick.
		TinyFit fixture;
		fixture.concentrations.push_back({"f2", 0.1});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f2'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("concentrations"));
	}

	SECTION("a repeated parameters sample_id") {
		TinyFit fixture;
		fixture.masses.push_back({"s1", 2.22});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
	}

	SECTION("a negative read count") {
		// pysyndna lets this become NaN under log10 and silently voids the whole
		// sample. A count below zero is structurally impossible rather than
		// degenerate-but-computable, so miint refuses it -- a deliberate
		// divergence, consistent with D23.
		TinyFit fixture;
		fixture.counts.push_back({"s1", "f4", -5.0});
		fixture.concentrations.push_back({"f4", 0.001});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f4'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("negative"));

		// Zero is NOT negative and must stay accepted -- it is the ordinary way
		// a dense count table spells "no reads", and dropping it is a filtering
		// step, not a rejection. Asserted as "the model is unchanged" rather
		// than merely "does not throw": the point is that the zero row is
		// filtered out of the regression, and only comparing the fit can say so.
		TinyFit with_zero;
		with_zero.counts.push_back({"s1", "f4", 0.0});
		with_zero.concentrations.push_back({"f4", 0.001});
		FitResult zero_result;
		REQUIRE_NOTHROW(zero_result = with_zero.Run());
		REQUIRE(zero_result.models.size() == 1);
		CHECK(zero_result.unfittable_sample_ids.empty());
		// f4 contributes no reads at all, so min_syndna_counts drops it too.
		CHECK(zero_result.dropped_syndna_ids == std::vector<std::string> {"f4"});
	}

	SECTION("a NULL or negative concentration") {
		// pandas' sum() would skip a NaN and carry it into that synDNA's
		// fraction, voiding every sample the synDNA appears in with no
		// explanation. The pool composition is configuration, not measurement.
		TinyFit null_conc;
		null_conc.concentrations[1].syndna_indiv_ng_ul = std::nan("");
		CHECK_THROWS_WITH(null_conc.Run(), ContainsSubstring("'f2'"));
		CHECK_THROWS_WITH(null_conc.Run(), ContainsSubstring("syndna_indiv_ng_ul"));

		TinyFit negative_conc;
		negative_conc.concentrations[2].syndna_indiv_ng_ul = -0.01;
		CHECK_THROWS_WITH(negative_conc.Run(), ContainsSubstring("'f3'"));
	}

	SECTION("concentrations that sum to zero") {
		// Every fraction would be 0/0. A zero for a single synDNA is legal (it
		// contributes no mass); all of them together is not a pool.
		TinyFit fixture;
		for (auto &row : fixture.concentrations) {
			row.syndna_indiv_ng_ul = 0.0;
		}
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("positive and finite"));

		// One zero among several is legal input, but it is not harmless, and
		// "does not throw" would hide what actually happens. f3 still carries
		// counts, so it survives min_syndna_counts and reaches the regression
		// with a mass of exactly zero -- log10(0) is -inf, whose difference from
		// the (also -inf) mean is NaN, which poisons the covariances and voids
		// the sample. That outcome is the assertion: the run is quiet, but the
		// sample lands in unfittable_sample_ids rather than producing a model
		// from an infinite point. A future change that "optimized" away the
		// zero-delta path could silently emit a finite model here instead.
		TinyFit one_zero;
		one_zero.concentrations[2].syndna_indiv_ng_ul = 0.0;
		FitResult zero_conc;
		REQUIRE_NOTHROW(zero_conc = one_zero.Run());
		CHECK(zero_conc.models.empty());
		CHECK(zero_conc.unfittable_sample_ids == std::vector<std::string> {"s1"});
	}
}

TEST_CASE("FitSyndnaModels enforces id consistency asymmetrically", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;

	SECTION("a synDNA with counts but no concentration is an error") {
		// Its mass is uncomputable. pysyndna raises here too.
		TinyFit fixture;
		fixture.counts.push_back({"s1", "f9", 500.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f9'"));
	}

	SECTION("a synDNA with a concentration but no counts is NOT an error") {
		// THE deliberate divergence. pysyndna raises, because its dense read
		// table still contains a synDNA that sequenced nowhere -- but a
		// long-form sparse relation cannot tell "zero in every sample" from
		// "absent", so enforcing this direction would reject legitimate input.
		//
		// The results are identical anyway: min_syndna_counts is required to be
		// >= 1, so a zero-count synDNA is dropped by pysyndna as well. Only the
		// error behavior differs, which is what makes the divergence safe.
		TinyFit fixture;
		fixture.concentrations.push_back({"f9", 0.001});
		const auto result = fixture.Run();
		REQUIRE(result.models.size() == 1);
		CHECK(result.dropped_syndna_ids == std::vector<std::string> {"f9"});
	}

	SECTION("a sample with counts but no parameter row is an error") {
		// Its reads cannot be converted to mass. Reporting it as
		// "parameter-filtered" instead would tell the user a value was bad when
		// the row was simply absent, which sends them looking in the wrong place.
		TinyFit fixture;
		fixture.counts.push_back({"s9", "f1", 900.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s9'"));
	}

	SECTION("a sample with a parameter row but no counts is NOT an error") {
		// A parameters relation covering a whole run, queried for one plate, is
		// ordinary usage -- but it is also exactly what a join on the wrong key
		// looks like, so it is reported rather than silently ignored. pysyndna
		// logs the same list.
		TinyFit fixture;
		fixture.masses.push_back({"s9", 1.11});
		const auto result = fixture.Run();
		REQUIRE(result.models.size() == 1);
		CHECK(result.samples_without_counts == std::vector<std::string> {"s9"});
		CHECK(result.filtered_sample_ids.empty());
		CHECK(result.unfittable_sample_ids.empty());
	}
}

TEST_CASE("FitSyndnaModels separates a bad parameter from an unusable one", "[absquant]") {
	// A NULL mass says "skip this sample"; a zero mass is a real measurement
	// that survives the filter and then dies on log10(0). Both end with no
	// model, and reporting them the same way would misdirect the user, so the
	// two buckets have to stay distinct.
	TinyFit null_mass;
	null_mass.masses[0].mass_syndna_input_ng = std::nan("");
	const auto filtered = null_mass.Run();
	CHECK(filtered.models.empty());
	CHECK(filtered.filtered_sample_ids == std::vector<std::string> {"s1"});
	CHECK(filtered.unfittable_sample_ids.empty());

	TinyFit zero_mass;
	zero_mass.masses[0].mass_syndna_input_ng = 0.0;
	const auto unfittable = zero_mass.Run();
	CHECK(unfittable.models.empty());
	CHECK(unfittable.filtered_sample_ids.empty());
	CHECK(unfittable.unfittable_sample_ids == std::vector<std::string> {"s1"});
}

TEST_CASE("FitSyndnaModels checks the call before the data", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// With several things wrong at once, the first error reported should be the
	// most basic one. A user who passed min_syndna_counts = 0 should hear about
	// that rather than about a duplicate row they may well have on purpose.
	TinyFit fixture;
	fixture.options.min_syndna_counts = 0;
	fixture.counts.push_back({"s1", "f2", 100.0});
	fixture.masses.push_back({"s1", 2.22});
	CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("min_syndna_counts"));
}

TEST_CASE("FitSyndnaModels rejects non-finite values in every relation", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	constexpr double kInf = std::numeric_limits<double>::infinity();

	// Infinity satisfies "not NaN" and "not negative", so it slips past every
	// check written the obvious way -- and then does more damage than any value
	// those checks do catch. None of these produces a wrong number; each
	// produces NO answer, silently, which is the outcome this whole layer
	// exists to replace with a sentence the user can act on.

	SECTION("an infinite concentration destroys the entire run, not one sample") {
		// It becomes the whole denominator, so finite/inf collapses EVERY other
		// synDNA's fraction to exactly 0.0, in every sample at once. y = log10(0)
		// is -inf for every point, the deviations are NaN, and the result is an
		// empty models list for a dataset that is otherwise perfectly good.
		TinyFit fixture;
		fixture.concentrations.push_back({"f4", kInf});
		fixture.counts.push_back({"s1", "f4", 10.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f4'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("finite"));
	}

	SECTION("concentrations that overflow to an infinite total") {
		// Every row is finite here; only the sum is not. Checking the rows alone
		// would leave the same collapse-every-fraction-to-zero outcome.
		TinyFit fixture;
		fixture.concentrations = {{"f1", 1e308}, {"f2", 1e308}, {"f3", 1e308}};
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("positive and finite"));
	}

	SECTION("an infinite read count") {
		TinyFit fixture;
		fixture.counts.push_back({"s1", "f4", kInf});
		fixture.concentrations.push_back({"f4", 0.001});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f4'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("finite"));
	}

	SECTION("a NaN read count, which would make its synDNA undroppable") {
		// The subtle one. The per-point filter excludes NaN, so the fit itself
		// looks fine -- but the min_syndna_counts total is a plain sum, and one
		// NaN makes it NaN, whereupon `NaN < min_syndna_counts` is FALSE and the
		// synDNA can never be dropped however few reads it truly has. pandas
		// sums with skipna=True and does not have this hole. Rejecting every
		// non-finite count closes it at the source rather than by remembering to
		// skip NaN at each later summation.
		TinyFit fixture;
		fixture.counts.push_back({"s1", "f4", std::nan("")});
		fixture.concentrations.push_back({"f4", 0.001});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f4'"));
	}

	SECTION("an infinite sample mass") {
		TinyFit fixture;
		fixture.masses[0].mass_syndna_input_ng = kInf;
		const auto result = fixture.Run();
		// Not an error: pysyndna filters a bad mass rather than raising, and an
		// infinite one is unusable in exactly the way NULL and negative are. It
		// must land in the bucket that says "the parameter was unusable", not
		// the one that says "the fit was attempted and failed".
		CHECK(result.models.empty());
		CHECK(result.filtered_sample_ids == std::vector<std::string> {"s1"});
		CHECK(result.unfittable_sample_ids.empty());
	}
}

TEST_CASE("FitSyndnaModels names a thrice-repeated cell once", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// The duplicate scan reports every occurrence past the first, so three
	// copies of one cell yield two entries. The message must still name it once
	// -- otherwise a handful of distinct problems can trip the ten-id cap and
	// the rest go unreported.
	TinyFit fixture;
	fixture.counts.push_back({"s1", "f2", 100.0});
	fixture.counts.push_back({"s1", "f2", 100.0});
	CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f2'"));
	CHECK_THROWS_WITH(fixture.Run(), !ContainsSubstring("'s1 / f2', 's1 / f2'"));
}

TEST_CASE("FitSyndnaModels reports a malformed relation before an id mismatch", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// The header claims three tiers -- bad call, then malformed relation, then
	// relations that disagree. The bad-call boundary is pinned above; this pins
	// the other one. A user whose counts relation has a duplicate row AND names
	// an unconfigured synDNA should hear about the duplicate, because a
	// malformed relation is often the cause of the mismatch rather than a
	// second, separate mistake.
	TinyFit fixture;
	fixture.counts.push_back({"s1", "f2", 100.0}); // duplicate cell
	fixture.counts.push_back({"s1", "f9", 500.0}); // no concentration row
	CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f2'"));
	CHECK_THROWS_WITH(fixture.Run(), !ContainsSubstring("'f9'"));
}

// ---------------------------------------------------------------------------
// Cell counts
// ---------------------------------------------------------------------------

namespace {

// The smallest input ComputeCellCounts accepts, for pinning one algebraic
// relationship at a time. slope 1 / intercept 0 makes the predicted mass equal
// the read count exactly, so every value below is traceable by hand.
struct TinyCells {
	std::vector<CountObservation> counts = {{"s1", "f1", 100.0}, {"s1", "f2", 1000.0}};
	std::vector<SampleRegression> models = {{"s1", 1.0, 0.0, 1.0}};
	std::vector<CoverageObservation> coverage = {{"s1", "f1", 0.5}, {"s1", "f2", 0.5}};
	std::vector<FeatureLength> lengths = {{"f1", 1000.0}, {"f2", 2000.0}};
	// 2.0 ng/ul over 100 ul is 2e-7 g of extracted gDNA, so against a 0.5-unit
	// denominator the sample-level ratio is exactly 4e-7 -- traceable by hand
	// like every other number in this fixture. Unread unless a test asks for a
	// sample-level metric.
	std::vector<SampleCellParams> params = {{"s1", 10.0, 2.0, 100.0, 0.5}};
	CellCountsOptions options;

	TinyCells() {
		options.min_coverage = 0.1;
		options.min_rsquared = 0.8;
	}
	CellCountsResult Run() const {
		return ComputeCellCounts(counts, models, coverage, lengths, params, options);
	}
	// Value of one cell, or NaN if it was dropped.
	double ValueOf(const std::string &sample_id, const std::string &feature_id) const {
		for (const auto &cell : Run().values) {
			if (cell.sample_id == sample_id && cell.feature_id == feature_id) {
				return cell.value;
			}
		}
		return std::nan("");
	}
};

} // namespace

TEST_CASE("The metric table is consistent in both directions", "[absquant]") {
	// MetricName is the SQL surface and DenominatorColumnName is the schema the
	// params relation must satisfy, so the two together ARE the user-visible
	// contract. Walking every enumerator is what makes a fifth metric that
	// forgets one of the switches fail here rather than at a user's query.
	// Taken from the header's own array rather than restated, so a fifth metric
	// added there is exercised here without anyone remembering to add it.
	const std::vector<CellCountsMetric> all(std::begin(miint::absquant::kAllCellCountsMetrics),
	                                        std::end(miint::absquant::kAllCellCountsMetrics));
	REQUIRE(all.size() == 4);

	SECTION("every name parses back to the metric it came from") {
		std::set<std::string> names;
		for (const auto metric : all) {
			CellCountsMetric parsed = CellCountsMetric::CellsPerGOfGdna;
			INFO("metric " << MetricName(metric));
			REQUIRE(ParseCellCountsMetric(MetricName(metric), parsed));
			CHECK(parsed == metric);
			names.insert(MetricName(metric));
		}
		// Two enumerators sharing a name would round-trip one of them to the
		// other and the loop above would not notice.
		CHECK(names.size() == all.size());
	}

	SECTION("parsing is case-insensitive, as the SQL argument is") {
		CellCountsMetric parsed = CellCountsMetric::CellsPerGOfGdna;
		REQUIRE(ParseCellCountsMetric("Cells_Per_UL_Of_Sample", parsed));
		CHECK(parsed == CellCountsMetric::CellsPerUlOfSample);
	}

	SECTION("anything else is rejected without touching the output") {
		CellCountsMetric parsed = CellCountsMetric::CellsPerCm2OfSample;
		CHECK_FALSE(ParseCellCountsMetric("cells_per_mg_of_sample", parsed));
		CHECK_FALSE(ParseCellCountsMetric("", parsed));
		CHECK_FALSE(ParseCellCountsMetric("cells_per_g_of_gdna ", parsed));
		CHECK(parsed == CellCountsMetric::CellsPerCm2OfSample);
	}

	SECTION("only cells_per_g_of_gdna has no denominator column") {
		CHECK(DenominatorColumnName(CellCountsMetric::CellsPerGOfGdna) == nullptr);
		std::set<std::string> columns;
		for (const auto metric : all) {
			if (metric == CellCountsMetric::CellsPerGOfGdna) {
				continue;
			}
			const char *column = DenominatorColumnName(metric);
			INFO("metric " << MetricName(metric));
			REQUIRE(column != nullptr);
			// Two metrics dividing by the same column would make them the same
			// metric under different names.
			CHECK(columns.insert(column).second);
		}
		// The names pysyndna's SAMPLE_LEVEL_METRICS_DICT keys on, which are also
		// the column names in data/syndna/cells_params.csv.
		CHECK(std::string(DenominatorColumnName(CellCountsMetric::CellsPerGOfSample)) ==
		      "calc_mass_sample_aliquot_input_g");
		CHECK(std::string(DenominatorColumnName(CellCountsMetric::CellsPerUlOfSample)) == "sample_volume_ul");
		CHECK(std::string(DenominatorColumnName(CellCountsMetric::CellsPerCm2OfSample)) == "sample_surface_area_cm2");
	}
}

TEST_CASE("ComputeCellCounts reproduces pysyndna's cells_per_g_of_gdna", "[absquant]") {
	// Set A: pysyndna's own fixture, traceable to the published spreadsheet and
	// the Zaramela notebook. min_coverage 0.1 and min_rsquared 0.8 are the
	// values data/syndna/README.md records the oracle was generated with.
	auto fixture = LoadCellsFixture("cells");
	fixture.options.min_coverage = 0.1;
	fixture.options.min_rsquared = 0.8;

	const auto result = fixture.Run();
	CheckCellsAgainstOracle(result, "data/syndna/cells_oracle.csv", "cells_per_g_of_gdna");

	// 13 OGUs x 2 samples, less the two that fail coverage in both samples,
	// less the (example2, Neisseria subflava) cell whose count is zero -- which
	// is one of those already-dropped four, so 26 - 4 = 22.
	CHECK(result.values.size() == 22);
}

TEST_CASE("ComputeCellCounts drops features below min_coverage", "[absquant]") {
	// pysyndna's fixture engineers these two to fail in every sample so the
	// filter is always exercised: N. subflava at 7.9% and H. influenzae at 1.4%
	// coverage, against a 10% threshold.
	auto fixture = LoadCellsFixture("cells");
	fixture.options.min_coverage = 0.1;
	const auto result = fixture.Run();

	CHECK(Contains(result.low_coverage_feature_ids, "Neisseria subflava"));
	CHECK(Contains(result.low_coverage_feature_ids, "Haemophilus influenzae"));
	CHECK(result.low_coverage_feature_ids.size() == 2);
	for (const auto &cell : result.values) {
		CHECK(cell.feature_id != "Neisseria subflava");
		CHECK(cell.feature_id != "Haemophilus influenzae");
	}

	// At a threshold below both, nothing is dropped and all 25 nonzero cells
	// survive -- so the 22 above really is the filter's doing, not an artifact
	// of the fixture.
	fixture.options.min_coverage = 0.0;
	const auto unfiltered = fixture.Run();
	CHECK(unfiltered.low_coverage_feature_ids.empty());
	CHECK(unfiltered.values.size() == 25);
}

TEST_CASE("ComputeCellCounts min_coverage is a strict less-than", "[absquant]") {
	// A feature sitting EXACTLY on the threshold is kept. pysyndna filters with
	// `coverage >= min_coverage` (calc_cell_counts.py:602), and every ordinary
	// fixture lands clearly on one side, so without this a `<=` bug would pass
	// everything else.
	TinyCells fixture;
	fixture.coverage[0].coverage = 0.25;
	fixture.options.min_coverage = 0.25;
	CHECK_FALSE(std::isnan(fixture.ValueOf("s1", "f1")));
	CHECK(fixture.Run().low_coverage_feature_ids.empty());

	// The smallest representable step above the threshold drops it.
	fixture.options.min_coverage = std::nextafter(0.25, 1.0);
	CHECK(std::isnan(fixture.ValueOf("s1", "f1")));
	CHECK(Contains(fixture.Run().low_coverage_feature_ids, "f1"));
}

TEST_CASE("ComputeCellCounts scales as the method requires", "[absquant]") {
	// The oracle pins the whole chain at once; these pin each factor's role
	// separately, so a wrong constant or an inverted ratio is localized rather
	// than showing up as "all 22 values differ".
	const TinyCells base;
	const double v = base.ValueOf("s1", "f1");
	REQUIRE(std::isfinite(v));
	REQUIRE(v > 0.0);

	// Genome length is a denominator: twice the genome, half the genomes.
	TinyCells longer;
	longer.lengths[0].ogu_len_in_bp *= 2.0;
	CHECK(CloseEnough(longer.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));

	// Normalizing by the sequenced gDNA mass: twice the gDNA, half the density.
	TinyCells heavier;
	heavier.params[0].sequenced_sample_gdna_mass_ng *= 2.0;
	CHECK(CloseEnough(heavier.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));

	// The model acts through log10, so +1 on the intercept is x10 on the mass.
	TinyCells shifted;
	shifted.models[0].intercept += 1.0;
	CHECK(CloseEnough(shifted.ValueOf("s1", "f1"), v * 10.0, 1e-14, 0.0));

	// ... and the slope acts on log10(count): at slope 1 the mass IS the count,
	// so f2's ten-fold count with f1's length would be ten-fold the value.
	TinyCells same_length;
	same_length.lengths[1].ogu_len_in_bp = same_length.lengths[0].ogu_len_in_bp;
	CHECK(CloseEnough(same_length.ValueOf("s1", "f2"), v * 10.0, 1e-14, 0.0));

	// Every cell is per-sample: a second sample with its own model and gDNA mass
	// must not disturb the first.
	TinyCells two_samples;
	two_samples.counts.push_back({"s2", "f1", 100.0});
	two_samples.models.push_back({"s2", 1.0, 1.0, 1.0});
	two_samples.coverage.push_back({"s2", "f1", 0.5});
	two_samples.params.push_back({"s2", 20.0, 2.0, 100.0});
	CHECK(CloseEnough(two_samples.ValueOf("s1", "f1"), v, 1e-14, 0.0));
	CHECK(CloseEnough(two_samples.ValueOf("s2", "f1"), v * 10.0 / 2.0, 1e-14, 0.0));
}

TEST_CASE("ComputeCellCounts reproduces pysyndna's sample-level metrics", "[absquant]") {
	// Same Set A, same thresholds, same 22 surviving cells: the sample-level
	// metrics change only the unit, so any one of them dropping or inventing a
	// cell would mean the metric had leaked into the filters.
	for (const auto metric : {CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
	                          CellCountsMetric::CellsPerCm2OfSample}) {
		INFO("metric " << MetricName(metric));
		auto fixture = LoadCellsFixture("cells", metric);
		fixture.options.min_coverage = 0.1;
		fixture.options.min_rsquared = 0.8;

		const auto result = fixture.Run();
		CheckCellsAgainstOracle(result, "data/syndna/cells_oracle.csv", MetricName(metric));
		CHECK(result.values.size() == 22);
	}
}

TEST_CASE("ComputeCellCounts converts to sample units by one per-sample ratio", "[absquant]") {
	// The oracle pins all four metrics at once; this pins the ratio's shape, so
	// an inverted division or a denominator read from the wrong column is
	// localized instead of showing up as "all 22 values differ".
	const TinyCells base;
	const double gdna = base.ValueOf("s1", "f1");
	REQUIRE(std::isfinite(gdna));
	REQUIRE(gdna > 0.0);

	// 2.0 ng/ul * 100 ul / 1e9 = 2e-7 g extracted, over a 0.5-unit denominator.
	// Written as the formula rather than as 4e-7 so the test cannot agree with
	// an arithmetic mistake it shares with the code.
	const double extracted_gdna_mass_g =
	    base.params[0].extracted_gdna_concentration_ng_ul * base.params[0].vol_extracted_elution_ul / 1e9;
	const double ratio = extracted_gdna_mass_g / base.params[0].sample_denominator;

	SECTION("each sample-level metric scales cells_per_g_of_gdna by that ratio") {
		for (const auto metric : {CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
		                          CellCountsMetric::CellsPerCm2OfSample}) {
			INFO("metric " << MetricName(metric));
			TinyCells fixture;
			fixture.options.metric = metric;
			// All three read `sample_denominator`, which the caller has already
			// resolved -- so with the same denominator they must agree exactly.
			// That is what makes the enum a column selector rather than three
			// different formulas.
			CHECK(CloseEnough(fixture.ValueOf("s1", "f1"), gdna * ratio, 1e-14, 0.0));
		}
	}

	SECTION("the denominator divides and the extraction terms multiply") {
		TinyCells bigger;
		bigger.options.metric = CellCountsMetric::CellsPerGOfSample;
		const double v = bigger.ValueOf("s1", "f1");
		REQUIRE(std::isfinite(v));

		// Twice the sample, half the cells per unit of it.
		TinyCells twice_the_sample = bigger;
		twice_the_sample.params[0].sample_denominator *= 2.0;
		CHECK(CloseEnough(twice_the_sample.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));

		// Twice the extracted gDNA, from either factor, is twice the cells --
		// varying them separately is what catches a ratio that reads only one.
		TinyCells twice_the_concentration = bigger;
		twice_the_concentration.params[0].extracted_gdna_concentration_ng_ul *= 2.0;
		CHECK(CloseEnough(twice_the_concentration.ValueOf("s1", "f1"), v * 2.0, 1e-14, 0.0));

		TinyCells twice_the_elution = bigger;
		twice_the_elution.params[0].vol_extracted_elution_ul *= 2.0;
		CHECK(CloseEnough(twice_the_elution.ValueOf("s1", "f1"), v * 2.0, 1e-14, 0.0));
	}

	SECTION("the ratio is per sample, not per feature") {
		// Two features of the same sample must move together: a ratio applied
		// per cell would still pass every single-feature check above.
		TinyCells fixture;
		fixture.options.metric = CellCountsMetric::CellsPerUlOfSample;
		CHECK(CloseEnough(fixture.ValueOf("s1", "f2"), base.ValueOf("s1", "f2") * ratio, 1e-14, 0.0));
	}

	SECTION("cells_per_g_of_gdna is untouched by any of the three inputs") {
		// Exact equality, not a tolerance: the gdna metric must not enter the
		// ratio path at all, so M3's values are reproduced bit for bit rather
		// than merely closely.
		TinyCells fixture;
		fixture.params[0].sample_denominator = 12345.0;
		fixture.params[0].extracted_gdna_concentration_ng_ul = 7.0;
		fixture.params[0].vol_extracted_elution_ul = 3.0;
		CHECK(fixture.ValueOf("s1", "f1") == gdna);
	}
}

TEST_CASE("ComputeCellCounts matches pysyndna once every filter is in play", "[absquant]") {
	// Set B exists to exercise the filters Set A cannot: a sample below the r^2
	// gate, a sample whose every feature fails coverage, a feature that passes
	// in one sample and not another, and a zero read count. 0.1 and 0.8 are the
	// values gen_syndna_oracle.py generated cellsb_oracle.csv with.
	auto fixture = LoadCellsFixture("cellsb");
	fixture.options.min_coverage = 0.1;
	fixture.options.min_rsquared = 0.8;

	const auto result = fixture.Run();
	CheckCellsAgainstOracle(result, "data/syndna/cellsb_oracle.csv", "cells_per_g_of_gdna");

	// 9 dense golden cells, of which two are the zeros D10 says we omit:
	// (s2, f_only1) failed coverage and (s2, f_zero) had no reads.
	CHECK(result.values.size() == 7);
}

TEST_CASE("ComputeCellCounts pins both strict-< thresholds at once", "[absquant]") {
	// The boundary oracle sets min_coverage and min_rsquared EXACTLY to values
	// present in the data: f_only1 has coverage 0.02 in s2, and slowr2's
	// rvalue 0.5 gives r^2 of exactly 0.25. Both are exactly representable, so
	// this is a real equality test rather than a lucky rounding. Everything else
	// in the suite lands clearly on one side of each threshold, which is exactly
	// how a `<=` for either would otherwise ship.
	auto fixture = LoadCellsFixture("cellsb");
	fixture.options.min_coverage = 0.02;
	fixture.options.min_rsquared = 0.25;

	const auto result = fixture.Run();
	CheckCellsAgainstOracle(result, "data/syndna/cellsb_boundary_oracle.csv", "cells_per_g_of_gdna");

	// 12 dense golden cells less the one zero, (s2, f_zero). slowr2 is back
	// (3 cells) and so is (s2, f_only1), against the 7 above.
	CHECK(result.values.size() == 11);
	CHECK(result.low_rsquared_sample_ids.empty());
	CHECK(Contains(result.uncovered_sample_ids, "sallcov"));

	// One representable step up on either threshold and the sample sitting on
	// it is gone -- which is what makes the equality above load-bearing.
	auto rsquared = fixture;
	rsquared.options.min_rsquared = std::nextafter(0.25, 1.0);
	const auto tightened = rsquared.Run();
	CHECK(tightened.low_rsquared_sample_ids == std::vector<std::string> {"slowr2"});
	CHECK(tightened.values.size() == 8);

	auto coverage = fixture;
	coverage.options.min_coverage = std::nextafter(0.02, 1.0);
	const auto narrowed = coverage.Run();
	CHECK(narrowed.values.size() == 10);
	for (const auto &cell : narrowed.values) {
		CHECK_FALSE((cell.sample_id == "s2" && cell.feature_id == "f_only1"));
	}
}

TEST_CASE("ComputeCellCounts screens only the requested metric's denominator", "[absquant]") {
	// pysyndna's filter set is sequenced_sample_gdna_mass_ng, plus
	// REQUIRED_DNA_PREP_INFO_KEYS, plus the requested metric's OWN denominator
	// and no other (calc_ogu_cell_counts_biom:959-965). Set B's `snull` has a
	// NULL calc_mass_sample_aliquot_input_g and nothing else missing, so it is
	// the one sample that tells the three candidate readings apart:
	//
	//   screen all three sample columns -> snull is dropped from every metric
	//   screen none of them             -> snull yields NaN cells for per-g
	//   screen only the requested one   -> the goldens
	//
	// All three reproduce every VALUE of the samples that survive, so only the
	// membership below can distinguish them.
	SECTION("the metric whose denominator is missing loses the sample") {
		auto fixture = LoadCellsFixture("cellsb", CellCountsMetric::CellsPerGOfSample);
		fixture.options.min_coverage = 0.1;
		fixture.options.min_rsquared = 0.8;

		const auto result = fixture.Run();
		CheckCellsAgainstOracle(result, "data/syndna/cellsb_oracle.csv", "cells_per_g_of_sample");

		CHECK(result.filtered_sample_ids == std::vector<std::string> {"snull"});
		for (const auto &cell : result.values) {
			CHECK(cell.sample_id != "snull");
		}
		// 6 dense golden cells less the two zeros, against 9 less two for the
		// metrics that keep snull.
		CHECK(result.values.size() == 4);
	}

	SECTION("the metrics that do not read that column keep it") {
		for (const auto metric : {CellCountsMetric::CellsPerUlOfSample, CellCountsMetric::CellsPerCm2OfSample}) {
			INFO("metric " << MetricName(metric));
			auto fixture = LoadCellsFixture("cellsb", metric);
			fixture.options.min_coverage = 0.1;
			fixture.options.min_rsquared = 0.8;

			const auto result = fixture.Run();
			CheckCellsAgainstOracle(result, "data/syndna/cellsb_oracle.csv", MetricName(metric));

			CHECK(result.filtered_sample_ids.empty());
			CHECK(result.values.size() == 7);
			bool has_snull = false;
			for (const auto &cell : result.values) {
				has_snull = has_snull || cell.sample_id == "snull";
			}
			CHECK(has_snull);
		}
	}

	SECTION("a bad denominator is filtered, not thrown, and only for its metric") {
		// The same three unusable values the other required parameters get, on
		// the one column whose relevance depends on the metric.
		for (const double bad : {std::nan(""), -1.0, std::numeric_limits<double>::infinity()}) {
			INFO("denominator " << bad);
			TinyCells fixture;
			fixture.params[0].sample_denominator = bad;

			fixture.options.metric = CellCountsMetric::CellsPerGOfSample;
			CHECK(fixture.Run().filtered_sample_ids == std::vector<std::string> {"s1"});
			CHECK(fixture.Run().values.empty());

			// cells_per_g_of_gdna has no denominator column at all, so the same
			// row is perfectly usable for it.
			fixture.options.metric = CellCountsMetric::CellsPerGOfGdna;
			CHECK(fixture.Run().filtered_sample_ids.empty());
			CHECK(fixture.Run().values.size() == 2);
		}
	}
}

TEST_CASE("ComputeCellCounts pins the strict-< thresholds for every metric", "[absquant]") {
	// The boundary oracle at min_coverage 0.02 / min_rsquared 0.25, where both
	// thresholds sit exactly on a value in the data. Running it per metric is
	// what keeps the metric out of the filters: a threshold comparison that
	// somehow depended on the unit would show up here and nowhere else.
	for (const auto metric : {CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
	                          CellCountsMetric::CellsPerCm2OfSample}) {
		INFO("metric " << MetricName(metric));
		auto fixture = LoadCellsFixture("cellsb", metric);
		fixture.options.min_coverage = 0.02;
		fixture.options.min_rsquared = 0.25;

		const auto result = fixture.Run();
		CheckCellsAgainstOracle(result, "data/syndna/cellsb_boundary_oracle.csv", MetricName(metric));

		// slowr2 is admitted at exactly 0.25, and sallcov is still wholly
		// uncovered, whichever unit is asked for.
		CHECK(result.low_rsquared_sample_ids.empty());
		CHECK(result.uncovered_sample_ids == std::vector<std::string> {"sallcov"});
		// 9 golden cells less one zero for the metric that drops snull, 12 less
		// one for the two that keep it.
		CHECK(result.values.size() == (metric == CellCountsMetric::CellsPerGOfSample ? 8u : 11u));
	}
}

TEST_CASE("ComputeCellCounts reports a sample whose every cell came out zero", "[absquant]") {
	// A zero extraction ratio is the one way a sample can pass every filter and
	// still contribute nothing. Both factors of extracted_gdna_mass_g may
	// legitimately be zero -- pysyndna's parameter screen tests `< 0`, so zero
	// passes it, and a blank really can extract no DNA -- and then every cell of
	// the sample is exactly 0.0 and the sparse invariant omits all of them.
	//
	// This cannot arise for cells_per_g_of_gdna, whose values are strictly
	// positive, which is why M3 had no such list. Without one the sample would
	// vanish from the output with nothing said about it, and "no rows" would be
	// indistinguishable from "no such sample".
	for (size_t which = 0; which < 2; ++which) {
		INFO("zeroed extraction factor " << which);
		TinyCells fixture;
		fixture.options.metric = CellCountsMetric::CellsPerGOfSample;
		double *factors[] = {&fixture.params[0].extracted_gdna_concentration_ng_ul,
		                     &fixture.params[0].vol_extracted_elution_ul};
		*factors[which] = 0.0;

		const auto result = fixture.Run();
		CHECK(result.values.empty());
		CHECK(result.zero_valued_sample_ids == std::vector<std::string> {"s1"});
		// Named ONCE. A zero parameter is not a missing one, so this must not
		// also report the sample as filtered.
		CHECK(result.filtered_sample_ids.empty());
		CHECK(result.samples_without_models.empty());
		CHECK(result.low_rsquared_sample_ids.empty());
		CHECK(result.uncovered_sample_ids.empty());

		// The identical parameters under cells_per_g_of_gdna are fine: it never
		// forms the ratio, so there is nothing to be zero.
		fixture.options.metric = CellCountsMetric::CellsPerGOfGdna;
		const auto gdna = fixture.Run();
		CHECK(gdna.values.size() == 2);
		CHECK(gdna.zero_valued_sample_ids.empty());
	}

	SECTION("only the affected sample is lost") {
		TinyCells fixture;
		fixture.options.metric = CellCountsMetric::CellsPerUlOfSample;
		fixture.counts.push_back({"s2", "f1", 100.0});
		fixture.models.push_back({"s2", 1.0, 0.0, 1.0});
		fixture.coverage.push_back({"s2", "f1", 0.5});
		fixture.params.push_back({"s2", 10.0, 2.0, 100.0, 0.5});
		fixture.params[0].vol_extracted_elution_ul = 0.0;

		const auto result = fixture.Run();
		CHECK(result.zero_valued_sample_ids == std::vector<std::string> {"s1"});
		REQUIRE(result.values.size() == 1);
		CHECK(result.values[0].sample_id == "s2");
	}
}

TEST_CASE("ComputeCellCounts reports every sample it discards", "[absquant]") {
	// Checked against pysyndna's own log for this exact run, captured when the
	// oracle was generated:
	//
	//   The following items have coverage lower than the minimum of 10.0:
	//       ['f_neither', 'f_only1', 'f_both', 'f_zero']
	//   R^2 of linear regression for sample slowr2 is 0.25, which is less than
	//       the minimum allowed value of 0.8.
	//   No cell counts calculated for sample slowr2
	//
	// Note what is NOT there: nothing about sallcov, whose four features are all
	// below coverage and which therefore produces no output at all. That is the
	// one place miint says more than pysyndna does.
	auto fixture = LoadCellsFixture("cellsb");
	fixture.options.min_coverage = 0.1;
	fixture.options.min_rsquared = 0.8;
	const auto result = fixture.Run();

	CHECK(result.low_coverage_feature_ids == std::vector<std::string> {"f_both", "f_neither", "f_only1", "f_zero"});
	CHECK(result.low_rsquared_sample_ids == std::vector<std::string> {"slowr2"});
	CHECK(result.uncovered_sample_ids == std::vector<std::string> {"sallcov"});
	CHECK(result.filtered_sample_ids.empty());
	CHECK(result.samples_without_models.empty());

	// f_both and f_zero reach that list only through sallcov: they pass coverage
	// in every other sample. So the list really is "failed somewhere", not
	// "failed everywhere", which is what pysyndna's is.
	CHECK(Contains(result.low_coverage_feature_ids, "f_both"));
	for (const auto &cell : result.values) {
		CHECK(cell.sample_id != "sallcov");
		CHECK(cell.sample_id != "slowr2");
	}

	// Total accounting: five samples in, three producing cells and two named.
	// Every list is in the sum, so a sample that starts landing in none of them
	// -- or in two -- fails here rather than silently vanishing.
	std::set<std::string> with_values;
	for (const auto &cell : result.values) {
		with_values.insert(cell.sample_id);
	}
	CHECK(with_values == std::set<std::string> {"s1", "s2", "snull"});
	CHECK(with_values.size() + result.filtered_sample_ids.size() + result.samples_without_models.size() +
	          result.low_rsquared_sample_ids.size() + result.uncovered_sample_ids.size() +
	          result.zero_valued_sample_ids.size() ==
	      5);
}

TEST_CASE("ComputeCellCounts filters bad parameters before it filters coverage", "[absquant]") {
	// The order is invisible in the values and decides the diagnostics. pysyndna
	// drops bad-parameter samples out of the counts table itself before the
	// coverage filter runs, so a feature that is poorly covered ONLY in such a
	// sample is never reported -- rightly, since no surviving sample was going
	// to use it. Get the order wrong and the user is sent after a coverage
	// problem in data that was already discarded for another reason.
	TinyCells fixture;
	fixture.counts.push_back({"s2", "f3", 100.0});
	fixture.models.push_back({"s2", 1.0, 0.0, 1.0});
	fixture.coverage.push_back({"s2", "f3", 0.01}); // below min_coverage 0.1
	fixture.lengths.push_back({"f3", 1000.0});
	fixture.params.push_back({"s2", 10.0, 2.0, 100.0});

	// With s2's parameters intact, f3 fails coverage and is named.
	const auto reported = fixture.Run();
	CHECK(reported.low_coverage_feature_ids == std::vector<std::string> {"f3"});
	CHECK(reported.filtered_sample_ids.empty());

	// The only change is s2's gDNA mass. f3 must now go unmentioned.
	fixture.params[1].sequenced_sample_gdna_mass_ng = std::nan("");
	const auto silent = fixture.Run();
	CHECK(silent.low_coverage_feature_ids.empty());
	CHECK(silent.filtered_sample_ids == std::vector<std::string> {"s2"});
	CHECK(silent.values.size() == 2); // s1's two cells, untouched
}

TEST_CASE("ComputeCellCounts requires the parameters it never divides by", "[absquant]") {
	// Only sequenced_sample_gdna_mass_ng enters this metric's arithmetic. The
	// other two are pysyndna's REQUIRED_DNA_PREP_INFO_KEYS: they feed
	// extracted_gdna_mass_g, which only the sample-level metrics need, and yet
	// filter_data_by_sample_info screens on them for every metric. Dropping the
	// requirement would change which samples appear in the output while every
	// value stayed identical -- a divergence no parity test on values could see,
	// which is exactly why it is pinned here.
	for (size_t which = 0; which < 3; ++which) {
		TinyCells fixture;
		double *fields[] = {&fixture.params[0].sequenced_sample_gdna_mass_ng,
		                    &fixture.params[0].extracted_gdna_concentration_ng_ul,
		                    &fixture.params[0].vol_extracted_elution_ul};
		INFO("required parameter " << which);

		*fields[which] = std::nan("");
		CHECK(fixture.Run().filtered_sample_ids == std::vector<std::string> {"s1"});
		CHECK(fixture.Run().values.empty());

		*fields[which] = -1.0;
		CHECK(fixture.Run().filtered_sample_ids == std::vector<std::string> {"s1"});

		// Infinity passes both of pysyndna's tests and then quietly voids the
		// sample; IsUsableSampleParameter rejects it, same as everywhere else in
		// this port.
		*fields[which] = std::numeric_limits<double>::infinity();
		CHECK(fixture.Run().filtered_sample_ids == std::vector<std::string> {"s1"});
	}

	// Zero is a real measurement, not a missing one: pysyndna tests with a
	// strict `< 0` and this filter must not exceed it. It is only a problem for
	// the one column that is a denominator, and that is the validation layer's
	// business -- refusing it here would report a structurally impossible number
	// as though the metadata were merely absent.
	TinyCells zeroed;
	zeroed.params[0].extracted_gdna_concentration_ng_ul = 0.0;
	zeroed.params[0].vol_extracted_elution_ul = 0.0;
	CHECK(zeroed.Run().filtered_sample_ids.empty());
	CHECK(zeroed.Run().values.size() == 2);
}

TEST_CASE("ComputeCellCounts reports a sample with no model", "[absquant]") {
	// absquant_fit_models returns fewer models than it was given samples
	// whenever a standard curve cannot be fit, so its output fed straight back
	// in will routinely be missing samples the counts relation has. That has to
	// be a warning, not an error, or the two functions would not compose.
	TinyCells fixture;
	fixture.counts.push_back({"s2", "f1", 100.0});
	fixture.coverage.push_back({"s2", "f1", 0.5});
	fixture.params.push_back({"s2", 10.0, 2.0, 100.0});

	const auto result = fixture.Run();
	CHECK(result.samples_without_models == std::vector<std::string> {"s2"});
	CHECK(result.values.size() == 2); // s1's, unaffected
	CHECK(result.low_rsquared_sample_ids.empty());

	// A sample that is BOTH unmodelled and wholly uncovered is reported once,
	// as uncovered. pysyndna reaches the model lookup only for samples still in
	// the working table after the coverage filter, so it says nothing about the
	// model either; the point of the single bucket is that the first reason the
	// data was discarded is the one worth acting on.
	fixture.coverage.back().coverage = 0.01;
	const auto uncovered = fixture.Run();
	CHECK(uncovered.uncovered_sample_ids == std::vector<std::string> {"s2"});
	CHECK(uncovered.samples_without_models.empty());
}

TEST_CASE("ComputeCellCounts accepts its own minimal fixture", "[absquant]") {
	// Guards every negative test below: each breaks one field of TinyCells and
	// expects a throw, which proves nothing if the untouched fixture throws too.
	TinyCells fixture;
	CHECK_NOTHROW(fixture.Run());
	CHECK(fixture.Run().values.size() == 2);
}

TEST_CASE("ComputeCellCounts rejects out-of-range options", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	constexpr double kInf = std::numeric_limits<double>::infinity();

	SECTION("min_coverage outside [0, 1]") {
		// D9: coverage in miint is a FRACTION, so the threshold is one too.
		// pysyndna takes either and explicitly refuses to guard the distinction
		// ("IT IS UP TO THE USER..."), which makes min_coverage = 10 against
		// fractional coverages a filter that drops everything, silently.
		TinyCells fixture;
		for (const double bad : {-0.001, 1.001, 10.0, kInf, std::nan("")}) {
			fixture.options.min_coverage = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("min_coverage"));
		}
		// Both endpoints are legal: 0 keeps everything with a coverage row, and
		// 1 keeps only fully covered features.
		fixture.options.min_coverage = 0.0;
		CHECK_NOTHROW(fixture.Run());
		fixture.options.min_coverage = 1.0;
		CHECK_NOTHROW(fixture.Run());
	}

	SECTION("min_rsquared outside [0, 1]") {
		// r^2 is a squared correlation and cannot leave [0, 1], so a threshold
		// outside it is not a strict filter -- it is one that admits every
		// sample or none, with no indication that is what happened.
		TinyCells fixture;
		for (const double bad : {-0.001, 1.001, kInf, std::nan("")}) {
			fixture.options.min_rsquared = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("min_rsquared"));
		}
		fixture.options.min_rsquared = 0.0;
		CHECK_NOTHROW(fixture.Run());
		fixture.options.min_rsquared = 1.0;
		CHECK_NOTHROW(fixture.Run());
	}
}

TEST_CASE("ComputeCellCounts rejects malformed relations", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	constexpr double kInf = std::numeric_limits<double>::infinity();

	SECTION("a duplicated counts cell") {
		// Every relation here is keyed, and a repeat means the caller's join fanned
		// out. Taking either row silently would give an answer to a question the
		// data cannot actually answer.
		TinyCells fixture;
		fixture.counts.push_back({"s1", "f1", 100.0});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
	}

	SECTION("a duplicated coverage cell") {
		TinyCells fixture;
		fixture.coverage.push_back({"s1", "f2", 0.5});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f2'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("coverage"));
	}

	SECTION("a duplicated feature length") {
		TinyCells fixture;
		fixture.lengths.push_back({"f1", 1000.0});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f1'"));
	}

	SECTION("a duplicated model") {
		TinyCells fixture;
		fixture.models.push_back({"s1", 2.0, 0.0, 1.0});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
	}

	SECTION("a duplicated parameter row") {
		TinyCells fixture;
		fixture.params.push_back({"s1", 10.0, 2.0, 100.0});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
	}

	SECTION("a read count that is not a positive finite number") {
		// Stricter than the fit's rule, which admits zero because it drops those
		// points before taking the log. There is no such step here: log10(0) is
		// -inf, and under a NEGATIVE slope that comes back as 10^(+inf), an
		// INFINITE cell count -- which is what pysyndna emits for that input.
		// The DuckDB reader drops zero-valued cells, so this is unreachable in
		// practice; it is the guard that keeps it that way if the reader ever
		// changes.
		for (const double bad : {0.0, -1.0, kInf, -kInf, std::nan("")}) {
			TinyCells fixture;
			fixture.counts[0].count = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
			// Naming the message, not just the cell: an infinite count reaches
			// the arithmetic as 10^(+inf) and the overflow guard refuses it too,
			// pointing at the model coefficients and the extraction parameters
			// instead. Only the zero case distinguishes the two checks without
			// this, so the non-finite ones would otherwise pass either way.
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("read counts"));
		}
	}

	SECTION("a coverage outside [0, 1]") {
		// The percent-versus-fraction trap, made loud. pysyndna's own fixtures
		// carry coverages like 92.597 and compare them against min_coverage = 10;
		// handing those to miint's fractional threshold would keep every feature
		// and quietly disable the filter. Rejecting > 1 catches it at the door.
		TinyCells fixture;
		fixture.coverage[0].coverage = 92.597;
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("fraction"));

		for (const double bad : {-0.001, 1.001, kInf, std::nan("")}) {
			TinyCells each;
			each.coverage[0].coverage = bad;
			CHECK_THROWS_AS(each.Run(), std::invalid_argument);
		}
		// The endpoints are ordinary data: nothing aligned, or every base did.
		TinyCells edges;
		edges.coverage[0].coverage = 0.0;
		edges.coverage[1].coverage = 1.0;
		edges.options.min_coverage = 0.0;
		CHECK_NOTHROW(edges.Run());
	}

	SECTION("a genome length that is not positive and finite") {
		// A denominator. pysyndna divides by it with no check at all, so a zero
		// length yields inf cells and a NaN one yields NaN -- in both cases for
		// that OGU in every sample, and in both cases silently.
		for (const double bad : {0.0, -1.0, kInf, std::nan("")}) {
			TinyCells fixture;
			fixture.lengths[0].ogu_len_in_bp = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f1'"));
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("ogu_len_in_bp"));
		}
	}

	SECTION("a zero sequenced gDNA mass, which the parameter filter lets past") {
		// The other denominator, and the one case pysyndna's NaN/negative screen
		// cannot see: zero is neither. It divides through to inf cells for every
		// feature in the sample, with no log message -- the same situation D23
		// already refuses for total_biological_reads_r1r2.
		TinyCells fixture;
		fixture.params[0].sequenced_sample_gdna_mass_ng = 0.0;
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("sequenced_sample_gdna_mass_ng"));
	}

	SECTION("a zero sample denominator, for the metric that divides by it") {
		// Same reasoning one level out: the sample-level metrics divide by their
		// own column, pysyndna does so unguarded, and a zero there yields inf
		// cells for every feature in the sample. Zero passes the NaN/negative
		// screen, so this is the only thing standing between a blank metadata
		// cell and an infinite answer.
		for (const auto metric : {CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
		                          CellCountsMetric::CellsPerCm2OfSample}) {
			INFO("metric " << MetricName(metric));
			TinyCells fixture;
			fixture.options.metric = metric;
			fixture.params[0].sample_denominator = 0.0;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
			// The message must name the column THIS metric reads, not a generic
			// "denominator" and not one of the other two. Sending someone to
			// sample_volume_ul to fix a blank calc_mass_sample_aliquot_input_g
			// is worse than saying nothing.
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring(DenominatorColumnName(metric)));
			for (const auto other : {CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
			                         CellCountsMetric::CellsPerCm2OfSample}) {
				if (other != metric) {
					CHECK_THROWS_WITH(fixture.Run(), !ContainsSubstring(DenominatorColumnName(other)));
				}
			}
		}

		// cells_per_g_of_gdna has no denominator column, and its callers leave
		// the field at zero -- so an unconditional check here would make every
		// query for M3's metric throw.
		TinyCells gdna;
		gdna.params[0].sample_denominator = 0.0;
		CHECK_NOTHROW(gdna.Run());
		CHECK(gdna.Run().values.size() == 2);
	}

	SECTION("a cell count that overflows to a non-finite value") {
		// Every input here is individually finite and passes every screen; it is
		// their combination that leaves the representable range. pysyndna emits
		// the inf. Refusing is the same call the zero denominators get, and for
		// the same reason -- an infinite cell count is not a measurement.
		//
		// Reachable two ways, so both are pinned. First an extreme model, which
		// is NOT specific to the new metrics: cells_per_g_of_gdna overflows the
		// same way, and did so silently before this check existed.
		TinyCells huge_model;
		huge_model.models[0].intercept = 400.0;
		CHECK_THROWS_AS(huge_model.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(huge_model.Run(), ContainsSubstring("non-finite"));
		CHECK_THROWS_WITH(huge_model.Run(), ContainsSubstring("'s1 / f1'"));

		// Second the extraction product, which is new surface: the multiply that
		// makes a sample-level metric can overflow before the divide undoes it.
		TinyCells huge_extraction;
		huge_extraction.options.metric = CellCountsMetric::CellsPerUlOfSample;
		huge_extraction.params[0].extracted_gdna_concentration_ng_ul = 1e200;
		huge_extraction.params[0].vol_extracted_elution_ul = 1e200;
		CHECK_THROWS_WITH(huge_extraction.Run(), ContainsSubstring("non-finite"));

		// The identical extraction parameters are harmless to the metric that
		// never forms the ratio.
		TinyCells gdna = huge_extraction;
		gdna.options.metric = CellCountsMetric::CellsPerGOfGdna;
		CHECK_NOTHROW(gdna.Run());
		CHECK(gdna.Run().values.size() == 2);
	}

	SECTION("a predicted mass that underflows is omitted, not reported per cell") {
		// The other end of the same range. An extreme NEGATIVE intercept drives
		// 10^x to exactly zero, and under D10 a zero cell is an absent cell --
		// so the cell is omitted and nothing is said, because the sample is
		// still in the output and still means what it says.
		//
		// This is why the zero-valued list documents itself as "not reachable
		// through a realistic model" rather than "impossible": here is the
		// unrealistic one, and it behaves.
		TinyCells fixture;
		fixture.models[0].intercept = -400.0;
		const auto result = fixture.Run();
		CHECK(result.values.empty());
		// Every cell of the sample went that way, so the whole-sample case DOES
		// report -- for cells_per_g_of_gdna, which the header says is not
		// reachable in practice and does not claim is impossible.
		CHECK(result.zero_valued_sample_ids == std::vector<std::string> {"s1"});

		// One cell underflowing while others survive is silent, and correctly
		// so. Splitting the sample takes a wide spread of counts and an
		// intercept between them: at -330 with slope 1, f1's 10^-328 flushes to
		// zero while f2's 10^-300 is merely subnormal and comes out 9.3e-284.
		// (-320 does NOT split them: 10^-318 is subnormal, not zero.)
		TinyCells split;
		split.lengths[1].ogu_len_in_bp = split.lengths[0].ogu_len_in_bp;
		split.models[0].intercept = -330.0;
		split.counts[1].count = 1e30;
		const auto partial = split.Run();
		REQUIRE(partial.values.size() == 1);
		CHECK(partial.values[0].feature_id == "f2");
		CHECK(partial.zero_valued_sample_ids.empty());
	}

	SECTION("an rvalue outside [-1, 1]") {
		// A correlation coefficient cannot be 1.5, so this is a malformed models
		// relation rather than a weak model -- and squaring it would turn it into
		// an r^2 above 1 that passes any legal min_rsquared.
		for (const double bad : {-1.001, 1.001, 2.0}) {
			TinyCells fixture;
			fixture.models[0].rvalue = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("rvalue"));
		}
		// Both endpoints are what a perfect fit actually produces.
		for (const double good : {-1.0, 1.0}) {
			TinyCells fixture;
			fixture.models[0].rvalue = good;
			CHECK_NOTHROW(fixture.Run());
		}
	}
}

TEST_CASE("ComputeCellCounts treats an unusable model as no model", "[absquant]") {
	constexpr double kInf = std::numeric_limits<double>::infinity();
	// NOT an error, unlike an out-of-range rvalue. pysyndna records None for a
	// sample whose fit had any NaN field (_convert_linregressresults_to_dict
	// :329-334), so "absent" and "present but unusable" are the same thing to it
	// and reach the same log line. Keeping them one bucket here is what lets
	// absquant_fit_models' output be fed straight back in.
	for (size_t which = 0; which < 3; ++which) {
		for (const double bad : {std::nan(""), kInf, -kInf}) {
			TinyCells fixture;
			double *fields[] = {&fixture.models[0].slope, &fixture.models[0].intercept, &fixture.models[0].rvalue};
			INFO("model field " << which);
			*fields[which] = bad;

			const auto result = fixture.Run();
			CHECK(result.samples_without_models == std::vector<std::string> {"s1"});
			CHECK(result.values.empty());
			CHECK(result.low_rsquared_sample_ids.empty());
		}
	}
}

TEST_CASE("ComputeCellCounts enforces id consistency asymmetrically", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// Counts are the subject: every cell must be describable. Reference
	// relations covering more than the counts do is ordinary -- a lengths table
	// is a whole reference database, a parameters table a whole study -- so the
	// other direction is not checked at all.

	SECTION("a counted feature with no length") {
		TinyCells fixture;
		fixture.counts.push_back({"s1", "f9", 100.0});
		fixture.coverage.push_back({"s1", "f9", 0.5});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f9'"));
	}

	SECTION("a counted sample with no parameters") {
		TinyCells fixture;
		fixture.counts.push_back({"s9", "f1", 100.0});
		fixture.coverage.push_back({"s9", "f1", 0.5});
		fixture.models.push_back({"s9", 1.0, 0.0, 1.0});
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s9'"));
	}

	SECTION("a counted cell with no coverage") {
		// The upstream defect this closes: pysyndna validates only the two AXIS
		// id sets, so a sample and a feature can each be present while their
		// pair is not. Its left join then yields NaN, `NaN >= min_coverage` is
		// false, and the cell is dropped -- without even joining the
		// too-low-coverage list, because `NaN < min_coverage` is false too. The
		// cell vanishes with no log line of any kind.
		TinyCells fixture;
		fixture.counts.push_back({"s1", "f3", 100.0});
		fixture.lengths.push_back({"f3", 1000.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		// s1 already has cells and f3 now has a length, so both AXIS id sets are
		// satisfied; only the pair is missing. That is precisely the input
		// pysyndna's two set checks wave through, and why this check has to be
		// on the pair rather than on either axis.
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f3'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("coverage"));
	}

	SECTION("reference relations may describe more than the counts do") {
		TinyCells fixture;
		fixture.lengths.push_back({"f_unused", 1234.0});
		fixture.params.push_back({"s_unused", 10.0, 2.0, 100.0});
		fixture.models.push_back({"s_unused", 1.0, 0.0, 1.0});
		fixture.coverage.push_back({"s_unused", "f_unused", 0.5});
		CHECK_NOTHROW(fixture.Run());
		CHECK(fixture.Run().values.size() == 2);
	}

	SECTION("and their unused rows are not validated") {
		// Deliberately narrower than FitSyndnaModels, which screens every
		// concentration because it sums the whole relation into a denominator.
		// Nothing here is summed: a lengths table of ten thousand genomes should
		// not fail a query over twelve of them because of one bad row elsewhere.
		TinyCells fixture;
		fixture.lengths.push_back({"f_unused", -1.0});
		fixture.params.push_back({"s_unused", 0.0, 2.0, 100.0});
		fixture.models.push_back({"s_unused", 1.0, 0.0, 42.0});
		CHECK_NOTHROW(fixture.Run());

		// Including the denominator: a parameters relation is a whole study, and
		// a zero aliquot mass on a sample this query never counts is not this
		// query's problem.
		fixture.options.metric = CellCountsMetric::CellsPerGOfSample;
		fixture.params.back().sample_denominator = 0.0;
		CHECK_NOTHROW(fixture.Run());
		CHECK(fixture.Run().values.size() == 2);
	}
}

TEST_CASE("ComputeCellCounts checks the call before the data", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// Same three tiers as FitSyndnaModels: a bad call, then a malformed
	// relation, then relations that disagree. With several things wrong at once
	// the user should hear the most basic one first, because it is often the
	// cause of the others rather than a separate mistake.
	TinyCells worst;
	worst.options.min_coverage = 10.0;
	worst.counts.push_back({"s1", "f1", 100.0}); // duplicate cell
	worst.counts.push_back({"s1", "f9", 100.0}); // no length, no coverage
	CHECK_THROWS_WITH(worst.Run(), ContainsSubstring("min_coverage"));

	TinyCells malformed;
	malformed.counts.push_back({"s1", "f1", 100.0});
	malformed.counts.push_back({"s1", "f9", 100.0});
	CHECK_THROWS_WITH(malformed.Run(), ContainsSubstring("'s1 / f1'"));
	CHECK_THROWS_WITH(malformed.Run(), !ContainsSubstring("'f9'"));
}

namespace {

struct OrfFixture {
	std::vector<CountObservation> counts;
	std::vector<OrfCoords> coords;
	std::vector<SampleOrfParams> params;

	OrfCopiesResult Run() const {
		return ComputeOrfCopies(counts, coords, params);
	}
};

// Loads one ORF fixture family ("orf" = Set A, pysyndna's own; "orfb" = Set B,
// the synthetic filter case).
OrfFixture LoadOrfFixture(const std::string &prefix) {
	OrfFixture fixture;
	const std::string dir = "data/syndna/";

	const auto counts = ReadCsv(dir + prefix + "_counts.csv");
	const size_t c_sample = counts.Col("sample_id");
	const size_t c_feature = counts.Col("feature_id");
	const size_t c_value = counts.Col("value");
	for (const auto &row : counts.rows) {
		// The DuckDB reader drops zero-valued cells before the core ever sees
		// them (the sparse invariant), so the fixture loader must too -- these
		// CSVs are stored as dense as pysyndna's own input.
		const double value = AsDouble(row[c_value]);
		if (value == 0.0) {
			continue;
		}
		fixture.counts.push_back({row[c_sample], row[c_feature], value});
	}

	const auto coords = ReadCsv(dir + prefix + "_coords.csv");
	const size_t o_feature = coords.Col("feature_id");
	const size_t o_start = coords.Col("ogu_orf_start");
	const size_t o_end = coords.Col("ogu_orf_end");
	for (const auto &row : coords.rows) {
		fixture.coords.push_back({row[o_feature], AsDouble(row[o_start]), AsDouble(row[o_end])});
	}

	// An empty cell is pysyndna's NaN, which is what SQL NULL will arrive as;
	// AsDouble turns it into one. Set B's `rnan` sample is exactly that.
	const auto params = ReadCsv(dir + prefix + "_params.csv");
	const size_t p_sample = params.Col("sample_id");
	const size_t p_mass = params.Col("calc_mass_sample_aliquot_input_g");
	const size_t p_conc = params.Col("total_rna_concentration_ng_ul");
	const size_t p_vol = params.Col("vol_extracted_elution_ul");
	const size_t p_reads = params.Col("total_biological_reads_r1r2");
	for (const auto &row : params.rows) {
		fixture.params.push_back({row[p_sample], AsDouble(row[p_mass]), AsDouble(row[p_conc]), AsDouble(row[p_vol]),
		                          AsDouble(row[p_reads])});
	}
	return fixture;
}

// Compares ORF copy output against a committed oracle. Every row of an ORF
// oracle belongs to the one metric this function computes, so unlike the cells
// case there is nothing to select -- the whole file is the golden.
void CheckOrfAgainstOracle(const OrfCopiesResult &result, const std::string &oracle_path) {
	const auto gold = ReadCsv(oracle_path);
	const size_t g_sample = gold.Col("sample_id");
	const size_t g_feature = gold.Col("feature_id");
	const size_t g_value = gold.Col("value");

	std::vector<GoldenCell> cells;
	for (const auto &row : gold.rows) {
		cells.push_back({row[g_sample], row[g_feature], AsDouble(row[g_value])});
	}
	CheckValuesAgainstGolden(result.values, cells, oracle_path);
}

// The smallest input ComputeOrfCopies accepts, for pinning one factor of the
// chain at a time. Every number is round: f1 spans 100 bases and f2 spans 200,
// the sample has 1000 biological reads, 2 ng/uL over 50 uL is 1e-7 g of ssRNA,
// and half a gram of sample went into the extraction.
struct TinyOrf {
	std::vector<CountObservation> counts = {{"s1", "f1", 100.0}, {"s1", "f2", 400.0}};
	std::vector<OrfCoords> coords = {{"f1", 1.0, 100.0}, {"f2", 201.0, 400.0}};
	std::vector<SampleOrfParams> params = {{"s1", 0.5, 2.0, 50.0, 1000.0}};

	OrfCopiesResult Run() const {
		return ComputeOrfCopies(counts, coords, params);
	}
	// Value of one cell, or NaN if it was dropped.
	double ValueOf(const std::string &sample_id, const std::string &feature_id) const {
		for (const auto &cell : Run().values) {
			if (cell.sample_id == sample_id && cell.feature_id == feature_id) {
				return cell.value;
			}
		}
		return std::nan("");
	}
};

} // namespace

TEST_CASE("ComputeOrfCopies reproduces pysyndna's ORF copy goldens", "[absquant]") {
	// Set A: pysyndna's own RNA fixture, two metatranscriptomes against ten
	// OGU+ORFs. No options to set -- this workflow has no standard curve, no
	// coverage filter and no r^2 gate, so the oracle pins the entire chain.
	const auto result = LoadOrfFixture("orf").Run();
	CheckOrfAgainstOracle(result, "data/syndna/orf_oracle.csv");

	// 20 dense golden cells, of which 9 are zero because their read count is --
	// dropped by the reader, and matched above against an absent cell.
	CHECK(result.values.size() == 11);
	CHECK(result.filtered_sample_ids.empty());
}

TEST_CASE("ComputeOrfCopies reproduces the Set B goldens and filters `rnan`", "[absquant]") {
	// Set B exists for one reason: `rnan` has a NULL calc_mass_sample_aliquot_input_g
	// and is therefore absent from the oracle. pysyndna's REQUIRED_PARAM_KEYS
	// screen drops it (quant_orfs.py:313-317) rather than dividing by NaN, and
	// the golden is what pins that -- a port that kept the sample would emit
	// three NaN cells the parity check would catch in both directions.
	const auto fixture = LoadOrfFixture("orfb");
	const auto result = fixture.Run();
	CheckOrfAgainstOracle(result, "data/syndna/orfb_oracle.csv");

	// Two samples x three ORFs, less the two cells whose count is zero.
	CHECK(result.values.size() == 4);
	CHECK(result.filtered_sample_ids == std::vector<std::string> {"rnan"});
	for (const auto &cell : result.values) {
		CHECK(cell.sample_id != "rnan");
	}

	// ... and `rnan` really was in the input, so the assertion above is about
	// the filter rather than about a fixture that never had the sample.
	bool counted = false;
	for (const auto &row : fixture.counts) {
		counted = counted || row.sample_id == "rnan";
	}
	CHECK(counted);
}

TEST_CASE("ComputeOrfCopies measures an ORF inclusively and without direction", "[absquant]") {
	// Both halves of `|end - start| + 1`, which pysyndna spells out in three
	// statements (quant_orfs.py:50-57) and which four of Set A's ten features
	// depend on.
	const TinyOrf base;
	const double v = base.ValueOf("s1", "f1");
	REQUIRE(std::isfinite(v));
	REQUIRE(v > 0.0);

	// Reverse strand: start > end is legal and spans the same bases, so it must
	// give the same answer bit for bit rather than a negative length, which
	// would flip the sign of every copy count.
	TinyOrf reversed;
	reversed.coords[0] = {"f1", 100.0, 1.0};
	CHECK(reversed.ValueOf("s1", "f1") == v);

	// Inclusive: a single-base ORF has length 1, not 0. Without the +1 this
	// divides by zero and returns inf, which is the failure the constant guards
	// against; with a wrong +2 it would merely be quietly off.
	TinyOrf single;
	single.coords[0] = {"f1", 42.0, 42.0};
	const double one_base = single.ValueOf("s1", "f1");
	REQUIRE(std::isfinite(one_base));
	CHECK(CloseEnough(one_base, v * 100.0, 1e-14, 0.0));

	// Length is a denominator: an endpoint at 200 rather than 100 doubles the
	// span to 200 bases and halves the value. Written as an endpoint rather than
	// as a length so an off-by-one in the span is visible here.
	TinyOrf longer;
	longer.coords[0].ogu_orf_end = 200.0;
	CHECK(CloseEnough(longer.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));
}

TEST_CASE("ComputeOrfCopies matches pysyndna's published copies-per-gram table", "[absquant]") {
	// The one factor that is pure arithmetic on the coordinates, pinned against
	// pysyndna's OWN fixture table (tests/test_quant_orfs.py:39-43) rather than
	// against this port's oracle. That makes it an independent check on both the
	// 340 g/mole constant and the inclusive length -- the oracle would agree
	// with a port that had the same two mistakes it does.
	//
	// The rest of the chain is arranged to be exactly 1: the read count equals
	// the sample's total biological reads, 1e7 ng/uL over 100 uL is exactly one
	// gram of ssRNA, and one gram of sample went in. What is left IS
	// copies_per_g_ogu_orf_ss_rna.
	struct Published {
		const char *feature_id;
		double start;
		double end;
		double copies_per_g;
	};
	const Published table[] = {
	    {"G000005825_1", 816, 2168, 1.3091041e18},           {"G000005825_2", 2348, 3490, 1.5496219e18},
	    {"G000005825_3", 3744, 3959, 8.2000827e18},          {"G000005825_4", 3971, 5086, 1.5871128e18},
	    {"G000005825_5", 5098, 5373, 6.4174561e18},          {"G900163845_3247", 3392209, 3390413, 9.8565268e17},
	    {"G900163845_3248", 3393051, 3392206, 2.0936381e18}, {"G900163845_3249", 3393938, 3393048, 1.9878988e18},
	    {"G900163845_3250", 3394702, 3393935, 2.3062733e18}, {"G900163845_3251", 3395077, 3395721, 2.7460742e18}};

	for (const auto &row : table) {
		OrfFixture fixture;
		fixture.counts = {{"s1", row.feature_id, 1000.0}};
		fixture.coords = {{row.feature_id, row.start, row.end}};
		fixture.params = {{"s1", 1.0, 1e7, 100.0, 1000.0}};

		const auto result = fixture.Run();
		INFO("feature " << row.feature_id);
		REQUIRE(result.values.size() == 1);
		// The table is printed to eight significant figures, so that is the
		// most this can assert; the committed oracle carries the full precision.
		CHECK(CloseEnough(result.values[0].value, row.copies_per_g, 1e-7, 0.0));
	}
}

TEST_CASE("ComputeOrfCopies scales as the method requires", "[absquant]") {
	// The oracles pin the whole chain at once; these pin each factor's role
	// separately, so an inverted division or a swapped parameter column is
	// localized rather than showing up as "all 11 values differ".
	const TinyOrf base;
	const double v = base.ValueOf("s1", "f1");
	REQUIRE(std::isfinite(v));
	REQUIRE(v > 0.0);

	// The absolute value first, so that everything below is anchored rather than
	// merely self-consistent: a chain that is wrong by a constant factor
	// satisfies every ratio in this test. Worked out from the fixture's own
	// numbers, away from this implementation --
	//
	//   fraction     = 100 / 1000                          = 0.1
	//   g_total      = 2 * 50 / 1e9                         = 1e-7
	//   copies_per_g = 6.02214076e23 / (100 * 340)          = 1.7712178705882353e19
	//   value        = 0.1 * 1e-7 * 1.7712178705882353e19 / 0.5
	//
	// -- and f2 is 4x the reads over 2x the length, so exactly twice f1.
	CHECK(CloseEnough(v, 354243574117.64703, 1e-12, 0.0));
	CHECK(CloseEnough(base.ValueOf("s1", "f2"), 708487148235.2941, 1e-12, 0.0));

	// The read count is a numerator: twice the reads on this ORF, twice the
	// copies of it.
	TinyOrf counted;
	counted.counts[0].count *= 2.0;
	CHECK(CloseEnough(counted.ValueOf("s1", "f1"), v * 2.0, 1e-14, 0.0));

	// ... but only because it is a FRACTION of the sample's reads: twice the
	// biological reads with the same count on this ORF is half the share, so
	// half the copies. This is the factor a port is likeliest to drop entirely,
	// since omitting it changes nothing about the shape of the answer.
	TinyOrf deeper;
	deeper.params[0].total_biological_reads_r1r2 *= 2.0;
	CHECK(CloseEnough(deeper.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));

	// The extracted ssRNA mass is a numerator, and reaches the value through
	// two separate columns whose product it is -- so both are checked, and a
	// port that read the volume where it meant the concentration would still
	// pass one of them.
	TinyOrf concentrated;
	concentrated.params[0].total_rna_concentration_ng_ul *= 2.0;
	CHECK(CloseEnough(concentrated.ValueOf("s1", "f1"), v * 2.0, 1e-14, 0.0));

	TinyOrf eluted;
	eluted.params[0].vol_extracted_elution_ul *= 2.0;
	CHECK(CloseEnough(eluted.ValueOf("s1", "f1"), v * 2.0, 1e-14, 0.0));

	// The sample aliquot is the final denominator, and the one that makes the
	// answer per-gram-of-sample rather than per-extraction: twice the sample in,
	// half the copies per gram.
	TinyOrf heavier;
	heavier.params[0].calc_mass_sample_aliquot_input_g *= 2.0;
	CHECK(CloseEnough(heavier.ValueOf("s1", "f1"), v / 2.0, 1e-14, 0.0));

	// Every cell is per-sample: a second sample with its own parameters must not
	// disturb the first.
	TinyOrf two_samples;
	two_samples.counts.push_back({"s2", "f1", 100.0});
	two_samples.params.push_back({"s2", 1.0, 2.0, 50.0, 1000.0});
	CHECK(CloseEnough(two_samples.ValueOf("s1", "f1"), v, 1e-14, 0.0));
	CHECK(CloseEnough(two_samples.ValueOf("s2", "f1"), v / 2.0, 1e-14, 0.0));
}

TEST_CASE("ComputeOrfCopies screens all four parameter columns", "[absquant]") {
	// pysyndna filters on REQUIRED_PARAM_KEYS, which is the union of its
	// sample-info and RNA-prep key lists (quant_orfs.py:21-22) -- all four
	// columns, not just the ones that happen to be denominators. Walking them by
	// pointer-to-member is what makes a fifth column, or a forgotten fourth,
	// fail here rather than at a user's query.
	double SampleOrfParams::*const columns[] = {
	    &SampleOrfParams::calc_mass_sample_aliquot_input_g, &SampleOrfParams::total_rna_concentration_ng_ul,
	    &SampleOrfParams::vol_extracted_elution_ul, &SampleOrfParams::total_biological_reads_r1r2};
	REQUIRE(std::end(columns) - std::begin(columns) == 4);

	// NaN is SQL NULL, and negative is a value that cannot be what it claims to
	// be. Infinity is miint's own addition (see IsUsableSampleParameter), and
	// what pysyndna does with it here depends on WHICH column carries it: an
	// infinite denominator (the biological reads, the aliquot mass) collapses the
	// sample to zeros, while an infinite numerator (the concentration, the
	// elution volume) writes +inf into every cell of it. Two different wrong
	// answers from one unusable input, and reporting it as a bad parameter is the
	// difference between a diagnosis and a mystery.
	for (const double bad : {std::nan(""), -1.0, std::numeric_limits<double>::infinity()}) {
		for (const auto column : columns) {
			TinyOrf fixture;
			fixture.params[0].*column = bad;
			const auto result = fixture.Run();
			INFO("bad value " << bad);
			CHECK(result.filtered_sample_ids == std::vector<std::string> {"s1"});
			CHECK(result.values.empty());
		}
	}

	// Zero passes this filter, exactly as pysyndna's strict `< 0` does. Checked
	// on the two columns for which zero stays a legal input -- a blank
	// extraction really can yield no ssRNA -- because for the two denominators
	// zero is refused outright, which is the validation pass's business rather
	// than the filter's.
	for (const auto column :
	     {&SampleOrfParams::total_rna_concentration_ng_ul, &SampleOrfParams::vol_extracted_elution_ul}) {
		TinyOrf fixture;
		fixture.params[0].*column = 0.0;
		CHECK(fixture.Run().filtered_sample_ids.empty());
	}
}

TEST_CASE("ComputeOrfCopies rejects malformed relations", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	constexpr double kInf = std::numeric_limits<double>::infinity();

	SECTION("a duplicated counts cell") {
		// Every relation here is keyed, and a repeat means the caller's join
		// fanned out. pysyndna cannot express one -- its input is a biom table,
		// unique by construction -- so this is a check the long-form shape needs
		// and the reference implementation has no equivalent of.
		TinyOrf fixture;
		fixture.counts.push_back({"s1", "f1", 100.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
	}

	SECTION("a duplicated coordinate row") {
		TinyOrf fixture;
		fixture.coords.push_back({"f1", 1.0, 100.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f1'"));
	}

	SECTION("a duplicated parameter row") {
		TinyOrf fixture;
		fixture.params.push_back({"s1", 0.5, 2.0, 50.0, 1000.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
	}

	SECTION("a read count that is negative or not finite") {
		for (const double bad : {-1.0, kInf, -kInf, std::nan("")}) {
			TinyOrf fixture;
			fixture.counts[0].count = bad;
			CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
			// The MESSAGE matters here, not just that something threw. An
			// infinite or NaN count sails through to the arithmetic and comes
			// back as a non-finite value, so the overflow guard would refuse it
			// too -- with a message sending the user to check the ORF
			// coordinates and the extraction parameters, neither of which is
			// wrong. Asserting only that it throws lets this check be deleted
			// with nothing but the diagnosis changing, which is exactly how a
			// mutation survived here.
			CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("read counts"));
		}
	}

	SECTION("a zero read count is legal, unlike in the cell-count chain") {
		// The difference is the log10. ComputeCellCounts refuses a zero count
		// because log10(0) under a negative slope returns an INFINITE cell count;
		// nothing here takes a logarithm, so a zero count simply means no copies
		// of that ORF. Refusing it would be a rule invented for no reason.
		TinyOrf fixture;
		fixture.counts[0].count = 0.0;
		CHECK_NOTHROW(fixture.Run());

		// ... and the cell is omitted rather than written as a zero, which is the
		// sparse invariant (D10) doing exactly what it does for a zero VALUE.
		const auto result = fixture.Run();
		REQUIRE(result.values.size() == 1);
		CHECK(result.values[0].feature_id == "f2");
		CHECK(result.zero_valued_sample_ids.empty());
	}

	SECTION("a coordinate that is not finite") {
		for (const double bad : {kInf, -kInf, std::nan("")}) {
			TinyOrf start;
			start.coords[0].ogu_orf_start = bad;
			CHECK_THROWS_AS(start.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(start.Run(), ContainsSubstring("'f1'"));

			TinyOrf end;
			end.coords[0].ogu_orf_end = bad;
			CHECK_THROWS_AS(end.Run(), std::invalid_argument);
			CHECK_THROWS_WITH(end.Run(), ContainsSubstring("'f1'"));
		}
	}

	SECTION("a coordinate that is not a whole number") {
		// pysyndna's cast_cols would take 100.5 and hand back a fractional ORF
		// length, which then divides into Avogadro's number and produces a
		// perfectly plausible-looking copy count. A genome coordinate that is not
		// a whole base is not a measurement -- the same call D9 makes on
		// percent-versus-fraction coverage.
		TinyOrf fixture;
		fixture.coords[0].ogu_orf_end = 100.5;
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f1'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("ogu_orf"));

		// Negative coordinates are NOT refused: |end - start| + 1 is well defined
		// for any pair of whole numbers, and woltka's format says nothing about
		// the origin. Inventing a rule the reference does not have would reject
		// legitimate input for the sake of tidiness.
		TinyOrf negative;
		negative.coords[0] = {"f1", -100.0, -1.0};
		CHECK_NOTHROW(negative.Run());
		CHECK(negative.ValueOf("s1", "f1") == TinyOrf().ValueOf("s1", "f1"));
	}

	SECTION("a zero total_biological_reads_r1r2, which the filter lets past") {
		// The denominator D23 is about. Zero is neither NaN nor negative, so it
		// passes pysyndna's screen and then divides through to an INFINITE copy
		// count for every ORF in the sample, silently. Refusing it is the same
		// call ComputeCellCounts already makes for its two denominators.
		TinyOrf fixture;
		fixture.params[0].total_biological_reads_r1r2 = 0.0;
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("total_biological_reads_r1r2"));
	}

	SECTION("a zero calc_mass_sample_aliquot_input_g, likewise") {
		TinyOrf fixture;
		fixture.params[0].calc_mass_sample_aliquot_input_g = 0.0;
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1'"));
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("calc_mass_sample_aliquot_input_g"));
	}

	SECTION("the two zero denominators are refused, the other two columns are not") {
		// The split this function draws, stated in one place: zero is an ERROR
		// where it would divide, and ordinary data where it would not. A blank
		// extraction genuinely yields no ssRNA, and saying so is a measurement.
		TinyOrf concentration;
		concentration.params[0].total_rna_concentration_ng_ul = 0.0;
		CHECK_NOTHROW(concentration.Run());

		TinyOrf volume;
		volume.params[0].vol_extracted_elution_ul = 0.0;
		CHECK_NOTHROW(volume.Run());
	}
}

TEST_CASE("ComputeOrfCopies screens the coordinates only where the counts reach them", "[absquant]") {
	// A coords relation is a whole annotation -- typically every ORF of every
	// genome in a reference database -- while the counts name the handful that
	// sequenced. Failing a query over two ORFs because a third nobody asked
	// about has a malformed coordinate would be its own kind of wrong answer.
	// Same rule ComputeCellCounts applies to its lengths and parameters.
	TinyOrf fixture;
	fixture.coords.push_back({"unused", std::nan(""), 100.5});
	CHECK_NOTHROW(fixture.Run());
	CHECK(fixture.Run().values.size() == 2);

	// The same coordinate becomes an error the moment something counts it.
	TinyOrf counted;
	counted.coords.push_back({"unused", std::nan(""), 100.5});
	counted.counts.push_back({"s1", "unused", 10.0});
	CHECK_THROWS_AS(counted.Run(), std::invalid_argument);
}

TEST_CASE("ComputeOrfCopies enforces id consistency asymmetrically", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// The counts relation is the subject: every cell it names must have
	// coordinates and its sample must have parameters, or nothing defined can be
	// said about that cell. The reverse never has to hold. pysyndna splits it
	// the same way at quant_orfs.py:291, whose own comment says extra parameter
	// samples "could just be samples that failed sequencing/etc.".
	SECTION("a counted feature with no coordinates") {
		// pysyndna does not validate this direction at all -- it reaches the
		// coordinates through a bare `.at[]` inside a biom transform, so a
		// missing ORF surfaces as a pandas KeyError from inside a lambda.
		TinyOrf fixture;
		fixture.counts.push_back({"s1", "f9", 10.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'f9'"));
	}

	SECTION("a counted sample with no parameters") {
		TinyOrf fixture;
		fixture.counts.push_back({"s9", "f1", 10.0});
		CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
		CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s9'"));
	}

	SECTION("coordinates for ORFs nobody counted") {
		TinyOrf fixture;
		fixture.coords.push_back({"f9", 1.0, 500.0});
		CHECK_NOTHROW(fixture.Run());
		CHECK(fixture.Run().values.size() == 2);
	}

	SECTION("parameters for samples that were never sequenced") {
		TinyOrf fixture;
		fixture.params.push_back({"s9", 0.5, 2.0, 50.0, 1000.0});
		CHECK_NOTHROW(fixture.Run());
		CHECK(fixture.Run().values.size() == 2);
		CHECK(fixture.Run().filtered_sample_ids.empty());
		CHECK(fixture.Run().zero_valued_sample_ids.empty());
	}
}

TEST_CASE("ComputeOrfCopies reports a malformed relation before an id mismatch", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// With several things wrong at once the user should hear the most basic one
	// first, because it is often the cause of the others rather than a separate
	// mistake -- a fanned-out join produces both symptoms at the same time.
	TinyOrf fixture;
	fixture.counts.push_back({"s1", "f1", 100.0}); // duplicate cell
	fixture.counts.push_back({"s1", "f9", 100.0}); // no coordinates
	CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));
	CHECK_THROWS_WITH(fixture.Run(), !ContainsSubstring("'f9'"));
}

TEST_CASE("ComputeOrfCopies refuses a copy count that overflowed", "[absquant]") {
	using Catch::Matchers::ContainsSubstring;
	// Every input finite, their combination not: the extracted ssRNA mass is a
	// product formed before the ng->g divide, so two large-but-legal values
	// overflow it to +inf and carry that through the whole chain. pysyndna emits
	// the inf. An infinite copy count is not a measurement, and a DOUBLE column
	// full of inf is worse than a failed query.
	TinyOrf fixture;
	fixture.params[0].total_rna_concentration_ng_ul = 1e300;
	fixture.params[0].vol_extracted_elution_ul = 1e300;
	CHECK_THROWS_AS(fixture.Run(), std::invalid_argument);
	CHECK_THROWS_WITH(fixture.Run(), ContainsSubstring("'s1 / f1'"));

	// The guard is not simply refusing anything large. At 1e140 the same two
	// columns give a copy count around 1e289 -- absurd as biology, entirely
	// representable as a double -- and it is computed rather than refused.
	//
	// Note where the ceiling actually sits. At 1e150 the PRODUCT is finite
	// (1e300) and it is the later multiply by copies-per-gram that overflows, so
	// guarding the two columns, or their product, would still have let an inf
	// through. Guarding the computed value is what closes it.
	TinyOrf large;
	large.params[0].total_rna_concentration_ng_ul = 1e140;
	large.params[0].vol_extracted_elution_ul = 1e140;
	CHECK_NOTHROW(large.Run());
	CHECK(std::isfinite(large.ValueOf("s1", "f1")));
	CHECK(large.ValueOf("s1", "f1") > 1e280);

	TinyOrf midway;
	midway.params[0].total_rna_concentration_ng_ul = 1e150;
	midway.params[0].vol_extracted_elution_ul = 1e150;
	CHECK(std::isfinite(midway.params[0].total_rna_concentration_ng_ul * midway.params[0].vol_extracted_elution_ul /
	                    1e9));
	CHECK_THROWS_AS(midway.Run(), std::invalid_argument);
}

TEST_CASE("ComputeOrfCopies reports a sample whose every cell came out zero", "[absquant]") {
	// The sparse form cannot distinguish "every cell of this sample is zero"
	// from "this sample was never here", so the all-zero case is named. A zero
	// extraction reaches it legitimately: both columns pass pysyndna's `< 0`
	// screen, and a blank really can yield no ssRNA.
	for (const auto column :
	     {&SampleOrfParams::total_rna_concentration_ng_ul, &SampleOrfParams::vol_extracted_elution_ul}) {
		TinyOrf fixture;
		fixture.params[0].*column = 0.0;
		const auto result = fixture.Run();
		CHECK(result.values.empty());
		CHECK(result.zero_valued_sample_ids == std::vector<std::string> {"s1"});
		// Reported here and NOT as a bad parameter: zero is a value, not a gap,
		// and the two lists answer different questions.
		CHECK(result.filtered_sample_ids.empty());
	}

	// Reachable by underflow too, which is the honest boundary: 1e-300 squared
	// is below the smallest subnormal, so the product is exactly 0.0 before the
	// ng->g divide ever runs.
	TinyOrf underflowed;
	underflowed.params[0].total_rna_concentration_ng_ul = 1e-300;
	underflowed.params[0].vol_extracted_elution_ul = 1e-300;
	const auto result = underflowed.Run();
	CHECK(result.values.empty());
	CHECK(result.zero_valued_sample_ids == std::vector<std::string> {"s1"});

	// A single zero cell among others needs no diagnostic: under D10 an omitted
	// cell and a dense 0.0 are the same claim, and the sample is still in the
	// output saying what it is. Only the WHOLE sample going that way is
	// ambiguous.
	TinyOrf one_zero;
	one_zero.counts[0].count = 0.0;
	CHECK(one_zero.Run().values.size() == 1);
	CHECK(one_zero.Run().zero_valued_sample_ids.empty());

	// A filtered sample is not also reported as zero-valued -- every sample
	// lands in exactly one place.
	TinyOrf filtered;
	filtered.params[0].total_rna_concentration_ng_ul = std::nan("");
	CHECK(filtered.Run().filtered_sample_ids == std::vector<std::string> {"s1"});
	CHECK(filtered.Run().zero_valued_sample_ids.empty());
}
