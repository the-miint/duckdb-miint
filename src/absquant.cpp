#include "absquant.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>
#include <utility>

namespace miint {
namespace absquant {

namespace {

constexpr double kNaN = std::numeric_limits<double>::quiet_NaN();

// Continued-fraction expansion of the incomplete beta (DLMF 8.17.22), evaluated
// by the modified Lentz algorithm. Valid where it converges quickly, i.e. for
// x < (a+1)/(a+b+2); the caller reflects otherwise.
//
// Lentz's method carries the fraction as a running product of correction
// factors, guarding each denominator against zero with kFpMin rather than
// building a numerator/denominator pair that can overflow. kEps sets the
// convergence target; at 3e-16 the result sits at the limit of double
// precision, far tighter than the ~6 significant figures a pvalue needs.
double BetaContinuedFraction(double a, double b, double x) {
	constexpr int kMaxIter = 300;
	constexpr double kEps = 3.0e-16;
	constexpr double kFpMin = 1.0e-300;

	const double qab = a + b;
	const double qap = a + 1.0;
	const double qam = a - 1.0;

	double c = 1.0;
	double d = 1.0 - qab * x / qap;
	if (std::fabs(d) < kFpMin) {
		d = kFpMin;
	}
	d = 1.0 / d;
	double h = d;

	for (int m = 1; m <= kMaxIter; ++m) {
		const double dm = static_cast<double>(m);
		const double m2 = 2.0 * dm;

		// Even step.
		double aa = dm * (b - dm) * x / ((qam + m2) * (a + m2));
		d = 1.0 + aa * d;
		if (std::fabs(d) < kFpMin) {
			d = kFpMin;
		}
		c = 1.0 + aa / c;
		if (std::fabs(c) < kFpMin) {
			c = kFpMin;
		}
		d = 1.0 / d;
		h *= d * c;

		// Odd step.
		aa = -(a + dm) * (qab + dm) * x / ((a + m2) * (qap + m2));
		d = 1.0 + aa * d;
		if (std::fabs(d) < kFpMin) {
			d = kFpMin;
		}
		c = 1.0 + aa / c;
		if (std::fabs(c) < kFpMin) {
			c = kFpMin;
		}
		d = 1.0 / d;
		const double del = d * c;
		h *= del;

		if (std::fabs(del - 1.0) < kEps) {
			return h;
		}
	}
	// Non-convergence is unreachable for the (a, b, x) this port generates --
	// a = df/2, b = 1/2, with the reflection applied -- but returning NaN rather
	// than a half-converged value keeps a silent wrong answer impossible should
	// a future caller widen the domain.
	return kNaN;
}

} // namespace

double RegularizedIncompleteBetaC(double a, double b, double x, double one_minus_x) {
	if (std::isnan(a) || std::isnan(b) || std::isnan(x) || std::isnan(one_minus_x)) {
		return kNaN;
	}
	if (a <= 0.0 || b <= 0.0 || x < 0.0 || x > 1.0 || one_minus_x < 0.0 || one_minus_x > 1.0) {
		return kNaN;
	}
	// Exact at the boundaries; the continued fraction is neither needed nor
	// well-behaved there. Keyed on the complement too, so a genuinely-zero
	// complement is treated as the boundary even when x rounded to 1.0.
	if (x == 0.0) {
		return 0.0;
	}
	if (one_minus_x == 0.0) {
		return 1.0;
	}

	// Front factor x^a (1-x)^b / B(a,b), in log space so the large-df cases
	// (a = df/2, up to 500 on the committed grid) cannot overflow. The two log
	// terms take x and its complement separately -- that is the whole point of
	// this entry point, since log(x) is harmless when x has rounded to 1.0 but
	// log(1-x) computed from that same rounded x would be log(0).
	const double log_front =
	    std::lgamma(a + b) - std::lgamma(a) - std::lgamma(b) + a * std::log(x) + b * std::log(one_minus_x);
	const double front = std::exp(log_front);

	// Reflect when x is past the continued fraction's fast-convergence region,
	// using I_x(a,b) = 1 - I_{1-x}(b,a). The subtraction is safe here precisely
	// because this branch is only taken when the result is near 1: the term
	// being subtracted is at most O(1), never a value whose cancellation would
	// wipe out the answer. Getting that backwards is exactly how an "obvious"
	// small-x-first branch destroys the far tail -- I_y(0.5, 500) at y = 0.0909
	// is 1 - 1.7e-22, so forming 1 - I_y there yields 0 instead of 8.3e-23.
	if (x < (a + 1.0) / (a + b + 2.0)) {
		return front * BetaContinuedFraction(a, b, x) / a;
	}
	return 1.0 - front * BetaContinuedFraction(b, a, one_minus_x) / b;
}

double RegularizedIncompleteBeta(double a, double b, double x) {
	// Safe whenever the caller has no better handle on 1-x than this; callers
	// that do (StudentTSurvival) must use the C form.
	return RegularizedIncompleteBetaC(a, b, x, 1.0 - x);
}

double StudentTSurvival(double t, double df) {
	if (std::isnan(t) || std::isnan(df) || df <= 0.0) {
		return kNaN;
	}
	// Exact, and the value a caller is most likely to eyeball. Going through the
	// incomplete beta would give 0.5 only to within rounding.
	if (t == 0.0) {
		return 0.5;
	}
	if (std::isinf(t)) {
		return t > 0.0 ? 0.0 : 1.0;
	}

	// sf(t, df) = 0.5 * I_{df/(df+t^2)}(df/2, 1/2) for t >= 0, mirrored below.
	// Reproduces scipy.stats.t.sf on data/syndna/studentt_sf_oracle.csv.
	//
	// w = t^2/df is formed as (|t|/df)*|t| rather than t*t/df so the enormous t
	// the TINY guard produces at a perfect fit (up to ~2.2e11 on that grid) keeps
	// its ratio: t*t alone reaches ~5e22, which is representable, but dividing
	// first avoids ever forming it.
	const double abs_t = std::fabs(t);
	const double w = (abs_t / df) * abs_t;

	// x and its complement are both derived from w, so both stay accurate in the
	// regime where they matter: x -> 1 as t -> 0 and rounds away, but the
	// complement is then ~w and exact; x is tiny for large t and exact there.
	// Handing both to the incomplete beta means neither regime has to reconstruct
	// the other by subtracting from 1.
	const double x = 1.0 / (1.0 + w);
	const double one_minus_x = w / (1.0 + w);
	const double half = 0.5 * RegularizedIncompleteBetaC(0.5 * df, 0.5, x, one_minus_x);
	return t > 0.0 ? half : 1.0 - half;
}

namespace {

// scipy.stats.linregress's TINY, added inside the t statistic so that a perfect
// fit produces a huge but finite t instead of dividing by zero. Not a p-value
// floor: at df = 1 it yields pvalue 9.0e-11, a genuine number.
constexpr double kTiny = 1.0e-20;

} // namespace

// pysyndna removes a sample whose required per-sample parameter is missing or
// negative (util.py:258-266, isna().any() | lt(0).any()). Zero passes -- the
// test is strict.
bool IsUsableSampleParameter(double value) {
	// isfinite, not !isnan: +inf passes every naive test (it is not NaN and it
	// is >= 0) and then poisons the arithmetic downstream with no diagnostic --
	// an infinite pool mass makes every log10 in the sample +inf, whose
	// differences are NaN, so the sample is reported "unfittable" as though the
	// fit had been attempted and failed.
	return std::isfinite(value) && value >= 0.0;
}

// The metric table: one switch per mapping the rest of the port needs, the SQL
// spelling and the parameters column. Both are exhaustive with no `default:`,
// but this build does not turn on -Wswitch, so a fifth enumerator would fall
// through to the CellsPerGOfGdna answer rather than fail to compile. What
// actually keeps a metric from being half-added is that every caller iterates
// kAllCellCountsMetrics instead of restating the set.
const char *MetricName(CellCountsMetric metric) {
	switch (metric) {
	case CellCountsMetric::CellsPerGOfSample:
		return "cells_per_g_of_sample";
	case CellCountsMetric::CellsPerUlOfSample:
		return "cells_per_ul_of_sample";
	case CellCountsMetric::CellsPerCm2OfSample:
		return "cells_per_cm2_of_sample";
	case CellCountsMetric::CellsPerGOfGdna:
		break;
	}
	return "cells_per_g_of_gdna";
}

const char *DenominatorColumnName(CellCountsMetric metric) {
	switch (metric) {
	case CellCountsMetric::CellsPerGOfSample:
		return "calc_mass_sample_aliquot_input_g";
	case CellCountsMetric::CellsPerUlOfSample:
		return "sample_volume_ul";
	case CellCountsMetric::CellsPerCm2OfSample:
		return "sample_surface_area_cm2";
	case CellCountsMetric::CellsPerGOfGdna:
		break;
	}
	// Not "no such column": this metric normalizes by a parameter it already
	// has, so requiring one of the three sample columns to exist would refuse
	// relations that can answer the question perfectly well.
	return nullptr;
}

bool ParseCellCountsMetric(const std::string &name, CellCountsMetric &out) {
	// Folded here rather than assumed of the caller. The DuckDB wrapper does
	// lowercase its input, but a precondition nothing enforces is one refactor
	// away from being a silently rejected metric.
	std::string lowered(name);
	std::transform(lowered.begin(), lowered.end(), lowered.begin(),
	               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	for (const auto metric : kAllCellCountsMetrics) {
		if (lowered == MetricName(metric)) {
			out = metric;
			return true;
		}
	}
	return false;
}

std::vector<std::string> IdsMissingFrom(const std::vector<std::string> &subject,
                                        const std::vector<std::string> &reference) {
	const std::unordered_set<std::string> known(reference.begin(), reference.end());
	std::vector<std::string> missing;
	for (const auto &id : subject) {
		if (known.find(id) == known.end()) {
			missing.push_back(id);
		}
	}
	std::sort(missing.begin(), missing.end());
	missing.erase(std::unique(missing.begin(), missing.end()), missing.end());
	return missing;
}

std::vector<std::string> DuplicatedIds(const std::vector<std::string> &ids) {
	std::vector<std::string> sorted(ids);
	std::sort(sorted.begin(), sorted.end());
	std::vector<std::string> repeated;
	for (size_t i = 1; i < sorted.size(); ++i) {
		if (sorted[i] == sorted[i - 1] && (repeated.empty() || repeated.back() != sorted[i])) {
			repeated.push_back(sorted[i]);
		}
	}
	return repeated;
}

std::string FormatIdList(std::vector<std::string> ids) {
	std::sort(ids.begin(), ids.end());
	// Deduplicate here rather than trusting callers to. Most arrive from
	// DuplicatedIds or IdsMissingFrom, which already dedupe, but a caller that
	// collects offenders row by row will hand over one entry per offending ROW
	// -- and then a cell repeated three times renders as "'x', 'x'" and eats two
	// of the ten slots. Cheap, and it removes a contract nobody would remember.
	ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
	constexpr size_t kCap = 10;
	std::string out;
	for (size_t i = 0; i < ids.size() && i < kCap; ++i) {
		if (i > 0) {
			out += ", ";
		}
		out += "'" + ids[i] + "'";
	}
	if (ids.size() > kCap) {
		out += ", ... (" + std::to_string(ids.size() - kCap) + " more)";
	}
	return out;
}

namespace {

// Reject every input FitSyndnaModels cannot compute a defined answer from,
// before any arithmetic runs. Ordered from "the call itself is wrong" through
// "a relation is malformed" to "the relations disagree", so that the first
// error a user sees is the most basic thing wrong with their input rather than
// whichever check happened to be written first.
//
// Messages carry no function-name prefix: the DuckDB wrapper adds one when it
// re-throws as InvalidInputException, the same contract BuildDenseDistanceMatrix
// has with ReadDistanceTable.
void ValidateFitInputs(const std::vector<CountObservation> &counts,
                       const std::vector<SyndnaConcentration> &concentrations, const std::vector<SampleMass> &masses,
                       const FitOptions &options) {
	// Written as !(x > 0) rather than x <= 0 so a NaN option is rejected too.
	if (!(options.syndna_contributing_fraction > 0.0) || !(options.syndna_contributing_fraction <= 1.0)) {
		throw std::invalid_argument("syndna_contributing_fraction must be > 0 and <= 1 (got " +
		                            std::to_string(options.syndna_contributing_fraction) + ")");
	}
	// pysyndna requires >= 1 (fit_syndna_models.py:468-471). It is what makes
	// the concentrations-only synDNA rule below safe: a synDNA with zero total
	// reads is always dropped, so miint never has to distinguish it from an
	// absent one.
	if (options.min_syndna_counts < 1) {
		throw std::invalid_argument("min_syndna_counts must be >= 1 (got " + std::to_string(options.min_syndna_counts) +
		                            ")");
	}

	std::vector<std::string> concentration_ids;
	concentration_ids.reserve(concentrations.size());
	for (const auto &row : concentrations) {
		concentration_ids.push_back(row.feature_id);
	}
	const auto repeated_concentrations = DuplicatedIds(concentration_ids);
	if (!repeated_concentrations.empty()) {
		throw std::invalid_argument("the synDNA concentrations relation has more than one row for feature_id " +
		                            FormatIdList(repeated_concentrations));
	}

	std::vector<std::string> mass_ids;
	mass_ids.reserve(masses.size());
	for (const auto &row : masses) {
		mass_ids.push_back(row.sample_id);
	}
	const auto repeated_masses = DuplicatedIds(mass_ids);
	if (!repeated_masses.empty()) {
		throw std::invalid_argument("the sample parameters relation has more than one row for sample_id " +
		                            FormatIdList(repeated_masses));
	}

	// (sample_id, feature_id) must be unique. Compared as index-free string
	// pairs rather than a packed key so no separator character can collide with
	// an id that legitimately contains it.
	std::set<std::pair<std::string, std::string>> seen;
	std::vector<std::string> repeated_cells;
	for (const auto &row : counts) {
		if (!seen.emplace(row.sample_id, row.feature_id).second) {
			repeated_cells.push_back(row.sample_id + " / " + row.feature_id);
		}
	}
	if (!repeated_cells.empty()) {
		throw std::invalid_argument(
		    "the synDNA counts relation has more than one row for the same (sample_id, feature_id): " +
		    FormatIdList(repeated_cells));
	}

	// A read count must be a finite number that is not negative. pysyndna lets
	// anything else become NaN under log10 and silently voids the whole sample;
	// none of these is a degenerate-but-computable case, they are structurally
	// impossible, so miint refuses them -- a deliberate divergence consistent
	// with D23.
	//
	// One rule covering NaN, +/-inf and negatives rather than three, because the
	// alternative is remembering to skip the bad values again at every later
	// summation. `feature_totals` below is exactly that trap: a raw `+=` over a
	// NaN makes the total NaN, `NaN < min_syndna_counts` is false, and the
	// synDNA then becomes undroppable no matter how few reads it really has.
	std::vector<std::string> unusable_cells;
	for (const auto &row : counts) {
		if (!std::isfinite(row.count) || row.count < 0.0) {
			unusable_cells.push_back(row.sample_id + " / " + row.feature_id);
		}
	}
	if (!unusable_cells.empty()) {
		throw std::invalid_argument(
		    "read counts must be finite and not negative; offending (sample_id / feature_id): " +
		    FormatIdList(unusable_cells));
	}

	std::vector<std::string> unusable_concentrations;
	double total_ng_ul = 0.0;
	for (const auto &row : concentrations) {
		if (!std::isfinite(row.syndna_indiv_ng_ul) || row.syndna_indiv_ng_ul < 0.0) {
			unusable_concentrations.push_back(row.feature_id);
			continue;
		}
		total_ng_ul += row.syndna_indiv_ng_ul;
	}
	if (!unusable_concentrations.empty()) {
		// pandas' sum() would skip a NaN here and carry it into that synDNA's
		// fraction, voiding every sample the synDNA appears in with no
		// explanation. An infinite one is worse still: it becomes the whole
		// denominator, so EVERY other synDNA's fraction collapses to exactly
		// zero and the entire run returns no models at all. The pool
		// composition is configuration, not measurement -- a hole in it is a
		// mistake to report, not a value to propagate.
		throw std::invalid_argument(
		    "syndna_indiv_ng_ul must be present, finite, and not negative; offending feature_id " +
		    FormatIdList(unusable_concentrations));
	}
	// Finite as well as positive: every row above is finite, but enough large
	// ones still sum to +inf, and that would zero every fraction just as surely.
	if (!(total_ng_ul > 0.0) || !std::isfinite(total_ng_ul)) {
		throw std::invalid_argument("the synDNA concentrations sum to " + std::to_string(total_ng_ul) +
		                            "; the total must be positive and finite");
	}

	std::vector<std::string> count_features;
	std::vector<std::string> count_samples;
	count_features.reserve(counts.size());
	count_samples.reserve(counts.size());
	for (const auto &row : counts) {
		count_features.push_back(row.feature_id);
		count_samples.push_back(row.sample_id);
	}

	const auto unconfigured = IdsMissingFrom(count_features, concentration_ids);
	if (!unconfigured.empty()) {
		throw std::invalid_argument("the synDNA counts relation names feature_id " + FormatIdList(unconfigured) +
		                            ", which the concentrations relation does not describe, so their mass cannot be "
		                            "computed");
	}

	const auto unparameterized = IdsMissingFrom(count_samples, mass_ids);
	if (!unparameterized.empty()) {
		throw std::invalid_argument("the synDNA counts relation names sample_id " + FormatIdList(unparameterized) +
		                            ", which the sample parameters relation has no mass_syndna_input_ng row for");
	}
}

// Lookup maps over the four reference relations, built ONCE and shared by the
// validation pass and the computation that follows it.
//
// Shared rather than built twice because the two passes need exactly the same
// four maps, and a second construction is not just wasted work on relations the
// header describes as "a whole reference database" / "a whole study" -- it is
// two places that have to agree on how each map is keyed, with nothing to
// notice if they drift apart.
//
// Duplicate keys are recorded here rather than rediscovered, since inserting is
// what detects them. Each list holds one entry per occurrence past the first;
// FormatIdList deduplicates when rendering.
struct CellCountsIndex {
	std::unordered_map<std::string, double> length_by_feature;
	std::unordered_map<std::string, const SampleRegression *> model_by_sample;
	std::unordered_map<std::string, const SampleCellParams *> params_by_sample;
	std::map<std::pair<std::string, std::string>, double> coverage_by_cell;

	std::vector<std::string> repeated_lengths;
	std::vector<std::string> repeated_models;
	std::vector<std::string> repeated_params;
	std::vector<std::string> repeated_coverage;
};

// Borrows from `models` and `params`: the returned index holds pointers into
// them, so both must outlive it. Both are const references held by
// ComputeCellCounts for the duration, so this is safe by construction there.
CellCountsIndex BuildCellCountsIndex(const std::vector<SampleRegression> &models,
                                     const std::vector<CoverageObservation> &coverage,
                                     const std::vector<FeatureLength> &lengths,
                                     const std::vector<SampleCellParams> &params) {
	CellCountsIndex index;

	index.length_by_feature.reserve(lengths.size());
	for (const auto &row : lengths) {
		if (!index.length_by_feature.emplace(row.feature_id, row.ogu_len_in_bp).second) {
			index.repeated_lengths.push_back(row.feature_id);
		}
	}
	index.model_by_sample.reserve(models.size());
	for (const auto &row : models) {
		if (!index.model_by_sample.emplace(row.sample_id, &row).second) {
			index.repeated_models.push_back(row.sample_id);
		}
	}
	index.params_by_sample.reserve(params.size());
	for (const auto &row : params) {
		if (!index.params_by_sample.emplace(row.sample_id, &row).second) {
			index.repeated_params.push_back(row.sample_id);
		}
	}
	for (const auto &row : coverage) {
		if (!index.coverage_by_cell.emplace(std::make_pair(row.sample_id, row.feature_id), row.coverage).second) {
			index.repeated_coverage.push_back(row.sample_id + " / " + row.feature_id);
		}
	}
	return index;
}

// The factor that turns "per gram of sequenced gDNA" into "per unit of the
// material actually collected".
//
// extracted_gdna_mass_g is the gDNA recovered from the WHOLE extraction, not
// the portion that reached the sequencer, which is what makes the quotient a
// sample-side quantity rather than an instrument-side one.
//
// Exactly 1.0 for cells_per_g_of_gdna -- multiplying by it is the identity on
// every double, so that metric's values are M3's bit for bit. A branch in the
// per-cell loop would say the same thing less obviously.
double SampleUnitRatio(const SampleCellParams &params, CellCountsMetric metric) {
	if (metric == CellCountsMetric::CellsPerGOfGdna) {
		return 1.0;
	}
	const double extracted_gdna_mass_g =
	    params.extracted_gdna_concentration_ng_ul * params.vol_extracted_elution_ul / 1e9;
	return extracted_gdna_mass_g / params.sample_denominator;
}

// Reject every input ComputeCellCounts cannot compute a defined answer from,
// before any arithmetic runs. Same three tiers and the same no-prefix message
// contract as ValidateFitInputs above.
//
// Narrower than that function in one respect, deliberately: the four reference
// relations are screened only where the counts relation actually reaches them.
// ValidateFitInputs screens every concentration because it sums the whole
// relation into a denominator, so one bad row there corrupts every sample.
// Nothing here is summed, and these relations are routinely far wider than the
// query -- a lengths table is a whole reference database, a parameters table a
// whole study. Failing a query over twelve genomes because a thirteenth nobody
// asked about has a bad length would be its own kind of wrong answer.
void ValidateCellCountsInputs(const std::vector<CountObservation> &counts, const CellCountsIndex &index,
                              const CellCountsOptions &options) {
	// Written as !(x >= 0) rather than x < 0 so a NaN option is rejected too.
	if (!(options.min_coverage >= 0.0) || !(options.min_coverage <= 1.0)) {
		throw std::invalid_argument("min_coverage must be between 0 and 1 inclusive (got " +
		                            std::to_string(options.min_coverage) +
		                            "); coverage is a fraction here, not a percent");
	}
	// r^2 is a squared correlation and cannot leave [0, 1]. A threshold outside
	// that is not a strict filter, it is one that keeps every sample or none --
	// and pysyndna, which does not check, gives no sign which happened.
	if (!(options.min_rsquared >= 0.0) || !(options.min_rsquared <= 1.0)) {
		throw std::invalid_argument("min_rsquared must be between 0 and 1 inclusive (got " +
		                            std::to_string(options.min_rsquared) + ")");
	}

	const auto &length_by_feature = index.length_by_feature;
	const auto &model_by_sample = index.model_by_sample;
	const auto &params_by_sample = index.params_by_sample;
	const auto &coverage_by_cell = index.coverage_by_cell;

	if (!index.repeated_lengths.empty()) {
		throw std::invalid_argument("the feature lengths relation has more than one row for feature_id " +
		                            FormatIdList(index.repeated_lengths));
	}
	if (!index.repeated_models.empty()) {
		throw std::invalid_argument("the models relation has more than one row for sample_id " +
		                            FormatIdList(index.repeated_models));
	}
	if (!index.repeated_params.empty()) {
		throw std::invalid_argument("the sample parameters relation has more than one row for sample_id " +
		                            FormatIdList(index.repeated_params));
	}

	// The counts relation is keyed on the pair, like the coverage relation.
	// Compared as index-free string pairs rather than a packed key so no
	// separator character can collide with an id that legitimately contains it.
	std::set<std::pair<std::string, std::string>> seen_counts;
	std::vector<std::string> repeated_cells;
	for (const auto &row : counts) {
		if (!seen_counts.emplace(row.sample_id, row.feature_id).second) {
			repeated_cells.push_back(row.sample_id + " / " + row.feature_id);
		}
	}
	if (!repeated_cells.empty()) {
		throw std::invalid_argument("the counts relation has more than one row for the same (sample_id, feature_id): " +
		                            FormatIdList(repeated_cells));
	}

	if (!index.repeated_coverage.empty()) {
		throw std::invalid_argument(
		    "the coverage relation has more than one row for the same (sample_id, feature_id): " +
		    FormatIdList(index.repeated_coverage));
	}

	// One pass over the counts relation collects everything that depends on it,
	// so the reference relations are consulted exactly where they are used. The
	// lists are reported afterwards in tier order -- a malformed value before a
	// relation that is merely incomplete -- rather than in the order found.
	std::vector<std::string> unusable_cells;
	std::vector<std::string> percent_coverage;
	std::vector<std::string> bad_lengths;
	std::vector<std::string> zero_masses;
	std::vector<std::string> zero_denominators;
	std::vector<std::string> bad_rvalues;
	// nullptr for cells_per_g_of_gdna, whose callers leave sample_denominator at
	// zero because there is no column to fill it from. Checking it regardless
	// would make every query for that metric throw.
	const char *denominator_column = DenominatorColumnName(options.metric);
	std::vector<std::string> uncovered_cells;
	std::vector<std::string> unmeasured_features;
	std::vector<std::string> unparameterized_samples;
	for (const auto &row : counts) {
		// Stricter than the fit's rule, which admits a zero count because it
		// drops those points before taking the log. There is no such step here:
		// log10(0) is -inf, and under a NEGATIVE slope that returns as
		// 10^(+inf), an infinite cell count. The reader drops zero-valued cells
		// long before this, so the check is unreachable in practice -- it is
		// what keeps it that way if the reader ever changes.
		if (!std::isfinite(row.count) || !(row.count > 0.0)) {
			unusable_cells.push_back(row.sample_id + " / " + row.feature_id);
		}

		const auto covered = coverage_by_cell.find({row.sample_id, row.feature_id});
		if (covered == coverage_by_cell.end()) {
			uncovered_cells.push_back(row.sample_id + " / " + row.feature_id);
		} else if (!(covered->second >= 0.0) || !(covered->second <= 1.0)) {
			percent_coverage.push_back(row.sample_id + " / " + row.feature_id);
		}

		const auto length = length_by_feature.find(row.feature_id);
		if (length == length_by_feature.end()) {
			unmeasured_features.push_back(row.feature_id);
		} else if (!std::isfinite(length->second) || !(length->second > 0.0)) {
			bad_lengths.push_back(row.feature_id);
		}

		const auto sample_params = params_by_sample.find(row.sample_id);
		if (sample_params == params_by_sample.end()) {
			unparameterized_samples.push_back(row.sample_id);
		} else {
			if (sample_params->second->sequenced_sample_gdna_mass_ng == 0.0) {
				zero_masses.push_back(row.sample_id);
			}
			if (denominator_column != nullptr && sample_params->second->sample_denominator == 0.0) {
				zero_denominators.push_back(row.sample_id);
			}
		}

		// A sample with no model is a warning, not an error -- see the header --
		// so only a model that is present and impossible is caught here.
		const auto model = model_by_sample.find(row.sample_id);
		if (model != model_by_sample.end() && std::isfinite(model->second->rvalue) &&
		    !(model->second->rvalue >= -1.0 && model->second->rvalue <= 1.0)) {
			bad_rvalues.push_back(row.sample_id);
		}
	}

	if (!unusable_cells.empty()) {
		throw std::invalid_argument("read counts must be finite and positive; offending (sample_id / feature_id): " +
		                            FormatIdList(unusable_cells));
	}
	if (!percent_coverage.empty()) {
		// D9: miint is fraction-only. pysyndna takes fractions or percents and
		// explicitly leaves it to the caller to pass a matching min_coverage,
		// which turns the commonest mistake -- percents against a fractional
		// threshold -- into a filter that silently keeps everything.
		// "present" as well as in range: a SQL NULL arrives here as NaN and fails
		// the same test, and a message that only mentioned percents would send
		// that user looking for a unit problem they do not have.
		throw std::invalid_argument(
		    "coverage must be present and a fraction between 0 and 1 inclusive, not a percent; offending "
		    "(sample_id / feature_id): " +
		    FormatIdList(percent_coverage));
	}
	if (!bad_lengths.empty()) {
		// A denominator, divided by with no check at all on pysyndna's side: a
		// zero length yields inf cells for that feature in every sample and a
		// NaN one yields NaN, both silently.
		throw std::invalid_argument("ogu_len_in_bp must be finite and positive; offending feature_id " +
		                            FormatIdList(bad_lengths));
	}
	if (!zero_masses.empty()) {
		// The other denominator, and the one case pysyndna's parameter screen
		// cannot see: zero is neither NaN nor negative, so it passes the filter
		// and then divides through to inf cells with no log message. Same
		// treatment as D23 gives total_biological_reads_r1r2, for the same
		// reason. NaN, negative and infinite masses are NOT errors -- those the
		// filter drops, matching pysyndna.
		throw std::invalid_argument("sequenced_sample_gdna_mass_ng must not be zero; offending sample_id " +
		                            FormatIdList(zero_masses));
	}
	if (!zero_denominators.empty()) {
		// The sample-level metrics' denominator, and the same story a third
		// time. Named from DenominatorColumnName rather than spelled here, so
		// the column the message blames is the column the wrapper read -- a
		// message pointing at sample_volume_ul when the missing value was in
		// calc_mass_sample_aliquot_input_g is worse than no message at all.
		throw std::invalid_argument(std::string(denominator_column) + " must not be zero; offending sample_id " +
		                            FormatIdList(zero_denominators));
	}
	if (!bad_rvalues.empty()) {
		// Malformed models relation rather than a weak model: a correlation
		// coefficient cannot be 1.5, and squaring one gives an r^2 above 1 that
		// passes every legal min_rsquared. A non-finite rvalue is a different
		// thing -- an unusable model, reported as such and not an error.
		throw std::invalid_argument("rvalue must be between -1 and 1 inclusive; offending sample_id " +
		                            FormatIdList(bad_rvalues));
	}

	if (!unmeasured_features.empty()) {
		throw std::invalid_argument("the counts relation names feature_id " + FormatIdList(unmeasured_features) +
		                            ", which the feature lengths relation has no ogu_len_in_bp row for");
	}
	if (!unparameterized_samples.empty()) {
		throw std::invalid_argument("the counts relation names sample_id " + FormatIdList(unparameterized_samples) +
		                            ", which the sample parameters relation does not describe");
	}
	if (!uncovered_cells.empty()) {
		// Checked on the PAIR, which is where pysyndna's equivalent has a hole:
		// it validates the two axis id sets separately, so a sample and a
		// feature can each be present while their cell is not. Its left join
		// then yields NaN, `NaN >= min_coverage` is false so the cell is
		// dropped, and `NaN < min_coverage` is false too so it does not even
		// reach the too-low-coverage log. The cell disappears in silence.
		throw std::invalid_argument(
		    "the counts relation names (sample_id, feature_id) " + FormatIdList(uncovered_cells) +
		    ", which the coverage relation has no row for; every counted cell needs a coverage to be filtered on");
	}
}

} // namespace

LinregressResult Linregress(const std::vector<double> &x, const std::vector<double> &y) {
	LinregressResult out;
	const size_t n = x.size();
	// scipy raises on empty input; a length mismatch is a caller bug. Both are
	// "no model" here rather than exceptions, because every failure mode of this
	// function ends the same way -- the sample is omitted with a warning.
	if (n == 0 || y.size() != n) {
		return out;
	}

	// scipy raises when x is constant and n > 1 (np.amax(x) == np.amin(x)).
	// NaN propagates through amax/amin and NaN == NaN is false, so a NaN in x
	// SUPPRESSES that check on scipy's side; mirrored here so the two agree even
	// though log10 of a positive count can never actually produce one.
	if (n > 1) {
		bool has_nan = false;
		double lo = x[0];
		double hi = x[0];
		for (size_t i = 0; i < n; ++i) {
			if (std::isnan(x[i])) {
				has_nan = true;
				break;
			}
			lo = std::min(lo, x[i]);
			hi = std::max(hi, x[i]);
		}
		if (!has_nan && lo == hi) {
			return out;
		}
	}

	const double dn = static_cast<double>(n);
	double sum_x = 0.0;
	double sum_y = 0.0;
	for (size_t i = 0; i < n; ++i) {
		sum_x += x[i];
		sum_y += y[i];
	}
	const double xmean = sum_x / dn;
	const double ymean = sum_y / dn;

	double ssxm = 0.0;
	double ssxym = 0.0;
	double ssym = 0.0;
	for (size_t i = 0; i < n; ++i) {
		const double dx = x[i] - xmean;
		const double dy = y[i] - ymean;
		ssxm += dx * dx;
		ssxym += dx * dy;
		ssym += dy * dy;
	}
	// BIASED (population) covariances -- np.cov(x, y, bias=1) divides by n, not
	// n-1. Every other use below is a ratio in which the choice cancels, but
	// intercept_stderr multiplies by sqrt(ssxm + xmean^2), where ssxm stands
	// alone; the goldens therefore distinguish biased from unbiased.
	ssxm /= dn;
	ssxym /= dn;
	ssym /= dn;

	if (ssxm == 0.0 || ssym == 0.0) {
		// A degenerate spread has no correlation to report. scipy's `else 0.0`
		// arm needs ssxym != 0 with a zero variance, which exact arithmetic
		// cannot produce; it is reproduced anyway rather than assumed dead.
		out.rvalue = (ssxym == 0.0) ? kNaN : 0.0;
	} else {
		double r = ssxym / std::sqrt(ssxm * ssym);
		// Explicit clamp: rounding can push |r| a hair past 1, and everything
		// downstream takes sqrt(1 - r^2).
		if (r > 1.0) {
			r = 1.0;
		} else if (r < -1.0) {
			r = -1.0;
		}
		out.rvalue = r;
	}

	out.slope = ssxym / ssxm;
	out.intercept = ymean - out.slope * xmean;

	if (n == 2) {
		// scipy short-circuits two points: the line is exact, so both standard
		// errors are zero by construction and the p-value is degenerate. Note
		// the pvalue == 1.0 arm always coexists with a NaN rvalue (equal y means
		// ssym == 0), so it can never reach output -- it is reproduced for
		// fidelity, not because a caller will see it.
		out.pvalue = (y[0] == y[1]) ? 1.0 : 0.0;
		out.stderr_ = 0.0;
		out.intercept_stderr = 0.0;
	} else {
		const double df = dn - 2.0;
		const double t = out.rvalue * std::sqrt(df / ((1.0 - out.rvalue + kTiny) * (1.0 + out.rvalue + kTiny)));
		out.pvalue = 2.0 * StudentTSurvival(std::fabs(t), df);

		// Algebraically scipy's sqrt((1 - r^2)*ssym/ssxm/df): n*ssym*(1 - r^2) is
		// the residual sum of squares, so accumulating the residuals directly
		// computes the same quantity without the cancellation. See the header for
		// why that matters -- on a perfect fit the transcribed form is both wrong
		// against the golden and platform-dependent. Residuals are formed as
		// dy - slope*dx rather than y - (slope*x + intercept) so the sample means
		// never have to cancel back out.
		double resid_var = 0.0;
		for (size_t i = 0; i < n; ++i) {
			const double e = (y[i] - ymean) - out.slope * (x[i] - xmean);
			resid_var += e * e;
		}
		resid_var /= dn;
		if (std::isnan(out.rvalue)) {
			// scipy's form carries a NaN r straight into stderr; the residual
			// form would not (constant y gives exact zero residuals), so the
			// propagation is restored explicitly.
			resid_var = kNaN;
		}
		out.stderr_ = std::sqrt(resid_var / ssxm / df);
		out.intercept_stderr = out.stderr_ * std::sqrt(ssxm + xmean * xmean);
	}

	// pysyndna discards the whole sample if ANY field is NaN -- so a NaN rvalue
	// throws away an otherwise good slope (_convert_linregressresults_to_dict
	// :329-334 breaks out of the field loop and records None).
	out.ok = !(std::isnan(out.slope) || std::isnan(out.intercept) || std::isnan(out.rvalue) || std::isnan(out.pvalue) ||
	           std::isnan(out.stderr_) || std::isnan(out.intercept_stderr));
	return out;
}

FitResult FitSyndnaModels(const std::vector<CountObservation> &counts,
                          const std::vector<SyndnaConcentration> &concentrations, const std::vector<SampleMass> &masses,
                          const FitOptions &options) {
	ValidateFitInputs(counts, concentrations, masses, options);

	FitResult result;

	// The feature universe is every synDNA named by either relation, sorted.
	// Sorted because the order features are visited fixes the summation order
	// inside the fit and so its last few ULP, and SQL hands rows over in
	// whatever order the scan produced -- without this the same query could
	// return answers differing in the last bits.
	std::vector<std::string> feature_ids;
	feature_ids.reserve(counts.size() + concentrations.size());
	for (const auto &row : counts) {
		feature_ids.push_back(row.feature_id);
	}
	for (const auto &row : concentrations) {
		feature_ids.push_back(row.feature_id);
	}
	std::sort(feature_ids.begin(), feature_ids.end());
	feature_ids.erase(std::unique(feature_ids.begin(), feature_ids.end()), feature_ids.end());

	std::unordered_map<std::string, size_t> feature_pos;
	feature_pos.reserve(feature_ids.size());
	for (size_t i = 0; i < feature_ids.size(); ++i) {
		feature_pos[feature_ids[i]] = i;
	}

	// Totals across EVERY sample, including the ones the parameter filter is
	// about to remove -- see the ordering note in the header. A plain sum is
	// safe only because validation has already rejected every non-finite count;
	// pandas' equivalent skips NaN, and a raw `+=` that met one would make the
	// total NaN and the comparison below permanently false.
	std::vector<double> feature_totals(feature_ids.size(), 0.0);
	for (const auto &row : counts) {
		feature_totals[feature_pos.at(row.feature_id)] += row.count;
	}

	const double min_counts = static_cast<double>(options.min_syndna_counts);
	std::vector<bool> dropped(feature_ids.size(), false);
	for (size_t i = 0; i < feature_ids.size(); ++i) {
		// Strict `<`, so a synDNA sitting exactly on the threshold survives. A
		// synDNA named only by the concentrations relation has a total of zero
		// and lands here too, which is the outcome pysyndna reaches from an
		// all-zero row of its dense table.
		if (feature_totals[i] < min_counts) {
			dropped[i] = true;
			result.dropped_syndna_ids.push_back(feature_ids[i]);
		}
	}

	// Mass fraction of the pool contributed by each synDNA. The denominator sums
	// the FULL concentrations relation and is never rescaled after the drop
	// above. Validation has already established that every value present is a
	// non-negative number and that they do not all sum to zero; a feature with
	// no concentration row at all keeps its NaN, and validation has established
	// that no such feature carries counts, so the NaN never reaches a fit.
	std::vector<double> concentration(feature_ids.size(), kNaN);
	for (const auto &row : concentrations) {
		concentration[feature_pos.at(row.feature_id)] = row.syndna_indiv_ng_ul;
	}
	double total_ng_ul = 0.0;
	for (size_t i = 0; i < concentration.size(); ++i) {
		if (!std::isnan(concentration[i])) {
			total_ng_ul += concentration[i];
		}
	}
	std::vector<double> fraction_of_pool(feature_ids.size(), kNaN);
	for (size_t i = 0; i < concentration.size(); ++i) {
		fraction_of_pool[i] = concentration[i] / total_ng_ul;
	}

	std::unordered_map<std::string, double> mass_by_sample;
	mass_by_sample.reserve(masses.size());
	for (const auto &row : masses) {
		mass_by_sample[row.sample_id] = row.mass_syndna_input_ng;
	}

	// Group observations by sample in first-appearance order. Unlike the feature
	// order this one is arithmetically inert -- it decides only which output row
	// comes first -- so it is left as the caller presented it.
	std::vector<std::string> sample_ids;
	std::unordered_map<std::string, size_t> sample_pos;
	std::vector<std::vector<std::pair<size_t, double>>> observations;
	for (const auto &row : counts) {
		auto found = sample_pos.find(row.sample_id);
		size_t s;
		if (found == sample_pos.end()) {
			s = sample_ids.size();
			sample_pos.emplace(row.sample_id, s);
			sample_ids.push_back(row.sample_id);
			observations.emplace_back();
		} else {
			s = found->second;
		}
		observations[s].emplace_back(feature_pos.at(row.feature_id), row.count);
	}

	std::vector<double> log_counts;
	std::vector<double> log_mass_ng;
	for (size_t s = 0; s < sample_ids.size(); ++s) {
		// Validation has established that every sample in `counts` has a row
		// here, so a miss is impossible; the NaN keeps the lookup total rather
		// than relying on that.
		const auto found = mass_by_sample.find(sample_ids[s]);
		const double mass_ng = (found == mass_by_sample.end()) ? kNaN : found->second;
		if (!IsUsableSampleParameter(mass_ng)) {
			result.filtered_sample_ids.push_back(sample_ids[s]);
			continue;
		}
		const double pool_mass_ng = mass_ng * options.syndna_contributing_fraction;

		auto &sample_rows = observations[s];
		std::sort(sample_rows.begin(), sample_rows.end());
		log_counts.clear();
		log_mass_ng.clear();
		for (const auto &observation : sample_rows) {
			if (dropped[observation.first]) {
				continue;
			}
			// Zero counts must go before the log10, not after. Written as
			// !(count > 0) rather than count == 0 to match pandas' mask exactly;
			// validation has already rejected NaN, so the two spellings agree
			// here, but the mask form stays correct if that ever loosens.
			if (!(observation.second > 0.0)) {
				continue;
			}
			log_counts.push_back(std::log10(observation.second));
			log_mass_ng.push_back(std::log10(pool_mass_ng * fraction_of_pool[observation.first]));
		}
		if (log_counts.empty()) {
			result.unfittable_sample_ids.push_back(sample_ids[s]);
			continue;
		}
		const LinregressResult fit = Linregress(log_counts, log_mass_ng);
		if (!fit.ok) {
			result.unfittable_sample_ids.push_back(sample_ids[s]);
			continue;
		}
		result.models.push_back({sample_ids[s], fit});
	}

	// Metadata describing samples that carry no reads. Not an error -- a
	// parameters relation covering a whole run, queried for one plate, is
	// ordinary usage -- but reported, because the same shape is what a join on
	// the wrong key looks like. Reported in the order `masses` presented them,
	// which is the order the user wrote.
	for (const auto &row : masses) {
		if (sample_pos.find(row.sample_id) == sample_pos.end()) {
			result.samples_without_counts.push_back(row.sample_id);
		}
	}
	return result;
}

CellCountsResult ComputeCellCounts(const std::vector<CountObservation> &counts,
                                   const std::vector<SampleRegression> &models,
                                   const std::vector<CoverageObservation> &coverage,
                                   const std::vector<FeatureLength> &lengths,
                                   const std::vector<SampleCellParams> &params, const CellCountsOptions &options) {
	// Built once and used by both the validation pass and everything below it;
	// `index` borrows from `models` and `params`, which are alive for the whole
	// call. See CellCountsIndex.
	const CellCountsIndex index = BuildCellCountsIndex(models, coverage, lengths, params);
	ValidateCellCountsInputs(counts, index, options);

	CellCountsResult result;

	const auto &model_by_sample = index.model_by_sample;
	const auto &length_by_feature = index.length_by_feature;
	const auto &params_by_sample = index.params_by_sample;
	const auto &coverage_by_cell = index.coverage_by_cell;

	// Sorted and deduplicated, so the diagnostic lists below come out in an
	// order that does not depend on the order the counts rows arrived in.
	std::vector<std::string> sample_ids;
	sample_ids.reserve(counts.size());
	for (const auto &row : counts) {
		sample_ids.push_back(row.sample_id);
	}
	std::sort(sample_ids.begin(), sample_ids.end());
	sample_ids.erase(std::unique(sample_ids.begin(), sample_ids.end()), sample_ids.end());

	// Samples that will produce nothing, whatever the reason. Each one is also
	// named in exactly one of the result's lists as it is added here.
	std::unordered_set<std::string> discarded;

	// Filter 1 of 3: unusable parameters. This runs FIRST, and the order is
	// load-bearing for the diagnostics rather than for the values -- see the
	// header. pysyndna drops these samples out of the counts table itself
	// (calc_cell_counts.py:1013) before the coverage filter ever sees them.
	//
	// All three base columns are screened even though only the gDNA mass is
	// divided by, because pysyndna screens on REQUIRED_DNA_PREP_INFO_KEYS
	// regardless of the metric requested. Relaxing that would change which
	// samples appear while leaving every value identical.
	//
	// The metric's own denominator joins them, and no other denominator does:
	// pysyndna's filter set is those keys plus the REQUESTED metric's column
	// (calc_ogu_cell_counts_biom:959-965). A sample missing a column this query
	// never reads is not a sample this query cannot answer for.
	const bool screens_denominator = DenominatorColumnName(options.metric) != nullptr;
	for (const auto &sample_id : sample_ids) {
		// Validation has established that every sample in `counts` has a row
		// here; `at` rather than `find` so a future edit that breaks that throws
		// instead of quietly dropping the sample into no bucket at all.
		const SampleCellParams &sample_params = *params_by_sample.at(sample_id);
		if (!IsUsableSampleParameter(sample_params.sequenced_sample_gdna_mass_ng) ||
		    !IsUsableSampleParameter(sample_params.extracted_gdna_concentration_ng_ul) ||
		    !IsUsableSampleParameter(sample_params.vol_extracted_elution_ul) ||
		    (screens_denominator && !IsUsableSampleParameter(sample_params.sample_denominator))) {
			result.filtered_sample_ids.push_back(sample_id);
			discarded.insert(sample_id);
		}
	}

	// Filter 2 of 3: coverage, over every (sample, feature) still standing, and
	// before anything is computed. The comparison is pysyndna's
	// `coverage >= min_coverage` (calc_cell_counts.py:602) written as its
	// negation, so a feature sitting exactly on the threshold survives.
	//
	// Offenders are reported by feature rather than by cell, matching pysyndna's
	// own log line: a genome too poorly covered to trust is a property of the
	// assembly far more often than of one sample.
	std::vector<std::string> low_coverage;
	std::vector<const CountObservation *> surviving;
	std::unordered_set<std::string> samples_with_covered_cells;
	surviving.reserve(counts.size());
	for (const auto &row : counts) {
		if (discarded.count(row.sample_id) != 0) {
			continue;
		}
		const auto found = coverage_by_cell.find({row.sample_id, row.feature_id});
		const double cell_coverage = (found == coverage_by_cell.end()) ? kNaN : found->second;
		if (!(cell_coverage >= options.min_coverage)) {
			low_coverage.push_back(row.feature_id);
			continue;
		}
		surviving.push_back(&row);
		samples_with_covered_cells.insert(row.sample_id);
	}
	std::sort(low_coverage.begin(), low_coverage.end());
	low_coverage.erase(std::unique(low_coverage.begin(), low_coverage.end()), low_coverage.end());
	result.low_coverage_feature_ids = std::move(low_coverage);

	// Filter 3 of 3: the per-sample model gates. pysyndna reaches these only for
	// the sample ids still present in the working table after the coverage
	// filter (calc_cell_counts.py:505), so a sample with nothing left to compute
	// is never asked about its model -- which is why the wholly uncovered case
	// is settled first here and reported on its own rather than as a missing or
	// weak model it may well also be.
	for (const auto &sample_id : sample_ids) {
		if (discarded.count(sample_id) != 0) {
			continue;
		}
		if (samples_with_covered_cells.count(sample_id) == 0) {
			result.uncovered_sample_ids.push_back(sample_id);
			discarded.insert(sample_id);
			continue;
		}
		// "No model" and "a model none of whose fields is a number" are one
		// bucket, because they are one thing to pysyndna: a fit with any NaN
		// field is recorded as None (_convert_linregressresults_to_dict:329-334)
		// and reaches the same log line as a sample that was never fitted.
		// Keeping them together is what lets absquant_fit_models' output be fed
		// straight back in without projection.
		const auto found = model_by_sample.find(sample_id);
		if (found == model_by_sample.end() || !std::isfinite(found->second->slope) ||
		    !std::isfinite(found->second->intercept) || !std::isfinite(found->second->rvalue)) {
			result.samples_without_models.push_back(sample_id);
			discarded.insert(sample_id);
			continue;
		}
		// r^2, not |r|: pysyndna squares before comparing (calc_cell_counts.py
		// :680), so min_rsquared 0.25 admits a correlation of 0.5, and the two
		// spellings differ on every threshold except 0 and 1.
		const double rvalue = found->second->rvalue;
		if (rvalue * rvalue < options.min_rsquared) {
			result.low_rsquared_sample_ids.push_back(sample_id);
			discarded.insert(sample_id);
		}
	}

	std::unordered_set<std::string> samples_with_values;
	std::vector<std::string> overflowed_cells;
	for (const auto *row : surviving) {
		if (discarded.count(row->sample_id) != 0) {
			continue;
		}
		// Every lookup is guaranteed to hit: validation requires a length and a
		// parameter row for each counted cell, and a sample with no usable model
		// was discarded above. `at` rather than `find` so that a future edit
		// which breaks that invariant throws instead of computing with garbage.
		const SampleRegression &model = *model_by_sample.at(row->sample_id);
		const SampleCellParams &sample_params = *params_by_sample.at(row->sample_id);
		const double ogu_len_in_bp = length_by_feature.at(row->feature_id);

		// Predict this feature's gDNA mass from its read count, convert that
		// mass to genome copies, and normalize by the gDNA actually sequenced.
		// Read counts are guaranteed positive here -- see the header comment on
		// ComputeCellCounts for why that matters to the log10.
		const double gdna_mass_ng = std::pow(10.0, model.slope * std::log10(row->count) + model.intercept);
		const double genomes = gdna_mass_ng * kAvogadro / (ogu_len_in_bp * kGramsPerMolePerBasePair * 1e9);
		const double genomes_per_g_gdna = genomes / sample_params.sequenced_sample_gdna_mass_ng * 1e9;

		// One genome per cell: pysyndna's own explicit simplifying assumption.
		// The ratio then restates that per unit of collected sample, and is
		// exactly 1 for cells_per_g_of_gdna.
		const double value = genomes_per_g_gdna * SampleUnitRatio(sample_params, options.metric);
		if (!std::isfinite(value)) {
			// Every input reaching here is finite, but their combination need
			// not be: a model with a large positive intercept overflows the
			// 10^x, and an extraction concentration times an elution volume can
			// overflow before the divide. pysyndna emits the inf. Refusing is
			// the same call D23 and the zero-denominator check already make --
			// an infinite cell count is not a measurement, and a DOUBLE column
			// full of inf is worse than a failed query.
			overflowed_cells.push_back(row->sample_id + " / " + row->feature_id);
			continue;
		}
		if (value == 0.0) {
			// The sparse invariant (D10): a zero cell IS an absent cell, which
			// is exactly what pysyndna's dense output spells as an explicit
			// 0.0. So omitting one loses nothing and needs no diagnostic --
			// the sample is still in the output and still says what it is.
			// Only when EVERY cell of a sample goes this way does the omission
			// become ambiguous, and that case is reported below.
			continue;
		}
		samples_with_values.insert(row->sample_id);
		result.values.push_back({row->sample_id, row->feature_id, value});
	}
	if (!overflowed_cells.empty()) {
		throw std::invalid_argument(
		    "the predicted cell count overflowed to a non-finite value; check the model coefficients and the "
		    "extraction parameters for (sample_id / feature_id) " +
		    FormatIdList(overflowed_cells));
	}

	// A sample can pass every filter above and still emit nothing, if every one
	// of its cells came out zero. Reported rather than left to be inferred from
	// an absence: the sparse form has no way to distinguish "this sample is all
	// zeros" from "this sample was never here". Walking `sample_ids` keeps the
	// list sorted and independent of the order the counts arrived in.
	for (const auto &sample_id : sample_ids) {
		if (discarded.count(sample_id) == 0 && samples_with_values.count(sample_id) == 0) {
			result.zero_valued_sample_ids.push_back(sample_id);
		}
	}
	return result;
}

namespace {

// The number of bases an ORF spans, counting both endpoints.
//
// The absolute value is not defensive: woltka writes a reverse-strand ORF with
// start > end and pysyndna says so explicitly (quant_orfs.py:248), so a signed
// difference would hand four of the ten features in data/syndna/orf_coords.csv
// a negative length and a negative copy count. The +1 is inclusivity, and it is
// also what keeps a single-base ORF from dividing by zero.
double OrfLength(const OrfCoords &coords) {
	return std::fabs(coords.ogu_orf_end - coords.ogu_orf_start) + 1.0;
}

// A genome coordinate has to be a whole number of bases.
//
// pysyndna's cast_cols would accept 100.5 and hand back a fractional ORF
// length, which divides into Avogadro's number and yields a copy count nothing
// downstream can tell from a real one. Same call D9 makes on percent-versus-
// fraction coverage: the input that cannot mean what it says is refused rather
// than silently reinterpreted.
bool IsWholeNumber(double value) {
	return std::isfinite(value) && std::floor(value) == value;
}

// Index over the two reference relations, built once and shared by the
// validation pass and the arithmetic below it. Borrows from `coords` and
// `params`, so both must outlive it; ComputeOrfCopies holds them as const
// references for the whole call, so this is safe by construction there.
//
// Duplicate keys are recorded here rather than rediscovered, since inserting is
// what detects them. Each list holds one entry per occurrence past the first;
// FormatIdList deduplicates when rendering.
struct OrfCopiesIndex {
	std::unordered_map<std::string, const OrfCoords *> coords_by_feature;
	std::unordered_map<std::string, const SampleOrfParams *> params_by_sample;

	std::vector<std::string> repeated_coords;
	std::vector<std::string> repeated_params;
};

OrfCopiesIndex BuildOrfCopiesIndex(const std::vector<OrfCoords> &coords, const std::vector<SampleOrfParams> &params) {
	OrfCopiesIndex index;

	index.coords_by_feature.reserve(coords.size());
	for (const auto &row : coords) {
		if (!index.coords_by_feature.emplace(row.feature_id, &row).second) {
			index.repeated_coords.push_back(row.feature_id);
		}
	}
	index.params_by_sample.reserve(params.size());
	for (const auto &row : params) {
		if (!index.params_by_sample.emplace(row.sample_id, &row).second) {
			index.repeated_params.push_back(row.sample_id);
		}
	}
	return index;
}

// Reject every input ComputeOrfCopies cannot compute a defined answer from,
// before any arithmetic runs. Same tier order and the same no-prefix message
// contract as ValidateCellCountsInputs above: a malformed relation is reported
// before relations that merely disagree, because a fanned-out join produces
// both symptoms at once and only the first is the cause.
//
// There is no options tier here, this function having no options at all.
//
// The coordinates are screened only where the counts relation reaches them,
// for the reason ValidateCellCountsInputs gives about its own reference
// relations, and more sharply: a coords relation is a whole annotation, often
// every ORF of every genome in a reference database, against a counts relation
// naming the handful that sequenced.
void ValidateOrfCopiesInputs(const std::vector<CountObservation> &counts, const OrfCopiesIndex &index) {
	if (!index.repeated_coords.empty()) {
		throw std::invalid_argument("the ORF coordinates relation has more than one row for feature_id " +
		                            FormatIdList(index.repeated_coords));
	}
	if (!index.repeated_params.empty()) {
		throw std::invalid_argument("the sample parameters relation has more than one row for sample_id " +
		                            FormatIdList(index.repeated_params));
	}

	// Keyed on the pair, and compared as index-free string pairs rather than a
	// packed key so no separator character can collide with an id that
	// legitimately contains it. pysyndna has no equivalent check because its
	// input is a biom table, unique by construction.
	std::set<std::pair<std::string, std::string>> seen_counts;
	std::vector<std::string> repeated_cells;
	for (const auto &row : counts) {
		if (!seen_counts.emplace(row.sample_id, row.feature_id).second) {
			repeated_cells.push_back(row.sample_id + " / " + row.feature_id);
		}
	}
	if (!repeated_cells.empty()) {
		throw std::invalid_argument("the counts relation has more than one row for the same (sample_id, feature_id): " +
		                            FormatIdList(repeated_cells));
	}

	// One pass over the counts relation collects everything that depends on it,
	// so the reference relations are consulted exactly where they are used. The
	// lists are reported afterwards in tier order rather than in the order found.
	std::vector<std::string> unusable_cells;
	std::vector<std::string> bad_coords;
	std::vector<std::string> zero_reads;
	std::vector<std::string> zero_masses;
	std::vector<std::string> uncharted_features;
	std::vector<std::string> unparameterized_samples;
	for (const auto &row : counts) {
		// Looser than the cell-count rule, which refuses a zero count because
		// log10(0) under a negative slope comes back as an infinite answer.
		// Nothing here takes a logarithm, so a zero count simply means no copies
		// of that ORF -- refusing it would be a rule invented for no reason.
		if (!std::isfinite(row.count) || row.count < 0.0) {
			unusable_cells.push_back(row.sample_id + " / " + row.feature_id);
		}

		const auto orf = index.coords_by_feature.find(row.feature_id);
		if (orf == index.coords_by_feature.end()) {
			uncharted_features.push_back(row.feature_id);
		} else if (!IsWholeNumber(orf->second->ogu_orf_start) || !IsWholeNumber(orf->second->ogu_orf_end)) {
			// Negative coordinates are deliberately NOT refused: |end - start| + 1
			// is well defined for any pair of whole numbers, and woltka's format
			// says nothing about the origin.
			bad_coords.push_back(row.feature_id);
		}

		const auto sample_params = index.params_by_sample.find(row.sample_id);
		if (sample_params == index.params_by_sample.end()) {
			unparameterized_samples.push_back(row.sample_id);
		} else {
			if (sample_params->second->total_biological_reads_r1r2 == 0.0) {
				zero_reads.push_back(row.sample_id);
			}
			if (sample_params->second->calc_mass_sample_aliquot_input_g == 0.0) {
				zero_masses.push_back(row.sample_id);
			}
		}
	}

	if (!unusable_cells.empty()) {
		throw std::invalid_argument(
		    "read counts must be finite and not negative; offending (sample_id / feature_id): " +
		    FormatIdList(unusable_cells));
	}
	if (!bad_coords.empty()) {
		throw std::invalid_argument(
		    "ogu_orf_start and ogu_orf_end must be finite whole numbers of bases; offending feature_id " +
		    FormatIdList(bad_coords));
	}
	if (!zero_reads.empty()) {
		// The denominator D23 is about, and the one case pysyndna's parameter
		// screen cannot see: zero is neither NaN nor negative, so it passes the
		// filter and then divides through to an infinite copy count for every ORF
		// in the sample, with no log message. NaN, negative and infinite values
		// are NOT errors -- those the filter drops, matching pysyndna.
		throw std::invalid_argument("total_biological_reads_r1r2 must not be zero; offending sample_id " +
		                            FormatIdList(zero_reads));
	}
	if (!zero_masses.empty()) {
		// The other denominator, same story. Between them these are the whole of
		// the zero-is-an-error rule here: the remaining two parameter columns are
		// multiplied through, where zero is ordinary data.
		throw std::invalid_argument("calc_mass_sample_aliquot_input_g must not be zero; offending sample_id " +
		                            FormatIdList(zero_masses));
	}

	if (!uncharted_features.empty()) {
		// pysyndna does not validate this direction at all: it reaches the
		// coordinates through a bare `.at[]` inside a biom transform
		// (quant_orfs.py:158-161), so a counted ORF the annotation does not
		// describe surfaces as a pandas KeyError from inside a lambda.
		throw std::invalid_argument("the counts relation names feature_id " + FormatIdList(uncharted_features) +
		                            ", which the ORF coordinates relation has no ogu_orf_start/ogu_orf_end row for");
	}
	if (!unparameterized_samples.empty()) {
		throw std::invalid_argument("the counts relation names sample_id " + FormatIdList(unparameterized_samples) +
		                            ", which the sample parameters relation does not describe");
	}
}

} // namespace

OrfCopiesResult ComputeOrfCopies(const std::vector<CountObservation> &counts, const std::vector<OrfCoords> &coords,
                                 const std::vector<SampleOrfParams> &params) {
	// Built once and used by both the validation pass and everything below it;
	// `index` borrows from `coords` and `params`, which are alive for the whole
	// call. See OrfCopiesIndex.
	const OrfCopiesIndex index = BuildOrfCopiesIndex(coords, params);
	ValidateOrfCopiesInputs(counts, index);

	OrfCopiesResult result;

	const auto &coords_by_feature = index.coords_by_feature;
	const auto &params_by_sample = index.params_by_sample;

	// Sorted and deduplicated, so filtered_sample_ids comes out in an order that
	// does not depend on the order the counts rows arrived in.
	std::vector<std::string> sample_ids;
	sample_ids.reserve(counts.size());
	for (const auto &row : counts) {
		sample_ids.push_back(row.sample_id);
	}
	std::sort(sample_ids.begin(), sample_ids.end());
	sample_ids.erase(std::unique(sample_ids.begin(), sample_ids.end()), sample_ids.end());

	// The only filter this function has. All four columns are screened, which is
	// pysyndna's REQUIRED_PARAM_KEYS -- the union of its sample-info and RNA-prep
	// lists, passed whole to filter_data_by_sample_info (quant_orfs.py:313-317).
	// Unlike the cell-count parameters there is no column here carried purely for
	// parity: every one of the four is multiplied or divided through below.
	std::unordered_set<std::string> discarded;
	for (const auto &sample_id : sample_ids) {
		// Validation has established that every sample in `counts` has a row
		// here; `at` rather than `find` so a future edit that breaks that throws
		// instead of quietly dropping the sample into no bucket at all.
		const SampleOrfParams &sample_params = *params_by_sample.at(sample_id);
		if (!IsUsableSampleParameter(sample_params.calc_mass_sample_aliquot_input_g) ||
		    !IsUsableSampleParameter(sample_params.total_rna_concentration_ng_ul) ||
		    !IsUsableSampleParameter(sample_params.vol_extracted_elution_ul) ||
		    !IsUsableSampleParameter(sample_params.total_biological_reads_r1r2)) {
			result.filtered_sample_ids.push_back(sample_id);
			discarded.insert(sample_id);
		}
	}

	std::unordered_set<std::string> samples_with_values;
	std::vector<std::string> overflowed_cells;
	for (const auto &row : counts) {
		if (discarded.count(row.sample_id) != 0) {
			continue;
		}
		// Every lookup is guaranteed to hit: validation requires coordinates for
		// each counted feature and a parameter row for each counted sample. `at`
		// rather than `find` so that a future edit which breaks that invariant
		// throws instead of computing with garbage.
		const SampleOrfParams &sample_params = *params_by_sample.at(row.sample_id);
		const OrfCoords &orf = *coords_by_feature.at(row.feature_id);

		// pysyndna applies these as four successive biom transforms
		// (quant_orfs.py:113-166), and the association order is preserved rather
		// than rearranged into something tidier: the goldens reproduce bit for
		// bit this way, and any regrouping would put them an ulp out and turn a
		// real regression later into a judgement call about tolerance.
		const double copies_per_g_ssrna = kAvogadro / (OrfLength(orf) * kGramsPerMolePerRnaBase);
		const double fraction_of_reads = row.count / sample_params.total_biological_reads_r1r2;
		const double g_total_ssrna =
		    sample_params.total_rna_concentration_ng_ul * sample_params.vol_extracted_elution_ul / 1e9;

		const double value =
		    fraction_of_reads * g_total_ssrna * copies_per_g_ssrna / sample_params.calc_mass_sample_aliquot_input_g;
		if (!std::isfinite(value)) {
			// Every input reaching here is finite, but their combination need
			// not be: the extracted ssRNA mass is a product formed BEFORE the
			// ng->g divide, so two large-but-legal values overflow it to +inf and
			// carry that through the rest of the chain. pysyndna emits the inf.
			// Refusing is the same call the zero-denominator checks already make
			// -- an infinite copy count is not a measurement, and a DOUBLE column
			// full of inf is worse than a failed query.
			overflowed_cells.push_back(row.sample_id + " / " + row.feature_id);
			continue;
		}
		if (value == 0.0) {
			// The sparse invariant (D10): a zero cell IS an absent cell, which is
			// exactly what pysyndna's dense output spells as an explicit 0.0. So
			// omitting one loses nothing and needs no diagnostic. Only when EVERY
			// cell of a sample goes this way does the omission become ambiguous,
			// and that case is reported below.
			continue;
		}
		samples_with_values.insert(row.sample_id);
		result.values.push_back({row.sample_id, row.feature_id, value});
	}
	if (!overflowed_cells.empty()) {
		throw std::invalid_argument("the predicted copy count overflowed to a non-finite value; check the ORF "
		                            "coordinates and the extraction parameters for (sample_id / feature_id) " +
		                            FormatIdList(overflowed_cells));
	}

	// A sample can pass the filter and still emit nothing, if every one of its
	// cells came out zero. Reported rather than left to be inferred from an
	// absence. Walking `sample_ids` keeps the list sorted and independent of the
	// order the counts arrived in.
	for (const auto &sample_id : sample_ids) {
		if (discarded.count(sample_id) == 0 && samples_with_values.count(sample_id) == 0) {
			result.zero_valued_sample_ids.push_back(sample_id);
		}
	}
	return result;
}

} // namespace absquant
} // namespace miint
