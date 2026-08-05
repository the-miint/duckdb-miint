#include "absquant.hpp"

#include <algorithm>
#include <cmath>
#include <limits>
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

} // namespace absquant
} // namespace miint
