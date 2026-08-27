#include "KsTwoSample.hpp"
#include "alignment_functions_internal.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <numeric>
#include <string>
#include <utility>

namespace miint {

namespace {

std::string Lowered(const std::string &s) {
	std::string out(s);
	std::transform(out.begin(), out.end(), out.begin(),
	               [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
	return out;
}

// floor/ceil division with a POSITIVE divisor and a possibly-negative numerator.
// C++ integer division truncates toward zero, which is neither of these for
// negative numerators -- and the band's lower edge goes negative for small
// columns, so getting this wrong silently widens the band near (0,0).
int64_t FloorDiv(int64_t a, int64_t b) {
	const int64_t q = a / b;
	return (a % b != 0 && a < 0) ? q - 1 : q;
}

int64_t CeilDiv(int64_t a, int64_t b) {
	const int64_t q = a / b;
	return (a % b != 0 && a > 0) ? q + 1 : q;
}

void RejectNaNIn(const std::vector<double> &v, const char *which) {
	for (const double x : v) {
		if (std::isnan(x)) {
			throw InvalidInputException(std::string("ks_2samp: sample ") + which +
			                            " contains NaN. NaN has no position in the sort order, so the ECDF -- and "
			                            "therefore the statistic -- would be undefined. Filter or impute it first.");
		}
	}
}

} // namespace

void RejectKsNaN(const std::vector<double> &a, const std::vector<double> &b) {
	RejectNaNIn(a, "a");
	RejectNaNIn(b, "b");
}

void ValidateKsMethod(const std::string &method) {
	const std::string m = Lowered(method);
	if (m == "auto" || m == "exact") {
		return;
	}
	if (m == "asymp") {
		throw InvalidInputException(
		    "ks_2samp: method 'asymp' is not implemented. SciPy's asymptotic branch is not the Kolmogorov series -- it "
		    "is the exact one-sample KS distribution at an effective sample size, selected per region among the Durbin "
		    "matrix, Pomeranz and Pelz-Good expansions -- so an approximation here would silently disagree with "
		    "previously published values. Use 'auto' or 'exact' (both exact, for samples up to 10000); see issue "
		    "#218.");
	}
	throw InvalidInputException("ks_2samp: unknown method '" + method + "' (must be one of 'auto', 'exact')");
}

double KsStatistic(std::vector<double> &a, std::vector<double> &b) {
	std::sort(a.begin(), a.end());
	std::sort(b.begin(), b.end());

	const size_t na = a.size();
	const size_t nb = b.size();
	size_t i = 0;
	size_t j = 0;
	double best = 0.0;

	// Advance to the next distinct pooled value, consume EVERY observation equal to
	// it in both samples, and only then compare the ECDFs. Comparing mid-tie would
	// report a jump that neither ECDF actually has.
	while (i < na && j < nb) {
		const double v = std::min(a[i], b[j]);
		while (i < na && a[i] == v) {
			i++;
		}
		while (j < nb && b[j] == v) {
			j++;
		}
		const double diff =
		    static_cast<double>(i) / static_cast<double>(na) - static_cast<double>(j) / static_cast<double>(nb);
		best = std::max(best, std::fabs(diff));
	}

	// No tail pass is needed. The loop ends with one sample exhausted, so that ECDF
	// is already 1 and the other's can only rise toward 1 over the remaining values
	// -- the gap shrinks monotonically from the value just recorded.
	//
	// `best` is not yet D. The sweep rounds three times per step -- two divisions and
	// a subtraction -- so it lands within a few ULP of the answer rather than on it:
	// at na=1, nb=3 the peak 1.0/1.0 - 1.0/3.0 comes out one ULP above 2/3. Every
	// achievable D is an exact multiple of 1/lcm(na,nb), so recover WHICH multiple and
	// divide once. A single IEEE division of two exactly-representable integers is
	// correctly rounded by construction, which makes the result the true D to the bit
	// (issue #256). Exact representability holds while lcm <= 2^53; KS_MAX_EXACT_N
	// keeps it under 1e8.
	const int64_t n1 = static_cast<int64_t>(na);
	const int64_t n2 = static_cast<int64_t>(nb);
	const int64_t lcm = n1 / std::gcd(n1, n2) * n2;
	return static_cast<double>(std::llround(best * static_cast<double>(lcm))) / static_cast<double>(lcm);
}

double KsExactPValue(int64_t n1, int64_t n2, double d) {
	// Transposing the lattice maps every path to a path and leaves |i*ng - j*mg|
	// unchanged, so the answer is symmetric in the two sizes -- but only exactly so
	// if the arithmetic runs in the same order. Canonicalising to n1 >= n2 buys two
	// things: ks_2samp(a,b) and ks_2samp(b,a) agree to the bit rather than to 1 ULP,
	// and the working column becomes O(min(n1,n2)) instead of O(max(n1,n2)), which
	// for a 2-vs-10000 comparison is three doubles rather than ten thousand.
	if (n1 < n2) {
		std::swap(n1, n2);
	}

	const int64_t g = std::gcd(n1, n2);
	const int64_t mg = n1 / g;
	const int64_t ng = n2 / g;
	const int64_t lcm = mg * n2;

	// D is always an exact multiple of 1/lcm; recover which one. Rounding rather
	// than truncating matters because d arrives as a double ratio.
	const int64_t h = static_cast<int64_t>(std::llround(d * static_cast<double>(lcm)));
	if (h <= 0) {
		return 1.0;
	}

	// Walk the interleaving as a random sequence of steps: at (i,j), with n1-i
	// observations of A and n2-j of B still unplaced, the next one comes from A with
	// probability (n1-i)/(n1+n2-i-j). Make the OUTSIDE of the band absorbing and
	// accumulate the mass crossing into it, so the p-value is a sum of POSITIVE
	// terms. Computing it as 1 - P(stay inside) instead would cancel: at n1=n2=20
	// with D=1 the answer is 1.45e-11, and subtracting from 1 leaves only five
	// significant digits (measured 4.2e-06 relative error) exactly where the
	// p-value matters most.
	//
	// f is one column of the lattice, updated in place: reading f[j] before writing
	// it gives column i-1's value (the +x predecessor) while f[j-1] has already
	// become column i's (the +y predecessor).
	std::vector<double> f(static_cast<size_t>(n2 + 1), 0.0);
	// Seed the origin. h >= 1 here, so (0,0) is always inside the band, and the sweep
	// visits it exactly once -- so pre-seeding is identical to special-casing it in the
	// recurrence, and it keeps the inner loop to a single uniform expression.
	f[0] = 1.0;
	double escape = 0.0;
	int64_t jlo_prev = 0;

	for (int64_t i = 0; i <= n1; i++) {
		// Inside the band means |i*ng - j*mg| < h, i.e. j strictly between
		// (i*ng - h)/mg and (i*ng + h)/mg.
		const int64_t jlo = std::max<int64_t>(FloorDiv(i * ng - h, mg) + 1, 0);
		const int64_t jhi = std::min<int64_t>(CeilDiv(i * ng + h, mg) - 1, n2);

		// Rows worth visiting in this column: the band itself, the row just above it
		// (reachable by a +y step out of the band), and the rows the band has left
		// behind since the previous column (reachable by a +x step out of it). Both
		// edges rise monotonically with i, so that is one contiguous range.
		const int64_t vlo = std::max<int64_t>(jlo_prev, 0);
		const int64_t vhi = std::min<int64_t>(jhi + 1, n2);

		for (int64_t j = vlo; j <= vhi; j++) {
			if (i == 0 && j == 0) {
				continue; // seeded above
			}

			const double denom = static_cast<double>(n1 + n2 - i - j + 1);
			double inflow = 0.0;
			if (i > 0) {
				inflow += f[static_cast<size_t>(j)] * static_cast<double>(n1 - i + 1) / denom;
			}
			if (j > 0) {
				inflow += f[static_cast<size_t>(j - 1)] * static_cast<double>(n2 - j + 1) / denom;
			}

			if (j >= jlo && j <= jhi) {
				f[static_cast<size_t>(j)] = inflow;
			} else {
				escape += inflow;
				f[static_cast<size_t>(j)] = 0.0;
			}
		}

		jlo_prev = jlo;
	}

	// A sum of probabilities can land a rounding step above 1 -- SciPy hits exactly
	// this at n1=n2=5 with D=0.2, where its own exact routine returns
	// 1.0000000000000002, declares failure and falls back to its asymptotic branch
	// before clipping to 1.0. Clamping reaches the same answer without the detour.
	return std::min(std::max(escape, 0.0), 1.0);
}

KsTwoSampleResult KsTwoSample(std::vector<double> a, std::vector<double> b) {
	if (a.empty() || b.empty()) {
		throw InvalidInputException(
		    "ks_2samp: both samples must not be empty (got " + std::to_string(a.size()) + " and " +
		    std::to_string(b.size()) +
		    " observations). The statistic is undefined against an empty sample; ks_2samp returns NULL for this rather "
		    "than 0.0, which would read as 'the two samples are identical'.");
	}

	RejectKsNaN(a, b);

	const int64_t na = static_cast<int64_t>(a.size());
	const int64_t nb = static_cast<int64_t>(b.size());
	if (na > KS_MAX_EXACT_N || nb > KS_MAX_EXACT_N) {
		throw InvalidInputException(
		    "ks_2samp: sample sizes " + std::to_string(na) + " and " + std::to_string(nb) + " exceed the maximum of " +
		    std::to_string(KS_MAX_EXACT_N) +
		    ". miint computes only the exact p-value and caps it at this size to bound the cost of the lattice sweep. "
		    "Two different things meet at this number, and only one of them is SciPy's: for method 'auto' SciPy also "
		    "stops here (its MAX_AUTO_N) and switches to 'asymp', which miint does not implement because approximating "
		    "it would silently disagree with previously published values; for method 'exact' SciPy has NO such cap, so "
		    "there the ceiling is miint's own choice and could be raised. See issue #218.");
	}

	KsTwoSampleResult result;
	result.statistic = KsStatistic(a, b);
	result.pvalue = KsExactPValue(na, nb, result.statistic);
	return result;
}

} // namespace miint
