#pragma once

#include <cstdint>
#include <string>
#include <vector>

namespace miint {

struct KsTwoSampleResult {
	double statistic; // two-sided D = sup|F_a(x) - F_b(x)|, in [0, 1]
	double pvalue;    // P(D_{n1,n2} >= statistic) under H0, exact
};

// Largest sample size for which ks_2samp answers at all (issue #218).
//
// This is SciPy's MAX_AUTO_N. Above it, scipy.stats.ks_2samp switches to its
// asymptotic branch, and that branch is NOT the textbook Kolmogorov series
// 2*sum (-1)^(k-1) exp(-2 k^2 lambda^2). It is
//
//     kstwo.sf(d, round(n1*n2/(n1+n2)))
//
// -- the EXACT ONE-SAMPLE KS distribution at an effective sample size, which SciPy
// evaluates by choosing per (n, d) region among the Durbin matrix
// (Marsaglia-Tsang-Wang 2003), Pomeranz (1974), and the Pelz-Good (1976) asymptotic
// expansion. Measured on this machine against SciPy 1.18.0: the plain Kolmogorov
// series is 0.7% off in the most favourable reachable case (n1=n2=20000, d=0.01)
// and 47% off at small effective n, which method='auto' does reach (n1=20000,
// n2=5 gives an effective n of 5). So it is not a substitute, and reproducing
// SciPy here means reproducing its region selection -- including the regions where
// SciPy deliberately returns an approximation rather than the truth.
//
// Rather than return a p-value that silently disagrees with previously published
// micov numbers, ks_2samp raises above this size. micov never approaches it: its n
// is the number of SAMPLES in a metadata group, because cumulative_coverage_curve
// emits one rank per sample, i.e. n of 10-100 against a ceiling of 10000. If the
// asymptotic branch is ever needed, issue #218 records what it would take.
//
// One honest caveat about reusing the same number for both methods: MAX_AUTO_N gates
// only SciPy's method='auto'. SciPy's method='exact' has NO size cap -- it computes
// the exact p-value at any size and bails only when lcm(n1,n2) >= 2^31. So for
// 'exact' this ceiling is miint's own cost bound, not SciPy behaviour, and it is a
// conservative one: the lcm at 10000 x 9999 is 9.999e7, three orders below SciPy's
// own guard, and the sweep is O(n1 * band) rather than O(n1 * n2). Raising it would
// need nothing but a cost decision and fresh timings -- no approximation is involved.
constexpr int64_t KS_MAX_EXACT_N = 10000;

// Rejects a `method` that ks_2samp cannot honour. 'auto' and 'exact' are accepted
// (case-insensitively) and behave identically, because the exact distribution is
// the only one implemented -- see KS_MAX_EXACT_N. 'asymp' gets its own message
// saying so rather than being lumped in with a typo.
//
// Throws miint::InvalidInputException.
void ValidateKsMethod(const std::string &method);

// Two-sided KS statistic D = sup|F_a(x) - F_b(x)|. Sorts both vectors in place.
//
// The ECDFs are compared only AFTER every observation equal to the current pooled
// value has been consumed. That is what makes D == 0 for a sample compared against
// itself when it contains repeats; comparing mid-tie would report a spurious jump.
// Matches SciPy's searchsorted(..., side='right') formulation.
//
// Requires both vectors non-empty and NaN-free; the caller checks that.
double KsStatistic(std::vector<double> &a, std::vector<double> &b);

// P(D_{n1,n2} >= d) under H0 that both samples come from the same continuous
// distribution, computed exactly.
//
// Requires n1, n2 >= 1 and d in [0, 1]. Two further conditions the signature cannot
// express, and which the SQL path satisfies for free but a direct C++ caller must not
// assume:
//
//   * KS_MAX_EXACT_N is NOT enforced here -- KsTwoSample enforces it. This function
//     allocates min(n1,n2)+1 doubles and sweeps O(n1 * band), so calling it with
//     n1 = n2 = 1e6 will run for a very long time rather than reject the input.
//   * `d` is SNAPPED to the nearest achievable statistic, llround(d * lcm)/lcm, not
//     validated against it. D is always an exact multiple of 1/lcm(n1,n2), so any
//     d produced by KsStatistic is already one; but KsExactPValue(5, 5, 0.3) quietly
//     answers for D = 0.4 instead of complaining, because 0.3 is not achievable at
//     n1 = n2 = 5.
//
// Under H0 both samples are assumed drawn from a CONTINUOUS distribution, so that all
// C(n1+n2, n1) interleavings are equally likely. Ties violate that assumption, which
// makes the p-value conservative rather than wrong -- see the note in docs/diversity.md.
double KsExactPValue(int64_t n1, int64_t n2, double d);

// Throws miint::InvalidInputException if either sample contains a NaN, naming which.
//
// Exposed separately because the SQL wrapper turns an empty sample into NULL and must
// still report a NaN in the OTHER sample: ks_2samp([], ['nan']) is a data error, not an
// absence, and the empty-to-NULL short circuit would otherwise swallow it. KsTwoSample
// calls this too, so a direct C++ caller cannot skip it.
void RejectKsNaN(const std::vector<double> &a, const std::vector<double> &b);

// The whole test. NULL elements are dropped by the caller, so a and b here are the
// surviving observations.
//
// Throws miint::InvalidInputException if either sample is empty, contains a NaN, or
// is larger than KS_MAX_EXACT_N.
KsTwoSampleResult KsTwoSample(std::vector<double> a, std::vector<double> b);

} // namespace miint
