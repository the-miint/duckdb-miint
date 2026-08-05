#pragma once
//
// Pure (DuckDB-free) core for the absquant_* absolute-quantification functions.
//
// Reimplements the synDNA spike-in method of Zaramela et al. 2022 (mSystems 7(6)
// e0044722), as realized by pysyndna (biocore, BSD-3-Clause, (c) Amanda
// Birmingham). See THIRD_PARTY_LICENSES.md; duckdb-miint never invokes pysyndna,
// and the parity goldens under data/syndna/ are numbers only.
//
// Header/impl split with NO DuckDB dependency so the Catch2 unit-test target
// links it without libduckdb (mirrors community_distances.{hpp,cpp} and
// simulate_resemblance.{hpp,cpp}). The DuckDB table-function wrappers live in
// absquant_function.cpp.

#include <cstdint>
#include <string>
#include <vector>

namespace miint {
namespace absquant {

//! Regularized incomplete beta I_x(a, b), for a > 0, b > 0 and x in [0, 1].
//!
//! Exposed rather than kept file-local because it IS the numerically delicate
//! part -- StudentTSurvival is a thin change of variables over it -- and testing
//! it directly against known values is far sharper than only testing it through
//! the survival function.
//!
//! Evaluated by the classical continued fraction (DLMF 8.17.22) using Lentz's
//! method, with the standard x > (a+1)/(a+b+2) reflection to keep convergence
//! fast. Returns NaN for out-of-domain input.
double RegularizedIncompleteBeta(double a, double b, double x);

//! As above, but with 1-x supplied by the caller instead of computed here.
//!
//! Not a micro-optimization -- it is the difference between a right and a wrong
//! answer at the extremes. Near t = 0 the Student-t change of variables gives
//! x = 1/(1+t^2/df), which ROUNDS TO EXACTLY 1.0 once t^2/df falls below ~1e-16,
//! taking the whole result with it: sf(1e-8, 1) came back as exactly 0.5 rather
//! than 0.4999999952568131. The information survives only in the complement
//! (1e-16 here), which the caller can compute accurately but this function
//! could never recover from a rounded x.
//!
//! `x + one_minus_x` should be 1 to within rounding; the two are used in
//! separate log terms and are never added, so no accuracy is lost either way.
double RegularizedIncompleteBetaC(double a, double b, double x, double one_minus_x);

//! Survival function (upper tail, 1 - CDF) of Student's t with `df` degrees of
//! freedom: P(T > t).
//!
//! This is the one quantity in the whole port that DuckDB cannot supply -- the
//! regr_* aggregate family covers slope/intercept/r/stderr, but a p-value needs
//! an incomplete beta. It is reached via
//!
//!     sf(t, df) = 0.5 * I_{df/(df+t^2)}(df/2, 1/2)          for t >= 0
//!     sf(t, df) = 1 - sf(-t, df)                            for t < 0
//!
//! which reproduces scipy.stats.t.sf bit-for-bit on the committed grid
//! (data/syndna/studentt_sf_oracle.csv, generated from scipy 1.17.1). scipy is
//! BSD-3-Clause and is credited as the reference implementation; the continued
//! fraction here is written from the standard definition rather than
//! transcribed from scipy's vendored Cephes.
//!
//! `df` must be > 0. Returns NaN otherwise, and for NaN `t`.
double StudentTSurvival(double t, double df);

//! The six fields scipy.stats.linregress returns, in its order.
//!
//! `ok == false` means no usable model: either the fit raised (constant x with
//! n > 1, which scipy rejects) or some field came out NaN. pysyndna represents
//! those two cases differently (key present with None vs key absent); miint
//! collapses both to "omitted, with a warning" -- see the divergence ledger.
struct LinregressResult {
	double slope = 0.0;
	double intercept = 0.0;
	double rvalue = 0.0;
	double pvalue = 0.0;
	double stderr_ = 0.0;
	double intercept_stderr = 0.0;
	bool ok = false;
};

//! Ordinary least squares matching scipy.stats.linregress, including four
//! behaviors a textbook OLS gets wrong. Reverse-engineered from scipy 1.17.1
//! and pinned by data/syndna/fitb_*.csv:
//!
//!   - covariances are BIASED (divide by n, not n-1). Observable: every other
//!     use of ssxm is a ratio, but intercept_stderr multiplies by
//!     sqrt(ssxm + xmean^2), where it appears alone;
//!   - a TINY = 1e-20 term inside the t statistic keeps a perfect fit finite
//!     rather than dividing by zero (at r == 1 it yields a huge but real t);
//!   - n == 2 is special-cased: pvalue is exactly 1.0 or 0.0 and both standard
//!     errors are hard-set to 0.0;
//!   - intercept_stderr is stderr*sqrt(ssxm + xmean^2), NOT the textbook
//!     stderr*sqrt(1/n + xmean^2/Sxx).
//!
//! ONE deliberate numerical difference, in stderr. scipy writes it as
//! sqrt((1 - r^2)*ssym/ssxm/df); this accumulates the residual variance
//! directly. The two are algebraically identical -- mean((dy - slope*dx)^2) IS
//! (1 - r^2)*ssym -- and agree to ~1e-15 on well-conditioned data.
//!
//! They part company as r approaches 1, where 1 - r^2 is pure cancellation and
//! scipy's form loses a digit for every decade the residuals shrink. Measured
//! on y = 2x + 1 plus a shrinking wobble, against the true value the residual
//! form tracks exactly:
//!
//!     1 - r      scipy stderr            residual stderr        scipy is off by
//!     7.7e-05    0.010103717379028586    0.01010371737902853    5e-15
//!     7.7e-07    0.001010371737737878    0.0010103717379028527  1.6e-10
//!     7.7e-11    1.0103704639979164e-05  1.0103717379024465e-05 1.3e-06
//!     7.8e-15    1.017943078004628e-07   1.0103717376714997e-07 7.4e-03
//!     0          0.0                     1.0103717414876474e-08 everything
//!
//! So scipy is the inaccurate side for any fit cleaner than about r = 0.999999,
//! and a parity comparison against pysyndna on a very clean standard curve will
//! show a "mismatch" that is this port being right. Same situation as
//! StudentTSurvival at df = 1, and pinned the same way.
//!
//! At r == 1 exactly it is worse than inaccurate, it is unstable: whether r
//! lands on 1.0 is decided by summation order, and scipy reaches 0.0 only
//! because numpy's pairwise summation happens to win that coin flip. Straight
//! transcription gives 9.4e-9 against a golden of 0.0, and which side of the
//! flip we land on would vary with libm across x86, arm and WASM. The residual
//! form gives 8.5e-17 on every target.
//!
//! `x` and `y` must be the same length. n < 2, or constant x with n > 1, or any
//! NaN in the result yields `ok == false`.
LinregressResult Linregress(const std::vector<double> &x, const std::vector<double> &y);

//! True when a required per-sample parameter is usable: present, finite, and
//! not negative. NaN carries SQL NULL, and zero is deliberately allowed --
//! pysyndna tests with a strict `< 0` (util.py:262-264).
//!
//! Infinity is rejected even though pysyndna accepts it. It passes both of
//! pysyndna's tests and then makes every log10 in the sample infinite, whose
//! differences are NaN, so the sample is discarded anyway -- but as though a
//! fit had been attempted and failed rather than as a bad parameter.
//!
//! Shared rather than inlined because M3-M5 apply exactly this rule to
//! different parameter columns (calc_mass_sample_aliquot_input_g and friends),
//! and "NULL or negative means drop the sample" has to mean one thing across
//! all four functions.
bool IsUsableSampleParameter(double value);

//! Ids present in `subject` but absent from `reference`, deduplicated, sorted.
//!
//! The shape every absquant_* function needs: each takes a long-form counts
//! relation plus keyed reference relations, and must reject counts naming
//! something no reference can describe.
std::vector<std::string> IdsMissingFrom(const std::vector<std::string> &subject,
                                        const std::vector<std::string> &reference);

//! Ids occurring more than once in `ids`, deduplicated, sorted.
std::vector<std::string> DuplicatedIds(const std::vector<std::string> &ids);

//! Render ids for an error message: sorted, DEDUPLICATED, quoted,
//! comma-separated, and capped at ten with a ", ... (N more)" tail. Matches
//! NewickTree.cpp's `join_names`, which is the established rendering for "here
//! are the offending names" in this codebase -- a relation with ten thousand
//! bad ids must not produce a ten thousand id error string.
//!
//! (`join_names` is file-local to NewickTree.cpp, so this is a second copy
//! rather than a shared call. A third would be too many: whenever either file
//! is next touched for its own reasons, lift the pair into a shared header.)
std::string FormatIdList(std::vector<std::string> ids);

//! One long-form observation from the synDNA counts relation.
//!
//! Must be unique on (sample_id, feature_id); FitSyndnaModels rejects repeats.
//! pysyndna's input is a dense table where the cell is unique by construction,
//! so it has no equivalent of a repeat, and left unchecked a duplicate would
//! silently double-weight that synDNA in the regression.
//!
//! `count` must be finite and not negative. Zero is legal and means what it
//! says; the fit drops those points before taking the log.
struct CountObservation {
	std::string sample_id;
	std::string feature_id;
	double count = 0.0;
};

//! One row of the per-sample parameter relation. NaN carries SQL NULL.
//!
//! Must be unique on sample_id; FitSyndnaModels rejects repeats. Taking the
//! last row would not merely be a different tie-break from pysyndna: its
//! merge(how='left') fans the count rows out once per duplicate key, so the two
//! would disagree on the shape of the data rather than on one value.
struct SampleMass {
	std::string sample_id;
	double mass_syndna_input_ng = 0.0;
};

//! One row of the synDNA pool composition relation.
//!
//! Must be unique on feature_id, and syndna_indiv_ng_ul must be present,
//! finite, and not negative; FitSyndnaModels rejects all three. A duplicate
//! diverges from pysyndna twice over -- it fans out on the merge AND
//! double-counts inside the mass-fraction denominator. An infinite value is the
//! most destructive input this whole relation accepts if unchecked: it becomes
//! the entire denominator, collapsing every other synDNA's fraction to exactly
//! zero, so the run returns no models at all.
struct SyndnaConcentration {
	std::string feature_id;
	double syndna_indiv_ng_ul = 0.0;
};

struct FitOptions {
	//! Fraction of the input synDNA pool mass that actually contributed reads.
	double syndna_contributing_fraction = 1.0;
	//! A synDNA whose counts summed over all samples fall BELOW this is dropped.
	//! Strict `<`, so a synDNA sitting exactly on the threshold is kept.
	int64_t min_syndna_counts = 1;
};

struct SampleModel {
	std::string sample_id;
	LinregressResult fit;
};

//! Every sample present in `counts` lands in exactly one of `models`,
//! `filtered_sample_ids` and `unfittable_sample_ids`, so nothing disappears
//! without the caller being able to say why. The other two lists are on
//! different axes: `dropped_syndna_ids` is keyed by synDNA, and
//! `samples_without_counts` names samples that were never in `counts` at all.
struct FitResult {
	//! Samples that produced a usable model.
	std::vector<SampleModel> models;
	//! synDNAs whose total counts fell below `min_syndna_counts`.
	std::vector<std::string> dropped_syndna_ids;
	//! Samples removed for a NULL/NaN or negative mass_syndna_input_ng. A
	//! sample with no parameter row at all is an error, not a member here --
	//! "the value says to skip this one" and "the metadata is incomplete" are
	//! different situations and only the first is routine.
	std::vector<std::string> filtered_sample_ids;
	//! Samples that reached the fit but yielded no usable model -- no points
	//! left after filtering, scipy's constant-x rejection, or a NaN field.
	std::vector<std::string> unfittable_sample_ids;
	//! Samples named in `masses` with no rows in `counts`. Nothing is computed
	//! for them and they are not an error (pysyndna logs the same list), but
	//! silence would let a botched join look like a successful run.
	std::vector<std::string> samples_without_counts;
};

//! Fit log10(syndna_ng) ~ log10(read_count) per sample, reproducing pysyndna's
//! fit_linear_regression_models.
//!
//! Two orderings in pysyndna's pipeline are load-bearing and are preserved here
//! because the committed goldens distinguish them:
//!
//!   1. `min_syndna_counts` totals are summed across EVERY sample, including
//!      ones the parameter filter is about to remove. `fitb`'s b6 reaches
//!      exactly 20 reads only because the NULL-mass and negative-mass samples
//!      each contribute 2, and fitb_boundary_models_oracle.csv (generated at
//!      min_syndna_counts = 20) keeps b6 -- which is what lets three further
//!      samples fit at all.
//!   2. The mass-fraction denominator sums the FULL concentrations relation and
//!      is never rescaled after that drop. Summing over survivors instead is
//!      the obvious implementation and shifts every fitted intercept.
//!
//! Id-set mismatches are asymmetric, in both dimensions, and each asymmetry is
//! deliberate:
//!
//!   - A sample in `counts` with no `masses` row THROWS; a sample in `masses`
//!     with no counts is reported in `samples_without_counts`. pysyndna splits
//!     it the same way (util.py's validate_id_consistency_between_datasets
//!     raises on subset-not-in-superset and returns the other direction to be
//!     logged): reads you cannot convert to mass are unusable, metadata for a
//!     sample you did not sequence is merely untidy.
//!   - A synDNA in `counts` with no `concentrations` row THROWS; a synDNA in
//!     `concentrations` with no counts is treated as zero-count and so dropped
//!     by min_syndna_counts. pysyndna raises on BOTH, and this is a deliberate
//!     divergence: its read data is a dense table where a synDNA that sequenced
//!     nowhere is still present, whereas a long-form sparse relation cannot
//!     tell "zero everywhere" from "absent", so enforcing that direction would
//!     reject legitimate input. Because min_syndna_counts is required to be
//!     >= 1, a zero-count synDNA is dropped by pysyndna too -- the results are
//!     identical, only the error behavior differs.
//!
//! Throws std::invalid_argument on any of the above, on a malformed relation
//! (duplicate keys; a count that is negative or non-finite; a concentration
//! that is NULL, negative or non-finite; concentrations that do not sum to a
//! positive finite total) and on out-of-range options. The DuckDB wrapper
//! re-throws these as InvalidInputException prefixed with the function name,
//! per the BuildDenseDistanceMatrix / ReadDistanceTable precedent.
FitResult FitSyndnaModels(const std::vector<CountObservation> &counts,
                          const std::vector<SyndnaConcentration> &concentrations, const std::vector<SampleMass> &masses,
                          const FitOptions &options);

} // namespace absquant
} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers the absquant_* table functions into the extension catalog.
void RegisterAbsQuant(ExtensionLoader &loader);
} // namespace duckdb
