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
//! Infinity is rejected even though pysyndna accepts it, and what pysyndna
//! then does with it differs by function -- which is the argument for refusing
//! it rather than the argument against bothering:
//!
//!   - FitSyndnaModels: every log10 in the sample becomes infinite, their
//!     differences NaN, so the sample is discarded anyway -- but as though a fit
//!     had been attempted and failed rather than as a bad parameter;
//!   - ComputeOrfCopies, which takes no logarithm at all: an infinite
//!     DENOMINATOR (total_biological_reads_r1r2,
//!     calc_mass_sample_aliquot_input_g) collapses the sample to zeros, while an
//!     infinite NUMERATOR (total_rna_concentration_ng_ul,
//!     vol_extracted_elution_ul) writes +inf into every one of its cells.
//!
//! Three different wrong answers from one unusable input, and only one of them
//! even looks wrong. Refusing it up front is what names the problem.
//!
//! Shared rather than inlined because M2-M5 apply the same RULE to different
//! parameter columns (calc_mass_sample_aliquot_input_g and friends), and
//! "NULL, negative or infinite means drop the sample" has to mean one thing
//! across all four functions.
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

//! Throw if any (sample_id, feature_id) occurs more than once in `counts`.
//!
//! `relation` names the offending relation and is the ONLY thing that differed
//! between the three copies this replaces ("the synDNA counts relation" for the
//! fit, "the counts relation" for the other two); it is spliced in ahead of
//! " has more than one row for the same (sample_id, feature_id): ".
//!
//! Shared because all three absquant_* functions take a long-form counts
//! relation whose key SQL cannot enforce, and the check is a set, a loop and a
//! throw rather than a one-line message -- which is what makes a third copy
//! worth removing where the sibling id-mismatch messages, at two copies each,
//! are not. See FormatIdList's note for the rule.
void RejectDuplicateCells(const std::vector<CountObservation> &counts, const char *relation);

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

// ---------------------------------------------------------------------------
// Cell counts: applying a fitted model to turn OGU read counts into absolute
// quantities. Reproduces pysyndna's calc_ogu_cell_counts_biom.
// ---------------------------------------------------------------------------

//! Avogadro's number, in genomes per mole.
//!
//! pysyndna also carries a truncated 6.022e23 behind an `is_test` flag, to
//! reproduce the Zaramela notebook's own value. That path is deliberately NOT
//! ported (D2): it exists only to match a rounding in a published analysis, and
//! offering a knob that makes results *less* accurate invites someone to leave
//! it on.
constexpr double kAvogadro = 6.02214076e23;

//! Average molar mass of one base pair of double-stranded DNA, in g/mole.
//! pysyndna's `calc_copies_genomic_element_per_g_series` default.
constexpr double kGramsPerMolePerBasePair = 650.0;

//! The unit the predicted cell counts are expressed per.
//!
//! `CellsPerGOfGdna` normalizes by the gDNA that went into the sequencer, which
//! is an instrument quantity. The other three normalize by the material that
//! was actually collected, which is what makes them comparable across samples,
//! and each divides by its own column of the parameters relation.
enum class CellCountsMetric { CellsPerGOfGdna, CellsPerGOfSample, CellsPerUlOfSample, CellsPerCm2OfSample };

//! Every metric, in the order the user-facing documentation lists them.
//!
//! The ONE place the set is enumerated. MetricName and DenominatorColumnName
//! are exhaustive switches with no `default:`, but this build does not enable
//! -Wswitch (CMakeLists.txt has -Wall/-Wextra commented out for the extension
//! target), so a fifth enumerator would not be diagnosed there either. Deriving
//! the parser and the accepted-value list from this array is what actually
//! makes a half-added metric impossible, rather than restating the set in each
//! of them and hoping they stay in step.
inline constexpr CellCountsMetric kAllCellCountsMetrics[] = {
    CellCountsMetric::CellsPerGOfGdna, CellCountsMetric::CellsPerGOfSample, CellCountsMetric::CellsPerUlOfSample,
    CellCountsMetric::CellsPerCm2OfSample};

//! The SQL-facing spelling of a metric, e.g. "cells_per_g_of_gdna".
const char *MetricName(CellCountsMetric metric);

//! Parses one of those spellings, case-insensitively. Returns false and leaves
//! `out` untouched for anything else.
bool ParseCellCountsMetric(const std::string &name, CellCountsMetric &out);

//! The parameters-relation column a metric divides by, or nullptr for
//! `CellsPerGOfGdna`, which has no sample-side denominator.
//!
//! This is the single source of truth for that mapping: the DuckDB wrapper
//! calls it to decide which column to read, and the core calls it to name the
//! column in its diagnostics, so the column read and the column blamed cannot
//! drift apart.
const char *DenominatorColumnName(CellCountsMetric metric);

//! One row of the per-(sample, feature) coverage relation.
//!
//! Coverage is a FRACTION in [0, 1], never a percent (D9). pysyndna accepts
//! either and puts the burden on the caller to pass a matching `min_coverage`,
//! which makes handing it percents against a fractional threshold a silent
//! no-op that keeps every feature. miint rejects anything outside [0, 1], so
//! that mistake fails loudly instead.
//!
//! Must be unique on (sample_id, feature_id).
struct CoverageObservation {
	std::string sample_id;
	std::string feature_id;
	double coverage = 0.0;
};

//! One row of the feature-length relation. `ogu_len_in_bp` is a denominator, so
//! it must be finite and strictly positive; must be unique on feature_id.
struct FeatureLength {
	std::string feature_id;
	double ogu_len_in_bp = 0.0;
};

//! A fitted per-sample model, as produced by absquant_fit_models. Only these
//! three of its six fields are used: slope and intercept predict the mass, and
//! rvalue feeds the r^2 gate. Must be unique on sample_id.
struct SampleRegression {
	std::string sample_id;
	double slope = 0.0;
	double intercept = 0.0;
	double rvalue = 0.0;
};

//! Per-sample parameters for the cell-count calculation. Must be unique on
//! sample_id.
//!
//! Only `sequenced_sample_gdna_mass_ng` enters the cells_per_g_of_gdna
//! arithmetic. The other two are pysyndna's REQUIRED_DNA_PREP_INFO_KEYS
//! (calc_cell_counts.py:53): they compute extracted_gdna_mass_g, which only the
//! sample-level metrics need, yet they are required and NaN/negative-filtered
//! for every metric. Keeping that is parity -- dropping the requirement would
//! silently change which samples appear in the output.
struct SampleCellParams {
	std::string sample_id;
	double sequenced_sample_gdna_mass_ng = 0.0;
	double extracted_gdna_concentration_ng_ul = 0.0;
	double vol_extracted_elution_ul = 0.0;
	//! Whichever column DenominatorColumnName picked for the requested metric,
	//! already resolved by the caller. ONE field rather than one per metric,
	//! because only one is ever read: three named fields would leave two
	//! permanently unset, which a later reader would reasonably take for data.
	//!
	//! Unused, and left at zero, for cells_per_g_of_gdna -- which is also why it
	//! is screened for the other three metrics only. pysyndna filters on the
	//! REQUESTED metric's denominator and no other
	//! (calc_ogu_cell_counts_biom:959-965), so a sample missing a column the
	//! query never reads stays in.
	double sample_denominator = 0.0;
};

struct CellCountsOptions {
	//! Which unit to express the counts per. Selects `sample_denominator`'s
	//! meaning; see CellCountsMetric.
	CellCountsMetric metric = CellCountsMetric::CellsPerGOfGdna;
	//! A (sample, feature) whose coverage is BELOW this is dropped. Strict `<`,
	//! so a feature sitting exactly on the threshold is kept. A fraction in
	//! [0, 1].
	double min_coverage = 0.0;
	//! A sample whose model has r^2 BELOW this is dropped entirely. Strict `<`.
	double min_rsquared = 0.8;
};

//! One cell of a long-form output feature table.
//!
//! Shared by ComputeCellCounts and ComputeOrfCopies rather than spelled twice:
//! both return the same `(sample_id, feature_id, value)` sparse shape, and both
//! omit zero-valued cells under the same invariant (D10). Two identical structs
//! would let that shape drift apart between the two functions for no reason.
struct FeatureTableValue {
	std::string sample_id;
	std::string feature_id;
	double value = 0.0;
};

//! Every sample present in `counts` either contributes cells to `values` or
//! appears in exactly one of the five sample lists below, so no sample
//! disappears without the caller being able to say why. (Input that would break
//! that -- a counted sample with no parameter row, a counted cell with no
//! coverage row -- is rejected before any of this runs.)
//! `low_coverage_feature_ids` is on the other axis, keyed by feature rather than
//! by sample.
struct CellCountsResult {
	//! Surviving (sample, feature) cells. Zero-valued cells are omitted, per
	//! miint's long-form sparse invariant (D10); pysyndna's dense output spells
	//! the same thing as an explicit 0.0 after its final fillna.
	std::vector<FeatureTableValue> values;
	//! Features dropped somewhere for coverage below min_coverage. Keyed by
	//! feature, not by (sample, feature): pysyndna reports the same list this
	//! way, and a feature that fails in one sample is worth naming once.
	//!
	//! Sorted, where pysyndna's is in the order its melt happens to produce.
	//! That is presentation only -- the contents are identical -- and being
	//! independent of input row order is what makes it assertable.
	std::vector<std::string> low_coverage_feature_ids;
	//! Samples removed because a required parameter was unusable: NULL/NaN,
	//! negative, or infinite. Same rule and same name as FitResult's list.
	std::vector<std::string> filtered_sample_ids;
	//! Samples with counts but no USABLE model -- either no row in `models` at
	//! all, or one whose slope, intercept or rvalue is NULL/NaN or infinite.
	//! One bucket because they are one thing to pysyndna, which records None for
	//! a fit with any NaN field (_convert_linregressresults_to_dict:329-334) and
	//! reaches the same log line either way.
	//!
	//! A warning rather than an error, matching pysyndna (calc_cell_counts.py
	//! :674) and matching M2's treatment of a sample that could not be fit:
	//! absquant_fit_models legitimately returns fewer models than it was given
	//! samples, so feeding its output straight back in must not be an error.
	//! (An rvalue that is a number but outside [-1, 1] is a different matter --
	//! that is a malformed relation and it throws.)
	std::vector<std::string> samples_without_models;
	//! Samples whose model's r^2 fell below min_rsquared.
	std::vector<std::string> low_rsquared_sample_ids;
	//! Samples every one of whose cells fell below min_coverage, so nothing was
	//! left to compute. pysyndna emits NO message for these at all -- it loops
	//! over the sample ids remaining in the working table
	//! (calc_cell_counts.py:505), which a wholly uncovered sample has already
	//! left -- so the sample vanishes from its output unremarked. The values
	//! agree; miint just says so. Naming the features is not a substitute,
	//! because in a real run that list is capped at ten ids.
	std::vector<std::string> uncovered_sample_ids;
	//! Samples that passed every filter and still produced nothing, because
	//! every one of their cells came out exactly zero and the sparse invariant
	//! omits those. In practice that means a zero extraction ratio: either
	//! extracted_gdna_concentration_ng_ul or vol_extracted_elution_ul was zero,
	//! which pysyndna's `< 0` parameter screen admits and a blank extraction can
	//! genuinely produce.
	//!
	//! Not reachable for cells_per_g_of_gdna through any realistic model, whose
	//! values are strictly positive, which is why M3 had no such list. It is not
	//! *impossible* there: a hand-built model with an extreme negative intercept
	//! underflows 10^x to exactly zero, and this list would then name the
	//! sample. That is the honest boundary -- the metric does not make the case
	//! disappear, it makes it unreachable from absquant_fit_models' output.
	//!
	//! Reported only when the WHOLE sample goes this way. A single zero cell
	//! among others needs no diagnostic, because under D10 an omitted cell and a
	//! dense 0.0 are the same claim and the sample is still there to be read.
	//! It is the all-zero sample that the sparse form cannot express: "no rows
	//! for this sample" would otherwise be indistinguishable from "no such
	//! sample". pysyndna never faces this, its output being dense.
	std::vector<std::string> zero_valued_sample_ids;
};

//! Predict cells of each feature per unit of sample or of gDNA, reproducing
//! pysyndna's calc_ogu_cell_counts_biom.
//!
//! Per surviving (sample, feature):
//!
//!     ogu_gdna_mass_ng   = 10 ^ (slope * log10(read_count) + intercept)
//!     ogu_genomes        = ogu_gdna_mass_ng * kAvogadro
//!                          / (ogu_len_in_bp * kGramsPerMolePerBasePair * 1e9)
//!     genomes_per_g_gdna = ogu_genomes / sequenced_sample_gdna_mass_ng * 1e9
//!     cells_per_g_gdna   = genomes_per_g_gdna
//!
//! The last line is an explicit simplifying assumption -- one genome per cell --
//! which pysyndna's own comment acknowledges is wrong for dividing and
//! polyploid microbes. It is reproduced because it is the published method.
//!
//! The three sample-level metrics then scale that by one per-sample ratio:
//!
//!     extracted_gdna_mass_g = extracted_gdna_concentration_ng_ul
//!                             * vol_extracted_elution_ul / 1e9
//!     ratio                 = extracted_gdna_mass_g / sample_denominator
//!     cells_per_<unit>      = cells_per_g_gdna * ratio
//!
//! `extracted_gdna_mass_g` is the gDNA recovered from the whole extraction, not
//! the portion that was sequenced -- which is exactly what turns "per gram of
//! sequenced gDNA" into "per unit of collected material". The ratio is 1 for
//! cells_per_g_of_gdna, whose values are unaffected.
//!
//! Three filters run before that arithmetic, and their ORDER is load-bearing --
//! not for the values, which are the same whichever way round it goes, but for
//! the diagnostics:
//!
//!   1. samples with an unusable required parameter  -> filtered_sample_ids
//!   2. cells with coverage < min_coverage           -> low_coverage_feature_ids
//!   3. samples with no model, or with r^2 < min     -> samples_without_models,
//!                                                      low_rsquared_sample_ids
//!
//! Filter 1 screens sequenced_sample_gdna_mass_ng, both extraction parameters,
//! and -- for the three sample-level metrics only -- sample_denominator. That
//! last one is per metric, matching pysyndna: screening all three sample columns
//! or none of them reproduces every VALUE and still gets the sample membership
//! wrong.
//!
//! pysyndna removes the bad-parameter samples from the counts table itself
//! (calc_cell_counts.py:1013) before the coverage filter ever runs, so such a
//! sample never contributes its poorly covered features to the coverage list.
//! Running the coverage filter first would name features that no surviving
//! sample was ever going to use, sending the user after a data problem that
//! does not exist.
//!
//! Read counts reach here already positive: the DuckDB reader drops zero-valued
//! cells, and a count that is not finite and positive is refused. That matters
//! more than it looks. log10(0) is -inf, and with a NEGATIVE slope the mass
//! would come back as 10^(+inf) = inf, i.e. an infinite cell count -- which is
//! exactly what pysyndna emits for a zero count under a negative model. The
//! check is unreachable through the reader, and is what keeps it that way.
//!
//! Id-set enforcement is asymmetric, deliberately. The counts relation is the
//! subject: every cell it names must have a coverage row, its feature must have
//! a length and its sample must have parameters, or nothing defined can be said
//! about that cell. The reverse never has to hold -- a lengths table is a whole
//! reference database and a parameters table a whole study -- so the reference
//! relations may describe far more than the counts do, and their unused rows
//! are not screened at all. The one id relationship that is NOT enforced is
//! `models`: a sample with no usable model is reported, not refused, so that
//! absquant_fit_models' output composes with this function directly.
//!
//! Throws std::invalid_argument on out-of-range options (min_coverage or
//! min_rsquared outside [0, 1]), on a malformed relation (duplicate keys; a
//! count that is not finite and positive; a coverage outside [0, 1]; a
//! non-positive or non-finite ogu_len_in_bp; a zero sequenced_sample_gdna_mass_ng;
//! a zero sample_denominator when the metric divides by one; an rvalue outside
//! [-1, 1]), on the id mismatches above, and on a computed cell count that
//! overflowed to a non-finite value. The DuckDB wrapper re-throws these as
//! InvalidInputException prefixed with the function name.
//!
//! Note the split the two denominators share: zero is an ERROR, because dividing
//! by it is structurally impossible and pysyndna silently emits inf cells for
//! every feature in the sample, while NaN, negative and infinite values are
//! FILTERED with a diagnostic, matching pysyndna's own screen.
CellCountsResult ComputeCellCounts(const std::vector<CountObservation> &counts,
                                   const std::vector<SampleRegression> &models,
                                   const std::vector<CoverageObservation> &coverage,
                                   const std::vector<FeatureLength> &lengths,
                                   const std::vector<SampleCellParams> &params, const CellCountsOptions &options);

// ---------------------------------------------------------------------------
// ORF copies: turning metatranscriptomic OGU+ORF read counts into copies of
// each ORF's ssRNA per gram of sample. Reproduces pysyndna's
// calc_copies_of_ogu_orf_ssrna_per_g_sample.
// ---------------------------------------------------------------------------

//! Average molar mass of one base of single-stranded RNA, in g/mole.
//! pysyndna's RNA_BASE_G_PER_MOLE (util.py:9).
//!
//! Beside kGramsPerMolePerBasePair rather than instead of it: 650 is a base
//! PAIR of double-stranded DNA and 340 is a single RNA base, so which constant
//! applies is decided by the molecule, not by the function.
constexpr double kGramsPerMolePerRnaBase = 340.0;

//! One row of the ORF coordinate relation. Must be unique on feature_id.
//!
//! Both coordinates are 1-based and CLOSED, the convention woltka's coords.txt
//! uses, so the ORF spans |end - start| + 1 bases. `start > end` is legal and
//! marks a reverse-strand ORF -- pysyndna says so explicitly
//! (quant_orfs.py:248, "it is NOT required that start be less than end, per
//! woltka docs") and four of the ten features in data/syndna/orf_coords.csv
//! are that way round.
//!
//! Doubles rather than integers because that is the shape a SQL DOUBLE column
//! and a NULL (arriving as NaN) have. A coordinate that is not finite and
//! integral is refused rather than truncated: pysyndna's cast_cols would take
//! 100.5 and produce a fractional length, and a fractional genome coordinate is
//! not a measurement.
struct OrfCoords {
	std::string feature_id;
	double ogu_orf_start = 0.0;
	double ogu_orf_end = 0.0;
};

//! Per-sample parameters for the ORF copy calculation. Must be unique on
//! sample_id.
//!
//! All four are required, all four are screened, and all four are actually
//! divided by or multiplied through -- unlike the cell-count parameters, none
//! of these is carried only for parity. pysyndna screens exactly this set:
//! REQUIRED_PARAM_KEYS is the union of its sample-info and RNA-prep key lists
//! (quant_orfs.py:10-22) and goes to filter_data_by_sample_info whole.
struct SampleOrfParams {
	std::string sample_id;
	double calc_mass_sample_aliquot_input_g = 0.0;
	double total_rna_concentration_ng_ul = 0.0;
	double vol_extracted_elution_ul = 0.0;
	double total_biological_reads_r1r2 = 0.0;
};

//! Every sample present in `counts` either contributes cells to `values` or is
//! named in exactly one of the two lists below, so no sample disappears without
//! the caller being able to say why.
struct OrfCopiesResult {
	//! Surviving (sample, feature) cells, in the same sparse form and under the
	//! same invariant as CellCountsResult::values.
	std::vector<FeatureTableValue> values;
	//! Samples removed because a required parameter was unusable: NULL/NaN,
	//! negative, or infinite. Same rule and same name as FitResult's list.
	std::vector<std::string> filtered_sample_ids;
	//! Samples that passed the filter and still produced nothing, because every
	//! one of their cells came out exactly zero and the sparse invariant omits
	//! those. In practice a zero extraction: either
	//! total_rna_concentration_ng_ul or vol_extracted_elution_ul was zero, which
	//! pysyndna's `< 0` screen admits and a blank really can produce. Underflow
	//! reaches it too -- their product is formed before the ng->g divide, so two
	//! very small values can round it to exactly 0.0.
	//!
	//! Reported only when the WHOLE sample goes this way, for the same reason
	//! CellCountsResult gives: a single zero cell among others needs no
	//! diagnostic, because under D10 an omitted cell and a dense 0.0 are the
	//! same claim. It is the all-zero sample the sparse form cannot express --
	//! "no rows for this sample" would otherwise be indistinguishable from "no
	//! such sample". pysyndna never faces this, its output being dense.
	std::vector<std::string> zero_valued_sample_ids;
};

//! Copies of each ORF's ssRNA per gram of sample, reproducing pysyndna's
//! calc_copies_of_ogu_orf_ssrna_per_g_sample.
//!
//! Per (sample, feature):
//!
//!     ogu_orf_len        = |ogu_orf_end - ogu_orf_start| + 1
//!     copies_per_g_ssrna = kAvogadro / (ogu_orf_len * kGramsPerMolePerRnaBase)
//!     fraction_of_reads  = count / total_biological_reads_r1r2
//!     g_total_ssrna      = total_rna_concentration_ng_ul
//!                          * vol_extracted_elution_ul / 1e9
//!     value              = fraction_of_reads * g_total_ssrna
//!                          * copies_per_g_ssrna
//!                          / calc_mass_sample_aliquot_input_g
//!
//! The quantity is copies of the ORF's OWN ssRNA, not of the transcript
//! containing it -- a transcript may carry further ORFs and be heavier.
//! pysyndna's docstring is explicit about this and the distinction is the
//! difference between a per-gene and a per-message abundance.
//!
//! No standard curve, no coverage filter and no r^2 gate: this workflow does
//! not go through a spike-in at all. The read fraction, the extracted ssRNA
//! mass and the ORF's length are enough on their own, which is why there is no
//! options struct and no metric to choose -- pysyndna exposes exactly one.
//!
//! Read counts reach here already positive: the DuckDB reader drops
//! zero-valued cells. A zero count is accepted regardless, unlike in
//! ComputeCellCounts -- the difference is the log10. There is none here, so a
//! zero count simply means no copies of that ORF, and refusing it would be a
//! rule invented for no reason.
//!
//! One filter runs before the arithmetic: samples with an unusable value in any
//! of the four parameter columns are dropped into filtered_sample_ids.
//!
//! Id-set enforcement is asymmetric, deliberately. The counts relation is the
//! subject: every cell it names must have coordinates and its sample must have
//! parameters, or nothing defined can be said about that cell. The reverse
//! never has to hold -- a coords relation is a whole annotation and a params
//! relation a whole study -- so both may describe far more than the counts do,
//! and their unused rows are not screened at all. pysyndna splits the sample
//! axis the same way (quant_orfs.py:291, whose comment says extra parameter
//! samples "could just be samples that failed sequencing/etc.") and does not
//! check the feature axis at all: it reaches the coordinates through a bare
//! `.at[]` inside a biom transform, so a counted ORF the annotation does not
//! describe surfaces as a pandas KeyError from inside a lambda.
//!
//! Throws std::invalid_argument on a malformed relation (duplicate keys; a
//! count that is negative or non-finite; an ogu_orf_start or ogu_orf_end that
//! is not a finite whole number; a zero total_biological_reads_r1r2 or
//! calc_mass_sample_aliquot_input_g), on the id mismatches above, and on a
//! computed copy count that overflowed to a non-finite value. The DuckDB
//! wrapper re-throws these as InvalidInputException prefixed with the function
//! name.
//!
//! Note the split the two denominators share, the same one ComputeCellCounts
//! draws: zero is an ERROR, because dividing by it is structurally impossible
//! and pysyndna silently emits inf copies for every ORF in the sample, while
//! NaN, negative and infinite values are FILTERED with a diagnostic, matching
//! pysyndna's own screen. The other two parameter columns are multiplied
//! through rather than divided by, so zero there is ordinary data -- a blank
//! extraction genuinely yields no ssRNA -- and lands the sample in
//! zero_valued_sample_ids instead.
OrfCopiesResult ComputeOrfCopies(const std::vector<CountObservation> &counts, const std::vector<OrfCoords> &coords,
                                 const std::vector<SampleOrfParams> &params);

} // namespace absquant
} // namespace miint

namespace duckdb {
class ExtensionLoader;
//! Registers absquant_fit_models into the extension catalog.
void RegisterAbsQuant(ExtensionLoader &loader);
//! Registers absquant_cell_counts into the extension catalog. Separate from
//! RegisterAbsQuant, and separately sourced, because the two functions are ~350
//! lines of bind/read apiece and the established convention is one Register* per
//! function called from LoadInternal.
void RegisterAbsQuantCellCounts(ExtensionLoader &loader);
//! Registers absquant_orf_copies into the extension catalog, on the same
//! one-Register*-per-function convention.
void RegisterAbsQuantOrfCopies(ExtensionLoader &loader);
} // namespace duckdb
