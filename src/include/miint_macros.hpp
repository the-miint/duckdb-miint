#pragma once

#include "duckdb.hpp"
#include "duckdb/main/extension_helper.hpp"
#include "duckdb/parser/parser.hpp"
#include "duckdb/parser/statement/create_statement.hpp"
#include "duckdb/parser/parsed_data/create_macro_info.hpp"

namespace duckdb {

// woltka_ogu is now a C++ table function — see src/woltka_ogu_function.cpp

const std::string MIINT_WARNINGS = // NOLINT
    "CREATE OR REPLACE MACRO miint_warnings() AS TABLE "
    "SELECT timestamp, message "
    "FROM duckdb_logs() "
    "WHERE type = 'MiintWarning' "
    "ORDER BY timestamp;";

// GFF3 (col 9) reserves ; = & , and requires them percent-encoded inside keys and
// values, so decode every key/value with url_decode. The value is everything after
// the FIRST '=' (rejoined), not string_split(x,'=')[2] — a raw '=' in a value (e.g.
// an insert name "gc=0.46") must not be truncated, which would silently collapse
// distinct features onto one id. A token with no '=' keeps a NULL value.
const std::string PARSE_GFF_ATTRIBUTES = // NOLINT
    "CREATE OR REPLACE MACRO parse_gff_attributes(kvp_string) AS ( "
    "  map_from_entries( "
    "    list_transform( "
    "      string_split(kvp_string, ';'), "
    "      lambda x: struct_pack( "
    "        key := url_decode(string_split(x, '=')[1]), "
    "        value := CASE WHEN contains(x, '=') "
    "                      THEN url_decode(array_to_string(string_split(x, '=')[2:], '=')) "
    "                      ELSE NULL END "
    "      ) "
    "    ) "
    "  ) "
    "); ";

const std::string READ_GFF = // NOLINT
    "CREATE OR REPLACE MACRO read_gff(path) AS TABLE "
    "SELECT "
    "   column0 AS seqid, "
    "   column1 AS source, "
    "   column2 AS type, "
    "   column3::INTEGER AS position, "
    // GFF3's `end` is 1-based CLOSED; miint's project-wide convention is 1-based HALF-OPEN
    // [position, stop_position), as used by read_alignments, align_*, alignment_slice,
    // compute_coverage_depth and genome_coverage. Normalize on read (+1) so `stop_position`
    // means the same thing everywhere. Before this, read_gff was the outlier while sharing
    // the column NAME, so composing it with any of those silently dropped the interval's
    // last base -- genome_coverage's SUM(stop - start) under-counted every feature by one.
    // The raw GFF `end` is recoverable as `stop_position - 1` (issue #196).
    "   column4::INTEGER + 1 AS stop_position, "
    "   CASE  "
    "     WHEN column5 = '.' THEN NULL  "
    "     ELSE column5::DOUBLE  "
    "   END AS score, "
    "   CASE  "
    "     WHEN column6 = '.' THEN NULL  "
    "     ELSE column6  "
    "   END AS strand, "
    "   CASE  "
    "     WHEN column7 = '.' THEN NULL  "
    "     ELSE column7::INTEGER  "
    "   END AS phase, "
    "   parse_gff_attributes(column8) AS attributes "
    "FROM read_csv(path, "
    " delim = '\t', "
    "   header = false, "
    "   columns = { "
    "     'column0': 'VARCHAR', "
    "     'column1': 'VARCHAR', "
    "     'column2': 'VARCHAR', "
    "     'column3': 'VARCHAR', "
    "     'column4': 'VARCHAR', "
    "     'column5': 'VARCHAR', "
    "     'column6': 'VARCHAR', "
    "     'column7': 'VARCHAR', "
    "     'column8': 'VARCHAR' "
    "   }, auto_detect=false, "
    "   comment = '#', "
    "   skip = 0, "
    "   null_padding = true "
    " ) "
    "WHERE column0 NOT LIKE '##%' "
    // The GFF3 sequence section (issue #186). Everything after a ##FASTA directive is
    // FASTA, and prokka/bakta ALWAYS append the genome -- so this is the common case. Those
    // lines contain no tabs, so with null_padding every mandatory field after column0
    // parses as NULL. Dropped silently: appended FASTA is expected file content.
    //
    // Some-but-not-all mandatory fields present is a different thing -- a broken feature
    // line -- and is raised rather than dropped, so we don't trade one silent wrong answer
    // for another. column8 (attributes) is deliberately NOT required: a legitimate line may
    // end with an empty attributes field, which also parses as NULL, so keying off it would
    // reject real prokka output.
    //
    // This MUST live in the WHERE clause, not a projection: DuckDB prunes unused
    // projections, so a check inside a SELECT expression would not fire for
    // `SELECT type FROM read_gff(...)` -- the very shape the issue reported.
    "  AND CASE "
    // A FASTA header line, checked FIRST because its description may itself contain tabs
    // (e.g. ">ctg1<TAB>some description"), which would otherwise populate column1 and trip
    // the malformed-feature branch below on a perfectly valid file. '>' is not a legal GFF3
    // seqid character, so this prefix identifies FASTA unambiguously.
    "        WHEN column0 LIKE '>%' THEN FALSE "
    "        WHEN column1 IS NULL AND column2 IS NULL AND column3 IS NULL "
    "         AND column4 IS NULL AND column5 IS NULL AND column6 IS NULL "
    "         AND column7 IS NULL THEN FALSE "
    "        WHEN column1 IS NULL OR column2 IS NULL OR column3 IS NULL "
    "          OR column4 IS NULL OR column5 IS NULL OR column6 IS NULL "
    "          OR column7 IS NULL "
    "        THEN error(printf('read_gff: malformed GFF feature line "
    "(expected 9 tab-separated fields, found fewer): %s', column0)) "
    "        ELSE TRUE "
    "      END; ";

// read_jplace(path)
//
// Read jplace phylogenetic placement file(s) and return the best placement
// for each fragment. Supports glob patterns for multiple files.
// The jplace format is defined in:
// Matsen FA, Hoffman NG, Gallagher A, Stamatakis A (2012) A Format for
// Phylogenetic Placements. PLoS ONE 7(2): e31009.
//
// Parameters:
// path : VARCHAR, path to a jplace file (supports glob patterns)
//
// Returns a table with columns:
// - fragment: VARCHAR (fragment/sequence name)
// - edge_num: INTEGER (edge number in the reference tree)
// - likelihood: DOUBLE (log likelihood of placement)
// - like_weight_ratio: DOUBLE (likelihood weight ratio)
// - distal_length: DOUBLE (distance from distal end of edge)
// - pendant_length: DOUBLE (pendant branch length)
// - filepath: VARCHAR (source file path)
//
// Note: Only the best placement (first in 'p' array) is returned.
// Supports both 'nm' (named multiplicities) and 'n' (names) formats.
const std::string READ_JPLACE = // NOLINT
    "CREATE OR REPLACE MACRO read_jplace(path) AS TABLE "
    "SELECT "
    "    COALESCE( "
    "        json_extract_string(placement, '$.nm[0][0]'), "
    "        json_extract_string(placement, '$.n[0]') "
    "    ) AS fragment, "
    "    json_extract(placement, '$.p[0][0]')::INTEGER AS edge_num, "
    "    json_extract(placement, '$.p[0][1]')::DOUBLE AS likelihood, "
    "    json_extract(placement, '$.p[0][2]')::DOUBLE AS like_weight_ratio, "
    "    json_extract(placement, '$.p[0][3]')::DOUBLE AS distal_length, "
    "    json_extract(placement, '$.p[0][4]')::DOUBLE AS pendant_length, "
    "    filepath "
    "FROM ( "
    "    SELECT unnest(placements) AS placement, filename AS filepath "
    "    FROM read_json(path, filename := true, maximum_object_size=1000000000) "
    "); ";

// mz_within(observed, target, tolerance_da)
//
// Returns true if abs(observed - target) < tolerance_da (strict inequality).
// Matches MassQL TOLERANCEMZ semantics.
const std::string MZ_WITHIN = // NOLINT
    "CREATE OR REPLACE MACRO mz_within(observed, target, tolerance_da) AS "
    "(observed > target - tolerance_da AND observed < target + tolerance_da);";

// mz_within_ppm(observed, target, ppm)
//
// Returns true if abs(observed - target) < target * ppm * 1e-6 (strict inequality).
// Matches MassQL TOLERANCEPPM semantics.
const std::string MZ_WITHIN_PPM = // NOLINT
    "CREATE OR REPLACE MACRO mz_within_ppm(observed, target, ppm) AS "
    "(observed > target * (1.0 - ppm * 1e-6) AND observed < target * (1.0 + ppm * 1e-6));";

// mzml_peaks(relation)
//
// Unnests mz_array and intensity_array into per-peak rows.
// Takes any relation containing read_mzml output (table, view, subquery).
// Each output row is one (mz, intensity) peak with its parent spectrum's metadata.
// i_norm = intensity / base_peak_intensity (NULL when base_peak_intensity is NULL).
// i_tic = intensity / total_ion_current (NULL when TIC is 0 or NULL — division-by-zero safe).
// neutral_loss = precursor_mz - mz (NULL for MS1 spectra).
const std::string MZML_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_peaks(relation) AS TABLE "
    "SELECT spectrum_index, spectrum_id, scan_number, ms_level, retention_time, spectrum_type, polarity, "
    "       base_peak_intensity, total_ion_current, "
    "       precursor_mz, precursor_charge, precursor_intensity, ms1_scan_index, "
    "       mz, intensity, intensity / base_peak_intensity AS i_norm, "
    "       CASE WHEN total_ion_current > 0 THEN intensity / total_ion_current ELSE NULL END AS i_tic, "
    "       precursor_mz - mz AS neutral_loss "
    "FROM ( "
    "    SELECT spectrum_index, spectrum_id, scan_number, ms_level, retention_time, spectrum_type, polarity, "
    "           base_peak_intensity, total_ion_current, "
    "           precursor_mz, precursor_charge, precursor_intensity, ms1_scan_index, "
    "           UNNEST(mz_array) AS mz, UNNEST(intensity_array) AS intensity "
    "    FROM query_table(relation) "
    "); ";

// mzml_scaninfo(relation)
//
// Re-aggregates peak-level data back to one-row-per-scan with summary statistics.
// Takes any relation (typically filtered output of mzml_peaks).
// Uses query_table(relation) for the relation parameter.
const std::string MZML_SCANINFO = // NOLINT
    "CREATE OR REPLACE MACRO mzml_scaninfo(relation) AS TABLE "
    "SELECT "
    "    spectrum_index, "
    "    first(scan_number) AS scan_number, "
    "    first(ms_level) AS ms_level, "
    "    first(retention_time) AS retention_time, "
    "    first(spectrum_type) AS spectrum_type, "
    "    first(polarity) AS polarity, "
    "    SUM(intensity) AS total_ion_current, "
    "    MAX(intensity) AS base_peak_intensity, "
    "    MAX(i_norm) AS i_norm, "
    "    first(precursor_mz) AS precursor_mz, "
    "    first(precursor_charge) AS precursor_charge, "
    "    first(precursor_intensity) AS precursor_intensity, "
    "    first(ms1_scan_index) AS ms1_scan_index "
    "FROM query_table(relation) "
    "GROUP BY spectrum_index; ";

// mzml_scansum(relation)
//
// Aggregate peak-level data to one row per scan with total intensity.
const std::string MZML_SCANSUM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_scansum(relation) AS TABLE "
    "SELECT spectrum_index, SUM(intensity) AS total_intensity "
    "FROM query_table(relation) "
    "GROUP BY spectrum_index; ";

// mzml_scannum(relation)
//
// Return distinct spectrum indices from peak-level data.
const std::string MZML_SCANNUM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_scannum(relation) AS TABLE "
    "SELECT DISTINCT spectrum_index "
    "FROM query_table(relation); ";

// mzml_scanmz(relation)
//
// Return distinct non-NULL precursor_mz values from peak-level data.
const std::string MZML_SCANMZ = // NOLINT
    "CREATE OR REPLACE MACRO mzml_scanmz(relation) AS TABLE "
    "SELECT DISTINCT precursor_mz "
    "FROM query_table(relation) "
    "WHERE precursor_mz IS NOT NULL; ";

// mzml_scanmaxint(relation)
//
// Aggregate peak-level data to one row per scan with max intensity.
const std::string MZML_SCANMAXINT = // NOLINT
    "CREATE OR REPLACE MACRO mzml_scanmaxint(relation) AS TABLE "
    "SELECT spectrum_index, MAX(intensity) AS max_intensity "
    "FROM query_table(relation) "
    "GROUP BY spectrum_index; ";

// mzml_ms1_peaks(relation)
//
// Convenience wrapper: returns only MS1 peaks from mzml_peaks.
const std::string MZML_MS1_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms1_peaks(relation) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1; ";

// mzml_ms2_peaks(relation)
//
// Convenience wrapper: returns only MS2 peaks from mzml_peaks.
const std::string MZML_MS2_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms2_peaks(relation) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2; ";

// mzml_ms1_parent_peaks(spectra_rel, ms2_peaks_rel)
//
// Given filtered MS2 peaks, return MS1 peaks from their parent scans.
// Uses ms1_scan_index to link MS2 → MS1. Orphan MS2 (NULL ms1_scan_index) are ignored.
const std::string MZML_MS1_PARENT_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms1_parent_peaks(spectra_rel, ms2_peaks_rel) AS TABLE "
    "SELECT * FROM mzml_peaks(spectra_rel) WHERE ms_level = 1 "
    "AND spectrum_index IN ("
    "  SELECT DISTINCT ms1_scan_index FROM query_table(ms2_peaks_rel) "
    "  WHERE ms1_scan_index IS NOT NULL"
    "); ";

// mzml_ms2_child_peaks(spectra_rel, ms1_peaks_rel)
//
// Given filtered MS1 peaks, return MS2 peaks from child scans.
// Uses ms1_scan_index to link MS1 → MS2.
const std::string MZML_MS2_CHILD_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms2_child_peaks(spectra_rel, ms1_peaks_rel) AS TABLE "
    "SELECT * FROM mzml_peaks(spectra_rel) WHERE ms_level = 2 "
    "AND ms1_scan_index IN ("
    "  SELECT DISTINCT spectrum_index FROM query_table(ms1_peaks_rel)"
    "); ";

// mzml_ms1_where_ms2prod(relation, target_mz [, tolerance])
//
// MassQL: QUERY scansum(MS1DATA) WHERE MS2PROD=target_mz
// Returns MS1 peaks from parent scans of MS2 spectra containing a product ion near target_mz.
// Default tolerance: 0.1 Da.
// Note: inlines cross-level logic rather than composing mzml_ms1_parent_peaks/mzml_ms2_peaks
// because DuckDB table macros cannot appear in subquery positions.
const std::string MZML_MS1_WHERE_MS2PROD = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms1_where_ms2prod(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1 "
    "AND spectrum_index IN ("
    "  SELECT DISTINCT ms1_scan_index"
    "  FROM mzml_peaks(relation)"
    "  WHERE ms_level = 2"
    "    AND mz_within(mz, target_mz, tolerance)"
    "    AND ms1_scan_index IS NOT NULL"
    "); ";

// mzml_ms2_where_ms1mz(relation, target_mz [, tolerance])
//
// MassQL: QUERY scaninfo(MS2DATA) WHERE MS1MZ=target_mz
// Returns MS2 peaks from child scans of MS1 spectra containing a peak near target_mz.
// Default tolerance: 0.1 Da.
// Note: inlines cross-level logic because DuckDB table macros cannot appear in subquery positions.
const std::string MZML_MS2_WHERE_MS1MZ = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms2_where_ms1mz(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 "
    "AND ms1_scan_index IN ("
    "  SELECT DISTINCT spectrum_index"
    "  FROM mzml_peaks(relation)"
    "  WHERE ms_level = 1"
    "    AND mz_within(mz, target_mz, tolerance)"
    "); ";

// mzml_ms1_where_ms2prec(relation, target_mz [, tolerance])
//
// MassQL: QUERY scansum(MS1DATA) WHERE MS2PREC=target_mz
// Returns MS1 peaks from parent scans of MS2 spectra whose precursor m/z matches.
// Default tolerance: 0.1 Da.
const std::string MZML_MS1_WHERE_MS2PREC = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms1_where_ms2prec(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1 "
    "AND spectrum_index IN ("
    "  SELECT DISTINCT ms1_scan_index FROM ("
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2"
    "  ) sub WHERE mz_within(sub.precursor_mz, target_mz, tolerance)"
    "  AND sub.ms1_scan_index IS NOT NULL"
    "); ";

// mzml_ms2_where_ms2prod_and_ms1mz(relation, prod_mz, ms1_mz [, tolerance])
//
// MassQL: QUERY scaninfo(MS2DATA) WHERE MS2PROD=prod_mz AND MS1MZ=ms1_mz
// Combined cross-level: MS2 scans that contain a product ion near prod_mz
// AND whose parent MS1 scan contains a peak near ms1_mz.
// Default tolerance: 0.1 Da.
const std::string MZML_MS2_WHERE_MS2PROD_AND_MS1MZ = // NOLINT
    "CREATE OR REPLACE MACRO mzml_ms2_where_ms2prod_and_ms1mz("
    "  relation, prod_mz, ms1_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 "
    "AND spectrum_index IN ("
    "  SELECT DISTINCT spectrum_index FROM ("
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2"
    "  ) sub WHERE mz_within(sub.mz, prod_mz, tolerance)"
    ") AND ms1_scan_index IN ("
    "  SELECT DISTINCT spectrum_index FROM ("
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1"
    "  ) sub2 WHERE mz_within(sub2.mz, ms1_mz, tolerance)"
    "); ";

// mzml_filter_mz(relation, target_mz, tolerance)
//
// MassQL FILTER clause: select only peaks near target_mz from a peak-level relation.
// Composes with WHERE macros: first select scans, then filter which peaks to return.
const std::string MZML_FILTER_MZ = // NOLINT
    "CREATE OR REPLACE MACRO mzml_filter_mz(relation, target_mz, tolerance) AS TABLE "
    "SELECT * FROM query_table(relation) "
    "WHERE mz > target_mz - tolerance AND mz < target_mz + tolerance; ";

// mzml_filter_nl(relation, target_nl, tolerance)
//
// MassQL FILTER clause for neutral loss: select only peaks whose neutral_loss is near target_nl.
const std::string MZML_FILTER_NL = // NOLINT
    "CREATE OR REPLACE MACRO mzml_filter_nl(relation, target_nl, tolerance) AS TABLE "
    "SELECT * FROM query_table(relation) "
    "WHERE neutral_loss > target_nl - tolerance AND neutral_loss < target_nl + tolerance; ";

// massdefect(mz_value)
//
// Returns the mass defect (fractional part) of an m/z value.
// Mass defect = mz_value - floor(mz_value).
const std::string MASSDEFECT = // NOLINT
    "CREATE OR REPLACE MACRO massdefect(mz_value) AS "
    "(mz_value - FLOOR(mz_value));";

// mz_massdefect_within(mz_value, min_defect, max_defect)
//
// Returns true if the mass defect of mz_value is within [min_defect, max_defect].
// Inclusive on both ends (matches MassQL range semantics).
const std::string MZ_MASSDEFECT_WITHIN = // NOLINT
    "CREATE OR REPLACE MACRO mz_massdefect_within(mz_value, min_defect, max_defect) AS "
    "((mz_value - FLOOR(mz_value)) BETWEEN min_defect AND max_defect);";

// mzml_x_offset_pair(relation, delta [, tolerance])
//
// MassQL: MS2PROD=X AND MS2PROD=X+delta
// Find MS2 spectra where two peaks differ by delta Da.
// Delegates to mzml_x_offset_ntuple([0, delta]).
// Default tolerance: 0.1 Da.
const std::string MZML_X_OFFSET_PAIR = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_offset_pair(relation, delta, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_x_offset_ntuple(relation, [0, delta], tolerance); ";

// mzml_x_offset_triplet(relation, delta2, delta3 [, tolerance])
//
// MassQL: MS2PROD=X AND MS2PROD=X+delta2 AND MS2PROD=X+delta3
// Find MS2 spectra where three peaks match X, X+delta2, X+delta3.
// Delegates to mzml_x_offset_ntuple([0, delta2, delta3]).
// Default tolerance: 0.1 Da.
const std::string MZML_X_OFFSET_TRIPLET = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_offset_triplet(relation, delta2, delta3, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_x_offset_ntuple(relation, [0, delta2, delta3], tolerance); ";

// mzml_x_offset_ntuple(relation, offsets [, tolerance])
//
// Generalized N-tuple offset matching: find MS2 spectra where peaks exist at
// X+offset[0], X+offset[1], ..., X+offset[N-1] simultaneously (ALL must match).
// X candidates are drawn from all MS2 peaks across all scans (cross-scan matching)
// and deduplicated with greedy 0.05 Da step.
// offsets must contain distinct values; duplicate offsets will produce no matches.
// Default tolerance: 0.1 Da.
const std::string MZML_X_OFFSET_NTUPLE = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_offset_ntuple(relation, offsets, tolerance := 0.1) AS TABLE "
    "WITH RECURSIVE ms2 AS ( "
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 AND intensity > 0 "
    "), "
    "x_candidates(x_val, next_min) AS ( "
    "    (SELECT mz, mz + 0.05 FROM ms2 ORDER BY mz LIMIT 1) "
    "    UNION ALL "
    "    (SELECT s.mz, s.mz + 0.05 "
    "     FROM x_candidates g "
    "     JOIN (SELECT DISTINCT mz FROM ms2) s ON s.mz >= g.next_min "
    "     ORDER BY s.mz "
    "     LIMIT 1) "
    "), "
    "offset_targets AS ( "
    "    SELECT UNNEST(offsets) AS offset_val "
    "), "
    "offset_count AS (SELECT COUNT(*) AS n FROM offset_targets), "
    "matched AS ( "
    "    SELECT xc.x_val, p.spectrum_index, ot.offset_val "
    "    FROM x_candidates xc "
    "    CROSS JOIN offset_targets ot "
    "    JOIN ms2 p ON p.mz > xc.x_val + ot.offset_val - tolerance "
    "              AND p.mz < xc.x_val + ot.offset_val + tolerance "
    "), "
    "qualifying AS ( "
    "    SELECT spectrum_index "
    "    FROM matched, offset_count "
    "    GROUP BY spectrum_index, x_val, offset_count.n "
    "    HAVING COUNT(DISTINCT offset_val) = offset_count.n "
    ") "
    "SELECT * FROM ms2 "
    "WHERE spectrum_index IN (SELECT spectrum_index FROM qualifying); ";

// mzml_x_prec_prod(relation, delta [, tolerance, min_intensity_pct])
//
// MassQL: MS2PREC=X AND MS2PROD=X-delta[:INTENSITYPERCENT=N]
// Per-spectrum: X is the spectrum's own precursor_mz. Find MS2 spectra that
// contain a product ion at precursor_mz - delta within tolerance.
// min_intensity_pct: product ion intensity as % of base peak (0 = no filter).
// Default tolerance: 0.1 Da.
const std::string MZML_X_PREC_PROD = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_prec_prod("
    "  relation, delta, tolerance := 0.1, min_intensity_pct := 0) AS TABLE "
    "WITH ms2 AS ( "
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 AND intensity > 0 "
    ") "
    "SELECT * FROM ms2 "
    "WHERE precursor_mz IS NOT NULL "
    "AND spectrum_index IN ( "
    "    SELECT DISTINCT p.spectrum_index "
    "    FROM ms2 p "
    "    WHERE p.precursor_mz IS NOT NULL "
    "      AND p.mz > p.precursor_mz - delta - tolerance "
    "      AND p.mz < p.precursor_mz - delta + tolerance "
    "      AND (min_intensity_pct = 0 OR p.i_norm >= min_intensity_pct / 100.0) "
    "); ";

// mzml_x_prec_massdefect(relation, min_defect, max_defect)
//
// MassQL: MS2PREC=X AND X=massdefect(min=min_defect, max=max_defect)
// X is the precursor_mz. Returns MS2 peaks from scans whose precursor_mz
// has mass defect in [min_defect, max_defect].
// Convenience wrapper for MassQL parity.
const std::string MZML_X_PREC_MASSDEFECT = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_prec_massdefect("
    "  relation, min_defect, max_defect) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) "
    "WHERE ms_level = 2 "
    "  AND precursor_mz IS NOT NULL "
    "  AND mz_massdefect_within(precursor_mz, min_defect, max_defect); ";

// mzml_x_ms1_ms2_prec(relation [, tolerance])
//
// MassQL: MS1MZ=X AND MS2PREC=X
// MS2 spectra whose precursor matches an MS1 peak in the parent scan.
// Default tolerance: 0.1 Da.
const std::string MZML_X_MS1_MS2_PREC = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_ms1_ms2_prec(relation, tolerance := 0.1) AS TABLE "
    "SELECT ms2.* FROM ("
    "  SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2"
    ") ms2 "
    "JOIN ("
    "  SELECT DISTINCT spectrum_index, mz FROM mzml_peaks(relation) WHERE ms_level = 1"
    ") ms1 ON ms2.ms1_scan_index = ms1.spectrum_index "
    "AND mz_within(ms1.mz, ms2.precursor_mz, tolerance); ";

// mzml_x_offset_pair_range(relation, delta, x_min, x_max [, tolerance])
//
// MassQL: MS2PROD=X AND MS2PROD=X+delta WHERE X=range(min=x_min, max=x_max)
// Same as mzml_x_offset_pair but X candidates are constrained to [x_min, x_max].
// Default tolerance: 0.1 Da.
const std::string MZML_X_OFFSET_PAIR_RANGE = // NOLINT
    "CREATE OR REPLACE MACRO mzml_x_offset_pair_range("
    "  relation, delta, x_min, x_max, tolerance := 0.1) AS TABLE "
    "WITH RECURSIVE ms2 AS ( "
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 AND intensity > 0 "
    "), "
    "ms2_in_range AS ( "
    "    SELECT * FROM ms2 WHERE mz >= x_min AND mz <= x_max "
    "), "
    "x_candidates(x_val, next_min) AS ( "
    "    (SELECT mz, mz + 0.05 FROM ms2_in_range ORDER BY mz LIMIT 1) "
    "    UNION ALL "
    "    (SELECT s.mz, s.mz + 0.05 "
    "     FROM x_candidates g "
    "     JOIN (SELECT DISTINCT mz FROM ms2_in_range) s ON s.mz >= g.next_min "
    "     ORDER BY s.mz "
    "     LIMIT 1) "
    ") "
    "SELECT * FROM ms2 "
    "WHERE spectrum_index IN ( "
    "    SELECT DISTINCT p1.spectrum_index "
    "    FROM x_candidates xc "
    "    JOIN ms2 p1 ON p1.mz > xc.x_val - tolerance AND p1.mz < xc.x_val + tolerance "
    "    JOIN ms2 p2 ON p2.spectrum_index = p1.spectrum_index "
    "        AND p2.mz > xc.x_val + delta - tolerance "
    "        AND p2.mz < xc.x_val + delta + tolerance "
    "); ";

// mzml_or_cardinality(relation, target_list, tolerance, min_card, max_card)
//
// MassQL: MS2PROD=(t1 OR t2 OR ...):CARDINALITY=range(min=min_card, max=max_card)
// Find MS2 spectra where the number of distinct target m/z values matched
// falls within [min_card, max_card].
// target_list is a LIST of DOUBLE values.
const std::string MZML_OR_CARDINALITY = // NOLINT
    "CREATE OR REPLACE MACRO mzml_or_cardinality("
    "  relation, target_list, tolerance, min_card, max_card) AS TABLE "
    "WITH ms2 AS ( "
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 "
    "), "
    "targets AS ( "
    "    SELECT UNNEST(target_list) AS target_mz "
    "), "
    "matches AS ( "
    "    SELECT ms2.spectrum_index, t.target_mz "
    "    FROM ms2, targets t "
    "    WHERE mz_within(ms2.mz, t.target_mz, tolerance) "
    "), "
    "cardinality_check AS ( "
    "    SELECT spectrum_index, COUNT(DISTINCT target_mz) AS match_count "
    "    FROM matches "
    "    GROUP BY spectrum_index "
    "    HAVING match_count >= min_card AND match_count <= max_card "
    ") "
    "SELECT * FROM ms2 "
    "WHERE spectrum_index IN (SELECT spectrum_index FROM cardinality_check); ";

// mzml_peak_pair is now a C++ table function — see src/mzml_peak_pair_function.cpp

// mzml_i_norm(intensity_array, base_peak_intensity)
//
// Max-normalize an intensity array by dividing each value by the base peak intensity.
// Returns NULL if base_peak_intensity is NULL (natural NULL propagation).
const std::string MZML_I_NORM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_i_norm(intensity_array, base_peak_intensity) AS ("
    "  list_transform(intensity_array, lambda x: x / base_peak_intensity)"
    ");";

// mzml_i_tic_norm(intensity_array, total_ion_current)
//
// TIC-normalize an intensity array by dividing each value by the total ion current.
// Returns NULL if total_ion_current is NULL (natural NULL propagation).
const std::string MZML_I_TIC_NORM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_i_tic_norm(intensity_array, total_ion_current) AS ("
    "  list_transform(intensity_array, lambda x: x / total_ion_current)"
    ");";

// mzml_excluded_ms2prod(relation, target_mz, tolerance := 0.1)
//
// MassQL: MS2PROD=target_mz:EXCLUDED
// Returns MS2 peaks from scans that do NOT contain a product ion near target_mz.
const std::string MZML_EXCLUDED_MS2PROD = // NOLINT
    "CREATE OR REPLACE MACRO mzml_excluded_ms2prod(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 "
    "AND spectrum_index NOT IN ("
    "  SELECT DISTINCT spectrum_index"
    "  FROM mzml_peaks(relation)"
    "  WHERE ms_level = 2"
    "    AND mz_within(mz, target_mz, tolerance)"
    "); ";

// mzml_excluded_ms1mz(relation, target_mz, tolerance := 0.1)
//
// MassQL: MS1MZ=target_mz:EXCLUDED
// Returns MS1 peaks from scans that do NOT contain a peak near target_mz.
const std::string MZML_EXCLUDED_MS1MZ = // NOLINT
    "CREATE OR REPLACE MACRO mzml_excluded_ms1mz(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1 "
    "AND spectrum_index NOT IN ("
    "  SELECT DISTINCT spectrum_index"
    "  FROM mzml_peaks(relation)"
    "  WHERE ms_level = 1"
    "    AND mz_within(mz, target_mz, tolerance)"
    "); ";

// mzml_excluded_ms2prec(relation, target_mz, tolerance := 0.1)
//
// MassQL: MS2PREC=target_mz:EXCLUDED
// Returns MS2 peaks from scans whose precursor m/z is NOT near target_mz.
// MS2 scans with NULL precursor_mz are retained (they cannot match the target).
const std::string MZML_EXCLUDED_MS2PREC = // NOLINT
    "CREATE OR REPLACE MACRO mzml_excluded_ms2prec(relation, target_mz, tolerance := 0.1) AS TABLE "
    "SELECT * FROM mzml_peaks(relation) WHERE ms_level = 2 "
    "AND spectrum_index NOT IN ("
    "  SELECT DISTINCT spectrum_index"
    "  FROM mzml_peaks(relation)"
    "  WHERE ms_level = 2"
    "    AND mz_within(precursor_mz, target_mz, tolerance)"
    "); ";

// mzml_isotope_pattern(relation, target_offsets, target_ratios, target_tol_pcts, mz_tolerance := 0.1)
//
// MassQL Y variable: isotope pattern matching.
// Finds MS1 spectra where a reference peak (X) has satellite peaks at specified Da offsets
// with intensity ratios (relative to X) within specified percentage tolerances.
// Returns all MS1 peaks from matching spectra.
//
// Parameters:
// target_offsets : LIST of DOUBLE, Da offsets from reference peak (e.g. [1.0, 2.0])
// target_ratios  : LIST of DOUBLE, expected intensity ratios (e.g. [0.5, 0.1])
// target_tol_pcts: LIST of DOUBLE, percentage tolerances (e.g. [30, 30] for ±30%)
// mz_tolerance   : DOUBLE, m/z matching tolerance in Da (default 0.1)
//
// Intensity semantics: Y = SUM(intensity) for all peaks within mz_tolerance of reference
// or satellite m/z, per MassQL convention. In centroid data this is typically one peak;
// in profile data or with wide mz_tolerance, multiple peaks may be summed.
//
// Ratio window uses strict inequalities (consistent with mz_within).
//
// Performance note: the recursive CTE enumerates all distinct m/z values across all spectra
// (deduplicated with 0.05 Da step). For large datasets, pre-filter to a retention time
// window or m/z range before calling this macro.
const std::string MZML_ISOTOPE_PATTERN = // NOLINT
    "CREATE OR REPLACE MACRO mzml_isotope_pattern("
    "  relation, target_offsets, target_ratios, target_tol_pcts, "
    "  mz_tolerance := 0.1) AS TABLE "
    "WITH RECURSIVE ms1 AS ( "
    "    SELECT * FROM mzml_peaks(relation) WHERE ms_level = 1 AND intensity > 0 "
    "), "
    "x_candidates(x_val, next_min) AS ( "
    "    (SELECT mz, mz + 0.05 FROM ms1 ORDER BY mz LIMIT 1) "
    "    UNION ALL "
    "    (SELECT s.mz, s.mz + 0.05 "
    "     FROM x_candidates g "
    "     JOIN (SELECT DISTINCT mz FROM ms1) s ON s.mz >= g.next_min "
    "     ORDER BY s.mz "
    "     LIMIT 1) "
    "), "
    "targets AS ( "
    "    SELECT UNNEST(target_offsets) AS offset_val, "
    "           UNNEST(target_ratios) AS ratio, "
    "           UNNEST(target_tol_pcts) AS tol_pct "
    "), "
    "target_count AS (SELECT COUNT(*) AS n FROM targets), "
    "ref_intensity AS ( "
    "    SELECT xc.x_val, ms1.spectrum_index, SUM(ms1.intensity) AS y_ref "
    "    FROM x_candidates xc "
    "    JOIN ms1 ON ms1.mz > xc.x_val - mz_tolerance "
    "            AND ms1.mz < xc.x_val + mz_tolerance "
    "    GROUP BY xc.x_val, ms1.spectrum_index "
    "), "
    "target_intensity AS ( "
    "    SELECT ri.x_val, ri.spectrum_index, ri.y_ref, "
    "           t.offset_val, t.ratio, t.tol_pct, "
    "           SUM(ms1.intensity) AS t_int "
    "    FROM ref_intensity ri "
    "    CROSS JOIN targets t "
    "    JOIN ms1 ON ms1.spectrum_index = ri.spectrum_index "
    "            AND ms1.mz > ri.x_val + t.offset_val - mz_tolerance "
    "            AND ms1.mz < ri.x_val + t.offset_val + mz_tolerance "
    "    GROUP BY ri.x_val, ri.spectrum_index, ri.y_ref, "
    "             t.offset_val, t.ratio, t.tol_pct "
    "), "
    "matched_targets AS ( "
    "    SELECT x_val, spectrum_index, offset_val "
    "    FROM target_intensity "
    "    WHERE t_int > y_ref * ratio * (1.0 - tol_pct / 100.0) "
    "      AND t_int < y_ref * ratio * (1.0 + tol_pct / 100.0) "
    "), "
    "qualifying AS ( "
    "    SELECT spectrum_index "
    "    FROM matched_targets, target_count "
    "    GROUP BY spectrum_index, x_val, target_count.n "
    "    HAVING COUNT(DISTINCT offset_val) = target_count.n "
    ") "
    "SELECT * FROM ms1 "
    "WHERE spectrum_index IN (SELECT spectrum_index FROM qualifying); ";

// genome_coverage(alignments, subject_total_length, subject_genome_id)
//
// Compute genome coverage from alignment data by:
// 1. Compressing overlapping alignment intervals per reference contig
// 2. Mapping contigs to genomes and summing covered bases
// 3. Joining with total genome lengths to compute proportion covered
//
// Parameters:
// alignments : a relation with columns: reference (VARCHAR), position (BIGINT),
//     stop_position (BIGINT)
// subject_total_length : a relation with columns: genome_id, total_length (BIGINT)
// subject_genome_id : a relation with columns: contig_id (VARCHAR), genome_id
//
// genome_id is passed through without casting -- VARCHAR, BIGINT and UUID all
// work, and the output keeps the type given in subject_genome_id. See the
// genome_coverage_per_sample comment below on mixed types across the two
// reference relations; the behavior is identical here, same code path.
//
// Returns: genome_id (as supplied), covered (BIGINT), proportion_covered (DOUBLE)
//
// Raises, rather than returning a wrong or missing answer, when a covered genome
// has no usable denominator:
//   - no total_length row for a genome that has coverage
//   - total_length that is NULL, zero or negative
// A genome listed in subject_total_length but with no coverage is not an error;
// it simply produces no row.
//
// `covered` is cast back to BIGINT deliberately -- do not remove the cast.
// DuckDB widens the return type of SUM unconditionally (SUM(BIGINT) and even
// SUM(INTEGER) yield HUGEINT), which is pointless here: compress_intervals
// merges overlaps per reference, so each contig contributes at most its own
// length and `covered` is bounded by the genome length -- BIGINT holds ~2.96
// billion human genomes. The width was not free: HUGEINT is 16 bytes rather
// than 8, and COPY ... TO a Parquet file writes it as DOUBLE, so an integer
// base count silently arrived as a float.
const std::string GENOME_COVERAGE = // NOLINT
    "CREATE OR REPLACE MACRO genome_coverage(alignments, subject_total_length, subject_genome_id) AS TABLE "
    "WITH compressed_intervals AS ( "
    "    SELECT "
    "        reference, "
    "        UNNEST(compress_intervals(position, stop_position)) AS ci "
    "    FROM query_table(alignments) "
    "    GROUP BY reference "
    "), "
    "internal_coverage AS ( "
    "    SELECT "
    "        sg.genome_id, "
    "        SUM(ci.stop - ci.start) AS covered_internal "
    "    FROM compressed_intervals "
    // SELECT DISTINCT, not a bare query_table(): a duplicated (contig_id,
    // genome_id) row would otherwise multiply every interval for that contig
    // before any GROUP BY sees it, reporting 11 covered bases as 22. Mapping
    // tables assembled by concatenating per-genome files hit this routinely, and
    // the inflated result looks like an ordinary low-coverage genome. The dedupe
    // is on the PAIR, so a contig legitimately mapped to two different genomes
    // still counts toward both.
    "    JOIN (SELECT DISTINCT contig_id, genome_id FROM query_table(subject_genome_id)) sg "
    "      ON reference = sg.contig_id "
    "    GROUP BY sg.genome_id, reference "
    "), "
    "total_coverage AS ( "
    "    SELECT "
    "        genome_id, "
    "        SUM(covered_internal) AS covered "
    "    FROM internal_coverage "
    "    GROUP BY genome_id "
    ") "
    // SUM widens to INT128 in DuckDB regardless of how small the values are, so this
    // returned HUGEINT while the docs promised BIGINT. Narrowing is always safe:
    // covered is bounded by total_length, and the largest known genome is ~150 Gbp,
    // about 38 bits.
    "SELECT "
    "    tc.genome_id, "
    "    tc.covered::BIGINT AS covered, "
    "    tc.covered::DOUBLE / tl.total_length AS proportion_covered "
    "FROM total_coverage tc "
    // LEFT JOIN + an explicit guard, not an inner join: a covered genome with no
    // usable total_length must be reported, not silently dropped or divided into
    // garbage. The guard lives in WHERE rather than the projection because
    // DuckDB prunes unused projections -- inside SELECT it would not fire for
    // `SELECT covered FROM genome_coverage(...)`, the very query that would
    // otherwise hide the problem. Same reasoning as read_gff's malformed-line
    // check above. A genome present in subject_total_length but uncovered here
    // simply produces no row; only a genome we must divide for and cannot is
    // fatal, so an empty result stays empty rather than raising.
    "LEFT JOIN query_table(subject_total_length) tl "
    "  ON tc.genome_id = tl.genome_id "
    "WHERE CASE "
    "        WHEN tl.genome_id IS NULL "
    "          THEN error(printf('genome_coverage: no total_length entry for genome_id ''%s''', "
    "                            tc.genome_id::VARCHAR)) "
    "        WHEN tl.total_length IS NULL OR tl.total_length <= 0 "
    "          THEN error(printf('genome_coverage: total_length must be positive for genome_id ''%s'' "
    "(got %s)', tc.genome_id::VARCHAR, COALESCE(tl.total_length::VARCHAR, 'NULL'))) "
    "        ELSE TRUE "
    "      END;";

// genome_coverage_per_sample(alignments, subject_total_length, subject_genome_id)
//
// Per-sample sibling of genome_coverage. Identical arithmetic, but every step is
// partitioned by a `sample_id` column that `alignments` must additionally carry:
// intervals compress per (sample_id, reference) and covered bases sum per
// (sample_id, genome_id).
//
// Use this whenever `alignments` holds more than one sample. genome_coverage has
// no sample dimension, so it POOLS multi-sample input into a single row per
// genome -- two samples covering [10,21)+[30,41) and [10,21) of a 100 bp genome
// yield one row at 0.22 instead of 0.22 and 0.11. That is not a bug in
// genome_coverage (it cannot see samples), but it does mean every sample below
// the pooled maximum is silently over-reported if you feed it pooled input.
//
// The two macros are deliberate near-duplicates rather than one parameterized
// macro: a DuckDB macro cannot vary its output schema, the single-sample form
// has no sample_id column to delegate with, and a macro cannot inject a constant
// column into a query_table() reference. Any change to the coverage arithmetic
// here must be mirrored in GENOME_COVERAGE above, and vice versa.
//
// Parameters:
// alignments : a relation with columns: sample_id (any type), reference
//     (VARCHAR), position (BIGINT), stop_position (BIGINT)
// subject_total_length : a relation with columns: genome_id, total_length (BIGINT)
// subject_genome_id : a relation with columns: contig_id (VARCHAR), genome_id
//
// sample_id and genome_id are passed through without casting, so VARCHAR,
// BIGINT and UUID identifiers all survive end to end; the output keeps the type
// given in subject_genome_id. The two reference relations do NOT have to declare
// genome_id with the same type -- it is a join key, so DuckDB implicit-casts and
// a BIGINT 77 matches a VARCHAR '77' (pinned by Test 14 of
// test/sql/genome_coverage_per_sample.test). Matching types are still worth
// preferring: an unconvertible value fails as a cast error naming the column
// rather than as a clear diagnosis of the mismatch.
//
// A (sample, genome) pair with no qualifying alignments produces no row at all,
// rather than a zero-coverage row.
//
// Returns: sample_id (as supplied), genome_id (as supplied), covered (BIGINT),
//     proportion_covered (DOUBLE)
//
// Raises on a NULL sample_id, and on the same unusable-denominator cases as
// genome_coverage (missing total_length row; NULL, zero or negative
// total_length). The two macros are kept in lockstep on all of it.
//
// `covered` is cast back to BIGINT for the same reason as genome_coverage above
// -- DuckDB widens SUM unconditionally, and the value is bounded by the genome
// length. Keep both macros in agreement; a type mismatch between siblings would
// be a nasty surprise for anyone UNION-ing their results.
const std::string GENOME_COVERAGE_PER_SAMPLE = // NOLINT
    "CREATE OR REPLACE MACRO genome_coverage_per_sample(alignments, subject_total_length, subject_genome_id) AS "
    "TABLE "
    "WITH compressed_intervals AS ( "
    "    SELECT "
    "        sample_id, "
    "        reference, "
    "        UNNEST(compress_intervals(position, stop_position)) AS ci "
    "    FROM query_table(alignments) "
    // Reject NULL sample_id rather than letting it form its own group. A NULL
    // sample is not a meaningful partition, and grouping it would emit a
    // plausible-looking row for what is almost always a data error. Every other
    // per-sample function in the extension rejects it (see
    // docs/internals/per-sample-pattern.md); the message reuses that family's
    // substring so callers see one behavior. In WHERE, not the projection, so
    // projection pruning cannot silence it.
    "    WHERE CASE WHEN sample_id IS NULL "
    "               THEN error('genome_coverage_per_sample: NULL values in sample_id column "
    "''sample_id''') "
    "               ELSE TRUE END "
    "    GROUP BY sample_id, reference "
    "), "
    "internal_coverage AS ( "
    "    SELECT "
    "        sample_id, "
    "        sg.genome_id, "
    "        SUM(ci.stop - ci.start) AS covered_internal "
    "    FROM compressed_intervals "
    // Deduplicated for the same reason as genome_coverage above -- a repeated
    // (contig_id, genome_id) row would double every covered count.
    "    JOIN (SELECT DISTINCT contig_id, genome_id FROM query_table(subject_genome_id)) sg "
    "      ON reference = sg.contig_id "
    "    GROUP BY sample_id, sg.genome_id, reference "
    "), "
    "total_coverage AS ( "
    "    SELECT "
    "        sample_id, "
    "        genome_id, "
    "        SUM(covered_internal) AS covered "
    "    FROM internal_coverage "
    "    GROUP BY sample_id, genome_id "
    ") "
    "SELECT "
    "    tc.sample_id, "
    "    tc.genome_id, "
    "    tc.covered::BIGINT AS covered, "
    "    tc.covered::DOUBLE / tl.total_length AS proportion_covered "
    "FROM total_coverage tc "
    // Mirrors genome_coverage's guard; see the reasoning there.
    "LEFT JOIN query_table(subject_total_length) tl "
    "  ON tc.genome_id = tl.genome_id "
    "WHERE CASE "
    "        WHEN tl.genome_id IS NULL "
    "          THEN error(printf('genome_coverage_per_sample: no total_length entry for genome_id "
    "''%s''', tc.genome_id::VARCHAR)) "
    "        WHEN tl.total_length IS NULL OR tl.total_length <= 0 "
    "          THEN error(printf('genome_coverage_per_sample: total_length must be positive for "
    "genome_id ''%s'' (got %s)', tc.genome_id::VARCHAR, COALESCE(tl.total_length::VARCHAR, 'NULL'))) "
    "        ELSE TRUE "
    "      END;";

// bin_of(pos, n_bins, genome_length)
// bin_start(b, n_bins, genome_length)
// interval_bins(start, stop, n_bins, genome_length)
//
// Divide a genome into n_bins near-equal-width bins. bin_of maps a position to its bin,
// bin_start maps a bin back to its first position, and interval_bins reports which bins
// a half-open interval touches. Used to find genomic regions whose coverage differs
// between sample groups: bin, count samples/reads per bin per group, then rank bins by
// the standard deviation of the per-group counts -- all plain SQL once the fan-out
// exists (issue #215).
//
// Widths are near-equal, not equal: they differ by at most one base when n_bins does not
// divide genome_length. That is the entire reason bin_start cannot be computed by the
// obvious floor formula -- see below.
//
// Bin indices are 0-BASED. They are array positions, not genomic coordinates, so they
// deliberately cut against the 1-based coordinate convention used everywhere else here.
//
// Genomic positions stay 1-based half-open, matching read_alignments and the normalized
// stop_position rule in docs/internals/architecture.md: bin b spans
// [bin_start(b), bin_start(b+1)) and its width is the plain difference, no +1. That makes
// bin_start(0) = 1 and bin_start(n_bins) = genome_length + 1, so (bin_start(b),
// bin_start(b+1)) is directly usable as a half-open (start, stop) pair.
//
// bin_of is the definition; bin_start is its left inverse:
//
//   bin_of(p)    = ((p - 1) * n_bins) // genome_length
//   bin_start(b) = ceil(b * genome_length / n_bins) + 1
//
// The invariant is CONTAINMENT, not round-tripping in both directions:
//
//   bin_start(bin_of(p)) <= p < bin_start(bin_of(p) + 1)   for every p in [1, genome_length]
//
// bin_of(bin_start(b)) = b fails when n_bins > genome_length, because zero-width bins make
// bin_start non-injective (n_bins=5, genome_length=3: bin_start(2) = 3 but bin_of(3) = 3).
//
// bin_start must use ceil. With floor, genome_length=10 / n_bins=4 gives bounds
// [1,3) [3,6) [6,8) [8,11), which places position 3 in bin 1 while bin_of(3) is 0 --
// the two disagree at 4 of 10 positions. bin_start is the join key callers carry
// around, so a mismatch there shifts every downstream bin and the ranking still looks
// perfectly plausible. Pinned by the containment assertions in
// test/sql/interval_bins.test, which report 23 violations against the floor variant.
//
// Every operand is normalized with ::BIGINT before it is used in arithmetic OR in a
// bounds comparison, and the normalized value is what the error messages quote. Two
// reasons. read_gff and read_ncbi_annotation emit INTEGER positions, and an unwidened
// (position * n_bins) overflows INT32 at ~4.3 Mbp / 500 bins -- an ordinary bacterial
// genome. And a genome_length arriving as DOUBLE or DECIMAL (an inferred CSV/Parquet
// column) must compare against the same rounded value the arithmetic uses, or the clamp
// and the bounds check disagree and the caller sees an error about a position they never
// supplied.
//
// NULL in any argument yields NULL, following the SQL convention compress_intervals uses.
// Note what that means downstream: UNNEST(NULL) emits no rows, so a genome_length that
// arrived NULL from a missed join drops those reads silently. The positivity check below
// cannot catch that -- a join miss produces NULL, not 0. Guard the join if it matters.
const std::string BIN_OF = // NOLINT
    "CREATE OR REPLACE MACRO bin_of(pos, n_bins, genome_length) AS ( "
    "  CASE "
    "    WHEN pos IS NULL OR n_bins IS NULL OR genome_length IS NULL THEN NULL "
    "    WHEN n_bins::BIGINT < 1 OR genome_length::BIGINT < 1 "
    "      THEN error(printf('bin_of: n_bins and genome_length must be positive "
    "(got n_bins=%s, genome_length=%s)', "
    "                        n_bins::BIGINT::VARCHAR, genome_length::BIGINT::VARCHAR)) "
    // The hint is not decoration. stop_position is an EXCLUSIVE end, so a read reaching
    // the reference end legitimately carries genome_length + 1, and
    // `SELECT bin_of(position, ...), bin_of(stop_position, ...)` is the natural thing to
    // write. It passes on small fixtures and aborts on a real BAM, so the error has to
    // name the fix.
    "    WHEN pos::BIGINT < 1 OR pos::BIGINT > genome_length::BIGINT "
    "      THEN error(printf('bin_of: position %s outside genome [1, %s] "
    "(if this is an exclusive stop_position, pass stop_position - 1)', "
    "                        pos::BIGINT::VARCHAR, genome_length::BIGINT::VARCHAR)) "
    "    ELSE ((pos::BIGINT - 1) * n_bins::BIGINT) // genome_length::BIGINT "
    "  END);";

// b may equal n_bins: that is the exclusive end of the last bin, genome_length + 1, and
// callers need it to close the final region.
const std::string BIN_START = // NOLINT
    "CREATE OR REPLACE MACRO bin_start(b, n_bins, genome_length) AS ( "
    "  CASE "
    "    WHEN b IS NULL OR n_bins IS NULL OR genome_length IS NULL THEN NULL "
    "    WHEN n_bins::BIGINT < 1 OR genome_length::BIGINT < 1 "
    "      THEN error(printf('bin_start: n_bins and genome_length must be positive "
    "(got n_bins=%s, genome_length=%s)', "
    "                        n_bins::BIGINT::VARCHAR, genome_length::BIGINT::VARCHAR)) "
    "    WHEN b::BIGINT < 0 OR b::BIGINT > n_bins::BIGINT "
    "      THEN error(printf('bin_start: bin index %s outside [0, %s]', "
    "                        b::BIGINT::VARCHAR, n_bins::BIGINT::VARCHAR)) "
    "    ELSE ((b::BIGINT * genome_length::BIGINT + n_bins::BIGINT - 1) "
    "          // n_bins::BIGINT) + 1 "
    "  END);";

// Two regimes, because range() materializes every element before anything can prune it.
//
// n_bins <= genome_length: every bin is at least one base wide, so no bin in the spanned
// range can be empty and the range IS the answer. Cost is the number of bins touched.
//
// n_bins > genome_length: some bins are zero-width. That is legitimate -- a fixed bin
// count applied across genomes of very different sizes -- but those bins hold no bases,
// so an interval cannot touch them and they must not be emitted. Enumerating the range
// and filtering would cost O(n_bins) (measured 19 ms/row at n_bins=1e6, genome_length=3);
// enumerating positions instead costs O(interval length), which in this regime is bounded
// by genome_length < n_bins. Same answer, 269x faster on that case.
//
// The bin arithmetic is inlined rather than delegated to bin_of, duplicating the floor
// formula on purpose. The reason is COST, not safety:
//
//   - Delegating re-runs bin_of's NULL and bounds guards once per enumerated position,
//     measured at ~1.4x, when all of them are already established here and the clamp below
//     already guarantees every p is inside [1, genome_length].
//
//   - A BARE `bin_of(...)` reference would also resolve in the CALLER's catalog, so a user
//     macro of that name silently hijacks every result (verified: a user
//     `CREATE MACRO bin_of(...) AS 999` turns interval_bins(1,5,4,10) into [999]).
//
// An earlier version of this comment claimed the hijack made delegation impossible. That
// is FALSE and was corrected here. Schema-qualifying the callee defeats it completely, and
// the qualification cannot be spoofed:
//
//   CREATE MACRO delegator_q(p,nb,gl) AS system.main.bin_of(p,nb,gl);
//   CREATE MACRO delegator_u(p,nb,gl) AS bin_of(p,nb,gl);
//   CREATE MACRO bin_of(p,nb,gl) AS 999;      -- shadow created AFTER both
//   SELECT delegator_q(1,4,10), delegator_u(1,4,10);   -->   0   |   999
//
//   ATTACH ':memory:' AS system;               -->  "reserved name"
//   CREATE MACRO system.main.bin_of(...);      -->  "Cannot create entry in system catalog"
//
// (miint's macros do live in system.main -- confirmed via duckdb_functions().) So the
// file's mzml family is right to delegate: MZML_X_OFFSET_PAIR and MZML_X_OFFSET_TRIPLET
// call mzml_x_offset_ntuple, and READ_GFF calls parse_gff_attributes. There is no
// project-wide "never delegate" rule, and this comment should not be read as one.
//
// A shared, schema-qualified `_miint_bin_index` helper would remove the triplication at no
// measured cost and is the obvious follow-up; it is deliberately NOT done here to keep this
// change to correcting the record. Until then, any change to the formula must be mirrored
// in BIN_OF above and in both branches below -- test/sql/interval_bins.test Test 5
// cross-checks interval_bins against bin_of over both regimes, including a non-dividing
// (3,5) pair, so the copies cannot drift silently.
const std::string INTERVAL_BINS = // NOLINT
    "CREATE OR REPLACE MACRO interval_bins(start, stop, n_bins, genome_length) AS ( "
    "  CASE "
    "    WHEN start IS NULL OR stop IS NULL OR n_bins IS NULL OR genome_length IS NULL "
    "      THEN NULL "
    "    WHEN n_bins::BIGINT < 1 OR genome_length::BIGINT < 1 "
    "      THEN error(printf('interval_bins: n_bins and genome_length must be positive "
    "(got n_bins=%s, genome_length=%s)', "
    "                        n_bins::BIGINT::VARCHAR, genome_length::BIGINT::VARCHAR)) "
    // Clamp to the genome, then treat an empty or inverted interval as touching nothing
    // so UNNEST drops the row.
    "    WHEN LEAST(stop::BIGINT, genome_length::BIGINT + 1) "
    "         <= GREATEST(start::BIGINT, 1) THEN []::BIGINT[] "
    "    WHEN n_bins::BIGINT <= genome_length::BIGINT "
    "      THEN range(((GREATEST(start::BIGINT, 1) - 1) * n_bins::BIGINT) "
    "                 // genome_length::BIGINT, "
    "                 ((LEAST(stop::BIGINT, genome_length::BIGINT + 1) - 2) "
    "                  * n_bins::BIGINT) // genome_length::BIGINT + 1) "
    "    ELSE list_sort(list_distinct(list_transform( "
    "           range(GREATEST(start::BIGINT, 1), "
    "                 LEAST(stop::BIGINT, genome_length::BIGINT + 1)), "
    "           lambda p: ((p - 1) * n_bins::BIGINT) // genome_length::BIGINT))) "
    "  END);";

// region_presence(positions, regions, samples)
//
// Is a genomic region present in a sample? The answer is THREE-state, and the third
// state is what makes downstream statistics honest (issue #216):
//
//   present         the sample has a covered interval overlapping the region
//   absent          the sample covers the genome, but none of it overlaps the region
//   not applicable  the sample has no coverage of the genome at all
//
// Collapsing 'not applicable' into 'absent' conflates "the organism is here and this
// region is missing from it" -- a strain-content claim -- with "the organism is not
// here", a detection failure. Those must not be pooled, and the distinction is easy to
// lose when hand-rolling the SQL.
//
// Parameters:
// positions : a relation with columns: sample_id (any type), genome_id (any type),
//     start, stop (any integral type). Intervals are 1-based half-open, matching
//     compress_intervals output. Note these are NOT read_alignments' column names --
//     positions is a covered-interval relation, so the caller renames
//     reference/position/stop_position on the way in. The docs show the rename.
// regions : a relation with columns: genome_id (any type), region_start, region_stop
//     (any integral type), region_id (any type). Half-open on the same convention.
// samples : a relation with a sample_id column -- the full cohort roster.
//
// Returns: sample_id (from samples), genome_id / region_id (from regions, uncast),
//     region_start / region_stop (from regions, as BIGINT), state (VARCHAR:
//     'present' / 'absent' / 'not applicable')
//
// One row per region ROW per sample. The coordinates are emitted, not just region_id,
// because the presence join is keyed on them: two rows sharing a region_id but covering
// different spans are scored independently, and without the coordinates in the output
// those rows would be indistinguishable -- a PIVOT would silently collapse them and a
// join back to metadata would duplicate samples. When region_ids are unique, which is
// the normal case, (sample_id, region_id) still keys the output.
//
// Long form, deliberately, rather than a sample x region matrix: it is what read_biom,
// woltka_ogu and the diversity functions already consume, and callers who want the
// matrix can PIVOT. micov builds the wide form internally with two PIVOTs and a Polars
// coalesce because pivot columns cannot be named in advance; long form makes that
// workaround unnecessary.
//
// samples is REQUIRED, not optional. 'not applicable' cannot be derived from positions
// alone -- a sample with no coverage contributes no rows, so it is indistinguishable
// from a sample that is not in the study. Without the roster the function could only
// emit 'present' and 'absent', which is the very conflation it exists to prevent. A
// macro also cannot vary its body on whether an argument was supplied. region_id is
// likewise required rather than synthesized, because a macro body cannot conditionally
// reference a column that may not exist.
//
// positions does NOT need to be pre-compressed. Presence asks only whether ANY interval
// overlaps, so overlapping input intervals give the same answer as their union, and the
// DISTINCTs below absorb the duplication.
//
// EVERY coordinate operand is normalized with ::BIGINT before it is compared, not just
// before it is used arithmetically -- the same rule bin_of/bin_start/interval_bins
// follow. Without it, coordinates loaded as VARCHAR (a BED file, an inferred CSV column)
// compare lexicographically: '9' < '10' is FALSE, so real intervals get discarded as
// degenerate, and '9' <= '100' is FALSE, so an inverted region passes validation. The
// result is a silent, cohort-wide 'not applicable' -- the exact conflation this function
// exists to prevent.
//
// NULL handling in positions splits on what the column means. The interval columns have
// a defined NULL meaning -- a NULL start or stop marks a sample that exists with zero
// coverage -- so those rows drop out without counting as coverage. The identity columns
// do not: a NULL sample_id or genome_id is an unattributable interval, and silently
// dropping it would downgrade real coverage to a non-detection, so it raises. An
// interval with start >= stop covers no bases and drops out for the same reason a NULL
// one does. That is deliberately NOT symmetric with the hard error on an empty region:
// a region is the QUESTION, and an unanswerable question must not be silently answered
// 'absent', whereas a degenerate interval is just data that contributes nothing.
//
// Samples present in positions but absent from the roster are ignored: the roster
// defines the cohort.
//
// The overlap test is `start < region_stop AND stop > region_start`. micov's predicate
// (_view.py:330) is `pos.start <= fm.stop AND pos.stop > fm.start`, which admits a
// zero-width touch at the region's exclusive end and reports a read starting exactly
// where the region finishes as present. It shares no base with the region.
//
// Three structural details that are load-bearing, not stylistic:
//
// 1. The positions guard and the degeneracy filter are ONE CASE, not two predicates.
//    Split across two CTEs, DuckDB collapses them into a single FILTER and evaluates the
//    cheap interval test first, so a row that is BOTH unattributable and degenerate is
//    discarded before the error branch is ever reached -- the guard silently does
//    nothing on exactly the rows most likely to be malformed.
// 2. _rp_checked is MATERIALIZED. That is an optimizer barrier: it stops the roster
//    semi-join below from being reordered ahead of the guard (a NULL sample_id matches
//    no roster entry, so the semi-join would prune the offending row before it could
//    raise). It also pins positions to a single scan, which matters for a relation that
//    can only be consumed once, such as an Arrow RecordBatchReader.
// 3. CTE names are prefixed. query_table() resolves in the caller's scope, so a CTE can
//    shadow a user relation of the same name passed as an argument, and the failure is an
//    unattributable binder error about a missing column. `roster` is exactly the name the
//    docs use, so the hazard is real here.
//
//    The PRECISE rule -- established later, while testing region_coverage, and stated in
//    full at _rc_reg below -- is narrower than "an unprefixed CTE shadows an argument":
//    a CTE shadows an argument if and only if it is DECLARED BEFORE the query_table() call
//    that resolves that argument. Whether the argument is a CTE, a view or a base table is
//    irrelevant. An earlier version of this comment stated the broad form as fact; it is
//    wrong, and prefixing is kept here as a blanket precaution rather than because every
//    name listed above could actually collide.
//
// One thing MATERIALIZED does NOT fix: the validation guards are still ordinary SQL
// expressions, so the optimizer may push a CALLER's filter beneath them and a malformed
// row outside that filter will not raise (verified: materializing does not change this).
// That never affects the rows returned -- each output row depends only on its own region
// and sample -- but validation covers what the query reads, not the whole relation.
const std::string REGION_PRESENCE = // NOLINT
    "CREATE OR REPLACE MACRO region_presence(positions, regions, samples) AS TABLE "
    "WITH _rp_roster AS ( "
    "    SELECT DISTINCT sample_id FROM query_table(samples) "
    // The message reuses the per-sample family's substring
    // (docs/internals/per-sample-pattern.md) so callers see one behavior.
    "    WHERE CASE WHEN sample_id IS NULL "
    "               THEN error('region_presence: NULL values in sample_id column "
    "''sample_id''') "
    "               ELSE TRUE END "
    "), "
    "_rp_reg AS ( "
    "    SELECT DISTINCT genome_id, region_start::BIGINT AS region_start, "
    "           region_stop::BIGINT AS region_stop, region_id "
    "    FROM query_table(regions) "
    // An unusable region must not be scored: a NULL bound makes the overlap predicate
    // NULL and an empty or inverted region can never overlap anything, so every sample
    // would come back 'absent' -- indistinguishable from a real negative result. NULLs
    // are checked first so the comparisons below never see one. Zero-width and inverted
    // are separate branches because they have different causes and different remedies:
    // filtering is right for a bin expansion and wrong for transposed coordinates.
    "    WHERE CASE "
    "            WHEN genome_id IS NULL OR region_start IS NULL "
    "                 OR region_stop IS NULL OR region_id IS NULL "
    "              THEN error(printf('region_presence: genome_id, region_start, "
    "region_stop and region_id must not be NULL (got genome_id=%s, region_start=%s, "
    "region_stop=%s, region_id=%s)', "
    "                                COALESCE(genome_id::VARCHAR, 'NULL'), "
    "                                COALESCE(region_start::VARCHAR, 'NULL'), "
    "                                COALESCE(region_stop::VARCHAR, 'NULL'), "
    "                                COALESCE(region_id::VARCHAR, 'NULL'))) "
    "            WHEN region_stop::BIGINT = region_start::BIGINT "
    "              THEN error(printf('region_presence: region is zero-width "
    "(region_id ''%s'': [%s, %s)). Zero-width regions arise from bin_start when "
    "n_bins > genome_length; filter them with WHERE region_stop > region_start', "
    "                                region_id::VARCHAR, region_start::VARCHAR, "
    "                                region_stop::VARCHAR)) "
    "            WHEN region_stop::BIGINT < region_start::BIGINT "
    "              THEN error(printf('region_presence: region coordinates are inverted "
    "(region_id ''%s'': [%s, %s)); region_stop must be greater than region_start', "
    "                                region_id::VARCHAR, region_start::VARCHAR, "
    "                                region_stop::VARCHAR)) "
    "            ELSE TRUE "
    "          END "
    "), "
    // See note 1 and note 2 above: one CASE, and MATERIALIZED.
    "_rp_checked AS MATERIALIZED ( "
    "    SELECT sample_id, genome_id, start::BIGINT AS start, stop::BIGINT AS stop "
    "    FROM query_table(positions) "
    "    WHERE CASE "
    "            WHEN sample_id IS NULL OR genome_id IS NULL "
    "              THEN error('region_presence: NULL values in sample_id or genome_id "
    "column of positions (an interval that cannot be attributed would silently become a "
    "non-detection)') "
    "            ELSE start::BIGINT < stop::BIGINT "
    "          END "
    "), "
    // Drop non-cohort samples once, here, rather than scoring them against every region
    // and discarding them at the final join.
    "_rp_cov AS ( "
    "    SELECT c.* FROM _rp_checked c "
    "    SEMI JOIN _rp_roster s ON c.sample_id = s.sample_id "
    "), "
    "_rp_genome_cov AS (SELECT DISTINCT sample_id, genome_id FROM _rp_cov), "
    "_rp_hit AS ( "
    "    SELECT DISTINCT c.sample_id, r.genome_id, r.region_start, r.region_stop, "
    "           r.region_id "
    "    FROM _rp_cov c "
    "    JOIN _rp_reg r ON c.genome_id = r.genome_id "
    "                  AND c.start < r.region_stop AND c.stop > r.region_start "
    ") "
    "SELECT s.sample_id, r.genome_id, r.region_id, r.region_start, r.region_stop, "
    "       CASE WHEN h.sample_id IS NOT NULL THEN 'present' "
    "            WHEN g.sample_id IS NOT NULL THEN 'absent' "
    "            ELSE 'not applicable' END AS state "
    "FROM _rp_roster s "
    "CROSS JOIN _rp_reg r "
    "LEFT JOIN _rp_hit h ON h.sample_id = s.sample_id "
    "                   AND h.genome_id = r.genome_id "
    "                   AND h.region_start = r.region_start "
    "                   AND h.region_stop = r.region_stop "
    "                   AND h.region_id = r.region_id "
    "LEFT JOIN _rp_genome_cov g ON g.sample_id = s.sample_id "
    "                          AND g.genome_id = r.genome_id;";

// region_coverage(positions, regions)
//
// Breadth of coverage of a sub-genome region, per sample, with the REGION's length as
// the denominator. genome_coverage divides by the genome, so a fully covered 5 kb region
// inside a 5 Mb genome reports 0.1% there and 100% here. That denominator is the whole
// point of the function (issue #217).
//
// Parameters:
// positions : a relation with columns sample_id (any type), genome_id (any type),
//     start (BIGINT), stop (BIGINT). 1-based half-open intervals.
// regions : a relation with columns genome_id (any type), region_start (BIGINT),
//     region_stop (BIGINT), region_id (any type). Half-open on the same convention.
//
// Returns: sample_id, genome_id, region_id, region_start, region_stop, covered (BIGINT),
//     region_length (BIGINT), proportion_covered (DOUBLE)
//
// The column contract is region_presence's, NOT genome_coverage's. The two region_*
// functions are called on the same relation in the same pipeline -- typically
// region_coverage for the number and region_presence for whether the number means
// anything -- and making one of them demand `reference`/`position`/`stop_position` while
// its sibling demands `genome_id`/`start`/`stop` would force a rename between two
// consecutive lines of a query. genome_coverage keeps the alignment names because it
// predates the sample dimension and consumes read_alignments output directly; a caller
// coming from there renames once, in a view, and both region functions then take it.
// Uncompressed alignments are safe input either way -- compress_intervals runs inside.
//
// Clipping is `GREATEST(start, region_start)` / `LEAST(stop, region_stop)` on the rows
// the overlap test admits, and it happens before the merge. That ordering is NOT
// load-bearing, contrary to the claim in #217 that unclipped alignments overlapping
// outside the region would merge into one interval spanning it: intersection distributes
// over union, so (A u B) n R = (A n R) u (B n R) for intervals, and both orders return
// the same number on a fixture built to trigger it (test/sql/region_coverage.test,
// Test 5, asserts the equality against merge-then-clip written out in SQL). Clipping
// first is still the right order because it bounds what the merge has to sort.
//
// The MERGE, by contrast, is load-bearing. Two alignments overlapping INSIDE the region
// cover their union, not the sum of their lengths, and an unmerged sum inflates the
// proportion -- past 1.0 on a short region, invisibly on a long one.
//
// No roster, and no row for a (sample, region) pair with no overlap. That is deliberate:
// a fabricated `covered = 0` asserts the region was measured and found empty, which is
// true of a sample that covers the genome elsewhere and false of a sample where the
// organism is simply not present -- the exact conflation region_presence exists to
// prevent, and it cannot be undone downstream. Callers who need the full grid LEFT JOIN
// this onto region_presence, which supplies both the roster and the distinction.
//
// Every coordinate is normalized ::BIGINT before it is compared, not merely before it is
// subtracted. VARCHAR coordinates -- what a BED file or an inferred CSV column gives you
// -- otherwise compare lexicographically ('9' < '60' is FALSE), which drops real
// intervals from the overlap test and clips the survivors to the wrong end. Both failures
// produce a smaller, entirely plausible coverage number. INTEGER coordinates, which
// read_gff and read_ncbi_annotation emit, are widened for the same reason.
//
// covered is cast to BIGINT. DuckDB widens SUM to INT128 regardless of magnitude; the
// cast is always safe because covered is bounded by region_length. genome_coverage had
// the same widening and was returning HUGEINT against docs that promised BIGINT, so it
// now casts too -- the two are consistent.
//
// The two structural details from region_presence apply unchanged, and both were
// re-verified here by mutation rather than assumed:
//
// 1. The positions guard and the degeneracy filter are ONE CASE behind a MATERIALIZED
//    barrier. Split into two predicates, DuckDB collapses them and evaluates the cheap
//    interval test first, so a row that is both unattributable and degenerate is
//    discarded before the error branch is reached. Without MATERIALIZED, the region
//    join reorders ahead of the guard and an unattributable row on a genome no region
//    mentions is pruned before it can raise -- measured, and note it takes only three
//    well-formed rows alongside it for DuckDB to pick that plan, which is why the guard
//    fixtures in the test are deliberately not two-row VALUES lists.
// 2. CTE names are prefixed. An argument name is resolved with the macro's own WITH
//    list in scope, so an argument collides with an internal CTE if and only if that
//    CTE is declared BEFORE the query_table() call that resolves the argument. The
//    failure is a bare "Referenced column not found" naming neither the macro nor the
//    collision. Whether the caller's relation is a CTE, a view or a base table makes
//    no difference -- all three collide; an earlier version of this comment claimed
//    base tables were immune, which is simply false.
//
//    For this macro that is exactly ONE slot: _rc_reg precedes query_table(positions),
//    so a `positions` argument named _rc_reg collides. _rc_checked contains the
//    query_table(positions) call itself and resolves outward, and _rc_clipped and
//    _rc_merged are declared after both query_table() calls -- none of those three can
//    collide with anything, prefixed or not. So the prefix earns its keep on one name.
//    It makes the collision improbable, not impossible: `_rc_reg` in the positions slot
//    still breaks, and test/sql/region_coverage.test Test 16 asserts that it does.
//
//    The corollary for anyone adding a CTE here: putting it above _rc_checked adds a
//    new colliding name, putting it below does not.
//
// The DISTINCT on _rc_reg is NOT in that category: it bounds the join fan-out when a
// caller passes a regions relation consolidated from shards, but duplicate region rows
// cannot reach the output regardless, because the final GROUP BY absorbs them. It is
// load-bearing in region_presence, where a CROSS JOIN and no final aggregate mean a
// duplicated region row duplicates the answer.
//
// As in region_presence, the guards are still ordinary SQL expressions: a caller's
// filter may be pushed beneath them, so validation covers the rows the query reads,
// not the whole relation.
const std::string REGION_COVERAGE = // NOLINT
    "CREATE OR REPLACE MACRO region_coverage(positions, regions) AS TABLE "
    "WITH _rc_reg AS ( "
    "    SELECT DISTINCT genome_id, region_start::BIGINT AS region_start, "
    "           region_stop::BIGINT AS region_stop, region_id "
    "    FROM query_table(regions) "
    // A region is the QUESTION. A NULL bound makes the overlap test NULL, and a
    // zero-width region has no denominator at all -- answering either with a number
    // would be inventing a measurement. Zero-width and inverted are separate branches
    // because the causes and the remedies differ: filtering is right for a bin
    // expansion and wrong for transposed coordinates.
    "    WHERE CASE "
    "            WHEN genome_id IS NULL OR region_start IS NULL "
    "                 OR region_stop IS NULL OR region_id IS NULL "
    "              THEN error(printf('region_coverage: genome_id, region_start, "
    "region_stop and region_id must not be NULL (got genome_id=%s, region_start=%s, "
    "region_stop=%s, region_id=%s)', "
    "                                COALESCE(genome_id::VARCHAR, 'NULL'), "
    "                                COALESCE(region_start::VARCHAR, 'NULL'), "
    "                                COALESCE(region_stop::VARCHAR, 'NULL'), "
    "                                COALESCE(region_id::VARCHAR, 'NULL'))) "
    "            WHEN region_stop::BIGINT = region_start::BIGINT "
    "              THEN error(printf('region_coverage: region is zero-width "
    "(region_id ''%s'': [%s, %s)); it has no denominator. Zero-width regions arise from "
    "bin_start when n_bins > genome_length; filter them with "
    "WHERE region_stop > region_start', "
    "                                region_id::VARCHAR, region_start::VARCHAR, "
    "                                region_stop::VARCHAR)) "
    "            WHEN region_stop::BIGINT < region_start::BIGINT "
    "              THEN error(printf('region_coverage: region coordinates are inverted "
    "(region_id ''%s'': [%s, %s)); region_stop must be greater than region_start', "
    "                                region_id::VARCHAR, region_start::VARCHAR, "
    "                                region_stop::VARCHAR)) "
    "            ELSE TRUE "
    "          END "
    "), "
    // ONE CASE, and MATERIALIZED -- see the header note.
    "_rc_checked AS MATERIALIZED ( "
    "    SELECT sample_id, genome_id, start::BIGINT AS start, stop::BIGINT AS stop "
    "    FROM query_table(positions) "
    "    WHERE CASE "
    "            WHEN sample_id IS NULL OR genome_id IS NULL "
    "              THEN error('region_coverage: NULL values in sample_id or genome_id "
    "column of positions (an interval that cannot be attributed would silently vanish "
    "from the coverage it belongs to)') "
    "            ELSE start::BIGINT < stop::BIGINT "
    "          END "
    "), "
    "_rc_clipped AS ( "
    "    SELECT c.sample_id, r.genome_id, r.region_id, r.region_start, r.region_stop, "
    "           GREATEST(c.start, r.region_start) AS start, "
    "           LEAST(c.stop, r.region_stop) AS stop "
    "    FROM _rc_checked c "
    "    JOIN _rc_reg r ON c.genome_id = r.genome_id "
    "                  AND c.start < r.region_stop AND c.stop > r.region_start "
    "), "
    // Summed inside the aggregate, NOT by UNNESTing and re-grouping. The earlier version
    // exploded compress_intervals() into one row per merged interval and then re-hashed the
    // SAME five grouping keys to SUM them back up. The keys do not change across the UNNEST,
    // so the second hash aggregate was pure overhead. genome_coverage keeps its UNNEST
    // legitimately, because its grouping keys DO change across it for the contig->genome
    // join; these do not.
    //
    // Measured on 1M real alignments x 500 regions -> 3000 groups, threads=1, results
    // verified byte-identical (0 mismatches, EXCEPT both ways): the aggregate step alone
    // 0.076 s -> 0.067 s, end to end 0.686 s -> 0.651 s. So this is a ~10% step win, not a
    // step-change -- on that fixture compress_intervals merged nothing, so the UNNEST
    // re-grouped 1M rows and DuckDB sums 1M rows into 3000 groups cheaply. The win scales
    // with how little the intervals merge; correctness of the simplification does not.
    //
    // covered is computed in its own CTE so the aggregate is written -- and evaluated --
    // once, rather than repeated in the proportion_covered expression.
    "_rc_covered AS ( "
    "    SELECT sample_id, genome_id, region_id, region_start, region_stop, "
    "           list_sum(list_transform(compress_intervals(start, stop), "
    "                                   lambda iv: iv.stop - iv.start))::BIGINT AS covered "
    "    FROM _rc_clipped "
    "    GROUP BY sample_id, genome_id, region_id, region_start, region_stop "
    ") "
    "SELECT sample_id, genome_id, region_id, region_start, region_stop, covered, "
    "       (region_stop - region_start)::BIGINT AS region_length, "
    "       covered::DOUBLE / (region_stop - region_start) AS proportion_covered "
    "FROM _rc_covered;";

// cumulative_coverage_curve(positions, roster, genome_length)
//
// The rank-ordered cumulative breadth curve of issue #214, wrapped so the caller does
// not have to assemble the ranking, the zero-coverage backfill and the aggregate by
// hand. micov's framing is long-exposure astrophotography: samples that individually
// cover little of a genome stack into a detectable signal, which is what exposes a
// region present in one sample group and absent from another.
//
// Parameters:
// positions : a relation with columns sample_id (any type), start (BIGINT), stop
//     (BIGINT), PRE-FILTERED to a single genome. There is no genome_id parameter
//     because the curve is only meaningful against one target -- pooling genomes would
//     rank samples by their summed breadth across unrelated references.
// roster : a relation with columns sample_id (any type), group_id (any type). The
//     cohort, and the grouping the curves are computed within. Required for the same
//     reason region_presence requires one: a sample with no coverage of the target
//     contributes no position rows, so without the roster it would silently vanish and
//     the group would look smaller than it is.
// genome_length : the breadth denominator. Scalar, must be positive.
//
// Returns: group_id, rank (INTEGER, 0-based), sample_id, covered (BIGINT),
//     proportion_covered (DOUBLE)
//
// The ranking is a window function here, not a parameter of the aggregate: samples are
// ordered by their OWN breadth ascending, ties broken by sample_id. The tiebreak is
// load-bearing rather than tidy -- without it two equal-breadth samples swap ranks
// between runs and thread counts, and micov's Monte Carlo null calls this ~100x per
// genome per group, so an unstable curve would surface as noise in the null.
//
// Zero-coverage samples reach the aggregate as a row with NULL start/stop, which is how
// cumulative_coverage is told "this rank exists and covers nothing". That comes free
// from the roster LEFT JOIN, and it puts those samples at the low ranks, so the curve
// correctly starts flat.
//
// The aggregate returns only (rank, covered); sample_id is reattached by joining the
// curve back to the ranking on (group_id, rank), where it costs nothing. That is why
// the aggregate needs no id-type mirroring in its bind, and why VARCHAR, BIGINT and
// UUID sample identifiers all survive uncast.
//
// Conventions carried from region_presence/region_coverage, with one correction each
// time a claim was actually measured rather than assumed:
//
// - Every coordinate is ::BIGINT-normalized before comparison, so VARCHAR coordinates
//   from a BED or inferred CSV column are not compared as text.
// - Internal CTE names are prefixed. An argument name is resolved with the macro's own
//   WITH list in scope and binds to an internal CTE only when that CTE is declared
//   BEFORE the query_table() call resolving the argument. Measured for this macro: the
//   prefix on _cc_checked is load-bearing for the `roster` slot; the prefixes on
//   _cc_own/_cc_ranked/_cc_feed/_cc_curve protect nothing, because those are declared
//   after both query_table() calls. Declaration ORDER does more here than the prefix
//   does -- see the note on _cc_len below.
// - MATERIALIZED on _cc_checked is belt-and-braces, and the comment says so because the
//   measurement did not support more. Unlike region_coverage, where a pruning join sits
//   ahead of the guard and dropping the barrier lets a malformed row escape, here the
//   guard fires either way and EXPLAIN shows a byte-identical plan with and without it:
//   DuckDB already materializes a CTE referenced twice. It is kept for the explicit
//   single-scan guarantee, which matters for a relation that can only be consumed once
//   such as an Arrow RecordBatchReader, and for consistency with the sibling macros --
//   not because it is what makes validation work.
const std::string CUMULATIVE_COVERAGE_CURVE = // NOLINT
    "CREATE OR REPLACE MACRO cumulative_coverage_curve(positions, roster, genome_length) AS TABLE "
    "WITH _cc_checked AS MATERIALIZED ( "
    "    SELECT sample_id, start::BIGINT AS start, stop::BIGINT AS stop "
    "    FROM query_table(positions) "
    // The `start < stop` clause is load-bearing and was missing in the first version.
    // Without it an inverted interval reaches compress_intervals in _cc_own, which
    // SWAPS it -- so a sample with transposed coordinates was ranked as the best-covered
    // in its group while the aggregate (which drops inverted intervals) credited it with
    // nothing. Two conventions in one pipeline: the ranking said 90 bases, the curve said
    // 0, and both the curve shape and the per-rank sample_id were wrong. Filtering here
    // makes both halves agree, and the _cc_feed LEFT JOIN then supplies the NULL row that
    // registers the sample as zero-coverage, so it keeps its place at rank 0.
    //
    // A half-NULL interval raises rather than being swept into the zero-coverage case:
    // both-NULL is a sample with no coverage, but exactly one NULL is a broken join or a
    // blank field, and treating it as "covers nothing" would silently discard real
    // coverage.
    "    WHERE CASE "
    "            WHEN sample_id IS NULL "
    "              THEN error('cumulative_coverage_curve: NULL values in sample_id column of "
    "positions (coverage that cannot be attributed to a sample would silently lower "
    "somebody''s breadth and reorder the curve)') "
    "            WHEN (start IS NULL) <> (stop IS NULL) "
    "              THEN error(printf('cumulative_coverage_curve: start and stop must both be "
    "NULL or both be present (sample_id ''%s'' has start=%s, stop=%s). One NULL is a broken "
    "join or a blank field, not a zero-coverage sample.', "
    "                                COALESCE(sample_id::VARCHAR, 'NULL'), "
    "                                COALESCE(start::VARCHAR, 'NULL'), "
    "                                COALESCE(stop::VARCHAR, 'NULL'))) "
    "            WHEN start IS NULL THEN FALSE "
    "            ELSE start::BIGINT < stop::BIGINT "
    "          END "
    "), "
    // The DISTINCT bounds the join fan-out for a roster consolidated from shards; it is
    // NOT required for correctness, because _cc_own's GROUP BY absorbs duplicate roster
    // rows either way (measured). Same situation as the DISTINCT in region_coverage.
    "_cc_roster AS ( "
    "    SELECT DISTINCT sample_id, group_id FROM query_table(roster) "
    "    WHERE CASE "
    "            WHEN sample_id IS NULL OR group_id IS NULL "
    "              THEN error('cumulative_coverage_curve: NULL values in sample_id or group_id "
    "column of roster (a cohort member that belongs to no group cannot be ranked)') "
    "            ELSE TRUE "
    "          END "
    "), "
    // Each sample's OWN breadth, which is what the ranking sorts on. LEFT JOIN so a
    // roster sample with no coverage survives with 0 rather than disappearing;
    // compress_intervals returns NULL for an all-NULL group, hence the COALESCE.
    "_cc_own AS ( "
    "    SELECT s.group_id, s.sample_id, "
    "           COALESCE(list_sum(list_transform( "
    "               compress_intervals(c.start, c.stop), lambda iv: iv.stop - iv.start)), 0) "
    "             AS own_covered "
    "    FROM _cc_roster s "
    "    LEFT JOIN _cc_checked c ON c.sample_id = s.sample_id "
    "    GROUP BY s.group_id, s.sample_id "
    "), "
    "_cc_ranked AS ( "
    "    SELECT group_id, sample_id, "
    "           (ROW_NUMBER() OVER (PARTITION BY group_id "
    "                               ORDER BY own_covered, sample_id) - 1)::INTEGER AS rank "
    "    FROM _cc_own "
    "), "
    // One row per (rank, interval); a sample with no intervals yields exactly one row
    // with NULL start/stop, which registers the rank with zero coverage. That is also
    // what guarantees the aggregate sees a contiguous 0..n-1 rank set.
    "_cc_feed AS ( "
    "    SELECT k.group_id, k.rank, c.start, c.stop "
    "    FROM _cc_ranked k "
    "    LEFT JOIN _cc_checked c ON c.sample_id = k.sample_id "
    "), "
    "_cc_curve AS ( "
    "    SELECT group_id, UNNEST(cumulative_coverage(rank, start, stop)) AS pt "
    "    FROM _cc_feed "
    "    GROUP BY group_id "
    "), "
    // Declared LAST on purpose. An argument name is resolved with the macro's own WITH
    // list in scope, and it binds to an internal CTE only if that CTE was declared
    // BEFORE the query_table() call resolving the argument. Sitting first, this CTE
    // shadowed BOTH arguments; sitting last it shadows neither, which leaves _cc_checked
    // as the only name a caller relation can collide with (and only in the roster slot).
    // Verified by passing caller relations named after each internal CTE, in both slots,
    // as base tables and as CTEs.
    // Every test is against the DOUBLE value, before any narrowing. Casting to BIGINT
    // first rounds, so 1.6 became 2 and silently divided by the wrong denominator, while
    // 0.4 was rejected as "not positive" -- which it is. A genome length is a count of
    // bases, so a fractional value is a mistake worth naming rather than rounding.
    "_cc_len AS ( "
    "    SELECT CASE "
    "             WHEN genome_length IS NULL "
    "               THEN error('cumulative_coverage_curve: genome_length must be positive "
    "(got NULL)') "
    "             WHEN genome_length::DOUBLE <> FLOOR(genome_length::DOUBLE) "
    "               THEN error(printf('cumulative_coverage_curve: genome_length must be a whole "
    "number of bases (got %s)', genome_length::VARCHAR)) "
    "             WHEN genome_length::DOUBLE < 1 "
    "               THEN error(printf('cumulative_coverage_curve: genome_length must be positive "
    "(got %s)', genome_length::VARCHAR)) "
    "             ELSE genome_length::BIGINT "
    "           END AS genome_length "
    ") "
    // covered > genome_length is impossible for a correct denominator -- the union of
    // intervals on a genome cannot exceed its length -- so it means genome_length is
    // wrong (a contig length where a genome length belongs, or a mismatched reference
    // build) or positions carries coordinates off the end. Unclamped this returned
    // proportion_covered = 5.0 with no signal, which is worse than useless because the
    // docs teach that a proportion above 1.0 is the symptom of the summing mistake this
    // function exists to avoid.
    "SELECT v.group_id, v.pt.rank AS rank, k.sample_id, "
    "       v.pt.covered AS covered, "
    "       CASE WHEN v.pt.covered > l.genome_length "
    "              THEN error(printf('cumulative_coverage_curve: covered (%s) exceeds "
    "genome_length (%s) at rank %s -- genome_length is too small for these positions "
    "(a contig length used where a genome length belongs, or a mismatched reference)', "
    "                                v.pt.covered::VARCHAR, l.genome_length::VARCHAR, "
    "                                v.pt.rank::VARCHAR)) "
    "            ELSE v.pt.covered::DOUBLE / l.genome_length END AS proportion_covered "
    "FROM _cc_curve v "
    "JOIN _cc_ranked k ON k.group_id = v.group_id AND k.rank = v.pt.rank "
    "CROSS JOIN _cc_len l;";

// circular_query_coverage(alignments, reference_lengths)
//
// How much of each read is explained by each reference, pooling every alignment record
// the aligner split the read into.
//
// cigar_query_coverage answers a per-row question -- "how much of the read does THIS
// record explain" -- and cannot see sibling records. A read that spans the origin of a
// circular reference held as a linearised contig is emitted as two or more records, so
// no single row reports more than about half the read and a query-coverage floor
// discards it silently. That is a systematic loss localised at the origin of every
// assembled circular genome and plasmid.
//
// Coverage here is the UNION of the fragments' query footprints, not their sum. Summing
// overshoots: a junction read with a couple of bases deleted across the origin produces
// fragments that overlap by a base or two. The union is bounded by 1.0 for free, is
// unaffected by a read wrapping the reference more than once, and on a group of one
// fragment equals cigar_query_coverage exactly -- so a caller can use this
// unconditionally rather than branching on whether a reference is circular.
//
// Grouping is (read_id, read1 bit, reference). The read1 bit is load-bearing: R1 and R2
// share a read_id but are different molecules, and pooling them would manufacture
// coverage. Same reason woltka_ogu partitions on alignment_is_read1.
//
// Secondary records are excluded -- they re-place query bases a sibling already covers,
// so they cannot raise a union but would inflate n_fragments. Unmapped records and
// references absent from reference_lengths produce no rows.
//
// A read explained by two references yields one row per reference, each reporting how much
// of the read that reference explains. That is the question being answered, but it means the
// rows do not partition the read: summing coverage across references counts it more than
// once. Recruitment gates one reference at a time and is unaffected; abundance is not, and
// wants woltka_ogu, which distributes multi-mapped reads fractionally.
//
// Parameters:
// alignments : a relation with columns: read_id, flags, reference, position (BIGINT),
//     stop_position (BIGINT), cigar (VARCHAR)
// reference_lengths : a relation with columns: reference, length (BIGINT),
//     is_circular (BOOLEAN). read_alignment_header() supplies the first two; is_circular
//     has to be added, because circularity is not recorded in an alignment and cannot be
//     inferred from one. A reference left NULL is rejected rather than assumed circular.
// coverage_type : 'aligned' (default) counts M/=/X; 'mapped' also counts insertions.
//     Same vocabulary as cigar_query_coverage and cigar_query_intervals, so a pipeline
//     already filtering on cigar_query_coverage(cigar, 'mapped') can migrate unchanged.
//
// Returns: read_id, is_read1 (BOOLEAN), reference, coverage (DOUBLE), identity (DOUBLE),
//     query_length (BIGINT), n_fragments (BIGINT), mixed_strand (BOOLEAN),
//     max_ref_gap (BIGINT)
//
// mixed_strand and max_ref_gap are topology evidence, REPORTED not enforced: two fragments
// of one read can reach coverage 1.0 and identity 1.0 while not being an origin span at all
// (inverted repeat, chimera, inversion, misassembly), and identity does not catch those --
// a chimera the aligner splits scores 1.0 on every fragment. max_ref_gap is the reference
// gap modulo the reference length between query-adjacent same-strand fragments; a genuine
// wrap closes to 0, because linearising a circle is the only reason the aligner could not
// chain the fragments itself. Testing abutment against the contig ends instead would be
// wrong twice over -- see the docs.
//
// identity is NULL for legacy M CIGARs. That is a trap worth knowing before gating on it;
// the consequence is spelled out at the identity expression below and in the docs.
//
// The user-facing reference -- worked examples, the recommended gate, the measured gap
// separation, and the full list of preconditions -- lives in
// docs/alignment_analysis.md#circular-query-coverage. Keep behavioural changes in step with
// it; the comments from here down explain the implementation, not the contract.
const std::string CIRCULAR_QUERY_COVERAGE = // NOLINT
    "CREATE OR REPLACE MACRO circular_query_coverage(alignments, reference_lengths, "
    "                                                coverage_type := 'aligned') AS TABLE "
    // One row per reference, validated. Collapsing first matters because the documented
    // way to build this relation -- read_alignment_header() over a glob, or a UNION ALL of
    // headers -- yields one row per contig PER FILE. Joining those duplicates would double
    // n_fragments and fabricate a nonzero max_ref_gap for an ordinary non-wrapping read,
    // dropping it from the very gate below. A length that is missing or non-positive is
    // rejected rather than tolerated: `x % 0` and `x % NULL` are both NULL in DuckDB, which
    // would silently make max_ref_gap NULL, and COALESCE(max_ref_gap, 0) then reads
    // "unmeasured" as "perfectly abutting" -- turning the one column that distinguishes an
    // origin span from a chimera into a false accept.
    "WITH ref_len AS ( "
    "    SELECT "
    "        reference, "
    "        CASE WHEN MIN(length) <> MAX(length) "
    "             THEN error('circular_query_coverage: reference ' || reference::VARCHAR "
    "                        || ' has more than one recorded length in reference_lengths') "
    "             WHEN MAX(length) IS NULL OR MAX(length) <= 0 "
    "             THEN error('circular_query_coverage: reference ' || reference::VARCHAR "
    "                        || ' has a missing or non-positive length; a reference length is ' "
    "                        || 'required to recognise wrapping') "
    "             ELSE MAX(length) END AS ref_length, "
    // Circularity is the premise the whole gap measure rests on and it cannot be inferred
    // from an alignment: the same two fragments, one ending where the contig ends and the
    // next starting at its origin, are a wrap on a circular reference and an end-join on a
    // linear one. Assuming circular is what makes an adapter chimera or an assembly-graph
    // artifact look like a perfect origin span, and the assumption would be invisible at the
    // call site, so an undeclared reference is rejected rather than guessed at.
    // BOOL_OR <> BOOL_AND is "both values occur", the same idiom mixed_strand uses; both
    // ignore NULLs, so an all-NULL group falls through to the second branch.
    "        CASE WHEN BOOL_OR(is_circular) <> BOOL_AND(is_circular) "
    "             THEN error('circular_query_coverage: reference ' || reference::VARCHAR "
    "                        || ' is recorded as both circular and linear in reference_lengths') "
    "             WHEN BOOL_AND(is_circular) IS NULL "
    "             THEN error('circular_query_coverage: reference ' || reference::VARCHAR "
    "                        || ' does not say whether it is circular; set is_circular in ' "
    "                        || 'reference_lengths, because a reference gap that wraps is only ' "
    "                        || 'evidence of an origin span if the reference actually is circular') "
    "             ELSE BOOL_AND(is_circular) END AS is_circular "
    "    FROM query_table(reference_lengths) "
    "    WHERE reference IS NOT NULL "
    "    GROUP BY reference "
    "), "
    "frag AS ( "
    "    SELECT "
    "        a.read_id AS read_id, "
    "        a.reference AS reference, "
    "        r.ref_length AS ref_length, "
    "        r.is_circular AS is_circular, "
    // flags is cast because alignment relations that have been through Parquet widen it
    // to BIGINT, and every alignment_is_* function is USMALLINT-only.
    "        alignment_is_read1(a.flags::USMALLINT) AS is_read1, "
    "        alignment_is_reverse(a.flags::USMALLINT) AS is_reverse, "
    // Carried as-is: every consumer takes a difference of these two, and the window below
    // only orders by them, so converting to 0-based would cancel out.
    "        a.position AS ref_start, "
    "        a.stop_position AS ref_stop, "
    "        a.cigar AS cigar, "
    "        cigar_query_length(a.cigar, true) AS query_length, "
    "        cigar_query_intervals(a.cigar, a.flags::USMALLINT, coverage_type) AS intervals "
    "    FROM query_table(alignments) a "
    "    JOIN ref_len r "
    "      ON a.reference = r.reference "
    "    WHERE NOT alignment_is_unmapped(a.flags::USMALLINT) "
    "      AND NOT alignment_is_secondary(a.flags::USMALLINT) "
    // A record whose read_id is NULL cannot be grouped: read_alignments maps the '*' QNAME
    // sentinel to NULL, and pooling every unnamed record together would merge unrelated
    // molecules. Excluded explicitly so it is a stated rule rather than an accident of
    // NULL never matching in the final equi-join.
    "      AND a.read_id IS NOT NULL "
    // A record with no reference coordinates cannot be placed, so it cannot take part in the
    // reference-contiguity test. Leaving it in is worse than dropping it: it would count
    // toward n_fragments and coverage while its gap came back NULL, and MAX() ignores NULLs,
    // so a chimera whose fragments sit 148 kb apart would report max_ref_gap NULL and the gate
    // documented above -- COALESCE(max_ref_gap, 0) <= 100 -- would read "unmeasured" as
    // "perfectly abutting" and admit it. That is the same false accept the ref_len validation
    // prevents, arriving from the alignments side. Dropping the record instead surfaces the
    // problem where a caller will see it: the read's coverage falls by whatever that record
    // would have explained, so it fails a coverage floor rather than passing a gap check it
    // was never measured against. Well-formed mapped records always carry both columns, so
    // this fires only on malformed input.
    "      AND a.position IS NOT NULL "
    "      AND a.stop_position IS NOT NULL "
    "), "
    // A mapped record that covers no query positions (a clip-only CIGAR) would
    // contribute to n_fragments while contributing nothing to the union.
    "covering AS ( "
    "    SELECT * FROM frag WHERE len(intervals) > 0 "
    "), "
    // Fragment granularity. Ordered by where each fragment starts on the READ, which is
    // not record order -- the aligner is free to make any fragment the primary. On a
    // reverse-strand fragment travel along the read runs backwards along the reference,
    // so the successor's reference END is what the current fragment's START must meet.
    // A strand change is not a gap to measure but a different event entirely, so it
    // yields NULL here and is reported through mixed_strand instead.
    "chain AS ( "
    "    SELECT "
    "        read_id, is_read1, reference, query_length, is_reverse, "
    "        CASE WHEN LEAD(is_reverse) OVER w <> is_reverse THEN NULL "
    "             WHEN is_reverse THEN ref_start - LEAD(ref_stop) OVER w "
    "             ELSE LEAD(ref_start) OVER w - ref_stop END AS ref_delta, "
    // On a circular reference, reduce the signed delta into [0, ref_length) and take its
    // distance to zero, which is the absolute signed representative without forming it. A
    // delta of -ref_length -- one fragment ending where the contig ends and the next starting
    // at its origin -- lands on 0, which is what makes wrapping indistinguishable from
    // contiguity. On a linear reference that identification is exactly what must not happen,
    // so the plain distance is reported and an end-join shows up as the whole contig.
    "        (ref_delta % ref_length + ref_length) % ref_length AS delta_mod, "
    "        CASE WHEN is_circular THEN LEAST(delta_mod, ref_length - delta_mod) "
    "             ELSE abs(ref_delta) END AS ref_gap "
    "    FROM covering "
    // ref_start/ref_stop/is_reverse are tie-breakers, not decoration: two same-strand
    // fragments of one read can begin at the same query offset (a tandem-repeat
    // supplementary, or a duplicated record from merged BAMs). Ordering on the query start
    // alone leaves LEAD to pick an arbitrary peer, so max_ref_gap would vary with physical
    // scan order and the gate below would admit or reject the same read depending on it.
    "    WINDOW w AS (PARTITION BY read_id, is_read1, reference "
    "                 ORDER BY intervals[1].start, ref_start, ref_stop, is_reverse) "
    "), "
    // Union of the fragments' read-axis footprints, via the same compress_intervals
    // aggregate genome_coverage uses on the reference axis -- one definition of "union of
    // half-open intervals" for the whole extension, rather than a second one in SQL that
    // could drift from IntervalCompressor. A hand-rolled window sweep was tried here to
    // avoid per-group aggregate state; benchmarked at 500k and 2M read groups the
    // aggregate is equal or faster (0.224s vs 0.277s at 2M) with identical results, so the
    // second implementation bought nothing.
    "covered AS ( "
    "    SELECT read_id, is_read1, reference, SUM(ci.stop - ci.start) AS covered_bases "
    "    FROM ( "
    "        SELECT read_id, is_read1, reference, "
    "               UNNEST(compress_intervals(x.start, x.stop)) AS ci "
    "        FROM (SELECT read_id, is_read1, reference, UNNEST(intervals) AS x FROM covering) "
    "        GROUP BY read_id, is_read1, reference "
    "    ) "
    "    GROUP BY read_id, is_read1, reference "
    "), "
    // Identity is aggregated straight off `covering` rather than alongside the other evidence
    // columns, which read from `chain`. Everything selected in `chain` travels through the
    // window operator's sort payload -- it does not project down to the columns the window
    // functions actually use -- so leaving the CIGAR there copies every record's string
    // through the sort and back for the sake of an aggregate that needs no ordering at all.
    // Measured over 2M records in 1M read groups, hoisting it cuts this section from 0.37s to
    // 0.13s wall and drops the sort's system time by an order of magnitude. It cannot ride
    // along in `covered` either: that CTE unnests intervals, so a fragment contributing two
    // intervals would have its CIGAR counted twice.
    "pooled_identity AS ( "
    "    SELECT read_id, is_read1, reference, "
    // sum(=) / sum(alignment columns) over the fragments, which on a single fragment is
    // exactly cigar_sequence_identity of that fragment. Note it does NOT charge the junction
    // as a gap -- the fragments are one molecule, and the reference discontinuity between
    // them is what max_ref_gap is for.
    //
    // NULL for legacy M CIGARs, where identity is not recoverable from the CIGAR at all,
    // matching cigar_sequence_identity's own behaviour on the same input, and NULL when the
    // fragments disagree about the encoding. Raising instead was tried and rejected: making
    // the raise conditional on the column being selected depends on projection pruning, and
    // pruning is an optimisation -- with the optimiser disabled the expression is evaluated
    // and the error fires even for a query that only asked for coverage. Semantics must not
    // vary with the optimiser, so this stays NULL and the precondition is documented in the
    // header comment instead.
    "           cigar_pooled_identity(cigar) AS identity "
    "    FROM covering "
    "    GROUP BY read_id, is_read1, reference "
    "), "
    "evidence AS ( "
    "    SELECT "
    "        read_id, is_read1, reference, "
    // Every record of one read reports the same full query length, because clipping is
    // what the aligner uses to account for the bases a fragment does not align.
    // Disagreement means the relation grouped rows from more than one read under a
    // single read_id, which is a caller error rather than a coverage figure to guess at.
    // MIN <> MAX rather than COUNT(DISTINCT ...): a distinct aggregate builds a second
    // hash table keyed on (group, value) over every row, which over millions of read-level
    // groups measured ~30% of the whole macro's CPU. Both ignore NULLs, so "more than one
    // distinct non-null value" is exactly MIN <> MAX, and an all-NULL group falls through
    // to MAX (NULL) either way.
    "        CASE WHEN MIN(query_length) <> MAX(query_length) "
    "             THEN error('circular_query_coverage: fragments under read_id ' "
    "                        || read_id::VARCHAR || ' report different query lengths; the ' "
    "                        || 'alignments relation groups rows from more than one read') "
    "             ELSE MAX(query_length) END AS query_length, "
    "        COUNT(*) AS n_fragments, "
    // "some query-adjacent pair differs in strand" is just "both strand values occur in the
    // group" -- a transition has to happen somewhere -- so this needs no ordering and no
    // window at all.
    "        BOOL_OR(is_reverse) <> BOOL_AND(is_reverse) AS mixed_strand, "
    "        MAX(ref_gap) AS max_ref_gap "
    "    FROM chain "
    "    GROUP BY read_id, is_read1, reference "
    ") "
    "SELECT "
    "    e.read_id, "
    "    e.is_read1, "
    "    e.reference, "
    "    c.covered_bases::DOUBLE / e.query_length AS coverage, "
    "    p.identity, "
    "    e.query_length, "
    "    e.n_fragments, "
    "    e.mixed_strand, "
    "    e.max_ref_gap "
    "FROM evidence e "
    "JOIN covered c USING (read_id, is_read1, reference) "
    "JOIN pooled_identity p USING (read_id, is_read1, reference);";

// infer_trim(original_reads, qcd_reads)
//
// Recover per-read 5'/3' trim coordinates from reads that were quality-
// controlled by an *external* tool that only trims read ends (no base edits)
// and may omit reads entirely. LEFT JOINs original -> QC'd on sequence_index
// and locates the (contiguous) QC'd sequence within the original.
//
// Both relations must expose:
//   sequence_index : BIGINT  (join key)
//   sequence       : VARCHAR (read sequence)
//
// CONTRACT / preconditions (documented, not all enforced):
//   sequence_index must be GLOBALLY UNIQUE per read and identify the SAME read
//   on both sides -- it is the caller's job to carry a stable key through the
//   external tool. read_fastx assigns sequence_index positionally and RESETS it
//   per input file, so it is a valid key only for a single file whose numbering
//   the external tool round-trips; you cannot re-derive a matching key by
//   re-running read_fastx on the QC'd output (QC drops reads, so the positional
//   numbering diverges). A non-unique sequence_index on the qcd side fans the
//   LEFT JOIN out (>1 row per original); this precondition is documented rather
//   than enforced, to keep the macro a single linear scan. A qcd row whose
//   sequence_index is absent from the originals is dropped (one row per
//   original), by design.
//
// Returns one row per original read:
//   sequence_index : BIGINT
//   trimmed_5p     : UINTEGER  bases removed from the 5' end (NULL if omitted)
//   trimmed_3p     : UINTEGER  bases removed from the 3' end (NULL if omitted)
//
// A read with no QC'd counterpart yields NULL/NULL (omitted). Per-row data-
// integrity violations FAIL LOUD (throw) regardless of which output columns the
// caller projects, because the checks live in a row-level WHERE filter the
// optimizer cannot prune away: a present QC row with a NULL or empty sequence,
// a NULL original sequence, or a QC sequence that is not a contiguous substring
// of its original (base edits, or a mismatched join key). This macro is strictly
// for pure end-trimming. When the kept block occurs at multiple offsets
// (self-repetitive reads), the leftmost match wins. position() is computed once
// in a CTE, not per output column. A qcd `matched` marker distinguishes an
// omitted read (no row) from a present row with a NULL sequence.
const std::string INFER_TRIM = // NOLINT
    "CREATE OR REPLACE MACRO infer_trim(original_reads, qcd_reads) AS TABLE "
    "WITH joined AS ( "
    "    SELECT o.sequence_index AS sequence_index, "
    "           o.sequence       AS oseq, "
    "           q.sequence       AS qseq, "
    "           q.matched        AS matched "
    "    FROM query_table(original_reads) o "
    "    LEFT JOIN (SELECT sequence_index, sequence, TRUE AS matched "
    "               FROM query_table(qcd_reads)) q USING (sequence_index) "
    "), "
    "computed AS ( "
    "    SELECT sequence_index, matched, oseq, qseq, "
    "           CASE WHEN matched IS NULL THEN NULL ELSE position(qseq IN oseq) END AS p, "
    "           length(oseq) AS olen, "
    "           length(qseq) AS qlen "
    "    FROM joined "
    "), "
    "validated AS ( "
    "    SELECT * FROM computed "
    "    WHERE CASE "
    "        WHEN matched IS NULL THEN TRUE "
    "        WHEN oseq IS NULL THEN error(printf("
    "            'infer_trim: original sequence for sequence_index=%d is NULL', sequence_index)) "
    "        WHEN qseq IS NULL THEN error(printf("
    "            'infer_trim: QC sequence for sequence_index=%d is NULL (a kept read must have a sequence)', "
    "            sequence_index)) "
    "        WHEN qseq = '' THEN error(printf("
    "            'infer_trim: QC sequence for sequence_index=%d is empty (omit the row to mark a dropped read)', "
    "            sequence_index)) "
    "        WHEN p = 0 THEN error(printf("
    "            'infer_trim: QC sequence for sequence_index=%d is not a contiguous substring of its original', "
    "            sequence_index)) "
    "        ELSE TRUE END "
    ") "
    "SELECT sequence_index, "
    "       CASE WHEN matched IS NULL THEN NULL ELSE (p - 1)::UINTEGER END AS trimmed_5p, "
    "       CASE WHEN matched IS NULL THEN NULL ELSE (olen - (p - 1) - qlen)::UINTEGER END AS trimmed_3p "
    "FROM validated; ";

// taxonomy_lineage(taxids := NULL, source := NULL, refresh := false)
//
// Offline counterpart to read_ncbi_lineage: walks an NCBI taxdump tree (read via
// read_ncbi_taxdump(source, refresh)) from each query taxon up to the root and
// rank-collapses the path into the SHARED lineage schema, so the two functions are
// drop-in interchangeable. `taxids` is a BIGINT[] of query taxa; NULL (default)
// collapses every live node. `source`/`refresh` pass straight through to
// read_ncbi_taxdump (NULL source => auto-download + cache the canonical NCBI dump).
//
// Rank -> column mapping is the single documented source of truth mirrored from
// src/taxonomy_lineage.cpp: legacy 'superkingdom' and newer 'domain' both collapse
// to `domain`; phylum/class/order/family/genus/species/strain map to like-named
// columns; every other rank (clade, subspecies, no rank, ...) is dropped. The
// formatted `lineage` renders d__;p__;c__;o__;f__;g__;s__ (empty where absent) and
// appends ;t__<strain> only when a strain is present. Absent ranks are NULL in
// their column but empty in the formatted string, matching read_ncbi_lineage.
//
// `tree` is MATERIALIZED so the (expensive) taxdump parse runs once, not once per
// recursive iteration. A requested taxid absent from the tree yields no row (same
// omission semantics as read_ncbi_lineage).
const std::string TAXONOMY_LINEAGE = // NOLINT
    "CREATE OR REPLACE MACRO taxonomy_lineage(taxids := NULL, source := NULL, refresh := false) AS TABLE "
    "WITH RECURSIVE tree AS MATERIALIZED ( "
    "    SELECT node_index, parent_index, name, rank "
    "    FROM read_ncbi_taxdump(source, refresh := refresh) "
    "), "
    "seeds(q_taxid) AS ( "
    "    SELECT node_index FROM tree WHERE taxids IS NULL "
    "    UNION ALL "
    // A NULL element in the taxids list is a user error; reject it row-wise via
    // error() rather than silently dropping it (parity with read_ncbi_lineage,
    // which throws "taxid list cannot contain NULL"). The check lives in the
    // CASE THEN branch so it fires only for a NULL element, never when the list
    // is clean or when taxids itself is NULL (then UNNEST yields no rows).
    "    SELECT DISTINCT CASE WHEN v IS NULL "
    "                         THEN error('taxonomy_lineage: taxids list cannot contain NULL') "
    "                         ELSE v END "
    "    FROM (SELECT UNNEST(taxids)::BIGINT AS v) "
    "), "
    "walk(q_taxid, cur, w_name, w_rank, w_parent) AS ( "
    "    SELECT s.q_taxid, t.node_index, t.name, t.rank, t.parent_index "
    "    FROM seeds s JOIN tree t ON t.node_index = s.q_taxid "
    "    UNION ALL "
    "    SELECT w.q_taxid, t.node_index, t.name, t.rank, t.parent_index "
    "    FROM walk w JOIN tree t ON t.node_index = w.w_parent "
    "), "
    "collapsed AS ( "
    "    SELECT "
    "        q_taxid AS taxid, "
    // name/rank of the query taxon itself: always emitted (never NULLed), even when
    // empty, matching read_ncbi_lineage which always writes name/rank cells.
    "        MAX(CASE WHEN cur = q_taxid THEN w_name END) AS name, "
    "        MAX(CASE WHEN cur = q_taxid THEN w_rank END) AS rank, "
    // Rank-collapse columns: NULLIF(..., '') so a node with no scientific name yields
    // NULL, not '', matching read_ncbi_lineage's SetRankCell (empty -> NULL column).
    // The formatted lineage below COALESCEs these back to '' so its segments are
    // unaffected.
    "        NULLIF(MAX(CASE WHEN w_rank IN ('superkingdom', 'domain') THEN w_name END), '') AS domain, "
    "        NULLIF(MAX(CASE WHEN w_rank = 'phylum'  THEN w_name END), '') AS phylum, "
    "        NULLIF(MAX(CASE WHEN w_rank = 'class'   THEN w_name END), '') AS \"class\", "
    "        NULLIF(MAX(CASE WHEN w_rank = 'order'   THEN w_name END), '') AS \"order\", "
    "        NULLIF(MAX(CASE WHEN w_rank = 'family'  THEN w_name END), '') AS family, "
    "        NULLIF(MAX(CASE WHEN w_rank = 'genus'   THEN w_name END), '') AS genus, "
    "        NULLIF(MAX(CASE WHEN w_rank = 'species' THEN w_name END), '') AS species, "
    "        NULLIF(MAX(CASE WHEN w_rank = 'strain'  THEN w_name END), '') AS strain "
    "    FROM walk "
    "    GROUP BY q_taxid "
    ") "
    "SELECT "
    "    taxid, name, rank, "
    "    domain, phylum, \"class\", \"order\", family, genus, species, strain, "
    "    'd__' || COALESCE(domain, '') || "
    "    ';p__' || COALESCE(phylum, '') || "
    "    ';c__' || COALESCE(\"class\", '') || "
    "    ';o__' || COALESCE(\"order\", '') || "
    "    ';f__' || COALESCE(family, '') || "
    "    ';g__' || COALESCE(genus, '') || "
    "    ';s__' || COALESCE(species, '') || "
    "    CASE WHEN strain IS NOT NULL AND strain <> '' THEN ';t__' || strain ELSE '' END "
    "    AS lineage "
    "FROM collapsed; ";

// beta_group_distances(distances, groups)
//
// Labels each pair in a condensed distance table (e.g. unifrac_distances output:
// sample_a, sample_b, distance) as within- or between-group by joining each side
// to a (sample_id, grouping) relation. Returns the per-pair rows; aggregate with
// GROUP BY comparison + quantile_cont/avg for the within/between distribution.
//
// Preconditions (documented, not enforced — pass clean inputs):
//   - `distances` holds each unordered pair once. If it carries multiple
//     `iteration`s (unifrac_distances with n_subsamples > 1), pre-filter to a
//     single iteration first, or the distribution pools across bootstrap
//     replicates. query_table() rejects subqueries, so filter into a table/view.
//   - `sample_id` is unique in `groups`; a duplicate fans the two joins out and
//     silently inflates the counts.
//   - A pair whose either endpoint is absent from `groups` is dropped (inner join).
//   - NULL grouping: `NULL = NULL` is NULL in SQL, so a pair of two NULL-group
//     samples is labeled 'between', the same as genuinely different groups.
//   - `groups` is referenced twice; pass a materialized table, not a view over an
//     expensive query.
const std::string BETA_GROUP_DISTANCES = // NOLINT
    "CREATE OR REPLACE MACRO beta_group_distances(distances, groups) AS TABLE "
    "SELECT d.sample_a, d.sample_b, d.distance, "
    "       ga.grouping AS group_a, gb.grouping AS group_b, "
    "       CASE WHEN ga.grouping = gb.grouping THEN 'within' ELSE 'between' END AS comparison "
    "FROM query_table(distances) d "
    "JOIN query_table(groups) ga ON d.sample_a = ga.sample_id "
    "JOIN query_table(groups) gb ON d.sample_b = gb.sample_id; ";

// beta_knn(distances, k)
//
// The k nearest neighbors of every sample over a condensed distance table. The
// input holds each unordered pair once, so both orientations are unioned before
// ranking. Returns (sample_id, neighbor, distance, rank), rank in 1..k. Ties are
// broken by neighbor id (deterministic); k <= 0 yields no rows.
//
// If `distances` carries multiple `iteration`s (unifrac_distances with
// n_subsamples > 1), pre-filter to a single iteration first, or neighbors are
// ranked across all iterations together. `distances` is referenced twice, so pass
// a materialized table, not a view over an expensive query.
const std::string BETA_KNN = // NOLINT
    "CREATE OR REPLACE MACRO beta_knn(distances, k) AS TABLE "
    "SELECT sample_id, neighbor, distance, rank FROM ( "
    "    SELECT sample_id, neighbor, distance, "
    "           ROW_NUMBER() OVER (PARTITION BY sample_id ORDER BY distance, neighbor) AS rank "
    "    FROM ( "
    "        SELECT sample_a AS sample_id, sample_b AS neighbor, distance FROM query_table(distances) "
    "        UNION ALL "
    "        SELECT sample_b AS sample_id, sample_a AS neighbor, distance FROM query_table(distances) "
    "    ) "
    ") WHERE rank <= k; ";

// beta_knn_from_sample(distances, k, source)
//
// The k samples nearest a single `source` sample over a condensed distance table,
// checking both orientations. Returns (neighbor, distance) ordered nearest-first;
// ties broken by neighbor id. An absent source yields no rows. Same
// iteration/materialization caveats as beta_knn.
const std::string BETA_KNN_FROM_SAMPLE = // NOLINT
    "CREATE OR REPLACE MACRO beta_knn_from_sample(distances, k, source) AS TABLE "
    "SELECT neighbor, distance FROM ( "
    "    SELECT sample_b AS neighbor, distance FROM query_table(distances) WHERE sample_a = source "
    "    UNION ALL "
    "    SELECT sample_a AS neighbor, distance FROM query_table(distances) WHERE sample_b = source "
    ") ORDER BY distance, neighbor LIMIT k; ";

// mmvec_train_test_split(relation, test_fraction, seed)
//
// Split a long-form feature-table's SAMPLES into 'train' and 'test', for holding
// data back from mmvec_fit and scoring on it with mmvec_score. Returns one row per
// distinct sample_id: (sample_id, split), sample_id passed through with its own
// type so it joins back to the feature table without a cast.
//
// ONE relation, not two. mmvec_fit requires its X and Y tables to describe exactly
// the same samples and validates that itself, so a second relation would carry no
// information -- split either one and filter both by the result.
//
// Exactly round(n * test_fraction) samples are assigned 'test' (rounding half away
// from zero, so 10 samples at 0.35 gives 4). test_fraction outside [0, 1] is an
// error rather than a silent all-train or all-test.
//
// NULL is rejected before the range test, and has to be: SQL's three-valued logic
// makes `NULL < 0 OR NULL > 1` evaluate to NULL rather than true, so a NULL would
// fall past a range test written on its own. It would then make n_test NULL, make
// `rn <= n_test` NULL, and land every sample in the ELSE branch -- an all-train
// split, silently, which is the exact outcome the range check exists to prevent.
// A NULL seed is rejected for the same class of reason: md5(NULL || ...) is NULL
// for every row, so the ordering would collapse to plain alphabetical by sample_id
// and the split would look seeded without being seeded.
//
// The assignment is a deterministic function of (sample_id, seed) alone: samples
// are ordered by md5(seed || ':' || sample_id) and the first n_test taken. Ties are
// broken by sample_id, so the result depends on neither the input row order nor how
// the scan was parallelized, and re-running it -- in another session, on a permuted
// table -- gives the same split. md5 rather than DuckDB's hash() deliberately: hash()
// is an implementation detail that may change between versions, whereas md5 is
// specified, so a split recorded in a paper still reproduces after an upgrade. Not a
// cryptographic claim -- this is a permutation, not a secret.
//
// A sample-wise split routinely leaves test-only FEATURES in the held-out tables,
// which mmvec_predict and mmvec_score reject by design (there is no conditional
// probability for a feature the model never saw). Restrict the test tables to the
// model's own features first -- their error messages name the exact predicate.
const std::string MMVEC_TRAIN_TEST_SPLIT = // NOLINT
    "CREATE OR REPLACE MACRO mmvec_train_test_split(relation, test_fraction, seed) AS TABLE "
    "SELECT sample_id, "
    "       CASE WHEN test_fraction IS NULL "
    "              THEN error('mmvec_train_test_split: test_fraction must not be NULL') "
    "            WHEN seed IS NULL "
    "              THEN error('mmvec_train_test_split: seed must not be NULL') "
    "            WHEN test_fraction < 0 OR test_fraction > 1 "
    "              THEN error(printf('mmvec_train_test_split: test_fraction must be in [0, 1], got %s', "
    "                                CAST(test_fraction AS VARCHAR))) "
    "            WHEN rn <= n_test THEN 'test' ELSE 'train' END AS split "
    "FROM ( "
    "    SELECT sample_id, "
    "           ROW_NUMBER() OVER (ORDER BY md5(CAST(seed AS VARCHAR) || ':' || CAST(sample_id AS VARCHAR)), "
    "                                       sample_id) AS rn, "
    "           CAST(round(COUNT(*) OVER () * test_fraction) AS BIGINT) AS n_test "
    "    FROM (SELECT DISTINCT sample_id FROM query_table(relation)) "
    "); ";

class MIINTMacros {
public:
	static void Register(ExtensionLoader &loader) {
		// Register macros into the system catalog (like DuckDB's built-in macros)
		// rather than executing CREATE statements against the default catalog.
		// The default catalog may be a database the user attached read-only, in
		// which case a CREATE would abort extension load with an internal error.
		// The system catalog is always writable and resolves across all catalogs.
		auto register_macro = [&](const std::string &sql, const char *name) {
			Parser parser;
			parser.ParseQuery(sql);
			if (parser.statements.size() != 1 || parser.statements[0]->type != StatementType::CREATE_STATEMENT) {
				throw InternalException("Failed to register macro '%s': expected a single CREATE statement", name);
			}
			auto &create = parser.statements[0]->Cast<CreateStatement>();
			if (create.info->type != CatalogType::MACRO_ENTRY && create.info->type != CatalogType::TABLE_MACRO_ENTRY) {
				throw InternalException("Failed to register macro '%s': not a macro definition", name);
			}
			auto &macro_info = create.info->Cast<CreateMacroInfo>();
			// Pin to the system catalog's default schema. RegisterFunction targets
			// the system catalog directly; the parsed statement leaves the schema
			// empty, which GetSchema would otherwise fail to resolve there.
			macro_info.schema = DEFAULT_SCHEMA;
			// The system catalog only accepts internal entries (same as built-in
			// macros and RegisterType). Mark accordingly.
			macro_info.internal = true;
			macro_info.temporary = true;
			loader.RegisterFunction(macro_info);
		};

		register_macro(MIINT_WARNINGS, "miint_warnings");

		register_macro(PARSE_GFF_ATTRIBUTES, "parse_gff_attributes");
		register_macro(READ_GFF, "read_gff");
		register_macro(GENOME_COVERAGE, "genome_coverage");
		register_macro(GENOME_COVERAGE_PER_SAMPLE, "genome_coverage_per_sample");
		register_macro(BIN_OF, "bin_of");
		register_macro(BIN_START, "bin_start");
		register_macro(INTERVAL_BINS, "interval_bins");
		register_macro(REGION_PRESENCE, "region_presence");
		register_macro(REGION_COVERAGE, "region_coverage");
		register_macro(CUMULATIVE_COVERAGE_CURVE, "cumulative_coverage_curve");
		register_macro(CIRCULAR_QUERY_COVERAGE, "circular_query_coverage");
		register_macro(INFER_TRIM, "infer_trim");

		register_macro(READ_JPLACE, "read_jplace");

		register_macro(TAXONOMY_LINEAGE, "taxonomy_lineage");

		register_macro(MZ_WITHIN, "mz_within");
		register_macro(MZ_WITHIN_PPM, "mz_within_ppm");
		register_macro(MZML_PEAKS, "mzml_peaks");
		register_macro(MZML_SCANINFO, "mzml_scaninfo");
		register_macro(MZML_SCANSUM, "mzml_scansum");
		register_macro(MZML_SCANNUM, "mzml_scannum");
		register_macro(MZML_SCANMZ, "mzml_scanmz");
		register_macro(MZML_SCANMAXINT, "mzml_scanmaxint");
		register_macro(MZML_MS1_PEAKS, "mzml_ms1_peaks");
		register_macro(MZML_MS2_PEAKS, "mzml_ms2_peaks");
		register_macro(MZML_MS1_PARENT_PEAKS, "mzml_ms1_parent_peaks");
		register_macro(MZML_MS2_CHILD_PEAKS, "mzml_ms2_child_peaks");
		register_macro(MZML_MS1_WHERE_MS2PROD, "mzml_ms1_where_ms2prod");
		register_macro(MZML_MS2_WHERE_MS1MZ, "mzml_ms2_where_ms1mz");
		register_macro(MZML_MS1_WHERE_MS2PREC, "mzml_ms1_where_ms2prec");
		register_macro(MZML_MS2_WHERE_MS2PROD_AND_MS1MZ, "mzml_ms2_where_ms2prod_and_ms1mz");
		register_macro(MZML_FILTER_MZ, "mzml_filter_mz");
		register_macro(MZML_FILTER_NL, "mzml_filter_nl");
		register_macro(MASSDEFECT, "massdefect");
		register_macro(MZ_MASSDEFECT_WITHIN, "mz_massdefect_within");
		// ntuple must be registered before pair and triplet (they delegate to it)
		register_macro(MZML_X_OFFSET_NTUPLE, "mzml_x_offset_ntuple");
		register_macro(MZML_X_OFFSET_PAIR, "mzml_x_offset_pair");
		register_macro(MZML_X_OFFSET_TRIPLET, "mzml_x_offset_triplet");
		register_macro(MZML_X_PREC_PROD, "mzml_x_prec_prod");
		register_macro(MZML_X_PREC_MASSDEFECT, "mzml_x_prec_massdefect");
		register_macro(MZML_X_MS1_MS2_PREC, "mzml_x_ms1_ms2_prec");
		register_macro(MZML_X_OFFSET_PAIR_RANGE, "mzml_x_offset_pair_range");
		register_macro(MZML_OR_CARDINALITY, "mzml_or_cardinality");
		// mzml_peak_pair is now a C++ table function (MzmlPeakPairFunction) for performance.
		register_macro(MZML_I_NORM, "mzml_i_norm");
		register_macro(MZML_I_TIC_NORM, "mzml_i_tic_norm");
		register_macro(MZML_EXCLUDED_MS2PROD, "mzml_excluded_ms2prod");
		register_macro(MZML_EXCLUDED_MS1MZ, "mzml_excluded_ms1mz");
		register_macro(MZML_EXCLUDED_MS2PREC, "mzml_excluded_ms2prec");
		register_macro(MZML_ISOTOPE_PATTERN, "mzml_isotope_pattern");

		register_macro(BETA_GROUP_DISTANCES, "beta_group_distances");
		register_macro(BETA_KNN, "beta_knn");
		register_macro(BETA_KNN_FROM_SAMPLE, "beta_knn_from_sample");

		register_macro(MMVEC_TRAIN_TEST_SPLIT, "mmvec_train_test_split");
	}
};

} // namespace duckdb
