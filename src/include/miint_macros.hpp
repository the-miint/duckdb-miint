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
// subject_total_length : a relation with columns: genome_id (VARCHAR),
//     total_length (BIGINT)
// subject_genome_id : a relation with columns: contig_id (VARCHAR),
//     genome_id (VARCHAR)
//
// Returns: genome_id (VARCHAR), covered (BIGINT), proportion_covered (DOUBLE)
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
    "    JOIN query_table(subject_genome_id) sg "
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
    "SELECT "
    "    tc.genome_id, "
    "    tc.covered, "
    "    tc.covered::DOUBLE / tl.total_length AS proportion_covered "
    "FROM total_coverage tc "
    "JOIN query_table(subject_total_length) tl "
    "  USING (genome_id);";

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
