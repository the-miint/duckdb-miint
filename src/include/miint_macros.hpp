#pragma once

#include "duckdb.hpp"
#include "duckdb/main/extension_helper.hpp"

namespace duckdb {

// TODO: overloading didn't work, unsure why, it's arguably better
// to be explicit anyway.
//
// TODO: there is a lot of common logic with the macros.

// woltka_ogu_per_sample(relation, sample_id_field, sequence_id_field)
//
// Compute Woltka OGU counts over SAM-like data for multiple samples.
//
// Parameters:
// relation : a relation, view, etc to SAM-like data. The relation must contain
//     what "sequence_id_field" and "sample_id_field" resolve to as provided
//     by the caller. The relation must also contain `subject_id`, and `flags`
// sample_id_field : a varchar, the field within the relation which contains
// 	   sample IDs.
// sequence_id_field : a VARCHAR, the field within the relation to group
// 	   sequence identifiers by. Note: this is abstracted out to allow users to
//     specify either the SAM standard `read_id` which is VARCHAR or to use
//     a numeric index if it exists. A numeric index will be faster for the
//     internal GROUP BY operations.
//
// IMPORTANT: queries should not be quoted. As in, an operation should look like
//   SELECT * FROM woltka_ogu_per_sample(some_table, some_field, some_other_field);
// A query should _not_ look like
// 	 SELECT * FROM woltka_ogu_per_sample('some_table', 'some_field', 'some_other_field');
// The specific issue is `some_field` and `some_other_field` will be re-interpreted as
// a string literal.
const std::string WOLTKA_OGU_PER_SAMPLE =
    "CREATE OR REPLACE MACRO woltka_ogu_per_sample(relation, sample_id_field, sequence_id_field) AS TABLE "
    "WITH "
    "    base AS ( "
    "        SELECT DISTINCT "
    "            sequence_id_field AS query_local_id_field, "
    "            sample_id_field AS query_local_sample_id, "
    "            reference AS feature_id, "
    "            alignment_is_read1(flags::USMALLINT) AS is_fwd "
    "        FROM query_table(relation) "
    "    ), "
    "    with_counts AS ( "
    "        SELECT "
    "            query_local_sample_id, "
    "            feature_id, "
    "            1.0 / COUNT(*) OVER (PARTITION BY query_local_id_field, is_fwd) AS local_value "
    "        FROM base "
    "    ) "
    "SELECT "
    "    query_local_sample_id AS sample_id, "
    "    feature_id, "
    "    SUM(local_value) AS value "
    "FROM with_counts "
    "GROUP BY query_local_sample_id, feature_id; ";

// woltka_ogu(relation, sequence_id_field)
//
// Compute Woltka OGU counts over SAM-like data irrespective of sample.
//
// Parameters:
// relation : a relation, view, etc to SAM-like data. The relation must contain
//     what "sequence_id_field" and "sample_id_field" resolve to as provided
//     by the caller. The relation must also contain `subject_id`, and `flags`
// sequence_id_field : a VARCHAR, the field within the relation to group
// 	   sequence identifiers by. Note: this is abstracted out to allow users to
//     specify either the SAM standard `read_id` which is VARCHAR or to use
//     a numeric index if it exists. A numeric index will be faster for the
//     internal GROUP BY operations.
//
// IMPORTANT: queries should not be quoted. As in, an operation should look like
//   SELECT * FROM woltka_ogu(some_table, some_field);
// A query should _not_ look like
// 	 SELECT * FROM woltka_ogu('some_table', 'some_field');
// The specific issue is `some_field` will be re-interpreted as a string literal.

const std::string WOLTKA_OGU = // NOLINT
    "CREATE OR REPLACE MACRO woltka_ogu(relation, sequence_id_field) AS TABLE "
    "WITH "
    "    base AS ( "
    "        SELECT DISTINCT "
    "            sequence_id_field AS query_local_id_field, "
    "            reference AS feature_id, "
    "            alignment_is_read1(flags::USMALLINT) AS is_fwd "
    "        FROM query_table(relation) "
    "    ), "
    "    with_counts AS ( "
    "        SELECT "
    "            feature_id, "
    "            1.0 / COUNT(*) OVER (PARTITION BY query_local_id_field, is_fwd) AS local_value "
    "        FROM base "
    "    ) "
    "SELECT "
    "    feature_id, "
    "    SUM(local_value) AS value "
    "FROM with_counts "
    "GROUP BY feature_id;";

const std::string PARSE_GFF_ATTRIBUTES = // NOLINT
    "CREATE OR REPLACE MACRO parse_gff_attributes(kvp_string) AS ( "
    "  map_from_entries( "
    "    list_transform( "
    "      string_split(kvp_string, ';'), "
    "      x -> struct_pack( "
    "        key := string_split(x, '=')[1], "
    "        value := string_split(x, '=')[2] "
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
    "   column4::INTEGER AS stop_position, "
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
    "WHERE column0 NOT LIKE '##%'; ";

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
    "    FROM read_json(path, filename := true) "
    "); ";

// mzml_peaks(filepath)
//
// Unnests mz_array and intensity_array from read_mzml into per-peak rows.
// Each output row is one (mz, intensity) peak with its parent spectrum's metadata.
// i_norm = intensity / base_peak_intensity (NULL when base_peak_intensity is NULL).
const std::string MZML_PEAKS = // NOLINT
    "CREATE OR REPLACE MACRO mzml_peaks(filepath) AS TABLE "
    "SELECT spectrum_index, ms_level, retention_time, spectrum_type, polarity, "
    "       base_peak_intensity, total_ion_current, "
    "       precursor_mz, precursor_charge, precursor_intensity, ms1_scan_index, "
    "       mz, intensity, intensity / base_peak_intensity AS i_norm "
    "FROM ( "
    "    SELECT spectrum_index, ms_level, retention_time, spectrum_type, polarity, "
    "           base_peak_intensity, total_ion_current, "
    "           precursor_mz, precursor_charge, precursor_intensity, ms1_scan_index, "
    "           UNNEST(mz_array) AS mz, UNNEST(intensity_array) AS intensity "
    "    FROM read_mzml(filepath) "
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

// mzml_i_norm(intensity_array, base_peak_intensity)
//
// Max-normalize an intensity array by dividing each value by the base peak intensity.
// Returns NULL if base_peak_intensity is NULL (natural NULL propagation).
const std::string MZML_I_NORM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_i_norm(intensity_array, base_peak_intensity) AS ("
    "  list_transform(intensity_array, x -> x / base_peak_intensity)"
    ");";

// mzml_i_tic_norm(intensity_array, total_ion_current)
//
// TIC-normalize an intensity array by dividing each value by the total ion current.
// Returns NULL if total_ion_current is NULL (natural NULL propagation).
const std::string MZML_I_TIC_NORM = // NOLINT
    "CREATE OR REPLACE MACRO mzml_i_tic_norm(intensity_array, total_ion_current) AS ("
    "  list_transform(intensity_array, x -> x / total_ion_current)"
    ");";

class MIINTMacros {

public:
	static void Register(ExtensionLoader &loader) {
		auto &instance = loader.GetDatabaseInstance();
		Connection con(instance);

		con.Query(WOLTKA_OGU_PER_SAMPLE);
		con.Query(WOLTKA_OGU);
		con.Query(PARSE_GFF_ATTRIBUTES);
		con.Query(READ_GFF);

		// read_jplace requires the json extension for read_json function
		ExtensionInstallOptions options;
		ExtensionHelper::InstallExtension(*con.context, "json", options);
		ExtensionHelper::AutoLoadExtension(instance, "json");
		con.Query(READ_JPLACE);

		con.Query(MZML_PEAKS);
		con.Query(MZML_SCANINFO);
		con.Query(MZML_I_NORM);
		con.Query(MZML_I_TIC_NORM);
	}
};

} // namespace duckdb
