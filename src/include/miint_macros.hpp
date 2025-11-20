#pragma once

#include "duckdb.hpp"

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

class MIINTMacros {

public:
	static void Register(ExtensionLoader &loader) {
		auto &instance = loader.GetDatabaseInstance();
		Connection con(instance);

		con.Query(WOLTKA_OGU_PER_SAMPLE);
		con.Query(WOLTKA_OGU);
		con.Query(PARSE_GFF_ATTRIBUTES);
		con.Query(READ_GFF);
	}
};

} // namespace duckdb
