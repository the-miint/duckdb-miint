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
    "    oriented AS ( "
    "        SELECT sequence_id_field AS query_local_id_field, "
    "               sample_id_field AS query_local_sample_id, "
    "               reference AS feature_id, "
    "               alignment_is_read1(flags::USMALLINT) AS is_fwd "
    "        FROM query_table(relation)), "
    "    fcount AS ( "
    "        SELECT query_local_id_field, "
    "               is_fwd, "
    "               1 / array_unique(array_agg(feature_id)) AS local_value "
    "        FROM oriented "
    "        GROUP By query_local_id_field, is_fwd), "
    "    coo AS ( "
    "        SELECT DISTINCT ON (s.query_local_id_field, "
    "                            s.feature_id, "
    "                            s.is_fwd) "
    "               query_local_sample_id, "
    "               feature_id, "
    "               fc.local_value, "
    "               s.is_fwd "
    "        FROM fcount fc "
    "            JOIN oriented s ON fc.query_local_id_field=s.query_local_id_field "
    "                AND fc.is_fwd=s.is_fwd) "
    "SELECT query_local_sample_id AS sample_id, "
    "       feature_id, "
    "       list_sum(array_agg(local_value)) AS value "
    "FROM coo "
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
    "    oriented AS ( "
    "        SELECT sequence_id_field AS query_local_id_field, "
    "               reference AS feature_id, "
    "               alignment_is_read1(flags::USMALLINT) AS is_fwd "
    "        FROM query_table(relation)), "
    "    fcount AS ( "
    "        SELECT query_local_id_field, "
    "               is_fwd, "
    "               1 / array_unique(array_agg(feature_id)) AS local_value "
    "        FROM oriented "
    "        GROUP By query_local_id_field, is_fwd), "
    "   coo AS ( "
    "       SELECT DISTINCT ON (s.query_local_id_field, "
    "                            s.feature_id, "
    "                            s.is_fwd) "
    "                  s.feature_id, "
    "                  fc.local_value, "
    "                  s.is_fwd "
    "       FROM fcount fc "
    "           JOIN oriented s ON fc.query_local_id_field=s.query_local_id_field "
    "               AND fc.is_fwd=s.is_fwd) "
    "SELECT feature_id, "
    "       list_sum(array_agg(local_value)) AS value "
    "FROM coo "
    "GROUP BY feature_id; ";

class MIINTMacros {

public:
	static void Register(ExtensionLoader &loader) {
		auto &instance = loader.GetDatabaseInstance();
		Connection con(instance);

		con.Query(WOLTKA_OGU_PER_SAMPLE);
		con.Query(WOLTKA_OGU);
	}
};

} // namespace duckdb
