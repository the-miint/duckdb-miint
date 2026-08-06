#pragma once

#include <map>
#include <string>
#include <vector>

#include "duckdb/main/client_context.hpp"

namespace duckdb::absquant_internal {

// A keyed reference relation read into a key vector plus one parallel vector
// per requested value column. All vectors have the same length; row `r` of the
// relation is `keys[r]` with `values.at("col")[r]`.
//
// Keyed by column NAME rather than by position on purpose. Every column here is
// a double, so unpacking them positionally -- `values[0][i], values[1][i]` --
// type-checks whether or not that order still matches the list of names handed
// to ReadKeyedColumns twenty lines earlier. Nothing would catch the two drifting
// apart: a transposed `slope`/`intercept` compiles, runs, and returns wrong
// numbers. Naming the columns at the point of use removes the possibility
// rather than relying on the two lists being maintained in step.
struct KeyedColumns {
	std::vector<std::string> keys;
	std::map<std::string, std::vector<double>> values;
};

// Read `(key_column, value_columns...)` from a user-named relation, following
// the separate-connection recipe in docs/internals/reading-tables-views.md.
//
// A NULL value becomes NaN, which is how SQL NULL travels into the pure cores:
// for a per-sample parameter that is a legitimate "skip this sample", and for
// everything else the core rejects it with a message naming the row. A NULL KEY
// is refused outright -- an unnamed row cannot be joined to anything, so
// dropping it silently would change a computation (the synDNA pool total, say)
// with nothing to show for it.
//
// (ReadFeatureTable, which reads the long-form counts relations, instead drops
// rows with a NULL id -- an established contract shared with unifrac_* and
// community_distances. The asymmetry is defensible on its own terms, a NULL cell
// in a sparse matrix being absence rather than a broken config, but it is worth
// revisiting across all of them.)
//
// `table_name` is quoted before interpolation and may come from the user.
// `key_column` and `value_columns` are NOT quoted and MUST be trusted
// compile-time literals -- wiring a user-supplied column name (a named
// parameter, say) into either is a SQL injection. Keep the column names baked in
// at the call site and select among them; never pass one through.
//
// `value_columns` must be distinct; a repeat is a caller bug and throws.
KeyedColumns ReadKeyedColumns(ClientContext &context, const std::string &table_name, const char *key_column,
                              const std::vector<const char *> &value_columns, const char *entity,
                              const std::string &caller_name);

// One row of a long-form `(sample_id, feature_id, <value>)` relation.
struct LongFormRow {
	std::string sample_id;
	std::string feature_id;
	double value = 0.0;
};

// Read a long-form per-(sample, feature) relation, PRESERVING zero-valued rows.
//
// That is the whole reason this exists next to ReadFeatureTable rather than
// reusing it. ReadFeatureTable drops zeros to maintain the sparse-storage
// invariant, which is right for counts -- a zero count contributes nothing and
// its absence means the same thing. It is wrong for coverage: zero is a real
// measurement there, and with `min_coverage := 0` pysyndna keeps a
// zero-coverage feature, so a reader that dropped the row would silently
// disagree about which cells exist.
//
// Rows with a NULL sample_id or feature_id are dropped, matching
// ReadFeatureTable; a NULL value becomes NaN and is left for the core to
// reject, since "this cell has no coverage" is exactly the situation that must
// not pass unremarked.
//
// `value_column` is NOT quoted and MUST be a trusted compile-time literal.
std::vector<LongFormRow> ReadLongFormValues(ClientContext &context, const std::string &table_name,
                                            const char *value_column, const char *entity,
                                            const std::string &caller_name);

} // namespace duckdb::absquant_internal
