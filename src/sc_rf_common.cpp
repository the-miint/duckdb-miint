#include "sc_rf_common.hpp"

#include <string>
#include <string_view>
#include <vector>

namespace duckdb {
namespace sc_rf {

namespace {

void ReleaseNothing(ArrowArray *array) {
	array->release = nullptr;
}
void ReleaseNothingSchema(ArrowSchema *schema) {
	schema->release = nullptr;
}

} // namespace

void BuildTargets(TargetArray &t, bool classification) {
	t.schema.format = classification ? "u" : "g";
	t.schema.name = nullptr;
	t.schema.metadata = nullptr;
	t.schema.flags = 0;
	t.schema.n_children = 0;
	t.schema.children = nullptr;
	t.schema.dictionary = nullptr;
	t.schema.release = ReleaseNothingSchema;
	t.schema.private_data = nullptr;

	int64_t length = 0;
	int64_t n_buffers = 0;
	if (classification) {
		length = static_cast<int64_t>(t.labels.size());
		t.offsets.reserve(t.labels.size() + 1);
		int32_t cursor = 0;
		t.offsets.push_back(cursor);
		for (const auto &s : t.labels) {
			t.chars.insert(t.chars.end(), s.begin(), s.end());
			cursor += static_cast<int32_t>(s.size());
			t.offsets.push_back(cursor);
		}
		t.buffers[1] = t.offsets.data();
		t.buffers[2] = t.chars.empty() ? nullptr : t.chars.data();
		n_buffers = 3;
	} else {
		length = static_cast<int64_t>(t.numbers.size());
		t.buffers[1] = t.numbers.data();
		n_buffers = 2;
	}
	// buffers[0] stays NULL: no nulls, which sc requires (null_count must be 0).
	t.array.length = length;
	t.array.null_count = 0;
	t.array.offset = 0;
	t.array.n_buffers = n_buffers;
	t.array.n_children = 0;
	t.array.buffers = t.buffers;
	t.array.children = nullptr;
	t.array.dictionary = nullptr;
	t.array.release = ReleaseNothing;
	t.array.private_data = nullptr;
}

[[noreturn]] void ThrowNotCooTriplet(const string &relation, const string &engine_error, const char *caller) {
	throw InvalidInputException("%s: Data relation '%s' does not match the required COO triplet schema.\n"
	                            "  Expected columns : sample_id, feature_id, value\n"
	                            "  Engine error     : %s\n\n"
	                            "Remedy:\n"
	                            "  If your relation uses different column names (e.g. ASV/taxa names, read counts),\n"
	                            "  wrap it in an aliased view:\n"
	                            "    CREATE VIEW my_counts AS\n"
	                            "      SELECT your_sample_col  AS sample_id,\n"
	                            "             your_feature_col AS feature_id,\n"
	                            "             your_count_col   AS value\n"
	                            "      FROM %s;",
	                            caller, relation, engine_error, relation);
}

//! Scan the data relation into `builder`, rejecting NULLs.
void ScanCounts(Connection &conn, const ScTrainingInput &bind, miint::CooBuilder &builder) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(bind.data_relation);
	// The casts are load-bearing, not cosmetic. Reading a vector's buffer
	// directly assumes its physical type, and DuckDB infers `42.0` as
	// DECIMAL(3,1) (physical INT16), an INTEGER count column as INT32, and
	// woltka's feature ids as BIGINT or UUID. Vector::GetValue used to convert
	// on the way out; a raw buffer read cannot. Casting in SQL makes DuckDB do
	// the conversion and guarantees the layout this loop reads.

	// given a table with triplet columns as input to the duck db table function sc_fit_*, we are invoking the duckdb
	// SQL engine to go fetch the data in columnar format, n output chunks of <= 2048 rows each
	auto result = conn.Query("SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + q);
	if (result->HasError()) {
		ThrowNotCooTriplet(bind.data_relation, result->GetError(), bind.caller);
	}
	// Read through UnifiedVectorFormat rather than Vector::GetValue(row).
	// GetValue materialises a duckdb::Value per cell -- a heap-allocating
	// variant -- and ToString() allocates again on top, so a scan of n cells
	// costs ~5n allocations that are discarded immediately. At 13M cells that
	// dominates the scan. This path reads string_t views straight out of
	// DuckDB's own buffers, and CooBuilder::Append takes string_view, so a
	// cell now costs no allocation at all beyond interning a genuinely new id.
	//
	// Unified (not FlatVector) because the column may arrive constant- or
	// dictionary-encoded, in which case the selection vector maps row -> slot.
	// chunk is the <= 2048 record /rows coming from DuckDB in one go
	// duckdb::unique_ptr<duckdb::DataChunk> chunk is owned by
	while (duckdb::unique_ptr<duckdb::DataChunk> chunk = result->Fetch()) {
		const idx_t n = chunk->size();
		// alocate a unified vector format for each column on the stack, and fill it with the data from the
		// corresponding chunk buffer.  This is a view into the chunk's data, not a copy.
		UnifiedVectorFormat sf, ff, vf;
		chunk->data[0].ToUnifiedFormat(n, sf);
		chunk->data[1].ToUnifiedFormat(n, ff);
		chunk->data[2].ToUnifiedFormat(n, vf);
		const auto *samples = UnifiedVectorFormat::GetData<string_t>(sf);
		const auto *features = UnifiedVectorFormat::GetData<string_t>(ff);
		const auto *values = UnifiedVectorFormat::GetData<double>(vf);

		for (idx_t row = 0; row < n; row++) {
			// Step 1: Translate logical row -> physical buffer index
			const auto si = sf.sel->get_index(row);
			const auto fi = ff.sel->get_index(row);
			const auto vi = vf.sel->get_index(row);
			// A NULL here means a broken join upstream, not a zero count. An
			// absent cell is already zero in a sparse matrix, so a NULL cannot
			// be passed through and guessing at it would hide the mistake.
			if (!sf.validity.RowIsValid(si) || !ff.validity.RowIsValid(fi) || !vf.validity.RowIsValid(vi)) {
				throw InvalidInputException(
				    "%s: NULL in data relation '%s' (sample_id/feature_id/value must all be non-NULL)", bind.caller,
				    bind.data_relation);
			}
			// Step 2: Use si to fetch the actual data from the raw buffer
			// assign the result to a readonly alias
			// bracket operator is dereference operator
			// prevent stack copy when dereferencing the pointer via & aliasing
			// & refers to the address of the object
			// allocation no new variable storage
			const string_t &s = samples[si];
			const string_t &f = features[fi];
			const double v = values[vi];
			//  append borrowed views of the sample_id, feature_id as string_view
			//  +-------------------+-------------------+
			//	|  const char* ptr  |    size_t len     |
			//	+-------------------+-------------------+
			//  16 bytes 8 byte pointer, a length number
			//  copy values as double, which is 8 bytes on most platforms
			builder.Append(std::string_view(s.GetData(), s.GetSize()), std::string_view(f.GetData(), f.GetSize()), v);
		}
	} // duckdb::unique_ptr<duckdb::DataChunk> chunk freed, by now relevant view data is copied on the heap for the
	  // CooBuilder
}

//! Reject duplicate cells, showing the values so the caller can tell a join
//! fanout (identical values, deduplicate) from repeat measurements (differing
//! values, maybe sum). Summing a fanout silently inflates every count, so the
//! message deliberately does not prescribe a repair.
void RequireNoDuplicateCells(const miint::CooBuilder &builder, const ScTrainingInput &bind) {
	const auto report = builder.FindDuplicateCells();
	if (report.Empty()) {
		return;
	}
	string msg = StringUtil::Format("%s: %llu duplicate (sample_id, feature_id) cell(s) in '%s'.", bind.caller,
	                                (unsigned long long)report.duplicate_cells, bind.data_relation.c_str());
	for (const auto &c : report.examples) {
		string values;
		for (size_t i = 0; i < c.values.size(); i++) {
			values += (i ? ", " : "") + StringUtil::Format("%g", c.values[i]);
		}
		msg += StringUtil::Format("\n  %s / %s  x%llu  values: %s", c.sample_id.c_str(), c.feature_id.c_str(),
		                          (unsigned long long)c.count, values.c_str());
	}
	msg += "\nIdentical values suggest a join fanout (deduplicate); differing values suggest repeat "
	       "measurements (aggregate deliberately).";
	throw InvalidInputException(msg);
}

//! Fill `params` with sklearn's defaults for the task.
//!
//! Classifier: gini, max_features='sqrt'. Regressor: squared_error,
//! max_features=1.0 (all). Both: 100 trees, unbounded depth, min_samples_split=2,
//! min_samples_leaf=1, bootstrap=True, max_samples=None.
void ApplyDefaults(sc_rf_params_t &params, bool classification) {
	params.criterion = classification ? SC_CRITERION_GINI : SC_CRITERION_SQUARED_ERROR;
	params.max_depth = 0; // sklearn's None
	params.max_features.kind = classification ? SC_MAX_FEATURES_SQRT : SC_MAX_FEATURES_ALL;
	params.min_samples_split.kind = SC_MIN_SAMPLES_COUNT;
	params.min_samples_split.count = 2;
	params.min_samples_leaf.kind = SC_MIN_SAMPLES_COUNT;
	params.min_samples_leaf.count = 1;
	params.min_weight_fraction_leaf = 0.0;
	params.min_impurity_decrease = 0.0;
	params.n_estimators = 100;
	params.bootstrap = true;
	params.max_samples.kind = SC_MAX_SAMPLES_ALL;
	params.random_state = 0;
	params.optimize_feature_selection = false;
	params.rfe_step = 0.0;
	params.parameter_tuning = false;
	params.cv = 0;
	params.n_threads = 0;
}

//! Is this value an SQL integer rather than a float?
//!
//! sklearn's int/float distinction is load-bearing for the tagged unions --
//! `max_features=10` means ten features, `max_features=0.3` means thirty
//! percent -- and SQL's literal types carry exactly that information, so the
//! named parameters take ANY and dispatch on what the user actually wrote.
bool IsIntegerValue(const Value &v) {
	switch (v.type().id()) {
	case LogicalTypeId::TINYINT:
	case LogicalTypeId::SMALLINT:
	case LogicalTypeId::INTEGER:
	case LogicalTypeId::BIGINT:
	case LogicalTypeId::HUGEINT:
	case LogicalTypeId::UTINYINT:
	case LogicalTypeId::USMALLINT:
	case LogicalTypeId::UINTEGER:
	case LogicalTypeId::UBIGINT:
		return true;
	default:
		return false;
	}
}

void ParseMaxFeatures(const Value &v, const char *caller, sc_max_features_t &out) {
	if (v.type().id() == LogicalTypeId::VARCHAR) {
		const auto s = v.GetValue<string>();
		if (StringUtil::CIEquals(s, "sqrt")) {
			out.kind = SC_MAX_FEATURES_SQRT;
		} else if (StringUtil::CIEquals(s, "log2")) {
			out.kind = SC_MAX_FEATURES_LOG2;
		} else if (StringUtil::CIEquals(s, "all") || StringUtil::CIEquals(s, "none")) {
			out.kind = SC_MAX_FEATURES_ALL;
		} else {
			throw InvalidInputException("%s: max_features must be 'sqrt', 'log2', 'all', a fraction in (0, 1], "
			                            "or a positive integer (got '%s')",
			                            caller, s);
		}
		return;
	}
	if (IsIntegerValue(v)) {
		const auto n = v.GetValue<int64_t>();
		if (n <= 0) {
			throw InvalidInputException("%s: max_features count must be > 0 (got %lld)", caller, (long long)n);
		}
		out.kind = SC_MAX_FEATURES_COUNT;
		out.count = static_cast<uint64_t>(n);
		return;
	}
	const auto f = v.GetValue<double>();
	if (!(f > 0.0 && f <= 1.0)) {
		throw InvalidInputException("%s: max_features fraction must be in (0, 1] (got %g)", caller, f);
	}
	out.kind = SC_MAX_FEATURES_FRACTION;
	out.fraction = f;
}

void ParseMinSamples(const Value &v, const char *name, uint64_t min_count, const char *caller, sc_min_samples_t &out) {
	if (IsIntegerValue(v)) {
		const auto n = v.GetValue<int64_t>();
		if (n < static_cast<int64_t>(min_count)) {
			throw InvalidInputException("%s: %s count must be >= %llu (got %lld)", caller, name,
			                            (unsigned long long)min_count, (long long)n);
		}
		out.kind = SC_MIN_SAMPLES_COUNT;
		out.count = static_cast<uint64_t>(n);
		return;
	}
	const auto f = v.GetValue<double>();
	if (!(f > 0.0 && f <= 1.0)) {
		throw InvalidInputException("%s: %s fraction must be in (0, 1] (got %g)", caller, name, f);
	}
	out.kind = SC_MIN_SAMPLES_FRACTION;
	out.fraction = f;
}

void ParseMaxSamples(const Value &v, const char *caller, sc_max_samples_t &out) {
	if (v.type().id() == LogicalTypeId::VARCHAR) {
		const auto s = v.GetValue<string>();
		if (!StringUtil::CIEquals(s, "all") && !StringUtil::CIEquals(s, "none")) {
			throw InvalidInputException("%s: max_samples must be 'all', a fraction in (0, 1], or a positive "
			                            "integer (got '%s')",
			                            caller, s);
		}
		out.kind = SC_MAX_SAMPLES_ALL;
		return;
	}
	if (IsIntegerValue(v)) {
		const auto n = v.GetValue<int64_t>();
		if (n <= 0) {
			throw InvalidInputException("%s: max_samples count must be > 0 (got %lld)", caller, (long long)n);
		}
		out.kind = SC_MAX_SAMPLES_COUNT;
		out.count = static_cast<uint64_t>(n);
		return;
	}
	const auto f = v.GetValue<double>();
	if (!(f > 0.0 && f <= 1.0)) {
		throw InvalidInputException("%s: max_samples fraction must be in (0, 1] (got %g)", caller, f);
	}
	out.kind = SC_MAX_SAMPLES_FRACTION;
	out.fraction = f;
}

void ParseCriterion(const Value &v, bool classification, const char *caller, sc_criterion_t &out) {
	const auto s = v.GetValue<string>();
	// sc rejects a criterion that does not match the task, but catching it here
	// names both the criterion and the function the user called.
	if (classification) {
		if (StringUtil::CIEquals(s, "gini")) {
			out = SC_CRITERION_GINI;
		} else if (StringUtil::CIEquals(s, "entropy")) {
			out = SC_CRITERION_ENTROPY;
		} else {
			throw InvalidInputException("%s_classifier: criterion must be 'gini' or 'entropy' (got '%s')", caller, s);
		}
		return;
	}
	if (StringUtil::CIEquals(s, "squared_error")) {
		out = SC_CRITERION_SQUARED_ERROR;
		return;
	}
	throw InvalidInputException("%s_regressor: criterion must be 'squared_error' (got '%s')", caller, s);
}

//! Work out which metadata column holds the target.
//!
//! With one candidate there is nothing to choose, so choose it -- `(sample_id,
//! month)` needs no configuration. With several the answer is genuinely unknown,
//! and guessing would train on the wrong variable without any error: a
//! `(sample_id, age, bmi, delivery_mode)` table has three equally plausible
//! targets and picking one silently would be the worst outcome.
//!
//! This replaces a hardcoded default of "value", which happened to suit
//! long-format metadata and quietly failed for everything else.
std::string ResolveTargetColumn(Connection &conn, const std::string &relation, const char *caller) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(relation);
	// LIMIT 0 binds the relation and returns its schema without scanning it.
	auto probe = conn.Query("SELECT * FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: metadata relation '%s' could not be read: %s", caller, relation,
		                            probe->GetError());
	}

	std::vector<std::string> candidates;
	bool has_sample_id = false;
	for (const auto &name : probe->names) {
		if (StringUtil::CIEquals(name, kSampleIdColumn)) {
			has_sample_id = true;
			continue; // if sample_id, dont push to candidates
		}
		candidates.push_back(name);
	}
	if (!has_sample_id) {
		throw InvalidInputException("%s: metadata relation '%s' has no 'sample_id' column", caller, relation);
	}
	if (candidates.empty()) {
		throw InvalidInputException(
		    "%s: metadata relation '%s' has only a sample_id column; it needs a target column too", caller, relation);
	}
	if (candidates.size() == 1) {
		return candidates[0];
	}

	std::string list;
	for (size_t i = 0; i < candidates.size(); i++) {
		list += (i ? ", " : "") + ("'" + candidates[i] + "'");
	}
	throw InvalidInputException(
	    "%s: Metadata relation '%s' has more than one column that could be the target.\n"
	    "  Candidates : [%s]\n\n"
	    "Remedy:\n"
	    "  Name the variable you are modelling, so the wrong one cannot be picked silently:\n"
	    "    SELECT * FROM %s_regressor('counts', '%s',\n"
	    "                                   target_column := 'one_of_the_above',\n"
	    "                                   name := 'my_model');\n"
	    "  A metadata relation with exactly one non-sample_id column needs no target_column at all.",
	    caller, relation, list, caller, relation);
}

namespace {

//! The declared type of one column, via a bind-only probe.
LogicalType ProbeColumnType(Connection &conn, const string &relation, const string &column, const char *caller) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(relation);
	const auto col = KeywordHelper::WriteOptionallyQuoted(column);
	// LIMIT 0 binds the relation and returns its schema without scanning it.
	auto probe = conn.Query("SELECT " + col + " FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: relation '%s' has no '%s' column: %s", caller, relation, column,
		                            probe->GetError());
	}
	return probe->types[0];
}

} // namespace

namespace {

//! Reject a column type that cannot be an id, naming the cast that fixes it.
//!
//! Ids are matched as text, so any castable type would *work*; the restriction
//! to VARCHAR / BIGINT / UUID is the codebase's (id_column_utils.hpp), kept so
//! an id column behaves the same across miint.
LogicalType RequireIdType(const LogicalType &type, const string &relation, const string &column, const char *caller) {
	if (!IsAllowedIdType(type)) {
		throw InvalidInputException("%s: '%s' in relation '%s' is %s; an id column must be %s.\n\n"
		                            "Remedy:\n"
		                            "  Cast it in a view, keeping the type you want returned:\n"
		                            "    CREATE VIEW typed AS SELECT %s::VARCHAR AS %s, * EXCLUDE (%s) FROM %s;",
		                            caller, column, relation, type.ToString(), AllowedIdTypeList(), column, column,
		                            column, relation);
	}
	return type;
}

} // namespace

ScCooIdTypes DetectCooIdTypes(Connection &conn, const string &relation, const char *caller) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(relation);
	// LIMIT 0 binds the relation and returns its schema without scanning it.
	auto probe = conn.Query("SELECT sample_id, feature_id, value FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		// The same message the scan would give, raised before any work happens.
		ThrowNotCooTriplet(relation, probe->GetError(), caller);
	}
	ScCooIdTypes out;
	out.sample_id_type = RequireIdType(probe->types[0], relation, kSampleIdColumn, caller);
	out.feature_id_type = RequireIdType(probe->types[1], relation, "feature_id", caller);
	return out;
}

LogicalType DetectColumnType(Connection &conn, const string &relation, const string &column, const char *caller) {
	return ProbeColumnType(conn, relation, column, caller);
}

} // namespace sc_rf
} // namespace duckdb
