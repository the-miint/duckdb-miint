#include "sc_common.hpp"

#include "id_column_utils.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <cstdlib>
#include <cstring>
#include <string>

namespace miint {

using duckdb::InternalException;
using duckdb::InvalidInputException;

ScContext::~ScContext() {
	if (ptr) {
		sc_context_free(ptr);
	}
}

ScModel::~ScModel() {
	if (ptr) {
		sc_model_free(ptr);
	}
}

OwnedArrowArray::~OwnedArrowArray() {
	// A released array has a NULL callback, so this is a no-op on an out-slot
	// that was never filled (an sc call that failed before exporting).
	if (array_.release) {
		array_.release(&array_);
	}
	if (schema_.release) {
		schema_.release(&schema_);
	}
}

namespace {

void RequireFormat(const ArrowSchema &schema, const char *expected, const char *what) {
	if (!schema.format || std::strcmp(schema.format, expected) != 0) {
		throw InvalidInputException("%s: expected Arrow format '%s', got '%s'", what, expected,
		                            schema.format ? schema.format : "(null)");
	}
}

} // namespace

// Read Arrow utf8 arrays (format "u")
std::vector<std::string> OwnedArrowArray::ReadUtf8(const char *what) const {
	RequireFormat(schema_, "u", what);
	std::vector<std::string> out;
	if (array_.length == 0) {
		return out;
	}
	// offset int32_t array is the second buffer in the ArrowArray, and the char array is the third buffer.  The first
	// buffer is the validity bitmap, which we don't use since we don't have nulls.
	const int32_t *offsets = static_cast<const int32_t *>(array_.buffers[1]);
	const char *chars = static_cast<const char *>(array_.buffers[2]);
	out.reserve(static_cast<size_t>(array_.length));
	for (int64_t i = 0; i < array_.length; i++) {
		// start chunking u8 at i..i+1
		const int32_t start = offsets[i];
		const int32_t end = offsets[i + 1];
		// An all-empty-string column allocates no data buffer at all, so guard
		// the pointer rather than the length.
		// std::string contructor bc of emplace_back used on std::vector<std::string> out. It takes two arguments, the
		// starting memory position and the lenghth and then takes a bite out of the chars buffer to construct a
		// std::string out of the chunk.
		// empty string for empty chars buffer, otherwise it will segfault when trying to read from a nullptr.  The
		// empty string is a valid string in C++ and is represented by a std::string with length 0 and no data.
		out.emplace_back(chars ? chars + start : "", static_cast<size_t>(end - start));
	}
	return out;
}

std::vector<double> OwnedArrowArray::ReadFloat64(const char *what) const {
	RequireFormat(schema_, "g", what);
	if (array_.length == 0) {
		return {};
	}
	const double *values = static_cast<const double *>(array_.buffers[1]);
	return std::vector<double>(values, values + array_.length);
}

const double *OwnedArrowArray::FixedSizeListFloat64Data(const char *what, int64_t &width) const {
	// "+w:N" -- the width is part of the type, not a struct field.
	const std::string fmt = schema_.format ? schema_.format : "";
	if (fmt.rfind("+w:", 0) != 0) {
		throw InvalidInputException("%s: expected Arrow format '+w:N' (fixed size list), got '%s'", what,
		                            fmt.empty() ? "(null)" : fmt.c_str());
	}
	width = std::strtoll(fmt.c_str() + 3, nullptr, 10);
	if (width <= 0) {
		throw InvalidInputException("%s: fixed size list width must be >= 1, got '%s'", what, fmt.c_str());
	}
	if (array_.n_children != 1 || array_.children == nullptr || array_.children[0] == nullptr) {
		throw InternalException("%s: fixed size list must have exactly one child, got %lld", what,
		                        (int64_t)array_.n_children);
	}
	// The values are the child's, not the parent's: a fixed size list stores no
	// offsets, so the outer array has only a validity buffer and the child holds
	// length * width doubles back to back.
	const ArrowArray &child = *array_.children[0];
	const auto total = static_cast<size_t>(array_.length * width);
	if (static_cast<size_t>(child.length) < total) {
		throw InternalException("%s: fixed size list child holds %lld values, expected %llu", what,
		                        (int64_t)child.length, (uint64_t)total);
	}
	if (total == 0) {
		return nullptr;
	}
	const double *values = static_cast<const double *>(child.buffers[1]);
	// The child may be sliced relative to its parent; honour its offset.
	return values + child.offset;
}

std::vector<double> OwnedArrowArray::ReadFixedSizeListFloat64(const char *what, int64_t &width) const {
	const double *values = FixedSizeListFloat64Data(what, width);
	if (!values) {
		return {};
	}
	return std::vector<double>(values, values + static_cast<size_t>(array_.length * width));
}

void ThrowSc(const char *what, sc_context_t *ctx, sc_status_t status) {
	const char *msg = ctx ? sc_context_last_error(ctx) : nullptr;
	throw InvalidInputException("%s: sc error %d%s%s", what, static_cast<int>(status), msg ? ": " : "", msg ? msg : "");
}

//! One place for every "which model did you mean" failure, so the three cases
//! read the same wherever they are raised.
[[noreturn]] void RejectModelSelection(const char *caller, const std::string &relation, const std::string &name,
                                       int64_t rows) {
	if (!name.empty() && rows == 0) {
		throw InvalidInputException("%s: No model named '%s' in relation '%s'.\n\n"
		                            "Remedy:\n"
		                            "  List what the relation holds, then use one of those names:\n"
		                            "    SELECT name, task, n_trees FROM %s;",
		                            caller, name, relation, relation);
	}
	if (!name.empty()) {
		throw InvalidInputException(
		    "%s: Relation '%s' holds %lld models named '%s'; a name has to identify exactly one.\n\n"
		    "Remedy:\n"
		    "  Inspect the duplicates and drop or rename all but one:\n"
		    "    SELECT name, task, n_trees, random_state FROM %s WHERE name = '%s';",
		    caller, relation, (int64_t)rows, name, relation, name);
	}
	if (rows == 0) {
		throw InvalidInputException(
		    "%s: Model relation '%s' is empty.\n\n"
		    "Remedy:\n"
		    "  Fit a model into it first:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation);
	}
	throw InvalidInputException("%s: Relation '%s' holds %lld models, so this call is ambiguous.\n\n"
	                            "Remedy:\n"
	                            "  Say which one, by the name it was fit with:\n"
	                            "    SELECT name, task, n_trees FROM %s;              -- see what is there\n"
	                            "    ... FROM %s(..., '%s', name := 'my_model');      -- then pick one",
	                            caller, relation, (int64_t)rows, relation, caller, relation);
}

std::string ModelNameFilter(const std::string &name) {
	if (name.empty()) {
		return "";
	}
	return " WHERE name = " + duckdb::KeywordHelper::WriteQuoted(name, '\'');
}

std::string ReadModelTask(duckdb::Connection &conn, const std::string &relation, const std::string &name,
                          const char *caller) {
	const auto q = duckdb::KeywordHelper::WriteOptionallyQuoted(relation);
	auto result = conn.Query("SELECT task FROM " + q + ModelNameFilter(name));
	if (result->HasError()) {
		throw InvalidInputException(
		    "%s: Model relation '%s' has no 'task' column.\n"
		    "  Engine error : %s\n\n"
		    "A model row carries more than its bytes: 'task' says whether it predicts a label or a\n"
		    "number, which is how the prediction column gets its type.\n\n"
		    "Remedy:\n"
		    "  Store the whole row from a fit, not just the blob:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation, result->GetError());
	}
	std::string task;
	int64_t rows = 0;
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			if (rows++ == 0) {
				auto v = chunk->data[0].GetValue(row);
				if (!v.IsNull()) {
					task = v.ToString();
				}
			}
		}
	}
	if (rows != 1) {
		RejectModelSelection(caller, relation, name, rows);
	}
	if (task != "classification" && task != "regression") {
		throw InvalidInputException("%s: model relation '%s' has task '%s'; expected 'classification' or "
		                            "'regression'",
		                            caller, relation, task);
	}
	return task;
}

ScIdTypes ProbeIdTypes(duckdb::Connection &conn, const std::string &relation, const char *caller) {
	const auto q = duckdb::KeywordHelper::WriteOptionallyQuoted(relation);
	// Ids uncast so their own types come through; value cast so a value column
	// that cannot become DOUBLE fails here rather than mid-scan. LIMIT 0 binds
	// the relation without reading it.
	auto probe = conn.Query("SELECT sample_id, feature_id, value::DOUBLE FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException(
		    "%s: Data relation '%s' does not match the required COO triplet schema.\n"
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
		    caller, relation, probe->GetError(), relation);
	}
	ScIdTypes types;
	types.sample_id = probe->types[0];
	types.feature_id = probe->types[1];
	auto require_allowed = [&](const char *column, const duckdb::LogicalType &type, const std::string &view_select) {
		if (!duckdb::IsAllowedIdType(type)) {
			throw InvalidInputException("%s: '%s' field %s has %s, but ids must be %s. Cast it first: CREATE VIEW v "
			                            "AS SELECT %s FROM %s;",
			                            caller, relation, column, type.ToString(), duckdb::AllowedIdTypeList(),
			                            view_select, relation);
		}
	};
	// Integers widen to BIGINT and keep sorting numerically; anything else
	// becomes text.
	const auto cast_to = [](const duckdb::LogicalType &type) {
		return std::string(type.IsIntegral() ? "BIGINT" : "VARCHAR");
	};
	require_allowed("sample_id", types.sample_id,
	                "sample_id::" + cast_to(types.sample_id) + " AS sample_id, feature_id, value");
	require_allowed("feature_id", types.feature_id,
	                "sample_id, feature_id::" + cast_to(types.feature_id) + " AS feature_id, value");
	return types;
}

duckdb::LogicalType ReadModelFeatureIdType(duckdb::Connection &conn, const std::string &relation,
                                           const std::string &name, const char *caller) {
	const auto q = duckdb::KeywordHelper::WriteOptionallyQuoted(relation);
	auto result = conn.Query("SELECT feature_id_type FROM " + q + ModelNameFilter(name));
	if (result->HasError()) {
		throw InvalidInputException(
		    "%s: Could not read 'feature_id_type' from model relation '%s'.\n"
		    "  Engine error : %s\n\n"
		    "A model row records the feature_id type it was fit on, which data must match.\n\n"
		    "Remedy:\n"
		    "  Store the whole row from a fit, refitting if the model predates this column:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation, result->GetError());
	}
	std::string type;
	int64_t rows = 0;
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			if (rows++ == 0) {
				auto v = chunk->data[0].GetValue(row);
				if (!v.IsNull()) {
					type = v.ToString();
				}
			}
		}
	}
	if (rows != 1) {
		RejectModelSelection(caller, relation, name, rows);
	}
	for (const auto &allowed :
	     {duckdb::LogicalType(duckdb::LogicalType::VARCHAR), duckdb::LogicalType(duckdb::LogicalType::BIGINT),
	      duckdb::LogicalType(duckdb::LogicalType::UUID)}) {
		if (type == allowed.ToString()) {
			return allowed;
		}
	}
	throw InvalidInputException("%s: model relation '%s' has feature_id_type '%s'; expected %s", caller, relation, type,
	                            duckdb::AllowedIdTypeList());
}

void RequireModelFeatureIdType(const duckdb::LogicalType &model_type, const duckdb::LogicalType &data_type,
                               const char *caller, const std::string &model_relation, const std::string &model_name,
                               const std::string &data_relation) {
	if (model_type == data_type) {
		return;
	}
	const auto model = model_name.empty() ? "the model in '" + model_relation + "'" : "model '" + model_name + "'";
	throw InvalidInputException("%s: %s was trained on feature_id %s, but '%s' field feature_id has %s. Cast it "
	                            "first: CREATE VIEW v AS SELECT sample_id, feature_id::%s AS feature_id, value FROM "
	                            "%s;",
	                            caller, model, model_type.ToString(), data_relation, data_type.ToString(),
	                            model_type.ToString(), data_relation);
}

void LoadModelFromRelation(duckdb::Connection &conn, const std::string &relation, const std::string &name,
                           const char *caller, sc_context_t *ctx, ScModel &out) {
	const auto q = duckdb::KeywordHelper::WriteOptionallyQuoted(relation);
	auto result = conn.Query("SELECT model_blob FROM " + q + ModelNameFilter(name));
	if (result->HasError()) {
		throw InvalidInputException("%s: model relation '%s' must expose a 'model_blob' column: %s", caller, relation,
		                            result->GetError());
	}

	duckdb::Value blob;
	int64_t rows = 0;
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			if (rows++ == 0) {
				blob = chunk->data[0].GetValue(row);
			}
		}
	}
	// A model table holds exactly one row. More than one is ambiguous -- there
	// is no basis for choosing -- and zero usually means an upstream filter ate
	// it, which is worth saying rather than failing later inside sc.
	if (rows != 1) {
		RejectModelSelection(caller, relation, name, rows);
	}
	if (blob.IsNull()) {
		throw InvalidInputException("%s: model relation '%s' has a NULL model", caller, relation);
	}

	const auto bytes = blob.GetValueUnsafe<duckdb::string_t>();
	if (auto st =
	        sc_model_deserialize(ctx, reinterpret_cast<const uint8_t *>(bytes.GetData()), bytes.GetSize(), &out.ptr);
	    st != SC_OK) {
		ThrowSc(caller, ctx, st);
	}
}

} // namespace miint
