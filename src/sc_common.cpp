#include "sc_common.hpp"
#include <cctype>

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <cstdlib>
#include <cstring>
#include <string>
#include <unordered_map>

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
		                        (long long)array_.n_children);
	}
	// The values are the child's, not the parent's: a fixed size list stores no
	// offsets, so the outer array has only a validity buffer and the child holds
	// length * width doubles back to back.
	const ArrowArray &child = *array_.children[0];
	const auto total = static_cast<size_t>(array_.length * width);
	if (static_cast<size_t>(child.length) < total) {
		throw InternalException("%s: fixed size list child holds %lld values, expected %llu", what,
		                        (long long)child.length, (unsigned long long)total);
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

namespace {

std::string TrimAscii(const std::string &s) {
	const auto first = s.find_first_not_of(" \t\n\r\f\v");
	if (first == std::string::npos) {
		return {};
	}
	return s.substr(first, s.find_last_not_of(" \t\n\r\f\v") - first + 1);
}

std::string LowerAscii(const std::string &s) {
	std::string out = s;
	for (auto &c : out) {
		c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
	}
	return out;
}

//! The digits of an integer id, or empty when `s` is not one.
//!
//! Catches the two ways an integer id is respelled: zero padding ('042') and a
//! trip through DOUBLE ('42.0'). Both compare equal here to '42'.
std::string NumericKey(const std::string &s) {
	std::string body = TrimAscii(s);
	if (body.size() > 2 && body.compare(body.size() - 2, 2, ".0") == 0) {
		body.erase(body.size() - 2);
	}
	if (body.empty() || body.find_first_not_of("0123456789") != std::string::npos) {
		return {};
	}
	const auto first = body.find_first_not_of('0');
	return first == std::string::npos ? "0" : body.substr(first);
}

bool LooksLikeUuid(const std::string &s) {
	if (s.size() != 36) {
		return false;
	}
	for (size_t i = 0; i < s.size(); i++) {
		const bool hyphen = i == 8 || i == 13 || i == 18 || i == 23;
		if (hyphen != (s[i] == '-') || (!hyphen && !std::isxdigit(static_cast<unsigned char>(s[i])))) {
			return false;
		}
	}
	return true;
}

} // namespace

std::string VocabularyMismatchHint(const std::vector<std::string> &dropped, const std::vector<std::string> &vocab,
                                   const std::string &relation) {
	if (dropped.empty() || vocab.empty()) {
		return {};
	}
	// Index the vocabulary once under each normalisation. Only reached when a
	// call is already failing, so an O(vocab) pass costs nothing that matters.
	std::unordered_map<std::string, const std::string *> folded, numeric;
	for (const auto &v : vocab) {
		folded.emplace(LowerAscii(TrimAscii(v)), &v);
		const auto n = NumericKey(v);
		if (!n.empty()) {
			numeric.emplace(n, &v);
		}
	}

	for (const auto &d : dropped) {
		const auto trimmed = TrimAscii(d);
		const char *difference = nullptr;
		const char *fix = nullptr;
		const std::string *match = nullptr;

		if (const auto it = folded.find(LowerAscii(trimmed)); it != folded.end()) {
			match = it->second;
			const bool case_differs = LowerAscii(trimmed) != trimmed || LowerAscii(*match) != *match;
			if (trimmed != d && case_differs) {
				difference = "surrounding whitespace and case";
				fix = "lower(trim(feature_id))";
			} else if (trimmed != d) {
				difference = "surrounding whitespace";
				fix = "trim(feature_id)";
			} else {
				difference = "case";
				fix = LooksLikeUuid(trimmed) && LooksLikeUuid(*match) ? "feature_id::UUID" : "lower(feature_id)";
			}
		} else if (const auto n = NumericKey(d); !n.empty()) {
			if (const auto it2 = numeric.find(n); it2 != numeric.end() && *it2->second != d) {
				match = it2->second;
				// '042' and '42.0' both mean 42; the cast renders it one way.
				difference = "zero padding or a decimal point";
				fix = "feature_id::BIGINT::VARCHAR";
			}
		}
		if (!match) {
			continue;
		}
		return duckdb::StringUtil::Format(
		    "\n\nThe ids match the model's except for %s -- the data has '%s', the model has '%s'. Ids are matched by "
		    "their text, so these are different features.\n"
		    "Remedy:\n"
		    "  Normalise them into a view, then pass that instead:\n"
		    "    CREATE VIEW fixed AS SELECT sample_id, %s AS feature_id, value FROM %s;",
		    difference, d, *match, fix, relation);
	}
	return {};
}

sc_coo_table_t AsScTable(const CooTable &table) {
	const auto &a = table.arrays();
	sc_coo_table_t out {};
	out.rows = a.rows;
	out.rows_schema = a.rows_schema;
	out.cols = a.cols;
	out.cols_schema = a.cols_schema;
	out.vals = a.vals;
	out.vals_schema = a.vals_schema;
	out.sample_ids = a.sample_ids;
	out.sample_ids_schema = a.sample_ids_schema;
	out.feature_ids = a.feature_ids;
	out.feature_ids_schema = a.feature_ids_schema;
	out.n_samples = a.n_samples;
	out.n_features = a.n_features;
	return out;
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
		    caller, relation, (long long)rows, name, relation, name);
	}
	if (rows == 0) {
		throw InvalidInputException(
		    "%s: Model relation '%s' is empty.\n\n"
		    "Remedy:\n"
		    "  Fit a model into it first:\n"
		    "    CREATE TABLE models AS SELECT * FROM sc_fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation);
	}
	throw InvalidInputException("%s: Relation '%s' holds %lld models, so this call is ambiguous.\n\n"
	                            "Remedy:\n"
	                            "  Say which one, by the name it was fit with:\n"
	                            "    SELECT name, task, n_trees FROM %s;              -- see what is there\n"
	                            "    ... FROM %s(..., '%s', name := 'my_model');      -- then pick one",
	                            caller, relation, (long long)rows, relation, caller, relation);
}

std::string ModelNameFilter(const std::string &name) {
	if (name.empty()) {
		return "";
	}
	return " WHERE name = " + duckdb::KeywordHelper::WriteQuoted(name, '\'');
}

ScModelTypes ReadModelTypes(duckdb::Connection &conn, duckdb::ClientContext &context, const std::string &relation,
                            const std::string &name, const char *caller) {
	ScModelTypes out;
	const auto q = duckdb::KeywordHelper::WriteOptionallyQuoted(relation);
	auto result = conn.Query("SELECT feature_id_type, target_type FROM " + q + ModelNameFilter(name));
	// Absent columns are not an error: a model table from before these existed
	// still predicts, it just cannot say what its ids used to be.
	if (result->HasError()) {
		return out;
	}
	auto parse = [&](const duckdb::Value &v, duckdb::LogicalType &into) {
		if (v.IsNull()) {
			return;
		}
		try {
			into = duckdb::TransformStringToLogicalType(v.ToString(), context);
		} catch (const std::exception &e) {
			// A stored type nobody can parse is a corrupt row, not a reason to
			// refuse the prediction -- fall back to text and say so.
			throw InvalidInputException("%s: model relation '%s' stores an unreadable id type '%s': %s", caller,
			                            relation, v.ToString(), e.what());
		}
	};
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			parse(chunk->data[0].GetValue(row), out.feature_id_type);
			parse(chunk->data[1].GetValue(row), out.target_type);
			return out; // one row; RejectModelSelection already guards the rest
		}
	}
	return out;
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
		    "    CREATE TABLE models AS SELECT * FROM sc_fit_classifier('counts', 'meta', name := 'm1');",
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
