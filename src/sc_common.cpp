#include "sc_common.hpp"

#include "id_column_utils.hpp"
#include "miint_log.hpp"
#include "miint_streaming_query.hpp"

#include "duckdb/common/error_data.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <cmath>
#include <cstdlib>
#include <cstring>
#include <string>
#include <string_view>

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

//! The readers below take buffers at face value: no slice offset, no nulls.
//! sc's exports never use either; anything else is refused, not misread.
void RequirePlain(const ArrowArray &array, const char *what, bool allow_offset = false) {
	if ((!allow_offset && array.offset != 0) || array.null_count != 0) {
		throw InternalException("%s: Arrow array has offset %lld and %lld nulls; expected neither", what,
		                        (int64_t)array.offset, (int64_t)array.null_count);
	}
}

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
	RequirePlain(array_, what);
	std::vector<std::string> out;
	if (array_.length == 0) {
		return out;
	}
	// buffers: [validity (unused, no nulls), int32 offsets, chars].
	const int32_t *offsets = static_cast<const int32_t *>(array_.buffers[1]);
	const char *chars = static_cast<const char *>(array_.buffers[2]);
	out.reserve(static_cast<size_t>(array_.length));
	for (int64_t i = 0; i < array_.length; i++) {
		const int32_t start = offsets[i];
		const int32_t end = offsets[i + 1];
		// An all-empty-string column allocates no data buffer at all, so guard
		// the pointer rather than the length.
		out.emplace_back(chars ? chars + start : "", static_cast<size_t>(end - start));
	}
	return out;
}

std::vector<double> OwnedArrowArray::ReadFloat64(const char *what) const {
	RequireFormat(schema_, "g", what);
	RequirePlain(array_, what);
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
	RequirePlain(array_, what);
	if (schema_.n_children != 1 || schema_.children == nullptr || schema_.children[0] == nullptr) {
		throw InternalException("%s: fixed size list schema must have exactly one child, got %lld", what,
		                        (int64_t)schema_.n_children);
	}
	RequireFormat(*schema_.children[0], "g", what);
	const ArrowArray &child = *array_.children[0];
	if (child.n_buffers != 2) {
		throw InternalException("%s: Float64 child has %lld buffers, expected 2", what, (int64_t)child.n_buffers);
	}
	RequirePlain(child, what, true);
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

void TakeArray(ArrowArray &src, ArrowSchema &src_schema, OwnedArrowArray &dst) {
	*dst.array() = src;
	*dst.schema() = src_schema;
	src.release = nullptr;
	src_schema.release = nullptr;
}

void ThrowSc(const char *what, sc_context_t *ctx, sc_status_t status) {
	// A call with no context records its error on the calling thread instead.
	const char *msg = ctx ? sc_context_last_error(ctx) : sc_thread_last_error();
	throw InvalidInputException("%s: sc error %d%s%s", what, static_cast<int>(status), msg ? ": " : "", msg ? msg : "");
}

void RequireNotNull(const std::string &parameter, const duckdb::Value &v, const char *caller) {
	if (v.IsNull()) {
		throw InvalidInputException("%s: named parameter '%s' must not be NULL", caller, parameter);
	}
}

void OpenScContext(int32_t n_threads, const char *caller, ScContext &out) {
	sc_config_t config {};
#ifndef __EMSCRIPTEN__
	config.n_threads = n_threads;
#else
	(void)n_threads;
#endif
	if (auto st = sc_context_new(&config, &out.ptr); st != SC_OK) {
		// No context exists yet, so sc's message is on this thread.
		ThrowSc(caller, nullptr, st);
	}
}

//! Every "which model did you mean" failure: no such name, a duplicated name, an
//! empty relation, or several models and no name. `takes_data` says whether the
//! caller's first argument is a data relation, which shapes the call the remedy
//! shows.
[[noreturn]] void RejectModelSelection(const char *caller, const std::string &relation, const std::string &name,
                                       int64_t rows, bool takes_data) {
	// The remedies are SQL to run, so every value they splice is quoted: the
	// relation as an identifier after FROM, as a string literal in a call.
	const auto from = duckdb::SQLIdentifier::ToString(relation);
	const auto arg = duckdb::SQLString::ToString(relation);
	if (!name.empty() && rows == 0) {
		throw InvalidInputException("%s: No model named '%s' in relation '%s'.\n\n"
		                            "Remedy:\n"
		                            "  List what the relation holds, then use one of those names:\n"
		                            "    SELECT name, task, n_estimators FROM %s;",
		                            caller, name, relation, from);
	}
	if (!name.empty()) {
		throw InvalidInputException(
		    "%s: Relation '%s' holds %lld models named '%s'; a name has to identify exactly one.\n\n"
		    "Remedy:\n"
		    "  Inspect the duplicates and drop or rename all but one:\n"
		    "    SELECT name, task, n_estimators, random_state FROM %s WHERE name = %s;",
		    caller, relation, (int64_t)rows, name, from, duckdb::SQLString::ToString(name));
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
	                            "    SELECT name, task, n_estimators FROM %s;  -- see what is there\n"
	                            "    SELECT * FROM %s(%s%s, name := 'my_model');  -- then pick one",
	                            caller, relation, (int64_t)rows, from, caller, takes_data ? "'counts', " : "", arg);
}

std::string ModelNameFilter(const std::string &name) {
	if (name.empty()) {
		return "";
	}
	return " WHERE name = " + duckdb::SQLString::ToString(name);
}

namespace {

//! Name the column a failed model-relation read lacked. The engine reports the
//! first column it could not bind, which is not always the one to fix: without
//! a `name` column the WHERE clause fails before `task` is looked at.
[[noreturn]] void ThrowModelReadError(duckdb::Connection &conn, const std::string &relation, const std::string &name,
                                      const char *caller, bool with_task, const std::string &error) {
	const auto q = duckdb::SQLIdentifier::ToString(relation);
	auto probe = conn.Query("SELECT * FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: Could not read model relation '%s'.\n  Engine error : %s", caller, relation,
		                            probe->GetError());
	}
	const auto has = [&](const char *column) {
		for (const auto &n : probe->GetNames()) {
			if (n == column) {
				return true;
			}
		}
		return false;
	};
	if (!name.empty() && !has("name")) {
		throw InvalidInputException(
		    "%s: Model relation '%s' has no 'name' column, so name := '%s' selects nothing.\n\n"
		    "Remedy:\n"
		    "  Leave out name := for a relation holding one model, or store the whole row from a fit:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation, name);
	}
	if (with_task && !has("task")) {
		throw InvalidInputException(
		    "%s: Model relation '%s' has no 'task' column.\n\n"
		    "A model row carries more than its bytes: 'task' says whether it predicts a label or a\n"
		    "number, which is how the prediction column gets its type.\n\n"
		    "Remedy:\n"
		    "  Store the whole row from a fit, not just the blob:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation);
	}
	if (!has("feature_id_type")) {
		throw InvalidInputException(
		    "%s: Could not read 'feature_id_type' from model relation '%s': it has no such column.\n\n"
		    "A model row records the feature_id type it was fit on, which data must match.\n\n"
		    "Remedy:\n"
		    "  Store the whole row from a fit:\n"
		    "    CREATE TABLE models AS SELECT * FROM fit_classifier('counts', 'meta', name := 'm1');",
		    caller, relation);
	}
	throw InvalidInputException("%s: Could not read model relation '%s'.\n  Engine error : %s", caller, relation,
	                            error);
}

} // namespace

ScModelHeader ReadModelHeader(duckdb::Connection &conn, const std::string &relation, const std::string &name,
                              const char *caller, bool with_task) {
	const auto q = duckdb::SQLIdentifier::ToString(relation);
	auto result = conn.Query(std::string("SELECT feature_id_type::VARCHAR") + (with_task ? ", task::VARCHAR" : "") +
	                         " FROM " + q + ModelNameFilter(name));
	if (result->HasError()) {
		ThrowModelReadError(conn, relation, name, caller, with_task, result->GetError());
	}
	duckdb::Value type, task;
	int64_t rows = 0;
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			if (rows++ == 0) {
				type = chunk->data[0].GetValue(row);
				if (with_task) {
					task = chunk->data[1].GetValue(row);
				}
			}
		}
	}
	// The callers that read a task are exactly the ones whose first argument is
	// a data relation, which shapes the remedy's example call.
	if (rows != 1) {
		RejectModelSelection(caller, relation, name, rows, with_task);
	}

	ScModelHeader header;
	if (with_task) {
		header.task = task.IsNull() ? "" : task.ToString();
		if (header.task != kClassification && header.task != kRegression) {
			throw InvalidInputException("%s: model relation '%s' has task '%s'; expected '%s' or '%s'", caller,
			                            relation, header.task, kClassification, kRegression);
		}
	}
	const auto type_name = type.IsNull() ? std::string() : type.ToString();
	for (const auto &allowed :
	     {duckdb::LogicalType(duckdb::LogicalType::VARCHAR), duckdb::LogicalType(duckdb::LogicalType::BIGINT),
	      duckdb::LogicalType(duckdb::LogicalType::UUID)}) {
		if (type_name == allowed.ToString()) {
			header.feature_id_type = allowed;
			return header;
		}
	}
	throw InvalidInputException("%s: model relation '%s' has feature_id_type '%s'; expected %s", caller, relation,
	                            type_name, duckdb::AllowedIdTypeList());
}

namespace {

[[noreturn]] void ThrowCooSchemaError(const char *caller, const std::string &relation, const std::string &error) {
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
	                            caller, relation, error, duckdb::SQLIdentifier::ToString(relation));
}

} // namespace

ScIdTypes ProbeIdTypes(duckdb::Connection &conn, const std::string &relation, const char *caller) {
	const auto q = duckdb::SQLIdentifier::ToString(relation);
	// Ids uncast so their own types come through; value cast so a value column
	// that cannot become DOUBLE fails here rather than mid-scan. LIMIT 0 binds
	// the relation without reading it.
	auto probe = conn.Query("SELECT sample_id, feature_id, value::DOUBLE FROM " + q + " LIMIT 0");
	if (probe->HasError()) {
		ThrowCooSchemaError(caller, relation, probe->GetError());
	}
	ScIdTypes types;
	types.sample_id = probe->GetTypes()[0];
	types.feature_id = probe->GetTypes()[1];
	auto require_allowed = [&](const char *column, const duckdb::LogicalType &type, const std::string &view_select) {
		if (!duckdb::IsAllowedIdType(type)) {
			throw InvalidInputException("%s: '%s' field %s has %s, but ids must be %s. Cast it first: CREATE VIEW v "
			                            "AS SELECT %s FROM %s;",
			                            caller, relation, column, type.ToString(), duckdb::AllowedIdTypeList(),
			                            view_select, q);
		}
	};
	// Integers widen to BIGINT, which results return as BIGINT; anything else
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

void ScanCooRelation(duckdb::ClientContext &context, duckdb::Connection &conn, const std::string &relation,
                     const char *caller, ScCooBuilder &builder) {
	using duckdb::idx_t;
	using duckdb::string_t;
	using duckdb::UnifiedVectorFormat;
	const auto q = duckdb::SQLIdentifier::ToString(relation);
	// The casts are required, not cosmetic. Reading a vector's buffer directly
	// assumes its physical type, and DuckDB infers `42.0` as DECIMAL(3,1)
	// (physical INT16), an INTEGER count column as INT32, and woltka's feature
	// ids as BIGINT or UUID. Casting in SQL makes DuckDB do the conversion and
	// guarantees the layout this loop reads.
	//
	// Streamed, so one chunk is held at a time rather than the whole relation
	// next to the builder's own copy of it.
	auto result = duckdb::SubmitStream(conn, "SELECT sample_id::VARCHAR, feature_id::VARCHAR, value::DOUBLE FROM " + q);
	if (result->HasError()) {
		ThrowCooSchemaError(caller, relation, result->GetError());
	}
	// A stream surfaces an error raised in a later chunk -- a value that fails the
	// DOUBLE cast -- in one of two ways: Fetch returns no chunk with the error
	// recorded, or, when chunks were already buffered, the next Fetch throws.
	// Both must fail the scan. The first, unchecked, looks exactly like the end
	// of the relation and would fit a model on whatever preceded the bad row.
	const auto next_chunk = [&]() {
		duckdb::unique_ptr<duckdb::DataChunk> chunk;
		try {
			chunk = result->Fetch();
		} catch (std::exception &ex) {
			throw InvalidInputException("%s: reading data relation '%s' failed: %s", caller, relation,
			                            result->HasError() ? result->GetError() : duckdb::ErrorData(ex).Message());
		}
		if (!chunk && result->HasError()) {
			throw InvalidInputException("%s: reading data relation '%s' failed: %s", caller, relation,
			                            result->GetError());
		}
		return chunk;
	};
	// UnifiedVectorFormat rather than Vector::GetValue: GetValue builds a
	// heap-allocating Value per cell, while
	// these string_t views are read in place and ScCooBuilder::Append copies an
	// id only the first time it sees it. Unified rather than flat because a
	// column may arrive constant- or dictionary-encoded.
	while (auto chunk = next_chunk()) {
		context.BlockingInterruptCheck();
		const idx_t n = chunk->size();
		UnifiedVectorFormat sf, ff, vf;
		chunk->data[0].ToUnifiedFormat(sf);
		chunk->data[1].ToUnifiedFormat(ff);
		chunk->data[2].ToUnifiedFormat(vf);
		const auto *samples = UnifiedVectorFormat::GetData<string_t>(sf);
		const auto *features = UnifiedVectorFormat::GetData<string_t>(ff);
		const auto *values = UnifiedVectorFormat::GetData<double>(vf);
		for (idx_t row = 0; row < n; row++) {
			const auto si = sf.sel->get_index(row);
			const auto fi = ff.sel->get_index(row);
			const auto vi = vf.sel->get_index(row);
			// A NULL here means a broken join upstream, not a zero count. An
			// absent cell is already zero in a sparse matrix, so guessing at a
			// NULL would hide the mistake.
			if (!sf.validity.RowIsValid(si) || !ff.validity.RowIsValid(fi) || !vf.validity.RowIsValid(vi)) {
				throw InvalidInputException(
				    "%s: NULL in data relation '%s' (sample_id/feature_id/value must all be non-NULL)", caller,
				    relation);
			}
			const string_t &s = samples[si];
			const string_t &f = features[fi];
			// NaN and infinity fail every `x <= threshold`, so a tree routes them
			// right and predicts confidently; sc's fit panics on them instead.
			if (!std::isfinite(values[vi])) {
				throw InvalidInputException("%s: non-finite value %s in data relation '%s' (sample '%s', feature "
				                            "'%s'); values must be finite",
				                            caller, std::to_string(values[vi]), relation, s.GetString(), f.GetString());
			}
			builder.Append(std::string_view(s.GetData(), s.GetSize()), std::string_view(f.GetData(), f.GetSize()),
			               values[vi]);
		}
	}
}

void RequireNoDuplicateCells(const ScCooBuilder &builder, const char *caller, const std::string &relation) {
	const auto report = builder.FindDuplicateCells();
	if (report.Empty()) {
		return;
	}
	auto msg = duckdb::StringUtil::Format("%s: %llu duplicate (sample_id, feature_id) cell(s) in '%s'.", caller,
	                                      (uint64_t)report.duplicate_cells, relation);
	for (const auto &c : report.examples) {
		std::string values;
		for (size_t i = 0; i < c.values.size(); i++) {
			values += (i ? ", " : "") + duckdb::StringUtil::Format("%g", c.values[i]);
		}
		msg += duckdb::StringUtil::Format("\n  %s / %s  x%llu  values: %s", c.sample_id, c.feature_id,
		                                  (uint64_t)c.count, values);
	}
	msg += "\nIdentical values suggest a join that duplicated rows (deduplicate); differing values suggest repeat "
	       "measurements (aggregate deliberately).";
	throw InvalidInputException(msg);
}

std::unique_ptr<ScCooTable> ScanForModel(duckdb::ClientContext &context, duckdb::Connection &conn, const ScModel &model,
                                         const std::string &relation, const char *caller, const char *empty_effect) {
	// The vocabulary comes out of the model, in the model's column order.
	OwnedArrowArray vocab;
	if (auto st = sc_model_feature_ids(model.ptr, vocab.array(), vocab.schema()); st != SC_OK) {
		ThrowSc("sc_model_feature_ids", nullptr, st);
	}
	ScCooBuilder builder;
	builder.SetFeatureVocabulary(vocab.ReadUtf8("sc_model_feature_ids"));
	ScanCooRelation(context, conn, relation, caller, builder);
	RequireNoDuplicateCells(builder, caller, relation);

	const auto dropped = builder.DroppedCells();
	auto table = builder.Finalize();
	if (!table) {
		throw InvalidInputException("%s: data relation '%s' produced no samples", caller, relation);
	}
	// Every cell dropped means the model shares no vocabulary at all with this
	// data -- a different reference database or pipeline, not a borderline case.
	// Applying the model to an all-zero matrix would return confident nonsense.
	if (dropped > 0 && table->NumNonZeros() == 0) {
		throw InvalidInputException("%s: none of the %llu cells in '%s' use a feature this model was trained on; "
		                            "the data and the model do not share a feature vocabulary",
		                            caller, (uint64_t)dropped, relation);
	}
	if (dropped > 0) {
		size_t empty_samples = 0;
		for (auto f : table->KnownFeatureFraction()) {
			if (f == 0.0) {
				empty_samples++;
			}
		}
		const auto empty_note = empty_samples > 0 ? " -- " + std::to_string(empty_samples) +
		                                                " sample(s) retained no features at all" + empty_effect
		                                          : std::string();
		EmitWarning(context,
		            "%s: dropped %llu cell(s) from '%s' whose feature the model was not trained on%s. See the "
		            "known_feature_fraction column.",
		            caller, (uint64_t)dropped, relation.c_str(), empty_note.c_str());
	}
	return table;
}

std::vector<double> ReadProbaMatrix(const OwnedArrowArray &proba, size_t n_classes, size_t n_samples,
                                    const char *what) {
	int64_t width = 0;
	auto p = proba.ReadFixedSizeListFloat64(what, width);
	if (static_cast<size_t>(width) != n_classes || p.size() != n_samples * n_classes) {
		throw InternalException("%s: %llu probabilities of width %lld for %llu samples x %llu classes", what,
		                        (uint64_t)p.size(), (int64_t)width, (uint64_t)n_samples, (uint64_t)n_classes);
	}
	return p;
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
	                            model_type.ToString(), duckdb::SQLIdentifier::ToString(data_relation));
}

void LoadModelFromRelation(duckdb::Connection &conn, const std::string &relation, const std::string &name,
                           const char *caller, const std::string &task, const duckdb::LogicalType &feature_id_type,
                           sc_context_t *ctx, ScModel &out) {
	const auto q = duckdb::SQLIdentifier::ToString(relation);
	auto result = conn.Query("SELECT model_blob, feature_id_type::VARCHAR" +
	                         std::string(task.empty() ? "" : ", task::VARCHAR") + " FROM " + q + ModelNameFilter(name));
	if (result->HasError()) {
		throw InvalidInputException("%s: model relation '%s' must expose a 'model_blob' column: %s", caller, relation,
		                            result->GetError());
	}
	// Read below as raw bytes, which only a BLOB holds; any other type would be
	// an INTERNAL error, and that invalidates the whole database.
	if (result->GetTypes()[0].id() != duckdb::LogicalTypeId::BLOB) {
		throw InvalidInputException("%s: model_blob in '%s' is %s; it must be the BLOB a fit produced", caller,
		                            relation, result->GetTypes()[0].ToString());
	}

	duckdb::Value blob, read_type, read_task;
	int64_t rows = 0;
	while (auto chunk = result->Fetch()) {
		for (duckdb::idx_t row = 0; row < chunk->size(); row++) {
			if (rows++ == 0) {
				blob = chunk->data[0].GetValue(row);
				read_type = chunk->data[1].GetValue(row);
				if (!task.empty()) {
					read_task = chunk->data[2].GetValue(row);
				}
			}
		}
	}
	// Bind found exactly one row, so another count means the relation changed.
	if (rows != 1) {
		throw InvalidInputException("%s: model relation '%s' changed between planning and execution: it selected "
		                            "one model, now %lld",
		                            caller, relation, (int64_t)rows);
	}
	if (blob.IsNull()) {
		throw InvalidInputException("%s: model relation '%s' has a NULL model", caller, relation);
	}
	const auto changed = [&](const char *column, const std::string &was, const duckdb::Value &now) {
		throw InvalidInputException("%s: model relation '%s' changed between planning and execution: its %s was %s, "
		                            "now %s",
		                            caller, relation, column, was, now.IsNull() ? "NULL" : now.ToString());
	};
	if (!task.empty() && (read_task.IsNull() || read_task.ToString() != task)) {
		changed("task", task, read_task);
	}
	if (read_type.IsNull() || read_type.ToString() != feature_id_type.ToString()) {
		changed("feature_id_type", feature_id_type.ToString(), read_type);
	}

	const auto bytes = blob.GetValueUnsafe<duckdb::string_t>();
	if (auto st =
	        sc_model_deserialize(ctx, reinterpret_cast<const uint8_t *>(bytes.GetData()), bytes.GetSize(), &out.ptr);
	    st != SC_OK) {
		ThrowSc(caller, ctx, st);
	}
}

} // namespace miint
