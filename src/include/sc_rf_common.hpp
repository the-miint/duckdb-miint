#pragma once

#include "coo_builder.hpp"

#include "duckdb/common/string_util.hpp"
#include "id_column_utils.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/parser/keyword_helper.hpp"

#include <string>
#include <unordered_map>
#include <vector>

#include "sc.h"

namespace duckdb {
namespace sc_rf {

//! The metadata column holding the sample id. Everything else is a candidate
//! target.
constexpr const char *kSampleIdColumn = "sample_id";

//! What fitting and cross-validating both need to know about their inputs.
//!
//! Both read the same two relations -- COO triplets, and (sample_id, target) --
//! validate them the same way and fit the same forests, so the helpers below
//! take this instead of either function's own bind data. A table function's bind
//! data inherits it, which is why the moved code still reads `bind.data_relation`.
//!
//! `caller` is the SQL function's name, so an error names the function the user
//! actually called rather than whichever one the code was first written for.
struct ScTrainingInput {
	string data_relation;
	string metadata_relation;
	//! Resolved at bind time: the caller's `target_column :=`, or the metadata
	//! relation's single non-sample_id column.
	string target_column;
	const char *caller = "sc_fit";
	//! The data relation's id types, captured at bind so every id this call
	//! returns goes back out as the type it came in as. Ids travel through sc as
	//! text -- that is how a BIGINT 42 and a VARCHAR '42' are one feature -- but
	//! handing a VARCHAR back to a caller whose table is BIGINT changes what
	//! ORDER BY means and spreads the string type into everything downstream.
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
	//! The metadata target's type, stored on the model so a classifier's labels
	//! come back as they went in. A regressor predicts a continuous value, so its
	//! prediction is DOUBLE whatever the column was.
	LogicalType target_type = LogicalType::VARCHAR;
};

//! The data relation's `(sample_id, feature_id)` types, validating the triplet
//! schema in the same pass.
//!
//! One bind-time probe answers both questions, and a relation that is not a COO
//! triplet fails here -- with the same message the scan would have given -- so a
//! wrong relation name costs nothing rather than a full scan.
struct ScCooIdTypes {
	LogicalType sample_id_type = LogicalType::VARCHAR;
	LogicalType feature_id_type = LogicalType::VARCHAR;
};
ScCooIdTypes DetectCooIdTypes(Connection &conn, const string &relation, const char *caller);

//! The type of `column` in `relation`, whatever it is.
//!
//! Targets carry no such restriction: any type that renders to text can label a
//! sample, and classification hands the label back by casting the text home.
LogicalType DetectColumnType(Connection &conn, const string &relation, const string &column, const char *caller);

//! Export one Arrow array of targets. sc takes Utf8 labels for a classifier and
//! Float64 for a regressor; the two are otherwise identical at this boundary.
struct TargetArray {
	ArrowArray array {};
	ArrowSchema schema {};
	std::vector<std::string> labels; // classification
	std::vector<double> numbers;     // regression
	std::vector<int32_t> offsets;
	std::vector<char> chars;
	const void *buffers[3] = {nullptr, nullptr, nullptr};
};

//! Fill `t` as an Arrow array sc can borrow. `t` owns the buffers and must
//! outlive the sc call -- sc's import borrows in place and never releases.
void BuildTargets(TargetArray &t, bool classification);

//! Reject a relation that is not `(sample_id, feature_id, value)`.
//!
//! Raised from the bind probe and from every scan, so it lives here rather than
//! in five near-identical copies -- the remedy it prints is the whole point of
//! the message, and a copy that drifts is worse than no message.
[[noreturn]] void ThrowNotCooTriplet(const string &relation, const string &engine_error, const char *caller);

//! Scan the data relation into `builder`, rejecting NULLs.
void ScanCounts(Connection &conn, const ScTrainingInput &bind, miint::CooBuilder &builder);

//! Reject duplicate cells, showing the values so the caller can tell a join
//! fanout (identical values, deduplicate) from repeat measurements (differing
//! values, maybe sum).
void RequireNoDuplicateCells(const miint::CooBuilder &builder, const ScTrainingInput &bind);

//! Fill `params` with sklearn's defaults for the task.
void ApplyDefaults(sc_rf_params_t &params, bool classification);

// Parse one named parameter's SQL value into its sc tagged union. `caller` names
// the SQL function in any error.
void ParseMaxFeatures(const Value &v, const char *caller, sc_max_features_t &out);
void ParseMinSamples(const Value &v, const char *name, uint64_t min_count, const char *caller, sc_min_samples_t &out);
void ParseMaxSamples(const Value &v, const char *caller, sc_max_samples_t &out);
void ParseCriterion(const Value &v, bool classification, const char *caller, sc_criterion_t &out);

//! Work out which metadata column holds the target.
std::string ResolveTargetColumn(Connection &conn, const std::string &relation, const char *caller);

// Templates: defined here because the target type varies by task -- Utf8 labels
// for a classifier, Float64 values for a regressor.

//! Scan the metadata relation into a sample_id -> target map, rejecting NULLs
//! and duplicate labels.
template <class T>
std::unordered_map<std::string, T> ScanTargets(Connection &conn, const ScTrainingInput &bind, const char *cast) {
	const auto q = KeywordHelper::WriteOptionallyQuoted(bind.metadata_relation);
	const auto col = KeywordHelper::WriteOptionallyQuoted(bind.target_column);
	// sample_id is cast exactly as the data relation's is (ScanCounts), so both
	// sides of the join render an id the same way. Without it the key would come
	// from Value::ToString(), which agrees today but is a separate code path.
	auto result = conn.Query("SELECT sample_id::VARCHAR, " + col + "::" + cast + " FROM " + q);
	if (result->HasError()) {
		throw InvalidInputException("%s: metadata relation '%s' must expose (sample_id, %s) castable to %s: %s",
		                            bind.caller, bind.metadata_relation, bind.target_column, cast, result->GetError());
	}
	std::unordered_map<std::string, T> targets;
	while (auto chunk = result->Fetch()) {
		for (idx_t row = 0; row < chunk->size(); row++) {
			auto s = chunk->data[0].GetValue(row);
			auto t = chunk->data[1].GetValue(row);
			if (s.IsNull() || t.IsNull()) {
				throw InvalidInputException("%s: NULL in metadata relation '%s' (sample_id and %s must be "
				                            "non-NULL)",
				                            bind.caller, bind.metadata_relation, bind.target_column);
			}
			const auto sample = s.ToString();
			// Two labels for one sample cannot be reconciled -- unlike duplicate
			// counts, there is no defensible way to combine them.
			if (!targets.emplace(sample, t.GetValue<T>()).second) {
				throw InvalidInputException(
				    "%s: sample '%s' has more than one %s in metadata relation '%s'; deduplicate it first", bind.caller,
				    sample, bind.target_column, bind.metadata_relation);
			}
		}
	}
	return targets;
}

//! Both relations must describe exactly the same samples. Report both
//! directions at once so one round trip fixes the whole mismatch.
template <class T>
void RequireSameSamples(const std::vector<std::string> &data_samples, const std::unordered_map<std::string, T> &targets,
                        const ScTrainingInput &bind) {
	// unlabelled: Sample is in Data, but missing from Metadata
	std::vector<std::string> unlabelled;
	for (const auto &s : data_samples) {
		if (targets.find(s) == targets.end()) {
			unlabelled.push_back(s);
		}
	}
	// undated = un-data'd: Sample is in Metadata, but missing from Data
	// in other words, which targets were missing from the target mapping above
	std::vector<std::string> undated;
	if (targets.size() + unlabelled.size() != data_samples.size()) {
		std::unordered_map<std::string, bool> present;
		for (const auto &s : data_samples) {
			present.emplace(s, true);
		}
		for (const auto &kv : targets) {
			if (present.find(kv.first) == present.end()) {
				undated.push_back(kv.first);
			}
		}
	}
	if (unlabelled.empty() && undated.empty()) {
		return;
	}

	auto sample_list = [](const std::vector<std::string> &v) {
		string out;
		const size_t shown = v.size() < 5 ? v.size() : 5;
		for (size_t i = 0; i < shown; i++) {
			out += (i ? ", " : "") + v[i];
		}
		if (v.size() > shown) {
			out += StringUtil::Format(", ... (+%llu more)", (unsigned long long)(v.size() - shown));
		}
		return out;
	};
	string msg = StringUtil::Format("%s: The data and metadata relations describe different sample sets.\n"
	                                "  Every sample with counts must have exactly one label, and vice versa.",
	                                bind.caller);
	if (!unlabelled.empty()) {
		msg +=
		    StringUtil::Format("\n  %llu sample(s) in '%s' with no %s: %s", (unsigned long long)unlabelled.size(),
		                       bind.data_relation.c_str(), bind.target_column.c_str(), sample_list(unlabelled).c_str());
	}
	if (!undated.empty()) {
		msg += StringUtil::Format("\n  %llu sample(s) in '%s' with no data: %s", (unsigned long long)undated.size(),
		                          bind.metadata_relation.c_str(), sample_list(undated).c_str());
	}
	msg += StringUtil::Format(
	    "\n\nRemedy:\n"
	    "  Restrict both sides to the samples they share, then refit:\n"
	    "    CREATE VIEW shared AS SELECT DISTINCT sample_id FROM %s SEMI JOIN %s USING (sample_id);\n"
	    "    CREATE VIEW counts AS SELECT * FROM %s SEMI JOIN shared USING (sample_id);\n"
	    "    CREATE VIEW labels AS SELECT * FROM %s SEMI JOIN shared USING (sample_id);\n"
	    "  Or fix the upstream join if the mismatch is unexpected -- a truncated metadata\n"
	    "  export trains a model on fewer samples than you think.",
	    bind.data_relation.c_str(), bind.metadata_relation.c_str(), bind.data_relation.c_str(),
	    bind.metadata_relation.c_str());
	throw InvalidInputException(msg);
}

} // namespace sc_rf
} // namespace duckdb
