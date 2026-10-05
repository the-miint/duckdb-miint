#pragma once
/*
 * Shared parameter parsing for align_sortmerna / align_sortmerna_rrna.
 *
 * Both functions bind from the same parameter set; only the output schema
 * differs. Keep Bind/InitGlobal/Execute in each table function's .cpp — they
 * are simple enough that pulling them behind a second layer would hurt more
 * than it helps — but factor the parameter parsing here so each function
 * parses identically.
 */

#include "SortMeRNAAligner.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/function/table_function.hpp"
#include <string>
#include <vector>
#include "miint_named_parameter.hpp"

namespace duckdb {

// Parse a LIST(VARCHAR) ref_paths named parameter into std::vector<std::string>.
// `fn_name` is the caller's SQL-facing name, used to prefix error messages.
inline std::vector<std::string> ParseSortMeRNARefPaths(const named_argument_map_t &params, const char *fn_name) {
	auto it = params.find("ref_paths");
	if (it == params.end() || it->second.IsNull()) {
		throw BinderException("%s requires ref_paths (LIST of FASTA paths)", fn_name);
	}
	auto &children = ListValue::GetChildren(it->second);
	if (children.empty()) {
		throw BinderException("%s: ref_paths must be a non-empty list", fn_name);
	}
	std::vector<std::string> result;
	result.reserve(children.size());
	for (auto &child : children) {
		if (child.IsNull()) {
			throw BinderException("%s: ref_paths contains NULL entry", fn_name);
		}
		result.push_back(child.ToString());
	}
	return result;
}

inline void ParseSortMeRNAConfigParams(const named_argument_map_t &params, miint::SortMeRNAConfig &cfg) {
	auto set_i32 = [&](const char *name, int32_t &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<int32_t>();
		}
	};
	auto set_u32 = [&](const char *name, uint32_t &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<uint32_t>();
		}
	};
	auto set_bool = [&](const char *name, bool &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<bool>();
		}
	};
	auto set_double = [&](const char *name, double &dst) {
		auto it = params.find(name);
		if (it != params.end() && !it->second.IsNull()) {
			dst = it->second.GetValue<double>();
		}
	};

	set_i32("num_threads", cfg.num_threads);
	set_i32("match", cfg.match);
	set_i32("mismatch", cfg.mismatch);
	set_i32("gap_open", cfg.gap_open);
	set_i32("gap_ext", cfg.gap_ext);
	set_i32("score_N", cfg.score_N);
	set_double("evalue", cfg.evalue);
	set_u32("seed_win_len", cfg.seed_win_len);
	set_u32("num_alignments", cfg.num_alignments);
	set_bool("best", cfg.best);
	set_bool("paired", cfg.paired);
	set_bool("forward_only", cfg.forward_only);
	set_bool("reverse_only", cfg.reverse_only);
	set_bool("full_search", cfg.full_search);
}

// Register the shared LIST(VARCHAR)/INTEGER/DOUBLE/BOOLEAN named parameters on
// a sortmerna-family TableFunction.
inline void RegisterSortMeRNANamedParameters(TableFunction &tf) {
	AddNamedParameter(tf, "ref_paths", LogicalType::LIST(LogicalType::VARCHAR));
	AddNamedParameter(tf, "num_threads", LogicalType::INTEGER);
	AddNamedParameter(tf, "match", LogicalType::INTEGER);
	AddNamedParameter(tf, "mismatch", LogicalType::INTEGER);
	AddNamedParameter(tf, "gap_open", LogicalType::INTEGER);
	AddNamedParameter(tf, "gap_ext", LogicalType::INTEGER);
	AddNamedParameter(tf, "score_N", LogicalType::INTEGER);
	AddNamedParameter(tf, "evalue", LogicalType::DOUBLE);
	AddNamedParameter(tf, "seed_win_len", LogicalType::UINTEGER);
	AddNamedParameter(tf, "num_alignments", LogicalType::UINTEGER);
	AddNamedParameter(tf, "best", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "paired", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "forward_only", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "reverse_only", LogicalType::BOOLEAN);
	AddNamedParameter(tf, "full_search", LogicalType::BOOLEAN);
}

} // namespace duckdb
