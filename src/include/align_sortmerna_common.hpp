#pragma once
/*
 * What align_sortmerna and align_sortmerna_rrna share: their named
 * parameters, the Bind-time check of where the references come from, and the
 * InitGlobal-time construction of the aligner from them. Bind/InitGlobal/
 * Execute themselves stay in each table function's .cpp and call these.
 */

#include "SortMeRNAAligner.hpp"
#include "catalog_utils.hpp"
#include "sequence_table_reader.hpp"
#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/value.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include <memory>
#include <string>
#include <vector>
#include "miint_named_parameter.hpp"

namespace duckdb {

// Parse the LIST(VARCHAR) ref_paths named parameter into std::vector<std::string>.
// `fn_name` is the caller's SQL-facing name, used to prefix error messages.
// SortMeRNA builds one index from one file, so the list must hold exactly one
// path; the parameter stays a list so existing calls keep working.
inline std::vector<std::string> ParseSortMeRNARefPaths(const Value &ref_paths, const char *fn_name) {
	auto &children = ListValue::GetChildren(ref_paths);
	if (children.empty()) {
		throw BinderException("%s: ref_paths must be a non-empty list", fn_name);
	}
	if (children.size() > 1) {
		throw BinderException("%s: ref_paths accepts one FASTA path (got %llu); concatenate the files or pass "
		                      "subject_table",
		                      fn_name, static_cast<idx_t>(children.size()));
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

// Where the references come from: one FASTA file (ref_paths) or a relation
// with read_id and sequence1 columns (subject_table). Bind fills exactly one.
struct SortMeRNAReferences {
	std::vector<std::string> ref_paths;
	std::string subject_table;
	SequenceTableSchema subject_schema;

	// Type of the output reference column(s). FASTA header ids are VARCHAR; a
	// subject table's ids keep their read_id type so results join back to it.
	LogicalType IdType() const {
		return subject_table.empty() ? LogicalType::VARCHAR : subject_schema.id_type;
	}
};

// Bind-time half: name checks and schema validation only. The subject table is
// read in InitGlobal (MakeSortMeRNAAligner).
inline SortMeRNAReferences BindSortMeRNAReferences(ClientContext &context, TableFunctionBindInput &input,
                                                   const char *fn_name) {
	auto &params = input.named_parameters;
	auto paths_it = params.find("ref_paths");
	auto subject_it = params.find("subject_table");
	const bool has_paths = paths_it != params.end() && !paths_it->second.IsNull();
	const bool has_subject = subject_it != params.end() && !subject_it->second.IsNull();

	SortMeRNAReferences refs;
	if (has_subject) {
		refs.subject_table = subject_it->second.ToString();
		RejectCTERelationName(input, refs.subject_table);
	}
	if (has_paths == has_subject) {
		throw BinderException("%s requires exactly one of ref_paths or subject_table", fn_name);
	}
	if (has_subject) {
		refs.subject_schema = ValidateSequenceTableSchema(context, refs.subject_table, /*allow_bigint=*/true);
	} else {
		refs.ref_paths = ParseSortMeRNARefPaths(paths_it->second, fn_name);
	}
	return refs;
}

// InitGlobal half: build the index from whichever source Bind recorded.
inline std::unique_ptr<miint::SortMeRNAAligner>
MakeSortMeRNAAligner(ClientContext &context, const miint::SortMeRNAConfig &cfg, const SortMeRNAReferences &refs) {
	if (refs.subject_table.empty()) {
		return std::make_unique<miint::SortMeRNAAligner>(cfg, refs.ref_paths);
	}
	auto subjects = ReadSubjectTable(context, refs.subject_table, refs.subject_schema);
	return std::make_unique<miint::SortMeRNAAligner>(cfg, subjects);
}

inline void ParseSortMeRNAConfigParams(const named_argument_map_t &params, miint::SortMeRNAConfig &cfg,
                                       const char *fn_name) {
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
	// A negative evalue is the aligner's "keep every hit" sentinel; from SQL it
	// (or NaN) would be a filter that silently filters nothing.
	auto evalue_it = params.find("evalue");
	if (evalue_it != params.end() && !evalue_it->second.IsNull() && !(cfg.evalue >= 0.0)) {
		throw BinderException("%s: evalue must be a non-negative number (got %s)", fn_name,
		                      evalue_it->second.ToString());
	}
	set_u32("seed_win_len", cfg.seed_win_len);
	set_u32("num_alignments", cfg.num_alignments);
	set_bool("best", cfg.best);
	set_bool("paired", cfg.paired);
	set_bool("forward_only", cfg.forward_only);
	set_bool("reverse_only", cfg.reverse_only);
	set_bool("full_search", cfg.full_search);
}

// Register the shared LIST(VARCHAR)/VARCHAR/INTEGER/DOUBLE/BOOLEAN named
// parameters on a sortmerna-family TableFunction.
inline void RegisterSortMeRNANamedParameters(TableFunction &tf) {
	AddNamedParameter(tf, "ref_paths", LogicalType::LIST(LogicalType::VARCHAR));
	AddNamedParameter(tf, "subject_table", LogicalType::VARCHAR);
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
