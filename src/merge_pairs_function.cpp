#include "merge_pairs_function.hpp"
#include "QualScore.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/types/string_type.hpp"
#include "duckdb/common/vector_operations/generic_executor.hpp"
#include "duckdb/execution/expression_executor.hpp"
#include "duckdb/function/scalar_function.hpp"
#include "duckdb/main/extension/extension_loader.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"

#include "table_function_common.hpp"
#include "vsearch_api.h"
#include "fastq_mergepairs.h"

#include <mutex>

namespace duckdb {

// Maximum merged length supported by vsearch's fixed-size buffer
static constexpr int VSEARCH_MAX_MERGED_LEN = 9999;

// ---------------------------------------------------------------------------
// Session: init-once, shared across threads, released when last ref drops
// ---------------------------------------------------------------------------
struct MergePairsSession {
	std::once_flag init_flag;
	bool initialized = false;

	void EnsureInit(int minovlen, int maxdiffs, double maxdiffpct, double maxee, int minlen, int maxlen) {
		std::call_once(init_flag, [&]() {
			vsearch_init_defaults();
			initialized = true;

			opt_fastq_minovlen = minovlen;
			opt_fastq_maxdiffs = maxdiffs;
			opt_fastq_maxdiffpct = maxdiffpct;
			opt_fastq_maxee = maxee;
			opt_fastq_minlen = minlen;
			opt_fastq_maxlen = maxlen;
			opt_fastq_ascii = 33;
			opt_fastq_qmin = 0;
			opt_fastq_qmax = 41;
			opt_fastq_qmaxout = 41;
			opt_fastq_qminout = 0;

			vsearch_apply_defaults_fixups();
			mergepairs_init();
		});
	}

	~MergePairsSession() {
		if (initialized) {
			vsearch_session_end();
		}
	}

	MergePairsSession() = default;
	MergePairsSession(const MergePairsSession &) = delete;
	MergePairsSession &operator=(const MergePairsSession &) = delete;
};

// ---------------------------------------------------------------------------
// Bind data: stores merge parameters + shared session
// ---------------------------------------------------------------------------
struct MergePairsBindData : public FunctionData {
	int minovlen = 10;
	int maxdiffs = 10;
	double maxdiffpct = 100.0;
	double maxee = 1e18; // effectively unlimited
	int minlen = 1;
	int maxlen = 1000000; // effectively unlimited

	// Shared session: init-once, released when last BindData copy is destroyed
	std::shared_ptr<MergePairsSession> session;

	MergePairsBindData() : session(std::make_shared<MergePairsSession>()) {
	}

	unique_ptr<FunctionData> Copy() const override {
		auto copy = make_uniq<MergePairsBindData>();
		copy->minovlen = minovlen;
		copy->maxdiffs = maxdiffs;
		copy->maxdiffpct = maxdiffpct;
		copy->maxee = maxee;
		copy->minlen = minlen;
		copy->maxlen = maxlen;
		copy->session = session; // shared ownership
		return copy;
	}

	bool Equals(const FunctionData &other_p) const override {
		auto &other = other_p.Cast<MergePairsBindData>();
		return minovlen == other.minovlen && maxdiffs == other.maxdiffs && maxdiffpct == other.maxdiffpct &&
		       maxee == other.maxee && minlen == other.minlen && maxlen == other.maxlen;
	}

	static unique_ptr<MergePairsBindData> Defaults() {
		return make_uniq<MergePairsBindData>();
	}

	// Extract and validate parameters from the 10-arg overload (args 4-9)
	static unique_ptr<MergePairsBindData> FromArgs10(ClientContext &context,
	                                                 vector<unique_ptr<Expression>> &arguments) {
		for (idx_t i = 4; i < 10; i++) {
			if (!arguments[i]->IsFoldable()) {
				throw InvalidInputException(
				    "merge_pairs tuning parameters must be constant values, not column references");
			}
		}
		auto data = make_uniq<MergePairsBindData>();
		data->minovlen = ExpressionExecutor::EvaluateScalar(context, *arguments[4]).GetValue<int>();
		data->maxdiffs = ExpressionExecutor::EvaluateScalar(context, *arguments[5]).GetValue<int>();
		data->maxdiffpct = ExpressionExecutor::EvaluateScalar(context, *arguments[6]).GetValue<double>();
		data->maxee = ExpressionExecutor::EvaluateScalar(context, *arguments[7]).GetValue<double>();
		data->minlen = ExpressionExecutor::EvaluateScalar(context, *arguments[8]).GetValue<int>();
		data->maxlen = ExpressionExecutor::EvaluateScalar(context, *arguments[9]).GetValue<int>();

		if (data->minovlen < 1) {
			throw InvalidInputException("merge_pairs_vsearch: minovlen must be >= 1 (got %d)", data->minovlen);
		}
		if (data->maxdiffs < 0) {
			throw InvalidInputException("merge_pairs_vsearch: maxdiffs must be >= 0 (got %d)", data->maxdiffs);
		}
		if (data->maxdiffpct < 0.0 || data->maxdiffpct > 100.0) {
			throw InvalidInputException("merge_pairs_vsearch: maxdiffpct must be 0.0-100.0 (got %g)", data->maxdiffpct);
		}
		if (data->maxee <= 0.0) {
			throw InvalidInputException("merge_pairs_vsearch: maxee must be > 0.0 (got %g)", data->maxee);
		}
		if (data->minlen < 1) {
			throw InvalidInputException("merge_pairs_vsearch: minlen must be >= 1 (got %d)", data->minlen);
		}
		if (data->maxlen < 1) {
			throw InvalidInputException("merge_pairs_vsearch: maxlen must be >= 1 (got %d)", data->maxlen);
		}
		if (data->minlen > data->maxlen) {
			throw InvalidInputException("merge_pairs_vsearch: minlen (%d) must be <= maxlen (%d)", data->minlen,
			                            data->maxlen);
		}
		return data;
	}
};

// ---------------------------------------------------------------------------
// Local state: per-thread reusable buffers (no session management)
// ---------------------------------------------------------------------------
struct MergePairsLocalState : public FunctionLocalState {
	std::string fwd_qual_ascii;
	std::string rev_qual_ascii;
};

static unique_ptr<FunctionLocalState>
MergePairsInitLocalState(ExpressionState &state, const BoundFunctionExpression &expr, FunctionData *bind_data_p) {
	auto &data = bind_data_p->Cast<MergePairsBindData>();

	// Ensure vsearch session is initialized (init-once across all threads)
	data.session->EnsureInit(data.minovlen, data.maxdiffs, data.maxdiffpct, data.maxee, data.minlen, data.maxlen);

	return make_uniq<MergePairsLocalState>();
}

// ---------------------------------------------------------------------------
// Shared type helpers
// ---------------------------------------------------------------------------
static LogicalType QualType() {
	return LogicalType::LIST(LogicalType::UTINYINT);
}

// 4-arg: fwd_seq, fwd_qual, rev_seq, rev_qual
static vector<LogicalType> FourArgTypes() {
	return {LogicalType::VARCHAR, QualType(), LogicalType::VARCHAR, QualType()};
}

// 10-arg: 4 inputs + minovlen, maxdiffs, maxdiffpct, maxee, minlen, maxlen
static vector<LogicalType> TenArgTypes() {
	return {LogicalType::VARCHAR, QualType(),           LogicalType::VARCHAR, QualType(),
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::DOUBLE,  LogicalType::DOUBLE,
	        LogicalType::INTEGER, LogicalType::INTEGER};
}

static LogicalType MergePairsReturnType() {
	return LogicalType::STRUCT({{"merged", LogicalType::BOOLEAN},
	                            {"sequence", LogicalType::VARCHAR},
	                            {"quality", QualType()},
	                            {"ee_merged", LogicalType::DOUBLE},
	                            {"ee_fwd", LogicalType::DOUBLE},
	                            {"ee_rev", LogicalType::DOUBLE},
	                            {"fwd_errors", LogicalType::INTEGER},
	                            {"rev_errors", LogicalType::INTEGER},
	                            {"overlap", LogicalType::INTEGER}});
}

// ---------------------------------------------------------------------------
// Execute (shared by both overloads — first 4 args are the same)
// ---------------------------------------------------------------------------
static void MergePairsExecute(DataChunk &args, ExpressionState &state, Vector &result) {
	auto &lstate = ExecuteFunctionState::GetFunctionState(state)->Cast<MergePairsLocalState>();

	// Input vectors (first 4 args only)
	UnifiedVectorFormat fwd_seq_data, fwd_qual_data, rev_seq_data, rev_qual_data;
	args.data[0].ToUnifiedFormat(args.size(), fwd_seq_data);
	args.data[1].ToUnifiedFormat(args.size(), fwd_qual_data);
	args.data[2].ToUnifiedFormat(args.size(), rev_seq_data);
	args.data[3].ToUnifiedFormat(args.size(), rev_qual_data);

	auto fwd_seq_ptr = UnifiedVectorFormat::GetData<string_t>(fwd_seq_data);
	auto rev_seq_ptr = UnifiedVectorFormat::GetData<string_t>(rev_seq_data);

	// Output struct entries
	auto &entries = StructVector::GetEntries(result);
	auto merged_data = FlatVector::GetData<bool>(*entries[0]);
	auto &seq_vec = *entries[1];
	auto &qual_list_vec = *entries[2]; // LIST(UTINYINT)
	auto ee_merged_data = FlatVector::GetData<double>(*entries[3]);
	auto ee_fwd_data = FlatVector::GetData<double>(*entries[4]);
	auto ee_rev_data = FlatVector::GetData<double>(*entries[5]);
	auto fwd_errors_data = FlatVector::GetData<int32_t>(*entries[6]);
	auto rev_errors_data = FlatVector::GetData<int32_t>(*entries[7]);
	auto overlap_data = FlatVector::GetData<int32_t>(*entries[8]);
	auto &result_validity = FlatVector::Validity(result);

	// Quality output: LIST(UTINYINT) managed via ListVector
	auto qual_list_entries = FlatVector::GetData<list_entry_t>(qual_list_vec);
	idx_t qual_child_offset = ListVector::GetListSize(qual_list_vec);

	// Reusable buffers
	std::string fwd_seq_buf, rev_seq_buf;

	for (idx_t i = 0; i < args.size(); i++) {
		auto fi = fwd_seq_data.sel->get_index(i);
		auto fqi = fwd_qual_data.sel->get_index(i);
		auto rsi = rev_seq_data.sel->get_index(i);
		auto rqi = rev_qual_data.sel->get_index(i);

		// NULL if any input is NULL
		if (!fwd_seq_data.validity.RowIsValid(fi) || !fwd_qual_data.validity.RowIsValid(fqi) ||
		    !rev_seq_data.validity.RowIsValid(rsi) || !rev_qual_data.validity.RowIsValid(rqi)) {
			result_validity.SetInvalid(i);
			qual_list_entries[i].offset = qual_child_offset;
			qual_list_entries[i].length = 0;
			continue;
		}

		fwd_seq_buf.assign(fwd_seq_ptr[fi].GetData(), fwd_seq_ptr[fi].GetSize());
		rev_seq_buf.assign(rev_seq_ptr[rsi].GetData(), rev_seq_ptr[rsi].GetSize());

		// Guard against buffer overflow in vsearch's fixed-size merge_result_s
		if (static_cast<int>(fwd_seq_buf.size() + rev_seq_buf.size()) > VSEARCH_MAX_MERGED_LEN) {
			throw InvalidInputException(
			    "merge_pairs_vsearch: combined input length (%llu + %llu = %llu) exceeds maximum merged length (%d). "
			    "This function is designed for short-read paired-end data.",
			    fwd_seq_buf.size(), rev_seq_buf.size(), fwd_seq_buf.size() + rev_seq_buf.size(),
			    VSEARCH_MAX_MERGED_LEN);
		}

		// Convert LIST(UTINYINT) Phred scores to ASCII quality strings
		const uint8_t *fwd_qual_ptr;
		const uint8_t *rev_qual_ptr;
		idx_t fwd_qual_len, rev_qual_len;
		GetListUInt8Slice(args.data[1], fwd_qual_data, i, fwd_qual_ptr, fwd_qual_len);
		GetListUInt8Slice(args.data[3], rev_qual_data, i, rev_qual_ptr, rev_qual_len);

		lstate.fwd_qual_ascii = miint::encode_quality_ascii(fwd_qual_ptr, fwd_qual_len);
		lstate.rev_qual_ascii = miint::encode_quality_ascii(rev_qual_ptr, rev_qual_len);

		// Call vsearch merge
		merge_result_s mr {};
		mergepairs_single(fwd_seq_buf.c_str(), lstate.fwd_qual_ascii.c_str(), static_cast<int>(fwd_seq_buf.size()),
		                  rev_seq_buf.c_str(), lstate.rev_qual_ascii.c_str(), static_cast<int>(rev_seq_buf.size()), "",
		                  "", &mr);

		merged_data[i] = mr.merged;
		if (mr.merged) {
			FlatVector::GetData<string_t>(seq_vec)[i] =
			    StringVector::AddString(seq_vec, mr.merged_sequence, mr.merged_length);

			// Write merged quality as LIST(UTINYINT) using QualScore::write_decoded
			ListVector::Reserve(qual_list_vec, qual_child_offset + mr.merged_length);
			auto &qual_child = ListVector::GetEntry(qual_list_vec);
			auto qual_child_data = FlatVector::GetData<uint8_t>(qual_child);
			qual_list_entries[i].offset = qual_child_offset;
			qual_list_entries[i].length = mr.merged_length;
			miint::QualScore merged_qual(std::string(mr.merged_quality, mr.merged_length));
			merged_qual.write_decoded(qual_child_data + qual_child_offset);
			qual_child_offset += mr.merged_length;
			ListVector::SetListSize(qual_list_vec, qual_child_offset);

			ee_merged_data[i] = mr.ee_merged;
			ee_fwd_data[i] = mr.ee_fwd;
			ee_rev_data[i] = mr.ee_rev;
			fwd_errors_data[i] = mr.fwd_errors;
			rev_errors_data[i] = mr.rev_errors;
			overlap_data[i] = mr.overlap_length;
		} else {
			FlatVector::GetData<string_t>(seq_vec)[i] = StringVector::AddString(seq_vec, "", 0);
			qual_list_entries[i].offset = qual_child_offset;
			qual_list_entries[i].length = 0;
			ee_merged_data[i] = 0.0;
			ee_fwd_data[i] = 0.0;
			ee_rev_data[i] = 0.0;
			fwd_errors_data[i] = 0;
			rev_errors_data[i] = 0;
			overlap_data[i] = 0;
		}
	}
}

// ---------------------------------------------------------------------------
// Register
// ---------------------------------------------------------------------------
void MergePairsFunction::Register(ExtensionLoader &loader) {
	ScalarFunctionSet function_set("merge_pairs_vsearch");

	// 4-arg: merge_pairs(fwd_seq, fwd_qual, rev_seq, rev_qual)
	ScalarFunction merge_4arg("merge_pairs_vsearch", FourArgTypes(), MergePairsReturnType(), MergePairsExecute);
	merge_4arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	merge_4arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = MergePairsReturnType();
		return unique_ptr<FunctionData>(MergePairsBindData::Defaults().release());
	};
	merge_4arg.init_local_state = MergePairsInitLocalState;
	function_set.AddFunction(merge_4arg);

	// 10-arg: merge_pairs(fwd_seq, fwd_qual, rev_seq, rev_qual,
	//                      minovlen, maxdiffs, maxdiffpct, maxee, minlen, maxlen)
	ScalarFunction merge_10arg("merge_pairs_vsearch", TenArgTypes(), MergePairsReturnType(), MergePairsExecute);
	merge_10arg.null_handling = FunctionNullHandling::SPECIAL_HANDLING;
	merge_10arg.bind = [](ClientContext &ctx, ScalarFunction &fn, vector<unique_ptr<Expression>> &args) {
		fn.return_type = MergePairsReturnType();
		return unique_ptr<FunctionData>(MergePairsBindData::FromArgs10(ctx, args).release());
	};
	merge_10arg.init_local_state = MergePairsInitLocalState;
	function_set.AddFunction(merge_10arg);

	loader.RegisterFunction(function_set);
}

} // namespace duckdb
