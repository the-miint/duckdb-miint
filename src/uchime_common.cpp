#include "uchime_common.hpp"

#include <algorithm>

namespace duckdb {

std::vector<std::string> GetUchimeOutputNames() {
	return {"score",        "query",      "parent_a", "parent_b",      "closest_parent", "id_query_model",
	        "id_query_a",   "id_query_b", "id_a_b",   "id_query_top",  "left_yes",       "left_no",
	        "left_abstain", "right_yes",  "right_no", "right_abstain", "divergence",     "flag"};
}

std::vector<LogicalType> GetUchimeOutputTypes() {
	return {LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	        LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
	        LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::DOUBLE,  LogicalType::VARCHAR};
}

idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset,
                          idx_t count) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = 0;

	auto score_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		score_data[i] = results[offset + i].score;
	}

	auto &query_vec = output.data[col++];
	auto &parent_a_vec = output.data[col++];
	auto &parent_b_vec = output.data[col++];
	auto &closest_parent_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		FlatVector::GetData<string_t>(query_vec)[i] = StringVector::AddString(query_vec, r.query_label);
		FlatVector::GetData<string_t>(parent_a_vec)[i] = StringVector::AddString(parent_a_vec, r.parent_a_label);
		FlatVector::GetData<string_t>(parent_b_vec)[i] = StringVector::AddString(parent_b_vec, r.parent_b_label);
		FlatVector::GetData<string_t>(closest_parent_vec)[i] =
		    StringVector::AddString(closest_parent_vec, r.closest_parent_label);
	}

	auto id_qm = FlatVector::GetData<double>(output.data[col++]);
	auto id_qa = FlatVector::GetData<double>(output.data[col++]);
	auto id_qb = FlatVector::GetData<double>(output.data[col++]);
	auto id_ab = FlatVector::GetData<double>(output.data[col++]);
	auto id_qt = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		id_qm[i] = r.id_query_model;
		id_qa[i] = r.id_query_a;
		id_qb[i] = r.id_query_b;
		id_ab[i] = r.id_a_b;
		id_qt[i] = r.id_query_top;
	}

	auto ly = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ln = FlatVector::GetData<int32_t>(output.data[col++]);
	auto la = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ry = FlatVector::GetData<int32_t>(output.data[col++]);
	auto rn = FlatVector::GetData<int32_t>(output.data[col++]);
	auto ra = FlatVector::GetData<int32_t>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		ly[i] = r.left_yes;
		ln[i] = r.left_no;
		la[i] = r.left_abstain;
		ry[i] = r.right_yes;
		rn[i] = r.right_no;
		ra[i] = r.right_abstain;
	}

	auto div_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		div_data[i] = results[offset + i].divergence;
	}

	auto &flag_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(flag_vec)[i] = StringVector::AddString(flag_vec, results[offset + i].flag);
	}

	D_ASSERT(col == output.ColumnCount());
	output.SetCardinality(actual);
	return actual;
}

} // namespace duckdb
