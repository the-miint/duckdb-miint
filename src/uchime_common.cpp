#include "uchime_common.hpp"

#include <algorithm>

namespace duckdb {

std::vector<std::string> GetUchimeOutputNames() {
	return {"score",        "query_id",   "parent_a_id", "parent_b_id",   "closest_parent_id", "id_query_model",
	        "id_query_a",   "id_query_b", "id_a_b",      "id_query_top",  "left_yes",          "left_no",
	        "left_abstain", "right_yes",  "right_no",    "right_abstain", "divergence",        "flag"};
}

std::vector<LogicalType> GetUchimeOutputTypes() {
	return {LogicalType::DOUBLE,  LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR,
	        LogicalType::VARCHAR, LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::DOUBLE,
	        LogicalType::DOUBLE,  LogicalType::DOUBLE,  LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER, LogicalType::INTEGER,
	        LogicalType::DOUBLE,  LogicalType::VARCHAR};
}

idx_t OutputUchimeResults(DataChunk &output, const std::vector<miint::UchimeResult> &results, idx_t offset, idx_t count,
                          idx_t start_col) {
	idx_t actual = std::min(count, static_cast<idx_t>(results.size()) - offset);
	if (actual == 0) {
		output.SetCardinality(0);
		return 0;
	}

	idx_t col = start_col;

	// score — always populated
	auto score_data = FlatVector::GetData<double>(output.data[col++]);
	for (idx_t i = 0; i < actual; i++) {
		score_data[i] = results[offset + i].score;
	}

	// query_id — always populated
	auto &query_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(query_vec)[i] =
		    StringVector::AddString(query_vec, results[offset + i].query_label);
	}

	// parent_a_id, parent_b_id, closest_parent_id — NULL when non-chimeric (empty label)
	auto &parent_a_vec = output.data[col++];
	auto &parent_b_vec = output.data[col++];
	auto &closest_parent_vec = output.data[col++];
	auto &pa_validity = FlatVector::Validity(parent_a_vec);
	auto &pb_validity = FlatVector::Validity(parent_b_vec);
	auto &cp_validity = FlatVector::Validity(closest_parent_vec);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		if (r.parent_a_label.empty()) {
			pa_validity.SetInvalid(i);
			pb_validity.SetInvalid(i);
			cp_validity.SetInvalid(i);
		} else {
			FlatVector::GetData<string_t>(parent_a_vec)[i] = StringVector::AddString(parent_a_vec, r.parent_a_label);
			FlatVector::GetData<string_t>(parent_b_vec)[i] = StringVector::AddString(parent_b_vec, r.parent_b_label);
			FlatVector::GetData<string_t>(closest_parent_vec)[i] =
			    StringVector::AddString(closest_parent_vec, r.closest_parent_label);
		}
	}

	// id_query_model, id_query_a, id_query_b, id_a_b, id_query_top — NULL when non-chimeric
	idx_t id_qm_col = col++;
	idx_t id_qa_col = col++;
	idx_t id_qb_col = col++;
	idx_t id_ab_col = col++;
	idx_t id_qt_col = col++;
	auto id_qm = FlatVector::GetData<double>(output.data[id_qm_col]);
	auto id_qa = FlatVector::GetData<double>(output.data[id_qa_col]);
	auto id_qb = FlatVector::GetData<double>(output.data[id_qb_col]);
	auto id_ab = FlatVector::GetData<double>(output.data[id_ab_col]);
	auto id_qt = FlatVector::GetData<double>(output.data[id_qt_col]);
	auto &id_qm_v = FlatVector::Validity(output.data[id_qm_col]);
	auto &id_qa_v = FlatVector::Validity(output.data[id_qa_col]);
	auto &id_qb_v = FlatVector::Validity(output.data[id_qb_col]);
	auto &id_ab_v = FlatVector::Validity(output.data[id_ab_col]);
	auto &id_qt_v = FlatVector::Validity(output.data[id_qt_col]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		if (r.flag == "N") {
			id_qm_v.SetInvalid(i);
			id_qa_v.SetInvalid(i);
			id_qb_v.SetInvalid(i);
			id_ab_v.SetInvalid(i);
			id_qt_v.SetInvalid(i);
		} else {
			id_qm[i] = r.id_query_model;
			id_qa[i] = r.id_query_a;
			id_qb[i] = r.id_query_b;
			id_ab[i] = r.id_a_b;
			id_qt[i] = r.id_query_top;
		}
	}

	// left_yes, left_no, left_abstain, right_yes, right_no, right_abstain
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

	// divergence — NULL when non-chimeric
	idx_t div_col = col++;
	auto div_data = FlatVector::GetData<double>(output.data[div_col]);
	auto &div_v = FlatVector::Validity(output.data[div_col]);
	for (idx_t i = 0; i < actual; i++) {
		auto &r = results[offset + i];
		if (r.flag == "N") {
			div_v.SetInvalid(i);
		} else {
			div_data[i] = r.divergence;
		}
	}

	// flag — always populated
	auto &flag_vec = output.data[col++];
	for (idx_t i = 0; i < actual; i++) {
		FlatVector::GetData<string_t>(flag_vec)[i] = StringVector::AddString(flag_vec, results[offset + i].flag);
	}

	// Exactly 18 uchimeout columns should have been written, starting at start_col.
	// (output.ColumnCount() equals start_col + 18 when the caller prepends columns,
	// so the two forms are equivalent — but this form stays honest about what
	// this function owns regardless of the caller's extra columns.)
	D_ASSERT(col == start_col + 18);
	output.SetCardinality(actual);
	return actual;
}

} // namespace duckdb
