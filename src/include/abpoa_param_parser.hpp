#pragma once

#include "AbpoaAligner.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/function/table_function.hpp"

namespace duckdb {

inline miint::AbpoaAlignParams ParseAbpoaParams(const named_parameter_map_t &named_params) {
	miint::AbpoaAlignParams params;
	for (auto &kv : named_params) {
		if (kv.first == "sample_id") {
			continue;
		} else if (kv.first == "match") {
			params.match = kv.second.GetValue<int>();
		} else if (kv.first == "mismatch") {
			params.mismatch = kv.second.GetValue<int>();
		} else if (kv.first == "gap_open1") {
			params.gap_open1 = kv.second.GetValue<int>();
		} else if (kv.first == "gap_open2") {
			params.gap_open2 = kv.second.GetValue<int>();
		} else if (kv.first == "gap_ext1") {
			params.gap_ext1 = kv.second.GetValue<int>();
		} else if (kv.first == "gap_ext2") {
			params.gap_ext2 = kv.second.GetValue<int>();
		} else if (kv.first == "align_mode") {
			params.align_mode = kv.second.GetValue<string>();
		} else if (kv.first == "progressive") {
			params.progressive = kv.second.GetValue<bool>();
		} else if (kv.first == "disable_seeding") {
			params.disable_seeding = kv.second.GetValue<bool>();
		} else if (kv.first == "amb_strand") {
			params.amb_strand = kv.second.GetValue<bool>();
		} else if (kv.first == "k") {
			params.k = kv.second.GetValue<int>();
		} else if (kv.first == "w") {
			params.w = kv.second.GetValue<int>();
		} else if (kv.first == "min_w") {
			params.min_w = kv.second.GetValue<int>();
		} else if (kv.first == "bandwidth") {
			params.bandwidth = kv.second.GetValue<int>();
		} else if (kv.first == "bandwidth_frac") {
			params.bandwidth_frac = kv.second.GetValue<float>();
		} else if (kv.first == "max_num_cons") {
			params.max_num_cons = kv.second.GetValue<int>();
		} else if (kv.first == "min_freq") {
			params.min_freq = kv.second.GetValue<float>();
		} else if (kv.first == "algorithm") {
			params.algorithm = kv.second.GetValue<string>();
		}
	}
	return params;
}

} // namespace duckdb
