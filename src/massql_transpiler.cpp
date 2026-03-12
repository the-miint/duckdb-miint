#include "massql_transpiler.hpp"

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>

namespace miint {

// ── Helpers ──────────────────────────────────────────────────────────────────

static std::string to_upper(const std::string &s) {
	std::string result = s;
	std::transform(result.begin(), result.end(), result.begin(), ::toupper);
	return result;
}

static std::string fmt_double(double v) {
	std::ostringstream ss;
	ss << std::setprecision(std::numeric_limits<double>::digits10) << v;
	return ss.str();
}

bool MassQLTranspiler::is_file_path(const std::string &source) {
	// Check for .mzml extension at end of string (case insensitive)
	if (source.size() >= 5) {
		auto ext = to_upper(source.substr(source.size() - 5));
		if (ext == ".MZML") {
			return true;
		}
	}
	return source.find('/') != std::string::npos;
}

// ── Condition field classification ───────────────────────────────────────────

static bool is_peak_level_field(ConditionField field) {
	switch (field) {
	case ConditionField::MS1MZ:
	case ConditionField::MS2PROD:
	case ConditionField::MS2NL:
		return true;
	default:
		return false;
	}
}

static bool is_metadata_field(ConditionField field) {
	switch (field) {
	case ConditionField::RTMIN:
	case ConditionField::RTMAX:
	case ConditionField::SCANMIN:
	case ConditionField::SCANMAX:
	case ConditionField::CHARGE:
	case ConditionField::POLARITY:
		return true;
	default:
		return false;
	}
}

static std::string peak_column(ConditionField field) {
	switch (field) {
	case ConditionField::MS1MZ:
	case ConditionField::MS2PROD:
		return "mz";
	case ConditionField::MS2NL:
		return "neutral_loss";
	default:
		throw std::runtime_error("MassQL transpiler: not a peak-level field");
	}
}

// ── Qualifier helpers ────────────────────────────────────────────────────────

static bool has_qualifier(const std::vector<Qualifier> &qualifiers, const std::string &name) {
	for (const auto &q : qualifiers) {
		if (q.name == name) {
			return true;
		}
	}
	return false;
}

static double get_qualifier_value(const std::vector<Qualifier> &qualifiers, const std::string &name) {
	for (const auto &q : qualifiers) {
		if (q.name == name) {
			return q.value;
		}
	}
	return 0.0;
}

// ── Tolerance expression using mz_within / mz_within_ppm macros ─────────────

static std::string tolerance_expr(const std::string &col, double target, const std::vector<Qualifier> &qualifiers) {
	for (const auto &q : qualifiers) {
		if (q.name == "TOLERANCEPPM") {
			return "mz_within_ppm(" + col + ", " + fmt_double(target) + ", " + fmt_double(q.value) + ")";
		}
		if (q.name == "TOLERANCEMZ") {
			return "mz_within(" + col + ", " + fmt_double(target) + ", " + fmt_double(q.value) + ")";
		}
	}
	// Default: 0.1 Da (matches Python MassQL)
	return "mz_within(" + col + ", " + fmt_double(target) + ", 0.1)";
}

static QualifierOp get_qualifier_op(const std::vector<Qualifier> &qualifiers, const std::string &name) {
	for (const auto &q : qualifiers) {
		if (q.name == name) {
			return q.op;
		}
	}
	return QualifierOp::NONE;
}

// Peak match with optional intensity qualifiers
static std::string peak_match_expr(const std::string &col, double target, const Condition &cond) {
	std::string expr = tolerance_expr(col, target, cond.qualifiers);

	if (has_qualifier(cond.qualifiers, "INTENSITYPERCENT")) {
		double pct = get_qualifier_value(cond.qualifiers, "INTENSITYPERCENT");
		auto op = get_qualifier_op(cond.qualifiers, "INTENSITYPERCENT");
		std::string cmp = (op == QualifierOp::GREATER_THAN) ? " > " : (op == QualifierOp::LESS_THAN) ? " < " : " >= ";
		expr += " AND i_norm" + cmp + fmt_double(pct / 100.0);
	}
	if (has_qualifier(cond.qualifiers, "INTENSITYVALUE")) {
		double val = get_qualifier_value(cond.qualifiers, "INTENSITYVALUE");
		auto op = get_qualifier_op(cond.qualifiers, "INTENSITYVALUE");
		std::string cmp = (op == QualifierOp::GREATER_THAN) ? " > " : (op == QualifierOp::LESS_THAN) ? " < " : " >= ";
		expr += " AND intensity" + cmp + fmt_double(val);
	}

	return expr;
}

// ── X-variable helpers ───────────────────────────────────────────────────────

static bool condition_has_x_variable(const Condition &cond) {
	for (const auto &v : cond.values) {
		if (v.has_x_variable) {
			return true;
		}
	}
	return false;
}

// Extract Da tolerance from X-conditions (PPM not supported for offset patterns)
static std::string x_offset_tolerance(const std::vector<const Condition *> &x_conds) {
	for (const auto *cond : x_conds) {
		for (const auto &q : cond->qualifiers) {
			if (q.name == "TOLERANCEPPM") {
				throw std::runtime_error(
				    "MassQL transpiler: PPM tolerance not supported for X-variable offset patterns");
			}
			if (q.name == "TOLERANCEMZ") {
				return fmt_double(q.value);
			}
		}
	}
	return "0.1";
}

// Generate CTEs for X-variable offset pattern (mirrors mzml_x_offset_ntuple macro)
static std::string x_offset_ctes(const std::vector<const Condition *> &x_conds, const std::string &base_cte,
                                 std::string &qualifying_cte_name) {
	// Extract and validate offsets
	std::vector<double> offsets;
	for (const auto *cond : x_conds) {
		if (cond->values.size() != 1 || !cond->values[0].has_x_variable) {
			throw std::runtime_error("MassQL transpiler: X-variable conditions must have exactly one X value");
		}
		if (cond->values[0].x_coefficient != 1.0) {
			throw std::runtime_error("MassQL transpiler: non-linear X coefficients not yet supported");
		}
		offsets.push_back(cond->values[0].constant_value);
	}

	// Normalize: subtract minimum so smallest offset is 0
	double min_offset = *std::min_element(offsets.begin(), offsets.end());
	for (auto &o : offsets) {
		o -= min_offset;
	}

	std::string tol = x_offset_tolerance(x_conds);
	qualifying_cte_name = "__x_qualifying";

	std::ostringstream ss;
	// Recursive CTE: enumerate distinct mz values as X candidates
	ss << "__x_candidates(x_val, next_min) AS (\n"
	   << "  (SELECT mz, mz + 0.05 FROM " << base_cte << " WHERE intensity > 0 ORDER BY mz LIMIT 1)\n"
	   << "  UNION ALL\n"
	   << "  (SELECT s.mz, s.mz + 0.05\n"
	   << "   FROM __x_candidates g\n"
	   << "   JOIN (SELECT DISTINCT mz FROM " << base_cte << " WHERE intensity > 0) s ON s.mz >= g.next_min\n"
	   << "   ORDER BY s.mz\n"
	   << "   LIMIT 1)\n"
	   << "),\n";

	// Match candidates against all offsets
	ss << "__x_matched AS (\n"
	   << "  SELECT xc.x_val, p.spectrum_index, ot.offset_val\n"
	   << "  FROM __x_candidates xc\n"
	   << "  CROSS JOIN (VALUES ";
	for (size_t i = 0; i < offsets.size(); i++) {
		if (i > 0) {
			ss << ", ";
		}
		ss << "(" << fmt_double(offsets[i]) << ")";
	}
	ss << ") AS ot(offset_val)\n"
	   << "  JOIN " << base_cte << " p ON p.mz > xc.x_val + ot.offset_val - " << tol << "\n"
	   << "              AND p.mz < xc.x_val + ot.offset_val + " << tol << "\n"
	   << "              AND p.intensity > 0\n"
	   << "),\n";

	// Qualifying: spectra where ALL offsets matched for some X value
	ss << qualifying_cte_name << " AS (\n"
	   << "  SELECT spectrum_index FROM __x_matched\n"
	   << "  GROUP BY spectrum_index, x_val\n"
	   << "  HAVING COUNT(DISTINCT offset_val) = " << offsets.size() << "\n"
	   << ")";

	return ss.str();
}

// Generate CTEs for peak-pair X pattern (mirrors mzml_peak_pair macro)
// Pattern: MS2PROD=X AND MS2PROD=2*(X-formula(Fe)) — N-way self-join with arbitrary coefficients
static std::string x_peak_pair_ctes(const std::vector<const Condition *> &x_conds, const std::string &base_cte,
                                    std::string &qualifying_cte_name) {
	// Validate each condition
	for (const auto *cond : x_conds) {
		if (cond->values.size() != 1 || !cond->values[0].has_x_variable) {
			throw std::runtime_error("MassQL transpiler: X-variable conditions must have exactly one X value");
		}
	}

	std::string tol = x_offset_tolerance(x_conds);
	qualifying_cte_name = "__x_qualifying";

	std::ostringstream ss;
	// Recursive CTE: enumerate distinct mz values as X candidates
	ss << "__x_candidates(x_val, next_min) AS (\n"
	   << "  (SELECT mz, mz + 0.05 FROM " << base_cte << " WHERE intensity > 0 ORDER BY mz LIMIT 1)\n"
	   << "  UNION ALL\n"
	   << "  (SELECT s.mz, s.mz + 0.05\n"
	   << "   FROM __x_candidates g\n"
	   << "   JOIN (SELECT DISTINCT mz FROM " << base_cte << " WHERE intensity > 0) s ON s.mz >= g.next_min\n"
	   << "   ORDER BY s.mz\n"
	   << "   LIMIT 1)\n"
	   << "),\n";

	// N-way self-join: one JOIN per X-condition
	ss << qualifying_cte_name << " AS (\n"
	   << "  SELECT DISTINCT p1.spectrum_index\n"
	   << "  FROM __x_candidates xc\n";

	for (size_t i = 0; i < x_conds.size(); i++) {
		std::string alias = "p" + std::to_string(i + 1);
		double coeff = x_conds[i]->values[0].x_coefficient;
		double offset = x_conds[i]->values[0].constant_value;

		// Build target expression: coeff * xc.x_val + offset
		std::string target;
		if (coeff == 1.0) {
			target = "xc.x_val";
		} else {
			target = fmt_double(coeff) + " * xc.x_val";
		}
		if (offset > 0.0) {
			target += " + " + fmt_double(offset);
		} else if (offset < 0.0) {
			target += " - " + fmt_double(-offset);
		}

		ss << "  JOIN " << base_cte << " " << alias << " ON ";
		if (i > 0) {
			ss << alias << ".spectrum_index = p1.spectrum_index AND ";
		}
		ss << alias << ".mz > " << target << " - " << tol << "\n"
		   << "              AND " << alias << ".mz < " << target << " + " << tol << "\n"
		   << "              AND " << alias << ".intensity > 0\n";
	}

	ss << ")";
	return ss.str();
}

// Generate CTEs for cross-level X pattern (mirrors mzml_x_ms1_ms2_prec macro)
// Pattern: MS1MZ=X AND MS2PREC=X — MS2 spectra whose precursor matches an MS1 peak
static std::string x_cross_level_ctes(const std::vector<const Condition *> &x_conds, const std::string &source,
                                      std::string &qualifying_cte_name) {
	// Extract tolerance from any X-condition
	std::string tol_value;
	std::string tol_func = "mz_within";
	for (const auto *cond : x_conds) {
		for (const auto &q : cond->qualifiers) {
			if (q.name == "TOLERANCEPPM") {
				tol_value = fmt_double(q.value);
				tol_func = "mz_within_ppm";
				break;
			}
			if (q.name == "TOLERANCEMZ") {
				tol_value = fmt_double(q.value);
				break;
			}
		}
		if (!tol_value.empty()) {
			break;
		}
	}
	if (tol_value.empty()) {
		tol_value = "0.1";
	}

	qualifying_cte_name = "__x_cross_qualifying";

	std::ostringstream ss;
	ss << "__ms1_peaks AS (\n"
	   << "  SELECT DISTINCT spectrum_index, mz FROM mzml_peaks(" << source << ") WHERE ms_level = 1\n"
	   << "),\n";

	ss << qualifying_cte_name << " AS (\n"
	   << "  SELECT DISTINCT b.spectrum_index FROM __base b\n"
	   << "  JOIN __ms1_peaks ms1 ON b.ms1_scan_index = ms1.spectrum_index\n"
	   << "  AND " << tol_func << "(ms1.mz, b.precursor_mz, " << tol_value << ")\n"
	   << ")";

	return ss.str();
}

// ── Scan-level matching CTE for one condition ────────────────────────────────

static std::string scan_match_cte(const Condition &cond, const std::string &base_cte, int idx) {
	std::ostringstream ss;
	std::string cte_name = "__match_" + std::to_string(idx);

	if (cond.field == ConditionField::MS2PREC) {
		ss << cte_name << " AS (\n  SELECT DISTINCT spectrum_index FROM " << base_cte << "\n  WHERE "
		   << tolerance_expr("precursor_mz", cond.values[0].constant_value, cond.qualifiers) << "\n)";
	} else if (is_peak_level_field(cond.field)) {
		auto col = peak_column(cond.field);

		if (cond.values.size() == 1) {
			ss << cte_name << " AS (\n  SELECT DISTINCT spectrum_index FROM " << base_cte << "\n  WHERE "
			   << peak_match_expr(col, cond.values[0].constant_value, cond) << "\n)";
		} else {
			// OR list with optional CARDINALITY
			int min_card = 1;
			int max_card = 0; // 0 means no upper limit
			for (const auto &q : cond.qualifiers) {
				if (q.name == "CARDINALITY") {
					min_card = static_cast<int>(q.value);
					max_card = static_cast<int>(q.max_value);
				}
			}

			ss << cte_name << " AS (\n  SELECT spectrum_index FROM (\n";
			for (size_t i = 0; i < cond.values.size(); i++) {
				if (i > 0) {
					ss << "    UNION ALL\n";
				}
				ss << "    SELECT DISTINCT spectrum_index, " << i << " AS target_idx FROM " << base_cte << "\n"
				   << "    WHERE " << peak_match_expr(col, cond.values[i].constant_value, cond) << "\n";
			}
			ss << "  )\n  GROUP BY spectrum_index\n  HAVING COUNT(DISTINCT target_idx) >= " << min_card;
			if (max_card > 0) {
				ss << " AND COUNT(DISTINCT target_idx) <= " << max_card;
			}
			ss << "\n)";
		}
	} else {
		throw std::runtime_error("MassQL transpiler: metadata condition should not generate scan match CTE");
	}

	return ss.str();
}

// ── Metadata WHERE clause ────────────────────────────────────────────────────

static std::string metadata_where(const Condition &cond) {
	std::ostringstream ss;
	switch (cond.field) {
	case ConditionField::RTMIN:
		ss << "retention_time >= " << fmt_double(cond.values[0].constant_value);
		break;
	case ConditionField::RTMAX:
		ss << "retention_time <= " << fmt_double(cond.values[0].constant_value);
		break;
	case ConditionField::SCANMIN:
		ss << "spectrum_index >= " << static_cast<int64_t>(cond.values[0].constant_value);
		break;
	case ConditionField::SCANMAX:
		ss << "spectrum_index <= " << static_cast<int64_t>(cond.values[0].constant_value);
		break;
	case ConditionField::CHARGE:
		ss << "precursor_charge = " << static_cast<int64_t>(cond.values[0].constant_value);
		break;
	case ConditionField::POLARITY:
		ss << "polarity = '" << cond.string_value << "'";
		break;
	default:
		throw std::runtime_error("MassQL transpiler: not a metadata field");
	}
	return ss.str();
}

// ── Filter expression ────────────────────────────────────────────────────────

static std::string filter_peak_expr(const Condition &cond) {
	if (cond.field == ConditionField::MS2PREC) {
		return tolerance_expr("precursor_mz", cond.values[0].constant_value, cond.qualifiers);
	}
	auto col = peak_column(cond.field);
	return peak_match_expr(col, cond.values[0].constant_value, cond);
}

// ── Aggregation SQL ──────────────────────────────────────────────────────────
// MUST be kept in sync with the corresponding macros in miint_macros.hpp:
//   SCANINFO  → mzml_scaninfo()    SCANSUM    → mzml_scansum()
//   SCANNUM   → mzml_scannum()     SCANMZ     → mzml_scanmz()
//   SCANMAXINT → mzml_scanmaxint()
// The macros use query_table(relation) which requires a table name, not a CTE,
// so we inline the simple GROUP BY / DISTINCT patterns here.
// If the macro column list or semantics change, update this function too.

static std::string aggregation_sql(AggFunction agg, const std::string &source_cte) {
	std::ostringstream ss;
	switch (agg) {
	case AggFunction::NONE:
		ss << "SELECT * FROM " << source_cte;
		break;
	case AggFunction::SCANINFO:
		ss << "SELECT spectrum_index, first(ms_level) AS ms_level, "
		   << "first(retention_time) AS retention_time, first(spectrum_type) AS spectrum_type, "
		   << "first(polarity) AS polarity, SUM(intensity) AS total_ion_current, "
		   << "MAX(intensity) AS base_peak_intensity, MAX(i_norm) AS i_norm, "
		   << "first(precursor_mz) AS precursor_mz, first(precursor_charge) AS precursor_charge, "
		   << "first(precursor_intensity) AS precursor_intensity, first(ms1_scan_index) AS ms1_scan_index "
		   << "FROM " << source_cte << " GROUP BY spectrum_index";
		break;
	case AggFunction::SCANSUM:
		ss << "SELECT spectrum_index, SUM(intensity) AS total_intensity FROM " << source_cte
		   << " GROUP BY spectrum_index";
		break;
	case AggFunction::SCANNUM:
		ss << "SELECT DISTINCT spectrum_index FROM " << source_cte;
		break;
	case AggFunction::SCANMZ:
		ss << "SELECT DISTINCT precursor_mz FROM " << source_cte << " WHERE precursor_mz IS NOT NULL";
		break;
	case AggFunction::SCANMAXINT:
		ss << "SELECT spectrum_index, MAX(intensity) AS max_intensity FROM " << source_cte
		   << " GROUP BY spectrum_index";
		break;
	}
	return ss.str();
}

// ── Main transpiler ──────────────────────────────────────────────────────────

std::string MassQLTranspiler::to_sql(const MassQLQuery &query, const std::string &source) {
	int ms_level = (query.data_type == DataType::MS1DATA) ? 1 : 2;

	// Detect and classify X-variable conditions
	enum class XPattern { NONE, OFFSET, PEAK_PAIR, CROSS_LEVEL };
	std::vector<const Condition *> x_conds;
	for (const auto &cond : query.where_conditions) {
		if (condition_has_x_variable(cond)) {
			x_conds.push_back(&cond);
		}
	}

	XPattern x_pattern = XPattern::NONE;
	if (!x_conds.empty()) {
		// Check if all X-conditions share the same peak-level field
		bool all_same_field = true;
		ConditionField x_field = x_conds[0]->field;
		for (const auto *c : x_conds) {
			if (c->field != x_field) {
				all_same_field = false;
				break;
			}
		}

		if (all_same_field && is_peak_level_field(x_field)) {
			// Same field: check for non-unit coefficients → peak-pair vs offset
			bool has_nonunit_coeff = false;
			for (const auto *c : x_conds) {
				if (c->values[0].x_coefficient != 1.0) {
					has_nonunit_coeff = true;
					break;
				}
			}
			x_pattern = has_nonunit_coeff ? XPattern::PEAK_PAIR : XPattern::OFFSET;
		} else {
			// Check for cross-level: MS1MZ=X AND MS2PREC=X
			bool has_ms1mz = false, has_ms2prec = false;
			for (const auto *c : x_conds) {
				if (c->field == ConditionField::MS1MZ) {
					has_ms1mz = true;
				}
				if (c->field == ConditionField::MS2PREC) {
					has_ms2prec = true;
				}
			}
			if (has_ms1mz && has_ms2prec && x_conds.size() == 2) {
				x_pattern = XPattern::CROSS_LEVEL;
			} else {
				throw std::runtime_error("MassQL transpiler: unsupported X-variable pattern");
			}
		}
	}

	std::ostringstream sql;

	// Base CTE: reuse mzml_peaks() macro for unnesting + enrichment
	std::string base_cte = "__base";
	bool needs_recursive = (x_pattern == XPattern::OFFSET || x_pattern == XPattern::PEAK_PAIR);
	sql << (needs_recursive ? "WITH RECURSIVE " : "WITH ") << base_cte << " AS (\n"
	    << "  SELECT * FROM mzml_peaks(" << source << ") WHERE ms_level = " << ms_level;

	for (const auto &cond : query.where_conditions) {
		if (is_metadata_field(cond.field)) {
			sql << " AND " << metadata_where(cond);
		}
	}
	sql << "\n),\n";

	// X-variable CTEs
	std::vector<std::pair<std::string, bool>> match_ctes; // name, is_excluded
	if (x_pattern == XPattern::OFFSET) {
		std::string x_qual_name;
		sql << x_offset_ctes(x_conds, base_cte, x_qual_name) << ",\n";
		match_ctes.emplace_back(x_qual_name, false);
	} else if (x_pattern == XPattern::PEAK_PAIR) {
		std::string x_qual_name;
		sql << x_peak_pair_ctes(x_conds, base_cte, x_qual_name) << ",\n";
		match_ctes.emplace_back(x_qual_name, false);
	} else if (x_pattern == XPattern::CROSS_LEVEL) {
		std::string x_qual_name;
		sql << x_cross_level_ctes(x_conds, source, x_qual_name) << ",\n";
		match_ctes.emplace_back(x_qual_name, false);
	}

	// Scan-matching CTEs for non-X peak-level and precursor conditions
	int match_idx = 0;
	for (const auto &cond : query.where_conditions) {
		if (is_metadata_field(cond.field) || condition_has_x_variable(cond)) {
			continue;
		}

		sql << scan_match_cte(cond, base_cte, match_idx) << ",\n";
		match_ctes.emplace_back("__match_" + std::to_string(match_idx), has_qualifier(cond.qualifiers, "EXCLUDED"));
		match_idx++;
	}

	// Build filtered CTE: INTERSECT all match CTEs (AND semantics)
	std::string filtered_cte = "__filtered";
	if (match_ctes.empty()) {
		sql << filtered_cte << " AS (SELECT * FROM " << base_cte << ")";
	} else {
		sql << filtered_cte << " AS (SELECT * FROM " << base_cte << " WHERE ";
		for (size_t i = 0; i < match_ctes.size(); i++) {
			if (i > 0) {
				sql << " AND ";
			}
			sql << "spectrum_index " << (match_ctes[i].second ? "NOT IN" : "IN") << " (SELECT spectrum_index FROM "
			    << match_ctes[i].first << ")";
		}
		sql << ")";
	}

	// Apply FILTER clause (peak-level post-filtering)
	std::string final_cte = filtered_cte;
	if (!query.filter_conditions.empty()) {
		final_cte = "__post_filter";
		sql << ",\n" << final_cte << " AS (SELECT * FROM " << filtered_cte << " WHERE ";
		for (size_t i = 0; i < query.filter_conditions.size(); i++) {
			if (i > 0) {
				sql << " AND ";
			}
			sql << "(" << filter_peak_expr(query.filter_conditions[i]) << ")";
		}
		sql << ")";
	}

	sql << "\n" << aggregation_sql(query.agg_function, final_cte);

	return sql.str();
}

} // namespace miint
