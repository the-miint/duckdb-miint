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

// Append intensity and massdefect qualifier expressions
static void append_qualifier_exprs(std::string &expr, const Condition &cond) {
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
	if (has_qualifier(cond.qualifiers, "INTENSITYTICPERCENT")) {
		double pct = get_qualifier_value(cond.qualifiers, "INTENSITYTICPERCENT");
		auto op = get_qualifier_op(cond.qualifiers, "INTENSITYTICPERCENT");
		std::string cmp = (op == QualifierOp::GREATER_THAN) ? " > " : (op == QualifierOp::LESS_THAN) ? " < " : " >= ";
		expr += " AND i_tic" + cmp + fmt_double(pct / 100.0);
	}
	if (has_qualifier(cond.qualifiers, "MASSDEFECT")) {
		double md_min = get_qualifier_value(cond.qualifiers, "MASSDEFECT");
		double md_max = 0.0;
		for (const auto &q : cond.qualifiers) {
			if (q.name == "MASSDEFECT") {
				md_max = q.max_value;
				break;
			}
		}
		expr += " AND (mz - FLOOR(mz)) > " + fmt_double(md_min) + " AND (mz - FLOOR(mz)) < " + fmt_double(md_max);
	}
}

// Peak match with optional intensity/massdefect qualifiers
static std::string peak_match_expr(const std::string &col, double target, const Condition &cond) {
	std::string expr = tolerance_expr(col, target, cond.qualifiers);
	append_qualifier_exprs(expr, cond);
	return expr;
}

// ANY wildcard match: no m/z filter, just intensity/massdefect qualifiers
static std::string any_match_expr(const Condition &cond) {
	std::string expr = "intensity > 0";
	append_qualifier_exprs(expr, cond);
	return expr;
}

// ── Y-variable helpers ───────────────────────────────────────────────────────

// Generate SQL expression for Y intensity factor: "y_col * constant" or "y_col * (constant + x_coeff * x_col)"
static std::string y_expr_sql(const Qualifier &q, const std::string &y_col, const std::string &x_col) {
	if (q.y_expr_has_x) {
		return y_col + " * (" + fmt_double(q.y_expr_constant) + " + " + fmt_double(q.y_expr_x_coeff) + " * " + x_col +
		       ")";
	}
	return y_col + " * " + fmt_double(q.y_expr_constant);
}

static bool condition_has_y_variable(const Condition &cond) {
	return has_qualifier(cond.qualifiers, "INTENSITYMATCH") ||
	       has_qualifier(cond.qualifiers, "INTENSITYMATCHREFERENCE");
}

static bool condition_is_y_ref(const Condition &cond) {
	return has_qualifier(cond.qualifiers, "INTENSITYMATCHREFERENCE");
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
	// Enumerate distinct mz values as X candidates
	ss << "__x_candidates AS (\n"
	   << "  SELECT DISTINCT mz AS x_val FROM " << base_cte << " WHERE intensity > 0\n"
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
	// Enumerate distinct mz values as X candidates
	ss << "__x_candidates AS (\n"
	   << "  SELECT DISTINCT mz AS x_val FROM " << base_cte << " WHERE intensity > 0\n"
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
	if (cond.values.empty() && cond.field != ConditionField::POLARITY) {
		throw std::runtime_error("MassQL transpiler: condition has no values");
	}
	std::ostringstream ss;
	std::string cte_name = "__match_" + std::to_string(idx);

	// Check for OTHERSCAN qualifier — expands match to neighboring scans by retention time
	bool has_otherscan = has_qualifier(cond.qualifiers, "OTHERSCAN");
	double rt_left = 0.0, rt_right = 0.0;
	if (has_otherscan) {
		for (const auto &q : cond.qualifiers) {
			if (q.name == "OTHERSCAN") {
				rt_left = q.value;
				rt_right = q.max_value;
				break;
			}
		}
	}

	if (has_otherscan && cond.values.size() > 1) {
		throw std::runtime_error("MassQL transpiler: OTHERSCAN is not supported with OR value lists");
	}

	if (cond.field == ConditionField::MS2PREC) {
		if (has_otherscan) {
			ss << cte_name << " AS (\n"
			   << "  SELECT DISTINCT b2.spectrum_index\n"
			   << "  FROM (SELECT DISTINCT spectrum_index, retention_time FROM " << base_cte << " WHERE "
			   << tolerance_expr("precursor_mz", cond.values[0].constant_value, cond.qualifiers) << ") b1\n"
			   << "  JOIN " << base_cte << " b2 ON b2.retention_time >= b1.retention_time - " << fmt_double(rt_left)
			   << "\n"
			   << "                AND b2.retention_time <= b1.retention_time + " << fmt_double(rt_right) << "\n"
			   << ")";
		} else {
			ss << cte_name << " AS (\n  SELECT DISTINCT spectrum_index FROM " << base_cte << "\n  WHERE "
			   << tolerance_expr("precursor_mz", cond.values[0].constant_value, cond.qualifiers) << "\n)";
		}
	} else if (is_peak_level_field(cond.field)) {
		auto col = peak_column(cond.field);

		if (has_otherscan && cond.values.size() == 1) {
			// OTHERSCAN: find scans with peak match, then expand to neighboring scans by RT window
			// Use subquery to avoid column ambiguity between b1 and b2 aliases
			std::string inner_match;
			if (cond.values[0].has_any_wildcard) {
				inner_match = any_match_expr(cond);
			} else {
				inner_match = peak_match_expr(col, cond.values[0].constant_value, cond);
			}
			ss << cte_name << " AS (\n"
			   << "  SELECT DISTINCT b2.spectrum_index\n"
			   << "  FROM (SELECT DISTINCT spectrum_index, retention_time FROM " << base_cte << " WHERE " << inner_match
			   << ") b1\n"
			   << "  JOIN " << base_cte << " b2 ON b2.retention_time >= b1.retention_time - " << fmt_double(rt_left)
			   << "\n"
			   << "                AND b2.retention_time <= b1.retention_time + " << fmt_double(rt_right) << "\n"
			   << ")";
		} else if (cond.values.size() == 1 && cond.values[0].has_any_wildcard) {
			ss << cte_name << " AS (\n  SELECT DISTINCT spectrum_index FROM " << base_cte << "\n  WHERE "
			   << any_match_expr(cond) << "\n)";
		} else if (cond.values.size() == 1) {
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
	if (cond.values.empty() && cond.field != ConditionField::POLARITY) {
		throw std::runtime_error("MassQL transpiler: metadata condition has no values");
	}
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
		if (cond.string_value != "positive" && cond.string_value != "negative") {
			throw std::runtime_error("MassQL transpiler: invalid polarity value '" + cond.string_value + "'");
		}
		ss << "polarity = '" << cond.string_value << "'";
		break;
	default:
		throw std::runtime_error("MassQL transpiler: not a metadata field");
	}
	return ss.str();
}

// ── Filter expression ────────────────────────────────────────────────────────

static std::string filter_peak_expr(const Condition &cond) {
	if (cond.values.empty()) {
		throw std::runtime_error("MassQL transpiler: filter condition has no values");
	}
	if (cond.field == ConditionField::MS2PREC) {
		return tolerance_expr("precursor_mz", cond.values[0].constant_value, cond.qualifiers);
	}
	if (cond.values.size() == 1 && cond.values[0].has_any_wildcard) {
		return any_match_expr(cond);
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

static std::string aggregation_sql(AggFunction agg, const std::string &source_cte, double scanrangesum_tol = 0.0) {
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
	case AggFunction::SCANRANGESUM: {
		double bin_width = (scanrangesum_tol > 0.0) ? scanrangesum_tol : 0.1;
		ss << "SELECT spectrum_index, FLOOR(mz / " << fmt_double(bin_width) << ") * " << fmt_double(bin_width)
		   << " AS mz_bin, SUM(intensity) AS total_intensity FROM " << source_cte << " GROUP BY spectrum_index, mz_bin";
		break;
	}
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
	sql << "WITH " << base_cte << " AS (\n"
	    << "  SELECT * FROM mzml_peaks(" << source << ") WHERE ms_level = " << ms_level;

	for (const auto &cond : query.where_conditions) {
		if (is_metadata_field(cond.field)) {
			sql << " AND " << metadata_where(cond);
		}
	}
	sql << "\n),\n";

	// Check if X-conditions also have Y-variable qualifiers
	bool x_has_y = false;
	for (const auto *c : x_conds) {
		if (condition_has_y_variable(*c)) {
			x_has_y = true;
			break;
		}
	}

	// X-variable CTEs
	std::vector<std::pair<std::string, bool>> match_ctes; // name, is_excluded
	if (x_has_y && x_pattern == XPattern::OFFSET) {
		// Combined X+Y pattern: X-enumeration with intensity matching
		// Find ref and match conditions among X-conditions
		const Condition *xy_ref = nullptr;
		std::vector<const Condition *> xy_matches;
		for (const auto *c : x_conds) {
			if (condition_is_y_ref(*c)) {
				xy_ref = c;
			} else {
				xy_matches.push_back(c);
			}
		}
		if (!xy_ref || xy_matches.empty()) {
			throw std::runtime_error("MassQL transpiler: X+Y pattern requires one INTENSITYMATCHREFERENCE condition");
		}

		// Extract offsets, normalize so ref offset is 0
		double ref_offset = xy_ref->values[0].constant_value;
		std::string tol = x_offset_tolerance(x_conds);

		// Enumerate distinct mz values as X candidates
		sql << "__x_candidates AS (\n"
		    << "  SELECT DISTINCT mz AS x_val FROM " << base_cte << " WHERE intensity > 0\n"
		    << "),\n";

		// Y-ref CTE: for each X candidate, find intensity at ref offset
		sql << "__y_ref AS (\n"
		    << "  SELECT xc.x_val, p.spectrum_index, SUM(p.intensity) AS y_val\n"
		    << "  FROM __x_candidates xc\n"
		    << "  JOIN " << base_cte << " p ON p.mz > xc.x_val + " << fmt_double(ref_offset) << " - " << tol << "\n"
		    << "              AND p.mz < xc.x_val + " << fmt_double(ref_offset) << " + " << tol << "\n"
		    << "              AND p.intensity > 0\n"
		    << "  GROUP BY xc.x_val, p.spectrum_index\n"
		    << "),\n";

		// Y-match CTEs: for each match condition, check intensity at offset
		std::vector<std::string> xy_match_names;
		int xy_idx = 0;
		for (const auto *mcond : xy_matches) {
			double match_offset = mcond->values[0].constant_value;
			double y_match_pct = 20.0;
			const Qualifier *y_match_qual = nullptr;
			for (const auto &q : mcond->qualifiers) {
				if (q.name == "INTENSITYMATCH") {
					y_match_qual = &q;
				}
				if (q.name == "INTENSITYMATCHPERCENT") {
					y_match_pct = q.value;
				}
			}

			// Build Y-expression SQL (may reference r.x_val for X-dependent expressions)
			std::string y_expected;
			if (y_match_qual) {
				y_expected = y_expr_sql(*y_match_qual, "r.y_val", "r.x_val");
			} else {
				y_expected = "r.y_val";
			}

			std::string cte_name = "__xy_match_" + std::to_string(xy_idx);
			sql << cte_name << " AS (\n"
			    << "  SELECT DISTINCT r.x_val, r.spectrum_index\n"
			    << "  FROM __y_ref r\n"
			    << "  JOIN " << base_cte << " p ON p.spectrum_index = r.spectrum_index\n"
			    << "    AND p.mz > r.x_val + " << fmt_double(match_offset) << " - " << tol << "\n"
			    << "    AND p.mz < r.x_val + " << fmt_double(match_offset) << " + " << tol << "\n"
			    << "    AND p.intensity > " << y_expected << " * (1.0 - " << fmt_double(y_match_pct / 100.0) << ")\n"
			    << "    AND p.intensity < " << y_expected << " * (1.0 + " << fmt_double(y_match_pct / 100.0) << ")\n"
			    << "),\n";
			xy_match_names.push_back(cte_name);
			xy_idx++;
		}

		// X+Y qualifying: require all match conditions agree on the same (x_val, spectrum_index)
		std::string xy_qual_name = "__xy_qualifying";
		sql << xy_qual_name << " AS (\n  SELECT DISTINCT spectrum_index FROM (\n";
		for (size_t i = 0; i < xy_match_names.size(); i++) {
			if (i > 0) {
				sql << "    INTERSECT\n";
			}
			sql << "    SELECT x_val, spectrum_index FROM " << xy_match_names[i] << "\n";
		}
		sql << "  )\n),\n";
		match_ctes.emplace_back(xy_qual_name, false);
	} else if (x_pattern == XPattern::OFFSET) {
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

	// Y-variable CTEs (only for non-X Y-conditions; X+Y conditions handled above)
	const Condition *y_ref_cond = nullptr;
	std::vector<const Condition *> y_match_conds;
	for (const auto &cond : query.where_conditions) {
		if (condition_has_y_variable(cond) && !condition_has_x_variable(cond)) {
			if (condition_is_y_ref(cond)) {
				y_ref_cond = &cond;
			} else {
				y_match_conds.push_back(&cond);
			}
		}
	}
	if (y_ref_cond && !y_match_conds.empty()) {
		// __y_ref CTE: match reference m/z, get intensity as y_val per spectrum
		auto col = peak_column(y_ref_cond->field);
		double ref_target = y_ref_cond->values[0].constant_value;
		sql << "__y_ref AS (\n"
		    << "  SELECT spectrum_index, SUM(intensity) AS y_val\n"
		    << "  FROM " << base_cte << "\n"
		    << "  WHERE " << tolerance_expr(col, ref_target, y_ref_cond->qualifiers) << "\n"
		    << "  GROUP BY spectrum_index\n"
		    << "),\n";

		// __y_match_N CTEs: for each match condition
		int y_match_idx = 0;
		for (const auto *mcond : y_match_conds) {
			auto mcol = peak_column(mcond->field);
			double mtarget = mcond->values[0].constant_value;

			// Get Y-expression parameters
			double y_match_pct = 20.0; // default tolerance
			const Qualifier *y_match_qual = nullptr;
			for (const auto &q : mcond->qualifiers) {
				if (q.name == "INTENSITYMATCH") {
					y_match_qual = &q;
				}
				if (q.name == "INTENSITYMATCHPERCENT") {
					y_match_pct = q.value;
				}
			}

			std::string y_expected;
			if (y_match_qual) {
				y_expected = y_expr_sql(*y_match_qual, "r.y_val", "0"); // no X in Y-only
			} else {
				y_expected = "r.y_val";
			}

			std::string cte_name = "__y_match_" + std::to_string(y_match_idx);
			sql << cte_name << " AS (\n"
			    << "  SELECT DISTINCT b.spectrum_index\n"
			    << "  FROM " << base_cte << " b\n"
			    << "  JOIN __y_ref r ON b.spectrum_index = r.spectrum_index\n"
			    << "  WHERE " << tolerance_expr("b." + mcol, mtarget, mcond->qualifiers) << "\n"
			    << "    AND b.intensity > " << y_expected << " * (1.0 - " << fmt_double(y_match_pct / 100.0) << ")\n"
			    << "    AND b.intensity < " << y_expected << " * (1.0 + " << fmt_double(y_match_pct / 100.0) << ")\n"
			    << "),\n";

			match_ctes.emplace_back(cte_name, false);
			y_match_idx++;
		}
	}

	// Scan-matching CTEs for non-X, non-Y peak-level and precursor conditions
	int match_idx = 0;
	for (const auto &cond : query.where_conditions) {
		if (is_metadata_field(cond.field) || condition_has_x_variable(cond) || condition_has_y_variable(cond)) {
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

	sql << "\n" << aggregation_sql(query.agg_function, final_cte, query.scanrangesum_tolerance);

	return sql.str();
}

} // namespace miint
