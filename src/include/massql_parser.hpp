#pragma once

#include <string>
#include <vector>

namespace miint {

enum class DataType { MS1DATA, MS2DATA };

enum class AggFunction { NONE, SCANINFO, SCANSUM, SCANNUM, SCANMZ, SCANMAXINT, SCANRANGESUM };

enum class ConditionField {
	// Peak-level fields
	MS1MZ,
	MS2PROD,
	MS2PREC,
	MS2NL,
	// Metadata fields
	RTMIN,
	RTMAX,
	SCANMIN,
	SCANMAX,
	CHARGE,
	POLARITY
};

struct ConditionValue {
	double constant_value = 0.0;
	bool has_x_variable = false;
	bool has_any_wildcard = false;
	double x_coefficient = 1.0;
};

enum class QualifierOp { EQUALS, GREATER_THAN, LESS_THAN, NONE };

struct Qualifier {
	std::string name;
	double value = 0.0;
	double max_value = 0.0; // for range(min=a,max=b)
	QualifierOp op = QualifierOp::NONE;
	// Y-expression fields (for INTENSITYMATCH qualifier)
	double y_expr_constant = 1.0; // Y * constant (or Y * (constant + x_coeff*X))
	double y_expr_x_coeff = 0.0;  // coefficient of X in Y expression
	bool y_expr_has_x = false;
};

struct Condition {
	ConditionField field;
	std::vector<ConditionValue> values;
	std::string string_value; // for POLARITY
	std::vector<Qualifier> qualifiers;
};

struct MassQLQuery {
	DataType data_type;
	AggFunction agg_function = AggFunction::NONE;
	double scanrangesum_tolerance = 0.0; // 0 means use default (0.1 Da)
	std::vector<Condition> where_conditions;
	std::vector<Condition> filter_conditions;
};

class MassQLParser {
public:
	static MassQLQuery parse(const std::string &query);
};

} // namespace miint
