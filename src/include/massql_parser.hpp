#pragma once

#include <string>
#include <vector>

namespace miint {

enum class DataType { MS1DATA, MS2DATA };

enum class AggFunction { NONE, SCANINFO, SCANSUM, SCANNUM, SCANMZ, SCANMAXINT };

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
	double x_coefficient = 1.0;
};

enum class QualifierOp { EQUALS, GREATER_THAN, NONE };

struct Qualifier {
	std::string name;
	double value = 0.0;
	double max_value = 0.0; // for range(min=a,max=b)
	QualifierOp op = QualifierOp::NONE;
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
	std::vector<Condition> where_conditions;
	std::vector<Condition> filter_conditions;
};

class MassQLParser {
public:
	static MassQLQuery parse(const std::string &query);
};

} // namespace miint
