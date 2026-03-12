#include <massql_parser.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using miint::AggFunction;
using miint::ConditionField;
using miint::DataType;
using miint::MassQLParser;
using miint::QualifierOp;

// ===== Cycle 1: Parser skeleton — QUERY keyword + data type =====

TEST_CASE("empty string throws", "[massql][parser]") {
	REQUIRE_THROWS(MassQLParser::parse(""));
}

TEST_CASE("missing QUERY keyword throws", "[massql][parser]") {
	REQUIRE_THROWS(MassQLParser::parse("scaninfo(MS2DATA)"));
}

TEST_CASE("bare QUERY MS2DATA", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA");
	REQUIRE(q.data_type == DataType::MS2DATA);
	REQUIRE(q.agg_function == AggFunction::NONE);
	REQUIRE(q.where_conditions.empty());
}

TEST_CASE("QUERY MS1DATA", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS1DATA");
	REQUIRE(q.data_type == DataType::MS1DATA);
}

TEST_CASE("case insensitive", "[massql][parser]") {
	auto q = MassQLParser::parse("query ms2data");
	REQUIRE(q.data_type == DataType::MS2DATA);
}

// ===== Cycle 2: Aggregation functions =====

TEST_CASE("scaninfo(MS2DATA)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scaninfo(MS2DATA)");
	REQUIRE(q.agg_function == AggFunction::SCANINFO);
	REQUIRE(q.data_type == DataType::MS2DATA);
}

TEST_CASE("scansum(MS1DATA)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scansum(MS1DATA)");
	REQUIRE(q.agg_function == AggFunction::SCANSUM);
	REQUIRE(q.data_type == DataType::MS1DATA);
}

TEST_CASE("scannum(MS2DATA)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA)");
	REQUIRE(q.agg_function == AggFunction::SCANNUM);
}

TEST_CASE("scanmz(MS2DATA)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scanmz(MS2DATA)");
	REQUIRE(q.agg_function == AggFunction::SCANMZ);
}

TEST_CASE("scanmaxint(MS2DATA)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scanmaxint(MS2DATA)");
	REQUIRE(q.agg_function == AggFunction::SCANMAXINT);
}

TEST_CASE("unknown aggregation function throws", "[massql][parser]") {
	REQUIRE_THROWS(MassQLParser::parse("QUERY bogus(MS2DATA)"));
}

// ===== Cycle 3: Simple WHERE conditions =====

TEST_CASE("WHERE MS2PROD=226.18", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scaninfo(MS2DATA) WHERE MS2PROD=226.18");
	REQUIRE(q.where_conditions.size() == 1);
	REQUIRE(q.where_conditions[0].field == ConditionField::MS2PROD);
	REQUIRE(q.where_conditions[0].values.size() == 1);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(226.18, 0.001));
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
}

TEST_CASE("WHERE MS2PREC=200", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PREC=200");
	REQUIRE(q.where_conditions.size() == 1);
	REQUIRE(q.where_conditions[0].field == ConditionField::MS2PREC);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(200.0, 0.001));
}

TEST_CASE("WHERE MS1MZ=100", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS1DATA) WHERE MS1MZ=100");
	REQUIRE(q.where_conditions[0].field == ConditionField::MS1MZ);
}

TEST_CASE("WHERE MS2NL=80", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2NL=80");
	REQUIRE(q.where_conditions[0].field == ConditionField::MS2NL);
}

TEST_CASE("MS2MZ is alias for MS2PROD", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2MZ=150");
	REQUIRE(q.where_conditions[0].field == ConditionField::MS2PROD);
}

// ===== Cycle 4: Metadata conditions =====

TEST_CASE("WHERE RTMIN=1.7", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE RTMIN=1.7");
	REQUIRE(q.where_conditions[0].field == ConditionField::RTMIN);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(1.7, 0.001));
}

TEST_CASE("WHERE RTMAX=2.0", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE RTMAX=2.0");
	REQUIRE(q.where_conditions[0].field == ConditionField::RTMAX);
}

TEST_CASE("WHERE CHARGE=2", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE CHARGE=2");
	REQUIRE(q.where_conditions[0].field == ConditionField::CHARGE);
}

TEST_CASE("WHERE POLARITY=Positive", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE POLARITY=Positive");
	REQUIRE(q.where_conditions[0].field == ConditionField::POLARITY);
	REQUIRE(q.where_conditions[0].string_value == "positive");
}

// ===== Cycle 5: Qualifiers =====

TEST_CASE("TOLERANCEMZ qualifier", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=200:TOLERANCEMZ=0.1");
	REQUIRE(q.where_conditions[0].qualifiers.size() == 1);
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "TOLERANCEMZ");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(0.1, 0.0001));
}

TEST_CASE("TOLERANCEPPM qualifier", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=200:TOLERANCEPPM=10");
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "TOLERANCEPPM");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(10.0, 0.001));
}

TEST_CASE("INTENSITYPERCENT qualifier with >", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT>50");
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "INTENSITYPERCENT");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(50.0, 0.001));
}

TEST_CASE("INTENSITYVALUE qualifier with >", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYVALUE>3000");
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "INTENSITYVALUE");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(3000.0, 0.001));
}

TEST_CASE("EXCLUDED qualifier", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=150:EXCLUDED");
	REQUIRE(q.where_conditions[0].qualifiers.size() == 1);
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "EXCLUDED");
}

TEST_CASE("multiple qualifiers chained", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=200:TOLERANCEMZ=0.1:INTENSITYPERCENT>50");
	REQUIRE(q.where_conditions[0].qualifiers.size() == 2);
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "TOLERANCEMZ");
	REQUIRE(q.where_conditions[0].qualifiers[1].name == "INTENSITYPERCENT");
}

// ===== Cycle 6: AND conjunction + FILTER clause =====

TEST_CASE("AND conjunction", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=150 AND MS2PREC=200");
	REQUIRE(q.where_conditions.size() == 2);
	REQUIRE(q.where_conditions[0].field == ConditionField::MS2PROD);
	REQUIRE(q.where_conditions[1].field == ConditionField::MS2PREC);
}

TEST_CASE("FILTER clause", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=220 FILTER MS2PROD=220");
	REQUIRE(q.where_conditions.size() == 1);
	REQUIRE(q.filter_conditions.size() == 1);
	REQUIRE(q.filter_conditions[0].field == ConditionField::MS2PROD);
}

// ===== Cycle 7: X-variable expressions =====

TEST_CASE("MS2PROD=X", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=X");
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].x_coefficient, Catch::Matchers::WithinAbs(1.0, 0.001));
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(0.0, 0.001));
}

TEST_CASE("MS2PROD=X+164.9", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=X+164.9");
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].x_coefficient, Catch::Matchers::WithinAbs(1.0, 0.001));
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(164.9, 0.001));
}

TEST_CASE("MS2PROD=X-50.5", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=X-50.5");
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(-50.5, 0.001));
}

TEST_CASE("X offset pair: MS2PROD=X AND MS2PROD=X+100", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=X AND MS2PROD=X+100");
	REQUIRE(q.where_conditions.size() == 2);
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(0.0, 0.001));
	REQUIRE(q.where_conditions[1].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[1].values[0].constant_value, Catch::Matchers::WithinAbs(100.0, 0.001));
}

// ===== Cycle 8: OR expressions + CARDINALITY =====

TEST_CASE("OR expression: MS2PROD=(150 OR 250)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250)");
	REQUIRE(q.where_conditions[0].values.size() == 2);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(150.0, 0.001));
	REQUIRE_THAT(q.where_conditions[0].values[1].constant_value, Catch::Matchers::WithinAbs(250.0, 0.001));
}

TEST_CASE("OR expression with three values", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=(100 OR 200 OR 300)");
	REQUIRE(q.where_conditions[0].values.size() == 3);
}

TEST_CASE("CARDINALITY qualifier", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250):CARDINALITY=range(min=2,max=5)");
	REQUIRE(q.where_conditions[0].qualifiers.size() == 1);
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "CARDINALITY");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(2.0, 0.001));
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].max_value, Catch::Matchers::WithinAbs(5.0, 0.001));
}

// ===== Cycle 9: Arithmetic constant folding =====

TEST_CASE("arithmetic in value: 157.0857+10", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=157.0857+10");
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(167.0857, 0.001));
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
}

// ===== Fix validations =====

TEST_CASE("lone dot is not a valid number", "[massql][parser]") {
	REQUIRE_THROWS(MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=."));
}

TEST_CASE("INTENSITYPERCENT= stores EQUALS op", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT=50");
	REQUIRE(q.where_conditions[0].qualifiers[0].op == QualifierOp::EQUALS);
}

TEST_CASE("INTENSITYPERCENT> stores GREATER_THAN op", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS1DATA) WHERE MS1MZ=200:INTENSITYPERCENT>50");
	REQUIRE(q.where_conditions[0].qualifiers[0].op == QualifierOp::GREATER_THAN);
}

TEST_CASE("EXCLUDED qualifier has NONE op", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=150:EXCLUDED");
	REQUIRE(q.where_conditions[0].qualifiers[0].op == QualifierOp::NONE);
}

// ===== Cycle 22: formula() and coefficient expressions =====

TEST_CASE("formula(Fe) as constant", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=formula(Fe)");
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(55.9349375, 1e-4));
}

TEST_CASE("formula(H2O) as constant", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=formula(H2O)");
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
	// H2O = 2*1.00782503207 + 15.99491461956 = 18.01056468
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(18.01056, 0.001));
}

TEST_CASE("X-formula(Fe)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=X-formula(Fe)");
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].x_coefficient, Catch::Matchers::WithinAbs(1.0, 1e-9));
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(-55.9349375, 1e-4));
}

TEST_CASE("2*(X-formula(Fe))", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=2*(X-formula(Fe))");
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].x_coefficient, Catch::Matchers::WithinAbs(2.0, 1e-9));
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(-2.0 * 55.9349375, 1e-4));
}

TEST_CASE("iron pattern full query parse", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scaninfo(MS2DATA) WHERE MS2PROD=X AND MS2PROD=2*(X-formula(Fe))");
	REQUIRE(q.where_conditions.size() == 2);
	REQUIRE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].x_coefficient, Catch::Matchers::WithinAbs(1.0, 1e-9));
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(0.0, 1e-9));
	REQUIRE(q.where_conditions[1].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[1].values[0].x_coefficient, Catch::Matchers::WithinAbs(2.0, 1e-9));
}

TEST_CASE("constant folding with formula: 100+formula(Fe)", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=100+formula(Fe)");
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(100.0 + 55.9349375, 1e-4));
}

// ===== Cycle 2 (Phase 7): MATCHCOUNT alias, < operator, division =====

TEST_CASE("MATCHCOUNT is alias for CARDINALITY", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=(150 OR 250):MATCHCOUNT=range(min=1,max=2)");
	REQUIRE(q.where_conditions[0].qualifiers.size() == 1);
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "CARDINALITY");
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(1.0, 0.001));
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].max_value, Catch::Matchers::WithinAbs(2.0, 0.001));
}

TEST_CASE("LESS_THAN operator: INTENSITYPERCENT<50", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY scannum(MS2DATA) WHERE MS2PROD=220:INTENSITYPERCENT<50");
	REQUIRE(q.where_conditions[0].qualifiers[0].name == "INTENSITYPERCENT");
	REQUIRE(q.where_conditions[0].qualifiers[0].op == QualifierOp::LESS_THAN);
	REQUIRE_THAT(q.where_conditions[0].qualifiers[0].value, Catch::Matchers::WithinAbs(50.0, 0.001));
}

TEST_CASE("division constant folding: MS2PROD=400/2", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=400/2");
	REQUIRE_FALSE(q.where_conditions[0].values[0].has_x_variable);
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(200.0, 0.001));
}

TEST_CASE("division with addition: MS2PROD=400/2+10", "[massql][parser]") {
	auto q = MassQLParser::parse("QUERY MS2DATA WHERE MS2PROD=400/2+10");
	REQUIRE_THAT(q.where_conditions[0].values[0].constant_value, Catch::Matchers::WithinAbs(210.0, 0.001));
}
