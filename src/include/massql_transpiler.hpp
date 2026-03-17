#pragma once

#include "massql_parser.hpp"

#include <string>

namespace miint {

struct MaterializationPlan {
	bool needs_ms1; // true for cross-level queries: MS1MZ=X AND MS2PREC=X
};

class MassQLTranspiler {
public:
	static std::string to_sql(const MassQLQuery &query, const std::string &source);
	static std::string to_sql_materialized(const MassQLQuery &query, const std::string &source,
	                                       const std::string &base_table, const std::string &ms1_table = "");
	static std::string materialize_base_sql(const MassQLQuery &query, const std::string &source,
	                                        const std::string &table_name);
	static std::string materialize_ms1_sql(const std::string &source, const std::string &table_name);
	static MaterializationPlan get_materialization_plan(const MassQLQuery &query);
	static bool is_file_path(const std::string &source);
};

} // namespace miint
