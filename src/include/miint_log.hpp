#pragma once

#include "duckdb/common/string_util.hpp"
#include "duckdb/logging/log_type.hpp"
#include "duckdb/logging/logging.hpp"

#include <string>

namespace duckdb {
class ClientContext;
class DatabaseInstance;
} // namespace duckdb

namespace miint {

class MiintWarningLogType : public duckdb::LogType {
public:
	static constexpr const char *NAME = "MiintWarning";
	static constexpr duckdb::LogLevel LEVEL = duckdb::LogLevel::LOG_WARNING;

	MiintWarningLogType();
};

// Register the miint warning log type with the database's LogManager. Call once
// from the extension loader; subsequent loads on the same instance are no-ops.
void RegisterMiintLogType(duckdb::DatabaseInstance &db);

// Emit a user-facing warning. Writes to BOTH:
//   1. The global MutableLogger (queryable via duckdb_logs() / miint_warnings())
//   2. Printer::Print (stderr, matches today's interactive visibility)
//
// The MutableLogger writes unconditionally, so miint_warnings() returns rows
// regardless of the user's enable_logging / logging_storage settings (as long
// as the default in-memory storage is kept).
void EmitWarning(duckdb::ClientContext &ctx, const std::string &msg);

template <typename... Args>
void EmitWarning(duckdb::ClientContext &ctx, const char *fmt, Args &&...args) {
	EmitWarning(ctx, duckdb::StringUtil::Format(fmt, std::forward<Args>(args)...));
}

} // namespace miint
