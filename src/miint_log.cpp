#include "miint_log.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/printer.hpp"
#include "duckdb/logging/log_manager.hpp"
#include "duckdb/logging/logger.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/database.hpp"

namespace miint {

MiintWarningLogType::MiintWarningLogType() : duckdb::LogType(NAME, LEVEL) {
}

void RegisterMiintLogType(duckdb::DatabaseInstance &db) {
	auto &log_manager = db.GetLogManager();
	// LookupLogType + RegisterLogType is check-then-act. DuckDB's extension
	// loader serializes Load() calls per DatabaseInstance, so concurrent
	// registration is not a concern in practice; this guard just makes a
	// second LOAD of the extension on the same instance a no-op.
	if (log_manager.LookupLogType(MiintWarningLogType::NAME)) {
		return;
	}
	log_manager.RegisterLogType(duckdb::make_uniq<MiintWarningLogType>());
}

void EmitWarning(duckdb::ClientContext &ctx, const std::string &msg) {
	// Why GlobalLogger instead of Logger::Get(ctx)?
	//   Logger::Get(ctx) returns ClientContext::logger, which is a NopLogger
	//   when enable_logging=false (DuckDB's default). Routing through it would
	//   make miint_warnings() silently empty for any user who hasn't run
	//   `SET enable_logging=true`. The plan requires that miint_warnings()
	//   works out of the box, so we write through the database-scoped
	//   GlobalLogger instead — LogManager::Initialize() always constructs it
	//   as a MutableLogger (mutable_settings=true in CreateLogger), whose
	//   WriteLogInternal skips the ShouldLog gate and writes to storage
	//   unconditionally.
	//
	// This relies on a DuckDB internal detail (GlobalLogger is always a
	// MutableLogger). If a future DuckDB refactor returns a NopLogger here,
	// miint_warnings() would silently return zero rows; the dedicated
	// emit-assertion SQL test exists to catch that regression.
	//
	// WarningsAsErrors interaction: LogManager::WriteLogEntry throws
	// InvalidInputException when the user has set warnings_as_errors=true.
	// That would convert every skip-warning into a hard query error that
	// aborts retries and surfaces partial-data skips as fatal — a behavior
	// change from the pre-logging Printer::Print path. We catch it here and
	// still emit to stderr so the user sees the warning in both modes; they
	// keep getting the "continue after skipping" behavior the calling code
	// was designed around.
	try {
		duckdb::LogManager::Get(ctx).GlobalLogger().WriteLog(MiintWarningLogType::NAME, MiintWarningLogType::LEVEL,
		                                                     msg);
	} catch (...) {
		// Fall through to stderr only. miint_warnings() will not see this row
		// under warnings_as_errors=true; that's acceptable because the
		// interactive stderr channel still delivers the message.
	}
	duckdb::Printer::Print(msg);
}

} // namespace miint
