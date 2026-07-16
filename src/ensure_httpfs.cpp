#include "ensure_httpfs.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace miint {

void EnsureHttpfsLoaded(duckdb::ClientContext &context) {
	duckdb::DatabaseInstance &instance = *context.db;
	if (instance.ExtensionIsLoaded("httpfs")) {
		return;
	}
	duckdb::Connection con(instance);
	auto result = con.Query("LOAD httpfs;");
	if (result->HasError()) {
		// httpfs may not be installed yet; try to install then load.
		con.Query("INSTALL httpfs;");
		con.Query("LOAD httpfs;");
	}
	if (!instance.ExtensionIsLoaded("httpfs")) {
		throw duckdb::IOException("remote NCBI access requires the httpfs extension, which could not be loaded "
		                          "automatically; run INSTALL httpfs; LOAD httpfs; first");
	}
}

} // namespace miint
