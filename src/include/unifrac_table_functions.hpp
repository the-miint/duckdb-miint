#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

void RegisterUnifracPcoa(ExtensionLoader &loader);
void RegisterUnifracPermanova(ExtensionLoader &loader);
void RegisterUnifracFaithPD(ExtensionLoader &loader);

} // namespace duckdb
