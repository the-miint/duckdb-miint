#pragma once

#include "duckdb/main/extension/extension_loader.hpp"

namespace duckdb {

void RegisterUnifracPcoa(ExtensionLoader &loader);
void RegisterUnifracPermanova(ExtensionLoader &loader);
void RegisterUnifracFaithPD(ExtensionLoader &loader);
void RegisterUnifracDistances(ExtensionLoader &loader);
void RegisterPcoaFromDistances(ExtensionLoader &loader);
void RegisterProgressivePcoaFromDistances(ExtensionLoader &loader);
void RegisterProgressivePcoaFromUnifrac(ExtensionLoader &loader);
void RegisterPermanovaFromDistances(ExtensionLoader &loader);
void RegisterRarefyFeatureTable(ExtensionLoader &loader);

} // namespace duckdb
