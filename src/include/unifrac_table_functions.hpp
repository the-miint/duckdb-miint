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
void RegisterProgressivePcoaFromFeatures(ExtensionLoader &loader);
void RegisterPermanovaFromDistances(ExtensionLoader &loader);
void RegisterRarefyFeatureTable(ExtensionLoader &loader);
// RegisterPickAnchors lives in pick_anchors.hpp: it consumes ordination
// coordinates, not a feature table or a tree, so it carries no UniFrac dependency
// and is registered unconditionally.

} // namespace duckdb
