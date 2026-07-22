#pragma once

#include "duckdb/main/client_context.hpp"
#include <map>
#include <string>
#include <unordered_map>

namespace duckdb {

// Validate that a long-form traits table/view has the required columns
// (name, trait, value). Shared by phylo_independent_contrasts and the
// phylo_ancestral_* functions (continuous and discrete traits share this
// long-form schema; only the `value` cast differs between readers).
void ValidateTraitsTableSchema(ClientContext &context, const std::string &table_name);

// Read long-form continuous traits into trait -> (tip name -> value). Uses a
// separate Connection (see docs/internals/reading-tables-views.md). `name` is
// cast to VARCHAR so UUID/BIGINT/VARCHAR tip keys all match tip labels by
// canonical text. A std::map keeps trait order stable (deterministic output).
// NULLs and duplicate (name, trait) pairs are hard errors.
std::map<std::string, std::unordered_map<std::string, double>> ReadContinuousTraits(ClientContext &context,
                                                                                    const std::string &table_name);

// Read long-form discrete traits into trait -> (tip name -> state string). Like
// ReadContinuousTraits, but `value` is cast to VARCHAR so integer- or text-coded
// states both work. The per-trait state alphabet (sorted distinct states mapped to
// 0..k-1) is built by the caller — it may instead come from a supplied cost matrix,
// which the reader has no knowledge of. A std::map keeps trait order stable. NULLs
// and duplicate (name, trait) pairs are hard errors.
std::map<std::string, std::unordered_map<std::string, std::string>> ReadDiscreteTraits(ClientContext &context,
                                                                                       const std::string &table_name);

} // namespace duckdb
