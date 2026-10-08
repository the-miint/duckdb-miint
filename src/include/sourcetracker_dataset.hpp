#pragma once
//
// The relation-to-array translation for sourcetracker(): long-form feature-table
// cells plus (sample_id, source_sink, env) metadata rows in, the index-based
// arrays st3_table_from_arrow takes out. DuckDB-free for the same reason
// mmvec_relation.hpp is: st3 addresses features and samples by INDEX, so a cell
// filed under the wrong sample, or a role attached to the wrong id, does not
// fail -- the sampler runs and returns well-formed proportions for the wrong
// sink. That has to be asserted against independently known values, in a unit
// test that does not need libduckdb.
//
// Every precondition st3 does not check itself is checked here, before any
// Arrow marshaling, and each error names the offending id. Errors are
// std::invalid_argument; the SQL layer wraps them into the user-facing
// exception with the function name prefixed.

#include <cstdint>
#include <optional>
#include <string>
#include <vector>

#include "feature_table_row.hpp"
#include "unifrac_metadata.hpp"

namespace miint::sourcetracker {

// Exactly the shape st3_table_from_arrow consumes: a COO count table over
// feature and sample indices, feature ids in feature order, and per-sample
// role/env in sample order. Both dictionaries are lexicographic, so the sink
// order of the SQL output and the key order of the assignments MAP are fixed
// by the ids, not by the scan.
struct Dataset {
	std::vector<std::string> feature_ids;
	std::vector<std::string> sample_ids;
	std::vector<bool> is_source;                  // aligned to sample_ids
	std::vector<std::optional<std::string>> envs; // aligned to sample_ids; nullopt for sinks
	std::vector<int32_t> rows;                    // feature index per cell
	std::vector<int32_t> cols;                    // sample index per cell
	std::vector<double> vals;                     // count per cell
	std::vector<double> sample_totals;            // aligned to sample_ids
};

// Build a Dataset from what ReadFeatureTable and ReadWideMetadata return.
//
// `metadata` is the unpivoted long form: one row per (sample_id, variable) with
// the variables `source_sink` and `env` (names matched case-insensitively, as
// are the role values `source`/`sink`). A NULL env arrives as "". A sink's env
// is ignored.
//
// Throws std::invalid_argument, naming the id, when:
//   - a feature-table sample has no metadata row, or a metadata sample has no
//     cell in the feature table (the reader has already dropped zero cells, so
//     an all-zero sample and an absent one are the same event);
//   - a (sample, feature) cell is duplicated;
//   - a count is negative or not finite;
//   - a sample is listed twice in the metadata, has no source_sink row, or has
//     a source_sink value other than source/sink;
//   - a source has an empty env;
//   - there is no source at all.
Dataset IngestDataset(const std::vector<miint::unifrac::CooRow> &table,
                      const std::vector<miint::unifrac::MetadataRow> &metadata);

// The per-sample rarefaction precondition SourceTracker2 states up front: in
// sink mode every sink must hold at least `sink_depth` sequences; in
// leave-one-out every source sample must hold at least `source_depth` (the sink
// depth is unused there). A depth of 0 disables the check for that role. The
// source depth in sink mode is NOT checked here: sources are collapsed by
// environment first, and st3 reports a shallow collapsed environment itself.
//
// Throws std::invalid_argument in SourceTracker2's words, naming the shallowest
// sample and its total.
void CheckDepths(const Dataset &dataset, int32_t source_depth, int32_t sink_depth, bool loo);

} // namespace miint::sourcetracker
