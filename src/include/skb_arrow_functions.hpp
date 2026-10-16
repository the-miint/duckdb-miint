#pragma once

#include "duckdb.hpp"

namespace duckdb {

/// Registers the table functions that call skb-arrow (skb_arrow/client.hpp):
/// - `skb_mantel(x, y, method, permutations, alternative, seed, threads)`:
///   the Mantel test between two condensed `(sample_a, sample_b, distance)`
///   relations. Returns one row: `statistic DOUBLE, pvalue DOUBLE, n BIGINT`.
/// - `skb_ancombc(table, metadata, formula, grouping, posthoc, pseudocount,
///   max_iter, tol, alpha, p_adjust, bootstraps, seed, threads)`: ANCOM-BC
///   over a long `(sample_id, feature_id, value)` feature table and a wide
///   metadata relation. Returns `feature_id, test, term, lfc, se, w, pvalue,
///   qvalue, signif`, in skb-arrow's order.
void RegisterSkbArrowFunctions(ExtensionLoader &loader);

} // namespace duckdb
