// SPDX-License-Identifier: MIT
//
// `DELETE FROM ena.X WHERE col=val` dispatch into the ENA Webin V2 targeted
// CANCEL lifecycle path.
//
// `DELETE` is a thin SQL alias over `ena_cancel`: the WHERE predicate must be
// a single column-equality on `alias` or the per-table accession column,
// e.g. `samea_accession`, `prjeb_accession`, `err_accession`, …
// Anything else is rejected at plan time with an actionable BinderException.
//
// The HTTP POST and submission_log write happen in the produced
// PhysicalOperator's source pipeline — never at plan time. EXPLAIN of a DELETE
// never causes a CANCEL to be issued.

#pragma once

namespace duckdb {
class ClientContext;
class ENACatalog;
class LogicalDelete;
class PhysicalOperator;
class PhysicalPlanGenerator;
} // namespace duckdb

namespace miint {

// Build the physical operator for `DELETE FROM ena.X WHERE col=val`.
// Throws BinderException on unsupported tables, missing/compound/non-equality
// WHERE clauses, RETURNING, and predicates on columns other than alias or
// the table's accession column.
duckdb::PhysicalOperator &PlanENALifecycleDelete(duckdb::ClientContext &context, duckdb::PhysicalPlanGenerator &planner,
                                                 duckdb::ENACatalog &catalog, duckdb::LogicalDelete &op);

} // namespace miint
