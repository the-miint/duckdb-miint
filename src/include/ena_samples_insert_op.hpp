// SPDX-License-Identifier: MIT
//
// PhysicalOperator subclass that implements `INSERT INTO ena.samples ... [RETURNING ...]`
// by buffering rows in Sink, POSTing the resulting envelope in Finalize, and
// scanning a per-statement ColumnDataCollection in GetDataInternal (for the
// RETURNING projection).
//
// Mirrors ENAProjectsInsert; the duplication is extracted into a shared base
// in the Phase 5 refactor step.

#pragma once

#include "duckdb/common/index_vector.hpp"
#include "duckdb/execution/physical_operator.hpp"

namespace duckdb {

class ENATableEntry;

class ENASamplesInsert : public PhysicalOperator {
public:
	static constexpr const PhysicalOperatorType TYPE = PhysicalOperatorType::EXTENSION;

public:
	ENASamplesInsert(PhysicalPlan &physical_plan, LogicalOperator &op, ENATableEntry &table,
	                 physical_index_vector_t<idx_t> column_index_map, bool return_chunk);

	ENATableEntry &table;
	physical_index_vector_t<idx_t> column_index_map;
	bool return_chunk;

public:
	// Source
	unique_ptr<GlobalSourceState> GetGlobalSourceState(ClientContext &context) const override;
	SourceResultType GetDataInternal(ExecutionContext &context, DataChunk &chunk,
	                                 OperatorSourceInput &input) const override;

	bool IsSource() const override {
		return true;
	}

	// Sink
	unique_ptr<GlobalSinkState> GetGlobalSinkState(ClientContext &context) const override;
	SinkResultType Sink(ExecutionContext &context, DataChunk &chunk, OperatorSinkInput &input) const override;
	SinkFinalizeType Finalize(Pipeline &pipeline, Event &event, ClientContext &context,
	                          OperatorSinkFinalizeInput &input) const override;

	bool IsSink() const override {
		return true;
	}

	bool ParallelSink() const override {
		return false;
	}

	string GetName() const override {
		return "ENA_SAMPLES_INSERT";
	}
};

} // namespace duckdb
