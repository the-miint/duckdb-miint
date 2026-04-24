#include "read_ena_attributes.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/planner/expression/bound_columnref_expression.hpp"
#include "duckdb/planner/expression/bound_comparison_expression.hpp"
#include "duckdb/planner/expression/bound_conjunction_expression.hpp"
#include "duckdb/planner/expression/bound_constant_expression.hpp"
#include "duckdb/planner/expression/bound_operator_expression.hpp"

namespace duckdb {

namespace {

// Walk a DuckDB bound filter expression into an ENAFilterNode. Unknown shapes
// map to ENAFilterNode::Kind::UNSUPPORTED, which the extractor then rejects
// (forcing XML fallback).
//
// Column-name resolution invariant (relied upon here):
//   When DuckDB calls `pushdown_complex_filter`, every
//   `BoundColumnRefExpression::binding.column_index` in the supplied filter
//   tree is a LOCAL index into `LogicalGet::GetColumnIds()` — i.e. an index
//   into the projected-columns subset, not into the full schema `names`.
//   To recover a name we do: `names[GetColumnIds()[local_idx].GetPrimaryIndex()]`.
//
//   Source of truth: `duckdb/src/common/multi_file/multi_file_list.cpp:68-69`,
//   where DuckDB itself constructs filter column-refs with exactly this
//   binding shape (`ColumnBinding(table_index, entry.first)` where `entry.first`
//   is the local index into column_ids). If a future DuckDB version breaks
//   this invariant, the out-of-bounds checks below cause us to emit
//   `MakeUnsupported` and silently fall back to the XML path — degraded but
//   still correct.
std::string ResolveColumnName(const BoundColumnRefExpression &ref, const LogicalGet &get) {
	const auto &col_ids = get.GetColumnIds();
	idx_t local_idx = ref.binding.column_index;
	if (local_idx >= col_ids.size()) {
		return "";
	}
	idx_t actual_col = col_ids[local_idx].GetPrimaryIndex();
	if (actual_col >= get.names.size()) {
		return "";
	}
	return get.names[actual_col];
}

std::unique_ptr<miint::ENAFilterNode> ExpressionToFilterNode(const Expression &expr, const LogicalGet &get);

// Handle `col = const` or `const = col`. Only VARCHAR constants are mapped;
// anything else (NULL, non-string) => unsupported.
std::unique_ptr<miint::ENAFilterNode> TranslateEqual(const BoundComparisonExpression &cmp, const LogicalGet &get) {
	const Expression *col_side = nullptr;
	const Expression *const_side = nullptr;
	if (cmp.left->GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF &&
	    cmp.right->GetExpressionClass() == ExpressionClass::BOUND_CONSTANT) {
		col_side = cmp.left.get();
		const_side = cmp.right.get();
	} else if (cmp.right->GetExpressionClass() == ExpressionClass::BOUND_COLUMN_REF &&
	           cmp.left->GetExpressionClass() == ExpressionClass::BOUND_CONSTANT) {
		col_side = cmp.right.get();
		const_side = cmp.left.get();
	} else {
		return miint::ENAFilterNode::MakeUnsupported();
	}

	const auto &col_ref = col_side->Cast<BoundColumnRefExpression>();
	const auto &const_expr = const_side->Cast<BoundConstantExpression>();
	if (const_expr.value.IsNull()) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	if (const_expr.return_type.id() != LogicalTypeId::VARCHAR) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	auto col_name = ResolveColumnName(col_ref, get);
	if (col_name.empty()) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	return miint::ENAFilterNode::MakeEqual(col_name, const_expr.value.ToString());
}

// Handle `col IN (c1, c2, ...)`. Children[0] is the column ref; remaining
// children are the values. Anything else in the IN list (NULL, non-string
// constant, subexpression) => unsupported.
std::unique_ptr<miint::ENAFilterNode> TranslateIn(const BoundOperatorExpression &op, const LogicalGet &get) {
	if (op.children.size() < 2) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	if (op.children[0]->GetExpressionClass() != ExpressionClass::BOUND_COLUMN_REF) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	const auto &col_ref = op.children[0]->Cast<BoundColumnRefExpression>();
	auto col_name = ResolveColumnName(col_ref, get);
	if (col_name.empty()) {
		return miint::ENAFilterNode::MakeUnsupported();
	}
	std::vector<std::string> values;
	values.reserve(op.children.size() - 1);
	for (idx_t i = 1; i < op.children.size(); i++) {
		if (op.children[i]->GetExpressionClass() != ExpressionClass::BOUND_CONSTANT) {
			return miint::ENAFilterNode::MakeUnsupported();
		}
		const auto &c = op.children[i]->Cast<BoundConstantExpression>();
		if (c.value.IsNull() || c.return_type.id() != LogicalTypeId::VARCHAR) {
			return miint::ENAFilterNode::MakeUnsupported();
		}
		values.push_back(c.value.ToString());
	}
	return miint::ENAFilterNode::MakeIn(col_name, values);
}

std::unique_ptr<miint::ENAFilterNode> TranslateConjunction(const BoundConjunctionExpression &conj,
                                                           const LogicalGet &get, bool is_and) {
	std::vector<std::unique_ptr<miint::ENAFilterNode>> children;
	children.reserve(conj.children.size());
	for (const auto &child : conj.children) {
		children.push_back(ExpressionToFilterNode(*child, get));
	}
	return is_and ? miint::ENAFilterNode::MakeAnd(std::move(children))
	              : miint::ENAFilterNode::MakeOr(std::move(children));
}

std::unique_ptr<miint::ENAFilterNode> ExpressionToFilterNode(const Expression &expr, const LogicalGet &get) {
	switch (expr.GetExpressionClass()) {
	case ExpressionClass::BOUND_COMPARISON: {
		const auto &cmp = expr.Cast<BoundComparisonExpression>();
		if (cmp.GetExpressionType() != ExpressionType::COMPARE_EQUAL) {
			return miint::ENAFilterNode::MakeUnsupported();
		}
		return TranslateEqual(cmp, get);
	}
	case ExpressionClass::BOUND_OPERATOR: {
		const auto &op = expr.Cast<BoundOperatorExpression>();
		// COMPARE_NOT_IN is a distinct expression type and correctly lands
		// here as UNSUPPORTED. BOUND_OPERATOR also covers other operators
		// (IS NULL, IS NOT NULL, COALESCE, ...) that we reject uniformly.
		if (op.GetExpressionType() != ExpressionType::COMPARE_IN) {
			return miint::ENAFilterNode::MakeUnsupported();
		}
		// `COMPARE_IN` with a non-constant RHS (e.g., a correlated subquery
		// that hasn't been decorrelated yet) is handled inside TranslateIn:
		// any non-`BOUND_CONSTANT` child returns UNSUPPORTED.
		return TranslateIn(op, get);
	}
	case ExpressionClass::BOUND_CONJUNCTION: {
		const auto &conj = expr.Cast<BoundConjunctionExpression>();
		if (conj.GetExpressionType() == ExpressionType::CONJUNCTION_AND) {
			return TranslateConjunction(conj, get, /*is_and=*/true);
		}
		if (conj.GetExpressionType() == ExpressionType::CONJUNCTION_OR) {
			return TranslateConjunction(conj, get, /*is_and=*/false);
		}
		return miint::ENAFilterNode::MakeUnsupported();
	}
	default:
		return miint::ENAFilterNode::MakeUnsupported();
	}
}

} // namespace

// ---- Data ----

ReadENAAttributesTableFunction::Data::Data(std::vector<std::string> accessions)
    : accessions(std::move(accessions)), names({"sample_accession", "tag", "value"}),
      types({LogicalType::VARCHAR, LogicalType::VARCHAR, LogicalType::VARCHAR}) {
}

// ---- GlobalState ----

ReadENAAttributesTableFunction::GlobalState::GlobalState(DatabaseInstance &db)
    : client(make_uniq<miint::ENAClient>(db)) {
}

std::vector<std::string>
ReadENAAttributesTableFunction::GlobalState::ResolveSampleAccessions(const std::vector<std::string> &accessions) {
	std::vector<std::string> sample_accs;

	for (const auto &acc : accessions) {
		auto acc_type = miint::ENAParser::DetectAccessionType(acc);

		if (acc_type == miint::ENAAccessionType::SAMPLE) {
			sample_accs.push_back(acc);
		} else if (acc_type == miint::ENAAccessionType::STUDY) {
			// Get all sample accessions for the study
			auto tsv = client->Search(acc, "sample", "sample_accession");
			auto parsed = miint::ENAParser::ParseTSV(tsv);
			for (size_t i = 0; i < parsed.column_names.size(); i++) {
				if (parsed.column_names[i] == "sample_accession") {
					for (const auto &row : parsed.rows) {
						if (i < row.size() && !row[i].empty()) {
							sample_accs.push_back(row[i]);
						}
					}
					break;
				}
			}
		} else {
			// Run or experiment: resolve to sample_accession via read_run
			auto tsv = client->Search(acc, "read_run", "sample_accession");
			auto parsed = miint::ENAParser::ParseTSV(tsv);
			for (size_t i = 0; i < parsed.column_names.size(); i++) {
				if (parsed.column_names[i] == "sample_accession") {
					for (const auto &row : parsed.rows) {
						if (i < row.size() && !row[i].empty()) {
							sample_accs.push_back(row[i]);
						}
					}
					break;
				}
			}
		}
	}

	return sample_accs;
}

void ReadENAAttributesTableFunction::GlobalState::ResolveAccessions(const std::vector<std::string> &accessions,
                                                                    bool has_pushdown) {
	// Mark as resolved immediately so that a network error during the
	// Search call below doesn't cause successive Execute callbacks to keep
	// retrying the same failing request in a loop.
	resolved = true;
	if (has_pushdown) {
		// Split inputs by type. Studies stay as studies — the structured query
		// runs directly against `study_accession="..."` and never enumerates
		// the study's samples on the client side. Samples pass through unchanged.
		// Runs / experiments still need one /search?result=read_run call each to
		// resolve to sample_accession (same as the legacy resolver), because the
		// structured search endpoint accepts only sample/study/run columns as
		// top-level filters on `result=sample`.
		for (const auto &acc : accessions) {
			auto acc_type = miint::ENAParser::DetectAccessionType(acc);
			if (acc_type == miint::ENAAccessionType::STUDY) {
				pending_study_accs.push_back(acc);
			} else if (acc_type == miint::ENAAccessionType::SAMPLE) {
				pending_sample_accs.push_back(acc);
			} else {
				auto tsv = client->Search(acc, "read_run", "sample_accession");
				auto parsed = miint::ENAParser::ParseTSV(tsv);
				for (size_t i = 0; i < parsed.column_names.size(); i++) {
					if (parsed.column_names[i] == "sample_accession") {
						for (const auto &row : parsed.rows) {
							if (i < row.size() && !row[i].empty()) {
								pending_sample_accs.push_back(row[i]);
							}
						}
						break;
					}
				}
			}
		}
	} else {
		pending_sample_accs = ResolveSampleAccessions(accessions);
	}
	// Progress denominator: studies count as 1 unit each (their batch is a
	// single HTTP call returning the filtered set); samples contribute their
	// raw count (advanced by BATCH_SIZE per HTTP call). This isn't perfectly
	// linear but terminates at 100%.
	total_samples_expected.store(pending_sample_accs.size() + pending_study_accs.size(), std::memory_order_relaxed);
	samples_fetched.store(0, std::memory_order_relaxed);
}

bool ReadENAAttributesTableFunction::GlobalState::FetchNextBatch() {
	if (pending_sample_offset >= pending_sample_accs.size()) {
		return false;
	}
	size_t remaining = pending_sample_accs.size() - pending_sample_offset;
	size_t n = std::min<size_t>(BATCH_SIZE, remaining);
	std::vector<std::string> batch(pending_sample_accs.begin() + pending_sample_offset,
	                               pending_sample_accs.begin() + pending_sample_offset + n);
	pending_sample_offset += n;

	auto xml = client->FetchXML(batch);
	current_batch = miint::ENAParser::ParseSampleAttributesXML(xml);
	current_batch_offset = 0;
	samples_fetched.fetch_add(n, std::memory_order_relaxed);
	return true;
}

bool ReadENAAttributesTableFunction::GlobalState::FetchNextStructuredBatch(
    const miint::ENAAttributePushdown &pushdown) {
	// Drain study-direct batches first: each one runs
	// `study_accession IN (...) AND <pushdown>` directly against the /search
	// endpoint, returning only the matching samples without client-side
	// enumeration. Batched so multiple study inputs share one request.
	if (pending_study_offset < pending_study_accs.size()) {
		size_t remaining = pending_study_accs.size() - pending_study_offset;
		size_t n = std::min<size_t>(BATCH_SIZE, remaining);
		std::vector<std::string> batch(pending_study_accs.begin() + pending_study_offset,
		                               pending_study_accs.begin() + pending_study_offset + n);
		pending_study_offset += n;

		auto url = miint::BuildStudyDirectSearchURL(batch, pushdown.tags, pushdown.tag_value_pairs);
		auto tsv = client->FetchURL(url);
		auto parsed = miint::ENAParser::ParseTSV(tsv);
		current_batch = miint::UnpivotStructuredTSV(parsed, pushdown.tags);
		current_batch_offset = 0;
		samples_fetched.fetch_add(n, std::memory_order_relaxed);
		return true;
	}

	if (pending_sample_offset >= pending_sample_accs.size()) {
		return false;
	}
	size_t remaining = pending_sample_accs.size() - pending_sample_offset;
	size_t n = std::min<size_t>(BATCH_SIZE, remaining);
	std::vector<std::string> batch(pending_sample_accs.begin() + pending_sample_offset,
	                               pending_sample_accs.begin() + pending_sample_offset + n);
	pending_sample_offset += n;

	auto url = miint::BuildStructuredSearchURL(batch, pushdown.tags, pushdown.tag_value_pairs);
	auto tsv = client->FetchURL(url);
	auto parsed = miint::ENAParser::ParseTSV(tsv);
	current_batch = miint::UnpivotStructuredTSV(parsed, pushdown.tags);
	current_batch_offset = 0;
	samples_fetched.fetch_add(n, std::memory_order_relaxed);
	return true;
}

// ---- Bind ----

unique_ptr<FunctionData> ReadENAAttributesTableFunction::Bind(ClientContext &context, TableFunctionBindInput &input,
                                                              vector<LogicalType> &return_types,
                                                              vector<std::string> &names) {
	std::vector<std::string> accessions;
	if (input.inputs[0].type().id() == LogicalTypeId::VARCHAR) {
		accessions.push_back(input.inputs[0].ToString());
	} else if (input.inputs[0].type().id() == LogicalTypeId::LIST) {
		auto &list_children = ListValue::GetChildren(input.inputs[0]);
		for (const auto &child : list_children) {
			accessions.push_back(child.ToString());
		}
	} else {
		throw InvalidInputException("read_ena_attributes: first argument must be VARCHAR or VARCHAR[]");
	}

	if (accessions.empty()) {
		throw InvalidInputException("read_ena_attributes: at least one accession must be provided");
	}
	for (auto &acc : accessions) {
		// Trim whitespace
		auto start = acc.find_first_not_of(" \t\n\r");
		if (start == std::string::npos) {
			acc.clear();
		} else {
			acc = acc.substr(start, acc.find_last_not_of(" \t\n\r") - start + 1);
		}
		if (acc.empty()) {
			throw InvalidInputException("read_ena_attributes: accession cannot be empty");
		}
	}

	auto data = make_uniq<Data>(std::move(accessions));

	for (auto &name : data->names) {
		names.emplace_back(name);
	}
	for (auto &type : data->types) {
		return_types.emplace_back(type);
	}

	return data;
}

// ---- InitGlobal / InitLocal ----

unique_ptr<GlobalTableFunctionState> ReadENAAttributesTableFunction::InitGlobal(ClientContext &context,
                                                                                TableFunctionInitInput &input) {
	auto &db = DatabaseInstance::GetDatabase(context);
	return make_uniq<GlobalState>(db);
}

unique_ptr<LocalTableFunctionState> ReadENAAttributesTableFunction::InitLocal(ExecutionContext &context,
                                                                              TableFunctionInitInput &input,
                                                                              GlobalTableFunctionState *global_state) {
	return make_uniq<LocalState>();
}

// ---- Execute ----

void ReadENAAttributesTableFunction::Execute(ClientContext &context, TableFunctionInput &data_p, DataChunk &output) {
	auto &bind_data = data_p.bind_data->Cast<Data>();
	auto &global_state = data_p.global_state->Cast<GlobalState>();

	lock_guard<mutex> guard(global_state.lock);

	// First call: resolve input accessions into the right pending queues.
	// When pushdown is active, study inputs stay as studies (drained via the
	// study-direct path in FetchNextStructuredBatch); otherwise studies get
	// enumerated into sample accessions up front. Subsequent fetches happen
	// lazily one batch at a time below.
	const bool use_structured = !bind_data.pushdown.tags.empty();
	if (!global_state.resolved) {
		global_state.ResolveAccessions(bind_data.accessions, use_structured);
	}

	// Advance to the next fetched batch if we've drained the current one.
	// Loop to skip empty batches (a 50-sample batch whose XML has no sample
	// with any attributes would otherwise end the scan prematurely). When
	// pushdown_complex_filter extracted a structured predicate, use the TSV
	// search path; otherwise fall back to per-sample XML.
	while (global_state.current_batch_offset >= global_state.current_batch.size()) {
		bool ok =
		    use_structured ? global_state.FetchNextStructuredBatch(bind_data.pushdown) : global_state.FetchNextBatch();
		if (!ok) {
			output.SetCardinality(0);
			return;
		}
	}

	auto &attrs = global_state.current_batch;
	size_t offset = global_state.current_batch_offset;
	idx_t remaining = attrs.size() - offset;
	idx_t count = MinValue<idx_t>(remaining, STANDARD_VECTOR_SIZE);

	// sample_accession (column 0)
	auto acc_data = FlatVector::GetData<string_t>(output.data[0]);
	for (idx_t i = 0; i < count; i++) {
		acc_data[i] = StringVector::AddString(output.data[0], attrs[offset + i].sample_accession);
	}

	// tag (column 1)
	auto tag_data = FlatVector::GetData<string_t>(output.data[1]);
	for (idx_t i = 0; i < count; i++) {
		tag_data[i] = StringVector::AddString(output.data[1], attrs[offset + i].tag);
	}

	// value (column 2)
	auto val_data = FlatVector::GetData<string_t>(output.data[2]);
	auto &val_validity = FlatVector::Validity(output.data[2]);
	for (idx_t i = 0; i < count; i++) {
		const auto &val = attrs[offset + i].value;
		if (val.empty()) {
			val_validity.SetInvalid(i);
		} else {
			val_data[i] = StringVector::AddString(output.data[2], val);
		}
	}

	global_state.current_batch_offset += count;
	output.SetCardinality(count);
}

// ---- Pushdown complex filter ----

void ReadENAAttributesTableFunction::PushdownComplexFilter(ClientContext &context, LogicalGet &get,
                                                           FunctionData *bind_data_p,
                                                           vector<unique_ptr<Expression>> &filters) {
	auto &bind_data = bind_data_p->Cast<Data>();

	std::vector<std::unique_ptr<miint::ENAFilterNode>> ast;
	ast.reserve(filters.size());
	for (const auto &expr : filters) {
		ast.push_back(ExpressionToFilterNode(*expr, get));
	}
	bind_data.pushdown = miint::ExtractPushdownPredicates(ast);

	// Per plan decision #5 (localdocs/PLAN-ena-predicate-maxseqs.md): we do
	// NOT erase any entries from `filters`. DuckDB re-applies the predicates
	// as a LogicalFilter above the scan, so any semantic divergence between
	// ENA's /search matcher and SQL equality degrades to a mild inefficiency,
	// not a correctness bug.
}

// ---- Progress ----

double ReadENAAttributesTableFunction::Progress(ClientContext &context, const FunctionData *bind_data,
                                                const GlobalTableFunctionState *global_state) {
	if (!global_state) {
		return -1.0;
	}
	auto &state = global_state->Cast<GlobalState>();
	// Before ResolveAccessions runs we don't know the scale of the work. Once
	// resolved, report samples-fetched / total. This gives users visible
	// progress for large studies like PRJEB11419 (33k+ samples) where a full
	// scan is long-running by nature.
	size_t total = state.total_samples_expected.load(std::memory_order_relaxed);
	if (total == 0) {
		return -1.0;
	}
	size_t fetched = state.samples_fetched.load(std::memory_order_relaxed);
	double pct = 100.0 * static_cast<double>(fetched) / static_cast<double>(total);
	return std::min(100.0, pct);
}

// ---- Registration ----

TableFunction ReadENAAttributesTableFunction::GetFunction() {
	auto tf = TableFunction("read_ena_attributes", {LogicalType::ANY}, Execute, Bind, InitGlobal, InitLocal);
	tf.table_scan_progress = Progress;
	tf.pushdown_complex_filter = PushdownComplexFilter;
	return tf;
}

void ReadENAAttributesTableFunction::Register(ExtensionLoader &loader) {
	loader.RegisterFunction(GetFunction());
}

} // namespace duckdb
