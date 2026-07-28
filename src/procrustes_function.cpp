#include "procrustes_table_function.hpp"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include "catalog_utils.hpp"
#include "id_column_utils.hpp"
#include "procrustes_core.hpp"

#include "duckdb/common/exception.hpp"
#include "duckdb/common/string_util.hpp"
#include "duckdb/common/types.hpp"
#include "duckdb/common/types/data_chunk.hpp"
#include "duckdb/common/vector_size.hpp"
#include "duckdb/function/table_function.hpp"
#include "duckdb/main/client_context.hpp"
#include "duckdb/main/connection.hpp"
#include "duckdb/main/database.hpp"

namespace duckdb {
namespace {

using ::miint::procrustes::ApplyToOther;
using ::miint::procrustes::ApplyToReference;
using ::miint::procrustes::Disparity;
using ::miint::procrustes::FitProcrustes;
using ::miint::procrustes::MonteCarloPValue;
using ::miint::procrustes::ProcrustesFit;

// ── id-type resolution ──────────────────────────────────────────────────────
// Mirror BIGINT/UUID onto the output; everything else collapses to VARCHAR. Kept
// local (rather than pulling in the UniFrac-gated ResolveSampleIdOutputType) so
// procrustes stays UniFrac-free. sample_id handling supports VARCHAR/BIGINT/UUID.
LogicalType ResolveIdType(const LogicalType &t) {
	if (t.id() == LogicalTypeId::BIGINT || t.id() == LogicalTypeId::UUID) {
		return t;
	}
	return LogicalType::VARCHAR;
}

// ── input readers ────────────────────────────────────────────────────────────
// A long-format ordination (sample_id, axis, coordinate) reshaped into a dense
// n×d row-major matrix. `sample_ids` is the sorted, distinct id dictionary (row
// order of `matrix`); `sample_id_type` is the resolved output id type.
struct OrdinationTable {
	std::vector<std::string> sample_ids;
	std::unordered_map<std::string, uint32_t> index;
	uint32_t n = 0;
	uint32_t d = 0;
	std::vector<double> matrix; // n*d row-major
	LogicalType sample_id_type = LogicalType::VARCHAR;

	const double *Row(uint32_t r) const {
		return matrix.data() + static_cast<size_t>(r) * d;
	}
};

OrdinationTable ReadOrdinationTable(ClientContext &context, const std::string &table_name, const std::string &role,
                                    const std::string &caller_name) {
	if (table_name.empty()) {
		throw BinderException("%s: %s ordination-table name must not be empty", caller_name, role);
	}
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const std::string select = "SELECT sample_id::VARCHAR, axis::INTEGER, coordinate::DOUBLE FROM " + qname;

	// LIMIT 0 probe — surface a missing column / bad cast as a bind-time error
	// before materializing (mirrors ReadDistanceTable).
	auto probe = conn.Query(select + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException(
		    "%s: %s ordination-table '%s' must expose (sample_id, axis INTEGER, coordinate DOUBLE): %s", caller_name,
		    role, table_name, probe->GetError());
	}

	// Capture sample_id's original SQL type (before the ::VARCHAR cast) so the
	// output column can mirror BIGINT/UUID ids. The probe proved it exists.
	LogicalType raw_id_type = LogicalType::VARCHAR;
	{
		auto cols = GetTableOrViewColumns(context, table_name, "ordination-table");
		for (idx_t i = 0; i < cols.names.size(); ++i) {
			if (StringUtil::Lower(cols.names[i]) == "sample_id") {
				raw_id_type = cols.types[i];
			}
		}
	}

	auto result = conn.Query(select);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read %s ordination-table '%s': %s", caller_name, role, table_name,
		                            result->GetError());
	}

	struct Entry {
		std::string id;
		int32_t axis;
		double coord;
	};
	std::vector<Entry> entries;
	std::vector<std::string> ids_raw;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat id_u, ax_u, co_u;
		chunk->data[0].ToUnifiedFormat(rn, id_u);
		chunk->data[1].ToUnifiedFormat(rn, ax_u);
		chunk->data[2].ToUnifiedFormat(rn, co_u);
		auto id_data = UnifiedVectorFormat::GetData<string_t>(id_u);
		auto ax_data = UnifiedVectorFormat::GetData<int32_t>(ax_u);
		auto co_data = UnifiedVectorFormat::GetData<double>(co_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ii = id_u.sel->get_index(i);
			const auto ai = ax_u.sel->get_index(i);
			const auto ci = co_u.sel->get_index(i);
			// Any NULL in a (sample_id, axis, coordinate) row makes the cell
			// undefined — drop the row; the completeness check below turns a
			// resulting hole into a named error rather than a silent zero.
			if (!id_u.validity.RowIsValid(ii) || !ax_u.validity.RowIsValid(ai) || !co_u.validity.RowIsValid(ci)) {
				continue;
			}
			const int32_t axis = ax_data[ai];
			if (axis < 0) {
				throw InvalidInputException("%s: %s ordination-table '%s' has a negative axis (%d); axes are 0-indexed",
				                            caller_name, role, table_name, axis);
			}
			std::string id = id_data[ii].GetString();
			ids_raw.push_back(id);
			entries.push_back({std::move(id), axis, co_data[ci]});
		}
	}

	OrdinationTable out;
	out.sample_id_type = ResolveIdType(raw_id_type);

	// Sorted, distinct id dictionary (deterministic row order; also lets the
	// full-overlap caller compare two id sets by ==).
	std::vector<std::string> sample_ids = std::move(ids_raw);
	std::sort(sample_ids.begin(), sample_ids.end());
	sample_ids.erase(std::unique(sample_ids.begin(), sample_ids.end()), sample_ids.end());
	const auto n = static_cast<uint32_t>(sample_ids.size());
	if (n == 0) {
		throw InvalidInputException("%s: %s ordination-table '%s' is empty after dropping NULL rows", caller_name, role,
		                            table_name);
	}
	out.index.reserve(sample_ids.size());
	for (uint32_t i = 0; i < n; ++i) {
		out.index.emplace(sample_ids[i], i);
	}

	// d = max axis + 1; then require every sample to carry a full, gap-free
	// 0..d-1 axis set exactly once (fail loud on ragged / duplicated input).
	// int64 accumulator: e.axis can be INT32_MAX (e.g. the wrong column mapped to
	// `axis`), and `max_axis + 1` in int32 would be signed overflow (UB).
	int64_t max_axis = -1;
	for (const auto &e : entries) {
		max_axis = std::max<int64_t>(max_axis, e.axis);
	}
	// An ordination of n points has at most n axes; a larger max axis means the
	// input is malformed (usually a non-axis column mapped to `axis`). Reject here,
	// before sizing an n*d buffer, rather than attempting a wild allocation.
	if (max_axis >= static_cast<int64_t>(n)) {
		throw InvalidInputException(
		    "%s: %s ordination-table '%s' has axis %lld but only %u sample(s); an ordination of n points has at most n "
		    "axes (is the wrong column mapped to `axis`?)",
		    caller_name, role, table_name, static_cast<long long>(max_axis), n);
	}
	const auto d = static_cast<uint32_t>(max_axis + 1);
	std::vector<double> matrix(static_cast<size_t>(n) * d, 0.0);
	std::vector<bool> filled(static_cast<size_t>(n) * d, false);
	for (const auto &e : entries) {
		const uint32_t r = out.index.at(e.id);
		const size_t pos = static_cast<size_t>(r) * d + static_cast<uint32_t>(e.axis);
		if (filled[pos]) {
			throw InvalidInputException("%s: %s ordination-table '%s' has duplicate (sample_id='%s', axis=%d) entries",
			                            caller_name, role, table_name, e.id, e.axis);
		}
		filled[pos] = true;
		matrix[pos] = e.coord;
	}
	for (uint32_t r = 0; r < n; ++r) {
		for (uint32_t a = 0; a < d; ++a) {
			if (!filled[static_cast<size_t>(r) * d + a]) {
				throw InvalidInputException(
				    "%s: %s ordination-table '%s' is ragged — sample '%s' is missing axis %u (expected a full 0..%u "
				    "axis set)",
				    caller_name, role, table_name, sample_ids[r], a, d - 1);
			}
		}
	}

	out.sample_ids = std::move(sample_ids);
	out.n = n;
	out.d = d;
	out.matrix = std::move(matrix);
	return out;
}

// Read a pairing table (reference_id, other_id) as string pairs; NULL rows are
// dropped. Deduplicated on reference_id and sorted by it (mirrors q2's
// `sorted(set(pairing.index) & ...)` ordering).
std::vector<std::pair<std::string, std::string>> ReadPairing(ClientContext &context, const std::string &table_name,
                                                             const std::string &caller_name) {
	auto &db = DatabaseInstance::GetDatabase(context);
	Connection conn(db);
	const auto qname = KeywordHelper::WriteOptionallyQuoted(table_name);
	const std::string select = "SELECT reference_id::VARCHAR, other_id::VARCHAR FROM " + qname;

	auto probe = conn.Query(select + " LIMIT 0");
	if (probe->HasError()) {
		throw InvalidInputException("%s: pairing-table '%s' must expose (reference_id, other_id): %s", caller_name,
		                            table_name, probe->GetError());
	}
	auto result = conn.Query(select);
	if (result->HasError()) {
		throw InvalidInputException("%s: failed to read pairing-table '%s': %s", caller_name, table_name,
		                            result->GetError());
	}

	std::vector<std::pair<std::string, std::string>> pairs;
	auto &materialized = result->Cast<MaterializedQueryResult>();
	while (auto chunk = materialized.Fetch()) {
		const idx_t rn = chunk->size();
		if (rn == 0) {
			break;
		}
		UnifiedVectorFormat r_u, o_u;
		chunk->data[0].ToUnifiedFormat(rn, r_u);
		chunk->data[1].ToUnifiedFormat(rn, o_u);
		auto r_data = UnifiedVectorFormat::GetData<string_t>(r_u);
		auto o_data = UnifiedVectorFormat::GetData<string_t>(o_u);
		for (idx_t i = 0; i < rn; ++i) {
			const auto ri = r_u.sel->get_index(i);
			const auto oi = o_u.sel->get_index(i);
			if (!r_u.validity.RowIsValid(ri) || !o_u.validity.RowIsValid(oi)) {
				continue; // an unpaired row carries no anchor
			}
			pairs.emplace_back(r_data[ri].GetString(), o_data[oi].GetString());
		}
	}
	// Enforce a 1:1 anchor pairing. A repeated reference_id would silently drop a
	// row; a repeated other_id would feed one physical `other` sample into the fit
	// as several independent anchors, silently corrupting the transform. Reject
	// either, mirroring ReadOrdinationTable's duplicate-(sample_id, axis) throw.
	std::unordered_set<std::string> seen_ref, seen_other;
	for (const auto &p : pairs) {
		if (!seen_ref.insert(p.first).second) {
			throw InvalidInputException(
			    "%s: pairing-table '%s' has a duplicate reference_id '%s'; the anchor pairing must be 1:1", caller_name,
			    table_name, p.first);
		}
		if (!seen_other.insert(p.second).second) {
			throw InvalidInputException(
			    "%s: pairing-table '%s' has a duplicate other_id '%s'; the anchor pairing must be 1:1", caller_name,
			    table_name, p.second);
		}
	}
	std::sort(pairs.begin(), pairs.end());
	return pairs;
}

// ── carrier / states ─────────────────────────────────────────────────────────
struct ProcrustesRow {
	bool is_other; // false = reference cloud, true = other cloud
	std::string sample_id;
	int32_t axis;
	double coordinate;
};

struct ProcrustesData : public TableFunctionData {
	std::vector<ProcrustesRow> rows;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	double m2 = 0.0;
	double pvalue = std::numeric_limits<double>::quiet_NaN(); // NaN -> NULL
};

struct ProcrustesGlobalState : public GlobalTableFunctionState {
	std::vector<ProcrustesRow> rows;
	size_t cursor = 0;
	LogicalType sample_id_type = LogicalType::VARCHAR;
	double m2 = 0.0;
	double pvalue = std::numeric_limits<double>::quiet_NaN();
	idx_t MaxThreads() const override {
		return 1;
	}
};

void DeclareProcrustesOutputSchema(const LogicalType &sample_id_type, vector<LogicalType> &return_types,
                                   vector<string> &names) {
	names.emplace_back("matrix");
	return_types.emplace_back(LogicalType::VARCHAR);
	names.emplace_back("sample_id");
	return_types.emplace_back(sample_id_type);
	names.emplace_back("axis");
	return_types.emplace_back(LogicalType::INTEGER);
	names.emplace_back("coordinate");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("m2");
	return_types.emplace_back(LogicalType::DOUBLE);
	names.emplace_back("pvalue");
	return_types.emplace_back(LogicalType::DOUBLE);
}

// First `d` columns of ordination row `r`, appended to `dst`.
void AppendUsedRow(const OrdinationTable &ord, uint32_t r, uint32_t d, std::vector<double> &dst) {
	const double *row = ord.Row(r);
	dst.insert(dst.end(), row, row + d);
}

// Emit one transformed cloud (`is_other`) as rows: `ids[i]` gets `coords[i*d+a]`.
void EmitCloud(bool is_other, const std::vector<std::string> &ids, const std::vector<double> &coords, uint32_t d,
               std::vector<ProcrustesRow> &out) {
	for (size_t i = 0; i < ids.size(); ++i) {
		for (uint32_t a = 0; a < d; ++a) {
			out.push_back({is_other, ids[i], static_cast<int32_t>(a), coords[i * d + a]});
		}
	}
}

uint64_t ResolveSeed(int32_t seed_param) {
	if (seed_param >= 0) {
		return static_cast<uint64_t>(seed_param);
	}
	std::random_device rd;
	return (static_cast<uint64_t>(rd()) << 32) ^ static_cast<uint64_t>(rd());
}

unique_ptr<FunctionData> ProcrustesBind(ClientContext &context, TableFunctionBindInput &input,
                                        vector<LogicalType> &return_types, vector<string> &names) {
	const std::string ref_name = input.inputs[0].GetValue<string>();
	const std::string other_name = input.inputs[1].GetValue<string>();

	std::string pairing_name;
	int32_t n_dims_param = 0; // 0 = use all available axes
	int32_t permutations = 999;
	int32_t seed_param = -1; // <0 = nondeterministic (mirrors pcoa's seed convention)
	for (const auto &kv : input.named_parameters) {
		const auto key = StringUtil::Lower(kv.first);
		if (key == "pairing") {
			pairing_name = kv.second.GetValue<string>();
		} else if (key == "n_dims") {
			n_dims_param = kv.second.GetValue<int32_t>();
		} else if (key == "permutations") {
			permutations = kv.second.GetValue<int32_t>();
		} else if (key == "seed") {
			seed_param = kv.second.GetValue<int32_t>();
		}
	}
	if (n_dims_param < 0) {
		throw BinderException("procrustes: n_dims must be >= 1 (got %d)", n_dims_param);
	}
	if (permutations < 0) {
		throw BinderException("procrustes: permutations must be >= 0 (got %d)", permutations);
	}

	auto ref = ReadOrdinationTable(context, ref_name, "reference", "procrustes");
	auto other = ReadOrdinationTable(context, other_name, "other", "procrustes");

	auto data = make_uniq<ProcrustesData>();
	// Output id type: shared BIGINT/UUID only if both inputs agree, else VARCHAR.
	// This deliberately DIVERGES from ReadDistanceTable's fail-loud on mismatched
	// ids: there sample_a/sample_b share one id space (a mismatch is a bug), but a
	// procrustes reference and other are legitimately different id spaces joined by
	// string equality (full) or an explicit pairing (partial), so collapsing a
	// mixed pair to VARCHAR is the correct, non-lossy behavior, not an error.
	data->sample_id_type = (ref.sample_id_type == other.sample_id_type) ? ref.sample_id_type : LogicalType::VARCHAR;

	const bool partial = !pairing_name.empty();
	uint32_t d;
	if (partial) {
		d = n_dims_param > 0 ? static_cast<uint32_t>(n_dims_param) : std::min(ref.d, other.d);
		if (d > ref.d || d > other.d) {
			throw BinderException("procrustes: n_dims (%u) exceeds available axes (reference=%u, other=%u)", d, ref.d,
			                      other.d);
		}
	} else {
		if (ref.d != other.d) {
			throw BinderException(
			    "procrustes: reference and other have different axis counts (%u vs %u); a full procrustes needs the "
			    "same ordination shape (supply a pairing for anchored alignment)",
			    ref.d, other.d);
		}
		d = n_dims_param > 0 ? static_cast<uint32_t>(n_dims_param) : ref.d;
		if (d > ref.d) {
			throw BinderException("procrustes: n_dims (%u) exceeds available axes (%u)", d, ref.d);
		}
	}
	if (d < 1) {
		throw BinderException("procrustes: need at least 1 axis");
	}

	if (partial) {
		auto pairs = ReadPairing(context, pairing_name, "procrustes");
		std::vector<double> ref_anchor;
		std::vector<double> other_anchor;
		uint32_t n_anchor = 0;
		for (const auto &p : pairs) {
			auto rit = ref.index.find(p.first);
			auto oit = other.index.find(p.second);
			if (rit == ref.index.end() || oit == other.index.end()) {
				continue; // a pair whose ids aren't both present is not an anchor
			}
			AppendUsedRow(ref, rit->second, d, ref_anchor);
			AppendUsedRow(other, oit->second, d, other_anchor);
			++n_anchor;
		}
		if (n_anchor < d + 1) {
			throw BinderException(
			    "procrustes: only %u usable anchor pair(s) after matching the pairing to both ordinations; a "
			    "%u-dimensional fit needs at least %u",
			    n_anchor, d, d + 1);
		}

		ProcrustesFit fit;
		try {
			fit = FitProcrustes(ref_anchor.data(), other_anchor.data(), n_anchor, d);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("procrustes: %s", e.what());
		}

		// Apply the anchor-fit transform to every row of each ordination.
		std::vector<double> ref_all;
		std::vector<double> other_all;
		for (uint32_t r = 0; r < ref.n; ++r) {
			AppendUsedRow(ref, r, d, ref_all);
		}
		for (uint32_t r = 0; r < other.n; ++r) {
			AppendUsedRow(other, r, d, other_all);
		}
		std::vector<double> ref_std(ref_all.size());
		std::vector<double> other_fit(other_all.size());
		ApplyToReference(fit, ref_all.data(), ref.n, ref_std.data());
		ApplyToOther(fit, other_all.data(), other.n, other_fit.data());

		// M^2 is the disparity over the anchor pairs; the p-value is undefined for
		// the anchored case (q2's partial_procrustes runs no Monte Carlo test).
		std::vector<double> ref_anchor_std(ref_anchor.size());
		std::vector<double> other_anchor_fit(other_anchor.size());
		ApplyToReference(fit, ref_anchor.data(), n_anchor, ref_anchor_std.data());
		ApplyToOther(fit, other_anchor.data(), n_anchor, other_anchor_fit.data());
		data->m2 = Disparity(ref_anchor_std.data(), other_anchor_fit.data(), n_anchor, d);
		data->pvalue = std::numeric_limits<double>::quiet_NaN();

		data->rows.reserve(static_cast<size_t>(ref.n + other.n) * d);
		EmitCloud(/*is_other=*/false, ref.sample_ids, ref_std, d, data->rows);
		EmitCloud(/*is_other=*/true, other.sample_ids, other_fit, d, data->rows);
	} else {
		// Full overlap: both ordinations must describe the same samples.
		if (ref.sample_ids != other.sample_ids) {
			throw BinderException("procrustes: reference and other describe different sample sets; a full procrustes "
			                      "requires identical samples (supply a pairing for anchored alignment)");
		}
		const uint32_t n = ref.n;
		if (n < d + 1) {
			throw BinderException(
			    "procrustes: %u sample(s) is too few for a %u-dimensional fit; at least %u are required", n, d, d + 1);
		}
		// Same sorted id order in both, so row i lines up by construction.
		std::vector<double> ref_use;
		std::vector<double> other_use;
		ref_use.reserve(static_cast<size_t>(n) * d);
		other_use.reserve(static_cast<size_t>(n) * d);
		for (uint32_t r = 0; r < n; ++r) {
			AppendUsedRow(ref, r, d, ref_use);
			AppendUsedRow(other, r, d, other_use);
		}

		ProcrustesFit fit;
		try {
			fit = FitProcrustes(ref_use.data(), other_use.data(), n, d);
		} catch (const std::invalid_argument &e) {
			throw InvalidInputException("procrustes: %s", e.what());
		}
		std::vector<double> ref_std(ref_use.size());
		std::vector<double> other_fit(other_use.size());
		ApplyToReference(fit, ref_use.data(), n, ref_std.data());
		ApplyToOther(fit, other_use.data(), n, other_fit.data());
		data->m2 = Disparity(ref_std.data(), other_fit.data(), n, d);
		data->pvalue = permutations > 0 ? MonteCarloPValue(ref_use.data(), other_use.data(), n, d, data->m2,
		                                                   static_cast<uint32_t>(permutations), ResolveSeed(seed_param))
		                                : std::numeric_limits<double>::quiet_NaN();

		data->rows.reserve(static_cast<size_t>(2 * n) * d);
		EmitCloud(/*is_other=*/false, ref.sample_ids, ref_std, d, data->rows);
		EmitCloud(/*is_other=*/true, other.sample_ids, other_fit, d, data->rows);
	}

	DeclareProcrustesOutputSchema(data->sample_id_type, return_types, names);
	return std::move(data);
}

unique_ptr<GlobalTableFunctionState> ProcrustesInitGlobal(ClientContext &, TableFunctionInitInput &input) {
	// CastNoConst (not Cast): bind_data is optional_ptr<const FunctionData>, so
	// Cast<> yields a const ref and the std::move below would silently deep-copy
	// the whole row buffer via copy-assignment. Every other table function here
	// (unifrac_pcoa, woltka_ogu, massql, ...) uses CastNoConst for exactly this.
	auto &data = input.bind_data->CastNoConst<ProcrustesData>();
	auto state = make_uniq<ProcrustesGlobalState>();
	state->rows = std::move(data.rows);
	state->sample_id_type = data.sample_id_type;
	state->m2 = data.m2;
	state->pvalue = data.pvalue;
	return std::move(state);
}

void ProcrustesExecute(ClientContext &, TableFunctionInput &data_p, DataChunk &output) {
	auto &gstate = data_p.global_state->Cast<ProcrustesGlobalState>();
	const size_t total = gstate.rows.size();
	const idx_t n = MinValue<idx_t>(STANDARD_VECTOR_SIZE, total - gstate.cursor);

	auto matrix_data = FlatVector::GetData<string_t>(output.data[0]);
	auto &sample_id_vec = output.data[1];
	auto axis_data = FlatVector::GetData<int32_t>(output.data[2]);
	auto coord_data = FlatVector::GetData<double>(output.data[3]);
	auto m2_data = FlatVector::GetData<double>(output.data[4]);
	auto pvalue_data = FlatVector::GetData<double>(output.data[5]);
	auto &pvalue_validity = FlatVector::Validity(output.data[5]);
	const bool pvalue_is_null = std::isnan(gstate.pvalue);

	for (idx_t i = 0; i < n; ++i) {
		const auto &r = gstate.rows[gstate.cursor + i];
		matrix_data[i] = StringVector::AddString(output.data[0], r.is_other ? "other" : "reference");
		EmitIdCell(sample_id_vec, i, r.sample_id, gstate.sample_id_type);
		axis_data[i] = r.axis;
		coord_data[i] = r.coordinate;
		m2_data[i] = gstate.m2;
		if (pvalue_is_null) {
			pvalue_data[i] = 0.0;
			pvalue_validity.SetInvalid(i);
		} else {
			pvalue_data[i] = gstate.pvalue;
		}
	}
	gstate.cursor += n;
	output.SetCardinality(n);
}

} // namespace

void RegisterProcrustes(ExtensionLoader &loader) {
	TableFunction fn("procrustes", {LogicalType::VARCHAR, LogicalType::VARCHAR}, ProcrustesExecute, ProcrustesBind,
	                 ProcrustesInitGlobal);
	fn.named_parameters["pairing"] = LogicalType::VARCHAR;
	fn.named_parameters["n_dims"] = LogicalType::INTEGER;
	fn.named_parameters["permutations"] = LogicalType::INTEGER;
	fn.named_parameters["seed"] = LogicalType::INTEGER;
	loader.RegisterFunction(fn);
}

} // namespace duckdb
