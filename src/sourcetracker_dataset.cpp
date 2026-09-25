#include "sourcetracker_dataset.hpp"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <unordered_map>
#include <unordered_set>

namespace miint::sourcetracker {

using miint::unifrac::CooRow;
using miint::unifrac::MetadataRow;

namespace {

std::string Lower(std::string s) {
	for (auto &c : s) {
		c = static_cast<char>(std::tolower(static_cast<unsigned char>(c)));
	}
	return s;
}

// Counts are doubles (a feature table may hold non-integer values), but a whole
// number must read as one: "7 sequences", not "7.000000".
std::string FormatCount(double v) {
	std::ostringstream os;
	os << std::setprecision(15) << v;
	return os.str();
}

[[noreturn]] void Fail(const std::string &message) {
	throw std::invalid_argument(message);
}

struct SampleMeta {
	std::optional<bool> is_source; // set by the sample's source_sink row
	std::string env;
};

} // namespace

Dataset IngestDataset(const std::vector<CooRow> &table, const std::vector<MetadataRow> &metadata) {
	// Metadata first: a role or environment problem is a property of the study
	// design and should be reported before any question about the counts.
	std::unordered_map<std::string, SampleMeta> meta;
	for (const auto &row : metadata) {
		const std::string variable = Lower(row.variable);
		auto &m = meta[row.sample_id];
		if (variable == "source_sink") {
			if (m.is_source.has_value()) {
				Fail("duplicate metadata rows for sample '" + row.sample_id +
				     "': source_sink is given more than once (is the sample listed twice?)");
			}
			const std::string role = Lower(row.value);
			if (role == "source") {
				m.is_source = true;
			} else if (role == "sink") {
				m.is_source = false;
			} else {
				Fail("sample '" + row.sample_id + "' has source_sink = '" + row.value +
				     "'; expected 'source' or 'sink'");
			}
		} else if (variable == "env") {
			m.env = row.value;
		}
	}

	Dataset ds;
	ds.sample_ids.reserve(meta.size());
	for (const auto &entry : meta) {
		ds.sample_ids.push_back(entry.first);
	}
	std::sort(ds.sample_ids.begin(), ds.sample_ids.end());
	if (ds.sample_ids.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
		Fail("too many samples for st3's 32-bit indices");
	}

	size_t n_sources = 0;
	ds.is_source.reserve(ds.sample_ids.size());
	ds.envs.reserve(ds.sample_ids.size());
	for (const auto &sample : ds.sample_ids) {
		const auto &m = meta.at(sample);
		if (!m.is_source.has_value()) {
			Fail("sample '" + sample + "' has no source_sink value in the metadata");
		}
		const bool is_source = *m.is_source;
		if (is_source && m.env.empty()) {
			Fail("source sample '" + sample + "' has an empty env; every source needs an environment");
		}
		ds.is_source.push_back(is_source);
		ds.envs.push_back(is_source ? std::optional<std::string> {m.env} : std::nullopt);
		n_sources += is_source ? 1 : 0;
	}
	if (n_sources == 0) {
		Fail("no source samples in the metadata; at least one sample must have source_sink = 'source'");
	}

	// Feature dictionary from the cells, then every cell to its indices.
	{
		std::unordered_set<std::string> seen;
		for (const auto &cell : table) {
			if (seen.insert(cell.feature_id).second) {
				ds.feature_ids.push_back(cell.feature_id);
			}
		}
	}
	std::sort(ds.feature_ids.begin(), ds.feature_ids.end());
	if (ds.feature_ids.size() > static_cast<size_t>(std::numeric_limits<int32_t>::max())) {
		Fail("too many features for st3's 32-bit indices");
	}
	std::unordered_map<std::string, int32_t> sample_index;
	for (size_t i = 0; i < ds.sample_ids.size(); ++i) {
		sample_index[ds.sample_ids[i]] = static_cast<int32_t>(i);
	}
	std::unordered_map<std::string, int32_t> feature_index;
	for (size_t i = 0; i < ds.feature_ids.size(); ++i) {
		feature_index[ds.feature_ids[i]] = static_cast<int32_t>(i);
	}

	ds.rows.reserve(table.size());
	ds.cols.reserve(table.size());
	ds.vals.reserve(table.size());
	ds.sample_totals.assign(ds.sample_ids.size(), 0.0);
	std::unordered_set<int64_t> cells_seen;
	cells_seen.reserve(table.size());
	const auto n_samples = static_cast<int64_t>(ds.sample_ids.size());
	for (const auto &cell : table) {
		const auto s = sample_index.find(cell.sample_id);
		if (s == sample_index.end()) {
			Fail("sample '" + cell.sample_id + "' is in the feature table but has no metadata row");
		}
		if (!std::isfinite(cell.count) || cell.count < 0) {
			Fail("feature table cell (sample '" + cell.sample_id + "', feature '" + cell.feature_id + "') has count " +
			     FormatCount(cell.count) + "; counts must be finite and >= 0");
		}
		// st3 stores integer counts and floors what it is given, so a table of
		// relative abundances would be silently truncated toward zero.
		// SourceTracker2 refuses non-integer tables too.
		if (cell.count != std::floor(cell.count)) {
			Fail("feature table cell (sample '" + cell.sample_id + "', feature '" + cell.feature_id + "') has count " +
			     FormatCount(cell.count) + "; counts must be whole numbers");
		}
		const int32_t fi = feature_index.at(cell.feature_id);
		const int32_t si = s->second;
		if (!cells_seen.insert(static_cast<int64_t>(fi) * n_samples + si).second) {
			Fail("duplicate feature table cell (sample '" + cell.sample_id + "', feature '" + cell.feature_id +
			     "'); each (sample_id, feature_id) may appear once");
		}
		ds.rows.push_back(fi);
		ds.cols.push_back(si);
		ds.vals.push_back(cell.count);
		ds.sample_totals[static_cast<size_t>(si)] += cell.count;
	}

	// Every metadata sample must have contributed at least one cell. The reader
	// drops zero cells, so a sample with none is either absent from the table
	// or all zeros -- an empty sink has nothing to attribute and an empty source
	// silently thins its environment, and neither is what the user meant.
	std::vector<bool> has_cell(ds.sample_ids.size(), false);
	for (const auto si : ds.cols) {
		has_cell[static_cast<size_t>(si)] = true;
	}
	for (size_t i = 0; i < ds.sample_ids.size(); ++i) {
		if (!has_cell[i]) {
			Fail("sample '" + ds.sample_ids[i] + "' is in the metadata but has no non-zero cells in the feature table");
		}
	}
	return ds;
}

void CheckDepths(const Dataset &dataset, int32_t source_depth, int32_t sink_depth, bool loo) {
	// Sink mode rarefies sinks per sample and sources per collapsed environment
	// (st3 checks the latter); leave-one-out rarefies source samples and has no
	// sinks. So exactly one role is checked here, against its own depth.
	const bool check_sources = loo;
	const int32_t depth = loo ? source_depth : sink_depth;
	const char *role = loo ? "source" : "sink";
	if (depth == 0) {
		return;
	}
	if (depth < 0) {
		Fail(std::string(role) + " rarefaction depth must be >= 0");
	}

	size_t n_shallow = 0;
	size_t shallowest = 0;
	double shallowest_total = std::numeric_limits<double>::infinity();
	for (size_t i = 0; i < dataset.sample_ids.size(); ++i) {
		if (dataset.is_source[i] != check_sources) {
			continue;
		}
		const double total = dataset.sample_totals[i];
		if (total < static_cast<double>(depth)) {
			++n_shallow;
			if (total < shallowest_total) {
				shallowest_total = total;
				shallowest = i;
			}
		}
	}
	if (n_shallow == 0) {
		return;
	}
	std::ostringstream os;
	os << "You requested rarefaction of " << role << " samples at " << depth << ", but " << n_shallow << " " << role
	   << (n_shallow == 1 ? " sample has" : " samples have")
	   << " fewer sequences than that. The shallowest of these is '" << dataset.sample_ids[shallowest] << "' with "
	   << FormatCount(shallowest_total) << " sequences.";
	Fail(os.str());
}

} // namespace miint::sourcetracker
