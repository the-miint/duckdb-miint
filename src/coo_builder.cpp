#include "coo_builder.hpp"

#include "arrow_export.hpp"

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <numeric>
#include <stdexcept>

namespace miint {

// helper functions go in here so that nobody else using miint::<name> can see them
namespace {

//! How many distinct dropped feature ids to remember for diagnostics.
constexpr size_t kMaxDroppedExamples = 5;

} // namespace

// c++ destructor for CooTable, releases all the arrow arrays and schemas if they are live when CooTable goes out of
// scope and is removed from memory. This is important because the arrow arrays and schemas are allocated on the heap,
// and if they are not released, they will cause a memory leak. The destructor is called automatically when the CooTable
// object goes out of scope, so we don't have to worry about manually calling it. (as opposed to C++ constructor
// function or allocating an instance in memory)
CooTable::~CooTable() {
	ReleaseIfLive(arrays_.rows, arrays_.rows_schema);
	ReleaseIfLive(arrays_.cols, arrays_.cols_schema);
	ReleaseIfLive(arrays_.vals, arrays_.vals_schema);
	ReleaseIfLive(arrays_.sample_ids, arrays_.sample_ids_schema);
	ReleaseIfLive(arrays_.feature_ids, arrays_.feature_ids_schema);
}

/*
    params
    std::unordered_map<std::string, int64_t> &index: a reference to an unordered map that maps strings to integers. This
   is used to store the mapping of sample/feature ids to their corresponding indices.

    std::vector<std::string> &ids: a reference to a vector of strings that stores the unique sample/feature ids in the
   order they were first seen. This is used to maintain the order of the ids for later sorting.

    std::string_view id: a string view that represents the sample/feature id to be interned. This is the id that we want
   to encode into a unique integer index.

    dynamic unique number mapping to string id, if the string is already in the index, return the existing index,
   otherwise add it to the index and return the new index. This is used to encode the sample and feature ids into
   numeric indices for the COO matrix.

    Named Intern to reflect string interning.  The function is like a coat check. Give coat -> get numerical tag.
*/
int64_t CooBuilder::Intern(std::unordered_map<std::string, int64_t> &index, std::vector<std::string> &ids,
                           std::string_view id) {
	auto it = index.find(std::string(id));
	if (it != index.end()) {
		return it->second; // if found return the unique id int64
	}
	// if not found, add it to the index and return the new unique id int64
	const auto next = static_cast<int64_t>(ids.size());
	ids.emplace_back(id);
	index.emplace(ids.back(), next);
	return next;
}

// This is where we process each record triple (sample_id, feature_id, value) and store them in the COO format. We use
// the Intern function to get the unique indices for the sample and feature ids, and then we store the row index, column
// index, and value in their respective vectors. This allows us to build the COO matrix incrementally as we process each
// record.
void CooBuilder::SetFeatureVocabulary(std::vector<std::string> vocab) {
	feature_ids_ = std::move(vocab);
	feature_index_.clear();
	for (size_t i = 0; i < feature_ids_.size(); i++) {
		// First occurrence wins, so the index always points at the model's own
		// column for that feature.
		feature_index_.emplace(feature_ids_[i], static_cast<int64_t>(i));
	}
	has_fixed_features_ = true;
	dropped_cells_ = 0;
}

// intake a triplet (sample_id, feature_id, value) and store it in the COO format. We use the Intern function to get the
// unique indices for the sample and feature ids, and then we store the row index, column index, and value in their
// respective vectors. If a fixed vocabulary is set, we check if the feature_id is in the vocabulary, and if not, we
// drop the cell and increment the dropped_cells_ counter. This allows us to build the COO matrix incrementally as we
// process each record.
void CooBuilder::Append(std::string_view sample_id, std::string_view feature_id, double value) {
	// ------------------ samples ------------------
	// Intern the sample first, unconditionally. A sample every one of whose
	// features is unknown to the model still belongs in the output: it gets an
	// all-zero row and a prediction, rather than silently vanishing.
	const auto row = Intern(sample_index_, sample_ids_, sample_id);
	// Intern hands out sequential indices, so the per-sample counters only ever
	// need to grow by one.
	if (static_cast<size_t>(row) == sample_cells_.size()) {
		sample_cells_.push_back(0);
		sample_matched_.push_back(0);
	}
	sample_cells_[static_cast<size_t>(row)]++;

	// ------------------ features ------------------
	int64_t col;
	if (has_fixed_features_) {
		// O(1) hash lookup, not a scan or a binary search. The model's
		// vocabulary happens to be sorted today (our fit sorts it), but binary
		// search would bake that in -- a model fit through sc's C API directly,
		// or a future RFE-reordered bundle, would then mis-resolve silently.
		const auto it = feature_index_.find(std::string(feature_id));
		if (it == feature_index_.end()) {
			dropped_cells_++; // drop it, its not in the model's vocabulary for prediction, and count it for reporting
			// Keep the first few distinct ones. A linear scan over <= 5 entries
			// beats a set, and the cap keeps a wholly-foreign table from
			// accumulating one string per feature.
			if (dropped_examples_.size() < kMaxDroppedExamples &&
			    std::find(dropped_examples_.begin(), dropped_examples_.end(), feature_id) == dropped_examples_.end()) {
				dropped_examples_.emplace_back(feature_id);
			}
			return;
		}
		col = it->second;
	} else {
		col = Intern(feature_index_, feature_ids_, feature_id);
	}
	sample_matched_[static_cast<size_t>(row)]++;

	rows_.push_back(row);
	cols_.push_back(col);
	vals_.push_back(value);
}

namespace {

// index translation pipeline for canonical ordering of the COO matrix. This is where we sort the sample and feature
// ids, and remap the row and column indices to match the sorted order. This ensures that the COO matrix is in a
// consistent order regardless of the order in which the records were appended (db scanned). The SortDictionary function
// is used to sort the ids and produce a remapping of the indices.
/*
    Intake a vector of strings std::vector<std::string>& ids, argsort them to get ordered set of vocab for this
   dimension.

    Original IDs scanned:    ['Zebra', 'Apple', 'Mango']
    Original Rows (int64 aranged):   [0, 1, 0, 2] -> ['Zebra', 'Apple', 'Zebra', 'Mango']

    --- ARGSORT --- Deterministic ordering for the dimension, so that the COO matrix is always in the same order
   regardless of the order in which the records were appended (db scanned). order:           [1, 2, 0] -> ['Apple',
   'Mango', 'Zebra']  # indices of the original ids that would sort them

    --- INVERSION map ---
    remap:           [2, 0, 1] - map original aranged int64_t to new sorted canonical indices. Deterministic ordering
   for each COO dimension. remap[original[0]] = 2 = 'Zebra' remap[original[1]] = 0 = 'Apple' remap[original[2]] = 1 =
   'Mango'


*/
std::vector<int64_t> SortDictionary(std::vector<std::string> &ids) {
	// allocate a vector of int64_t with the same size as ids, and fill it with the values 0, 1, 2, ..., ids.size() - 1.
	// This will be used to keep track of the original indices of the ids before sorting.
	std::vector<int64_t> order(ids.size());
	std::iota(order.begin(), order.end(), 0);

	// order is an int64 list argsorted on input strings byte for byte
	// provides a deterministic ordering of the input dim where the elements are indices of the original ids vector
	// indicating a sorted order
	std::sort(order.begin(), order.end(), [&ids](int64_t a, int64_t b) {
		return ids[static_cast<size_t>(a)] < ids[static_cast<size_t>(b)];
	}); // whenever you see size_t type think indexing something

	// order[new] = old, so invert it into remap[old] = new.
	std::vector<int64_t> remap(ids.size());
	// std::vector is a 24 byte struct containing 3 pointers: pointer to the data, size, and capacity.
	std::vector<std::string> sorted;
	sorted.reserve(ids.size());
	for (size_t newpos = 0; newpos < order.size(); newpos++) {
		const auto oldpos = static_cast<size_t>(order[newpos]);
		remap[oldpos] = static_cast<int64_t>(newpos);
		sorted.push_back(std::move(ids[oldpos]));
	}
	// std::move is zero copy pointer swaps
	// 1) deallocates old buffer ids
	// 2) steals pointers: ids copies the three pointers from sorted
	// 3) nulls out sorted's pointers so it doesn't free the buffer when it goes out of scope
	// the ids struct is not equivalent to the sorted struct. Sorted struct is nullified so buffer is not freed as soon
	// as this function goes out of scope which happens in the next few lines. this way ids does not become a dangling
	// pointer and we get segfault when we try to access it later.
	ids = std::move(sorted);
	return remap;
}

} // namespace

/*
    First CooBuilder is initialized and BIOM triples are appended to the <string, int> mapping

    Then CooBuilder::Finalize() is called.
*/
namespace {
/*
bitwise packing to reduce memory overhead instead of using a hash table

Step 1: static_cast<uint64_t>(row=2)  -> 64-bit container:
[ 0000 0000 ... 0000 0000 ] [ 0000 0000 ... 0000 0010 ]
   Upper 32 bits (32..63)       Lower 32 bits (0..31)

Step 2: << 32  (Shift into the upper half):
[ 0000 0000 ... 0000 0010 ] [ 0000 0000 ... 0000 0000 ]
   Upper 32 bits (row=2)        Lower 32 bits (empty zeros)

Step 3: static_cast<uint32_t>(col=1)  (Lives in lower half):
[ 0000 0000 ... 0000 0000 ] [ 0000 0000 ... 0000 0001 ]

Step 4: Combine with | :
[ 0000 0000 ... 0000 0010 ] [ 0000 0000 ... 0000 0001 ]
   Upper 32 bits = 2            Lower 32 bits = 1
*/
//! Pack a (row, col) index pair into one sortable key. Both are dictionary
//! positions, so they are bounded by the number of distinct samples / features
//! and fit in 32 bits for any table that could be built in memory.
inline uint64_t PackCell(int64_t row, int64_t col) {
	return (static_cast<uint64_t>(row) << 32) | static_cast<uint32_t>(col);
}

} // namespace
/*
 * Bitwise packs (row, col) into a uint64_t key to detect duplicate cells via
 * sorting (O(N log N)) rather than a hash set (O(1)).
 *
 * A contiguous std::vector avoids the node-allocation overhead (~40+ bytes/pair)
 * and cache misses of std::unordered_set, using exactly 8 bytes per cell.
 * Minimizing peak memory overhead is critical in WebAssembly environments,
 * where linear memory is constrained and shared between DuckDB and this extension.
 */
DuplicateReport CooBuilder::FindDuplicateCells(size_t max_examples) const {
	DuplicateReport report;
	if (rows_.size() < 2) {
		return report;
	}

	// 8 bytes per cell, freed on return. A hash set of pairs would cost several
	// times this in node overhead alone.
	std::vector<uint64_t> keys;
	keys.reserve(rows_.size());
	for (size_t i = 0; i < rows_.size(); i++) {
		keys.push_back(PackCell(rows_[i], cols_[i])); // register a row and col pair
	}
	std::sort(keys.begin(), keys.end());

	// Duplicates are adjacent once sorted, so walk the runs of equal keys. Count
	// every run longer than one; keep only the first few keys for the report.
	std::vector<uint64_t> wanted;
	for (size_t i = 0; i < keys.size();) {
		size_t j = i + 1;
		while (j < keys.size() && keys[j] == keys[i]) {
			j++;
		}
		if (j - i > 1) {
			report.duplicate_cells++;
			if (wanted.size() < max_examples) {
				wanted.push_back(keys[i]);
			}
		}
		i = j;
	}
	if (wanted.empty()) {
		return report;
	}

	// Second pass over the triples, collecting values for the sampled keys only,
	// so the report can show whether they are identical (a join fanout) or
	// different (genuine repeat measurements).
	std::vector<DuplicateCell> cells(wanted.size());
	for (size_t i = 0; i < rows_.size(); i++) {
		const auto key = PackCell(rows_[i], cols_[i]);
		for (size_t w = 0; w < wanted.size(); w++) {
			if (wanted[w] != key) {
				continue;
			}
			auto &cell = cells[w];
			if (cell.count == 0) {
				cell.sample_id = sample_ids_[static_cast<size_t>(rows_[i])];
				cell.feature_id = feature_ids_[static_cast<size_t>(cols_[i])];
			}
			cell.count++;
			cell.values.push_back(vals_[i]);
			break;
		}
	}
	report.examples = std::move(cells);
	return report;
}

std::unique_ptr<CooTable> CooBuilder::Finalize() {
	// sc rejects a 0 x N or N x 0 matrix (`from_coo`: "dimensions must be > 0"),
	// so an empty input has no representable table. Say so here rather than
	// building one sc will refuse.
	if (sample_ids_.empty() || feature_ids_.empty()) {
		return nullptr;
	}

	// sort and produce remappings
	const auto sample_remap = SortDictionary(sample_ids_);

	// Canonicalize row (sample) IDs to guarantee deterministic training & CV.
	// Parallel DuckDB scans deliver chunks in arbitrary order. Because downstream
	// RNG operations (bootstrap sampling, chunk-based CV fold assignment) sample by
	// positional index, an unstable row order causes identical random seeds to produce
	// different models and scores.
	// take old id and map to new id sorted along row strings
	for (auto &r : rows_) {
		r = sample_remap[static_cast<size_t>(r)];
	}
	// canonicalize features - extremely important. Data from duckdb can stream in any order. If we assigned feature ids
	// in the order they were seen, then the model would be trained on one permutation of a set of features and then we
	// would try to predict on another permutation of the a set of features, the feature ids would be different and the
	// model would be invalid. So we need to sort the feature ids and remap the feature ids to the sorted order so that
	// the model is trained on a consistent set of features.
	// this provides a deterministic ordering
	// [a,b,c], [c,b,a], [b,c,a] all map to [a,b,c] and the model is trained on the same feature ids regardless of the
	// order they were seen in the input data.
	// take old id and map to new id sorted along col strings
	// ...unless the vocabulary was fixed to a model's. Then the column order IS
	// the model's definition of what each column means, and re-sorting it here
	// would silently re-point every learned split at a different feature.
	if (!has_fixed_features_) {
		const auto feature_remap = SortDictionary(feature_ids_);
		for (auto &c : cols_) {
			c = feature_remap[static_cast<size_t>(c)];
		}
	}

	// exception safety, atomic creation and wrapping / automatic destruction of the CooTable object even though its on
	// the heap.  If any of the Export* functions throw an exception, the partially constructed CooTable will be
	// destroyed and its destructor will release any allocated memory.
	// Coverage rides along the same remap as the dictionary, so it stays aligned
	// with SampleIds().
	std::vector<double> coverage(sample_ids_.size(), 1.0);
	for (size_t old_pos = 0; old_pos < sample_cells_.size(); old_pos++) {
		const auto seen = sample_cells_[old_pos];
		const auto matched = sample_matched_[old_pos];
		const auto at = static_cast<size_t>(sample_remap[old_pos]);
		coverage[at] = seen > 0 ? static_cast<double>(matched) / static_cast<double>(seen) : 0.0;
	}

	auto out = std::make_unique<CooTable>();
	out->sample_coverage_ = std::move(coverage);
	// n_something is always some kind of boundary in cpp or a length
	out->arrays_.n_samples = static_cast<int64_t>(sample_ids_.size());
	out->arrays_.n_features = static_cast<int64_t>(feature_ids_.size());
	ExportInt64(out->arrays_.rows, out->arrays_.rows_schema, std::move(rows_));
	ExportInt64(out->arrays_.cols, out->arrays_.cols_schema, std::move(cols_));
	ExportFloat64(out->arrays_.vals, out->arrays_.vals_schema, std::move(vals_));
	ExportUtf8(out->arrays_.sample_ids, out->arrays_.sample_ids_schema, sample_ids_);
	ExportUtf8(out->arrays_.feature_ids, out->arrays_.feature_ids_schema, feature_ids_);
	out->sample_ids_ = std::move(sample_ids_);
	out->feature_ids_ = std::move(feature_ids_);

	sample_index_.clear();
	feature_index_.clear();
	has_fixed_features_ = false;
	dropped_cells_ = 0;
	dropped_examples_.clear();
	sample_cells_.clear();
	sample_matched_.clear();
	sample_ids_.clear();
	feature_ids_.clear();
	rows_.clear();
	cols_.clear();
	vals_.clear();
	return out;
}

CooBatcher::CooBatcher(const CooTable &table) : table_(table) {
	const auto n_samples = static_cast<size_t>(table.arrays_.n_samples);
	const auto nnz = static_cast<size_t>(table.arrays_.rows.length);
	const auto *rows = static_cast<const int64_t *>(table.arrays_.rows.buffers[1]);

	// Counting sort of cell indices by sample: count, prefix-sum, place. Stable,
	// so a sample's cells keep the order they had in the full table.
	sample_start_.assign(n_samples + 1, 0);
	for (size_t i = 0; i < nnz; i++) {
		sample_start_[static_cast<size_t>(rows[i]) + 1]++;
	}
	for (size_t s = 0; s < n_samples; s++) {
		sample_start_[s + 1] += sample_start_[s];
	}
	cell_order_.resize(nnz);
	std::vector<size_t> next(sample_start_.begin(), sample_start_.end() - 1);
	for (size_t i = 0; i < nnz; i++) {
		cell_order_[next[static_cast<size_t>(rows[i])]++] = i;
	}
}

std::unique_ptr<CooTable> CooBatcher::Batch(size_t first, size_t count) const {
	const size_t n_samples = sample_start_.size() - 1;
	if (count == 0 || first > n_samples || count > n_samples - first) {
		throw std::out_of_range("CooBatcher::Batch: samples [" + std::to_string(first) + ", " +
		                        std::to_string(first + count) + ") are outside a table of " +
		                        std::to_string(n_samples));
	}
	const auto &src = table_.arrays_;
	const auto *src_rows = static_cast<const int64_t *>(src.rows.buffers[1]);
	const auto *src_cols = static_cast<const int64_t *>(src.cols.buffers[1]);
	const auto *src_vals = static_cast<const double *>(src.vals.buffers[1]);

	const size_t begin = sample_start_[first];
	const size_t end = sample_start_[first + count];
	std::vector<int64_t> rows, cols;
	std::vector<double> vals;
	rows.reserve(end - begin);
	cols.reserve(end - begin);
	vals.reserve(end - begin);
	for (size_t k = begin; k < end; k++) {
		const size_t i = cell_order_[k];
		rows.push_back(src_rows[i] - static_cast<int64_t>(first));
		cols.push_back(src_cols[i]);
		vals.push_back(src_vals[i]);
	}

	auto out = std::make_unique<CooTable>();
	out->sample_ids_.assign(table_.sample_ids_.begin() + first, table_.sample_ids_.begin() + first + count);
	out->sample_coverage_.assign(table_.sample_coverage_.begin() + first,
	                             table_.sample_coverage_.begin() + first + count);
	out->feature_ids_ = table_.feature_ids_;
	out->arrays_.n_samples = static_cast<int64_t>(count);
	out->arrays_.n_features = src.n_features;
	ExportInt64(out->arrays_.rows, out->arrays_.rows_schema, std::move(rows));
	ExportInt64(out->arrays_.cols, out->arrays_.cols_schema, std::move(cols));
	ExportFloat64(out->arrays_.vals, out->arrays_.vals_schema, std::move(vals));
	ExportUtf8(out->arrays_.sample_ids, out->arrays_.sample_ids_schema, out->sample_ids_);
	ExportUtf8(out->arrays_.feature_ids, out->arrays_.feature_ids_schema, out->feature_ids_);
	return out;
}

} // namespace miint
