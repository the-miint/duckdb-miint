#pragma once

#include <cstdint>
#include <memory>
#include <string>
#include <string_view>
#include <unordered_map>
#include <vector>

#include "duckdb/common/arrow/arrow.hpp"

namespace miint {

//! One COO matrix as Arrow C Data Interface arrays: `(rows, cols, vals)` plus
//! the id dictionaries the indices refer to.
//!
//! Deliberately not any consumer's struct. A consumer that wants its own layout
//! copies these fields into it -- see `AsScTable` in sc_common.hpp, which is the
//! whole of the sc-specific part.
struct CooArrays {
	ArrowArray rows {}, cols {}, vals {}, sample_ids {}, feature_ids {};
	ArrowSchema rows_schema {}, cols_schema {}, vals_schema {}, sample_ids_schema {}, feature_ids_schema {};
	int64_t n_samples = 0;
	int64_t n_features = 0;
};

//! Owns every buffer behind a `CooArrays` and releases them on destruction.
//!
//! A consumer that borrows across the C Data Interface -- reading the buffers in
//! place and never invoking `release`, which is what sc's importer does -- needs
//! this object to outlive the call that reads it.
class CooTable {
public:
	CooTable() = default;
	~CooTable();
	CooTable(const CooTable &) = delete;
	CooTable &operator=(const CooTable &) = delete;

	//! Borrowed; valid while this object lives.
	const CooArrays &arrays() const {
		return arrays_;
	}
	int64_t NumSamples() const {
		return arrays_.n_samples;
	}
	int64_t NumFeatures() const {
		return arrays_.n_features;
	}
	//! Number of stored triples. Duplicates are NOT collapsed here — sc sums
	//! them (scipy COO semantics, sc-core `matrix.rs` `from_coo`).
	int64_t NumNonZeros() const {
		return arrays_.rows.length;
	}

	//! Per sample, the fraction of its observed cells whose feature the model
	//! knows: `matched / observed`, in `SampleIds()` order.
	//!
	//! The denominator is the *sample's* features, not the model's. Coverage
	//! against the model's vocabulary is always tiny in sparse data -- a sample
	//! legitimately carries a handful of 200k features -- so it would flag
	//! everything. This ratio instead answers "how much of what I observed here
	//! can the model actually use", where a low value means the model is being
	//! applied to different data: another reference database, another pipeline,
	//! another 16S region.
	//!
	//! Always 1.0 without a fixed vocabulary, where every feature is known by
	//! construction.
	const std::vector<double> &SampleCoverage() const {
		return sample_coverage_;
	}

	//! Dictionaries, in the sorted order the `cols` / `rows` indices refer to.
	const std::vector<std::string> &SampleIds() const {
		return sample_ids_;
	}
	const std::vector<std::string> &FeatureIds() const {
		return feature_ids_;
	}

private:
	friend class CooBuilder;
	friend class CooBatcher;
	CooArrays arrays_ {};
	std::vector<double> sample_coverage_;
	std::vector<std::string> sample_ids_;
	std::vector<std::string> feature_ids_;
};

//! Cuts a finalized table into standalone tables of consecutive samples.
//!
//! For a caller whose sc output grows with the samples in one call -- sc_shap's
//! dense attribution matrix -- and so must bound it. sc scores each sample
//! independently, so how samples are grouped changes no result.
//!
//! The cells are indexed by sample once, one `size_t` per cell; each batch then
//! costs only its own cells rather than a scan of the whole table. `table` must
//! outlive this object.
class CooBatcher {
public:
	explicit CooBatcher(const CooTable &table);

	//! Samples `[first, first + count)` as their own table: rows renumbered from
	//! 0, the same feature columns and vocabulary, coverage carried along. A
	//! batch whose samples have no cells is valid -- all-zero rows.
	std::unique_ptr<CooTable> Batch(size_t first, size_t count) const;

private:
	const CooTable &table_;
	//! Sample s's cells are cell_order_[sample_start_[s] .. sample_start_[s + 1]].
	std::vector<size_t> sample_start_;
	std::vector<size_t> cell_order_;
};

//! One `(sample_id, feature_id)` pair that was appended more than once, with the
//! values that were seen for it.
//!
//! The values matter as much as the keys: identical values point at a join
//! fanout (a duplicated key upstream, where one row is spurious), differing
//! values at genuine repeat measurements. The correct repair is opposite in the
//! two cases -- deduplicating versus summing -- and summing a fanout silently
//! inflates every affected count, so a caller reporting this should show them.
struct DuplicateCell {
	std::string sample_id;
	std::string feature_id;
	size_t count = 0;
	std::vector<double> values;
};

//! What [`CooBuilder::FindDuplicateCells`] found.
struct DuplicateReport {
	//! Distinct `(sample_id, feature_id)` pairs appearing more than once.
	size_t duplicate_cells = 0;
	//! A bounded sample of them, for an error message.
	std::vector<DuplicateCell> examples;

	bool Empty() const {
		return duplicate_cells == 0;
	}
};

//! Accumulates long-format `(sample_id, feature_id, value)` cells and emits the
//! COO feature table sc expects.
//!
//! This is the bridge between the shape DuckDB has and the shape sc wants.
//! `read_biom`, `woltka_ogu_per_sample`, and anything else projected to those
//! three columns all arrive here — sc is coupled to the triple, not to a file
//! format.
//!
//! **Both dictionaries are sorted**, so a given set of ids always produces the
//! same integer encoding regardless of the order rows arrive in. That matters
//! beyond tidiness: a model trained when `OTU_7` was column 6 returns nonsense
//! if prediction data encodes it as column 2, and sc only validates the feature
//! *count* (sc-core `model.rs` `check_features`), never the identity. It also
//! matches the reference pipeline sc's own fixtures were generated with
//! (`oracle/generator/gen_forest.py` `load_aligned`, which sorts both axes).
class CooBuilder {
public:
	//! Ids are copied on first sight and interned thereafter.
	//!
	//! With a fixed vocabulary (see [`SetFeatureVocabulary`]) a feature the
	//! model never saw is dropped and counted. The sample is interned either
	//! way, so a sample whose every feature was dropped still gets a row -- an
	//! all-zero one, which is the truthful representation of "none of the
	//! model's features were observed here" and is what lets it still receive a
	//! prediction rather than silently vanishing from the output.
	void Append(std::string_view sample_id, std::string_view feature_id, double value);

	//! Encode features against a model's training vocabulary instead of
	//! deriving one from the data.
	//!
	//! This is what makes a model reusable across datasets. Column `i` means
	//! `vocab[i]` because that is what the model was trained on -- so the
	//! vocabulary is used in the model's own order and is NOT re-sorted. Two
	//! datasets that differ by one feature in each direction have the same
	//! width, and sc validates only the width (`RandomForest::check_features`),
	//! so re-deriving an encoding here would produce confident, silent
	//! nonsense.
	//!
	//! Consequences, all of which are correct for a sparse matrix:
	//!   * a feature not in `vocab` is dropped -- the model has no column for it
	//!   * a feature in `vocab` but absent from the data is simply not stored,
	//!     and absent already means zero
	//!   * `n_features` is always `vocab.size()`, never what the data happened
	//!     to contain
	//!
	//! Must be called before the first [`Append`].
	void SetFeatureVocabulary(std::vector<std::string> vocab);

	//! A few distinct feature ids that were dropped, for diagnosing a mismatch.
	//!
	//! Bounded and only collected when a vocabulary is fixed, so this costs
	//! nothing on the fit path and at most a handful of small strings on the
	//! predict path -- enough to tell "different database" from "same ids,
	//! different spelling".
	const std::vector<std::string> &DroppedExamples() const {
		return dropped_examples_;
	}

	//! Cells dropped because their feature is not in the fixed vocabulary.
	//!
	//! Worth surfacing: a prediction table sharing almost no features with the
	//! model still yields confident predictions from a near-empty matrix, and
	//! nothing else in the pipeline will say so.
	size_t DroppedCells() const {
		return dropped_cells_;
	}

	//! Report `(sample_id, feature_id)` pairs appended more than once.
	//!
	//! Policy-free by design: sc's `from_coo` *sums* duplicates (scipy COO
	//! semantics), which is right for genuine repeat measurements and wrong for
	//! a join fanout. The builder cannot tell those apart, so it reports and
	//! lets the caller decide.
	//!
	//! Costs one temporary `uint64` per stored cell (8 bytes, freed on return)
	//! plus an O(n log n) sort -- deliberately not a hash set, whose per-node
	//! overhead would be several times larger. That matters for the wasm build,
	//! where DuckDB and this extension share one linear memory.
	//!
	//! Safe to call before [`Finalize`]; it reads the interned indices and does
	//! not modify them.
	DuplicateReport FindDuplicateCells(size_t max_examples = 5) const;

	//! Sorts the dictionaries, remaps the interned indices, and builds the
	//! Arrow arrays. The builder is left empty and reusable.
	//!
	//! With a fixed vocabulary only the sample dictionary is sorted; the
	//! feature columns already mean what the model says they mean.
	//!
	//! Returns nullptr if nothing was appended: sc rejects a 0-row or 0-column
	//! matrix outright, so an empty table has no valid representation.
	std::unique_ptr<CooTable> Finalize();

	size_t NumNonZeros() const {
		return vals_.size();
	}
	size_t NumSamples() const {
		return sample_index_.size();
	}
	size_t NumFeatures() const {
		return feature_index_.size();
	}

private:
	//! Intern `id`, returning its provisional (insertion-order) index.
	static int64_t Intern(std::unordered_map<std::string, int64_t> &index, std::vector<std::string> &ids,
	                      std::string_view id);

	bool has_fixed_features_ = false;
	size_t dropped_cells_ = 0;
	//! At most kMaxDroppedExamples distinct ids; see DroppedExamples().
	std::vector<std::string> dropped_examples_;

	std::unordered_map<std::string, int64_t> sample_index_;
	std::unordered_map<std::string, int64_t> feature_index_;
	std::vector<std::string> sample_ids_;
	std::vector<std::string> feature_ids_;
	//! Indexed by provisional sample index; remapped alongside the dictionary.
	std::vector<int64_t> sample_cells_;
	std::vector<int64_t> sample_matched_;
	std::vector<int64_t> rows_;
	std::vector<int64_t> cols_;
	std::vector<double> vals_;
};

} // namespace miint
