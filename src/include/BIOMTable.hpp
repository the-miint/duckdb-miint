#include <vector>
#include <hdf5.h>
#include <string>

namespace miint {

enum class BIOMTableField { SAMPLE_ID = 0, FEATURE_ID, VALUE };

// Sparse matrix representation for CSR/CSC formats
struct SparseMatrix {
	std::vector<double> data;
	std::vector<int32_t> indices;
	std::vector<int32_t> indptr;

	SparseMatrix() = default;
	SparseMatrix(std::vector<double> data_, std::vector<int32_t> indices_, std::vector<int32_t> indptr_)
	    : data(std::move(data_)), indices(std::move(indices_)), indptr(std::move(indptr_)) {
	}
};

class BIOMTable {
public:
	BIOMTable(hid_t indices, hid_t indptr, hid_t data, hid_t obs_ids, hid_t samp_ids);
	BIOMTable(const std::vector<std::string> &feature_ids, const std::vector<std::string> &sample_ids,
	          const std::vector<double> &values);
	// Constructor accepting pre-computed integer indices and ordered ID lists
	// Optimized for performance - skips string hashing that happens in string-based constructor
	BIOMTable(const std::vector<size_t> &feature_indices, const std::vector<size_t> &sample_indices,
	          const std::vector<double> &values, std::vector<std::string> feature_ids_ordered,
	          std::vector<std::string> sample_ids_ordered);
	BIOMTable();
	uint32_t nnz() const;
	const std::vector<size_t> &COOFeatureIndices() const;
	const std::vector<size_t> &COOSampleIndices() const;
	const std::vector<double> &COOValues() const;

	// Get ordered feature and sample IDs
	const std::vector<std::string> &FeatureIDs() const {
		return feature_ids_ordered;
	}
	const std::vector<std::string> &SampleIDs() const {
		return sample_ids_ordered;
	}

	// Get matrix dimensions
	size_t NumFeatures() const {
		return feature_ids_ordered.size();
	}
	size_t NumSamples() const {
		return sample_ids_ordered.size();
	}

	// Convert COO to CSR (Compressed Sparse Row) format
	// Used for observation/matrix in BIOM (features = rows, samples = columns)
	SparseMatrix ToCSR() const;

	// Convert COO to CSC (Compressed Sparse Column) format
	// Used for sample/matrix in BIOM (features = rows, samples = columns)
	SparseMatrix ToCSC() const;

private:
	std::vector<size_t> coo_feature_indices;
	std::vector<size_t> coo_sample_indices;
	std::vector<double> coo_values;

	std::vector<std::string> feature_ids_ordered;
	std::vector<std::string> sample_ids_ordered;

	void InitCOOFromCSC(const std::vector<int32_t> &indptr, const std::vector<int32_t> &indices);
	void InitCOOFromCOO(const std::vector<std::string> &feature_ids, const std::vector<std::string> &sample_ids);

	// accumulate duplicate sample / feature pairs
	// remove zeros
	void compress_coo();
	std::vector<std::string> indices_to_ids(const std::vector<size_t> &indices, const std::vector<std::string> &names);

	// Helper for ToCSR/ToCSC: convert COO to compressed format
	// by_row=true: CSR, by_row=false: CSC
	SparseMatrix ConvertCOOToCompressed(bool by_row) const;

	std::vector<std::string> load_dataset_1D_str(hid_t ds_id);

	template <typename T>
	std::vector<T> load_dataset_1D(hid_t ds_id, hid_t exp_dtype);

	template <typename T>
	hid_t get_hdf5_type();
};

std::vector<std::string> unique_ids_in_order(const std::vector<std::string> &ids);
std::vector<size_t> ids_to_indices(const std::vector<std::string> &ids, const std::vector<std::string> &ordered);
void sort_by_row_then_column(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values);
void apply_permutation(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values,
                       const std::vector<size_t> &indices);
} // namespace miint
