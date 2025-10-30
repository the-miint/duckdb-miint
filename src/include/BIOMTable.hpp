#include <vector>
#include "H5Cpp.h"

namespace miint {

enum class BIOMTableField { SAMPLE_ID = 0, FEATURE_ID, VALUE };

class BIOMTable {
public:
	BIOMTable(const H5::DataSet &indices, const H5::DataSet &indptr, const H5::DataSet &data,
	          const H5::DataSet &obs_ids, const H5::DataSet &samp_ids);
	BIOMTable(const std::vector<std::string> &feature_ids, const std::vector<std::string> &sample_ids,
	          const std::vector<double> &values);
	BIOMTable();
	uint32_t nnz(); // not const as it implicitly calls compress_coo
	const std::vector<std::string> &COOFeatures() const;
	const std::vector<size_t> &COOFeatureIndices() const;
	const std::vector<std::string> &COOSamples() const;
	const std::vector<size_t> &COOSampleIndices() const;
	const std::vector<double> &COOValues() const;

private:
	std::vector<size_t> coo_feature_indices;
	std::vector<size_t> coo_sample_indices;
	std::vector<double> coo_values;

	std::vector<std::string> feature_ids_ordered;
	std::vector<std::string> sample_ids_ordered;
	std::vector<std::string> coo_sample_indices_as_ids;
	std::vector<std::string> coo_feature_indices_as_ids;

	void InitCOOFromCSC(const std::vector<int32_t> &indptr, const std::vector<int32_t> &indices);
	void InitCOOFromCOO(const std::vector<std::string> &feature_ids, const std::vector<std::string> &sample_ids);

	// accumulate duplicate sample / feature pairs
	// remove zeros
	void compress_coo();
	std::vector<std::string> indices_to_ids(const std::vector<size_t> &indices, const std::vector<std::string> &names);

	std::vector<std::string> load_dataset_1D_str(const H5::DataSet &ds_ids);

	template <typename T>
	std::vector<T> load_dataset_1D(const H5::DataSet &ds, const H5::PredType &exp_dtype);

	template <typename T>
	H5::DataType get_hdf5_type();
};

std::vector<std::string> unique_ids_in_order(const std::vector<std::string> &ids);
std::vector<size_t> ids_to_indices(const std::vector<std::string> &ids, const std::vector<std::string> &ordered);
void sort_by_row_then_column(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values);
void apply_permutation(std::vector<size_t> &rows, std::vector<size_t> &cols, std::vector<double> &values,
                       const std::vector<size_t> &indices);
} // namespace miint
