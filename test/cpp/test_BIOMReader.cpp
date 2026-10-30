#define CATCH_CONFIG_MAIN
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <chrono>
#include <filesystem>
#include <string>
#include <thread>
#include <vector>
#include "BIOMReader.hpp"
#ifndef _WIN32
#include <unistd.h>
#endif

TEST_CASE("BIOM empty table", "[BIOMReader]") {
	auto path = "data/biom/empty.biom";
	miint::BIOMReader reader(path);
	auto table = reader.read();

	REQUIRE((table.COOFeatureIndices().size() == 0));
	REQUIRE((table.COOSampleIndices().size() == 0));
	REQUIRE((table.COOValues().size() == 0));
	REQUIRE((table.nnz() == 0));
}

TEST_CASE("BIOM basic table", "[BIOMReader]") {
	auto path = "data/biom/test.biom";
	miint::BIOMReader reader(path);
	auto table = reader.read();

	REQUIRE((table.COOFeatureIndices().size() == 15));
	REQUIRE((table.COOSampleIndices().size() == 15));
	REQUIRE((table.COOValues().size() == 15));

	REQUIRE((table.nnz() == 15));

	std::vector<std::string> exp_sample_ids = {"Sample1", "Sample1", "Sample2", "Sample2", "Sample2",
	                                           "Sample3", "Sample3", "Sample3", "Sample3", "Sample4",
	                                           "Sample4", "Sample5", "Sample6", "Sample6", "Sample6"};

	std::vector<std::string> exp_feature_ids = {"GG_OTU_2", "GG_OTU_4", "GG_OTU_2", "GG_OTU_4", "GG_OTU_5",
	                                            "GG_OTU_1", "GG_OTU_3", "GG_OTU_4", "GG_OTU_5", "GG_OTU_2",
	                                            "GG_OTU_3", "GG_OTU_2", "GG_OTU_2", "GG_OTU_3", "GG_OTU_4"};

	std::vector<double> exp_values = {5.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 2.0, 4.0, 3.0, 1.0, 2.0, 1.0};

	// Convert indices to IDs for comparison
	const auto &sample_indices = table.COOSampleIndices();
	const auto &feature_indices = table.COOFeatureIndices();
	const auto &sample_ids_ordered = table.SampleIDs();
	const auto &feature_ids_ordered = table.FeatureIDs();

	std::vector<std::string> actual_sample_ids;
	std::vector<std::string> actual_feature_ids;
	for (size_t i = 0; i < sample_indices.size(); i++) {
		actual_sample_ids.push_back(sample_ids_ordered[sample_indices[i]]);
		actual_feature_ids.push_back(feature_ids_ordered[feature_indices[i]]);
	}

	REQUIRE((actual_sample_ids == exp_sample_ids));
	REQUIRE((actual_feature_ids == exp_feature_ids));
	REQUIRE((table.COOValues() == exp_values));
}

TEST_CASE("BIOM IsBIOM", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error printing
	auto path1 = "data/biom/test.biom";
	auto path2 = "data/biom/empty.biom";
	auto path3 = "data/sam/foo_no_header.sam";
	auto path4 = "data/biom/doesntexist";
	auto path5 = "data/biom/notbiom.h5";
	REQUIRE(miint::BIOMReader::IsBIOM(path1));
	REQUIRE(miint::BIOMReader::IsBIOM(path2));
	REQUIRE(!miint::BIOMReader::IsBIOM(path3));
	REQUIRE(!miint::BIOMReader::IsBIOM(path4));
}

TEST_CASE("BIOM file error recovery", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error printing
	// Try to open invalid file (should throw)
	REQUIRE_THROWS_AS(miint::BIOMReader("data/biom/notbiom.h5"), std::runtime_error);

	// Verify we can still open a valid file after the error
	REQUIRE_NOTHROW(miint::BIOMReader("data/biom/test.biom"));

	// And use it
	miint::BIOMReader reader("data/biom/test.biom");
	auto table = reader.read();
	REQUIRE((table.nnz() == 15));
}

TEST_CASE("BIOM IsBIOM refuses a file HDF5 cannot open, and the recorder keeps why", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr); // Disable HDF5 error printing, as the extension does at load
	SECTION("truncated") {
		miint::Hdf5FailureRecorder failure;
		REQUIRE_THROWS_AS(miint::BIOMReader::IsBIOM("data/biom/truncated.biom"), std::runtime_error);
		REQUIRE_THAT(failure.Cause(), Catch::Matchers::StartsWith("truncated file"));
	}
#ifndef _WIN32
	SECTION("unreadable") {
		if (geteuid() == 0) {
			SKIP("root reads a file whatever its mode");
		}
		// A copy of test.biom no one may read, given back its mode and removed however the section ends
		struct UnreadableCopy {
			std::filesystem::path path = std::filesystem::temp_directory_path() /
			                             ("miint_isbiom_unreadable_" + std::to_string(getpid()) + ".biom");
			UnreadableCopy() {
				std::filesystem::copy_file("data/biom/test.biom", path);
				std::filesystem::permissions(path, std::filesystem::perms::none);
			}
			~UnreadableCopy() {
				std::error_code ignored;
				std::filesystem::permissions(path, std::filesystem::perms::owner_read, ignored);
				std::filesystem::remove(path, ignored);
			}
		} unreadable;
		miint::Hdf5FailureRecorder failure;
		REQUIRE_THROWS_AS(miint::BIOMReader::IsBIOM(unreadable.path.string()), std::runtime_error);
		REQUIRE(failure.Cause() == "Permission denied");
	}
#endif
}

namespace {
herr_t SilentHdf5Handler(hid_t, void *) {
	return 0;
}
} // namespace

TEST_CASE("Hdf5FailureRecorder keeps the first failure's cause and restores HDF5's handler", "[BIOMReader]") {
	// A handler and data of the test's own, so a recorder that leaves any other handler behind is caught
	int handler_data = 0;
	H5Eset_auto2(H5E_DEFAULT, SilentHdf5Handler, &handler_data);
	{
		miint::Hdf5FailureRecorder failure;
		REQUIRE(failure.Cause().empty());
		// An HDF5 file without BIOM's datasets: the first failure is the missing dataset, and neither the cleanup
		// closes that follow it nor a later failure replaces the recorded cause
		REQUIRE_THROWS_AS(miint::BIOMReader("data/biom/notbiom.h5"), std::runtime_error);
		REQUIRE_THAT(failure.Cause(), Catch::Matchers::ContainsSubstring("component not found"));
		REQUIRE_THROWS_AS(miint::BIOMReader::IsBIOM("data/biom/truncated.biom"), std::runtime_error);
		REQUIRE_THAT(failure.Cause(), Catch::Matchers::ContainsSubstring("component not found"));
	}
	H5E_auto2_t func = nullptr;
	void *data = nullptr;
	H5Eget_auto2(H5E_DEFAULT, &func, &data);
	H5Eset_auto2(H5E_DEFAULT, nullptr, nullptr);
	REQUIRE((func == SilentHdf5Handler));
	REQUIRE((data == &handler_data));
}

TEST_CASE("Hdf5FailureRecorder holds the HDF5 mutex while it is alive", "[BIOMReader]") {
	// Asked from another thread: try_lock on a mutex the calling thread holds is undefined
	auto free_elsewhere = [] {
		bool got = false;
		std::thread([&got] {
			got = miint::g_hdf5_mutex.try_lock();
			if (got) {
				miint::g_hdf5_mutex.unlock();
			}
		}).join();
		return got;
	};
	{
		miint::Hdf5FailureRecorder failure;
		REQUIRE_FALSE(free_elsewhere());
	}
	REQUIRE(free_elsewhere());
}

// An HDF5 file whose root has a format-version attribute of `n` integers, the first `first`
namespace {
std::filesystem::path WriteFormatVersionFile(const std::string &name, hsize_t n, int first) {
	auto path =
	    std::filesystem::temp_directory_path() /
	    ("miint_" + name + "_" + std::to_string(std::chrono::steady_clock::now().time_since_epoch().count()) + ".h5");
	hid_t file = H5Fcreate(path.string().c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	REQUIRE(file >= 0);
	hid_t space = H5Screate_simple(1, &n, nullptr);
	hid_t attr =
	    H5Acreate_by_name(file, "/", "format-version", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	std::vector<int> values(n, 0);
	values[0] = first;
	REQUIRE(H5Awrite(attr, H5T_NATIVE_INT, values.data()) >= 0);
	H5Aclose(attr);
	H5Sclose(space);
	H5Fclose(file);
	return path;
}
} // namespace

TEST_CASE("BIOM IsBIOM takes only a format-version of two integers starting with 2", "[BIOMReader]") {
	H5Eset_auto(H5E_DEFAULT, nullptr, nullptr);
	struct Case {
		std::string name;
		hsize_t n;
		int first;
		bool is_biom;
	};
	// Three integers would overrun a reader that reads the attribute into two
	for (const auto &c : {Case {"v2", 2, 2, true}, Case {"v1", 2, 1, false}, Case {"v2_three_ints", 3, 2, false},
	                      Case {"v2_one_int", 1, 2, false}}) {
		auto path = WriteFormatVersionFile(c.name, c.n, c.first);
		INFO(c.name);
		bool is_biom = miint::BIOMReader::IsBIOM(path.string());
		std::filesystem::remove(path);
		REQUIRE(is_biom == c.is_biom);
	}
}
