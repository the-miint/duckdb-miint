#include <MzMLReader.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>

using miint::MzMLReader;

static const std::string DATA_DIR = "data/mzml/";

// ===== Error handling =====

TEST_CASE("MzMLReader rejects nonexistent file", "[mzml][reader]") {
	CHECK_THROWS_WITH(MzMLReader("nonexistent.mzML"), Catch::Matchers::ContainsSubstring("cannot open"));
}

TEST_CASE("MzMLReader rejects non-XML file", "[mzml][reader]") {
	// FASTQ file: constructor succeeds (lazy parsing), but read_spectra triggers XML parse error
	MzMLReader reader("data/fastq/small_a.fq");
	CHECK_THROWS_WITH(reader.read_spectra(1), Catch::Matchers::ContainsSubstring("parse error"));
}

// ===== Spectrum metadata (basic_3spectra.mzML) =====

TEST_CASE("MzMLReader parses MS1 spectrum metadata", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	CHECK(batch.spectrum_index[0] == 0);
	CHECK(batch.spectrum_id[0] == "scan=1");
	CHECK(batch.ms_level[0] == 1);
	CHECK(batch.retention_time_valid[0]);
	CHECK_THAT(batch.retention_time[0], Catch::Matchers::WithinRel(1.5, 1e-9));
	CHECK(batch.spectrum_type[0] == "centroid");
	CHECK(batch.polarity[0] == "positive");
	CHECK(batch.base_peak_mz_valid[0]);
	CHECK(batch.base_peak_mz[0] == 200.0);
	CHECK(batch.base_peak_intensity_valid[0]);
	CHECK(batch.base_peak_intensity[0] == 5000.0);
	CHECK(batch.total_ion_current_valid[0]);
	CHECK(batch.total_ion_current[0] == 8000.0);
	CHECK(batch.lowest_mz_valid[0]);
	CHECK(batch.lowest_mz[0] == 100.0);
	CHECK(batch.highest_mz_valid[0]);
	CHECK(batch.highest_mz[0] == 300.0);
	CHECK(batch.default_array_length[0] == 3);
	CHECK(batch.filter_string[0] == "FTMS + p NSI Full ms [100.00-500.00]");
	CHECK(batch.scan_window_lower_valid[0]);
	CHECK(batch.scan_window_lower[0] == 100.0);
	CHECK(batch.scan_window_upper_valid[0]);
	CHECK(batch.scan_window_upper[0] == 500.0);

	// MS1 has no precursor info
	CHECK_FALSE(batch.precursor_mz_valid[0]);
	// ms1_scan_index is NULL for MS1
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);
}

TEST_CASE("MzMLReader parses MS2 precursor info", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(3);

	REQUIRE(batch.size() == 3);

	// Second spectrum (index 1): MS2 with CID
	CHECK(batch.ms_level[1] == 2);
	CHECK(batch.precursor_mz_valid[1]);
	CHECK(batch.precursor_mz[1] == 200.0);
	CHECK(batch.precursor_charge_valid[1]);
	CHECK(batch.precursor_charge[1] == 2);
	CHECK(batch.precursor_intensity_valid[1]);
	CHECK(batch.precursor_intensity[1] == 5000.0);
	CHECK(batch.isolation_window_target_valid[1]);
	CHECK(batch.isolation_window_target[1] == 200.0);
	CHECK(batch.isolation_window_lower_valid[1]);
	CHECK(batch.isolation_window_lower[1] == 1.5);
	CHECK(batch.isolation_window_upper_valid[1]);
	CHECK(batch.isolation_window_upper[1] == 1.5);
	CHECK(batch.activation_method[1] == "CID");
	CHECK(batch.collision_energy_valid[1]);
	CHECK(batch.collision_energy[1] == 35.0);
}

TEST_CASE("MzMLReader converts RT seconds to minutes", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(3);

	REQUIRE(batch.size() == 3);
	// Spectrum 2: RT=108 seconds = 1.8 minutes
	CHECK(batch.retention_time_valid[2]);
	CHECK_THAT(batch.retention_time[2], Catch::Matchers::WithinRel(1.8, 1e-9));
}

TEST_CASE("MzMLReader handles negative polarity and profile type", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(3);

	REQUIRE(batch.size() == 3);
	CHECK(batch.spectrum_type[2] == "profile");
	CHECK(batch.polarity[2] == "negative");
	CHECK(batch.activation_method[2] == "HCD");
}

// ===== Binary array decoding =====

TEST_CASE("MzMLReader decodes mz and intensity arrays", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	CHECK(batch.mz_array[0][0] == 100.0);
	CHECK(batch.mz_array[0][1] == 200.0);
	CHECK(batch.mz_array[0][2] == 300.0);

	REQUIRE(batch.intensity_array[0].size() == 3);
	CHECK(batch.intensity_array[0][0] == 1000.0);
	CHECK(batch.intensity_array[0][1] == 5000.0);
	CHECK(batch.intensity_array[0][2] == 2000.0);
}

// ===== ms1_scan_index tracking =====

TEST_CASE("MzMLReader tracks ms1_scan_index", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(3);

	REQUIRE(batch.size() == 3);

	// MS1 at index 0: ms1_scan_index is NULL
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);

	// MS2 at index 1: ms1_scan_index = 0 (preceding MS1)
	CHECK(batch.ms1_scan_index_valid[1]);
	CHECK(batch.ms1_scan_index[1] == 0);

	// MS2 at index 2: ms1_scan_index = 0 (same preceding MS1)
	CHECK(batch.ms1_scan_index_valid[2]);
	CHECK(batch.ms1_scan_index[2] == 0);
}

// ===== Missing optional metadata =====

TEST_CASE("MzMLReader handles missing optional metadata", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "missing_optional.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	CHECK(batch.ms_level[0] == 1);
	CHECK(batch.spectrum_type[0].empty());
	CHECK(batch.polarity[0].empty());
	CHECK_FALSE(batch.base_peak_mz_valid[0]);
	CHECK_FALSE(batch.base_peak_intensity_valid[0]);
	CHECK_FALSE(batch.total_ion_current_valid[0]);
	CHECK_FALSE(batch.retention_time_valid[0]);
	CHECK_FALSE(batch.precursor_mz_valid[0]);
	CHECK(batch.filter_string[0].empty());
}

// ===== Zlib compressed data =====

TEST_CASE("MzMLReader decodes zlib compressed data", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "compressed.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	CHECK(batch.mz_array[0][0] == 100.0);
	CHECK(batch.mz_array[0][1] == 200.0);
	CHECK(batch.mz_array[0][2] == 300.0);
	REQUIRE(batch.intensity_array[0].size() == 3);
	CHECK(batch.intensity_array[0][0] == 1000.0);
	CHECK(batch.intensity_array[0][1] == 5000.0);
	CHECK(batch.intensity_array[0][2] == 2000.0);
}

// ===== 32-bit float data =====

TEST_CASE("MzMLReader decodes 32-bit float data", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "float32.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	CHECK_THAT(batch.mz_array[0][0], Catch::Matchers::WithinRel(100.0, 1e-5));
	CHECK_THAT(batch.mz_array[0][1], Catch::Matchers::WithinRel(200.0, 1e-5));
	CHECK_THAT(batch.mz_array[0][2], Catch::Matchers::WithinRel(300.0, 1e-5));
}

// ===== Referenceable param group resolution =====

TEST_CASE("MzMLReader resolves referenceableParamGroupRef", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "param_groups.mzML");
	auto batch = reader.read_spectra(1);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	CHECK(batch.mz_array[0][0] == 100.0);
	CHECK(batch.mz_array[0][1] == 200.0);
	CHECK(batch.mz_array[0][2] == 300.0);
	REQUIRE(batch.intensity_array[0].size() == 3);
	CHECK(batch.intensity_array[0][0] == 1000.0);
}

// ===== Empty spectrumList =====

TEST_CASE("MzMLReader handles empty spectrumList", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "empty_spectra.mzML");
	auto batch = reader.read_spectra(10);
	CHECK(batch.empty());
}

// ===== Sequential batch reading =====

TEST_CASE("MzMLReader supports sequential batch reading", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");

	auto batch1 = reader.read_spectra(2);
	REQUIRE(batch1.size() == 2);
	CHECK(batch1.spectrum_index[0] == 0);
	CHECK(batch1.spectrum_index[1] == 1);

	auto batch2 = reader.read_spectra(2);
	REQUIRE(batch2.size() == 1);
	CHECK(batch2.spectrum_index[0] == 2);

	auto batch3 = reader.read_spectra(2);
	CHECK(batch3.empty());
}

// ===== Chromatogram reading =====

TEST_CASE("MzMLReader reads chromatograms", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "with_chromatograms.mzML");

	auto chrom_batch = reader.read_chromatograms(10);
	REQUIRE(chrom_batch.size() == 3);

	// TIC
	CHECK(chrom_batch.chromatogram_index[0] == 0);
	CHECK(chrom_batch.chromatogram_id[0] == "TIC");
	CHECK(chrom_batch.chromatogram_type[0] == "TIC");
	CHECK_FALSE(chrom_batch.precursor_mz_valid[0]);
	CHECK_FALSE(chrom_batch.product_mz_valid[0]);
	REQUIRE(chrom_batch.time_array[0].size() == 4);
	CHECK(chrom_batch.time_array[0][0] == 0.5);
	CHECK(chrom_batch.time_array[0][3] == 2.0);
	REQUIRE(chrom_batch.intensity_array[0].size() == 4);
	CHECK(chrom_batch.intensity_array[0][0] == 10000.0);

	// SRM with precursor and product
	CHECK(chrom_batch.chromatogram_type[1] == "SRM");
	CHECK(chrom_batch.precursor_mz_valid[1]);
	CHECK(chrom_batch.precursor_mz[1] == 500.0);
	CHECK(chrom_batch.product_mz_valid[1]);
	CHECK(chrom_batch.product_mz[1] == 250.0);

	// BPC
	CHECK(chrom_batch.chromatogram_type[2] == "BPC");
}

TEST_CASE("MzMLReader returns empty chromatograms for spectrum-only file", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto chrom_batch = reader.read_chromatograms(10);
	CHECK(chrom_batch.empty());
}

// ===== Orphan MS2 (no preceding MS1) =====

TEST_CASE("MzMLReader handles orphan MS2 with no preceding MS1", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "orphan_ms2.mzML");
	auto batch = reader.read_spectra(10);

	REQUIRE(batch.size() == 2);

	// Both are MS2 with no preceding MS1 — ms1_scan_index should be NULL
	CHECK(batch.ms_level[0] == 2);
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);

	CHECK(batch.ms_level[1] == 2);
	CHECK_FALSE(batch.ms1_scan_index_valid[1]);
}

// ===== ms1_scan_index across batch boundaries =====

TEST_CASE("MzMLReader tracks ms1_scan_index across batch boundaries", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");

	// Read MS1 (index 0) in first batch
	auto batch1 = reader.read_spectra(1);
	REQUIRE(batch1.size() == 1);
	CHECK(batch1.ms_level[0] == 1);
	CHECK_FALSE(batch1.ms1_scan_index_valid[0]); // MS1 has no ms1_scan_index

	// Read MS2 (index 1) in second batch — should reference MS1 at index 0
	auto batch2 = reader.read_spectra(1);
	REQUIRE(batch2.size() == 1);
	CHECK(batch2.ms_level[0] == 2);
	CHECK(batch2.ms1_scan_index_valid[0]);
	CHECK(batch2.ms1_scan_index[0] == 0);

	// Read MS2 (index 2) in third batch — should still reference MS1 at index 0
	auto batch3 = reader.read_spectra(1);
	REQUIRE(batch3.size() == 1);
	CHECK(batch3.ms_level[0] == 2);
	CHECK(batch3.ms1_scan_index_valid[0]);
	CHECK(batch3.ms1_scan_index[0] == 0);
}

// ===== Bare mzML root (no indexedmzML wrapper) =====

TEST_CASE("MzMLReader reads bare mzML without indexedmzML wrapper", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "no_index.mzML");
	auto batch = reader.read_spectra(10);

	REQUIRE(batch.size() == 1);
	CHECK(batch.spectrum_index[0] == 0);
	CHECK(batch.ms_level[0] == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	CHECK(batch.mz_array[0][0] == 100.0);
	CHECK(batch.mz_array[0][1] == 200.0);
	CHECK(batch.mz_array[0][2] == 300.0);
}

// ===== Mixed precision (32-bit mz + 64-bit intensity) =====

TEST_CASE("MzMLReader handles mixed precision arrays", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "mixed_precision.mzML");
	auto batch = reader.read_spectra(10);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 3);
	// 32-bit floats have ~6 decimal digits of precision
	CHECK_THAT(batch.mz_array[0][0], Catch::Matchers::WithinRel(100.0, 1e-5));
	CHECK_THAT(batch.mz_array[0][1], Catch::Matchers::WithinRel(200.0, 1e-5));
	CHECK_THAT(batch.mz_array[0][2], Catch::Matchers::WithinRel(300.0, 1e-5));

	// 64-bit intensity should be exact
	REQUIRE(batch.intensity_array[0].size() == 3);
	CHECK(batch.intensity_array[0][0] == 1000.0);
	CHECK(batch.intensity_array[0][1] == 5000.0);
	CHECK(batch.intensity_array[0][2] == 2000.0);
}

// ===== Large spectrum (10K data points) =====

TEST_CASE("MzMLReader handles large spectrum with 10K data points", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "large_spectrum.mzML");
	auto batch = reader.read_spectra(10);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.mz_array[0].size() == 10000);
	REQUIRE(batch.intensity_array[0].size() == 10000);

	// First mz value: 100.0
	CHECK_THAT(batch.mz_array[0][0], Catch::Matchers::WithinRel(100.0, 1e-9));
	// Last mz value: 100.0 + 9999 * 0.1 = 1099.9
	CHECK_THAT(batch.mz_array[0][9999], Catch::Matchers::WithinRel(1099.9, 1e-9));
}

// ===== Zero intensity spectrum =====

TEST_CASE("MzMLReader handles zero intensity spectrum", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "zero_intensity.mzML");
	auto batch = reader.read_spectra(10);

	REQUIRE(batch.size() == 1);
	REQUIRE(batch.intensity_array[0].size() == 3);
	CHECK(batch.intensity_array[0][0] == 0.0);
	CHECK(batch.intensity_array[0][1] == 0.0);
	CHECK(batch.intensity_array[0][2] == 0.0);
	CHECK(batch.base_peak_intensity_valid[0]);
	CHECK(batch.base_peak_intensity[0] == 0.0);
	CHECK(batch.total_ion_current_valid[0]);
	CHECK(batch.total_ion_current[0] == 0.0);
}

// ===== Total count =====

TEST_CASE("MzMLReader reads all 3 spectra from basic file", "[mzml][reader]") {
	MzMLReader reader(DATA_DIR + "basic_3spectra.mzML");
	auto batch = reader.read_spectra(100);
	CHECK(batch.size() == 3);
}
