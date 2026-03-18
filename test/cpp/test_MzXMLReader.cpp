#include <MzXMLReader.hpp>
#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <cmath>
#include <cstdlib>
#include <fstream>

using miint::MzXMLReader;

// ===== Constructor errors (Cycle 2.1) =====

TEST_CASE("MzXMLReader: nonexistent file throws", "[mzxml]") {
	CHECK_THROWS_WITH(MzXMLReader("nonexistent.mzXML"), Catch::Matchers::ContainsSubstring("cannot open"));
}

TEST_CASE("MzXMLReader: non-XML file triggers parse error on read", "[mzxml]") {
	MzXMLReader reader("data/fastq/small_a.fq");
	CHECK_THROWS_WITH(reader.read_spectra(1), Catch::Matchers::ContainsSubstring("parse error"));
}

// ===== Empty mzXML (Cycle 2.2) =====

TEST_CASE("MzXMLReader: empty file returns empty batch", "[mzxml]") {
	MzXMLReader reader("data/mzxml/empty.mzXML");
	auto batch = reader.read_spectra(10);
	CHECK(batch.empty());
}

// ===== MS1 scan attributes (Cycle 2.3) =====

TEST_CASE("MzXMLReader: MS1 metadata from basic_3spectra", "[mzxml]") {
	MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
	auto batch = reader.read_spectra(1);
	REQUIRE(batch.size() == 1);

	CHECK(batch.spectrum_index[0] == 0);
	CHECK(batch.spectrum_id[0] == "scan=1");
	CHECK(batch.scan_number[0] == 1);
	CHECK(batch.scan_number_valid[0]);
	CHECK(batch.ms_level[0] == 1);
	CHECK(batch.spectrum_type[0] == "centroid");
	CHECK(batch.polarity[0] == "positive");
	CHECK(batch.base_peak_mz[0] == 200.0);
	CHECK(batch.base_peak_intensity[0] == 5000.0);
	CHECK(batch.total_ion_current[0] == 8000.0);
	CHECK(batch.lowest_mz[0] == 100.0);
	CHECK(batch.highest_mz[0] == 300.0);
	CHECK(batch.default_array_length[0] == 3);
}

// ===== Retention time (Cycle 2.4) =====

TEST_CASE("MzXMLReader: retention time parsing", "[mzxml]") {
	SECTION("PT90S = 1.5 min") {
		MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
		auto batch = reader.read_spectra(1);
		REQUIRE(batch.size() == 1);
		CHECK(batch.retention_time_valid[0]);
		CHECK_THAT(batch.retention_time[0], Catch::Matchers::WithinRel(1.5, 1e-9));
	}

	SECTION("PT108S = 1.8 min") {
		MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
		auto batch = reader.read_spectra(3);
		REQUIRE(batch.size() == 3);
		CHECK(batch.retention_time_valid[2]);
		CHECK_THAT(batch.retention_time[2], Catch::Matchers::WithinRel(1.8, 1e-9));
	}

	SECTION("PT1M30S = 1.5 min") {
		MzXMLReader reader("data/mzxml/combined_rt.mzXML");
		auto batch = reader.read_spectra(1);
		REQUIRE(batch.size() == 1);
		CHECK(batch.retention_time_valid[0]);
		CHECK_THAT(batch.retention_time[0], Catch::Matchers::WithinRel(1.5, 1e-9));
	}

	SECTION("PT1H30M = 90 min") {
		MzXMLReader reader("data/mzxml/hours_rt.mzXML");
		auto batch = reader.read_spectra(1);
		REQUIRE(batch.size() == 1);
		CHECK(batch.retention_time_valid[0]);
		CHECK_THAT(batch.retention_time[0], Catch::Matchers::WithinRel(90.0, 1e-9));
	}
}

// ===== Binary peak data (Cycle 2.5) =====

TEST_CASE("MzXMLReader: binary peak arrays", "[mzxml]") {
	MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
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

// ===== MS2 precursor info + isolation window (Cycle 2.6) =====

TEST_CASE("MzXMLReader: MS2 precursor info", "[mzxml]") {
	MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
	auto batch = reader.read_spectra(3);
	REQUIRE(batch.size() == 3);

	// Spectrum 1 (index 1): MS2, CID, precursor=200.0
	CHECK(batch.ms_level[1] == 2);
	CHECK(batch.precursor_mz_valid[1]);
	CHECK(batch.precursor_mz[1] == 200.0);
	CHECK(batch.precursor_charge_valid[1]);
	CHECK(batch.precursor_charge[1] == 2);
	CHECK(batch.precursor_intensity_valid[1]);
	CHECK(batch.precursor_intensity[1] == 5000.0);
	CHECK(batch.activation_method[1] == "CID");

	// Isolation window: target=200.0, lower=upper=1.5 (from windowWideness=3.0)
	CHECK(batch.isolation_window_target_valid[1]);
	CHECK(batch.isolation_window_target[1] == 200.0);
	CHECK(batch.isolation_window_lower_valid[1]);
	CHECK(batch.isolation_window_lower[1] == 1.5);
	CHECK(batch.isolation_window_upper_valid[1]);
	CHECK(batch.isolation_window_upper[1] == 1.5);
}

// ===== Filter string and scan window (Cycle 2.7) =====

TEST_CASE("MzXMLReader: filter string and scan window", "[mzxml]") {
	MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
	auto batch = reader.read_spectra(1);
	REQUIRE(batch.size() == 1);

	CHECK(batch.filter_string[0] == "FTMS + p NSI Full ms [100.00-500.00]");
	CHECK(batch.scan_window_lower_valid[0]);
	CHECK(batch.scan_window_lower[0] == 100.0);
	CHECK(batch.scan_window_upper_valid[0]);
	CHECK(batch.scan_window_upper[0] == 500.0);
}

// ===== ms1_scan_index flat layout (Cycle 2.8) =====

TEST_CASE("MzXMLReader: ms1_scan_index flat layout", "[mzxml]") {
	MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
	auto batch = reader.read_spectra(3);
	REQUIRE(batch.size() == 3);

	// MS1 at index 0: NULL ms1_scan_index
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);

	// MS2 at index 1: ms1_scan_index = 0
	CHECK(batch.ms1_scan_index_valid[1]);
	CHECK(batch.ms1_scan_index[1] == 0);

	// MS2 at index 2: ms1_scan_index = 0
	CHECK(batch.ms1_scan_index_valid[2]);
	CHECK(batch.ms1_scan_index[2] == 0);
}

// ===== ms1_scan_index nested layout (Cycle 2.9) =====

TEST_CASE("MzXMLReader: ms1_scan_index nested scans", "[mzxml]") {
	MzXMLReader reader("data/mzxml/nested_scans.mzXML");
	auto batch = reader.read_spectra(10);
	REQUIRE(batch.size() == 3);

	// Emission order: child MS2s complete before parent MS1
	// scan 2 (MS2) -> spectrum_index 0
	// scan 3 (MS2) -> spectrum_index 1
	// scan 1 (MS1) -> spectrum_index 2

	CHECK(batch.spectrum_id[0] == "scan=2");
	CHECK(batch.ms_level[0] == 2);
	CHECK(batch.spectrum_index[0] == 0);

	CHECK(batch.spectrum_id[1] == "scan=3");
	CHECK(batch.ms_level[1] == 2);
	CHECK(batch.spectrum_index[1] == 1);

	CHECK(batch.spectrum_id[2] == "scan=1");
	CHECK(batch.ms_level[2] == 1);
	CHECK(batch.spectrum_index[2] == 2);

	// Both MS2s should have ms1_scan_index = 2 (parent MS1)
	CHECK(batch.ms1_scan_index_valid[0]);
	CHECK(batch.ms1_scan_index[0] == 2);
	CHECK(batch.ms1_scan_index_valid[1]);
	CHECK(batch.ms1_scan_index[1] == 2);

	// MS1 has no ms1_scan_index
	CHECK_FALSE(batch.ms1_scan_index_valid[2]);
}

// ===== ms1_scan_index via precursorScanNum (Cycle 2.10) =====

TEST_CASE("MzXMLReader: ms1_scan_index via precursorScanNum", "[mzxml]") {
	MzXMLReader reader("data/mzxml/precursorscannum.mzXML");
	auto batch = reader.read_spectra(10);
	REQUIRE(batch.size() == 2);

	// Flat: MS1 (scan 1, index 0), MS2 (scan 2, index 1, precursorScanNum=1)
	CHECK(batch.ms_level[0] == 1);
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);

	CHECK(batch.ms_level[1] == 2);
	CHECK(batch.ms1_scan_index_valid[1]);
	CHECK(batch.ms1_scan_index[1] == 0);
}

// ===== Orphan MS2 (Cycle 2.11) =====

TEST_CASE("MzXMLReader: orphan MS2 has NULL ms1_scan_index", "[mzxml]") {
	MzXMLReader reader("data/mzxml/orphan_ms2.mzXML");
	auto batch = reader.read_spectra(10);
	REQUIRE(batch.size() == 2);

	CHECK(batch.ms_level[0] == 2);
	CHECK_FALSE(batch.ms1_scan_index_valid[0]);

	CHECK(batch.ms_level[1] == 2);
	CHECK_FALSE(batch.ms1_scan_index_valid[1]);
}

// ===== Missing optional metadata (Cycle 2.12) =====

TEST_CASE("MzXMLReader: missing optional metadata", "[mzxml]") {
	MzXMLReader reader("data/mzxml/missing_optional.mzXML");
	auto batch = reader.read_spectra(1);
	REQUIRE(batch.size() == 1);

	CHECK(batch.spectrum_type[0].empty());
	CHECK(batch.polarity[0].empty());
	CHECK_FALSE(batch.base_peak_mz_valid[0]);
	CHECK_FALSE(batch.base_peak_intensity_valid[0]);
	CHECK_FALSE(batch.total_ion_current_valid[0]);
	CHECK_FALSE(batch.lowest_mz_valid[0]);
	CHECK_FALSE(batch.highest_mz_valid[0]);
	CHECK_FALSE(batch.retention_time_valid[0]);
	CHECK(batch.filter_string[0].empty());
	CHECK_FALSE(batch.scan_window_lower_valid[0]);
	CHECK_FALSE(batch.scan_window_upper_valid[0]);
	CHECK_FALSE(batch.precursor_mz_valid[0]);

	// This file also omits compressionType on <peaks> — verify peaks still decode
	REQUIRE(batch.mz_array[0].size() == 1);
	CHECK(batch.mz_array[0][0] == 150.0);
	CHECK(batch.intensity_array[0][0] == 999.0);
}

// ===== Zero peaks (Cycle 2.13) =====

TEST_CASE("MzXMLReader: zero peaks", "[mzxml]") {
	MzXMLReader reader("data/mzxml/zero_peaks.mzXML");
	auto batch = reader.read_spectra(1);
	REQUIRE(batch.size() == 1);

	CHECK(batch.mz_array[0].empty());
	CHECK(batch.intensity_array[0].empty());
	CHECK(batch.default_array_length[0] == 0);
}

// ===== Format variants (Cycle 2.14) =====

TEST_CASE("MzXMLReader: format variants", "[mzxml]") {
	SECTION("compressed") {
		MzXMLReader reader("data/mzxml/compressed.mzXML");
		auto batch = reader.read_spectra(1);
		REQUIRE(batch.size() == 1);
		REQUIRE(batch.mz_array[0].size() == 3);
		CHECK(batch.mz_array[0][0] == 100.0);
		CHECK(batch.mz_array[0][1] == 200.0);
		CHECK(batch.mz_array[0][2] == 300.0);
		CHECK(batch.intensity_array[0][0] == 1000.0);
		CHECK(batch.intensity_array[0][1] == 5000.0);
		CHECK(batch.intensity_array[0][2] == 2000.0);
	}

	SECTION("float32") {
		MzXMLReader reader("data/mzxml/float32.mzXML");
		auto batch = reader.read_spectra(1);
		REQUIRE(batch.size() == 1);
		REQUIRE(batch.mz_array[0].size() == 3);
		CHECK_THAT(batch.mz_array[0][0], Catch::Matchers::WithinRel(100.0, 1e-5));
		CHECK_THAT(batch.mz_array[0][1], Catch::Matchers::WithinRel(200.0, 1e-5));
		CHECK_THAT(batch.mz_array[0][2], Catch::Matchers::WithinRel(300.0, 1e-5));
		CHECK_THAT(batch.intensity_array[0][0], Catch::Matchers::WithinRel(1000.0, 1e-5));
		CHECK_THAT(batch.intensity_array[0][1], Catch::Matchers::WithinRel(5000.0, 1e-5));
		CHECK_THAT(batch.intensity_array[0][2], Catch::Matchers::WithinRel(2000.0, 1e-5));
	}

	SECTION("profile, negative, HCD") {
		MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
		auto batch = reader.read_spectra(3);
		REQUIRE(batch.size() == 3);
		CHECK(batch.spectrum_type[2] == "profile");
		CHECK(batch.polarity[2] == "negative");
		CHECK(batch.activation_method[2] == "HCD");
	}
}

// ===== Batch reading (Cycle 2.15) =====

TEST_CASE("MzXMLReader: batch reading", "[mzxml]") {
	SECTION("batches of 2") {
		MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
		auto b1 = reader.read_spectra(2);
		CHECK(b1.size() == 2);
		auto b2 = reader.read_spectra(2);
		CHECK(b2.size() == 1);
		auto b3 = reader.read_spectra(2);
		CHECK(b3.size() == 0);
	}

	SECTION("ms1_scan_index persists across batches") {
		MzXMLReader reader("data/mzxml/basic_3spectra.mzXML");
		auto b1 = reader.read_spectra(1);
		REQUIRE(b1.size() == 1);
		CHECK(b1.ms_level[0] == 1);
		CHECK_FALSE(b1.ms1_scan_index_valid[0]);

		auto b2 = reader.read_spectra(1);
		REQUIRE(b2.size() == 1);
		CHECK(b2.ms_level[0] == 2);
		CHECK(b2.ms1_scan_index_valid[0]);
		CHECK(b2.ms1_scan_index[0] == 0);
	}
}

// ===== Real data validation (MassIVE MSV000080179) =====
// Ground truth extracted independently via Python xml.etree + struct.unpack:
//   total=5027, ms1=1010, ms2=4017
//   scan 32: first MS1 with peaks, 1 peak, mz[0]=322.4377136230, int[0]=89.0, RT=PT11.232S
//   scan 93: first MS2 with peaks, 304 peaks, precursor_mz=349.18398482, mz[0]=79.0552

TEST_CASE("MzXMLReader: real MassIVE data validation", "[mzxml][real]") {
	const char *path = std::getenv("MZXML_REAL_DATA");
	if (!path) {
		SKIP("MZXML_REAL_DATA not set (run via run_tests.sh to download)");
	}
	if (!std::ifstream(path).good()) {
		SKIP("MZXML_REAL_DATA file not found");
	}

	MzXMLReader reader(path);

	int total = 0, ms1_count = 0, ms2_count = 0;
	bool found_scan32 = false, found_scan93 = false;

	while (true) {
		auto batch = reader.read_spectra(2048);
		if (batch.empty()) {
			break;
		}

		for (size_t i = 0; i < batch.size(); i++) {
			total++;
			if (batch.ms_level[i] == 1) {
				ms1_count++;
			}
			if (batch.ms_level[i] == 2) {
				ms2_count++;
			}

			if (batch.scan_number[i] == 32 && batch.scan_number_valid[i]) {
				found_scan32 = true;
				CHECK(batch.ms_level[i] == 1);
				CHECK(batch.default_array_length[i] == 1);
				CHECK(batch.retention_time_valid[i]);
				CHECK_THAT(batch.retention_time[i], Catch::Matchers::WithinRel(11.232 / 60.0, 1e-6));
				REQUIRE(batch.mz_array[i].size() == 1);
				CHECK_THAT(batch.mz_array[i][0], Catch::Matchers::WithinRel(322.4377136230469, 1e-5));
				CHECK_THAT(batch.intensity_array[i][0], Catch::Matchers::WithinRel(89.0, 1e-5));
			}

			if (batch.scan_number[i] == 93 && batch.scan_number_valid[i]) {
				found_scan93 = true;
				CHECK(batch.ms_level[i] == 2);
				CHECK(batch.default_array_length[i] == 304);
				CHECK(batch.precursor_mz_valid[i]);
				CHECK_THAT(batch.precursor_mz[i], Catch::Matchers::WithinRel(349.18398482, 1e-6));
				CHECK(batch.precursor_intensity_valid[i]);
				CHECK_THAT(batch.precursor_intensity[i], Catch::Matchers::WithinRel(8924.0, 1e-5));
				REQUIRE(batch.mz_array[i].size() == 304);
				CHECK_THAT(batch.mz_array[i][0], Catch::Matchers::WithinRel(79.0551757812, 1e-4));
				CHECK_THAT(batch.intensity_array[i][0], Catch::Matchers::WithinRel(33.0, 1e-5));
			}
		}
	}

	CHECK(total == 5027);
	CHECK(ms1_count == 1010);
	CHECK(ms2_count == 4017);
	CHECK(found_scan32);
	CHECK(found_scan93);
}
