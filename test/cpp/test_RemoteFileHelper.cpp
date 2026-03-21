#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_string.hpp>
#include <fstream>
#include <remote_file_helper.hpp>

static const std::string TEST_FASTQ = "data/fastq/small_a.fq";

// ===== IsRemotePath =====

TEST_CASE("IsRemotePath detects https://", "[remote]") {
	CHECK(miint::RemoteFileHelper::IsRemotePath("https://example.com/file.fq"));
}

TEST_CASE("IsRemotePath detects http://", "[remote]") {
	CHECK(miint::RemoteFileHelper::IsRemotePath("http://example.com/file.fq"));
}

TEST_CASE("IsRemotePath detects s3://", "[remote]") {
	CHECK(miint::RemoteFileHelper::IsRemotePath("s3://bucket/file.fq"));
}

TEST_CASE("IsRemotePath returns false for local paths", "[remote]") {
	CHECK_FALSE(miint::RemoteFileHelper::IsRemotePath("/tmp/file.fq"));
	CHECK_FALSE(miint::RemoteFileHelper::IsRemotePath("relative/path.fq"));
	CHECK_FALSE(miint::RemoteFileHelper::IsRemotePath("file.fq"));
	CHECK_FALSE(miint::RemoteFileHelper::IsRemotePath(""));
}

// ===== ResolvedFile / ResolvedFileSet =====

TEST_CASE("ResolvedFileSet deletes temp files on destruction", "[remote]") {
	std::string temp_path = "/tmp/miint_test_cleanup.tmp";
	{
		std::ofstream ofs(temp_path);
		ofs << "test content";
	}
	CHECK(std::ifstream(temp_path).good());

	{
		miint::ResolvedFileSet rfs;
		rfs.Add(miint::ResolvedFile(temp_path, true));
	}

	CHECK_FALSE(std::ifstream(temp_path).good());
}

TEST_CASE("ResolvedFileSet does not delete non-temp files", "[remote]") {
	{
		miint::ResolvedFileSet rfs;
		rfs.Add(miint::ResolvedFile(TEST_FASTQ, false));
	}
	CHECK(std::ifstream(TEST_FASTQ).good());
}

TEST_CASE("ResolvedFileSet handles empty set", "[remote]") {
	miint::ResolvedFileSet rfs;
	CHECK(rfs.Files().empty());
}

// ===== ResolveToLocal (simple overload) =====

TEST_CASE("ResolveToLocal returns local path unchanged", "[remote]") {
	auto resolved = miint::RemoteFileHelper::ResolveToLocal(TEST_FASTQ);
	CHECK(resolved.local_path == TEST_FASTQ);
	CHECK_FALSE(resolved.is_temp);
}

TEST_CASE("ResolveToLocal throws for remote path without context", "[remote]") {
	CHECK_THROWS(miint::RemoteFileHelper::ResolveToLocal("https://example.com/file.fq"));
}
