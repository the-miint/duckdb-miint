#include "ena_secret.hpp"

#include <catch2/catch_all.hpp>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <string>

namespace {

// Helper: write a file with the given content, return its path.
// Caller is responsible for cleanup via std::remove.
std::string WriteTempPwFile(const std::string &content) {
	std::string path = std::string("/tmp/miint_test_pwfile_") + std::to_string(::getpid()) + "_" +
	                   std::to_string(rand());
	std::ofstream f(path);
	f << content;
	f.close();
	return path;
}

// RAII guard to set/restore an env var.
class EnvGuard {
public:
	EnvGuard(std::string name, const char *value) : name_(std::move(name)) {
		const char *prior = std::getenv(name_.c_str());
		if (prior) {
			prior_ = prior;
			had_prior_ = true;
		}
		if (value) {
			::setenv(name_.c_str(), value, 1);
		} else {
			::unsetenv(name_.c_str());
		}
	}
	~EnvGuard() {
		if (had_prior_) {
			::setenv(name_.c_str(), prior_.c_str(), 1);
		} else {
			::unsetenv(name_.c_str());
		}
	}

private:
	std::string name_;
	std::string prior_;
	bool had_prior_ {false};
};

} // namespace

TEST_CASE("ENA password indirection: literal", "[ena_secret]") {
	auto pw = miint::ResolvePasswordIndirection("hunter2", "", "");
	CHECK(pw == "hunter2");
}

TEST_CASE("ENA password indirection: env var found", "[ena_secret]") {
	EnvGuard g("MIINT_TEST_ENA_PW", "envpassword");
	auto pw = miint::ResolvePasswordIndirection("", "MIINT_TEST_ENA_PW", "");
	CHECK(pw == "envpassword");
}

TEST_CASE("ENA password indirection: env var unset throws", "[ena_secret]") {
	EnvGuard g("MIINT_TEST_ENA_UNSET", nullptr);
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("", "MIINT_TEST_ENA_UNSET", ""),
	                  Catch::Matchers::ContainsSubstring("MIINT_TEST_ENA_UNSET"));
}

TEST_CASE("ENA password indirection: file found", "[ena_secret]") {
	auto path = WriteTempPwFile("filepassword\n");
	auto pw = miint::ResolvePasswordIndirection("", "", path);
	std::remove(path.c_str());
	CHECK(pw == "filepassword");
}

TEST_CASE("ENA password indirection: file with no trailing newline", "[ena_secret]") {
	auto path = WriteTempPwFile("nonewline");
	auto pw = miint::ResolvePasswordIndirection("", "", path);
	std::remove(path.c_str());
	CHECK(pw == "nonewline");
}

TEST_CASE("ENA password indirection: file CRLF stripped", "[ena_secret]") {
	auto path = WriteTempPwFile("crlf\r\nignored");
	auto pw = miint::ResolvePasswordIndirection("", "", path);
	std::remove(path.c_str());
	CHECK(pw == "crlf");
}

TEST_CASE("ENA password indirection: missing file throws", "[ena_secret]") {
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("", "", "/tmp/__miint_test_no_such_file__"),
	                  Catch::Matchers::ContainsSubstring("/tmp/__miint_test_no_such_file__"));
}

TEST_CASE("ENA password indirection: empty file throws", "[ena_secret]") {
	auto path = WriteTempPwFile("");
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("", "", path),
	                  Catch::Matchers::ContainsSubstring(path));
	std::remove(path.c_str());
}

TEST_CASE("ENA password indirection: none of the three throws", "[ena_secret]") {
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("", "", ""),
	                  Catch::Matchers::ContainsSubstring("password"));
}

TEST_CASE("ENA password indirection: more than one throws", "[ena_secret]") {
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("a", "B", ""),
	                  Catch::Matchers::ContainsSubstring("exactly one"));
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("a", "", "/tmp/x"),
	                  Catch::Matchers::ContainsSubstring("exactly one"));
	CHECK_THROWS_WITH(miint::ResolvePasswordIndirection("", "B", "/tmp/x"),
	                  Catch::Matchers::ContainsSubstring("exactly one"));
}
