// SPDX-License-Identifier: MIT
//
// ENA Webin secret type for DuckDB-miint.
//
// Registers a custom DuckDB Secret with type='ena' so users can write
//   CREATE SECRET my_ena (TYPE ENA, USER 'Webin-12345', PASSWORD ..., ENDPOINT 'test');
// Passwords (literal, env var, or file) are resolved at create time and
// stored in a redacted KeyValueSecret.

#pragma once

#include <cstdlib>
#include <fstream>
#include <stdexcept>
#include <string>

namespace duckdb {
class ExtensionLoader;
}

namespace miint {

// Resolve password indirection. Exactly one of literal/env_var_name/file_path
// must be non-empty. Throws std::runtime_error if zero or more than one
// are provided, if an env var is unset, or if a file is missing/empty.
// Only the first line of a password file is read; getline strips its trailing
// LF and we additionally strip a trailing CR (Windows-authored CRLF files).
//
// Inline so the test executable can call this without linking ena_secret.cpp,
// which depends on DuckDB symbols not present in the unit-test target.
inline std::string ResolvePasswordIndirection(const std::string &literal,
                                               const std::string &env_var_name,
                                               const std::string &file_path) {
	const int set_count = (!literal.empty() ? 1 : 0) + (!env_var_name.empty() ? 1 : 0) +
	                      (!file_path.empty() ? 1 : 0);
	if (set_count == 0) {
		throw std::runtime_error("ENA secret requires one of 'password', 'password_env', or 'password_file'");
	}
	if (set_count > 1) {
		throw std::runtime_error("ENA secret: specify exactly one of 'password', 'password_env', 'password_file'");
	}
	if (!literal.empty()) {
		return literal;
	}
	if (!env_var_name.empty()) {
		const char *v = std::getenv(env_var_name.c_str());
		if (!v) {
			throw std::runtime_error("ENA secret: env var '" + env_var_name + "' is not set");
		}
		return std::string(v);
	}
	// file_path
	std::ifstream f(file_path);
	if (!f.is_open()) {
		throw std::runtime_error("ENA secret: cannot open password file '" + file_path + "'");
	}
	std::string line;
	std::getline(f, line);
	// strip trailing CR (CRLF on Windows-authored files)
	if (!line.empty() && line.back() == '\r') {
		line.pop_back();
	}
	if (line.empty()) {
		throw std::runtime_error("ENA secret: password file '" + file_path + "' is empty");
	}
	return line;
}

// Register the ENA secret type and the 'config' provider.
void RegisterENASecretType(duckdb::ExtensionLoader &loader);

} // namespace miint
