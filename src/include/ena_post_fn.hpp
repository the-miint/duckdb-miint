// SPDX-License-Identifier: MIT
//
// Tiny shared header for the ENA submit-side post functor type. Kept separate
// so per-table insert headers (ena_projects_insert.hpp, ena_samples_insert.hpp,
// ...) don't have to depend on each other to import a shared typedef.

#pragma once

#include <functional>
#include <string>

namespace miint {

// Functor signature: (url, body, user, password, content_type) -> response_body.
// Production wires this to ENAClient::PostJSON; tests inject a fake.
using ENAPostFn = std::function<std::string(const std::string &url, const std::string &body, const std::string &user,
                                            const std::string &password, const std::string &content_type)>;

} // namespace miint
