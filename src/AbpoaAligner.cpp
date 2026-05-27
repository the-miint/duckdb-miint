#include "AbpoaAligner.hpp"
#include <stdexcept>

extern "C" {
#include "abpoa.h"
}

namespace miint {

AbpoaAligner::AbpoaAligner() {
}

AbpoaAligner::~AbpoaAligner() {
}

AbpoaMsaResult AbpoaAligner::align(const std::vector<std::string> &names,
                                    const std::vector<std::string> &sequences,
                                    const AbpoaAlignParams &params) {
	(void)names;
	(void)sequences;
	(void)params;
	throw std::runtime_error("AbpoaAligner::align not yet implemented");
}

AbpoaConsensusResult AbpoaAligner::consensus(const std::vector<std::string> &names,
                                              const std::vector<std::string> &sequences,
                                              const AbpoaAlignParams &params) {
	(void)names;
	(void)sequences;
	(void)params;
	throw std::runtime_error("AbpoaAligner::consensus not yet implemented");
}

} // namespace miint
