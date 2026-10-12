//
// See OutputFailure.h.
//

#include "OutputFailure.h"

#include <fmt/core.h>

#include <cstdio>
#include <cstdlib>

namespace transmission_nets::core::io {
    void failOutput(const std::string& action, const std::filesystem::path& path, const std::string& reason) {
        fmt::print(stderr, "ERROR: cannot {} output file {}: {}\n", action, path.string(), reason);
        std::fflush(stderr);
        std::_Exit(OUTPUT_FAILURE);
    }
}// namespace transmission_nets::core::io
