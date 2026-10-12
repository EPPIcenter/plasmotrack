//
// Created by Maxwell Murphy on 5/28/20.
//

#include "utils.h"


#include <fmt/core.h>
#include <sstream>

namespace transmission_nets::core::io {

    std::basic_regex<char> getInvalidPathChars() {
        static std::basic_regex invalid_path_chars = std::basic_regex(R"([\\\\'|/:*?\"<>|])");
        return invalid_path_chars;
    }

    fs::path getPathFromEnvVar(const char* envVar) {
        const auto dirStr = getenv(envVar);
        if (dirStr == nullptr) {
            fmt::print(stderr, "Please set {} environment variable.\n", envVar);
            exit(1);
        }

        fs::path dir{dirStr};

        fmt::print("Loading Path: {}\n", std::string(dir));
        fmt::print("Path Exists: {}\n", fs::exists(dir));

        if (!fs::is_directory(dir)) {
            fmt::print(stderr, "Please set {} to an appropriate directory.\n", envVar);
            exit(1);
        }

        return dir;
    }

    std::string makePathValid(const std::string& input) {
        return std::regex_replace(input, getInvalidPathChars(), "_");
    }

    std::vector<char> decompressGzipFile(const fs::path& filePath) {
        auto file = gzopen(filePath.c_str(), "rb");
        if (file == nullptr) {
            fmt::print(stderr, "Failed to open file\n");
            exit(1);
        }

        auto output = std::vector<char>();
        do {
            char buffer[1024 * 1024]; // 1MB buffer
            const int numRead = gzread(file, buffer, 1024 * 1024);
            output.insert(output.end(), buffer, buffer + numRead);
        } while (!gzeof(file));

        return output;
    }

}// namespace transmission_nets::core::io
