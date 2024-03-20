//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_UTILS_H
#define TRANSMISSION_NETWORKS_APP_UTILS_H

#include "core/io/serialize.h"

#include <boost/algorithm/string.hpp>

#include <fmt/core.h>
#include <zlib.h>

#include <filesystem>
#include <fstream>
#include <regex>
#include <vector>

namespace transmission_nets::core::io {

    namespace fs = std::filesystem;
    fs::path getPathFromEnvVar(const char* envVar);

    std::string getLastLine(std::istream& in);
    std::string getCompressedLastLine(gzFile in);
    double hotloadDouble(const fs::path& filePath);
    std::vector<double> hotloadVector(const fs::path& filePath);
    std::string hotloadString(const fs::path& filePath);
    std::string makePathValid(const std::string& input);
    std::vector<char> decompressGzipFile(const fs::path& filePath);

    /*
     * Load a vector of doubles from a file.
     * Will concatenate all lines in the file into a single vector.
     * @param filePath The path to the file.
     * @return A vector of doubles.
     */
    template<typename T>
    std::vector<T> loadVectorFromFile(const fs::path& filePath) {
        std::vector<std::string> tokens;
        std::vector<T> values;
        std::ifstream input(filePath);

        if (input) {
            std::string line;
            while (std::getline(input, line)) {
                boost::split(tokens, line, boost::is_any_of(","));
                values.reserve(tokens.size());
                std::transform(tokens.begin(), tokens.end(), std::back_inserter(values), [](const std::string& val) { return std::stold(val); });
            }
        } else {
            fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
            exit(0);
        }

        return values;
    }

}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_UTILS_H
