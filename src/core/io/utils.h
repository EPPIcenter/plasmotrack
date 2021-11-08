//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_UTILS_H
#define TRANSMISSION_NETWORKS_APP_UTILS_H

#include "core/io/serialize.h"

#include <boost/algorithm/string.hpp>
#include <fmt/core.h>

#include <fstream>
#include <iostream>
#include <regex>
#include <filesystem>
#include <vector>

namespace transmission_nets::core::io {

    namespace fs = std::filesystem;
    fs::path getPathFromEnvVar(const char *envVar);

    std::string getLastLine(std::ifstream& in);
    double hotloadParameter(const fs::path& input);
    std::vector<double> hotloadVector(const fs::path& input);
    std::string hotloadString(const fs::path& input);
    std::string makePathValid(const std::string& input);

}


#endif//TRANSMISSION_NETWORKS_APP_UTILS_H
