//
// Created by Maxwell Murphy on 5/28/20.
//

#include "path_parsing.h"

namespace fs = boost::filesystem;

fs::path getPathFromEnvVar(const char *envVar) {
    const auto dirStr = getenv(envVar);
    if (dirStr == NULL) {
        std::cerr << "Please set " << envVar << " environment variable.\n";
        exit(1);
    }

    fs::path dir{dirStr};

    if(!fs::exists(dir)) {
        std::cerr << "Please set " << envVar << " to an appropriate directory.\n";
        exit(1);
    }

    return dir;
}