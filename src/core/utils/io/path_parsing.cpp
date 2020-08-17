//
// Created by Maxwell Murphy on 5/28/20.
//

#include "path_parsing.h"


fs::path getPathFromEnvVar(const char *envVar) {
    const auto dirStr = getenv(envVar);
    if (dirStr == NULL) {
        std::cerr << "Please set " << envVar << " environment variable.\n";
        exit(1);
    }

    fs::path dir{dirStr};

    std::cout << "Loading Path: " << dir << std::endl;
    if(!fs::is_directory(dir)) {
        std::cerr << "Please set " << envVar << " to an appropriate directory.\n";
        exit(1);
    }

    return dir;
}