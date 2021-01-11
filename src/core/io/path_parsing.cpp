//
// Created by Maxwell Murphy on 5/28/20.
//
#include <iostream>

#include "path_parsing.h"

namespace transmission_nets::core::io {

    fs::path getPathFromEnvVar(const char *envVar) {
        const auto dirStr = getenv(envVar);
        if (dirStr == nullptr) {
            std::cerr << "Please set " << envVar << " environment variable.\n";
            exit(1);
        }

        fs::path dir{dirStr};

        std::cout << "Loading Path: " << dir << std::endl;

        std::cout<< "Path Exists: " <<  fs::exists(dir) << std::endl;

        if(!fs::is_directory(dir)) {
            std::cerr << "Please set " << envVar << " to an appropriate directory.\n";
            exit(1);
        }

        return dir;
    }

}

