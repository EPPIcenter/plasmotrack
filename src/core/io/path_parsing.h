//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_PATH_PARSING_H
#define TRANSMISSION_NETWORKS_APP_PATH_PARSING_H


#include <filesystem>

namespace transmission_nets::core::io {

    namespace fs = std::filesystem;
    fs::path getPathFromEnvVar(const char *envVar);

}


#endif //TRANSMISSION_NETWORKS_APP_PATH_PARSING_H
