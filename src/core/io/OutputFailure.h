//
// How a sink ends the Run when its output cannot be written.
//

#ifndef TRANSMISSION_NETWORKS_APP_OUTPUTFAILURE_H
#define TRANSMISSION_NETWORKS_APP_OUTPUTFAILURE_H

#include <filesystem>
#include <string>

namespace transmission_nets::core::io {
    // The exit code of a Run whose output could not be written.
    inline constexpr int OUTPUT_FAILURE = 3;

    // Prints "ERROR: cannot <action> output file <path>: <reason>" to stderr and ends the process with OUTPUT_FAILURE.
    // Safe from any thread: the process ends without exit handlers or static destructors another thread could race.
    [[noreturn]] void failOutput(const std::string& action, const std::filesystem::path& path, const std::string& reason);
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_OUTPUTFAILURE_H
