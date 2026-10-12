//
// SHA-256 of a file's bytes, as the Run manifest records each input's identity.
//

#ifndef TRANSMISSION_NETWORKS_APP_SHA256_H
#define TRANSMISSION_NETWORKS_APP_SHA256_H

#include <filesystem>
#include <string>
#include <string_view>

namespace transmission_nets::core::utils {
    // The digest as 64 lowercase hex digits.
    [[nodiscard]] std::string sha256Hex(std::string_view bytes);

    // Throws std::runtime_error if the file can't be read.
    [[nodiscard]] std::string sha256HexOfFile(const std::filesystem::path& file);
}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_SHA256_H
