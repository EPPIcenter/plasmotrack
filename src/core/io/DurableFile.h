//
// Writes output files so that a crash leaves either the old file or the whole new one under a file's name. Every
// failure ends the process through failOutput, naming the path; see OutputFailure.h.
//

#ifndef TRANSMISSION_NETWORKS_APP_DURABLEFILE_H
#define TRANSMISSION_NETWORKS_APP_DURABLEFILE_H

#include <filesystem>
#include <functional>
#include <string_view>

namespace transmission_nets::core::io {
    // `write` fills a temporary file next to `path`, named `.<name>.tmp`; then moveDurably moves it to `path`. Creates
    // the directory if it does not exist.
    void writeDurably(const std::filesystem::path& path, const std::function<void(const std::filesystem::path& temporary)>& write);

    // Syncs the file written at `temporary` to disk, renames it to `path`, then syncs the directory, so the rename
    // never makes a file visible before its contents are on disk. `temporary` must be in the directory of `path`.
    void moveDurably(const std::filesystem::path& temporary, const std::filesystem::path& path);

    // Writes `bytes` to `path`, replacing any file there. Not durable on its own: use it inside writeDurably.
    void writeFile(const std::filesystem::path& path, std::string_view bytes);

    // Creates `dir` and any missing parents.
    void createDirectories(const std::filesystem::path& dir);

    // Flushes a file, or a directory's entries, to disk.
    void sync(const std::filesystem::path& path);
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_DURABLEFILE_H
