//
// See DurableFile.h.
//

#include "DurableFile.h"
#include "OutputFailure.h"

#include <fcntl.h>
#include <unistd.h>

#include <cerrno>
#include <cstring>
#include <system_error>

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;

    void createDirectories(const fs::path& dir) {
        std::error_code ec;
        fs::create_directories(dir, ec);
        if (ec) {
            failOutput("create the directory of", dir, ec.message());
        }
    }

    void sync(const fs::path& path) {
        const int fd = ::open(path.c_str(), O_RDONLY);
        if (fd < 0) {
            failOutput("open", path, std::strerror(errno));
        }
        if (::fsync(fd) != 0) {
            const int err = errno;
            ::close(fd);
            failOutput("sync", path, std::strerror(err));
        }
        ::close(fd);
    }

    void writeFile(const fs::path& path, const std::string_view bytes) {
        const int fd = ::open(path.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644);
        if (fd < 0) {
            failOutput("open", path, std::strerror(errno));
        }
        std::size_t written = 0;
        while (written < bytes.size()) {
            const ssize_t n = ::write(fd, bytes.data() + written, bytes.size() - written);
            if (n < 0) {
                if (errno == EINTR) {
                    continue;
                }
                const int err = errno;
                ::close(fd);
                failOutput("write", path, std::strerror(err));
            }
            written += static_cast<std::size_t>(n);
        }
        if (::close(fd) != 0) {
            failOutput("write", path, std::strerror(errno));
        }
    }

    void moveDurably(const fs::path& temporary, const fs::path& path) {
        sync(temporary);
        std::error_code ec;
        fs::rename(temporary, path, ec);
        if (ec) {
            failOutput("rename into place", path, ec.message());
        }
        sync(path.has_parent_path() ? path.parent_path() : fs::path("."));
    }

    void writeDurably(const fs::path& path, const std::function<void(const fs::path&)>& write) {
        const fs::path dir = path.has_parent_path() ? path.parent_path() : fs::path(".");
        createDirectories(dir);
        const auto temporary = dir / ("." + path.filename().string() + ".tmp");
        write(temporary);
        moveDurably(temporary, path);
    }
}// namespace transmission_nets::core::io
