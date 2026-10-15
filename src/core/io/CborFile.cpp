//
// CborFile. See CborFile.h.
//

#include "CborFile.h"

#include "core/io/OutputFailure.h"

#include <fcntl.h>
#include <unistd.h>

#include <bit>
#include <cerrno>
#include <cstring>

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;
    using nlohmann::json;

    CborFile::CborFile(fs::path path) : path_(std::move(path)), fd_(::open(path_.c_str(), O_WRONLY | O_CREAT | O_TRUNC, 0644)) {
        if (fd_ < 0) {
            failOutput("open", path_, std::strerror(errno));
        }
        buffer_.reserve(BUFFER_BYTES);
    }

    CborFile::~CborFile() {
        if (fd_ >= 0) {
            ::close(fd_);
        }
    }

    void CborFile::map(const std::size_t entries) { head(5, entries); }
    void CborFile::array(const std::size_t items) { head(4, items); }

    void CborFile::text(const std::string_view value) {
        head(3, value.size());
        append(value);
    }

    void CborFile::bytes(const std::vector<std::uint8_t>& value) {
        head(2, value.size());
        append({reinterpret_cast<const char*>(value.data()), value.size()});
    }

    void CborFile::number(const std::uint64_t value) { head(0, value); }

    void CborFile::number(const std::int64_t value) {
        if (value < 0) {
            head(1, static_cast<std::uint64_t>(-(value + 1)));
        } else {
            head(0, static_cast<std::uint64_t>(value));
        }
    }

    void CborFile::number(const double value) {
        const auto bits = std::bit_cast<std::uint64_t>(value);
        buffer_.push_back(static_cast<char>(0xfb));
        for (int shift = 56; shift >= 0; shift -= 8) {
            buffer_.push_back(static_cast<char>((bits >> shift) & 0xff));
        }
        flushIfFull();
    }

    void CborFile::texts(const std::vector<std::string>& values) {
        array(values.size());
        for (const auto& value : values) {
            text(value);
        }
    }

    void CborFile::value(const json& small) {
        const auto encoded = json::to_cbor(small);
        append({reinterpret_cast<const char*>(encoded.data()), encoded.size()});
    }

    void CborFile::close() {
        flush();
        const int fd = fd_;
        fd_          = -1;
        if (::close(fd) != 0) {
            failOutput("write", path_, std::strerror(errno));
        }
    }

    void CborFile::head(const std::uint8_t major, const std::uint64_t value) {
        const auto type = static_cast<std::uint8_t>(major << 5);
        if (value < 24) {
            buffer_.push_back(static_cast<char>(type | value));
        } else {
            const int width = value <= 0xff ? 1 : value <= 0xffff ? 2 : value <= 0xffffffff ? 4 : 8;
            buffer_.push_back(static_cast<char>(type | (width == 1 ? 24 : width == 2 ? 25 : width == 4 ? 26 : 27)));
            for (int shift = 8 * (width - 1); shift >= 0; shift -= 8) {
                buffer_.push_back(static_cast<char>((value >> shift) & 0xff));
            }
        }
        flushIfFull();
    }

    void CborFile::append(const std::string_view raw) {
        buffer_.append(raw);
        flushIfFull();
    }

    void CborFile::flushIfFull() {
        if (buffer_.size() >= BUFFER_BYTES) {
            flush();
        }
    }

    void CborFile::flush() {
        std::size_t written = 0;
        while (written < buffer_.size()) {
            const ssize_t n = ::write(fd_, buffer_.data() + written, buffer_.size() - written);
            if (n < 0) {
                if (errno == EINTR) {
                    continue;
                }
                failOutput("write", path_, std::strerror(errno));
            }
            written += static_cast<std::size_t>(n);
        }
        buffer_.clear();
    }
}// namespace transmission_nets::core::io
