//
// CborFile: writes CBOR straight to a file through a small buffer, so a value is never held whole in memory. The
// caller writes the structure item by item: a map or array head with its count, then that many items. Small parts
// can go through nlohmann_json whole. A write that fails ends the process with OUTPUT_FAILURE, naming the path (see
// OutputFailure.h).
//

#ifndef TRANSMISSION_NETWORKS_APP_CBORFILE_H
#define TRANSMISSION_NETWORKS_APP_CBORFILE_H

#include <nlohmann/json.hpp>

#include <cstdint>
#include <filesystem>
#include <string>
#include <string_view>
#include <vector>

namespace transmission_nets::core::io {
    class CborFile {
    public:
        // Opens `path` for writing, truncating it.
        explicit CborFile(std::filesystem::path path);
        CborFile(const CborFile&)            = delete;
        CborFile& operator=(const CborFile&) = delete;
        ~CborFile();

        // A map of `entries` key-value pairs, or an array of `items`, whose contents follow.
        void map(std::size_t entries);
        void array(std::size_t items);
        void text(std::string_view value);
        void bytes(const std::vector<std::uint8_t>& value);
        void number(std::uint64_t value);
        void number(std::int64_t value);
        void number(double value);// always as a 64-bit float
        template<typename Number>
        void numbers(const std::vector<Number>& values) {
            array(values.size());
            for (const auto value : values) {
                number(value);
            }
        }
        void texts(const std::vector<std::string>& values);
        // `small` encoded whole, as one item.
        void value(const nlohmann::json& small);

        // Flushes and closes; after this the file is complete. The destructor closes without flushing what remains.
        void close();

    private:
        static constexpr std::size_t BUFFER_BYTES = std::size_t{1} << 20;

        void head(std::uint8_t major, std::uint64_t value);
        void append(std::string_view raw);
        void flushIfFull();
        void flush();

        std::filesystem::path path_;
        int fd_;
        std::string buffer_;
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_CBORFILE_H
