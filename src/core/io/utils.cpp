//
// Created by Maxwell Murphy on 5/28/20.
//

#include "utils.h"


#include <fmt/core.h>

namespace transmission_nets::core::io {

    std::basic_regex<char> getInvalidPathChars() {
        static std::basic_regex invalid_path_chars = std::basic_regex(R"([\\\\'|/:*?\"<>|])");
        return invalid_path_chars;
    }

    fs::path getPathFromEnvVar(const char* envVar) {
        const auto dirStr = getenv(envVar);
        if (dirStr == nullptr) {
            fmt::print(stderr, "Please set {} environment variable.\n", envVar);
            exit(1);
        }

        fs::path dir{dirStr};

        fmt::print("Loading Path: {}\n", std::string(dir));
        fmt::print("Path Exists: {}\n", fs::exists(dir));

        if (!fs::is_directory(dir)) {
            fmt::print(stderr, "Please set {} to an appropriate directory.\n", envVar);
            exit(1);
        }

        return dir;
    }

    std::string makePathValid(const std::string& input) {
        return std::regex_replace(input, getInvalidPathChars(), "_");
    }

    std::string getLastLine(std::istream& in) {
        in.seekg(-1, std::ios_base::end);
        if (in.peek() == '\n') {
            in.seekg(-1, std::ios_base::cur);
            for (long i = in.tellg(); i >= 0; --i) {
                if (in.peek() == '\n') {
                    in.get();
                    break;
                }
                in.seekg(i, std::ios_base::beg);
            }
        }
        std::string lastLine;
        getline(in, lastLine);
        in.seekg(0, std::ios_base::beg);
        in.clear();
        return lastLine;
    }

    std::string getCompressedLastLine(gzFile in) {
        std::stringstream output;
        std::string lastLine;
        do {
            char buffer[1024 * 1024];
            const int numRead = gzread(in, buffer, 1024 * 1024);
            if (gzeof(in)) {
                output.write(buffer, numRead);
                lastLine = getLastLine(output);
            }
        } while (!gzeof(in));
        return lastLine;
    }

    double hotloadDouble(const fs::path& filePath) {
        std::string lastLine;
        if (filePath.extension() == ".gz") {
            lastLine = getCompressedLastLine(gzopen(filePath.c_str(), "rb"));
        } else {
            if (std::ifstream input(filePath); input) {
                lastLine = getLastLine(input);
            } else {
                fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
                exit(1);
            }
        }
        return std::stod(lastLine);
   }

    std::vector<double> hotloadVector(const fs::path& filePath) {
        std::vector<std::string> tokens;
        std::vector<double> values;

        std::string lastLine;
        if (filePath.extension() == ".gz") {
            lastLine = getCompressedLastLine(gzopen(filePath.c_str(), "rb"));
        } else {
            if (std::ifstream input(filePath); input) {
                lastLine = getLastLine(input);
            } else {
                fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
                exit(1);
            }
        }
        split(tokens, lastLine, boost::is_any_of(","));
        values.reserve(tokens.size());
        std::ranges::transform(tokens, std::back_inserter(values), [](const std::string& val) { return std::stod(val); });
        return values;
    }

    std::string hotloadString(const fs::path& filePath) {
        std::string value;

        if (filePath.extension() == ".gz") {
            value = getCompressedLastLine(gzopen(filePath.c_str(), "rb"));
        } else {
            if (std::ifstream input(filePath); input) {
                value = getLastLine(input);
            } else {
                fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
                exit(1);
            }
        }

        return value;
    }

    std::vector<char> decompressGzipFile(const fs::path& filePath) {
        std::ifstream file(filePath, std::ios::binary);
        std::vector<char> buffer((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
        file.close();

        std::string decompressed;
        z_stream zs;
        zs.zalloc = Z_NULL;
        zs.zfree = Z_NULL;
        zs.opaque = Z_NULL;
        zs.avail_in = 0;
        zs.next_in = Z_NULL;
        if (inflateInit2(&zs, MAX_WBITS + 16) != Z_OK) {
            fmt::print(stderr, "inflateInit2 failed while reading file {}\n", filePath.c_str());
            exit(1);
        }

        zs.avail_in = buffer.size();
        zs.next_in = reinterpret_cast<Bytef*>(buffer.data());

        // create a temporary file to write the decompressed string to
        // std::string tempPath = filePath.string() + ".decompressed";
        std::vector<char> decompressedData;
        do {
            std::vector<char> buffer(4096);
            zs.avail_out = buffer.size();
            zs.next_out = reinterpret_cast<Bytef*>(buffer.data());

            int ret = inflate(&zs, Z_NO_FLUSH);
            if (ret == Z_NEED_DICT || ret == Z_DATA_ERROR || ret == Z_MEM_ERROR) {
                std::cerr << "Failed to decompress gzip file" << std::endl;
                inflateEnd(&zs);
                exit(1);
            }

            decompressedData.insert(decompressedData.end(), buffer.begin(), buffer.begin() + buffer.size() - zs.avail_out);
        } while (zs.avail_out == 0);

        inflateEnd(&zs);

        return decompressedData;
    }

}// namespace transmission_nets::core::io
