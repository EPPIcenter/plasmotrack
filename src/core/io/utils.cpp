//
// Created by Maxwell Murphy on 5/28/20.
//

#include "utils.h"

namespace transmission_nets::core::io {
    auto invalid_path_chars = std::basic_regex(R"([\\\\'|/:*?\"<>|])");

    fs::path getPathFromEnvVar(const char *envVar) {
        const auto dirStr = getenv(envVar);
        if (dirStr == nullptr) {
            std::cerr << "Please set " << envVar << " environment variable.\n";
            exit(1);
        }

        fs::path dir{dirStr};

        std::cout << "Loading Path: " << dir << std::endl;

        std::cout << "Path Exists: " <<  fs::exists(dir) << std::endl;

        if(!fs::is_directory(dir)) {
            std::cerr << "Please set " << envVar << " to an appropriate directory.\n";
            exit(1);
        }

        return dir;
    }

    std::string makePathValid(const std::string& input) {
        return std::regex_replace(input, invalid_path_chars, "_");
    }

    std::string getLastLine(std::ifstream& in) {
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
        return lastLine;
    }

    double hotloadParameter(const fs::path& filePath) {

        double value;
        std::ifstream input(filePath);

        if (input) {
            std::string lastLine = getLastLine(input);
            value = std::stod(lastLine);
        } else {
            fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
            exit(1);
        }

//        fmt::print("Value: {}, file: {}\n", value, filePath.c_str());
        return value;
    }

    std::vector<double> hotloadVector(const fs::path& filePath) {
        std::vector<std::string> tokens;
        std::vector<double> values;
        std::ifstream input(filePath);

        if (input) {
            std::string lastLine = getLastLine(input);
            boost::split(tokens, lastLine, boost::is_any_of(","));
            values.reserve(tokens.size());
            std::transform(tokens.begin(), tokens.end(), std::back_inserter(values), [](const std::string& val) {return std::stod(val);});
        } else {
            fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
            exit(1);
        }

//        fmt::print("Value: {}, file: {}\n", serialize(values), filePath.c_str());
        return values;
    }

    std::string hotloadString(const fs::path& filePath) {
        std::string value;
        std::ifstream input(filePath);

        if (input) {
            value = getLastLine(input);
        } else {
            fmt::print(stderr, "Problem opening file {}\n", filePath.c_str());
            exit(1);
        }

        return value;
    }

}

