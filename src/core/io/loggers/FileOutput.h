//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
#define TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H

#include <filesystem>
#include <map>
#include <ostream>
#include "LambdaLogger.h"
#include "AbstractOutput.h"

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;
    class FileOutput : public AbstractOutput {
    public:
        explicit FileOutput(fs::path outputPath, std::string header = "");

        void reset() override;

        void initialize() override;

        void write(const std::string&) const override;

    protected:
        static std::map<std::string, std::shared_ptr<std::ofstream>>* outputFiles_;
        fs::path outputPath_;
        std::string header_;
        std::shared_ptr<std::ofstream> outputFile_;
    };
}



#endif//TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
