//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
#define TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H

#include <filesystem>
#include <map>
#include <ostream>
#include <sstream>
#include "LambdaLogger.h"
#include "AbstractOutput.h"

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;
    class FileOutput : public AbstractOutput {
    public:
        explicit FileOutput(fs::path outputPath, std::string header = "", bool resetOutput = true);
        ~FileOutput() override;

        FileOutput(const FileOutput&) = delete;
        FileOutput& operator=(FileOutput&) = delete;


        void reset() override;

        void initialize(bool reset = true) final;

        void finalize() override;

        void write(const std::string&) override;

        void write_buffer();


    protected:
//        static std::map<std::string, std::shared_ptr<std::ofstream>> outputFiles_;
        static std::map<std::string, std::shared_ptr<std::ofstream>>& getOutputFiles();
        long buffer_size_ = 4096;
        fs::path outputPath_;
        std::string header_;
        std::stringstream buffer_;
        std::shared_ptr<std::ofstream> outputFile_;
    };
}



#endif//TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
