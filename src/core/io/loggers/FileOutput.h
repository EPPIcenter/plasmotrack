//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
#define TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H

#include "AbstractOutput.h"
#include "LambdaLogger.h"
#include <filesystem>
#include <map>
#include <ostream>
#include <sstream>

namespace transmission_nets::core::io {
    namespace fs = std::filesystem;
    class FileOutput : public AbstractOutput {
    public:
        explicit FileOutput(fs::path outputPath, std::string header = "", bool resetOutput = true);
        ~FileOutput() override;

        FileOutput(const FileOutput&) = delete;
        FileOutput& operator=(FileOutput&) = delete;

        void reset() final;

        void initialize(bool reset) final;

        void finalize() final;

        void write(const std::string&) final;

        void write_buffer();


    protected:
        //        static std::map<std::string, std::shared_ptr<std::ofstream>> outputFiles_;
        static std::map<std::string, std::shared_ptr<std::ofstream>>& getOutputFiles();
        static std::map<std::string, std::shared_ptr<std::stringstream>>& getOutputBuffers();
        long buffer_size_ = 4096;
        fs::path outputPath_;
        std::string header_;
        std::shared_ptr<std::stringstream> pBuffer_;
        std::shared_ptr<std::ofstream> outputFile_;
    };
}// namespace transmission_nets::core::io


#endif//TRANSMISSION_NETWORKS_APP_FILEOUTPUT_H
