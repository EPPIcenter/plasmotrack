//
// Created by Maxwell Murphy on 5/27/20.
//


#include "FileOutput.h"

#include <fmt/core.h>

#include <fstream>
#include <iostream>
#include <utility>

namespace transmission_nets::core::io {

    /*
     * Takes as arguments a path and optionally a header to be written immediately.
     * Creates a new output to a file. If the path has already been used to create a file output, the same file handle will be used and the header will not be
     * written even if passed in, unless resetOutput is true.
     */
    FileOutput::FileOutput(fs::path outputPath, std::string header, bool resetOutput) : outputPath_(std::move(outputPath)), header_(std::move(header)) {
        auto &outputFiles = getOutputFiles();
        auto &outputBuffers = getOutputBuffers();
        if (outputFiles.contains(outputPath_)) {
            outputFile_ = outputFiles.at(outputPath_);
            pBuffer_ = outputBuffers.at(outputPath_);
            if (resetOutput) {
                reset();
                if (!header_.empty()) {
                    write(header_);
                }
            }
        } else {
            initialize(resetOutput);
        }
    }

    std::map<std::string, std::shared_ptr<std::ofstream>> &FileOutput::getOutputFiles() {
        static std::map<std::string, std::shared_ptr<std::ofstream>> outputFiles_{};
        return outputFiles_;
    }

    std::map<std::string, std::shared_ptr<std::stringstream>> &FileOutput::getOutputBuffers() {
        static std::map<std::string, std::shared_ptr<std::stringstream>> outputBuffers_{};
        return outputBuffers_;
    }

    void FileOutput::initialize(bool resetOutput) {
        auto &outputFiles = getOutputFiles();
        auto &outputBuffers = getOutputBuffers();

        if (!fs::exists(outputPath_.parent_path())) {
            fs::create_directory(outputPath_.parent_path());
        }

        outputFile_ = std::make_shared<std::ofstream>();
        pBuffer_ = std::make_shared<std::stringstream>();

        if (resetOutput) {
            outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::trunc);
        } else {
            outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::app);
        }

        outputFiles.insert({outputPath_, outputFile_});
        outputBuffers.insert({outputPath_, pBuffer_});

        if (!header_.empty() and resetOutput) {
            write(header_);
        }

        outputFile_->close();
    }

    void FileOutput::finalize() {
        write_buffer();
    }

    void FileOutput::reset() {
        outputFile_->close();
        outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::trunc);
        outputFile_->close();
        pBuffer_->str(std::string());
    }

    void FileOutput::write_buffer() {
        outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::app);
        (*outputFile_) << pBuffer_->rdbuf();
        outputFile_->close();
        pBuffer_->str(std::string());// clear the buffer
    }

    void FileOutput::write(const std::string &val) {
        *pBuffer_ << val << "\n";
        if (pBuffer_->tellp() > buffer_size_) {
            write_buffer();
        }
    }

    FileOutput::~FileOutput() {
        write_buffer();
    }

}// namespace transmission_nets::core::io
