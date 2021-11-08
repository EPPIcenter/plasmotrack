//
// Created by Maxwell Murphy on 5/27/20.
//

#include <fstream>
#include <iostream>
#include <utility>

#include "FileOutput.h"

namespace transmission_nets::core::io {

//    std::map<std::string, std::shared_ptr<std::ofstream>>*FileOutput::outputFiles_ = new std::map<std::string, std::shared_ptr<std::ofstream>>();

    /*
     * Takes as arguments a path and optionally a header to be written immediately.
     * Creates a new output to a file. If the path has already been used to create a file output, the same file handle will be used and the header will not be
     * written even if passed in.
     */
    FileOutput::FileOutput(fs::path outputPath, std::string header, bool resetOutput): outputPath_(std::move(outputPath)), header_(std::move(header)) {
        auto& outputFiles = getOutputFiles();
        if (outputFiles.contains(outputPath_)) {
            outputFile_ = outputFiles.at(outputPath_);
        } else {
            initialize(resetOutput);
        }
    }

    std::map<std::string, std::shared_ptr<std::ofstream>>& FileOutput::getOutputFiles() {
        static std::map<std::string, std::shared_ptr<std::ofstream>> outputFiles_{};
        return outputFiles_;
    }

    void FileOutput::initialize(bool resetOutput) {
        auto& outputFiles = getOutputFiles();
        if(!fs::exists(outputPath_.parent_path())) {
            fs::create_directory(outputPath_.parent_path());
        }

        outputFile_ = std::make_shared<std::ofstream>();
        if(resetOutput) {
            outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::trunc);
        } else {
            outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::app);
        }
        outputFiles.insert({outputPath_, outputFile_});
        if(!header_.empty()) {
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
    }

    void FileOutput::write_buffer() {
        outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::app);
        (*outputFile_) << buffer_.rdbuf();
        outputFile_->close();
        buffer_.str(std::string()); // clear the buffer
    }

    void FileOutput::write(const std::string& val) {
        buffer_ << val << "\n";
        if (buffer_.tellp() > buffer_size_) {
            write_buffer();
        }
    }

    FileOutput::~FileOutput() {
        write_buffer();
    }

}

