//
// Created by Maxwell Murphy on 5/27/20.
//

#include <fstream>
#include <iostream>
#include <utility>

#include "FileOutput.h"

namespace transmission_nets::core::io {

    std::map<std::string, std::shared_ptr<std::ofstream>>*FileOutput::outputFiles_ = new std::map<std::string, std::shared_ptr<std::ofstream>>();

    /*
     * Takes as arguments a path and optionally a header to be written immediately.
     * Creates a new output to a file. If the path has already been used to create a file output, the same file handle will be used and the header will not be
     * written even if passed in.
     */
    FileOutput::FileOutput(fs::path outputPath, std::string header) : outputPath_(std::move(outputPath)), header_(std::move(header)) {
        if (outputFiles_->contains(outputPath_)) {
            outputFile_ = outputFiles_->at(outputPath_);
            std::cout << "Using old path!" << std::endl;
        } else {
            initialize();
        }
    }

    void FileOutput::initialize() {
        if(!fs::exists(outputPath_.parent_path())) {
            fs::create_directory(outputPath_.parent_path());
        }

        outputFile_ = std::make_shared<std::ofstream>();
        outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::trunc);
        outputFiles_->insert({outputPath_, outputFile_});
        if(!header_.empty()) {
            write(header_);
        }
    }

    void FileOutput::reset() {
        outputFile_->close();
        outputFile_->open(outputPath_, std::ofstream::out | std::ofstream::trunc);
    }

    void FileOutput::write(const std::string& val) const {
        std::cout << "writing" << std::endl;
        (*outputFile_) << val << "\n";
    }

}

