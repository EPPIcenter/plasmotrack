//
// Created by mmurphy on 2/6/24.
//

#include "CompressedFileOutput.h"

namespace transmission_nets::core::io {

    CompressedFileOutput::CompressedFileOutput(fs::path outputPath, std::string header, bool resetOutput) : outputPath_(std::move(outputPath)), header_(std::move(header)) {
        const auto& outputBuffers = getOutputBuffers();
        if (outputBuffers.contains(outputPath_)) {
            buffer_ = outputBuffers.at(outputPath_);
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

    std::map<std::string, std::shared_ptr<std::stringstream>>& CompressedFileOutput::getOutputBuffers() {
        static std::map<std::string, std::shared_ptr<std::stringstream>> outputBuffers_{};
        return outputBuffers_;
    }

    void CompressedFileOutput::initialize(const bool reset) {
        auto& outputBuffers = getOutputBuffers();
        if (!exists(outputPath_.parent_path())) {
            create_directory(outputPath_.parent_path());
        }

        buffer_     = std::make_shared<std::stringstream>();
        outputBuffers.insert({outputPath_, buffer_});

        if (!header_.empty() and reset) {
            this->reset();
            write(header_);
        }
    }

    void CompressedFileOutput::finalize() {
        write_buffer();
    }

    void CompressedFileOutput::reset() {
        const auto outputFile_ = gzopen(outputPath_.c_str(), "wb");
        if (outputFile_ == nullptr) {
            fmt::print(stderr, "Failed to open file\n");
            exit(0);
        }
        gzclose(outputFile_);
    }

    void CompressedFileOutput::write(const std::string& data) {
        *buffer_ << data << "\n";
        if (buffer_->tellp() > COMPRESSED_CHUNK_SIZE) {
            write_buffer();
        }
    }

    void CompressedFileOutput::write_buffer() const {
        const auto outputFile_ = gzopen(outputPath_.c_str(), "ab");
        if (outputFile_ == nullptr) {
            fmt::print(stderr, "Failed to open file\n");
            exit(0);
        }
        gzwrite(outputFile_, buffer_->str().c_str(), buffer_->str().size());
        gzclose(outputFile_);
        buffer_->str(std::string());
    }

    CompressedFileOutput::~CompressedFileOutput() {
        write_buffer();
    }

}