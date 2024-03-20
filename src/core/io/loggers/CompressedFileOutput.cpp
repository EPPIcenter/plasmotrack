//
// Created by mmurphy on 2/6/24.
//

#include "CompressedFileOutput.h"

namespace transmission_nets::core::io {

    CompressedFileOutput::CompressedFileOutput(fs::path outputPath, std::string header, bool resetOutput) : outputPath_(std::move(outputPath)), header_(std::move(header)) {
        const auto& outputFiles   = getOutputFiles();
        const auto& outputBuffers = getOutputBuffers();
        if (outputFiles.contains(outputPath_)) {
            outputFile_ = outputFiles.at(outputPath_);
            buffer_     = outputBuffers.at(outputPath_);
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

    std::map<std::string, std::shared_ptr<std::ofstream>>& CompressedFileOutput::getOutputFiles() {
        static std::map<std::string, std::shared_ptr<std::ofstream>> outputFiles_{};
        return outputFiles_;
    }

    std::map<std::string, std::shared_ptr<std::stringstream>>& CompressedFileOutput::getOutputBuffers() {
        static std::map<std::string, std::shared_ptr<std::stringstream>> outputBuffers_{};
        return outputBuffers_;
    }

    void CompressedFileOutput::initialize(const bool reset) {
        auto& outputFiles   = getOutputFiles();
        auto& outputBuffers = getOutputBuffers();

        if (!exists(outputPath_.parent_path())) {
            create_directory(outputPath_.parent_path());
        }

        outputFile_ = std::make_shared<std::ofstream>();
        buffer_     = std::make_shared<std::stringstream>();

        if (reset) {
            outputFile_->open(outputPath_, std::ofstream::binary | std::ofstream::trunc);
        } else {
            outputFile_->open(outputPath_, std::ofstream::binary | std::ofstream::app);
        }

        outputFiles.insert({outputPath_, outputFile_});
        outputBuffers.insert({outputPath_, buffer_});

        if (!header_.empty() and reset) {
            write(header_);
        }

        outputFile_->close();
    }

    void CompressedFileOutput::finalize() {
        write_buffer();
        outputFile_->close();
    }

    void CompressedFileOutput::reset() {
        outputFile_->close();
        outputFile_->open(outputPath_, std::ofstream::binary | std::ofstream::trunc);
    }

    void CompressedFileOutput::write(const std::string& data) {
        *buffer_ << data << "\n";
        if (buffer_->tellp() > COMPRESSED_CHUNK_SIZE - 256) {
            write_buffer();
        }
    }

    void CompressedFileOutput::write_buffer() {

        z_stream zStream_;

        zStream_.zalloc = Z_NULL;
        zStream_.zfree  = Z_NULL;
        zStream_.opaque = Z_NULL;
        zStream_.avail_in = -1;
        zStream_.next_in = Z_NULL;

        if (const int ret = deflateInit2(&zStream_, compressionLevel_, Z_DEFLATED, 15 | 16, 8, Z_DEFAULT_STRATEGY); ret != Z_OK) {
            fmt::print(stderr, "Failed to initialize zlib stream: {}\n", zStream_.msg);
            exit(0);
        }

        int flush;

        outputFile_->open(outputPath_, std::ofstream::binary | std::ofstream::app);
        if (outputFile_->fail()) {
            fmt::print(stderr, "Failed to open file\n");
            exit(0);
        }

        do {
            char inBuffer[COMPRESSED_CHUNK_SIZE];
            buffer_->read(inBuffer, COMPRESSED_CHUNK_SIZE);
            if (buffer_->bad()) {
                fmt::print(stderr, "Failed to read from buffer\n");
                exit(0);
            }
            zStream_.avail_in = buffer_->gcount();
            flush = buffer_->eof() ? Z_FINISH : Z_NO_FLUSH;
            zStream_.next_in = reinterpret_cast<Bytef*>(inBuffer);

            do {
                char outBuffer[COMPRESSED_CHUNK_SIZE];
                zStream_.avail_out = COMPRESSED_CHUNK_SIZE;
                zStream_.next_out = reinterpret_cast<Bytef*>(outBuffer);
                if (const int ret = deflate(&zStream_, flush); ret != Z_OK && ret != Z_STREAM_END) {
                    fmt::print(stderr, "Failed to compress data: {}\n", ret);
                    outputFile_->close();
                    exit(0);
                }
                const long have = COMPRESSED_CHUNK_SIZE - zStream_.avail_out;

                outputFile_->write(outBuffer, have);
                if (outputFile_->fail()) {
                    fmt::print(stderr, "Failed to write to file\n");
                    outputFile_->close();
                    exit(0);
                }
            } while (zStream_.avail_out == 0);
        } while (flush != Z_FINISH);

        deflateEnd(&zStream_);
        outputFile_->close();

        buffer_->clear();
        buffer_->str(std::string());
        buffer_->seekg(0, std::ios::beg);
    }

    CompressedFileOutput::~CompressedFileOutput() {
        write_buffer();
    }

}