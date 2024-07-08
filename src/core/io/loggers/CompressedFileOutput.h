//
// Created by mmurphy on 2/6/24.
//

#ifndef COMPRESSEDFILEOUTPUT_H
#define COMPRESSEDFILEOUTPUT_H

#include "AbstractOutput.h"

#include <zlib.h>
#include <fmt/core.h>

#include <filesystem>
#include <map>
#include <sstream>
#include <vector>


namespace transmission_nets::core::io {
    namespace fs = std::filesystem;
    class CompressedFileOutput : public AbstractOutput {
    public:
        explicit CompressedFileOutput(fs::path outputPath, std::string header = "", bool resetOutput = true);

        ~CompressedFileOutput() override;
        CompressedFileOutput(const CompressedFileOutput&) = delete;
        CompressedFileOutput& operator=(CompressedFileOutput&) = delete;

        void reset() final;

        void initialize(bool reset) final;

        void finalize() override;

        void write(const std::string&) final;

        void write_buffer() const;

    protected:
        static std::map<std::string, std::shared_ptr<std::stringstream>>& getOutputBuffers();
        const long COMPRESSED_CHUNK_SIZE = 1024 * 8; // 8KB
        fs::path outputPath_;
        std::string header_;
        std::shared_ptr<std::stringstream> buffer_;
    };
};



#endif //COMPRESSEDFILEOUTPUT_H
