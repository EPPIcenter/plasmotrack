//
// Created by Maxwell Murphy on 5/27/20.
//

#include <fstream>

#include "AbstractLogger.h"

namespace transmission_nets::core::io {

    AbstractLogger::AbstractLogger(fs::path outputPath) : outputPath_(std::move(outputPath)) {}

    void AbstractLogger::clearFile() {
        outputFile_.close();
        outputFile_.open(outputPath_, std::ofstream::out | std::ofstream::trunc);
    }

}

