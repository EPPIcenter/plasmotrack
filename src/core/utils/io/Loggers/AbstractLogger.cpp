//
// Created by Maxwell Murphy on 5/27/20.
//

#include "AbstractLogger.h"

AbstractLogger::AbstractLogger(fs::path outputPath) : outputPath_(std::move(outputPath)) {}

void AbstractLogger::clearFile() {
    outputFile_.close();
    outputFile_.open(outputPath_, fs::ofstream::out | fs::ofstream::trunc);
}
