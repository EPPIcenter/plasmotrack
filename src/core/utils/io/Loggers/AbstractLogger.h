//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H

//#include <boost/filesystem.hpp>
#include <utility>
#include <filesystem>
#include <fstream>


//namespace fs = boost::filesystem;
namespace fs = std::filesystem;

class AbstractLogger {

public:

    explicit AbstractLogger(fs::path outputPath);

    void clearFile();

    virtual void initializeFile() {};

    virtual void logValue() = 0;

protected:
    fs::path outputPath_;
    std::ofstream outputFile_;
    bool initialized{false};
};

#endif //TRANSMISSION_NETWORKS_APP_ABSTRACTLOGGER_H
