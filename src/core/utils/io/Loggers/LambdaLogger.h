//
// Created by Maxwell Murphy on 6/22/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H
#define TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H

#include "AbstractLogger.h"

template<typename Loggable>
class LambdaLogger : public AbstractLogger {

public:
    LambdaLogger(const fs::path& outputPath, Loggable f);

    void logValue() noexcept override;

private:

    Loggable f_;
};

template<typename Loggable>
void LambdaLogger<Loggable>::logValue() noexcept {
    outputFile_ << f_() << "\n";
}

template<typename Loggable>
LambdaLogger<Loggable>::LambdaLogger(const fs::path &outputPath, Loggable f) : AbstractLogger(outputPath), f_(f) {
    if(!fs::exists(outputPath_.parent_path())) {
        fs::create_directory(outputPath_.parent_path());
    }

    if(fs::exists(outputPath_)) {
        initialized = true;
    }

    outputFile_.open(outputPath, fs::ofstream::app);
}


#endif //TRANSMISSION_NETWORKS_APP_LAMBDALOGGER_H
