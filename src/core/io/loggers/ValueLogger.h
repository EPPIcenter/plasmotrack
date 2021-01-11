//
// Created by Maxwell Murphy on 5/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_VALUELOGGER_H
#define TRANSMISSION_NETWORKS_APP_VALUELOGGER_H

#include <fstream>

#include "core/io/serialize.h"
#include "core/io/loggers/AbstractLogger.h"


namespace transmission_nets::core::io {

    template<typename T>
    class ValueLogger : public AbstractLogger {

    public:
        ValueLogger(fs::path outputPath, T& target);

        void logValue() noexcept override;

    private:
        T& target_;
        bool initialized{false};

    };

    template<typename T>
    ValueLogger<T>::ValueLogger(fs::path outputPath, T &target) : AbstractLogger(outputPath), target_(target) {
        if(!fs::exists(outputPath_.parent_path())) {
            fs::create_directory(outputPath_.parent_path());
        }

        if(fs::exists(outputPath_)) {
            initialized = true;
        }

        outputFile_.open(outputPath, std::ofstream::app);
    }

    template<typename T>
    void ValueLogger<T>::logValue() noexcept {
        outputFile_ << serialize(target_.value()) << "\n";
    }

}



#endif //TRANSMISSION_NETWORKS_APP_VALUELOGGER_H
