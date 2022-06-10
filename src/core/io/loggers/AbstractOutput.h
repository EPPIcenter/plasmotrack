//
// Created by Maxwell Murphy on 1/13/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_ABSTRACTOUTPUT_H
#define TRANSMISSION_NETWORKS_APP_ABSTRACTOUTPUT_H

#include <string>

namespace transmission_nets::core::io {
    class AbstractOutput {
    public:
        virtual void write(const std::string& val) = 0;
        virtual void initialize([[maybe_unused]] bool resetOutput = true){};
        virtual void finalize(){};
        virtual void reset(){};
        virtual ~AbstractOutput() = default;
    };
}// namespace transmission_nets::core::io

#endif//TRANSMISSION_NETWORKS_APP_ABSTRACTOUTPUT_H
