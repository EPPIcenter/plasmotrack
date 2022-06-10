//
// Created by Maxwell Murphy on 5/17/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOGPQ_H
#define TRANSMISSION_NETWORKS_APP_LOGPQ_H

#include <vector>

namespace transmission_nets::core::utils {
    struct LogPQ {
        LogPQ(const std::vector<double>& x);
        std::vector<double> logP{};
        std::vector<double> logQ{};
    };
}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_LOGPQ_H
