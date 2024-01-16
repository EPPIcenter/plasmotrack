//
// Created by Maxwell Murphy on 5/17/21.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOGPQ_H
#define TRANSMISSION_NETWORKS_APP_LOGPQ_H

#include <vector>

namespace transmission_nets::core::utils {
    struct LogPQ {
        explicit LogPQ(const std::vector<float>& x);
        std::vector<float> logP{};
        std::vector<float> logQ{};
    };
}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_LOGPQ_H
