//
// Created by Maxwell Murphy on 5/17/21.
//

#include "LogPQ.h"

#include <cmath>

namespace transmission_nets::core::utils {
    LogPQ::LogPQ(const std::vector<double>& x) {
        logP.reserve(x.size());
        logQ.reserve(x.size());
        for (const auto el : x) {
            double ex = std::exp(el);
            if (el < 0) {
                logQ.push_back(-std::log1p(ex));
                logP.push_back(logQ.back() + el);
            } else {
                logP.push_back(-std::log1p(1 / ex));
                logQ.push_back(logP.back() - el);
            }
        }
    }
}// namespace transmission_nets::core::utils
