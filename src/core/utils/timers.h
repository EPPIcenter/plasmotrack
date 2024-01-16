//
// Created by Maxwell Murphy on 10/6/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_TIMERS_H
#define TRANSMISSION_NETWORKS_APP_TIMERS_H

#include <chrono>

namespace transmission_nets::core::utils::timers {

    using Time = std::chrono::high_resolution_clock;
    using dsec = std::chrono::duration<float>;
    using isec = std::chrono::duration<int>;
    using fsec = std::chrono::duration<float>;
    using lsec = std::chrono::duration<long>;

    auto time() {
        return Time::now();
    }

}

#endif//TRANSMISSION_NETWORKS_APP_TIMERS_H
