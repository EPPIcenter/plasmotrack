//
// Created by Maxwell Murphy on 10/6/22.
//

#ifndef TRANSMISSION_NETWORKS_APP_TIMERS_H
#define TRANSMISSION_NETWORKS_APP_TIMERS_H

#include <fmt/core.h>

#include <chrono>
#include <string>

namespace transmission_nets::core::utils::timers {

    using Time = std::chrono::high_resolution_clock;
    using dsec = std::chrono::duration<double>;
    using isec = std::chrono::duration<int>;
    using fsec = std::chrono::duration<double>;
    using lsec = std::chrono::duration<long>;

    inline auto time() {
        return Time::now();
    }

    // UTC, to the millisecond: 2026-10-11T12:00:00.000Z.
    inline std::string wallTime(const std::chrono::system_clock::time_point t) {
        const auto ms  = std::chrono::floor<std::chrono::milliseconds>(t);
        const auto day = std::chrono::floor<std::chrono::days>(ms);
        const std::chrono::year_month_day ymd{day};
        const std::chrono::hh_mm_ss time{ms - day};
        return fmt::format("{:04}-{:02}-{:02}T{:02}:{:02}:{:02}.{:03}Z", static_cast<int>(ymd.year()), static_cast<unsigned>(ymd.month()), static_cast<unsigned>(ymd.day()),
                           time.hours().count(), time.minutes().count(), time.seconds().count(), time.subseconds().count());
    }

}

#endif//TRANSMISSION_NETWORKS_APP_TIMERS_H
