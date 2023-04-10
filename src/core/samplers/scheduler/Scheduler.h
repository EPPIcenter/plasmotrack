//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SCHEDULER_H

#include <memory>
#include <numeric>
#include <vector>

#include "core/samplers/AbstractSampler.h"

namespace transmission_nets::core::samplers {

    struct ScheduledSampler {
        std::unique_ptr<AbstractSampler> sampler;
        int adaptation_start_ = -1;
        int adaptation_end_   = -1;
        bool scaled_adaptation_{false};
        int update_start_   = -1;
        int update_end_     = std::numeric_limits<int>::infinity();
        int update_frequency_ = 0;

        void update() const;
        void adapt() const;
        void adapt(unsigned int step) const;
    };


    inline bool isBetween(unsigned long val, unsigned long lower, unsigned long upper) {
        return val >= lower and val < upper;
    }

    inline bool isUpdateStep(int frequency, int current_step) {
        return ((current_step + 1) % frequency) == 0;
    }

    class Scheduler {

    public:
        void registerSampler(std::unique_ptr<AbstractSampler> sampler);
        void registerSampler(ScheduledSampler sampler);

        void update(ScheduledSampler& sampler) const;

        void adapt(ScheduledSampler& sampler) const;

        void step();

        [[nodiscard]] const std::vector<ScheduledSampler>& samplers() const noexcept;

    protected:
        std::vector<ScheduledSampler> samplers_{};
        int totalSteps = 0;
    };

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SCHEDULER_H
