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

        std::string id = "Unknown";
        int adaptationStart = 0;
        int adaptationEnd = 0;
        bool scaledAdaptation  = false;
        long unsigned int updateStart = 0;
        long unsigned int updateEnd = std::numeric_limits<int>::max();
        int updateFrequency = 1;
        double weight = 1.0;
        bool debug = false;

        void update() const;
        void adapt() const;
        void adapt(unsigned int step) const;
    };


    inline bool isBetween(unsigned long val, unsigned long lower, unsigned long upper) {
        return (val >= lower) && (val < upper);
    }

    inline bool isUpdateStep(int frequency, int current_step) {
        return ((current_step + 1) % frequency) == 0;
    }

    class Scheduler {

    public:
        explicit Scheduler(int numSamples);
        void registerSampler(std::unique_ptr<AbstractSampler> sampler);
        void registerSampler(ScheduledSampler sampler);

        void update(ScheduledSampler& sampler) const;

        void adapt(ScheduledSampler& sampler) const;

        void step();

        [[nodiscard]] const std::vector<ScheduledSampler>& samplers() const noexcept;

    protected:
        std::vector<ScheduledSampler> samplers_{};
        int numSamples = 0;
        int totalSteps = 0;
    };

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SCHEDULER_H
