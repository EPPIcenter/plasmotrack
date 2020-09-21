//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SCHEDULER_H

#include <vector>

#include "AbstractSampler.h"

namespace transmission_nets::core::samplers {

    struct ScheduledSampler {

        AbstractSampler* sampler{};

        int adaptationStart = 0;
        int adaptationEnd = 0;
        bool scaledAdaptation{false};

        int updateStart = 0;
        int updateEnd = std::numeric_limits<int>::infinity();

        int updateFrequency = 1;

        void update() const;

        void adapt() const;

        void adapt(int step) const;
    };


    inline bool isBetween(int val, int lower, int upper) {
        return val >= lower and val < upper;
    }

    inline bool isUpdateStep(int frequency, int currentStep) {
        return ((currentStep + 1) % frequency) == 0;
    }

    class Scheduler {

    public:
        void registerSampler(AbstractSampler* sampler);

        void registerSampler(ScheduledSampler sampler);

        void update(const ScheduledSampler & sampler) const;

        void adapt(const ScheduledSampler & sampler) const;

        void step();

        [[nodiscard]] const std::vector<ScheduledSampler>& samplers() const noexcept;

    protected:

        std::vector<ScheduledSampler> samplers_{};
        int totalSteps = 0;
    };

}



#endif //TRANSMISSION_NETWORKS_APP_SCHEDULER_H
