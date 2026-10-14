//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_SCHEDULER_H

#include <memory>
#include <numeric>
#include <vector>

#include "core/samplers/AbstractSampler.h"
#include "core/samplers/scheduler/SamplerIndex.h"

namespace transmission_nets::core::samplers {

    struct ScheduledSampler {
        std::unique_ptr<AbstractSampler> sampler;

        std::string id;// `<type>/<target>`, unique within the scheduler
        // Adapts on the Sampling iterations in [adaptationStart, adaptationEnd), and updates every `updateFrequency` iterations.
        int adaptationStart = 0;
        int adaptationEnd   = 0;
        int updateFrequency = 1;

        void update() const;
        void adapt() const;
    };


    inline bool isBetween(const int val, const int lower, const int upper) {
        return (val >= lower) && (val < upper);
    }

    inline bool isUpdateStep(const int frequency, const int current_step) {
        return ((current_step + 1) % frequency) == 0;
    }

    class Scheduler {

    public:
        Scheduler() = default;
        // Throws std::invalid_argument when the sampler's id is already registered.
        void registerSampler(ScheduledSampler sampler);

        void update(ScheduledSampler& sampler) const;

        void adapt(ScheduledSampler& sampler) const;

        // Updates every sampler; samplers inside their adaptation window also adapt, unless `adapt` is false.
        void step(bool adapt = true);

        [[nodiscard]] const std::vector<ScheduledSampler>& samplers() const noexcept;

        // Sampler ids in registration order.
        [[nodiscard]] const std::vector<std::string>& samplerIds() const noexcept;

        // Throws std::out_of_range for an id that isn't registered.
        [[nodiscard]] AbstractSampler& sampler(const std::string& id) const;

        // Sampling iterations run so far, which decide when each sampler updates and adapts.
        [[nodiscard]] int iterations() const noexcept { return totalSteps; }

        // Continue from `iterations`, as when resuming from a checkpoint.
        void restoreIterations(const int iterations) noexcept { totalSteps = iterations; }

    protected:
        std::vector<ScheduledSampler> samplers_{};
        SamplerIndex index_{};
        int totalSteps = 0;
    };

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_SCHEDULER_H
