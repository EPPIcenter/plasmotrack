//
// Created by Maxwell Murphy on 4/19/20.
//

#include "Scheduler.h"

namespace transmission_nets::core::samplers {
    
    void Scheduler::registerSampler(AbstractSampler *sampler) {
        samplers_.push_back({.sampler=sampler});
    }

    void Scheduler::registerSampler(ScheduledSampler sampler) {
        samplers_.push_back(sampler);
    }

    void Scheduler::update(const ScheduledSampler &sampler) const {
        if (isUpdateStep(sampler.updateFrequency, totalSteps) and
            isBetween(totalSteps, sampler.updateStart, sampler.updateEnd)) {
            sampler.update();
        }
    }

    void Scheduler::adapt(const ScheduledSampler &sampler) const {
        if (isBetween(totalSteps, sampler.adaptationStart, sampler.adaptationEnd)) {
            if (sampler.scaledAdaptation) {
                sampler.adapt(totalSteps - sampler.updateStart);
            } else {
                sampler.adapt();
            }
        }
    }

    void Scheduler::step() {
        for (const auto &sampler : samplers_) {
            update(sampler);
            adapt(sampler);
        }
    }

    const std::vector<ScheduledSampler> &Scheduler::samplers() const noexcept {
        return samplers_;
    }

    void ScheduledSampler::update() const {
        sampler->update();
    }

    void ScheduledSampler::adapt() const {
        sampler->adapt();
    }

    void ScheduledSampler::adapt(int step) const {
        sampler->adapt(step);
    }

}

