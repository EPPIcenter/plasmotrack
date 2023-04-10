//
// Created by Maxwell Murphy on 4/19/20.
//

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    void Scheduler::registerSampler(std::unique_ptr<AbstractSampler> sampler) {
        samplers_.push_back(ScheduledSampler{.sampler = std::move(sampler)});
    }

    void Scheduler::registerSampler(ScheduledSampler sampler) {
        samplers_.push_back(std::move(sampler));
    }

    void Scheduler::update(ScheduledSampler& sampler) const {
        if (isUpdateStep(sampler.update_frequency_, totalSteps) and
            isBetween(totalSteps, sampler.update_start_, sampler.update_end_)) {
            sampler.update();
        }
    }

    void Scheduler::adapt(ScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.adaptation_start_, sampler.adaptation_end_)) {
            if (sampler.scaled_adaptation_) {
                sampler.adapt(totalSteps - sampler.update_start_);
            } else {
                sampler.adapt();
            }
        }
    }

    void Scheduler::step() {
        for (auto& sampler : samplers_) {
            update(sampler);
            adapt(sampler);
        }
    }

    const std::vector<ScheduledSampler>& Scheduler::samplers() const noexcept {
        return samplers_;
    }

    void ScheduledSampler::update() const {
        sampler->update();
    }

    void ScheduledSampler::adapt() const {
        sampler->adapt();
    }

    void ScheduledSampler::adapt(unsigned int step) const {
        sampler->adapt(step);
    }

}// namespace transmission_nets::core::samplers
