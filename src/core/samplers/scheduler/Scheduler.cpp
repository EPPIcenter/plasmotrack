//
// Created by Maxwell Murphy on 4/19/20.
//

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    Scheduler::Scheduler(int numSamples) : numSamples(numSamples) {}

    void Scheduler::registerSampler(std::unique_ptr<AbstractSampler> sampler) {
        samplers_.push_back(ScheduledSampler{.sampler = std::move(sampler)});
    }

    void Scheduler::registerSampler(ScheduledSampler sampler) {
        if (sampler.debug) {
            sampler.sampler->setDebug();
        }
        samplers_.push_back(std::move(sampler));
    }

    void Scheduler::update(ScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.updateStart, sampler.updateEnd)) {
            sampler.update();
        }
    }

    void Scheduler::adapt(ScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.adaptationStart, sampler.adaptationEnd)) {
            if (sampler.scaledAdaptation) {
                sampler.adapt(totalSteps - sampler.updateStart);
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
        ++totalSteps;
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
