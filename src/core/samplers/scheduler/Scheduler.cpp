//
// Created by Maxwell Murphy on 4/19/20.
//

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    void Scheduler::registerSampler(ScheduledSampler sampler) {
        index_.add(sampler.id, *sampler.sampler);
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

    void Scheduler::step(const bool adapt) {
        for (auto& sampler : samplers_) {
            update(sampler);
            if (adapt) {
                this->adapt(sampler);
            }
        }
        ++totalSteps;
    }

    const std::vector<ScheduledSampler>& Scheduler::samplers() const noexcept {
        return samplers_;
    }

    const std::vector<std::string>& Scheduler::samplerIds() const noexcept {
        return index_.ids();
    }

    AbstractSampler& Scheduler::sampler(const std::string& id) const {
        return index_.at(id);
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
