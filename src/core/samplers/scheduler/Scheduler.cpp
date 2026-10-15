//
// Created by Maxwell Murphy on 4/19/20.
//

#include "Scheduler.h"

#include <fmt/core.h>

#include <stdexcept>

namespace transmission_nets::core::samplers {

    void Scheduler::registerSampler(ScheduledSampler sampler) {
        index_.add(sampler.id, *sampler.sampler);
        samplers_.push_back(std::move(sampler));
    }

    void Scheduler::update(ScheduledSampler& sampler) const {
        if (isUpdateStep(sampler.updateFrequency, totalSteps)) {
            sampler.update();
        }
    }

    void Scheduler::adapt(ScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.adaptationStart, sampler.adaptationEnd)) {
            sampler.adapt();
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

    SamplerSetRecord Scheduler::record() const {
        SamplerSetRecord saved{.iterations = totalSteps};
        saved.counts.reserve(samplers_.size());
        for (std::uint32_t ii = 0; ii < samplers_.size(); ++ii) {
            const auto& sampler = *samplers_[ii].sampler;
            saved.counts.push_back(sampler.counts());
            if (auto adaptation = sampler.adaptationState(); !adaptation.scales.empty() or adaptation.updates != 0) {
                saved.adaptation.emplace(ii, std::move(adaptation));
            }
        }
        return saved;
    }

    void Scheduler::restore(const SamplerSetRecord& saved) {
        if (saved.counts.size() != samplers_.size()) {
            throw std::invalid_argument(fmt::format("record has counts for {} samplers, not {}", saved.counts.size(), samplers_.size()));
        }
        for (const auto& [position, adaptation] : saved.adaptation) {
            if (position >= samplers_.size()) {
                throw std::invalid_argument(fmt::format("record has adaptation state for sampler {} of {}", position, samplers_.size()));
            }
        }
        totalSteps = saved.iterations;
        for (std::uint32_t ii = 0; ii < samplers_.size(); ++ii) {
            auto& sampler = *samplers_[ii].sampler;
            sampler.restoreCounts(saved.counts[ii]);
            if (const auto adaptation = saved.adaptation.find(ii); adaptation != saved.adaptation.end()) {
                sampler.restoreAdaptationState(adaptation->second);
            }
        }
    }

    void ScheduledSampler::update() const {
        sampler->update();
    }

    void ScheduledSampler::adapt() const {
        sampler->adapt();
    }

}// namespace transmission_nets::core::samplers
