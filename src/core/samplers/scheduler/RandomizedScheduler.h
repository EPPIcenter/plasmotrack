//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H

#include <boost/random.hpp>

#include <optional>

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    struct WeightedScheduledSampler {
        std::unique_ptr<AbstractSampler> sampler;

        std::string id;// `<type>/<target>`, unique within the scheduler
        int adaptationStart           = 0;
        int adaptationEnd             = 0;
        bool scaledAdaptation         = false;
        long unsigned int updateStart = 0;
        long unsigned int updateEnd   = std::numeric_limits<int>::max();
        double weight            = 1.0;

        bool operator<(const WeightedScheduledSampler& rhs) const {
            return weight < rhs.weight;
        }

        bool operator>(const WeightedScheduledSampler& rhs) const {
            return rhs < *this;
        }

        bool operator<=(const WeightedScheduledSampler& rhs) const {
            return !(rhs < *this);
        }

        bool operator>=(const WeightedScheduledSampler& rhs) const {
            return !(*this < rhs);
        }


        void update() const {
            sampler->update();
        };

        void adapt() const {
            sampler->adapt();
        };

        void adapt(unsigned int step) const {
            sampler->adapt(step);
        };
    };


    template<typename Engine = boost::random::mt19937>
    class RandomizedScheduler {

    public:
        // Each Sampling iteration makes `updatesPerIteration` updates, each by a sampler drawn by weight. Without a
        // count it makes one per registered sampler, so every sampler gets one update opportunity in expectation when
        // the weights are equal, as with the sequential Scheduler.
        explicit RandomizedScheduler(std::shared_ptr<Engine> rng, std::optional<int> updatesPerIteration = std::nullopt);

        // Throws std::invalid_argument when the sampler's id is already registered.
        void registerSampler(WeightedScheduledSampler sampler);

        // Sampler ids in registration order.
        [[nodiscard]] const std::vector<std::string>& samplerIds() const noexcept { return index_.ids(); }

        // Throws std::out_of_range for an id that isn't registered.
        [[nodiscard]] AbstractSampler& sampler(const std::string& id) const { return index_.at(id); }

        // Updates randomly chosen samplers; samplers inside their adaptation window also adapt, unless `adapt` is false.
        void step(bool adapt = true);

        void update(const WeightedScheduledSampler& sampler) const;

        void adapt(const WeightedScheduledSampler& sampler) const;

    private:
        std::shared_ptr<Engine> rng_;
        std::optional<int> updates_per_iteration_;

        std::vector<WeightedScheduledSampler> samplers_{};
        SamplerIndex index_{};
        std::vector<double> cumulative_weight_{};
        int total_steps_         = 0;
        double total_weight_ = 0;
        boost::random::uniform_int_distribution<> dist_;

        bool cumulative_weights_calculated_{false};
        void calculateCumulativeWeights();
    };

    template<typename Engine>
    RandomizedScheduler<Engine>::RandomizedScheduler(std::shared_ptr<Engine> rng, std::optional<int> updatesPerIteration) : rng_(rng), updates_per_iteration_(updatesPerIteration) {}

    template<typename Engine>
    void RandomizedScheduler<Engine>::step(const bool adapt) {
        if (!cumulative_weights_calculated_) {
            calculateCumulativeWeights();
        }

        const auto updates = updates_per_iteration_.value_or(static_cast<int>(samplers_.size()));
        for (int i = 0; i < updates; ++i) {
            long int v = dist_(*rng_);
            auto it = std::ranges::lower_bound(cumulative_weight_, v);
            const int idx = it - cumulative_weight_.begin();
            update(samplers_[idx]);
            if (adapt) {
                this->adapt(samplers_[idx]);
            }
        }
        ++total_steps_;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(WeightedScheduledSampler sampler) {
        index_.add(sampler.id, *sampler.sampler);
        samplers_.push_back(std::move(sampler));

        total_weight_ += samplers_.back().weight;
        cumulative_weights_calculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::calculateCumulativeWeights() {
        std::ranges::sort(samplers_, std::greater<>());
        double total = 0;
        cumulative_weight_.clear();
        cumulative_weight_.reserve(samplers_.size());
        for (auto& sampler : samplers_) {
            total += sampler.weight;
            cumulative_weight_.push_back(total);
        }
        cumulative_weights_calculated_ = true;

        using bounds = boost::random::uniform_int_distribution<>::param_type;
        dist_.param(bounds(0, total_weight_));
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::update(const WeightedScheduledSampler& sampler) const {
        if (isBetween(total_steps_, sampler.updateStart, sampler.updateEnd)) {
            sampler.update();
        }
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::adapt(const WeightedScheduledSampler& sampler) const {
        if (isBetween(total_steps_, sampler.adaptationStart, sampler.adaptationEnd)) {
            if (sampler.scaledAdaptation) {
                sampler.adapt(total_steps_ - sampler.updateStart);
            } else {
                sampler.adapt();
            }
        }
    }

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
