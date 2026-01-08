//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H

#include <boost/random.hpp>

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    struct WeightedScheduledSampler {
        std::unique_ptr<AbstractSampler> sampler;

        std::string id                = "Unknown";
        int adaptationStart           = 0;
        int adaptationEnd             = 0;
        bool scaledAdaptation         = false;
        long unsigned int updateStart = 0;
        long unsigned int updateEnd   = std::numeric_limits<int>::max();
        double weight            = 1.0;
        bool debug                    = false;

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
        explicit RandomizedScheduler(std::shared_ptr<Engine> rng, int numSamples);

        void registerSampler(std::unique_ptr<AbstractSampler> sampler);

        void registerSampler(WeightedScheduledSampler sampler);

        void step();

        void update(const WeightedScheduledSampler& sampler) const;

        void adapt(const WeightedScheduledSampler& sampler) const;

    private:
        std::shared_ptr<Engine> rng_;
        int num_samples_;

        std::vector<WeightedScheduledSampler> samplers_{};
        std::vector<double> cumulative_weight_{};
        int total_steps_         = 0;
        double total_weight_ = 0;
        boost::random::uniform_int_distribution<> dist_;

        bool cumulative_weights_calculated_{false};
        void calculateCumulativeWeights();
    };

    template<typename Engine>
    RandomizedScheduler<Engine>::RandomizedScheduler(std::shared_ptr<Engine> rng, int numSamples) : rng_(rng), num_samples_(numSamples) {}

    template<typename Engine>
    void RandomizedScheduler<Engine>::step() {
        if (!cumulative_weights_calculated_) {
            calculateCumulativeWeights();
        }

        for (int i = 0; i < num_samples_; ++i) {
            long int v = dist_(*rng_);
            auto it = std::ranges::lower_bound(cumulative_weight_, v);
            const int idx = it - cumulative_weight_.begin();
            update(samplers_[idx]);
            adapt(samplers_[idx]);
        }
        ++total_steps_;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(std::unique_ptr<AbstractSampler> sampler) {
        samplers_.push_back(WeightedScheduledSampler{.sampler = std::move(sampler)});
        total_weight_ += samplers_.back().weight;
        cumulative_weights_calculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(WeightedScheduledSampler sampler) {
        if (sampler.debug) {
            sampler.sampler->setDebug();
        }
        sampler.sampler->setIdentifier(sampler.id);
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
