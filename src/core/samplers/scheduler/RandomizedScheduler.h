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
        long double weight            = 1.0;

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
        int numSamples_;

        std::vector<WeightedScheduledSampler> samplers_{};
        std::vector<long double> cumulativeWeight_{};
        int totalSteps           = 0;
        long double totalWeight_ = 0;
        boost::random::uniform_int_distribution<> dist_;

        bool cumulativeWeightsCalculated_{false};
        void calculateCumulativeWeights();
    };

    template<typename Engine>
    RandomizedScheduler<Engine>::RandomizedScheduler(std::shared_ptr<Engine> rng, int numSamples) : rng_(rng), numSamples_(numSamples) {}

    template<typename Engine>
    void RandomizedScheduler<Engine>::step() {
        if (!cumulativeWeightsCalculated_) {
            calculateCumulativeWeights();
        }

        long int v;
        for (int i = 0; i < numSamples_; ++i) {
            //            fmt::print("---------------------STEP-------------------------\n");
            v       = dist_(*rng_);
            auto it = std::lower_bound(cumulativeWeight_.begin(), cumulativeWeight_.end(), v);
            int idx = it - cumulativeWeight_.begin();
            //            fmt::print("Sampler: {}\n", samplers_[idx].id);
            update(samplers_[idx]);
            adapt(samplers_[idx]);
        }
        ++totalSteps;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(std::unique_ptr<AbstractSampler> sampler) {
        samplers_.push_back(WeightedScheduledSampler{.sampler = std::move(sampler)});
        totalWeight_ += samplers_.back().weight;
        cumulativeWeightsCalculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(WeightedScheduledSampler sampler) {
        samplers_.push_back(std::move(sampler));
        totalWeight_ += samplers_.back().weight;
        cumulativeWeightsCalculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::calculateCumulativeWeights() {
        std::sort(samplers_.begin(), samplers_.end(), std::greater<>());
        long double total = 0;
        cumulativeWeight_.clear();
        cumulativeWeight_.reserve(samplers_.size());
        for (auto& sampler : samplers_) {
            total += sampler.weight;
            cumulativeWeight_.push_back(total);
        }
        cumulativeWeightsCalculated_ = true;

        using bounds = boost::random::uniform_int_distribution<>::param_type;
        dist_.param(bounds(0, totalWeight_));
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::update(const WeightedScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.updateStart, sampler.updateEnd)) {
            sampler.update();
        }
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::adapt(const WeightedScheduledSampler& sampler) const {
        if (isBetween(totalSteps, sampler.adaptationStart, sampler.adaptationEnd)) {
            if (sampler.scaledAdaptation) {
                sampler.adapt(totalSteps - sampler.updateStart);
            } else {
                sampler.adapt();
            }
        }
    }

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
