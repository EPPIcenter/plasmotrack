//
// Created by Maxwell Murphy on 5/28/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
#define TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H

#include <boost/random.hpp>

#include "Scheduler.h"

namespace transmission_nets::core::samplers {

    struct WeightedScheduledSampler {
        AbstractSampler* sampler{};

        int adaptationStart = 0;
        int adaptationEnd = 0;
        bool scaledAdaptation{false};

        int updateStart = 0;
        int updateEnd = std::numeric_limits<int>::max();

        long double weight = 1.0;

        bool operator<(const WeightedScheduledSampler &rhs) const {
            return weight < rhs.weight;
        }

        bool operator>(const WeightedScheduledSampler &rhs) const {
            return rhs < *this;
        }

        bool operator<=(const WeightedScheduledSampler &rhs) const {
            return !(rhs < *this);
        }

        bool operator>=(const WeightedScheduledSampler &rhs) const {
            return !(*this < rhs);
        }


        void update() const {
            sampler->update();
        };

        void adapt() const {
            sampler->adapt();
        };

        void adapt(int step) const {
            sampler->adapt(step);
        };
    };


    template<typename Engine=boost::random::mt19937>
    class RandomizedScheduler {

    public:
        explicit RandomizedScheduler(Engine* rng, int numSamples);

        void registerSampler(AbstractSampler* sampler);

        void registerSampler(WeightedScheduledSampler sampler);

        void step();

        void update(const WeightedScheduledSampler& sampler) const;

        void adapt(const WeightedScheduledSampler& sampler) const;

    private:
        Engine* rng_;
        int numSamples_;

        std::vector<WeightedScheduledSampler> samplers_{};
        std::vector<long double> cumulativeWeight_{};
        int totalSteps = 0;
        long double totalWeight_ = 0;
        boost::random::uniform_int_distribution<> dist_;

        bool cumulativeWeightsCalculated_{false};
        void calculateCumulativeWeights();
    };

    template<typename Engine>
    RandomizedScheduler<Engine>::RandomizedScheduler(Engine* rng, int numSamples) : rng_(rng), numSamples_(numSamples) {}

    template<typename Engine>
    void RandomizedScheduler<Engine>::step() {
        if (!cumulativeWeightsCalculated_) {
            calculateCumulativeWeights();
        }

        long int v;
        for (int i = 0; i < numSamples_; ++i) {
            v = dist_(*rng_);
            auto it = std::lower_bound(cumulativeWeight_.begin(), cumulativeWeight_.end(), v);
            int idx = it - cumulativeWeight_.begin();
            update(samplers_[idx]);
            adapt(samplers_[idx]);
        }
        ++totalSteps;


        //    auto indices = randomSequence(0, samplers_.size(), rng_);
        //    for (const auto idx : indices) {
        //        update(samplers_[idx]);
        //        adapt(samplers_[idx]);
        //    }
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(AbstractSampler* sampler) {
        samplers_.push_back(WeightedScheduledSampler{.sampler = sampler});
        totalWeight_ += samplers_.back().weight;
        cumulativeWeightsCalculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::registerSampler(WeightedScheduledSampler sampler) {
        samplers_.push_back(sampler);
        totalWeight_ += samplers_.back().weight;
        cumulativeWeightsCalculated_ = false;
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::calculateCumulativeWeights() {
        std::sort(samplers_.begin(), samplers_.end(), std::greater<>());
        long double cumulant = 0;
        cumulativeWeight_.clear();
        cumulativeWeight_.reserve(samplers_.size());
        for (auto & sampler : samplers_) {
            cumulativeWeight_.push_back(cumulant + sampler.weight);
            cumulant += sampler.weight;
        }
        //    std::function<long double(const WeightedScheduledSampler&, const WeightedScheduledSampler&)> binOp = [](auto left, auto right) {return left.weight + right.weight};
        //    std::partial_sum(
        //            samplers_.begin(),
        //            samplers_.end(),
        //            cumulativeWeight_.begin(),
        //            [](const WeightedScheduledSampler& left, const WeightedScheduledSampler& right) -> long double
        //            {
        //                return left.weight + right.weight;
        //            });
        cumulativeWeightsCalculated_ = true;

        using bounds = boost::random::uniform_int_distribution<>::param_type;
        dist_.param(bounds(0, totalWeight_));
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::update(const WeightedScheduledSampler &sampler) const {
        if (isBetween(totalSteps, sampler.updateStart, sampler.updateEnd)) {
            sampler.update();
        }
    }

    template<typename Engine>
    void RandomizedScheduler<Engine>::adapt(const WeightedScheduledSampler &sampler) const {
        if (isBetween(totalSteps, sampler.adaptationStart, sampler.adaptationEnd)) {
            if (sampler.scaledAdaptation) {
                sampler.adapt(totalSteps - sampler.updateStart);
            } else {
                sampler.adapt();
            }
        }
    }

}





#endif //TRANSMISSION_NETWORKS_APP_RANDOMIZEDSCHEDULER_H
