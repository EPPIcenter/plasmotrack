//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_DISCRETERANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_DISCRETERANDOMWALK_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <utility>

#include "core/parameters/Parameter.h"

#include "core/samplers/AbstractSampler.h"


namespace transmission_nets::core::samplers {

    template<typename T, typename Engine = boost::random::mt19937>
    class DiscreteRandomWalk : public AbstractSampler {

    public:
        DiscreteRandomWalk(std::shared_ptr<parameters::Parameter<int>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        DiscreteRandomWalk(std::shared_ptr<parameters::Parameter<int>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, unsigned int maxDistance) noexcept;


        [[nodiscard]] unsigned int acceptances() const noexcept;

        [[nodiscard]] unsigned int rejections() const noexcept;

        [[nodiscard]] double acceptanceRate() const noexcept;

        virtual double sampleStride([[maybe_unused]] int current) noexcept;

        virtual Likelihood logMetropolisHastingsAdjustment([[maybe_unused]] int curr, [[maybe_unused]] int proposed) noexcept;

        void update() noexcept override;


    protected:
        std::shared_ptr<parameters::Parameter<int>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        unsigned int max_distance_ = 1;
        boost::random::normal_distribution<> normal_dist_{0, 1};
        boost::random::uniform_01<> uniform_dist_{};
        boost::random::uniform_int_distribution<> stride_sampling_dist_;

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;
    };

    template<typename T, typename Engine>
    DiscreteRandomWalk<T, Engine>::DiscreteRandomWalk(std::shared_ptr<parameters::Parameter<int>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng) {
        stride_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, max_distance_));
    }

    template<typename T, typename Engine>
    DiscreteRandomWalk<T, Engine>::DiscreteRandomWalk(std::shared_ptr<parameters::Parameter<int>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                      unsigned int maxDistance) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), max_distance_(maxDistance) {
        stride_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, max_distance_));
    }

    template<typename T, typename Engine>
    unsigned int DiscreteRandomWalk<T, Engine>::acceptances() const noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine>
    unsigned int DiscreteRandomWalk<T, Engine>::rejections() const noexcept {
        return rejections_;
    }

    template<typename T, typename Engine>
    double DiscreteRandomWalk<T, Engine>::acceptanceRate() const noexcept {
        return double(acceptances_) / double(rejections_ + acceptances_);
    }

    template<typename T, typename Engine>
    double DiscreteRandomWalk<T, Engine>::sampleStride([[maybe_unused]] int current) noexcept {
        auto stride = stride_sampling_dist_(*rng_);
        return uniform_dist_(*rng_) > .5 ? stride : -stride;
    }

    template<typename T, typename Engine>
    Likelihood DiscreteRandomWalk<T, Engine>::logMetropolisHastingsAdjustment(
            [[maybe_unused]] int curr, [[maybe_unused]] int proposed) noexcept {
        return 0;
    }

    template<typename T, typename Engine>
    void DiscreteRandomWalk<T, Engine>::update() noexcept {
        const std::string stateId = "State1";
        Likelihood curLik         = target_->value();
        parameter_->saveState(stateId);

        const int stride        = sampleStride(parameter_->value());
        Likelihood mhAdjustment = logMetropolisHastingsAdjustment(parameter_->value(), parameter_->value() + stride);

        parameter_->setValue(parameter_->value() + stride);

        const Likelihood acceptanceRatio = target_->value() - curLik + mhAdjustment;
        const bool accept                = log(uniform_dist_(*rng_)) <= acceptanceRatio;

        if (accept) {
            acceptances_ += 1;
            parameter_->acceptState();
        } else {
            rejections_ += 1;
            parameter_->restoreState(stateId);
            assert(curLik == target_->value());
        }

        total_updates_++;
    }

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_DISCRETERANDOMWALK_H
