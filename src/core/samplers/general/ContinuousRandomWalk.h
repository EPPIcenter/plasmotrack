//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H


#include <boost/random.hpp>
#include <cmath>
#include <utility>

#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
// Random Walk Metropolis Hastings using a gaussian proposal distribution centered at the current value

namespace transmission_nets::core::samplers {
    template<typename T, typename Engine = boost::random::mt19937, typename U = double>
    class ContinuousRandomWalk : public AbstractSampler {

    public:
        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance) noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, double minVariance,
                             double maxVariance) noexcept;


        [[nodiscard]] double variance() const noexcept;

        [[nodiscard]] double acceptanceRate() const noexcept;

        void setAdaptationRate(double rate) noexcept;

        virtual double sampleProposal() noexcept;

        [[nodiscard]] virtual Likelihood
        logMetropolisHastingsAdjustment([[maybe_unused]] U curr, [[maybe_unused]] U proposed) const noexcept;

        void update() noexcept override;

        void adapt() noexcept override;

        [[nodiscard]] AdaptationState adaptationState() const override;

        void restoreAdaptationState(const AdaptationState& state) override;

    protected:
        std::shared_ptr<parameters::Parameter<double>> parameter_{};
        std::shared_ptr<T> target_{};
        std::shared_ptr<Engine> rng_{};
        double variance_     = 1;
        double min_variance_ = 1e-12;
        double max_variance_ = 1e6;

        boost::random::normal_distribution<> normal_dist_{0, 1};
        double adaptation_rate_        = .5;
        double target_acceptance_rate_ = .23;
        unsigned int total_updates_ = 0;


    };

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept
        : parameter_(std::move(parameter)),
          target_(target), rng_(rng) {}

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             double variance) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance) {
        assert(variance > 0);
    }

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             double variance, double minVariance,
                                                             double maxVariance) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance), min_variance_(minVariance),
                                                                                                                    max_variance_(maxVariance) {}

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::update() noexcept {
        auto proposal           = this->propose(*target_, *rng_);
        const double currentVal = parameter_->value();
        const double proposed   = sampleProposal();
        proposal.change(*parameter_, proposed);
        proposal.metropolis(logMetropolisHastingsAdjustment(currentVal, proposed));
        total_updates_++;
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::adapt() noexcept {
        variance_ += (acceptanceRate() - target_acceptance_rate_) / std::pow(total_updates_ + 1, adaptation_rate_);
        variance_ = std::clamp(variance_, min_variance_, max_variance_);
        if (std::isnan(variance_)) {
            variance_ = min_variance_;
        }
    }

    template<typename T, typename Engine, typename U>
    AdaptationState ContinuousRandomWalk<T, Engine, U>::adaptationState() const {
        return {.scales  = {{.variance = variance_, .counts = this->counts()}},
                .updates = total_updates_};
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::restoreAdaptationState(const AdaptationState& state) {
        requireScales(state, 1);
        const auto& scale = state.scales.front();
        variance_         = scale.variance;
        total_updates_    = static_cast<unsigned int>(state.updates);
        // The walk's one scale counts every proposal it makes, so its counts are the sampler's.
        this->restoreCounts(scale.counts);
    }

    template<typename T, typename Engine, typename U>
    double ContinuousRandomWalk<T, Engine, U>::acceptanceRate() const noexcept {
        const auto counts = this->counts();
        return static_cast<double>(counts.accepted) / static_cast<double>(counts.proposed);
    }

    template<typename T, typename Engine, typename U>
    double ContinuousRandomWalk<T, Engine, U>::sampleProposal() noexcept {
        return parameter_->value() + normal_dist_(*rng_) * variance_;
    }

    template<typename T, typename Engine, typename U>
    Likelihood
    ContinuousRandomWalk<T, Engine, U>::logMetropolisHastingsAdjustment([[maybe_unused]] U curr,
                                                                        [[maybe_unused]] U proposed) const noexcept {
        return 0;
    }

    template<typename T, typename Engine, typename U>
    double ContinuousRandomWalk<T, Engine, U>::variance() const noexcept {
        return variance_;
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::setAdaptationRate(double rate) noexcept {
        adaptation_rate_ = rate;
    }
}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
