//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H


#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <cmath>
#include <utility>

#include "core/parameters/Parameter.h"
#include "core/samplers/AbstractSampler.h"
#include "core/samplers/SamplerEnum.h"


// Random Walk Metropolis Hastings using a gaussian proposal distribution centered at the current value

namespace transmission_nets::core::samplers {
    template<typename T, typename Engine = boost::random::mt19937, typename U = float>
    class ContinuousRandomWalk : public AbstractSampler {

    public:
        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, float variance) noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, float variance, float minVariance,
                             float maxVariance) noexcept;


        [[nodiscard]] unsigned int acceptances() const noexcept;

        [[nodiscard]] unsigned int rejections() const noexcept;

        [[nodiscard]] float variance() const noexcept;

        [[nodiscard]] float acceptanceRate() const noexcept;

        void setTargetAcceptanceRate(float target) noexcept;

        void setAdaptationRate(float rate) noexcept;

        virtual float sampleProposal() noexcept;

        [[nodiscard]] virtual Likelihood
        logMetropolisHastingsAdjustment([[maybe_unused]] U curr, [[maybe_unused]] U proposed) const noexcept;

        void update() noexcept override;

        void adapt() noexcept override;

        void adapt(unsigned int idx) noexcept override;


    protected:
        std::shared_ptr<parameters::Parameter<float>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        float variance_     = 1;
        float min_variance_ = 1e-12;
        float max_variance_ = 1e6;

        boost::random::normal_distribution<> normal_dist_{0, 1};
        boost::random::uniform_01<> uniform_dist_{};

        float adaptation_rate_        = .5;
        float target_acceptance_rate_ = .23;

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;

        float timeSpentSampling_ = 0;


    };

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng) noexcept
        : parameter_(std::move(parameter)),
          target_(target), rng_(rng) {}

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             float variance) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance) {
        assert(variance > 0);
    }

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<float>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             float variance, float minVariance,
                                                             float maxVariance) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance), min_variance_(minVariance),
                                                                                                                    max_variance_(maxVariance) {}

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::update() noexcept {
        SAMPLER_STATE_ID stateId = SAMPLER_STATE_ID::ConstrainedContinuousRandomWalkID;
        Likelihood curLik         = target_->value();
        parameter_->saveState(stateId);

        const float currentVal = parameter_->value();
        const float proposal   = sampleProposal();

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const Likelihood adj = logMetropolisHastingsAdjustment(currentVal, proposal);

        if (debug_) {
            std::cout << "ID:" << identifier_ << std::endl;
            std::cout << "Curr: " << currentVal << std::endl;
            std::cout << "Prop: " << proposal << std::endl;
            std::cout << "Var: " << variance() << std::endl;
            std::cout << "Acceptance Rate: " << acceptanceRate() << std::endl;
            std::cout << "Total Samples: " << acceptances_ + rejections_ << std::endl;
//            std::cout << "LLik: " << curLik << " " << target_->value() << " " << adj << " " << acceptanceRatio << std::endl;
            std::cout << "LLik: " << curLik << " " << target_->value() << " " << adj << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
        }

        const Likelihood acceptanceRatio = target_->value() - curLik + adj;
        assert(!target_->isDirty());
//        if (debug_) {
//            std::cout << "ID:" << identifier_ << std::endl;
//            std::cout << "Curr: " << currentVal << std::endl;
//            std::cout << "Prop: " << proposal << std::endl;
//            std::cout << "Var: " << variance() << std::endl;
//            std::cout << "LLik: " << curLik << " " << target_->value() << " " << adj << " " << acceptanceRatio << std::endl;
//            std::cout << "----------------------------------------------------" << std::endl;
//        }
        const bool accept = log(uniform_dist_(*rng_)) <= acceptanceRatio;

        if (accept) {
            acceptances_ += 1;
            parameter_->acceptState();
        } else {
            rejections_ += 1;
            parameter_->restoreState(stateId);
            assert(curLik == target_->value());
        }
        assert(!std::isinf(target_->value()));
        assert(!target_->isDirty());

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
    void ContinuousRandomWalk<T, Engine, U>::adapt(unsigned int idx) noexcept {
        variance_ += (acceptanceRate() - target_acceptance_rate_) / std::pow(idx, adaptation_rate_);
        variance_ = std::clamp(variance_, min_variance_, max_variance_);
        if (std::isnan(variance_)) {
            variance_ = min_variance_;
        }
    }

    template<typename T, typename Engine, typename U>
    unsigned int ContinuousRandomWalk<T, Engine, U>::acceptances() const noexcept {
        return acceptances_;
    }

    template<typename T, typename Engine, typename U>
    unsigned int ContinuousRandomWalk<T, Engine, U>::rejections() const noexcept {
        return rejections_;
    }

    template<typename T, typename Engine, typename U>
    float ContinuousRandomWalk<T, Engine, U>::acceptanceRate() const noexcept {
        return float(acceptances_) / float(rejections_ + acceptances_);
    }

    template<typename T, typename Engine, typename U>
    float ContinuousRandomWalk<T, Engine, U>::sampleProposal() noexcept {
        return parameter_->value() + normal_dist_(*rng_) * variance_;
    }

    template<typename T, typename Engine, typename U>
    Likelihood
    ContinuousRandomWalk<T, Engine, U>::logMetropolisHastingsAdjustment([[maybe_unused]] U curr,
                                                                        [[maybe_unused]] U proposed) const noexcept {
        return 0;
    }

    template<typename T, typename Engine, typename U>
    float ContinuousRandomWalk<T, Engine, U>::variance() const noexcept {
        return variance_;
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::setTargetAcceptanceRate(float target) noexcept {
        target_acceptance_rate_ = target;
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::setAdaptationRate(float rate) noexcept {
        adaptation_rate_ = rate;
    }
}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
