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


// Random Walk Metropolis Hastings using a gaussian proposal distribution centered at the current value

namespace transmission_nets::core::samplers {
    template<typename T, typename Engine = boost::random::mt19937, typename U = double>
    class ContinuousRandomWalk : public AbstractSampler {

    public:
        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, std::string identifier = "") noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, std::string identifier = "") noexcept;

        ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, double variance, double minVariance,
                             double maxVariance, std::string identifier = "") noexcept;


        [[nodiscard]] unsigned int acceptances() const noexcept;

        [[nodiscard]] unsigned int rejections() const noexcept;

        [[nodiscard]] double variance() const noexcept;

        [[nodiscard]] double acceptanceRate() const noexcept;

        void setTargetAcceptanceRate(double target) noexcept;

        void setAdaptationRate(double rate) noexcept;

        virtual double sampleProposal() noexcept;

        [[nodiscard]] virtual Likelihood
        logMetropolisHastingsAdjustment([[maybe_unused]] U curr, [[maybe_unused]] U proposed) const noexcept;

        void update() noexcept override;

        void adapt() noexcept override;

        void adapt(unsigned int idx) noexcept override;

        void setVerbose(bool verbose = true) {
            verbose_ = verbose;
        }


    protected:
        std::shared_ptr<parameters::Parameter<double>> parameter_;
        std::shared_ptr<T> target_;
        std::shared_ptr<Engine> rng_;
        double variance_     = 1;
        double min_variance_ = 1e-12;
        double max_variance_ = 1e6;

        boost::random::normal_distribution<> normal_dist_{0, 1};
        boost::random::uniform_01<> uniform_dist_{};

        double adaptation_rate_        = .5;
        double target_acceptance_rate_ = .23;

        unsigned int acceptances_   = 0;
        unsigned int rejections_    = 0;
        unsigned int total_updates_ = 0;

        bool verbose_ = false;
        std::string identifier_;
    };

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng, std::string identifier) noexcept
        : parameter_(std::move(parameter)),
          target_(target), rng_(rng), identifier_(std::move(identifier)) {}

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             double variance, std::string identifier) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance), identifier_(std::move(identifier)) {
        assert(variance > 0);
    }

    template<typename T, typename Engine, typename U>
    ContinuousRandomWalk<T, Engine, U>::ContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, std::shared_ptr<Engine> rng,
                                                             double variance, double minVariance,
                                                             double maxVariance, std::string identifier) noexcept : parameter_(std::move(parameter)), target_(target), rng_(rng), variance_(variance), min_variance_(minVariance),
                                                                                                                    max_variance_(maxVariance), identifier_(std::move(identifier)) {}

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::update() noexcept {
        const std::string stateId = "ContinuousRW";
        Likelihood curLik         = target_->value();
        parameter_->saveState(stateId);

        const double currentVal = parameter_->value();
        const double proposal   = sampleProposal();

        assert(!target_->isDirty());
        parameter_->setValue(proposal);
        const Likelihood adj = logMetropolisHastingsAdjustment(currentVal, proposal);


        const Likelihood acceptanceRatio = target_->value() - curLik + adj;
        assert(!target_->isDirty());
        if (verbose_) {
            std::cout << "ID:" << identifier_ << std::endl;
            std::cout << "Curr: " << currentVal << std::endl;
            std::cout << "Prop: " << proposal << std::endl;
            std::cout << "Var: " << variance() << std::endl;
            std::cout << "LLik: " << curLik << " " << target_->value() << " " << adj << " " << acceptanceRatio << std::endl;
            std::cout << "----------------------------------------------------" << std::endl;
        }
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
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::adapt(unsigned int idx) noexcept {
        variance_ += (acceptanceRate() - target_acceptance_rate_) / std::pow(idx, adaptation_rate_);
        variance_ = std::clamp(variance_, min_variance_, max_variance_);
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
    double ContinuousRandomWalk<T, Engine, U>::acceptanceRate() const noexcept {
        return double(acceptances_) / double(rejections_ + acceptances_);
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
    void ContinuousRandomWalk<T, Engine, U>::setTargetAcceptanceRate(double target) noexcept {
        target_acceptance_rate_ = target;
    }

    template<typename T, typename Engine, typename U>
    void ContinuousRandomWalk<T, Engine, U>::setAdaptationRate(double rate) noexcept {
        adaptation_rate_ = rate;
    }
}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
