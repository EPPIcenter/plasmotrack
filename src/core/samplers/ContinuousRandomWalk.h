//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H


#include <boost/random.hpp>
#include <boost/math/distributions.hpp>
#include <cmath>

#include "core/parameters/Parameter.h"

#include "core/samplers/AbstractSampler.h"


// Random Walk Metropolis Hastings using a gaussian proposal distribution centered at the current value

template<typename T, typename Engine=boost::random::mt19937>
class ContinuousRandomWalk : public AbstractSampler {

public:
    ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng) noexcept;

    ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng, double variance) noexcept;

    ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng, double variance, double minVariance,
                         double maxVariance) noexcept;


    [[nodiscard]] unsigned int acceptances() const noexcept;

    [[nodiscard]] unsigned int rejections() const noexcept;

    [[nodiscard]] double variance() const noexcept;

    [[nodiscard]] double acceptanceRate() const noexcept;

    void setTargetAcceptanceRate(double target) noexcept;

    void setAdaptationRate(double rate) noexcept;

    virtual double sampleProposal() noexcept;

    virtual double
    logMetropolisHastingsAdjustment([[maybe_unused]] double curr, [[maybe_unused]] double proposed) noexcept;

    void update() noexcept override;

    void adapt() noexcept override;

    void adapt(unsigned int idx) noexcept override;

protected:
    Parameter<double> &parameter_;
    T &target_;
    Engine *rng_;
    double variance_ = 1;
    double min_variance_ = 1e-12;
    double max_variance_ = 1e6;

    boost::random::normal_distribution<> normal_dist_{0, 1};
    boost::random::uniform_01<> uniform_dist_{};

    double adaptation_rate_ = .5;
    double target_acceptance_rate_ = .23;

    unsigned int acceptances_ = 0;
    unsigned int rejections_ = 0;
    unsigned int total_updates_ = 0;

};

template<typename T, typename Engine>
ContinuousRandomWalk<T, Engine>::ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng) noexcept
        : parameter_(parameter),
          target_(target), rng_(rng) {}

template<typename T, typename Engine>
ContinuousRandomWalk<T, Engine>::ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng,
                                                      double variance) noexcept : parameter_(
        parameter), target_(target), rng_(rng), variance_(variance) {
    assert(variance > 0);
}

template<typename T, typename Engine>
ContinuousRandomWalk<T, Engine>::ContinuousRandomWalk(Parameter<double> &parameter, T &target, Engine *rng,
                                                      double variance, double minVariance,
                                                      double maxVariance) noexcept :
        parameter_(parameter), target_(target), rng_(rng), variance_(variance), min_variance_(minVariance),
        max_variance_(maxVariance) {}

template<typename T, typename Engine>
void ContinuousRandomWalk<T, Engine>::update() noexcept {
    double curLik = target_.value();
    parameter_.saveState();

    const double currentVal = parameter_.value();
    const double proposal = sampleProposal();

    assert(!target_.isDirty());
    parameter_.setValue(proposal);

    assert(target_.isDirty());
    const double acceptanceRatio = target_.value() - curLik + logMetropolisHastingsAdjustment(currentVal, proposal);

    const bool accept = log(uniform_dist_(*rng_)) <= acceptanceRatio;

    if (accept) {
        acceptances_ += 1;
        parameter_.acceptState();
    } else {
        rejections_ += 1;
        parameter_.restoreState();
        assert(curLik == target_.value());
    }
    assert(!target_.isDirty());

    total_updates_++;
}

template<typename T, typename Engine>
void ContinuousRandomWalk<T, Engine>::adapt() noexcept {
    variance_ += (acceptanceRate() - target_acceptance_rate_) / std::pow(total_updates_ + 1, adaptation_rate_);
    variance_ = std::clamp(variance_, min_variance_, max_variance_);
}

template<typename T, typename Engine>
void ContinuousRandomWalk<T, Engine>::adapt(unsigned int idx) noexcept {
    variance_ += (acceptanceRate() - target_acceptance_rate_) / std::pow(idx, adaptation_rate_);
    variance_ = std::clamp(variance_, min_variance_, max_variance_);
}

template<typename T, typename Engine>
unsigned int ContinuousRandomWalk<T, Engine>::acceptances() const noexcept {
    return acceptances_;
}

template<typename T, typename Engine>
unsigned int ContinuousRandomWalk<T, Engine>::rejections() const noexcept {
    return rejections_;
}

template<typename T, typename Engine>
double ContinuousRandomWalk<T, Engine>::acceptanceRate() const noexcept {
    return double(acceptances_) / double(rejections_ + acceptances_);
}

template<typename T, typename Engine>
double ContinuousRandomWalk<T, Engine>::sampleProposal() noexcept {
    return parameter_.value() + normal_dist_(*rng_) * variance_;
}

template<typename T, typename Engine>
double
ContinuousRandomWalk<T, Engine>::logMetropolisHastingsAdjustment([[maybe_unused]] double curr,
                                                                 [[maybe_unused]] double proposed) noexcept {
    return 0;
}

template<typename T, typename Engine>
double ContinuousRandomWalk<T, Engine>::variance() const noexcept {
    return variance_;
}

template<typename T, typename Engine>
void ContinuousRandomWalk<T, Engine>::setTargetAcceptanceRate(double target) noexcept {
    target_acceptance_rate_ = target;
}

template<typename T, typename Engine>
void ContinuousRandomWalk<T, Engine>::setAdaptationRate(double rate) noexcept {
    adaptation_rate_ = rate;
}

#endif //TRANSMISSION_NETWORKS_APP_CONTINUOUSRANDOMWALK_H
