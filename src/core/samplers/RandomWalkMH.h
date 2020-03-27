//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_RANDOMWALKMH_H
#define TRANSMISSION_NETWORKS_APP_RANDOMWALKMH_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <cmath>

#include "AbstractSampler.h"


// Random Walk Metropolis Hastings using a gaussian proposal distribution centered at the current value

template<typename T>
class RandomWalkMH : AbstractSampler {

public:
    RandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng);

    RandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng, double variance);

    [[nodiscard]] unsigned int acceptances() const;

    [[nodiscard]] unsigned int rejections() const;

    [[nodiscard]] double variance() const;

    [[nodiscard]] double acceptanceRate() const;

    void setTargetAcceptanceRate(double target);

    virtual double sampleProposal();

    virtual double logMetropolisHastingsAdjustment([[maybe_unused]] double curr, [[maybe_unused]] double proposed);

    void update() override;

    void adapt();

protected:
    Parameter<double> &parameter_;
    T &target_;
    gsl_rng *rng;
    double variance_ = 1;
//    unsigned int max_adaptation_steps_ = 0;
    double min_variance_ = 1e-6;
    double max_variance_ = 1e6;
    double target_acceptance_rate_ = .23;

//    bool adaptive_ = false;
    unsigned int acceptances_ = 0;
    unsigned int rejections_ = 0;
    unsigned int total_updates_ = 0;

};

template<typename T>
void RandomWalkMH<T>::update() {
    double curLik = target_.value();
    parameter_.saveState();

    const double currentVal = parameter_.value();
    const double proposal = sampleProposal();

    parameter_.setValue(proposal);
    const double acceptanceRatio = target_.value() - curLik + logMetropolisHastingsAdjustment(currentVal, proposal);

    const bool accept = log(gsl_ran_flat(rng, 0, 1)) <= acceptanceRatio;

    if (accept) {
        acceptances_ += 1;
        parameter_.acceptState();
    } else {
        rejections_ += 1;
        parameter_.restoreState();
    }

    total_updates_++;
}

template<typename T>
void RandomWalkMH<T>::adapt() {
    variance_ += (acceptanceRate() - target_acceptance_rate_) / std::sqrt(total_updates_ + 1);
    variance_ = std::clamp(variance_, min_variance_, max_variance_);
}

template<typename T>
RandomWalkMH<T>::RandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng) : parameter_(parameter),
                                                                                       target_(target), rng(rng) {}

template<typename T>
RandomWalkMH<T>::RandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng, double variance) : parameter_(
        parameter), target_(target), rng(rng), variance_(variance) {
    assert(variance > 0);
}

template<typename T>
unsigned int RandomWalkMH<T>::acceptances() const {
    return acceptances_;
}

template<typename T>
unsigned int RandomWalkMH<T>::rejections() const {
    return rejections_;
}

template<typename T>
double RandomWalkMH<T>::acceptanceRate() const {
    return double(acceptances_) / double(rejections_ + acceptances_);
}

template<typename T>
double RandomWalkMH<T>::sampleProposal() {
    return parameter_.value() + gsl_ran_gaussian_ziggurat(rng, variance_);
}

template<typename T>
double
RandomWalkMH<T>::logMetropolisHastingsAdjustment([[maybe_unused]] double curr, [[maybe_unused]] double proposed) {
    return 0;
}

template<typename T>
double RandomWalkMH<T>::variance() const {
    return variance_;
}

template<typename T>
void RandomWalkMH<T>::setTargetAcceptanceRate(double target) {
    target_acceptance_rate_ = target;
}


#endif //TRANSMISSION_NETWORKS_APP_RANDOMWALKMH_H
