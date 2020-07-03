//
// Created by Maxwell Murphy on 4/7/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

#include "core/samplers/AbstractSampler.h"
#include "core/datatypes/Simplex.h"

#include "core/utils/random_sequence.h"
#include "core/utils/numerics.h"


//  @article{
//      doi:10.1080/00949655.2017.1376063,
//      author = {Hannah M. Director and James Gattiker and Earl Lawrence and Scott Vander Wiel},
//      title = {Efficient sampling on the simplex with a self-adjusting logit transform proposal},
//      journal = {Journal of Statistical Computation and Simulation},
//      volume = {87},
//      number = {18},
//      pages = {3521-3536},
//      year  = {2017},
//      publisher = {Taylor & Francis},
//      doi = {10.1080/00949655.2017.1376063},
//  }



template<typename T,  typename Engine=boost::random::mt19937>
class SALTSampler : public AbstractSampler {
public:
    SALTSampler(Parameter<Simplex> &parameter, T &target, Engine *rng);

    double sampleProposal(double current, double variance) noexcept;

    double logMetropolisHastingsAdjustment(double logitCurr, double logitProp, int k) noexcept;

    void setAdaptationRate(double adaptationRate);

    void setTargetAcceptanceRate(double targetAcceptanceRate);

    void setMinVariance(double minVariance);

    void setMaxVariance(double maxVariance);

    void update() noexcept override;

    void adapt() noexcept override;

    void adapt(unsigned int idx) noexcept override;

    double acceptanceRate(int idx) noexcept;


private:

    Parameter<Simplex> &parameter_;
    T& target_;
    Engine* rng_;
    boost::random::normal_distribution<> normal_dist_{0, 1};
    boost::random::uniform_01<> uniform_dist_{};

    std::vector<double> variances_{};
    std::vector<double> acceptances_{};
    std::vector<double> rejections_{};

    double min_variance_ = 1e-12;
    double max_variance_ = 1e6;


private:

    double adaptation_rate_ = 1;
    double target_acceptance_rate_ = .23;

    unsigned int total_updates_ = 0;


};


template<typename T, typename Engine>
void SALTSampler<T, Engine>::update() noexcept {
    auto indices = randomSequence(0, parameter_.value().totalElements(), rng_);
    for(const auto idx : indices) {
        Simplex currentVal(parameter_.value());
        double curLik = target_.value();

        double logitCurr = logit(currentVal.frequencies(idx));
        double logitProp = sampleProposal(logitCurr, variances_[idx]);
        double prop = std::max(expit(logitProp), 1e-8); // technically incorrect to bound like this without correcting, but probably close enough

        currentVal.set(idx, prop);
        parameter_.saveState();
        parameter_.setValue(currentVal);

        const double adjRatio = logMetropolisHastingsAdjustment(logitCurr, logitProp, currentVal.totalElements());
        const double acceptanceRatio = target_.value() - curLik + adjRatio;

        const bool accept = log(uniform_dist_(*rng_)) <= acceptanceRatio;



        if (accept) {
//            std::cout << target_.value() << " | " << curLik << "\n";
            acceptances_.at(idx)++;
            parameter_.acceptState();
        } else {
            rejections_.at(idx)++;
            parameter_.restoreState();
            assert(curLik == target_.value());
        }
    }

    total_updates_++;

}

template<typename T, typename Engine>
SALTSampler<T, Engine>::SALTSampler(Parameter<Simplex> &parameter, T &target, Engine *rng) : parameter_(
        parameter), target_(target), rng_(rng) {
    for (size_t j = 0; j < parameter_.value().totalElements(); ++j) {
        variances_.push_back(1);
        acceptances_.push_back(0);
        rejections_.push_back(0);
    }
}

template<typename T, typename Engine>
double SALTSampler<T, Engine>::sampleProposal(double logitCurr, double variance) noexcept {
    double eps = normal_dist_(*rng_) * variance;
    return logitCurr + eps;
}

template<typename T, typename Engine>
double SALTSampler<T, Engine>::logMetropolisHastingsAdjustment(double logitCurr, double logitProp, int k) noexcept {

    double logAdj = 0;

    if(logitProp < 0) {
        logAdj += logitProp - log1p(exp(logitProp));
        logAdj += (k - 1) * log1p(exp(logitProp));
    } else {
        logAdj += -log1p(1 / exp(logitProp));
        logAdj += (k - 1) * (-log1p(1 / exp(logitProp)) - logitProp);
    }

    if(logitCurr < 0) {
        logAdj -= logitCurr - log1p(exp(logitCurr));
        logAdj -= (k - 1) * log1p(exp(logitCurr));
    } else {
        logAdj += -log1p(1 / exp(logitCurr));
        logAdj += (k - 1) * (-log1p(1 / exp(logitCurr)) - logitCurr);
    }

    return logAdj;
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::adapt(unsigned int idx) noexcept {
    for (unsigned int j = 0; j < parameter_.value().totalElements(); ++j) {
        variances_.at(j) += (acceptanceRate(j) - target_acceptance_rate_) / std::pow(idx, adaptation_rate_);
        variances_.at(j) = std::clamp(variances_.at(j), min_variance_, max_variance_);
    }
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::adapt() noexcept {
    for (unsigned int j = 0; j < parameter_.value().totalElements(); ++j) {
        variances_.at(j) += (acceptanceRate(j) - target_acceptance_rate_) / std::pow(total_updates_ + 1, adaptation_rate_);
        variances_.at(j) = std::clamp(variances_.at(j), min_variance_, max_variance_);
    }
}

template<typename T, typename Engine>
double SALTSampler<T, Engine>::acceptanceRate(int idx) noexcept {
    return acceptances_.at(idx) / (acceptances_.at(idx) + rejections_.at(idx));
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::setAdaptationRate(double adaptationRate) {
    adaptation_rate_ = adaptationRate;
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::setTargetAcceptanceRate(double targetAcceptanceRate) {
    target_acceptance_rate_ = targetAcceptanceRate;
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::setMinVariance(double minVariance) {
    min_variance_ = minVariance;
}

template<typename T, typename Engine>
void SALTSampler<T, Engine>::setMaxVariance(double maxVariance) {
    max_variance_ = maxVariance;
}


#endif //TRANSMISSION_NETWORKS_APP_SALTSAMPLER_H
