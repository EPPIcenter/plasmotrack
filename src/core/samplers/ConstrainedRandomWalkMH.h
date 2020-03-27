//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H

#include <cmath>

#include "RandomWalkMH.h"

template<typename T>
class ConstrainedRandomWalkMH : public RandomWalkMH<T> {
public:
    ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng, double lowerBound,
                            double upperBound);

    ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng, double variance, double lowerBound,
                            double upperBound);

private:
    double logMetropolisHastingsAdjustment(double curr, double proposed) override;

    double sampleProposal() override;

private:
    double lower_bound_;
    double upper_bound_;
};

template<typename T>
ConstrainedRandomWalkMH<T>::ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng,
                                                    double lowerBound, double upperBound)
        : RandomWalkMH<T>(parameter, target, rng), lower_bound_(lowerBound), upper_bound_(upperBound) {}

template<typename T>
ConstrainedRandomWalkMH<T>::ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, gsl_rng *rng,
                                                    double variance, double lowerBound, double upperBound)
        :RandomWalkMH<T>(
        parameter, target, rng, variance), lower_bound_(lowerBound), upper_bound_(upperBound) {}


template<typename T>
double ConstrainedRandomWalkMH<T>::logMetropolisHastingsAdjustment(double curr, double proposed) {
    return log(proposed - lower_bound_) +
           log(upper_bound_ - proposed) -
           log(curr - lower_bound_) -
           log(upper_bound_ - curr);
}

template<typename T>
double ConstrainedRandomWalkMH<T>::sampleProposal() {
    double eps = gsl_ran_gaussian_ziggurat(this->rng, this->variance_);
    double unconstrained = log(this->parameter_.value() - lower_bound_) - log(upper_bound_ - this->parameter_.value());
    double exp_prop = exp(eps + unconstrained);
    double prop = (upper_bound_ * exp_prop + lower_bound_) / (exp_prop + 1);
    assert(prop <= upper_bound_);
    assert(prop >= lower_bound_);
    return prop;
}

#endif //TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H
