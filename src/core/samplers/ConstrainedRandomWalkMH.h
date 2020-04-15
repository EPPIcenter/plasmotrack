//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H

#include <cmath>

#include "RandomWalkMH.h"

// Random Walk constrained to a range (lower, upper)

template<typename T, typename Engine>
class ConstrainedRandomWalkMH : public RandomWalkMH<T, Engine> {
public:
    ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng, double lowerBound,
                            double upperBound) noexcept;

    ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng, double variance, double lowerBound,
                            double upperBound) noexcept;

    ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng, double variance, double minVariance,
                            double maxVariance, double lowerBound, double upperBound);

private:
    using RandomWalkMH<T, Engine>::rng_;
    using RandomWalkMH<T, Engine>::variance_;
    using RandomWalkMH<T, Engine>::normal_dist_;
    using RandomWalkMH<T, Engine>::parameter_;
    double logMetropolisHastingsAdjustment(double curr, double proposed) noexcept override;
    double sampleProposal() noexcept override;
    double lower_bound_;
    double upper_bound_;
};

template<typename T, typename Engine>
ConstrainedRandomWalkMH<T, Engine>::ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng,
                                                    double lowerBound, double upperBound) noexcept
        : RandomWalkMH<T, Engine>(parameter, target, rng), lower_bound_(lowerBound), upper_bound_(upperBound) {}

template<typename T, typename Engine>
ConstrainedRandomWalkMH<T, Engine>::ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng,
                                                    double variance, double lowerBound, double upperBound) noexcept
        :RandomWalkMH<T, Engine>(
        parameter, target, rng, variance), lower_bound_(lowerBound), upper_bound_(upperBound) {}


template<typename T, typename Engine>
ConstrainedRandomWalkMH<T, Engine>::ConstrainedRandomWalkMH(Parameter<double> &parameter, T &target, Engine *rng,
                                                            double variance, double minVariance, double maxVariance,
                                                            double lowerBound, double upperBound):RandomWalkMH<T, Engine>(
        parameter, target, rng, variance, minVariance, maxVariance), lower_bound_(lowerBound), upper_bound_(
        upperBound) {}

template<typename T, typename Engine>
double ConstrainedRandomWalkMH<T, Engine>::logMetropolisHastingsAdjustment(double curr, double proposed) noexcept {
    return log(proposed - lower_bound_) +
           log(upper_bound_ - proposed) -
           log(curr - lower_bound_) -
           log(upper_bound_ - curr);
}

template<typename T, typename Engine>
double ConstrainedRandomWalkMH<T, Engine>::sampleProposal() noexcept {
    double eps = normal_dist_(*rng_) * variance_;
    double unconstrained = log(parameter_.value() - lower_bound_) - log(upper_bound_ - parameter_.value());
    double exp_prop = exp(eps + unconstrained);
    double prop = (upper_bound_ * exp_prop + lower_bound_) / (exp_prop + 1);
    assert(prop <= upper_bound_);
    assert(prop >= lower_bound_);
    return prop;
}

#endif //TRANSMISSION_NETWORKS_APP_CONSTRAINEDRANDOMWALKMH_H
