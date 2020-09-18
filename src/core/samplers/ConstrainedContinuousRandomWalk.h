//
// Created by Maxwell Murphy on 3/25/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H

#include <cmath>

#include "ContinuousRandomWalk.h"

namespace transmission_nets::core::samplers {

    // Continuous Random Walk constrained to a range (lower, upper)

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine=boost::random::mt19937>
    class ConstrainedContinuousRandomWalk : public ContinuousRandomWalk<T, Engine> {
    public:
        ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, Engine *rng);

        ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, Engine *rng, double variance);

        ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, Engine *rng, double variance, double minVariance,
                                        double maxVariance);

    private:
        using ContinuousRandomWalk<T, Engine>::rng_;
        using ContinuousRandomWalk<T, Engine>::variance_;
        using ContinuousRandomWalk<T, Engine>::normal_dist_;
        using ContinuousRandomWalk<T, Engine>::parameter_;
        double logMetropolisHastingsAdjustment(double curr, double proposed) noexcept override;
        double sampleProposal() noexcept override;

    };

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    double ConstrainedContinuousRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::logMetropolisHastingsAdjustment(double curr, double proposed) noexcept {
        return log(proposed - LOWER_BOUND) +
               log(UPPER_BOUND - proposed) -
               log(curr - LOWER_BOUND) -
               log(UPPER_BOUND - curr);
    }

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    double ConstrainedContinuousRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::sampleProposal() noexcept {
        double eps = normal_dist_(*rng_) * variance_;
        double unconstrained = log(parameter_.value() - LOWER_BOUND) - log(UPPER_BOUND - parameter_.value());
        double exp_prop = exp(eps + unconstrained);
        double prop = (UPPER_BOUND * exp_prop + LOWER_BOUND) / (exp_prop + 1);
        assert(prop <= UPPER_BOUND);
        assert(prop >= LOWER_BOUND);
        return prop;
    }

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    ConstrainedContinuousRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter,
                                                                                                          T &target, Engine *rng)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng) {}


    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    ConstrainedContinuousRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter,
                                                                                                          T &target, Engine *rng,
                                                                                                          double variance)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance) {}


    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    ConstrainedContinuousRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::ConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter,
                                                                                                          T &target, Engine *rng,
                                                                                                          double variance,
                                                                                                          double minVariance,
                                                                                                          double maxVariance)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance, minVariance, maxVariance) {}

}

#endif //TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H
