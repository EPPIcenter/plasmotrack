//
// Created by Maxwell Murphy on 7/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_DOUBLECONSTRAINEDCONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_DOUBLECONSTRAINEDCONTINUOUSRANDOMWALK_H

#include <cmath>
#include <algorithm>

#include "ContinuousRandomWalk.h"

namespace transmission_nets::core::samplers {

    // Continuous Random Walk constrained to a range (lower, upper)

    template<typename T, typename Engine=boost::random::mt19937>
    class DoubleConstrainedContinuousRandomWalk : public ContinuousRandomWalk<T, Engine> {
    public:
        DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, double lower_bound, double upper_bound, Engine *rng);

        DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target,  double lower_bound, double upper_bound, Engine *rng, double variance);

        DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, double lower_bound, double upper_bound, Engine *rng, double variance, double minVariance, double maxVariance);

    private:
        using ContinuousRandomWalk<T, Engine>::rng_;
        using ContinuousRandomWalk<T, Engine>::variance_;
        using ContinuousRandomWalk<T, Engine>::normal_dist_;
        using ContinuousRandomWalk<T, Engine>::parameter_;
        double lower_bound_;
        double upper_bound_;
        double logMetropolisHastingsAdjustment(double curr, double proposed) noexcept override;
        double sampleProposal() noexcept override;

    };


    template<typename T, typename Engine>
    double DoubleConstrainedContinuousRandomWalk<T, Engine>::logMetropolisHastingsAdjustment(double curr, double proposed) noexcept {
        return log(proposed - lower_bound_) +
               log(upper_bound_ - proposed) -
               log(curr - lower_bound_) -
               log(upper_bound_ - curr);
    }

    template<typename T, typename Engine>
    double DoubleConstrainedContinuousRandomWalk<T, Engine>::sampleProposal() noexcept {
        double eps = normal_dist_(*rng_) * variance_;
        double unconstrained = log(parameter_.value() - lower_bound_) - log(upper_bound_ - parameter_.value());
        double exp_prop = exp(eps + unconstrained);
        double prop = (upper_bound_ * exp_prop + lower_bound_) / (exp_prop + 1);
        prop = std::clamp(prop, lower_bound_, upper_bound_);
#ifndef NDEBUG
        if(prop > upper_bound_) {
            std::cerr << "Proposal " << prop << " exceeds upper bound " << upper_bound_ << std::endl;
        }
        if(prop < lower_bound_) {
            std::cerr << "Proposal " << prop << " below lower bound " << lower_bound_ << std::endl;
        }
#endif
        return prop;
    }

    template<typename T, typename Engine>
    DoubleConstrainedContinuousRandomWalk<T, Engine>::DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter, T &target, double lower_bound, double upper_bound, Engine *rng)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng), lower_bound_(lower_bound), upper_bound_(upper_bound) {}


    template<typename T, typename Engine>
    DoubleConstrainedContinuousRandomWalk<T, Engine>::DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter,
                                                                                            T &target, double lower_bound, double upper_bound, Engine *rng,
                                                                                            double variance)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance), lower_bound_(lower_bound), upper_bound_(upper_bound) {}


    template<typename T, typename Engine>
    DoubleConstrainedContinuousRandomWalk<T, Engine>::DoubleConstrainedContinuousRandomWalk(parameters::Parameter<double> &parameter,
                                                                                            T &target, double lower_bound, double upper_bound, Engine *rng,
                                                                                            double variance,
                                                                                            double minVariance,
                                                                                            double maxVariance)
            : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance, minVariance, maxVariance), lower_bound_(lower_bound), upper_bound_(upper_bound) {}

}


#endif//TRANSMISSION_NETWORKS_APP_DOUBLECONSTRAINEDCONTINUOUSRANDOMWALK_H
