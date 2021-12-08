//
// Created by Maxwell Murphy on 7/29/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H

#include <algorithm>
#include <cmath>

#include "ContinuousRandomWalk.h"

namespace transmission_nets::core::samplers {

    // Continuous Random Walk constrained to a range (lower, upper)

    template<typename T, typename Engine = boost::random::mt19937, typename U = double>
    class ConstrainedContinuousRandomWalk : public ContinuousRandomWalk<T, Engine> {
    public:
        ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng);

        ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng, double variance);

        ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng, double variance, double minVariance, double maxVariance);

    private:
        using ContinuousRandomWalk<T, Engine>::rng_;
        using ContinuousRandomWalk<T, Engine>::variance_;
        using ContinuousRandomWalk<T, Engine>::normal_dist_;
        using ContinuousRandomWalk<T, Engine>::parameter_;
        U lower_bound_;
        U upper_bound_;
        [[nodiscard]] Likelihood logMetropolisHastingsAdjustment(U curr, U proposed) const noexcept override;
        double sampleProposal() noexcept override;
    };


    template<typename T, typename Engine, typename U>
    Likelihood ConstrainedContinuousRandomWalk<T, Engine, U>::logMetropolisHastingsAdjustment(U curr, U proposed) const noexcept {
        auto adj = log(proposed - lower_bound_) +
                   log(upper_bound_ - proposed) -
                   log(curr - lower_bound_) -
                   log(upper_bound_ - curr);
        return adj;
    }

    template<typename T, typename Engine, typename U>
    double ConstrainedContinuousRandomWalk<T, Engine, U>::sampleProposal() noexcept {
        double eps = normal_dist_(*rng_) * variance_;
        double unconstrained = std::log(parameter_->value() - lower_bound_) - std::log(upper_bound_ - parameter_->value());
        double exp_prop = std::exp(eps + unconstrained);
        double prop = (upper_bound_ * exp_prop + lower_bound_) / (exp_prop + 1);
        prop = std::clamp(prop, lower_bound_, upper_bound_);
#ifndef NDEBUG
        if (prop > upper_bound_) {
            std::cerr << "Proposal " << prop << " exceeds upper bound " << upper_bound_ << std::endl;
        }
        if (prop < lower_bound_) {
            std::cerr << "Proposal " << prop << " below lower bound " << lower_bound_ << std::endl;
        }
#endif
        return prop;
    }

    template<typename T, typename Engine, typename U>
    ConstrainedContinuousRandomWalk<T, Engine, U>::ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter, std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng)
        : ContinuousRandomWalk<T, Engine>(parameter, target, rng), lower_bound_(lower_bound), upper_bound_(upper_bound) {}


    template<typename T, typename Engine, typename U>
    ConstrainedContinuousRandomWalk<T, Engine, U>::ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter,
                                                                                   std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng,
                                                                                   double variance)
        : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance), lower_bound_(lower_bound), upper_bound_(upper_bound) {}


    template<typename T, typename Engine, typename U>
    ConstrainedContinuousRandomWalk<T, Engine, U>::ConstrainedContinuousRandomWalk(std::shared_ptr<parameters::Parameter<double>> parameter,
                                                                                   std::shared_ptr<T> target, U lower_bound, U upper_bound, std::shared_ptr<Engine> rng,
                                                                                   double variance,
                                                                                   double minVariance,
                                                                                   double maxVariance)
        : ContinuousRandomWalk<T, Engine>(parameter, target, rng, variance, minVariance, maxVariance), lower_bound_(lower_bound), upper_bound_(upper_bound) {}

}// namespace transmission_nets::core::samplers


#endif//TRANSMISSION_NETWORKS_APP_CONSTRAINEDCONTINUOUSRANDOMWALK_H
