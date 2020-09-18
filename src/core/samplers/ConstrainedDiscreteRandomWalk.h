//
// Created by Maxwell Murphy on 4/16/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_CONSTRAINEDDISCRETERANDOMWALK_H
#define TRANSMISSION_NETWORKS_APP_CONSTRAINEDDISCRETERANDOMWALK_H

#include <algorithm>

#include "DiscreteRandomWalk.h"


namespace transmission_nets::core::samplers {

// Discrete Random Walk constrained to a range [lower, upper]

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    class ConstrainedDiscreteRandomWalk : public DiscreteRandomWalk<T, Engine> {
    public:
        ConstrainedDiscreteRandomWalk(parameters::Parameter<int> &parameter, T &target, Engine *rng)
                : DiscreteRandomWalk<T, Engine>(parameter, target, rng) {
            assert(parameter.value() > LOWER_BOUND);
            assert(parameter.value() < UPPER_BOUND);
        }


        ConstrainedDiscreteRandomWalk(parameters::Parameter<int> &parameter, T &target, Engine *rng, unsigned int maxDistance)
                : DiscreteRandomWalk<T, Engine>(parameter, target, rng, maxDistance) {
            assert(parameter.value() > LOWER_BOUND);
            assert(parameter.value() < UPPER_BOUND);
        }

        double sampleStride(int current) noexcept override;

        double logMetropolisHastingsAdjustment(int current, int proposed) noexcept override;

    protected:
        using DiscreteRandomWalk<T, Engine>::uniform_dist_;
        using DiscreteRandomWalk<T, Engine>::stride_sampling_dist_;
        using DiscreteRandomWalk<T, Engine>::max_distance_;
        using DiscreteRandomWalk<T, Engine>::rng_;
    };

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    double ConstrainedDiscreteRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::sampleStride(int current) noexcept {
        auto left_range = std::min(current - LOWER_BOUND, int(max_distance_));
        auto right_range = std::min(UPPER_BOUND - current, int(max_distance_));
        int step;

        if(uniform_dist_(*rng_) <= (double(left_range) / double(right_range)) ) {
            stride_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, left_range));
            step = -stride_sampling_dist_(*rng_);
        } else {
            stride_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, right_range));
            step = stride_sampling_dist_(*rng_);
        }
        return step;
    }

    template<int LOWER_BOUND, int UPPER_BOUND, typename T, typename Engine>
    double
    ConstrainedDiscreteRandomWalk<LOWER_BOUND, UPPER_BOUND, T, Engine>::logMetropolisHastingsAdjustment(int current,
                                                                                                        int proposed) noexcept {
        auto curr_range = std::min(current - LOWER_BOUND, int(max_distance_)) + std::min(UPPER_BOUND - current, int(max_distance_));
        auto proposed_range = std::min(proposed - LOWER_BOUND, int(max_distance_)) + std::min(UPPER_BOUND - proposed, int(max_distance_));

        return log(curr_range) - log(proposed_range);

    }

}

#endif //TRANSMISSION_NETWORKS_APP_CONSTRAINEDDISCRETERANDOMWALK_H
