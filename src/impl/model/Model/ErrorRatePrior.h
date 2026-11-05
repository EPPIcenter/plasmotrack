//
// The hierarchical error model's prior on one Infection's rate: Gamma(shape, mean / shape) truncated to the range the
// rate's walk keeps it in, [ERROR_RATE_LOWER, ERROR_RATE_UPPER]. Without the truncation's normalizer the walk's bounds
// would weigh the mean and shape by the Gamma's mass inside them, once per Infection.
//

#ifndef TRANSMISSION_NETWORKS_APP_ERRORRATEPRIOR_H
#define TRANSMISSION_NETWORKS_APP_ERRORRATEPRIOR_H

#include "config.h"

#include <boost/math/special_functions/gamma.hpp>

#include <cmath>
#include <limits>

namespace transmission_nets::impl::Model {

    // log of the Gamma(shape, mean / shape) density at x, truncated to [lower, upper].
    inline double truncatedGammaLogDensity(const double x, const double mean, const double shape, const double lower = ERROR_RATE_LOWER, const double upper = ERROR_RATE_UPPER) {
        if (not(x >= lower and x <= upper)) {
            return -std::numeric_limits<double>::infinity();
        }
        const double scale  = mean / shape;
        const double inside = boost::math::gamma_p(shape, upper / scale) - boost::math::gamma_p(shape, lower / scale);
        return (shape - 1) * std::log(x) - x / scale - std::lgamma(shape) - shape * std::log(scale) - std::log(inside);
    }

}// namespace transmission_nets::impl::Model

#endif//TRANSMISSION_NETWORKS_APP_ERRORRATEPRIOR_H
