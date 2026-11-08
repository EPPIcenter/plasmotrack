//
// A draw among options in proportion to exp(value), the values given on the log scale. The guided samplers' walks and
// Gibbs Proposals draw this way.
//

#ifndef TRANSMISSION_NETWORKS_APP_LOGWEIGHTEDDRAW_H
#define TRANSMISSION_NETWORKS_APP_LOGWEIGHTEDDRAW_H

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <vector>

namespace transmission_nets::core::samplers {

    // The largest value, and the sum of exp(value - top) in order. top is -inf or not finite when no option is.
    struct LogWeights {
        double top;
        double total;
    };

    inline LogWeights logWeights(const std::vector<double>& values) {
        LogWeights w{.top = *std::ranges::max_element(values), .total = 0};
        for (const double v : values) w.total += std::exp(v - w.top);
        return w;
    }

    // An index drawn in proportion to exp(values[i]): one uniform draw, walked through the options in order.
    template<typename Engine>
    size_t drawLogWeighted(const std::vector<double>& values, const LogWeights& w, Engine& rng) {
        boost::random::uniform_01<> unit;
        double r    = unit(rng) * w.total;
        size_t pick = 0;
        while (pick + 1 < values.size() and (r -= std::exp(values[pick] - w.top)) > 0) ++pick;
        return pick;
    }

}// namespace transmission_nets::core::samplers

#endif//TRANSMISSION_NETWORKS_APP_LOGWEIGHTEDDRAW_H
