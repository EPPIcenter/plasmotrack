//
// An infection-duration prior, for proposals that draw durations from it within a range.
//

#ifndef TRANSMISSION_NETWORKS_APP_DURATIONSHAPE_H
#define TRANSMISSION_NETWORKS_APP_DURATIONSHAPE_H

#include <boost/random.hpp>

#include <algorithm>
#include <cmath>
#include <utility>
#include <vector>

namespace transmission_nets::core::samplers::specialized {

    /*
     * A duration density given as a weight per whole day: a duration d falls in day round(d), so day k covers
     * [k - 0.5, k + 0.5) and the density at d is the weight of its day.
     */
    class DurationShape {
    public:
        explicit DurationShape(std::vector<double> weights) : weights_(std::move(weights)) {}

        // The density at a duration; 0 outside the days.
        [[nodiscard]] double density(const double duration) const {
            const auto day = static_cast<long>(std::round(duration));
            return day < 0 or day >= static_cast<long>(weights_.size()) ? 0.0 : weights_[day];
        }

        // The mass within [first, second).
        [[nodiscard]] double mass(const double first, const double second) const {
            double mass = 0;
            for (size_t day = 0; day < weights_.size(); ++day) {
                const double overlap = std::min(second, day + 0.5) - std::max(first, day - 0.5);
                if (overlap > 0) {
                    mass += weights_[day] * overlap;
                }
            }
            return mass;
        }

        // A duration drawn from the density restricted to [first, second), whose mass is `mass` (positive): a day by
        // its mass within the range, then a duration uniformly within the day and the range. Should rounding leave
        // mass unspent, the last day with any is used.
        template<typename Engine>
        double draw(const double first, const double second, const double mass, Engine& rng) const {
            boost::random::uniform_01<> uniform{};
            double remaining = uniform(rng) * mass;
            double low = first, high = first;
            for (size_t day = 0; day < weights_.size(); ++day) {
                const double dayLow  = std::max(first, day - 0.5);
                const double dayHigh = std::min(second, day + 0.5);
                const double weight  = dayHigh > dayLow ? weights_[day] * (dayHigh - dayLow) : 0.0;
                if (weight <= 0) {
                    continue;
                }
                low  = dayLow;
                high = dayHigh;
                if (remaining < weight) {
                    break;
                }
                remaining -= weight;
            }
            return low + uniform(rng) * (high - low);
        }

    private:
        std::vector<double> weights_;
    };

}// namespace transmission_nets::core::samplers::specialized

#endif//TRANSMISSION_NETWORKS_APP_DURATIONSHAPE_H
