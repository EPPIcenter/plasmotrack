//
// Created by Maxwell Murphy on 4/14/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_ORDERSAMPLER_H
#define TRANSMISSION_NETWORKS_APP_ORDERSAMPLER_H

#include <boost/math/distributions.hpp>
#include <boost/random.hpp>
#include <boost/range/algorithm.hpp>

#include "core/samplers/AbstractSampler.h"
#include "core/parameters/Ordering.h"

template<typename T,  typename OrderingElement, typename Engine=boost::random::mt19937>
class OrderSampler : public AbstractSampler {
public:
    OrderSampler(Ordering<OrderingElement> &parameter, T &target, Engine *rng) noexcept;

    OrderSampler(Ordering<OrderingElement> &parameter, T &target, Engine *rng, unsigned int max_distance) noexcept;

    void update() noexcept override;

    std::tuple<int, int> sampleProposal() noexcept;

    [[nodiscard]] unsigned int acceptances() noexcept;

    [[nodiscard]] unsigned int rejections() noexcept;

    [[nodiscard]] double acceptanceRate() noexcept;

private:
    Ordering<OrderingElement> &parameter_;
    T &target_;
    Engine *rng_;
    unsigned int max_distance_;
    unsigned int num_elements_;
    boost::random::uniform_01<> uniform_dist_{};
    boost::random::uniform_int_distribution<> pivot_sampling_dist_;
    boost::random::uniform_int_distribution<> offset_sampling_dist_;

    unsigned int acceptances_ = 0;
    unsigned int rejections_ = 0;
    unsigned int total_updates_ = 0;
};


template<typename T, typename OrderingElement, typename Engine>
void OrderSampler<T, OrderingElement, Engine>::update() noexcept {
    double curLik = target_.value();
    parameter_.saveState();

    auto [pivot, offset] = sampleProposal();

    parameter_.swap(pivot, pivot + offset);

    const double acceptanceRatio = target_.value() - curLik;

    const bool accept = log(uniform_dist_(*rng_)) <= acceptanceRatio;

    if (accept) {
        acceptances_++;
        parameter_.acceptState();
    } else {
        rejections_++;
        parameter_.restoreState();
    }

    total_updates_++;
}

template<typename T, typename OrderingElement, typename Engine>
OrderSampler<T, OrderingElement, Engine>::OrderSampler(Ordering<OrderingElement> &parameter, T &target, Engine *rng) noexcept
        :parameter_(parameter), target_(target), rng_(rng), max_distance_(1) {
    offset_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, max_distance_));
    num_elements_ = parameter.value().size();
    assert(max_distance_ <= (num_elements_ / 2));
}

template<typename T, typename OrderingElement, typename Engine>
OrderSampler<T, OrderingElement, Engine>::OrderSampler(Ordering<OrderingElement> &parameter, T &target, Engine *rng, unsigned int max_distance) noexcept
        :parameter_(parameter), target_(target), rng_(rng), max_distance_(max_distance) {
    offset_sampling_dist_.param(boost::random::uniform_int_distribution<>::param_type(1, max_distance_));
    num_elements_ = parameter.value().size();
    assert(max_distance_ <= (num_elements_ / 2));
}

template<typename T, typename OrderingElement, typename Engine>
std::tuple<int, int> OrderSampler<T, OrderingElement, Engine>::sampleProposal() noexcept {
    auto offset = offset_sampling_dist_(*rng_);
    offset = uniform_dist_(*rng_) > .5 ? offset : -offset;

    if (offset > 0) {
        pivot_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type{0, static_cast<int>(num_elements_) - 1 - offset}
        );
    } else {
        pivot_sampling_dist_.param(
                boost::random::uniform_int_distribution<>::param_type{-offset, static_cast<int>(num_elements_) - 1}
        );
    }

    auto pivot = pivot_sampling_dist_(*rng_);

    return std::tuple<int, int>(pivot, offset);
}

template<typename T, typename OrderingElement, typename Engine>
unsigned int OrderSampler<T, OrderingElement, Engine>::acceptances() noexcept {
    return acceptances_;
}

template<typename T, typename OrderingElement, typename Engine>
unsigned int OrderSampler<T, OrderingElement, Engine>::rejections() noexcept {
    return rejections_;
}

template<typename T, typename OrderingElement, typename Engine>
double OrderSampler<T, OrderingElement, Engine>::acceptanceRate() noexcept {
    return double(acceptances_) / (acceptances_ + rejections_);
}

#endif //TRANSMISSION_NETWORKS_APP_ORDERSAMPLER_H
