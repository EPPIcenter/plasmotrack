//
// Created by Maxwell Murphy on 2/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NUMERICS_H
#define TRANSMISSION_NETWORKS_APP_NUMERICS_H

#include <numeric>
#include <cmath>

/**
 * Does not assume the iterable is sorted.
 * @tparam Iter
 * @param begin
 * @param end
 * @return
 */
template<typename Iter>
typename std::iterator_traits<Iter>::value_type logSumExp(Iter begin, Iter end) {
    using ValueType = typename std::iterator_traits<Iter>::value_type;

    if (begin == end) {
        return ValueType{};
    }

    auto max_el = *std::max_element(begin, end);
    auto sum = std::accumulate(
            begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); }
            );
    return max_el + log(sum);
}

/**
 * Assumes the provided iterable is sorted largest to smallest.
 * @tparam Iter
 * @param begin
 * @param end
 * @return
 */
template<typename Iter>
typename std::iterator_traits<Iter>::value_type sortedLogSumExp(Iter begin, Iter end) {
    using ValueType = typename std::iterator_traits<Iter>::value_type;

    if (begin == end) {
        return ValueType{};
    }

    auto max_el = *begin;
    auto sum = std::accumulate(
            begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); }
    );
    return max_el + std::log(sum);
}

/**
 * Assumes the max value is known
 * @tparam Iter
 * @param begin
 * @param end
 * @return
 */
template<typename Iter>
typename std::iterator_traits<Iter>::value_type logSumExpKnownMax(Iter begin, Iter end, typename std::iterator_traits<Iter>::value_type max_el) {
    using ValueType = typename std::iterator_traits<Iter>::value_type;

    if (begin == end) {
        return ValueType{};
    }

    auto sum = std::accumulate(
            begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); }
    );
    return max_el + std::log(sum);
}

#endif //TRANSMISSION_NETWORKS_APP_NUMERICS_H
