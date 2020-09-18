//
// Created by Maxwell Murphy on 2/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NUMERICS_H
#define TRANSMISSION_NETWORKS_APP_NUMERICS_H

#include <algorithm>
#include <cmath>
#include <numeric>


namespace transmission_nets::core::utils {
    /**
     * Does not assume the iterable is sorted.
     * @tparam Iter
     * @param begin
     * @param end
     * @return
     */
    template<typename Iter>
    typename std::iterator_traits<Iter>::value_type logSumExp(const Iter& begin, const Iter& end) {
        using ValueType = typename std::iterator_traits<Iter>::value_type;

        if (begin == end) {
            return ValueType{};
        }

        auto max_el = *std::max_element(begin, end);

        if (max_el == -std::numeric_limits<double>::infinity()) {
            return -std::numeric_limits<double>::infinity();
        }

        auto sum = std::accumulate(
                begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); }
        );
        return max_el + std::log(sum);
    }

    template <typename T>
    double logSumExp(const T& iterable) {
        return logSumExp(iterable.begin(), iterable.end());
    }


    inline double logSumExp(const double a, const double b) {
        double max_el = std::max(a, b);
        if (max_el == -std::numeric_limits<double>::infinity()) {
            return -std::numeric_limits<double>::infinity();
        }

        double sum = std::exp(a - max_el) + std::exp(b - max_el);
        return max_el + std::log(sum);
    }

    /**
     * Assumes the provided iterable is sorted largest to smallest.
     * @tparam Iter
     * @param begin
     * @param end
     * @return
     */
    template<typename Iter>
    typename std::iterator_traits<Iter>::value_type sortedLogSumExp(const Iter& begin, const Iter& end) {
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
    typename std::iterator_traits<Iter>::value_type logSumExpKnownMax(const Iter& begin, const Iter& end, const typename std::iterator_traits<Iter>::value_type& max_el) {
        using ValueType = typename std::iterator_traits<Iter>::value_type;

        if (begin == end) {
            return ValueType{};
        }

        auto sum = std::accumulate(
                begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); }
        );
        return max_el + std::log(sum);
    }


    template<typename T>
    inline double logit(const T x) {
        if (x < .5) {
            return log(x) - log1p(x);
        } else {
            return log(x / (1 - x));
        }
    }

    template<typename T>
    inline double expit(const T x) {
        return 1 / (1 + exp(-x));
    }

    template<typename T>
    inline double logLogit(const T x) {
        if (x < log(.5)) {
            return x - log1p(-exp(x));
        } else {
            return x - log(-expm1(x));
        }
    }
}


#endif //TRANSMISSION_NETWORKS_APP_NUMERICS_H
