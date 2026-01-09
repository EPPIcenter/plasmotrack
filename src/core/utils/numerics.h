//
// Created by Maxwell Murphy on 2/26/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_NUMERICS_H
#define TRANSMISSION_NETWORKS_APP_NUMERICS_H

#include "LogPQ.h"


#include <fmt/core.h>

#include <algorithm>
#include <cmath>
#include <numeric>

// Check if std::execution::unseq is available
// libc++ on macOS may not support unseq even with C++20
#if __has_include(<execution>)
#include <execution>
// Check if unseq is actually available (some implementations don't support it)
#if defined(__cpp_lib_parallel_algorithm) && __cpp_lib_parallel_algorithm >= 201603L
    #if defined(_LIBCPP_VERSION) && _LIBCPP_VERSION < 14000
        // libc++ before version 14 doesn't support unseq
        #define HAS_STD_EXECUTION_UNSEQ 0
    #else
        #define HAS_STD_EXECUTION_UNSEQ 1
    #endif
#else
    #define HAS_STD_EXECUTION_UNSEQ 0
#endif
#else
#define HAS_STD_EXECUTION_UNSEQ 0
#endif

namespace transmission_nets::core::utils {
    template<typename T>
    inline double fast_exp(T x) {
        x = 1.0 + x / 1024;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        x *= x;
        return x;
    }

    inline double _int_as_double(int32_t a) {
        double r;
        memcpy(&r, &a, sizeof r);
        return r;
    }
    inline int32_t _double_as_int(double a) {
        int32_t r;
        memcpy(&r, &a, sizeof r);
        return r;
    }
    // https://stackoverflow.com/a/39822314/2755374
    inline double fast_log(double a) {
        double m, r, s, t, i, f;
        int32_t e;

        e = (_double_as_int(a) - 0x3f2aaaab) & 0xff800000;
        m = _int_as_double(_double_as_int(a) - e);
        i = (double) e * 1.19209290e-7f;// 0x1.0p-23
        /* m in [2/3, 4/3] */
        f = m - 1.0;
        s = f * f;
        /* Compute log1p(f) for f in [-1/3, 1/3] */
        r = fmaf(0.230836749f, f, -0.279208571f);// 0x1.d8c00p-3, -0x1.1de8dap-2
        t = fmaf(0.331826031f, f, -0.498910338f);// 0x1.53ca34p-2, -0x1.fee25ap-2
        r = fmaf(r, s, t);
        r = fmaf(r, s, f);
        r = fmaf(i, 0.693147182f, r);// 0x1.62e430p-1 // log(2)
        return r;
    }

    inline double logit(double x) {
        if (x < .5) {
            return std::log(x) - std::log1p(-x);
        } else {
            return std::log(x / (1 - x));
        }
    }

    template<typename Iter>
    typename std::vector<double> logit(const Iter& begin, const Iter& end) {
        std::vector<double> out{};
        std::transform(begin, end, std::back_inserter(out), [](double c) -> double { return logit(c); });
        return out;
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    std::vector<double> logit(const It& x) {
        return logit(x.begin(), x.end());
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    double logitSum(const It& logx) {
        auto logx_sorted = logx;
        std::sort(logx_sorted.rbegin(), logx_sorted.rend());
        auto lpq = LogPQ(logx_sorted);
        return logitSum(logx_sorted, lpq);
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    double logitSum(const It& logx_sorted, const LogPQ& lpq) {
        double out;
        double cumsum = 0.0;
        if (logx_sorted[0] < 0) {
            double lp1 = lpq.logP[0];
            for (unsigned int ii = 1; ii < lpq.logP.size(); ++ii) {
                cumsum += std::exp(lpq.logP[ii] - lp1);
            }
            out = lp1 + std::log1p(cumsum);
        } else {
            double lq1 = lpq.logQ[0];
            for (unsigned int ii = 1; ii < lpq.logQ.size(); ++ii) {
                cumsum += std::exp(lpq.logP[ii]);
            }
            out = std::log1p(-std::exp(lq1) + cumsum);
        }
        return out;
    }

    inline double expit(double x) {
        return 1 / (1 + std::exp(-x));
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    std::vector<double> logitScale(const It& x, double scale) {
        std::vector<double> out{};
        out.reserve(x.size());
        double l2;
        double u;
        double v;
        double ev;
        double eumo;
        bool ok;
        for (auto const& el : x) {
            ok = (scale < std::log(2)) and (std::abs(scale) < std::abs(el + scale));
            u = -scale - (!ok) * el;
            v = -scale - ok * el;
            ev = std::exp(v);
            eumo = std::expm1(u);

            if (std::isinf(eumo)) {
                l2 = std::max(u, v) + std::log1p(std::exp(-std::abs(u - v)));
            } else {
                l2 = std::log(eumo + ev);
            }

            if (v > std::log(2 * std::abs(eumo))) {
                out.push_back(-(v + std::log1p(eumo / ev)));
            } else {
                out.push_back(l2);
            }
        }
        return out;
    }

    template<typename Iter>
    typename std::vector<double> expit(const Iter& begin, const Iter& end) {
        std::vector<double> out{};
        std::transform(begin, end, std::back_inserter(out), [](double c) -> double { return expit(c); });
        return out;
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    std::vector<double> expit(const It& x) {
        return expit(x.begin(), x.end());
    }

    inline double logLogit(double x) {
        if (x < std::log(.5)) {
            return x - std::log1p(-std::exp(x));
        } else {
            return x - log(-expm1(x));
        }
    }

    /**
     * @brief Numerically stable log(exp(a) + exp(b))
     * @param a First argument
     * @param b Second argument
     * @return The log of the sum of exp(a) and exp(b)
     */
    inline double logSumExp(const double a, const double b) {
        double max_el = std::max(a, b);
        if (max_el == -std::numeric_limits<double>::infinity()) {
            return -std::numeric_limits<double>::infinity();
        }
        double sum = std::exp(a - max_el) + std::exp(b - max_el);
        return max_el + std::log(sum);
    }

    /**
     * Numerically stable log(∑(exp(a)))
     * Does not assume the iterable is sorted.
     * @tparam Iter implements iterable
     * @param begin iterator pointer to beginning
     * @param end iterator pointer to end
     * @return
     */
    template<typename Iter>
    typename std::iterator_traits<Iter>::value_type logSumExp(const Iter& begin, const Iter& end) {
        using ValueType = typename std::iterator_traits<Iter>::value_type;

        if (begin == end) {
            return ValueType{};
        }

        auto max_el = *std::max_element(begin, end);

        if (max_el == -std::numeric_limits<ValueType>::infinity()) {
            return -std::numeric_limits<ValueType>::infinity();
        }

#if HAS_STD_EXECUTION_UNSEQ
        auto sum = std::transform_reduce(
            std::execution::unseq,
            begin,
            end,
            ValueType{},
            std::plus<ValueType>{},
            [max_el](ValueType a) { return std::exp(a - max_el); }
            );
#else
        auto sum = std::transform_reduce(
            begin,
            end,
            ValueType{},
            std::plus<ValueType>{},
            [max_el](ValueType a) { return std::exp(a - max_el); }
            );
#endif

        // auto sum = std::accumulate(
        //         begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); });
        return max_el + std::log(sum);
    }

    template<typename It, typename T = std::decay_t<decltype(*begin(std::declval<It>()))>>
    double logSumExp(const It& x) {
        return logSumExp(x.begin(), x.end());
    }

    template<typename T>
    double logDiffExp(const T a, const T b, const bool verbose = false) {
        /**
         * return the numerically stable log(exp(a) - exp(b))
         * if b > a, returns -inf
         */
        if (b > a) {
            //#ifndef NDEBUG
            //            fmt::print(stderr, "Warning -- b > a : {} {} {}\n", b, a, b + std::log(1 - std::exp(a - b)));
            //            std::cerr << "Warning -- b > a : " << b << " " << a << " " << b + std::log(1 - std::exp(a - b)) << std::endl;
            //#endif
            return -std::numeric_limits<T>::infinity();
        }
        if (verbose) {
            fmt::print("Log Diff: {}, {}\n", std::exp(b - a), std::log(1.0 - std::exp(b - a)));
        }
        //        mpfr_t a_;
        //        mpfr_t b_;
        //        mpfr_t out_;
        //        mpfr_t tmp1_;
        //        mpfr_t one_;
        //        mpfr_inits2(256, a_, b_, out_, tmp1_, one_, (mpfr_ptr) 0);
        //        mpfr_set_ld(a_, a, mpfr_get_default_rounding_mode());
        //        mpfr_set_ld(b_, b, mpfr_get_default_rounding_mode());
        //        mpfr_set_ld(one_, (double) 1.0, mpfr_get_default_rounding_mode());
        //
        //        mpfr_sub(out_, b_, a_, mpfr_get_default_rounding_mode());
        //        mpfr_exp(out_, out_, mpfr_get_default_rounding_mode());
        //        mpfr_sub(out_, one_, out_, mpfr_get_default_rounding_mode());
        //        mpfr_log(out_, out_, mpfr_get_default_rounding_mode());
        //        mpfr_add(out_, a_, out_, mpfr_get_default_rounding_mode());
        //        return mpfr_get_ld(out_, mpfr_get_default_rounding_mode());

        return a + std::log(1.0 - std::exp(b - a));
    }

    template<typename T>
    T absLogDiff(const T a, const T b) {
        return b > a ? b + std::log(1 - std::exp(a - b)) : a + std::log(1 - std::exp(b - a));
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
                begin, end, ValueType{}, [max_el](ValueType a, ValueType b) { return a + std::exp(b - max_el); });
        return max_el + std::log(sum);
    }


    /**
     * Assumes the max value is known
     * @tparam Iter
     * @param begin
     * @tparam Iter
     * @param end
     * @param max_el
     * @return
     */
    template<typename Iter>
    typename std::iterator_traits<Iter>::value_type logSumExpKnownMax(const Iter& begin, const Iter& end, const typename std::iterator_traits<Iter>::value_type& max_el) {
        using ValueType = typename std::iterator_traits<Iter>::value_type;

        if (begin == end) {
            return ValueType{};
        }

        if (max_el == -std::numeric_limits<ValueType>::infinity()) {
            return -std::numeric_limits<ValueType>::infinity();
        }

        // auto sum = std::transform_reduce(std::execution::unseq,
        //         begin,
        //         end,
        //         ValueType{},
        //         std::plus<ValueType>{},
        //         [max_el](ValueType a) { return std::exp(a - max_el); });

        // auto sum = std::accumulate(
                // begin, end, ValueType{}, [max_el](const ValueType a, const ValueType b) { return a + std::exp(b - max_el); });

#if HAS_STD_EXECUTION_UNSEQ
        auto red_sum = std::transform_reduce(
                std::execution::unseq,
                begin,
                end,
                ValueType{},
                std::plus<ValueType>{},
                [max_el](const ValueType a) { return std::exp(a - max_el); });
#else
        auto red_sum = std::transform_reduce(
                begin,
                end,
                ValueType{},
                std::plus<ValueType>{},
                [max_el](const ValueType a) { return std::exp(a - max_el); });
#endif
        // if (std::abs(sum - red_sum) > 1e-6) {
            // fmt::print("Sum: {}, Red Sum: {}\n", sum, red_sum);
        // }
        return max_el + std::log(red_sum);
    }

    template<typename T>
    T logSumExpKnownMax(const std::vector<T>& x, const T& max_el) {
        return logSumExpKnownMax(x.begin(), x.end(), max_el);
    }

    /**
     * Normalize the exp(x) vector
     * @tparam T
     * @param x
     * @return
     */
    template<typename T>
    std::vector<double> expNormalize(const T& x) {
        std::vector<double> out{};
        out.reserve(x.size());

        double sum = 0.0;
        auto max_el = *std::max_element(x.begin(), x.end());

        for (auto& el : x) {
            out.push_back(std::exp(el - max_el));
            sum += out.back();
        }

        for (auto& el : out) {
            el = el / sum;
        }

        return out;
    }

    inline constexpr unsigned long const_pow_helper(unsigned long base, unsigned long exp) {
        return (exp == 0) ? 1 : base * const_pow_helper(base, exp - 1);
    }

    inline constexpr unsigned long const_pow(unsigned int a, unsigned int b) {
        return const_pow_helper(a, b);
    }


    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE>
    constexpr std::array<unsigned int, MAX_PARENTSET_SIZE + 1> idxToKvec(unsigned int idx) {
        std::array<unsigned int, MAX_PARENTSET_SIZE + 1> indxs;
        indxs.fill(1);

        for (unsigned int i = 0; i < MAX_PARENTSET_SIZE; ++i) {
            indxs[i] += idx % MAX_TRANSMISSIONS;
            idx /= MAX_TRANSMISSIONS;
        }

        return indxs;
    }

    template<unsigned int MAX_TRANSMISSIONS, unsigned int MAX_PARENTSET_SIZE>
    constexpr std::array<std::array<unsigned int, MAX_PARENTSET_SIZE + 1>, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> initKvecs() {
        std::array<std::array<unsigned int, MAX_PARENTSET_SIZE + 1>, core::utils::const_pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE)> kvecs;
        for (unsigned int i = 0; i < pow(MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE); i++) {
            kvecs[i] = idxToKvec<MAX_TRANSMISSIONS, MAX_PARENTSET_SIZE>(i);
        }
        return kvecs;
    }

}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_NUMERICS_H
