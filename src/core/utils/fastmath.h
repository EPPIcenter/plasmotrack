//
// exp and log for loops over many cells, written so the compiler can vectorise them (no library call per element).
// Accurate to a few units in the last place over the ranges they are used on; not a general replacement.
//

#ifndef TRANSMISSION_NETWORKS_APP_FASTMATH_H
#define TRANSMISSION_NETWORKS_APP_FASTMATH_H

#include <bit>
#include <cmath>
#include <cstdint>
#include <limits>

namespace transmission_nets::core::utils {

    // exp(x) for x in [-700, 0].
    inline double expNonPositive(const double x) noexcept {
        constexpr double log2e = 1.4426950408889634;
        constexpr double ln2hi = 6.93147180369123816490e-01;
        constexpr double ln2lo = 1.90821492927058770002e-10;
        const double n         = std::nearbyint(x * log2e);
        const double r         = (x - n * ln2hi) - n * ln2lo;// |r| <= ln2 / 2
        // Taylor series to degree 13: the remainder is below 1e-17 relative on |r| <= ln2 / 2.
        double p = 1.0 / 6227020800.0;
        p        = p * r + 1.0 / 479001600.0;
        p        = p * r + 1.0 / 39916800.0;
        p        = p * r + 1.0 / 3628800.0;
        p        = p * r + 1.0 / 362880.0;
        p        = p * r + 1.0 / 40320.0;
        p        = p * r + 1.0 / 5040.0;
        p        = p * r + 1.0 / 720.0;
        p        = p * r + 1.0 / 120.0;
        p        = p * r + 1.0 / 24.0;
        p        = p * r + 1.0 / 6.0;
        p        = p * r + 0.5;
        p        = p * r + 1.0;
        p        = p * r + 1.0;
        const auto e = static_cast<std::int64_t>(n) + 1023;
        return p * std::bit_cast<double>(static_cast<std::uint64_t>(e) << 52);
    }

    // log(x) for normal positive x.
    inline double logPositive(const double x) noexcept {
        constexpr double sqrt2 = 1.4142135623730951;
        constexpr double ln2   = 0.6931471805599453;
        const auto bits        = std::bit_cast<std::uint64_t>(x);
        double e               = static_cast<double>(static_cast<std::int64_t>((bits >> 52) & 0x7ff) - 1023);
        double m               = std::bit_cast<double>((bits & 0x000fffffffffffffULL) | 0x3ff0000000000000ULL);// [1, 2)
        if (m > sqrt2) {
            m *= 0.5;
            e += 1.0;
        }
        // log(m) = 2 atanh(z), z = (m - 1) / (m + 1), |z| <= 0.1716; the series to z^23 is below 1e-17.
        const double z  = (m - 1.0) / (m + 1.0);
        const double z2 = z * z;
        double p        = 1.0 / 23.0;
        p               = p * z2 + 1.0 / 21.0;
        p               = p * z2 + 1.0 / 19.0;
        p               = p * z2 + 1.0 / 17.0;
        p               = p * z2 + 1.0 / 15.0;
        p               = p * z2 + 1.0 / 13.0;
        p               = p * z2 + 1.0 / 11.0;
        p               = p * z2 + 1.0 / 9.0;
        p               = p * z2 + 1.0 / 7.0;
        p               = p * z2 + 1.0 / 5.0;
        p               = p * z2 + 1.0 / 3.0;
        p               = p * z2 + 1.0;
        return e * ln2 + 2.0 * z * p;
    }

}// namespace transmission_nets::core::utils

#endif//TRANSMISSION_NETWORKS_APP_FASTMATH_H
