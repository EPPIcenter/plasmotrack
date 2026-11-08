//
// The per-Allele rate a genotype proposal reads from an Infection's sampled expected false positives or negatives.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONRATE_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONRATE_H

#include <algorithm>

namespace transmission_nets::core::samplers::genetics {

    // An expected count spread over `alleles` Alleles, kept in [1e-4, 0.49] so that the proposal rules no outcome out
    // and still favours the Observed genotype.
    inline double perAlleleRate(const double expectedCount, const double alleles) noexcept {
        return std::clamp(expectedCount / alleles, 1e-4, 0.49);
    }

}// namespace transmission_nets::core::samplers::genetics

#endif//TRANSMISSION_NETWORKS_APP_OBSERVATIONRATE_H
