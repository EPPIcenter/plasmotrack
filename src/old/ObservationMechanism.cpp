//
// Created by Maxwell Murphy on 10/16/19.
//

#include "ObservationMechanism.h"


LogLik ObservationMechanism::calc_prob_observed(const ObservedNode& node, const std::vector<int>& total_alleles,
                                                const int total_loci, const ObservationMechanismParameters params) {
    assert(total_loci > 0);
    assert(params.epsilon_pos < 1 && params.epsilon_pos > 0);
    assert(params.epsilon_neg < 1 && params.epsilon_neg > 0);

    LogLik llik = 0.0;
    unsigned int true_positive = 0;
    unsigned int true_negative = 0;
    unsigned int false_positive = 0;
    unsigned int false_negative = 0;
    for (int i = 0; i < total_loci; ++i) {
        if (node.has_data[i]) {
            true_positive += (node.observed_genotype[i] & node.latent_genotype[i]).count();
            true_negative += (~node.observed_genotype[i] & ~node.latent_genotype[i]).count() -
                             (MAX_ALLELES - total_alleles[i]); // Excess 0s off the end of the bitset
            false_positive += (node.observed_genotype[i] & ~node.latent_genotype[i]).count();
            false_negative += (~node.observed_genotype[i] & node.latent_genotype[i]).count();
        }
    }

    llik += (LogLik) true_positive * log(1 - params.epsilon_pos) +
            (LogLik) true_negative * log(1 - params.epsilon_neg) +
            (LogLik) false_positive * log(params.epsilon_pos) +
            (LogLik) false_negative * log(params.epsilon_neg);

    return llik;
}

ObservationMechanism::ObservationMechanism(const ObservationMechanismParameters &observationMechParameters)
        : observation_mech_parameters(observationMechParameters) {}
