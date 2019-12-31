//
// Created by Maxwell Murphy on 11/25/19.
//

#ifndef TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISMPARAMETERS_H
#define TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISMPARAMETERS_H

#include <memory>
#include <bitset>
#include "utils.h"
#include "TransmissionNetworkModelConfig.h"

class ObservationMechanismParameters {
public:
    std::shared_ptr<Probability> false_positive_rate;
    std::shared_ptr<Probability> false_negative_rate;
    std::shared_ptr<std::bitset<MAX_ALLELES>> latent_locus;
    std::shared_ptr<std::bitset<MAX_ALLELES>> observed_locus;
};


#endif //TRANSMISSION_NETWORKS_APP_OBSERVATIONMECHANISMPARAMETERS_H
