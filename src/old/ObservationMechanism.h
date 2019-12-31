//
// Created by Maxwell Murphy on 10/16/19.
//

#ifndef MALARIA_NET_OBSERVATIONMECHANISM_H
#define MALARIA_NET_OBSERVATIONMECHANISM_H

#include <cmath>
#include "utils.h"
#include "Node.h"

struct ObservationMechanismParameters {
    Probability epsilon_pos;
    Probability epsilon_neg;
};

class ObservationMechanism {
public:
    ObservationMechanism(const ObservationMechanismParameters &observation_mech_parameters);

private:
    ObservationMechanismParameters observation_mech_parameters;

public:
    static LogLik calc_prob_observed(const ObservedNode &node, const std::vector<int>& total_alleles,
            int total_loci, ObservationMechanismParameters params);
};


#endif //MALARIA_NET_OBSERVATIONMECHANISM_H
