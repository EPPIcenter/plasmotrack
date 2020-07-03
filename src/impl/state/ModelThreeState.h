//
// Created by Maxwell Murphy on 6/18/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H
#define TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H

#include <vector>
#include <map>

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/AlleleFrequencyContainer.h"
#include "core/containers/TransmissionNetwork.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Parameter.h"


struct ModelThreeState {
    static constexpr int MAX_ALLELES = 32;
    static constexpr int MAX_PARENT_SET = 1;
    using LocusImpl = Locus;
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = Infection<GeneticsImpl, LocusImpl>;
    using AlleleFrequencyImpl = Simplex;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;


    std::map<std::string, LocusImpl *> loci{};
    std::vector<InfectionEvent *> infections{};

    AlleleFrequencyContainer alleleFrequencies;

    // Network Structure
    TransmissionNetwork<InfectionEvent> transmissionNetwork;

    // Observation Process
    Parameter<double> observationFalsePositiveRate;
    Parameter<double> observationFalseNegativeRate;

    // Node Transmission Process
    Parameter<double> geometricGenerationProb;
    Parameter<double> lossProb;
    Parameter<double> mutationProb;

    // Source Transmission Process
    Parameter<double> meanCOI;
};


#endif //TRANSMISSION_NETWORKS_APP_MODELTHREESTATE_H
