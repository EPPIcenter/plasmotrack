//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELONESTATE_H
#define TRANSMISSION_NETWORKS_APP_MODELONESTATE_H

#include <vector>
#include <map>

#include "core/containers/Infection.h"
#include "core/containers/Locus.h"
#include "core/containers/AlleleFrequencyContainer.h"

#include "core/datatypes/Alleles.h"
#include "core/datatypes/Simplex.h"

#include "core/parameters/Parameter.h"
#include "core/parameters/Ordering.h"


/// State implementation of a simple model
struct ModelOneState {
    static constexpr int MAX_ALLELES = 32;
    using LocusImpl = Locus;
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = Infection<GeneticsImpl, LocusImpl>;
    using AlleleFrequencyImpl = Simplex;
    using AlleleFrequencyContainerImpl = AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;


    std::map<std::string, LocusImpl *> loci{};
    std::vector<InfectionEvent *> infections{};

    AlleleFrequencyContainerImpl alleleFrequencies;

    // Network Structure
    Ordering<InfectionEvent> infectionEventOrdering;

    // Observation Process
    Parameter<double> observationFalsePositiveRate;
    Parameter<double> observationFalseNegativeRate;

    // Node Transmission Process
    Parameter<double> geometricGenerationProb;
    Parameter<double> ztMultiplicativeBinomialProb;
    Parameter<double> ztMultiplicativeBinomialAssoc;

    // Source Transmission Process
    Parameter<double> geometricCOIProb;

};


#endif //TRANSMISSION_NETWORKS_APP_MODELONESTATE_H
