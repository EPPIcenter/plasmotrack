//
// Created by Maxwell Murphy on 4/19/20.
//

#ifndef TRANSMISSION_NETWORKS_APP_MODELONESTATE_H
#define TRANSMISSION_NETWORKS_APP_MODELONESTATE_H

#include <vector>

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
    using GeneticsImpl = AllelesBitSet<MAX_ALLELES>;
    using InfectionEvent = Infection<GeneticsImpl, Locus>;
    using AlleleFrequencyImpl = Simplex;
    using AlleleFrequencyContainer = AlleleFrequencyContainer<AlleleFrequencyImpl, Locus>;


    std::vector<Locus *> loci{};
    std::vector<InfectionEvent *> infections{};

    AlleleFrequencyContainer alleleFrequencies;

    Ordering<InfectionEvent> infectionEventOrdering;

    Parameter<double> observationFalsePositiveRate;
    Parameter<double> observationFalseNegativeRate;

    Parameter<double> geometricGenerationProb;
    Parameter<double> ztMultiplicativeBinomialProb;
    Parameter<double> ztMultiplicativeBinomialAssoc;

    Parameter<double> geometricCOIProb;

};


#endif //TRANSMISSION_NETWORKS_APP_MODELONESTATE_H
