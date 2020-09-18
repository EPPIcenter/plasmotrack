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


namespace transmission_nets::impl {

    /// State implementation of a simple model
    struct ModelOneState {
        static constexpr int MAX_ALLELES = 32;
        using LocusImpl = core::containers::Locus;
        using GeneticsImpl = core::datatypes::AllelesBitSet<MAX_ALLELES>;
        using InfectionEvent = core::containers::Infection<GeneticsImpl, LocusImpl>;
        using AlleleFrequencyImpl = core::datatypes::Simplex;
        using AlleleFrequencyContainerImpl = core::containers::AlleleFrequencyContainer<AlleleFrequencyImpl, LocusImpl>;


        std::map<std::string, LocusImpl *> loci{};
        std::vector<InfectionEvent *> infections{};

        AlleleFrequencyContainerImpl alleleFrequencies;

        // Network Structure
        core::parameters::Ordering<InfectionEvent> infectionEventOrdering;

        // Observation Process
        core::parameters::Parameter<double> observationFalsePositiveRate;
        core::parameters::Parameter<double> observationFalseNegativeRate;

        // Node Transmission Process
        core::parameters::Parameter<double> geometricGenerationProb;
        core::parameters::Parameter<double> ztMultiplicativeBinomialProb;
        core::parameters::Parameter<double> ztMultiplicativeBinomialAssoc;

        // Source Transmission Process
        core::parameters::Parameter<double> geometricCOIProb;

    };


}


#endif //TRANSMISSION_NETWORKS_APP_MODELONESTATE_H
